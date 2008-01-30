/*
  This file is part of p4est.
  p4est is a C library to manage a parallel collection of quadtrees and/or
  octrees.

  Copyright (C) 2007,2008 Carsten Burstedde, Lucas Wilcox.

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/*
 * Usage: p4est_timings <level>
 */

#include <p4est_algorithms.h>
#include <p4est_base.h>
#include <p4est_vtk.h>

#ifdef HAVE_MPI

typedef struct
{
  MPI_Comm            mpicomm;
  int                 mpisize;
  int                 mpirank;
}
mpi_context_t;

static int          refine_level = 0;
static int          level_shift = 0;

static int
refine_fractal (p4est_t * p4est, int32_t which_tree, p4est_quadrant_t * q)
{
  int8_t              qid;

  if (q->level >= refine_level) {
    return 0;
  }
  if (q->level < refine_level - level_shift) {
    return 1;
  }

  qid = p4est_quadrant_child_id (q);
  return (qid == 0 || qid == 3);
}

static void
abort_fn (void *data)
{
  int                 mpiret;
  mpi_context_t      *mpi = data;

  fprintf (stderr, "[%d] p4est_timings abort handler\n", mpi->mpirank);

  mpiret = MPI_Abort (mpi->mpicomm, 1);
  P4EST_CHECK_MPI (mpiret);
}

#endif /* HAVE_MPI */

int
main (int argc, char **argv)
{
#ifdef HAVE_MPI
  int                 mpiret;
  unsigned            crc;
  int32_t             count_refined, count_balanced;
  p4est_t            *p4est;
  p4est_connectivity_t *connectivity;
  double              start, elapsed_refine;
  double              elapsed_balance, elapsed_rebalance;
  double              elapsed_partition, elapsed_repartition;
  mpi_context_t       mpi_context, *mpi = &mpi_context;

  /* initialize MPI and p4est internals */
  mpiret = MPI_Init (&argc, &argv);
  P4EST_CHECK_MPI (mpiret);
  mpi->mpicomm = MPI_COMM_WORLD;
  mpiret = MPI_Comm_size (mpi->mpicomm, &mpi->mpisize);
  P4EST_CHECK_MPI (mpiret);
  mpiret = MPI_Comm_rank (mpi->mpicomm, &mpi->mpirank);
  P4EST_CHECK_MPI (mpiret);
  p4est_init (stdout, mpi->mpirank, abort_fn, mpi);

  /* get command line argument: maximum refinement level */
  P4EST_CHECK_ABORT (argc == 2, "Give level");
  refine_level = atoi (argv[1]);
  level_shift = 4;

  /* print general setup information */
  if (mpi->mpirank == 0) {
    printf ("Processors %d level %d shift %d\n",
            mpi->mpisize, refine_level, level_shift);
  }

  /* create connectivity and forest structures */
  connectivity = p4est_connectivity_new_unitsquare ();
  p4est = p4est_new (mpi->mpicomm, connectivity, 0, NULL);

  /* time refine */
  mpiret = MPI_Barrier (mpi->mpicomm);
  P4EST_CHECK_MPI (mpiret);
  start = -MPI_Wtime ();
  p4est_refine (p4est, refine_fractal, NULL);
  mpiret = MPI_Barrier (mpi->mpicomm);
  P4EST_CHECK_MPI (mpiret);
  elapsed_refine = start + MPI_Wtime ();
  if (refine_level <= 8) {
    p4est_vtk_write_file (p4est, "mesh_timings_refined");
  }
  count_refined = p4est->global_num_quadrants;

  /* time balance */
  mpiret = MPI_Barrier (mpi->mpicomm);
  P4EST_CHECK_MPI (mpiret);
  start = -MPI_Wtime ();
  p4est_balance (p4est, NULL);
  mpiret = MPI_Barrier (mpi->mpicomm);
  P4EST_CHECK_MPI (mpiret);
  elapsed_balance = start + MPI_Wtime ();
  if (refine_level <= 8) {
    p4est_vtk_write_file (p4est, "mesh_timings_balanced");
  }
  count_balanced = p4est->global_num_quadrants;
  crc = p4est_checksum (p4est);

  /* time rebalance - is a noop on the tree */
  mpiret = MPI_Barrier (mpi->mpicomm);
  P4EST_CHECK_MPI (mpiret);
  start = -MPI_Wtime ();
  p4est_balance (p4est, NULL);
  mpiret = MPI_Barrier (mpi->mpicomm);
  P4EST_CHECK_MPI (mpiret);
  elapsed_rebalance = start + MPI_Wtime ();
  P4EST_ASSERT (count_balanced == p4est->global_num_quadrants);
  P4EST_ASSERT (crc == p4est_checksum (p4est));

  /* time partition */
  mpiret = MPI_Barrier (mpi->mpicomm);
  P4EST_CHECK_MPI (mpiret);
  start = -MPI_Wtime ();
  p4est_partition (p4est, NULL);
  mpiret = MPI_Barrier (mpi->mpicomm);
  P4EST_CHECK_MPI (mpiret);
  elapsed_partition = start + MPI_Wtime ();
  if (refine_level <= 8) {
    p4est_vtk_write_file (p4est, "mesh_timings_partitioned");
  }
  P4EST_ASSERT (crc == p4est_checksum (p4est));

  /* time repartition - is a noop on the tree */
  mpiret = MPI_Barrier (mpi->mpicomm);
  P4EST_CHECK_MPI (mpiret);
  start = -MPI_Wtime ();
  p4est_partition (p4est, NULL);
  mpiret = MPI_Barrier (mpi->mpicomm);
  P4EST_CHECK_MPI (mpiret);
  elapsed_repartition = start + MPI_Wtime ();
  P4EST_ASSERT (crc == p4est_checksum (p4est));

  /* print checksum and timings */
  P4EST_GLOBAL_STATISTICSF ("Processors %d level %d shift %d"
                            " tree checksum 0x%x\n",
                            mpi->mpisize, refine_level, level_shift, crc);
  P4EST_GLOBAL_STATISTICSF ("Level %d refined to %lld balanced to %lld\n",
                            refine_level, (long long) count_refined,
                            (long long) count_balanced);
  P4EST_GLOBAL_STATISTICSF ("Level %d refinement %.3gs"
                            " balance %.3gs rebalance %.3gs\n",
                            refine_level, elapsed_refine,
                            elapsed_balance, elapsed_rebalance);
  P4EST_GLOBAL_STATISTICSF ("Level %d partition %.3gs repartition %.3gs\n",
                            refine_level, elapsed_partition,
                            elapsed_repartition);

  /* destroy the p4est and its connectivity structure */
  p4est_destroy (p4est);
  p4est_connectivity_destroy (connectivity);

  /* clean up and exit */
  p4est_memory_check ();

  mpiret = MPI_Finalize ();
  P4EST_CHECK_MPI (mpiret);
#else
  P4EST_CHECK_ABORT (0, "This example requires --enable-mpi to run.");
#endif /* HAVE_MPI */

  return 0;
}

/* EOF timings.c */
