/*
  This file is part of p4est.
  p4est is a C library to manage a parallel collection of quadtrees and/or
  octrees.

  Copyright (C) 2007 Carsten Burstedde, Lucas Wilcox.

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
#endif /* HAVE_MPI */

int
main (int argc, char **argv)
{
#ifdef HAVE_MPI
  int                 mpiret;
  int32_t             count_refined, count_balanced;
  p4est_t            *p4est;
  p4est_connectivity_t *connectivity;
  double              start, elapsed_refine;
  double              elapsed_balance, elapsed_rebalance;
  MPI_Comm            mpicomm;

  /* initialize MPI */
  mpiret = MPI_Init (&argc, &argv);
  P4EST_CHECK_MPI (mpiret);
  mpicomm = MPI_COMM_WORLD;

  /* get command line argument: maximum refinement level */
  P4EST_CHECK_ABORT (argc == 2, "Give level");
  refine_level = atoi (argv[1]);
  level_shift = 4;

  /* create connectivity and forest structures */
  connectivity = p4est_connectivity_new_unitsquare ();
  p4est = p4est_new (mpicomm, NULL, connectivity, 0, NULL);

  /* time refine */
  start = -MPI_Wtime ();
  p4est_refine (p4est, refine_fractal, NULL);
  elapsed_refine = start + MPI_Wtime ();
  if (refine_level <= 6) {
    p4est_vtk_write_file (p4est, "mesh_timings_refined");
  }
  count_refined = p4est->global_num_quadrants;

  /* time balance */
  start = -MPI_Wtime ();
  p4est_balance (p4est, NULL);
  elapsed_balance = start + MPI_Wtime ();
  if (refine_level <= 6) {
    p4est_vtk_write_file (p4est, "mesh_timings_balanced");
  }
  count_balanced = p4est->global_num_quadrants;

  /* time rebalance - is a noop on the tree */
  start = -MPI_Wtime ();
  p4est_balance (p4est, NULL);
  elapsed_rebalance = start + MPI_Wtime ();
  P4EST_ASSERT (count_balanced == p4est->global_num_quadrants);

  /* print timings */
  if (p4est->mpirank == 0) {
    printf ("Level %d refined to %lld balanced to %lld\n", refine_level,
            (long long) count_refined, (long long) count_balanced);
    printf ("Level %d refinement %.3gs balance %.3gs rebalance %.3gs\n",
            refine_level, elapsed_refine, elapsed_balance, elapsed_rebalance);
  }

  /* destroy the p4est and its connectivity structure */
  p4est_destroy (p4est);
  p4est_connectivity_destroy (connectivity);

  /* clean up and exit */
  p4est_memory_check ();

  mpiret = MPI_Finalize ();
  P4EST_CHECK_MPI (mpiret);
#else
  P4EST_CHECK_ABORT (0, "This example requires the --enable-mpi flag to run.");
#endif

  return 0;
}

/* EOF timings.c */
