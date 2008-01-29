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
 * This program needs to be run with 5 processors.
 * Usage: p4est_second
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

static int
refine_fn (p4est_t * p4est, int32_t which_tree, p4est_quadrant_t * q)
{
  if ((which_tree == 0 && q->x != 0 && q->y != 0) ||
      (which_tree == 2 && q->x >= P4EST_LD (2) &&
       q->y >= P4EST_QH (1) && q->y < P4EST_LD (2))) {
    return (q->level < 5);
  }

  return 0;
}

static void
abort_fn (void *data)
{
  int                 mpiret;
  mpi_context_t      *mpi = data;

  fprintf (stderr, "[%d] p4est_second abort handler\n", mpi->mpirank);

  mpiret = MPI_Abort (mpi->mpicomm, 1);
  P4EST_CHECK_MPI (mpiret);
}

#endif /* HAVE_MPI */

int
main (int argc, char **argv)
{
#ifdef HAVE_MPI
  int                 mpiret;
  p4est_t            *p4est;
  p4est_connectivity_t *connectivity;
  mpi_context_t       mpi_context, *mpi = &mpi_context;
  const int32_t       given[5] = { 3, 7, 36, 10, 10 };

  /* initialize MPI and p4est internals */
  mpiret = MPI_Init (&argc, &argv);
  P4EST_CHECK_MPI (mpiret);
  mpi->mpicomm = MPI_COMM_WORLD;
  mpiret = MPI_Comm_size (mpi->mpicomm, &mpi->mpisize);
  P4EST_CHECK_MPI (mpiret);
  mpiret = MPI_Comm_rank (mpi->mpicomm, &mpi->mpirank);
  P4EST_CHECK_MPI (mpiret);
  p4est_init (stdout, mpi->mpirank, abort_fn, mpi);
  P4EST_CHECK_ABORT (mpi->mpisize == 5, "This example requires np=5.");

  /* create connectivity and forest structures */
  connectivity = p4est_connectivity_new_star ();
  p4est = p4est_new (mpi->mpicomm, connectivity, 0, NULL);

  /* partition and refine the mesh */
  p4est_partition_given (p4est, given);
  p4est_refine (p4est, refine_fn, NULL);
  p4est_vtk_write_file (p4est, "mesh_second_refined");

  /* balance the mesh */
  p4est_balance (p4est, NULL);
  p4est_vtk_write_file (p4est, "mesh_second_balanced");

  /* print forest checksum */
  P4EST_GLOBAL_INFOF ("Tree checksum 0x%x\n", p4est_checksum (p4est));

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

/* EOF second.c */
