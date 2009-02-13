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
#include <p4est_vtk.h>

typedef struct
{
  MPI_Comm            mpicomm;
  int                 mpisize;
  int                 mpirank;
}
mpi_context_t;

static int
refine_fn (p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t * q)
{
  if ((which_tree == 0 && q->x != 0 && q->y != 0) ||
      (which_tree == 2 && q->x >= P4EST_LAST_OFFSET (2) &&
       q->y >= P4EST_QUADRANT_LEN (1) && q->y < P4EST_LAST_OFFSET (2))) {
    return (q->level < 5);
  }

  return 0;
}

int
main (int argc, char **argv)
{
  int                 mpiret;
  unsigned            crc;
  p4est_t            *p4est;
  p4est_connectivity_t *connectivity;
  mpi_context_t       mpi_context, *mpi = &mpi_context;
  const p4est_locidx_t given[5] = { 3, 7, 36, 10, 10 };

  /* initialize MPI and p4est internals */
  mpiret = MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  mpi->mpicomm = MPI_COMM_WORLD;
  mpiret = MPI_Comm_size (mpi->mpicomm, &mpi->mpisize);
  SC_CHECK_MPI (mpiret);
  mpiret = MPI_Comm_rank (mpi->mpicomm, &mpi->mpirank);
  SC_CHECK_MPI (mpiret);

  sc_init (mpi->mpicomm, true, true, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  SC_CHECK_ABORT (mpi->mpisize == 5, "This example requires MPI with np=5.");

  /* create connectivity and forest structures */
  connectivity = p4est_connectivity_new_star ();
  p4est = p4est_new (mpi->mpicomm, connectivity, 15, 0, NULL, NULL);

  /* partition and refine the mesh */
  (void) p4est_partition_given (p4est, given);
  p4est_refine (p4est, true, refine_fn, NULL);
  p4est_vtk_write_file (p4est, "mesh_second_refined");

  /* balance the mesh */
  p4est_balance (p4est, P4EST_BALANCE_FULL, NULL);
  p4est_vtk_write_file (p4est, "mesh_second_balanced");
  crc = p4est_checksum (p4est);

  /* print and verify forest checksum */
  P4EST_GLOBAL_STATISTICSF ("Tree checksum 0x%x\n", crc);
  if (mpi->mpirank == 0)
    SC_CHECK_ABORT (crc == 0x324eb631U, "Checksum mismatch");

  /* destroy the p4est and its connectivity structure */
  p4est_destroy (p4est);
  p4est_connectivity_destroy (connectivity);

  /* clean up and exit */
  sc_finalize ();

  mpiret = MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}

/* EOF second.c */
