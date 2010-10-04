/*
  This file is part of p4est.
  p4est is a C library to manage a collection (a forest) of multiple
  connected adaptive quadtrees or octrees in parallel.

  Copyright (C) 2010 The University of Texas System
  Written by Carsten Burstedde, Lucas C. Wilcox, and Tobin Isaac

  p4est is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  p4est is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with p4est; if not, write to the Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/*
 * This example creates uniformly refined forests in both 2D and 3D.
 * We'll figure out elsewhere how to interface this to PETSc/Sieve.
 */

#include <p4est_bits.h>
#include <p4est_extended.h>
#include <p4est_ghost.h>
#include <p4est_vtk.h>

/*
 * This program triggers an assertion in p4est_ghost_new.
 * mpirun -np 10 p4est_ptest2
 */

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
  return ((int) q->level < 2 && p4est_quadrant_child_id (q) == 0);
}

int
main (int argc, char **argv)
{
  int                 mpiret;
  mpi_context_t       mpi_context, *mpi = &mpi_context;
  p4est_t            *p4est;
  p4est_connectivity_t *connectivity;
  p4est_ghost_t      *ghost;

  /* initialize MPI and p4est internals */
  mpiret = MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  mpi->mpicomm = MPI_COMM_WORLD;        /* your favourite comm here */
  mpiret = MPI_Comm_size (mpi->mpicomm, &mpi->mpisize);
  SC_CHECK_MPI (mpiret);
  mpiret = MPI_Comm_rank (mpi->mpicomm, &mpi->mpirank);
  SC_CHECK_MPI (mpiret);

  /* this should alwaps be MPI_COMM_WORLD (no effect on p4est) */
  sc_init (MPI_COMM_WORLD, 0, 0, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  /* create 2D connectivity and forest structures */
  connectivity = p4est_connectivity_new_unitsquare ();
  p4est = p4est_new_ext (mpi->mpicomm, connectivity, 0, 0, 1, 0, NULL, NULL);

  /* refine and partition */
  p4est_refine (p4est, 1, refine_fn, NULL);
  p4est_partition (p4est, NULL);

  /* write vtk output */
  p4est_vtk_write_file (p4est, NULL, "p4est_ptest2");

  /* create and destroy ghost layer */
  ghost = p4est_ghost_new (p4est, P4EST_BALANCE_FULL);
  p4est_ghost_destroy (ghost);

  /* destroy the p4est and its connectivity structure */
  p4est_destroy (p4est);
  p4est_connectivity_destroy (connectivity);

  /* clean up and exit */
  sc_finalize ();

  mpiret = MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
