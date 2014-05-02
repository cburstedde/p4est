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

#include <p4est_extended.h>
#include <p4est_vtk.h>
#include <p8est_extended.h>
#include <p8est_vtk.h>

typedef struct
{
  sc_MPI_Comm         mpicomm;
  int                 mpisize;
  int                 mpirank;
}
mpi_context_t;

int
main (int argc, char **argv)
{
  int                 mpiret;
  int                 refine_level;
  size_t              size2, size3;
  long long           lsize[2], gsize[2];
  mpi_context_t       mpi_context, *mpi = &mpi_context;
  p4est_t            *p4est;
  p4est_connectivity_t *connectivity2;
  p8est_t            *p8est;
  p8est_connectivity_t *connectivity3;

  /* initialize MPI and p4est internals */
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  mpi->mpicomm = sc_MPI_COMM_WORLD;     /* your favourite comm here */
  mpiret = sc_MPI_Comm_size (mpi->mpicomm, &mpi->mpisize);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (mpi->mpicomm, &mpi->mpirank);
  SC_CHECK_MPI (mpiret);

  /* this should alwaps be sc_MPI_COMM_WORLD (no effect on p4est) */
  sc_init (sc_MPI_COMM_WORLD, 0, 0, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  /* initial level of refinement */
  refine_level = 2;

  /* create 2D connectivity and forest structures */
  connectivity2 = p4est_connectivity_new_corner ();
  size2 = p4est_connectivity_memory_used (connectivity2);
  p4est = p4est_new_ext (mpi->mpicomm, connectivity2,
                         0, refine_level, 1, 0, NULL, NULL);
  size2 += p4est_memory_used (p4est);

  /* create 3D connectivity and forest structures */
  connectivity3 = p8est_connectivity_new_rotcubes ();
  size3 = p8est_connectivity_memory_used (connectivity3);
  p8est = p8est_new_ext (mpi->mpicomm, connectivity3,
                         0, refine_level, 1, 0, NULL, NULL);
  size3 += p8est_memory_used (p8est);

  /* compute total size of forest storage */
  lsize[0] = (long long) size2;
  lsize[1] = (long long) size3;
  mpiret =
    sc_MPI_Reduce (lsize, gsize, 2, sc_MPI_LONG_LONG_INT, sc_MPI_SUM, 0,
                   mpi->mpicomm);
  SC_CHECK_MPI (mpiret);
  P4EST_GLOBAL_INFOF ("Total forest byte sizes: %lld (2D), %lld (3D)\n",
                      gsize[0], gsize[1]);

  /* write vtk output files */
  p4est_vtk_write_file (p4est, NULL, "p4est_sieve");
  p8est_vtk_write_file (p8est, NULL, "p8est_sieve");

  /* destroy the p4est and its connectivity structure */
  p4est_destroy (p4est);
  p4est_connectivity_destroy (connectivity2);
  p8est_destroy (p8est);
  p8est_connectivity_destroy (connectivity3);

  /* clean up and exit */
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
