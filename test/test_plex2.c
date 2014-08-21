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

#ifndef P4_TO_P8
#include <p4est_plex.h>
#include <p4est_extended.h>
#else
#include <p8est_plex.h>
#include <p8est_extended.h>
#endif

#ifdef P4EST_PETSC
#include <petsc.h>
#include <petscdmplex.h>
const char help[] = "test creating DMPlex from " P4EST_STRING "\n";
#endif

int
main (int argc, char **argv)
{
  sc_MPI_Comm         mpicomm;
  int                 mpiret;
  int                 mpisize, mpirank;
  p4est_t            *p4est;
  p4est_connectivity_t *conn;
  sc_array_t         *points_per_dim, *cone_sizes, *cones,
                     *cone_orientations, *coords,
                     *children, *parents,
                     *childids, *leaves, *remotes;
  p4est_locidx_t      first_local_quad = -1;

  /* initialize MPI */
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  mpicomm = sc_MPI_COMM_WORLD;
  mpiret = sc_MPI_Comm_size (mpicomm, &mpisize);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (mpicomm, &mpirank);
  SC_CHECK_MPI (mpiret);

  sc_init (mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

#ifndef P4_TO_P8
  conn = p4est_connectivity_new_moebius ();
#else
  conn = p8est_connectivity_new_rotcubes ();
#endif
  p4est = p4est_new_ext (mpicomm, conn, 0, 0, 1, 0, NULL, NULL);

  points_per_dim = sc_array_new (sizeof (p4est_locidx_t));
  cone_sizes = sc_array_new (sizeof (p4est_locidx_t));
  cones = sc_array_new (sizeof (p4est_locidx_t));
  cone_orientations = sc_array_new (sizeof (p4est_locidx_t));
  coords = sc_array_new (3 * sizeof (double));
  children = sc_array_new (sizeof (p4est_locidx_t));
  parents = sc_array_new (sizeof (p4est_locidx_t));
  childids = sc_array_new (sizeof (p4est_locidx_t));
  leaves = sc_array_new (sizeof (p4est_locidx_t));
  remotes = sc_array_new (sizeof (p4est_locidx_t));

  p4est_get_plex_data (p4est, P4EST_CONNECT_FULL, (mpisize > 1) ? 2 : 0,
                       &first_local_quad, points_per_dim, cone_sizes, cones,
                       cone_orientations, coords, children, parents, childids,
                       leaves, remotes);

#ifdef P4EST_PETSC
  {
    PetscErrorCode ierr;
    DM plex;

    ierr = PetscInitialize (&argc, &argv, 0, help);CHKERRQ(ierr);

    ierr = DMPlexCreate (mpicomm, &plex);CHKERRQ(ierr);
    ierr = DMSetDimension (plex, P4EST_DIM);CHKERRQ(ierr);
    ierr = DMPlexCreateFromDAG (plex, P4EST_DIM,
                                (PetscInt *) points_per_dim->array,
                                (PetscInt *) cone_sizes->array,
                                (PetscInt *) cones->array,
                                (PetscInt *) cone_orientations->array,
                                (PetscScalar *) coords->array);CHKERRQ(ierr);
    ierr = DMViewFromOptions(plex,NULL,"-dm_view");CHKERRQ(ierr);
    ierr = DMDestroy (&plex);CHKERRQ(ierr);

    ierr = PetscFinalize();
  }
#endif

  sc_array_destroy (points_per_dim);
  sc_array_destroy (cone_sizes);
  sc_array_destroy (cones);
  sc_array_destroy (cone_orientations);
  sc_array_destroy (coords);
  sc_array_destroy (children);
  sc_array_destroy (parents);
  sc_array_destroy (childids);
  sc_array_destroy (leaves);
  sc_array_destroy (remotes);

  p4est_destroy (p4est);
  p4est_connectivity_destroy (conn);

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}

/* EOF test_plex2.c */
