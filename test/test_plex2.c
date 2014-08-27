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
#include <p4est_bits.h>
#else
#include <p8est_plex.h>
#include <p8est_extended.h>
#include <p8est_bits.h>
#endif

#ifdef P4EST_WITH_PETSC
#include <petsc.h>
#include <petscdmplex.h>
#include <petscsf.h>

const char          help[] = "test creating DMPlex from " P4EST_STRING "\n";

static void
locidx_to_PetscInt (sc_array_t * array)
{
  sc_array_t         *newarray;
  size_t              zz, count = array->elem_count;
  P4EST_ASSERT (array->elem_size == sizeof (p4est_locidx_t));

  if (sizeof (p4est_locidx_t) == sizeof (PetscInt)) {
    return;
  }

  newarray = sc_array_new_size (sizeof (PetscInt), array->elem_count);
  for (zz = 0; zz < count; zz++) {
    p4est_locidx_t      il = *((p4est_locidx_t *) sc_array_index (array, zz));
    PetscInt           *ip = (PetscInt *) sc_array_index (newarray, zz);

    *ip = (PetscInt) il;
  }

  sc_array_reset (array);
  sc_array_init_size (array, sizeof (PetscInt), count);
  sc_array_copy (array, newarray);
  sc_array_destroy (newarray);
}

static void
coords_double_to_PetscScalar (sc_array_t * array)
{
  sc_array_t         *newarray;
  size_t              zz, count = array->elem_count;
  P4EST_ASSERT (array->elem_size == 3 * sizeof (double));

  if (sizeof (double) == sizeof (PetscScalar)) {
    return;
  }

  newarray = sc_array_new_size (3 * sizeof (PetscScalar), array->elem_count);
  for (zz = 0; zz < count; zz++) {
    double             *id = (double *) sc_array_index (array, zz);
    PetscScalar        *ip = (PetscScalar *) sc_array_index (newarray, zz);

    ip[0] = (PetscScalar) id[0];
    ip[1] = (PetscScalar) id[1];
    ip[2] = (PetscScalar) id[2];
  }

  sc_array_reset (array);
  sc_array_init_size (array, 3 * sizeof (PetscScalar), count);
  sc_array_copy (array, newarray);
  sc_array_destroy (newarray);
}

static void
locidx_pair_to_PetscSFNode (sc_array_t * array)
{
  sc_array_t         *newarray;
  size_t              zz, count = array->elem_count;
  P4EST_ASSERT (array->elem_size == 2 * sizeof (p4est_locidx_t));

  newarray = sc_array_new_size (sizeof (PetscSFNode), array->elem_count);
  for (zz = 0; zz < count; zz++) {
    p4est_locidx_t     *il = (p4est_locidx_t *) sc_array_index (array, zz);
    PetscSFNode        *ip = (PetscSFNode *) sc_array_index (newarray, zz);

    ip->rank = (PetscInt) il[0];
    ip->index = (PetscInt) il[1];
  }

  sc_array_reset (array);
  sc_array_init_size (array, sizeof (PetscSFNode), count);
  sc_array_copy (array, newarray);
  sc_array_destroy (newarray);
}
#endif

#ifndef P4_TO_P8
static int          refine_level = 5;
#else
static int          refine_level = 3;
#endif

static int
refine_fn (p4est_t * p4est, p4est_topidx_t which_tree,
           p4est_quadrant_t * quadrant)
{
  int                 cid;

  if (which_tree == 2 || which_tree == 3) {
    return 0;
  }

  cid = p4est_quadrant_child_id (quadrant);

  if (cid == P4EST_CHILDREN - 1 ||
      (quadrant->x >= P4EST_LAST_OFFSET (P4EST_MAXLEVEL - 2) &&
       quadrant->y >= P4EST_LAST_OFFSET (P4EST_MAXLEVEL - 2)
#ifdef P4_TO_P8
       && quadrant->z >= P4EST_LAST_OFFSET (P4EST_MAXLEVEL - 2)
#endif
      )) {
    return 1;
  }
  if ((int) quadrant->level >= (refine_level - (int) (which_tree % 3))) {
    return 0;
  }
  if (quadrant->level == 1 && cid == 2) {
    return 1;
  }
  if (quadrant->x == P4EST_QUADRANT_LEN (2) &&
      quadrant->y == P4EST_LAST_OFFSET (2)) {
    return 1;
  }
  if (quadrant->y >= P4EST_QUADRANT_LEN (2)) {
    return 0;
  }

  return 1;
}

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
    *children, *parents, *childids, *leaves, *remotes;
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
  p4est = p4est_new_ext (mpicomm, conn, 0, 1, 1, 0, NULL, NULL);
  p4est_refine (p4est, 1, refine_fn, NULL);
  p4est_balance (p4est, P4EST_CONNECT_FULL, NULL);
  p4est_partition (p4est, 0, NULL);

  points_per_dim = sc_array_new (sizeof (p4est_locidx_t));
  cone_sizes = sc_array_new (sizeof (p4est_locidx_t));
  cones = sc_array_new (sizeof (p4est_locidx_t));
  cone_orientations = sc_array_new (sizeof (p4est_locidx_t));
  coords = sc_array_new (3 * sizeof (double));
  children = sc_array_new (sizeof (p4est_locidx_t));
  parents = sc_array_new (sizeof (p4est_locidx_t));
  childids = sc_array_new (sizeof (p4est_locidx_t));
  leaves = sc_array_new (sizeof (p4est_locidx_t));
  remotes = sc_array_new (2 * sizeof (p4est_locidx_t));

  p4est_get_plex_data (p4est, P4EST_CONNECT_FULL, (mpisize > 1) ? 2 : 0,
                       &first_local_quad, points_per_dim, cone_sizes, cones,
                       cone_orientations, coords, children, parents, childids,
                       leaves, remotes);

#ifdef P4EST_WITH_PETSC
  {
    PetscErrorCode      ierr;
    DM                  plex, refTree;
    PetscInt            pStart, pEnd;
    PetscSection        parentSection;
    PetscSF             pointSF;
    size_t              zz, count;

    locidx_to_PetscInt (points_per_dim);
    locidx_to_PetscInt (cone_sizes);
    locidx_to_PetscInt (cones);
    locidx_to_PetscInt (cone_orientations);
    coords_double_to_PetscScalar (coords);
    locidx_to_PetscInt (children);
    locidx_to_PetscInt (parents);
    locidx_to_PetscInt (childids);
    locidx_to_PetscInt (leaves);
    locidx_pair_to_PetscSFNode (remotes);

    P4EST_GLOBAL_PRODUCTION ("Begin PETSc routines\n");
    ierr = PetscInitialize (&argc, &argv, 0, help);
    CHKERRQ (ierr);

    ierr = DMPlexCreate (mpicomm, &plex);
    CHKERRQ (ierr);
    ierr = DMSetDimension (plex, P4EST_DIM);
    CHKERRQ (ierr);
    ierr = DMSetCoordinateDim (plex, 3);
    CHKERRQ (ierr);
    ierr = DMPlexCreateFromDAG (plex, P4EST_DIM,
                                (PetscInt *) points_per_dim->array,
                                (PetscInt *) cone_sizes->array,
                                (PetscInt *) cones->array,
                                (PetscInt *) cone_orientations->array,
                                (PetscScalar *) coords->array);
    CHKERRQ (ierr);
    ierr = PetscSFCreate (mpicomm, &pointSF);
    CHKERRQ (ierr);
    ierr =
      DMPlexCreateDefaultReferenceTree (mpicomm, P4EST_DIM, PETSC_FALSE,
                                        &refTree);
    CHKERRQ (ierr);
    ierr = DMPlexSetReferenceTree (plex, refTree);
    CHKERRQ (ierr);
    ierr = DMDestroy (&refTree);
    CHKERRQ (ierr);
    ierr = PetscSectionCreate (mpicomm, &parentSection);
    CHKERRQ (ierr);
    ierr = DMPlexGetChart (plex, &pStart, &pEnd);
    CHKERRQ (ierr);
    ierr = PetscSectionSetChart (parentSection, pStart, pEnd);
    CHKERRQ (ierr);
    count = children->elem_count;
    for (zz = 0; zz < count; zz++) {
      PetscInt            child =
        *((PetscInt *) sc_array_index (children, zz));

      ierr = PetscSectionSetDof (parentSection, child, 1);
      CHKERRQ (ierr);
    }
    ierr = PetscSectionSetUp (parentSection);
    CHKERRQ (ierr);
    ierr =
      DMPlexSetTree (plex, parentSection, (PetscInt *) parents->array,
                     (PetscInt *) childids->array);
    CHKERRQ (ierr);
    ierr = PetscSectionDestroy (&parentSection);
    CHKERRQ (ierr);
    ierr =
      PetscSFSetGraph (pointSF, pEnd - pStart, (PetscInt) leaves->elem_count,
                       (PetscInt *) leaves->array, PETSC_COPY_VALUES,
                       (PetscSFNode *) remotes->array, PETSC_COPY_VALUES);
    CHKERRQ (ierr);
    ierr = DMViewFromOptions (plex, NULL, "-dm_view");
    CHKERRQ (ierr);
    /* TODO: test with rigid body modes as in plex ex3 */
    ierr = DMDestroy (&plex);
    CHKERRQ (ierr);

    ierr = PetscFinalize ();
    P4EST_GLOBAL_PRODUCTION ("End   PETSc routines\n");
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
