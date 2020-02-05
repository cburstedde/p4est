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

#include <inttypes.h>
#include <unistd.h>
#ifndef P4_TO_P8
#include <p4est_extended.h>
#include <p4est_ghost.h>
#include <p4est_mesh.h>
#else /* !P4_TO_P8 */
#include <p8est_extended.h>
#include <p8est_ghost.h>
#include <p8est_mesh.h>
#endif /* !P4_TO_P8 */

/** Set constants needed to decode an encoding obtained from \ref
 * p4est_mesh_get_neighbors
 * CAUTION: Encoding needs to be normalized to 0 (done by (abs(enc) - 1)).
 * \param [in]      direction             Encoded direction for neighbor search
 *                                        analog to direction used in
 *                                        \ref p4est_mesh_get_neighbors
 * \param [in]      btype                 Connectivity type (used for asserting
 *                                        that direction is valid).
 * \param     [out] l_same_size           Minimum value of an encoding
 *                                        indicating a same-sized neighbor.
 * \param     [out] u_same_size           Maximum value of an encoding
 *                                        indicating a same-sized neighbor.
 * \param     [out] l_double_size         Minimum value of an encoding
 *                                        indicating a double-sized neighbor.
 * \param     [out] u_double_size         Maximum value of an encoding
 *                                        indicating a double-sized neighbor.
 * \param     [out] l_half_size           Minimum value of an encoding
 *                                        indicating a half-sized neighbor.
 * \param     [out] u_half_size           Maximum value of an encoding
 *                                        indicating a half-sized neighbor.
 * \param     [out] n_neighbor_entities   Number of different entities of the
 *                                        respective neighbor entity, i.e.
 *                                        number of faces, (edges), or corners.
 * \param     [out] n_hanging_quads       Number of hanging quads across face
 *                                        (or edge).
 */
static int
set_limits (int direction, p4est_connect_type_t btype,
            int *l_same_size, int *u_same_size, int *l_double_size,
            int *u_double_size, int *l_half_size, int *u_half_size,
            int *n_neighbor_entities, int *n_hanging_quads)
{
#ifndef P4_TO_P8
  if (0 <= direction && direction < P4EST_FACES) {
    *l_same_size = 0;
    *u_same_size = *l_double_size = 8;
    *u_double_size = *l_half_size = 24;
    *u_half_size = 32;
    *n_neighbor_entities = P4EST_FACES;
    *n_hanging_quads = P4EST_HALF;
  }
  else if (P4EST_FACES <= direction &&
           direction < (P4EST_FACES + P4EST_CHILDREN)) {
    P4EST_ASSERT (btype == P4EST_CONNECT_CORNER);
    *l_same_size = 0;
    *u_same_size = 4;
    *l_double_size = *u_double_size = *l_half_size = *u_half_size = 4;
    *n_neighbor_entities = P4EST_CHILDREN;
    *n_hanging_quads = 1;
  }
  else {
    SC_ABORT_NOT_REACHED ();
  }
#else /* !P4_TO_P8 */
  if (0 <= direction && direction < P4EST_FACES) {
    *l_same_size = 0;
    *u_same_size = *l_double_size = 24;
    *u_double_size = *l_half_size = 120;
    *u_half_size = 144;
    *n_neighbor_entities = P4EST_FACES;
    *n_hanging_quads = P4EST_HALF;
  }
  else if (P4EST_FACES <= direction &&
           direction < (P4EST_FACES + P8EST_EDGES)) {
    P4EST_ASSERT (btype >= P8EST_CONNECT_EDGE);
    *l_same_size = 0;
    *u_same_size = *l_double_size = 24;
    *u_double_size = *l_half_size = 72;
    *u_half_size = 96;
    *n_neighbor_entities = P8EST_EDGES;
    *n_hanging_quads = 2;
  }
  else if ((P4EST_FACES + P8EST_EDGES) <= direction &&
           direction < (P4EST_FACES + P8EST_EDGES + P4EST_CHILDREN)) {
    P4EST_ASSERT (btype == P4EST_CONNECT_CORNER);
    *l_same_size = 0;
    *u_same_size = 8;
    *l_double_size = *u_double_size = *l_half_size = *u_half_size = 8;
    *n_neighbor_entities = P4EST_CHILDREN;
    *n_hanging_quads = 1;
  }
  else {
    SC_ABORT_NOT_REACHED ();
  }
#endif /* !P4_TO_P8 */

  return 0;
}

static int
encode_direction (int dir, int entity)
{
  switch (entity) {
  case 0:
    return dir;
  case 1:
    return dir + P4EST_FACES;
#ifdef P4_TO_P8
  case 2:
    return dir + P4EST_FACES + P8EST_EDGES;
#endif /* P4_TO_P8 */
  default:
    SC_ABORT_NOT_REACHED ();
  }
}

/** Decode encoding obtained in neighbor search
 * \param[in]      enc           The normalized encoding, i.e. 0 based and no
 *                               longer containing ghost status, i.e. 0 <= enc
 * \param[in]      n_entities    Number of faces, edges, or corners, depending
 *                               on the direction the neighbor has been looked
 *                               up.
 * \param[in]      l_same_size   Lower bound for encoding a neighbor of same
 *                               size.
 * \param[in]      u_same_size   Upper bound for encoding a neighbor of same
 *                               size.
 * \param[in]      l_double_size Lower bound for encoding a neighbor of double
 *                               size.
 * \param[in]      u_double_size Upper bound for encoding a neighbor of double
 *                               size.
 * \param[in]      l_half_size   Lower bound for encoding a neighbor of half
 *                               size.
 * \param[in]      u_half_size   Upper bound for encoding a neighbor of half
 *                               size.
 */
static int
decode_encoding (int enc, int n_entities, int l_same_size,
                 int u_same_size, int l_double_size,
                 int u_double_size, int l_half_size, int u_half_size,
                 int *subquad, int *orientation, int *entity)
{
  P4EST_ASSERT (u_same_size == l_double_size && u_double_size == l_half_size);
  P4EST_ASSERT (l_same_size == 0);
  P4EST_ASSERT (l_same_size < l_double_size);
  P4EST_ASSERT (l_double_size <= l_half_size);
  P4EST_ASSERT (0 <= enc);

#ifdef P4EST_ENABLE_DEBUG
  int8_t              upper_bnd;
  if (n_entities == P4EST_FACES) {
    upper_bnd = P4EST_HALF;
  }
#ifdef P4_TO_P8
  else if (n_entities == P8EST_EDGES) {
    upper_bnd = 2;
  }
#endif /* P4_TO_P8 */
  else if (n_entities == P4EST_CHILDREN) {
    upper_bnd = 1;
  }
  else {
    SC_ABORT_NOT_REACHED ();
  }
#endif /* P4EST_ENABLE_DEBUG */

  if (l_same_size <= enc && enc < u_same_size) {
    *orientation = enc / n_entities;
    *entity = enc % n_entities;
    *subquad = -1;
  }
  else if (l_double_size <= enc && enc < u_double_size) {
    int                 e = enc;
    e -= l_double_size;
    *subquad = e / l_double_size;
    e -= (l_double_size * *subquad);
    *orientation = e / n_entities;
    *entity = e % n_entities;

#ifdef P4EST_ENABLE_DEBUG
    P4EST_ASSERT (0 <= *subquad && *subquad < upper_bnd);
#endif /* P4EST_ENABLE_DEBUG */
  }
  else if (l_half_size <= enc && enc < u_half_size) {
    int                 e = enc;
    e -= l_half_size;
    *orientation = e / n_entities;
    *entity = e % n_entities;
    *subquad = -1;
  }
  else {
    SC_ABORT_NOT_REACHED ();
  }

#ifdef P4EST_ENABLE_DEBUG
  P4EST_ASSERT (0 <= *orientation && *orientation < upper_bnd);
#endif /* P4EST_ENABLE_DEBUG */
  P4EST_ASSERT (0 <= *entity && *entity < n_entities);

  return 0;
}

void
check_bijectivity (p4est_t * p4est, p4est_ghost_t * ghost,
                   p4est_mesh_t * mesh)
{
  int                 quad, norm_quad;
  int                 i, imax, iinv, j, k;
  int                 dir;
  int                 n_neighbor_entities, n_hanging_quads;
  int                 l_same_size, u_same_size;
  int                 l_double_size, u_double_size;
  int                 l_half_size, u_half_size;
  sc_array_t         *neighboring_quads;
  sc_array_t         *neighboring_encs;
  sc_array_t         *neighboring_qids;
  int                 neighbor_qid, neighbor_enc, neighbor_sub_ctr;
  int                 neighbor_subquad, neighbor_entity, neighbor_orientation;
  sc_array_t         *found_quads;
  sc_array_t         *found_encs;
  sc_array_t         *found_qids;
  int                 found_qid, found_enc, found_sub_ctr;
  int                 found_subquad, found_entity, found_orientation;
  int                 entity;
  int                 success;

  switch (ghost->btype) {
  case P4EST_CONNECT_FACE:
    imax = P4EST_FACES;
    break;
#ifdef P4_TO_P8
  case P8EST_CONNECT_EDGE:
    imax = P4EST_FACES + P8EST_EDGES;
    break;
#endif /* P4_TO_P8 */
  case P4EST_CONNECT_FULL:
#ifdef P4_TO_P8
    imax = P4EST_FACES + P8EST_EDGES + P4EST_CHILDREN;
#else /* P4_TO_P8 */
    imax = P4EST_FACES + P4EST_CHILDREN;
#endif /* P4_TO_P8 */
    break;
  default:
    SC_ABORT_NOT_REACHED ();
  }

  for (quad = 0; quad < p4est->global_num_quadrants; ++quad) {
    sc_MPI_Barrier (p4est->mpicomm);
    /** loop over all quads, verify only on the processor owning the respective
        quadrant. */
    if (p4est->global_first_quadrant[p4est->mpirank] <= quad &&
        quad < p4est->global_first_quadrant[p4est->mpirank + 1]) {
      /** norm global quad index to local index */
      norm_quad = quad - p4est->global_first_quadrant[p4est->mpirank];

      for (i = 0; i < imax; ++i) {
        /** allocate containers */
        neighboring_quads = sc_array_new (sizeof (p4est_quadrant_t *));
        neighboring_encs = sc_array_new (sizeof (int));
        neighboring_qids = sc_array_new (sizeof (int));

        /** set constants for decoding */
        set_limits (i, ghost->btype, &l_same_size, &u_same_size,
                    &l_double_size, &u_double_size, &l_half_size,
                    &u_half_size, &n_neighbor_entities, &n_hanging_quads);
        if (0 <= i && i < P4EST_FACES) {
          dir = i;
          entity = 0;
        }
#ifdef P4_TO_P8
        else if (P4EST_FACES <= i && i < (P4EST_FACES + P8EST_EDGES)) {
          dir = i - P4EST_FACES;
          entity = 1;
        }
        else if ((P4EST_FACES + P8EST_EDGES) <= i &&
                 i < (P4EST_FACES + P8EST_EDGES + P4EST_CHILDREN)) {
          dir = i - (P4EST_FACES + P8EST_EDGES);
          entity = 2;
        }
#else /* P4_TO_P8 */
        else if (P4EST_FACES <= i && i < (P4EST_FACES + P4EST_CHILDREN)) {
          dir = i - P4EST_FACES;
          entity = 1;
        }
#endif /* P4_TO_P8 */
        else {
          SC_ABORT_NOT_REACHED ();
        }

        success = 0;
        neighbor_sub_ctr = 0;
        found_sub_ctr = 0;

        /** search neighbors */
        p4est_mesh_get_neighbors (p4est, ghost, mesh, norm_quad, i,
                                  neighboring_quads, neighboring_encs,
                                  neighboring_qids);
        P4EST_ASSERT (neighboring_encs->elem_count ==
                      neighboring_quads->elem_count);
        P4EST_ASSERT (neighboring_encs->elem_count ==
                      neighboring_qids->elem_count);
        if (0 == neighboring_quads->elem_count) {
          success = 1;
        }

        /** inspect obtained neighbors */
        for (j = 0; j < (int) neighboring_quads->elem_count; ++j) {
          neighbor_qid = *(int *) sc_array_index (neighboring_qids, j);
          neighbor_enc = *(int *) sc_array_index (neighboring_encs, j);

          /** we can only search neighbors of local quadrants */
          if (neighbor_enc < 0) {
            if (l_half_size <= -1 - neighbor_enc &&
                -1 - neighbor_enc < u_half_size) {
              neighbor_sub_ctr = (neighbor_sub_ctr + 1) % n_hanging_quads;
            }
            continue;
          }

          /** normalize encoding before decoding (abs not necessary, because
              quadrant must be local) */
          --neighbor_enc;
          decode_encoding (neighbor_enc, n_neighbor_entities, l_same_size,
                           u_same_size, l_double_size, u_double_size,
                           l_half_size, u_half_size, &neighbor_subquad,
                           &neighbor_orientation, &neighbor_entity);
          iinv = encode_direction (neighbor_entity, entity);

          if (l_same_size <= neighbor_enc && neighbor_enc < u_same_size) {
            P4EST_ASSERT (0 == neighbor_sub_ctr);
            P4EST_ASSERT (0 == found_sub_ctr);
            /** allocate containers */
            found_quads = sc_array_new (sizeof (p4est_quadrant_t *));
            found_encs = sc_array_new (sizeof (int));
            found_qids = sc_array_new (sizeof (int));

            /** search neighbors */
            p4est_mesh_get_neighbors (p4est, ghost, mesh, neighbor_qid, iinv,
                                      found_quads, found_encs, found_qids);
            P4EST_ASSERT (found_encs->elem_count == found_quads->elem_count);
            P4EST_ASSERT (found_encs->elem_count == found_qids->elem_count);
            P4EST_ASSERT (0 < found_quads->elem_count);

            for (k = 0; k < (int) found_quads->elem_count; ++k) {
              found_qid = *(int *) sc_array_index (found_qids, k);
              found_enc = *(int *) sc_array_index (found_encs, k);
              if (found_qid == norm_quad && 0 < found_enc) {
                /** normalize encoding before decoding (abs not necessary,
                    because quadrant must be local by design) */
                --found_enc;
                decode_encoding (found_enc, n_neighbor_entities, l_same_size,
                                 u_same_size, l_double_size, u_double_size,
                                 l_half_size, u_half_size, &found_subquad,
                                 &found_orientation, &found_entity);
                /** the original quadrant must have the same size */
                P4EST_ASSERT (l_same_size <= found_enc &&
                              found_enc < u_same_size);
                P4EST_ASSERT (found_subquad == neighbor_subquad);
                SC_CHECK_ABORT (found_entity == dir, "Direction test");
                P4EST_ASSERT (found_orientation == neighbor_orientation);

                success = 1;
              }
            }

            /** de-allocate containers */
            sc_array_destroy (found_quads);
            sc_array_destroy (found_encs);
            sc_array_destroy (found_qids);
          }
          else if (l_double_size <= neighbor_enc &&
                   neighbor_enc < u_double_size) {
            P4EST_ASSERT (0 == neighbor_sub_ctr);
            P4EST_ASSERT (0 <= found_sub_ctr
                          && found_sub_ctr < n_hanging_quads);

            /** allocate containers */
            found_quads = sc_array_new (sizeof (p4est_quadrant_t *));
            found_encs = sc_array_new (sizeof (int));
            found_qids = sc_array_new (sizeof (int));

            /** search neighbors */
            p4est_mesh_get_neighbors (p4est, ghost, mesh, neighbor_qid, iinv,
                                      found_quads, found_encs, found_qids);
            P4EST_ASSERT (found_encs->elem_count == found_quads->elem_count);
            P4EST_ASSERT (found_encs->elem_count == found_qids->elem_count);
            P4EST_ASSERT (0 < found_quads->elem_count);

            for (k = 0; k < (int) found_quads->elem_count; ++k) {
              /** inspect quad id and encoding beforehand: */
              found_qid = *(int *) sc_array_index (found_qids, k);
              found_enc = *(int *) sc_array_index (found_encs, k);

              if (found_qid == norm_quad && 0 < found_enc) {
                /** normalize encoding before decoding (abs not necessary,
                    because quadrant must be local by design) */
                --found_enc;
                decode_encoding (found_enc, n_neighbor_entities, l_same_size,
                                 u_same_size, l_double_size, u_double_size,
                                 l_half_size, u_half_size, &found_subquad,
                                 &found_orientation, &found_entity);
                /** the original quadrant must be half size */
                P4EST_ASSERT (l_half_size <= found_enc &&
                              found_enc < u_half_size);
                P4EST_ASSERT (found_sub_ctr == neighbor_subquad);
                SC_CHECK_ABORT (found_entity == dir, "Direction test");
                P4EST_ASSERT (found_orientation == neighbor_orientation);
                found_sub_ctr = (found_sub_ctr + 1) % n_hanging_quads;

                success = 1;
              }
              else {
                /** hanging quads are stored contiguous in memory, such that we
                 * can cope with multiple adjacent trees */
                int                 temp_enc =
                  *(int *) sc_array_index (found_encs, k);
                temp_enc = abs (temp_enc) - 1;

                if (l_half_size <= temp_enc && temp_enc < u_half_size) {
                  found_sub_ctr = (found_sub_ctr + 1) % n_hanging_quads;
                }
                else {
                  P4EST_ASSERT (0 == found_sub_ctr);
                }
              }
            }

            /** de-allocate containers */
            sc_array_destroy (found_quads);
            sc_array_destroy (found_encs);
            sc_array_destroy (found_qids);
          }
          else if (l_half_size <= neighbor_enc && neighbor_enc < u_half_size) {
            P4EST_ASSERT (0 <= neighbor_sub_ctr &&
                          neighbor_sub_ctr < n_hanging_quads);
            P4EST_ASSERT (0 == found_sub_ctr);

            /** allocate containers */
            found_quads = sc_array_new (sizeof (p4est_quadrant_t *));
            found_encs = sc_array_new (sizeof (int));
            found_qids = sc_array_new (sizeof (int));

            /** search neighbors */
            p4est_mesh_get_neighbors (p4est, ghost, mesh, neighbor_qid, iinv,
                                      found_quads, found_encs, found_qids);
            P4EST_ASSERT (found_encs->elem_count == found_quads->elem_count);
            P4EST_ASSERT (found_encs->elem_count == found_qids->elem_count);
            P4EST_ASSERT (0 < found_quads->elem_count);

            for (k = 0; k < (int) found_quads->elem_count; ++k) {
              /** inspect quad id and encoding beforehand: */
              found_qid = *(int *) sc_array_index (found_qids, k);
              found_enc = *(int *) sc_array_index (found_encs, k);

              if (found_qid == norm_quad && 0 < found_enc) {
                /** normalize encoding before decoding (abs not necessary,
                    because quadrant must be local by design) */
                --found_enc;
                decode_encoding (found_enc, n_neighbor_entities, l_same_size,
                                 u_same_size, l_double_size, u_double_size,
                                 l_half_size, u_half_size, &found_subquad,
                                 &found_orientation, &found_entity);
                /** the original quadrant must be double size */
                P4EST_ASSERT (l_double_size <= found_enc &&
                              found_enc < u_double_size);
                P4EST_ASSERT (neighbor_sub_ctr == found_subquad);
                SC_CHECK_ABORT (found_entity == dir, "Direction test");
                P4EST_ASSERT (found_orientation == neighbor_orientation);

                success = 1;
                neighbor_sub_ctr = (neighbor_sub_ctr + 1) % n_hanging_quads;
              }
            }

            /** de-allocate containers */
            sc_array_destroy (found_quads);
            sc_array_destroy (found_encs);
            sc_array_destroy (found_qids);
          }
          else {
            SC_ABORT_NOT_REACHED ();
          }
          SC_CHECK_ABORT (1 == success, "Overall success test");
        }

        /** de-allocate containers */
        sc_array_destroy (neighboring_quads);
        sc_array_destroy (neighboring_encs);
        sc_array_destroy (neighboring_qids);
      }
    }
  }
}

/* Function for testing p4est-mesh for a single tree scenario
 *
 * \param [in] p4est     The forest.
 * \param [in] conn      The connectivity structure
 * \param [in] periodic  Flag for checking if we have periodic boundaries
 * \param [in] mpicomm   MPI communicator
 * \returns 0 for success, -1 for failure
 */
int
test_mesh_one_tree (p4est_t * p4est, p4est_connectivity_t * conn,
                    int8_t periodic, sc_MPI_Comm mpicomm)
{
  /* ensure that we have null pointers at beginning and end of function */
  P4EST_ASSERT (p4est == NULL);
  P4EST_ASSERT (conn == NULL);

  if (periodic) {
    P4EST_VERBOSE ("Check if get_neighbors is bijective for single tree,"
                   " periodic\n");
  }
  else {
    P4EST_VERBOSE ("Check if get_neighbors is bijective for single tree,"
                   " non-periodic\n");
  }

  /* create connectivity */
#ifndef P4_TO_P8
  conn = periodic == 1 ? p4est_connectivity_new_periodic ()
    : p4est_connectivity_new_unitsquare ();
#else /* !P4_TO_P8 */
  conn = periodic == 1 ? p8est_connectivity_new_periodic ()
    : p8est_connectivity_new_unitcube ();
#endif /* !P4_TO_P8 */
  /* setup p4est */
  int                 minLevel = 3;
  p4est = p4est_new_ext (mpicomm, conn, 0, minLevel, 0, 0, NULL, NULL);

  /* create mesh */
  p4est_connect_type_t btype;
#ifdef P4_TO_P8
  btype = P8EST_CONNECT_EDGE;
#else /* P4_TO_P8 */
  btype = P4EST_CONNECT_FULL;
#endif /* P4_TO_P8 */

  p4est_balance (p4est, btype, NULL);

  p4est_ghost_t      *ghost = p4est_ghost_new (p4est, btype);
  p4est_mesh_t       *mesh = p4est_mesh_new_ext (p4est, ghost, 1, 1, btype);

  /* check mesh */
  check_bijectivity (p4est, ghost, mesh);

  /* cleanup */
  p4est_ghost_destroy (ghost);
  p4est_mesh_destroy (mesh);
  p4est_destroy (p4est);
  p4est_connectivity_destroy (conn);

  conn = 0;
  p4est = 0;
  P4EST_ASSERT (p4est == NULL);
  P4EST_ASSERT (conn == NULL);

  return 0;
}

/** Function for testing p4est_mesh for all existing two_trees scenarios
 * \param [in] p4est     The forest, NULL.
 * \param [in] conn      The connectivity structure, NULL.
 * \param [in] mpicomm   MPI communicator
 */
int
test_mesh_two_trees (p4est_t * p4est, p4est_connectivity_t * conn,
                     sc_MPI_Comm mpicomm)
{
  /* ensure that we have null pointers at beginning and end of function */
  P4EST_ASSERT (p4est == NULL);
  P4EST_ASSERT (conn == NULL);
  int                 conn_face_tree1, conn_face_tree2, orientation;
  for (conn_face_tree1 = 0; conn_face_tree1 < P4EST_FACES; ++conn_face_tree1) {
    for (conn_face_tree2 = 0; conn_face_tree2 < P4EST_FACES;
         ++conn_face_tree2) {
      for (orientation = 0; orientation < P4EST_HALF; ++orientation) {
        P4EST_VERBOSEF ("Check if get_neighbors is bijective for two trees,"
                        " left face: %i, right face: %i, orientation: %i \n",
                        conn_face_tree1, conn_face_tree2, orientation);
        /* create connectivity */
        conn =
          p4est_connectivity_new_twotrees (conn_face_tree1, conn_face_tree2,
                                           orientation);
        /* setup p4est */
        int                 minLevel = 3;
        p4est = p4est_new_ext (mpicomm, conn, 0, minLevel, 0, 0, NULL, NULL);

        p4est_connect_type_t btype;
#ifdef P4_TO_P8
        btype = P8EST_CONNECT_EDGE;
#else /* P4_TO_P8 */
        btype = P4EST_CONNECT_FULL;
#endif /* P4_TO_P8 */

        p4est_balance (p4est, btype, NULL);

        /* create mesh */
        p4est_ghost_t      *ghost = p4est_ghost_new (p4est, btype);
        p4est_mesh_t       *mesh =
          p4est_mesh_new_ext (p4est, ghost, 1, 1, btype);

        /* check mesh */
        check_bijectivity (p4est, ghost, mesh);

        /* cleanup */
        p4est_ghost_destroy (ghost);
        p4est_mesh_destroy (mesh);
        p4est_destroy (p4est);
        p4est_connectivity_destroy (conn);

        conn = 0;
        p4est = 0;
        P4EST_ASSERT (p4est == NULL);
        P4EST_ASSERT (conn == NULL);
      }
    }
  }

  return 0;
}

/* Function for testing p4est-mesh for multiple trees in a brick scenario
 *
 * \param [in] p4est     The forest.
 * \param [in] conn      The connectivity structure
 * \param [in] periodic  Flag for checking if we have periodic boundaries
 * \param [in] mpicomm   MPI communicator
 * \returns 0 for success, -1 for failure
 */
int
test_mesh_multiple_trees_brick (p4est_t * p4est, p4est_connectivity_t * conn,
                                int8_t periodic, sc_MPI_Comm mpicomm)
{
  /* ensure that we have null pointers at beginning and end of function */
  P4EST_ASSERT (p4est == NULL);
  P4EST_ASSERT (conn == NULL);

  if (periodic) {
    P4EST_VERBOSE ("Check if get_neighbors is bijective for brick of trees,"
                   " periodic\n");
  }
  else {
    P4EST_VERBOSE ("Check if get_neighbors is bijective for brick of trees,"
                   " non-periodic\n");
  }

  /* create connectivity */
#ifndef P4_TO_P8
  conn = p4est_connectivity_new_brick (2, 1, periodic, periodic);
#else /* !P4_TO_P8 */
  conn = p8est_connectivity_new_brick (2, 1, 1, periodic, periodic, periodic);
#endif /* !P4_TO_P8 */
  /* setup p4est */
  int                 minLevel = 3;
  p4est = p4est_new_ext (mpicomm, conn, 0, minLevel, 0, 0, NULL, NULL);

  p4est_connect_type_t btype;
#ifdef P4_TO_P8
  btype = P8EST_CONNECT_EDGE;
#else /* P4_TO_P8 */
  btype = P4EST_CONNECT_FULL;
#endif /* P4_TO_P8 */

  p4est_balance (p4est, btype, NULL);

  /* create mesh */
  p4est_ghost_t      *ghost = p4est_ghost_new (p4est, btype);
  p4est_mesh_t       *mesh = p4est_mesh_new_ext (p4est, ghost, 1, 1, btype);

  /* check mesh */
  check_bijectivity (p4est, ghost, mesh);

  /* cleanup */
  p4est_ghost_destroy (ghost);
  p4est_mesh_destroy (mesh);
  p4est_destroy (p4est);
  p4est_connectivity_destroy (conn);

  conn = 0;
  p4est = 0;
  P4EST_ASSERT (p4est == NULL);
  P4EST_ASSERT (conn == NULL);

  return 0;
}

/* Function for testing p4est-mesh for multiple trees in a non-brick scenario
 *
 * \param [in] p4est     The forest.
 * \param [in] conn      The connectivity structure
 * \param [in] periodic  Flag for checking if we have periodic boundaries
 * \returns 0 for success, -1 for failure
 */
int
test_mesh_multiple_trees_nonbrick (p4est_t * p4est,
                                   p4est_connectivity_t * conn,
                                   int8_t periodic, sc_MPI_Comm mpicomm)
{
  return 0;
}

int
main (int argc, char **argv)
{
  sc_MPI_Comm         mpicomm;
  int                 mpiret;
  int                 mpisize, mpirank;
  p4est_t            *p4est;
  p4est_connectivity_t *conn;
  int8_t              periodic_boundaries;
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
  p4est = 0;
  conn = 0;

  int                 test_single, test_two_trees;
  int                 test_multi_brick, test_multi_non_brick;
  test_single = 1;
  test_two_trees = 1;
  test_multi_brick = 1;
  test_multi_non_brick = 0;

  /* test both periodic and non-periodic boundaries */
  if (test_single) {
    /* test one tree */
    periodic_boundaries = 0;
    test_mesh_one_tree (p4est, conn, periodic_boundaries, mpicomm);

    periodic_boundaries = 1;
    test_mesh_one_tree (p4est, conn, periodic_boundaries, mpicomm);
  }

  if (test_two_trees) {
    /* test two trees
     * (using all combinations of faces connected and all orientations) */
    test_mesh_two_trees (p4est, conn, mpicomm);
  }

  if (test_multi_brick) {
    /* test multiple trees; brick */
    periodic_boundaries = 0;
    test_mesh_multiple_trees_brick (p4est, conn, periodic_boundaries,
                                    mpicomm);
    periodic_boundaries = 1;
    test_mesh_multiple_trees_brick (p4est, conn, periodic_boundaries,
                                    mpicomm);
  }

  if (test_multi_non_brick) {
    /* test multiple trees; non-brick */
    periodic_boundaries = 0;
    test_mesh_multiple_trees_nonbrick (p4est, conn, periodic_boundaries,
                                       mpicomm);
    periodic_boundaries = 1;
    test_mesh_multiple_trees_nonbrick (p4est, conn, periodic_boundaries,
                                       mpicomm);
  }
  /* exit */
  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;
}
