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

#include <p4est_bits.h>
#include <p4est_extended.h>
#include <p4est_nodes.h>

typedef struct
{
  p4est_topidx_t      a;
  int64_t             sum;
}
user_data_t;

typedef struct p4est_vert
{
  double              x, y, z;
  p4est_topidx_t      treeid;
}
p4est_vert_t;

static int          refine_level = 6;

static void
init_fn (p4est_t * p4est, p4est_topidx_t which_tree,
         p4est_quadrant_t * quadrant)
{
  user_data_t        *data = (user_data_t *) quadrant->p.user_data;

  data->a = which_tree;
  data->sum = quadrant->x + quadrant->y + quadrant->level;
}

static int
refine_fn (p4est_t * p4est, p4est_topidx_t which_tree,
           p4est_quadrant_t * quadrant)
{
  if (quadrant->level >= (refine_level - (which_tree % 3))) {
    return 0;
  }
  if (quadrant->level == 1 && p4est_quadrant_child_id (quadrant) == 3) {
    return 1;
  }
  if (quadrant->x == P4EST_LAST_OFFSET (2) &&
      quadrant->y == P4EST_LAST_OFFSET (2)) {
    return 1;
  }
  if (quadrant->x >= P4EST_QUADRANT_LEN (2)) {
    return 0;
  }

  return 1;
}

static int
weight_one (p4est_t * p4est, p4est_topidx_t which_tree,
            p4est_quadrant_t * quadrant)
{
  return 1;
}

static int
p4est_vert_compare (const void *a, const void *b)
{
  const p4est_vert_t *v1 = (const p4est_vert_t *) a;
  const p4est_vert_t *v2 = (const p4est_vert_t *) b;
  const double        eps = 1e-15;
  double              xdiff, ydiff, zdiff;
  int                 retval = 0;

  if (v1->treeid != v2->treeid) {
    return (int) v1->treeid - (int) v2->treeid;
  }

  xdiff = fabs (v1->x - v2->x);
  if (xdiff < eps) {
    ydiff = fabs (v1->y - v2->y);
    if (ydiff < eps) {
      zdiff = fabs (v1->z - v2->z);
      if (zdiff < eps) {
        retval = 0;
      }
      else {
        retval = (v1->z < v2->z) ? -1 : 1;
      }
    }
    else {
      retval = (v1->y < v2->y) ? -1 : 1;
    }
  }
  else {
    retval = (v1->x < v2->x) ? -1 : 1;
  }

  return retval;
}

static void
p4est_check_local_order (p4est_t * p4est, p4est_connectivity_t * connectivity)
{
  const double        intsize = 1.0 / P4EST_ROOT_LEN;
  double             *vertices;
  double              h, eta1, eta2;
  double              v0x, v0y, v0z, v1x, v1y, v1z;
  double              v2x, v2y, v2z, v3x, v3y, v3z;
  double              w0x, w0y, w0z, w1x, w1y, w1z;
  double              w2x, w2y, w2z, w3x, w3y, w3z;
  size_t              iz;
  size_t              num_quads;
  size_t              quad_count;
  p4est_topidx_t      jt;
  p4est_topidx_t     *tree_to_vertex;
  p4est_topidx_t      first_local_tree;
  p4est_topidx_t      last_local_tree;
  p4est_topidx_t      v0, v1, v2, v3;
  p4est_locidx_t      kl;
  p4est_locidx_t      lv0, lv1, lv2, lv3;
  p4est_locidx_t      num_uniq_local_vertices;
  p4est_locidx_t     *quadrant_to_local_vertex;
  p4est_qcoord_t      inth;
  p4est_tree_t       *tree;
  p4est_quadrant_t   *quad;
  p4est_vert_t       *vert_locations;
  p4est_nodes_t      *nodes;
  sc_array_t         *trees;
  sc_array_t         *quadrants;

  nodes = p4est_nodes_new (p4est, NULL);
  quadrant_to_local_vertex = nodes->local_nodes;
  num_uniq_local_vertices = nodes->num_owned_indeps;
  SC_CHECK_ABORT ((size_t) num_uniq_local_vertices ==
                  nodes->indep_nodes.elem_count, "Node count mismatch");

  P4EST_INFOF ("Unique local vertices %lld\n",
               (long long) num_uniq_local_vertices);

  vert_locations = P4EST_ALLOC (p4est_vert_t, num_uniq_local_vertices);
  for (kl = 0; kl < num_uniq_local_vertices; ++kl) {
    vert_locations[kl].treeid = -1;
  }

  tree_to_vertex = connectivity->tree_to_vertex;
  vertices = connectivity->vertices;
  first_local_tree = p4est->first_local_tree;
  last_local_tree = p4est->last_local_tree;
  trees = p4est->trees;

  for (jt = first_local_tree, quad_count = 0; jt <= last_local_tree; ++jt) {
    tree = p4est_tree_array_index (trees, jt);

    P4EST_ASSERT (0 <= jt && jt < connectivity->num_trees);

    v0 = tree_to_vertex[jt * 4 + 0];
    v1 = tree_to_vertex[jt * 4 + 1];
    v2 = tree_to_vertex[jt * 4 + 2];
    v3 = tree_to_vertex[jt * 4 + 3];

    P4EST_ASSERT (0 <= v0 && v0 < connectivity->num_vertices);
    P4EST_ASSERT (0 <= v1 && v1 < connectivity->num_vertices);
    P4EST_ASSERT (0 <= v2 && v2 < connectivity->num_vertices);
    P4EST_ASSERT (0 <= v3 && v3 < connectivity->num_vertices);

    v0x = vertices[v0 * 3 + 0];
    v0y = vertices[v0 * 3 + 1];
    v0z = vertices[v0 * 3 + 2];
    v1x = vertices[v1 * 3 + 0];
    v1y = vertices[v1 * 3 + 1];
    v1z = vertices[v1 * 3 + 2];
    v2x = vertices[v2 * 3 + 0];
    v2y = vertices[v2 * 3 + 1];
    v2z = vertices[v2 * 3 + 2];
    v3x = vertices[v3 * 3 + 0];
    v3y = vertices[v3 * 3 + 1];
    v3z = vertices[v3 * 3 + 2];

    quadrants = &tree->quadrants;
    num_quads = quadrants->elem_count;

    /* loop over the elements in the tree */
    for (iz = 0; iz < num_quads; ++iz, ++quad_count) {
      quad = p4est_quadrant_array_index (quadrants, iz);
      inth = P4EST_QUADRANT_LEN (quad->level);
      h = intsize * inth;
      eta1 = intsize * quad->x;
      eta2 = intsize * quad->y;

      w0x = v0x * (1.0 - eta1) * (1.0 - eta2)
        + v1x * (eta1) * (1.0 - eta2)
        + v2x * (1.0 - eta1) * (eta2)
        + v3x * (eta1) * (eta2);

      w0y = v0y * (1.0 - eta1) * (1.0 - eta2)
        + v1y * (eta1) * (1.0 - eta2)
        + v2y * (1.0 - eta1) * (eta2)
        + v3y * (eta1) * (eta2);

      w0z = v0z * (1.0 - eta1) * (1.0 - eta2)
        + v1z * (eta1) * (1.0 - eta2)
        + v2z * (1.0 - eta1) * (eta2)
        + v3z * (eta1) * (eta2);

      w1x = v0x * (1.0 - eta1 - h) * (1.0 - eta2)
        + v1x * (eta1 + h) * (1.0 - eta2)
        + v2x * (1.0 - eta1 - h) * (eta2)
        + v3x * (eta1 + h) * (eta2);

      w1y = v0y * (1.0 - eta1 - h) * (1.0 - eta2)
        + v1y * (eta1 + h) * (1.0 - eta2)
        + v2y * (1.0 - eta1 - h) * (eta2)
        + v3y * (eta1 + h) * (eta2);

      w1z = v0z * (1.0 - eta1 - h) * (1.0 - eta2)
        + v1z * (eta1 + h) * (1.0 - eta2)
        + v2z * (1.0 - eta1 - h) * (eta2)
        + v3z * (eta1 + h) * (eta2);

      w2x = v0x * (1.0 - eta1) * (1.0 - eta2 - h)
        + v1x * (eta1) * (1.0 - eta2 - h)
        + v2x * (1.0 - eta1) * (eta2 + h)
        + v3x * (eta1) * (eta2 + h);

      w2y = v0y * (1.0 - eta1) * (1.0 - eta2 - h)
        + v1y * (eta1) * (1.0 - eta2 - h)
        + v2y * (1.0 - eta1) * (eta2 + h)
        + v3y * (eta1) * (eta2 + h);

      w2z = v0z * (1.0 - eta1) * (1.0 - eta2 - h)
        + v1z * (eta1) * (1.0 - eta2 - h)
        + v2z * (1.0 - eta1) * (eta2 + h)
        + v3z * (eta1) * (eta2 + h);

      w3x = v0x * (1.0 - eta1 - h) * (1.0 - eta2 - h)
        + v1x * (eta1 + h) * (1.0 - eta2 - h)
        + v2x * (1.0 - eta1 - h) * (eta2 + h)
        + v3x * (eta1 + h) * (eta2 + h);

      w3y = v0y * (1.0 - eta1 - h) * (1.0 - eta2 - h)
        + v1y * (eta1 + h) * (1.0 - eta2 - h)
        + v2y * (1.0 - eta1 - h) * (eta2 + h)
        + v3y * (eta1 + h) * (eta2 + h);

      w3z = v0z * (1.0 - eta1 - h) * (1.0 - eta2 - h)
        + v1z * (eta1 + h) * (1.0 - eta2 - h)
        + v2z * (1.0 - eta1 - h) * (eta2 + h)
        + v3z * (eta1 + h) * (eta2 + h);

      P4EST_ASSERT ((p4est_locidx_t) quad_count < p4est->local_num_quadrants);

      lv0 = quadrant_to_local_vertex[4 * quad_count + 0];
      lv1 = quadrant_to_local_vertex[4 * quad_count + 1];
      lv2 = quadrant_to_local_vertex[4 * quad_count + 2];
      lv3 = quadrant_to_local_vertex[4 * quad_count + 3];

      P4EST_ASSERT (0 <= lv0 && lv0 < num_uniq_local_vertices);
      P4EST_ASSERT (0 <= lv1 && lv1 < num_uniq_local_vertices);
      P4EST_ASSERT (0 <= lv2 && lv2 < num_uniq_local_vertices);
      P4EST_ASSERT (0 <= lv3 && lv3 < num_uniq_local_vertices);

      vert_locations[lv0].x = w0x;
      vert_locations[lv0].y = w0y;
      vert_locations[lv0].z = w0z;
      P4EST_ASSERT (vert_locations[lv0].treeid == -1 ||
                    vert_locations[lv0].treeid == jt);
      vert_locations[lv0].treeid = jt;

      vert_locations[lv1].x = w1x;
      vert_locations[lv1].y = w1y;
      vert_locations[lv1].z = w1z;
      P4EST_ASSERT (vert_locations[lv1].treeid == -1 ||
                    vert_locations[lv1].treeid == jt);
      vert_locations[lv1].treeid = jt;

      vert_locations[lv2].x = w2x;
      vert_locations[lv2].y = w2y;
      vert_locations[lv2].z = w2z;
      P4EST_ASSERT (vert_locations[lv2].treeid == -1 ||
                    vert_locations[lv2].treeid == jt);
      vert_locations[lv2].treeid = jt;

      vert_locations[lv3].x = w3x;
      vert_locations[lv3].y = w3y;
      vert_locations[lv3].z = w3z;
      P4EST_ASSERT (vert_locations[lv3].treeid == -1 ||
                    vert_locations[lv3].treeid == jt);
      vert_locations[lv3].treeid = jt;
    }
  }

  qsort (vert_locations, num_uniq_local_vertices, sizeof (p4est_vert_t),
         p4est_vert_compare);

  /* Check to make sure that we don't have any duplicates in the list */
  for (kl = 0; kl < num_uniq_local_vertices - 1; ++kl) {
    SC_CHECK_ABORT (p4est_vert_compare (vert_locations + kl,
                                        vert_locations + kl + 1) != 0,
                    "local ordering not unique");
  }

  P4EST_FREE (vert_locations);
  p4est_nodes_destroy (nodes);
}

static int          weight_counter;
static int          weight_index;

static int
weight_once (p4est_t * p4est, p4est_topidx_t which_tree,
             p4est_quadrant_t * quadrant)
{
  if (weight_counter++ == weight_index) {
    return 1;
  }

  return 0;
}

int
main (int argc, char **argv)
{
  int                 rank;
  int                 mpiret;
  sc_MPI_Comm         mpicomm;
  p4est_t            *p4est;
  p4est_connectivity_t *connectivity;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  mpicomm = sc_MPI_COMM_WORLD;
  mpiret = sc_MPI_Comm_rank (mpicomm, &rank);
  SC_CHECK_MPI (mpiret);

  sc_init (mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  /* create connectivity and forest structures */
  connectivity = p4est_connectivity_new_star ();
  p4est = p4est_new_ext (mpicomm, connectivity, 15, 0, 0,
                         sizeof (user_data_t), init_fn, NULL);

  /* refine to make the number of elements interesting */
  p4est_refine (p4est, 1, refine_fn, init_fn);

  /* balance the forest */
  p4est_balance (p4est, P4EST_CONNECT_FULL, init_fn);

  /* do a uniform partition, include the weight function for testing */
  p4est_partition (p4est, 0, weight_one);

  p4est_check_local_order (p4est, connectivity);

  /* do a weighted partition with many zero weights */
  weight_counter = 0;
  weight_index = (rank == 1) ? 1342 : 0;
  p4est_partition (p4est, 0, weight_once);

  p4est_check_local_order (p4est, connectivity);

  /* clean up */
  p4est_destroy (p4est);
  p4est_connectivity_destroy (connectivity);

  /* create connectivity and forest structures */
  connectivity = p4est_connectivity_new_periodic ();
  p4est = p4est_new_ext (mpicomm, connectivity, 15, 0, 0,
                         sizeof (user_data_t), init_fn, NULL);

  /* refine to make the number of elements interesting */
  p4est_refine (p4est, 1, refine_fn, init_fn);

  /* balance the forest */
  p4est_balance (p4est, P4EST_CONNECT_FULL, init_fn);

  /* do a uniform partition, include the weight function for testing */
  p4est_partition (p4est, 0, weight_one);

  p4est_check_local_order (p4est, connectivity);

  /* clean up and exit */
  p4est_destroy (p4est);
  p4est_connectivity_destroy (connectivity);
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
