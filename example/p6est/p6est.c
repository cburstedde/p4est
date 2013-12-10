/*
  This file is part of p4est.
  p4est is a C library to manage a collection (a forest) of multiple
  connected adaptive quadtrees or octrees in parallel.

  Copyright (C) 2013 The University of Texas System
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

#include <p6est.h>
#include <p8est.h>
#include <p4est_extended.h>
#include <sc_containers.h>

void
p6est_qcoord_to_vertex (p8est_connectivity_t * connectivity,
                        p4est_topidx_t treeid,
                        p4est_qcoord_t x, p4est_qcoord_t y,
                        p6est_zcoord_t z,
                        double vxyz[3])
{
  p4est_qcoord_t      qz = (p4est_qcoord_t) z;

  p8est_qcoord_to_vertex (connectivity, treeid, x, y, qz, vxyz);
}

size_t
p6est_memory_used (p6est_t * p6est)
{
  size_t              size;

  size = p4est_memory_used (p6est->p4est);
  size += p4est_connectivity_memory_used (p6est->conn4);
  size += sc_array_memory_used (p6est->quads, 1);
  if (p6est->data_size > 0) {
    size += sc_mempool_memory_used (p6est->user_data_pool);
  }
  size += sc_mempool_memory_used (p6est->quadrant_pool);

  return size;
}

static int
p8est_connectivity_is_flat (p8est_connectivity_t *conn8)
{
  p4est_topidx_t num_trees = conn8->num_trees;
  p4est_topidx_t ti;
  p4est_topidx_t j;

  P4EST_ASSERT (p8est_connectivity_is_valid (conn8));
  for (ti = 0; ti < num_trees; ti++) {
    for (j = 4; j < 6; j++) {
      if (conn8->tree_to_tree[P8EST_FACES * ti + j] != ti) {
        return 0;
      }
      if (conn8->tree_to_face[P8EST_FACES * ti + j] != j) {
        return 0;
      }
    }
    for (j = 0; j < 8; j++) {
      if (conn8->tree_to_edge[P8EST_EDGES * ti + j] != -1) {
        return 0;
      }
    }
  }

  return 1;
}

static p4est_connectivity_t *
p8est_connectivity_flatten (p8est_connectivity_t *conn8)
{
  p4est_connectivity_t *conn4;
  p4est_topidx_t num_trees = conn8->num_trees;
  p4est_topidx_t num_edges = conn8->num_edges;
  p4est_topidx_t ti, tj;
  p4est_topidx_t *ett = conn8->ett_offset;
  int j;

  P4EST_ASSERT (p8est_connectivity_is_flat (conn8));
  conn4 = p4est_connectivity_new (0, num_trees, conn8->num_edges,
                                  conn8->ett_offset[conn8->num_edges]);

  for (ti = 0; ti < num_trees; ti++) {
    for (j = 0; j < P4EST_FACES; j++) {
      conn4->tree_to_tree[P4EST_FACES * ti + j] =
        conn8->tree_to_tree[P8EST_FACES * ti + j];
      conn4->tree_to_face[P4EST_FACES * ti + j] =
        (conn8->tree_to_face[P8EST_FACES * ti + j] % P8EST_FACES);
    }
    for (j = 0; j < P4EST_CHILDREN; j++) {
      conn4->tree_to_corner[P4EST_CHILDREN * ti + j] =
        conn8->tree_to_edge[P8EST_EDGES * ti + 8 + j];
    }
  }

  for (ti = 0; ti < num_edges; ti++) {
    conn4->ctt_offset[ti] = ett[ti];
    conn4->ctt_offset[ti + 1] = ett[ti + 1];
    for (tj = ett[ti]; tj < ett[ti + 1]; ti++) {
      conn4->corner_to_tree[tj] = conn8->edge_to_tree[tj];
      conn8->corner_to_corner[tj] = (conn8->edge_to_edge[tj] % P8EST_EDGES);
    }
  }

  return conn4;
}

typedef struct p6est_init_data
{
  int         min_zlevel;
  sc_array_t *quads;
  p6est_init_t init_fn;
  void       *user_pointer;
}
p6est_init_data_t;

static void
p6est_init_fn (p4est_t *p4est, p4est_topidx_t which_tree,
               p4est_quadrant_t *quadrant)
{
  p6est_t *p6est = (p6est_t *) p4est->user_pointer;
  p6est_init_data_t *init_data = (p6est_init_data_t *) p6est->user_pointer;
  int nquads = 1 << (init_data->min_zlevel);
  sc_array_t *quads = init_data->quads;
  size_t incount = quads->elem_count, zz;

  p6est->user_pointer = init_data->user_pointer;

  quadrant->p.piggy3.local_num = nquads;
  quadrant->p.piggy3.which_tree = quads->elem_count;

  sc_array_resize (quads, quads->elem_count + nquads);

  for (zz = incount; zz < quads->elem_count; zz++) {
    p6est_quadrant_t * colquad = p6est_quadrant_array_index (quads, zz);
    P6EST_QUADRANT_INIT (colquad);

    colquad->x = quadrant->x;
    colquad->y = quadrant->y;
    colquad->level = quadrant->level;
    colquad->zlevel = init_data->min_zlevel;
    colquad->z = (zz - incount) * P6EST_QUADRANT_LEN (colquad->zlevel);

    if (init_data->init_fn) {
      init_data->init_fn (p6est, which_tree, colquad);
    }
  }

  p6est->user_pointer = (void *) init_data;
}

p6est_t            *
p6est_new_ext (MPI_Comm mpicomm, p8est_connectivity_t * connectivity,
               p4est_locidx_t min_quadrants, int min_level, int min_zlevel,
               int fill_uniform, size_t data_size, p6est_init_t init_fn,
               void *user_pointer)
{
  p6est_t             *p6est = P4EST_ALLOC (p6est_t, 1);
  p4est_t             *p4est;
  p4est_connectivity_t *conn4;
  sc_array_t          *quads;
  sc_mempool_t        *user_data_pool;
  sc_mempool_t        *quadrant_pool;
  p6est_init_data_t    init_data;
  int                  mpiret, num_procs, rank;

  mpiret = MPI_Comm_size (mpicomm, &num_procs);
  SC_CHECK_MPI (mpiret);
  mpiret = MPI_Comm_rank (mpicomm, &rank);
  SC_CHECK_MPI (mpiret);
  /* get the 2D connectivity for the forest */
  conn4 = p8est_connectivity_flatten (connectivity);

  quads = sc_array_new (sizeof (p6est_quadrant_t));

  if (data_size > 0) {
    user_data_pool = sc_mempool_new (data_size);
  }
  else {
    user_data_pool = NULL;
  }

  quadrant_pool = sc_mempool_new (sizeof (p6est_quadrant_t));

  p6est->mpicomm = mpicomm;
  p6est->mpisize = num_procs;
  p6est->mpirank = rank;
  p6est->data_size = data_size;
  p6est->user_pointer = user_pointer;
  p6est->connectivity = connectivity;
  p6est->conn4 = conn4;
  p6est->quads = quads;
  p6est->user_data_pool = user_data_pool;

  P4EST_ASSERT (min_zlevel <= P6EST_QMAXLEVEL);
  init_data.min_zlevel = min_zlevel;
  init_data.quads      = quads;
  init_data.init_fn    = init_fn;
  p6est->user_pointer = &init_data;
  p4est = p4est_new_ext (mpicomm, conn4, min_quadrants / (1 << min_zlevel),
                         min_level, fill_uniform, 0, p6est_init_fn,
                         (void *) p6est);

  p6est->user_pointer = user_pointer;
  p6est->p4est = p4est;

  /* order the columns? */

  return p6est;
}

p6est_t            *
p6est_new (MPI_Comm mpicomm, p8est_connectivity_t * connectivity,
           size_t data_size, p6est_init_t init_fn, void *user_pointer)
{
  return p6est_new_ext (mpicomm, connectivity, 0, 0, 0, 1,
                        data_size, init_fn, user_pointer);
}

