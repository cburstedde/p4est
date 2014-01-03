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

void
p6est_quadrant_init_data (p6est_t * p6est, p4est_topidx_t which_tree,
                          p6est_quadrant_t * quad, p6est_init_t init_fn)
{
  if (p6est->data_size > 0) {
    quad->p.user_data = sc_mempool_alloc (p6est->user_data_pool);
  }
  else {
    quad->p.user_data = NULL;
  }
  if (init_fn != NULL) {
    init_fn (p6est, which_tree, quad);
  }
}

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

  quadrant->p.piggy3.local_num = (p4est_locidx_t) quads->elem_count;
  quadrant->p.piggy3.which_tree = (p4est_topidx_t) nquads;

  sc_array_resize (quads, quads->elem_count + nquads);

  for (zz = incount; zz < quads->elem_count; zz++) {
    p6est_quadrant_t * colquad = p6est_quadrant_array_index (quads, zz);
    P6EST_QUADRANT_INIT (colquad);

    colquad->x = quadrant->x;
    colquad->y = quadrant->y;
    colquad->level = quadrant->level;
    colquad->zlevel = init_data->min_zlevel;
    colquad->z = (zz - incount) * P6EST_QUADRANT_LEN (colquad->zlevel);

    p6est_quadrant_init_data (p6est, which_tree, colquad, init_data->init_fn);
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
  int                  quadpercol = (1 << min_zlevel);
  int                  i;

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
  p4est = p4est_new_ext (mpicomm, conn4, min_quadrants / quadpercol,
                         min_level, fill_uniform, 0, p6est_init_fn,
                         (void *) p6est);

  p6est->user_pointer = user_pointer;
  p6est->p4est = p4est;
  p6est->global_first_quadrant = P4EST_ALLOC (p4est_gloidx_t, num_procs + 1);
  for (i = 0; i <= num_procs; i++) {
    p6est->global_first_quadrant[i] = quadpercol * p4est->global_first_quadrant[i];
  }

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

void
p6est_quadrant_free_data (p6est_t * p6est, p6est_quadrant_t * quad)
{
  if (p6est->data_size > 0) {
    sc_mempool_free (p6est->user_data_pool, quad->p.user_data);
  }
  quad->p.user_data = NULL;
}

void
p6est_destroy (p6est_t *p6est)
{
  sc_array_t *quads = p6est->quads;
  size_t quadcount = quads->elem_count;
  size_t zz;

  for (zz = 0; zz < quadcount; zz++) {
    p6est_quadrant_t *colquad = p6est_quadrant_array_index (quads, zz);

    p6est_quadrant_free_data (p6est, colquad);
  }
  sc_array_destroy (p6est->quads);

  p4est_destroy (p6est->p4est);
  p4est_connectivity_destroy (p6est->conn4);
  if (p6est->user_data_pool != NULL) {
    sc_mempool_destroy (p6est->user_data_pool);
  }
  sc_mempool_destroy (p6est->quadrant_pool);
  P4EST_FREE (p6est->global_first_quadrant);
  P4EST_FREE (p6est);
}

p6est_t *
p6est_copy (p6est_t *input, int copy_data)
{
  p6est_t *p6est = P4EST_ALLOC (p6est_t, 1);
  size_t zz, qcount = input->quads->elem_count;

  memcpy (p6est, input, sizeof (p4est_t));
  p6est->conn4 = p4est_connectivity_new_copy (input->conn4->num_vertices,
                                              input->conn4->num_trees,
                                              input->conn4->num_corners,
                                              input->conn4->vertices,
                                              input->conn4->tree_to_vertex,
                                              input->conn4->tree_to_tree,
                                              input->conn4->tree_to_face,
                                              input->conn4->tree_to_corner,
                                              input->conn4->ctt_offset,
                                              input->conn4->corner_to_tree,
                                              input->conn4->corner_to_corner);
  p6est->quads = sc_array_new_size (input->quads->elem_size, input->quads->elem_count);
  sc_array_copy (input->quads, p6est->quads);
  p6est->p4est = p4est_copy (input->p4est, 0);
  p6est->p4est->user_pointer = p6est;
  if (copy_data && p6est->data_size > 0) {
    p6est->user_data_pool = sc_mempool_new (p6est->data_size);
  }
  else {
    p6est->data_size = 0;
  }
  p6est->quadrant_pool = sc_mempool_new (sizeof (p6est_quadrant_t));

  if (p6est->data_size > 0) {
    P4EST_ASSERT (copy_data);
    for (zz = 0; zz < qcount; zz++) {
      p6est_quadrant_t *inquad = p6est_quadrant_array_index (input->quads, zz);
      p6est_quadrant_t *outquad = p6est_quadrant_array_index (p6est->quads, zz);

      outquad->p.user_data = sc_mempool_alloc (p6est->user_data_pool);
      memcpy (outquad->p.user_data, inquad->p.user_data, p6est->data_size);
    }
  }
  p6est->global_first_quadrant = P4EST_ALLOC (p4est_gloidx_t, p6est->mpisize + 1);
  memcpy (p6est->global_first_quadrant,
          input->global_first_quadrant,
          (p6est->mpisize + 1) * sizeof (p4est_gloidx_t));
  return p6est;
}

void
p6est_reset_data (p6est_t *p6est, size_t data_size, p6est_init_t init_fn,
                  void *user_pointer)
{
  int                 doresize;
  size_t              zz, zy, first, last;
  p4est_topidx_t      jt;
  p4est_quadrant_t   *col;
  p6est_quadrant_t   *q;
  p4est_tree_t       *tree;
  sc_array_t         *tquadrants;

  doresize = (p6est->data_size != data_size);

  p6est->data_size = data_size;
  p6est->user_pointer = user_pointer;

  if (doresize) {
    if (p6est->user_data_pool != NULL) {
      sc_mempool_destroy (p6est->user_data_pool);
    }
    if (p6est->data_size > 0) {
      p6est->user_data_pool = sc_mempool_new (p6est->data_size);
    }
    else {
      p6est->user_data_pool = NULL;
    }
  }

  for (jt = p6est->p4est->first_local_tree; jt <= p6est->p4est->last_local_tree; ++jt) {
    tree = p4est_tree_array_index (p6est->p4est->trees, jt);
    tquadrants = &tree->quadrants;
    for (zz = 0; zz < tquadrants->elem_count; ++zz) {
      col = p4est_quadrant_array_index (tquadrants, zz);
      first = col->p.piggy3.local_num;
      last = first + col->p.piggy3.which_tree;
      for (zy = first; zy < last; zy++) {
        q = p6est_quadrant_array_index (p6est->quads, zy);
        if (doresize) {
          if (p6est->data_size > 0) {
            q->p.user_data = sc_mempool_alloc (p6est->user_data_pool);
          }
          else {
            q->p.user_data = NULL;
          }
        }
        if (init_fn != NULL) {
          init_fn (p6est, jt, q);
        }
      }
    }
  }
}

typedef void (*p6est_replace_t) (p6est_t *p6est,
                                 p4est_topidx_t which_tree,
                                 int num_outgoing,
                                 p6est_quadrant_t *outgoing[],
                                 int num_incoming,
                                 p6est_quadrant_t *incoming[]);

typedef struct p6est_refine_col_data
{
  p6est_init_t init_fn;
  p6est_replace_t replace_fn;
  void *user_pointer;
}
p6est_refine_col_data_t;

void
p6est_refine_column (p4est_t *p4est, p4est_topidx_t which_tree,
                     int num_outgoing,
                     p4est_quadrant_t * outgoing[],
                     int num_incoming,
                     p4est_quadrant_t * incoming[])
{
  p6est_t *p6est = (p6est_t *) p4est->user_pointer;
  p6est_refine_col_data_t *refine_col = (p6est_refine_col_data_t *) p6est->user_pointer;
  int nquads;
  size_t first, ifirst;
  int i, j;
  p6est_quadrant_t *oq, *q;

  p6est->user_pointer = refine_col->user_pointer;
  P4EST_ASSERT (num_outgoing == 1);
  P4EST_ASSERT (num_incoming == P4EST_CHILDREN);

  nquads = (int) outgoing[0]->p.piggy3.which_tree;
  first = (size_t) outgoing[0]->p.piggy3.local_num;

  for (i = 0; i < num_incoming; i++) {
    incoming[i]->p.piggy3.which_tree = (p4est_topidx_t) nquads;
    ifirst = p6est->quads->elem_count;
    incoming[i]->p.piggy3.local_num = (p4est_locidx_t) ifirst;
    sc_array_resize (p6est->quads, ifirst + nquads);
    for (j = 0; j < nquads; j++) {
      oq = p6est_quadrant_array_index (p6est->quads, first + j);
      q = p6est_quadrant_array_index (p6est->quads, ifirst + j);
      P6EST_QUADRANT_INIT (q);
      q->x = incoming[i]->x;
      q->y = incoming[i]->y;
      q->level = incoming[i]->level;
      q->z = oq->z;
      q->zlevel = oq->zlevel;
      p6est_quadrant_init_data (p6est, which_tree, q, refine_col->init_fn);
    }
  }

  if (refine_col->replace_fn != NULL) {
    for (j = 0; j < nquads; j++) {
      p6est_quadrant_t *inq[P4EST_CHILDREN];
      oq = p6est_quadrant_array_index (p6est->quads, first + j);

      ifirst = incoming[i]->p.piggy3.local_num;
      for (i = 0; i < P4EST_CHILDREN; i++) {
        inq[i] = p6est_quadrant_array_index (p6est->quads, ifirst + j);
      }

      refine_col->replace_fn (p6est, which_tree, 1, &oq, P4EST_CHILDREN,
                              inq);
    }
  }

  for (j = 0; j < nquads; j++) {
    oq = p6est_quadrant_array_index (p6est->quads, first + j);
    p6est_quadrant_free_data (p6est, oq);
  }
  p6est->user_pointer = (void *) refine_col;
}

void
p6est_compress_columns (p6est_t *p6est)
{
  size_t              zz, first;
  p4est_topidx_t      jt;
  p4est_quadrant_t   *col;
  p4est_tree_t       *tree;
  sc_array_t         *tquadrants;
  size_t              offset = 0;
  int                 count;

  for (jt = p6est->p4est->first_local_tree; jt <= p6est->p4est->last_local_tree; ++jt) {
    tree = p4est_tree_array_index (p6est->p4est->trees, jt);
    tquadrants = &tree->quadrants;
    for (zz = 0; zz < tquadrants->elem_count; ++zz) {
      col = p4est_quadrant_array_index (tquadrants, zz);
      first = col->p.piggy3.local_num;
      count = col->p.piggy3.which_tree;
      P4EST_ASSERT (first >= offset);
      if (first != offset) {
        memcpy (sc_array_index (p6est->quads, offset),
                sc_array_index (p6est->quads, first),
                count * sizeof (p6est_quadrant_t));
      }
      offset += count;
    }
  }
  sc_array_resize (p6est->quads, offset);
}

void
p6est_refine_ext (p6est_t *p6est, int refine_recursive, int allowed_level,
                  p4est_refine_t refine_fn, int allowed_zlevel,
                  p6est_refine_t zrefine_fn, p6est_init_t zinit_fn,
                  p6est_replace_t zreplace_fn)
{
#ifdef P4EST_DEBUG
  size_t              quadrant_pool_size, data_pool_size;
#endif
  int                 firsttime;
  int                 i, maxlevel;
  p4est_topidx_t      nt;
  size_t              incount, current, restpos, movecount;
  sc_list_t          *list;
  p4est_tree_t       *tree;
  p4est_quadrant_t   *q, *qalloc, *qpop;
  p4est_quadrant_t   *c0, *c1, *c2, *c3;
#ifdef P4_TO_P8
  p4est_quadrant_t   *c4, *c5, *c6, *c7;
#endif
  sc_array_t         *tquadrants;
  p4est_quadrant_t   *family[8];
  p4est_quadrant_t    parent, *pp = &parent;
  p6est_refine_col_data_t refine_col;
  void               *orig_user_pointer = p6est->user_pointer;
  size_t              zz, first;
  p4est_topidx_t      jt;
  p4est_quadrant_t   *col;
  size_t              offset = 0;
  int                 count;

  if (refine_fn) {
    refine_col.init_fn = zinit_fn;
    refine_col.replace_fn = zreplace_fn;
    refine_col.user_pointer = orig_user_pointer;

    p6est->user_pointer = (void *) &refine_col;
    p4est_refine_ext (p6est->p4est, refine_recursive, allowed_level, refine_fn,
                      NULL, p6est_refine_column);
    p6est->user_pointer = orig_user_pointer;

    p6est_compress_columns (p6est);
  }
  if (zrefine_fn) {
    for (jt = p6est->p4est->first_local_tree; jt <= p6est->p4est->last_local_tree; ++jt) {
      tree = p4est_tree_array_index (p6est->p4est->trees, jt);
      tquadrants = &tree->quadrants;
      for (zz = 0; zz < tquadrants->elem_count; ++zz) {
        col = p4est_quadrant_array_index (tquadrants, zz);
        first = col->p.piggy3.local_num;
        count = col->p.piggy3.which_tree;
      }
    }
  }
  /* update counts */
}
