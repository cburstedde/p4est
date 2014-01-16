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

p6est_connectivity_t *
p6est_connectivity_new (p4est_connectivity_t * conn4,
                        double *top_to_vertex, double height[3])
{
  p6est_connectivity_t *conn = P4EST_ALLOC (p6est_connectivity_t, 1);

  conn->conn4 = p4est_connectivity_new_copy (conn4->num_vertices,
                                             conn4->num_trees,
                                             conn4->num_corners,
                                             conn4->vertices,
                                             conn4->tree_to_vertex,
                                             conn4->tree_to_tree,
                                             conn4->tree_to_face,
                                             conn4->tree_to_corner,
                                             conn4->ctt_offset,
                                             conn4->corner_to_tree,
                                             conn4->corner_to_corner);

  if (top_to_vertex != NULL) {
    conn->top_to_vertex = P4EST_ALLOC (double, 3 * conn4->num_vertices);
    memcpy (conn->top_to_vertex, top_to_vertex,
            3 * conn4->num_vertices * sizeof (double));
  }
  else {
    conn->top_to_vertex = NULL;
    P4EST_ASSERT (height != NULL);
    conn->height[0] = height[0];
    conn->height[1] = height[1];
    conn->height[2] = height[2];
  }

  return conn;
}

void
p6est_connectivity_destroy (p6est_connectivity_t * conn)
{
  p4est_connectivity_destroy (conn->conn4);
  if (conn->top_to_vertex != NULL) {
    P4EST_FREE (conn->top_to_vertex);
  }
  P4EST_FREE (conn);
}

void
p6est_tree_get_vertices (p6est_connectivity_t * conn,
                         p4est_topidx_t which_tree, double vertices[24])
{
  const double       *btv = conn->conn4->vertices;
  const double       *ttv = conn->top_to_vertex;
  int                 i, j;

  P4EST_ASSERT (conn->conn4->num_vertices > 0);
  P4EST_ASSERT (btv != NULL);
  P4EST_ASSERT (which_tree >= 0 && which_tree < conn->conn4->num_trees);
  P4EST_ASSERT (vertices != NULL);

  memcpy (vertices, btv + which_tree * 12, 12 * sizeof (double));
  if (ttv != NULL) {
    memcpy (vertices + 12, ttv + which_tree * 12, 12 * sizeof (double));
  }
  else {
    memcpy (vertices + 12, btv + which_tree * 12, 12 * sizeof (double));

    for (i = 0; i < P4EST_CHILDREN; i++) {
      for (j = 0; j < 3; j++) {
        vertices[12 + 3 * i + j] += conn->height[j];
      }
    }
  }
}

void
p6est_qcoord_to_vertex (p6est_connectivity_t * conn,
                        p4est_topidx_t treeid,
                        p4est_qcoord_t x, p4est_qcoord_t y,
                        p4est_qcoord_t z, double vxyz[3])
{
  double              bottom[3], top[3];
  double              eta = (double) z / (double) P4EST_ROOT_LEN;

  p4est_qcoord_to_vertex (conn->conn4, treeid, x, y, bottom);
  if (conn->top_to_vertex != NULL) {
    double             *orig = conn->conn4->vertices;

    conn->conn4->vertices = conn->top_to_vertex;
    p4est_qcoord_to_vertex (conn->conn4, treeid, x, y, top);
    conn->conn4->vertices = orig;
  }
  else {
    top[0] = bottom[0] + conn->height[0];
    top[1] = bottom[1] + conn->height[1];
    top[2] = bottom[2] + conn->height[2];
  }
  vxyz[0] = (1. - eta) * bottom[0] + eta * top[0];
  vxyz[1] = (1. - eta) * bottom[1] + eta * top[1];
  vxyz[2] = (1. - eta) * bottom[2] + eta * top[2];
}

size_t
p6est_connectivity_memory_used (p6est_connectivity_t * conn)
{
  return
    p4est_connectivity_memory_used (conn->conn4) +
    conn->top_to_vertex == NULL ? 0 :
    (conn->conn4->num_vertices * 3 * sizeof (double));
}

size_t
p6est_memory_used (p6est_t * p6est)
{
  size_t              size;

  size = p4est_memory_used (p6est->columns);
  size += sc_array_memory_used (p6est->layers, 1);
  if (p6est->data_size > 0) {
    size += sc_mempool_memory_used (p6est->user_data_pool);
  }
  size += sc_mempool_memory_used (p6est->layer_pool);

  return size;
}

static int
p8est_connectivity_is_flat (p8est_connectivity_t * conn8)
{
  p4est_topidx_t      num_trees = conn8->num_trees;
  p4est_topidx_t      ti;
  p4est_topidx_t      j;

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

typedef struct p6est_init_data
{
  int                 min_zlevel;
  sc_array_t         *layers;
  p6est_init_t        init_fn;
  void               *user_pointer;
}
p6est_init_data_t;

void
p6est_layer_init_data (p6est_t * p6est, p4est_topidx_t which_tree,
                       p4est_quadrant_t * column,
                       p2est_quadrant_t * layer, p6est_init_t init_fn)
{
  if (p6est->data_size > 0) {
    layer->p.user_data = sc_mempool_alloc (p6est->user_data_pool);
  }
  else {
    layer->p.user_data = NULL;
  }
  if (init_fn != NULL) {
    init_fn (p6est, which_tree, column, layer);
  }
}

static void
p6est_init_fn (p4est_t * p4est, p4est_topidx_t which_tree,
               p4est_quadrant_t * col)
{
  p6est_t            *p6est = (p6est_t *) p4est->user_pointer;
  p6est_init_data_t  *init_data = (p6est_init_data_t *) p6est->user_pointer;
  int                 nlayers = 1 << (init_data->min_zlevel);
  sc_array_t         *layers = init_data->layers;
  size_t              incount = layers->elem_count, zz;
  size_t              last = incount + nlayers;

  /* we have to sneak the user_pointer in */
  p6est->user_pointer = init_data->user_pointer;

  P6EST_COLUMN_SET_RANGE (col, layers->elem_count, last)

    sc_array_resize (layers, last);

  for (zz = incount; zz < last; zz++) {
    p2est_quadrant_t   *layer = p2est_quadrant_array_index (layers, zz);
    P2EST_QUADRANT_INIT (layer);

    layer->level = init_data->min_zlevel;
    layer->z = (zz - incount) * P4EST_QUADRANT_LEN (layer->level);

    p6est_layer_init_data (p6est, which_tree, col, layer, init_data->init_fn);
  }

  /* we have to sneak the user_pointer out */
  p6est->user_pointer = (void *) init_data;
}

p6est_t            *
p6est_new_ext (MPI_Comm mpicomm, p6est_connectivity_t * connectivity,
               p4est_locidx_t min_quadrants, int min_level, int min_zlevel,
               int fill_uniform, size_t data_size, p6est_init_t init_fn,
               void *user_pointer)
{
  p6est_t            *p6est = P4EST_ALLOC (p6est_t, 1);
  p4est_t            *p4est;
  sc_array_t         *layers;
  sc_mempool_t       *user_data_pool;
  p6est_init_data_t   init_data;
  int                 mpiret, num_procs, rank;
  int                 quadpercol = (1 << min_zlevel);
  int                 i;

  mpiret = MPI_Comm_size (mpicomm, &num_procs);
  SC_CHECK_MPI (mpiret);
  mpiret = MPI_Comm_rank (mpicomm, &rank);
  SC_CHECK_MPI (mpiret);

  layers = sc_array_new (sizeof (p2est_quadrant_t));

  if (data_size > 0) {
    user_data_pool = sc_mempool_new (data_size);
  }
  else {
    user_data_pool = NULL;
  }

  p6est->layer_pool = sc_mempool_new (sizeof (p2est_quadrant_t));

  p6est->mpicomm = mpicomm;
  p6est->mpisize = num_procs;
  p6est->mpirank = rank;
  p6est->data_size = data_size;
  p6est->user_pointer = user_pointer;
  p6est->connectivity = connectivity;
  p6est->layers = layers;
  p6est->user_data_pool = user_data_pool;

  P4EST_ASSERT (min_zlevel <= P4EST_QMAXLEVEL);

  init_data.min_zlevel = min_zlevel;
  init_data.layers = layers;
  init_data.init_fn = init_fn;
  init_data.user_pointer = user_pointer;
  p6est->user_pointer = &init_data;

  p4est =
    p4est_new_ext (mpicomm, connectivity->conn4, min_quadrants / quadpercol,
                   min_level, fill_uniform, 0, p6est_init_fn, (void *) p6est);

  p6est->user_pointer = user_pointer;
  p6est->columns = p4est;
  p6est->global_first_layer = P4EST_ALLOC (p4est_gloidx_t, num_procs + 1);
  for (i = 0; i <= num_procs; i++) {
    p6est->global_first_layer[i] =
      quadpercol * p4est->global_first_quadrant[i];
  }

  return p6est;
}

p6est_t            *
p6est_new (MPI_Comm mpicomm, p6est_connectivity_t * connectivity,
           size_t data_size, p6est_init_t init_fn, void *user_pointer)
{
  return p6est_new_ext (mpicomm, connectivity, 0, 0, 0, 1,
                        data_size, init_fn, user_pointer);
}

void
p6est_layer_free_data (p6est_t * p6est, p2est_quadrant_t * layer)
{
  if (p6est->data_size > 0) {
    sc_mempool_free (p6est->user_data_pool, layer->p.user_data);
  }
  layer->p.user_data = NULL;
}

void
p6est_destroy (p6est_t * p6est)
{
  sc_array_t         *layers = p6est->layers;
  size_t              layercount = layers->elem_count;
  size_t              zz;

  for (zz = 0; zz < layercount; zz++) {
    p2est_quadrant_t   *layer = p2est_quadrant_array_index (layers, zz);

    p6est_layer_free_data (p6est, layer);
  }
  sc_array_destroy (p6est->layers);

  p4est_destroy (p6est->columns);
  if (p6est->user_data_pool != NULL) {
    sc_mempool_destroy (p6est->user_data_pool);
  }
  sc_mempool_destroy (p6est->layer_pool);
  P4EST_FREE (p6est->global_first_layer);
  P4EST_FREE (p6est);
}

p6est_t            *
p6est_copy (p6est_t * input, int copy_data)
{
  p6est_t            *p6est = P4EST_ALLOC (p6est_t, 1);
  size_t              zz, qcount = input->layers->elem_count;

  memcpy (p6est, input, sizeof (p6est_t));
  p6est->layers =
    sc_array_new_size (input->layers->elem_size, input->layers->elem_count);
  sc_array_copy (p6est->layers, input->layers);
  p6est->columns = p4est_copy (input->columns, 0);
  p6est->columns->user_pointer = p6est;
  if (copy_data && p6est->data_size > 0) {
    p6est->user_data_pool = sc_mempool_new (p6est->data_size);
  }
  else {
    p6est->data_size = 0;
  }
  p6est->layer_pool = sc_mempool_new (sizeof (p2est_quadrant_t));

  if (p6est->data_size > 0) {
    P4EST_ASSERT (copy_data);
    for (zz = 0; zz < qcount; zz++) {
      p2est_quadrant_t   *inlayer =
        p2est_quadrant_array_index (input->layers, zz);
      p2est_quadrant_t   *outlayer =
        p2est_quadrant_array_index (p6est->layers, zz);

      outlayer->p.user_data = sc_mempool_alloc (p6est->user_data_pool);
      memcpy (outlayer->p.user_data, inlayer->p.user_data, p6est->data_size);
    }
  }
  p6est->global_first_layer =
    P4EST_ALLOC (p4est_gloidx_t, p6est->mpisize + 1);
  memcpy (p6est->global_first_layer, input->global_first_layer,
          (p6est->mpisize + 1) * sizeof (p4est_gloidx_t));
  return p6est;
}

void
p6est_reset_data (p6est_t * p6est, size_t data_size, p6est_init_t init_fn,
                  void *user_pointer)
{
  int                 doresize;
  size_t              zz, zy, first, last;
  p4est_topidx_t      jt;
  p4est_quadrant_t   *col;
  p2est_quadrant_t   *q;
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

  for (jt = p6est->columns->first_local_tree;
       jt <= p6est->columns->last_local_tree; ++jt) {
    tree = p4est_tree_array_index (p6est->columns->trees, jt);
    tquadrants = &tree->quadrants;
    for (zz = 0; zz < tquadrants->elem_count; ++zz) {
      col = p4est_quadrant_array_index (tquadrants, zz);
      P6EST_COLUMN_GET_RANGE (col, &first, &last);
      for (zy = first; zy < last; zy++) {
        q = p2est_quadrant_array_index (p6est->layers, zy);
        if (doresize) {
          if (p6est->data_size > 0) {
            q->p.user_data = sc_mempool_alloc (p6est->user_data_pool);
          }
          else {
            q->p.user_data = NULL;
          }
        }
        if (init_fn != NULL) {
          init_fn (p6est, jt, col, q);
        }
      }
    }
  }
}

typedef void        (*p6est_replace_t) (p6est_t * p6est,
                                        p4est_topidx_t which_tree,
                                        int num_outcolumns,
                                        int num_outlayers,
                                        p4est_quadrant_t * outcolumns[],
                                        p2est_quadrant_t * outlayers[],
                                        int num_incolumns,
                                        int num_inlayers,
                                        p4est_quadrant_t * incolumns[],
                                        p2est_quadrant_t * inlayers[]);

typedef struct p6est_refine_col_data
{
  p6est_refine_column_t refine_col_fn;
  p6est_init_t        init_fn;
  p6est_replace_t     replace_fn;
  void               *user_pointer;
}
p6est_refine_col_data_t;

static int
p6est_refine_column_int (p4est_t * p4est, p4est_topidx_t which_tree,
                         p4est_quadrant_t * quadrant)
{
  p6est_t            *p6est = (p6est_t *) p4est->user_pointer;
  p6est_refine_col_data_t *refine_col =
    (p6est_refine_col_data_t *) p6est->user_pointer;
  int                 retval;

  /* sneak the user pointer in */
  p6est->user_pointer = refine_col->user_pointer;
  retval = refine_col->refine_col_fn (p6est, which_tree, quadrant);

  /* sneak the user pointer_out */
  p6est->user_pointer = (void *) refine_col;

  return retval;
}

static void
p6est_replace_column_split (p4est_t * p4est, p4est_topidx_t which_tree,
                            int num_outgoing, p4est_quadrant_t * outgoing[],
                            int num_incoming, p4est_quadrant_t * incoming[])
{
  p6est_t            *p6est = (p6est_t *) p4est->user_pointer;
  p6est_refine_col_data_t *refine_col =
    (p6est_refine_col_data_t *) p6est->user_pointer;
  int                 nlayers;
  size_t              first, last, ifirst, ilast;
  int                 i, j;
  p2est_quadrant_t   *oq, *q;

  /* sneak the user pointer in */
  p6est->user_pointer = refine_col->user_pointer;
  P4EST_ASSERT (num_outgoing == 1);
  P4EST_ASSERT (num_incoming == P4EST_CHILDREN);

  P6EST_COLUMN_GET_RANGE (outgoing[0], &first, &last);
  nlayers = last - first;

  for (i = 0; i < num_incoming; i++) {
    ifirst = p6est->layers->elem_count;
    ilast = ifirst + nlayers;
    sc_array_resize (p6est->layers, ilast);
    P6EST_COLUMN_SET_RANGE (incoming[i], ifirst, ilast);
    for (j = 0; j < nlayers; j++) {
      oq = p2est_quadrant_array_index (p6est->layers, first + j);
      q = p2est_quadrant_array_index (p6est->layers, ifirst + j);
      P2EST_QUADRANT_INIT (q);
      q->z = oq->z;
      q->level = oq->level;
      p6est_layer_init_data (p6est, which_tree, incoming[i], q,
                             refine_col->init_fn);
    }
  }

  if (refine_col->replace_fn != NULL) {
    for (j = 0; j < nlayers; j++) {
      p2est_quadrant_t   *inq[P4EST_CHILDREN];

      oq = p2est_quadrant_array_index (p6est->layers, first + j);
      for (i = 0; i < P4EST_CHILDREN; i++) {
        P6EST_COLUMN_GET_RANGE (incoming[i], &ifirst, &ilast);
        inq[i] = p2est_quadrant_array_index (p6est->layers, ifirst + j);
      }

      refine_col->replace_fn (p6est, which_tree,
                              1, 1,
                              outgoing, &oq,
                              P4EST_CHILDREN, P4EST_CHILDREN, incoming, inq);
    }
  }

  for (j = 0; j < nlayers; j++) {
    oq = p2est_quadrant_array_index (p6est->layers, first + j);
    p6est_layer_free_data (p6est, oq);
  }
  /* sneak the user pointer out */
  p6est->user_pointer = (void *) refine_col;
}

void
p6est_compress_columns (p6est_t * p6est)
{
  size_t              zz, zy, first, last;
  p4est_topidx_t      jt;
  p4est_quadrant_t   *col;
  p4est_tree_t       *tree;
  sc_array_t         *tquadrants;
  size_t              offset, nkeep;
  int                 count;
  size_t              this_col;
  p4est_t            *columns = p6est->columns;
  sc_array_t         *layers = p6est->layers;
  size_t             *newindex;
  sc_array_t         *na;
  size_t              nlayers = layers->elem_count;

  na = sc_array_new_size (sizeof (size_t), nlayers);
  newindex = (size_t *) na->array;
  for (zy = 0; zy < nlayers; zy++) {
    newindex[zy] = nlayers;
  }

  offset = 0;
  for (this_col = 0, jt = columns->first_local_tree;
       jt <= columns->last_local_tree; ++jt) {
    tree = p4est_tree_array_index (columns->trees, jt);
    tquadrants = &tree->quadrants;
    for (zz = 0; zz < tquadrants->elem_count; ++zz, this_col++) {
      col = p4est_quadrant_array_index (tquadrants, zz);
      P6EST_COLUMN_GET_RANGE (col, &first, &last);
      count = last - first;
      P6EST_COLUMN_SET_RANGE (col, offset, offset + count);
      for (zy = first; zy < last; zy++) {
        newindex[zy] = offset++;
      }
    }
  }
  nkeep = offset;
  P4EST_ASSERT (nkeep <= nlayers);

  for (zy = 0; zy < nlayers; zy++) {
    if (newindex[zy] == nlayers) {
      newindex[zy] = offset++;
    }
  }
  P4EST_ASSERT (offset == nlayers);

  sc_array_permute (layers, na, 0);
  sc_array_resize (p6est->layers, nkeep);
  sc_array_destroy (na);
}

void
p6est_update_offsets (p6est_t * p6est)
{
  int                 p, mpiret;
  p4est_gloidx_t     *gfl = p6est->global_first_layer;
  p4est_gloidx_t      mycount = p6est->layers->elem_count;
  p4est_gloidx_t      psum = 0, thiscount;

  mpiret = MPI_Allgather (&mycount, 1, P4EST_MPI_GLOIDX, gfl, 1,
                          P4EST_MPI_GLOIDX, p6est->mpicomm);
  SC_CHECK_MPI (mpiret);

  for (p = 0; p < p6est->mpisize; p++) {
    thiscount = gfl[p];
    gfl[p] += psum;
    psum += thiscount;
  }
  gfl[p6est->mpisize] = psum;
}

void
p6est_refine_columns_ext (p6est_t * p6est, int refine_recursive,
                          int allowed_level, p6est_refine_column_t refine_fn,
                          p6est_init_t init_fn, p6est_replace_t replace_fn)
{
  p6est_refine_col_data_t refine_col;
  void               *orig_user_pointer = p6est->user_pointer;

  refine_col.refine_col_fn = refine_fn;
  refine_col.init_fn = init_fn;
  refine_col.replace_fn = replace_fn;
  refine_col.user_pointer = orig_user_pointer;

  p6est->user_pointer = (void *) &refine_col;
  p4est_refine_ext (p6est->columns, refine_recursive, allowed_level,
                    p6est_refine_column_int, NULL,
                    p6est_replace_column_split);
  p6est->user_pointer = orig_user_pointer;

  p6est_compress_columns (p6est);
  p6est_update_offsets (p6est);
}

void
p6est_refine_layers_ext (p6est_t * p6est, int refine_recursive,
                         int allowed_level, p6est_refine_layer_t refine_fn,
                         p6est_init_t init_fn, p6est_replace_t replace_fn)
{
  p4est_t            *columns = p6est->columns;
  sc_array_t         *layers = p6est->layers;
  sc_array_t         *newcol = sc_array_new (sizeof (p2est_quadrant_t));
  p4est_topidx_t      jt;
  p4est_tree_t       *tree;
  sc_array_t         *tquadrants;
  p4est_quadrant_t   *col;
  p2est_quadrant_t   *q, *newq;
  p2est_quadrant_t    nextq[P4EST_MAXLEVEL];
  p2est_quadrant_t    c[2];
  p2est_quadrant_t    p, *parent = &p;
  p2est_quadrant_t   *child[2];
  size_t              first, last, zz, current, old_count;
  int                 any_change;
  int                 level;
  int                 stop_recurse;

  for (jt = columns->first_local_tree; jt <= columns->last_local_tree; ++jt) {
    tree = p4est_tree_array_index (columns->trees, jt);
    tquadrants = &tree->quadrants;

    for (zz = 0; zz < tquadrants->elem_count; ++zz) {
      col = p4est_quadrant_array_index (tquadrants, zz);
      P6EST_COLUMN_GET_RANGE (col, &first, &last);

      any_change = 0;

      for (current = first; current < last; current++) {
        q = p2est_quadrant_array_index (layers, current);
        stop_recurse = 0;
        level = q->level;
        parent = q;
        for (;;) {
          if (!stop_recurse && refine_fn (p6est, jt, col, parent) &&
              (allowed_level < 0 || (int) parent->level < allowed_level)) {
            level++;
            any_change = 1;
            c[0] = *parent;
            c[0].level = level;
            c[1] = *parent;
            c[1].level = level;
            c[1].z += P4EST_QUADRANT_LEN (level);
            child[0] = &c[0];
            child[1] = &c[1];
            p6est_layer_init_data (p6est, jt, col, child[0], init_fn);
            p6est_layer_init_data (p6est, jt, col, child[1], init_fn);
            if (replace_fn != NULL) {
              replace_fn (p6est, jt, 1, 1, &col, &parent, 1, 2, &col, child);
            }
            P4EST_ASSERT (parent->p.user_data != NULL);
            p6est_layer_free_data (p6est, parent);
            p = c[0];
            parent = &p;
            nextq[level] = c[1];
            stop_recurse = !refine_recursive;
          }
          else {
            /* parent is accepted */
            newq = (p2est_quadrant_t *) sc_array_push (newcol);
            *newq = *parent;
            if (parent == &p) {
              parent = &nextq[level];
            }
            else {
              while (--level > q->level && parent->z > nextq[level].z) {
              }
              if (level <= q->level) {
                break;
              }
              parent = &(nextq[level]);
            }
          }
        }
      }
      if (any_change) {
        old_count = layers->elem_count;
        sc_array_resize (layers, old_count + newcol->elem_count);
        memcpy (sc_array_index (layers, old_count),
                sc_array_index (newcol, 0),
                newcol->elem_size * newcol->elem_count);
        P6EST_COLUMN_SET_RANGE (col, old_count,
                                old_count + newcol->elem_count);
      }
      sc_array_truncate (newcol);
    }
  }
  sc_array_destroy (newcol);
  p6est_compress_columns (p6est);
  p6est_update_offsets (p6est);
}

void
p6est_refine_columns (p6est_t * p6est, int refine_recursive,
                      p6est_refine_column_t refine_fn, p6est_init_t init_fn)
{
  p6est_refine_columns_ext (p6est, refine_recursive, -1, refine_fn, init_fn,
                            NULL);
}

void
p6est_refine_layers (p6est_t * p6est, int refine_recursive,
                     p6est_refine_layer_t refine_fn, p6est_init_t init_fn)
{
  p6est_refine_layers_ext (p6est, refine_recursive, -1, refine_fn, init_fn,
                           NULL);
}

typedef struct p6est_coarsen_col_data
{
  p6est_coarsen_column_t coarsen_col_fn;
  p6est_init_t        init_fn;
  p6est_replace_t     replace_fn;
  void               *user_pointer;
  sc_array_t         *work_array;       /* for merging columns */
}
p6est_coarsen_col_data_t;

static int
p6est_coarsen_column_int (p4est_t * p4est, p4est_topidx_t which_tree,
                          p4est_quadrant_t * quadrants[])
{
  p6est_t            *p6est = (p6est_t *) p4est->user_pointer;
  p6est_coarsen_col_data_t *coarsen_col =
    (p6est_coarsen_col_data_t *) p6est->user_pointer;
  int                 retval;

  /* sneak the user pointer in */
  p6est->user_pointer = coarsen_col->user_pointer;
  retval = coarsen_col->coarsen_col_fn (p6est, which_tree, quadrants);

  /* sneak the user pointer_out */
  p6est->user_pointer = (void *) coarsen_col;

  return retval;
}

static int
p2est_quadrant_compare (const void *A, const void *B)
{
  const p2est_quadrant_t *a = (const p2est_quadrant_t *) A;
  const p2est_quadrant_t *b = (const p2est_quadrant_t *) B;
  p4est_qcoord_t      diff = a->z - b->z;

  if (!diff) {
    return (int) a->level - b->level;
  }

  return (int) diff;
}

static int
p2est_quadrant_is_equal (p2est_quadrant_t * a, p2est_quadrant_t * b)
{
  return (a->z == b->z && a->level == b->level);
}

static int
p2est_quadrant_is_ancestor (p2est_quadrant_t * a, p2est_quadrant_t * b)
{
  if (a->level >= b->level) {
    return 0;
  }

  return (b->z >= a->z && b->z < a->z + P4EST_QUADRANT_LEN (a->level));
}

static void
p6est_coarsen_all_layers (p6est_t * p6est, p4est_topidx_t which_tree,
                          p4est_quadrant_t * column, int ancestor_level,
                          sc_array_t * descendants,
                          int coarsen_recursive,
                          int callback_orphan,
                          p6est_coarsen_layer_t coarsen_fn,
                          p6est_init_t init_fn, p6est_replace_t replace_fn)
{
  p2est_quadrant_t    prevq[P4EST_QMAXLEVEL];
  size_t              old_count = descendants->elem_count;
  size_t              zz, new_count = 0;
  p2est_quadrant_t   *q, *r, a, s, *sib[2];
  int                 i, stackheight;
#ifdef P4EST_DEBUG
  size_t              mcount = p6est->user_data_pool->elem_count;
  p4est_qcoord_t      startpos, endpos;
#endif

  P4EST_ASSERT (old_count > 0);

  q = p2est_quadrant_array_index (descendants, 0);
#ifdef P4EST_DEBUG
  startpos = q->z;
#endif

  if (q->level == ancestor_level) {
    P4EST_ASSERT (old_count == 1);
    if (coarsen_fn != NULL && callback_orphan) {
      sib[0] = q;
      sib[1] = NULL;
      coarsen_fn (p6est, which_tree, column, sib);
    }
    return;
  }

  P4EST_ASSERT (old_count > 1);

  /* the stack contains only first siblings until it is time to move the stack
   * back to the array */
  P4EST_ASSERT (!(q->z & P4EST_QUADRANT_LEN (q->level)));

  stackheight = 0;
  zz = 1;
  for (;;) {
    /* if this is a second sibling */
    if (q->z & P4EST_QUADRANT_LEN (q->level)) {
      if (!stackheight) {
        prevq[stackheight++] = *q;
      }
      else {
        /* pop the first sibling off the stack */
        sib[0] = &prevq[stackheight - 1];
        s = *q;
        sib[1] = &s;
        P4EST_ASSERT (sib[0]->level == sib[1]->level);
        P4EST_ASSERT (sib[0]->z + P4EST_QUADRANT_LEN (sib[0]->level) ==
                      sib[1]->z);
        /* test if we will coarsen this pair */
        if (coarsen_fn == NULL || coarsen_fn (p6est, which_tree, column, sib)) {
          /* coarsen */
          a = *sib[0];
          q = &a;
          q->level--;
          P4EST_ASSERT (p2est_quadrant_is_ancestor (q, sib[0]));
          P4EST_ASSERT (p2est_quadrant_is_ancestor (q, sib[1]));
          p6est_layer_init_data (p6est, which_tree, column, q, init_fn);
          if (replace_fn) {
            replace_fn (p6est, which_tree, 1, 2, &column, sib, 1, 1, &column,
                        &q);
          }
          p6est_layer_free_data (p6est, sib[0]);
          p6est_layer_free_data (p6est, sib[1]);
          /* if we are doing recursive coarsening, do not advance */
          if (coarsen_recursive) {
            stackheight--;
          }
          else {
            /* put the parent on the end of the stack */
            prevq[stackheight - 1] = *q;
          }
        }
        else {
          prevq[stackheight++] = *q;
        }
      }
    }
    else {
      prevq[stackheight++] = *q;
      if (q->level > ancestor_level) {
        P4EST_ASSERT (zz < old_count);
        q = p2est_quadrant_array_index (descendants, zz++);
      }
    }

    if (stackheight && p2est_quadrant_is_equal (&prevq[stackheight - 1], q)) {
      /* clear the stack into the array from the bottom up */
      for (i = 0; i < stackheight; i++) {
        r = p2est_quadrant_array_index (descendants, new_count++);
        *r = prevq[i];
        if (coarsen_fn != NULL && callback_orphan) {
          sib[0] = r;
          sib[1] = NULL;
          coarsen_fn (p6est, which_tree, column, sib);
        }
      }
      stackheight = 0;
      if (zz == old_count) {
        break;
      }
      else {
        q = p2est_quadrant_array_index (descendants, zz++);
      }
    }
  }

  sc_array_resize (descendants, new_count);

#ifdef P4EST_DEBUG
  P4EST_ASSERT (mcount - p6est->user_data_pool->elem_count ==
                (old_count - new_count));

  q = p2est_quadrant_array_index (descendants, new_count - 1);
  endpos = q->z + P4EST_QUADRANT_LEN (q->level);
  P4EST_ASSERT (endpos - startpos == P4EST_QUADRANT_LEN (ancestor_level));

#endif
}

static void
p6est_replace_column_join (p4est_t * p4est, p4est_topidx_t which_tree,
                           int num_outgoing, p4est_quadrant_t * outgoing[],
                           int num_incoming, p4est_quadrant_t * incoming[])
{
  p6est_t            *p6est = (p6est_t *) p4est->user_pointer;
  p6est_coarsen_col_data_t *coarsen_col =
    (p6est_coarsen_col_data_t *) p6est->user_pointer;
  size_t              nlayers[P4EST_CHILDREN];
  size_t              first[P4EST_CHILDREN], last[P4EST_CHILDREN];
  size_t              zw[P4EST_CHILDREN];
  size_t              view_first, view_count;
  size_t              new_first, new_count, new_last;
  int                 j;
  p2est_quadrant_t   *q[P4EST_CHILDREN], *p;
  sc_array_t         *layers = p6est->layers;
  sc_array_t         *work_array = coarsen_col->work_array;
  sc_array_t          view;
  p6est_init_t        init_fn = coarsen_col->init_fn;
  p6est_replace_t     replace_fn = coarsen_col->replace_fn;

  /* sneak the user pointer in */
  p6est->user_pointer = coarsen_col->user_pointer;

  P4EST_ASSERT (num_outgoing == P4EST_CHILDREN);
  P4EST_ASSERT (num_incoming == 1);
  P4EST_ASSERT (!work_array->elem_count);

  new_first = layers->elem_count;
  new_count = 0;

  /* zero counters for each column */
  for (j = 0; j < num_outgoing; j++) {
    zw[j] = 0;
    P6EST_COLUMN_GET_RANGE (outgoing[j], &first[j], &last[j]);
    nlayers[j] = last[j] - first[j];
  }

  while (zw[0] < nlayers[0]) {
    for (j = 0; j < num_outgoing; j++) {
      P4EST_ASSERT (zw[j] < nlayers[j]);
      q[j] = p2est_quadrant_array_index (layers, first[j] + zw[j]);
    }
    p = (p2est_quadrant_t *) sc_array_push (work_array);
    *p = *q[0];
    p6est_layer_init_data (p6est, which_tree, incoming[0], p, init_fn);
    for (j = 1; j < num_outgoing; j++) {
      P4EST_ASSERT (q[j]->z == p->z);
      if (q[j]->level < p->level) {
        *p = *q[j];
      }
    }
    for (j = 0; j < num_outgoing; j++) {
      if (q[j]->level > p->level) {
        P4EST_ASSERT (p2est_quadrant_is_ancestor (p, q[j]));
        view_count = 1;
        view_first = first[j] + zw[j];
        P4EST_ASSERT (zw[j] < nlayers[j] - 1);
        while (++zw[j] < nlayers[j] &&
               p2est_quadrant_is_ancestor (p,
                                           p2est_quadrant_array_index (layers,
                                                                       first
                                                                       [j] +
                                                                       zw
                                                                       [j])))
        {
          view_count++;
        }
        sc_array_init_view (&view, layers, view_first, view_count);
        /* coarsen within this column */
        p6est_coarsen_all_layers (p6est, which_tree, outgoing[j], p->level,
                                  &view, 1, 0, NULL, init_fn, replace_fn);
        q[j] = p2est_quadrant_array_index (&view, 0);
      }
      else {
        zw[j]++;
      }
    }
    if (replace_fn != NULL) {
      replace_fn (p6est, which_tree, P4EST_CHILDREN, 1, outgoing, q, 1, 1,
                  incoming, &p);
    }
    for (j = 0; j < num_outgoing; j++) {
      p6est_layer_free_data (p6est, q[j]);
    }
  }

  new_count = work_array->elem_count;
  new_last = new_first + new_count;
  P6EST_COLUMN_SET_RANGE (incoming[0], new_first, new_last);

  /* create the new column */
  sc_array_resize (layers, new_last);
  memcpy (sc_array_index (layers, new_first),
          sc_array_index (work_array, 0), new_count * work_array->elem_size);

  sc_array_truncate (work_array);
  /* sneak the user pointer out */
  p6est->user_pointer = (void *) coarsen_col;
}

void
p6est_coarsen_columns_ext (p6est_t * p6est, int coarsen_recursive,
                           int callback_orphans,
                           p6est_coarsen_column_t coarsen_fn,
                           p6est_init_t init_fn, p6est_replace_t replace_fn)
{
  p6est_coarsen_col_data_t coarsen_col;
  void               *orig_user_pointer = p6est->user_pointer;

  coarsen_col.coarsen_col_fn = coarsen_fn;
  coarsen_col.init_fn = init_fn;
  coarsen_col.replace_fn = replace_fn;
  coarsen_col.user_pointer = orig_user_pointer;
  coarsen_col.work_array = sc_array_new (sizeof (p2est_quadrant_t));

  p6est->user_pointer = (void *) &coarsen_col;
  p4est_coarsen_ext (p6est->columns, coarsen_recursive, callback_orphans,
                     p6est_coarsen_column_int, NULL,
                     p6est_replace_column_join);
  p6est->user_pointer = orig_user_pointer;

  sc_array_destroy (coarsen_col.work_array);
  p6est_compress_columns (p6est);
  p6est_update_offsets (p6est);
}

void
p6est_coarsen_columns (p6est_t * p6est, int coarsen_recursive,
                       p6est_coarsen_column_t coarsen_fn,
                       p6est_init_t init_fn)
{
  p6est_coarsen_columns_ext (p6est, coarsen_recursive, 0,
                             coarsen_fn, init_fn, NULL);
}

void
p6est_coarsen_layers_ext (p6est_t * p6est, int coarsen_recursive,
                          int callback_orphans,
                          p6est_coarsen_layer_t coarsen_fn,
                          p6est_init_t init_fn, p6est_replace_t replace_fn)
{
  p4est_t            *columns = p6est->columns;
  sc_array_t         *layers = p6est->layers;
  sc_array_t          view;
  p4est_topidx_t      jt;
  p4est_tree_t       *tree;
  sc_array_t         *tquadrants;
  p4est_quadrant_t   *col;
  size_t              first, last, zz, count;

  for (jt = columns->first_local_tree; jt <= columns->last_local_tree; ++jt) {
    tree = p4est_tree_array_index (columns->trees, jt);
    tquadrants = &tree->quadrants;

    for (zz = 0; zz < tquadrants->elem_count; ++zz) {
      col = p4est_quadrant_array_index (tquadrants, zz);
      P6EST_COLUMN_GET_RANGE (col, &first, &last);

      count = last - first;
      sc_array_init_view (&view, layers, first, count);
      p6est_coarsen_all_layers (p6est, jt, col, 0, &view,
                                coarsen_recursive, callback_orphans,
                                coarsen_fn, init_fn, replace_fn);
      P4EST_ASSERT (view.elem_count > 0);
      P4EST_ASSERT (view.elem_count <= count);
      last = first + view.elem_count;
      P6EST_COLUMN_SET_RANGE (col, first, last);
    }
  }
  p6est_compress_columns (p6est);
  p6est_update_offsets (p6est);
  P4EST_ASSERT (p6est->user_data_pool->elem_count == layers->elem_count);
}

void
p6est_coarsen_layers (p6est_t * p6est, int coarsen_recursive,
                      p6est_coarsen_layer_t coarsen_fn, p6est_init_t init_fn)
{
  p6est_coarsen_layers_ext (p6est, coarsen_recursive, 0, coarsen_fn, init_fn,
                            NULL);
}

typedef struct p6est_weight_column
{
  p6est_weight_t      layer_weight_fn;
  void               *user_pointer;
}
p6est_weight_column_t;

static int
p6est_weight_fn (p4est_t * p4est, p4est_topidx_t which_tree,
                 p4est_quadrant_t * q)
{
  p6est_t            *p6est = (p6est_t *) p4est->user_pointer;
  p6est_weight_column_t *wc = (p6est_weight_column_t *) p6est->user_pointer;
  void               *orig_pointer = wc->user_pointer;
  size_t              first, last, zz;
  int                 weight = 0;

  p6est->user_pointer = orig_pointer;

  P6EST_COLUMN_GET_RANGE (q, &first, &last);

  if (wc->layer_weight_fn == NULL) {
    weight = last - first;
  }
  else {
    for (zz = first; zz < last; zz++) {
      p2est_quadrant_t   *layer =
        p2est_quadrant_array_index (p6est->layers, zz);

      weight += wc->layer_weight_fn (p6est, which_tree, q, layer);
    }
  }

  p6est->user_pointer = (void *) wc;

  return weight;
}

static int
int_compare_overlap (const void *key, const void *array)
{
  const int           search = *((const int *) key);
  const int          *range = (const int *) array;

  if (search >= range[0]) {
    if (search < range[1]) {
      return 0;
    }
    else {
      return 1;
    }
  }
  else {
    return -1;
  }
}

p4est_gloidx_t
p6est_partition_ext (p6est_t * p6est, int partition_for_coarsening,
                     p6est_weight_t weight_fn)
{
  p6est_weight_column_t wc;
  size_t              zz, offset, count, first, last;
  p4est_gloidx_t      my_count;
  p4est_gloidx_t     *new_gfl, *old_gfl = p6est->global_first_layer;
  p4est_topidx_t      jt;
  p4est_tree_t       *tree;
  p4est_t            *columns = p6est->columns;
  p4est_quadrant_t   *col;
  sc_array_t         *tquadrants;
  sc_array_t          new_gfl_bsearch, old_gfl_bsearch;
  int                 mpiret, p, rank = p6est->mpirank;
  p4est_gloidx_t      psum, thiscount;
  ssize_t             search;
  sc_array_t         *recv, *send;
  size_t              data_size = p6est->data_size;
  int                 overlap;
  int                 mpisize = p6est->mpisize;
  int                 send_count, recv_count;
  sc_array_t         *send_requests;
  sc_array_t         *recv_requests;
  MPI_Request        *req;
  void               *layer_data;
  p2est_quadrant_t   *layer;
  p2est_quadrant_t   *packedlayer;
  int                 i, nrecv, nsend, nleft;
  int                 outcount, *array_of_indices;
  size_t              self_offset, self_count;
  sc_array_t         *new_layers;
  int                 local_offset, local_last;
  sc_array_t         *old_layers = p6est->layers;
  p4est_gloidx_t      shipped;

  /* wrap the p6est_weight_t in a p4est_weight_t */
  wc.layer_weight_fn = weight_fn;
  wc.user_pointer = p6est->user_pointer;
  p6est->user_pointer = &wc;
  /* repartition the columns */
  p4est_partition_ext (p6est->columns, partition_for_coarsening,
                       p6est_weight_fn);

  /* the column counts (last - first) are correct, but not the offsets. update
   * the offsets */
  offset = 0;
  for (jt = columns->first_local_tree; jt <= columns->last_local_tree; ++jt) {
    tree = p4est_tree_array_index (columns->trees, jt);
    tquadrants = &tree->quadrants;

    for (zz = 0; zz < tquadrants->elem_count; ++zz) {
      col = p4est_quadrant_array_index (tquadrants, zz);
      P6EST_COLUMN_GET_RANGE (col, &first, &last);
      count = last - first;
      first = offset;
      last = offset + count;
      offset += count;
      P6EST_COLUMN_SET_RANGE (col, first, last);
    }
  }
  my_count = (p4est_gloidx_t) offset;

  /* calculate the new global_first_layer */
  new_gfl = P4EST_ALLOC (p4est_gloidx_t, p6est->mpisize + 1);
  mpiret = MPI_Allgather (&my_count, 1, P4EST_MPI_GLOIDX, new_gfl, 1,
                          P4EST_MPI_GLOIDX, p6est->mpicomm);
  SC_CHECK_MPI (mpiret);

  psum = 0;
  for (p = 0; p < p6est->mpisize; p++) {
    thiscount = new_gfl[p];
    new_gfl[p] += psum;
    psum += thiscount;
  }
  new_gfl[p6est->mpisize] = psum;
  P4EST_ASSERT (new_gfl[p6est->mpisize] == old_gfl[p6est->mpisize]);

  if (old_gfl[rank] == new_gfl[rank] &&
      old_gfl[rank + 1] == new_gfl[rank + 1]) {
    /* if my range is unchanged, I neither send nor receive to anyone */
    p6est->global_first_layer = new_gfl;
    P4EST_FREE (old_gfl);
    return 0;
  }

  /* calculate the global number of shipped (received) layers */
  shipped = 0;
  for (i = 0; i < mpisize; i++) {
    p4est_locidx_t      same, diff;

    same =
      SC_MIN (new_gfl[i + 1], old_gfl[i + 1]) - SC_MAX (new_gfl[i],
                                                        old_gfl[i]);
    diff = my_count - same;
    shipped += diff;
  }

  new_layers = sc_array_new_size (sizeof (p2est_quadrant_t), my_count);

  /* initialize views that we can use to bsearch */
  sc_array_init_data (&new_gfl_bsearch, new_gfl, 0, (size_t) mpisize);
  sc_array_init_data (&old_gfl_bsearch, old_gfl, 0, (size_t) mpisize);

  /* create the mpi recv requests */
  recv_requests = sc_array_new (sizeof (MPI_Request));

  /* create the receive buffer */
  if (!data_size) {
    recv = sc_array_new_view (new_layers, 0, my_count);
  }
  else {
    recv =
      sc_array_new_size (sizeof (p2est_quadrant_t) + data_size, my_count);
  }

  /* find the first proc that owns layers to send to me */
  search = sc_array_bsearch (&old_gfl_bsearch, new_gfl + rank,
                             int_compare_overlap);
  P4EST_ASSERT (search >= 0 && search < mpisize);
  overlap = search;
  P4EST_ASSERT (old_gfl[overlap] <= new_gfl[rank] &&
                new_gfl[rank] < old_gfl[overlap + 1]);

  offset = new_gfl[rank];
  self_offset = my_count;
  /* for every proc whose old range overlaps with this rank's new range */
  while (overlap < mpisize &&
         old_gfl[overlap] < new_gfl[rank + 1] &&
         new_gfl[rank] < old_gfl[overlap + 1]) {
    /* get the count that this proc will send this rank */
    recv_count = SC_MIN (new_gfl[rank + 1], old_gfl[overlap + 1]) - offset;
    local_offset = offset - new_gfl[rank];
    if (overlap != rank) {
      if (recv_count) {
        /* post the receive */
        req = (MPI_Request *) sc_array_push (recv_requests);
        mpiret = MPI_Irecv (sc_array_index (recv, local_offset),
                            (int) (recv_count * recv->elem_size), MPI_BYTE,
                            overlap, P6EST_COMM_PARTITION, p6est->mpicomm,
                            req);
        SC_CHECK_MPI (mpiret);
      }
    }
    else {
      self_offset = local_offset;
      self_count = recv_count;
    }

    offset += recv_count;

    overlap++;
  }
  P4EST_ASSERT (offset == new_gfl[rank + 1]);
  nrecv = (int) recv_requests->elem_count;

  /* create the mpi send requests */
  send_requests = sc_array_new (sizeof (MPI_Request));

  /* create the send buffer */
  if (!data_size) {
    send = sc_array_new_view (old_layers, 0, old_layers->elem_count);
  }
  else {
    send = sc_array_new_size (sizeof (p2est_quadrant_t) + data_size,
                              old_layers->elem_count);
  }

  /* find the first proc that owns layers to send to me */
  search = sc_array_bsearch (&new_gfl_bsearch, old_gfl + rank,
                             int_compare_overlap);
  P4EST_ASSERT (search >= 0 && search < mpisize);
  overlap = search;
  P4EST_ASSERT (new_gfl[overlap] <= old_gfl[rank] &&
                old_gfl[rank] < new_gfl[overlap + 1]);

  offset = old_gfl[rank];
  /* for every proc whose new range overlaps with this rank's old range */
  while (overlap < mpisize &&
         new_gfl[overlap] < old_gfl[rank + 1] &&
         old_gfl[rank] < new_gfl[overlap + 1]) {
    /* get the count that this proc will send this rank */
    send_count = SC_MIN (old_gfl[rank + 1], new_gfl[overlap + 1]) - offset;
    local_offset = offset - old_gfl[rank];
    if (overlap != rank) {
      if (data_size) {
        /* pack data for shipping */
        for (zz = 0; zz < send_count; zz++) {
          layer = p2est_quadrant_array_index (old_layers, zz + local_offset);
          packedlayer =
            (p2est_quadrant_t *) sc_array_index (send, zz + local_offset);
          *packedlayer = *layer;
          packedlayer->p.user_data = NULL;
          layer_data = (void *) (packedlayer + 1);
          memcpy (layer_data, layer->p.user_data, data_size);
          p6est_layer_free_data (p6est, layer);
        }
      }
      if (send_count) {
        /* post the send */
        req = (MPI_Request *) sc_array_push (send_requests);
        mpiret = MPI_Isend (sc_array_index (send, local_offset),
                            (int) (send_count * send->elem_size), MPI_BYTE,
                            overlap, P6EST_COMM_PARTITION, p6est->mpicomm,
                            req);
        SC_CHECK_MPI (mpiret);
      }
    }
    else {
      P4EST_ASSERT (send_count == self_count);
      /* copy locally */
      memcpy (sc_array_index (new_layers, self_offset),
              sc_array_index (old_layers, offset - old_gfl[rank]),
              self_count * new_layers->elem_size);
    }
    offset += send_count;
    overlap++;
  }
  P4EST_ASSERT (offset == old_gfl[rank + 1]);
  nsend = (int) send_requests->elem_count;

  nleft = nrecv;
  array_of_indices = P4EST_ALLOC (int, nrecv);
  while (nleft > 0) {
    /* finalize some receives */
    mpiret = MPI_Waitsome (nrecv, (MPI_Request *) (recv_requests->array),
                           &outcount, array_of_indices, MPI_STATUSES_IGNORE);
    SC_CHECK_MPI (mpiret);
    P4EST_ASSERT (outcount != MPI_UNDEFINED);
    P4EST_ASSERT (outcount >= 0);
    if (!outcount) {
      continue;
    }
    for (i = 0; i < outcount; i++) {
      p = array_of_indices[i];
      P4EST_ASSERT (p >= 0 && p < mpisize);
      P4EST_ASSERT (p != rank);
      if (data_size) {
        /* get the range in the array of this section */
        local_offset = SC_MAX (0, old_gfl[p] - new_gfl[rank]);
        local_last =
          SC_MIN (new_gfl[rank + 1], old_gfl[p + 1]) - new_gfl[rank];
        for (zz = local_offset; zz < local_last; zz++) {
          /* unpack */
          layer = p2est_quadrant_array_index (new_layers, zz);
          packedlayer = (p2est_quadrant_t *) sc_array_index (recv, zz);
          *layer = *packedlayer;
          layer_data = (void *) (packedlayer + 1);
          layer->p.user_data = sc_mempool_alloc (p6est->user_data_pool);
          memcpy (layer->p.user_data, layer_data, data_size);
        }
      }
    }
    nleft -= outcount;
  }
  P4EST_ASSERT (!nleft);
  P4EST_FREE (array_of_indices);

  /* clean up recv */
  sc_array_destroy (recv);
  sc_array_destroy (recv_requests);

  /* switch to the new layers */
  p6est->layers = new_layers;
  p6est->global_first_layer = new_gfl;
  P4EST_FREE (old_gfl);

  /* wait for sends to complete */
  mpiret = MPI_Waitall (nsend, (MPI_Request *) send->array,
                        MPI_STATUSES_IGNORE);
  SC_CHECK_MPI (mpiret);

  /* clean up send and old_layers */
  sc_array_destroy (send);
  sc_array_destroy (send_requests);
  sc_array_destroy (old_layers);

  return shipped;
}
