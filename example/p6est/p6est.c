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
p6est_connecitivty_new (p4est_connectivity_t *conn4,
                        double *top_to_vertex,
                        double height[3])
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
p6est_connecitivity_destroy (p6est_connectivity_t *conn)
{
  p4est_connectivity_destroy (conn->conn4);
  if (conn->top_to_vertex != NULL) {
    P4EST_FREE (conn->top_to_vertex);
  }
  P4EST_FREE (conn);
}

void
p6est_tree_get_vertices (p6est_connectivity_t * conn,
                         p4est_topidx_t which_tree,
                         double vertices[24])
{
  const double       *btv = conn->conn4->vertices;
  const double       *ttv = conn->top_to_vertex;
  int i, j;

  P4EST_ASSERT (conn->conn4->num_vertices > 0);
  P4EST_ASSERT (btv != NULL);
  P4EST_ASSERT (which_tree >= 0 && which_tree < conn->conn4->num_trees);
  P4EST_ASSERT (vertices != NULL);

  memcpy (vertices, btv,
          which_tree * 12 * sizeof (double));
  if (ttv != NULL) {
    memcpy (vertices + 12, ttv,
            which_tree * 12 * sizeof (double));
  }
  else {
    memcpy (vertices + 12, btv,
            which_tree * 12 * sizeof (double));

    for (i = 0; i < P4EST_HALF; i++) {
      for (j = 0; j < 3; j++) {
        vertices[12 + 3 * i + 0] += conn->height[0];
        vertices[12 + 3 * i + 1] += conn->height[1];
        vertices[12 + 3 * i + 2] += conn->height[2];
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
    double *orig = conn->conn4->vertices;

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
p6est_connectivity_memory_used (p6est_connectivity_t *conn)
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

typedef struct p6est_init_data
{
  int         min_zlevel;
  sc_array_t *layers;
  p6est_init_t init_fn;
  void       *user_pointer;
}
p6est_init_data_t;

void
p6est_layer_init_data (p6est_t * p6est, p4est_topidx_t which_tree,
                       p4est_quadrant_t * column,
                       p2est_quadrant_t * layer,
                       p6est_init_t init_fn)
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
p6est_init_fn (p4est_t *p4est, p4est_topidx_t which_tree,
               p4est_quadrant_t *col)
{
  p6est_t *p6est               = (p6est_t *) p4est->user_pointer;
  p6est_init_data_t *init_data = (p6est_init_data_t *) p6est->user_pointer;
  int nlayers                   = 1 << (init_data->min_zlevel);
  sc_array_t *layers           = init_data->layers;
  size_t incount               = layers->elem_count, zz;
  size_t last                  = incount + nlayers;

  /* we have to sneak the user_pointer in */
  p6est->user_pointer = init_data->user_pointer;

  P6EST_COLUMN_SET_RANGE (col, layers->elem_count, last)

  sc_array_resize (layers, last);

  for (zz = incount; zz < last; zz++) {
    p2est_quadrant_t * layer = p2est_quadrant_array_index (layers, zz);
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
  p6est_t             *p6est = P4EST_ALLOC (p6est_t, 1);
  p4est_t             *p4est;
  p4est_connectivity_t *conn4;
  sc_array_t          *layers;
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

  layers = sc_array_new (sizeof (p2est_quadrant_t));

  if (data_size > 0) {
    user_data_pool = sc_mempool_new (data_size);
  }
  else {
    user_data_pool = NULL;
  }

  quadrant_pool = sc_mempool_new (sizeof (p2est_quadrant_t));

  p6est->mpicomm = mpicomm;
  p6est->mpisize = num_procs;
  p6est->mpirank = rank;
  p6est->data_size = data_size;
  p6est->user_pointer = user_pointer;
  p6est->connectivity = connectivity;
  p6est->layers = layers;
  p6est->user_data_pool = user_data_pool;

  P4EST_ASSERT (min_zlevel <= P4EST_QMAXLEVEL);

  init_data.min_zlevel   = min_zlevel;
  init_data.layers       = layers;
  init_data.init_fn      = init_fn;
  init_data.user_pointer = user_pointer;
  p6est->user_pointer    = &init_data;

  p4est = p4est_new_ext (mpicomm, conn4, min_quadrants / quadpercol,
                         min_level, fill_uniform, 0, p6est_init_fn,
                         (void *) p6est);

  p6est->user_pointer = user_pointer;
  p6est->columns = p4est;
  p6est->global_first_quadrant = P4EST_ALLOC (p4est_gloidx_t, num_procs + 1);
  for (i = 0; i <= num_procs; i++) {
    p6est->global_first_quadrant[i] = quadpercol * p4est->global_first_quadrant[i];
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
p6est_destroy (p6est_t *p6est)
{
  sc_array_t *layers = p6est->layers;
  size_t layercount = layers->elem_count;
  size_t zz;

  for (zz = 0; zz < layercount; zz++) {
    p2est_quadrant_t *layer = p2est_quadrant_array_index (layers, zz);

    p6est_layer_free_data (p6est, layer);
  }
  sc_array_destroy (p6est->layers);

  p4est_destroy (p6est->columns);
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
  size_t zz, qcount = input->layers->elem_count;

  memcpy (p6est, input, sizeof (p4est_t));
  p6est->layers = sc_array_new_size (input->layers->elem_size, input->layers->elem_count);
  sc_array_copy (input->layers, p6est->layers);
  p6est->columns = p4est_copy (input->columns, 0);
  p6est->columns->user_pointer = p6est;
  if (copy_data && p6est->data_size > 0) {
    p6est->user_data_pool = sc_mempool_new (p6est->data_size);
  }
  else {
    p6est->data_size = 0;
  }
  p6est->quadrant_pool = sc_mempool_new (sizeof (p2est_quadrant_t));

  if (p6est->data_size > 0) {
    P4EST_ASSERT (copy_data);
    for (zz = 0; zz < qcount; zz++) {
      p2est_quadrant_t *inlayer = p2est_quadrant_array_index (input->layers, zz);
      p2est_quadrant_t *outlayer = p2est_quadrant_array_index (p6est->layers, zz);

      outlayer->p.user_data = sc_mempool_alloc (p6est->user_data_pool);
      memcpy (outlayer->p.user_data, inlayer->p.user_data, p6est->data_size);
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

  for (jt = p6est->columns->first_local_tree; jt <= p6est->columns->last_local_tree; ++jt) {
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

typedef void (*p6est_replace_t) (p6est_t *p6est,
                                 p4est_topidx_t which_tree,
                                 int num_outcolumns,
                                 int num_outlayers,
                                 p4est_quadrant_t *outcolumns[],
                                 p2est_quadrant_t *outlayers[],
                                 int num_incolumns,
                                 int num_inlayers,
                                 p4est_quadrant_t *incolumns[],
                                 p2est_quadrant_t *inlayers[]);

typedef struct p6est_refine_col_data
{
  p6est_init_t init_fn;
  p6est_replace_t replace_fn;
  void *user_pointer;
}
p6est_refine_col_data_t;

static void
p6est_refine_column (p4est_t *p4est, p4est_topidx_t which_tree,
                     int num_outgoing,
                     p4est_quadrant_t * outgoing[],
                     int num_incoming,
                     p4est_quadrant_t * incoming[])
{
  p6est_t *p6est = (p6est_t *) p4est->user_pointer;
  p6est_refine_col_data_t *refine_col = (p6est_refine_col_data_t *) p6est->user_pointer;
  int nlayers;
  size_t first, last, ifirst, ilast;
  int i, j;
  p2est_quadrant_t *oq, *q;

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
      p6est_layer_init_data (p6est, which_tree, incoming[i], q, refine_col->init_fn);
    }
  }

  if (refine_col->replace_fn != NULL) {
    for (j = 0; j < nlayers; j++) {
      p2est_quadrant_t *inq[P4EST_CHILDREN];

      oq = p2est_quadrant_array_index (p6est->layers, first + j);
      ifirst = incoming[i]->p.piggy3.local_num;
      for (i = 0; i < P4EST_CHILDREN; i++) {
        inq[i] = p2est_quadrant_array_index (p6est->layers, ifirst + j);
      }

      refine_col->replace_fn (p6est, which_tree,
                              1, 1,
                              outgoing, &oq,
                              P4EST_CHILDREN, P4EST_CHILDREN,
                              incoming, inq);
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
p6est_compress_columns (p6est_t *p6est)
{
  size_t              zz, first;
  p4est_topidx_t      jt;
  p4est_quadrant_t   *col;
  p4est_tree_t       *tree;
  sc_array_t         *tquadrants;
  size_t              offset = 0;
  int                 count;

  for (jt = p6est->columns->first_local_tree; jt <= p6est->columns->last_local_tree; ++jt) {
    tree = p4est_tree_array_index (p6est->columns->trees, jt);
    tquadrants = &tree->quadrants;
    for (zz = 0; zz < tquadrants->elem_count; ++zz) {
      col = p4est_quadrant_array_index (tquadrants, zz);
      first = col->p.piggy3.local_num;
      count = col->p.piggy3.which_tree;
      P4EST_ASSERT (first >= offset);
      if (first != offset) {
        memcpy (sc_array_index (p6est->layers, offset),
                sc_array_index (p6est->layers, first),
                count * sizeof (p2est_quadrant_t));
      }
      offset += count;
    }
  }
  sc_array_resize (p6est->layers, offset);
}

#if 0
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
  p6est_quadrant_t   *q, *newq, *qalloc, *qpop;
  p4est_quadrant_t   *c0, *c1, *c2, *c3;
#ifdef P4_TO_P8
  p4est_quadrant_t   *c4, *c5, *c6, *c7;
#endif
  sc_array_t         *quads = p6est->quads;
  sc_array_t         *tquadrants;
  p6est_quadrant_t   *family[2];
  p6est_quadrant_t    parent, *pp = &parent;
  p6est_refine_col_data_t refine_col;
  void               *orig_user_pointer = p6est->user_pointer;
  size_t              zz, first;
  p4est_topidx_t      jt;
  p4est_quadrant_t   *col;
  size_t              offset = 0;
  int                 count;
  sc_array_t         *newcol;
  int                 any_change, this_change;


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
    newcol = sc_array_new (sizeof (p6est_quadrant_t));

    for (jt = p6est->p4est->first_local_tree; jt <= p6est->p4est->last_local_tree; ++jt) {
      tree = p4est_tree_array_index (p6est->p4est->trees, jt);
      tquadrants = &tree->quadrants;

      for (zz = 0; zz < tquadrants->elem_count; ++zz) {
        col = p4est_quadrant_array_index (tquadrants, zz);
        first = (size_t) col->p.piggy3.local_num;
        count = (int) col->p.piggy3.which_tree;

        any_change = 0;
        parent = NULL;

        for (current = first; current < first + count; current++) {
          q = p6est_quadrant_array_index (quads, current);
          parent = q;
          newq = sc_array_push (newcol);
          *newq = *q;
          stop_recurse = 0;
          this_change = 0;
          for (;;) {
            if (!stop_recurse && zrefine_fn (p4est, nt, newq) && (int) newq->zlevel < allowed_zlevel) {
              this_change = 1;
              any_change = 1;
              newq->zlevel++;
              stop_recurse = !refine_recursive;
            }
            else {
              q = newq;
              newq = sc_array_push (newcol);
              *newq = *q;
              newq->z += P6EST_QUADRANT_LEN (newq->zlevel);
              if (current + 1 >= first + count ||
                  (p6est_quadrant_array_index (quads, current + 1))->z <= newq->z) {
                /* pop */
                P4EST_ASSERT (newcol->elem_count > 1);
                newcol->elem_count--;
                break;
              }
              while (newq->zlevel > 0 && !(newq->z & P6EST_QUADRANT_LEN (newq->zlevel))) {
                newq->zlevel--;
              }
            }
          }
          if (this_change) {
          }
        }
        if (any_change) {
          old_count = quads->elem_count;
          sc_array_resize (quads, old_count + newcol->elem_count);
          memcpy (sc_array_index (quads, old_count),
                  sc_array_index (newcol, 0), newcol->elem_size * newcol->elem_count);
          col->p.piggy3.local_num = (p4est_locidx_t) old_count;
          col->p.piggy3.which_tree = (p4est_topidx_t) newcol->elem_count;
        }
        sc_array_truncate (newcol);
      }
    }

    if (current != quads->elem_count) {
      P4EST_ASSERT (q != NULL);

      /* now we have a quadrant to refine, prepend it to the list */
      qalloc = p6est_quadrant_mempool_alloc (p4est->quadrant_pool);
      *qalloc = *q;               /* never prepend array members directly */
      sc_list_prepend (list, qalloc);     /* only newly allocated quadrants */

      P6EST_QUADRANT_INIT (&parent);

      /*
         current points to the next array member to write
         restpos points to the next array member to read
         */
      restpos = current + 1;
    }
  }
  /* update counts */
}
#endif
