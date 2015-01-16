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
#include <p6est_ghost.h>
#include <p6est_profile.h>
#include <p4est_lnodes.h>
#include <p8est.h>
#include <p4est_extended.h>
#include <sc_containers.h>
#include <p4est_communication.h>
#include <sc_io.h>
#include <p6est_extended.h>

/* htonl is in either of these two */
#ifdef P4EST_HAVE_ARPA_NET_H
#include <arpa/inet.h>
#endif
#ifdef P4EST_HAVE_NETINET_IN_H
#include <netinet/in.h>
#endif
#ifdef P4EST_HAVE_ZLIB
#include <zlib.h>
#endif

p6est_connectivity_t *
p6est_connectivity_new (p4est_connectivity_t * conn4,
                        double *top_vertices, double height[3])
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

  if (top_vertices != NULL) {
    conn->top_vertices = P4EST_ALLOC (double, 3 * conn4->num_vertices);
    memcpy (conn->top_vertices, top_vertices,
            3 * conn4->num_vertices * sizeof (double));
  }
  else {
    conn->top_vertices = NULL;
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
  if (conn->top_vertices != NULL) {
    P4EST_FREE (conn->top_vertices);
  }
  P4EST_FREE (conn);
}

void
p6est_tree_get_vertices (p6est_connectivity_t * conn,
                         p4est_topidx_t which_tree, double vertices[24])
{
  p4est_connectivity_t *conn4 = conn->conn4;
  const double       *btv = conn4->vertices;
  const double       *ttv = conn->top_vertices;
  int                 i, j, k;
  double             *height = conn->height;
  double              zerooff[3] = { 0., 0., 0. };

  P4EST_ASSERT (conn4->num_vertices > 0);
  P4EST_ASSERT (btv != NULL);
  P4EST_ASSERT (which_tree >= 0 && which_tree < conn4->num_trees);
  P4EST_ASSERT (vertices != NULL);

  for (k = 0; k < 2; k++) {
    const double       *verts = k ? ttv : btv;
    const double       *off = (verts != NULL) ? zerooff : height;

    if (verts == NULL) {
      verts = btv;
    }

    for (i = 0; i < P4EST_CHILDREN; i++) {
      p4est_topidx_t      vert;

      vert = conn4->tree_to_vertex[P4EST_CHILDREN * which_tree + i];
      for (j = 0; j < 3; j++) {
        vertices[3 * (P4EST_CHILDREN * k + i) + j] = verts[3 * vert + j] +
          off[j];
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
  if (conn->top_vertices != NULL) {
    double             *orig = conn->conn4->vertices;

    conn->conn4->vertices = conn->top_vertices;
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
    conn->top_vertices == NULL ? 0 :
    (conn->conn4->num_vertices * 3 * sizeof (double));
}

int
p6est_connectivity_extra_sink (p6est_connectivity_t * conn,
                               sc_io_sink_t * sink)
{
  int                 retval;
  size_t              u64z, tcount;
  double             *v;
  uint64_t            num_vertices;

  u64z = sizeof (uint64_t);
  num_vertices = (conn->top_vertices != NULL) ?
    (uint64_t) conn->conn4->num_vertices : 0;

  retval = sc_io_sink_write (sink, &num_vertices, 1 * u64z);

  if (conn->top_vertices != NULL) {
    tcount = (size_t) (3 * conn->conn4->num_vertices);
    v = conn->top_vertices;
  }
  else {
    tcount = 3;
    v = conn->height;
  }
  retval = retval || sc_io_sink_write (sink, v, tcount * sizeof (double));

  return retval;
}

p6est_connectivity_t *
p6est_connectivity_extra_source (p4est_connectivity_t * conn4,
                                 sc_io_source_t * source)
{
  p6est_connectivity_t *conn;
  int                 retval;
  size_t              u64z;
  uint64_t            num_vertices;
  double             *top_vertices;
  double              height[3];

  u64z = sizeof (uint64_t);

  retval = sc_io_source_read (source, &num_vertices, 1 * u64z, NULL);
  if (retval) {
    return NULL;
  }
  if (num_vertices > 0) {
    height[0] = height[1] = height[2] = 0.;

    if (num_vertices != (uint64_t) conn4->num_vertices) {
      return NULL;
    }
    top_vertices = P4EST_ALLOC (double, 3 * num_vertices);
    retval = sc_io_source_read (source, top_vertices,
                                3 * num_vertices * sizeof (double), NULL);

    if (retval) {
      P4EST_FREE (top_vertices);
      return NULL;
    }
  }
  else {
    top_vertices = NULL;
    retval = sc_io_source_read (source, height, 3 * sizeof (double), NULL);
    if (retval) {
      return NULL;
    }
  }
  conn = P4EST_ALLOC (p6est_connectivity_t, 1);

  conn->conn4 = conn4;
  conn->top_vertices = top_vertices;
  conn->height[0] = height[0];
  conn->height[1] = height[1];
  conn->height[2] = height[2];

  return conn;
}

int
p6est_connectivity_sink (p6est_connectivity_t * conn, sc_io_sink_t * sink)
{
  int                 retval;

  retval = p4est_connectivity_sink (conn->conn4, sink);

  retval = retval || p6est_connectivity_extra_sink (conn, sink);

  return retval;
}

p6est_connectivity_t *
p6est_connectivity_source (sc_io_source_t * source)
{
  p4est_connectivity_t *conn4;

  conn4 = p4est_connectivity_source (source);
  if (conn4 == NULL) {
    return NULL;
  }

  return p6est_connectivity_extra_source (conn4, source);
}

int
p6est_connectivity_save (const char *filename, p6est_connectivity_t * conn)
{
  int                 retval;
  sc_io_sink_t       *sink;

  sink = sc_io_sink_new (SC_IO_TYPE_FILENAME, SC_IO_MODE_WRITE,
                         SC_IO_ENCODE_NONE, filename);
  if (sink == NULL) {
    return -1;
  }

  /* Close file even on earlier write error */
  retval = p6est_connectivity_sink (conn, sink);
  retval = sc_io_sink_destroy (sink) || retval;

  return retval;
}

p6est_connectivity_t *
p6est_connectivity_load (const char *filename, size_t * bytes)
{
  int                 retval;
  size_t              bytes_in;
  sc_io_source_t     *source;
  p6est_connectivity_t *conn;

  source = sc_io_source_new (SC_IO_TYPE_FILENAME,
                             SC_IO_ENCODE_NONE, filename);
  if (source == NULL) {
    return NULL;
  }

  /* Get byte length and close file even on earlier read error */
  conn = p6est_connectivity_source (source);
  retval = sc_io_source_complete (source, &bytes_in, NULL) || conn == NULL;
  retval = sc_io_source_destroy (source) || retval;
  if (retval) {
    if (conn != NULL) {
      p6est_connectivity_destroy (conn);
    }
    return NULL;
  }

  if (bytes != NULL) {
    *bytes = bytes_in;
  }
  return conn;
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

typedef struct p6est_init_data
{
  int                 min_zlevel;
  int                 num_zroot;
  sc_array_t         *layers;
  p6est_init_t        init_fn;
  void               *user_pointer;
}
p6est_init_data_t;

static void
p6est_init_fn (p4est_t * p4est, p4est_topidx_t which_tree,
               p4est_quadrant_t * col)
{
  p6est_t            *p6est = (p6est_t *) p4est->user_pointer;
  p6est_init_data_t  *init_data = (p6est_init_data_t *) p6est->user_pointer;
  int                 log_zroot = SC_LOG2_32 (init_data->num_zroot - 1) + 1;
  int                 nlayers =
    (1 << (init_data->min_zlevel - log_zroot)) * init_data->num_zroot;
  sc_array_t         *layers = init_data->layers;
  size_t              incount = layers->elem_count, zz;
  size_t              last = incount + nlayers;
  p2est_quadrant_t   *layer;

  /* we have to sneak the user_pointer in */
  p6est->user_pointer = init_data->user_pointer;

  P6EST_COLUMN_SET_RANGE (col, layers->elem_count, last);

  layer = (p2est_quadrant_t *) sc_array_push_count (layers, nlayers);

  for (zz = incount; zz < last; zz++, layer++) {
    P2EST_QUADRANT_INIT (layer);

    layer->level = init_data->min_zlevel;
    layer->z = (zz - incount) * P4EST_QUADRANT_LEN (layer->level);

    p6est_layer_init_data (p6est, which_tree, col, layer, init_data->init_fn);
  }

  /* we have to sneak the user_pointer out */
  p6est->user_pointer = (void *) init_data;
}

p6est_t            *
p6est_new_ext (sc_MPI_Comm mpicomm, p6est_connectivity_t * connectivity,
               p4est_locidx_t min_quadrants, int min_level, int min_zlevel,
               int num_zroot,
               int fill_uniform, size_t data_size, p6est_init_t init_fn,
               void *user_pointer)
{
  p6est_t            *p6est = P4EST_ALLOC (p6est_t, 1);
  p4est_t            *p4est;
  sc_array_t         *layers;
  sc_mempool_t       *user_data_pool;
  p6est_init_data_t   init_data;
  int                 mpiret, num_procs, rank;
  int                 log_zroot = SC_LOG2_32 (num_zroot - 1) + 1;
  int                 quadpercol =
    (1 << (min_zlevel - log_zroot)) * num_zroot;
  int                 i;

  P4EST_GLOBAL_PRODUCTIONF
    ("Into p6est_new with min layers %lld z-level %d\n",
     (long long) min_quadrants, SC_MAX (min_zlevel, 0));
  p4est_log_indent_push ();

  mpiret = sc_MPI_Comm_size (mpicomm, &num_procs);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (mpicomm, &rank);
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
  p6est->root_len = num_zroot * P4EST_QUADRANT_LEN (log_zroot);

  P4EST_ASSERT (min_zlevel <= P4EST_QMAXLEVEL);
  P4EST_ASSERT (num_zroot >= 1);
  P4EST_ASSERT (min_zlevel >= log_zroot);

  init_data.min_zlevel = min_zlevel;
  init_data.num_zroot = num_zroot;
  init_data.layers = layers;
  init_data.init_fn = init_fn;
  init_data.user_pointer = user_pointer;
  p6est->user_pointer = &init_data;

  P4EST_GLOBAL_VERBOSE ("Creating p4est for columns\n");
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

  /* print more statistics */
  P4EST_VERBOSEF ("total local layers %lld\n",
                  (long long) layers->elem_count);

  p4est_log_indent_pop ();
  P4EST_GLOBAL_PRODUCTIONF
    ("Done p6est_new with %lld total layers in %lld total columns\n",
     (long long) p6est->global_first_layer[p6est->mpisize],
     (long long) p6est->columns->global_num_quadrants);

  return p6est;
}

p6est_t            *
p6est_new (sc_MPI_Comm mpicomm, p6est_connectivity_t * connectivity,
           size_t data_size, p6est_init_t init_fn, void *user_pointer)
{
  return p6est_new_ext (mpicomm, connectivity, 0, 0, 0, 1, 1,
                        data_size, init_fn, user_pointer);
}

p6est_t            *
p6est_new_from_p4est (p4est_t * p4est, double *top_vertices, double height[3],
                      int min_zlevel, size_t data_size, p6est_init_t init_fn,
                      void *user_pointer)
{
  p6est_t            *p6est = P4EST_ALLOC (p6est_t, 1);
  p6est_connectivity_t *conn;
  sc_array_t         *layers;
  sc_mempool_t       *user_data_pool;
  p6est_init_data_t   init_data;
  int                 num_procs;
  int                 quadpercol = (1 << min_zlevel);
  int                 i;

  P4EST_GLOBAL_PRODUCTIONF
    ("Into p6est_new_from_p4est with z-level %d\n", SC_MAX (min_zlevel, 0));
  p4est_log_indent_push ();

  layers = sc_array_new (sizeof (p2est_quadrant_t));

  if (data_size > 0) {
    user_data_pool = sc_mempool_new (data_size);
  }
  else {
    user_data_pool = NULL;
  }

  conn = p6est_connectivity_new (p4est->connectivity, top_vertices, height);

  p6est->layer_pool = sc_mempool_new (sizeof (p2est_quadrant_t));

  p6est->mpicomm = p4est->mpicomm;
  p6est->mpisize = num_procs = p4est->mpisize;
  p6est->mpirank = p4est->mpirank;
  p6est->data_size = data_size;
  p6est->user_pointer = user_pointer;
  p6est->connectivity = conn;
  p6est->layers = layers;
  p6est->user_data_pool = user_data_pool;
  p6est->columns = p4est_copy (p4est, 0);
  p6est->columns->connectivity = conn->conn4;

  P4EST_ASSERT (min_zlevel <= P4EST_QMAXLEVEL);

  init_data.min_zlevel = min_zlevel;
  init_data.layers = layers;
  init_data.init_fn = init_fn;
  init_data.user_pointer = user_pointer;
  p6est->user_pointer = &init_data;

  p4est_reset_data (p6est->columns, 0, p6est_init_fn, (void *) p6est);

  p6est->user_pointer = user_pointer;
  p6est->global_first_layer = P4EST_ALLOC (p4est_gloidx_t, num_procs + 1);
  for (i = 0; i <= num_procs; i++) {
    p6est->global_first_layer[i] =
      quadpercol * p4est->global_first_quadrant[i];
  }

  /* print more statistics */
  P4EST_VERBOSEF ("total local layers %lld\n",
                  (long long) layers->elem_count);

  p4est_log_indent_pop ();
  P4EST_GLOBAL_PRODUCTIONF
    ("Done p6est_new with %lld total layers in %lld total columns\n",
     (long long) p6est->global_first_layer[p6est->mpisize],
     (long long) p6est->columns->global_num_quadrants);

  return p6est;
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

void
p6est_save (const char *filename, p6est_t * p6est, int save_data)
{
  p6est_save_ext (filename, p6est, save_data, 1);
}

void
p6est_save_ext (const char *filename, p6est_t * p6est,
                int save_data, int save_partition)
{
  FILE               *file;
  int                 rank = p6est->mpirank;
  long                fpos = -1, foffset;
  int                 align = 32;
  uint64_t            u64a;
  sc_io_sink_t       *sink;
  int                 retval, mpiret;
  p4est_t            *savecolumns;
#ifdef P4EST_MPIIO_WRITE
  MPI_File            mpifile;
  MPI_Offset          mpipos;
  MPI_Offset          mpithis;
#else
  long                fthis;
  sc_MPI_Status       mpistatus;
#endif
  int                 num_procs = p6est->mpisize;

  P4EST_GLOBAL_PRODUCTION ("Into p6est_save\n");
  p4est_log_indent_push ();

  savecolumns = p4est_copy (p6est->columns, 0);
  p4est_reset_data (savecolumns, 2 * sizeof (p4est_locidx_t), NULL, NULL);
  size_t              comb_size, data_size = p6est->data_size;
  size_t              zz, nlayers = p6est->layers->elem_count;
  char               *lbuf, *bp;

  if (!data_size) {
    save_data = 0;
  }
  if (!save_data) {
    data_size = 0;
  }

  comb_size = 2 * sizeof (p4est_qcoord_t) + data_size;

  {
    p4est_topidx_t      jt;
    p4est_tree_t       *tree, *savetree;
    sc_array_t         *tquadrants, *savetquadrants;
    p4est_quadrant_t   *col, *savecol;
    size_t              first, last;
    p4est_locidx_t      lfirst, llast;

    for (jt = p6est->columns->first_local_tree;
         jt <= p6est->columns->last_local_tree; ++jt) {
      tree = p4est_tree_array_index (p6est->columns->trees, jt);
      savetree = p4est_tree_array_index (savecolumns->trees, jt);
      tquadrants = &tree->quadrants;
      savetquadrants = &savetree->quadrants;
      for (zz = 0; zz < tquadrants->elem_count; ++zz) {
        p4est_locidx_t     *savedata;
        col = p4est_quadrant_array_index (tquadrants, zz);
        savecol = p4est_quadrant_array_index (savetquadrants, zz);
        P6EST_COLUMN_GET_RANGE (col, &first, &last);
        lfirst = (p4est_locidx_t) first;
        llast = (p4est_locidx_t) last;
        savedata = (p4est_locidx_t *) savecol->p.user_data;
        savedata[0] = lfirst;
        savedata[1] = llast;
      }
    }
  }

  /* save the columns */
  p4est_save_ext (filename, savecolumns, 1, save_partition);

  p4est_destroy (savecolumns);

  /* we need a barrier so that all files have finished writing in
   * p4est_save_ext before we start writing additional data to the file */
  mpiret = sc_MPI_Barrier (p6est->mpicomm);
  SC_CHECK_MPI (mpiret);

  if (rank == 0) {
    file = fopen (filename, "ab");
    SC_CHECK_ABORT (file != NULL, "file open");

    /* align */
    fpos = ftell (file);
    SC_CHECK_ABORT (fpos > 0, "first file tell");
    while (fpos % align != 0) {
      retval = fputc ('\0', file);
      SC_CHECK_ABORT (retval == 0, "first file align");
      ++fpos;
    }

    /* write the p6est_connectivity extra data */
    sink = sc_io_sink_new (SC_IO_TYPE_FILEFILE, SC_IO_MODE_APPEND,
                           SC_IO_ENCODE_NONE, file);
    SC_CHECK_ABORT (sink != NULL, "file sink");
    retval = p6est_connectivity_extra_sink (p6est->connectivity, sink);
    SC_CHECK_ABORT (retval == 0, "sink connectivity");
    retval = sc_io_sink_destroy (sink);
    SC_CHECK_ABORT (retval == 0, "destroy sink");

    fpos = ftell (file);
    SC_CHECK_ABORT (fpos > 0, "second file tell");
    while (fpos % align != 0) {
      retval = fputc ('\0', file);
      SC_CHECK_ABORT (retval == 0, "second file align");
      ++fpos;
    }

    /* write the data size */
    u64a = (uint64_t) data_size;
    sc_fwrite (&u64a, sizeof (uint64_t), 1, file, "write data size");

    fpos = ftell (file);
    SC_CHECK_ABORT (fpos > 0, "third file tell");
    while (fpos % align != 0) {
      retval = fputc ('\0', file);
      SC_CHECK_ABORT (retval == 0, "third file align");
      ++fpos;
    }

#ifdef P4EST_MPIIO_WRITE
    /* We will close the sequential access to the file */
    /* best attempt to flush file to disk */
    retval = fflush (file);
    SC_CHECK_ABORT (retval == 0, "file flush");
#ifdef P4EST_HAVE_FSYNC
    retval = fsync (fileno (file));
    SC_CHECK_ABORT (retval == 0, "file fsync");
#endif
    retval = fclose (file);
    SC_CHECK_ABORT (retval == 0, "file close");
    file = NULL;
#else
    /* file is still open for sequential write mode */
#endif
  }
  else {
    file = NULL;
  }

#ifndef P4EST_MPIIO_WRITE
  if (rank > 0) {
    /* wait for sequential synchronization */
    mpiret = sc_MPI_Recv (&fpos, 1, sc_MPI_LONG, rank - 1, P4EST_COMM_SAVE,
                          p6est->mpicomm, &mpistatus);
    SC_CHECK_MPI (mpiret);

    /* open file after all previous processors have written to it */
    file = fopen (filename, "rb+");
    SC_CHECK_ABORT (file != NULL, "file open");
  }
#else
  /* Every core opens the file in append mode */
  mpiret = MPI_File_open (p6est->mpicomm, (char *) filename,
                          MPI_MODE_WRONLY | MPI_MODE_APPEND |
                          MPI_MODE_UNIQUE_OPEN, MPI_INFO_NULL, &mpifile);
  SC_CHECK_MPI (mpiret);
  mpiret = MPI_File_get_position (mpifile, &mpipos);
  SC_CHECK_MPI (mpiret);
#endif

  if (rank > 0) {
    /* seek to the beginning of this processor's storage */
    foffset = (long) (p6est->global_first_layer[rank] * comb_size);

#ifndef P4EST_MPIIO_WRITE
    fthis = fpos + foffset;
    retval = fseek (file, fthis, SEEK_SET);
    SC_CHECK_ABORT (retval == 0, "seek data");
#else
    mpithis = mpipos + (MPI_Offset) foffset;
    mpiret = MPI_File_seek (mpifile, mpithis, MPI_SEEK_SET);
    SC_CHECK_MPI (mpiret);
#endif
  }

  /* write layers and data interleaved */
  bp = lbuf = P4EST_ALLOC (char, comb_size * nlayers);
  for (zz = 0; zz < nlayers; zz++) {
    p2est_quadrant_t   *layer;
    p4est_qcoord_t     *qpos;

    qpos = (p4est_locidx_t *) bp;

    layer = p2est_quadrant_array_index (p6est->layers, zz);
    *qpos++ = layer->z;
    *qpos++ = (p4est_qcoord_t) layer->level;
    if (save_data) {
      memcpy (qpos, layer->p.user_data, data_size);
    }
    bp += comb_size;
  }
#ifndef P4EST_MPIIO_WRITE
  sc_fwrite (lbuf, comb_size, nlayers, file, "write quadrants");
#else
  sc_mpi_write (mpifile, lbuf, comb_size * nlayers, MPI_BYTE,
                "write quadrants");
#endif
  P4EST_FREE (lbuf);

#ifndef P4EST_MPIIO_WRITE
  /* best attempt to flush file to disk */
  retval = fflush (file);
  SC_CHECK_ABORT (retval == 0, "file flush");
#ifdef P4EST_HAVE_FSYNC
  retval = fsync (fileno (file));
  SC_CHECK_ABORT (retval == 0, "file fsync");
#endif
  retval = fclose (file);
  SC_CHECK_ABORT (retval == 0, "file close");
  file = NULL;

  /* initiate sequential synchronization */
  if (rank < num_procs - 1) {
    mpiret = sc_MPI_Send (&fpos, 1, sc_MPI_LONG, rank + 1, P4EST_COMM_SAVE,
                          p6est->mpicomm);
    SC_CHECK_MPI (mpiret);
  }
#else
  mpiret = MPI_File_close (&mpifile);
  SC_CHECK_MPI (mpiret);
#endif

  p4est_log_indent_pop ();
  P4EST_GLOBAL_PRODUCTION ("Done p6est_save\n");
}

p6est_t            *
p6est_load (const char *filename, sc_MPI_Comm mpicomm, size_t data_size,
            int load_data, void *user_pointer,
            p6est_connectivity_t ** connectivity)
{
  return p6est_load_ext (filename, mpicomm, data_size, load_data,
                         0, 0, user_pointer, connectivity);
}

p6est_t            *
p6est_load_ext (const char *filename, sc_MPI_Comm mpicomm, size_t data_size,
                int load_data, int autopartition, int broadcasthead,
                void *user_pointer, p6est_connectivity_t ** connectivity)
{
  int                 align = 32;
  int                 retval;
  sc_io_source_t     *src;
  p4est_t            *loadcolumns, *columns;
  p4est_connectivity_t *conn4;
  p6est_connectivity_t *conn;
  p4est_topidx_t      jt;
  p4est_tree_t       *tree, *loadtree;
  p4est_quadrant_t   *col, *loadcol;
  sc_array_t         *tquadrants, *loadtquadrants;
  size_t              zz, first, last, zcount, zpadding;
  p4est_locidx_t     *loaddata, lfirst, llast, nlayers;
  p6est_t            *p6est;
  size_t              save_data_size, comb_size;
  uint64_t            u64a;
  p4est_gloidx_t     *gfl;
  int                 rank, mpisize;

  P4EST_GLOBAL_PRODUCTIONF ("Into p6est_load %s\n", filename);
  p4est_log_indent_push ();

  /* open file on all processors */
  src = sc_io_source_new (SC_IO_TYPE_FILENAME, SC_IO_ENCODE_NONE, filename);
  SC_CHECK_ABORT (src != NULL, "file source");

  loadcolumns = p4est_source_ext (src, mpicomm, 2 * sizeof (p4est_locidx_t),
                                  1, autopartition, broadcasthead, NULL,
                                  &conn4);
  columns = p4est_copy (loadcolumns, 0);

  nlayers = 0;
  for (jt = columns->first_local_tree; jt <= columns->last_local_tree; ++jt) {
    tree = p4est_tree_array_index (columns->trees, jt);
    loadtree = p4est_tree_array_index (loadcolumns->trees, jt);
    tquadrants = &tree->quadrants;
    loadtquadrants = &loadtree->quadrants;
    for (zz = 0; zz < tquadrants->elem_count; ++zz) {
      col = p4est_quadrant_array_index (tquadrants, zz);
      loadcol = p4est_quadrant_array_index (loadtquadrants, zz);
      loaddata = (p4est_locidx_t *) loadcol->p.user_data;
      lfirst = loaddata[0];
      llast = loaddata[1];
      first = (size_t) nlayers;
      nlayers += llast - lfirst;
      last = first + (size_t) (llast - lfirst);

      P6EST_COLUMN_SET_RANGE (col, first, last);
    }
  }
  columns->connectivity = conn4;
  p4est_destroy (loadcolumns);

  /* padding */
  zcount = src->bytes_out;
  zpadding = (align - zcount % align) % align;
  retval = sc_io_source_read (src, NULL, zpadding, NULL);
  SC_CHECK_ABORT (!retval, "source padding");

  conn = p6est_connectivity_extra_source (columns->connectivity, src);
  if (connectivity != NULL) {
    *connectivity = conn;
  }

  /* padding */
  zcount = src->bytes_out;
  zpadding = (align - zcount % align) % align;
  retval = sc_io_source_read (src, NULL, zpadding, NULL);
  SC_CHECK_ABORT (!retval, "source padding");

  /* read data size */
  retval = sc_io_source_read (src, &u64a, sizeof (uint64_t), NULL);
  SC_CHECK_ABORT (!retval, "source data size");
  save_data_size = (size_t) u64a;

  if (load_data) {
    SC_CHECK_ABORT (save_data_size == data_size, "data size mismatch");
  }

  /* padding */
  zcount = src->bytes_out;
  zpadding = (align - zcount % align) % align;
  retval = sc_io_source_read (src, NULL, zpadding, NULL);
  SC_CHECK_ABORT (!retval, "source padding");

  p6est = P4EST_ALLOC (p6est_t, 1);
  columns->user_pointer = p6est;

  p6est->columns = columns;
  p6est->connectivity = conn;
  p6est->data_size = data_size;
  p6est->mpicomm = mpicomm;
  p6est->mpisize = mpisize = columns->mpisize;
  p6est->mpirank = rank = columns->mpirank;
  p6est->global_first_layer = gfl = P4EST_ALLOC (p4est_gloidx_t,
                                                 p6est->mpisize + 1);
  p6est->layers =
    sc_array_new_size (sizeof (p2est_quadrant_t), (size_t) nlayers);
  p6est->layer_pool = sc_mempool_new (sizeof (p2est_quadrant_t));
  p6est->user_pointer = user_pointer;
  p6est->user_data_pool = data_size ? sc_mempool_new (data_size) : NULL;

  p6est_update_offsets (p6est);

  comb_size = 2 * sizeof (p4est_qcoord_t) + save_data_size;

  if (gfl[rank]) {
    /* seek to the start of my layers */
    retval = sc_io_source_read (src, NULL, (long) (gfl[rank] * comb_size),
                                NULL);
    SC_CHECK_ABORT (!retval, "seek to layers");
  }

  if (nlayers) {
    sc_array_t         *loaddata;

    loaddata = sc_array_new_size (comb_size, (size_t) nlayers);

    retval =
      sc_io_source_read (src, loaddata->array, comb_size * nlayers, NULL);
    SC_CHECK_ABORT (!retval, "read layers");

    for (zz = 0; zz < (size_t) nlayers; zz++) {
      p2est_quadrant_t   *layer =
        p2est_quadrant_array_index (p6est->layers, zz);
      p4est_qcoord_t     *loadlayer =
        (p4est_qcoord_t *) sc_array_index (loaddata, zz);

      P2EST_QUADRANT_INIT (layer);
      layer->z = loadlayer[0];
      layer->level = (int8_t) loadlayer[1];
      if (load_data) {
        layer->p.user_data = sc_mempool_alloc (p6est->user_data_pool);
        memcpy (layer->p.user_data, &loadlayer[2], save_data_size);
      }
    }
    sc_array_destroy (loaddata);
  }

  retval = sc_io_source_destroy (src);
  SC_CHECK_ABORT (!retval, "destroy source");

  p4est_log_indent_pop ();
  P4EST_GLOBAL_PRODUCTIONF ("Done p6est_load %s\n", filename);

  return p6est;
}

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
    q = (p2est_quadrant_t *) sc_array_push_count (p6est->layers, nlayers);
    oq = (p2est_quadrant_t *) sc_array_index (p6est->layers, first);
    P6EST_COLUMN_SET_RANGE (incoming[i], ifirst, ilast);
    for (j = 0; j < nlayers; j++, oq++, q++) {
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
      P4EST_ASSERT (count > 0);
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

  mpiret = sc_MPI_Allgather (&mycount, 1, P4EST_MPI_GLOIDX, gfl, 1,
                             P4EST_MPI_GLOIDX, p6est->mpicomm);
  SC_CHECK_MPI (mpiret);

  for (p = 0; p < p6est->mpisize; p++) {
    thiscount = gfl[p];
    gfl[p] = psum;
    psum += thiscount;
  }
  gfl[p6est->mpisize] = psum;
  P4EST_ASSERT ((size_t) (gfl[p6est->mpirank + 1] - gfl[p6est->mpirank]) ==
                p6est->layers->elem_count);
}

void
p6est_refine_columns_ext (p6est_t * p6est, int refine_recursive,
                          int allowed_level, p6est_refine_column_t refine_fn,
                          p6est_init_t init_fn, p6est_replace_t replace_fn)
{
  p6est_refine_col_data_t refine_col;
  void               *orig_user_pointer = p6est->user_pointer;

  P4EST_GLOBAL_PRODUCTIONF ("Into p6est_refine_columns with %lld total layers"
                            " in %lld total columns\n",
                            (long long) p6est->global_first_layer[p6est->
                                                                  mpisize],
                            (long long) p6est->columns->global_num_quadrants);
  p4est_log_indent_push ();
  refine_col.refine_col_fn = refine_fn;
  refine_col.init_fn = init_fn;
  refine_col.replace_fn = replace_fn;
  refine_col.user_pointer = orig_user_pointer;

  p6est->user_pointer = (void *) &refine_col;
  P4EST_GLOBAL_VERBOSE ("Refining p4est for columns\n");
  p4est_refine_ext (p6est->columns, refine_recursive, allowed_level,
                    p6est_refine_column_int, NULL,
                    p6est_replace_column_split);
  p6est->user_pointer = orig_user_pointer;

  p6est_compress_columns (p6est);
  p6est_update_offsets (p6est);
  p4est_log_indent_pop ();
  P4EST_GLOBAL_PRODUCTIONF ("Done p6est_refine_columns with %lld total layers"
                            " in %lld total columns\n",
                            (long long) p6est->global_first_layer[p6est->
                                                                  mpisize],
                            (long long) p6est->columns->global_num_quadrants);
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

  P4EST_GLOBAL_PRODUCTIONF ("Into p6est_refine_layers with %lld total layers"
                            " in %lld total columns, allowed level %d\n",
                            (long long) p6est->global_first_layer[p6est->
                                                                  mpisize],
                            (long long) p6est->columns->global_num_quadrants,
                            allowed_level);
  p4est_log_indent_push ();
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
        newq = (p2est_quadrant_t *) sc_array_push_count (layers,
                                                         newcol->elem_count);
        memcpy (newq, sc_array_index (newcol, 0),
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
  p4est_log_indent_pop ();
  P4EST_GLOBAL_PRODUCTIONF ("Done p6est_refine_layers with %lld total layers "
                            " in %lld total columns\n",
                            (long long) p6est->global_first_layer[p6est->
                                                                  mpisize],
                            (long long) p6est->columns->global_num_quadrants);
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
      /* stop if a) we are as big as the root ancestor, or
       * b) if we have reached the end of the root domain */
      if (q->level > ancestor_level &&
          ((q->z + P4EST_QUADRANT_LEN (q->level)) < p6est->root_len)) {
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
  P4EST_ASSERT ((endpos - startpos == P4EST_QUADRANT_LEN (ancestor_level)) ||
                (startpos == 0 && endpos == p6est->root_len));

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
  p = (p2est_quadrant_t *) sc_array_push_count (layers, new_count);
  memcpy (p, sc_array_index (work_array, 0),
          new_count * work_array->elem_size);

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

  P4EST_GLOBAL_PRODUCTIONF
    ("Into p6est_coarsen_columns with %lld total layers"
     " in %lld total columns\n",
     (long long) p6est->global_first_layer[p6est->mpisize],
     (long long) p6est->columns->global_num_quadrants);
  p4est_log_indent_push ();

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
  p4est_log_indent_pop ();
  P4EST_GLOBAL_PRODUCTIONF
    ("Done p6est_coarsen_columns with %lld total layers "
     " in %lld total columns\n",
     (long long) p6est->global_first_layer[p6est->mpisize],
     (long long) p6est->columns->global_num_quadrants);
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

  P4EST_GLOBAL_PRODUCTIONF ("Into p6est_coarsen_layers with %lld total layers"
                            " in %lld total columns\n",
                            (long long) p6est->global_first_layer[p6est->
                                                                  mpisize],
                            (long long) p6est->columns->global_num_quadrants);
  p4est_log_indent_push ();

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

  p4est_log_indent_pop ();
  P4EST_GLOBAL_PRODUCTIONF
    ("Done p6est_coarsen_layers with %lld total layers "
     " in %lld total columns\n",
     (long long) p6est->global_first_layer[p6est->mpisize],
     (long long) p6est->columns->global_num_quadrants);
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
gloidx_compare_overlap (const void *key, const void *array)
{
  const p4est_gloidx_t search = *((const p4est_gloidx_t *) key);
  const p4est_gloidx_t *range = (const p4est_gloidx_t *) array;

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
  sc_MPI_Request     *req;
  void               *layer_data;
  p2est_quadrant_t   *layer;
  p2est_quadrant_t   *packedlayer;
  int                 i, nrecv, nsend, nleft;
  int                 outcount, *array_of_indices;
  size_t              self_offset, self_count;
  sc_array_t         *new_layers;
  p4est_gloidx_t      local_offset, local_last;
  sc_array_t         *old_layers = p6est->layers;
  sc_array_t         *recv_procs;
  p4est_gloidx_t      shipped;
  int                *ip;
  void               *orig_user_pointer = p6est->user_pointer;

  P4EST_GLOBAL_PRODUCTIONF ("Into p6est_parition with %lld total layers"
                            " in %lld total columns\n",
                            (long long) p6est->global_first_layer[p6est->
                                                                  mpisize],
                            (long long) p6est->columns->global_num_quadrants);
  p4est_log_indent_push ();
  /* wrap the p6est_weight_t in a p4est_weight_t */
  wc.layer_weight_fn = weight_fn;
  wc.user_pointer = orig_user_pointer;
  p6est->user_pointer = &wc;
  p6est_compress_columns (p6est);
  /* repartition the columns */
  p4est_partition_ext (p6est->columns, partition_for_coarsening,
                       p6est_weight_fn);
  p6est->user_pointer = orig_user_pointer;

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
      P4EST_ASSERT (count > 0);
      first = offset;
      last = offset + count;
      offset += count;
      P6EST_COLUMN_SET_RANGE (col, first, last);
    }
  }
  my_count = (p4est_gloidx_t) offset;

  /* calculate the new global_first_layer */
  new_gfl = P4EST_ALLOC (p4est_gloidx_t, p6est->mpisize + 1);
  mpiret = sc_MPI_Allgather (&my_count, 1, P4EST_MPI_GLOIDX, new_gfl, 1,
                             P4EST_MPI_GLOIDX, p6est->mpicomm);
  SC_CHECK_MPI (mpiret);

  psum = 0;
  for (p = 0; p < p6est->mpisize; p++) {
    thiscount = new_gfl[p];
    new_gfl[p] = psum;
    psum += thiscount;
  }
  new_gfl[p6est->mpisize] = psum;
  P4EST_ASSERT (new_gfl[p6est->mpisize] == old_gfl[p6est->mpisize]);

  /* calculate the global number of shipped (received) layers */
  shipped = 0;
  for (i = 0; i < mpisize; i++) {
    p4est_locidx_t      new_total, same, diff;

    new_total = new_gfl[i + 1] - new_gfl[i];
    same =
      SC_MIN (new_gfl[i + 1], old_gfl[i + 1]) - SC_MAX (new_gfl[i],
                                                        old_gfl[i]);
    diff = new_total - SC_MAX (0, same);
    shipped += diff;
  }

  if (old_gfl[rank] == new_gfl[rank] &&
      old_gfl[rank + 1] == new_gfl[rank + 1]) {
    /* if my range is unchanged, I neither send nor receive to anyone */
    p6est->global_first_layer = new_gfl;
    P4EST_FREE (old_gfl);
    p4est_log_indent_pop ();
    P4EST_GLOBAL_PRODUCTIONF
      ("Done p6est_partition shipped %lld layers %.3g%%\n",
       (long long) shipped,
       shipped * 100. / p6est->global_first_layer[p6est->mpisize]);

    return shipped;
  }

  new_layers = sc_array_new_size (sizeof (p2est_quadrant_t), my_count);

  /* initialize views that we can use to bsearch */
  sc_array_init_data (&new_gfl_bsearch, new_gfl, sizeof (p4est_gloidx_t),
                      (size_t) mpisize);
  sc_array_init_data (&old_gfl_bsearch, old_gfl, sizeof (p4est_gloidx_t),
                      (size_t) mpisize);

  /* create the mpi recv requests */
  recv_requests = sc_array_new (sizeof (sc_MPI_Request));
  /* create the list of procs from which to receive */
  recv_procs = sc_array_new (sizeof (int));

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
                             gloidx_compare_overlap);
  P4EST_ASSERT (search >= 0 && search < mpisize);
  overlap = search;
  P4EST_ASSERT (old_gfl[overlap] <= new_gfl[rank] &&
                new_gfl[rank] < old_gfl[overlap + 1]);

  offset = new_gfl[rank];
  self_offset = my_count;
  self_count = 0;
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
        req = (sc_MPI_Request *) sc_array_push (recv_requests);
        ip = (int *) sc_array_push (recv_procs);
        *ip = overlap;
        mpiret = sc_MPI_Irecv (sc_array_index (recv, local_offset),
                               (int) (recv_count * recv->elem_size),
                               sc_MPI_BYTE, overlap, P6EST_COMM_PARTITION,
                               p6est->mpicomm, req);
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
  P4EST_ASSERT (offset == (size_t) new_gfl[rank + 1]);
  nrecv = (int) recv_requests->elem_count;

  /* create the mpi send requests */
  send_requests = sc_array_new (sizeof (sc_MPI_Request));

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
                             gloidx_compare_overlap);
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
    if (send_count) {
      if (overlap != rank) {
        if (data_size) {
          /* pack data for shipping */
          for (zz = 0; zz < (size_t) send_count; zz++) {
            layer =
              p2est_quadrant_array_index (old_layers, zz + local_offset);
            packedlayer =
              (p2est_quadrant_t *) sc_array_index (send, zz + local_offset);
            *packedlayer = *layer;
            packedlayer->p.user_data = NULL;
            if (data_size) {
              layer_data = (void *) (packedlayer + 1);
              memcpy (layer_data, layer->p.user_data, data_size);
              p6est_layer_free_data (p6est, layer);
            }
          }
        }
        /* post the send */
        req = (sc_MPI_Request *) sc_array_push (send_requests);
        mpiret = sc_MPI_Isend (sc_array_index (send, local_offset),
                               (int) (send_count * send->elem_size),
                               sc_MPI_BYTE, overlap, P6EST_COMM_PARTITION,
                               p6est->mpicomm, req);
        SC_CHECK_MPI (mpiret);
      }
      else {
        P4EST_ASSERT ((size_t) send_count == self_count);
        /* copy locally */
        memcpy (sc_array_index (new_layers, self_offset),
                sc_array_index (old_layers, offset - old_gfl[rank]),
                self_count * new_layers->elem_size);
      }
    }
    offset += send_count;
    overlap++;
  }
  P4EST_ASSERT (offset == (size_t) old_gfl[rank + 1]);
  nsend = (int) send_requests->elem_count;

  nleft = nrecv;
  array_of_indices = P4EST_ALLOC (int, nrecv);
  while (nleft > 0) {
    /* finalize some receives */
    mpiret =
      sc_MPI_Waitsome (nrecv, (sc_MPI_Request *) (recv_requests->array),
                       &outcount, array_of_indices, sc_MPI_STATUSES_IGNORE);
    SC_CHECK_MPI (mpiret);
    P4EST_ASSERT (outcount != sc_MPI_UNDEFINED);
    P4EST_ASSERT (outcount >= 0);
    if (!outcount) {
      continue;
    }
    for (i = 0; i < outcount; i++) {
      int                 j;
      j = array_of_indices[i];
      p = *((int *) sc_array_index_int (recv_procs, j));
      P4EST_ASSERT (p >= 0 && p < mpisize);
      P4EST_ASSERT (p != rank);
      if (data_size) {
        /* get the range in the array of this section */
        local_offset = SC_MAX (0, old_gfl[p] - new_gfl[rank]);
        local_last =
          SC_MIN (new_gfl[rank + 1], old_gfl[p + 1]) - new_gfl[rank];
        P4EST_ASSERT (local_last > local_offset);
        for (zz = local_offset; zz < (size_t) local_last; zz++) {
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
  sc_array_destroy (recv_procs);

  /* switch to the new layers */
  p6est->layers = new_layers;
  p6est->global_first_layer = new_gfl;
  P4EST_FREE (old_gfl);

  /* wait for sends to complete */
  mpiret = sc_MPI_Waitall (nsend, (sc_MPI_Request *) (send_requests->array),
                           sc_MPI_STATUSES_IGNORE);
  SC_CHECK_MPI (mpiret);

  /* clean up send and old_layers */
  sc_array_destroy (send);
  sc_array_destroy (send_requests);
  sc_array_destroy (old_layers);

  p4est_log_indent_pop ();
  P4EST_GLOBAL_PRODUCTIONF
    ("Done p6est_partition shipped %lld layers %.3g%%\n",
     (long long) shipped,
     shipped * 100. / p6est->global_first_layer[p6est->mpisize]);

  return shipped;
}

p4est_gloidx_t
p6est_partition (p6est_t * p6est, p6est_weight_t weight_fn)
{
  return p6est_partition_ext (p6est, 0, weight_fn);
}

static int
p6est_column_refine_thin_layer (p6est_t * p6est,
                                p4est_topidx_t which_tree,
                                p4est_quadrant_t * column)
{
  int                 max_diff = *((int *) p6est->user_pointer);
  int8_t              horz_level = column->level;
  size_t              first, last, zz;
  sc_array_t         *layers = p6est->layers;

  P6EST_COLUMN_GET_RANGE (column, &first, &last);

  for (zz = first; zz < last; zz++) {
    p2est_quadrant_t   *layer = p2est_quadrant_array_index (layers, zz);
    int8_t              vert_level = layer->level;

    if (vert_level - horz_level > max_diff) {
      return 1;
    }
  }

  return 0;
}

static int
p6est_layer_refine_thick_layer (p6est_t * p6est,
                                p4est_topidx_t which_tree,
                                p4est_quadrant_t * column,
                                p2est_quadrant_t * layer)
{
  int                 min_diff = *((int *) p6est->user_pointer);
  int8_t              horz_level = column->level;
  int8_t              vert_level = layer->level;

  if (vert_level - horz_level < min_diff) {
    return 1;
  }

  return 0;
}

void
p6est_balance_ext (p6est_t * p6est, p8est_connect_type_t btype,
                   int max_diff, int min_diff,
                   p6est_init_t init_fn, p6est_replace_t replace_fn)
{
  p4est_connect_type_t hbtype;
  p6est_refine_col_data_t refine_col;
  void               *orig_user_pointer = p6est->user_pointer;
  p6est_profile_t    *profile;
  int                 any_change;
  int                 niter;

  P4EST_GLOBAL_PRODUCTIONF ("Into p6est_balance with %lld total layers"
                            " in %lld total columns\n",
                            (long long) p6est->global_first_layer[p6est->
                                                                  mpisize],
                            (long long) p6est->columns->global_num_quadrants);
  p4est_log_indent_push ();

  /* first refine columns whose layers are too thin */
  if (max_diff >= min_diff) {
    int                 print_max_diff = SC_MIN (max_diff, P4EST_QMAXLEVEL);

    print_max_diff = SC_MAX (print_max_diff, -P4EST_QMAXLEVEL);

    P4EST_GLOBAL_PRODUCTIONF
      ("Enforcing maximum layer width:height ratio 2^%d:1\n", print_max_diff);
    p6est->user_pointer = (void *) &max_diff;
    p6est_refine_columns_ext (p6est, 1, -1, p6est_column_refine_thin_layer,
                              init_fn, replace_fn);
    p6est->user_pointer = orig_user_pointer;
  }

  /* next perform the horizontal balancing */
  if (btype == P8EST_CONNECT_FACE) {
    hbtype = P4EST_CONNECT_FACE;
  }
  else {
    hbtype = P4EST_CONNECT_FULL;
  }
  refine_col.refine_col_fn = NULL;
  refine_col.init_fn = init_fn;
  refine_col.replace_fn = replace_fn;
  refine_col.user_pointer = orig_user_pointer;
  p6est->user_pointer = (void *) &refine_col;
  P4EST_GLOBAL_VERBOSE ("Balancing p4est for columns\n");
  p4est_balance_ext (p6est->columns, hbtype, NULL,
                     p6est_replace_column_split);
  p6est->user_pointer = orig_user_pointer;
  p6est_compress_columns (p6est);
  p6est_update_offsets (p6est);

  /* next refine layers that are too thick */
  if (max_diff >= min_diff) {
    int                 print_min_diff = SC_MIN (min_diff, P4EST_QMAXLEVEL);

    print_min_diff = SC_MAX (print_min_diff, -P4EST_QMAXLEVEL);

    P4EST_GLOBAL_PRODUCTIONF
      ("Enforcing minimum layer width:height ratio 2^%d:1\n", print_min_diff);
    p6est->user_pointer = (void *) &min_diff;
    p6est_refine_layers_ext (p6est, 1, -1, p6est_layer_refine_thick_layer,
                             init_fn, replace_fn);
    p6est->user_pointer = orig_user_pointer;
  }

  /* finally, the real work: balance neighboring layers */

  /* initialize the lnodes profile and balance locally */
  profile =
    p6est_profile_new_local (p6est, NULL, P6EST_PROFILE_UNION, btype, 2);
  niter = 0;
  do {
    P4EST_GLOBAL_VERBOSE ("p6est_balance layer balancing iteration begin\n");
    any_change = 0;
    p6est_profile_balance_local (profile);
    any_change = p6est_profile_sync (profile);
    P4EST_GLOBAL_VERBOSE ("p6est_balance layer balancing iteration end\n");
    niter++;
  } while (any_change);
  P4EST_GLOBAL_STATISTICSF ("p6est layers balanced in %d iterations\n",
                            niter);

  p6est_refine_to_profile (p6est, profile, init_fn, replace_fn);

  p6est_profile_destroy (profile);

  p4est_log_indent_pop ();
  P4EST_GLOBAL_PRODUCTIONF ("Done p6est_balance with %lld total layers "
                            "in %lld total columns\n",
                            (long long) p6est->global_first_layer[p6est->
                                                                  mpisize],
                            (long long) p6est->columns->global_num_quadrants);
}

void
p6est_balance (p6est_t * p6est, p8est_connect_type_t btype,
               p6est_init_t init_fn)
{
  p6est_balance_ext (p6est, btype, 0, 1, init_fn, NULL);
}

unsigned
p2est_quadrant_checksum (sc_array_t * quadrants,
                         sc_array_t * checkarray, size_t first_quadrant)
{
#ifdef P4EST_HAVE_ZLIB
  int                 own_check;
  size_t              kz, qcount;
  unsigned            crc;
  uint32_t           *check;
  p2est_quadrant_t   *q;

  qcount = quadrants->elem_count;

  P4EST_ASSERT (quadrants->elem_size == sizeof (p2est_quadrant_t));
  P4EST_ASSERT (first_quadrant <= qcount);

  if (checkarray == NULL) {
    checkarray = sc_array_new (4);
    own_check = 1;
  }
  else {
    P4EST_ASSERT (checkarray->elem_size == 4);
    own_check = 0;
  }

  sc_array_resize (checkarray, (qcount - first_quadrant) * 2);
  for (kz = first_quadrant; kz < qcount; ++kz) {
    q = p2est_quadrant_array_index (quadrants, kz);
    check =
      (uint32_t *) sc_array_index (checkarray, (kz - first_quadrant) * 2);
    check[0] = htonl ((uint32_t) q->z);
    check[1] = htonl ((uint32_t) q->level);
  }
  crc = sc_array_checksum (checkarray);

  if (own_check) {
    sc_array_destroy (checkarray);
  }

  return crc;
#else
  SC_ABORT ("Configure did not find a recent enough zlib.  Abort.\n");
  return 0;
#endif
}

unsigned
p6est_checksum (p6est_t * p6est)
{
#ifdef P4EST_HAVE_ZLIB
  uLong               columncrc, locallayercrc, layercrc;
  sc_array_t          checkarray;
  size_t              scount, globalscount;

  columncrc = p4est_checksum (p6est->columns);

  sc_array_init (&checkarray, 4);
  locallayercrc =
    (uLong) p2est_quadrant_checksum (p6est->layers, &checkarray, 0);
  scount = 4 * checkarray.elem_count;
  P4EST_ASSERT (scount == p6est->layers->elem_count * 8);
  sc_array_reset (&checkarray);

  layercrc = p4est_comm_checksum (p6est->columns, locallayercrc, scount);
  /* 4 bytes per field, two fields */
  globalscount = p6est->global_first_layer[p6est->mpisize] * 8;

  return adler32_combine (columncrc, layercrc, globalscount);
#else
  sc_abort_collective
    ("Configure did not find a recent enough zlib.  Abort.\n");
  return 0;
#endif
}
