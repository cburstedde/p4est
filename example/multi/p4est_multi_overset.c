/*
  This file is part of p4est.
  p4est is a C library to manage a collection (a forest) of multiple
  connected adaptive quadtrees or octrees in parallel.

  Copyright (C) 2010 The University of Texas System
  Additional copyright (C) 2011 individual authors
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
#include <p4est_search.h>
#include "p4est_multi_overset.h"
#else
#include <p8est_search.h>
#include "p8est_multi_overset.h"
#endif

typedef struct p4est_overset_background
{
  p4est_t            *bgp4est;
  p4est_interpolate_point_t intpl_fn;
}
p4est_overset_background_t;

typedef struct p4est_overset_nearbody
{
  sc_array_t         *qpoints;
  int                 point_size;
  sc_array_t         *send_buffer;
  sc_array_t         *send_indices;
  sc_array_t         *intpl_data;
  sc_array_t         *intpl_indices;
}
p4est_overset_nearbody_t;

typedef struct p4est_overset
{
  /* global communication */
  sc_MPI_Comm         glocomm;
  int                 glosize;
  int                 glorank;

  /* mesh communication */
  sc_MPI_Comm         rolecomm;
  int                 rolerank;
  int                 rolesize;

  /* inter-mesh communication */
  sc_MPI_Comm         headcomm;
  int                 headsize;
  int                 headrank;
  int                 myrole;   /* zero for background mesh */
  int                 ishead;   /* first rank of mesh */
  int                 num_meshes;       /* background plus overset meshes */
  int                 num_overset;      /* number of overset meshes */
  const int          *mesh_offsets;     /* mesh offsets in glocomm */

  /* role specific overset information */
  p4est_intersect_t   intsc_fn;
  union role
  {
    p4est_overset_background_t bg;
    p4est_overset_nearbody_t nb;
  }
  r;
  void               *user;
}
p4est_overset_t;

typedef struct overset_point
{
  void               *point;
  int                 bgrank;
  int                 nbrank;
}
overset_point_t;

typedef struct overset_send_buf
{
  int                 rank;
  sc_array_t          ops;
}
overset_send_buf_t;

typedef struct overset_send_ind
{
  int                 rank;
  sc_array_t          oqs;
}
overset_send_ind_t;

static void
overset_add_point_to_buffer (p4est_overset_nearbody_t *nb,
                             overset_point_t *op, int rank)
{
  size_t              nbuffer;
  overset_send_buf_t *sb;
  overset_send_ind_t *si;

  P4EST_ASSERT (nb != NULL);
  P4EST_ASSERT (0 <= rank);
  P4EST_ASSERT (nb->send_buffer != NULL && nb->send_indices != NULL);

  op->bgrank = rank;            /* mark, that we added this point to the process buffer */

  /* We store the point in the send_buffer of nb. This buffer contains one array
   * of query points for every background process that it needs to send points
   * to. We append the point to the array of query points that corresponds to
   * the background process with rank op->bgrank. For every entry of the send
   * buffer we store its index in the nb->qpoints array in nb->send_indices */
  nbuffer = nb->send_buffer->elem_count;
  sb = NULL;
  si = NULL;
  /* We search in the background mesh partition one process after another.
   * So, op either belongs into the last array of the current send buffer or we
   * need to initialize a new array for a higher background rank */
  if (nbuffer > 0) {
    sb = (overset_send_buf_t *) sc_array_index (nb->send_buffer, nbuffer - 1);
    si =
      (overset_send_ind_t *) sc_array_index (nb->send_indices, nbuffer - 1);
    P4EST_ASSERT (sb->rank == si->rank);
    P4EST_ASSERT (sb->rank <= rank);
    P4EST_ASSERT (sb->ops.elem_count == si->oqs.elem_count);
    P4EST_ASSERT (sb->ops.elem_count > 0);
  }
  if (nbuffer == 0 || sb->rank < rank) {
    sb = (overset_send_buf_t *) sc_array_push (nb->send_buffer);
    si = (overset_send_ind_t *) sc_array_push (nb->send_indices);
    sb->rank = si->rank = rank;
    sc_array_init (&sb->ops, nb->point_size);
    sc_array_init (&si->oqs, sizeof (size_t));
  }
  memcpy (sc_array_push (&sb->ops), op->point, nb->point_size);
  memcpy (sc_array_push (&si->oqs), &op->nbrank, si->oqs.elem_size);
}

static int
overset_intersect_partition_fn (p4est_t *p4est, p4est_topidx_t which_tree,
                                p4est_quadrant_t *quadrant, int pfirst,
                                int plast, void *point)
{
  int                 intersects;

  P4EST_ASSERT (p4est != NULL);
  P4EST_ASSERT (p4est->user_pointer != NULL);
  P4EST_ASSERT (quadrant != NULL);
  P4EST_ASSERT (point != NULL);

  /* verify that point is a valid overset_point_t */
  overset_point_t    *op = (overset_point_t *) point;
  P4EST_ASSERT (op->point != NULL);

  if (op->bgrank >= 0) {
    return 0;                   /* avoid sending the same point twice */
  }

  /* verify that p4est user pointer contains a valid p4est_overset_t */
  p4est_overset_t    *o = (p4est_overset_t *) p4est->user_pointer;
  P4EST_ASSERT (o->myrole > 0);
  P4EST_ASSERT (o->intsc_fn != NULL);

  /* set rank of p4est, if all quadrant ancestors belong to the same process */
  p4est->mpirank = (pfirst == plast) ? pfirst : -1;

  /* actual intersection test */
  intersects =
    o->intsc_fn (p4est, which_tree, quadrant, -1, op->point, o->user);

  if (!intersects) {
    return 0;
  }

  /* we have located the point in the intersection quadrant */
  if (pfirst == plast) {
    /* the point intersects a leaf quadrant of the partition search tree */
    overset_add_point_to_buffer (&o->r.nb, op, pfirst);
  }

  return 1;
}

static void
overset_search_partition (p4est_overset_t *o)
{
  int                 mpiret;
  int                 bgsize;
  p4est_gloidx_t     *gfq = NULL;
  p4est_quadrant_t   *gfp = NULL;

  /* distribute partition information of background to overset meshes */
  bgsize = o->mesh_offsets[1] - o->mesh_offsets[0];
  P4EST_ASSERT (o->myrole != 0 || bgsize == o->r.bg.bgp4est->mpisize);
  if (o->myrole > 0) {
    gfq = P4EST_ALLOC (p4est_gloidx_t, bgsize + 1);
    gfp = P4EST_ALLOC (p4est_quadrant_t, bgsize + 1);
  }
  if (o->ishead) {
    P4EST_LDEBUGF ("Head process of mesh %d is global %d\n", o->myrole, o->glorank);

    /* communicate partition from background to each mesh head */
    if (o->myrole == 0) {
      P4EST_ASSERT (o->glorank == 0);
      P4EST_ASSERT (gfq == NULL && gfp == NULL);
      gfq = o->r.bg.bgp4est->global_first_quadrant;
      gfp = o->r.bg.bgp4est->global_first_position;
    }
    mpiret = sc_MPI_Bcast (gfq, bgsize + 1, P4EST_MPI_GLOIDX, 0, o->headcomm);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Bcast (gfp, (bgsize + 1) * sizeof (p4est_quadrant_t),
                           sc_MPI_BYTE, 0, o->headcomm);
    SC_CHECK_MPI (mpiret);
    if (o->myrole == 0) {
      gfq = NULL;
      gfp = NULL;
    }
  }
  if (o->myrole > 0) {
    /* broadcast partition encoding within each overset mesh */
    /* Todo: switch to p4est_search partition_gfx? */
    mpiret = sc_MPI_Bcast (gfq, bgsize + 1, P4EST_MPI_GLOIDX, 0, o->rolecomm);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Bcast (gfp, (bgsize + 1) * sizeof (p4est_quadrant_t),
                           sc_MPI_BYTE, 0, o->rolecomm);
    SC_CHECK_MPI (mpiret);

    /* search the partition given by gfp and gfq */
    p4est_search_partition_gfx (gfq, gfp, bgsize, gfp[bgsize].p.which_tree, 0,
                                o, NULL, overset_intersect_partition_fn,
                                o->r.nb.qpoints);
  }

  P4EST_FREE (gfq);
  P4EST_FREE (gfp);
}

void
p4est_multi_overset (sc_MPI_Comm glocomm, sc_MPI_Comm headcomm,
                     sc_MPI_Comm rolecomm, int myrole,
                     int num_meshes, const int *mesh_offsets,
                     p4est_t *bgp4est, sc_array_t *qpoints,
                     p4est_intersect_t intsc_fn, sc_array_t *intpl_data,
                     sc_array_t *intpl_indices,
                     p4est_interpolate_point_t intpl_fn, void *user)
{
  int                 mpiret;
  size_t              iz, nqpz;
  overset_point_t    *op;
  overset_send_buf_t *sb;
  overset_send_ind_t *si;
  p4est_overset_t     overset, *o = &overset;

  P4EST_ASSERT (0 <= myrole);
  P4EST_ASSERT (myrole < num_meshes);
  P4EST_ASSERT (mesh_offsets != NULL);
  P4EST_ASSERT ((bgp4est != NULL) == (myrole == 0));
  P4EST_ASSERT ((qpoints != NULL) == (myrole > 0));
  P4EST_ASSERT (myrole == 0 || qpoints->elem_size == 4 * sizeof (double));
  P4EST_ASSERT (intsc_fn != NULL);
  P4EST_ASSERT ((intpl_data != NULL) == (myrole > 0));
  P4EST_ASSERT ((intpl_indices != NULL) == (myrole > 0));
  P4EST_ASSERT (myrole == 0 || intpl_data->elem_count == 0);
  P4EST_ASSERT (myrole == 0 || intpl_indices->elem_count == 0);
  P4EST_ASSERT ((intpl_fn != NULL) == (myrole == 0));
  P4EST_ASSERT (mesh_offsets[0] == 0);

  /* run a barrier to double-check collective calling */
  mpiret = sc_MPI_Barrier (glocomm);
  SC_CHECK_MPI (mpiret);

  /* initialize overset data struct */
  memset (o, -1, sizeof (p4est_overset_t));
  o->glocomm = glocomm;
  mpiret = sc_MPI_Comm_size (glocomm, &o->glosize);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (glocomm, &o->glorank);
  SC_CHECK_MPI (mpiret);
  P4EST_ASSERT (o->glosize == mesh_offsets[num_meshes]);
  P4EST_ASSERT (0 <= o->glorank && o->glorank < o->glosize);
  o->rolecomm = rolecomm;
  mpiret = sc_MPI_Comm_size (rolecomm, &o->rolesize);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (rolecomm, &o->rolerank);
  SC_CHECK_MPI (mpiret);
  o->headcomm = headcomm;
  o->myrole = myrole;
  o->ishead = (o->glorank == mesh_offsets[myrole]);        /* first rank of mesh */
  if (o->ishead) {
    mpiret = sc_MPI_Comm_size (headcomm, &o->headsize);
    SC_CHECK_MPI (mpiret);
    P4EST_ASSERT (o->headsize == num_meshes);
    mpiret = sc_MPI_Comm_rank (headcomm, &o->headrank);
    SC_CHECK_MPI (mpiret);
    P4EST_ASSERT (o->headrank == myrole);
  } else {
    o->headsize = 0;
    o->headrank = -1;
  }
  o->num_meshes = num_meshes;
  o->num_overset = num_meshes - 1;
  o->mesh_offsets = mesh_offsets;
  o->intsc_fn = intsc_fn;
  if (myrole == 0) {
    o->r.bg.bgp4est = bgp4est;
    o->r.bg.intpl_fn = intpl_fn;
  }
  else {
    /* store the query points as overlap_point_t to be able to mark their last
     * appearance in the search */
    nqpz = qpoints->elem_count;
    o->r.nb.qpoints = sc_array_new_count (sizeof (overset_point_t), nqpz);
    for (iz = 0; iz < nqpz; ++iz) {
      op = (overset_point_t *) sc_array_index (o->r.nb.qpoints, iz);
      memset (op, -1, sizeof (overset_point_t));
      op->point = sc_array_index (qpoints, iz);
      op->bgrank = -1;
      op->nbrank = o->rolerank;
    }
    o->r.nb.point_size = qpoints->elem_size;
    o->r.nb.send_buffer = sc_array_new (sizeof (overset_send_buf_t));
    o->r.nb.send_indices = sc_array_new (sizeof (overset_send_ind_t));
    o->r.nb.intpl_data = intpl_data;
    o->r.nb.intpl_indices = intpl_indices;
  }
  o->user = user;

  P4EST_LDEBUGF ("Multi overset global rank %d/%d role %d/%d points %ld\n",
                 o->glorank, o->glosize, myrole, num_meshes,
                 myrole > 0 ? (long) qpoints->elem_count : 0);

  /* search query points in background partition */
  overset_search_partition (o);

  /* free overset struct */
  if (myrole > 0) {
    sc_array_destroy (o->r.nb.qpoints);
    for (iz = 0; iz < o->r.nb.send_buffer->elem_count; iz++) {
      sb = (overset_send_buf_t *) sc_array_index (o->r.nb.send_buffer, iz);
      P4EST_ASSERT (sb->ops.elem_count > 0);
      sc_array_reset (&sb->ops);
    }
    sc_array_destroy (o->r.nb.send_buffer);
    for (iz = 0; iz < o->r.nb.send_indices->elem_count; iz++) {
      si = (overset_send_ind_t *) sc_array_index (o->r.nb.send_indices, iz);
      P4EST_ASSERT (si->oqs.elem_count > 0);
      sc_array_reset (&si->oqs);
    }
    sc_array_destroy (o->r.nb.send_indices);
  }
}
