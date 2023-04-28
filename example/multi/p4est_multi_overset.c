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
    o->r.nb.qpoints = qpoints;
    o->r.nb.intpl_data = intpl_data;
    o->r.nb.intpl_indices = intpl_indices;
  }
  o->user = user;

  P4EST_LDEBUGF ("Multi overset global rank %d/%d role %d/%d points %ld\n",
                 o->glorank, o->glosize, myrole, num_meshes,
                 myrole > 0 ? (long) qpoints->elem_count : 0);

  /* search query points in background partition */
  overset_search_partition (o);
}
