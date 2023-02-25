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

static const int tag_gfq = P4EST_DIM * 373 + 151;
static const int tag_gfp = P4EST_DIM * 373 + 152;

void
p4est_multi_overset (sc_MPI_Comm glocomm, int myrole,
                     int num_meshes, const int *mesh_offsets,
                     p4est_t *bgp4est, sc_array_t *qpoints)
{
  int                 mpiret;
  int                 glosize, glorank, bgsize;
  int                 ishead;
  int                 i;
  p4est_gloidx_t     *gfq = NULL;
  p4est_quadrant_t   *gfp = NULL;
  sc_MPI_Request     *pr = NULL;

  P4EST_ASSERT (0 <= myrole);
  P4EST_ASSERT (myrole < num_meshes);
  P4EST_ASSERT (mesh_offsets != NULL);
  P4EST_ASSERT ((bgp4est != NULL) == (myrole == 0));
  P4EST_ASSERT ((qpoints != NULL) == (myrole > 0));
  P4EST_ASSERT (myrole == 0 || qpoints->elem_size == 4 * sizeof (double));

  mpiret = sc_MPI_Comm_size (glocomm, &glosize);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (glocomm, &glorank);
  SC_CHECK_MPI (mpiret);
  bgsize = mesh_offsets[1] - mesh_offsets[0];
  P4EST_ASSERT (myrole != 0 || bgsize == bgp4est->mpisize);

  P4EST_ASSERT (mesh_offsets[0] == 0);
  P4EST_ASSERT (glosize == mesh_offsets[num_meshes]);
  P4EST_ASSERT (0 <= glorank && glorank < glosize);

  P4EST_LDEBUGF ("Multi overset global rank %d/%d role %d/%d points %ld\n",
                 glorank, glosize, myrole, num_meshes,
                 myrole > 0 ? (long) qpoints->elem_count : 0);

  /* running a barrier to double-check collective calling */
  mpiret = sc_MPI_Barrier (glocomm);
  SC_CHECK_MPI (mpiret);

  /* distribute partition information of background to overset meshes */
  ishead = 0;
  if (glorank == mesh_offsets[myrole]) {
    /* this is the head process of a mesh */
    ishead = 1;
    P4EST_LDEBUGF ("Head process of mesh %d is global %d\n",
                   myrole, glorank);

    /* communicate partition encoding to overset heads */
    if (myrole == 0) {
      pr = P4EST_ALLOC (sc_MPI_Request, 2 * num_meshes);
      pr[0] = pr[num_meshes] = sc_MPI_REQUEST_NULL;
      for (i = 1; i < num_meshes; ++i) {
        mpiret = sc_MPI_Isend (bgp4est->global_first_quadrant,
                               bgsize + 1, P4EST_MPI_GLOIDX,
                               mesh_offsets[i], tag_gfq, glocomm, pr + i);
        SC_CHECK_MPI (mpiret);
        mpiret = sc_MPI_Isend (bgp4est->global_first_position,
                               (bgsize + 1) * sizeof (p4est_quadrant_t),
                               sc_MPI_BYTE, mesh_offsets[i],
                               tag_gfp, glocomm, pr + (i + num_meshes));
        SC_CHECK_MPI (mpiret);
      }
      mpiret = sc_MPI_Waitall (2 * num_meshes, pr, sc_MPI_STATUSES_IGNORE);
      SC_CHECK_MPI (mpiret);
    }
    else {
      pr = P4EST_ALLOC (sc_MPI_Request, 4);
      pr[0] = pr[1] = pr[2] = pr[3] = sc_MPI_REQUEST_NULL;
      gfq = P4EST_ALLOC (p4est_gloidx_t, bgsize + 1);
      mpiret = sc_MPI_Irecv (gfq, bgsize + 1, P4EST_MPI_GLOIDX,
                             mesh_offsets[0], tag_gfq, glocomm, pr + 0);
      SC_CHECK_MPI (mpiret);
      gfp = P4EST_ALLOC (p4est_quadrant_t, bgsize + 1);
      mpiret = sc_MPI_Irecv (gfp, (bgsize + 1) * sizeof (p4est_quadrant_t),
                             sc_MPI_BYTE, mesh_offsets[0],
                             tag_gfp, glocomm, pr + 1);
      SC_CHECK_MPI (mpiret);

      mpiret = sc_MPI_Waitall (4, pr, sc_MPI_STATUSES_IGNORE);
      SC_CHECK_MPI (mpiret);
    }
    P4EST_FREE (pr);
  }

  /* p4est_search_partition_gfx */

  P4EST_FREE (gfq);
  P4EST_FREE (gfp);
}
