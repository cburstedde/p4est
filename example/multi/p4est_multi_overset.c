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
  int                 glosize, glorank, headsize, headrank;
  int                 bgsize;
  p4est_gloidx_t     *gfq = NULL;
  p4est_quadrant_t   *gfp = NULL;

  P4EST_ASSERT (0 <= myrole);
  P4EST_ASSERT (myrole < num_meshes);
  P4EST_ASSERT (mesh_offsets != NULL);
  P4EST_ASSERT ((bgp4est != NULL) == (myrole == 0));
  P4EST_ASSERT ((qpoints != NULL) == (myrole > 0));
  P4EST_ASSERT (intsc_fn != NULL);
  P4EST_ASSERT ((intpl_data != NULL) == (myrole > 0));
  P4EST_ASSERT ((intpl_indices != NULL) == (myrole > 0));
  P4EST_ASSERT ((intpl_fn != NULL) == (myrole == 0));
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

  /* run a barrier to double-check collective calling */
  mpiret = sc_MPI_Barrier (glocomm);
  SC_CHECK_MPI (mpiret);

  /* distribute partition information of background to overset meshes */
  gfq = NULL;
  gfp = NULL;
  headsize = 0;
  headrank = -1;
  if (myrole > 0) {
    gfq = P4EST_ALLOC (p4est_gloidx_t, bgsize + 1);
    gfp = P4EST_ALLOC (p4est_quadrant_t, bgsize + 1);
  }
  if (glorank == mesh_offsets[myrole]) {
    /* this is the head process of a mesh */
    P4EST_LDEBUGF ("Head process of mesh %d is global %d\n", myrole, glorank);

    /* query head communicator */
    mpiret = sc_MPI_Comm_size (headcomm, &headsize);
    SC_CHECK_MPI (mpiret);
    P4EST_ASSERT (headsize == num_meshes);
    mpiret = sc_MPI_Comm_rank (headcomm, &headrank);
    SC_CHECK_MPI (mpiret);
    P4EST_ASSERT (headrank == myrole);

    /* communicate partition from background to each mesh head */
    if (myrole == 0) {
      P4EST_ASSERT (glorank == 0);
      P4EST_ASSERT (gfq == NULL && gfp == NULL);
      gfq = bgp4est->global_first_quadrant;
      gfp = bgp4est->global_first_position;
    }
    mpiret = sc_MPI_Bcast (gfq, bgsize + 1, P4EST_MPI_GLOIDX, 0, headcomm);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Bcast (gfp, (bgsize + 1) * sizeof (p4est_quadrant_t),
                           sc_MPI_BYTE, 0, headcomm);
    SC_CHECK_MPI (mpiret);
    if (myrole == 0) {
      gfq = NULL;
      gfp = NULL;
    }
  }
  if (myrole > 0) {
    /* broadcast partition encoding within each overset mesh */
    mpiret = sc_MPI_Bcast (gfq, bgsize + 1, P4EST_MPI_GLOIDX, 0, rolecomm);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Bcast (gfp, (bgsize + 1) * sizeof (p4est_quadrant_t),
                           sc_MPI_BYTE, 0, rolecomm);
    SC_CHECK_MPI (mpiret);
  }

  /* p4est_search_partition_gfx */

  P4EST_FREE (gfq);
  P4EST_FREE (gfp);
}
