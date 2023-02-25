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
p4est_multi_overset (sc_MPI_Comm glocomm, int myrole,
                     int num_meshes, const int *mesh_offsets,
                     sc_array_t *qpoints)
{
  int                 mpiret;
  int                 glosize, glorank;

  P4EST_ASSERT (0 <= myrole);
  P4EST_ASSERT (myrole < num_meshes);
  P4EST_ASSERT (mesh_offsets != NULL);
  P4EST_ASSERT (qpoints != NULL);
  P4EST_ASSERT (qpoints->elem_size == 4 * sizeof (double));

  mpiret = sc_MPI_Comm_size (glocomm, &glosize);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (glocomm, &glorank);
  SC_CHECK_MPI (mpiret);

  P4EST_ASSERT (mesh_offsets[0] == 0);
  P4EST_ASSERT (glosize == mesh_offsets[num_meshes]);
  P4EST_ASSERT (0 <= glorank && glorank < glosize);

  P4EST_LDEBUGF ("Hello multi overset global rank %d/%d role %d/%d\n",
                 glorank, glosize, myrole, num_meshes);
}
