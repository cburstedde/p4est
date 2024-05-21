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

#include <p4est_to_p8est.h>

#define P4EST_TNODES_MAXNE 125  /* yet excluding face nodes */

/* *INDENT-OFF* */

/* cube corners */
static const int    n_cornr[ 8] = {  0,  1,  2,  3,  4,  5,  6,  7 };

#if 0 /* not yet implemented */

/* cube center */
static const int    n_center =       8;

/* face midpoints */
static const int    n_mface[ 6] = {  9, 10, 11, 12, 13, 14 };

/* edge midpoints */
static const int    n_medge[12] = { 15, 16, 17, 18, 19, 20,
                                    21, 22, 23, 24, 25, 26 };

/* edges between center and corners */
static const int    n_cedge[ 8] = { 27, 28, 29, 30, 31, 32, 33, 34 };

/* edges between center and face midpoints */
static const int    n_fedge[ 6] = { 35, 36, 37, 38, 39, 40 };

/* quarter cube faces */
static const int    n_hface[6][4] =
  {{ 41, 42, 43, 44 }, { 45, 46, 47, 48 }, { 49, 50, 51, 52 },
   { 53, 54, 55, 56 }, { 57, 58, 59, 60 }, { 61, 62, 63, 63 }};

#endif

/* nodes for all 90 tetrahedra: 4 corners, 6 edges, 4 faces, 1 volume */
const int p8est_tnodes_tet_nodes[90][15] =
  /* six original tetrahedra that split the cube */
  {{  0, 1, 3, 7 },
   {  0, 2, 6, 7 },
   {  0, 3, 2, 7 },
   {  0, 4, 5, 7 },
   {  0, 5, 1, 7 },
   {  0, 6, 4, 7 },
   { -1 }};

/* *INDENT-OF* */

#include "p4est_tnodes.c"
