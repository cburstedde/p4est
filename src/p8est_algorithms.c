/*
  This file is part of p4est.
  p4est is a C library to manage a parallel collection of quadtrees and/or
  octrees.

  Copyright (C) 2008,2009 Carsten Burstedde, Lucas Wilcox.

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <p4est_to_p8est.h>
#include <p8est_connectivity.h>

/* *INDENT-OFF* */

/** Store the number of quadrants to add for complete and balance stages. */
static const int    p4est_balance_count[P4EST_DIM + 1] =
{ 9, 12, 15, 16 };

/** Store coordinates of quadrants to add for balancing. */
static const p4est_qcoord_t p4est_balance_coord[26][P4EST_DIM] =
{ /* faces */
  { -1,  1,  1 },
  {  2,  0,  0 },
  {  1, -1,  1 },
  {  0,  2,  0 },
  {  1,  1, -1 },
  {  0,  0,  2 },
  /* edges */
  {  1, -1, -1 },
  {  1,  2, -1 },
  {  0, -1,  2 },
  {  0,  2,  2 },
  { -1,  1, -1 },
  {  2,  1, -1 },
  { -1,  0,  2 },
  {  2,  0,  2 },
  { -1, -1,  1 },
  {  2, -1,  1 },
  { -1,  2,  0 },
  {  2,  2,  0 },
  /* corners */
  { -1, -1, -1 },
  {  2, -1, -1 },
  { -1,  2, -1 },
  {  2,  2, -1 },
  { -1, -1,  2 },
  {  2, -1,  2 },
  { -1,  2,  2 },
  {  2,  2,  2 }};

/** Offset for edges into p4est_balance_coord */
static const int    pbeo = P4EST_FACES;

/** Offset for corners into p4est_balance_coord */
static const int    pbco = P4EST_FACES + P8EST_EDGES;

/* *INDENT-ON* */

#include "p4est_algorithms.c"
