/*
  This file is part of p4est.
  p4est is a C library to manage a collection (a forest) of multiple
  connected adaptive quadtrees or octrees in parallel.

  Copyright (C) 2010 The University of Texas System
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

/** \file p8est_step3.c
 *
 * This 3D example program uses p4est to solve a simple advection problem.  It
 * is numerically very simple, and intended to demonstrate several methods of
 * interacting with the p4est data after it has been refined and partitioned.
 * It demonstrates the construction of ghost layers (see p8est_ghost.h) and
 * communication of ghost-layer data, and it demonstrates iteracting with the
 * quadrants and quadrant boundaries through the p8est_iterate routine (see
 * p8est_iterate.h).
 *
 * the header file p4est_to_p8est.h defines preprocessor macros that map
 * 2D p4est routines and objects to their 3D p8est counterparts.  By including
 * this file and then including the source for the 2D example p4est_step3, we
 * convert the 2D example to a 3D example.
 *
 * It is entirely possible to write a 3D-only program without relying on this
 * mechanism.  In this case use the p8est* header files, functions, and data
 * structures.
 * */
#include <p4est_to_p8est.h>
#include "p4est_step3.c"
