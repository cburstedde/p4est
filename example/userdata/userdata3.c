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

/*
 * This example program demonstrates how to manage application data.
 *
 *   p8est_userdata <OPTIONS> [<configuration> [<level>]]
 *
 * with possible configurations (default is "unit"):
 *   o unit          Refinement on the unit cube.
 *   o periodic      Unit cube with all-periodic boundary conditions.
 *   o brick         Refinement on a 2x3x5 brick of octrees.
 *   o rotcubes      A collection of six connected rotated cubes.
 *   o sphere        A 13-tree geometry of a solid sphere.
 *   o shell         A 24-tree geometry of a hollow sphere.
 *   o torus         A multi-tree discretization of a torus.
 *
 * The maximum refinement level may be appended (default is 4).
 *
 * The following options are recognized:
 *   --help          Display a usage and help message and exit successfully.
 *   --level         The level may alternatively be specified as an option.
 *
 * Invalid options or arguments result in an error message and exit status.
 */

#include <p4est_to_p8est.h>
#include "userdata2.c"
