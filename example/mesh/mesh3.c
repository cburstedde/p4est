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

/*
 * Usage: p8est_mesh <configuration> <level>
 *        possible configurations:
 *        o unit      Refinement on the unit cube.
 *        o periodic  Refinement on the all-periodic unit cube.
 *        o rotwrap   Refinement on a funny some-periodic unit cube.
 *        o twocubes  Refinement on two connected cubes.
 *        o twowrap   Refinement on two connected and wrapped cubes.
 *        o rotcubes  Refinement on a 6-tree configuration.
 *        o shell     Refinement on a 24-tree spherical shell.
 *        o sphere    Refinement on a 13-tree solid sphere.
 */

#include <p4est_to_p8est.h>
#include "mesh2.c"
