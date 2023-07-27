/*
  This file is part of p4est, version 3.
  p4est is a C library to manage a collection (a forest) of multiple
  connected adaptive quadtrees or octrees in parallel.

  Copyright (C) 2019 individual authors
  Originally written by Carsten Burstedde, Lucas C. Wilcox, and Tobin Isaac

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

#ifndef P4EST3_TO_P8EST3_H
#define P4EST3_TO_P8EST3_H

#ifdef P4EST3_H
#error "The include files p4est3.h and p4est3_to_p8est3.h cannot be combined"
#endif
#define P4_TO_P8

/* functions in p4est3_p4est */
#define p4est3_new_p4est                    p4est3_new_p8est
#define p4est3_connectivity_new_p4est       p4est3_connectivity_new_p8est
#define p4est3_connectivity_new_p4est_brick     p4est3_connectivity_new_p8est_brick
#define p4est3_quadrant_vtable_p4est        p4est3_quadrant_vtable_p8est

/* functions in p4est3_quadrant_yx */
#define p4est3_quadrant_yx_vtable           p4est3_quadrant_zyx_vtable

/* functions in p4est3_quadrant_mort */
#define p4est3_quadrant_mort2d_vtable       p4est3_quadrant_mort3d_vtable

#endif /* !P4EST3_TO_P8EST3_H */
