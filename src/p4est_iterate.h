/*
  This file is part of p4est.
  p4est is a C library to manage a parallel collection of quadtrees and/or
  octrees.

  Copyright (C) 2009 Carsten Burstedde, Lucas Wilcox,
                     Toby Isaac.

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

#ifndef P4EST_ITERATE_H
#define P4EST_ITERATE_H

#include <p4est.h>

SC_EXTERN_C_BEGIN;

typedef struct p4est_iter_volume_info
{
  p4est_t            *p4est;
  sc_array_t         *ghost_layer;
  p4est_quadrant_t   *quad;
  size_t              quadid;
  p4est_topidx_t      treeid;
}
p4est_iter_volume_info_t;

typedef void        (*p4est_iter_volume_t) (p4est_iter_volume_info_t * info,
                                            void *user_data);

typedef struct p4est_iter_face_info
{
  p4est_t            *p4est;
  sc_array_t         *ghost_layer;
  p4est_quadrant_t   *left_quad;
  ssize_t             left_quadid;
  p4est_topidx_t      left_treeid;
  int                 left_outgoing_face;
  int                 left_corner;
  p4est_quadrant_t   *right_quad;
  ssize_t             right_quadid;
  p4est_topidx_t      right_treeid;
  int                 right_outgoing_face;
  int                 right_corner;
  int                 orientation;
  bool                is_hanging;
}
p4est_iter_face_info_t;

typedef void        (*p4est_iter_face_t) (p4est_iter_face_info_t * info,
                                          void *user_data);

typedef struct p4est_iter_corner_info
{
  p4est_t            *p4est;
  sc_array_t         *ghost_layer;
  sc_array_t         *quads;
  sc_array_t         *quadids;
  sc_array_t         *treeids;
  sc_array_t         *corners;
}
p4est_iter_corner_info_t;

typedef void        (*p4est_iter_corner_t) (p4est_iter_corner_info_t * info,
                                            void *user_data);

void                p4est_iterate (p4est_t * p4est, sc_array_t * ghost_layer,
                                   void *user_data,
                                   p4est_iter_volume_t iter_volume,
                                   p4est_iter_face_t iter_face,
                                   p4est_iter_corner_t iter_corner);

SC_EXTERN_C_END;

#endif /* !P4EST_ITERATE_H */
