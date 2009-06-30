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

#ifndef P8EST_ITERATOR_H
#define P8EST_ITERATOR_H

#include <p8est.h>

SC_EXTERN_C_BEGIN;

typedef struct p8est_qcb_info
{
  p8est_t            *p4est;
  sc_array_t         *ghost_layer;
  p8est_quadrant_t   *quad;
  p4est_locidx_t      tree_local_num;
  p4est_topidx_t      tree;
}
p8est_qcb_info_t;

typedef void        (*p8est_qcb_func_t) (p8est_qcb_info_t * info,
                                         void *user_data);

typedef struct p8est_fcb_info
{
  p8est_t            *p4est;
  sc_array_t         *ghost_layer;
  p8est_quadrant_t   *left_quad;
  p4est_locidx_t      left_tree_local_num;
  p4est_topidx_t      left_tree;
  int                 left_outgoing_face;
  int                 left_corner;
  p8est_quadrant_t   *right_quad;
  p4est_locidx_t      right_tree_local_num;
  p4est_topidx_t      right_tree;
  int                 right_outgoing_face;
  int                 right_corner;
  int                 orientation;
  bool                hanging_flag;
}
p8est_fcb_info_t;

typedef void        (*p8est_fcb_func_t) (p8est_fcb_info_t * info,
                                         void *user_data);

typedef struct p8est_ecb_info
{
  p8est_t            *p4est;
  sc_array_t         *ghost_layer;
  sc_array_t         *quads;
  sc_array_t         *tree_local_num;
  sc_array_t         *tree;
  sc_array_t         *edge_in_zorder;
  sc_array_t         *common_corner;
  bool                hanging_flag;
}
p8est_ecb_info_t;

typedef void        (*p8est_ecb_func_t) (p8est_ecb_info_t * info,
                                         void *user_data);

typedef struct p8est_vcb_info
{
  p8est_t            *p4est;
  sc_array_t         *ghost_layer;
  sc_array_t         *quads;
  sc_array_t         *tree_local_num;
  sc_array_t         *tree;
  sc_array_t         *corner_in_zorder;
}
p8est_vcb_info_t;

typedef void        (*p8est_vcb_func_t) (p8est_vcb_info_t * info,
                                         void *user_data);

void                p8est_iterator (p8est_t * p4est, sc_array_t * ghost_layer,
                                    void *user_data,
                                    p8est_qcb_func_t qcb_func,
                                    p8est_fcb_func_t fcb_func,
                                    p8est_ecb_func_t ecb_func,
                                    p8est_vcb_func_t vcb_func);

SC_EXTERN_C_END;

#endif /* !P8EST_ITERATOR_H */
