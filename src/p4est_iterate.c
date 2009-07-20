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

#ifdef P4_TO_P8
#include <p8est_algorithms.h>
#include <p8est_bits.h>
#include <p8est_iterate.h>
#else
#include <p4est_algorithms.h>
#include <p4est_bits.h>
#include <p4est_iterate.h>
#endif

typedef struct p4est_iter_corner_args
{
  p4est_t            *p4est;
  int                 level;
  sc_array_t         *corners;
  int                 num_sides;
  int                *start_idx2;
  int                *quad_idx2;
  sc_array_t         *ghost_layer;
  sc_array_t        **quadrants;
  sc_array_t         *quads;
  sc_array_t         *quadids;
  sc_array_t         *treeids;
  size_t            **index;
  size_t             *first_index;
  size_t             *count;
  p4est_quadrant_t  **test;
}
p4est_iter_corner_args_t;

static void
p4est_corner_iterate (p4est_iter_corner_args_t * args, void *user_data,
                      p4est_iter_corner_t iter_corner)
{
  const int           local = 0;
  const int           ghost = 1;
  const int           idx2_stride = P4EST_CHILDREN + 1;

  int                 side, sidetype;
  int                 level = args->level;
  sc_array_t         *corners = args->corners;
  int                 num_sides = args->num_sides;
  size_t              num_ghosts = args->ghost_layer->elem_count;
  const int          *start_idx2 = args->start_idx2;
  int                *quad_idx2 = args->quad_idx2;
  int                 this_corner;
  sc_array_t        **quadrants = args->quadrants;
  size_t            **index = args->index;
  size_t             *first_index = args->first_index;
  size_t             *count = args->count;
  p4est_quadrant_t  **test = args->test;
  p4est_quadrant_t  **ptemp;
  p4est_quadrant_t    temp;
  p4est_qcoord_t      mask =
    ((p4est_qcoord_t) - 1) << (P4EST_MAXLEVEL - level);
  sc_array_t          test_view;
  p4est_iter_corner_info_t info;
  ssize_t             temp_idx;
  ssize_t            *ploc;
  int                 level_idx2;
  int                 type;
  bool                all_empty, has_local;

  info.p4est = args->p4est;
  info.ghost_layer = args->ghost_layer;
  info.quads = args->quads;
  info.quadids = args->quadids;
  info.treeids = args->treeids;
  info.corners = args->corners;

  test_view.elem_size = sizeof (p4est_quadrant_t);
  test_view.byte_alloc = -1;

  level_idx2 = level * idx2_stride;

  for (side = 0; side < num_sides; side++) {
    quad_idx2[side] = level_idx2 + start_idx2[side];
    ptemp = sc_array_index_int (info.quads, side);
    *ptemp = NULL;
    for (type = local; type <= ghost; type++) {
      sidetype = side * 2 + type;
      first_index[sidetype] = index[sidetype][quad_idx2[side]];
      count[sidetype] = (index[sidetype][quad_idx2[side] + 1] -
                         first_index[sidetype]);
    }
  }

  all_empty = true;
  for (side = 0; side < num_sides; side++) {
    if (count[side * 2 + local]) {
      all_empty = false;
    }
  }
  if (all_empty) {
    return;
  }

  has_local = false;
  for (side = 0; side < num_sides; side++) {
    P4EST_ASSERT ((size_t) side < corners->elem_count);
    this_corner = *((int *) sc_array_index_int (corners, side));
    for (type = local; type <= ghost; type++) {
      sidetype = side * 2 + type;
      if (count[sidetype]) {
        if (count[sidetype] == 1) {
          test[sidetype] = sc_array_index (quadrants[sidetype],
                                           first_index[sidetype]);
          temp_idx = 0;
        }
        else {
          switch (this_corner) {
          case (0):
            P4EST_ASSERT (first_index[sidetype] <
                          quadrants[sidetype]->elem_count);
            test[sidetype] = sc_array_index (quadrants[sidetype],
                                             first_index[sidetype]);
            temp_idx = 0;
            break;
          case (P4EST_CHILDREN - 1):
            P4EST_ASSERT (first_index[sidetype] + count[sidetype] - 1 <
                          quadrants[sidetype]->elem_count);
            test[sidetype] = sc_array_index (quadrants[sidetype],
                                             first_index[sidetype] +
                                             count[sidetype] - 1);
            temp_idx = ((ssize_t) count[sidetype]) - 1;
            break;
          default:
            P4EST_ASSERT (first_index[sidetype] <
                          quadrants[sidetype]->elem_count);
            test[sidetype] = sc_array_index (quadrants[sidetype],
                                             first_index[sidetype]);
            temp = *(test[sidetype]);
            temp.x &= mask;
            temp.y &= mask;
#ifdef P4_TO_P8
            temp.z &= mask;
#endif
            temp.level = P4EST_QMAXLEVEL;
            P4EST_ASSERT (p4est_quadrant_is_valid (&temp));
            temp.x += (this_corner % 2) ?
              P4EST_QUADRANT_LEN (level) -
              P4EST_QUADRANT_LEN (P4EST_QMAXLEVEL) : 0;
            temp.y += ((this_corner % 4) / 2) ?
              P4EST_QUADRANT_LEN (level) -
              P4EST_QUADRANT_LEN (P4EST_QMAXLEVEL) : 0;
#ifdef P4_TO_P8
            temp.z += (this_corner / 4) ?
              P4EST_QUADRANT_LEN (level) -
              P4EST_QUADRANT_LEN (P4EST_QMAXLEVEL) : 0;
#endif
            P4EST_ASSERT (p4est_quadrant_is_valid (&temp));
            test_view.elem_count = count[sidetype];
            test_view.array = quadrants[sidetype]->array +
              sizeof (p4est_quadrant_t) * first_index[sidetype];
            temp_idx = p4est_find_higher_bound (&test_view, &temp, 0);
            if (temp_idx == -1) {
              test[sidetype] = NULL;
            }
            else {
              P4EST_ASSERT ((size_t) temp_idx < test_view.elem_count);
              test[sidetype] = sc_array_index_ssize_t (&test_view, temp_idx);
            }
            break;
          }
        }
        if (test[sidetype] != NULL) {
          temp_idx += first_index[sidetype];
          temp = *(test[sidetype]);
          temp.x += (this_corner % 2) ? P4EST_QUADRANT_LEN (temp.level) : 0;
          temp.y += ((this_corner % 4) / 2) ?
            P4EST_QUADRANT_LEN (temp.level) : 0;
#ifdef P4_TO_P8
          temp.z += (this_corner / 4) ? P4EST_QUADRANT_LEN (temp.level) : 0;
#endif
          P4EST_ASSERT (p4est_quadrant_is_extended (&temp));
          if (temp.x & ~mask || temp.y & ~mask
#ifdef P4_TO_P8
              || temp.z & ~mask
#endif
            ) {
            test[sidetype] = NULL;
          }
          else {
            P4EST_ASSERT ((size_t) side < info.quads->elem_count);
            ptemp = sc_array_index_int (info.quads, side);
            *ptemp = test[sidetype];
            P4EST_ASSERT ((size_t) side < info.quadids->elem_count);
            ploc = sc_array_index_int (info.quadids, side);
            *ploc = temp_idx - (ssize_t) ((type == local) ? 0 : num_ghosts);
            if (type == local) {
              has_local = true;
            }
          }
        }
      }
    }
  }
  if (!has_local) {
    return;
  }
#ifdef P4EST_DEBUG
  for (side = 0; side < num_sides; side++) {
    ptemp = sc_array_index_int (info.quads, side);
    if (*ptemp == NULL) {
      P4EST_PRODUCTIONF ("side %d not filled\n", side);
    }
  }
#endif
  iter_corner (&info, user_data);
}

#ifdef P4_TO_P8
typedef struct p8est_iter_edge_args
{
  int                 num_sides;
  p8est_t            *p8est;
  sc_array_t         *ghost_layer;
  int                 level;
  int                *start_idx2;
  int                *level_num;
  sc_array_t         *common_corners[2];
  sc_array_t        **quadrants;
  sc_array_t         *edges;
  sc_array_t         *treeids;
  sc_array_t         *quadids;
  sc_array_t         *quads;
  size_t            **index;
  size_t             *first_index;
  p8est_quadrant_t  **test;
  size_t             *count;
  int                *test_level;
  int                *quad_idx2;
  bool               *refine;
  p4est_iter_corner_args_t *corner_args;
}
p8est_iter_edge_args_t;

static void
p8est_edge_iterate (p8est_iter_edge_args_t * args, void *user_data,
                    p8est_iter_edge_t iter_edge,
                    p8est_iter_corner_t iter_corner)
{
  const int           local = 0;
  const int           ghost = 1;
  const int           idx2_stride = P4EST_CHILDREN + 1;

  int                 num_sides = args->num_sides;
  int                 start_level = args->level;
  size_t              num_ghosts = args->ghost_layer->elem_count;
  int                *start_idx2 = args->start_idx2;
  int                *level_num = args->level_num;
  sc_array_t        **quadrants = args->quadrants;
  size_t            **index = args->index;
  size_t             *first_index = args->first_index;
  sc_array_t        **common_corners = args->common_corners;
  p8est_quadrant_t  **test = args->test;
  size_t             *count = args->count;
  int                *test_level = args->test_level;
  int                *quad_idx2 = args->quad_idx2;
  bool               *refine = args->refine;
  p8est_quadrant_t   *test2;
  p8est_quadrant_t  **temp;
  ssize_t            *temp_loc;
  int                *temp_int;
  int                 i;
  int                 level;
  int                 side;
  int                 type;
  int                 sidetype;
  int                 level_idx2;
  p8est_iter_edge_info_t info;
  sc_array_t          test_view;
  bool                all_empty, stop_refine;
  int                 temp_idx2;
  p4est_iter_corner_args_t *corner_args = args->corner_args;
  sc_array_t         *corners = corner_args->corners;
  int                *c;

  level = start_level;
  level_idx2 = level * idx2_stride;
  for (side = 0; side < num_sides; side++) {
    quad_idx2[side] = level_idx2 + start_idx2[side];
    for (type = local; type <= ghost; type++) {
      sidetype = side * 2 + type;
      first_index[sidetype] = index[sidetype][quad_idx2[side]];
      count[sidetype] = (index[sidetype][quad_idx2[side] + 1] -
                         first_index[sidetype]);
    }
  }
  all_empty = true;
  for (side = 0; side < num_sides; side++) {
    if (count[side * 2 + local]) {
      all_empty = false;
      break;
    }
  }
  if (all_empty) {
    return;
  }

  info.p4est = args->p8est;
  info.ghost_layer = args->ghost_layer;
  info.quads = args->quads;
  info.quadids = args->quadids;
  info.treeids = args->treeids;
  info.edges = args->edges;

  test_view.elem_size = sizeof (p4est_quadrant_t);
  test_view.byte_alloc = -1;

  level_num[start_level] = 0;
  for (;;) {
    for (side = 0; side < num_sides; side++) {
      for (type = local; type <= ghost; type++) {
        sidetype = side * 2 + type;
        if (count[sidetype]) {
          P4EST_ASSERT (first_index[sidetype] <
                        quadrants[sidetype]->elem_count);
          test[sidetype] = sc_array_index (quadrants[sidetype],
                                           first_index[sidetype]);
          test_level[sidetype] = (int) test[sidetype]->level;
        }
        else {
          test_level[sidetype] = -1;
        }
      }
      refine[side] = true;
    }
    stop_refine = false;
    for (side = 0; side < num_sides; side++) {
      for (type = local; type <= ghost; type++) {
        sidetype = side * 2 + type;
        if (test_level[sidetype] == level) {
          P4EST_ASSERT ((size_t) side < info.quads->elem_count);
          temp = sc_array_index_int (info.quads, side);
          *temp = test[sidetype];
          P4EST_ASSERT ((size_t) side < info.quadids->elem_count);
          temp_loc = sc_array_index_int (info.quadids, side);
          *temp_loc = (ssize_t) first_index[sidetype] -
            (ssize_t) ((type == local) ? 0 : num_ghosts);
          refine[side] = false;
          stop_refine = true;
        }
      }
    }
    all_empty = true;
    for (side = 0; side < num_sides; side++) {
      if (refine[side]) {
        all_empty = false;
        break;
      }
    }
    if (all_empty) {
      info.is_hanging = false;
      info.common_corners = common_corners[0];
      iter_edge (&info, user_data);
      level_num[level]++;
      goto change_search_area;
    }
    for (side = 0; side < num_sides; side++) {
      if (refine[side]) {
        quad_idx2[side] = level_idx2 + idx2_stride;
        for (type = local; type <= ghost; type++) {
          sidetype = side * 2 + type;
          test_view.elem_count = count[sidetype];
          test_view.array = quadrants[sidetype]->array +
            sizeof (p4est_quadrant_t) * first_index[sidetype];
          p4est_split_array (&test_view, level,
                             index[sidetype] + quad_idx2[side]);
          for (i = 0; i < idx2_stride; i++) {
            index[sidetype][quad_idx2[side] + i] += first_index[sidetype];
          }
        }
      }
    }
    if (stop_refine) {
      info.is_hanging = true;
      for (i = 0; i < 2; i++) {
        all_empty = true;
        info.common_corners = common_corners[i];
        for (side = 0; side < num_sides; side++) {
          if (refine[side]) {
            P4EST_ASSERT ((size_t) side < common_corners[i]->elem_count);
            temp_int = sc_array_index_int (common_corners[i], side);
            quad_idx2[side] = level_idx2 + idx2_stride + *temp_int;
            for (type = local; type <= ghost; type++) {
              sidetype = side * 2 + type;
              first_index[sidetype] = index[sidetype][quad_idx2[side]];
              count[sidetype] =
                (size_t) index[sidetype][quad_idx2[side] + 1] -
                first_index[sidetype];
              if (count[sidetype]) {
                P4EST_ASSERT (count[sidetype] == 1);
                P4EST_ASSERT (first_index[sidetype] <
                              quadrants[sidetype]->elem_count);
                test2 =
                  sc_array_index (quadrants[sidetype], first_index[sidetype]);
                P4EST_ASSERT ((int) test2->level == level + 1);
                P4EST_ASSERT ((size_t) side < info.quads->elem_count);
                temp = sc_array_index_int (info.quads, side);
                *temp = test2;
                P4EST_ASSERT ((size_t) side < info.quadids->elem_count);
                temp_loc = sc_array_index_int (info.quadids, side);
                *temp_loc = (ssize_t) first_index[sidetype] -
                  (ssize_t) ((type == local) ? 0 : num_ghosts);
                if (type == local) {
                  all_empty = false;
                }
              }
            }
          }
          else {
            temp_loc = sc_array_index_int (info.quadids, side);
            if (*temp_loc >= 0) {
              all_empty = false;
            }
          }
        }
        if (!all_empty) {
          iter_edge (&info, user_data);
        }
      }
      level_num[level]++;
      goto change_search_area;
    }
    level_num[++level] = 0;
    level_idx2 += idx2_stride;
  change_search_area:
    if (level_num[start_level] > 0) {
      break;
    }
    if (level_num[level] == 2) {
      // corner iterate goes here
      if (iter_corner != NULL) {
        corner_args->level = level;
        for (i = 0; i < 2 * num_sides; i++) {
          c = sc_array_index_int (corners, i);
          temp_int = sc_array_index_int (common_corners[1 - i / num_sides],
                                         i % num_sides);
          *c = *temp_int;
          temp_int = sc_array_index_int (common_corners[i / num_sides],
                                         i % num_sides);
          start_idx2[i] = *temp_int;
          temp_idx2 = level_idx2 + start_idx2[i];
          index[i * 2 + local][temp_idx2] =
            index[(i % num_sides) * 2 + local][temp_idx2];
          index[i * 2 + local][temp_idx2 + 1] =
            index[(i % num_sides) * 2 + local][temp_idx2 + 1];
          index[i * 2 + ghost][temp_idx2] =
            index[(i % num_sides) * 2 + ghost][temp_idx2];
          index[i * 2 + ghost][temp_idx2 + 1] =
            index[(i % num_sides) * 2 + ghost][temp_idx2 + 1];
        }
        p4est_corner_iterate (corner_args, user_data, iter_corner);
      }
      level_num[--level]++;
      level_idx2 -= idx2_stride;
      goto change_search_area;
    }
    all_empty = true;
    for (side = 0; side < num_sides; side++) {
      P4EST_ASSERT ((size_t) side <
                    common_corners[level_num[level]]->elem_count);
      temp_int = sc_array_index_int (common_corners[level_num[level]], side);
      quad_idx2[side] = level_idx2 + *temp_int;
      for (type = local; type <= ghost; type++) {
        sidetype = side * 2 + type;
        first_index[sidetype] = index[sidetype][quad_idx2[side]];
        count[sidetype] = (index[sidetype][quad_idx2[side] + 1] -
                           first_index[sidetype]);
        if (type == local && count[sidetype]) {
          all_empty = false;
        }
      }
    }
    if (all_empty) {
      level_num[level]++;
      goto change_search_area;
    }
  }
}
#endif

typedef struct p4est_iter_face_args
{
  p4est_t            *p4est;
  sc_array_t         *ghost_layer;
  int                 level;
  int                *start_idx2;
  int                 face[2];                /** z-order */
  int                *level_num;
  int                 orientation;
  int                *num_to_child;
  sc_array_t        **quadrants;
  size_t            **index;
  size_t             *first_index;
  p4est_topidx_t      treeids[2];
  p4est_quadrant_t  **test;
  size_t             *count;
  int                *test_level;
  int                *quad_idx2;
  bool               *refine;
  bool                outside_face;
#ifdef P4_TO_P8
  p8est_iter_edge_args_t *edge_args;
#endif
  p4est_iter_corner_args_t *corner_args;
}
p4est_iter_face_args_t;

static void
p4est_face_iterate (p4est_iter_face_args_t * args, void *user_data,
                    p4est_iter_face_t iter_face,
#ifdef P4_TO_P8
                    p8est_iter_edge_t iter_edge,
#endif
                    p4est_iter_corner_t iter_corner)
{

  const int           left = 0;
  const int           right = 1;
  const int           local = 0;
  const int           ghost = 1;
  const int           idx2_stride = P4EST_CHILDREN + 1;
  const int           ntc_str = P4EST_CHILDREN / 2;

  int                 start_level = args->level;
  int                *start_idx2 = args->start_idx2;
  const int          *face = args->face;
  int                *level_num = args->level_num;
  sc_array_t        **quadrants = args->quadrants;
  size_t            **index = args->index;
  size_t             *first_index = args->first_index;
  int                *num_to_child = args->num_to_child;
  p4est_quadrant_t  **test = args->test;
  size_t             *count = args->count;
  int                *test_level = args->test_level;
  int                *quad_idx2 = args->quad_idx2;
  bool               *refine = args->refine;
  int                 limit;

  p4est_quadrant_t   *test2;
  int                 i, j;
  int                 level;
  int                 side, n_side;
  int                 type, n_type;
  int                 sidetype, nsidentype;
  int                 level_idx2;
  p4est_iter_face_info_t info;
  sc_array_t          test_view;
  size_t              num_ghosts = args->ghost_layer->elem_count;
  int                 temp_idx2;
  p4est_iter_corner_args_t *corner_args = args->corner_args;
  sc_array_t         *corners = corner_args->corners;
  int                *c;
#ifdef P4_TO_P8
  int                 v, true_dir, dir, k;
  p8est_iter_edge_args_t *edge_args = args->edge_args;
  sc_array_t         *common_corners[2];

  sc_array_t         *edges = edge_args->edges;
  int                *e_common_corners[2];
  int                *e_edges;

  common_corners[0] = edge_args->common_corners[0];
  common_corners[1] = edge_args->common_corners[1];
  e_common_corners[0] = (int *) common_corners[0]->array;
  e_common_corners[1] = (int *) common_corners[1]->array;
  e_edges = (int *) edges->array;
#endif

  limit = args->outside_face ? left : right;
  level = start_level;
  level_idx2 = level * idx2_stride;
  for (side = left; side <= limit; side++) {
    quad_idx2[side] = level_idx2 + start_idx2[side];
    for (type = local; type <= ghost; type++) {
      sidetype = side * 2 + type;
      first_index[sidetype] = index[sidetype][quad_idx2[side]];
      count[sidetype] = (index[sidetype][quad_idx2[side] + 1] -
                         first_index[sidetype]);
    }
  }
  if (!args->outside_face) {
    if (!count[left * 2 + local] && !count[right * 2 + local]) {
      return;
    }
  }
  else {
    if (!count[left * 2 + local]) {
      return;
    }
  }

  info.p4est = args->p4est;
  info.ghost_layer = args->ghost_layer;
  info.orientation = args->orientation;

  test_view.elem_size = sizeof (p4est_quadrant_t);
  test_view.byte_alloc = -1;

  level_num[start_level] = 0;
  for (;;) {
    for (side = left; side <= limit; side++) {
      for (type = local; type <= ghost; type++) {
        sidetype = side * 2 + type;
        if (count[sidetype]) {
          P4EST_ASSERT (first_index[sidetype] <
                        quadrants[sidetype]->elem_count);
          test[sidetype] = sc_array_index (quadrants[sidetype],
                                           first_index[sidetype]);
          test_level[sidetype] = (int) test[sidetype]->level;
        }
        else {
          test_level[sidetype] = -1;
        }
      }
    }
    refine[left] = refine[right] = true;
    for (side = left; side <= limit; side++) {
      for (type = local; type <= ghost; type++) {
        sidetype = side * 2 + type;
        if (test_level[sidetype] == level) {
          if (iter_face != NULL) {
            P4EST_ASSERT (count[sidetype] == 1);
            P4EST_ASSERT (count[side * 2 + (type ^ 1)] == 0);
            info.left_quad = test[sidetype];
            info.left_quadid = (ssize_t) first_index[sidetype] -
              (ssize_t) ((type == local) ? 0 : num_ghosts);
            info.left_outgoing_face = face[side];
            info.left_treeid = args->treeids[side];
            refine[side] = false;
            if (!args->outside_face) {
              n_side = side ^ 1;
              info.right_outgoing_face = face[n_side];
              info.right_treeid = args->treeids[n_side];
              for (n_type = type; n_type <= ghost; n_type++) {
                nsidentype = n_side * 2 + n_type;
                if ((n_type > type || n_side > side) &&
                    test_level[nsidentype] == level) {
                  P4EST_ASSERT (count[nsidentype] == 1);
                  P4EST_ASSERT (count[n_side * 2 + (n_type ^ 1)] == 0);
                  info.right_quad = test[nsidentype];
                  info.right_quadid =
                    (ssize_t) first_index[nsidentype] -
                    (ssize_t) ((n_type == local) ? 0 : num_ghosts);
                  info.left_corner = num_to_child[side * ntc_str];
                  info.right_corner = num_to_child[n_side * ntc_str];
                  info.is_hanging = false;
                  P4EST_ASSERT (!(type == ghost && n_type == ghost));
                  iter_face (&info, user_data);
                  level_num[level]++;
                  goto change_search_area;
                }
              }
#ifdef P4EST_DEBUG
              if (type == local) {
                P4EST_ASSERT (count[n_side * 2 + type] > 0 ||
                              count[n_side * 2 + (type ^ 1)] > 0);
              }
#endif
            }
            else {
              info.right_outgoing_face = info.left_outgoing_face;
              info.right_treeid = info.left_treeid;
              info.right_quad = info.left_quad;
              info.right_quadid = info.left_quadid;
              info.left_corner = num_to_child[side * ntc_str];
              info.right_corner = info.left_corner;
              info.is_hanging = false;
              iter_face (&info, user_data);
              level_num[level]++;
              goto change_search_area;
            }
          }
          else {
            level_num[level]++;
            goto change_search_area;
          }
        }
      }
    }
    for (side = left; side <= limit; side++) {
      if (refine[side]) {
        quad_idx2[side] = level_idx2 + idx2_stride;
        for (type = local; type <= ghost; type++) {
          sidetype = side * 2 + type;
          test_view.elem_count = count[sidetype];
          test_view.array = quadrants[sidetype]->array +
            sizeof (p4est_quadrant_t) * first_index[sidetype];
          p4est_split_array (&test_view, level,
                             index[sidetype] + quad_idx2[side]);
          for (i = 0; i < idx2_stride; i++) {
            index[sidetype][quad_idx2[side] + i] += first_index[sidetype];
          }
        }
      }
    }
    for (side = left; side <= limit; side++) {
      if (!refine[side]) {
        n_side = side ^ 1;
        info.is_hanging = true;
        for (type = local; type <= ghost; type++) {
          if (test_level[side * 2 + type] == level) {
            for (i = 0; i < P4EST_CHILDREN / 2; i++) {
              info.left_corner = num_to_child[side * ntc_str + i];
              info.right_corner = num_to_child[n_side * ntc_str + i];
              quad_idx2[n_side] = level_idx2 + idx2_stride +
                num_to_child[n_side * ntc_str + i];
              for (n_type = local; n_type <= ghost; n_type++) {
                nsidentype = n_side * 2 + n_type;
                first_index[nsidentype] =
                  index[nsidentype][quad_idx2[n_side]];
                count[nsidentype] =
                  (index[nsidentype][quad_idx2[n_side] + 1] -
                   first_index[nsidentype]);
                if (count[nsidentype] && (type + n_type < 2)) {
                  P4EST_ASSERT (count[nsidentype] == 1);
                  P4EST_ASSERT (index[n_side * 2 + (n_type ^ 1)]
                                [quad_idx2[n_side]]
                                == index[n_side * 2 + (n_type ^ 1)]
                                [quad_idx2[n_side] + 1]);
                  test2 = sc_array_index (quadrants[nsidentype],
                                          first_index[nsidentype]);
                  info.right_quad = test2;
                  info.right_quadid =
                    (ssize_t) first_index[nsidentype] -
                    (ssize_t) ((n_type == local) ? 0 : num_ghosts);
                  iter_face (&info, user_data);
                  P4EST_ASSERT ((int) test2->level == level + 1);
                }
              }
#ifdef P4EST_DEBUG
              if (type == local) {
                P4EST_ASSERT (count[n_side * 2 + local] == 1 ||
                              count[n_side * 2 + ghost] == 1);
              }
#endif
            }
            level_num[level]++;
            goto change_search_area;
          }
        }
      }
    }
    level_num[++level] = 0;
    level_idx2 += idx2_stride;
  change_search_area:
    if (level_num[start_level] > 0) {
      break;
    }
    if (level_num[level] == P4EST_CHILDREN / 2) {
#ifdef P4_TO_P8
      // edge iterator goes here
      if (iter_edge != NULL || iter_corner != NULL) {
        edge_args->level = level;
        if (!args->outside_face) {
          for (dir = 0; dir < 2; dir++) {
            for (side = 0; side < 2; side++) {
              for (j = 0; j < 4; j++) {
                k = j >> 1;
                if (dir == 0) {
                  e_common_corners[side][j] =
                    num_to_child[(j % 2) * ntc_str + (2 * (1 - k) + side)];
                }
                else {
                  e_common_corners[side][j] =
                    num_to_child[(j % 2) * ntc_str + ((1 - k) + 2 * side)];
                }
              }
            }
            for (j = 0; j < 4; j++) {
              v = e_common_corners[0][j];
              true_dir = e_common_corners[1][j] - v;
              true_dir = (true_dir > 0) ? true_dir : -true_dir;
              e_edges[j] = p8est_corner_edges[v][true_dir >> 1];
            }
            for (side = 0; side < 2; side++) {
              for (j = 0; j < 4; j++) {
                k = j >> 1;
                if (dir == 0) {
                  start_idx2[j] =
                    num_to_child[(j % 2) * ntc_str + side + 2 * k];
                }
                else {
                  start_idx2[j] =
                    num_to_child[(j % 2) * ntc_str + 2 * side + k];
                }
                temp_idx2 = level_idx2 + start_idx2[j];
                index[j * 2 + local][temp_idx2] = index[(j % 2) * 2 + local]
                  [temp_idx2];
                index[j * 2 + local][temp_idx2 + 1] =
                  index[(j % 2) * 2 + local]
                  [temp_idx2 + 1];
                index[j * 2 + ghost][temp_idx2] = index[(j % 2) * 2 + ghost]
                  [temp_idx2];
                index[j * 2 + ghost][temp_idx2 + 1] =
                  index[(j % 2) * 2 + ghost]
                  [temp_idx2 + 1];
              }
              p8est_edge_iterate (edge_args, user_data, iter_edge,
                                  iter_corner);
            }
          }
        }
        else {
          for (dir = 0; dir < 2; dir++) {
            for (side = 0; side < 2; side++) {
              for (j = 0; j < 2; j++) {
                if (dir == 0) {
                  e_common_corners[side][j] =
                    num_to_child[(2 * (1 - j) + side)];
                }
                else {
                  e_common_corners[side][j] =
                    num_to_child[((1 - j) + 2 * side)];
                }
              }
            }
            for (j = 0; j < 2; j++) {
              v = e_common_corners[0][j];
              true_dir = e_common_corners[1][j] - v;
              true_dir = (true_dir > 0) ? true_dir : -true_dir;
              e_edges[j] = p8est_corner_edges[v][true_dir >> 1];
            }
            for (side = 0; side < 2; side++) {
              for (j = 0; j < 2; j++) {
                if (dir == 0) {
                  start_idx2[j] = num_to_child[side + 2 * j];
                }
                else {
                  start_idx2[j] = num_to_child[2 * side + j];
                }
                temp_idx2 = level_idx2 + start_idx2[j];
                index[j * 2 + local][temp_idx2] = index[local][temp_idx2];
                index[j * 2 + local][temp_idx2 + 1] = index[local]
                  [temp_idx2 + 1];
                index[j * 2 + ghost][temp_idx2] = index[ghost][temp_idx2];
                index[j * 2 + ghost][temp_idx2 + 1] = index[ghost]
                  [temp_idx2 + 1];
              }
              p8est_edge_iterate (edge_args, user_data, iter_edge,
                                  iter_corner);
            }
          }
        }
      }
#endif
      // corner iterate goes here
      if (iter_corner != NULL) {
        corner_args->level = level;
        if (!args->outside_face) {
          for (j = 0; j < P4EST_CHILDREN; j++) {
            c = sc_array_index_int (corners, j);
            *c = num_to_child
              [(j % 2) * ntc_str + P4EST_CHILDREN / 2 - 1 - (j / 2)];
            start_idx2[j] = num_to_child[(j % 2) * ntc_str + (j / 2)];
            temp_idx2 = level_idx2 + start_idx2[j];
            index[j * 2 + local][temp_idx2] = index[(j % 2) * 2 + local]
              [temp_idx2];
            index[j * 2 + local][temp_idx2 + 1] = index[(j % 2) * 2 + local]
              [temp_idx2 + 1];
            index[j * 2 + ghost][temp_idx2] = index[(j % 2) * 2 + ghost]
              [temp_idx2];
            index[j * 2 + ghost][temp_idx2 + 1] = index[(j % 2) * 2 + ghost]
              [temp_idx2 + 1];
          }
        }
        else {
          for (j = 0; j < P4EST_CHILDREN / 2; j++) {
            c = sc_array_index_int (corners, j);
            *c = num_to_child[P4EST_CHILDREN / 2 - 1 - j];
            start_idx2[j] = num_to_child[j];
            temp_idx2 = level_idx2 + start_idx2[j];
            index[j * 2 + local][temp_idx2] = index[local][temp_idx2];
            index[j * 2 + local][temp_idx2 + 1] = index[local]
              [temp_idx2 + 1];
            index[j * 2 + ghost][temp_idx2] = index[ghost][temp_idx2];
            index[j * 2 + ghost][temp_idx2 + 1] = index[ghost]
              [temp_idx2 + 1];
          }
        }
        p4est_corner_iterate (corner_args, user_data, iter_corner);
      }
      level_num[--level]++;
      level_idx2 -= idx2_stride;
      goto change_search_area;
    }
    for (side = left; side <= limit; side++) {
      quad_idx2[side] =
        level_idx2 + num_to_child[side * ntc_str + level_num[level]];
    }
    for (side = left; side <= limit; side++) {
      for (type = local; type <= ghost; type++) {
        sidetype = side * 2 + type;
        first_index[sidetype] = index[sidetype][quad_idx2[side]];
        count[sidetype] = (index[sidetype][quad_idx2[side] + 1] -
                           first_index[sidetype]);
      }
    }
    if (!args->outside_face) {
      if (!count[left * 2 + local] && !count[right * 2 + local]) {
        level_num[level]++;
        goto change_search_area;
      }
    }
    else {
      if (!count[left * 2 + local]) {
        level_num[level]++;
        goto change_search_area;
      }
    }
  }
}

void static
p4est_volume_iterate ()
{
};

void
p4est_iterate (p4est_t * p4est, sc_array_t * ghost_layer, void *user_data,
               p4est_iter_volume_t iter_volume, p4est_iter_face_t iter_face,
#ifdef P4_TO_P8
               p8est_iter_edge_t iter_edge,
#endif
               p4est_iter_corner_t iter_corner)
{
  /* constants */
  const int           left = 0;
  const int           right = 1;
  const int           local = 0;
  const int           ghost = 1;
  const int           idx2_stride = P4EST_CHILDREN + 1;
  const int           ntc_str = P4EST_CHILDREN / 2;

  sc_array_t         *trees = p4est->trees, **quadrants, test_view;
  p4est_tree_t       *tree;
  p4est_connectivity_t *conn = p4est->connectivity;
  size_t              num_ghosts, num_local;
  p4est_quadrant_t  **test;
#ifndef P4_TO_P8
  int                 rface[2];
#endif
  int                 level, type, left_side, right_side,
    *test_level, *level_num, modulus, i, j, k, *face,
    num_to_child[P4EST_CHILDREN], *start_idx2;
  int8_t             *ttf = conn->tree_to_face;
  p4est_topidx_t     *ttt = conn->tree_to_tree, nt, t, c;
  size_t              global_num_trees, guess_tree;
  int                 alloc_size = 4;
  size_t            **index, *this_index, *tree_first_ghost, *first_index;
  size_t              guess, guess_low, *count;
#ifdef P4_TO_P8
  p4est_topidx_t     *tte = conn->tree_to_edge;
  p4est_topidx_t      num_edges = conn->num_edges;
  p4est_topidx_t     *ett_offset = conn->ett_offset;
  p4est_topidx_t     *edge_to_tree = conn->edge_to_tree;
  p4est_topidx_t      this_edge, e;
  int8_t             *edge_to_edge = conn->edge_to_edge;
  int                 max_edge_size = 4;
  p8est_iter_edge_args_t edge_args;
  int                 edge_size;
  sc_array_t          common_corners[2];
  sc_array_t          edges;
  sc_array_t          edge_treeids;
  sc_array_t          edge_quadids;
  sc_array_t          edge_quads;
  int                *e_common_corners[2], *e_edges;
  int                 perm, corner0, corner1, ne, m, l, dir, side;
  ssize_t            *e_quadids;
  p4est_quadrant_t  **e_quads;
  p4est_topidx_t     *e_treeids;
#endif
  int                 orientation, nc;
  p4est_topidx_t      this_corner;
#ifndef P4_TO_P8
  p4est_topidx_t     *ttc = conn->tree_to_vertex;
  p4est_topidx_t      num_corners = conn->num_vertices;
  p4est_topidx_t     *ctt_offset = conn->vtt_offset;
  p4est_corner_transform_t *ct;
  sc_array_t          ctransforms, *cta = &ctransforms;
#else
  p4est_topidx_t     *ttc = conn->tree_to_corner;
  p4est_topidx_t      num_corners = conn->num_corners;
  p4est_topidx_t     *ctt_offset = conn->ctt_offset;
  p4est_topidx_t     *corner_to_tree = conn->corner_to_tree;
  int8_t             *corner_to_corner = conn->corner_to_corner;
  p4est_topidx_t      jt, kt;
#endif
  int                 max_corner_size = P4EST_CHILDREN;
  int                 corner_size;
  sc_array_t          corners;
  sc_array_t          corner_treeids;
  sc_array_t          corner_quadids;
  sc_array_t          corner_quads;
  int                *c_corners;
  ssize_t            *c_quadids;
  p4est_quadrant_t  **c_quads;
  p4est_topidx_t     *c_treeids;
  int                *quad_idx2;
  bool               *refine;
  p4est_iter_face_args_t face_args;
  p4est_iter_corner_args_t corner_args;
  int                 level_idx2;
  p4est_iter_volume_info_t info;

  P4EST_ASSERT (p4est_is_valid (p4est));
  P4EST_ASSERT (ghost_layer != NULL);
  P4EST_ASSERT (sc_array_is_sorted
                (ghost_layer, p4est_quadrant_compare_piggy));

  if (p4est->first_local_tree < 0) {
    return;
  }

  if (iter_face == NULL && iter_corner == NULL
#ifdef P4_TO_P8
      && iter_edge == NULL
#endif
    ) {
    p4est_volume_iterate ();
    return;
  }

#ifdef P4_TO_P8
  if (iter_edge != NULL || iter_corner != NULL) {
    for (e = 0; e < num_edges; e++) {
      edge_size = (int) (ett_offset[e + 1] - ett_offset[e]);
      max_edge_size = (edge_size > max_edge_size) ? edge_size : max_edge_size;
    }
  }
  alloc_size = (2 * max_edge_size > alloc_size) ?
    2 * max_edge_size : alloc_size;
  max_corner_size = (max_edge_size * 2 > max_corner_size) ?
    max_edge_size * 2 : max_corner_size;
#endif

  if (iter_corner != NULL) {
    for (c = 0; c < num_corners; c++) {
      corner_size = (int) (ctt_offset[c + 1] - ctt_offset[c]);
      max_corner_size = (corner_size > max_corner_size) ? corner_size :
        max_corner_size;
    }
  }
  alloc_size = (2 * max_corner_size > alloc_size) ?
    2 * max_corner_size : alloc_size;

  /** initialize arrays that keep track of where we are in the search */
  index = P4EST_ALLOC (size_t *, alloc_size);
  quadrants = P4EST_ALLOC (sc_array_t *, alloc_size);
  start_idx2 = P4EST_ALLOC (int, alloc_size / 2);
  first_index = P4EST_ALLOC (size_t, alloc_size);
  test = P4EST_ALLOC (p4est_quadrant_t *, alloc_size);
  count = P4EST_ALLOC (size_t, alloc_size);
  test_level = P4EST_ALLOC (int, alloc_size);
  quad_idx2 = P4EST_ALLOC (int, alloc_size / 2);
  refine = P4EST_ALLOC (bool, alloc_size / 2);
  for (i = 0; i < alloc_size; i++) {
    index[i] = P4EST_ALLOC (size_t, P4EST_MAXLEVEL * idx2_stride);
    if (i % 2)
      quadrants[i] = ghost_layer;
  }
  level_num = P4EST_ALLOC (int, P4EST_MAXLEVEL);

  test_view.elem_size = sizeof (p4est_quadrant_t);
  test_view.byte_alloc = -1;

  /** Divide the ghost_layer by p.which_tree */
  global_num_trees = trees->elem_count;
  tree_first_ghost = P4EST_ALLOC (size_t, global_num_trees + 1);
  guess_low = 0;
  tree_first_ghost[0] = 0;
  num_ghosts = ghost_layer->elem_count;
  for (t = 1; t <= global_num_trees; t++) {
    tree_first_ghost[t] = num_ghosts;
  }
  if (num_ghosts) {
    guess = num_ghosts / 2;
    t = 1;
    for (;;) {
      P4EST_ASSERT (guess < ghost_layer->elem_count);
      test[0] = sc_array_index (ghost_layer, guess);
      guess_tree = test[0]->p.which_tree;
      if (guess_tree < t) {
        guess_low = guess + 1;
      }
      else {
        for (nt = t; nt <= guess_tree; nt++) {
          tree_first_ghost[nt] = guess;
        }
      }
      while (tree_first_ghost[t] == guess_low) {
        t++;
      }
      if (t == global_num_trees || guess_low == num_ghosts)
        break;
      guess = guess_low + (tree_first_ghost[t] - guess_low) / 2;
    }
  }

  /** set up the arguments passed to face_iterate that are invariant */
  face_args.p4est = p4est;
  face_args.ghost_layer = ghost_layer;
  face_args.level_num = level_num;
  face_args.index = index;
  face_args.quadrants = quadrants;
  face_args.num_to_child = num_to_child;
  face_args.start_idx2 = start_idx2;
  face_args.first_index = first_index;
  face_args.test = test;
  face_args.count = count;
  face_args.test_level = test_level;
  face_args.quad_idx2 = quad_idx2;
  face_args.refine = refine;
#ifdef P4_TO_P8
  face_args.edge_args = &edge_args;
#endif
  face_args.corner_args = &corner_args;
  face = face_args.face;

  /** set up the arguments passed to corner_iterate that are invariant */
  if (iter_corner != NULL) {
    corner_args.p4est = p4est;
    corner_args.ghost_layer = ghost_layer;
    corner_args.quadrants = quadrants;
    corner_args.index = index;
    corner_args.start_idx2 = start_idx2;
    corner_args.first_index = first_index;
    corner_args.test = test;
    corner_args.count = count;
    corner_args.quad_idx2 = quad_idx2;
    c_corners = P4EST_ALLOC (int, max_corner_size);
    corners.elem_size = sizeof (int);
    corners.byte_alloc = -1;
    corners.array = (char *) (c_corners);
    corner_args.corners = &corners;
    c_treeids = P4EST_ALLOC (p4est_topidx_t, max_corner_size);
    c_quadids = P4EST_ALLOC (ssize_t, max_corner_size);
    c_quads = P4EST_ALLOC (p4est_quadrant_t *, max_corner_size);
    corner_treeids.elem_size = sizeof (p4est_topidx_t);
    corner_treeids.byte_alloc = -1;
    corner_treeids.array = (char *) c_treeids;
    corner_quadids.elem_size = sizeof (ssize_t);
    corner_quadids.byte_alloc = -1;
    corner_quadids.array = (char *) c_quadids;
    corner_quads.elem_size = sizeof (p4est_quadrant_t *);
    corner_quads.byte_alloc = -1;
    corner_quads.array = (char *) c_quads;
    corner_args.treeids = &corner_treeids;
    corner_args.quadids = &corner_quadids;
    corner_args.quads = &corner_quads;
  }
  else {
    c_treeids = NULL;
    c_quads = NULL;
    c_quadids = NULL;
    c_corners = NULL;
  }

#ifdef P4_TO_P8
  /** set up the arguments passed to edge_iterate that are invariant */
  if (iter_edge != NULL || iter_corner != NULL) {
    edge_args.p8est = p4est;
    edge_args.ghost_layer = ghost_layer;
    edge_args.level_num = level_num;
    edge_args.quadrants = quadrants;
    edge_args.index = index;
    edge_args.start_idx2 = start_idx2;
    edge_args.first_index = first_index;
    edge_args.test = test;
    edge_args.count = count;
    edge_args.test_level = test_level;
    edge_args.quad_idx2 = quad_idx2;
    edge_args.refine = refine;
    e_common_corners[0] = P4EST_ALLOC (int, max_edge_size);
    e_common_corners[1] = P4EST_ALLOC (int, max_edge_size);
    e_edges = P4EST_ALLOC (int, max_edge_size);
    common_corners[0].elem_size = sizeof (int);
    common_corners[0].byte_alloc = -1;
    common_corners[0].array = (char *) (e_common_corners[0]);
    common_corners[1].elem_size = sizeof (int);
    common_corners[1].byte_alloc = -1;
    common_corners[1].array = (char *) (e_common_corners[1]);
    edges.elem_size = sizeof (int);
    edges.byte_alloc = -1;
    edges.array = (char *) e_edges;
    edge_args.common_corners[0] = &(common_corners[0]);
    edge_args.common_corners[1] = &(common_corners[1]);
    edge_args.edges = &edges;
    edge_args.corner_args = &corner_args;
    e_treeids = P4EST_ALLOC (p4est_topidx_t, max_edge_size);
    e_quadids = P4EST_ALLOC (ssize_t, max_edge_size);
    e_quads = P4EST_ALLOC (p4est_quadrant_t *, max_edge_size);
    edge_treeids.elem_size = sizeof (p4est_topidx_t);
    edge_treeids.byte_alloc = -1;
    edge_treeids.array = (char *) e_treeids;
    edge_quadids.elem_size = sizeof (ssize_t);
    edge_quadids.byte_alloc = -1;
    edge_quadids.array = (char *) e_quadids;
    edge_quads.elem_size = sizeof (p4est_quadrant_t *);
    edge_quads.byte_alloc = -1;
    edge_quads.array = (char *) e_quads;
    edge_args.treeids = &edge_treeids;
    edge_args.quadids = &edge_quadids;
    edge_args.quads = &edge_quads;
  }
  else {
    e_treeids = NULL;
    e_quads = NULL;
    e_quadids = NULL;
    e_edges = NULL;
  }
#endif

  /** set up the arguments passed to iter_volume that are invariant */
  info.p4est = p4est;
  info.ghost_layer = ghost_layer;

  for (t = 0; t < global_num_trees; t++) {

    tree = sc_array_index (trees, t);
    quadrants[local] = &(tree->quadrants);

    P4EST_ASSERT (p4est_tree_is_sorted (tree));
    P4EST_ASSERT (tree->maxlevel <= P4EST_QMAXLEVEL);

    this_index = index[left * 2 + local];
    this_index[0] = 0;
    num_local = quadrants[local]->elem_count;
    this_index[1] = (p4est_locidx_t) num_local;
    this_index = index[left * 2 + ghost];
    this_index[0] = tree_first_ghost[t];
    this_index[1] = tree_first_ghost[t + 1];
    num_ghosts = (this_index[1] - this_index[0]);

#ifdef P4EST_DEBUG
    if (num_local) {
      this_index = index[left * 2 + local];
      P4EST_ASSERT (this_index[0] < quadrants[local]->elem_count);
      test[0] = sc_array_index (quadrants[local], this_index[0]);
      P4EST_ASSERT (p4est_quadrant_is_inside_root (test[0]));
      P4EST_ASSERT ((this_index[1] - 1)
                    < quadrants[local]->elem_count);
      test[0] = sc_array_index (quadrants[local], this_index[1] - 1);
      P4EST_ASSERT (p4est_quadrant_is_inside_root (test[0]));
    }
    if (num_ghosts) {
      this_index = index[left * 2 + ghost];
      P4EST_ASSERT (this_index[0] < quadrants[ghost]->elem_count);
      test[0] = sc_array_index (quadrants[ghost], this_index[0]);
      P4EST_ASSERT (p4est_quadrant_is_inside_root (test[0]));
      P4EST_ASSERT ((this_index[1] - 1) < quadrants[ghost]->elem_count);
      test[0] = sc_array_index (quadrants[ghost], this_index[1] - 1);
      P4EST_ASSERT (p4est_quadrant_is_inside_root (test[0]));
    }
#endif

    face_args.treeids[left] = t;
    face_args.treeids[right] = t;
    face_args.outside_face = false;
    face_args.quadrants[left * 2 + local] = quadrants[local];
    face_args.quadrants[right * 2 + local] = quadrants[local];
    face_args.orientation = 0;

#ifdef P4_TO_P8
    if (iter_edge != NULL || iter_corner != NULL) {
      edge_args.num_sides = 4;
      common_corners[0].elem_count = 4;
      common_corners[1].elem_count = 4;
      edges.elem_count = 4;
      edge_treeids.elem_count = 4;
      e_treeids[0] = e_treeids[1] = e_treeids[2] = e_treeids[3] = t;
      edge_quadids.elem_count = 4;
      edge_quads.elem_count = 4;
      quadrants[2] = quadrants[local];
      quadrants[4] = quadrants[local];
      quadrants[6] = quadrants[local];
    }
#endif

    if (iter_corner != NULL) {
      corner_args.num_sides = P4EST_CHILDREN;
      corners.elem_count = P4EST_CHILDREN;
      for (i = 0; i < P4EST_CHILDREN; i++) {
        c_treeids[i] = t;
        quadrants[2 * i] = quadrants[local];
      }
      corner_treeids.elem_count = P4EST_CHILDREN;
      corner_quadids.elem_count = P4EST_CHILDREN;
      corner_quads.elem_count = P4EST_CHILDREN;
    }

    info.treeid = t;

    level = 0;
    level_idx2 = 0;
    level_num[0] = 0;

    if (num_local) {
      while (level_num[0] == 0) {
        if (level_num[level] == P4EST_CHILDREN) {
          for (left_side = 0, right_side = 1; right_side < 2 * P4EST_DIM;
               left_side += 2, right_side += 2) {
            face[left] = right_side;
            face[right] = left_side;
            for (i = 0; i < P4EST_CHILDREN / 2; i++) {
#ifndef P4_TO_P8
              num_to_child[left * ntc_str + i] = p4est_face_corners
                [p4est_zface_to_rface[right_side]][i];
              num_to_child[right * ntc_str + i] = p4est_face_corners
                [p4est_zface_to_rface[left_side]][i];
#else
              num_to_child[left * ntc_str + i] =
                p8est_face_corners[right_side][i];
              num_to_child[right * ntc_str + i] =
                p8est_face_corners[left_side][i];
#endif
            }
            face_args.level = level;
            for (i = 0; i < P4EST_CHILDREN / 2; i++) {
#ifndef P4_TO_P8
              start_idx2[left] =
                p4est_face_corners[p4est_zface_to_rface[left_side]][i];
              start_idx2[right] =
                p4est_face_corners[p4est_zface_to_rface[right_side]][i];
#else
              start_idx2[left] = p8est_face_corners[left_side][i];
              start_idx2[right] = p8est_face_corners[right_side][i];
#endif
              /** initialize index[right] for use in face_adjacency */
              quad_idx2[0] = level_idx2 + start_idx2[right];
              index[right * 2 + local][quad_idx2[0]] =
                index[left * 2 + local][quad_idx2[0]];
              index[right * 2 + local][quad_idx2[0] + 1] =
                index[left * 2 + local][quad_idx2[0] + 1];
              index[right * 2 + ghost][quad_idx2[0]] =
                index[left * 2 + ghost][quad_idx2[0]];
              index[right * 2 + ghost][quad_idx2[0] + 1] =
                index[left * 2 + ghost][quad_idx2[0] + 1];

              p4est_face_iterate (&face_args, user_data, iter_face,
#ifdef P4_TO_P8
                                  iter_edge,
#endif
                                  iter_corner);
            }
          }
#ifdef P4_TO_P8
          if (iter_edge != NULL || iter_corner != NULL) {
            edge_args.level = level;
            for (dir = 0; dir < P4EST_DIM; dir++) {
              for (side = 0; side < 2; side++) {
                for (j = 0; j < 4; j++) {
                  e_common_corners[side][j] =
                    p8est_face_corners[dir * 2 + side][3 - j];
                }
              }
              for (j = 0; j < 4; j++) {
                e_edges[j] = 4 * dir + (3 - j);
              }
              for (side = 0; side < 2; side++) {
                for (j = 0; j < 4; j++) {
                  start_idx2[j] = p8est_face_corners[dir * 2 + side][j];
                  if (j != 0) {
                    quad_idx2[0] = level_idx2 + start_idx2[j];
                    index[j * 2 + local][quad_idx2[0]] =
                      index[local][quad_idx2[0]];
                    index[j * 2 + local][quad_idx2[0] + 1] =
                      index[local][quad_idx2[0] + 1];
                    index[j * 2 + ghost][quad_idx2[0]] =
                      index[ghost][quad_idx2[0]];
                    index[j * 2 + ghost][quad_idx2[0] + 1] =
                      index[ghost][quad_idx2[0] + 1];
                  }
                }
                p8est_edge_iterate (&edge_args, user_data, iter_edge,
                                    iter_corner);
              }
            }
          }
#endif
          if (iter_corner != NULL) {
            corner_args.level = level;
            for (i = 0; i < P4EST_CHILDREN; i++) {
              c_corners[i] = P4EST_CHILDREN - i - 1;
              start_idx2[i] = i;
              quad_idx2[0] = level_idx2 + start_idx2[i];
              index[i * 2 + local][quad_idx2[0]] = index[local][quad_idx2[0]];
              index[i * 2 + local][quad_idx2[0] + 1] =
                index[local][quad_idx2[0] + 1];
              index[i * 2 + ghost][quad_idx2[0]] = index[ghost][quad_idx2[0]];
              index[i * 2 + ghost][quad_idx2[0] + 1] =
                index[ghost][quad_idx2[0] + 1];
            }
            p4est_corner_iterate (&corner_args, user_data, iter_corner);
          }
          level_num[--level]++;
          level_idx2 -= idx2_stride;
          continue;
        }
        quad_idx2[0] = level_idx2 + level_num[level];
        for (type = local; type <= ghost; type++) {
          this_index = index[left * 2 + type];
          first_index[type] = this_index[quad_idx2[0]];
          count[type] = (this_index[quad_idx2[0] + 1] - first_index[type]);
        }
        if (!count[local]) {
          level_num[level]++;
          continue;
        }
        for (type = local; type <= ghost; type++) {
          if (count[type]) {
            P4EST_ASSERT (first_index[type] < quadrants[type]->elem_count);
            test[type] = sc_array_index (quadrants[type], first_index[type]);
            test_level[type] = (int) test[type]->level;
            if (test_level[type] == level) {
              if (type == local) {
                info.quad = test[type];
                info.quadid = first_index[type];
                if (iter_volume != NULL)
                  iter_volume (&info, user_data);
              }
              break;
            }
          }
          else {
            test_level[type] = -1;
          }
        }
        if (type <= ghost) {
          level_num[level]++;
          continue;
        }
        quad_idx2[0] = level_idx2 + idx2_stride;
        for (type = local; type <= ghost; type++) {
          test_view.elem_count = count[type];
          P4EST_ASSERT (first_index[type] + count[type]
                        <= quadrants[type]->elem_count);
          test_view.array = quadrants[type]->array +
            sizeof (p4est_quadrant_t) * first_index[type];
          this_index = index[left * 2 + type];
          p4est_split_array (&test_view, level, this_index + quad_idx2[0]);

          for (i = 0; i < idx2_stride; i++) {
            this_index[quad_idx2[0] + i] += first_index[type];
          }
        }
        level_num[++level] = 0;
        level_idx2 += idx2_stride;
      }
    }

    for (i = 0; i < 2 * P4EST_DIM; i++) {
      face_args.level = 0;
      face_args.start_idx2[left] = face_args.start_idx2[right] = 0;
#ifndef P4_TO_P8
      rface[left] = p4est_zface_to_rface[i];
#endif
      for (j = 0; j < P4EST_CHILDREN / 2; j++) {
#ifndef P4_TO_P8
        num_to_child[left * ntc_str + j] = p4est_face_corners[rface[left]][j];
#else
        num_to_child[left * ntc_str + j] = p4est_face_corners[i][j];
#endif
      }
      face[left] = i;
#ifndef P4_TO_P8
      nt = ttt[t * 2 * P4EST_DIM + rface[left]];
      rface[right] = (int) ttf[t * 2 * P4EST_DIM + rface[left]];
      face_args.orientation = rface[right] / 4;
      rface[right] %= 4;
      P4EST_ASSERT (ttt[nt * 2 * P4EST_DIM + rface[right]] == t);
      face[right] = p4est_rface_to_zface[rface[right]];
#else
      nt = ttt[t * 2 * P4EST_DIM + face[left]];
      face[right] = (int) ttf[t * 2 * P4EST_DIM + face[left]];
      face_args.orientation = face[right] / 6;
      face[right] %= 6;
#endif
      if ((t > nt) || (t == nt && face[left] > face[right])) {
        face_args.outside_face = false;
#ifndef P4_TO_P8
        modulus = face[left] ^ face[right];
        modulus = ((modulus & 2) >> 1) ^ (modulus & 1) ^ 1;
        if ((modulus ^ face_args.orientation) == 0) {
          num_to_child[right * ntc_str] = p4est_face_corners[rface[right]][0];
          num_to_child[right * ntc_str + 1] =
            p4est_face_corners[rface[right]][1];
        }
        else {
          num_to_child[right * ntc_str] = p4est_face_corners[rface[right]][1];
          num_to_child[right * ntc_str + 1] =
            p4est_face_corners[rface[right]][0];
        }
#else
        modulus = p8est_face_permutation_refs[face[left]][face[right]];
        modulus = p8est_face_permutation_sets[modulus][face_args.orientation];
        for (j = 0; j < P4EST_CHILDREN / 2; j++) {
          k = p8est_face_permutations[modulus][j];
          num_to_child[right * ntc_str + j] =
            p8est_face_corners[face[right]][k];
        }
#endif
        level_num[0] = 0;
        P4EST_ASSERT (nt < trees->elem_count);
        tree = sc_array_index (trees, nt);
        face_args.quadrants[right * 2 + local] = &(tree->quadrants);
        face_args.treeids[right] = nt;
        index[right * 2 + local][0] = 0;
        index[right * 2 + local][1] = num_local = tree->quadrants.elem_count;
        index[right * 2 + ghost][0] = tree_first_ghost[nt];
        index[right * 2 + ghost][1] = tree_first_ghost[nt + 1];

#ifdef P4_TO_P8
        if (iter_edge != NULL || iter_corner != NULL) {
          edge_args.num_sides = 4;
          common_corners[0].elem_count = 4;
          common_corners[1].elem_count = 4;
          edges.elem_count = 4;
          edge_treeids.elem_count = 4;
          e_treeids[0] = e_treeids[2] = t;
          e_treeids[1] = e_treeids[3] = nt;
          edge_quadids.elem_count = 4;
          edge_quads.elem_count = 4;
          quadrants[2] = &(tree->quadrants);
          quadrants[4] = quadrants[local];
          quadrants[6] = &(tree->quadrants);
        }
#endif

        if (iter_corner != NULL) {
          corner_args.num_sides = P4EST_CHILDREN;
          corners.elem_count = P4EST_CHILDREN;
          for (j = 0; j < P4EST_CHILDREN; j++) {
            if (j % 2) {
              c_treeids[j] = nt;
              quadrants[2 * j] = &(tree->quadrants);
            }
            else {
              c_treeids[j] = t;
              quadrants[2 * j] = quadrants[local];
            }
          }
          corner_treeids.elem_count = P4EST_CHILDREN;
          corner_quadids.elem_count = P4EST_CHILDREN;
          corner_quads.elem_count = P4EST_CHILDREN;
        }

        p4est_face_iterate (&face_args, user_data, iter_face,
#ifdef P4_TO_P8
                            iter_edge,
#endif
                            iter_corner);
      }
      else if (t == nt && face[left] == face[right]) {
        face_args.outside_face = true;
        level_num[0] = 0;
#ifdef P4_TO_P8
        if (iter_edge != NULL || iter_corner != NULL) {
          edge_args.num_sides = 2;
          common_corners[0].elem_count = 2;
          common_corners[1].elem_count = 2;
          edges.elem_count = 2;
          edge_treeids.elem_count = 2;
          e_treeids[0] = e_treeids[1] = t;
          edge_quadids.elem_count = 2;
          edge_quads.elem_count = 2;
          quadrants[2] = quadrants[local];
        }
#endif

        if (iter_corner != NULL) {
          corner_args.num_sides = P4EST_CHILDREN / 2;
          corners.elem_count = P4EST_CHILDREN / 2;
          for (j = 0; j < P4EST_CHILDREN / 2; j++) {
            c_treeids[j] = t;
            quadrants[2 * j] = quadrants[local];
          }
          corner_treeids.elem_count = P4EST_CHILDREN / 2;
          corner_quadids.elem_count = P4EST_CHILDREN / 2;
          corner_quads.elem_count = P4EST_CHILDREN / 2;
        }

        p4est_face_iterate (&face_args, user_data, iter_face,
#ifdef P4_TO_P8
                            iter_edge,
#endif
                            iter_corner);
      }
    }

#ifdef P4_TO_P8
    if (iter_edge != NULL || iter_corner != NULL) {
      for (i = 0; i < 12; i++) {
        edge_args.level = 0;
        e_edges[0] = i;
        if (num_edges > 0) {
          this_edge = tte[t * 12 + i];
        }
        else {
          this_edge = -1;
        }
        modulus = 0;
        if (this_edge > -1) {
          for (jt = ett_offset[this_edge];
               jt < ett_offset[this_edge + 1]; jt++) {
            if (edge_to_tree[jt] == t) {
              modulus = (int) edge_to_edge[jt];
              if (modulus % 12 == i) {
                modulus /= 12;
                break;
              }
            }
          }
        }
        edge_args.num_sides = 1;
        e_treeids[0] = t;
        if (modulus == 0) {
          e_common_corners[0][0] = p8est_edge_corners[i][0];
          e_common_corners[1][0] = p8est_edge_corners[i][1];
        }
        else {
          e_common_corners[0][0] = p8est_edge_corners[i][1];
          e_common_corners[1][0] = p8est_edge_corners[i][0];
        }
        for (side = left; side <= right; side++) {
          face[0] = p8est_edge_faces[i][side];
          nt = ttt[t * 6 + face[0]];
          face[1] = (int) ttf[t * 6 + face[0]];
          orientation = face[1] / 6;
          face[1] %= 6;
          if (nt != t || face[1] != face[0]) {
            perm = p8est_face_permutation_refs[face[0]][face[1]];
            perm = p8est_face_permutation_sets[perm][orientation];
            corner0 = p8est_face_permutations[perm]
              [p8est_edge_face_corners[i][face[0]][0]];
            corner1 = p8est_face_permutations[perm]
              [p8est_edge_face_corners[i][face[0]][1]];
            k = edge_args.num_sides;
            if (modulus == 0) {
              e_common_corners[0][k] = p8est_face_corners[face[1]][corner0];
              e_common_corners[1][k] = p8est_face_corners[face[1]][corner1];
            }
            else {
              e_common_corners[0][k] = p8est_face_corners[face[1]][corner1];
              e_common_corners[1][k] = p8est_face_corners[face[1]][corner0];
            }
            ne = e_common_corners[1][k] - e_common_corners[0][k];
            ne = (ne > 0) ? ne : -ne;
            ne = p8est_corner_edges[e_common_corners[0][k]][ne >> 1];
            e_edges[k] = ne;
            e_treeids[k] = nt;
            edge_args.num_sides++;
          }
        }
        if (this_edge > -1) {
          for (jt = ett_offset[this_edge];
               jt < ett_offset[this_edge + 1]; jt++) {
            nt = edge_to_tree[jt];
            ne = (int) edge_to_edge[jt];
            orientation = ne / 12;
            ne %= 12;
            if (orientation == 0) {
              corner0 = p8est_edge_corners[ne][0];
              corner1 = p8est_edge_corners[ne][1];
            }
            else {
              corner0 = p8est_edge_corners[ne][1];
              corner1 = p8est_edge_corners[ne][0];
            }
            for (k = 0; k < edge_args.num_sides; k++) {
              if (nt == e_treeids[k]) {
                if (corner0 == e_common_corners[0][k] &&
                    corner1 == e_common_corners[1][k]) {
                  break;
                }
              }
            }
            if (k == edge_args.num_sides) {
              e_treeids[k] = nt;
              e_edges[k] = ne;
              e_common_corners[0][k] = corner0;
              e_common_corners[1][k] = corner1;
              edge_args.num_sides++;
            }
          }
        }
        for (j = 1; j < edge_args.num_sides; j++) {
          if (e_treeids[j] > t) {
            break;
          }
          if (e_treeids[j] == t && e_edges[j] > i) {
            break;
          }
        }
        if (j == edge_args.num_sides) {
          for (j = 0; j < edge_args.num_sides; j++) {
            edge_args.start_idx2[j] = 0;
            index[j * 2 + local][0] = 0;
            tree = p4est_array_index_topidx (trees, e_treeids[j]);
            quadrants[j * 2 + local] = &(tree->quadrants);
            index[j * 2 + local][1] = quadrants[j * 2 + local]->elem_count;
            index[j * 2 + ghost][0] = tree_first_ghost[e_treeids[j]];
            index[j * 2 + ghost][1] = tree_first_ghost[e_treeids[j] + 1];
          }
          common_corners[0].elem_count = (size_t) edge_args.num_sides;
          common_corners[1].elem_count = (size_t) edge_args.num_sides;
          edges.elem_count = (size_t) edge_args.num_sides;
          edge_treeids.elem_count = (size_t) edge_args.num_sides;
          edge_quadids.elem_count = (size_t) edge_args.num_sides;
          edge_quads.elem_count = (size_t) edge_args.num_sides;

          if (iter_corner != NULL) {
            corner_args.num_sides = edge_args.num_sides * 2;
            corners.elem_count = (size_t) corner_args.num_sides;
            for (j = 0; j < corner_args.num_sides; j++) {
              c_treeids[j] = e_treeids[j % edge_args.num_sides];
              quadrants[2 * j] = quadrants[2 * (j % edge_args.num_sides)];
            }
            corner_treeids.elem_count = (size_t) corner_args.num_sides;
            corner_quadids.elem_count = (size_t) corner_args.num_sides;
            corner_quads.elem_count = (size_t) corner_args.num_sides;
          }

          p8est_edge_iterate (&edge_args, user_data, iter_edge, iter_corner);
        }
      }
    }
#endif

    if (iter_corner != NULL) {
#ifndef P4_TO_P8
      sc_array_init (cta, sizeof (p4est_corner_transform_t));
#endif
      for (i = 0; i < P4EST_CHILDREN; i++) {
        corner_args.level = 0;
        c_treeids[0] = t;
        c_corners[0] = i;
        corner_args.num_sides = 1;

        for (j = 0; j < P4EST_DIM; j++) {
          face[0] = p4est_corner_faces[i][j];
          nt = ttt[t * P4EST_DIM * 2 + face[0]];
          face[1] = (int) ttf[t * P4EST_DIM * 2 + face[0]];
          orientation = face[1] / (P4EST_DIM * 2);
          face[1] = face[1] % (P4EST_DIM * 2);
          if (nt == t && face[0] == face[1]) {
            continue;
          }
          c_treeids[corner_args.num_sides] = nt;
          for (k = 0; k < P4EST_CHILDREN / 2; k++) {
            if (p4est_face_corners[face[0]][k] == i) {
#ifndef P4_TO_P8
              rface[0] = face[0];
              rface[1] = face[1];
              face[0] = p4est_rface_to_zface[rface[0]];
              face[1] = p4est_rface_to_zface[rface[1]];
              modulus = face[0] ^ face[1];
              modulus = ((modulus & 2) >> 1) ^ (modulus & 1) ^ 1;
              face[0] = rface[0];
              face[1] = rface[1];
              c_corners[corner_args.num_sides] =
                p4est_face_corners[face[1]][modulus ^ orientation ^ k];
#else
              modulus = p8est_face_permutation_refs[face[0]][face[1]];
              modulus = p8est_face_permutation_sets[modulus][orientation];
              c_corners[corner_args.num_sides] = p8est_face_corners[face[1]]
                [p8est_face_permutations[modulus][k]];
#endif
            }
          }
          corner_args.num_sides++;
        }

#ifdef P4_TO_P8
        if (num_edges > 0) {
          for (j = 0; j < P4EST_DIM; j++) {
            m = p8est_corner_edges[i][j];
            this_edge = tte[t * 12 + m];
            if (this_edge == -1) {
              continue;
            }
            modulus = 0;
            for (kt = ett_offset[this_edge];
                 kt < ett_offset[this_edge + 1]; kt++) {
              if (edge_to_tree[kt] == t) {
                modulus = (int) edge_to_edge[kt];
                if (modulus % 12 == m) {
                  modulus /= 12;
                  break;
                }
              }
            }
            if (p8est_edge_corners[m][modulus] == i) {
              modulus = 0;
            }
            else {
              modulus = 1;
            }
            for (kt = ett_offset[this_edge];
                 kt < ett_offset[this_edge + 1]; kt++) {
              nt = edge_to_tree[kt];
              ne = (int) edge_to_edge[kt];
              orientation = ne / 12;
              ne %= 12;
              nc = p8est_edge_corners[ne][orientation ^ modulus];
              for (l = 0; l < corner_args.num_sides; l++) {
                if (nt == c_treeids[l] && nc == c_corners[l]) {
                  break;
                }
              }
              if (l == corner_args.num_sides) {
                c_treeids[l] = nt;
                c_corners[l] = nc;
                corner_args.num_sides++;
              }
            }
          }
        }
#endif

        if (num_corners > 0) {
#ifndef P4_TO_P8
          this_corner = ttc[t * P4EST_CHILDREN + p4est_corner_to_zorder[i]];
#else
          this_corner = ttc[t * P4EST_CHILDREN + i];
#endif
        }
        else {
          this_corner = -1;
        }

        if (this_corner != -1) {
#ifndef P4_TO_P8
          p4est_find_corner_transform (conn, t, p4est_corner_to_zorder[i],
                                       cta);
          for (j = 0; j < (int) cta->elem_count; j++) {
            ct = sc_array_index_int (cta, j);
            nt = ct->ntree;
            nc = (int) ct->ncorner;
            for (k = 0; k < corner_args.num_sides; k++) {
              if (nt == c_treeids[k] && nc == c_corners[k]) {
                break;
              }
            }
            if (k == corner_args.num_sides) {
              c_treeids[k] = nt;
              c_corners[k] = nc;
              corner_args.num_sides++;
            }
          }
#else
          for (jt = ctt_offset[this_corner];
               jt < ctt_offset[this_corner + 1]; jt++) {
            nt = corner_to_tree[jt];
            nc = (int) corner_to_corner[jt];
            for (k = 0; k < corner_args.num_sides; k++) {
              if (nt == c_treeids[k] && nc == c_corners[k]) {
                break;
              }
            }
            if (k == corner_args.num_sides) {
              c_treeids[k] = nt;
              c_corners[k] = nc;
              corner_args.num_sides++;
            }
          }
#endif
        }

        for (j = 0; j < corner_args.num_sides; j++) {
          if (c_treeids[j] > t) {
            break;
          }
          if (c_treeids[j] == t && c_corners[j] > i) {
            break;
          }
        }
        if (j == corner_args.num_sides) {
          for (j = 0; j < corner_args.num_sides; j++) {
            corner_args.start_idx2[j] = 0;
            index[j * 2 + local][0] = 0;
            tree = p4est_array_index_topidx (trees, c_treeids[j]);
            quadrants[j * 2 + local] = &(tree->quadrants);
            index[j * 2 + local][1] = quadrants[j * 2 + local]->elem_count;
            index[j * 2 + ghost][0] = tree_first_ghost[c_treeids[j]];
            index[j * 2 + ghost][1] = tree_first_ghost[c_treeids[j] + 1];
          }
          corners.elem_count = (size_t) corner_args.num_sides;
          corner_treeids.elem_count = (size_t) corner_args.num_sides;
          corner_quadids.elem_count = (size_t) corner_args.num_sides;
          corner_quads.elem_count = (size_t) corner_args.num_sides;

          p4est_corner_iterate (&corner_args, user_data, iter_corner);
        }
      }
#ifndef P4_TO_P8
      sc_array_reset (cta);
#endif
    }

  }

  if (iter_corner != NULL) {
    P4EST_FREE (c_corners);
    P4EST_FREE (c_treeids);
    P4EST_FREE (c_quadids);
    P4EST_FREE (c_quads);
  }

#ifdef P4_TO_P8
  if (iter_edge != NULL || iter_corner != NULL) {
    P4EST_FREE (e_common_corners[0]);
    P4EST_FREE (e_common_corners[1]);
    P4EST_FREE (e_treeids);
    P4EST_FREE (e_quads);
    P4EST_FREE (e_quadids);
    P4EST_FREE (e_edges);
  }
#endif

  P4EST_FREE (tree_first_ghost);
  P4EST_FREE (level_num);
  for (i = 0; i < alloc_size; i++) {
    P4EST_FREE (index[i]);
  }
  P4EST_FREE (index);
  P4EST_FREE (refine);
  P4EST_FREE (quad_idx2);
  P4EST_FREE (test_level);
  P4EST_FREE (count);
  P4EST_FREE (test);
  P4EST_FREE (first_index);
  P4EST_FREE (start_idx2);
  P4EST_FREE (quadrants);
}
