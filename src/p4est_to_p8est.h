/*
  This file is part of p4est.
  p4est is a C library to manage a parallel collection of quadtrees and/or
  octrees.

  Copyright (C) 2008 Carsten Burstedde, Lucas Wilcox.

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

#ifndef P4EST_TO_P8EST_H
#define P4EST_TO_P8EST_H

#define P4_TO_P8

/* redefine macros */
#define P4EST_DIM                       P8EST_DIM
#define P4EST_MAXLEVEL                  P8EST_MAXLEVEL
#define P4EST_ROOT_LEN                  P8EST_ROOT_LEN
#define P4EST_QUADRANT_LEN              P8EST_QUADRANT_LEN
#define P4EST_QUADRANT_INIT             P8EST_QUADRANT_INIT

/* redefine types */
#define p4est_t                         p8est_t
#define p4est_tree_t                    p8est_tree_t
#define p4est_quadrant_t                p8est_quadrant_t
#define p4est_connectivity_t            p8est_connectivity_t
#define p4est_position_t                p8est_position_t
#define p4est_init_t                    p8est_init_t

/* redefine functions */
#define p4est_new                       p8est_new
#define p4est_is_valid                  p8est_is_valid
#define p4est_connectivity_is_valid     p8est_connectivity_is_valid
#define p4est_quadrant_compare          p8est_quadrant_compare
#define p4est_quadrant_is_valid         p8est_quadrant_is_valid
#define p4est_quadrant_is_inside        p8est_quadrant_is_inside
#define p4est_quadrant_is_extended      p8est_quadrant_is_extended
#define p4est_quadrant_is_equal         p8est_quadrant_is_equal
#define p4est_quadrant_hash             p8est_quadrant_hash
#define p4est_quadrant_child_id         p8est_quadrant_child_id
#define p4est_quadrant_linear_id        p8est_quadrant_linear_id
#define p4est_quadrant_set_morton       p8est_quadrant_set_morton
#define p4est_quadrant_init_data        p8est_quadrant_init_data
#define p4est_quadrant_free_data        p8est_quadrant_free_data
#define p4est_complete_region           p8est_complete_region
#define p4est_comm_count_quadrants      p8est_comm_count_quadrants
#define p4est_comm_global_partition     p8est_comm_global_partition
#define p4est_comm_find_owner           p8est_comm_find_owner

#endif /* !__P4EST_TO_P8EST_H */
