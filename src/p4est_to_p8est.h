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
#define P4EST_CHILDREN                  P8EST_CHILDREN
#define P4EST_MAXLEVEL                  P8EST_MAXLEVEL
#define P4EST_ROOT_LEN                  P8EST_ROOT_LEN
#define P4EST_QUADRANT_LEN              P8EST_QUADRANT_LEN
#define P4EST_LAST_OFFSET               P8EST_LAST_OFFSET
#define P4EST_QUADRANT_INIT             P8EST_QUADRANT_INIT
#define P4EST_COMM_BALANCE_FIRST_COUNT  P8EST_COMM_BALANCE_FIRST_COUNT
#define P4EST_COMM_BALANCE_FIRST_LOAD   P8EST_COMM_BALANCE_FIRST_LOAD
#define P4EST_COMM_BALANCE_SECOND_COUNT P8EST_COMM_BALANCE_SECOND_COUNT
#define P4EST_COMM_BALANCE_SECOND_LOAD  P8EST_COMM_BALANCE_SECOND_LOAD
#define P4EST_COMM_PARTITION_GIVEN      P8EST_COMM_PARTITION_GIVEN
#define P4EST_COMM_PARTITION_WEIGHTED_LOW P8EST_COMM_PARTITION_WEIGHTED_LOW
#define P4EST_COMM_PARTITION_WEIGHTED_HIGH P8EST_COMM_PARTITION_WEIGHTED_HIGH

/* redefine types */
#define p4est_t                         p8est_t
#define p4est_tree_t                    p8est_tree_t
#define p4est_quadrant_t                p8est_quadrant_t
#define p4est_connectivity_t            p8est_connectivity_t
#define p4est_position_t                p8est_position_t
#define p4est_init_t                    p8est_init_t
#define p4est_refine_t                  p8est_refine_t
#define p4est_coarsen_t                 p8est_coarsen_t
#define p4est_weight_t                  p8est_weight_t

/* redefine functions */
#define p4est_connectivity_destroy      p8est_connectivity_destroy
#define p4est_connectivity_is_valid     p8est_connectivity_is_valid
#define p4est_new                       p8est_new
#define p4est_destroy                   p8est_destroy
#define p4est_copy                      p8est_copy
#define p4est_refine                    p8est_refine
#define p4est_coarsen                   p8est_coarsen
#define p4est_partition                 p8est_partition
#define p4est_checksum                  p8est_checksum
#define p4est_quadrant_compare          p8est_quadrant_compare
#define p4est_quadrant_compare_piggy    p8est_quadrant_compare_piggy
#define p4est_quadrant_is_equal         p8est_quadrant_is_equal
#define p4est_quadrant_hash             p8est_quadrant_hash
#define p4est_quadrant_child_id         p8est_quadrant_child_id
#define p4est_quadrant_is_inside        p8est_quadrant_is_inside
#define p4est_quadrant_is_valid         p8est_quadrant_is_valid
#define p4est_quadrant_is_extended      p8est_quadrant_is_extended
#define p4est_quadrant_is_sibling       p8est_quadrant_is_sibling
#define p4est_quadrant_is_sibling_D     p8est_quadrant_is_sibling_D
/* p8est_quadrant_is_family is implemented in p8est_algorithms.c */
#define p4est_quadrant_is_parent        p8est_quadrant_is_parent
#define p4est_quadrant_is_parent_D      p8est_quadrant_is_parent_D
#define p4est_quadrant_is_ancestor      p8est_quadrant_is_ancestor
#define p4est_quadrant_is_ancestor_D    p8est_quadrant_is_ancestor_D
#define p4est_quadrant_is_next          p8est_quadrant_is_next
#define p4est_quadrant_is_next_D        p8est_quadrant_is_next_D
#define p4est_quadrant_parent           p8est_quadrant_parent
#define p4est_quadrant_sibling          p8est_quadrant_sibling
/* p8est_quadrant_children is implemented in p8est_algorithms.c */
#define p4est_quadrant_first_descendent p8est_quadrant_first_descendent
#define p4est_quadrant_last_descendent p8est_quadrant_last_descendent
#define p4est_nearest_common_ancestor   p8est_nearest_common_ancestor
#define p4est_nearest_common_ancestor_D p8est_nearest_common_ancestor_D
#define p4est_quadrant_linear_id        p8est_quadrant_linear_id
#define p4est_quadrant_set_morton       p8est_quadrant_set_morton
#define p4est_quadrant_init_data        p8est_quadrant_init_data
#define p4est_quadrant_free_data        p8est_quadrant_free_data
#define p4est_quadrant_print            p8est_quadrant_print
#define p4est_quadrant_checksum         p8est_quadrant_checksum
#define p4est_tree_is_sorted            p8est_tree_is_sorted
#define p4est_tree_is_linear            p8est_tree_is_linear
/* p4est_tree_is_almost_sorted not implemented yet */
#define p4est_tree_is_complete          p8est_tree_is_complete
#define p4est_tree_print                p8est_tree_print
#define p4est_is_valid                  p8est_is_valid
#define p4est_find_lower_bound          p8est_find_lower_bound
#define p4est_find_higher_bound         p8est_find_higher_bound
#define p4est_tree_uniqify_overlap      p8est_tree_uniqify_overlap
#define p4est_complete_region           p8est_complete_region
#define p4est_partition_given           p8est_partition_given
#define p4est_comm_count_quadrants      p8est_comm_count_quadrants
#define p4est_comm_global_partition     p8est_comm_global_partition
#define p4est_comm_find_owner           p8est_comm_find_owner

#endif /* !__P4EST_TO_P8EST_H */
