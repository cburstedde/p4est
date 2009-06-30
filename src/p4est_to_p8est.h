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

#ifdef P4EST_H
#error "The include files p4est.h and p4est_to_p8est.h cannot be combined"
#endif
#define P4_TO_P8

/* redefine macros */
#define P4EST_ONDISK_FORMAT             P8EST_ONDISK_FORMAT
#define P4EST_DIM                       P8EST_DIM
#define P4EST_CHILDREN                  P8EST_CHILDREN
#define P4EST_INSUL                     P8EST_INSUL
#define P4EST_STRING                    P8EST_STRING
#define P4EST_MAXLEVEL                  P8EST_MAXLEVEL
#define P4EST_QMAXLEVEL                 P8EST_QMAXLEVEL
#define P4EST_ROOT_LEN                  P8EST_ROOT_LEN
#define P4EST_QUADRANT_LEN              P8EST_QUADRANT_LEN
#define P4EST_LAST_OFFSET               P8EST_LAST_OFFSET
#define P4EST_QUADRANT_INIT             P8EST_QUADRANT_INIT

/* redefine enums */
#define P4EST_COMM_BALANCE_FIRST_COUNT  P8EST_COMM_BALANCE_FIRST_COUNT
#define P4EST_COMM_BALANCE_FIRST_LOAD   P8EST_COMM_BALANCE_FIRST_LOAD
#define P4EST_COMM_BALANCE_SECOND_COUNT P8EST_COMM_BALANCE_SECOND_COUNT
#define P4EST_COMM_BALANCE_SECOND_LOAD  P8EST_COMM_BALANCE_SECOND_LOAD
#define P4EST_COMM_PARTITION_GIVEN      P8EST_COMM_PARTITION_GIVEN
#define P4EST_COMM_PARTITION_WEIGHTED_LOW P8EST_COMM_PARTITION_WEIGHTED_LOW
#define P4EST_COMM_PARTITION_WEIGHTED_HIGH P8EST_COMM_PARTITION_WEIGHTED_HIGH
#define P4EST_COMM_GHOST_COUNT          P8EST_COMM_GHOST_COUNT
#define P4EST_COMM_GHOST_LOAD           P8EST_COMM_GHOST_LOAD
#define P4EST_COMM_NODES_QUERY          P8EST_COMM_NODES_QUERY
#define P4EST_COMM_NODES_REPLY          P8EST_COMM_NODES_REPLY
#define P4EST_COMM_SAVE                 P8EST_COMM_SAVE
#define P4EST_BALANCE_DEFAULT           P8EST_BALANCE_DEFAULT
#define P4EST_BALANCE_FULL              P8EST_BALANCE_FULL

/* redefine types */
#define p4est_connectivity_t            p8est_connectivity_t
#define p4est_geometry_t                p8est_geometry_t
#define p4est_balance_type_t            p8est_balance_type_t
#define p4est_t                         p8est_t
#define p4est_tree_t                    p8est_tree_t
#define p4est_quadrant_t                p8est_quadrant_t
#define p4est_position_t                p8est_position_t
#define p4est_init_t                    p8est_init_t
#define p4est_refine_t                  p8est_refine_t
#define p4est_coarsen_t                 p8est_coarsen_t
#define p4est_weight_t                  p8est_weight_t
#define p4est_indep_t                   p8est_indep_t
#define p4est_nodes_t                   p8est_nodes_t

/* redefine external variables */
#define p4est_face_corners              p8est_face_corners
#define p4est_face_dual                 p8est_face_dual
#define p4est_face_child_hang           p8est_face_child_hang
#define P4EST_DATA_UNINITIALIZED        P8EST_DATA_UNINITIALIZED
#define p4est_num_ranges                p8est_num_ranges

/* functions in p4est_connectivity */
#define p4est_connectivity_destroy      p8est_connectivity_destroy
#define p4est_connectivity_is_valid     p8est_connectivity_is_valid
#define p4est_connectivity_is_equal     p8est_connectivity_is_equal
#define p4est_connectivity_save         p8est_connectivity_save
#define p4est_connectivity_load         p8est_connectivity_load

/* functions in p4est */
#define p4est_new                       p8est_new
#define p4est_destroy                   p8est_destroy
#define p4est_copy                      p8est_copy
#define p4est_reset_data                p8est_reset_data
#define p4est_refine                    p8est_refine
#define p4est_refine_level              p8est_refine_level
#define p4est_coarsen                   p8est_coarsen
#define p4est_balance                   p8est_balance
#define p4est_partition                 p8est_partition
#define p4est_checksum                  p8est_checksum
#define p4est_save                      p8est_save
#define p4est_load                      p8est_load
#define p4est_balance_type_int          p8est_balance_type_int
#define p4est_balance_type_string       p8est_balance_type_string

/* functions in p4est_points */
#define p4est_new_points                p8est_new_points

/* functions in p4est_bits */
#define p4est_quadrant_print            p8est_quadrant_print
#define p4est_quadrant_is_equal         p8est_quadrant_is_equal
#define p4est_quadrant_is_equal_piggy   p8est_quadrant_is_equal_piggy
#define p4est_quadrant_compare          p8est_quadrant_compare
#define p4est_quadrant_compare_piggy    p8est_quadrant_compare_piggy
#define p4est_quadrant_equal_fn         p8est_quadrant_equal_fn
#define p4est_quadrant_hash_fn          p8est_quadrant_hash_fn
#define p4est_node_equal_piggy_fn       p8est_node_equal_piggy_fn
#define p4est_node_hash_piggy_fn        p8est_node_hash_piggy_fn
#define p4est_node_clamp_inside         p8est_node_clamp_inside
#define p4est_node_unclamp              p8est_node_unclamp
#define p4est_node_to_quadrant          p8est_node_to_quadrant
#define p4est_quadrant_contains_node    p8est_quadrant_contains_node
#define p4est_quadrant_higher_child_id  p8est_quadrant_higher_child_id
#define p4est_quadrant_child_id         p8est_quadrant_child_id
#define p4est_quadrant_is_inside_root   p8est_quadrant_is_inside_root
#define p4est_quadrant_is_inside_3x3    p8est_quadrant_is_inside_3x3
#define p4est_quadrant_is_outside_face  p8est_quadrant_is_outside_face
#define p4est_quadrant_is_outside_corner p8est_quadrant_is_outside_corner
#define p4est_quadrant_is_node          p8est_quadrant_is_node
#define p4est_quadrant_is_valid         p8est_quadrant_is_valid
#define p4est_quadrant_is_extended      p8est_quadrant_is_extended
#define p4est_quadrant_is_sibling       p8est_quadrant_is_sibling
#define p4est_quadrant_is_sibling_D     p8est_quadrant_is_sibling_D
#define p4est_quadrant_is_familyv       p8est_quadrant_is_familyv
#define p4est_quadrant_is_familypv      p8est_quadrant_is_familypv
#define p4est_quadrant_is_parent        p8est_quadrant_is_parent
#define p4est_quadrant_is_parent_D      p8est_quadrant_is_parent_D
#define p4est_quadrant_is_ancestor      p8est_quadrant_is_ancestor
#define p4est_quadrant_is_ancestor_D    p8est_quadrant_is_ancestor_D
#define p4est_quadrant_is_next          p8est_quadrant_is_next
#define p4est_quadrant_is_next_D        p8est_quadrant_is_next_D
#define p4est_quadrant_overlaps_tree    p8est_quadrant_overlaps_tree
#define p4est_quadrant_is_inside_tree   p8est_quadrant_is_inside_tree
#define p4est_quadrant_parent           p8est_quadrant_parent
#define p4est_quadrant_sibling          p8est_quadrant_sibling
#define p4est_quadrant_face_neighbor    p8est_quadrant_face_neighbor
#define p4est_quadrant_half_face_neighbors p8est_quadrant_half_face_neighbors
#define p4est_quadrant_all_face_neighbors p8est_quadrant_all_face_neighbors
#define p4est_quadrant_corner_neighbor  p8est_quadrant_corner_neighbor
#define p4est_quadrant_corner_node      p8est_quadrant_corner_node
#define p4est_quadrant_childrenv        p8est_quadrant_childrenv
#define p4est_quadrant_first_descendent p8est_quadrant_first_descendent
#define p4est_quadrant_last_descendent  p8est_quadrant_last_descendent
#define p4est_nearest_common_ancestor   p8est_nearest_common_ancestor
#define p4est_nearest_common_ancestor_D p8est_nearest_common_ancestor_D
#define p4est_quadrant_touches_corner   p8est_quadrant_touches_corner
#define p4est_quadrant_transform_corner p8est_quadrant_transform_corner
#define p4est_quadrant_shift_corner     p8est_quadrant_shift_corner
#define p4est_quadrant_linear_id        p8est_quadrant_linear_id
#define p4est_quadrant_set_morton       p8est_quadrant_set_morton

/* functions in p4est_algorithms */
#define p4est_quadrant_init_data        p8est_quadrant_init_data
#define p4est_quadrant_free_data        p8est_quadrant_free_data
#define p4est_quadrant_checksum         p8est_quadrant_checksum
#define p4est_tree_is_sorted            p8est_tree_is_sorted
#define p4est_tree_is_linear            p8est_tree_is_linear
#define p4est_tree_is_almost_sorted     p8est_tree_is_almost_sorted
#define p4est_tree_is_complete          p8est_tree_is_complete
#define p4est_tree_print                p8est_tree_print
#define p4est_is_equal                  p8est_is_equal
#define p4est_is_valid                  p8est_is_valid
#define p4est_split_array               p8est_split_array
#define p4est_find_lower_bound          p8est_find_lower_bound
#define p4est_find_higher_bound         p8est_find_higher_bound
#define p4est_tree_compute_overlap      p8est_tree_compute_overlap
#define p4est_tree_uniqify_overlap      p8est_tree_uniqify_overlap
#define p4est_tree_remove_nonowned      p8est_tree_remove_nonowned
#define p4est_complete_region           p8est_complete_region
#define p4est_complete_subtree          p8est_complete_subtree
#define p4est_balance_subtree           p8est_balance_subtree
#define p4est_linearize_tree            p8est_linearize_tree
#define p4est_partition_given           p8est_partition_given

/* functions in p4est_communication */
#define p4est_comm_count_quadrants      p8est_comm_count_quadrants
#define p4est_comm_global_partition     p8est_comm_global_partition
#define p4est_comm_find_owner           p8est_comm_find_owner
#define p4est_comm_tree_info            p8est_comm_tree_info
#define p4est_comm_neighborhood_owned   p8est_comm_neighborhood_owned
#define p4est_comm_sync_flag            p8est_comm_sync_flag

/* functions in p4est_vtk */
#define p4est_vtk_write_file            p8est_vtk_write_file
#define p4est_vtk_write_footer          p8est_vtk_write_footer

/* functions in p4est_ghost */
#define p4est_quadrant_find_owner       p8est_quadrant_find_owner
#define p4est_build_ghost_layer         p8est_build_ghost_layer
#define p4est_face_quadrant_exists      p8est_face_quadrant_exists
#define p4est_quadrant_exists           p8est_quadrant_exists
#define p4est_is_balanced               p8est_is_balanced

/* functions in p4est_nodes */
#define p4est_nodes_new                 p8est_nodes_new
#define p4est_nodes_destroy             p8est_nodes_destroy
#define p4est_nodes_is_valid            p8est_nodes_is_valid

#endif /* !P4EST_TO_P8EST_H */
