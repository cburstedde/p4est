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

#ifndef P4EST_TO_P8EST_H
#define P4EST_TO_P8EST_H

#ifdef P4EST_H
#error "The include files p4est.h and p4est_to_p8est.h cannot be combined"
#endif
#define P4_TO_P8

/* redefine macros */
#define P4EST_ONDISK_FORMAT             P8EST_ONDISK_FORMAT
#define P4EST_DIM                       P8EST_DIM
#define P4EST_FACES                     P8EST_FACES
#define P4EST_CHILDREN                  P8EST_CHILDREN
#define P4EST_HALF                      P8EST_HALF
#define P4EST_FTRANSFORM                P8EST_FTRANSFORM
#define P4EST_INSUL                     P8EST_INSUL
#define P4EST_STRING                    P8EST_STRING
#define P4EST_MAXLEVEL                  P8EST_MAXLEVEL
#define P4EST_QMAXLEVEL                 P8EST_QMAXLEVEL
#define P4EST_ROOT_LEN                  P8EST_ROOT_LEN
#define P4EST_QUADRANT_LEN              P8EST_QUADRANT_LEN
#define P4EST_LAST_OFFSET               P8EST_LAST_OFFSET
#define P4EST_QUADRANT_INIT             P8EST_QUADRANT_INIT

/* redefine enums */
#define P4EST_COMM_COUNT_PERTREE        P8EST_COMM_COUNT_PERTREE
#define P4EST_COMM_BALANCE_FIRST_COUNT  P8EST_COMM_BALANCE_FIRST_COUNT
#define P4EST_COMM_BALANCE_FIRST_LOAD   P8EST_COMM_BALANCE_FIRST_LOAD
#define P4EST_COMM_BALANCE_SECOND_COUNT P8EST_COMM_BALANCE_SECOND_COUNT
#define P4EST_COMM_BALANCE_SECOND_LOAD  P8EST_COMM_BALANCE_SECOND_LOAD
#define P4EST_COMM_PARTITION_GIVEN      P8EST_COMM_PARTITION_GIVEN
#define P4EST_COMM_PARTITION_WEIGHTED_LOW P8EST_COMM_PARTITION_WEIGHTED_LOW
#define P4EST_COMM_PARTITION_WEIGHTED_HIGH P8EST_COMM_PARTITION_WEIGHTED_HIGH
#define P4EST_COMM_PARTITION_CORRECTION P8EST_COMM_PARTITION_CORRECTION
#define P4EST_COMM_GHOST_COUNT          P8EST_COMM_GHOST_COUNT
#define P4EST_COMM_GHOST_LOAD           P8EST_COMM_GHOST_LOAD
#define P4EST_COMM_GHOST_EXCHANGE       P8EST_COMM_GHOST_EXCHANGE
#define P4EST_COMM_GHOST_EXPAND_COUNT   P8EST_COMM_GHOST_EXPAND_COUNT
#define P4EST_COMM_GHOST_EXPAND_LOAD    P8EST_COMM_GHOST_EXPAND_LOAD
#define P4EST_COMM_NODES_QUERY          P8EST_COMM_NODES_QUERY
#define P4EST_COMM_NODES_REPLY          P8EST_COMM_NODES_REPLY
#define P4EST_COMM_SAVE                 P8EST_COMM_SAVE
#define P4EST_COMM_LNODES_TEST          P8EST_COMM_LNODES_TEST
#define P4EST_COMM_LNODES_PASS          P8EST_COMM_LNODES_PASS
#define P4EST_COMM_LNODES_OWNED         P8EST_COMM_LNODES_OWNED
#define P4EST_COMM_LNODES_ALL           P8EST_COMM_LNODES_ALL
#define P4EST_CONNECT_FACE              P8EST_CONNECT_FACE
#define P4EST_CONNECT_CORNER            P8EST_CONNECT_CORNER
#define P4EST_CONNECT_DEFAULT           P8EST_CONNECT_DEFAULT
#define P4EST_CONNECT_FULL              P8EST_CONNECT_FULL
#define P4EST_CONN_ENCODE_NONE          P8EST_CONN_ENCODE_NONE
#define P4EST_WRAP_NONE                 P8EST_WRAP_NONE
#define P4EST_WRAP_REFINE               P8EST_WRAP_REFINE
#define P4EST_WRAP_COARSEN              P8EST_WRAP_COARSEN

/* redefine types */
#define p4est_connect_type_t            p8est_connect_type_t
#define p4est_connectivity_encode_t     p8est_connectivity_encode_t
#define p4est_connectivity_t            p8est_connectivity_t
#define p4est_corner_transform_t        p8est_corner_transform_t
#define p4est_corner_info_t             p8est_corner_info_t
#define p4est_geometry_t                p8est_geometry_t
#define p4est_t                         p8est_t
#define p4est_tree_t                    p8est_tree_t
#define p4est_quadrant_t                p8est_quadrant_t
#define p4est_inspect_t                 p8est_inspect_t
#define p4est_position_t                p8est_position_t
#define p4est_init_t                    p8est_init_t
#define p4est_refine_t                  p8est_refine_t
#define p4est_coarsen_t                 p8est_coarsen_t
#define p4est_weight_t                  p8est_weight_t
#define p4est_ghost_t                   p8est_ghost_t
#define p4est_indep_t                   p8est_indep_t
#define p4est_nodes_t                   p8est_nodes_t
#define p4est_lnodes_t                  p8est_lnodes_t
#define p4est_lnodes_code_t             p8est_lnodes_code_t
#define p4est_lnodes_rank_t             p8est_lnodes_rank_t
#define p4est_lnodes_buffer_t           p8est_lnodes_buffer_t
#define p4est_iter_volume_t             p8est_iter_volume_t
#define p4est_iter_volume_info_t        p8est_iter_volume_info_t
#define p4est_iter_face_t               p8est_iter_face_t
#define p4est_iter_face_info_t          p8est_iter_face_info_t
#define p4est_iter_face_side_t          p8est_iter_face_side_t
#define p4est_iter_corner_t             p8est_iter_corner_t
#define p4est_iter_corner_side_t        p8est_iter_corner_side_t
#define p4est_iter_corner_info_t        p8est_iter_corner_info_t
#define p4est_search_query_t            p8est_search_query_t
#define p4est_mesh_t                    p8est_mesh_t
#define p4est_mesh_face_neighbor_t      p8est_mesh_face_neighbor_t
#define p4est_wrap_t                    p8est_wrap_t
#define p4est_wrap_leaf_t               p8est_wrap_leaf_t
#define p4est_wrap_flags_t              p8est_wrap_flags_t

/* redefine external variables */
#define p4est_face_corners              p8est_face_corners
#define p4est_face_dual                 p8est_face_dual
#define p4est_corner_faces              p8est_corner_faces
#define p4est_corner_face_corners       p8est_corner_face_corners
#define p4est_child_corner_faces        p8est_child_corner_faces
#define P4EST_DATA_UNINITIALIZED        P8EST_DATA_UNINITIALIZED
#define p4est_num_ranges                p8est_num_ranges

/* functions in p4est_connectivity */
#define p4est_connectivity_memory_used  p8est_connectivity_memory_used
#define p4est_connectivity_new          p8est_connectivity_new
#define p4est_connectivity_new_brick    p8est_connectivity_new_brick
#define p4est_connectivity_new_byname   p8est_connectivity_new_byname
#define p4est_connectivity_new_copy     p8est_connectivity_new_copy
#define p4est_connectivity_destroy      p8est_connectivity_destroy
#define p4est_connectivity_set_attr     p8est_connectivity_set_attr
#define p4est_connectivity_is_valid     p8est_connectivity_is_valid
#define p4est_connectivity_is_equal     p8est_connectivity_is_equal
#define p4est_connectivity_sink         p8est_connectivity_sink
#define p4est_connectivity_deflate      p8est_connectivity_deflate
#define p4est_connectivity_save         p8est_connectivity_save
#define p4est_connectivity_source       p8est_connectivity_source
#define p4est_connectivity_inflate      p8est_connectivity_inflate
#define p4est_connectivity_load         p8est_connectivity_load
#define p4est_connectivity_complete     p8est_connectivity_complete
#define p4est_expand_face_transform     p8est_expand_face_transform
#define p4est_find_face_transform       p8est_find_face_transform
#define p4est_find_corner_transform     p8est_find_corner_transform
#define p4est_corner_array_index        p8est_corner_array_index
#define p4est_connectivity_reorder      p8est_connectivity_reorder
#define p4est_connectivity_permute      p8est_connectivity_permute
#define p4est_connectivity_join_faces   p8est_connectivity_join_faces
#define p4est_connectivity_is_equivalent p8est_connectivity_is_equivalent

/* functions in p4est */
#define p4est_qcoord_to_vertex          p8est_qcoord_to_vertex
#define p4est_memory_used               p8est_memory_used
#define p4est_new                       p8est_new
#define p4est_destroy                   p8est_destroy
#define p4est_copy                      p8est_copy
#define p4est_reset_data                p8est_reset_data
#define p4est_refine                    p8est_refine
#define p4est_coarsen                   p8est_coarsen
#define p4est_balance                   p8est_balance
#define p4est_partition                 p8est_partition
#define p4est_checksum                  p8est_checksum
#define p4est_save                      p8est_save
#define p4est_load                      p8est_load
#define p4est_connect_type_int          p8est_connect_type_int
#define p4est_connect_type_string       p8est_connect_type_string
#define p4est_tree_array_index          p8est_tree_array_index
#define p4est_quadrant_array_index      p8est_quadrant_array_index
#define p4est_quadrant_array_push       p8est_quadrant_array_push
#define p4est_quadrant_mempool_alloc    p8est_quadrant_mempool_alloc
#define p4est_quadrant_list_pop         p8est_quadrant_list_pop

/* functions in p4est_extended */
#define p4est_replace_t                 p8est_replace_t
#define p4est_new_ext                   p8est_new_ext
#define p4est_refine_ext                p8est_refine_ext
#define p4est_coarsen_ext               p8est_coarsen_ext
#define p4est_balance_ext               p8est_balance_ext
#define p4est_balance_subtree_ext       p8est_balance_subtree_ext
#define p4est_partition_ext             p8est_partition_ext
#define p4est_save_ext                  p8est_save_ext
#define p4est_load_ext                  p8est_load_ext

/* functions in p4est_iterate */
#define p4est_iterate                   p8est_iterate
#define p4est_iterate_ext               p8est_iterate_ext
#define p4est_iter_fside_array_index    p8est_iter_fside_array_index
#define p4est_iter_fside_array_index_int p8est_iter_fside_array_index_int
#define p4est_iter_cside_array_index    p8est_iter_cside_array_index
#define p4est_iter_cside_array_index_int p8est_iter_cside_array_index_int

/* functions in p4est_points */
#define p4est_new_points                p8est_new_points

/* functions in p4est_bits */
#define p4est_quadrant_print            p8est_quadrant_print
#define p4est_quadrant_is_equal         p8est_quadrant_is_equal
#define p4est_quadrant_overlaps         p8est_quadrant_overlaps
#define p4est_quadrant_is_equal_piggy   p8est_quadrant_is_equal_piggy
#define p4est_quadrant_compare          p8est_quadrant_compare
#define p4est_quadrant_disjoint         p8est_quadrant_disjoint
#define p4est_quadrant_compare_piggy    p8est_quadrant_compare_piggy
#define p4est_quadrant_equal_fn         p8est_quadrant_equal_fn
#define p4est_quadrant_hash_fn          p8est_quadrant_hash_fn
#define p4est_node_equal_piggy_fn       p8est_node_equal_piggy_fn
#define p4est_node_hash_piggy_fn        p8est_node_hash_piggy_fn
#define p4est_node_clamp_inside         p8est_node_clamp_inside
#define p4est_node_unclamp              p8est_node_unclamp
#define p4est_node_to_quadrant          p8est_node_to_quadrant
#define p4est_quadrant_contains_node    p8est_quadrant_contains_node
#define p4est_quadrant_ancestor_id      p8est_quadrant_ancestor_id
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
#define p4est_quadrant_is_family        p8est_quadrant_is_family
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
#define p4est_quadrant_ancestor         p8est_quadrant_ancestor
#define p4est_quadrant_parent           p8est_quadrant_parent
#define p4est_quadrant_sibling          p8est_quadrant_sibling
#define p4est_quadrant_face_neighbor    p8est_quadrant_face_neighbor
#define p4est_quadrant_face_neighbor_extra p8est_quadrant_face_neighbor_extra
#define p4est_quadrant_half_face_neighbors p8est_quadrant_half_face_neighbors
#define p4est_quadrant_all_face_neighbors p8est_quadrant_all_face_neighbors
#define p4est_quadrant_corner_neighbor  p8est_quadrant_corner_neighbor
#define p4est_quadrant_corner_neighbor_extra \
                                        p8est_quadrant_corner_neighbor_extra
#define p4est_quadrant_half_corner_neighbor \
                                        p8est_quadrant_half_corner_neighbor
#define p4est_quadrant_corner_node      p8est_quadrant_corner_node
#define p4est_quadrant_children         p8est_quadrant_children
#define p4est_quadrant_childrenv        p8est_quadrant_childrenv
#define p4est_quadrant_childrenpv       p8est_quadrant_childrenpv
#define p4est_quadrant_first_descendant p8est_quadrant_first_descendant
#define p4est_quadrant_last_descendant  p8est_quadrant_last_descendant
#define p4est_quadrant_corner_descendant p8est_quadrant_corner_descendant
#define p4est_nearest_common_ancestor   p8est_nearest_common_ancestor
#define p4est_nearest_common_ancestor_D p8est_nearest_common_ancestor_D
#define p4est_quadrant_transform_face   p8est_quadrant_transform_face
#define p4est_quadrant_touches_corner   p8est_quadrant_touches_corner
#define p4est_quadrant_transform_corner p8est_quadrant_transform_corner
#define p4est_quadrant_shift_corner     p8est_quadrant_shift_corner
#define p4est_quadrant_linear_id        p8est_quadrant_linear_id
#define p4est_quadrant_set_morton       p8est_quadrant_set_morton

/* functions in p4est_search */
#define p4est_find_lower_bound          p8est_find_lower_bound
#define p4est_find_higher_bound         p8est_find_higher_bound
#define p4est_split_array               p8est_split_array
#define p4est_find_range_boundaries     p8est_find_range_boundaries
#define p4est_search                    p8est_search

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
#define p4est_tree_compute_overlap      p8est_tree_compute_overlap
#define p4est_tree_uniqify_overlap      p8est_tree_uniqify_overlap
#define p4est_tree_remove_nonowned      p8est_tree_remove_nonowned
#define p4est_complete_region           p8est_complete_region
#define p4est_complete_subtree          p8est_complete_subtree
#define p4est_balance_subtree           p8est_balance_subtree
#define p4est_balance_border            p8est_balance_border
#define p4est_linearize_tree            p8est_linearize_tree
#define p4est_next_nonempty_process     p8est_next_nonempty_process
#define p4est_partition_correction      p8est_partition_correction
#define p4est_partition_for_coarsening  p8est_partition_for_coarsening
#define p4est_partition_given           p8est_partition_given

/* functions in p4est_communication */
#define p4est_comm_count_quadrants      p8est_comm_count_quadrants
#define p4est_comm_global_partition     p8est_comm_global_partition
#define p4est_comm_count_pertree        p8est_comm_count_pertree
#define p4est_comm_is_owner             p8est_comm_is_owner
#define p4est_comm_find_owner           p8est_comm_find_owner
#define p4est_comm_tree_info            p8est_comm_tree_info
#define p4est_comm_neighborhood_owned   p8est_comm_neighborhood_owned
#define p4est_comm_sync_flag            p8est_comm_sync_flag
#define p4est_comm_checksum             p8est_comm_checksum

/* functions in p4est_io */
#define p4est_deflate_quadrants         p8est_deflate_quadrants
#define p4est_inflate                   p8est_inflate

/* functions in p4est_vtk */
#define p4est_vtk_write_file            p8est_vtk_write_file
#define p4est_vtk_write_all             p8est_vtk_write_all
#define p4est_vtk_write_header          p8est_vtk_write_header
#define p4est_vtk_write_point_scalar    p8est_vtk_write_point_scalar
#define p4est_vtk_write_point_vector    p8est_vtk_write_point_vector
#define p4est_vtk_write_footer          p8est_vtk_write_footer

/* functions in p4est_ghost */
#define p4est_quadrant_find_owner       p8est_quadrant_find_owner
#define p4est_ghost_memory_used         p8est_ghost_memory_used
#define p4est_ghost_new                 p8est_ghost_new
#define p4est_ghost_destroy             p8est_ghost_destroy
#define p4est_ghost_exchange_data       p8est_ghost_exchange_data
#define p4est_ghost_exchange_custom     p8est_ghost_exchange_custom
#define p4est_ghost_exchange_custom_levels p8est_ghost_exchange_custom_levels
#define p4est_ghost_bsearch             p8est_ghost_bsearch
#define p4est_ghost_contains            p8est_ghost_contains
#define p4est_face_quadrant_exists      p8est_face_quadrant_exists
#define p4est_quadrant_exists           p8est_quadrant_exists
#define p4est_is_balanced               p8est_is_balanced
#define p4est_ghost_checksum            p8est_ghost_checksum
#define p4est_ghost_expand              p8est_ghost_expand

/* functions in p4est_nodes */
#define p4est_nodes_new                 p8est_nodes_new
#define p4est_nodes_destroy             p8est_nodes_destroy
#define p4est_nodes_is_valid            p8est_nodes_is_valid

/* functions in p4est_lnodes */
#define p4est_lnodes_new                p8est_lnodes_new
#define p4est_lnodes_destroy            p8est_lnodes_destroy
#define p4est_lnodes_decode             p8est_lnodes_decode
#define p4est_lnodes_share_owned_begin  p8est_lnodes_share_owned_begin
#define p4est_lnodes_share_owned_end    p8est_lnodes_share_owned_end
#define p4est_lnodes_share_owned        p8est_lnodes_share_owned
#define p4est_lnodes_share_all_begin    p8est_lnodes_share_all_begin
#define p4est_lnodes_share_all_end      p8est_lnodes_share_all_end
#define p4est_lnodes_share_all          p8est_lnodes_share_all
#define p4est_lnodes_buffer_destroy     p8est_lnodes_buffer_destroy
#define p4est_lnodes_rank_array_index   p8est_lnodes_rank_array_index
#define p4est_lnodes_rank_array_index_int p8est_lnodes_rank_array_index_int
#define p4est_lnodes_global_index       p8est_lnodes_global_index

/* functions in p4est_mesh */
#define p4est_mesh_memory_used          p8est_mesh_memory_used
#define p4est_mesh_new                  p8est_mesh_new
#define p4est_mesh_destroy              p8est_mesh_destroy
#define p4est_mesh_quadrant_cumulative  p8est_mesh_quadrant_cumulative
#define p4est_mesh_face_neighbor_init   p8est_mesh_face_neighbor_init
#define p4est_mesh_face_neighbor_init2  p8est_mesh_face_neighbor_init2
#define p4est_mesh_face_neighbor_next   p8est_mesh_face_neighbor_next

/* functions in p4est_balance */
#define p4est_balance_seeds_face        p8est_balance_seeds_face
#define p4est_balance_seeds_corner      p8est_balance_seeds_corner
#define p4est_balance_seeds             p8est_balance_seeds

/* functions in p4est_wrap */
#define p4est_wrap_new_conn             p8est_wrap_new_conn
#define p4est_wrap_new_world            p8est_wrap_new_world
#define p4est_wrap_destroy              p8est_wrap_destroy
#define p4est_wrap_get_ghost            p8est_wrap_get_ghost
#define p4est_wrap_get_mesh             p8est_wrap_get_mesh
#define p4est_wrap_mark_refine          p8est_wrap_mark_refine
#define p4est_wrap_mark_coarsen         p8est_wrap_mark_coarsen
#define p4est_wrap_adapt                p8est_wrap_adapt
#define p4est_wrap_partition            p8est_wrap_partition
#define p4est_wrap_complete             p8est_wrap_complete
#define p4est_wrap_leaf_next            p8est_wrap_leaf_next
#define p4est_wrap_leaf_first           p8est_wrap_leaf_first

#endif /* !P4EST_TO_P8EST_H */
