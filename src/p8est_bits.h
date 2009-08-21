/*
  This file is part of p4est.
  p4est is a C library to manage a parallel collection of quadtrees and/or
  octrees.

  Copyright (C) 2007,2008 Carsten Burstedde, Lucas Wilcox.

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

#ifndef P8EST_BITS_H
#define P8EST_BITS_H

#include <p8est.h>

SC_EXTERN_C_BEGIN;

/** Prints one line with quadrant's x, y, z and level.
 */
void                p8est_quadrant_print (int log_priority,
                                          const p8est_quadrant_t * q);

/** Test if two quadrants have equal Morton indices.
 * \return true if \a q1 describes the same quadrant as \a q2.
 */
bool                p8est_quadrant_is_equal (const p8est_quadrant_t * q1,
                                             const p8est_quadrant_t * q2);

/** Test if two quadrants have equal Morton indices and the same tree id.
 * \return          true if \a q1 describes the same quadrant as \a q2
 *                  and the p.which_tree fields are equal.
 */
bool                p8est_quadrant_is_equal_piggy (const p8est_quadrant_t *
                                                   q1,
                                                   const p8est_quadrant_t *
                                                   q2);

/** Compare two quadrants in their Morton ordering.
 * Both quadrants must be valid.
 * \return Returns < 0 if \a v1 < \a v2,
 *                   0 if \a v1 == \a v2,
 *                 > 0 if \a v1 > \a v2
 */
int                 p8est_quadrant_compare (const void *v1, const void *v2);

/** Compare two quadrants in their Morton ordering and the which_tree member.
 * Both quadrants must be extended (superset of valid, see below).
 * \return Returns < 0 if \a v1 < \a v2,
 *                   0 if \a v1 == \a v2,
 *                 > 0 if \a v1 > \a v2
 */
int                 p8est_quadrant_compare_piggy (const void *v1,
                                                  const void *v2);

/** Test if two quadrants have equal Morton indices, callback version.
 * \return true if \a v1 describes the same quadrant as \a v2.
 */
bool                p8est_quadrant_equal_fn (const void *v1, const void *v2,
                                             const void *u);

/** Computes a hash value for a quadrant by the lookup3 method.
 */
unsigned            p8est_quadrant_hash_fn (const void *v, const void *u);

/** Test if two nodes are in the same tree and have equal Morton indices.
 * \param [in] v1   Pointer to a clamped or unclamped node, depending on u.
 * \param [in] v2   Pointer to a clamped or unclamped node, depending on u.
 * \param [in] u    User data, points to a bool holding the clamped flag.
 */
bool                p8est_node_equal_piggy_fn (const void *v1,
                                               const void *v2, const void *u);

/** Compute hash value of a node based on its tree and Morton index.
 * \param [in] v    Pointer to a clamped or unclamped node, depending on u.
 * \param [in] u    User data, points to a bool holding the clamped flag.
 */
unsigned            p8est_node_hash_piggy_fn (const void *v, const void *u);

/** Clamp a node inside the unit tree if it sits on a high border.
 * \param [in] n    Node to be clamped. Must not yet be clamped.
 * \param [out] r   Existing node overwritten by the clamped result.
 */
void                p8est_node_clamp_inside (const p8est_quadrant_t * n,
                                             p8est_quadrant_t * r);

/** Move a clamped node out on the border.
 * \param [in] n    Node to be unclamped in-place.
 */
void                p8est_node_unclamp (p8est_quadrant_t * n);

/** Find the enclosing quadrant of a given node at a given level.
 * \param [in] n        Clamped node.
 * \param [in] level    Level of the quadrant to be created.
 * \param [out] q       Output quadrant, n == q is permitted.
 */
void                p8est_node_to_quadrant (const p8est_quadrant_t * n,
                                            int level, p8est_quadrant_t * q);

/** Decide if a node is completely contained within a quadrant.
 * \param [in] q        Valid quadrant.
 * \param [in] n        Clamped node.
 */
bool                p8est_quadrant_contains_node (const p8est_quadrant_t * q,
                                                  const p8est_quadrant_t * n);

/** Compute the position of the ancestor of this child at level \a level within
 * its siblings.
 * \return Returns its child id in 0..7
 */
int                 p8est_quadrant_ancestor_id (const p8est_quadrant_t * q,
                                                const int level);

/** Compute the position of this child within its siblings.
 * \return Returns its child id in 0..7
 */
int                 p8est_quadrant_child_id (const p8est_quadrant_t * q);

/** Test if a quadrant is inside the unit tree.
 * \param [in] q Quadrant to be tested.
 * \return Returns true if \a q is inside the unit tree.
 */
bool                p8est_quadrant_is_inside_root (const p8est_quadrant_t *
                                                   q);

/** Test if a quadrant is inside the 3x3 box around the root tree.
 * \param [in] q Quadrant to be tested.
 * \return Returns true if \a q is inside the unit tree.
 */
bool                p8est_quadrant_is_inside_3x3 (const p8est_quadrant_t * q);

/** Test if a quadrant is outside a tree face boundary (no edge or corner).
 * \param [in] q Quadrant to be tested.
 * \return Returns true if \a q is outside across a unit tree face.
 */
bool                p8est_quadrant_is_outside_face (const p8est_quadrant_t *
                                                    q);

/** Test if a quadrant is outside a tree edge boundary (no corner).
 * \param [in] q Quadrant to be tested.
 * \return Returns true if \a q is outside across a unit tree edge.
 */
bool                p8est_quadrant_is_outside_edge (const p8est_quadrant_t *
                                                    q);

/** Test if a quadrant is outside a tree edge boundary (no corner).
 * \param [in] q       Quadrant to be tested.
 * \param [out] edge   The tree edge number is computed if outside edge.
 *                     This pointer may be NULL.
 * \return Returns true if \a q is outside across a unit tree edge.
 */
bool                p8est_quadrant_is_outside_edge_extra (const
                                                          p8est_quadrant_t *
                                                          q, int *edge);

/** Test if a quadrant is outside a tree corner boundary.
 * \param [in] q Quadrant to be tested.
 * \return Returns true if \a q is outside across a unit tree corner.
 */
bool                p8est_quadrant_is_outside_corner (const p8est_quadrant_t *
                                                      q);

/** Test if a quadrant is used to represent a mesh node.
 * \param [in] q        Quadrant to be tested.
 * \param [in] inside   If true, boundary nodes must be clamped inside.
 *                      If false, nodes must align with the quadrant grid.
 * \return Returns true if \a q is a node.
 */
bool                p8est_quadrant_is_node (const p8est_quadrant_t * q,
                                            bool inside);

/** Test if a quadrant has valid Morton indices and is inside the unit tree.
 * \param [in] q Quadrant to be tested.
 * \return Returns true if \a q is valid.
 */
bool                p8est_quadrant_is_valid (const p8est_quadrant_t * q);

/** Test if a quadrant has valid Morton indices in the 3x3 box around root.
 * \param [in] q Quadrant to be tested.
 * \return Returns true if \a q is extended.
 */
bool                p8est_quadrant_is_extended (const p8est_quadrant_t * q);

/** Test if two quadrants are siblings.
 * \param [in] q1 First quadrant to be tested.
 * \param [in] q2 Second quadrant to be tested.
 * \return true if \a q1 is unequal to and a sibling of \a q2.
 */
bool                p8est_quadrant_is_sibling (const p8est_quadrant_t * q1,
                                               const p8est_quadrant_t * q2);

/** Test if two quadrants are siblings.
 * Descriptive, slower version of \a p8est_quadrant_is_sibling.
 * For debugging and educational purposes only.
 */
bool                p8est_quadrant_is_sibling_D (const p8est_quadrant_t * q1,
                                                 const p8est_quadrant_t * q2);

/** Test if 8 quadrants are siblings in Morton ordering.
 */
bool                p8est_quadrant_is_family (const p8est_quadrant_t * q0,
                                              const p8est_quadrant_t * q1,
                                              const p8est_quadrant_t * q2,
                                              const p8est_quadrant_t * q3,
                                              const p8est_quadrant_t * q4,
                                              const p8est_quadrant_t * q5,
                                              const p8est_quadrant_t * q6,
                                              const p8est_quadrant_t * q7);

/** Test if 8 quadrants are siblings in Morton ordering, array version.
 * \param [in] q   Array of 8 quadrants.
 */
bool                p8est_quadrant_is_familyv (const p8est_quadrant_t q[]);

/** Test if 8 quadrants are siblings in Morton ordering, array version.
 * \param [in] q   Array of 8 pointers to quadrants.
 */
bool                p8est_quadrant_is_familypv (p8est_quadrant_t * q[]);

/** Test if a quadrant it the parent of another quadrant.
 * \param [in] q Quadrant to be tested.
 * \param [in] r Possible child quadrant.
 * \return true if \a q is the parent of \a r.
 */
bool                p8est_quadrant_is_parent (const p8est_quadrant_t * q,
                                              const p8est_quadrant_t * r);

/** Test if a quadrant it the parent of another quadrant.
 * Descriptive, slower version of \a p8est_quadrant_is_parent.
 * For debugging and educational purposes only.
 */
bool                p8est_quadrant_is_parent_D (const p8est_quadrant_t * q,
                                                const p8est_quadrant_t * r);

/** Test if a quadrant is an ancestor of another quadrant.
 * \param [in] q Quadrant to be tested.
 * \param [in] r Descendent quadrant.
 * \return true if \a q is unequal to and an ancestor of \a r.
 */
bool                p8est_quadrant_is_ancestor (const p8est_quadrant_t * q,
                                                const p8est_quadrant_t * r);

/** Test if a quadrant is an ancestor of another quadrant.
 * Descriptive, slower version of \a p8est_quadrant_is_ancestor.
 * Contrary to \a p8est_quadrant_is_ancestor, it aborts for inter-tree q, r.
 * For debugging and educational purposes only.
 */
bool                p8est_quadrant_is_ancestor_D (const p8est_quadrant_t * q,
                                                  const p8est_quadrant_t * r);

/** Test if two quadrants follow each other in the tree with no holes.
 * \param [in] q A quadrant
 * \param [in] r Another quadrant
 * \return true if \a q is immediately before \a r in the tree.
 * \note for every \a q there are between 0 and P8EST_MAXLEVEL+1 possible nexts.
 */
bool                p8est_quadrant_is_next (const p8est_quadrant_t * q,
                                            const p8est_quadrant_t * r);

/** Test if two quadrants follow each other in the tree with no holes.
 * Descriptive, slower version of \a p8est_quadrant_is_next.
 * For debugging and educational purposes only.
 */
bool                p8est_quadrant_is_next_D (const p8est_quadrant_t * q,
                                              const p8est_quadrant_t * r);

/** Test if a quadrant has at least partial overlap with a tree.
 */
bool                p8est_quadrant_overlaps_tree (p8est_tree_t * tree,
                                                  const p8est_quadrant_t * q);

/** Test if a quadrant is completely contained within a tree.
 */
bool                p8est_quadrant_is_inside_tree (p8est_tree_t * tree,
                                                   const p8est_quadrant_t *
                                                   q);

/** Compute the parent of a quadrant.
 * \param [in]  q Input quadrant.
 * \param [in,out] r Existing quadrant whose Morton index will be filled
 *                   with the Morton index of the parent of \a q.
 *                   Its user_data will be untouched.
 * \note \a q may point to the same quadrant as \a r.
         The user_data of \a r is never modified.
 */
void                p8est_quadrant_parent (const p8est_quadrant_t * q,
                                           p8est_quadrant_t * r);

/** Compute a specific sibling of a quadrant.
 * \param [in]     q  Input quadrant.
 * \param [in,out] r  Existing quadrant whose Morton index will be filled
 *                    with the coordinates of sibling no. sibling_id of q.
 * \param [in]     sibling_id The id of the sibling computed, 0..3.
 */
void                p8est_quadrant_sibling (const p8est_quadrant_t * q,
                                            p8est_quadrant_t * r,
                                            int sibling_id);

/** Compute the face neighbor of a quadrant.
 * \param [in]     q      Input quadrant, must be valid.
 * \param [in]     face   The face across which to generate the neighbor.
 * \param [in,out] r      Existing quadrant whose Morton index will be filled.
 * \note \a q may point to the same quadrant as \a r.
 */
void                p8est_quadrant_face_neighbor (const p8est_quadrant_t * q,
                                                  int face,
                                                  p8est_quadrant_t * r);

/** Compute the face neighbor of a quadrant, transforming across tree
 * boundaries if necessary.
 * \param [in]     q      Input quadrant, must be valid.
 * \param [in]     t      Tree that contains \q.
 * \param [in]     face   The face across which to generate the neighbor.
 * \param [in,out] r      Existing quadrant whose Morton index will be filled.
 *                        By convention, if there is no tree across \face,
 *                        \r has the same Morton index as \q.
 * \param [in]     conn   The connectivity structure for the forest.
 * \return Returns the tree that contains \r.  By convention, if there is no
 * tree across \face, then -1 is returned.
 */
p4est_locidx_t      p8est_quadrant_face_neighbor_extra (const p8est_quadrant_t
                                                        * q, p4est_locidx_t t,
                                                        int face,
                                                        p8est_quadrant_t * r,
                                                        p8est_connectivity_t *
                                                        conn);

/** Get the smaller face neighbors of \a q.
 *
 * Gets the smaller face neighbors, which are half of the size assuming the
 * 2-1 constant.
 *
 * The order of the \a n[i] is given in the Morton ordering.
 *
 * \param [in]  q      The quadrant whose face neighbors will be constructed.
 * \param [in]  face   The face across which to generate the neighbors.
 * \param [out] n[0]..n[3] Filled with the four smaller face neighbors.
 * \param [out] nur[0]..nur[3] If not NULL, filled with smallest quadrants
 *                     that fit in the upper right corners of \a n.
 */
void                p8est_quadrant_half_face_neighbors (const p8est_quadrant_t
                                                        * q, int face,
                                                        p8est_quadrant_t n[],
                                                        p8est_quadrant_t
                                                        nur[]);

/** Create all possible face neighbors of \a q.
 *
 * Gets the face neighbors, possible assuming the 2-1 constraint.
 * If the larger or smaller quadrants do not exist than they are returned
 * as initialized by P4EST_QUADRANT_INIT.
 *
 * The order of \a n[0] through \a n[3] are given in Morton ordering.
 *
 * \param [in]  q      The quadrant whose face neighbors will be constructed.
 * \param [in]  face   The face across which to generate the neighbors.
 * \param [out] n[0]..n[3] Filled with the smaller possible face neighbors,
 *                     which are half of the size if they exist
 *                     or initialized to P4EST_QUADRANT_INIT.
 * \param [out] n[4]   Filled with the face neighbor, which is the same size.
 * \param [out] n[5]   Filled with the face neighbor, which is twice the size
 *                     if it exists or initialized to P4EST_QUADRANT_INIT.
 */
void                p8est_quadrant_all_face_neighbors (const p8est_quadrant_t
                                                       * q, int face,
                                                       p8est_quadrant_t n[]);

/** Compute the edge neighbor of a quadrant.
 * \param [in]     q      Input quadrant, must be valid.
 * \param [in]     edge   The edge across which to generate the neighbor.
 * \param [in,out] r      Existing quadrant whose Morton index will be filled.
 * \note \a q may point to the same quadrant as \a r.
 */
void                p8est_quadrant_edge_neighbor (const p8est_quadrant_t * q,
                                                  int edge,
                                                  p8est_quadrant_t * r);

/** Compute the edge neighbors of a quadrant, transforming across tree
 * boundaries if necessary.  Only computes neighbors that are not face
 * neighbors.
 * \param [in]     q      Input quadrant, must be valid.
 * \param [in]     t      Tree that contains \q.
 * \param [in]     edge   The edge across which to generate the neighbor.
 * \param [in,out] quads  An initialized but empty array where the edge
 *                        neighbors will be placed.
 * \param [in,out] treeids An initialized but empty array where the ids of the
 *                        trees containing the edge neighbors will be placed.
 * \param [in]     conn   The connectivity structure for the forest.
 */
void                p8est_quadrant_edge_neighbor_extra (const p8est_quadrant_t
                                                        * q, p4est_locidx_t t,
                                                        int edge, sc_array_t *
                                                        quads, sc_array_t *
                                                        treeids,
                                                        p8est_connectivity_t *
                                                        conn);

/** Compute the corner neighbor of a quadrant.
 * \param [in]     q      Input quadrant, must be valid.
 * \param [in]     corner The corner across which to generate the neighbor.
 * \param [in,out] r      Existing quadrant whose Morton index will be filled.
 * \note \a q may point to the same quadrant as \a r.
 */
void                p8est_quadrant_corner_neighbor (const p8est_quadrant_t *
                                                    q, int corner,
                                                    p8est_quadrant_t * r);

/** Compute the corner neighbors of a quadrant, transforming across tree
 * boundaries if necessary.  Only computes neighbors that are not face or edge
 * neighbors.
 * \param [in]     q      Input quadrant, must be valid.
 * \param [in]     t      Tree that contains \q.
 * \param [in]     corner The corner across which to generate the neighbor.
 * \param [in,out] quads  An initialized but empty array where the corner
 *                        neighbors will be placed.
 * \param [in,out] treeids An initialized but empty array where the ids of the
 *                        trees containing the corner neighbors will be placed.
 * \param [in]     conn   The connectivity structure for the forest.
 */
void                p8est_quadrant_corner_neighbor_extra (const
                                                          p8est_quadrant_t *
                                                          q, p4est_locidx_t t,
                                                          int corner,
                                                          sc_array_t * quads,
                                                          sc_array_t *
                                                          treeids,
                                                          p8est_connectivity_t
                                                          * conn);

/** Compute the corner node of a quadrant.
 * \param [in]     q      Input quadrant, must be valid.
 * \param [in]     corner The corner across which to generate the neighbor.
 * \param [in,out] r      Node that will not be clamped inside.
 * \note \a q may point to the same quadrant as \a r.
 */
void                p8est_quadrant_corner_node (const p8est_quadrant_t * q,
                                                int corner,
                                                p8est_quadrant_t * r);

/** Compute the 8 children of a quadrant.
 * \param [in]     q  Input quadrant.
 * \param [in,out] c0 First computed child.
 *                    \a q may point to the same quadrant as \a c0.
 * \note The user_data of \a c0, c1, c2, c3, c4, c5, c6, c7 is never modified.
 */
void                p8est_quadrant_children (const p8est_quadrant_t * q,
                                             p8est_quadrant_t * c0,
                                             p8est_quadrant_t * c1,
                                             p8est_quadrant_t * c2,
                                             p8est_quadrant_t * c3,
                                             p8est_quadrant_t * c4,
                                             p8est_quadrant_t * c5,
                                             p8est_quadrant_t * c6,
                                             p8est_quadrant_t * c7);

/** Compute the 8 children of a quadrant, array version.
 * \param [in]     q  Input quadrant.
 * \param [in,out] c  The 8 computed children in z-order.
 *                    q may point to the same quadrant as c[0].
 * \note The user_data of c[i] is never modified.
 */
void                p8est_quadrant_childrenv (const p8est_quadrant_t * q,
                                              p8est_quadrant_t c[]);

/** Compute the first descendent of a quadrant on a given level.
 * \param [in]  q      Input quadrant.
 * \param [out] fd     First descendent of \a q on level \a level.
 * \param [in]  level  Level must be greater equal than q's level.
 */
void                p8est_quadrant_first_descendent (const p8est_quadrant_t *
                                                     q, p8est_quadrant_t * fd,
                                                     int level);

/** Compute the last descendent of a quadrant on a given level.
 * \param [in]  q      Input quadrant.
 * \param [out] ld     Last descendent of \a q on level \a level.
 * \param [in]  level  Level must be greater equal than q's level.
 */
void                p8est_quadrant_last_descendent (const p8est_quadrant_t *
                                                    q, p8est_quadrant_t * ld,
                                                    int level);

/** Computes the nearest common ancestor of two quadrants in the same tree.
 * \param [in]     q1 First input quadrant.
 * \param [in]     q2 Second input quadrant.
 * \param [in,out] r Existing quadrant whose Morton index will be filled.
 *                   Its user_data will be untouched.
 * \note \a q1, \a q2, \a r may point to the same quadrant.
 *       The user_data of \a r is never modified.
 */
void                p8est_nearest_common_ancestor (const p8est_quadrant_t *
                                                   q1,
                                                   const p8est_quadrant_t *
                                                   q2, p8est_quadrant_t * r);

/** Computes the nearest common ancestor of two quadrants in the same tree.
 * Descriptive, slower version of \a p8est_nearest_common_ancestor.
 * For debugging and educationial purposes only.
 */
void                p8est_nearest_common_ancestor_D (const p8est_quadrant_t *
                                                     q1,
                                                     const p8est_quadrant_t *
                                                     q2,
                                                     p8est_quadrant_t * r);

/** Transforms a quadrant/node across a face between trees.
 * \param [in]     q        Input quadrant/non-clamped node.
 * \param [in,out] r        Quadrant/node whose Morton index will be filled.
 * \param [in] ftransform   This array holds 9 integers.
 *             [0]..[2]     The coordinate axis sequence of the origin face.
 *             [3]..[5]     The coordinate axis sequence of the target face.
 *             [6]..[8]     Edge reverse flag for axes 0, 1; face code for 2.
 * \note \a q and \q r may NOT point to the same quadrant structure.
 */
void                p8est_quadrant_transform_face (const p8est_quadrant_t * q,
                                                   p8est_quadrant_t * r,
                                                   const int ftransform[]);

/** Checks if a quadrant touches an edge (diagonally inside or outside).
 */
bool                p8est_quadrant_touches_edge (const p8est_quadrant_t * q,
                                                 int edge, bool inside);

/** Transforms a quadrant across an edge between trees.
 * \param [in]     q          Input quadrant.
 * \param [in,out] r          Quadrant whose Morton index will be filled.
 * \param [in]     edge       Edge index of the originating quadrant.
 * \param [in]     ei         Edge information computed previously.
 * \param [in]     inside     The quadrant will be placed inside or outside.
 */
void                p8est_quadrant_transform_edge (const p8est_quadrant_t * q,
                                                   p8est_quadrant_t * r,
                                                   const p8est_edge_info_t *
                                                   ei,
                                                   const
                                                   p8est_edge_transform_t *
                                                   et, bool inside);

/** Shifts a quadrant until it touches the specified edge from the inside.
 * \param [in]     q          Valid input quadrant.
 * \param [in,out] r          Quadrant whose Morton index will be filled.
 * \param [in]     edge       Edge index.
 */
void                p8est_quadrant_shift_edge (const p8est_quadrant_t * q,
                                               p8est_quadrant_t * r,
                                               int edge);

/** Checks if a quadrant touches a corner (diagonally inside or outside).
 */
bool                p8est_quadrant_touches_corner (const p8est_quadrant_t * q,
                                                   int corner, bool inside);

/** Move a quadrant inside or diagonally outside a corner position.
 * \param [in,out] q        This quadrant only requires a valid level.
 * \param [in]     icorner  Number of the corner in 0..7.
 * \param [int]    inside   Boolean flag for inside or diagonally outside.
 */
void                p8est_quadrant_transform_corner (p8est_quadrant_t * r,
                                                     int corner, bool inside);

/** Shifts a quadrant until it touches the specified corner from the inside.
 * \param [in]     q          Valid input quadrant.
 * \param [in,out] r          Quadrant whose Morton index will be filled.
 * \param [in]     corner     Corner index.
 */
void                p8est_quadrant_shift_corner (const p8est_quadrant_t * q,
                                                 p8est_quadrant_t * r,
                                                 int corner);

/** Computes the linear position of a quadrant in a uniform grid.
 * \param [in] quadrant  Quadrant whose id will be computed.
 * \return Returns the linear position of this quadrant on a grid.
 * \note This is the inverse operation of p8est_quadrant_set_morton.
 *       The user_data of \a quadrant is never modified.
 */
uint64_t            p8est_quadrant_linear_id (const p8est_quadrant_t *
                                              quadrant, int level);

/** Set quadrant Morton indices based on linear position in uniform grid.
 * \param [in,out] quadrant  Quadrant whose Morton indices will be set.
 * \param [in]     id        The linear position of this quadrant on a grid.
 * \note This is the inverse operation of p8est_quadrant_linear_id.
 *       The user_data of \a quadrant is never modified.
 */
void                p8est_quadrant_set_morton (p8est_quadrant_t * quadrant,
                                               int level, uint64_t id);

SC_EXTERN_C_END;

#endif /* !P8EST_BITS_H */
