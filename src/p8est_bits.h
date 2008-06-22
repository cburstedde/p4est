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

/** Test if two quadrants have equal Morton indices.
 * \return true if \a v1 describes the same quadrant as \a v2.
 */
bool                p8est_quadrant_is_equal (const void *v1, const void *v2);

/** Computes a hash value for a quadrant in 0..2^30-1.
 */
unsigned            p8est_quadrant_hash (const void *v);

/** Compute the position of this child within its siblings.
 * \return Returns its child id in 0..3
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

/** Test if 4 quadrants are siblings in Morton ordering, array version.
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

#if 0
/** Compute the level of balance needed at a specified corner.
 * \param [in]  zcorner  Corner index in z-ordering.
 * \return  Returns the maximum of level and this quadrants' corner level.
 */
int                 p8est_quadrant_corner_level (const p8est_quadrant_t * q,
                                                 int zcorner, int level);

/** Move a quadrant inside or diagonally outside a corner position.
 * \param [in,out] q        This quadrant only requires a valid level.
 * \param [in]     zcorner  Number of the corner in z-order, in 0..3.
 * \param [int]    inside   Boolean flag for inside or diagonally outside.
 */
void                p8est_quadrant_corner (p8est_quadrant_t * q,
                                           int zcorner, int inside);

/** Shift a quadrant by the size of a tree depending on the face.
 * \param [in,out] q     The quadrant to be modified.
 * \param [in]     face  Number of the face to move across, in 0..3.
 */
void                p8est_quadrant_translate (p8est_quadrant_t * q, int face);

/** Transforms the node of quadrant between trees.
 *
 * This gives the node of the transformed quadrant cooresponding to
 * the node passed in.
 *
 * \param [in]     node The node in z-order of the quadrant that is transformed.
 * \param [in] transform_type Transformation as in p8est_connectivity.h.
 * \return The transformed node number coresponding to \a node.
 */
int                 p8est_node_transform (int node, int transform_type);
#endif

/** Transforms a quadrant between trees.
 * \param [in]     q          Input quadrant.
 * \param [in,out] r          Quadrant whose Morton index will be filled.
 * \param [in] my_axis        The coordinate axis sequence of the origin face.
 * \param [in] target_axis    The coordinate axis sequence of the target face.
 * \param [in] edge_reverse   Edge reverse flag for axes 0, 1; face code for 2.
 * \note \a q and \q r may NOT point to the same quadrant structure.
 */
void                p8est_quadrant_transform (const p8est_quadrant_t * q,
                                              p8est_quadrant_t * r,
                                              int my_axis[3],
                                              int target_axis[3],
                                              int edge_reverse[3]);

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

#endif /* !P8EST_BITS_H */
