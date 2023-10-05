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

/** \file p4est3_quadrant_vtable.h
 *
 * Virtual function mechanism for operations on quadrants.
 *
 * We use the term quadrants for the nodes of an adaptive quadtree or octree.
 * There are multiple operations on quadrants required for a working forest.
 * We abstract these functions using a virtual table to allow to exchange
 * various implementations of the quadrant functionality.
 *
 * For convenience, we provide both a 2D and a 3D implementation of the
 * quadrant definition used up to version 2 to make it available to version 3.
 *
 * \ingroup p4est3
 */

#ifndef P4EST3_QUADRANT_VTABLE_H
#define P4EST3_QUADRANT_VTABLE_H

#include <sc3_array.h>
#include <p4est3_base.h>

/** Normalize coordinates of a quadrant to maxlevel of P4EST3_REF_MAXLEVEL.
 * This is done to make the coordinates of various virtual quadrant
 * implementations comparable with each other.
 */
#define P4EST3_REF_MAXLEVEL 31

#ifdef __cplusplus
extern              "C"
{
#if 0
}
#endif
#endif

/* *INDENT-OFF* */

/*** Generic prototypes for quadrant functions -- debugging ***/

/** Generic prototype according to the \c sc3_object_is* query format.
 * The _is_ functions carry a \a reason parameter to facilitate debugging.
 */
typedef int         (*p4est3_quadrant_is_t) (const void * q, char *reason);

/** Generic prototype according to the \c sc3_object_is2* query format.
 * The _is2_ functions carry a \a reason parameter to facilitate debugging.
 */
typedef int         (*p4est3_quadrant_is2_t) (const void * q1, const void *q2,
                                              char *reason);

/*** Generic prototypes for quadrant functions -- production ***/

/** Generic prototype to compute a global index from an integer input. */
typedef p4est3_gloidx (*p4est3_quadrant_glo_t) (int level);

/** Generic prototype to extract one int from a quadrant. */
typedef sc3_error_t *(*p4est3_quadrant_in_j_t) (const void *q, int *j);
/** Generic prototype to take and extract one int from a quadrant. */
typedef sc3_error_t *(*p4est3_quadrant_in_i_j_t) (const void *q, int i, int *j);
/** Generic prototype to take two quadrants inputs and output an int. */
typedef sc3_error_t *(*p4est3_quadrant_in2_j_t) (const void * q1,
                                                 const void * q2, int *j);
/** Prototype to take a quadrant and output an array. */
typedef sc3_error_t *(*p4est3_quadrant_in_arr_t) (const void *q,
                                                  sc3_array_t * a);

/** Generic prototype to set a quadrant. */
typedef sc3_error_t *(*p4est3_quadrant_out_t) (void *r);
/** Generic prototype to take a quadrant and output another. */
typedef sc3_error_t *(*p4est3_quadrant_in_out_t) (const void *q, void *r);
/** Generic prototype to take a quadrant and int and output another value. */
typedef sc3_error_t *(*p4est3_quadrant_in_i_out_t) (const void *q, int i,
                                                    void *r);
/** Prototype to set a quadrant based on a linear index. */
typedef sc3_error_t *(*p4est3_quadrant_morton_t) (int l, p4est3_gloidx i,
                                                  void *r);

/** Prototype to set a common nearest ancestor of two quadrants. */
typedef sc3_error_t *(*p4est3_nearest_common_ancestor_t) (const void * q1,
                                                          const void * q2,
                                                          void *r);
/** Prototype to set a linear index of a quadrant based on its Morton index. */
typedef sc3_error_t *(*p4est3_quadrant_linear_id_t) (const void *q, int l,
                                                     p4est3_gloidx *i);

/** Prototype to construct the face neighbor of a quadrant accross a tree. */
typedef sc3_error_t *(*p4est3_quadrant_tree_face_neighbor_t) (const void *q,
                                                              sc3_array_t * t,
                                                              int i, void *r);

/*** Specific prototypes for quadrant query functions ***/

/** Prototype to query the number of quadrants for a uniform refinement. */
typedef p4est3_quadrant_glo_t p4est3_quadrant_num_uniform_t;
/** Prototype to query the level of a quadrant. */
typedef p4est3_quadrant_in_j_t p4est3_quadrant_level_t;
/** Prototype to query the child id of a quadrant. */
typedef p4est3_quadrant_in_j_t p4est3_quadrant_child_id_t;
/** Prototype to query the ancestor id of a quadrant. */
typedef p4est3_quadrant_in_i_j_t p4est3_quadrant_ancestor_id_t;
/** Prototype to compare two quadrants by linear index. */
typedef p4est3_quadrant_in2_j_t p4est3_quadrant_compare_t;
/** Prototype to query the one quadrant is ancetor of another. */
typedef p4est3_quadrant_in2_j_t p4est3_quadrant_is_ancestor_t;
/** Prototype to query the one quadrant is parent of another. */
typedef p4est3_quadrant_in2_j_t p4est3_quadrant_is_parent_t;

/*** Specific prototypes for quadrant creation functions ***/

/** Prototype to constract the root quadrant of a tree. */
typedef p4est3_quadrant_out_t p4est3_quadrant_root_t;
/** Prototype to copy one quadrant to another. */
typedef p4est3_quadrant_in_out_t p4est3_quadrant_copy_t;
/** Prototype to construct the parent of a quadrant. */
typedef p4est3_quadrant_in_out_t p4est3_quadrant_parent_t;
/** Prototype to construct the sibling of a quadrant. */
typedef p4est3_quadrant_in_i_out_t p4est3_quadrant_sibling_t;
/** Prototype to compute the predecessor of a quadrant according to linear index. */
typedef p4est3_quadrant_in_out_t p4est3_quadrant_predecessor_t;
/** Prototype to compute the successor of a quadrant according to linear index. */
typedef p4est3_quadrant_in_out_t p4est3_quadrant_successor_t;
/** Prototype to compute an ancestor of a quadrant (up the tree). */
typedef p4est3_quadrant_in_i_out_t p4est3_quadrant_ancestor_t;
/** Prototype to compute a child of a quadrant (down the tree). */
typedef p4est3_quadrant_in_i_out_t p4est3_quadrant_child_t;
/** Prototype to construct the first descendant at maximum level of a quadrant. */
typedef p4est3_quadrant_in_i_out_t p4est3_quadrant_first_descendant_t;
/** Prototype to construct the last descendant at maximum level of a quadrant. */
typedef p4est3_quadrant_in_i_out_t p4est3_quadrant_last_descendant_t;
/** Prototype to construct one quadrant from coordinates and level */
typedef p4est3_quadrant_in_i_out_t p4est3_quadrant_quadrant_t;

/** Prototype to query if a quadrant touches a given tree face boundary */
typedef p4est3_quadrant_in_i_j_t p4est3_quadrant_get_tree_boundary_t;
/** Prototype to query if a quadrant touches tree face boundaries and which */
typedef p4est3_quadrant_in_arr_t p4est3_quadrant_tree_boundaries_t;
/** Prototype to construct the face neighbor of a quadrant. */
typedef p4est3_quadrant_in_i_out_t p4est3_quadrant_face_neighbor_t;

/* *INDENT-ON* */

/** We abstract the operations on quadrants to allow for multiple implementations.
 * The virtual table must be populated by the members defining the implementation.
 * It is then passed to \ref p4est3_set_quadrant_vtable and its contents copied.
 *
 * The functions placed in this table are supposted to return NULL on
 * successful operation and an \c sc3_error_t object otherwise.
 * They are expected to assert the validity of their inputs.
 */
typedef struct p4est3_quadrant_vtable
{
  /* Constant numbers are defined as integers for simplicity. */
  const char         *name;     /**< String identifier of implementation. */
  int                 dim;      /**< Spatial dimension. */
  int                 max_level;        /**< Maximum level to be reached. */
  size_t              quadrant_size;    /**< Quadrant object size in bytes. */

  /* functions that are independent of any quadrant input argument */

  p4est3_quadrant_num_uniform_t quadrant_num_uniform;   /**< Quadrants per level. */

  /* functions that create a quadrant from equivalent information */

  p4est3_quadrant_root_t quadrant_root; /**< Generate root quadrant. */
  p4est3_quadrant_morton_t quadrant_morton;     /**< Generate by linear index. */
  p4est3_quadrant_quadrant_t quadrant_quadrant; /**< Generate from coordinates. */

  /* functions that take one (constant) quadrant input argument */

  /** Examine validity.  Pointer may be NULL, in which case validity is true. */
  p4est3_quadrant_is_t quadrant_is_valid;

  /* information about quadrants */
  p4est3_quadrant_level_t quadrant_level;               /**< Query the level. */
  p4est3_quadrant_child_id_t quadrant_child_id;         /**< Query child id. */
  p4est3_quadrant_ancestor_id_t quadrant_ancestor_id;   /**< Query ancestor id. */
  p4est3_quadrant_in_out_t quadrant_coordinates;        /**< Query coordinates. */
  p4est3_quadrant_linear_id_t quadrant_linear_id;       /**< Query Morton index. */

  /** Query tree boundary in a specific direction. */
  p4est3_quadrant_get_tree_boundary_t quadrant_get_tree_boundary;
  p4est3_quadrant_in_arr_t quadrant_tree_boundaries;    /**< Query tree boundary. */

  /** Deep copy one quadrant to another.
   * This pointer may be NULL, in which case we memcpy (3) the quadrant. */
  p4est3_quadrant_copy_t quadrant_copy;
  p4est3_quadrant_child_t quadrant_child;       /**< Generate a child by number. */
  p4est3_quadrant_sibling_t quadrant_sibling;   /**< Generate sibling */
  p4est3_quadrant_parent_t quadrant_parent;     /**< Generate parent. */
  p4est3_quadrant_ancestor_t quadrant_ancestor; /**< Generate ancestor by level. */
  p4est3_quadrant_predecessor_t quadrant_predecessor;   /**< Generate predecessor. */
  p4est3_quadrant_successor_t quadrant_successor;       /**< Generate successor. */

  /** Generate first smallest descendant of a quadrant at \a max_level. */
  p4est3_quadrant_first_descendant_t quadrant_first_descendant;
  /** Generate last smallest descendant of a quadrant at \a max_level. */
  p4est3_quadrant_last_descendant_t quadrant_last_descendant;

  p4est3_quadrant_face_neighbor_t quadrant_face_neighbor; /**< Generate face neighbor */
  /** Generate face neighbor across a tree boundary*/
  p4est3_quadrant_tree_face_neighbor_t quadrant_tree_face_neighbor;

  /* functions that take two (constant) quadrant input arguments */

  p4est3_quadrant_compare_t quadrant_compare;   /**< Compare linear indices. */

  /** Examine equality.
   * Pointer may be NULL, in which case we memcmp (3) the contents. */
  p4est3_quadrant_is2_t quadrant_is_equal;

  /**< Query if a quadrant is a ancestor of another. */
  p4est3_quadrant_is_ancestor_t quadrant_is_ancestor;

  /**< Query if a quadrant is a parent of another. */
  p4est3_quadrant_is_parent_t quadrant_is_parent;

  p4est3_nearest_common_ancestor_t nearest_common_ancestor;
}
const p4est3_quadrant_vtable_t;

/** Query validity of a quadrant virtual table.
 * \param [in] qvt      NULL is allowed and considered not valid.
 * \param [out] reason  May be NULL.  Otherwise, set to "" on output when valid
 *                      or, as applicable, the reason for not being valid.
 * \return              True if valid, false otherwise.
 */
int                 p4est3_quadrant_vtable_is_valid (const p4est3_quadrant_vtable_t * qvt, char *reason);

/** Return user-defined id of the implementation.
 * \param [in] qvt  Valid virtual quadrant table.
 * \return          User-defined id if \a qvt valid, negative number otherwise.
 */
int                 p4est3_quadrant_id (const p4est3_quadrant_vtable_t * qvt);

/** Return spatial dimension of the implementation.
 * \param [in] qvt  Valid virtual quadrant table.
 * \return          Dimension if \a qvt valid, negative number otherwise.
 */
int                 p4est3_quadrant_dim (const p4est3_quadrant_vtable_t * qvt);

/** Return maximum refinement level of the implementation.
 * \param [in] qvt  Valid virtual quadrant table.
 * \return          Maximum level if \a qvt valid, negative number otherwise.
 */
int                 p4est3_quadrant_max_level (const p4est3_quadrant_vtable_t *
                                               qvt);

/** Return number of children of a quadrant in this implementation.
 * \param [in] qvt  Valid virtual quadrant table.
 * \return          Number of children if \a qvt valid,
 *                  negative number otherwise.
 */
int                 p4est3_quadrant_num_children (const p4est3_quadrant_vtable_t *
                                                  qvt);

/** Return memory size in bytes of a quadrant in this implementation.
 * \param [in] qvt  Valid virtual quadrant table.
 * \return          Memory size if \a qvt valid, zero otherwise.
 */
size_t              p4est3_quadrant_size (const p4est3_quadrant_vtable_t * qvt);

/** Query number of quadrants for a uniform refinement at a given level.
 * \param [in] qvt      Valid virtual quadrant table.
 * \param [in] level    Valid level for this implementation.
 * \return              0 if \a qvt or \a level invalid, positive result otherwise.
 */
p4est3_gloidx       p4est3_quadrant_num_uniform (const p4est3_quadrant_vtable_t *
                                                 qvt, int level);

/** Query validity of a quadrant in the style of \c sc3_<object>_is2_valid.
 * \param [in] qvt      NULL is allowed and considered not valid.
 * \param [in] q        NULL is allowed and considered not valid.
 * \param [out] reason  May be NULL.  Otherwise, set to "" on output when valid
 *                      or, as applicable, the reason for not being valid.
 * \return              True if the quadrant \a q is valid, false otherwise.
 */
int                 p4est3_quadrant_vtable_is2_valid
  (const p4est3_quadrant_vtable_t * qvt, const void *q, char *reason);

/** Query equality of two quadrants in the style of \c sc3_<object>_is3_valid.
 * \param [in] qvt      NULL is allowed and considered not equal.
 * \param [in] q1       NULL is allowed and considered not equal.
 * \param [in] q2       NULL is allowed and considered not equal.
 * \param [out] reason  May be NULL.  Otherwise, set to "" on output when valid
 *                      or, as applicable, the reason for not being valid.
 * \return              True if the quadrant objects pointed to by \a q1 and
 *                      \a q2 are equal, false otherwise.
 */
int                 p4est3_quadrant_is3_equal (const p4est3_quadrant_vtable_t * qvt,
                                               const void *q1, const void *q2,
                                               char *reason);

/** Query is a quadrant touches a specific boundary of a tree it is contained in
 * in the style of \c sc3_<object>_is3_valid.
 * \param [in] qvt      Valid virtual quadrant table.
 * \param [in] q        Valid quadrant in this implementation.
 * \param [in] i        A index of a tree's face for that the intersection
 *                      with the quadrant \a q will be checked.
 *                      For 3D: 0..5. For 2D: 0..3.
 * \param [out] j       True if the quadrant \a q touches the \a i-th boundary
 *                      face of the tree, false otherwise.
* \return               NULL on success, error object otherwise.
 */
sc3_error_t
  * p4est3_quadrant_get_tree_boundary (const p4est3_quadrant_vtable_t * qvt,
                                       const void *q, int i, int *j);

/** Query if a quadrant touches a tree face boundaries and which if so 
 * \param [in] qvt      Valid virtual quadrant table.
 * \param [in] q        Valid quadrant in this implementation.
 * \param [out] nf      Array of size dimension of \a qvt. Every element
 *                      corresponds to a spatial direction and contains
 *                      a face number of a quadrant touching the tree boundary
 *                      in this direction.  If quadrant touches no boundaries,
 *                      then the element of array is filled by -1.  If quadrant
 *                      touches all the boundaries (i.e., the quadrant is
 *                      the whole tree), then the entry is filled by -2.
 * \return              NULL on success, error object otherwise.
*/
sc3_error_t        *p4est3_quadrant_tree_boundaries (const p4est3_quadrant_vtable_t
                                                     * qvt, const void *q,
                                                     sc3_array_t * nf);

/** Query the refinement level of a quadrant.
 * \param [in] qvt      Valid virtual quadrant table.
 * \param [in] q        Valid quadrant in this implementation.
 * \param [out] l       Non-NULL reference assigned level on output.
 * \return              NULL on success, error object otherwise.
 */
sc3_error_t        *p4est3_quadrant_level (const p4est3_quadrant_vtable_t * qvt,
                                           const void *q, int *l);

/** Query the child id of a quadrant relative to its parent.
 * \param [in] qvt      Valid virtual quadrant table.
 * \param [in] q        Valid quadrant in this implementation.
 * \param [out] j       Non-NULL reference assigned child id on output.
 * \return              NULL on success, error object otherwise.
 */
sc3_error_t        *p4est3_quadrant_child_id (const p4est3_quadrant_vtable_t * qvt,
                                              const void *q, int *j);

/** Query the child id of a quadrant relative to a given ancestor.
 * \param [in] qvt      Valid virtual quadrant table.
 * \param [in] q        Valid quadrant in this implementation.
 * \param [in] l        Non-negative level less or equal that of \a q.
 * \param [out] j       Non-NULL reference assigned ancestor id on output.
 * \return              NULL on success, error object otherwise.
 */
sc3_error_t        *p4est3_quadrant_ancestor_id (const p4est3_quadrant_vtable_t *
                                                 qvt, const void *q, int l,
                                                 int *j);

/** Query the integer coordinates of a quadrant with respect to the unit tree
 * with a maximum level equal to P4EST_MAXLEVEL.
 * \param [in] qvt      Valid virtual quadrant table.
 * \param [in] q        Valid quadrant in this implementation.
 * \param [out] j       Output array of \a P4EST_DIM coordinates.
 *                      We assume that the number of bits in an int suffices.
 * \return              NULL on success, error object otherwise.
 */
sc3_error_t        *p4est3_quadrant_coordinates (const p4est3_quadrant_vtable_t *
                                                 qvt, const void *q, void *j);

/** Construct a quadrant assigning its coordinates and level.
 * \param [in] qvt      Valid virtual quadrant table.
 * \param [in] c        Array of size dimension of \a qvt storing values
 *                      for quadrant coordinates.
 * \param [in] l        Desired level of the quadrant.
 * \param [in] q        Constructed quadrant is placed here.
 * \return              NULL on success, error object otherwise.
*/
sc3_error_t        *p4est3_quadrant_quadrant (const p4est3_quadrant_vtable_t * qvt,
                                              const void *c, int l, void *q);

/** Translate a quadrant from one representation to another.
 * \param [in] vtold    Valid virtual table of the quadrant to translate.
 * \param [in] qin      Valid quadrant to translate.
 * \param [in] vtnew    Valid virtual table of target quadrant representation.
 * \param [out] qout    Translated quadrant is plased here.
 * \return              NULL on success, error object otherwise.
 */
sc3_error_t        *p4est3_quadrant_translate (const p4est3_quadrant_vtable_t * vtold,
                                               const void *qin,
                                               const p4est3_quadrant_vtable_t * vtnew,
                                               void *qout);

/** Compare two quadrants by linear index.
 * \param [in] qvt      Valid virtual quadrant table.
 * \param [in] q1       Valid quadrant in this implementation.
 * \param [in] q2       Valid quadrant in this implementation.
 * \param [out] j       Non-NULL reference assigned comparison on output:
 *                      negative if \a q1 is less than \a q2,
 *                      positive if \a q1 is greater than \a q2,
 *                      zero otherwise (on equality).
 * \return              NULL on success, error object otherwise.
 */
sc3_error_t        *p4est3_quadrant_compare (const p4est3_quadrant_vtable_t *
                                             qvt, const void *q1,
                                             const void *q2, int *j);

/** Generate the root quadrant of a tree.
 * \param [in] qvt      Valid virtual quadrant table.
 * \param [out] r       The root quadrant is placed here.
 * \return              NULL on success, error object otherwise.
 */
sc3_error_t        *p4est3_quadrant_root (const p4est3_quadrant_vtable_t * qvt,
                                          void *r);

/** Deep copy one quadrant into another.
 * \param [in] qvt      Valid virtual quadrant table.
 * \param [in] q        Valid quadrant in this implementation.
 * \param [out] r       The copied quadrant is placed here.
 * \return              NULL on success, error object otherwise.
 */
sc3_error_t        *p4est3_quadrant_copy (const p4est3_quadrant_vtable_t * qvt,
                                          const void *q, void *r);

/** Generate the parent quadrant.
 * \param [in] qvt      Valid virtual quadrant table.
 * \param [in] q        Valid quadrant in this implementation.
 * \param [out] r       The parent quadrant is placed here.
 * \return              NULL on success, error object otherwise.
 */
sc3_error_t        *p4est3_quadrant_parent (const p4est3_quadrant_vtable_t * qvt,
                                            const void *q, void *r);

/** Generate the sibling quadrant.
 * \param [in] qvt      Valid virtual quadrant table.
 * \param [in] q        Valid quadrant in this implementation.
 * \param [in] i        The id of the sibling computed.
 * \param [in] r        The i-th sibling of the quadrant q will be placed here.
 * \return              NULL on success, error object otherwise.
 */
sc3_error_t        *p4est3_quadrant_sibling (const p4est3_quadrant_vtable_t * qvt,
                                             const void *q, int i, void *r);

/** Generate the face neighbor quadrant within the same tree.
 * \param [in] qvt      Valid virtual quadrant table.
 * \param [in] q        Valid quadrant in this implementation.
 * \param [in] i        The face across which to generate the neighbor.
 * \param [out] r       The neighbor quadrant is placed here in existing memory.
 * \return              NULL on success, error object otherwise.
 */
sc3_error_t        *p4est3_quadrant_face_neighbor (const p4est3_quadrant_vtable_t *
                                                   qvt, const void *q, int i,
                                                   void *r);

/** Compute the face neighbor of a quadrant across a tree boundary.
 * \param [in] qvt      Valid virtual quadrant table.
 * \param [in] q        Valid quadrant in this implementation.
 * \param [in] transform  This array holds 9 integers.
 *                        For 3D:
 *            [0]..[2]    The coordinate axis sequence of the origin face.
 *            [3]..[5]    The coordinate axis sequence of the target face.
 *            [6]..[8]    Edge reverse flag for axes t1, t2; face code for n.
 *                        For 2D:
 *            [0,2]       The coordinate axis sequence of the origin face.
 *            [3,5]       The coordinate axis sequence of the target face.
 *            [6,8]       Edge reverse flag for axis t; face code for axis n.
 *            [1,4,7]     0 (unused for compatibility with 3D).
 * \param [in] i        The face across which to generate the neighbor.
 * \param [out] r       The neighbor quadrant is placed here in existing memory.
 * \return              NULL on success, error object otherwise.
 */
sc3_error_t
  * p4est3_quadrant_tree_face_neighbor (const p4est3_quadrant_vtable_t * qvt,
                                        const void *q,
                                        sc3_array_t * transform, int i,
                                        void *r);

/** Generate the predecessor quadrant.
 * \param [in] qvt      Valid virtual quadrant table.
 * \param [in] q        Valid quadrant in this implementation,
 *                      must not be the first on its level in this tree.
 * \param [out] r       The predecessor by linear index is placed here.
 * \return              NULL on success, error object otherwise.
 */
sc3_error_t        *p4est3_quadrant_predecessor (const p4est3_quadrant_vtable_t *
                                                 qvt, const void *q, void *r);

/** Generate the successor quadrant.
 * \param [in] qvt      Valid virtual quadrant table.
 * \param [in] q        Valid quadrant in this implementation,
 *                      must not be the last on its level in this tree.
 * \param [out] r       The successor by linear index is placed here.
 * \return              NULL on success, error object otherwise.
 */
sc3_error_t        *p4est3_quadrant_successor (const p4est3_quadrant_vtable_t * qvt,
                                               const void *q, void *r);

/** Generate a child quadrant specified by child id.
 * \param [in] qvt      Valid virtual quadrant table.
 * \param [in] q        Valid quadrant in this implementation,
 *                      must not be of maximum level.
 * \param [in] i        Child id valid with respect to \a q.
 * \param [out] r       The i-th child of \a q is placed here.
 * \return              NULL on success, error object otherwise.
 */
sc3_error_t        *p4est3_quadrant_child (const p4est3_quadrant_vtable_t * qvt,
                                           const void *q, int i, void *r);

/** Generate an ancestor quadrant specified by ancestor level.
 * \param [in] qvt      Valid virtual quadrant table.
 * \param [in] q        Valid quadrant in this implementation.
 * \param [in] l        Non-negative ancestor level less or equal than \a q's.
 * \param [out] r       The ancestor of \a q on level \a l is placed here.
 * \return              NULL on success, error object otherwise.
 */
sc3_error_t        *p4est3_quadrant_ancestor (const p4est3_quadrant_vtable_t * qvt,
                                              const void *q, int l, void *r);

/** Generate the first descendant quadrant on a specified level.
 * \param [in] qvt      Valid virtual quadrant table.
 * \param [in] q        Valid quadrant in this implementation.
 * \param [in] l        Non-negative level greater or equal than \a q's.
 * \param [out] r       First descendant of \a q on level \a l is placed here.
 * \return              NULL on success, error object otherwise.
 */
sc3_error_t        *p4est3_quadrant_first_descendant (const p4est3_quadrant_vtable_t
                                                      * qvt, const void *q,
                                                      int l, void *r);

/** Generate the last descendant quadrant on a specified level.
 * \param [in] qvt      Valid virtual quadrant table.
 * \param [in] q        Valid quadrant in this implementation.
 * \param [in] l        Non-negative level greater or equal than \a q's.
 * \param [out] r       Last descendant of \a q on level \a l is placed here.
 * \return              NULL on success, error object otherwise.
 */
sc3_error_t        *p4est3_quadrant_last_descendant (const p4est3_quadrant_vtable_t
                                                     * qvt, const void *q,
                                                     int l, void *r);

/** Generate a quadrant by linear index on a given level.
 * This function works up to level l=31 in 2D and l=21 in 3D.
 * \param [in] qvt      Valid virtual quadrant table.
 * \param [in] l        Valid level available with this prototype.
 * \param [in] g        Valid linear index for level \a l.
 * \param [out] r       Quadrant of given linear index is placed here.
 * \return              NULL on success, error object otherwise.
 */
sc3_error_t        *p4est3_quadrant_morton (const p4est3_quadrant_vtable_t * qvt,
                                            int l, p4est3_gloidx id, void *r);

/** Generate a common nearest ancestor of two quadrants. 
 * \param [in] qvt      Valid virtual quadrant table.
 * \param [in] q1       Valid quadrant in this implementation.
 * \param [in] q2       Valid quadrant in this implementation.
 * \param [out] r       Quadrant, common nearest ancestor of q1 and q2.
 * \return              NULL on success, error object otherwise.
*/
sc3_error_t        *p4est3_nearest_common_ancestor (const p4est3_quadrant_vtable_t
                                                    * qvt, const void *q1,
                                                    const void *q2, void *r);

/** Generate a linear index by quadrant on a given level.
 * This function works up to level l=31 in 2D and l=21 in 3D.
 * This function is reverse for p4est3_quadrant_morton.
 * \param [in] qvt      Valid virtual quadrant table.
 * \param [in] q        Valid quadrant in this implementation.
 * \param [in] l        Valid level available with this prototype.
 * \param [out] id      Linear index of given quadrant is placed here.
 * \return              NULL on success, error object otherwise.
*/
sc3_error_t        *p4est3_quadrant_linear_id (const p4est3_quadrant_vtable_t * qvt,
                                               const void *q, int l,
                                               p4est3_gloidx * id);

/** Query if quadrant q1 is an ancestor of quadrant q2.
 * \param [in] qvt      Valid virtual quadrant table.
 * \param [in] q1       Valid quadrant in this implementation.
 * \param [in] q2       Valid quadrant in this implementation.
 * \param [out] j       True if q1 is ancestor of q2, false otherwise.
 * \return              NULL on success, error object otherwise.
*/
sc3_error_t        *p4est3_quadrant_is_ancestor (const p4est3_quadrant_vtable_t
                                                 * qvt, const void *q1,
                                                 const void *q2, int *j);

/** Query if quadrant q1 is a parent of quadrant q2.
 * \param [in] qvt      Valid virtual quadrant table.
 * \param [in] q1       Valid quadrant in this implementation.
 * \param [in] q2       Valid quadrant in this implementation.
 * \param [out] j       True if q1 is parent of q2, false otherwise.
 * \return              NULL on success, error object otherwise.
*/
sc3_error_t        *p4est3_quadrant_is_parent (const p4est3_quadrant_vtable_t
                                               * qvt, const void *q1,
                                               const void *q2, int *j);

/************************ convenience functions **************************/

/** Create an array of given length intended to store quadrants.
 * \param [in,out] alloc    Allocator must be setup.
 * \param [in] qvt          Valid virtual quadrant table.
 * \param [in] n            Number of quadrants to allocate.
 * \param [out] arr         This array is of correct size and setup.
 * \return                  NULL on success, error object otherwise.
 */
sc3_error_t        *p4est3_quadrant_array_new (sc3_allocator_t * alloc,
                                               const p4est3_quadrant_vtable_t * qvt,
                                               p4est3_locidx n,
                                               sc3_array_t ** arr);

#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif /* !P4EST3_QUADRANT_VTABLE_H */
