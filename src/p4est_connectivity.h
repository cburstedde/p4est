/*
  This file is part of p4est.
  p4est is a C library to manage a collection (a forest) of multiple
  connected adaptive quadtrees or octrees in parallel.

  Copyright (C) 2010 The University of Texas System
  Additional copyright (C) 2011 individual authors
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

/** \file p4est_connectivity.h
 *
 * The coarse topological description of the forest.
 *
 * \ingroup p4est
 */

#ifndef P4EST_CONNECTIVITY_H
#define P4EST_CONNECTIVITY_H

#ifdef P4_TO_P8
#error "Including a p4est header with P4_TO_P8 defined"
#endif

#include <sc_io.h>
#include <p4est_base.h>

SC_EXTERN_C_BEGIN;

/** The spatial dimension */
#define P4EST_DIM 2
/** The number of faces of a quadrant */
#define P4EST_FACES (2 * P4EST_DIM)
/** The number of children of a quadrant
 *
 * also the number of corners */
#define P4EST_CHILDREN 4
/** The number of children/corners touching one face */
#define P4EST_HALF (P4EST_CHILDREN / 2)
/** The size of insulation layer */
#define P4EST_INSUL 9

/* size of face transformation encoding */
#define P4EST_FTRANSFORM 9

/** p4est identification string */
#define P4EST_STRING "p4est"

/* Increase this number whenever the on-disk format for
 * p4est_connectivity, p4est, or any other 2D data structure changes.
 * The format for reading and writing must be the same.
 */
#define P4EST_ONDISK_FORMAT 0x2000009

/** Characterize a type of adjacency.
 *
 * Several functions involve relationships between neighboring trees and/or
 * quadrants, and their behavior depends on how one defines adjacency:
 * 1) entities are adjacent if they share a face, or
 * 2) entities are adjacent if they share a face or corner.
 * p4est_connect_type_t is used to choose the desired behavior.
 * This enum must fit into an int8_t.
 */
typedef enum
{
  /* make sure to have different values 2D and 3D */
  P4EST_CONNECT_FACE = 21,
  P4EST_CONNECT_CORNER = 22,
  P4EST_CONNECT_FULL = P4EST_CONNECT_CORNER
}
p4est_connect_type_t;

#ifdef P4EST_BACKWARD_DEALII
typedef p4est_connect_type_t p4est_balance_type_t;
#endif

/** Typedef for serialization method. */
typedef enum
{
  P4EST_CONN_ENCODE_NONE = SC_IO_ENCODE_NONE,
  P4EST_CONN_ENCODE_LAST        /**< Invalid entry to close the list. */
}
p4est_connectivity_encode_t;

/** Convert the p4est_connect_type_t into a number.
 * \param [in] btype    The balance type to convert.
 * \return              Returns 1 or 2.
 */
int                 p4est_connect_type_int (p4est_connect_type_t btype);

/** Convert the p4est_connect_type_t into a const string.
 * \param [in] btype    The balance type to convert.
 * \return              Returns a pointer to a constant string.
 */
const char         *p4est_connect_type_string (p4est_connect_type_t btype);

/** This structure holds the 2D inter-tree connectivity information.
 * Identification of arbitrary faces and corners is possible.
 *
 * The arrays tree_to_* are stored in z ordering.
 * For corners the order wrt. yx is 00 01 10 11.
 * For faces the order is -x +x -y +y.
 * They are allocated [0][0]..[0][3]..[num_trees-1][0]..[num_trees-1][3].
 *
 * The values for tree_to_face are 0..7
 * where ttf % 4 gives the face number and ttf / 4 the face orientation code.
 * The orientation is 0 for edges that are aligned in z-order,
 * and 1 for edges that are running opposite in z-order.
 *
 * It is valid to specify num_vertices as 0.
 * In this case vertices and tree_to_vertex are set to NULL.
 * Otherwise the vertex coordinates are stored in the array vertices as
 * [0][0]..[0][2]..[num_vertices-1][0]..[num_vertices-1][2].
 *
 * The corners are only stored when they connect trees.
 * In this case tree_to_corner indexes into \a ctt_offset.
 * Otherwise the tree_to_corner entry must be -1 and this corner is ignored.
 * If num_corners == 0, tree_to_corner and corner_to_* arrays are set to NULL.
 *
 * The arrays corner_to_* store a variable number of entries per corner.
 * For corner c these are at position [ctt_offset[c]]..[ctt_offset[c+1]-1].
 * Their number for corner c is ctt_offset[c+1] - ctt_offset[c].
 * The entries encode all trees adjacent to corner c.
 * The size of the corner_to_* arrays is num_ctt = ctt_offset[num_corners].
 *
 * The *_to_attr arrays may have arbitrary contents defined by the user.
 */
typedef struct p4est_connectivity
{
  p4est_topidx_t      num_vertices; /**< the number of vertices that define
                                         the \a embedding of the forest (not
                                         the topology) */
  p4est_topidx_t      num_trees;    /**< the number of trees */
  p4est_topidx_t      num_corners;  /**< the number of corners that help
                                         define topology */
  double             *vertices;     /**< an array of size
                                         (3 * \a num_vertices) */
  p4est_topidx_t     *tree_to_vertex; /**< embed each tree into \f$R^3\f$ for
                                           e.g. visualization (see
                                           p4est_vtk.h) */

  size_t              tree_attr_bytes;  /**< bytes per tree in tree_to_attr */
  char               *tree_to_attr;     /**< not touched by p4est */

  p4est_topidx_t     *tree_to_tree; /**< (4 * \a num_trees) neighbors across
                                         faces */
  int8_t             *tree_to_face; /**< (4 * \a num_trees) face to
                                         face+orientation (see description) */

  p4est_topidx_t     *tree_to_corner; /**< (4 * \a num_trees) or NULL (see
                                           description) */
  p4est_topidx_t     *ctt_offset; /**< corner to offset in \a corner_to_tree
                                       and \a corner_to_corner */
  p4est_topidx_t     *corner_to_tree; /**< list of trees that meet at a corner */
  int8_t             *corner_to_corner; /**< list of tree-corners that meet at
                                             a corner */
}
p4est_connectivity_t;

/** Calculate memory usage of a connectivity structure.
 * \param [in] conn   Connectivity structure.
 * \return            Memory used in bytes.
 */
size_t              p4est_connectivity_memory_used (p4est_connectivity_t *
                                                    conn);

typedef struct
{
  p4est_topidx_t      ntree;
  int8_t              ncorner;
}
p4est_corner_transform_t;

typedef struct
{
  p4est_topidx_t      icorner;
  sc_array_t          corner_transforms;
}
p4est_corner_info_t;

/** Store the corner numbers 0..4 for each tree face. */
extern const int    p4est_face_corners[4][2];

/** Store the face numbers in the face neighbor's system. */
extern const int    p4est_face_dual[4];

/** Store the face numbers 0..3 for each tree corner. */
extern const int    p4est_corner_faces[4][2];

/** Store the face corner numbers for the faces touching a tree corner. */
extern const int    p4est_corner_face_corners[4][4];

/** Store the faces for each child and corner, can be -1. */
extern const int    p4est_child_corner_faces[4][4];

/** Transform a face corner across one of the adjacent faces into a neighbor tree.
 * This version expects the neighbor face and orientation separately.
 * \param [in] fc   A face corner number in 0..1.
 * \param [in] f    A face that the face corner number \a fc is relative to.
 * \param [in] nf   A neighbor face that is on the other side of \f.
 * \param [in] o    The orientation between tree boundary faces \a f and \nf.
 * \return          The face corner number relative to the neighbor's face.
 */
int                 p4est_connectivity_face_neighbor_face_corner
  (int fc, int f, int nf, int o);

/** Transform a corner across one of the adjacent faces into a neighbor tree.
 * This version expects the neighbor face and orientation separately.
 * \param [in] c    A corner number in 0..3.
 * \param [in] f    A face number that touches the corner \a c.
 * \param [in] nf   A neighbor face that is on the other side of \f.
 * \param [in] o    The orientation between tree boundary faces \a f and \nf.
 * \return          The number of the corner seen from the neighbor tree.
 */
int                 p4est_connectivity_face_neighbor_corner
  (int c, int f, int nf, int o);

/** Allocate a connectivity structure.
 * The attribute fields are initialized to NULL.
 * \param [in] num_vertices   Number of total vertices (i.e. geometric points).
 * \param [in] num_trees      Number of trees in the forest.
 * \param [in] num_corners    Number of tree-connecting corners.
 * \param [in] num_ctt        Number of total trees in corner_to_tree array.
 * \return                    A connectivity structure with allocated arrays.
 */
p4est_connectivity_t *p4est_connectivity_new (p4est_topidx_t num_vertices,
                                              p4est_topidx_t num_trees,
                                              p4est_topidx_t num_corners,
                                              p4est_topidx_t num_ctt);

/** Allocate a connectivity structure and populate from constants.
 * The attribute fields are initialized to NULL.
 * \param [in] num_vertices   Number of total vertices (i.e. geometric points).
 * \param [in] num_trees      Number of trees in the forest.
 * \param [in] num_corners    Number of tree-connecting corners.
 * \param [in] coff           Corner-to-tree offsets (num_corners + 1 values).
 *                            This must always be non-NULL; in trivial cases
 *                            it is just a pointer to a p4est_topix value of 0.
 * \return                    The connectivity is checked for validity.
 */
p4est_connectivity_t *p4est_connectivity_new_copy (p4est_topidx_t
                                                   num_vertices,
                                                   p4est_topidx_t num_trees,
                                                   p4est_topidx_t num_corners,
                                                   const double *vertices,
                                                   const p4est_topidx_t * ttv,
                                                   const p4est_topidx_t * ttt,
                                                   const int8_t * ttf,
                                                   const p4est_topidx_t * ttc,
                                                   const p4est_topidx_t *
                                                   coff,
                                                   const p4est_topidx_t * ctt,
                                                   const int8_t * ctc);

/** Broadcast a connectivity structure that exists only on one process to all.
 *  On the other processors, it will be allocated using p4est_connectivity_new.
 *  \param [in] conn_in For the root process the connectivity to be broadcast,
 *                      for the other processes it must be NULL.
 *  \param [in] root    The rank of the process that provides the connectivity.
 *  \param [in] comm    The MPI communicator.
 *  \return             For the root process this is a pointer to \a conn_in.
 *                      Else, a pointer to a newly allocated connectivity
 *                      structure with the same values as \a conn_in on the
 *                      root process.
 */
p4est_connectivity_t *p4est_connectivity_bcast (p4est_connectivity_t *
                                                conn_in, int root,
                                                sc_MPI_Comm comm);

/** Destroy a connectivity structure.  Also destroy all attributes.
 */
void                p4est_connectivity_destroy (p4est_connectivity_t *
                                                connectivity);

/** Allocate or free the attribute fields in a connectivity.
 * \param [in,out] conn         The conn->*_to_attr fields must either be NULL
 *                              or previously be allocated by this function.
 * \param [in] bytes_per_tree   If 0, tree_to_attr is freed (being NULL is ok).
 *                              If positive, requested space is allocated.
 */
void                p4est_connectivity_set_attr (p4est_connectivity_t * conn,
                                                 size_t bytes_per_tree);

/** Examine a connectivity structure.
 * \return          Returns true if structure is valid, false otherwise.
 */
int                 p4est_connectivity_is_valid (p4est_connectivity_t *
                                                 connectivity);

/** Check two connectivity structures for equality.
 * \return          Returns true if structures are equal, false otherwise.
 */
int                 p4est_connectivity_is_equal (p4est_connectivity_t * conn1,
                                                 p4est_connectivity_t *
                                                 conn2);

/** Write connectivity to a sink object.
 * \param [in] conn     The connectivity to be written.
 * \param [in,out] sink The connectivity is written into this sink.
 * \return              0 on success, nonzero on error.
 */
int                 p4est_connectivity_sink (p4est_connectivity_t * conn,
                                             sc_io_sink_t * sink);

/** Allocate memory and store the connectivity information there.
 * \param [in] conn     The connectivity structure to be exported to memory.
 * \param [in] code     Encoding and compression method for serialization.
 * \return              Newly created array that contains the information.
 */
sc_array_t         *p4est_connectivity_deflate (p4est_connectivity_t * conn,
                                                p4est_connectivity_encode_t
                                                code);

/** Save a connectivity structure to disk.
 * \param [in] filename         Name of the file to write.
 * \param [in] connectivity     Valid connectivity structure.
 * \return                      Returns 0 on success, nonzero on file error.
 */
int                 p4est_connectivity_save (const char *filename,
                                             p4est_connectivity_t *
                                             connectivity);

/** Read connectivity from a source object.
 * \param [in,out] source       The connectivity is read from this source.
 * \return              The newly created connectivity, or NULL on error.
 */
p4est_connectivity_t *p4est_connectivity_source (sc_io_source_t * source);

/** Create new connectivity from a memory buffer.
 * \param [in] buffer   The connectivity is created from this memory buffer.
 * \return              The newly created connectivity, or NULL on error.
 */
p4est_connectivity_t *p4est_connectivity_inflate (sc_array_t * buffer);

/** Load a connectivity structure from disk.
 * \param [in] filename         Name of the file to read.
 * \param [in,out] bytes        Size in bytes of connectivity on disk or NULL.
 * \return              Returns valid connectivity, or NULL on file error.
 */
p4est_connectivity_t *p4est_connectivity_load (const char *filename,
                                               size_t * bytes);

/** Create a connectivity structure for the unit square.
 */
p4est_connectivity_t *p4est_connectivity_new_unitsquare (void);

/** Create a connectivity structure for an all-periodic unit square.
 */
p4est_connectivity_t *p4est_connectivity_new_periodic (void);

/** Create a connectivity structure for a periodic unit square.
 * The left and right faces are identified, and bottom and top opposite.
 */
p4est_connectivity_t *p4est_connectivity_new_rotwrap (void);

/** Create a connectivity structure for two trees being rotated
 * w.r.t. each other in a user-defined way
 * \param[in] l_face      index of left face
 * \param[in] r_face      index of right face
 * \param[in] orientation orientation of trees w.r.t. each other
 */
p4est_connectivity_t *p4est_connectivity_new_twotrees (int l_face,
                                                       int r_face,
                                                       int orientation);

/** Create a connectivity structure for a three-tree mesh around a corner.
 */
p4est_connectivity_t *p4est_connectivity_new_corner (void);

/** Create a connectivity structure for two trees on top of each other.
 */
p4est_connectivity_t *p4est_connectivity_new_pillow (void);

/** Create a connectivity structure for a five-tree moebius band.
 */
p4est_connectivity_t *p4est_connectivity_new_moebius (void);

/** Create a connectivity structure for a six-tree star.
 */
p4est_connectivity_t *p4est_connectivity_new_star (void);

/** Create a connectivity structure for the six sides of a unit cube.
 * The ordering of the trees is as follows: 0 1
 *                                            2 3 <-- 3: axis-aligned top side
 *                                              4 5.
 * This choice has been made for maximum symmetry (see tree_to_* in .c file).
 */
p4est_connectivity_t *p4est_connectivity_new_cubed (void);

/** Create a connectivity structure for a five-tree flat spherical disk.
 * This disk can just as well be used as a square to test non-Cartesian maps.
 * Without any mapping this connectivity covers the square [-3, 3]**2.
 * \note The API of this function has changed to accept two arguments.
 *       You can query the #define P4EST_CONN_DISK_PERIODIC to check
 *       whether the new version with the argument is in effect.
 *
 * The ordering of the trees is as follows:   4
 *                                          1 2 3
 *                                            0.
 *
 * The outside x faces may be identified topologically.
 * The outside y faces may be identified topologically.
 * Both identifications may be specified simultaneously.
 * The general shape and periodicity are the same as those obtained with
 * \ref p4est_connectivity_new_brick (1, 1, periodic_a, periodic_b).
 *
 * \param [in] periodic_a       Bool to make disk periodic in x direction.
 * \param [in] periodic_b       Bool to make disk periodic in y direction.
 */
p4est_connectivity_t *p4est_connectivity_new_disk (int periodic_a,
                                                   int periodic_b);

/** A rectangular m by n array of trees with configurable periodicity.
 * The brick is periodic in x and y if periodic_a and periodic_b are true,
 * respectively.
 */
p4est_connectivity_t *p4est_connectivity_new_brick (int mi, int ni,
                                                    int periodic_a,
                                                    int periodic_b);

/** Create connectivity structure from predefined catalogue.
 * \param [in]  name            Invokes connectivity_new_* function.
 *              brick23         brick (2, 3, 0, 0)
 *              corner          corner
 *              cubed           cubed
 *              disk            disk
 *              moebius         moebius
 *              periodic        periodic
 *              pillow          pillow
 *              rotwrap         rotwrap
 *              star            star
 *              unit            unitsquare
 * \return      An initialized connectivity if name is defined, NULL else.
 */
p4est_connectivity_t *p4est_connectivity_new_byname (const char *name);

/** Uniformly refine a connectivity.
 * This is useful if you would like to uniformly refine by something other
 * than a power of 2.
 *
 * \param [in] conn         a valid connectivity
 * \param [in] num_per_edge the number of new trees in each direction
 *
 * \return a refined connectivity.
 */
p4est_connectivity_t *p4est_connectivity_refine (p4est_connectivity_t * conn,
                                                 int num_per_edge);

/** Fill an array with the axis combination of a face neighbor transform.
 * \param [in]  iface       The number of the originating face.
 * \param [in]  nface       Encoded as nface = r * 4 + nf, where nf = 0..3 is
 *                          the neigbbor's connecting face number and r = 0..1
 *                          is the relative orientation to the neighbor's face.
 *                          This encoding matches p4est_connectivity_t.
 * \param [out] ftransform  This array holds 9 integers.
 *              [0,2]       The coordinate axis sequence of the origin face,
 *                          the first referring to the tangential and the second
 *                          to the normal.  A permutation of (0, 1).
 *              [3,5]       The coordinate axis sequence of the target face.
 *              [6,8]       Edge reversal flag for tangential axis (boolean);
 *                          face code in [0, 3] for the normal coordinate q:
 *                          0: q' = -q
 *                          1: q' = q + 1
 *                          2: q' = q - 1
 *                          3: q' = 2 - q
 *              [1,4,7]     0 (unused for compatibility with 3D).
 */
void                p4est_expand_face_transform (int iface, int nface,
                                                 int ftransform[]);

/** Fill an array with the axis combinations of a tree neighbor transform.
 * \param [in]  itree       The number of the originating tree.
 * \param [in]  iface       The number of the originating tree's face.
 * \param [out] ftransform  This array holds 9 integers.
 *              [0,2]       The coordinate axis sequence of the origin face.
 *              [3,5]       The coordinate axis sequence of the target face.
 *              [6,8]       Edge reverse flag for axis t; face code for axis n.
 *              [1,4,7]     0 (unused for compatibility with 3D).
 * \return                  The face neighbor tree if it exists, -1 otherwise.
 */
p4est_topidx_t      p4est_find_face_transform (p4est_connectivity_t *
                                               connectivity,
                                               p4est_topidx_t itree,
                                               int iface, int ftransform[]);

/** Fills an array with information about corner neighbors.
 * \param [in] itree    The number of the originating tree.
 * \param [in] icorner  The number of the originating corner.
 * \param [in,out] ci   A p4est_corner_info_t structure with initialized array.
 */
void                p4est_find_corner_transform (p4est_connectivity_t *
                                                 connectivity,
                                                 p4est_topidx_t itree,
                                                 int icorner,
                                                 p4est_corner_info_t * ci);

/** Internally connect a connectivity based on tree_to_vertex information.
 * Periodicity that is not inherent in the list of vertices will be lost.
 * \param [in,out] conn     The connectivity needs to have proper vertices
 *                          and tree_to_vertex fields.  The tree_to_tree
 *                          and tree_to_face fields must be allocated
 *                          and satisfy p4est_connectivity_is_valid (conn)
 *                          but will be overwritten.  The corner
 *                          fields will be freed and allocated anew.
 */
void                p4est_connectivity_complete (p4est_connectivity_t * conn);

/** Removes corner information of a connectivity
 *  such that enough information is left to run p4est_connectivity_complete successfully.
 *  The reduced connectivity still passes p4est_connectivity_is_valid.
 * \param [in,out] conn     The connectivity to be reduced.
 */
void                p4est_connectivity_reduce (p4est_connectivity_t * conn);

/** p4est_connectivity_permute
 * Given a permutation \a perm of the trees in a connectivity \a conn,
 * permute the trees of \a conn in place and update \a conn to match.
 * \param [in,out] conn                The connectivity whose trees are
 *                                     permuted.
 * \param [in] perm                    A permutation array, whose elements are
 *                                     size_t's.
 * \param [in] is_current_to_new       if true, the jth entry of perm is the
 *                                     new index for the entry whose current
 *                                     index is j, otherwise the jth entry of
 *                                     perm is the current index of the tree
 *                                     whose index will be j after the
 *                                     permutation.
 */
void                p4est_connectivity_permute (p4est_connectivity_t * conn,
                                                sc_array_t * perm,
                                                int is_current_to_new);
#ifdef P4EST_WITH_METIS

/** p4est_connectivity_reorder
 * This function takes a connectivity \a conn and a parameter \a k,
 * which will typically be the number of processes, and reorders the trees
 * such that if every processes is assigned (num_trees / k) trees, the
 * communication volume will be minimized.  This is intended for use with
 * connectivities that contain a large number of trees.  This should be done
 * BEFORE a p4est is created using the connectivity.  This is done in place:
 * any data structures that use indices to refer to trees before this
 * procedure will be invalid.  Note that this routine calls metis and not
 * parmetis because the connectivity is copied on every process.
 * A communicator is required because I'm not positive that metis is
 * deterministic. \a ctype determines when an edge exist between two trees in
 * the dual graph used by metis in the reordering.
 * \param [in]     comm       MPI communicator.
 * \param [in]     k          if k > 0, the number of pieces metis will use to
 *                            guide the reordering; if k = 0, the number of
 *                            pieces will be determined from the MPI
 *                            communicator.
 * \param [in,out] conn       connectivity that will be reordered.
 * \param [in]     ctype      determines when an edge exists in the dual graph
 *                            of the connectivity structure.
 */
void                p4est_connectivity_reorder (sc_MPI_Comm comm, int k,
                                                p4est_connectivity_t * conn,
                                                p4est_connect_type_t ctype);

#endif /* P4EST_WITH_METIS */

/** p4est_connectivity_join_faces
 * This function takes an existing valid connectivity \a conn and modifies it
 * by joining two tree faces that are currently boundary faces.
 * \param [in,out] conn        connectivity that will be altered.
 * \param [in]     tree_left   tree that will be on the left side of the joined
 *                             faces.
 * \param [in]     tree_right  tree that will be on the right side of the
 *                             joined faces.
 * \param [in]     face_left   face of \a tree_left that will be joined.
 * \param [in]     face_right  face of \a tree_right that will be joined.
 * \param [in]     orientation the orientation of \a face_left and
 *                             \a face_right once joined (see the description
 *                             of p4est_connectivity_t to understand
 *                             orientation).
 */
void                p4est_connectivity_join_faces (p4est_connectivity_t *
                                                   conn,
                                                   p4est_topidx_t tree_left,
                                                   p4est_topidx_t tree_right,
                                                   int face_left,
                                                   int face_right,
                                                   int orientation);

/** p4est_connectivity_is_equivalent
 * This function compares two connectivities for equivalence: it returns
 * \a true if they are the same connectivity, or if they have the same
 * topology.  The definition of topological sameness is strict: there is no
 * attempt made to determine whether permutation and/or rotation of the trees
 * makes the connectivities equivalent.
 *
 * \param[in]      conn1    a valid connectivity
 * \param[out]     conn2    a valid connectivity
 */
int                 p4est_connectivity_is_equivalent (p4est_connectivity_t *
                                                      conn1,
                                                      p4est_connectivity_t *
                                                      conn2);

/** Return a pointer to a p4est_corner_transform_t array element. */
/*@unused@*/
static inline p4est_corner_transform_t *
p4est_corner_array_index (sc_array_t * array, size_t it)
{
  P4EST_ASSERT (array->elem_size == sizeof (p4est_corner_transform_t));
  P4EST_ASSERT (it < array->elem_count);

  return
    (p4est_corner_transform_t *) (array->array +
                                  sizeof (p4est_corner_transform_t) * it);
}

/** Read an ABAQUS input file from a file stream.
 *
 * This utility function reads a basic ABAQUS file supporting element type with
 * the prefix C2D4, CPS4, and S4 in 2D and of type C3D8 reading them as
 * bilinear quadrilateral and trilinear hexahedral trees respectively.
 *
 * A basic 2D mesh is given below.  The \c *Node section gives the vertex
 * number and x, y, and z components for each vertex.  The \c *Element section
 * gives the 4 vertices in 2D (8 vertices in 3D) of each element in counter
 * clockwise order. So in 2D the nodes are given as:
 *
 *   4                     3
 *   +-------------------+
 *   |                   |
 *   |                   |
 *   |                   |
 *   |                   |
 *   |                   |
 *   |                   |
 *   +-------------------+
 *   1                   2
 *
 * and in 3D they are given as:
 *
 * 8                     7
 *  +---------------------+
 *  |\                    |\
 *  | \                   | \
 *  |  \                  |  \
 *  |   \                 |   \
 *  |   5+---------------------+6
 *  |    |                |    |
 *  +----|----------------+    |
 *  4\   |               3 \   |
 *    \  |                  \  |
 *     \ |                   \ |
 *      \|                    \|
 *       +---------------------+
 *       1                     2
 *
 * \code
 * *Heading
 *  box.inp
 * *Node
 * 1,  -5, -5, 0
 * 2,   5, -5, 0
 * 3,   5,  5, 0
 * 4,  -5,  5, 0
 * 5,   0, -5, 0
 * 6,   5,  0, 0
 * 7,   0,  5, 0
 * 8,  -5,  0, 0
 * 9,   1, -1, 0
 * 10,  0,  0, 0
 * 11, -2,  1, 0
 * *Element, type=CPS4, ELSET=Surface1
 * 1,  1, 10, 11, 8
 * 2,  3, 10, 9,  6
 * 3,  9, 10, 1,  5
 * 4,  7,  4, 8, 11
 * 5, 11, 10, 3,  7
 * 6,  2,  6, 9,  5
 * \endcode
 *
 * This code can be called two ways.  The first, when \c vertex==NULL and \c
 * tree_to_vertex==NULL, is used to count the number of trees and vertices in
 * the connectivity to be generated by the \c .inp mesh in the \a stream.  The
 * second, when \c vertices!=NULL and \c tree_to_vertex!=NULL, fill \c vertices
 * and \c tree_to_vertex.  In this case \c num_vertices and \c num_trees need
 * to be set to the maximum number of entries allocated in \c vertices and \c
 * tree_to_vertex.
 *
 * \param[in,out]  stream         file stream to read the connectivity from
 * \param[in,out]  num_vertices   the number of vertices in the connectivity
 * \param[in,out]  num_trees      the number of trees in the connectivity
 * \param[out]     vertices       the list of \c vertices of the connectivity
 * \param[out]     tree_to_vertex the \c tree_to_vertex map of the connectivity
 *
 * \returns 0 if successful and nonzero if not
 */
int                 p4est_connectivity_read_inp_stream (FILE * stream,
                                                        p4est_topidx_t *
                                                        num_vertices,
                                                        p4est_topidx_t *
                                                        num_trees,
                                                        double *vertices,
                                                        p4est_topidx_t *
                                                        tree_to_vertex);

/** Create a p4est connectivity from an ABAQUS input file.
 *
 * This utility function reads a basic ABAQUS file supporting element type with
 * the prefix C2D4, CPS4, and S4 in 2D and of type C3D8 reading them as
 * bilinear quadrilateral and trilinear hexahedral trees respectively.
 *
 * A basic 2D mesh is given below.  The \c *Node section gives the vertex
 * number and x, y, and z components for each vertex.  The \c *Element section
 * gives the 4 vertices in 2D (8 vertices in 3D) of each element in counter
 * clockwise order. So in 2D the nodes are given as:
 *
 *   4                     3
 *   +-------------------+
 *   |                   |
 *   |                   |
 *   |                   |
 *   |                   |
 *   |                   |
 *   |                   |
 *   +-------------------+
 *   1                   2
 *
 * and in 3D they are given as:
 *
 * 8                     7
 *  +---------------------+
 *  |\                    |\
 *  | \                   | \
 *  |  \                  |  \
 *  |   \                 |   \
 *  |   5+---------------------+6
 *  |    |                |    |
 *  +----|----------------+    |
 *  4\   |               3 \   |
 *    \  |                  \  |
 *     \ |                   \ |
 *      \|                    \|
 *       +---------------------+
 *       1                     2
 *
 * \code
 * *Heading
 *  box.inp
 * *Node
 * 1,  -5, -5, 0
 * 2,   5, -5, 0
 * 3,   5,  5, 0
 * 4,  -5,  5, 0
 * 5,   0, -5, 0
 * 6,   5,  0, 0
 * 7,   0,  5, 0
 * 8,  -5,  0, 0
 * 9,   1, -1, 0
 * 10,  0,  0, 0
 * 11, -2,  1, 0
 * *Element, type=CPS4, ELSET=Surface1
 * 1,  1, 10, 11, 8
 * 2,  3, 10, 9,  6
 * 3,  9, 10, 1,  5
 * 4,  7,  4, 8, 11
 * 5, 11, 10, 3,  7
 * 6,  2,  6, 9,  5
 * \endcode
 *
 * This function reads a mesh from \a filename and returns an associated p4est
 * connectivity.
 *
 * \param[in]  filename         file to read the connectivity from
 *
 * \return  an allocated connectivity associated with the mesh in \a filename
 *          or NULL if an error occurred.
 */
p4est_connectivity_t *p4est_connectivity_read_inp (const char *filename);

SC_EXTERN_C_END;

#endif /* !P4EST_CONNECTIVITY_H */
