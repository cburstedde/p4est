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

#ifndef P8EST_ITERATE_H
#define P8EST_ITERATE_H

#include <p8est.h>
#include <p8est_ghost.h>

SC_EXTERN_C_BEGIN;

/** The information that is available to the user-defined p8est_iter_volume_t
 * callback function.
 *
 * \a treeid gives the index in \a p4est->trees of the tree to which
 *    \a quad belongs.
 * \a quadid gives the index of \a quad within \a tree's quadrants array.
 */
typedef struct p8est_iter_volume_info
{
  p8est_t            *p4est;
  p8est_ghost_t      *ghost_layer;
  p8est_quadrant_t   *quad;
  p4est_locidx_t      quadid;
  p4est_topidx_t      treeid;
}
p8est_iter_volume_info_t;

/** The prototype for a function that p8est_iterate will execute at every
 * quadrant local to the current process.
 */
typedef void        (*p8est_iter_volume_t) (p8est_iter_volume_info_t * info,
                                            void *user_data);

/* Information about one side of a face in the forest.  If a \a quad is local
 * (\a is_ghost is false), then its \a quadid indexes the tree's quadrant array;
 * otherwise, it indexes the ghosts array. If the face is hanging, then the
 * quadrants are listed in z-order. If a quadrant should be present, but it is
 * not included in the ghost layer, then quad = NULL, is_ghost is true, and
 * quadid = -1.
 */
typedef struct p8est_iter_face_side
{
  p4est_topidx_t      treeid;
  int8_t              face;
  int8_t              is_hanging;
  union p8est_iter_face_side_data
  {
    struct
    {
      int8_t              is_ghost;
      p8est_quadrant_t   *quad;
      p4est_locidx_t      quadid;
    }
    full;
    struct
    {
      int8_t              is_ghost[4];
      p8est_quadrant_t   *quad[4];
      p4est_locidx_t      quadid[4];
    }
    hanging;
  }
  is;
}
p8est_iter_face_side_t;

/** The information about both sides of a face in the forest.  The orientation
 * is 0 if the face is within one tree; otherwise, it is the same as the
 * orientation value between the two trees given in the connectivity.  If the
 * face is on the outside of the forest, then there is only one side.
 * If tree_boundary is false, the face is on the interior of a tree.
 * When tree_boundary false, sides[0] contains the lowest z-order quadrant that
 * touches the face.
 * When tree_boundary is true, its value is P8EST_CONNECT_FACE.
 */
typedef struct p8est_iter_face_info
{
  p8est_t            *p4est;
  p8est_ghost_t      *ghost_layer;
  int8_t              orientation;
  int8_t              tree_boundary;
  sc_array_t          sides;    /* p8est_iter_face_side_t */
}
p8est_iter_face_info_t;

/** The prototype for a function that p8est_iterate will execute wherever two
 * quadrants share a face: the face can be a 2:1 hanging face, it does not have
 * to be conformal.
 *
 * Note: the forest must be face balanced for p8est_iterate to execute a
 * callback function on faces.
 */
typedef void        (*p8est_iter_face_t) (p8est_iter_face_info_t * info,
                                          void *user_data);

/* Information about one side of an edge in the forest.  If a \a quad is local
 * (\a is_ghost is false), then its \a quadid indexes the tree's quadrant
 * array; otherwise, it indexes the ghosts array. If the edge is hanging, then
 * the quadrants are listed in z-order. If an edge is in the interior of a
 * tree, orientation is 0; if an edge is between trees, orientation is the same
 * as edge orientation in the connectivity. If a quadrant should be present,
 * but it is not included in the ghost layer, then quad = NULL, is_ghost is
 * true, and quadid = -1.
                      *
 * the \a faces field provides some additional information about the local
 * neighborhood: if side[i]->faces[j] == side[k]->faces[l], this indicates that
 * there is a common face between these two sides of the edge.
 */
typedef struct p8est_iter_edge_side
{
  p4est_topidx_t      treeid;
  int8_t              edge;
  int8_t              orientation;
  int8_t              is_hanging;
  union p8est_iter_edge_side_data
  {
    struct
    {
      int8_t              is_ghost;
      p8est_quadrant_t   *quad;
      p4est_locidx_t      quadid;
    }
    full;
    struct
    {
      int8_t              is_ghost[2];
      p8est_quadrant_t   *quad[2];
      p4est_locidx_t      quadid[2];
    }
    hanging;
  }
  is;
  int8_t              faces[2];
}
p8est_iter_edge_side_t;

/** The information about all sides of an edge in the forest.
 * If tree_boundary is false, the edge is on the interior of a tree.
 * When tree_boundary is false, sides[0] contains the lowest z-order quadrant
 * that touches the edge.
 * When tree_boundary is true, its value is P8EST_CONNECT_FACE/EDGE
 * depending on the location of the edge relative to the tree.
 */
typedef struct p8est_iter_edge_info
{
  p8est_t            *p4est;
  p8est_ghost_t      *ghost_layer;
  int8_t              tree_boundary;
  sc_array_t          sides;    /* p8est_iter_edge_side_t */
}
p8est_iter_edge_info_t;

/** The prototype for a function that p8est_iterate will execute wherever
 * the edge is an edge of all quadrants that touch it i.e. the callback will
 * not execute on an edge the sits on a hanging face.
 *
 * Note: the forest must be edge balanced for p8est_iterate to execute a
 * callback function on edges.
 */
typedef void        (*p8est_iter_edge_t) (p8est_iter_edge_info_t * info,
                                          void *user_data);

/* Information about one side of a corner in the forest.  If a \a quad is local,
 * then its \a quadid indexes the tree's quadrant array; otherwise, it indexes
 * the ghosts array.
 *
 * the \a faces and \a edges field provides some additional information about
 * the local neighborhood: if side[i]->faces[j] == side[k]->faces[l], this
 * indicates that there is a common face between these two sides of the
 * corner.
 */
typedef struct p8est_iter_corner_side
{
  p4est_topidx_t      treeid;
  int8_t              corner;
  int8_t              is_ghost;
  p8est_quadrant_t   *quad;
  p4est_locidx_t      quadid;
  int8_t              faces[3];
  int8_t              edges[3];
}
p8est_iter_corner_side_t;

/** The information about all sides of a face in the forest.
 * If tree_boundary is false, the corner is on the interior of a tree.
 * When tree_boundary is false, sides[0] contains the lowest z-order quadrant
 * that touches the corner.
 * When tree_boundary is true, its value is P8EST_CONNECT_FACE/EDGE/CORNER
 * depending on the location of the corner relative to the tree.
 */
typedef struct p8est_iter_corner_info
{
  p8est_t            *p4est;
  p8est_ghost_t      *ghost_layer;
  int8_t              tree_boundary;
  sc_array_t          sides;    /* p8est_iter_corner_side_t */
}
p8est_iter_corner_info_t;

/** The prototype for a function that p8est_iterate will execute wherever
 * the corner is a corner for all quadrants that touch it i.e. the callback
 * will not execute on a corner that sits on a hanging face or edge.
 *
 * Note: the forest does not need to be corner balanced for p8est_iterate to
 * execute a callback function at corners, only face and edge balanced.
 */
typedef void        (*p8est_iter_corner_t) (p8est_iter_corner_info_t * info,
                                            void *user_data);

/** p8est_iterate executes the user-supplied callback functions at every
 * volume, face, edge and corner in the local forest. The ghost_layer may be
 * NULL. The \a user_data pointer is not touched by p8est_iterate, but is
 * passed to each of the callbacks. Any of the callback functions may be NULL.
 * The callback functions are interspersed with each other, i.e. some face
 * callbacks will occur between volume callbacks, and some edge callbacks will
 * occur between face callbacks, etc.:
 *
 * 1) volume callbacks occur in the sorted Morton-index order.
 * 2) a face callback is not executed until after the volume callbacks have
 *    been executed for the quadrants that share it.
 * 3) an edge callback is not executed until the face callbacks have been
 *    executed for all faces that touch the edge.
 * 4) a corner callback is not executed until the edge callbacks have been
 *    executed for all edges that touch the corner.
 * 5) it is not always the case that every face callback for a given quadrant
 *    is executed before any of the edge or corner callbacks, and it is not
 *    always the case that every edge callback for a given quadrant is executed
 *    before any of the corner callbacks.
 * 6) callbacks are not executed at faces, edges or corners that only involve
 *    ghost quadrants, i.e. that are not adjacent in the local section of the
 *    forest.
 */
void                p8est_iterate (p8est_t * p4est,
                                   p8est_ghost_t * ghost_layer,
                                   void *user_data,
                                   p8est_iter_volume_t iter_volume,
                                   p8est_iter_face_t iter_face,
                                   p8est_iter_edge_t iter_edge,
                                   p8est_iter_corner_t iter_corner);

/** Return a pointer to a iter_corner_side array element indexed by a int.
 */
/*@unused@*/
static inline p8est_iter_corner_side_t *
p8est_iter_cside_array_index_int (sc_array_t * array, int it)
{
  P4EST_ASSERT (array->elem_size == sizeof (p8est_iter_corner_side_t));
  P4EST_ASSERT (it >= 0 && (size_t) it < array->elem_count);

  return (p8est_iter_corner_side_t *)
    (array->array + sizeof (p8est_iter_corner_side_t) * it);
}

/** Return a pointer to a iter_corner_side array element indexed by a size_t.
 */
/*@unused@*/
static inline p8est_iter_corner_side_t *
p8est_iter_cside_array_index (sc_array_t * array, size_t it)
{
  P4EST_ASSERT (array->elem_size == sizeof (p8est_iter_corner_side_t));
  P4EST_ASSERT (it < array->elem_count);

  return (p8est_iter_corner_side_t *)
    (array->array + sizeof (p8est_iter_corner_side_t) * it);
}

/** Return a pointer to a iter_edge_side array element indexed by a int.
 */
/*@unused@*/
static inline p8est_iter_edge_side_t *
p8est_iter_eside_array_index_int (sc_array_t * array, int it)
{
  P4EST_ASSERT (array->elem_size == sizeof (p8est_iter_edge_side_t));
  P4EST_ASSERT (it >= 0 && (size_t) it < array->elem_count);

  return (p8est_iter_edge_side_t *)
    (array->array + sizeof (p8est_iter_edge_side_t) * it);
}

/** Return a pointer to a iter_edge_side array element indexed by a size_t.
 */
/*@unused@*/
static inline p8est_iter_edge_side_t *
p8est_iter_eside_array_index (sc_array_t * array, size_t it)
{
  P4EST_ASSERT (array->elem_size == sizeof (p8est_iter_edge_side_t));
  P4EST_ASSERT (it < array->elem_count);

  return (p8est_iter_edge_side_t *)
    (array->array + sizeof (p8est_iter_edge_side_t) * it);
}

/** Return a pointer to a iter_face_side array element indexed by a int.
 */
/*@unused@*/
static inline p8est_iter_face_side_t *
p8est_iter_fside_array_index_int (sc_array_t * array, int it)
{
  P4EST_ASSERT (array->elem_size == sizeof (p8est_iter_face_side_t));
  P4EST_ASSERT (it >= 0 && (size_t) it < array->elem_count);

  return (p8est_iter_face_side_t *)
    (array->array + sizeof (p8est_iter_face_side_t) * it);
}

/** Return a pointer to a iter_face_side array element indexed by a size_t.
 */
/*@unused@*/
static inline p8est_iter_face_side_t *
p8est_iter_fside_array_index (sc_array_t * array, size_t it)
{
  P4EST_ASSERT (array->elem_size == sizeof (p8est_iter_face_side_t));
  P4EST_ASSERT (it < array->elem_count);

  return (p8est_iter_face_side_t *)
    (array->array + sizeof (p8est_iter_face_side_t) * it);
}

SC_EXTERN_C_END;

#endif /* !P8EST_ITERATE_H */
