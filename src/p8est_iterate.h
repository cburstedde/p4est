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

#ifndef P8EST_ITERATE_H
#define P8EST_ITERATE_H

#include <p8est.h>

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
  sc_array_t         *ghost_layer;
  p8est_quadrant_t   *quad;
  size_t              quadid;
  p4est_topidx_t      treeid;
}
p8est_iter_volume_info_t;

/** The prototype for a function that p8est_iterate will execute at every
 * quadrant local to the current process.
 */
typedef void        (*p8est_iter_volume_t) (p8est_iter_volume_info_t * info,
                                            void *user_data);

/** The information that is available to the user defined p8est_iter_face_t
 * callback function about the quadrants on either side of a face.
 * \a left_treeid and \a right_treeid are as above for the volume case.
 * \a left_quadid and/or \a right_quadid may be negative, which indicates that
 *    the quadrant is not local, and resides in the \a ghost_layer at the index
 *    quadid + ghost_layer->elem_count;
 * \a left_outgoing_face is the **z-order** face of \a left_quad that touches
 *    \a right_quad, and vice versa.
 * if \a is_hanging is true, then that indicates that \a left_quad is larger
 *    than right quad.
 * if \a is_hanging is true, then \a left_corner is the **z-order** corner that
 *    touches \a right_corner.  Said another way, of \a left_quad's children,
 *    \a left_corner is adjacent to right_quad, and \a right_corner is
 *    \a right_quad's child_id. When \a is_hanging is false, \a left_corner and
 *    \a right_corner are just one of the two sets of matching corners on the
 *    face.
 *
 * Note: at the outside boundary of the p4est, if it is not-periodic, a
 * quadrants face neighbor is considered to be itself.
 */
typedef struct p8est_iter_face_info
{
  p8est_t            *p4est;
  sc_array_t         *ghost_layer;
  p8est_quadrant_t   *left_quad;
  ssize_t             left_quadid;
  p4est_topidx_t      left_treeid;
  int                 left_outgoing_face;
  int                 left_corner;
  p8est_quadrant_t   *right_quad;
  ssize_t             right_quadid;
  p4est_topidx_t      right_treeid;
  int                 right_outgoing_face;
  int                 right_corner;
  int                 orientation;
  bool                is_hanging;
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

/** The information that is available to the user defined p8est_iter_edge_t
 * callback function about the quadrants surrounding an edge.
 * There may be a variable number of quadrants.  If \a is_hanging is true, then
 * some of the quadrants around the edge are half the size of the the others.
 * However, if this is the case, then they all meet at a corner, and
 * \a common_corners gives the corner id for the shared corner from each
 * quadrant.  As above, a\ quadids may be positive or negative.
 */
typedef struct p8est_iter_edge_info
{
  p8est_t            *p4est;
  sc_array_t         *ghost_layer;
  sc_array_t         *quads;          /** elements are (p8est_quadrant_t *) */
  sc_array_t         *quadids;        /** elements are (ssize_t)            */
  sc_array_t         *treeids;        /** elements are (p4est_locidx_t)     */
  sc_array_t         *edges;          /** elements are (int)                */
  sc_array_t         *common_corners; /** elements are (int)                */
  bool                is_hanging;
}
p8est_iter_edge_info_t;

/** The prototype for a function that p8est_iterate will execute wherever
 * quadrants meet at a conformal edge i.e. the callback will not execute on
 * an edge the sits on a hanging face.
 *
 * Note: the forest must be edge balanced for p8est_iterate to execute a
 * callback function on edges.
 */
typedef void        (*p8est_iter_edge_t) (p8est_iter_edge_info_t * info,
                                          void *user_data);

/** The information that is available to the user defined p8est_iter_corner_t
 * callback function about the quadrants surrounding a corner.
 * There may be a variable number of quadrants. \a corners gives the **z-order**
 * corner that touches the other shared corners, all other information is as
 * above.
 */
typedef struct p8est_iter_corner_info
{
  p8est_t            *p4est;
  sc_array_t         *ghost_layer;
  sc_array_t         *quads;      /** elements are (p8est_quadrant_t *) */
  sc_array_t         *quadids;    /** elements are (ssize_t)            */
  sc_array_t         *treeids;    /** elements are (p4est_locidx_t)     */
  sc_array_t         *corners;    /** elements are (int)                */
}
p8est_iter_corner_info_t;

/** The prototype for a function that p8est_iterate will execute wherever
 * quadrants meet at a conformal corner i.e. the callback will not execute on
 * a corner that sits on a hanging face or edge.
 *
 * Note: the forest does not need to be corner balanced for p8est_iterate to
 * execute a callback function at corners, only face and edge balance.
 * However, the ghost_layer must be created with the P4EST_BALANCE_FULL option.
 */
typedef void        (*p8est_iter_corner_t) (p8est_iter_corner_info_t * info,
                                            void *user_data);

/** p8est_iterate executes the user-supplied callback functions at every
 * volume, face, edge and corner in the local forest. The \a user_data pointer
 * is not touched by p4est_iterate, but is passed to each of the callbacks.
 * The callback functions are interspersed with each other, i.e. some face
 * callbacks will occur between volume callbacks, and some edge callbacks
 * will occur between face callbacks, etc.:
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
void                p8est_iterate (p8est_t * p4est, sc_array_t * ghost_layer,
                                   void *user_data,
                                   p8est_iter_volume_t iter_volume,
                                   p8est_iter_face_t iter_face,
                                   p8est_iter_edge_t iter_edge,
                                   p8est_iter_corner_t iter_corner);

SC_EXTERN_C_END;

#endif /* !P8EST_ITERATE_H */
