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

#ifndef P4EST_ITERATE_H
#define P4EST_ITERATE_H

#include <p4est.h>

SC_EXTERN_C_BEGIN;

/** The information that is available to the user defined p4est_iter_volume_t
 * callback function.
 *
 * \a treeid gives the index in \a p4est->trees of the the tree to which
 *    \a quad belongs.
 * \a quadid gives the index of \a quad within \a tree's quadrants array.
 */
typedef struct p4est_iter_volume_info
{
  p4est_t            *p4est;
  sc_array_t         *ghost_layer;
  p4est_quadrant_t   *quad;
  size_t              quadid;
  p4est_topidx_t      treeid;
}
p4est_iter_volume_info_t;

/** The prototype for a function that p4est_iterate will execute at every
 * quadrant local to the current process.
 */
typedef void        (*p4est_iter_volume_t) (p4est_iter_volume_info_t * info,
                                            void *user_data);

/** The information that is available to the user defined p4est_iter_face_t
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
typedef struct p4est_iter_face_info
{
  p4est_t            *p4est;
  sc_array_t         *ghost_layer;
  p4est_quadrant_t   *left_quad;
  ssize_t             left_quadid;
  p4est_topidx_t      left_treeid;
  int                 left_outgoing_face;
  int                 left_corner;
  p4est_quadrant_t   *right_quad;
  ssize_t             right_quadid;
  p4est_topidx_t      right_treeid;
  int                 right_outgoing_face;
  int                 right_corner;
  int                 orientation;
  bool                is_hanging;
}
p4est_iter_face_info_t;

/** The prototype for a function that p4est_iterate will execute wherever two
 * quadrants share a face: the face can be a 2:1 hanging face, it does not have
 * to be conformal.
 *
 * Note: the forest must be face balanced for p4est_iterate to execute a
 * callback function on faces.
 */
typedef void        (*p4est_iter_face_t) (p4est_iter_face_info_t * info,
                                          void *user_data);

/** The information that is available to the user defined p4est_iter_corner_t
 * callback function about the quadrants surrounding a corner.
 * There may be a variable number of quadrants. \a corners gives the **z-order**
 * corner that touches the other shared corners, all other information is as
 * above.
 */
typedef struct p4est_iter_corner_info
{
  p4est_t            *p4est;
  sc_array_t         *ghost_layer;
  sc_array_t         *quads;      /** elements are (p4est_quadrant_t *) */
  sc_array_t         *quadids;    /** elements are (ssize_t)            */
  sc_array_t         *treeids;    /** elements are (p4est_locidx_t)     */
  sc_array_t         *corners;    /** elements are (int)                */
}
p4est_iter_corner_info_t;

/** The prototype for a function that p4est_iterate will execute wherever
 * quadrants meet at a conformal corner i.e. the callback will not execute on
 * a hanging corner.
 *
 * Note: the forest does not need to be corner balanced for p4est_iterate to
 * execute a callback function at corners, only face balanced. However, the
 * ghost_layer must be created with the P4EST_BALANCE_FULL option if there is a
 * corner callback function.
 */
typedef void        (*p4est_iter_corner_t) (p4est_iter_corner_info_t * info,
                                            void *user_data);

/** p4est_iterate executes the user-supplied call back functions at every
 * volume, face, and corner in the local forest. The \a user_data pointer is
 * not touched by p4est_iterate, but is passed to each of the callbakcs.  The
 * callback functions are interspersed with each other, i.e. some face
 * callbacks will occur between volume call backs, and some corner callbacks
 * will occur between face callbacks:
 *
 * 1) volume callbacks occur in the sorted morton-index order.
 * 2) a face callback is not executed until after the volume callbacks have
 *    been executed for the quadrants that share it.
 * 3) a corner callback is not executed until the face callbacks have been
 *    executed for all faces that touch the corner.
 * 4) it is not always the case that every face callback for a given quadrant
 *    is executed before any of the corner callbacks.
 * 5) callbacks are not executed at faces or corners that only involve ghost
 *    quadrants, i.e. that are not adjacency the the local section of the
 *    forest.
 */
void                p4est_iterate (p4est_t * p4est, sc_array_t * ghost_layer,
                                   void *user_data,
                                   p4est_iter_volume_t iter_volume,
                                   p4est_iter_face_t iter_face,
                                   p4est_iter_corner_t iter_corner);

SC_EXTERN_C_END;

#endif /* !P4EST_ITERATE_H */
