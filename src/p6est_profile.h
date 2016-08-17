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

#ifndef P6EST_PROFILE_H
#define P6EST_PROFILE_H

#include <p6est.h>
#include <p6est_ghost.h>
#include <p4est_lnodes.h>
#include <p6est_lnodes.h>

SC_EXTERN_C_BEGIN;

/** A p6est profile is used to (a) balance a p6est, and (b) generate a
 * p6est_lnodes.  In every case, layers in a column are compressed to one
 * int8_t each.  The resulting column profiles and be quickly intersected,
 * unioned and balanced.  We use a p4est_lnodes_t to communicate between
 * neighboring columns, and a p4est_ghost_t to communicate between neighboring
 * processes.
 */

typedef enum
{
  P6EST_PROFILE_UNION,
  P6EST_PROFILE_INTERSECTION
}
p6est_profile_type_t;

typedef struct p6est_profile
{
  p6est_profile_type_t ptype;
  p8est_connect_type_t btype;
  p4est_lnodes_t     *lnodes;
  p4est_ghost_t      *cghost;
  int                 ghost_owned;
  p4est_locidx_t     *lnode_ranges;
  sc_array_t         *lnode_columns;
  int                *lnode_changed[2];
  p4est_locidx_t     *enode_counts;
  int                 evenodd;
  p4est_qcoord_t      diff;
}
p6est_profile_t;

/** Create a new profile.
 * \param[in] p6est
 * \param[in] ghost optional, can be NULL
 * \param[in] ptype P6EST_PROFILE_UNION if we are balancing,
 *                  P6EST_PROFILE_INTERSECTION if we are generating lnodes
 * \param[in] btype Type of 3D balance law desired.
 * \param[in] degree degree of underlying lnodes, should be 2 if used for
 *                   balancing
 */

p6est_profile_t    *p6est_profile_new_local (p6est_t * p6est,
                                             p6est_ghost_t * ghost,
                                             p6est_profile_type_t ptype,
                                             p8est_connect_type_t btype,
                                             int degree);

/** Destroy a profile */
void                p6est_profile_destroy (p6est_profile_t * profile);

/** Enforce balance between the column profiles locally: no communication */
void                p6est_profile_balance_local (p6est_profile_t * profile);

/** Synchronize the data from other processors, taking unions or
 * intersections, as determined at profile creation in \a
 * p6est_profile_new_local
 *
 * \return whether any change has occured.
 * */
int                 p6est_profile_sync (p6est_profile_t * profile);

/** modify a p6est, whose columns match the column profiles, to match the
 * column profiles */
void                p6est_refine_to_profile (p6est_t * p6est,
                                             p6est_profile_t * profile,
                                             p6est_init_t init_fn,
                                             p6est_replace_t replace_fn);

void                p6est_profile_element_to_node (p6est_t * p6est,
                                                   p6est_profile_t * profile,
                                                   p4est_locidx_t * offsets,
                                                   p4est_locidx_t *
                                                   elem_to_node,
                                                   p6est_lnodes_code_t * fc);

SC_EXTERN_C_END;

#endif /* !P6EST_PROFILE_H */
