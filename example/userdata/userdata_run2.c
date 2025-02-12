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

/* the core functionality of the exmple program, for both 2D and 3D */
#include "userdata_global.h"
#ifndef P4_TO_P8
#include <p4est_bits.h>
#include <p4est_extended.h>
#else
#include <p8est_bits.h>
#include <p8est_extended.h>
#endif

/* demonstration data for each quadrant */
typedef struct userdata_quadrant
{
  /* The tree number is implicit in the forest, no need to store it.
     We do this here just for demonstration purposes. */
  p4est_topidx_t      which_tree;

  /* The quadrant number is implicit in the forest iterators, no need to
     store it.  We do this here just for demonstration purposes. */
  p4est_locidx_t      quadid;

  /* Store a piecewise constant variable field. */
  double              value;
}
userdata_quadrant_t;

/* analytic solution function for demonstration purposes */
static double
userdata_analytic (p4est_userdata_global_t *g, const double coords[3])
{
  return
    .5 * sin (M_PI * coords[0]) +
    .25 * exp (-.5 * pow ((coords[1] - .5) / 2.5, 2.)) +
    .25 * cos (4 * M_PI * coords[2]);
}

/* compute the value of the given function at quadrant midpoint */
static double
userdata_value (p4est_userdata_global_t *g,
                p4est_topidx_t which_tree, p4est_quadrant_t *quadrant)
{
  p4est_qcoord_t      coords_in[P4EST_DIM];
  double              coords_out[3];

  /* transform quadrant midpoint into the geometry system */
  p4est_quadrant_volume_coordinates (quadrant, coords_in);
  p4est_geometry_transform_coordinates (g->geom, which_tree,
                                        coords_in, coords_out);

  /* evaluate function at this point */
  return userdata_analytic (g, coords_out);
}

/* callback to initialize internal quadrant data */
static void
userdata_init_internal (p4est_t *p4est,
                        p4est_topidx_t which_tree,
                        p4est_quadrant_t * quadrant)
{
  /* the global data structure is stashed into the forest's user pointer */
  p4est_userdata_global_t *g =
    (p4est_userdata_global_t *) p4est->user_pointer;

  /* p4est is agnostic to the quadrant user data */
  userdata_quadrant_t *qdat = (userdata_quadrant_t *) quadrant->p.user_data;

  /* exemplarily populate some quadrant data */
  qdat->which_tree = which_tree;

  /* we know that the callback is executed in order for every quadrant */
  qdat->quadid = g->qcount++;

  /* compute field value from analytic expression */
  qdat->value = userdata_value (g, which_tree, quadrant);
}

/* core demo with quadrant data stored internal to p4est */
static int
userdata_run_internal (p4est_userdata_global_t *g)
{
  P4EST_ASSERT (g->p4est == NULL);

  /* create initial forest and populate quadrant data by callback */
  P4EST_ASSERT (g->qcount == 0);
  g->p4est = p4est_new_ext
    (g->mpicomm, g->conn, 0, SC_MAX (g->maxlevel - 1, 0), 1,
     sizeof (userdata_quadrant_t), userdata_init_internal, g);
  P4EST_ASSERT (g->qcount == g->p4est->local_num_quadrants);
  g->qcount = 0;

  /* destroy forest */
  p4est_destroy (g->p4est);
  g->p4est = NULL;
  return 0;
}

/* core demo with quadrant data stored external to p4est */
static int
userdata_run_external (p4est_userdata_global_t *g)
{
  return 0;
}

/* execute the demonstration */
int
p4est_userdata_run (p4est_userdata_global_t *g)
{
  int                 erres;

  P4EST_ASSERT (g != NULL);
  P4EST_ASSERT (g->options != NULL);
  P4EST_ASSERT (g->conn != NULL);
  P4EST_ASSERT (p4est_connectivity_is_valid (g->conn));

  /* for consistency, track error status the same way */
  erres = 0;

  /* run the example once with p4est-allocated application data */
  if (!erres && (erres = userdata_run_internal (g))) {
    P4EST_GLOBAL_LERROR ("ERROR: run with internal data\n");
  }

  /* run the example another time with user-allocated application data */
  if (!erres && (erres = userdata_run_external (g))) {
    P4EST_GLOBAL_LERROR ("ERROR: run with external data\n");
  }

  /* return error status */
  return erres;
}
