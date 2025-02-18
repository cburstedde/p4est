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
#include <p4est_vtk.h>
#else
#include <p8est_bits.h>
#include <p8est_extended.h>
#include <p8est_vtk.h>
#endif

/************ code used regardless of internal or external data **************/

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
    .25 * cos (4. * M_PI * coords[2]);
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

/* wrapper for an iteration over every local element */
static void
userdata_iterate_volume (p4est_userdata_global_t *g,
                         p4est_iter_volume_t volume_callback)
{
  /* iterate over local quadrants passing an arbitrary callback */
  P4EST_ASSERT (g->qcount == 0);
  p4est_iterate (g->p4est, NULL, g, volume_callback, NULL,
#ifdef P4_TO_P8
                 NULL,
#endif
                 NULL);
  P4EST_ASSERT (g->qcount == g->p4est->local_num_quadrants);
  g->qcount = 0;
}

/************ functions that expect user data internal to p4est **************/

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

/* callback for local elements' consistency check */
static void
userdata_verify_internal_volume (p4est_iter_volume_info_t *v,
                                 void *user_data)
{
  /* the global data structure is passed by the iterator */
  p4est_userdata_global_t *g = (p4est_userdata_global_t *) user_data;
#ifdef P4EST_ENABLE_DEBUG
  p4est_tree_t       *tree;
  userdata_quadrant_t *qdat;

  /* check call consistency */
  P4EST_ASSERT (v != NULL);
  P4EST_ASSERT (g != NULL);
  P4EST_ASSERT (v->p4est == g->p4est);

  /* verify quadrant user data */
  qdat = (userdata_quadrant_t *) v->quad->p.user_data;
  P4EST_ASSERT (qdat != NULL);
  P4EST_ASSERT (qdat->which_tree == v->treeid); /* iterator # is per-tree */
  tree = p4est_tree_array_index (g->p4est->trees, qdat->which_tree);
  P4EST_ASSERT (qdat->quadid == tree->quadrants_offset + v->quadid);
  P4EST_ASSERT (qdat->quadid == g->qcount);
#endif
  ++g->qcount;
}

/* run consistency check over all local elements' data */
static void
userdata_verify_internal (p4est_userdata_global_t *g)
{
  /* iterate over local quadrants and verify their data */
  P4EST_ASSERT (g != NULL);
  userdata_iterate_volume (g, userdata_verify_internal_volume);
}

/* callback for placing the element data in the VTK file */
static void
userdata_vtk_internal_volume (p4est_iter_volume_info_t *v, void *user_data)
{
  /* the global data structure is passed by the iterator */
  p4est_userdata_global_t *g = (p4est_userdata_global_t *) user_data;
  userdata_quadrant_t *qdat;

  /* check call consistency */
  P4EST_ASSERT (v != NULL);
  P4EST_ASSERT (g != NULL);
  P4EST_ASSERT (v->p4est == g->p4est);

  /* access quadrant user data */
  qdat = (userdata_quadrant_t *) v->quad->p.user_data;
  P4EST_ASSERT (qdat != NULL);
  P4EST_ASSERT (qdat->which_tree == v->treeid);
  P4EST_ASSERT (qdat->quadid == g->qcount);

  /* write quadrant value into output array and advance counter */
  *(double *) sc_array_index (g->qarray, g->qcount) = qdat->value;
  ++g->qcount;
}

/* provide function for consistent deallocation */
static int
userdata_vtk_internal_return (int retval, sc_array_t *farray)
{
  /* cleanup function context and return */
  if (farray != NULL) {
    sc_array_destroy (farray);
  }
  return retval;
}

/* write a VTK file with the mesh and the element data */
static int
userdata_vtk_internal (p4est_userdata_global_t *g, const char *filename)
{
  const char         *fnames[1] = { "value" };
  sc_array_t         *fvalues[1] = { NULL };
  p4est_vtk_context_t *vtk, *rvtk;

  P4EST_ASSERT (g != NULL);
  if (g->novtk) {
    /* output disabled */
    return userdata_vtk_internal_return (0, fvalues[0]);
  }

  /* ensure consistent error cleanup */
  fvalues[0] = sc_array_new (sizeof (double));

  /* write file in multpile steps */
  P4EST_ASSERT (g->p4est != NULL);
  vtk = p4est_vtk_context_new (g->p4est, filename);
  p4est_vtk_context_set_geom (vtk, g->geom);
  if ((rvtk = p4est_vtk_write_header (vtk)) == NULL) {
    P4EST_GLOBAL_LERRORF ("ERROR: write VTK header for %s\n", filename);
    return userdata_vtk_internal_return (-1, fvalues[0]);
  }
  P4EST_ASSERT (rvtk == vtk);

  /* the cell data must be gathered in a contiguous array of values */
  sc_array_resize (fvalues[0], g->p4est->local_num_quadrants);
  P4EST_ASSERT (g->qarray == NULL);
  g->qarray = fvalues[0];
  userdata_iterate_volume (g, userdata_vtk_internal_volume);
  g->qarray = NULL;
  if ((rvtk = p4est_vtk_write_cell_data (vtk, 1, 1, 1, 0, 1, 0,
                                         fnames, fvalues)) == NULL) {
    P4EST_GLOBAL_LERRORF ("ERROR: write VTK data for %s\n", filename);
    return userdata_vtk_internal_return (-1, fvalues[0]);
  }
  P4EST_ASSERT (rvtk == vtk);

  /* finalize the output files */
  if (p4est_vtk_write_footer (vtk)) {
    P4EST_GLOBAL_LERRORF ("ERROR: write VTK footer for %s\n", filename);
    return userdata_vtk_internal_return (-1, fvalues[0]);
  }
  vtk = NULL;

  /* return memory neutral */
  return userdata_vtk_internal_return (0, fvalues[0]);
}

/* callback to tell p4est which quadrants shall be refined */
static int
userdata_refine_internal (p4est_t *p4est, p4est_topidx_t which_tree,
                          p4est_quadrant_t *quadrant)
{
  int                 refine;

  /* the global data structure is stashed into the forest's user pointer */
  p4est_userdata_global_t *g =
    (p4est_userdata_global_t *) p4est->user_pointer;

  /* update the quadrant user data contents */
  userdata_quadrant_t *qdat = (userdata_quadrant_t *) quadrant->p.user_data;
  P4EST_ASSERT (qdat != NULL);

  /* refinement does not change the tree index */
  P4EST_ASSERT (qdat->which_tree == which_tree);

  /* placeholder for a proper refinement criterion */
  refine = ((which_tree % 3) == 0);
  if (!refine || quadrant->level >= g->maxlevel) {
    /* quadrant is unchanged, keep its data consistent */
    qdat->quadid = g->qcount++;
    return 0;
  }
  else {
    /* calculate new quadrant data in the replacement callback */
    return 1;
  }
}

/* helper function to update quadrant data when not coarsening */
static int
userdata_coarsen_internal_dont (p4est_userdata_global_t *g,
                                p4est_topidx_t which_tree,
                                p4est_quadrant_t *quadrant)
{
  /* we do not coarsen: this call is for proper counting */
  userdata_quadrant_t *qdat = (userdata_quadrant_t *) quadrant->p.user_data;
  P4EST_ASSERT (qdat != NULL);

  /* coarsening does not change the tree index */
  P4EST_ASSERT (qdat->which_tree == which_tree);

  /* the local quadrant index generally changes on adaptation */
  qdat->quadid = g->qcount++;

  /* save a line of code below */
  return 0;
}

/* callback to tell p4est which quadrants shall be coarsened */
static int
userdata_coarsen_internal (p4est_t *p4est, p4est_topidx_t which_tree,
                           p4est_quadrant_t *quadrant[])
{
  int                 coarsen;

  /* the global data structure is stashed into the forest's user pointer */
  p4est_userdata_global_t *g =
    (p4est_userdata_global_t *) p4est->user_pointer;

  /* we do not coarsen, we are only passed one quadrant */
  if (quadrant[1] == NULL) {
    return userdata_coarsen_internal_dont (g, which_tree, quadrant[0]);
  }

  /* really should determine proper coarsening criterion */
  coarsen = ((which_tree % 3) == 1);
  if (!coarsen) {
    /* we decided not to coarsen, just update the first quadrant */
    return userdata_coarsen_internal_dont (g, which_tree, quadrant[0]);
  }
  else {
    /* calculate new quadrant data in the replacement callback */
    return 1;
  }
}

/* update element data for a family of quadrants during adaptation */
static void
userdata_replace_internal (p4est_t *p4est, p4est_topidx_t which_tree,
                           int num_outgoing, p4est_quadrant_t *outgoing[],
                           int num_incoming, p4est_quadrant_t *incoming[])
{
  int                 i;
  double              sum;

  /* the global data structure is stashed into the forest's user pointer */
  p4est_userdata_global_t *g =
    (p4est_userdata_global_t *) p4est->user_pointer;

  /* check context invariant */
  P4EST_ASSERT (g->in_balance || g->bcount == 0);

  if (num_incoming > 1) {
    p4est_locidx_t      addcount;

    /* we are refining */
    P4EST_ASSERT (num_outgoing == 1);
    P4EST_ASSERT (num_incoming == P4EST_CHILDREN);

    /* access old (larger) quadrant's data */
    userdata_quadrant_t *qold =
      (userdata_quadrant_t *) outgoing[0]->p.user_data;
    P4EST_ASSERT (qold != NULL);
    P4EST_ASSERT (qold->which_tree == which_tree);

    /* within 2:1 balance, we do not have an iterator invariant */
    if (!g->in_balance) {
      addcount = g->qcount;
      g->qcount += P4EST_CHILDREN;
    }
    else {
      /* determine the count by accessing the outgoing quadrant */
      addcount = qold->quadid + g->bcount;
      g->bcount += P4EST_CHILDREN - 1;
    }

    /* access new (smaller) quadrants' data */
    for (i = 0; i < P4EST_CHILDREN; ++i) {
      userdata_quadrant_t *qnew =
        (userdata_quadrant_t *) incoming[i]->p.user_data;
      P4EST_ASSERT (qnew != NULL);
      qnew->which_tree = which_tree;
      qnew->quadid = addcount + i;

      /* we just copy the old value into the refined elements */
      qnew->value = qold->value;
    }
  }
  else {
    /* we are coarsening */
    P4EST_ASSERT (!g->in_balance);
    P4EST_ASSERT (num_outgoing == P4EST_CHILDREN);
    P4EST_ASSERT (num_incoming == 1);

    /* access new (larger) quadrant's data */
    userdata_quadrant_t *qnew =
      (userdata_quadrant_t *) incoming[0]->p.user_data;
    P4EST_ASSERT (qnew != NULL);
    qnew->which_tree = which_tree;
    qnew->quadid = g->qcount++;

    /* access old (smaller) quadrants' data */
    sum = 0.;
    for (i = 0; i < P4EST_CHILDREN; ++i) {
      userdata_quadrant_t *qold =
        (userdata_quadrant_t *) outgoing[i]->p.user_data;
      P4EST_ASSERT (qold != NULL);
      P4EST_ASSERT (qold->which_tree == which_tree);

      /* just copy the old value into the refined elements */
      sum += qold->value;
    }

    /* we just overage the old values into the coarsened element */
    qnew->value = sum / P4EST_CHILDREN;
  }
}

/* provide function for consistent deallocation */
static int
userdata_run_internal_return (int retval, p4est_userdata_global_t *g)
{
  P4EST_ASSERT (g != NULL);

  /* clean up what has been allocated in the same function */
  if (g->p4est != NULL) {
    /* destroy forest */
    p4est_destroy (g->p4est);
    g->p4est = NULL;
  }
  return retval;
}

/* post-balance callback to update the local elements count */
static void
userdata_balance_internal_volume (p4est_iter_volume_info_t *v,
                                  void *user_data)
{
  /* the global data structure is passed by the iterator */
  p4est_userdata_global_t *g = (p4est_userdata_global_t *) user_data;
  userdata_quadrant_t *qdat;

  /* check call consistency */
  P4EST_ASSERT (v != NULL);
  P4EST_ASSERT (g != NULL);
  P4EST_ASSERT (g->in_balance);
  P4EST_ASSERT (v->p4est == g->p4est);

  /* verify quadrant user data */
  qdat = (userdata_quadrant_t *) v->quad->p.user_data;
  P4EST_ASSERT (qdat != NULL);
  P4EST_ASSERT (qdat->which_tree == v->treeid);

  /* the quadrant id is only correct for newly created quadrants,
     otherwise it may be short since it has the pre-balance value */
  P4EST_ASSERT (qdat->quadid <= g->qcount);
  qdat->quadid = g->qcount++;
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
  userdata_verify_internal (g);

  /* write VTK files to visualize geometry and data */
  if (userdata_vtk_internal (g, P4EST_STRING "_userdata_internal_new")) {
    P4EST_GLOBAL_LERROR ("ERROR: write VTK output for forest_new\n");
    return userdata_run_internal_return (-1, g);
  }

  /* refine the mesh adaptively and non-recursively */
  P4EST_ASSERT (g->qcount == 0);
  p4est_refine_ext (g->p4est, 0, g->maxlevel, userdata_refine_internal,
                    NULL, userdata_replace_internal);
  P4EST_ASSERT (g->qcount == g->p4est->local_num_quadrants);
  g->qcount = 0;
  userdata_verify_internal (g);

  /* coarsen the mesh adaptively and non-recursively */
  P4EST_ASSERT (g->qcount == 0);
  p4est_coarsen_ext (g->p4est, 0, 1, userdata_coarsen_internal,
                     NULL, userdata_replace_internal);
  P4EST_ASSERT (g->qcount == g->p4est->local_num_quadrants);
  g->qcount = 0;
  userdata_verify_internal (g);

  /* execute 2:1 balance on the mesh */
  P4EST_ASSERT (g->qcount == 0);
  P4EST_ASSERT (g->bcount == 0);
  P4EST_ASSERT (!g->in_balance);
  g->in_balance = 1;
  p4est_balance_ext (g->p4est, P4EST_CONNECT_FULL,
                     NULL, userdata_replace_internal);

  /* we invoke the volume iteration since we keep the local index
     in the quadrant's data for demonstration purposes */
  userdata_iterate_volume (g, userdata_balance_internal_volume);
  g->bcount = 0;
  g->in_balance = 0;
  userdata_verify_internal (g);

  /* return memory neutral */
  return userdata_run_internal_return (0, g);
}

/************ functions that expect user data external to p4est **************/

/* core demo with quadrant data stored external to p4est */
static int
userdata_run_external (p4est_userdata_global_t *g)
{
  /* allocating and maintaining userdata outside of p4est works, too */
  /* this part of the demo is still to be written */
  return 0;
}

/* *********** this function is the entry point into this file ***************/

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
