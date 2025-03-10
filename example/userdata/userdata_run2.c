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
#include <p4est_communication.h>
#include <p4est_extended.h>
#include <p4est_vtk.h>
#else
#include <p8est_bits.h>
#include <p8est_communication.h>
#include <p8est_extended.h>
#include <p8est_vtk.h>
#endif

/************ code used regardless of internal or external data **************/

/* demonstration data for each quadrant */
typedef struct userdata_quadrant_internal
{
  /* Store a piecewise constant variable field. */
  double              value;
}
userdata_quadrant_internal_t;

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

/* write a VTK file with the mesh and the element data */
static int
userdata_vtk_general (p4est_userdata_global_t *g, const char *filename)
{
  const char         *fnames[1] = { "value" };
  sc_array_t         *fvalues[1] = { NULL };
  p4est_vtk_context_t *vtk, *rvtk;

  /* check preconditions */
  P4EST_ASSERT (g != NULL);
  P4EST_ASSERT (g->p4est != NULL);
  P4EST_ASSERT (!g->novtk);

  /* verify structure of VTK data */
  P4EST_ASSERT (g->qarray != NULL);
  P4EST_ASSERT (g->qarray->elem_size == sizeof (double));
  P4EST_ASSERT (g->qarray->elem_count ==
                (size_t) g->p4est->local_num_quadrants);
  fvalues[0] = g->qarray;

  /* write file in multpile steps */
  P4EST_ASSERT (g->p4est != NULL);
  vtk = p4est_vtk_context_new (g->p4est, filename);
  p4est_vtk_context_set_geom (vtk, g->geom);
  if ((rvtk = p4est_vtk_write_header (vtk)) == NULL) {
    P4EST_GLOBAL_LERRORF ("ERROR: write VTK header for %s\n", filename);
    return -1;
  }
  P4EST_ASSERT (rvtk == vtk);

  /* the cell data is written from a contiguous array of values */
  if ((rvtk = p4est_vtk_write_cell_data (vtk, 1, 1, 1, 0, 1, 0,
                                         fnames, fvalues)) == NULL) {
    P4EST_GLOBAL_LERRORF ("ERROR: write VTK data for %s\n", filename);
    return -1;
  }
  P4EST_ASSERT (rvtk == vtk);

  /* finalize the output files and deallocate context */
  if (p4est_vtk_write_footer (vtk)) {
    P4EST_GLOBAL_LERRORF ("ERROR: write VTK footer for %s\n", filename);
    return -1;
  }

  /* return success */
  return 0;
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
  userdata_quadrant_internal_t *qdat =
    (userdata_quadrant_internal_t *) quadrant->p.user_data;

  /* compute field value from analytic expression */
  qdat->value = userdata_value (g, which_tree, quadrant);

  /* we know that the callback is executed in order for every quadrant */
  g->qcount++;
}

/* callback for placing the element data in the VTK file */
static void
userdata_vtk_internal_volume (p4est_iter_volume_info_t *v, void *user_data)
{
  /* the global data structure is passed by the iterator */
  p4est_userdata_global_t *g = (p4est_userdata_global_t *) user_data;
  userdata_quadrant_internal_t *qdat;

  /* check call consistency */
  P4EST_ASSERT (v != NULL);
  P4EST_ASSERT (g != NULL);
  P4EST_ASSERT (v->p4est == g->p4est);
  P4EST_ASSERT (g->qarray != NULL);

  /* access quadrant user data */
  qdat = (userdata_quadrant_internal_t *) v->quad->p.user_data;
  P4EST_ASSERT (qdat != NULL);

  /* write quadrant value into output array and advance counter */
  *(double *) sc_array_index (g->qarray, g->qcount) = qdat->value;
  ++g->qcount;
}

/* write a VTK file with the mesh and the element data */
static int
userdata_vtk_internal (p4est_userdata_global_t *g, const char *filename)
{
  P4EST_ASSERT (g != NULL);
  P4EST_ASSERT (g->p4est != NULL);
  P4EST_ASSERT (g->qarray == NULL);
  P4EST_ASSERT (g->in_internal);
  if (g->novtk) {
    /* output disabled */
    return 0;
  }

  /* populate temporary array for output data */
  g->qarray = sc_array_new_count (sizeof (double), (size_t)
                                  g->p4est->local_num_quadrants);
  userdata_iterate_volume (g, userdata_vtk_internal_volume);
  if (userdata_vtk_general (g, filename)) {
    P4EST_GLOBAL_LERRORF ("ERROR: write VTK file %s\n", filename);
    sc_array_destroy_null (&g->qarray);
    return -1;
  }

  /* return memory neutral */
  sc_array_destroy_null (&g->qarray);
  return 0;
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
  userdata_quadrant_internal_t *qdat =
    (userdata_quadrant_internal_t *) quadrant->p.user_data;
  P4EST_ASSERT (qdat != NULL);

  /* placeholder for a proper refinement criterion */
  refine = ((which_tree % 3) == 0) || (fabs (qdat->value) > .8);

  if (!refine || quadrant->level >= g->maxlevel) {
    /* quadrant is unchanged, keep counting */
    g->qcount++;
    return 0;
  }
  else {
    /* calculate new quadrant data in the replacement callback */
    return 1;
  }
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
    ++g->qcount;
    return 0;
  }

  /* really should determine proper coarsening criterion */
  coarsen = ((which_tree % 3) == 1);
  if (!coarsen) {
    /* we decided not to coarsen this family of quadrants */
    ++g->qcount;
    return 0;
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
    /* we are refining */
    P4EST_ASSERT (num_outgoing == 1);
    P4EST_ASSERT (num_incoming == P4EST_CHILDREN);

    /* access old (larger) quadrant's data */
    userdata_quadrant_internal_t *qold =
      (userdata_quadrant_internal_t *) outgoing[0]->p.user_data;
    P4EST_ASSERT (qold != NULL);

    /* determine offset for new quadrants' indices */
    if (!g->in_balance) {
      /* we rely on the iterator property of the refinement algorithm */
      g->qcount += P4EST_CHILDREN;
    }
    else {
      /* within 2:1 balance, we do not have an iterator property;
         determine the new index by accessing the outgoing quadrant */
      g->bcount += P4EST_CHILDREN - 1;
    }

    /* access new (smaller) quadrants' data */
    for (i = 0; i < P4EST_CHILDREN; ++i) {
      userdata_quadrant_internal_t *qnew =
        (userdata_quadrant_internal_t *) incoming[i]->p.user_data;
      P4EST_ASSERT (qnew != NULL);

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
    userdata_quadrant_internal_t *qnew =
      (userdata_quadrant_internal_t *) incoming[0]->p.user_data;
    P4EST_ASSERT (qnew != NULL);
    g->qcount++;

    /* access old (smaller) quadrants' data */
    sum = 0.;
    for (i = 0; i < P4EST_CHILDREN; ++i) {
      userdata_quadrant_internal_t *qold =
        (userdata_quadrant_internal_t *) outgoing[i]->p.user_data;
      P4EST_ASSERT (qold != NULL);

      /* just copy the old value into the refined elements */
      sum += qold->value;
    }

    /* we just overage the old values into the coarsened element */
    qnew->value = sum * (1. / P4EST_CHILDREN);
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

/* core demo with quadrant data stored internal to p4est */
static int
userdata_run_internal (p4est_userdata_global_t *g)
{
  P4EST_ASSERT (g != NULL);
  P4EST_ASSERT (g->p4est == NULL);
  P4EST_ASSERT (g->in_internal);

  /* this is the demo for userdata allocated by p4est and stored internally */
  P4EST_GLOBAL_PRODUCTION (P4EST_STRING
                           "_userdata: application data INTERNAL\n");

  /* create initial forest and populate quadrant data by callback */
  P4EST_ASSERT (g->qcount == 0);
  g->p4est = p4est_new_ext
    (g->mpicomm, g->conn, 0, SC_MAX (g->maxlevel - 1, 0), 1,
     sizeof (userdata_quadrant_internal_t), userdata_init_internal, g);
  P4EST_ASSERT (g->qcount == g->p4est->local_num_quadrants);
  g->qcount = 0;

  /* We like the invariant that after partitioning, partition-independent
     coarsening is always possible since every family of siblings is placed
     on a single process; this is not guaranteed by p4est_new_ext. */
  p4est_partition (g->p4est, 1, NULL);

  /* write VTK files to visualize geometry and data */
  if (userdata_vtk_internal (g, P4EST_STRING "_userdata_internal_new")) {
    P4EST_GLOBAL_LERROR ("ERROR: write VTK output after forest_new\n");
    return userdata_run_internal_return (-1, g);
  }

  /* refine the mesh adaptively and non-recursively */
  P4EST_ASSERT (g->qcount == 0);
  p4est_refine_ext (g->p4est, 0, g->maxlevel, userdata_refine_internal,
                    NULL, userdata_replace_internal);
  P4EST_ASSERT (g->qcount == g->p4est->local_num_quadrants);
  g->qcount = 0;

  /* coarsen the mesh adaptively and non-recursively */
  P4EST_ASSERT (g->qcount == 0);
  p4est_coarsen_ext (g->p4est, 0, 1, userdata_coarsen_internal,
                     NULL, userdata_replace_internal);
  P4EST_ASSERT (g->qcount == g->p4est->local_num_quadrants);
  g->qcount = 0;

  /* execute 2:1 balance on the mesh */
  P4EST_ASSERT (g->qcount == 0);
  P4EST_ASSERT (g->bcount == 0);
  P4EST_ASSERT (!g->in_balance);
  g->in_balance = 1;
  p4est_balance_ext (g->p4est, P4EST_CONNECT_FULL,
                     NULL, userdata_replace_internal);
  g->bcount = 0;
  g->in_balance = 0;

  /* write VTK files to visualize geometry and data after adaptation */
  if (userdata_vtk_internal (g, P4EST_STRING "_userdata_internal_adapt")) {
    P4EST_GLOBAL_LERROR ("ERROR: write VTK output after adaptation\n");
    return userdata_run_internal_return (-1, g);
  }

  /* repartition the mesh and the quadrant data */
  p4est_partition (g->p4est, 1, NULL);

  /* write VTK files to visualize geometry and data after partitioning */
  if (userdata_vtk_internal (g, P4EST_STRING "_userdata_internal_partition")) {
    P4EST_GLOBAL_LERROR ("ERROR: write VTK output after adaptation\n");
    return userdata_run_internal_return (-1, g);
  }

  /* return memory neutral */
  return userdata_run_internal_return (0, g);
}

/************ functions that expect user data external to p4est **************/

/* metadata for each quadrant */
typedef struct userdata_quadrant_external
{
  /* refinement flag */
  int8_t              refine;

  /* coarsening flag */
  int8_t              coarsen;
}
userdata_quadrant_external_t;

/* callback to initialize quadrant metadata */
static void
userdata_init_external (p4est_t *p4est,
                        p4est_topidx_t which_tree,
                        p4est_quadrant_t * quadrant)
{
  /* the global data structure is stashed into the forest's user pointer */
  p4est_userdata_global_t *g =
    (p4est_userdata_global_t *) p4est->user_pointer;

  /* p4est is agnostic to the quadrant user data */
  userdata_quadrant_external_t *qdat =
    (userdata_quadrant_external_t *) quadrant->p.user_data;

  /* zero all bytes to initialize compiler padding */
  memset (qdat, 0, sizeof (*qdat));

  /* advance debug counter */
  ++g->qcount;
}

/* callback for local elements' initialization */
static void
userdata_init_external_volume (p4est_iter_volume_info_t *v, void *user_data)
{
  /* the global data structure is passed into every iterator */
  p4est_userdata_global_t *g = (p4est_userdata_global_t *) user_data;

  /* check call consistency */
  P4EST_ASSERT (v != NULL);
  P4EST_ASSERT (g != NULL);
  P4EST_ASSERT (v->p4est == g->p4est);

  /* compute value in external userdata array */
  P4EST_ASSERT (g->qarray != NULL);
  *(double *) sc_array_index (g->qarray, g->qcount++) =
    userdata_value (g, v->treeid, v->quad);
}

/* write a VTK file with the mesh and the element data */
static int
userdata_vtk_external (p4est_userdata_global_t *g, const char *filename)
{
  P4EST_ASSERT (g != NULL);
  P4EST_ASSERT (g->p4est != NULL);
  P4EST_ASSERT (g->qarray != NULL);
  P4EST_ASSERT (g->in_external);
  if (g->novtk) {
    /* output disabled */
    return 0;
  }

  /* the data array is already in the required form */
  if (userdata_vtk_general (g, filename)) {
    P4EST_GLOBAL_LERRORF ("ERROR: write VTK file %s\n", filename);
    return -1;
  }

  /* return success */
  return 0;
}

/* callback to evaluate refinement/coarsening criterio */
static void
userdata_criterion_external_volume (p4est_iter_volume_info_t *v,
                                    void *user_data)
{
  /* the global data structure is passed into every iterator */
  p4est_userdata_global_t *g = (p4est_userdata_global_t *) user_data;

  /* check call consistency */
  P4EST_ASSERT (v != NULL);
  P4EST_ASSERT (g != NULL);
  P4EST_ASSERT (v->p4est == g->p4est);

  /* p4est is agnostic to the quadrant user data */
  userdata_quadrant_external_t *qdat =
    (userdata_quadrant_external_t *) v->quad->p.user_data;

  /* refinement criteria for demonstration */
  qdat->refine = ((v->treeid % 3) == 0);

  /* coarsening criteria for demonstration */
  qdat->coarsen = ((v->treeid % 3) == 1);

  /* refine also if value is near the extremes */
  P4EST_ASSERT (g->qarray != NULL);
  if (fabs (*(double *) sc_array_index (g->qarray, g->qcount++)) > .8) {
    qdat->refine = 1;
  }

  /* at this point, refinement and coarsening may both be requested */
}

/* callback to tell p4est which quadrants shall be refined */
static int
userdata_refine_external (p4est_t *p4est, p4est_topidx_t which_tree,
                          p4est_quadrant_t *quadrant)
{
  /* the global data structure is stashed into the forest's user pointer */
  p4est_userdata_global_t *g =
    (p4est_userdata_global_t *) p4est->user_pointer;

  /* access the quadrant user data contents */
  userdata_quadrant_external_t *qdat =
    (userdata_quadrant_external_t *) quadrant->p.user_data;
  P4EST_ASSERT (qdat != NULL);

  /* access refinement criterion */
  if (!qdat->refine || quadrant->level >= g->maxlevel) {
    return 0;
  }
  return 1;
}

/* callback to tell p4est which quadrants shall be coarsened */
static int
userdata_coarsen_external (p4est_t *p4est, p4est_topidx_t which_tree,
                           p4est_quadrant_t *quadrant[])
{
  int                 i;

  /* with external user data, we do not use p4est_coarsen_ext */
  P4EST_ASSERT (quadrant[1] != NULL);
  for (i = 0; i < P4EST_CHILDREN; ++i) {
    /* access the quadrant user data contents */
    userdata_quadrant_external_t *qdat =
      (userdata_quadrant_external_t *) quadrant[i]->p.user_data;
    P4EST_ASSERT (qdat != NULL);

    /* access coarsening criterion */
    if (!qdat->coarsen || qdat->refine) {
      return 0;
    }
  }
  return 1;
}

/* interpolate/project data after refinement/ */
static void
userdata_project_external (p4est_userdata_global_t *g)
{
  int                 i;
  double              value;
  size_t              pdindex, ndindex;
  size_t              pdbound;
#ifdef P4EST_ENABLE_DEBUG
  size_t              ndbound;
#endif
  sc_array_t         *parray, *narray;
  p4est_topidx_t      tt;
  p4est_tree_t       *ptree, *ntree;
  p4est_quadrant_t   *pquad, *nquad;

  /* check invariants */
  P4EST_ASSERT (g != NULL);
  P4EST_ASSERT (g->p4est != NULL);
  P4EST_ASSERT (g->n4est != NULL);
  P4EST_ASSERT (g->qarray != NULL);
  P4EST_ASSERT (g->in_external);

  /* allocate a second data array and deallocate the first afterwards */
  pdindex = ndindex = 0;
  parray = g->qarray;
  P4EST_ASSERT (parray->elem_count == (size_t) g->p4est->local_num_quadrants);
  narray = sc_array_new_count (sizeof (double),
                               (size_t) g->n4est->local_num_quadrants);
  P4EST_ASSERT (g->p4est->first_local_tree == g->n4est->first_local_tree);
  P4EST_ASSERT (g->p4est->last_local_tree == g->n4est->last_local_tree);
  for (tt = g->n4est->first_local_tree; tt <= g->n4est->last_local_tree; ++tt) {
    /* simultaneous loop over old and new mesh with same partition */
    ptree = p4est_tree_array_index (g->p4est->trees, tt);
    pquad = p4est_quadrant_array_index (&ptree->quadrants, 0);
    pdbound = (size_t) ptree->quadrants_offset + ptree->quadrants.elem_count;
    ntree = p4est_tree_array_index (g->n4est->trees, tt);
    nquad = p4est_quadrant_array_index (&ntree->quadrants, 0);
#ifdef P4EST_ENABLE_DEBUG
    ndbound = (size_t) ntree->quadrants_offset + ntree->quadrants.elem_count;
#endif
    for (;;) {
      /* simultaneous loop over old and new quadrants in this tree */
      P4EST_ASSERT (abs (pquad->level - nquad->level) <= 1);
      if (pquad->level == nquad->level) {
        /* old and new quadrant are same size: copy value */
        *(double *) sc_array_index (narray, ndindex++) =
          *(double *) sc_array_index (parray, pdindex++);
        pquad++;
        nquad++;
      }
      else if (pquad->level + 1 == nquad->level) {
        /* new quadrants are refined from the old one: copy value */
        value = *(double *) sc_array_index (parray, pdindex++);
        for (i = 0; i < P4EST_CHILDREN; ++i) {
          /* all new quadrants must be siblings of a family */
          *(double *) sc_array_index (narray, ndindex++) = value;
          P4EST_ASSERT (p4est_quadrant_child_id (nquad) == i);
          P4EST_ASSERT (nquad->level == pquad->level + 1);
          ++nquad;
        }
        ++pquad;
      }
      else if (nquad->level + 1 == pquad->level) {
        /* new quadrant is coarsened from the old ones */
        value = 0.;
        for (i = 0; i < P4EST_CHILDREN; ++i) {
          /* all old quadrants must be siblings of a family */
          value += *(double *) sc_array_index (parray, pdindex++);
          P4EST_ASSERT (p4est_quadrant_child_id (pquad) == i);
          P4EST_ASSERT (pquad->level == nquad->level + 1);
          ++pquad;
        }
        /* for demonstration just compute the average value */
        *(double *) sc_array_index (narray, ndindex++) =
          value * (1. / P4EST_CHILDREN);
        ++nquad;
      }
      else {
        /* no other situation can arise */
        SC_ABORT_NOT_REACHED ();
      }

      /* break quadrant loop if the local tree is traversed */
      if (pdindex == pdbound) {
        P4EST_ASSERT (ndindex == ndbound);
        break;
      }
      P4EST_ASSERT (ndindex < ndbound);
    }

    /* this tree is now processed */
    P4EST_ASSERT (pquad - p4est_quadrant_array_index (&ptree->quadrants, 0) ==
                  (ptrdiff_t) ptree->quadrants.elem_count);
    P4EST_ASSERT (nquad - p4est_quadrant_array_index (&ntree->quadrants, 0) ==
                  (ptrdiff_t) ntree->quadrants.elem_count);
  }

  /* by construction we have traversed both meshes simultaneously */
  P4EST_ASSERT (pdindex == (size_t) g->p4est->local_num_quadrants);
  P4EST_ASSERT (ndindex == (size_t) g->n4est->local_num_quadrants);

  /* replace the old array with the new one */
  sc_array_destroy (g->qarray);
  g->qarray = narray;
}

/* adapt the mesh according to refinement criteria and project application data */
static void
userdata_adapt_external (p4est_userdata_global_t *g)
{
  P4EST_ASSERT (g != NULL);
  P4EST_ASSERT (g->in_external);

  /* execute refinement on a copy of the forest */
  P4EST_ASSERT (g->n4est == NULL);
  g->n4est = p4est_copy (g->p4est, 1);

  /* adapt the new forest non-recursively */
  p4est_refine (g->n4est, 0,
                userdata_refine_external, userdata_init_external);
  p4est_coarsen (g->n4est, 0,
                 userdata_coarsen_external, userdata_init_external);
  p4est_balance (g->n4est, P4EST_CONNECT_FULL, userdata_init_external);

  /* interpolate/project values from old to adapted forest */
  userdata_project_external (g);
  p4est_destroy (g->p4est);
  g->p4est = g->n4est;
  g->n4est = NULL;
}

/* execute partitioning with associated data transfer */
static void
userdata_partition_external (p4est_userdata_global_t *g)
{
  P4EST_ASSERT (g != NULL);
  P4EST_ASSERT (g->qarray != NULL);
  P4EST_ASSERT (g->n4est == NULL);
  P4EST_ASSERT (g->in_external);

  /* partition moves the quadrant user data around */
  g->n4est = p4est_copy (g->p4est, 1);
  if (p4est_partition_ext (g->n4est, 1, NULL) == 0) {
    /* the partition is unchanged, do nothing */
    P4EST_ASSERT (g->p4est->local_num_quadrants ==
                  g->n4est->local_num_quadrants);
    p4est_destroy (g->n4est);
  }
  else {
    sc_array_t         *parray, *narray;

    /* transfer quadrant data from old to new partition */
    parray = g->qarray;
    P4EST_ASSERT (parray->elem_count ==
                  (size_t) g->p4est->local_num_quadrants);
    narray = sc_array_new_count (sizeof (double),
                                 (size_t) g->n4est->local_num_quadrants);
    p4est_transfer_fixed (g->n4est->global_first_quadrant,
                          g->p4est->global_first_quadrant,
                          /* the tag can really be anything */
                          g->p4est->mpicomm, P4EST_COMM_TAG_LAST + 0,
                          /* don't assert invalid index on empty array */
                          sc_array_index_null (narray, 0),
                          sc_array_index_null (parray, 0), parray->elem_size);

    /* free old forest and application data */
    sc_array_destroy (g->qarray);
    g->qarray = narray;
    p4est_destroy (g->p4est);
    g->p4est = g->n4est;
  }

  /* maintain invariant */
  g->n4est = NULL;
}

/* provide function for consistent deallocation */
static int
userdata_run_external_return (int retval, p4est_userdata_global_t *g)
{
  P4EST_ASSERT (g != NULL);

  /* clean up what has been allocated in the same function */
  if (g->qarray != NULL) {
    /* destroy application userdata */
    sc_array_destroy_null (&g->qarray);
  }
  if (g->n4est != NULL) {
    /* destroy forest */
    p4est_destroy (g->n4est);
    g->n4est = NULL;
  }
  if (g->p4est != NULL) {
    /* destroy forest */
    p4est_destroy (g->p4est);
    g->p4est = NULL;
  }
  return retval;
}

/* core demo with quadrant data stored external to p4est */
static int
userdata_run_external (p4est_userdata_global_t *g)
{
  P4EST_ASSERT (g != NULL);
  P4EST_ASSERT (g->p4est == NULL);
  P4EST_ASSERT (g->qarray == NULL);
  P4EST_ASSERT (g->in_external);

  /* this is the demo for userdata allocated by the application developer */
  P4EST_GLOBAL_PRODUCTION (P4EST_STRING
                           "_userdata: application data EXTERNAL\n");

  /* create initial forest and populate metadata by callback */
  P4EST_ASSERT (g->qcount == 0);
  g->p4est = p4est_new_ext
    (g->mpicomm, g->conn, 0, SC_MAX (g->maxlevel - 1, 0), 1,
     sizeof (userdata_quadrant_external_t), userdata_init_external, g);
  P4EST_ASSERT (g->qcount == g->p4est->local_num_quadrants);
  g->qcount = 0;

  /* We like the invariant that after partitioning, partition-independent
     coarsening is always possible since every family of siblings is placed
     on a single process; this is not guaranteed by p4est_new_ext. */
  p4est_partition (g->p4est, 1, NULL);

  /* populate quadrant data by volume iterator */
  P4EST_ASSERT (g->qarray == NULL);
  g->qarray = sc_array_new_count (sizeof (double),
                                  (size_t) g->p4est->local_num_quadrants);
  userdata_iterate_volume (g, userdata_init_external_volume);
  g->qcount = 0;

  /* write VTK files to visualize geometry and data */
  if (userdata_vtk_external (g, P4EST_STRING "_userdata_external_new")) {
    P4EST_GLOBAL_LERROR ("ERROR: write VTK output after forest_new\n");
    return userdata_run_external_return (-1, g);
  }

  /* evaluate refinement criteria */
  userdata_iterate_volume (g, userdata_criterion_external_volume);

  /* refine, coarsen, and balance the mesh, then project data once */
  userdata_adapt_external (g);

  /* write VTK files to visualize geometry and data */
  if (userdata_vtk_external (g, P4EST_STRING "_userdata_external_adapt")) {
    P4EST_GLOBAL_LERROR ("ERROR: write VTK output after forest_adapt\n");
    return userdata_run_external_return (-1, g);
  }

  /* repartition the mesh and the application data */
  userdata_partition_external (g);

  /* write VTK files to visualize geometry and data */
  if (userdata_vtk_external (g, P4EST_STRING "_userdata_external_partition")) {
    P4EST_GLOBAL_LERROR ("ERROR: write VTK output after forest_partition\n");
    return userdata_run_external_return (-1, g);
  }

  /* return memory neutral */
  return userdata_run_external_return (0, g);
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
  P4EST_ASSERT (!g->in_internal);
  P4EST_ASSERT (!g->in_external);
  erres = 0;

  /* run the example once with p4est-allocated application data */
  g->in_internal = 1;
  if (!erres && !g->noint && (erres = userdata_run_internal (g))) {
    P4EST_GLOBAL_LERROR ("ERROR: run with internal data\n");
  }
  g->in_internal = 0;

  /* run the example another time with user-allocated application data */
  g->in_external = 1;
  if (!erres && !g->noext && (erres = userdata_run_external (g))) {
    P4EST_GLOBAL_LERROR ("ERROR: run with external data\n");
  }
  g->in_external = 0;

  /* return error status */
  P4EST_ASSERT (!g->in_internal);
  P4EST_ASSERT (!g->in_external);
  return erres;
}
