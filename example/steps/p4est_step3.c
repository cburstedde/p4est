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

/* p4est has two separate interfaces for 2D and 3D, p4est*.h and p8est*.h.
 * Most API functions are available for both dimensions.  The header file
 * p4est_to_p8est.h #define's the 2D names to the 3D names such that most code
 * only needs to be written once.  In this example, we rely on this. */
#ifndef P4_TO_P8
#include <p4est_vtk.h>
#include <p4est_bits.h>
#include <p4est_extended.h>
#include <p4est_iterate.h>
#else
#include <p8est_vtk.h>
#include <p8est_bits.h>
#include <p8est_extended.h>
#include <p8est_iterate.h>
#endif

/* In this example we store data with each quadrant/octant. */

typedef struct step3_data
{
  double              u;
  double              du[P4EST_DIM];
}
step3_data_t;

typedef struct step3_ctx
{
  double              center[P4EST_DIM];        /* center of the initial condition Gaussian bump */
  double              bump_width;       /* width of the initial condition Gaussian bump */
  double              max_err;  /* maximum global interpolation error */
  double              v[P4EST_DIM];     /* advection direction */
}
step3_ctx_t;

/* The initial condition: a Gaussian bump.
 * Returns the value of the initial condition at x, and optionally the
 * gradient at that point as well */
static double
initial_condition (double x[], double du[], step3_ctx_t * ctx)
{
  int                 i;
  double             *c = ctx->center;
  double              bump_width = ctx->bump_width;
  double              r2, d[P4EST_DIM];
  double              arg, retval;

  r2 = 0.;
  for (i = 0; i < P4EST_DIM; i++) {
    d[i] = x[i] - c[i];
    r2 += d[i] * d[i];
  }

  arg = -(1. / 2.) * r2 / bump_width / bump_width;
  retval = exp (arg);

  if (du) {
    for (i = 0; i < P4EST_DIM; i++) {
      du[i] = -(1. / bump_width / bump_width) * d[i] * retval;
    }
  }

  return retval;
}

/* This function converts the quadrant data to the coordinates of the center
 * of a quadrant */
static void
quad_get_midpoint (p4est_t * p4est, p4est_topidx_t which_tree,
                   p4est_quadrant_t * q, double xyz[3])
{
  p4est_qcoord_t      half_length = P4EST_QUADRANT_LEN (q->level) / 2;

  p4est_qcoord_to_vertex (p4est->connectivity, which_tree,
                          q->x + half_length, q->y + half_length,
#ifdef P4_TO_P8
                          q->z + half_length,
#endif
                          xyz);
}

/* Initialize the initial condition value and derivative for each newly
 * created quadrant */
static void
init_initial_condition (p4est_t * p4est, p4est_topidx_t which_tree,
                        p4est_quadrant_t * q)
{
  step3_ctx_t        *ctx = (step3_ctx_t *) p4est->user_pointer;
  step3_data_t       *data = (step3_data_t *) q->p.user_data;
  double              midpoint[3];

  quad_get_midpoint (p4est, which_tree, q, midpoint);
  data->u = initial_condition (midpoint, data->du, ctx);
}

/* Estimate the square of the L2 error on the quadrant by assuming that the
 * function is linear and that the value and derivative are correct at the
 * midpoint. */
static double
error_sqr_estimate (p4est_quadrant_t * q)
{
  step3_data_t       *data = (step3_data_t *) q->p.user_data;
  int                 i;
  double              diff2;
  double             *du = data->du;
  double              h =
    (double) P4EST_QUADRANT_LEN (q->level) / (double) P4EST_ROOT_LEN;
  double              vol;

#ifdef P4_TO_P8
  vol = h * h * h;
#else
  vol = h * h;
#endif

  diff2 = 0.;
  /* use the approximate derivative to esimate the L2 error */
  for (i = 0; i < P4EST_DIM; i++) {
    diff2 += du[i] * du[i] * (1. / 12.) * h * h * vol;
  }

  return diff2;
}

/* Refine by error estimate above.  Given the maximum global error, we
 * enforce that each quadrant's portion of the error must not exceed is
 * fraction of the total volume of the domain (which is 1). */
static int
refine_err_estimate (p4est_t * p4est, p4est_topidx_t which_tree,
                     p4est_quadrant_t * q)
{
  step3_ctx_t        *ctx = (step3_ctx_t *) p4est->user_pointer;
  double              global_err = ctx->max_err;
  double              global_err2 = global_err * global_err;
  double              h =
    (double) P4EST_QUADRANT_LEN (q->level) / (double) P4EST_ROOT_LEN;
  double              vol, err2;

  /* the quadrant's volume is also its volume fraction */
#ifdef P4_TO_P8
  vol = h * h * h;
#else
  vol = h * h;
#endif

  err2 = error_sqr_estimate (q);
  if (err2 > global_err2 * vol) {
    return 1;
  }
  else {
    return 0;
  }
}

/* Coarsen by error estimate above.  Given the maximum global error, we
 * enforce that each quadrant's portion of the error must not exceed is
 * fraction of the total volume of the domain (which is 1). */
static int
coarsen_err_estimate_initial_condition (p4est_t * p4est,
                                        p4est_topidx_t which_tree,
                                        p4est_quadrant_t * children[])
{
  p4est_quadrant_t    parent;
  step3_ctx_t        *ctx = (step3_ctx_t *) p4est->user_pointer;
  double              global_err = ctx->max_err;
  double              global_err2 = global_err * global_err;
  double              h;
  step3_data_t        parentdata;
  double              parentmidpoint[P4EST_DIM];
  double              vol, err2;

  /* get the parent of the first child (the parent of all children) */
  p4est_quadrant_parent (children[0], &parent);
  quad_get_midpoint (p4est, which_tree, &parent, parentmidpoint);
  parentdata.u = initial_condition (parentmidpoint, parentdata.du, ctx);
  h = (double) P4EST_QUADRANT_LEN (parent.level) / (double) P4EST_ROOT_LEN;
  /* the quadrant's volume is also its volume fraction */
#ifdef P4_TO_P8
  vol = h * h * h;
#else
  vol = h * h;
#endif
  parent.p.user_data = (void *) (&parentdata);

  err2 = error_sqr_estimate (&parent);
  if (err2 < global_err2 * vol) {
    return 1;
  }
  else {
    return 0;
  }
}

/* Callback function for interpolating the solution from quadrant midpoints to
 * corners.  This callback is passed to p4est_iterate, which calls it for
 * every quadrant. The info object passes information about the current
 * quadrant to the callback: it is populated by p4est_iterate.  (See
 * p4est_iterate.h) */
static void
interpolate_solution (p4est_iter_volume_info_t * info, void *user_data)
{
  double             *u_interp = (double *) user_data;  /* we passed the array of values to fill as the user_data in the call to p4est_iterate */
  p4est_t            *p4est = info->p4est;
  p4est_quadrant_t   *q = info->quad;
  p4est_topidx_t      which_tree = info->treeid;
  p4est_locidx_t      local_id = info->quadid;  /* this is the index of q *within its tree's numbering*.  We want to convert it its index for all the quadrants on this process, which we do below */
  p4est_tree_t       *tree;
  step3_data_t       *data = (step3_data_t *) q->p.user_data;
  double              h;
  p4est_locidx_t      arrayoffset;
  double              this_u;
  int                 i, j;

  tree = p4est_tree_array_index (p4est->trees, which_tree);
  local_id += tree->quadrants_offset;   /* now the id is relative to the MPI process */
  arrayoffset = P4EST_CHILDREN * local_id;      /* each local quadrant has 2^d (P4EST_CHILDREN) values in u_interp */
  h = (double) P4EST_QUADRANT_LEN (q->level) / (double) P4EST_ROOT_LEN;

  for (i = 0; i < P4EST_CHILDREN; i++) {
    this_u = data->u;
    /* loop over the derivative components and linearly interpolate from the
     * midpoint to the corners */
    for (j = 0; j < P4EST_DIM; j++) {
      /* In order to know whether the direction from the midpoint to the corner is
       * negative or positive, we take advantage of the fact that the corners
       * are in z-order.  If i is an odd number, it is on the +x side; if it
       * is even, it is on the -x side.  If (i / 2) is an odd number, it is on
       * the +y side, etc. */
      this_u += (h / 2) * data->du[j] * ((i & (1 << j)) ? 1. : -1.);
    }
    u_interp[arrayoffset + i] = this_u;
  }

}

static void
step3_write_solution (p4est_t * p4est, int timestep)
{
  char                filename[BUFSIZ] = { '\0' };
  double             *u_interp;
  p4est_locidx_t      numquads;

  snprintf (filename, 16, P4EST_STRING "_step3_%4.4d", timestep);

  numquads = p4est->local_num_quadrants;

  /* create a vector with one value for the corner of every local quadrant
   * (the number of children is always the same as the number of corners) */
  u_interp = P4EST_ALLOC (double, numquads * P4EST_CHILDREN);

  /* Use the iterator to visit every cell and fill in the solution values.
   * Using the iterator is not absolutely necessary in this case: we could
   * also loop over every tree (there is only one tree in this case) and loop
   * over every quadrant within every tree, but we are trying to demonstrate
   * the usage of p4est_iterate in this example */
  p4est_iterate (p4est, NULL,   /* we don't need any ghost quadrants for this loop */
                 (void *) u_interp,     /* pass in u_interp so that we can fill it */
                 interpolate_solution,  /* callback function that interpolate from the cell center to the cell corners, defined above */
                 NULL,          /* there is no callback for the faces between quadrants */
#ifdef P4_TO_P8
                 NULL,          /* there is no callback for the edges between quadrants */
#endif
                 NULL);         /* there is no callback for the corners between quadrants */

  p4est_vtk_write_all (p4est, NULL,     /* we do not need to transform from the vertex space into physical space, so we do not need a p4est_geometry_t * pointer */
                       0.99,    /* draw each quadrant at almost full scale */
                       0,       /* do not write the tree id's of each quadrant (there is only one tree in this example) */
                       1,       /* do write the refinement level of each quadrant */
                       1,       /* do write the mpi process id of each quadrant */
                       0,       /* do not wrap the mpi rank (if this were > 0, the modulus of the rank relative to this number would be written instead of the rank) */
                       1,       /* write one scalar field: the solution value */
                       0,       /* write no vector fields */
                       filename, "solution", u_interp);

  P4EST_FREE (u_interp);
}

int
main (int argc, char **argv)
{
  int                 mpiret;
  int                 recursive, partforcoarsen, balance;
  sc_MPI_Comm         mpicomm;
  p4est_t            *p4est;
  p4est_connectivity_t *conn;
  step3_ctx_t         ctx;

  /* Initialize MPI; see sc_mpi.h.
   * If configure --enable-mpi is given these are true MPI calls.
   * Else these are dummy functions that simulate a single-processor run. */
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  mpicomm = sc_MPI_COMM_WORLD;

  /* These functions are optional.  If called they store the MPI rank as a
   * static variable so subsequent global p4est log messages are only issued
   * from processor zero.  Here we turn off most of the logging; see sc.h. */
  sc_init (mpicomm, 1, 1, NULL, SC_LP_ESSENTIAL);
  p4est_init (NULL, SC_LP_PRODUCTION);
  P4EST_GLOBAL_PRODUCTIONF
    ("This is the p4est %dD demo example/steps/%s_step3\n",
     P4EST_DIM, P4EST_STRING);

  ctx.bump_width = 0.1;
  ctx.max_err = 1.e-2;
  ctx.center[0] = 0.5;
  ctx.center[1] = 0.5;
#ifdef P4_TO_P8
  ctx.center[2] = 0.5;
#endif

  /* Create a forest that consists of just one periodic quadtree/octree. */
#ifndef P4_TO_P8
  conn = p4est_connectivity_new_periodic ();
#else
  conn = p8est_connectivity_new_periodic ();
#endif

  /* Create a forest that has 16 cells in each direction */
  p4est = p4est_new_ext (mpicomm, conn, 0,      /* minimum quadrants per mpi process */
                         4,     /* minimum level of refinement */
                         1,     /* fill uniform */
                         sizeof (step3_data_t), /* data size */
                         init_initial_condition,        /* data initializiation */
                         (void *) (&ctx));      /* user data */

  recursive = 1;
  p4est_refine (p4est, recursive, refine_err_estimate,
                init_initial_condition);
  p4est_coarsen (p4est, recursive, coarsen_err_estimate_initial_condition,
                 init_initial_condition);

  /* Partition: The quadrants are redistributed for equal element count.  The
   * partition can optionally be modified such that a family of octants, which
   * are possibly ready for coarsening, are never split between processors. */
  partforcoarsen = 1;
  p4est_partition (p4est, partforcoarsen, NULL);

  /* If we call the 2:1 balance we ensure that neighbors do not differ in size
   * by more than a factor of 2.  This can optionally include diagonal
   * neighbors across edges or corners as well; see p4est.h. */
  balance = 1;
  if (balance) {
    p4est_balance (p4est, P4EST_CONNECT_FACE, init_initial_condition);
    p4est_partition (p4est, partforcoarsen, NULL);
  }

  /* Write the forest and the solution to disk for visualization, one file per
   * processor. */
  step3_write_solution (p4est, 0);

  /* Destroy the p4est and the connectivity structure. */
  p4est_destroy (p4est);
  p4est_connectivity_destroy (conn);

  /* Verify that allocations internal to p4est and sc do not leak memory.
   * This should be called if sc_init () has been called earlier. */
  sc_finalize ();

  /* This is standard MPI programs.  Without --enable-mpi, this is a dummy. */
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;
}
