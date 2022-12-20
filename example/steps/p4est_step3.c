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

/** \file p4est_step3.c
 *
 * This 2D example program uses p4est to solve a simple advection problem.  It
 * is numerically very simple, and intended to demonstrate several methods of
 * interacting with the p4est data after it has been refined and partitioned.
 * It demonstrates the construction of ghost layers (see p4est_ghost_t in
 * p4est_ghost.h) and communication of ghost-layer data, and it demonstrates
 * interacting with the quadrants and quadrant boundaries through the
 * p4est_iterate() routine (see p4est_iterate.h).
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
#include <p4est_io.h>
#include <p4est_communication.h>
#else
#include <p8est_vtk.h>
#include <p8est_bits.h>
#include <p8est_extended.h>
#include <p8est_iterate.h>
#include <p8est_io.h>
#include <p8est_communication.h>
#endif
#include <sc_options.h>

#ifndef P4_TO_P8
#define P4EST_DATA_FILE_EXT "p4d" /**< file extension of p4est data files */
#else
#define P4EST_DATA_FILE_EXT               P8EST_DATA_FILE_EXT
#define P8EST_DATA_FILE_EXT "p8d" /**< file extension of p8est data files */
#endif

#define STEP3_BLOCK_SIZE (sizeof (step3_ctx_t)) /**< number of bytes of
                                                      the simulation context */
#define STEP3_ENDIAN_CHECK 0x30415062   /* check endian; bPA0 in ASCII */

/** We had 1. / 0. here to create a NaN but that is not portable. */
static const double step3_invalid = -1.;

/** A Boolean to decide if checkpoint files are written to disk. */
static int          step3_checkpoint = 0;

/* In this example we store data with each quadrant/octant. */

/** Per-quadrant data for this example.
 *
 * In this problem, we keep track of the state variable u, its
 * derivatives in space, and in time.
 */
typedef struct step3_data
{
  double              u;             /**< the state variable */
  double              du[P4EST_DIM]; /**< the spatial derivatives */
  double              dudt;          /**< the time derivative */
}
step3_data_t;

/** The example parameters.
 *
 * This describes the advection problem and time-stepping used in this
 * example.
 */
typedef struct step3_ctx
{
  double              center[P4EST_DIM];  /**< coordinates of the center of
                                               the initial condition Gaussian
                                               bump */
  double              bump_width;         /**< width of the initial condition
                                               Gaussian bump */
  double              max_err;            /**< maximum allowed global
                                               interpolation error */
  double              v[P4EST_DIM];       /**< the advection velocity */
  int                 refine_period;      /**< the number of time steps
                                               between mesh refinement */
  int                 repartition_period; /**< the number of time steps
                                               between repartitioning */
  int                 write_period;       /**< the number of time steps
                                               between writing vtk files */
  double              current_time;       /**< the current time */
  int                 time_step;          /**< current time step
                                               counted from the first start. */
}
step3_ctx_t;

/** Compute the value and derivatives of the initial condition.
 *
 * \param [in]  x   the coordinates
 * \param [out] du  the derivative at \a x
 * \param [in]  ctx the example parameters
 *
 * \return the initial condition at \a x
 */
static double
step3_initial_condition (double x[], double du[], step3_ctx_t * ctx)
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

/** Get the coordinates of the midpoint of a quadrant.
 *
 * \param [in]  p4est      the forest
 * \param [in]  which_tree the tree in the forest containing \a q
 * \param [in]  q          the quadrant
 * \param [out] xyz        the coordinates of the midpoint of \a q
 */
static void
step3_get_midpoint (p4est_t * p4est, p4est_topidx_t which_tree,
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

/** Initialize the initial condition data of a quadrant.
 *
 * This function matches the p4est_init_t prototype that is used by
 * p4est_new(), p4est_refine(), p4est_coarsen(), and p4est_balance().
 *
 * \param [in] p4est          the forest
 * \param [in] which_tree     the tree in the forest containing \a q
 * \param [in,out] q          the quadrant whose data gets initialized
 */
static void
step3_init_initial_condition (p4est_t * p4est, p4est_topidx_t which_tree,
                              p4est_quadrant_t * q)
{
  /* the data associated with a forest is accessible by user_pointer */
  step3_ctx_t        *ctx = (step3_ctx_t *) p4est->user_pointer;
  /* the data associated with a quadrant is accessible by p.user_data */
  step3_data_t       *data = (step3_data_t *) q->p.user_data;
  double              midpoint[3];

  step3_get_midpoint (p4est, which_tree, q, midpoint);
  /* initialize the data */
  data->u = step3_initial_condition (midpoint, data->du, ctx);
}

/** Estimate the square of the approximation error on a quadrant.
 *
 * We compute our estimate by integrating the difference of a constant
 * approximation at the midpoint and a linear approximation that interpolates
 * at the midpoint.
 *
 * \param [in] q a quadrant
 *
 * \return the square of the error estimate for the state variables contained
 * in \a q's data.
 */
static double
step3_error_sqr_estimate (p4est_quadrant_t * q)
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
  /* use the approximate derivative to estimate the L2 error */
  for (i = 0; i < P4EST_DIM; i++) {
    diff2 += du[i] * du[i] * (1. / 12.) * h * h * vol;
  }

  return diff2;
}

/** Refine by the L2 error estimate.
 *
 * Given the maximum global error, we enforce that each quadrant's portion of
 * the error must not exceed is fraction of the total volume of the domain
 * (which is 1).
 *
 * This function matches the p4est_refine_t prototype that is used by
 * p4est_refine() and p4est_refine_ext().
 *
 * \param [in] p4est          the forest
 * \param [in] which_tree     the tree in the forest containing \a q
 * \param [in] q              the quadrant
 *
 * \return 1 if \a q should be refined, 0 otherwise.
 */
static int
step3_refine_err_estimate (p4est_t * p4est, p4est_topidx_t which_tree,
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

  err2 = step3_error_sqr_estimate (q);
  if (err2 > global_err2 * vol) {
    return 1;
  }
  else {
    return 0;
  }
}

/** Coarsen by the L2 error estimate of the initial condition.
 *
 * Given the maximum global error, we enforce that each quadrant's portion of
 * the error must not exceed is fraction of the total volume of the domain
 * (which is 1).
 *
 * \param [in] p4est          the forest
 * \param [in] which_tree     the tree in the forest containing \a children
 * \param [in] children       a family of quadrants
 *
 * \return 1 if \a children should be coarsened, 0 otherwise.
 */
static int
step3_coarsen_initial_condition (p4est_t * p4est,
                                 p4est_topidx_t which_tree,
                                 p4est_quadrant_t * children[])
{
  p4est_quadrant_t    parent;
  step3_ctx_t        *ctx = (step3_ctx_t *) p4est->user_pointer;
  double              global_err = ctx->max_err;
  double              global_err2 = global_err * global_err;
  double              h;
  step3_data_t        parentdata;
  double              parentmidpoint[3];
  double              vol, err2;

  /* get the parent of the first child (the parent of all children) */
  p4est_quadrant_parent (children[0], &parent);
  step3_get_midpoint (p4est, which_tree, &parent, parentmidpoint);
  parentdata.u = step3_initial_condition (parentmidpoint, parentdata.du, ctx);
  h = (double) P4EST_QUADRANT_LEN (parent.level) / (double) P4EST_ROOT_LEN;
  /* the quadrant's volume is also its volume fraction */
#ifdef P4_TO_P8
  vol = h * h * h;
#else
  vol = h * h;
#endif
  parent.p.user_data = (void *) (&parentdata);

  err2 = step3_error_sqr_estimate (&parent);
  if (err2 < global_err2 * vol) {
    return 1;
  }
  else {
    return 0;
  }
}

/** Coarsen by the L2 error estimate of the current state approximation.
 *
 * Given the maximum global error, we enforce that each quadrant's portion of
 * the error must not exceed its fraction of the total volume of the domain
 * (which is 1).
 *
 * This function matches the p4est_coarsen_t prototype that is used by
 * p4est_coarsen() and p4est_coarsen_ext().
 *
 * \param [in] p4est          the forest
 * \param [in] which_tree     the tree in the forest containing \a children
 * \param [in] children       a family of quadrants
 *
 * \return 1 if \a children should be coarsened, 0 otherwise.
 */
static int
step3_coarsen_err_estimate (p4est_t * p4est,
                            p4est_topidx_t which_tree,
                            p4est_quadrant_t * children[])
{
  step3_ctx_t        *ctx = (step3_ctx_t *) p4est->user_pointer;
  double              global_err = ctx->max_err;
  double              global_err2 = global_err * global_err;
  double              h;
  step3_data_t       *data;
  double              vol, err2, childerr2;
  double              parentu;
  double              diff;
  int                 i;

  h =
    (double) P4EST_QUADRANT_LEN (children[0]->level) /
    (double) P4EST_ROOT_LEN;
  /* the quadrant's volume is also its volume fraction */
#ifdef P4_TO_P8
  vol = h * h * h;
#else
  vol = h * h;
#endif

  /* compute the average */
  parentu = 0.;
  for (i = 0; i < P4EST_CHILDREN; i++) {
    data = (step3_data_t *) children[i]->p.user_data;
    parentu += data->u / P4EST_CHILDREN;
  }

  err2 = 0.;
  for (i = 0; i < P4EST_CHILDREN; i++) {
    childerr2 = step3_error_sqr_estimate (children[i]);

    if (childerr2 > global_err2 * vol) {
      return 0;
    }
    err2 += step3_error_sqr_estimate (children[i]);
    diff = (parentu - data->u) * (parentu - data->u);
    err2 += diff * vol;
  }
  if (err2 < global_err2 * (vol * P4EST_CHILDREN)) {
    return 1;
  }
  else {
    return 0;
  }
}

/** Initialize the state variables of incoming quadrants from outgoing
 * quadrants.
 *
 * The functions p4est_refine_ext(), p4est_coarsen_ext(), and
 * p4est_balance_ext() take as an argument a p4est_replace_t callback function,
 * which allows one to setup the quadrant data of incoming quadrants from the
 * data of outgoing quadrants, before the outgoing data is destroyed.  This
 * function matches the p4est_replace_t prototype.
 *
 * In this example, we linearly interpolate the state variable of a quadrant
 * that is refined to its children, and we average the midpoints of children
 * that are being coarsened to the parent.
 *
 * \param [in] p4est          the forest
 * \param [in] which_tree     the tree in the forest containing \a children
 * \param [in] num_outgoing   the number of quadrants that are being replaced:
 *                            either 1 if a quadrant is being refined, or
 *                            P4EST_CHILDREN if a family of children are being
 *                            coarsened.
 * \param [in] outgoing       the outgoing quadrants
 * \param [in] num_incoming   the number of quadrants that are being added:
 *                            either P4EST_CHILDREN if a quadrant is being refined, or
 *                            1 if a family of children are being
 *                            coarsened.
 * \param [in,out] incoming   quadrants whose data are initialized.
 */
static void
step3_replace_quads (p4est_t * p4est, p4est_topidx_t which_tree,
                     int num_outgoing,
                     p4est_quadrant_t * outgoing[],
                     int num_incoming, p4est_quadrant_t * incoming[])
{
  step3_data_t       *parent_data, *child_data;
  int                 i, j;
  double              h;
  double              du_old, du_est;

  if (num_outgoing > 1) {
    /* this is coarsening */
    parent_data = (step3_data_t *) incoming[0]->p.user_data;
    parent_data->u = 0.;
    for (j = 0; j < P4EST_DIM; j++) {
      parent_data->du[j] = step3_invalid;

    }
    for (i = 0; i < P4EST_CHILDREN; i++) {
      child_data = (step3_data_t *) outgoing[i]->p.user_data;
      parent_data->u += child_data->u / P4EST_CHILDREN;
      for (j = 0; j < P4EST_DIM; j++) {
        du_old = parent_data->du[j];
        du_est = child_data->du[j];

        if (du_old == du_old) {
          if (du_est * du_old >= 0.) {
            if (fabs (du_est) < fabs (du_old)) {
              parent_data->du[j] = du_est;
            }
          }
          else {
            parent_data->du[j] = 0.;
          }
        }
        else {
          parent_data->du[j] = du_est;
        }
      }
    }
  }
  else {
    /* this is refinement */
    parent_data = (step3_data_t *) outgoing[0]->p.user_data;
    h =
      (double) P4EST_QUADRANT_LEN (outgoing[0]->level) /
      (double) P4EST_ROOT_LEN;

    for (i = 0; i < P4EST_CHILDREN; i++) {
      child_data = (step3_data_t *) incoming[i]->p.user_data;
      child_data->u = parent_data->u;
      for (j = 0; j < P4EST_DIM; j++) {
        child_data->du[j] = parent_data->du[j];
        child_data->u +=
          (h / 4.) * parent_data->du[j] * ((i & (1 << j)) ? 1. : -1);
      }
    }
  }
}

/** Callback function for interpolating the solution from quadrant midpoints to
 * corners.
 *
 * The function p4est_iterate() takes as an argument a p4est_iter_volume_t
 * callback function, which it executes at every local quadrant (see
 * p4est_iterate.h).  This function matches the p4est_iter_volume_t prototype.
 *
 * In this example, we use the callback function to interpolate the state
 * variable to the corners, and write those corners into an array so that they
 * can be written out.
 *
 * \param [in] info          the information about this quadrant that has been
 *                           populated by p4est_iterate()
 * \param [in,out] user_data the user_data that was given as an argument to
 *                           p4est_iterate: in this case, it points to the
 *                           array of corner values that we want to write.
 *                           The values for the corner of the quadrant
 *                           described by \a info are written during the
 *                           execution of the callback.
 */
static void
step3_interpolate_solution (p4est_iter_volume_info_t * info, void *user_data)
{
  sc_array_t         *u_interp = (sc_array_t *) user_data;      /* we passed the array of values to fill as the user_data in the call to p4est_iterate */
  p4est_t            *p4est = info->p4est;
  p4est_quadrant_t   *q = info->quad;
  p4est_topidx_t      which_tree = info->treeid;
  p4est_locidx_t      local_id = info->quadid;  /* this is the index of q *within its tree's numbering*.  We want to convert it its index for all the quadrants on this process, which we do below */
  p4est_tree_t       *tree;
  step3_data_t       *data = (step3_data_t *) q->p.user_data;
  double              h;
  p4est_locidx_t      arrayoffset;
  double              this_u;
  double             *this_u_ptr;
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
    this_u_ptr = (double *) sc_array_index (u_interp, arrayoffset + i);
    this_u_ptr[0] = this_u;
  }

}

/** Write the state variable to vtk format, one file per process.
 *
 * \param [in] p4est    the forest, whose quadrant data contains the state
 * \param [in] timestep the timestep number, used to name the output files
 */
static void
step3_write_solution (p4est_t * p4est, int timestep)
{
  char                filename[BUFSIZ] = "";
  int                 retval;
  sc_array_t         *u_interp;
  p4est_locidx_t      numquads;
  p4est_vtk_context_t *context;

  snprintf (filename, BUFSIZ, P4EST_STRING "_step3_%04d", timestep);

  numquads = p4est->local_num_quadrants;

  /* create a vector with one value for the corner of every local quadrant
   * (the number of children is always the same as the number of corners) */
  u_interp = sc_array_new_size (sizeof (double), numquads * P4EST_CHILDREN);

  /* Use the iterator to visit every cell and fill in the solution values.
   * Using the iterator is not absolutely necessary in this case: we could
   * also loop over every tree (there is only one tree in this case) and loop
   * over every quadrant within every tree, but we are trying to demonstrate
   * the usage of p4est_iterate in this example */
  p4est_iterate (p4est, NULL,   /* we don't need any ghost quadrants for this loop */
                 (void *) u_interp,     /* pass in u_interp so that we can fill it */
                 step3_interpolate_solution,    /* callback function that interpolates from the cell center to the cell corners, defined above */
                 NULL,          /* there is no callback for the faces between quadrants */
#ifdef P4_TO_P8
                 NULL,          /* there is no callback for the edges between quadrants */
#endif
                 NULL);         /* there is no callback for the corners between quadrants */

  /* create VTK output context and set its parameters */
  context = p4est_vtk_context_new (p4est, filename);
  p4est_vtk_context_set_scale (context, 0.99);  /* quadrant at almost full scale */

  /* begin writing the output files */
  context = p4est_vtk_write_header (context);
  SC_CHECK_ABORT (context != NULL,
                  P4EST_STRING "_vtk: Error writing vtk header");

  /* do not write the tree id's of each quadrant
   * (there is only one tree in this example) */
  context = p4est_vtk_write_cell_dataf (context, 0, 1,  /* do write the refinement level of each quadrant */
                                        1,      /* do write the mpi process id of each quadrant */
                                        0,      /* do not wrap the mpi rank (if this were > 0, the modulus of the rank relative to this number would be written instead of the rank) */
                                        0,      /* there is no custom cell scalar data. */
                                        0,      /* there is no custom cell vector data. */
                                        context);       /* mark the end of the variable cell data. */
  SC_CHECK_ABORT (context != NULL,
                  P4EST_STRING "_vtk: Error writing cell data");

  /* write one scalar field: the solution value */
  context = p4est_vtk_write_point_dataf (context, 1, 0, /* write no vector fields */
                                         "solution", u_interp, context);        /* mark the end of the variable cell data. */
  SC_CHECK_ABORT (context != NULL,
                  P4EST_STRING "_vtk: Error writing cell data");

  retval = p4est_vtk_write_footer (context);
  SC_CHECK_ABORT (!retval, P4EST_STRING "_vtk: Error writing footer");

  sc_array_destroy (u_interp);
}

/** Write a checkpoint file of the current simulation.
 * The file can be loaded using \ref step3_restart to
 * restart the simulation.
 * The checkpoint file is not compiler independent since structures
 * may be padded by the compiler (e. g., sizeof (step3_data_t)).
 * This can result in different section sizes.
 *
 * \param [in] p4est    the forest, whose quadrant data contains the state
 * \param [in] timestep the timestep number, used to name the output files
 */
static void
step3_write_checkpoint (p4est_t * p4est, int timestep)
{
  char                filename[BUFSIZ] = "";
  char                user_string[P4EST_FILE_USER_STRING_BYTES] = "";
  char                quad_data_user_string[P4EST_FILE_USER_STRING_BYTES] =
    "";
  int                 errcode;
  p4est_file_context_t *fc;
  uint32_t            check_endianness = STEP3_ENDIAN_CHECK;
  sc_array_t          block_arr;

  /** To write the data to the checkpoint file we need to store it
   * in a linear array. Therefore, we first need to create such a linear
   * array consisting of the quadrant data that is not stored in a linear
   * array in this example code. In this case one could also use
   * p4est_{load,save} to read and write the p4est including its
   * quadrant data but we use the p4est_file functions for
   * demonstration purposes. The p4est_file functions allow us to
   * store less data than p4est_{load,save} and the p4est_file
   * functions are required if an application uses external data
   * storage, that is there is data that is not stored as quadrant
   * data but for example in a linear array associated by its indexing
   * to quadrants of a p4est.
   */

  snprintf (filename, BUFSIZ,
            P4EST_STRING "_step3_checkpoint%04d." P4EST_DATA_FILE_EXT,
            timestep);

  fc = p4est_file_open_create (p4est, filename, "Checkpoint file", &errcode);
  /* One could use the error codes in \ref p4est_io.h and \ref p4est_file_error_string
   * for more sophisticated error handling.
   */
  SC_CHECK_ABORT (fc != NULL
                  && errcode == P4EST_FILE_ERR_SUCCESS,
                  P4EST_STRING "_file_open_create: Error creating file");

  snprintf (user_string, P4EST_FILE_USER_STRING_BYTES, "%s", "Endianness");

  sc_array_init_data (&block_arr, &check_endianness, sizeof (uint32_t), 1);
  /* write data to check endianness */
  fc =
    p4est_file_write_block (fc, sizeof (uint32_t), &block_arr,
                            user_string, &errcode);
  SC_CHECK_ABORT (fc != NULL
                  && errcode == P4EST_FILE_ERR_SUCCESS,
                  P4EST_STRING "_file_write_block: Error writing endianness");

  snprintf (user_string, P4EST_FILE_USER_STRING_BYTES, "%s",
            "Simulation context");

  sc_array_init_data (&block_arr, p4est->user_pointer, STEP3_BLOCK_SIZE, 1);
  /* write the simulation context */
  fc =
    p4est_file_write_block (fc, STEP3_BLOCK_SIZE, &block_arr,
                            user_string, &errcode);
  SC_CHECK_ABORT (fc != NULL
                  && errcode == P4EST_FILE_ERR_SUCCESS,
                  P4EST_STRING
                  "_file_write_block: Error writing simulation context");

  snprintf (user_string, P4EST_FILE_USER_STRING_BYTES,
            P4EST_STRING " connecitivity");

  /* write connectivity of the p4est */
  fc =
    p4est_file_write_connectivity (fc, p4est->connectivity, user_string,
                                   &errcode);
  SC_CHECK_ABORT (fc != NULL
                  && errcode == P4EST_FILE_ERR_SUCCESS,
                  P4EST_STRING
                  "_file_write_connectivity: Error writing connectivity");

  snprintf (user_string, P4EST_FILE_USER_STRING_BYTES,
            "Quadrants of time step %04d.", timestep);
  snprintf (quad_data_user_string, P4EST_FILE_USER_STRING_BYTES,
            "Quadrant data of time step %04d.", timestep);

  /** Write the current p4est to the checkpoint file; we do not write the
   * connectivity to disk because the connectivity is always the same in
   * this example and can be created again for each restart.
   */
  fc = p4est_file_write_p4est (fc, p4est,
                               user_string, quad_data_user_string, &errcode);
  SC_CHECK_ABORT (fc != NULL
                  && errcode == P4EST_FILE_ERR_SUCCESS,
                  P4EST_STRING "_file_write_field: Error writing p4est");

  p4est_file_close (fc, &errcode);
  SC_CHECK_ABORT (errcode == P4EST_FILE_ERR_SUCCESS,
                  P4EST_STRING "_file_close: Error closing file");
}

static void         step3_timestep (p4est_t * p4est, double start_time,
                                    double end_time);

/** Load a checkpoint file to restart the simulation.
 * The checkpoint file is not compiler independent since structures
 * may be padded by the compiler (e. g., sizeof (step3_data_t)).
 * This can result in different section sizes.
 *
 * \param [in] filename The file path to the checkpoint file
 *                      created using \ref step3_write_checkpoint.
 * \param [in] mpicomm  The MPI communciatior that is used for
 *                      the parallel simulation.
 * \param [in] time_inc The time increment added to the current time
 *                      of the restarted simulation to obtain the new
 *                      end time.
 */
static void
step3_restart (const char *filename, sc_MPI_Comm mpicomm, double time_inc)
{
  int                 errcode;
  char                user_string[P4EST_FILE_USER_STRING_BYTES],
    quad_string[P4EST_FILE_USER_STRING_BYTES],
    quad_data_string[P4EST_FILE_USER_STRING_BYTES];
  step3_ctx_t         ctx;
  p4est_gloidx_t      global_num_quadrants;
  p4est_file_context_t *fc;
  p4est_t            *loaded_p4est;
  p4est_connectivity_t *conn;
  uint32_t            check_endianness = STEP3_ENDIAN_CHECK;
  int32_t             read_check_endianness;
  sc_array_t          block_arr;

  fc =
    p4est_file_open_read_ext (mpicomm, filename, user_string,
                              &global_num_quadrants, &errcode);
  SC_CHECK_ABORT (fc != NULL
                  && errcode == P4EST_FILE_ERR_SUCCESS,
                  P4EST_STRING "_file_open_read: Error opening file");

  sc_array_init_data (&block_arr, &read_check_endianness, sizeof (uint32_t),
                      1);
  /* read data endianness */
  fc =
    p4est_file_read_block (fc, sizeof (uint32_t), &block_arr,
                           user_string, &errcode);
  SC_CHECK_ABORT (fc != NULL
                  && errcode == P4EST_FILE_ERR_SUCCESS,
                  P4EST_STRING
                  "_file_read_block: Error reading data endianness");
  P4EST_GLOBAL_PRODUCTIONF ("Read data with user string: %s\n", user_string);

  /* Instead of handling the wrong endianness, we just abort on it.
   * In a more sophisticated application the data should be converted.
   */
  SC_CHECK_ABORT (memcmp
                  (&read_check_endianness, &check_endianness,
                   sizeof (uint32_t)) == 0, "Wrong endianness");

  sc_array_init_data (&block_arr, &ctx, STEP3_BLOCK_SIZE, 1);
  /* read the simulation context */
  fc =
    p4est_file_read_block (fc, STEP3_BLOCK_SIZE, &block_arr, user_string,
                           &errcode);
  SC_CHECK_ABORT (fc != NULL
                  && errcode == P4EST_FILE_ERR_SUCCESS,
                  P4EST_STRING
                  "_file_read_block: Error reading simulation context");
  P4EST_GLOBAL_PRODUCTIONF ("Read data with user string: %s\n", user_string);

  /* read the connectivity */
  fc = p4est_file_read_connectivity (fc, &conn, user_string, &errcode);
  SC_CHECK_ABORT (fc != NULL
                  && errcode == P4EST_FILE_ERR_SUCCESS,
                  P4EST_STRING
                  "_file_read_connectivity: Error reading connectivity");
  P4EST_ASSERT (conn != NULL);
  P4EST_GLOBAL_PRODUCTIONF ("Read data with user string: %s\n", user_string);

  fc = p4est_file_read_p4est (fc, conn,
                              sizeof (step3_data_t),
                              &loaded_p4est, quad_string, quad_data_string,
                              &errcode);
  SC_CHECK_ABORT (fc != NULL
                  && errcode == P4EST_FILE_ERR_SUCCESS,
                  P4EST_STRING "_file_read_field_ext: Error reading p4est");
  P4EST_GLOBAL_PRODUCTIONF ("Read quadrants with user string: %s\n",
                            quad_string);
  P4EST_GLOBAL_PRODUCTIONF ("Read quadrant data with user string: %s\n",
                            quad_data_string);

  /* assign simulation context pointer */
  loaded_p4est->user_pointer = (void *) block_arr.array;

  /* close the file */
  p4est_file_close (fc, &errcode);
  SC_CHECK_ABORT (errcode == P4EST_FILE_ERR_SUCCESS,
                  P4EST_STRING "_file_close: Error closing file");

  step3_timestep (loaded_p4est, ctx.current_time,
                  ctx.current_time + time_inc);

  /* clean up */
  conn = loaded_p4est->connectivity;
  p4est_destroy (loaded_p4est);
  p4est_connectivity_destroy (conn);
}

/** Approximate the divergence of (vu) on each quadrant
 *
 * We use piecewise constant approximations on each quadrant, so the value is
 * always 0.
 *
 * Like step3_interpolate_solution(), this function matches the
 * p4est_iter_volume_t prototype used by p4est_iterate().
 *
 * \param [in] info          the information about the quadrant populated by
 *                           p4est_iterate()
 * \param [in] user_data     not used
 */
static void
step3_quad_divergence (p4est_iter_volume_info_t * info, void *user_data)
{
  p4est_quadrant_t   *q = info->quad;
  step3_data_t       *data = (step3_data_t *) q->p.user_data;

  data->dudt = 0.;
}

/** Approximate the flux across a boundary between quadrants.
 *
 * We use a very simple upwind numerical flux.
 *
 * This function matches the p4est_iter_face_t prototype used by
 * p4est_iterate().
 *
 * \param [in] info the information about the quadrants on either side of the
 *                  interface, populated by p4est_iterate()
 * \param [in] user_data the user_data given to p4est_iterate(): in this case,
 *                       it points to the ghost_data array, which contains the
 *                       step3_data_t data for all of the ghost cells, which
 *                       was populated by p4est_ghost_exchange_data()
 */
static void
step3_upwind_flux (p4est_iter_face_info_t * info, void *user_data)
{
  int                 i, j;
  p4est_t            *p4est = info->p4est;
  step3_ctx_t        *ctx = (step3_ctx_t *) p4est->user_pointer;
  step3_data_t       *ghost_data = (step3_data_t *) user_data;
  step3_data_t       *udata;
  p4est_quadrant_t   *quad;
  double              vdotn = 0.;
  double              uavg;
  double              q;
  double              h, facearea;
  int                 which_face;
  int                 upwindside;
  p4est_iter_face_side_t *side[2];
  sc_array_t         *sides = &(info->sides);

  /* because there are no boundaries, every face has two sides */
  P4EST_ASSERT (sides->elem_count == 2);

  side[0] = p4est_iter_fside_array_index_int (sides, 0);
  side[1] = p4est_iter_fside_array_index_int (sides, 1);

  /* which of the quadrant's faces the interface touches */
  which_face = side[0]->face;

  switch (which_face) {
  case 0:                      /* -x side */
    vdotn = -ctx->v[0];
    break;
  case 1:                      /* +x side */
    vdotn = ctx->v[0];
    break;
  case 2:                      /* -y side */
    vdotn = -ctx->v[1];
    break;
  case 3:                      /* +y side */
    vdotn = ctx->v[1];
    break;
#ifdef P4_TO_P8
  case 4:                      /* -z side */
    vdotn = -ctx->v[2];
    break;
  case 5:                      /* +z side */
    vdotn = ctx->v[2];
    break;
#endif
  }
  upwindside = vdotn >= 0. ? 0 : 1;

  /* Because we have non-conforming boundaries, one side of an interface can
   * either have one large ("full") quadrant or 2^(d-1) small ("hanging")
   * quadrants: we have to compute the average differently in each case.  The
   * info populated by p4est_iterate() gives us the context we need to
   * proceed. */
  uavg = 0;
  if (side[upwindside]->is_hanging) {
    /* there are 2^(d-1) (P4EST_HALF) subfaces */
    for (j = 0; j < P4EST_HALF; j++) {
      if (side[upwindside]->is.hanging.is_ghost[j]) {
        /* *INDENT-OFF* */
        udata =
          (step3_data_t *) &ghost_data[side[upwindside]->is.hanging.quadid[j]];
        /* *INDENT-ON* */
      }
      else {
        udata =
          (step3_data_t *) side[upwindside]->is.hanging.quad[j]->p.user_data;
      }
      uavg += udata->u;
    }
    uavg /= P4EST_HALF;
  }
  else {
    if (side[upwindside]->is.full.is_ghost) {
      udata = (step3_data_t *) & ghost_data[side[upwindside]->is.full.quadid];
    }
    else {
      udata = (step3_data_t *) side[upwindside]->is.full.quad->p.user_data;
    }
    uavg = udata->u;
  }
  /* flux from side 0 to side 1 */
  q = vdotn * uavg;
  for (i = 0; i < 2; i++) {
    if (side[i]->is_hanging) {
      /* there are 2^(d-1) (P4EST_HALF) subfaces */
      for (j = 0; j < P4EST_HALF; j++) {
        quad = side[i]->is.hanging.quad[j];
        h =
          (double) P4EST_QUADRANT_LEN (quad->level) / (double) P4EST_ROOT_LEN;
#ifndef P4_TO_P8
        facearea = h;
#else
        facearea = h * h;
#endif
        if (!side[i]->is.hanging.is_ghost[j]) {
          udata = (step3_data_t *) quad->p.user_data;
          if (i == upwindside) {
            udata->dudt += vdotn * udata->u * facearea * (i ? 1. : -1.);
          }
          else {
            udata->dudt += q * facearea * (i ? 1. : -1.);
          }
        }
      }
    }
    else {
      quad = side[i]->is.full.quad;
      h = (double) P4EST_QUADRANT_LEN (quad->level) / (double) P4EST_ROOT_LEN;
#ifndef P4_TO_P8
      facearea = h;
#else
      facearea = h * h;
#endif
      if (!side[i]->is.full.is_ghost) {
        udata = (step3_data_t *) quad->p.user_data;
        udata->dudt += q * facearea * (i ? 1. : -1.);
      }
    }
  }
}

/** Compute the new value of the state from the computed time derivative.
 *
 * We use a simple forward Euler scheme.
 *
 * The derivative was computed by a p4est_iterate() loop by the callbacks
 * step3_quad_divergence() and step3_upwind_flux(). Now we multiply this by
 * the timestep and add to the current solution.
 *
 * This function matches the p4est_iter_volume_t prototype used by
 * p4est_iterate().
 *
 * \param [in] info          the information about this quadrant that has been
 *                           populated by p4est_iterate()
 * \param [in] user_data the user_data given to p4est_iterate(): in this case,
 *                       it points to the timestep.
 */
static void
step3_timestep_update (p4est_iter_volume_info_t * info, void *user_data)
{
  p4est_quadrant_t   *q = info->quad;
  step3_data_t       *data = (step3_data_t *) q->p.user_data;
  double              dt = *((double *) user_data);
  double              vol;
  double              h =
    (double) P4EST_QUADRANT_LEN (q->level) / (double) P4EST_ROOT_LEN;

#ifdef P4_TO_P8
  vol = h * h * h;
#else
  vol = h * h;
#endif

  data->u += dt * data->dudt / vol;
}

/** Reset the approximate derivatives.
 *
 * p4est_iterate() has an invariant to the order of callback execution: the
 * p4est_iter_volume_t callback will be executed on a quadrant before the
 * p4est_iter_face_t callbacks are executed on its faces.  This function
 * resets the derivative stored in the quadrant's data before
 * step3_minmod_estimate() updates the derivative based on the face neighbors.
 *
 * This function matches the p4est_iter_volume_t prototype used by
 * p4est_iterate().
 *
 * \param [in] info          the information about this quadrant that has been
 *                           populated by p4est_iterate()
 * \param [in] user_data     not used
 */
static void
step3_reset_derivatives (p4est_iter_volume_info_t * info, void *user_data)
{
  p4est_quadrant_t   *q = info->quad;
  step3_data_t       *data = (step3_data_t *) q->p.user_data;
  int                 j;

  for (j = 0; j < P4EST_DIM; j++) {
    data->du[j] = step3_invalid;
  }
}

/** For two quadrants on either side of a face, estimate the derivative normal
 * to the face.
 *
 * This function matches the p4est_iter_face_t prototype used by
 * p4est_iterate().
 *
 * \param [in] info          the information about this quadrant that has been
 *                           populated by p4est_iterate()
 * \param [in] user_data the user_data given to p4est_iterate(): in this case,
 *                       it points to the ghost_data array, which contains the
 *                       step3_data_t data for all of the ghost cells, which
 *                       was populated by p4est_ghost_exchange_data()
 */
static void
step3_minmod_estimate (p4est_iter_face_info_t * info, void *user_data)
{
  int                 i, j;
  p4est_iter_face_side_t *side[2];
  sc_array_t         *sides = &(info->sides);
  step3_data_t       *ghost_data = (step3_data_t *) user_data;
  step3_data_t       *udata;
  p4est_quadrant_t   *quad;
  double              uavg[2];
  double              h[2];
  double              du_est, du_old;
  int                 which_dir;

  /* because there are no boundaries, every face has two sides */
  P4EST_ASSERT (sides->elem_count == 2);

  side[0] = p4est_iter_fside_array_index_int (sides, 0);
  side[1] = p4est_iter_fside_array_index_int (sides, 1);

  which_dir = side[0]->face / 2;        /* 0 == x, 1 == y, 2 == z */

  for (i = 0; i < 2; i++) {
    uavg[i] = 0;
    if (side[i]->is_hanging) {
      /* there are 2^(d-1) (P4EST_HALF) subfaces */
      for (j = 0; j < P4EST_HALF; j++) {
        quad = side[i]->is.hanging.quad[j];
        h[i] =
          (double) P4EST_QUADRANT_LEN (quad->level) / (double) P4EST_ROOT_LEN;
        if (side[i]->is.hanging.is_ghost[j]) {
          udata = &ghost_data[side[i]->is.hanging.quadid[j]];
        }
        else {
          udata = (step3_data_t *) side[i]->is.hanging.quad[j]->p.user_data;
        }
        uavg[i] += udata->u;
      }
      uavg[i] /= P4EST_HALF;
    }
    else {
      quad = side[i]->is.full.quad;
      h[i] =
        (double) P4EST_QUADRANT_LEN (quad->level) / (double) P4EST_ROOT_LEN;
      if (side[i]->is.full.is_ghost) {
        udata = &ghost_data[side[i]->is.full.quadid];
      }
      else {
        udata = (step3_data_t *) side[i]->is.full.quad->p.user_data;
      }
      uavg[i] = udata->u;
    }
  }
  du_est = (uavg[1] - uavg[0]) / ((h[0] + h[1]) / 2.);
  for (i = 0; i < 2; i++) {
    if (side[i]->is_hanging) {
      /* there are 2^(d-1) (P4EST_HALF) subfaces */
      for (j = 0; j < P4EST_HALF; j++) {
        quad = side[i]->is.hanging.quad[j];
        if (!side[i]->is.hanging.is_ghost[j]) {
          udata = (step3_data_t *) quad->p.user_data;
          du_old = udata->du[which_dir];
          if (du_old == du_old) {
            /* there has already been an update */
            if (du_est * du_old >= 0.) {
              if (fabs (du_est) < fabs (du_old)) {
                udata->du[which_dir] = du_est;
              }
            }
            else {
              udata->du[which_dir] = 0.;
            }
          }
          else {
            udata->du[which_dir] = du_est;
          }
        }
      }
    }
    else {
      quad = side[i]->is.full.quad;
      if (!side[i]->is.full.is_ghost) {
        udata = (step3_data_t *) quad->p.user_data;
        du_old = udata->du[which_dir];
        if (du_old == du_old) {
          /* there has already been an update */
          if (du_est * du_old >= 0.) {
            if (fabs (du_est) < fabs (du_old)) {
              udata->du[which_dir] = du_est;
            }
          }
          else {
            udata->du[which_dir] = 0.;
          }
        }
        else {
          udata->du[which_dir] = du_est;
        }
      }
    }
  }
}

/** Compute the maximum state value.
 *
 * This function updates the maximum value from the value of a single cell.
 *
 * This function matches the p4est_iter_volume_t prototype used by
 * p4est_iterate().
 *
 * \param [in] info              the information about this quadrant that has been
 *                               populated by p4est_iterate()
 * \param [in,out] user_data     the user_data given to p4est_iterate(): in this case,
 *                               it points to the maximum value that will be updated
 */
static void
step3_compute_max (p4est_iter_volume_info_t * info, void *user_data)
{
  p4est_quadrant_t   *q = info->quad;
  step3_data_t       *data = (step3_data_t *) q->p.user_data;
  double              umax = *((double *) user_data);

  umax = SC_MAX (data->u, umax);

  *((double *) user_data) = umax;
}

/** Compute the timestep.
 *
 * Find the smallest quadrant and scale the timestep based on that length and
 * the advection velocity.
 *
 * \param [in] p4est the forest
 * \return the timestep.
 */
static double
step3_get_timestep (p4est_t * p4est)
{
  step3_ctx_t        *ctx = (step3_ctx_t *) p4est->user_pointer;
  p4est_topidx_t      t, flt, llt;
  p4est_tree_t       *tree;
  int                 max_level, global_max_level;
  int                 mpiret, i;
  double              min_h, vnorm;
  double              dt;

  /* compute the timestep by finding the smallest quadrant */
  flt = p4est->first_local_tree;
  llt = p4est->last_local_tree;

  max_level = 0;
  for (t = flt; t <= llt; t++) {
    tree = p4est_tree_array_index (p4est->trees, t);
    max_level = SC_MAX (max_level, tree->maxlevel);

  }
  mpiret =
    sc_MPI_Allreduce (&max_level, &global_max_level, 1, sc_MPI_INT,
                      sc_MPI_MAX, p4est->mpicomm);
  SC_CHECK_MPI (mpiret);

  min_h =
    (double) P4EST_QUADRANT_LEN (global_max_level) / (double) P4EST_ROOT_LEN;

  vnorm = 0;
  for (i = 0; i < P4EST_DIM; i++) {
    vnorm += ctx->v[i] * ctx->v[i];
  }
  vnorm = sqrt (vnorm);

  dt = min_h / 2. / vnorm;

  return dt;
}

/** Timestep the advection problem.
 *
 * Update the state, refine, repartition, and write the solution to file.
 *
 * \param [in,out] p4est the forest, whose state is updated
 * \param [in] start_time    simulation start time
 * \param [in] end_time      simulation end time
 */
static void
step3_timestep (p4est_t * p4est, double start_time, double end_time)
{
  double              t = start_time;
  double              dt = 0.;
  int                 i;
  step3_data_t       *ghost_data;
  step3_ctx_t        *ctx = (step3_ctx_t *) p4est->user_pointer;
  int                 refine_period = ctx->refine_period;
  int                 repartition_period = ctx->repartition_period;
  int                 write_period = ctx->write_period;
  int                 recursive = 0;
  int                 allowed_level = P4EST_QMAXLEVEL;
  int                 allowcoarsening = 1;
  int                 callbackorphans = 0;
  int                 mpiret;
  double              orig_max_err = ctx->max_err;
  double              umax, global_umax;
  p4est_ghost_t      *ghost;

  /* create the ghost quadrants */
  ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FULL);
  /* create space for storing the ghost data */
  ghost_data = P4EST_ALLOC (step3_data_t, ghost->ghosts.elem_count);
  /* synchronize the ghost data */
  p4est_ghost_exchange_data (p4est, ghost, ghost_data);

  /* initialize du/dx estimates */
  p4est_iterate (p4est, ghost, (void *) ghost_data,     /* pass in ghost data that we just exchanged */
                 step3_reset_derivatives,       /* blank the previously calculated derivatives */
                 step3_minmod_estimate, /* compute the minmod estimate of each cell's derivative */
#ifdef P4_TO_P8
                 NULL,          /* there is no callback for the edges between quadrants */
#endif
                 NULL);         /* there is no callback for the corners between quadrants */

  for (t = start_time, i = ctx->time_step; t < end_time; t += dt, i++) {
    P4EST_GLOBAL_PRODUCTIONF ("time %f\n", t);

    /* refine */
    if (!(i % refine_period)) {
      if (i) {
        /* compute umax */
        umax = 0.;
        /* initialize derivative estimates */
        p4est_iterate (p4est, NULL, (void *) &umax,     /* pass in ghost data that we just exchanged */
                       step3_compute_max,       /* blank the previously calculated derivatives */
                       NULL,    /* there is no callback for the faces between quadrants */
#ifdef P4_TO_P8
                       NULL,    /* there is no callback for the edges between quadrants */
#endif
                       NULL);   /* there is no callback for the corners between quadrants */

        mpiret =
          sc_MPI_Allreduce (&umax, &global_umax, 1, sc_MPI_DOUBLE, sc_MPI_MAX,
                            p4est->mpicomm);
        SC_CHECK_MPI (mpiret);
        ctx->max_err = orig_max_err * global_umax;
        P4EST_GLOBAL_PRODUCTIONF ("u_max %f\n", global_umax);

        /* adapt */
        p4est_refine_ext (p4est, recursive, allowed_level,
                          step3_refine_err_estimate, NULL,
                          step3_replace_quads);
        p4est_coarsen_ext (p4est, recursive, callbackorphans,
                           step3_coarsen_err_estimate, NULL,
                           step3_replace_quads);
        p4est_balance_ext (p4est, P4EST_CONNECT_FACE, NULL,
                           step3_replace_quads);

        p4est_ghost_destroy (ghost);
        P4EST_FREE (ghost_data);
        ghost = NULL;
        ghost_data = NULL;
      }
      dt = step3_get_timestep (p4est);
    }

    /* repartition */
    if (i && !(i % repartition_period)) {
      p4est_partition (p4est, allowcoarsening, NULL);

      if (ghost) {
        p4est_ghost_destroy (ghost);
        P4EST_FREE (ghost_data);
        ghost = NULL;
        ghost_data = NULL;
      }
    }

    /* write out solution */
    if (!(i % write_period)) {
      step3_write_solution (p4est, i);
    }

    /* synchronize the ghost data */
    if (!ghost) {
      ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FULL);
      ghost_data = P4EST_ALLOC (step3_data_t, ghost->ghosts.elem_count);
      p4est_ghost_exchange_data (p4est, ghost, ghost_data);
    }

    /* compute du/dt */
    /* *INDENT-OFF* */
    p4est_iterate (p4est,                 /* the forest */
                   ghost,                 /* the ghost layer */
                   (void *) ghost_data,   /* the synchronized ghost data */
                   step3_quad_divergence, /* callback to compute each quad's
                                             interior contribution to du/dt */
                   step3_upwind_flux,     /* callback to compute each quads'
                                             faces' contributions to du/du */
#ifdef P4_TO_P8
                   NULL,                  /* there is no callback for the
                                             edges between quadrants */
#endif
                   NULL);                 /* there is no callback for the
                                             corners between quadrants */
    /* *INDENT-ON* */

    /* update u */
    p4est_iterate (p4est, NULL, /* ghosts are not needed for this loop */
                   (void *) &dt,        /* pass in dt */
                   step3_timestep_update,       /* update each cell */
                   NULL,        /* there is no callback for the faces between quadrants */
#ifdef P4_TO_P8
                   NULL,        /* there is no callback for the edges between quadrants */
#endif
                   NULL);       /* there is no callback for the corners between quadrants */

    /* synchronize the ghost data */
    p4est_ghost_exchange_data (p4est, ghost, ghost_data);

    /* update du/dx estimate */
    p4est_iterate (p4est, ghost, (void *) ghost_data,   /* pass in ghost data that we just exchanged */
                   step3_reset_derivatives,     /* blank the previously calculated derivatives */
                   step3_minmod_estimate,       /* compute the minmod estimate of each cell's derivative */
#ifdef P4_TO_P8
                   NULL,        /* there is no callback for the edges between quadrants */
#endif
                   NULL);       /* there is no callback for the corners between quadrants */

  }

  ctx->time_step = i - 1;
  ctx->current_time = t - dt;
  if (step3_checkpoint) {
    /* write checkpoint file */
    step3_write_checkpoint (p4est, i - 1);
  }

  P4EST_FREE (ghost_data);
  p4est_ghost_destroy (ghost);
}

/** The main step 3 program.
 *
 * Setup of the example parameters; create the forest, with the state variable
 * stored in the quadrant data; refine, balance, and partition the forest;
 * timestep; clean up, and exit.
 */
int
main (int argc, char **argv)
{
  int                 mpiret;
  int                 recursive, partforcoarsen;
  sc_MPI_Comm         mpicomm;
  p4est_t            *p4est;
  p4est_connectivity_t *conn;
  step3_ctx_t         ctx;
  sc_options_t       *opt;
  const char         *filename;

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

  /* read command line options */
  opt = sc_options_new (argv[0]);
  sc_options_add_bool (opt, 'C', "write-checkpoint", &step3_checkpoint, 0,
                       "Write checkpoint files to disk");
  sc_options_add_string (opt, 'L', "load-checkpoint", &filename,
                         NULL, "Load and start from a checkpoint file");
  sc_options_parse (p4est_package_id, SC_LP_DEFAULT, opt, argc, argv);
  sc_options_print_usage (p4est_package_id, SC_LP_ERROR, opt, NULL);
  sc_options_print_summary (p4est_package_id, SC_LP_ESSENTIAL, opt);

  if (filename != NULL) {
    /* load a checkpoint file and restart the simulation */
    step3_restart (filename, mpicomm, 0.1);

    sc_options_destroy (opt);
    sc_finalize ();

    mpiret = sc_MPI_Finalize ();
    SC_CHECK_MPI (mpiret);
    return 0;
  }

  /* Avoid write of uninitialized bytes (valgrind warning)
   * due to compiler padding.
   */
  memset (&ctx, -1, sizeof (ctx));

  ctx.bump_width = 0.1;
  ctx.max_err = 2.e-2;
  ctx.center[0] = 0.5;
  ctx.center[1] = 0.5;
#ifdef P4_TO_P8
  ctx.center[2] = 0.5;
#endif
#ifndef P4_TO_P8
  /* randomly chosen advection direction */
  ctx.v[0] = -0.445868402501118;
  ctx.v[1] = -0.895098523991131;
#else
  ctx.v[0] = 0.485191768970225;
  ctx.v[1] = -0.427996381877778;
  ctx.v[2] = 0.762501176669961;
#endif
  ctx.refine_period = 2;
  ctx.repartition_period = 4;
  ctx.write_period = 8;
  ctx.time_step = 0;

  /* Create a forest that consists of just one periodic quadtree/octree. */
#ifndef P4_TO_P8
  conn = p4est_connectivity_new_periodic ();
#else
  conn = p8est_connectivity_new_periodic ();
#endif

  /* *INDENT-OFF* */
  p4est = p4est_new_ext (mpicomm, /* communicator */
                         conn,    /* connectivity */
                         0,       /* minimum quadrants per MPI process */
                         4,       /* minimum level of refinement */
                         1,       /* fill uniform */
                         sizeof (step3_data_t),         /* data size */
                         step3_init_initial_condition,  /* initializes data */
                         (void *) (&ctx));              /* context */
  /* *INDENT-ON* */

  /* refine and coarsen based on an interpolation error estimate */
  recursive = 1;
  p4est_refine (p4est, recursive, step3_refine_err_estimate,
                step3_init_initial_condition);
  p4est_coarsen (p4est, recursive, step3_coarsen_initial_condition,
                 step3_init_initial_condition);

  /* Partition: The quadrants are redistributed for equal element count.  The
   * partition can optionally be modified such that a family of octants, which
   * are possibly ready for coarsening, are never split between processors. */
  partforcoarsen = 1;

  /* If we call the 2:1 balance we ensure that neighbors do not differ in size
   * by more than a factor of 2.  This can optionally include diagonal
   * neighbors across edges or corners as well; see p4est.h. */
  p4est_balance (p4est, P4EST_CONNECT_FACE, step3_init_initial_condition);
  p4est_partition (p4est, partforcoarsen, NULL);

  /* time step */
  step3_timestep (p4est, 0., 0.1);

  /* Destroy the p4est and the connectivity structure. */
  p4est_destroy (p4est);
  p4est_connectivity_destroy (conn);

  sc_options_destroy (opt);

  /* Verify that allocations internal to p4est and sc do not leak memory.
   * This should be called if sc_init () has been called earlier. */
  sc_finalize ();

  /* This is standard MPI programs.  Without --enable-mpi, this is a dummy. */
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;
}
