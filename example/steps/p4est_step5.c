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

/** \file p4est_step5.c
 *
 * This 2D example program illustrates high-order visualization.
 * Contributed by: Grant Seastream, ExxonMobil
 */

#ifndef P4_TO_P8
#include <p4est_extended.h>
#include <p4est_vtk.h>
#else
#include <p8est_extended.h>
#include <p8est_vtk.h>
#endif

/* change the following to #if 0 to revert to low-order visualization */
#if 1
#define STEP5_HO
#endif

#define STEP5_NNODE_1D 8
#ifndef P4_TO_P8
#define STEP5_NNODE (STEP5_NNODE_1D * STEP5_NNODE_1D)
#else
#define STEP5_NNODE (STEP5_NNODE_1D * STEP5_NNODE_1D * STEP5_NNODE_1D)
#endif

/** Per-quadrant data for this example.
 *
 * We store coordinates and a data value for visualization.
 */
typedef struct step5_data
{
  /* Element coordinates in the form of x_transpose = [ x_0 ... x_n,
   * y_0 ... y_n, z_0 ... z_n] */
  double              xt[P4EST_DIM * STEP5_NNODE];
  double              data[STEP5_NNODE];        /* data values for points */
}
step5_data_t;

/* --------------------------------------------------------------------------
 * Retrieve info about GLL integration for specific number of points (in 1D)
 * See the following link to generate x and w:
 * https://www.mathworks.com/matlabcentral/fileexchange/
 * 4775-legende-gauss-lobatto-nodes-and-weights
 * --------------------------------------------------------------------------*/
static void
step5_get_GLL_info (double *x, double *w, const int num_points)
{
  switch (num_points) {
  case 2:
    x[0] = -1.0;
    x[1] = +1.0;

    w[0] = 1.0;
    w[1] = 1.0;
    break;

  case 3:
    x[0] = -1.0;
    x[1] = 0.0;
    x[2] = +1.0;

    w[0] = 0.333333333333333;
    w[1] = 1.333333333333333;
    w[2] = 0.333333333333333;
    break;

  case 4:
    x[0] = -1.0;
    x[1] = -0.447213595499957;
    x[2] = +0.447213595499957;
    x[3] = +1.0;

    w[0] = 0.166666666666667;
    w[1] = 0.833333333333333;
    w[2] = 0.833333333333333;
    w[3] = 0.166666666666667;
    break;

  case 5:
    x[0] = -1.0;
    x[1] = -0.654653670707977;
    x[2] = 0.0;
    x[3] = +0.654653670707977;
    x[4] = +1.0;

    w[0] = 0.1;
    w[1] = 0.544444444444444;
    w[2] = 0.711111111111111;
    w[3] = 0.544444444444444;
    w[4] = 0.1;
    break;

  case 6:
    x[0] = -1.0;
    x[1] = -0.765055323929464;
    x[2] = -0.285231516480645;
    x[3] = +0.285231516480645;
    x[4] = +0.765055323929464;
    x[5] = +1.0;

    w[0] = 0.0666666666666667;
    w[1] = 0.378474956297846;
    w[2] = 0.554858377035486;
    w[3] = 0.554858377035486;
    w[4] = 0.378474956297846;
    w[5] = 0.0666666666666667;
    break;

  case 7:
    x[0] = -1.0;
    x[1] = -0.830223896278567;
    x[2] = -0.468848793470714;
    x[3] = 0.0;
    x[4] = +0.468848793470714;
    x[5] = +0.830223896278567;
    x[6] = +1.0;

    w[0] = 0.047619047619048;
    w[1] = 0.276826047361566;
    w[2] = 0.431745381209863;
    w[3] = 0.487619047619048;
    w[4] = 0.431745381209863;
    w[5] = 0.276826047361566;
    w[6] = 0.047619047619048;
    break;

  case 8:
    x[0] = -1.0;
    x[1] = -0.871740148509607;
    x[2] = -0.591700181433142;
    x[3] = -0.209299217902479;
    x[4] = +0.209299217902479;
    x[5] = +0.591700181433142;
    x[6] = +0.871740148509607;
    x[7] = +1.0;

    w[0] = 0.035714285714286;
    w[1] = 0.210704227143506;
    w[2] = 0.341122692483504;
    w[3] = 0.412458794658704;
    w[4] = 0.412458794658704;
    w[5] = 0.341122692483504;
    w[6] = 0.210704227143506;
    w[7] = 0.035714285714286;
    break;
  }
}

/* ----------------------------------------------------------------------------
 * convert p4est coordinate to physical coordinates for all 4 or 8 vertices (2D
 * or 3D)
 * --------------------------------------------------------------------------*/
static void
step5_qcoord_to_vertex_all (const p4est_t * p4est,
                            const p4est_topidx_t which_tree,
                            const p4est_quadrant_t * q, double vxy_all[])
{
  /* *INDENT-OFF* */
  /* Transform a quadrant coordinate into the space spanned by tree vertices.
     see p4est vertex convention:

       2----3             6----7
       |    |            /    /
       |    |           4----5
       0----1           |
                        | 2----3
                        |/    /
                        0----1
   */
  /* *INDENT-ON* */
#ifndef P4_TO_P8
  /* each call will write to vxy_all[i], vxy_all[i+1], and vxy_all[i+2].
   * we'll let each next call overwrite the previous 3rd index since we don't
   * need it. we pass a temporary array of 3 for the final call. */
  double              vxy_temp[3];

  /* p4est vertex 0 */
  p4est_qcoord_to_vertex (p4est->connectivity, which_tree, q->x, q->y,
                          &vxy_all[0]);
  /* p4est vertex 1 */
  p4est_qcoord_to_vertex (p4est->connectivity,
                          which_tree,
                          q->x + P4EST_QUADRANT_LEN (q->level),
                          q->y, &vxy_all[2]);
  /* p4est vertex 2 */
  p4est_qcoord_to_vertex (p4est->connectivity,
                          which_tree,
                          q->x,
                          q->y + P4EST_QUADRANT_LEN (q->level), &vxy_all[4]);
  /* p4est vertex 3 */
  p4est_qcoord_to_vertex (p4est->connectivity,
                          which_tree,
                          q->x + P4EST_QUADRANT_LEN (q->level),
                          q->y + P4EST_QUADRANT_LEN (q->level), vxy_temp);
  vxy_all[6] = vxy_temp[0];
  vxy_all[7] = vxy_temp[1];
#else
  /* p4est vertex 0 */
  p4est_qcoord_to_vertex (p4est->connectivity, which_tree, q->x, q->y, q->z,
                          &vxy_all[0]);
  /* p4est vertex 1 */
  p4est_qcoord_to_vertex (p4est->connectivity,
                          which_tree,
                          q->x + P4EST_QUADRANT_LEN (q->level),
                          q->y, q->z, &vxy_all[3]);
  /* p4est vertex 2 */
  p4est_qcoord_to_vertex (p4est->connectivity,
                          which_tree,
                          q->x,
                          q->y + P4EST_QUADRANT_LEN (q->level),
                          q->z, &vxy_all[6]);
  /* p4est vertex 3 */
  p4est_qcoord_to_vertex (p4est->connectivity,
                          which_tree,
                          q->x + P4EST_QUADRANT_LEN (q->level),
                          q->y + P4EST_QUADRANT_LEN (q->level),
                          q->z, &vxy_all[9]);
  /* p4est vertex 4 */
  p4est_qcoord_to_vertex (p4est->connectivity,
                          which_tree,
                          q->x,
                          q->y,
                          q->z + P4EST_QUADRANT_LEN (q->level), &vxy_all[12]);
  /* p4est vertex 5 */
  p4est_qcoord_to_vertex (p4est->connectivity,
                          which_tree,
                          q->x + P4EST_QUADRANT_LEN (q->level),
                          q->y,
                          q->z + P4EST_QUADRANT_LEN (q->level), &vxy_all[15]);
  /* p4est vertex 6 */
  p4est_qcoord_to_vertex (p4est->connectivity,
                          which_tree,
                          q->x,
                          q->y + P4EST_QUADRANT_LEN (q->level),
                          q->z + P4EST_QUADRANT_LEN (q->level), &vxy_all[18]);
  /* p4est vertex 7 */
  p4est_qcoord_to_vertex (p4est->connectivity,
                          which_tree,
                          q->x + P4EST_QUADRANT_LEN (q->level),
                          q->y + P4EST_QUADRANT_LEN (q->level),
                          q->z + P4EST_QUADRANT_LEN (q->level), &vxy_all[21]);
#endif
}

#ifndef P4_TO_P8
/* ----------------------------------------------------------------------------
 * 2D linear shape functions evaluated at given coordinates
 * --------------------------------------------------------------------------*/
static void
step5_shape_p1 (const double r, const double s, double *fn)
{
  fn[0] = 0.25 * (1.0 - r) * (1.0 - s);
  fn[1] = 0.25 * (1.0 + r) * (1.0 - s);
  fn[2] = 0.25 * (1.0 - r) * (1.0 + s);
  fn[3] = 0.25 * (1.0 + r) * (1.0 + s);

}
#else
/* ----------------------------------------------------------------------------
 * 3D linear shape functions evaluated at given coordinates
 * --------------------------------------------------------------------------*/
static void
step5_shape_p1 (const double r, const double s, const double t, double *fn)
{
  fn[0] = 0.125 * (1.0 - r) * (1.0 - s) * (1.0 - t);
  fn[1] = 0.125 * (1.0 + r) * (1.0 - s) * (1.0 - t);
  fn[2] = 0.125 * (1.0 - r) * (1.0 + s) * (1.0 - t);
  fn[3] = 0.125 * (1.0 + r) * (1.0 + s) * (1.0 - t);
  fn[4] = 0.125 * (1.0 - r) * (1.0 - s) * (1.0 + t);
  fn[5] = 0.125 * (1.0 + r) * (1.0 - s) * (1.0 + t);
  fn[6] = 0.125 * (1.0 - r) * (1.0 + s) * (1.0 + t);
  fn[7] = 0.125 * (1.0 + r) * (1.0 + s) * (1.0 + t);

}
#endif

/* ----------------------------------------------------------------------------
 * construct GLL xyz for high-order elements
 * --------------------------------------------------------------------------*/
static void
step5_construct_GLL_xyz (const double vxy_all[],
                         const int xt_nrows_dim,
                         const int xt_ncols_nnode, double *xt)
{
  /* *INDENT-OFF* */
  /* construct GLL points for p2 elements given p4est physical coordinates -
     see convention below.

       2----3             6----7               ---- 3 ----
       |    |            /    /               /          /
       |    |           4----5               0          1
       0----1           |                   /          /
                        | 2----3           ---- 2 ----
                        |/    /
                        0----1
   */
  /* *INDENT-ON* */

  int                 i, j, ii, jj, visited_node, i_dim;
  double              xi, xj;
#ifdef P4_TO_P8
  int                 k;
  double              xk;
#endif
  /* import vertex geometry corresponding to a p1 element: */
  double              xt_p1[P4EST_DIM][P4EST_CHILDREN];
  /* compute geometry for high-order elements: */
  double              xint[STEP5_NNODE_1D];
  double              wint[STEP5_NNODE_1D];
  double              fn[P4EST_CHILDREN];       /* array containing shape
                                                 * functions for 3D p1 
                                                 * elements */
  double              xy_computed[P4EST_DIM];

  step5_get_GLL_info (xint, wint, STEP5_NNODE_1D);

#ifndef P4_TO_P8
  /* fill the corners */
  xt_p1[0][0] = vxy_all[0];
  xt_p1[1][0] = vxy_all[1];
  xt_p1[0][1] = vxy_all[2];
  xt_p1[1][1] = vxy_all[3];
  xt_p1[0][2] = vxy_all[4];
  xt_p1[1][2] = vxy_all[5];
  xt_p1[0][3] = vxy_all[6];
  xt_p1[1][3] = vxy_all[7];

  for (j = 0; j < STEP5_NNODE_1D; ++j) {
    /* y-direction in reference coordinate {s} */
    xj = xint[j];

    for (i = 0; i < STEP5_NNODE_1D; ++i) {
      /* x-direction in reference coordinate {r} */
      xi = xint[i];

      /* p1 shape functions are evaluated at GLL points of the desired
       * polynomial degree */
      step5_shape_p1 (xi, xj, fn);

      for (ii = 0; ii < P4EST_DIM; ++ii) {
        xy_computed[ii] = 0.0;
        for (jj = 0; jj < P4EST_CHILDREN; ++jj) {
          xy_computed[ii] += xt_p1[ii][jj] * fn[jj];
        }
      }

      /* construct geometry for a high-order element:
       * find node that is being visited */
      visited_node = i + j * STEP5_NNODE_1D;
      for (i_dim = 0; i_dim < xt_nrows_dim; i_dim++) {
        xt[visited_node + i_dim * xt_ncols_nnode] = xy_computed[i_dim];
      }
    }
  }
#else
  /* fill the corners */
  xt_p1[0][0] = vxy_all[0];
  xt_p1[1][0] = vxy_all[1];
  xt_p1[2][0] = vxy_all[2];
  xt_p1[0][1] = vxy_all[3];
  xt_p1[1][1] = vxy_all[4];
  xt_p1[2][1] = vxy_all[5];
  xt_p1[0][2] = vxy_all[6];
  xt_p1[1][2] = vxy_all[7];
  xt_p1[2][2] = vxy_all[8];
  xt_p1[0][3] = vxy_all[9];
  xt_p1[1][3] = vxy_all[10];
  xt_p1[2][3] = vxy_all[11];
  xt_p1[0][4] = vxy_all[12];
  xt_p1[1][4] = vxy_all[13];
  xt_p1[2][4] = vxy_all[14];
  xt_p1[0][5] = vxy_all[15];
  xt_p1[1][5] = vxy_all[16];
  xt_p1[2][5] = vxy_all[17];
  xt_p1[0][6] = vxy_all[18];
  xt_p1[1][6] = vxy_all[19];
  xt_p1[2][6] = vxy_all[20];
  xt_p1[0][7] = vxy_all[21];
  xt_p1[1][7] = vxy_all[22];
  xt_p1[2][7] = vxy_all[23];

  for (k = 0; k < STEP5_NNODE_1D; ++k) {
    /* z-direction in reference coordinate {t} */
    xk = xint[k];

    for (j = 0; j < STEP5_NNODE_1D; ++j) {
      /* y-direction in reference coordinate {s} */
      xj = xint[j];

      for (i = 0; i < STEP5_NNODE_1D; ++i) {
        /* x-direction in reference coordinate {r} */
        xi = xint[i];

        /* p1 shape functions are evaluated at GLL points of the desired
         * polynomial degree */
        step5_shape_p1 (xi, xj, xk, fn);

        for (ii = 0; ii < P4EST_DIM; ++ii) {
          xy_computed[ii] = 0.0;
          for (jj = 0; jj < P4EST_CHILDREN; ++jj) {
            xy_computed[ii] += xt_p1[ii][jj] * fn[jj];
          }
        }

        /* construct geometry for a high-order element:
         * find node that is being visited */
        visited_node = i + j * STEP5_NNODE_1D +
          k * STEP5_NNODE_1D * STEP5_NNODE_1D;
        for (i_dim = 0; i_dim < xt_nrows_dim; i_dim++) {
          xt[visited_node + i_dim * xt_ncols_nnode] = xy_computed[i_dim];
        }
      }
    }
  }
#endif
}

static void
step5_init_initial_condition (p4est_t * p4est,
                              p4est_topidx_t which_tree, p4est_quadrant_t * q)
{
  step5_data_t       *data = (step5_data_t *) q->p.user_data;
  double              vxy_all[P4EST_DIM * P4EST_CHILDREN];
  int                 i;
  double              x, y;
#ifdef P4_TO_P8
  double              z;
#endif

  /* generate coordinates */
  step5_qcoord_to_vertex_all (p4est, which_tree, q, vxy_all);
  step5_construct_GLL_xyz (vxy_all, P4EST_DIM, STEP5_NNODE, data->xt);

  for (i = 0; i < STEP5_NNODE; ++i) {
    x = 10 * data->xt[i];
    y = 10 * data->xt[STEP5_NNODE + i];
    data->data[i] = sin (x) + cos (y);
#ifdef P4_TO_P8
    z = 10 * data->xt[2 * STEP5_NNODE + i];
    data->data[i] += sqrt (z);
#endif
  }
}

static void
step5_collect_info (p4est_iter_volume_info_t * info, void *user_data)
{
  int                 n, idx;
  double             *this_o_ptr;
#ifndef STEP5_HO
  int                 npoints = P4EST_CHILDREN;
  int                 corner_list[npoints];     /* only need corners */
#else
  int                 npoints = STEP5_NNODE;
  int                 ndatas = 1;       /* 1 piece of data to extract */
  int                 d;
  p4est_locidx_t      numquads;
#endif

  /* we passed the array of values to fill as
   * the user_data in the call to p4est_iterate */
  sc_array_t         *output = (sc_array_t *) user_data;
  p4est_t            *p4est = info->p4est;
  p4est_quadrant_t   *q = info->quad;
  p4est_topidx_t      which_tree = info->treeid;

  /* this is the index of q *within its tree's numbering*.
   * We want to convert its index for all the
   * quadrants on this process, which we do below */
  p4est_locidx_t      local_id = info->quadid;
  p4est_tree_t       *tree;
  step5_data_t       *data = (step5_data_t *) q->p.user_data;

  tree = p4est_tree_array_index (p4est->trees, which_tree);
  local_id += tree->quadrants_offset;   /* make relative to the MPI process */

#ifndef STEP5_HO
#ifndef P4_TO_P8
  corner_list[0] = 0;
  corner_list[1] = STEP5_NNODE_1D - 1;
  corner_list[2] = STEP5_NNODE_1D * (STEP5_NNODE_1D - 1);
  corner_list[3] = STEP5_NNODE_1D * STEP5_NNODE_1D - 1;
#else
  corner_list[0] = 0;
  corner_list[1] = STEP5_NNODE_1D - 1;
  corner_list[2] = STEP5_NNODE_1D * (STEP5_NNODE_1D - 1);
  corner_list[3] = STEP5_NNODE_1D * STEP5_NNODE_1D - 1;
  corner_list[4] = STEP5_NNODE_1D * STEP5_NNODE_1D * (STEP5_NNODE_1D - 1);
  corner_list[5] = STEP5_NNODE_1D * STEP5_NNODE_1D * (STEP5_NNODE_1D - 1) +
    STEP5_NNODE_1D - 1;
  corner_list[6] = STEP5_NNODE_1D * STEP5_NNODE_1D * STEP5_NNODE_1D - 1 -
    (STEP5_NNODE_1D - 1);
  corner_list[7] = STEP5_NNODE_1D * STEP5_NNODE_1D * STEP5_NNODE_1D - 1;
#endif
#endif

  /* extracting data */
  for (n = 0; n < npoints; ++n) {
    this_o_ptr = (double *) sc_array_index (output, n + npoints * local_id);
#ifndef STEP5_HO
    idx = corner_list[n];
#else
    idx = n;
#endif
    *this_o_ptr = data->data[idx];
  }

#ifdef STEP5_HO
  numquads = p4est->local_num_quadrants;
  /* copying coordinates of each node
   * ordering: [ x_0, y_0, z_0 ... x_n, y_n, z_n ] */
  for (n = 0; n < STEP5_NNODE; ++n) {
    for (d = 0; d < P4EST_DIM; ++d) {
      this_o_ptr =
        (double *) sc_array_index (output,
                                   n * P4EST_DIM + d +
                                   STEP5_NNODE * ((local_id * P4EST_DIM) +
                                                  ndatas * numquads));
      *this_o_ptr = data->xt[n + d * STEP5_NNODE];
    }
  }
#endif
}

int
main (int argc, char **argv)
{
  int                 mpiret, retval;
  int                 numquads, array_size;
  int                 ndatas = 1;
#ifndef STEP5_HO
  int                 npoints = P4EST_CHILDREN;
  int                 vals_per_node = ndatas;
#else
  int                 npoints = STEP5_NNODE;
  int                 vals_per_node = ndatas + P4EST_DIM;
  sc_array_t         *positions;
#endif
  sc_array_t         *data, *output;
  sc_MPI_Comm         mpicomm;
  p4est_t            *p4est;
  p4est_connectivity_t *conn;
  p4est_vtk_context_t *vtk_context;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  mpicomm = sc_MPI_COMM_WORLD;

  sc_init (mpicomm, 1, 1, NULL, SC_LP_ESSENTIAL);
  p4est_init (NULL, SC_LP_PRODUCTION);
  P4EST_GLOBAL_PRODUCTIONF
    ("This is the p4est %dD demo example/steps/%s_step5\n",
     P4EST_DIM, P4EST_STRING);

#ifndef P4_TO_P8
  conn = p4est_connectivity_new_unitsquare ();
#else
  conn = p8est_connectivity_new_unitcube ();
#endif

  /* 4 cells in 2x2 square (in 2D), each point with coords and data values */
  p4est = p4est_new_ext (mpicomm,       /* communicator */
                         conn,  /* connectivity */
                         0,     /* minimum quadrants per MPI process */
                         1,     /* minimum level of refinement */
                         1,     /* fill uniform */
                         sizeof (step5_data_t), /* data size */
                         step5_init_initial_condition,  /* initializes data */
                         NULL); /* context */

  numquads = p4est->local_num_quadrants;

  output = sc_array_new_count (sizeof (double),
                               numquads * vals_per_node * npoints);

  /* *INDENT-OFF* */
#ifndef P4_TO_P8
  p4est_iterate(p4est,
                NULL,
                (void *)output,
                step5_collect_info,
                NULL,
                NULL);
#else
  p4est_iterate(p4est,
                NULL,
                (void *)output,
                step5_collect_info,
                NULL,
                NULL,
                NULL);
#endif
  /* *INDENT-ON* */

  array_size = numquads * npoints;

  /* begin writing the output files */
  vtk_context = p4est_vtk_context_new (p4est, P4EST_STRING "_step5");

#ifdef STEP5_HO
  positions = sc_array_new_size (sizeof (double), P4EST_DIM * array_size);
  sc_array_move_part (positions,
                      0,
                      output,
                      (size_t) (ndatas * array_size),
                      (size_t) (P4EST_DIM * array_size));

  /* In this example we do not shrink down each individual quadrant */
  p4est_vtk_context_set_scale (vtk_context, 1.);

  /* The continuous setting only effects the low-order output */
  p4est_vtk_context_set_continuous (vtk_context, 1);

  vtk_context = p4est_vtk_write_header_ho (vtk_context,
                                           positions, STEP5_NNODE_1D);
#else
  vtk_context = p4est_vtk_write_header (vtk_context);
#endif
  SC_CHECK_ABORT (vtk_context != NULL,
                  P4EST_STRING "_vtk: Error writing header");

  data = sc_array_new_size (sizeof (double), array_size);
  sc_array_move_part (data, 0, output, 0, (size_t) array_size);

  vtk_context = p4est_vtk_write_point_dataf (vtk_context,
                                             1, 0, "data", data, vtk_context);
  SC_CHECK_ABORT (vtk_context != NULL,
                  P4EST_STRING "_vtk: Error writing point data");

  /* write_footer also calls vtk_context_destroy */
  retval = p4est_vtk_write_footer (vtk_context);
  SC_CHECK_ABORT (!retval, P4EST_STRING "_vtk: Error writing footer");

#ifdef STEP5_HO
  sc_array_destroy (positions);
#endif
  sc_array_destroy (output);
  sc_array_destroy (data);

  p4est_destroy (p4est);
  p4est_connectivity_destroy (conn);

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;
}
