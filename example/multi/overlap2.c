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

/*
 * Usage: p4est_overlap
 *
 * Create two forest workflow apps in the same main program.
 * One app requires data from the other, which we search in parallel.
 * In this example, both apps use a duplicate of MPI_COMM_WORLD.
 * Thus, they execute their respective program alternatingly.
 *
 * The two apps are labeled producer and consumer.
 *
 * In 2D, the producer map is extruded in the third dimension by 1.
 * In 3D, producer and consumer maps are fully 3D.
 * All coordinate transforms are affine-linear.
 */

#include <sc_notify.h>
#include <sc_options.h>
#include <sc_flops.h>
#ifndef P4_TO_P8
#include <p4est_extended.h>
#include <p4est_search.h>
#include <p4est_vtk.h>
#include <p4est_bits.h>
#else
#include <p8est_extended.h>
#include <p8est_search.h>
#include <p8est_vtk.h>
#include <p8est_bits.h>
#endif

#define P4EST_CON_TOLERANCE SC_1000_EPS
#define P4EST_PRO_TOLERANCE (2 * SC_1000_EPS)

enum
{
  OVERLAP_EXCHANGE,
  OVERLAP_SEARCH_PARTITION,
#ifdef P4EST_ENABLE_MPI
  OVERLAP_NOTIFY,
  OVERLAP_POST_MESSAGES,
  OVERLAP_INTERPOLATE,
  OVERLAP_UPDATE_QUERY_POINTS,
  OVERLAP_WAITALL,
#endif
  OVERLAP_UPDATE_LOCAL,
  OVERLAP_FREE_COMMUNICATION_DATA,
  OVERLAP_NUM_LOCAL_CONS_QUADRANTS,
  OVERLAP_NUM_LOCAL_PROD_QUADRANTS,
  OVERLAP_NUM_QP_SENT,
  OVERLAP_NUM_QP_RECEIVED,
  OVERLAP_NUM_QP_SENTRECVD,
  OVERLAP_NUM_PROD_SEARCH_OPS,
  OVERLAP_NUM_CONS_SEARCH_OPS,
  OVERLAP_NUM_SEARCH_OPS,
  OVERLAP_NUM_STATS
};

static int          overlap_stats_type[OVERLAP_NUM_STATS] = { 0, 0,
#ifdef P4EST_ENABLE_MPI
  0, 0, 0, 0, 0,
#endif
  0, 0, 1, 1, 1, 1, 1, 1, 1, 1
};

typedef struct overlap_tstats
{
  sc_flopinfo_t       fi;
  sc_statinfo_t       stats[OVERLAP_NUM_STATS];
  size_t              lnum_qp_sent;
  size_t              lnum_qp_recvd;
}
overlap_tstats_t;

typedef struct overlap_prodata
{
  double              myvalue;
  int                 isset;
}
overlap_prodata_t;

typedef struct overlap_point
{
  int                 rank;
  int                 isboundary;
  p4est_locidx_t      lnum;
  double              xyz[3];
  overlap_prodata_t   prodata;
}
overlap_point_t;

typedef struct overlap_send_buf
{
  int                 rank;
  sc_array_t          ops;
}
overlap_send_buf_t;

typedef overlap_send_buf_t overlap_recv_buf_t;

typedef enum comm_tag
{
  COMM_TAG_CONSDATA = P4EST_COMM_TAG_LAST,
  COMM_TAG_PRODATA
}
comm_tag_t;

typedef void        (*overlap_invmap_t) (p4est_connectivity_t *conn,
                                         p4est_topidx_t which_tree,
                                         const double abc[3], double xyz[3]);

typedef struct overlap_producer
{
  /* mesh constituents */
  sc_MPI_Comm         procomm;
  p4est_connectivity_t *proconn;
  p4est_t            *pro4est;
  p4est_geometry_t   *progeom, producer_geometry;
  overlap_invmap_t    invmap;

  /* parameters */
  int                 pminl;

  /* work variables */
  p4est_locidx_t      lquad_idx;

  /* communication */
  sc_MPI_Comm         glocomm;
  int                 prorank;
  int                 iprorank;
  sc_array_t         *recv_buffer;
  sc_array_t         *recv_reqs;
  sc_array_t         *send_reqs;

  /* vtk cell data */
  sc_array_t         *interpolation_data;
  sc_array_t         *xyz_data;

  /* adaptive refinement */
  int                 refining;

  /* timings */
  overlap_tstats_t   *tstats;
}
overlap_producer_t;

typedef struct overlap_condata
{
  int                 refine;
}
overlap_condata_t;

typedef struct overlap_consumer
{
  /* mesh constituents */
  sc_MPI_Comm         concomm;
  p4est_connectivity_t *conconn;
  p4est_t            *con4est;
  p4est_geometry_t   *congeom, consumer_geometry;
  overlap_invmap_t    invmap;

  /* parameters */
  int                 cminl;

  /* work variables */
  p4est_locidx_t      lquad_idx;
  sc_array_t         *query_xyz;

  /* minimal knowledge of the producer's mesh */
  p4est_connectivity_t *producer_conn;
  const p4est_gloidx_t *producer_gfq;
  const p4est_quadrant_t *producer_gfp;
  int                 pronum_procs;
  p4est_topidx_t      pronum_trees;

  /* communication */
  sc_MPI_Comm         glocomm;
  int                 conrank;
  int                 iconrank;
  sc_array_t         *send_buffer;
  sc_array_t         *send_reqs;
  sc_array_t         *recv_buffer;
  sc_array_t         *recv_reqs;

  /* vtk cell data */
  sc_array_t         *interpolation_data;
  sc_array_t         *isset_data;
  sc_array_t         *xyz_data;

  /* timings */
  overlap_tstats_t   *tstats;
}
overlap_consumer_t;

typedef struct overlap_global
{
  sc_MPI_Comm         glocomm;
  int                 rounds;
  int                 refinement_method;
  int                 example;
  overlap_tstats_t   *tstats;
  overlap_producer_t  pro, *p;
  overlap_consumer_t  con, *c;
}
overlap_global_t;

#define OVERLAP_IROOTLEN (1. / P4EST_ROOT_LEN)

/** Turn statistics collected so far into one value */
static void
sc_stats_collapse (sc_statinfo_t *stats)
{
  double              value;

  SC_ASSERT (stats->dirty);
  value = stats->sum_values;
  stats->sum_values = value;
  stats->sum_squares = value * value;
  stats->min = stats->max = value;
  stats->count = 1;
}

void
sc_stats_print_x (int package_id, int log_priority, int nvars,
                  sc_statinfo_t *stats, int *stats_type, int full,
                  int summary)
{
  int                 i, ti, count;
  sc_statinfo_t      *si;
  char                buffer[BUFSIZ];

  if (full) {
    for (i = 0; i < nvars; ++i) {
      si = &stats[i];
      ti = stats_type[i];
      /* begin printing */
      if (ti) {                 /* the stat is integer */
        if (si->variable != NULL) {
          SC_GEN_LOGF (package_id, SC_LC_GLOBAL, log_priority,
                       "Statistics for   %s\n", si->variable);
        }
        else {
          SC_GEN_LOGF (package_id, SC_LC_GLOBAL, log_priority,
                       "Statistics for %d\n", i);
        }
        SC_GEN_LOGF (package_id, SC_LC_GLOBAL, log_priority,
                     "   Global number of values: %7ld\n", si->count);
        if (!si->count) {
          continue;
        }
        if (si->average != 0.) {        /* ignore the comparison warning */
          SC_GEN_LOGF (package_id, SC_LC_GLOBAL, log_priority,
                       "   Mean value (std. dev.):           %.0f (%.3g = %.3g%%)\n",
                       si->average, si->standev,
                       100. * si->standev / fabs (si->average));
        }
        else {
          SC_GEN_LOGF (package_id, SC_LC_GLOBAL, log_priority,
                       "   Mean value (std. dev.):           %.0f (%.3g)\n",
                       si->average, si->standev);
        }
        SC_GEN_LOGF (package_id, SC_LC_GLOBAL, log_priority,
                     "   Minimum attained at rank %7d: %.0f\n",
                     si->min_at_rank, si->min);
        SC_GEN_LOGF (package_id, SC_LC_GLOBAL, log_priority,
                     "   Maximum attained at rank %7d: %.0f\n",
                     si->max_at_rank, si->max);
      }
      else {
        if (si->variable != NULL) {
          SC_GEN_LOGF (package_id, SC_LC_GLOBAL, log_priority,
                       "Statistics for   %s\n", si->variable);
        }
        else {
          SC_GEN_LOGF (package_id, SC_LC_GLOBAL, log_priority,
                       "Statistics for %d\n", i);
        }
        SC_GEN_LOGF (package_id, SC_LC_GLOBAL, log_priority,
                     "   Global number of values: %7ld\n", si->count);
        if (!si->count) {
          continue;
        }
        if (si->average != 0.) {        /* ignore the comparison warning */
          SC_GEN_LOGF (package_id, SC_LC_GLOBAL, log_priority,
                       "   Mean value (std. dev.):           %g (%.3g = %.3g%%)\n",
                       si->average, si->standev,
                       100. * si->standev / fabs (si->average));
        }
        else {
          SC_GEN_LOGF (package_id, SC_LC_GLOBAL, log_priority,
                       "   Mean value (std. dev.):           %g (%.3g)\n",
                       si->average, si->standev);
        }
        SC_GEN_LOGF (package_id, SC_LC_GLOBAL, log_priority,
                     "   Minimum attained at rank %7d: %g\n",
                     si->min_at_rank, si->min);
        SC_GEN_LOGF (package_id, SC_LC_GLOBAL, log_priority,
                     "   Maximum attained at rank %7d: %g\n",
                     si->max_at_rank, si->max);
      }
    }
  }
  else {
    for (i = 0; i < nvars; ++i) {
      si = &stats[i];
      ti = stats_type[i];
      if (ti) {
        /* print just the average */
        if (si->variable != NULL) {
          snprintf (buffer, BUFSIZ, "for %s:", si->variable);
        }
        else {
          snprintf (buffer, BUFSIZ, "for %3d:", i);
        }
        if (si->average != 0.) {        /* ignore the comparison warning */
          SC_GEN_LOGF (package_id, SC_LC_GLOBAL, log_priority,
                       "Mean (sigma) %-23s %.0f (%.3g = %.3g%%)\n",
                       buffer, si->average, si->standev,
                       100. * si->standev / fabs (si->average));
        }
        else {
          SC_GEN_LOGF (package_id, SC_LC_GLOBAL, log_priority,
                       "Mean (sigma) %-23s %.0f (%.3g)\n", buffer,
                       si->average, si->standev);
        }
      }
      else {
        /* print just the average */
        if (si->variable != NULL) {
          snprintf (buffer, BUFSIZ, "for %s:", si->variable);
        }
        else {
          snprintf (buffer, BUFSIZ, "for %3d:", i);
        }
        if (si->average != 0.) {        /* ignore the comparison warning */
          SC_GEN_LOGF (package_id, SC_LC_GLOBAL, log_priority,
                       "Mean (sigma) %-23s %g (%.3g = %.3g%%)\n",
                       buffer, si->average, si->standev,
                       100. * si->standev / fabs (si->average));
        }
        else {
          SC_GEN_LOGF (package_id, SC_LC_GLOBAL, log_priority,
                       "Mean (sigma) %-23s %g (%.3g)\n", buffer,
                       si->average, si->standev);
        }
      }
    }
  }

  /* the summary always contains all variables */
  if (summary) {
    count = snprintf (buffer, BUFSIZ, "Summary = ");
    for (i = 0; i < nvars && count >= 0 && (size_t) count < BUFSIZ; ++i) {
      si = &stats[i];
      ti = stats_type[i];
      if (ti) {
        count += snprintf (buffer + count, BUFSIZ - count, "%s%.0f",
                           i == 0 ? "[ " : " ", si->average);
      }
      else {
        count += snprintf (buffer + count, BUFSIZ - count, "%s%g",
                           i == 0 ? "[ " : " ", si->average);
      }
    }
    if (count >= 0 && (size_t) count < BUFSIZ) {
      snprintf (buffer + count, BUFSIZ - count, "%s", " ];\n");
      SC_GEN_LOG (package_id, SC_LC_GLOBAL, log_priority, buffer);
    }
    else {
      SC_GEN_LOG (package_id, SC_LC_GLOBAL, log_priority,
                  "Summary overflow\n");
    }
    count = snprintf (buffer, BUFSIZ, "Maximum = ");
    for (i = 0; i < nvars && count >= 0 && (size_t) count < BUFSIZ; ++i) {
      si = &stats[i];
      ti = stats_type[i];
      if (ti) {
        count += snprintf (buffer + count, BUFSIZ - count, "%s%.0f",
                           i == 0 ? "[ " : " ", si->max);
      }
      else {
        count += snprintf (buffer + count, BUFSIZ - count, "%s%g",
                           i == 0 ? "[ " : " ", si->max);
      }
    }
    if (count >= 0 && (size_t) count < BUFSIZ) {
      snprintf (buffer + count, BUFSIZ - count, "%s", " ];\n");
      SC_GEN_LOG (package_id, SC_LC_GLOBAL, log_priority, buffer);
    }
    else {
      SC_GEN_LOG (package_id, SC_LC_GLOBAL, log_priority,
                  "Maximum overflow\n");
    }
  }
}

static double
overlap_producer_evaluate (overlap_producer_t *p, double pxyz[3])
{
  double              r[3];

  P4EST_ASSERT (p != NULL);
  P4EST_ASSERT (pxyz != NULL);

  r[0] = (pxyz[0] - .4) / 1.6;
  r[1] = (pxyz[1] + .3) / 1.1;
  r[2] = (pxyz[2] - .2) / 0.5;
  return
    .1 + .9 * exp (-.5 * (SC_SQR (r[0]) + SC_SQR (r[1]) + SC_SQR (r[2])));
}

#ifdef OVERLAP_WITH_CUBE_MAP    /* cube map not currently used */
#ifdef P4EST_ENABLE_DEBUG

static void
overlap_cube_invmap (p4est_connectivity_t *conn, p4est_topidx_t which_tree,
                     const double xyz[3], double abc[3]);

#endif /* P4EST_ENABLE_DEBUG */

static void
overlap_cube_map (p4est_geometry_t *geom, p4est_topidx_t which_tree,
                  const double abc[3], double xyz[3])
{
  double              a, co, si, x;
#ifdef P4EST_ENABLE_DEBUG
  double              def[3];
#endif
  double             *vert;
  p4est_topidx_t      vind;
  p4est_connectivity_t *conn;

  P4EST_ASSERT (geom != NULL);
  P4EST_ASSERT (geom->X == overlap_cube_map);

  /* preimage domain is [0, 2]^3 */

  /* access origin vertex of given tree */
  P4EST_ASSERT (geom != NULL);
  P4EST_ASSERT (geom->X == overlap_cube_map);
  conn = (p4est_connectivity_t *) geom->user;
  P4EST_ASSERT (conn != NULL && conn->vertices != NULL);
  vind = conn->tree_to_vertex[P4EST_CHILDREN * which_tree + 0];
  vert = &conn->vertices[3 * vind + 0];

  /* move brick back towards origin and scale slightly */
  xyz[0] = (vert[0] + abc[0] - .5) * 1.1;
  xyz[1] = (vert[1] + abc[1] - .5) * 1.2;
  xyz[2] = (vert[2] + abc[2] - .5) * 1.3;

  /* rotate 20 degrees around the y axis */
  a = 20. * M_PI / 180.;
  co = cos (a);
  si = sin (a);
  x = xyz[0];
  xyz[0] = co * x - si * xyz[2];
  xyz[2] = si * x + co * xyz[2];

#ifdef P4EST_ENABLE_DEBUG
  /* verify inverse mapping */
  overlap_cube_invmap (conn, which_tree, xyz, def);
  P4EST_ASSERT (fabs (abc[0] - def[0]) < SC_1000_EPS &&
                fabs (abc[1] - def[1]) < SC_1000_EPS &&
                fabs (abc[2] - def[2]) < SC_1000_EPS);
#endif
}

static void
overlap_cube_invmap (p4est_connectivity_t *conn, p4est_topidx_t which_tree,
                     const double xyz[3], double abc[3])
{
  double              a, co, si;
  double             *vert;
  p4est_topidx_t      vind;

  P4EST_ASSERT (conn != NULL && conn->vertices != NULL);
  P4EST_ASSERT (0 <= which_tree && which_tree < conn->num_trees);
  vind = conn->tree_to_vertex[P4EST_CHILDREN * which_tree + 0];
  vert = &conn->vertices[3 * vind + 0];

  /* rotate 20 degrees backwards around the y axis */
  a = -20. * M_PI / 180.;
  co = cos (a);
  si = sin (a);
  abc[0] = co * xyz[0] - si * xyz[2];
  abc[2] = si * xyz[0] + co * xyz[2];
  abc[1] = xyz[1];

  /* scale slightly and move away from origin into brick */
  abc[0] = abc[0] / 1.1 + .5 - vert[0];
  abc[1] = abc[1] / 1.2 + .5 - vert[1];
  abc[2] = abc[2] / 1.3 + .5 - vert[2];
}

#endif /* 0 */

static void
overlap_curved_invmap
  (p4est_connectivity_t *conn, p4est_topidx_t which_tree,
   const double xyz[3], double abc[3]);

static int          xbricks = 10;
static int          ybricks = 10;

static void
overlap_curved_map (p4est_geometry_t *geom, p4est_topidx_t which_tree,
                    const double abc[3], double xyz[3])
{
  double             *vert;
  p4est_topidx_t      vind;
#ifdef P4EST_ENABLE_DEBUG
  double              def[3];
#endif
  p4est_connectivity_t *conn;

  /* map to [0,xbricks]x[0,ybricks]x[0,2] quadrant */
  P4EST_ASSERT (geom != NULL);
  P4EST_ASSERT (geom->X == overlap_curved_map);
  conn = (p4est_connectivity_t *) geom->user;
  P4EST_ASSERT (conn != NULL && conn->vertices != NULL);
  vind = conn->tree_to_vertex[P4EST_CHILDREN * which_tree + 0];
  vert = &conn->vertices[3 * vind + 0];

  xyz[0] = vert[0] + abc[0];
  xyz[1] = vert[1] + abc[1];
  xyz[2] = vert[2] + abc[2];

  /* scale down x and y to avoid octants getting to elongated */
  xyz[0] *= 2. / (double) xbricks;
  xyz[1] *= 2. / (double) ybricks;

  /* shift y and z and according to a curve of x */
  xyz[1] += 0.2;
  xyz[1] *= 1. / 4.;
  xyz[1] += (0.75 - xyz[0]) * (0.75 - xyz[0]);

  xyz[2] -= 0.2;
  xyz[2] *= 1. / 4.;
  xyz[2] += (0.75 - xyz[0]) * (0.75 - xyz[0]);

#ifdef P4EST_ENABLE_DEBUG
  overlap_curved_invmap (conn, which_tree, xyz, def);
  P4EST_ASSERT (fabs (abc[0] - def[0]) < SC_1000_EPS &&
                fabs (abc[1] - def[1]) < SC_1000_EPS &&
                fabs (abc[2] - def[2]) < SC_1000_EPS);
#endif
}

static void
overlap_curved_invmap (p4est_connectivity_t *conn, p4est_topidx_t which_tree,
                       const double xyz[3], double abc[3])
{
  double             *vert;
  p4est_topidx_t      vind;

  abc[0] = xyz[0];
  abc[1] = xyz[1];
  abc[2] = xyz[2];

  /* invert shifting of y and z according to a curve of x */
  abc[1] -= (0.75 - abc[0]) * (0.75 - abc[0]);
  abc[1] *= 4.;
  abc[1] -= 0.2;

  abc[2] -= (0.75 - abc[0]) * (0.75 - abc[0]);
  abc[2] *= 4.;
  abc[2] += 0.2;

  /* invert scaling of x and y */
  abc[0] *= (double) xbricks / 2.;
  abc[1] *= (double) ybricks / 2.;

  /* invert mapping to [0,xbricks]x[0,ybricks]x[0,2] quadrant */
  P4EST_ASSERT (conn != NULL && conn->vertices != NULL);
  P4EST_ASSERT (0 <= which_tree && which_tree < conn->num_trees);
  vind = conn->tree_to_vertex[P4EST_CHILDREN * which_tree + 0];
  vert = &conn->vertices[3 * vind + 0];

  abc[0] = abc[0] - vert[0];
  abc[1] = abc[1] - vert[1];
  abc[2] = abc[2] - vert[2];
}

static void
overlap_brick_invmap
  (p4est_connectivity_t *conn, p4est_topidx_t which_tree,
   const double xyz[3], double abc[3]);

static void
overlap_brick_map (p4est_geometry_t *geom, p4est_topidx_t which_tree,
                   const double abc[3], double xyz[3])
{
  double              a, co, si, x;
#ifdef P4EST_ENABLE_DEBUG
  double              def[3];
#endif
  double             *vert;
  p4est_topidx_t      vind;
  p4est_connectivity_t *conn;

  /* preimage domain is 3x2x1 with origin in the lower left front */

  /* access origin vertex of given tree */
  P4EST_ASSERT (geom != NULL);
  P4EST_ASSERT (geom->X == overlap_brick_map);
  conn = (p4est_connectivity_t *) geom->user;
  P4EST_ASSERT (conn != NULL && conn->vertices != NULL);
  vind = conn->tree_to_vertex[P4EST_CHILDREN * which_tree + 0];
  vert = &conn->vertices[3 * vind + 0];

  /* center brick around origin and scale */
  xyz[0] = (vert[0] + abc[0] - 1.5) * 0.7;
  xyz[1] = (vert[1] + abc[1] - 1.0) * 0.6;
  xyz[2] = (vert[2] + abc[2] - 0.5) * 0.5;

  /* rotate 30 degrees around the z axis and shift */
  a = 30. * M_PI / 180.;
  co = cos (a);
  si = sin (a);
  x = xyz[0];
  xyz[0] = co * x - si * xyz[1] + 0.8;
  xyz[1] = si * x + co * xyz[1] + 0.8;
  xyz[2] += 0.5;

#ifdef P4EST_ENABLE_DEBUG
  overlap_brick_invmap (conn, which_tree, xyz, def);
  P4EST_ASSERT (fabs (abc[0] - def[0]) < SC_1000_EPS &&
                fabs (abc[1] - def[1]) < SC_1000_EPS &&
                fabs (abc[2] - def[2]) < SC_1000_EPS);
#endif
}

static void
overlap_brick_invmap (p4est_connectivity_t *conn, p4est_topidx_t which_tree,
                      const double xyz[3], double abc[3])
{
  double              a, co, si, x;
  double             *vert;
  p4est_topidx_t      vind;

  P4EST_ASSERT (conn != NULL && conn->vertices != NULL);
  P4EST_ASSERT (0 <= which_tree && which_tree < conn->num_trees);
  vind = conn->tree_to_vertex[P4EST_CHILDREN * which_tree + 0];
  vert = &conn->vertices[3 * vind + 0];

  /* shift back */
  abc[0] = xyz[0] - 0.8;
  abc[1] = xyz[1] - 0.8;
  abc[2] = xyz[2] - 0.5;

  /* complete 360 degree rotation */
  a = 330. * M_PI / 180.;
  co = cos (a);
  si = sin (a);
  x = abc[0];
  abc[0] = co * x - si * abc[1];
  abc[1] = si * x + co * abc[1];

  /* invert scaling and centering */
  abc[0] = abc[0] / 0.7 + 1.5 - vert[0];
  abc[1] = abc[1] / 0.6 + 1.0 - vert[1];
  abc[2] = abc[2] / 0.5 + 0.5 - vert[2];
}

static void
get_quadrant_center (p4est_quadrant_t *q, double qxyz[3])
{
  p4est_qcoord_t      h2;

  /* get quadrant center reference coordinates and store them in qxyz */
  h2 = P4EST_QUADRANT_LEN (q->level) >> 1;
  qxyz[0] = OVERLAP_IROOTLEN * (q->x + h2);
  qxyz[1] = OVERLAP_IROOTLEN * (q->y + h2);
#ifndef P4_TO_P8
  qxyz[2] = 0.;
#else
  qxyz[2] = OVERLAP_IROOTLEN * (q->z + h2);
#endif
}

static void
overlap_producer_compute (p4est_iter_volume_info_t *info, void *user_data)
{
  p4est_quadrant_t   *q;
  overlap_prodata_t  *d;
  overlap_producer_t *p;
  double              qxyz[3], phys[3];
  p4est_locidx_t      lnum;

  P4EST_ASSERT (info != NULL && info->p4est != NULL);
  P4EST_ASSERT (info->p4est->user_pointer == user_data);

  p = (overlap_producer_t *) info->p4est->user_pointer;
  P4EST_ASSERT (p->pro4est == info->p4est);
  P4EST_ASSERT (info->quad != NULL);
  q = info->quad;
  d = (overlap_prodata_t *) q->p.user_data;
  P4EST_ASSERT (d != NULL);

  /* transform consumer quadrant center to physical using map */
  get_quadrant_center (q, qxyz);
  p->progeom->X (p->progeom, info->treeid, qxyz, phys);

  /* interpolate prescribed field at that point */
  d->myvalue = overlap_producer_evaluate (p, phys);
  d->isset = 1;

  /* store vtk cell data */
  lnum = p->lquad_idx++;
  *(double *) sc_array_index (p->interpolation_data, lnum) = d->myvalue;
  *(double *) sc_array_index (p->xyz_data, 3 * lnum) = phys[0];
  *(double *) sc_array_index (p->xyz_data, 3 * lnum + 1) = phys[1];
  *(double *) sc_array_index (p->xyz_data, 3 * lnum + 2) = phys[2];

  P4EST_LDEBUGF ("Producer input tree %d level %d quad %g %g %g\n",
                 (int) info->treeid, q->level, qxyz[0], qxyz[1], qxyz[2]);
  P4EST_LDEBUGF ("Producer compute %g %g %g\n", phys[0], phys[1], phys[2]);
}

static void
overlap_consumer_compute_center (p4est_iter_volume_info_t *info,
                                 void *user_data)
{
  p4est_quadrant_t   *q;
  overlap_consumer_t *c;
  overlap_point_t    *op;
  double              qxyz[3], *phys;

  P4EST_ASSERT (info != NULL && info->p4est != NULL);
  P4EST_ASSERT (info->p4est->user_pointer == user_data);

  /* access quadrant */
  c = (overlap_consumer_t *) info->p4est->user_pointer;
  P4EST_ASSERT (c->con4est == info->p4est);
  P4EST_ASSERT (c->lquad_idx >= 0 &&
                c->lquad_idx < c->con4est->local_num_quadrants);
  P4EST_ASSERT (info->quad != NULL);
  q = info->quad;

  /* transform consumer quadrant center to physical using map */
  get_quadrant_center (q, qxyz);
  phys = (op = (overlap_point_t *)
          sc_array_index (c->query_xyz, (size_t) c->lquad_idx))->xyz;
  memset (op, 0, sizeof (overlap_point_t));
  c->congeom->X (c->congeom, info->treeid, qxyz, phys);
  op->lnum = c->lquad_idx++;
  op->rank = -1;
  op->isboundary = -1;          /* our interpolation scheme ignores boundary info */
  op->prodata.myvalue = 0.;
  op->prodata.isset = 0;

  /* store quadrant center coordinates for vtk output */
  *(double *) sc_array_index (c->xyz_data, 3 * op->lnum) = op->xyz[0];
  *(double *) sc_array_index (c->xyz_data, 3 * op->lnum + 1) = op->xyz[1];
  *(double *) sc_array_index (c->xyz_data, 3 * op->lnum + 2) = op->xyz[2];

  P4EST_LDEBUGF ("Consumer input tree %d level %d quad %g %g %g\n",
                 (int) info->treeid, q->level, qxyz[0], qxyz[1], qxyz[2]);
  P4EST_LDEBUGF ("Consumer point %ld compute %g %g %g\n",
                 (long) op->lnum, phys[0], phys[1], phys[2]);

  /* optimize: ignore this point if not intersecting producer domain */
}

static void
overlap_consumer_compute_corners (p4est_iter_volume_info_t *info,
                                  void *user_data)
{
  p4est_quadrant_t   *q;
  overlap_consumer_t *c;
  overlap_point_t    *op;
  p4est_qcoord_t      h, qcoords[3];
  double              qxyz[3], *phys;
  int                 i, dim, lu;

  P4EST_ASSERT (info != NULL && info->p4est != NULL);
  P4EST_ASSERT (info->p4est->user_pointer == user_data);

  /* access quadrant */
  c = (overlap_consumer_t *) info->p4est->user_pointer;
  P4EST_ASSERT (c->con4est == info->p4est);
  P4EST_ASSERT (c->lquad_idx >= 0 &&
                c->lquad_idx <
                P4EST_CHILDREN * c->con4est->local_num_quadrants);
  P4EST_ASSERT (info->quad != NULL);
  q = info->quad;

  /* iterate over all children */
  h = P4EST_QUADRANT_LEN (q->level);
  for (i = 0; i < P4EST_CHILDREN; i++) {
    /* compute reference coordinates of corner */
    qcoords[0] = q->x + ((i % 2) ? h : 0);
    qcoords[1] = q->y + (((i % 4) / 2) ? h : 0);
#ifndef P4_TO_P8
    qcoords[2] = 0;
#else
    qcoords[2] = q->z + ((i / 4) ? h : 0);
#endif

    qxyz[0] = OVERLAP_IROOTLEN * qcoords[0];
    qxyz[1] = OVERLAP_IROOTLEN * qcoords[1];
    qxyz[2] = OVERLAP_IROOTLEN * qcoords[2];
    phys = (op = (overlap_point_t *)
            sc_array_index (c->query_xyz, (size_t) c->lquad_idx))->xyz;
    memset (op, 0, sizeof (overlap_point_t));
    c->congeom->X (c->congeom, info->treeid, qxyz, phys);

    op->lnum = c->lquad_idx++;
    op->rank = -1;
    op->prodata.myvalue = 0.;
    op->prodata.isset = 0;

    /* determine if we are at the forest boundary or not */
    op->isboundary = 0;
    /* loop over all faces by looping over all dimensions and lower/upper */
    for (dim = 0; dim < P4EST_DIM; dim++) {
      for (lu = 0; lu < 2; lu++) {
        if (qcoords[dim] != lu * P4EST_ROOT_LEN) {
          continue;             /* we are not at a tree face */
        }

        if (c->conconn->tree_to_tree[info->treeid * P4EST_FACES + dim * 2 +
                                     lu] != info->treeid) {
          continue;             /* the tree face is not on the forest boundary */
        }
        op->isboundary = 1;
      }
    }
  }
}

static void
overlap_consumer_evaluate_corners (p4est_iter_volume_info_t *info,
                                   void *user_data)
{
  p4est_quadrant_t   *q;
  overlap_consumer_t *c;
  overlap_point_t    *op;
  int                 i, npin, npout;
  overlap_condata_t  *d;

  P4EST_ASSERT (info != NULL && info->p4est != NULL);
  P4EST_ASSERT (info->p4est->user_pointer == user_data);

  /* access quadrant */
  c = (overlap_consumer_t *) info->p4est->user_pointer;
  P4EST_ASSERT (c->con4est == info->p4est);
  P4EST_ASSERT (c->lquad_idx >= 0 &&
                c->lquad_idx <
                P4EST_CHILDREN * c->con4est->local_num_quadrants);
  P4EST_ASSERT (info->quad != NULL);
  q = info->quad;

  /* iterate over all children */
  npin = npout = 0;
  for (i = 0; i < P4EST_CHILDREN; i++) {
    op = (overlap_point_t *) sc_array_index (c->query_xyz, c->lquad_idx++);
    if (op->prodata.isset) {
      npin++;
    }
    else {
      npout++;
    }
  }

  d = (overlap_condata_t *) q->p.user_data;
  if (npin && npout) {
    /* we have points inside and outside the producer domain, so, the boundary
     * crosses this quadrant */
    d->refine = 1;
  }
  else {
    d->refine = 0;
  }
}

static int          refine_level = 3;

static int
refine_producer_geometrical_fn (p4est_t *p4est, p4est_topidx_t which_tree,
                                p4est_quadrant_t *quadrant)
{
  overlap_producer_t *p;
  double              qxyz[3];
  double              phys[3] = { 0, 0, 0 };
  double              dist;
  P4EST_ASSERT (p4est != NULL);
  P4EST_ASSERT (p4est->user_pointer != NULL);
  p = (overlap_producer_t *) p4est->user_pointer;
  P4EST_ASSERT (p->progeom != NULL);
  P4EST_ASSERT (p->progeom->X != NULL);

  /* transform producer quadrant center to physical using map */
  get_quadrant_center (quadrant, qxyz);
  p->progeom->X (p->progeom, which_tree, qxyz, phys);

  /* compute distance from point of interest */
  phys[0] -= 0.7;
  phys[1] -= 0.2;
  phys[2] -= 0.45;
  dist = sqrt (phys[0] * phys[0] + phys[1] * phys[1] + phys[2] * phys[2]);

  /* refine quadrants that are close enough to point of interest */
  if (quadrant->level < refine_level - floor (dist / 0.2)) {
    return 1;
  }

  return 0;
}

static int
refine_consumer_geometrical_fn (p4est_t *p4est, p4est_topidx_t which_tree,
                                p4est_quadrant_t *quadrant)
{
  overlap_consumer_t *c;
  double              qxyz[3];
  double              phys[3] = { 0, 0, 0 };
  double              dist;

  P4EST_ASSERT (p4est != NULL);
  P4EST_ASSERT (p4est->user_pointer != NULL);
  c = (overlap_consumer_t *) p4est->user_pointer;
  P4EST_ASSERT (c->congeom != NULL);
  P4EST_ASSERT (c->congeom->X != NULL);

  /* transform producer quadrant center to physical using map */
  get_quadrant_center (quadrant, qxyz);
  c->congeom->X (c->congeom, which_tree, qxyz, phys);

  /* compute distance from point of interest */
  phys[0] -= 0.7;
  phys[1] -= 0.2;
  phys[2] -= 0.45;
  dist = sqrt (phys[0] * phys[0] + phys[1] * phys[1] + phys[2] * phys[2]);

  /* refine quadrants that are close enough to point of interest */
  if (quadrant->level < refine_level - 2 - floor (dist / 0.1)) {
    return 1;
  }

  return 0;
}

static int
refine_producer_fn (p4est_t *p4est, p4est_topidx_t which_tree,
                    p4est_quadrant_t *quadrant)
{
  if ((int) quadrant->level >= (refine_level - (int) (which_tree % 3))) {
    return 0;
  }
  if (quadrant->level == 1 && p4est_quadrant_child_id (quadrant) == 3) {
    return 1;
  }
  if (quadrant->x == P4EST_LAST_OFFSET (2) &&
      quadrant->y == P4EST_LAST_OFFSET (2)) {
    return 1;
  }
  if (quadrant->x >= P4EST_QUADRANT_LEN (2)) {
    return 0;
  }

  return 1;
}

static int
refine_consumer_fn (p4est_t *p4est, p4est_topidx_t which_tree,
                    p4est_quadrant_t *quadrant)
{
  if ((int) quadrant->level >= refine_level) {
    return 0;
  }
  if (p4est->mpirank >= 5 && p4est->mpirank <= 15) {
    return 1;
  }

  return 0;
}

static int
refine_producer_adaptive_fn (p4est_t *p4est, p4est_topidx_t which_tree,
                             p4est_quadrant_t *quadrant)
{
  overlap_prodata_t  *d = (overlap_prodata_t *) quadrant->p.user_data;

  if (d->isset == 1) {
    return 1;                   /* the producer quadrant contains a boundary consumer point */
  }

  return 0;
}

static int
refine_consumer_adaptive_fn (p4est_t *p4est, p4est_topidx_t which_tree,
                             p4est_quadrant_t *quadrant)
{
  overlap_condata_t  *d = (overlap_condata_t *) quadrant->p.user_data;

  if (d->refine == 1) {
    return 1;                   /* the quadrant contains both inside and outside query points */
  }
  return 0;
}

static void
overlap_query_centers (overlap_global_t *g)
{
  overlap_consumer_t *c = g->c;

  /* generate a query point for every local quadrant center */
  c->lquad_idx = 0;
  c->query_xyz = sc_array_new_count (sizeof (overlap_point_t),
                                     c->con4est->local_num_quadrants);
  c->interpolation_data =
    sc_array_new_count (sizeof (double), c->con4est->local_num_quadrants);
  c->isset_data =
    sc_array_new_count (sizeof (double), c->con4est->local_num_quadrants);
  c->xyz_data =
    sc_array_new_count (sizeof (double), 3 * c->con4est->local_num_quadrants);

  p4est_iterate (c->con4est, NULL, c, overlap_consumer_compute_center, NULL
#ifdef P4_TO_P8
                 , NULL
#endif
                 , NULL);
  P4EST_ASSERT (c->lquad_idx == c->con4est->local_num_quadrants);
}

static void
overlap_query_corners (overlap_global_t *g)
{
  overlap_consumer_t *c = g->c;

  /* generate a query point for every local quadrant center */
  c->lquad_idx = 0;
  c->query_xyz = sc_array_new_count (sizeof (overlap_point_t),
                                     P4EST_CHILDREN *
                                     c->con4est->local_num_quadrants);
  p4est_iterate (c->con4est, NULL, c, overlap_consumer_compute_corners, NULL
#ifdef P4_TO_P8
                 , NULL
#endif
                 , NULL);
  P4EST_ASSERT (c->lquad_idx ==
                P4EST_CHILDREN * c->con4est->local_num_quadrants);
}

static void
overlap_evaluate_corners (overlap_global_t *g)
{
  overlap_consumer_t *c = g->c;

  /* generate a query point for every local quadrant center */
  c->lquad_idx = 0;
  p4est_iterate (c->con4est, NULL, c, overlap_consumer_evaluate_corners, NULL
#ifdef P4_TO_P8
                 , NULL
#endif
                 , NULL);
  P4EST_ASSERT (c->lquad_idx ==
                P4EST_CHILDREN * c->con4est->local_num_quadrants);
}

static void
overlap_init_quadrant_prodata (p4est_iter_volume_info_t *info,
                               void *user_data)
{
  p4est_quadrant_t   *q;
  overlap_prodata_t  *d;

  P4EST_ASSERT (info->quad != NULL);
  q = info->quad;
  P4EST_ASSERT (q->p.user_data != NULL);
  d = (overlap_prodata_t *) q->p.user_data;
  d->isset = 0;
  d->myvalue = 0;
}

static void
overlap_init_prodata (overlap_producer_t *p)
{
  p4est_iterate (p->pro4est, NULL, p, overlap_init_quadrant_prodata, NULL
#ifdef P4_TO_P8
                 , NULL
#endif
                 , NULL);

}

static void         overlap_exchange (overlap_global_t *g);

static void
overlap_apps_init (overlap_global_t *g, sc_MPI_Comm mpicomm)
{
  overlap_producer_t *p = g->p = &g->pro;
  overlap_consumer_t *c = g->c = &g->con;
  int                 mpiret;
  int                 i, nrefine;
  p4est_connectivity_t *conns[2];

  /* initialization of global data */
  g->glocomm = mpicomm;
  g->rounds = 0;

  /* setup timing info */
  p->tstats = c->tstats = g->tstats;

  /* Create two connectivities. They will be assigned to the producer and the
   * consumer based on the value of g->example. */
  conns[0] = p4est_connectivity_new_brick (xbricks, ybricks
#ifdef P4_TO_P8
                                           , 2
#endif
                                           , 0, 0
#ifdef P4_TO_P8
                                           , 0
#endif
    );
  conns[1] = p4est_connectivity_new_brick (3, 2
#ifdef P4_TO_P8
                                           , 1
#endif
                                           , 0, 0
#ifdef P4_TO_P8
                                           , 0
#endif
    );

  /***************************** PRODUCER ****************************/

  P4EST_GLOBAL_PRODUCTION ("OVERLAP: init producer\n");

  /* setup producer geometry */
  p->progeom = &p->producer_geometry;
  p->progeom->name = "producer";
  p->progeom->destroy = (p4est_geometry_destroy_t) 0;

  /* setup producer app with communicator and mesh */
  p->glocomm = g->glocomm;
  mpiret = sc_MPI_Comm_dup (g->glocomm, &p->procomm);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (p->procomm, &p->prorank);
  SC_CHECK_MPI (mpiret);
  p->progeom->user = p->proconn = conns[!g->example];
  p->pro4est = p4est_new_ext (p->procomm, p->proconn, 0, p->pminl, 1,
                              sizeof (overlap_prodata_t), NULL, p);

  /* assign the geometry depending on the value of g->example */
  if (g->example) {
    p->progeom->X = overlap_curved_map;
    p->invmap = overlap_curved_invmap;
  }
  else {
    p->progeom->X = overlap_brick_map;
    p->invmap = overlap_brick_invmap;
  }

  /***************************** CONSUMER ****************************/

  P4EST_GLOBAL_PRODUCTION ("OVERLAP: init consumer\n");

  /* setup consumer geometry */
  c->congeom = &c->consumer_geometry;
  c->congeom->name = "consumer";
  c->congeom->destroy = (p4est_geometry_destroy_t) 0;

  /* setup consumer app with communicator and mesh */
  c->glocomm = g->glocomm;
  mpiret = sc_MPI_Comm_dup (g->glocomm, &c->concomm);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (c->concomm, &c->conrank);
  SC_CHECK_MPI (mpiret);
  c->congeom->user = c->conconn = conns[g->example];
  c->con4est = p4est_new_ext (c->concomm, c->conconn, 0, c->cminl, 1,
                              sizeof (overlap_condata_t), NULL, c);

  /* assign the geometry that was not assigned to the producer side */
  if (g->example) {
    c->invmap = overlap_curved_invmap;  /* the inverse producer mapping */
    c->congeom->X = overlap_brick_map;
  }
  else {
    c->invmap = overlap_brick_invmap;   /* the inverse producer mapping */
    c->congeom->X = overlap_curved_map;
  }

  /**************************** REFINEMENT ***************************/

  if (g->refinement_method == 0) {

    /* Adaptively refine the boundary of the mesh intersection area.
     * We query all corners of all consumer quadrants and refine all quadrants,
     * that contains at least one point that was found in the exchange and at
     * least one that was not found.
     * We mark all producer quadrants that contain a boundary query point for
     * refinement. */
    p->refining = 1;            /* set refinement flag to 1 */

    /* compute the maximum numbers of refinements to stay below refine_level */
    nrefine = refine_level - SC_MAX (c->cminl, p->pminl);

    for (i = 0; i < nrefine; i++) {
      overlap_init_prodata (p); /* overlap_exchange may not touch all quadrants */
      overlap_query_corners (g);

      /* query consumer corners and set p->refine_quadrant during the process */
      overlap_exchange (g);

      /* evaluate which consumer quadrants have to be refined */
      overlap_evaluate_corners (g);

      /* actual refinement based on the exchange results */
      p4est_refine (p->pro4est, 0, refine_producer_adaptive_fn, NULL);
      p4est_refine (c->con4est, 0, refine_consumer_adaptive_fn, NULL);

      /* cleanup */
      sc_array_destroy (g->c->query_xyz);
    }
  }
  else if (g->refinement_method == 1) {
    /* refine producer and consumer based on geometrical properties */
    p4est_refine (p->pro4est, 1, refine_producer_geometrical_fn, NULL);
    p4est_refine (c->con4est, 1, refine_consumer_geometrical_fn, NULL);
  }
  else {
    p4est_refine (p->pro4est, 1, refine_producer_fn, NULL);
    p4est_refine (c->con4est, 1, refine_consumer_fn, NULL);
  }

  p->refining = 0;              /* set refinement flag to 0 */
  overlap_query_centers (g);

  /* generate a local set of cell values by interpolating a function */
  p->lquad_idx = 0;
  p->interpolation_data =
    sc_array_new_count (sizeof (double), p->pro4est->local_num_quadrants);
  p->xyz_data =
    sc_array_new_count (sizeof (double), 3 * p->pro4est->local_num_quadrants);
  p4est_iterate (p->pro4est, NULL, p, overlap_producer_compute, NULL
#ifdef P4_TO_P8
                 , NULL
#endif
                 , NULL);

  /* here we would need to make the global partition encoding available to
     the consumer, if the communicators were not congruent */

  /* overlap works by intra-communication on the global communicator */
  P4EST_GLOBAL_PRODUCTION ("OVERLAP: init done\n");

}

static int
consumer_quadrant (p4est_t *p4est, p4est_topidx_t which_tree,
                   p4est_quadrant_t *quadrant, int pfirst, int plast,
                   void *point)
{
#ifdef P4EST_ENABLE_DEBUG
  overlap_consumer_t *c;

  /* Tree, quadrant, pfirst and plast refer to the producer. */

  P4EST_ASSERT (p4est != NULL);
  P4EST_ASSERT (quadrant != NULL);
  P4EST_ASSERT (point == NULL);

  c = (overlap_consumer_t *) p4est->user_pointer;
  P4EST_ASSERT (c != NULL);
#endif

  /* don't mess with the point search */
  return 1;
}

static void
overlap_consumer_add (overlap_consumer_t *c, overlap_point_t *op, int rank)
{
  size_t              bcount;
  overlap_send_buf_t *sb;

  P4EST_ASSERT (c != NULL && c->send_buffer != NULL);
  P4EST_ASSERT (op != NULL && op->rank == -1);
  P4EST_ASSERT (0 <= rank && rank < c->pronum_procs);
  op->rank = rank;

  /* if we have a new rank, push new send buffer */
  bcount = c->send_buffer->elem_count;
  sb = NULL;
  if (bcount > 0) {
    sb = (overlap_send_buf_t *)
      sc_array_index (c->send_buffer, bcount - 1);
    P4EST_ASSERT (sb->rank <= rank);
    P4EST_ASSERT (sb->ops.elem_count > 0);
  }
  if (bcount == 0 || sb->rank < rank) {
    sb = (overlap_send_buf_t *) sc_array_push (c->send_buffer);
    sb->rank = rank;
    sc_array_init (&sb->ops, sizeof (overlap_point_t));
  }
  memcpy (sc_array_push (&sb->ops), op, sizeof (overlap_point_t));
}

static int
producer_intersect (p4est_connectivity_t *pro_conn,
                    p4est_topidx_t which_tree, p4est_quadrant_t *quadrant,
                    overlap_point_t *op, double tol, overlap_invmap_t invmap)
{
  const double       *phys;
  double              abc[3], dh, dhz;
  double              qxyz[3];

  phys = op->xyz;

  /* transform point back to producer reference */
  invmap (pro_conn, which_tree, phys, abc);

  P4EST_LDEBUGF ("Point %ld is %g %g %g\n",
                 (long) op->lnum, phys[0], phys[1], phys[2]);
  P4EST_LDEBUGF ("Tree %d level %d invert to %g %g %g\n",
                 (int) which_tree, quadrant->level, abc[0], abc[1], abc[2]);

  /* check for tree intersection */
  if ((abc[0] <= -tol || abc[0] >= 1. + tol) ||
      (abc[1] <= -tol || abc[1] >= 1. + tol) ||
      (abc[2] <= -tol || abc[2] >= 1. + tol)) {
    return 0;
  }

  P4EST_LDEBUGF ("Point %ld survive tree\n", (long) op->lnum);

  /* check for quadrant intersection */
  dh = OVERLAP_IROOTLEN * P4EST_QUADRANT_LEN (quadrant->level);
  qxyz[0] = OVERLAP_IROOTLEN * quadrant->x;
  qxyz[1] = OVERLAP_IROOTLEN * quadrant->y;
#ifndef P4_TO_P8
  qxyz[2] = 0.;
  dhz = 1.;
#else
  qxyz[2] = OVERLAP_IROOTLEN * quadrant->z;
  dhz = dh;
#endif
  if ((abc[0] < qxyz[0] - tol || abc[0] > qxyz[0] + dh + tol) ||
      (abc[1] < qxyz[1] - tol || abc[1] > qxyz[1] + dh + tol) ||
      (abc[2] < qxyz[2] - tol || abc[2] > qxyz[2] + dhz + tol)) {
    return 0;
  }

  P4EST_LDEBUGF ("Point %ld survive quadrant\n", (long) op->lnum);
  return 1;
}

static int
consumer_point (p4est_t *p4est, p4est_topidx_t which_tree,
                p4est_quadrant_t *quadrant, int pfirst, int plast,
                void *point)
{
  overlap_point_t    *op = (overlap_point_t *) point;
  overlap_consumer_t *c;
  int                 intersects;

  /* The point is owned by the consumer.
     Tree, quadrant, pfirst and plast refer to the producer. */

  P4EST_ASSERT (p4est != NULL);
  P4EST_ASSERT (quadrant != NULL);
  P4EST_ASSERT (op != NULL);
  if (op->rank >= 0) {
    /* skip a point of multiple intersections */
    P4EST_INFOF ("Skip point %ld rank %d on multiple match\n",
                 (long) op->lnum, op->rank);
    return 0;
  }

  c = (overlap_consumer_t *) p4est->user_pointer;
  P4EST_ASSERT (c != NULL);
  P4EST_ASSERT (0 <= which_tree && which_tree < c->pronum_trees);

  /* check if the point intersects the quadrant */
  P4EST_LDEBUGF ("Consumer point %ld intersection test\n", (long) op->lnum);
  c->tstats->stats[OVERLAP_NUM_CONS_SEARCH_OPS].sum_values++;
  intersects =
    producer_intersect (c->producer_conn, which_tree, quadrant, op,
                        P4EST_CON_TOLERANCE, c->invmap);
  if (!intersects) {
    return 0;
  }

  /* we have located the point in the intersection quadrant */
  if (pfirst == plast) {
    /* we have intersected with a leaf quadrant */
    P4EST_INFOF ("Point %ld leaf intersect tree %d rank %d\n",
                 (long) op->lnum, (int) which_tree, pfirst);
    overlap_consumer_add (c, op, pfirst);
  }
  return 1;
}

static int
producer_point (p4est_t *p4est, p4est_topidx_t which_tree,
                p4est_quadrant_t *quadrant, p4est_locidx_t local_num,
                void *point)
{
  int                 isleaf, intersects;
  overlap_point_t    *op = (overlap_point_t *) point;
  overlap_producer_t *p = (overlap_producer_t *) p4est->user_pointer;

  P4EST_ASSERT (p4est != NULL);
  P4EST_ASSERT (quadrant != NULL);
  P4EST_ASSERT (op != NULL);
  P4EST_ASSERT (p != NULL);
  if (op->prodata.isset) {
    /* skip a point of multiple intersections */
    P4EST_INFOF ("Skip point %ld on multiple match\n", (long) op->lnum);
    return 0;
  }

  /* check if the point intersects the quadrant */
  P4EST_LDEBUGF ("Producer point %ld intersection test\n", (long) op->lnum);
  p->tstats->stats[OVERLAP_NUM_PROD_SEARCH_OPS].sum_values++;
  intersects =
    producer_intersect (p->proconn, which_tree, quadrant, op,
                        P4EST_PRO_TOLERANCE, p->invmap);
  if (!intersects) {
    return 0;
  }

  isleaf = local_num >= 0;
  if (isleaf) {
    op->prodata.isset = 1;
    overlap_prodata_t  *d = (overlap_prodata_t *) quadrant->p.user_data;
    P4EST_ASSERT (d != NULL);
    if (p->refining && op->isboundary) {
      /* mark that the producer intersects a consumer boundary point */
      d->isset = 1;
    }
    else {
      /* apply producer interpolation data to consumer point */
      op->prodata.myvalue = d->myvalue;
      P4EST_LDEBUGF ("Producer point %ld prodata set to %f.\n",
                     (long) op->lnum, op->prodata.myvalue);
    }
  }

  return 1;
}

static void
consumer_search_partition (overlap_consumer_t *c)
{
  c->send_buffer = sc_array_new (sizeof (overlap_send_buf_t));
  p4est_search_partition_gfx (c->producer_gfq, c->producer_gfp,
                              c->pronum_procs, c->pronum_trees, 0,
                              c, consumer_quadrant, consumer_point,
                              c->query_xyz);
}

#ifdef P4EST_ENABLE_MPI

static void
consumer_producer_notify (overlap_global_t *g)
{
  overlap_producer_t *p = g->p;
  overlap_consumer_t *c = g->c;
  size_t              bz, bcount;
  sc_array_t         *receivers, *senders;
  sc_array_t         *payload_in, *payload_out;
  int                 num_receivers, num_senders;
  overlap_send_buf_t *sb;
  overlap_recv_buf_t *rb;
  int                 same_rank, num_ops, i;
  int                 mpiret;
  sc_notify_t        *notifyc;

  /* assemble and execute receiver and payload query */
  num_receivers = (int) (bcount = c->send_buffer->elem_count);
  receivers = sc_array_new_count (sizeof (int), bcount);
  senders = sc_array_new (sizeof (p4est_locidx_t));
  payload_in = sc_array_new_count (sizeof (int), bcount);
  payload_out = sc_array_new (sizeof (p4est_locidx_t));
  for (bz = 0; bz < bcount; ++bz) {
    sb = (overlap_send_buf_t *) sc_array_index (c->send_buffer, bz);
    *(int *) sc_array_index (receivers, bz) = sb->rank;
    *(p4est_locidx_t *) sc_array_index (payload_in, bz) =
      (p4est_locidx_t) sb->ops.elem_count;
    g->tstats->stats[OVERLAP_NUM_QP_SENT].sum_values += sb->ops.elem_count;
  }
  notifyc = sc_notify_new (g->glocomm);
  sc_notify_set_type (notifyc, SC_NOTIFY_NBX);
  sc_notify_payload (receivers, senders, payload_in, payload_out, 1, notifyc);
  sc_notify_destroy (notifyc);
  num_senders = (int) senders->elem_count;
  P4EST_INFOF ("Overlap exchange receivers %d senders %d\n",
               num_receivers, num_senders);

  /* post nonblocking receives for the point data of the consumer side */
  p->recv_buffer = sc_array_new_count (sizeof (overlap_recv_buf_t),
                                       num_senders);
  p->recv_reqs = sc_array_new_count (sizeof (sc_MPI_Request), num_senders);
  p->iprorank = c->iconrank = -1;
  for (i = 0; i < num_senders; ++i) {
    /* initalize and allocate the buffer according to the payload */
    rb = (overlap_recv_buf_t *) sc_array_index_int (p->recv_buffer, i);
    rb->rank = *(int *) sc_array_index_int (senders, i);
    same_rank = (rb->rank == p->prorank);
    g->tstats->stats[OVERLAP_NUM_QP_RECEIVED].sum_values +=
      *(int *) sc_array_index_int (payload_out, i);
    num_ops = same_rank ? 0 : *(int *) sc_array_index_int (payload_out, i);
    sc_array_init_size (&(rb->ops), sizeof (overlap_point_t),
                        (size_t) num_ops);

    if (same_rank) {
      p->iprorank = i;          /* save the index in the producer buffer */
      *(sc_MPI_Request *) sc_array_index_int (p->recv_reqs, i) =
        sc_MPI_REQUEST_NULL;
      continue;
    }

    /* receive the array of overlap_point_t data and store it in the buffer */
    mpiret =
      sc_MPI_Irecv (rb->ops.array, num_ops * sizeof (overlap_point_t),
                    sc_MPI_BYTE, rb->rank, COMM_TAG_CONSDATA, p->glocomm,
                    (sc_MPI_Request *) sc_array_index_int (p->recv_reqs, i));
    SC_CHECK_MPI (mpiret);
  }

  sc_array_destroy (receivers);
  sc_array_destroy (senders);
  sc_array_destroy (payload_in);
  sc_array_destroy (payload_out);

}

static void
consumer_post_messages (overlap_consumer_t *c)
{
  overlap_send_buf_t *sb;
  overlap_recv_buf_t *rb;
  int                 num_receivers, same_rank, num_ops, i;
  int                 mpiret;

  /* send the point data to the producer side in a nonblocking way */
  num_receivers = (int) c->send_buffer->elem_count;
  c->send_reqs = sc_array_new_count (sizeof (sc_MPI_Request), num_receivers);
  for (i = 0; i < num_receivers; ++i) {
    sb = (overlap_send_buf_t *) sc_array_index_int (c->send_buffer, i);

    if (sb->rank == c->conrank) {
      c->iconrank = i;          /* save the index in the consumer buffer */
      *(sc_MPI_Request *) sc_array_index_int (c->send_reqs, i) =
        sc_MPI_REQUEST_NULL;
      continue;
    }

    mpiret =
      sc_MPI_Isend (sb->ops.array,
                    sb->ops.elem_count * sizeof (overlap_point_t),
                    sc_MPI_BYTE, sb->rank, COMM_TAG_CONSDATA, c->glocomm,
                    (sc_MPI_Request *) sc_array_index_int (c->send_reqs, i));
    SC_CHECK_MPI (mpiret);
  }

  /* recv the updated point data from the producer side in a nonblocking way */
  c->recv_reqs = sc_array_new_count (sizeof (sc_MPI_Request), num_receivers);
  c->recv_buffer = sc_array_new_size (sizeof (overlap_recv_buf_t),
                                      num_receivers);
  for (i = 0; i < num_receivers; ++i) {
    rb = (overlap_recv_buf_t *) sc_array_index_int (c->recv_buffer, i);
    sb = (overlap_send_buf_t *) sc_array_index_int (c->send_buffer, i);
    rb->rank = sb->rank;
    same_rank = (rb->rank == c->conrank);
    num_ops = same_rank ? 0 : (int) sb->ops.elem_count;
    sc_array_init_size (&(rb->ops), sizeof (overlap_point_t),
                        (size_t) num_ops);

    if (same_rank) {
      *(sc_MPI_Request *) sc_array_index_int (c->recv_reqs, c->iconrank) =
        sc_MPI_REQUEST_NULL;
      continue;
    }

    /* receive the array of overlap_point_t data and store it in the buffer */
    mpiret =
      sc_MPI_Irecv (rb->ops.array,
                    rb->ops.elem_count * sizeof (overlap_point_t),
                    sc_MPI_BYTE, rb->rank, COMM_TAG_PRODATA, c->glocomm,
                    (sc_MPI_Request *) sc_array_index_int (c->recv_reqs, i));
    SC_CHECK_MPI (mpiret);
  }
}

static void
producer_interpolate (overlap_producer_t *p)
{
  overlap_recv_buf_t *rb;
  int                 num_senders, i;
  int                 remaining, received;
  int                *prod_indices;
  int                 mpiret;

  /* compute producer data for all incoming messages as soon as they come in */
  num_senders = (int) p->recv_reqs->elem_count;
  prod_indices = P4EST_ALLOC (int, num_senders);
  p->send_reqs = sc_array_new_count (sizeof (sc_MPI_Request), num_senders);
  remaining = num_senders;
  if (p->iprorank >= 0) {
    *(sc_MPI_Request *) sc_array_index_int (p->send_reqs, p->iprorank) =
      sc_MPI_REQUEST_NULL;
    remaining--;                /* since we set the iprorank-th request to null earlier */
  }
  while (remaining > 0) {
    mpiret =
      sc_MPI_Waitsome (num_senders, (sc_MPI_Request *) p->recv_reqs->array,
                       &received, prod_indices, sc_MPI_STATUSES_IGNORE);
    SC_CHECK_MPI (mpiret);
    P4EST_ASSERT (received != sc_MPI_UNDEFINED);
    P4EST_ASSERT (received > 0);

    for (i = 0; i < received; ++i) {
      /* compute the prodata for the points sent by process prod_indices[i] */
      P4EST_ASSERT (0 <= prod_indices[i] && prod_indices[i] < num_senders);
      rb = (overlap_recv_buf_t *) sc_array_index_int (p->recv_buffer,
                                                      prod_indices[i]);
      p4est_search_local (p->pro4est, 0, NULL, producer_point, &(rb->ops));

      /* send the requested producer data back in a nonblocking way */
      mpiret =
        sc_MPI_Isend (rb->ops.array,
                      rb->ops.elem_count * sizeof (overlap_point_t),
                      sc_MPI_BYTE, rb->rank, COMM_TAG_PRODATA, p->glocomm,
                      (sc_MPI_Request *) sc_array_index_int (p->send_reqs,
                                                             prod_indices
                                                             [i]));
      SC_CHECK_MPI (mpiret);
    }

    remaining -= received;
  }

  P4EST_FREE (prod_indices);
}

#endif /* P4EST_ENABLE_MPI */

static void
consumer_update_from_buffer
  (sc_array_t *query_xyz, sc_array_t *buffer, int bi)
{
  overlap_recv_buf_t *rb;
  overlap_point_t    *op, *qp;
  int                 j;

  /* obtain the array of points we want to update query_xyz with */
  P4EST_ASSERT (0 <= bi && bi < (int) buffer->elem_size);
  rb = (overlap_recv_buf_t *) sc_array_index_int (buffer, bi);

  /* copy prodata into the query-point array */
  for (j = 0; j < (int) rb->ops.elem_count; ++j) {
    op = (overlap_point_t *) sc_array_index_int (&(rb->ops), j);
    qp = (overlap_point_t *) sc_array_index_int (query_xyz,
                                                 (size_t) op->lnum);
    qp->prodata.isset = op->prodata.isset;
    qp->prodata = op->prodata;
  }
}

#ifdef P4EST_ENABLE_MPI

static void
consumer_update_query_points (overlap_consumer_t *c)
{
  int                 num_receivers, i;
  int                 remaining, received;
  int                *cons_indices;
  int                 mpiret;

  /* compute producer data for all incoming messages as soon as they come in */
  num_receivers = (int) c->recv_reqs->elem_count;
  cons_indices = P4EST_ALLOC (int, num_receivers);
  remaining = (c->iconrank >= 0) ? num_receivers - 1 : num_receivers;
  while (remaining > 0) {
    mpiret =
      sc_MPI_Waitsome (num_receivers,
                       (sc_MPI_Request *) c->recv_reqs->array, &received,
                       cons_indices, sc_MPI_STATUSES_IGNORE);
    SC_CHECK_MPI (mpiret);
    P4EST_ASSERT (received != sc_MPI_UNDEFINED);
    P4EST_ASSERT (received > 0);

    for (i = 0; i < received; ++i) {
      consumer_update_from_buffer (c->query_xyz, c->recv_buffer,
                                   cons_indices[i]);
    }

    remaining -= received;
  }

  P4EST_FREE (cons_indices);
}

static void
consumer_waitall (overlap_consumer_t *c)
{
  int                 mpiret;
  int                 num_receivers;

  /* wait for the nonblocking sends to complete */
  num_receivers = (int) c->send_reqs->elem_count;
  mpiret =
    sc_MPI_Waitall (num_receivers, (sc_MPI_Request *) c->send_reqs->array,
                    sc_MPI_STATUSES_IGNORE);
  SC_CHECK_MPI (mpiret);
}

static void
producer_waitall (overlap_producer_t *p)
{
  int                 mpiret;
  int                 num_senders;

  /* wait for the nonblocking sends to complete */
  num_senders = (int) p->send_reqs->elem_count;
  mpiret =
    sc_MPI_Waitall (num_senders, (sc_MPI_Request *) p->send_reqs->array,
                    sc_MPI_STATUSES_IGNORE);
  SC_CHECK_MPI (mpiret);
}

#endif /* P4EST_ENABLE_MPI */

static void
consumer_producer_update_local (overlap_global_t *g)
{
  overlap_consumer_t *c = g->c;
  overlap_producer_t *p = g->p;
  overlap_send_buf_t *sb;

  if (c->iconrank >= 0 && c->send_buffer->elem_count) {
    /* Interpolate point-data of local points. Instead of copying to the
     * producer buffer, we update the points in-place. */
    sb =
      (overlap_send_buf_t *) sc_array_index_int (c->send_buffer, c->iconrank);
    p4est_search_local (p->pro4est, 0, NULL, producer_point, &(sb->ops));
    consumer_update_from_buffer (c->query_xyz, c->send_buffer, c->iconrank);
  }
}

#if 0
static void
consumer_print_interpolation_data (overlap_consumer_t *c)
{
  overlap_point_t    *qp;
  p4est_locidx_t      cind;

  for (cind = 0; cind < c->con4est->local_num_quadrants; cind++) {
    qp = (overlap_point_t *) sc_array_index (c->query_xyz, cind);
    if (qp->prodata.isset) {
      P4EST_INFOF ("Consumer query point %ld assigned prodata %f\n",
                   (long) qp->lnum, qp->prodata.myvalue);
    }
    else {
      P4EST_INFOF ("Consumer query point %ld not in producer domain.\n",
                   (long) qp->lnum);
    }

    /* store vtk cell data */
    *(double *) sc_array_index (c->interpolation_data, cind) =
      qp->prodata.myvalue;
    *(double *) sc_array_index (c->isset_data, cind) =
      (double) qp->prodata.isset;
  }
}
#endif

/* write consumer p4est with interpolation data into vtk */
void
consumer_write_vtk (overlap_consumer_t *c)
{
  int                 retval;
  p4est_vtk_context_t *con_context;

  /* allocate context and set parameters */
  con_context =
    p4est_vtk_context_new (c->con4est, P4EST_STRING "_consumer_new");
  p4est_vtk_context_set_geom (con_context, c->congeom);
  p4est_vtk_context_set_continuous (con_context, 1);

  /* write header */
  con_context = p4est_vtk_write_header (con_context);
  SC_CHECK_ABORT (con_context != NULL,
                  P4EST_STRING "_vtk: Error writing header");

  /* write cell_data */
  con_context =
    p4est_vtk_write_cell_dataf (con_context, 1, 1, 1, 0, 2, 1,
                                "interpolation", c->interpolation_data,
                                "is_set", c->isset_data,
                                "xyz", c->xyz_data, con_context);
  SC_CHECK_ABORT (con_context != NULL,
                  P4EST_STRING "_vtk: Error writing celldata");

  /* properly write rest of the files' contents */
  retval = p4est_vtk_write_footer (con_context);
  SC_CHECK_ABORT (!retval, P4EST_STRING "_vtk: Error writing footer");
}

/* write producer p4est with interpolation data into vtk */
void
producer_write_vtk (overlap_producer_t *p)
{
  int                 retval;
  p4est_vtk_context_t *pro_context;

  /* allocate context and set parameters */
  pro_context =
    p4est_vtk_context_new (p->pro4est, P4EST_STRING "_producer_new");
  p4est_vtk_context_set_geom (pro_context, p->progeom);
  p4est_vtk_context_set_continuous (pro_context, 1);

  /* write header */
  pro_context = p4est_vtk_write_header (pro_context);
  SC_CHECK_ABORT (pro_context != NULL,
                  P4EST_STRING "_vtk: Error writing header");

  /* write cell_data */
  pro_context =
    p4est_vtk_write_cell_dataf (pro_context, 1, 1, 1, 0, 1, 1,
                                "interpolation", p->interpolation_data,
                                "xyz", p->xyz_data, pro_context);
  SC_CHECK_ABORT (pro_context != NULL,
                  P4EST_STRING "_vtk: Error writing celldata");

  /* properly write rest of the files' contents */
  retval = p4est_vtk_write_footer (pro_context);
  SC_CHECK_ABORT (!retval, P4EST_STRING "_vtk: Error writing footer");
}

static void
consumer_free_communication_data (overlap_consumer_t *c)
{
  overlap_send_buf_t *sb;
#ifdef P4EST_ENABLE_MPI
  overlap_recv_buf_t *rb;
#ifdef P4EST_ENABLE_DEBUG
  int                 prev_rank;
#endif
  size_t              bz, bcount;

  sc_array_destroy (c->recv_reqs);
  sc_array_destroy (c->send_reqs);
#ifdef P4EST_ENABLE_DEBUG
  prev_rank = -1;
#endif
  bcount = c->send_buffer->elem_count;
  for (bz = 0; bz < bcount; ++bz) {
    sb = (overlap_send_buf_t *) sc_array_index (c->send_buffer, bz);
    rb = (overlap_recv_buf_t *) sc_array_index (c->recv_buffer, bz);
    SC_ASSERT (sb->rank == rb->rank);
    SC_ASSERT (prev_rank < sb->rank);
#ifdef P4EST_ENABLE_DEBUG
    prev_rank = sb->rank;
    if (bz == (size_t) c->iconrank) {
      P4EST_ASSERT (rb->ops.elem_count == 0);
    }
    else {
      P4EST_ASSERT (rb->ops.elem_count == sb->ops.elem_count);
    }
#endif
    P4EST_ASSERT (sb->ops.elem_count > 0);
    sc_array_reset (&sb->ops);
    sc_array_reset (&rb->ops);
  }
  sc_array_destroy_null (&c->recv_buffer);
#else /* !P4EST_ENABLE_MPI */
  if (c->send_buffer->elem_count) {
    sb = (overlap_send_buf_t *) sc_array_index_int (c->send_buffer, 0);
    sc_array_reset (&sb->ops);
  }
#endif
  sc_array_destroy_null (&c->send_buffer);
}

static void
producer_free_communication_data (overlap_producer_t *p)
{
#ifdef P4EST_ENABLE_MPI
  overlap_recv_buf_t *rb;
  int                 num_senders, i;

  sc_array_destroy (p->recv_reqs);
  sc_array_destroy (p->send_reqs);
  num_senders = (int) p->recv_buffer->elem_count;
  for (i = 0; i < num_senders; ++i) {
    rb = (overlap_recv_buf_t *) sc_array_index_int (p->recv_buffer, i);
#ifdef P4EST_ENABLE_DEBUG
    if (i == p->iprorank) {
      P4EST_ASSERT (rb->ops.elem_count == 0);
    }
    else {
      P4EST_ASSERT (rb->ops.elem_count > 0);
    }
#endif
    sc_array_reset (&rb->ops);
  }
  sc_array_destroy_null (&p->recv_buffer);
#endif
}

static void
overlap_exchange (overlap_global_t *g)
{
  overlap_producer_t *p = g->p;
  overlap_consumer_t *c = g->c;
  int                 istat;
  sc_flopinfo_t       snapshot, snapshot_total, *fi;
  sc_statinfo_t      *stats = g->tstats->stats;

  /* total time of the exchange function */
  fi = &g->tstats->fi;
  sc_flops_snap (fi, &snapshot_total);
  P4EST_GLOBAL_PRODUCTION ("OVERLAP: exchange partition\n");

  /* consumer receives global partition encoding from producer */
  /* since their communicators are congruent, this is a copy */
  c->producer_gfq = p->pro4est->global_first_quadrant;
  c->producer_gfp = p->pro4est->global_first_position;
  c->pronum_procs = p->pro4est->mpisize;
  c->pronum_trees = (c->producer_conn = p->proconn)->num_trees;

  P4EST_GLOBAL_PRODUCTION ("OVERLAP: customer partition search\n");

  /* search for the query points in the producer-partition and create a buffer
   * to send them to the respective producer ranks */
  sc_flops_snap (fi, &snapshot);
  consumer_search_partition (c);
  sc_flops_shot (fi, &snapshot);
  sc_stats_set1 (&stats[OVERLAP_SEARCH_PARTITION], snapshot.iwtime,
                 "Search partition");

#ifdef P4EST_ENABLE_MPI
  /* global, communication-based part of the interpolation */
  /* during this process we will mark messages that allow for a local, in-place
   * solution by setting c->iconrank */

  /* notify the producer about the point-array-messages it will receive,
   * allocate an receive buffer according to the transmitted payloads and
   * post Irecvs for the point-arrays */
  sc_flops_snap (fi, &snapshot);
  consumer_producer_notify (g);
  sc_flops_shot (fi, &snapshot);
  sc_stats_set1 (&stats[OVERLAP_NOTIFY], snapshot.iwtime,
                 "Consumer producer notify");

  /* post Isends for the point-arrays as well as Irecvs for the updated
   * point-arrays containing the interpolation prodata */
  sc_flops_snap (fi, &snapshot);
  consumer_post_messages (c);
  sc_flops_shot (fi, &snapshot);
  sc_stats_set1 (&stats[OVERLAP_POST_MESSAGES], snapshot.iwtime,
                 "Consumer post messages");

  P4EST_GLOBAL_PRODUCTION ("OVERLAP: producer local search\n");

  /* interpolate the point-arrays as soon as they arrive and send them back to
   * the consumer side in a non-blocking way */
  sc_flops_snap (fi, &snapshot);
  producer_interpolate (p);
  sc_flops_shot (fi, &snapshot);
  sc_stats_set1 (&stats[OVERLAP_INTERPOLATE], snapshot.iwtime,
                 "Producer interpolate");

  P4EST_GLOBAL_PRODUCTION ("OVERLAP: consumer query point update\n");

  /* compute the interpolation data of the query points based on the
   * updated point-arrays */
  sc_flops_snap (fi, &snapshot);
  consumer_update_query_points (c);
  sc_flops_shot (fi, &snapshot);
  sc_stats_set1 (&stats[OVERLAP_UPDATE_QUERY_POINTS],
                 snapshot.iwtime, "Consumer update query points");

  /* wait for the communication to complete */
  sc_flops_snap (fi, &snapshot);
  consumer_waitall (c);
  producer_waitall (p);
  sc_flops_shot (fi, &snapshot);
  sc_stats_set1 (&stats[OVERLAP_WAITALL], snapshot.iwtime,
                 "Consumer producer waitall");

#else /* !P4EST_ENABLE_MPI */
  c->iconrank = 0;              /* indicate that the send buffer can be updated directly */
#endif

  P4EST_GLOBAL_PRODUCTION ("OVERLAP: local interpolation\n");

  /* local, in-place part of the interpolation */
  sc_flops_snap (fi, &snapshot);
  consumer_producer_update_local (g);
  sc_flops_shot (fi, &snapshot);
  sc_stats_set1 (&stats[OVERLAP_UPDATE_LOCAL], snapshot.iwtime,
                 "Consumer producer update local");

#if 0                           /* we do not want to output result data during our measurements */
  if (!p->refining) {
    /* we are not in an adaptive refinement query, output the resulting
     * interpolation data of all query points */
    consumer_print_interpolation_data (c);

    consumer_write_vtk (c);
    producer_write_vtk (p);
  }
#endif

  /* free remaining communication data */
  sc_flops_snap (fi, &snapshot);
  consumer_free_communication_data (c);
  producer_free_communication_data (p);
  sc_flops_shot (fi, &snapshot);
  sc_stats_set1 (&stats[OVERLAP_FREE_COMMUNICATION_DATA],
                 snapshot.iwtime,
                 "Consumer producer free communication data");

  /* finish timings and stats */
  sc_flops_shot (&g->tstats->fi, &snapshot_total);
  sc_stats_set1 (&stats[OVERLAP_EXCHANGE], snapshot_total.iwtime, "Exchange");
  sc_stats_set1 (&stats[OVERLAP_NUM_LOCAL_CONS_QUADRANTS],
                 g->c->con4est->local_num_quadrants,
                 "Number local consumer quadrants");
  sc_stats_set1 (&stats[OVERLAP_NUM_LOCAL_PROD_QUADRANTS],
                 g->p->pro4est->local_num_quadrants,
                 "Number local producer quadrants");

  /* calculate and print timings */
  stats[OVERLAP_NUM_QP_SENTRECVD].sum_values =
    stats[OVERLAP_NUM_QP_SENT].sum_values
    + stats[OVERLAP_NUM_QP_RECEIVED].sum_values;
  sc_stats_collapse (&stats[OVERLAP_NUM_QP_SENT]);
  sc_stats_collapse (&stats[OVERLAP_NUM_QP_RECEIVED]);
  sc_stats_collapse (&stats[OVERLAP_NUM_QP_SENTRECVD]);
  stats[OVERLAP_NUM_SEARCH_OPS].sum_values =
    stats[OVERLAP_NUM_CONS_SEARCH_OPS].sum_values
    + stats[OVERLAP_NUM_PROD_SEARCH_OPS].sum_values;
  sc_stats_collapse (&stats[OVERLAP_NUM_PROD_SEARCH_OPS]);
  sc_stats_collapse (&stats[OVERLAP_NUM_CONS_SEARCH_OPS]);
  sc_stats_collapse (&stats[OVERLAP_NUM_SEARCH_OPS]);
  sc_stats_compute (g->glocomm, OVERLAP_NUM_STATS, stats);
  /* sc_stats_print_x works the same as sc_stats_print, but takes an array
   * that indicates, if the stat is a double or an integer, to decide between
   * %g and %f. We use a hardcoded overlap_stats_type for all OVERLAP_NUM_STATS
   * stats */
  sc_stats_print_x (p4est_package_id, SC_LP_ESSENTIAL, OVERLAP_NUM_STATS,
                    stats, overlap_stats_type, 1, 1);

  for (istat = 0; istat < OVERLAP_NUM_STATS; istat++) {
    sc_stats_reset (&stats[istat], 0);
  }
}

static void
overlap_update (overlap_global_t *g)
{
  /* possibly modify mesh and data for next round in main program */
}

static void
overlap_apps_reset (overlap_global_t *g)
{
  overlap_producer_t *p = g->p;
  overlap_consumer_t *c = g->c;
  int                 mpiret;

  /* destroy producer */
  sc_array_destroy (p->interpolation_data);
  sc_array_destroy (p->xyz_data);
  p4est_destroy (p->pro4est);
  p4est_connectivity_destroy (p->proconn);
  mpiret = sc_MPI_Comm_free (&p->procomm);
  SC_CHECK_MPI (mpiret);

  /* destroy consumer */
  sc_array_destroy (c->interpolation_data);
  sc_array_destroy (c->isset_data);
  sc_array_destroy (c->xyz_data);
  sc_array_destroy (c->query_xyz);
  p4est_destroy (c->con4est);
  p4est_connectivity_destroy (c->conconn);
  mpiret = sc_MPI_Comm_free (&c->concomm);
  SC_CHECK_MPI (mpiret);
}

int
main (int argc, char **argv)
{
  int                 i;
  int                 mpiret;
  int                 first_argc;
  sc_MPI_Comm         mpicomm;
  sc_options_t       *opt;
  overlap_tstats_t    tstats;
  overlap_global_t global, *g = &global;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  mpicomm = sc_MPI_COMM_WORLD;
  sc_init (mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
#ifndef P4EST_ENABLE_DEBUG
  sc_set_log_defaults (NULL, NULL, SC_LP_STATISTICS);
#endif
  p4est_init (NULL, SC_LP_DEFAULT);

  /* process command line arguments */
  opt = sc_options_new (argv[0]);
  sc_options_add_int (opt, 'c', "cons_minlevel", &g->con.cminl, 0,
                      "Lowest consumer level");
  sc_options_add_int (opt, 'p', "prod_minlevel", &g->pro.pminl, 0,
                      "Lowest producer level");
  sc_options_add_int (opt, 'e', "example", &g->example, 0,
                      "Example mapping index");
  sc_options_add_int (opt, 'r', "refine_option", &g->refinement_method, 0,
                      "Refinement pattern");
  sc_options_add_int (opt, 'm', "max_level", &refine_level, 3,
                      "Maximum refinement level");

  first_argc = sc_options_parse (p4est_package_id, SC_LP_DEFAULT,
                                 opt, argc, argv);
  if (first_argc < 0 || first_argc != argc) {
    sc_options_print_usage (p4est_package_id, SC_LP_ERROR, opt, NULL);
    return EXIT_FAILURE;
  }
  sc_options_print_summary (p4est_package_id, SC_LP_ESSENTIAL, opt);

  /* start overall timing */
  mpiret = sc_MPI_Barrier (mpicomm);
  SC_CHECK_MPI (mpiret);
  sc_flops_start (&tstats.fi);
  g->tstats = &tstats;
  sc_stats_init (&g->tstats->stats[OVERLAP_NUM_QP_SENT],
                 "Number query points sent");
  sc_stats_init (&g->tstats->stats[OVERLAP_NUM_QP_RECEIVED],
                 "Number query points received");
  sc_stats_init (&g->tstats->stats[OVERLAP_NUM_QP_SENTRECVD],
                 "Number query points sent and received");
  sc_stats_init (&g->tstats->stats[OVERLAP_NUM_PROD_SEARCH_OPS],
                 "Number producer intersection tests");
  sc_stats_init (&g->tstats->stats[OVERLAP_NUM_CONS_SEARCH_OPS],
                 "Number consumer intersection tests");
  sc_stats_init (&g->tstats->stats[OVERLAP_NUM_SEARCH_OPS],
                 "Number intersection tests");

  overlap_apps_init (g, mpicomm);

  overlap_exchange (g);

  for (i = 0; i < g->rounds; ++i) {
    P4EST_GLOBAL_PRODUCTIONF ("Into round %d/%d\n", i, g->rounds);
    overlap_update (g);
    overlap_exchange (g);
  }

  overlap_apps_reset (g);

  sc_options_destroy (opt);
  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return EXIT_SUCCESS;
}
