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

/*
 * Usage:
 * p8est_tsearch [-l <LEVEL>] [-s <LEVEL-SHIFT>] [-N <NUM-POINTS>] [-V]
 *
 * NUM-POINTS determines how many points are positioned randomly and searched
 * in the forest.  The domain is the spherical shell with radii 0.55 and 1 and
 * the points are uniformly distributed in the unit cube.
 */

#include <p4est_to_p8est.h>
#ifndef P4_TO_P8
#error "This program is currently intended for 3D only"
#else
#include <p8est_bits.h>
#include <p8est_extended.h>
#include <p8est_search.h>
#include <p8est_vtk.h>
#endif
#include <sc_flops.h>
#include <sc_options.h>
#include <sc_statistics.h>

static int          refine_level, level_shift;

typedef enum tsearch_stats
{
  TSEARCH_NEW,
  TSEARCH_REFINE,
  TSEARCH_BALANCE,
  TSEARCH_PARTITION,
  TSEARCH_SEARCH_1,
  TSEARCH_SEARCH_N,
  TSEARCH_NUM_STATS
}
tsearch_stats_t;

typedef struct
{
  /* MPI related data */
  sc_MPI_Comm         mpicomm;
  int                 mpisize;
  int                 mpirank;

  /* global data for the program */
  size_t              matches;
  double              rout, rin;
  double              rout2, rin2;
  double              routbyrin, logrbyr;
  double              expected;
  sc_array_t         *points;
  int                 test_rays;
  int                 skip_1;

  /* data for the currently active quadrant */
  p4est_locidx_t      which_tree;
  p4est_quadrant_t   *sq;
  int                 is_leaf;
  double              width, radius2, R2, r2;
  double              qref[P4EST_DIM], center[P4EST_DIM];
  double              qref2[P4EST_DIM];
  double              bbox[3][2];
  double              planes[2][3][3];  /* (x,y) * (min,mid,max) * dir(3) */
  double             *ray_lengths;
  double             *ray_lengths_expected;
}
tsearch_global_t;

static int
refine_local_fn (p4est_t * p4est, p4est_topidx_t which_tree,
                 p4est_quadrant_t * q)
{
  if ((int) q->level >= refine_level) {
    return 0;
  }
  if ((int) q->level < refine_level - level_shift) {
    return 1;
  }

  return which_tree == 0 && p4est_quadrant_child_id (q) == 0;
}

static int
refine_fractal_fn (p4est_t * p4est, p4est_topidx_t which_tree,
                   p4est_quadrant_t * q)
{
  int                 qid;

  if ((int) q->level >= refine_level) {
    return 0;
  }
  if ((int) q->level < refine_level - level_shift) {
    return 1;
  }

  qid = p4est_quadrant_child_id (q);
  return (qid == 0 || qid == 3
#ifdef P4_TO_P8
          || qid == 5 || qid == 6
#endif
    );
}

typedef struct tsearch_point
{
  double              xy[P4EST_DIM];
  p4est_locidx_t      lid;
}
tsearch_point_t;

/* *INDENT-OFF* */
/* lookup table for the permutations of the coordinate axes in the shell */
static const int    rtop[6][3] = {{ 2, 0, 1 },
                                  { 1, 0, 2 },
                                  { 2, 0, 1 },
                                  { 1, 0, 2 },
                                  { 0, 2, 1 },
                                  { 0, 2, 1 }};
static const double rsign[6][3] = {{  1., -1., -1. },
                                   { -1., -1., -1. },
                                   { -1., -1.,  1. },
                                   {  1., -1.,  1. },
                                   { -1.,  1.,  1. },
                                   {  1., -1.,  1. }};
/* *INDENT-ON* */

static void
map_planes (const double *ref, p4est_topidx_t which_tree,
            double planeeast[3], double planenorth[3])
{
  double              qpole[3];
  double              pole[3];
  double              qeast[3];
  double              east[3];
  double              qnorth[3];
  double              north[3];
  int                 j;
  const int          *pi;
  const double       *si;
  double              sinx, cosx;
  double              siny, cosy;

  /* assert that input coordinates are in the expected range */
  P4EST_ASSERT (ref[0] < 1.0 + 1e-12 && ref[0] > -1.0 - 1e-12);
  P4EST_ASSERT (ref[1] < 1.0 + 1e-12 && ref[1] > -1.0 - 1e-12);
  P4EST_ASSERT (ref[2] < 1.0 + 1e-12 && ref[2] > -1e-12);

  qpole[0] = 0.;
  qpole[1] = 0.;
  qpole[2] = 1.;
  qeast[0] = 1.;
  qeast[1] = 0.;
  qeast[2] = 0.;
  qnorth[0] = 0.;
  qnorth[1] = 1.;
  qnorth[2] = 0.;

  pi = rtop[which_tree / 4];
  si = rsign[which_tree / 4];
  for (j = 0; j < P4EST_DIM; ++j) {
    pole[j] = si[j] * qpole[pi[j]];
    east[j] = si[j] * qeast[pi[j]];
    north[j] = si[j] * qnorth[pi[j]];
  }

  sinx = sin (ref[0] * M_PI_4);
  cosx = cos (ref[0] * M_PI_4);
  siny = sin (ref[1] * M_PI_4);
  cosy = cos (ref[1] * M_PI_4);

  for (j = 0; j < P4EST_DIM; j++) {
    planeeast[j] = cosx * east[j] - sinx * pole[j];
    planenorth[j] = cosy * north[j] - siny * pole[j];
  }
}

static void
reference_to_physical (tsearch_global_t * tsg, const double *ref,
                       double *phys)
{
  const int          *pi;
  const double       *si;
  int                 j;
  double              qvalues[6];
  double              x, y;
  double              R, q;

  /* assert that input coordinates are in the expected range */
  P4EST_ASSERT (ref[0] < 1.0 + 1e-12 && ref[0] > -1.0 - 1e-12);
  P4EST_ASSERT (ref[1] < 1.0 + 1e-12 && ref[1] > -1.0 - 1e-12);
  P4EST_ASSERT (ref[2] < 1.0 + 1e-12 && ref[2] > -1e-12);

  /* transform x and y for nicer grading */
  x = tan (ref[0] * M_PI_4);
  y = tan (ref[1] * M_PI_4);

  /* compute transformation ingredients */
  R = tsg->rin * pow (tsg->routbyrin, ref[2]);
  q = R / sqrt (x * x + y * y + 1.);

  /* compute physical coordinate values */
  qvalues[0] = q * x;
  qvalues[1] = q * y;
  qvalues[2] = q;

  /* assign coordinates based on tree id */
#if 1
  pi = rtop[tsg->which_tree / 4];
  si = rsign[tsg->which_tree / 4];
  for (j = 0; j < P4EST_DIM; ++j) {
    phys[j] = si[j] * qvalues[pi[j]];
  }
#else
  switch (tsg->which_tree / 4) {
  case 3:                      /* top */
    phys[0] = +q * y;
    phys[1] = -q * x;
    phys[2] = +q;
    break;
  case 2:                      /* left */
    phys[0] = -q;
    phys[1] = -q * x;
    phys[2] = +q * y;
    break;
  case 1:                      /* bottom */
    phys[0] = -q * y;
    phys[1] = -q * x;
    phys[2] = -q;
    break;
  case 0:                      /* right */
    phys[0] = +q;
    phys[1] = -q * x;
    phys[2] = -q * y;
    break;
  case 4:                      /* back */
    phys[0] = -q * x;
    phys[1] = +q;
    phys[2] = +q * y;
    break;
  case 5:                      /* front */
    phys[0] = +q * x;
    phys[1] = -q;
    phys[2] = +q * y;
    break;
  default:
    SC_ABORT_NOT_REACHED ();
  }
#endif
}

/** Return true if and only if the reference coordinates match the tree. */
static int
physical_to_reference (tsearch_global_t * tsg, const double *phys,
                       double *ref)
{
  const int          *pi;
  const double       *si;
  int                 j;
  double              q, R2;
  double              qvalues[3];
  const p4est_topidx_t tid = tsg->which_tree;

  /* reverse permutation of coordinate axes and compute radius */
  R2 = 0.;
  pi = rtop[tid / 4];
  si = rsign[tid / 4];
  for (j = 0; j < P4EST_DIM; ++j) {
    q = qvalues[pi[j]] = si[j] * phys[j];
    R2 += q * q;
  }

  /* recover the reference coordinates */
  for (j = 0; j < 2; ++j) {
    q = ref[j] = atan2 (qvalues[j], qvalues[2]) / M_PI_4;
    /* check x and y coordinate range */
    if (!(tid & (1 << j))) {
      q += 1.;
    }
    if (q < 0. || q > 1.) {
      return 0;
    }
  }
  q = ref[2] = log (sqrt (R2) / tsg->rin) / tsg->logrbyr;

  return q >= 0. && q <= 1.;
}

/** Rays **/
typedef struct ts_ray
{
  double              dir[3];   /* normal vector */
  double              off[3];   /* point on the ray closest to the origin */
  /* a ray is parametrized as off + s * dir */
  double              idir[3];  /* 1. / dir */
  double              off2;     /* off' * off */
  /* by choosing off to be the closest point to the origin,
   * we guarantee off' * dir = 0 */
  /* by chooseing dir to be normal, dir' * dir = 1 */
  p4est_locidx_t      lidx;
}
ts_ray_t;

static int
ts_ray_sphere_test (ts_ray_t * ray, double R2, double minmax[2])
{
  double              o2 = ray->off2;
  double              val;

  if (o2 > R2) {
    return 0;
  }
  else {
    val = sqrt (R2 - o2);
    minmax[0] = -val;
    minmax[1] = val;

    return 1;
  }
}

static int
ts_ray_bounding_box_test (ts_ray_t * ray, double bounds[3][2],
                          double minmax[2])
{
  double              o[3] = { ray->off[0], ray->off[1], ray->off[2] };
  double              id[3] = { ray->idir[0], ray->idir[1], ray->idir[2] };
  double              b[3][2] = { {bounds[0][0], bounds[0][1]},
  {bounds[1][0], bounds[1][1]},
  {bounds[2][0], bounds[2][1]}
  };
  double              s[3][2];
  double              min, max;
  double              dmin, dmax;
  int                 i, j;

  min = minmax[0];
  max = minmax[1];

  for (i = 0; i < 3; i++) {
    /* intercepts of parameter s */
    for (j = 0; j < 2; j++) {
      s[i][j] = (b[i][j] - o[i]) * id[i];
    }
    /* overlap interval */
    if (s[i][0] < s[i][1]) {
      dmin = s[i][0];
      dmax = s[i][1];
    }
    else {
      dmin = s[i][1];
      dmax = s[i][0];
    }
    /* if overlap interval does not intersect
     * previous overlap interval */
    if (min > dmax || dmin > max) {
      return 0;
    }
    min = SC_MAX (min, dmin);
    max = SC_MIN (max, dmax);
  }

  minmax[0] = min;
  minmax[1] = max;

  return 1;
}

static void
tsearch_setup (tsearch_global_t * tsg)
{
  int                 j, k;
  double              mlen, hwidth, dist2;
  double              ref[P4EST_DIM], phys[P4EST_DIM];
  double             *qref = tsg->qref;
  double             *qref2 = tsg->qref2;
  double              cref[3];
  double             *center = tsg->center;
  p4est_quadrant_t   *q = tsg->sq, c;
#ifdef P4EST_ENABLE_DEBUG
  int                 retval;
  double              nref[P4EST_DIM];
#endif
  double              bbox[3][2] = { {1. / 0., -1. / 0.},
  {1. / 0., -1. / 0.},
  {1. / 0., -1. / 0.}
  };

  mlen = 1. / P4EST_ROOT_LEN;
  hwidth = .5 * (tsg->width = mlen * P4EST_QUADRANT_LEN (q->level));

  /* transform center of the quadrant to physical space */
  ref[0] = cref[0] = q->x * mlen + hwidth - 1. + (tsg->which_tree & 1);
  ref[1] = cref[1] = q->y * mlen + hwidth - 1. + (tsg->which_tree & 2) / 2;
  ref[2] = cref[2] = q->z * mlen + hwidth;
  reference_to_physical (tsg, ref, center);

#ifdef P4EST_ENABLE_DEBUG
  retval = physical_to_reference (tsg, center, nref);
  P4EST_ASSERT (retval);
  for (j = 0; j < P4EST_DIM; ++j) {
    P4EST_ASSERT (fabs (ref[j] - nref[j]) < 1e-10);
  }
#endif

  /* transform all corners of the quadrant and take max distance to center */
  if (!tsg->test_rays) {
    tsg->radius2 = 0.;
    for (k = 0; k < P4EST_CHILDREN; ++k) {
      p4est_quadrant_corner_node (q, k, &c);
      ref[0] = c.x * mlen - 1. + (tsg->which_tree & 1);
      ref[1] = c.y * mlen - 1. + (tsg->which_tree & 2) / 2;
      ref[2] = c.z * mlen;
      reference_to_physical (tsg, ref, phys);
      dist2 = 0.;
      for (j = 0; j < P4EST_DIM; ++j) {
        dist2 += (phys[j] - center[j]) * (phys[j] - center[j]);
      }
      if (dist2 > tsg->radius2) {
        tsg->radius2 = dist2;
      }
      if (k == 0) {
        qref[0] = ref[0];
        qref[1] = ref[1];
        qref[2] = ref[2];
      }
    }
  }
  else {
    tsg->R2 = 0.;
    tsg->r2 = 1. / 0.;
    for (k = 0; k < P4EST_CHILDREN; ++k) {
      p4est_quadrant_corner_node (q, k, &c);
      ref[0] = c.x * mlen - 1. + (tsg->which_tree & 1);
      ref[1] = c.y * mlen - 1. + (tsg->which_tree & 2) / 2;
      ref[2] = c.z * mlen;
      reference_to_physical (tsg, ref, phys);
      for (j = 0; j < P4EST_DIM; ++j) {
        bbox[j][0] = SC_MIN (bbox[j][0], phys[j]);
        bbox[j][1] = SC_MAX (bbox[j][1], phys[j]);
      }
      tsg->R2 =
        SC_MAX (tsg->R2,
                phys[0] * phys[0] + phys[1] * phys[1] + phys[2] * phys[2]);
      tsg->r2 =
        SC_MIN (tsg->r2,
                phys[0] * phys[0] + phys[1] * phys[1] + phys[2] * phys[2]);
      if (!k) {
        qref[0] = ref[0];
        qref[1] = ref[1];
        qref[2] = ref[2];
      }
      else if (k == P4EST_CHILDREN - 1) {
        qref2[0] = ref[0];
        qref2[1] = ref[1];
        qref2[2] = ref[2];
      }
    }
    for (j = 0; j < 3; j++) {
      tsg->bbox[j][0] = bbox[j][0];
      tsg->bbox[j][1] = bbox[j][1];
    }
    if (tsg->is_leaf) {
      map_planes (qref, tsg->which_tree, &tsg->planes[0][0][0],
                  &tsg->planes[1][0][0]);
      map_planes (cref, tsg->which_tree, &tsg->planes[0][1][0],
                  &tsg->planes[1][1][0]);
      map_planes (qref2, tsg->which_tree, &tsg->planes[0][2][0],
                  &tsg->planes[1][2][0]);
    }
  }
}

static int
time_search_fn (p4est_t * p4est, p4est_topidx_t which_tree,
                p4est_quadrant_t * q, p4est_locidx_t local_num, void *point)
{
  tsearch_global_t   *tsg = (tsearch_global_t *) p4est->user_pointer;
  int                 i, j;
  double              width, r2;
  double              ref[P4EST_DIM];
  const double       *qref;

  if (point == NULL) {
    /* per-quadrant setup function */
    tsg->which_tree = which_tree;
    tsg->sq = q;
    tsg->is_leaf = local_num >= 0;
    tsearch_setup (tsg);
    return 1;
  }
  P4EST_ASSERT (tsg->sq == q);
  P4EST_ASSERT (tsg->is_leaf == (local_num >= 0));

  if (!tsg->test_rays) {
    tsearch_point_t    *t = (tsearch_point_t *) point;

    /* root level check to see if the point is contained in the shell */
    if (q->level == 0) {
      r2 = t->xy[0] * t->xy[0] + t->xy[1] * t->xy[1] + t->xy[2] * t->xy[2];
      if (r2 >= tsg->rout2 || r2 <= tsg->rin2) {
        return 0;
      }
    }

    if (!tsg->is_leaf) {
      /* perform over-optimistic check with the quadrant's bounding sphere */
      r2 = 0.;
      for (j = 0; j < P4EST_DIM; ++j) {
        r2 += (t->xy[j] - tsg->center[j]) * (t->xy[j] - tsg->center[j]);
      }
      return r2 < (1. + 1e-12) * tsg->radius2;
    }
    else {
      /* perform strict check by inverse coordinate transformation */
      if (physical_to_reference (tsg, t->xy, ref)) {
        /* the point is contained in the correct tree, now check quadrant */
        width = tsg->width;
        qref = tsg->qref;
        for (j = 0; j < P4EST_DIM; ++j) {
          if (ref[j] < qref[j] || ref[j] > qref[j] + width) {
            return 0;
          }
        }

        /* we have found a matching quadrant for this point */
        ++tsg->matches;
        t->lid = local_num;
        return 1;
      }

      /* the return value is irrelevant for leaves */
      return 0;
    }
  }
  else {
    ts_ray_t           *ray = (ts_ray_t *) point;
    double              minmax[2];
    double              innerminmax[2];
    double              o[3] = { ray->off[0], ray->off[1], ray->off[2] };
    double              d[3] = { ray->dir[0], ray->dir[1], ray->dir[2] };
    double              length;

    /* quick tests */
    if (!ts_ray_sphere_test (ray, tsg->R2, minmax) ||
        !ts_ray_bounding_box_test (ray, tsg->bbox, minmax)) {
      return 0;
    }

    if (!tsg->is_leaf) {
      return 1;
    }

    /* intersect with the planes that define the vertical boundaries of the
     * physical element */
    for (i = 0; i < 2; i++) {
      double              s[3];
      double              dmin, dmax;

      for (j = 0; j < 3; j++) {
        double             *plane = &tsg->planes[i][j][0];
        double              offdotplane;
        double              dirdotplane;

        offdotplane = plane[0] * o[0] + plane[1] * o[1] + plane[2] * o[2];
        dirdotplane = plane[0] * d[0] + plane[1] * d[1] + plane[2] * d[2];

        s[j] = -offdotplane / dirdotplane;
      }
      if (s[0] < s[2]) {
        dmin = s[0];
        dmax = s[2];
      }
      else {
        dmin = s[2];
        dmax = s[0];
      }
      if (s[1] < dmin || s[1] > dmax) {
        if (dmin <= minmax[0] && minmax[1] <= dmax) {
          return 0;
        }

        if (dmax < minmax[1]) {
          minmax[0] = SC_MAX (minmax[0], dmax);
        }
        if (dmin > minmax[0]) {
          minmax[1] = SC_MIN (minmax[1], dmin);
        }
      }
      else {
        if (minmax[0] > dmax || dmin > minmax[1]) {
          return 0;
        }
        minmax[0] = SC_MAX (dmin, minmax[0]);
        minmax[1] = SC_MIN (dmax, minmax[1]);
      }
    }

    length = minmax[1] - minmax[0];

    if (length && ts_ray_sphere_test (ray, tsg->r2, innerminmax)) {
      innerminmax[0] = SC_MAX (innerminmax[0], minmax[0]);
      innerminmax[1] = SC_MIN (innerminmax[1], minmax[1]);

      length -= SC_MAX (0., innerminmax[1] - innerminmax[0]);
    }

    if (length > 0.) {
      tsg->ray_lengths[ray->lidx] += length;
      return 1;
    }

    /* the return value is irrelevant for leaves */
    return 0;
  }
}

static void
time_search_1 (tsearch_global_t * tsg, p4est_t * p4est, size_t znum_points,
               sc_flopinfo_t * fi, sc_statinfo_t * stats)
{
  const double        expected = tsg->expected;
  int                 mpiret;
  long long           ll, gg;
  size_t              zz;
  sc_array_t          pview;
  sc_flopinfo_t       snapshot;

  /*
   * Search all points separately.
   *
   * For each point, perform a search to compute the local number of the
   * containing quadrant and the reference coordinates in [0, 1]^d relative to
   * this quadrant in-place.  This only happens if the quadrant is
   * processor-local.
   *
   * The points are identical on all processors.  Due to points on quadrant
   * boundaries, we will find slightly more points than we expect.
   */

  if (tsg->test_rays) {
    memset (tsg->ray_lengths, 0, znum_points * sizeof (double));
  }

  ll = 0;
  sc_flops_snap (fi, &snapshot);
  for (zz = 0; zz < znum_points; ++zz) {
    tsg->matches = 0;
    sc_array_init_view (&pview, tsg->points, zz, 1);
    p4est_search (p4est, time_search_fn, time_search_fn, &pview);
    ll += (long long) tsg->matches;
  }
  sc_flops_shot (fi, &snapshot);
  sc_stats_set1 (&stats[TSEARCH_SEARCH_1], snapshot.iwtime, "Search_1");

  if (!tsg->test_rays) {
    mpiret = sc_MPI_Allreduce (&ll, &gg, 1, sc_MPI_LONG_LONG_INT, sc_MPI_SUM,
                               tsg->mpicomm);
    SC_CHECK_MPI (mpiret);

    P4EST_GLOBAL_STATISTICSF
      ("Search_1 expected %lld found %lld of %lld error %.3g%%\n",
       (long long) round (expected), gg, (long long) znum_points,
       100. * fabs ((gg - expected) / expected));
  }
  else {
    double             *ray_lengths_buf = NULL;
    sc_statinfo_t       raystat;

    if (!p4est->mpirank) {
      ray_lengths_buf = P4EST_ALLOC_ZERO (double, znum_points);
    }

    mpiret = sc_MPI_Reduce (tsg->ray_lengths, ray_lengths_buf,
                            (int) znum_points, sc_MPI_DOUBLE, sc_MPI_SUM, 0,
                            p4est->mpicomm);
    SC_CHECK_MPI (mpiret);

    sc_stats_init (&raystat, "Search_1 ray length errors");

    if (!p4est->mpirank) {
      for (zz = 0; zz < znum_points; zz++) {
        double              diff;

        diff = ray_lengths_buf[zz] - tsg->ray_lengths_expected[zz];
        sc_stats_accumulate (&raystat, diff);
      }

      P4EST_FREE (ray_lengths_buf);
    }

    sc_stats_compute (tsg->mpicomm, 1, &raystat);
    sc_stats_print (p4est_package_id, SC_LP_STATISTICS, 1, &raystat, 1, 0);
  }
}

static void
time_search_N (tsearch_global_t * tsg, p4est_t * p4est, size_t znum_points,
               sc_flopinfo_t * fi, sc_statinfo_t * stats)
{
  const double        expected = tsg->expected;
  int                 mpiret;
  long long           ll, gg;
  sc_flopinfo_t       snapshot;

  /*
   * Search all points in one pass through the forest
   *
   * For each point, perform a search to compute the local number of the
   * containing quadrant and the reference coordinates in [0, 1]^d relative to
   * this quadrant in-place.  This only happens if the quadrant is
   * processor-local.
   *
   * The points are identical on all processors.  Due to points on quadrant
   * boundaries, we will find slightly more points than we expect.
   */

  if (tsg->test_rays) {
    memset (tsg->ray_lengths, 0, znum_points * sizeof (double));
  }

  tsg->matches = 0;
  sc_flops_snap (fi, &snapshot);
  p4est_search (p4est, time_search_fn, time_search_fn, tsg->points);
  sc_flops_shot (fi, &snapshot);
  sc_stats_set1 (&stats[TSEARCH_SEARCH_N], snapshot.iwtime, "Search_N");
  ll = (long long) tsg->matches;

  if (!tsg->test_rays) {
    mpiret = sc_MPI_Allreduce (&ll, &gg, 1, sc_MPI_LONG_LONG_INT, sc_MPI_SUM,
                               tsg->mpicomm);
    SC_CHECK_MPI (mpiret);

    P4EST_GLOBAL_STATISTICSF
      ("Search_N expected %lld found %lld of %lld error %.3g%%\n",
       (long long) round (expected), gg, (long long) znum_points,
       100. * fabs ((gg - expected) / expected));
  }
  else {
    double             *ray_lengths_buf = NULL;
    sc_statinfo_t       raystat;

    if (!p4est->mpirank) {
      ray_lengths_buf = P4EST_ALLOC_ZERO (double, znum_points);
    }

    mpiret = sc_MPI_Reduce (tsg->ray_lengths, ray_lengths_buf,
                            (int) znum_points, sc_MPI_DOUBLE, sc_MPI_SUM, 0,
                            p4est->mpicomm);
    SC_CHECK_MPI (mpiret);

    sc_stats_init (&raystat, "Search_N ray length errors");

    if (!p4est->mpirank) {
      size_t              zz;

      for (zz = 0; zz < znum_points; zz++) {
        double              diff;

        diff = ray_lengths_buf[zz] - tsg->ray_lengths_expected[zz];
        sc_stats_accumulate (&raystat, diff);
      }

      P4EST_FREE (ray_lengths_buf);
    }

    sc_stats_compute (tsg->mpicomm, 1, &raystat);
    sc_stats_print (p4est_package_id, SC_LP_STATISTICS, 1, &raystat, 1, 0);
  }
}

static void
time_search_all (tsearch_global_t * tsg, p4est_t * p4est, size_t znum_points,
                 sc_flopinfo_t * fi, sc_statinfo_t * stats)
{
  int                 i;
  double              ratio;
  size_t              zz;
  tsearch_point_t    *point;
  ts_ray_t           *ray;

  if (!tsg->test_rays) {
    tsg->points = sc_array_new_size (sizeof (tsearch_point_t), znum_points);
    for (zz = 0; zz < znum_points; ++zz) {
      point = (tsearch_point_t *) sc_array_index (tsg->points, zz);
      for (i = 0; i < P4EST_DIM; ++i) {
        point->xy[i] = 2. * (rand () / (RAND_MAX + 1.)) - 1.;
      }
      point->lid = -1;
    }

    ratio = 4. / 3. * M_PI * (pow (tsg->rout, 3.) - pow (tsg->rin, 3.)) / 8.;
    tsg->expected = ratio * znum_points;
  }
  else {
    tsg->points = sc_array_new_size (sizeof (ts_ray_t), znum_points);
    tsg->ray_lengths = P4EST_ALLOC (double, znum_points);
    tsg->ray_lengths_expected = P4EST_ALLOC (double, znum_points);
    for (zz = 0; zz < znum_points; ++zz) {
      double              theta, sintheta, costheta, sinphi, cosphi, r;
      double              x[3], y[3], d[3];
      double              minmax[2];

      ray = (ts_ray_t *) sc_array_index (tsg->points, zz);

      /* generate a random direction */
      theta = 2. * M_PI * ((double) rand () / (double) (RAND_MAX));
      sintheta = sin (theta);
      costheta = cos (theta);

      sinphi = 2. * ((double) rand () / (double) (RAND_MAX)) - 1.;
      cosphi = sqrt (1. - sinphi * sinphi);

      ray->dir[0] = d[0] = costheta * cosphi;
      ray->dir[1] = d[1] = sintheta * cosphi;
      ray->dir[2] = d[2] = sinphi;
      ray->idir[0] = 1. / d[0];
      ray->idir[1] = 1. / d[1];
      ray->idir[2] = 1. / d[2];

      /* get a basis for the plane orthogonal to the direction */
      if (d[0] * d[0] > d[1] * d[1] && d[0] * d[0] > d[2] * d[2]) {
        x[0] = d[1];
        x[1] = -d[0];
        x[2] = 0.;
      }
      else {
        x[0] = 0.;
        x[1] = d[2];
        x[2] = -d[1];
      }

      r = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
      r = sqrt (r);
      x[0] /= r;
      x[1] /= r;
      x[2] /= r;

      y[0] = d[1] * x[2] - d[2] * x[1];
      y[1] = d[2] * x[0] - d[0] * x[2];
      y[2] = d[0] * x[1] - d[1] * x[0];

      /* generate a random point in the circle orthogonal to the direction */
      theta = 2. * M_PI * ((double) rand () / (double) (RAND_MAX));
      sintheta = sin (theta);
      costheta = cos (theta);
      r =
        ((double) rand () / (double) (RAND_MAX)) +
        ((double) rand () / (double) (RAND_MAX));
      if (r > 1.) {
        r = 2. - r;
      }
      r *= tsg->rout;

      ray->off[0] = x[0] * r * costheta + y[0] * r * sintheta;
      ray->off[1] = x[1] * r * costheta + y[1] * r * sintheta;
      ray->off[2] = x[2] * r * costheta + y[2] * r * sintheta;
      ray->off2 = r * r;

      ray->lidx = (p4est_locidx_t) zz;

      if (ts_ray_sphere_test (ray, tsg->rout * tsg->rout, minmax)) {
        double              length = minmax[1] - minmax[0];
        double              iminmax[2];

        if (ts_ray_sphere_test (ray, tsg->rin * tsg->rin, iminmax)) {
          length -= (iminmax[1] - iminmax[0]);
        }
        tsg->ray_lengths_expected[zz] = length;
      }
      else {
        tsg->ray_lengths_expected[zz] = 0.;
      }
    }
  }

  if (!tsg->skip_1) {
    time_search_1 (tsg, p4est, znum_points, fi, stats);
  }
  else {
    sc_stats_set1 (&stats[TSEARCH_SEARCH_1], 0., "Search_1");
  }
  time_search_N (tsg, p4est, znum_points, fi, stats);

  sc_array_destroy (tsg->points);
  if (tsg->test_rays) {
    P4EST_FREE (tsg->ray_lengths);
    P4EST_FREE (tsg->ray_lengths_expected);
  }
}

int
main (int argc, char **argv)
{
  int                 mpiret;
  int                 first_argc;
  int                 refine_local;
  int                 write_vtk;
  int                 test_rays;
  size_t              znum_points;
  unsigned            crc;
  p4est_gloidx_t      count_refined, count_balanced;
  p4est_connectivity_t *connectivity;
  p4est_t            *p4est;
  sc_statinfo_t       stats[TSEARCH_NUM_STATS];
  sc_flopinfo_t       fi, snapshot;
  sc_options_t       *opt;
  tsearch_global_t    tsgt, *tsg = &tsgt;
  int                 skip_1;

  /* initialize MPI */
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  tsg->mpicomm = sc_MPI_COMM_WORLD;
  mpiret = sc_MPI_Comm_size (tsg->mpicomm, &tsg->mpisize);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (tsg->mpicomm, &tsg->mpirank);
  SC_CHECK_MPI (mpiret);

  /* initialize p4est internals */
  sc_init (tsg->mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
#ifndef P4EST_ENABLE_DEBUG
  sc_set_log_defaults (NULL, NULL, SC_LP_STATISTICS);
#endif
  p4est_init (NULL, SC_LP_DEFAULT);

  /* initialize global data */
  srand (tsg->mpisize);
  tsg->rout = 1.;
  tsg->rin = .55;
  tsg->rout2 = tsg->rout * tsg->rout;
  tsg->rin2 = tsg->rin * tsg->rin;
  tsg->routbyrin = tsg->rout / tsg->rin;
  tsg->logrbyr = log (tsg->routbyrin);

  /* process command line arguments */
  opt = sc_options_new (argv[0]);
  sc_options_add_int (opt, 'l', "level", &refine_level, 0,
                      "Refinement level");
  sc_options_add_int (opt, 's', "level-shift", &level_shift, 4,
                      "Refinement shift");
  sc_options_add_bool (opt, 'L', "refine-local", &refine_local, 0,
                       "Refine around one point");
  sc_options_add_bool (opt, 'V', "write-vtk", &write_vtk, 0,
                       "Write VTK files");
  sc_options_add_size_t (opt, 'N', "num-points", &znum_points, 0,
                         "Number of points");
  sc_options_add_switch (opt, 'r', "rays", &test_rays,
                         "Test rays instead of points");
  sc_options_add_switch (opt, 'm', "skip-1", &skip_1,
                         "only test parallelized search");
  first_argc = sc_options_parse (p4est_package_id, SC_LP_ERROR,
                                 opt, argc, argv);
  if (first_argc < 0 || first_argc != argc) {
    sc_options_print_usage (p4est_package_id, SC_LP_ERROR, opt, NULL);
    return 1;
  }
  sc_options_print_summary (p4est_package_id, SC_LP_PRODUCTION, opt);

  tsg->test_rays = test_rays;
  tsg->skip_1 = skip_1;

  /* start overall timing */
  mpiret = sc_MPI_Barrier (tsg->mpicomm);
  SC_CHECK_MPI (mpiret);
  sc_flops_start (&fi);

  /* create connectivity and forest */
  sc_flops_snap (&fi, &snapshot);
  connectivity = p8est_connectivity_new_shell ();
  p4est = p4est_new_ext (tsg->mpicomm, connectivity,
                         0, refine_level - level_shift, 1, 0, NULL, tsg);
  sc_flops_shot (&fi, &snapshot);
  sc_stats_set1 (&stats[TSEARCH_NEW], snapshot.iwtime, "New");
  if (write_vtk) {
    p4est_vtk_write_file (p4est, NULL, "tsearch_new");
  }

  /* time refine */
  sc_flops_snap (&fi, &snapshot);
  p4est_refine (p4est, 1,
                refine_local ? refine_local_fn : refine_fractal_fn, NULL);
  sc_flops_shot (&fi, &snapshot);
  sc_stats_set1 (&stats[TSEARCH_REFINE], snapshot.iwtime, "Refine");
  if (write_vtk) {
    p4est_vtk_write_file (p4est, NULL, "tsearch_refined");
  }
  count_refined = p4est->global_num_quadrants;

  /* time balance */
  sc_flops_snap (&fi, &snapshot);
  p4est_balance (p4est, P4EST_CONNECT_FULL, NULL);
  sc_flops_shot (&fi, &snapshot);
  sc_stats_set1 (&stats[TSEARCH_BALANCE], snapshot.iwtime, "Balance");
  if (write_vtk) {
    p4est_vtk_write_file (p4est, NULL, "tsearch_balanced");
  }
  count_balanced = p4est->global_num_quadrants;
  crc = p4est_checksum (p4est);

  /* time a uniform partition */
  sc_flops_snap (&fi, &snapshot);
  p4est_partition (p4est, 0, NULL);
  sc_flops_shot (&fi, &snapshot);
  sc_stats_set1 (&stats[TSEARCH_PARTITION], snapshot.iwtime, "Partition");
  if (write_vtk) {
    p4est_vtk_write_file (p4est, NULL, "tsearch_partitioned");
  }
  P4EST_ASSERT (crc == p4est_checksum (p4est));

  /* run search timings */
  time_search_all (tsg, p4est, znum_points, &fi, stats);

  /* print status and checksum */
  P4EST_GLOBAL_STATISTICSF
    ("Processors %d level %d shift %d points %llu checksum 0x%08x\n",
     tsg->mpisize, refine_level, level_shift,
     (unsigned long long) znum_points, crc);
  P4EST_GLOBAL_STATISTICSF ("Level %d refined to %lld balanced to %lld\n",
                            refine_level, (long long) count_refined,
                            (long long) count_balanced);

  /* calculate and print timings */
  sc_stats_compute (tsg->mpicomm, TSEARCH_NUM_STATS, stats);
  sc_stats_print (p4est_package_id, SC_LP_STATISTICS,
                  TSEARCH_NUM_STATS, stats, 1, 1);

  /* destroy the p4est and its connectivity structure */
  p4est_destroy (p4est);
  p4est_connectivity_destroy (connectivity);

  /* clean up and exit */
  sc_options_destroy (opt);

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
