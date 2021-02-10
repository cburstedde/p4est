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

#ifndef P4_TO_P8
#include <p4est_bits.h>
#include <p4est_build.h>
#include <p4est_communication.h>
#include <p4est_extended.h>
#include <p4est_search.h>
#include <p4est_vtk.h>
#else
#include <p8est_bits.h>
#include <p8est_build.h>
#include <p8est_communication.h>
#include <p8est_extended.h>
#include <p8est_search.h>
#include <p8est_vtk.h>
#endif /* P4_TO_P8 */
#include <sc_notify.h>
#include <sc_options.h>
#include "particles_global.h"

#define PARTICLES_xstr(s) PARTICLES_str(s)
#define PARTICLES_str(s) #s
#define PARTICLES_48() PARTICLES_xstr(P4EST_CHILDREN)

/** Send full particle information in first message, comment out if not */
#define PART_SENDFULL

/** Context data to compute initial particle positions */
typedef struct pi_data
{
  double              sigma;
  double              invs2;
  double              gnorm;
  double              center[3];
}
pi_data_t;

/** Data type for payload data inside each quadrant */
typedef struct qu_data
{
  union
  {
    /** Offset into local array of all particles after this quadrant */
    p4est_locidx_t      lpend;
    double              d;
  } u;

  /** counts of local particles remaining on this quadrant and recieved ones */
  p4est_locidx_t      premain, preceive;
}
qu_data_t;

/** Property data stored in a flat array over all particles */
typedef struct pa_data
{
  double              xv[6];
  double              wo[6];
  double              up[6];
  double              rm[3];
  p4est_gloidx_t      id;
}
pa_data_t;

#if 0
static void
ppad (part_global_t * g, const char *lead)
{
  size_t              zz;
  pa_data_t          *pad;

  P4EST_ASSERT (g != NULL);
  P4EST_ASSERT (g->padata != NULL);

  if (lead != NULL) {
    P4EST_VERBOSEF ("PPAD %s\n", lead);
  }

  for (zz = 0; zz < g->padata->elem_count; ++zz) {
    pad = (pa_data_t *) sc_array_index (g->padata, zz);
    P4EST_VERBOSEF ("PPAD L %d I %d\n", (int) zz, (int) pad->id);
  }
}

static void
lrem (part_global_t * g, const char *lead)
{
  size_t              zz;
  p4est_locidx_t      lrem;

  P4EST_ASSERT (g != NULL);
  P4EST_ASSERT (g->iremain != NULL);

  if (lead != NULL) {
    P4EST_VERBOSEF ("LREM %s\n", lead);
  }

  for (zz = 0; zz < g->iremain->elem_count; ++zz) {
    lrem = *(p4est_locidx_t *) sc_array_index (g->iremain, zz);
    P4EST_VERBOSEF ("LREM Z %d N %d\n", (int) zz, (int) lrem);
  }
}
#endif

#ifdef PART_SENDFULL
#define PART_MSGSIZE (sizeof (pa_data_t))
#else
#define PART_MSGSIZE (3 * sizeof (double))
#endif

/** Hash table entry for a process that we send messages to */
typedef struct comm_psend
{
  int                 rank;
  sc_array_t          message;     /** Message data to send */
}
comm_psend_t;

/** Array entry for a process that we send messages to */
typedef struct comm_prank
{
  int                 rank;
  comm_psend_t       *psend;        /**< Points to hash table entry */
}
comm_prank_t;

typedef enum comm_tag
{
  COMM_TAG_PART = P4EST_COMM_TAG_LAST,
  COMM_TAG_FIXED,
  COMM_TAG_CUSTOM,
  COMM_TAG_LAST
}
comm_tag_t;

static const double planet_xyz[3][3] = {
  {.48, .58, .59},
  {.58, .41, .46},
  {.51, .52, .42}
};
static const double planet_mass[3] = { .049, .167, .06 };

static const double pidensy_sigma = .07;
static const double pidensy_center[3] = { .3, .4, .5 };

#define PART_NQPOINTS 3
static double       qpoints[PART_NQPOINTS];
static double       qweights[PART_NQPOINTS];

static const double rk1b[1] = { 0. };   /* avoid -pedantic warning for [0] */
static const double rk1g[1] = { 1. };
static const double rk2b[1] = { 1. };
static const double rk2g[2] = { .5, .5 };
static const double rk3b[2] = { 1. / 3., 2. / 3. };
static const double rk3g[3] = { .25, 0., .75 };
static const double rk4b[3] = { .5, .5, 1. };
static const double rk4g[4] = { 1. / 6., 1. / 3., 1. / 3., 1. / 6. };

static const double *prk[4][2] = {
  {rk1b, rk1g},
  {rk2b, rk2g},
  {rk3b, rk3g},
  {rk4b, rk4g}
};

static const char  *snames[PART_STATS_LAST] = {
  "Notify",
  "Binary",
  "Nary",
  "Comm",
  "Wait_A",
  "Wait_B",
  "Num_Peers",
  "Search_P",
  "Search_L",
  "Trans_F",
  "Trans_C",
  "Build",
  "Pertree_F",
  "Pertree_B",
  "RK"
};

#if 0

static void
p4est_free_int (int **pptr)
{
  P4EST_ASSERT (pptr != NULL);
  P4EST_FREE (*pptr);
  *pptr = NULL;
}

#endif

static void
part_string_to_int (const char *str, int n, ...)
{
  int                 i, j;
  int                *pi;
  char                buf[BUFSIZ];
  va_list             ap;

  P4EST_ASSERT (n >= 0);
  if (str == NULL || n == 0) {
    return;
  }

  /* use a once-loop for clean early return */
  va_start (ap, n);
  do {
    if (n == 1) {
      pi = va_arg (ap, int *);
      *pi = atoi (str);
      break;
    }

    j = snprintf (buf, BUFSIZ, "%s", "%d");
    if (j >= BUFSIZ) {
      break;
    }
    for (i = 1; i < n; ++i) {
      j += snprintf (buf + j, BUFSIZ - j, ":%s", "%d");
      if (j >= BUFSIZ) {
        break;
      }
    }

    /* we ignore return values and rely on defaults */
    vsscanf (str, buf, ap);
  }
  while (0);
  va_end (ap);
}

static void        *
sc_array_index_begin (sc_array_t * arr)
{
  P4EST_ASSERT (arr != NULL);

  if (arr->elem_count == 0) {
    return NULL;
  }

  P4EST_ASSERT (arr->array != NULL);
  return (void *) arr->array;
}

#ifdef P4EST_ENABLE_DEBUG

static void        *
sc_array_index_end (sc_array_t * arr)
{
  P4EST_ASSERT (arr != NULL);

  if (arr->elem_count == 0) {
    return NULL;
  }

  P4EST_ASSERT (arr->array != NULL);
  return (void *) (arr->array + arr->elem_count * arr->elem_size);
}

#endif

/** With two initialized arrays of same metadata, copy array data */
static void
sc_array_paste (sc_array_t * dest, sc_array_t * src)
{
  P4EST_ASSERT (dest->elem_size == src->elem_size);
  P4EST_ASSERT (dest->elem_count == src->elem_count);

  memcpy (dest->array, src->array, src->elem_count * src->elem_size);
}

/** Initialize one array with contents from other, reinit the other */
static void
sc_array_swap_init (sc_array_t * array, sc_array_t * from, size_t elem_size)
{
  *array = *from;
  sc_array_init (from, elem_size);
}

/** Turn statistics collected so far into one value */
static void
sc_stats_collapse (sc_statinfo_t * stats)
{
  double              value;

  SC_ASSERT (stats->dirty);
  if (stats->count) {
    value = stats->sum_values / (double) stats->count;
    stats->sum_values = value;
    stats->sum_squares = value * value;
    stats->min = stats->max = value;
    stats->count = 1;
  }
}

static int
comm_prank_compare (const void *v1, const void *v2)
{
  return sc_int_compare (&((const comm_prank_t *) v1)->rank,
                         &((const comm_prank_t *) v2)->rank);
}

static double
gaussnorm (double sigma)
{
  return pow (2. * M_PI * sigma * sigma, -.5 * P4EST_DIM);
}

static double
pidense (double x, double y, double z, void *data)
{
  pi_data_t          *piddata = (pi_data_t *) data;

  P4EST_ASSERT (piddata != NULL);
  P4EST_ASSERT (piddata->sigma > 0.);
  P4EST_ASSERT (piddata->invs2 > 0.);
  P4EST_ASSERT (piddata->gnorm > 0.);

  return piddata->gnorm * exp (-.5 * (SC_SQR (x - piddata->center[0]) +
                                      SC_SQR (y - piddata->center[1]) +
                                      SC_SQR (z - piddata->center[2])) *
                               piddata->invs2);
}

static void
loopquad (part_global_t * g, p4est_topidx_t tt, p4est_quadrant_t * quad,
          double lxyz[3], double hxyz[3], double dxyz[3])
{

  int                 i;
  p4est_qcoord_t      qh;

  qh = P4EST_QUADRANT_LEN (quad->level);
  p4est_qcoord_to_vertex (g->conn, tt, quad->x, quad->y,
#ifdef P4_TO_P8
                          quad->z,
#endif
                          lxyz);
  p4est_qcoord_to_vertex (g->conn, tt, quad->x + qh, quad->y + qh,
#ifdef P4_TO_P8
                          quad->z + qh,
#endif
                          hxyz);
  for (i = 0; i < 3; ++i) {
    lxyz[i] /= g->bricklength;
    hxyz[i] /= g->bricklength;
    dxyz[i] = hxyz[i] - lxyz[i];
  }
}

static void
run_pre (part_global_t * g, pi_data_t * piddata)
{
  int                 i;

  /* prepare parallel statistics */
  for (i = 0; i < PART_STATS_LAST; ++i) {
    sc_stats_init (g->si + i, snames[i]);
  }

  /* prepare quadrature rule */
  P4EST_ASSERT (PART_NQPOINTS == 3);
  qpoints[2] = sqrt (3. / 5.);
  qpoints[1] = 0.;
  qpoints[0] = -qpoints[2];
  qweights[0] = qweights[2] = 5. / 9.;
  qweights[1] = 8. / 9.;
  for (i = 0; i < 3; ++i) {
    qpoints[i] = .5 * (1. + qpoints[i]);
    qweights[i] *= .5;
  }

  /* the length of the brick is a power of 2 */
  g->bricklength = (1 << g->bricklev);

  /* initial particle density */
  piddata->sigma = pidensy_sigma;
  piddata->invs2 = 1. / SC_SQR (piddata->sigma);
  piddata->gnorm = gaussnorm (piddata->sigma);
  piddata->center[0] = pidensy_center[0];
  piddata->center[1] = pidensy_center[1];
#ifndef P4_TO_P8
  piddata->center[2] = 0.;
#else
  piddata->center[2] = pidensy_center[2];
#endif
  g->pidense = pidense;
  g->piddata = piddata;

  /* allocate arrays that are reused a lot */
  for (i = 0; i < 2; ++i) {
    g->ilh[i] = sc_array_new (sizeof (p4est_locidx_t));
    g->jlh[i] = sc_array_new (sizeof (p4est_locidx_t));
#ifdef P4_TO_P8
    g->klh[i] = sc_array_new (sizeof (p4est_locidx_t));
#else
    P4EST_ASSERT (g->klh[i] == NULL);
#endif
  }
}

static void
run_post (part_global_t * g)
{
  int                 i;

  /* clean up simulation memory */
  for (i = 0; i < 2; ++i) {
    sc_array_destroy_null (&g->ilh[i]);
    sc_array_destroy_null (&g->jlh[i]);
#ifdef P4_TO_P8
    sc_array_destroy_null (&g->klh[i]);
#else
    P4EST_ASSERT (g->klh[i] == NULL);
#endif
  }

  /* analyze parallel statistics */
  if (g->collapse) {
    for (i = 0; i < PART_STATS_LAST; ++i) {
      sc_stats_collapse (g->si + i);
    }
  }
  sc_stats_compute (g->mpicomm, PART_STATS_LAST, g->si);
  sc_stats_print (p4est_package_id, SC_LP_ESSENTIAL,
                  PART_STATS_LAST, g->si, 1, 1);
}

static double
integrate (part_global_t * g, const double lxyz[3], const double dxyz[3])
{
  int                 i, j, k;
  double              wk, wkj, wkji;
  double              qpi[3], qpj[3], qpk[3];
  double              d;

  /*** run quadrature rule ***/
  P4EST_ASSERT (PART_NQPOINTS == 3);
  for (k = 0; k < 3; ++k) {
    P4EST_ASSERT (qpoints[k] >= 0. && qpoints[k] <= 1.);
    qpi[k] = lxyz[0] + qpoints[k] * dxyz[0];
    qpj[k] = lxyz[1] + qpoints[k] * dxyz[1];
    qpk[k] = lxyz[2] + qpoints[k] * dxyz[2];
  }
  d = 0.;
#ifdef P4_TO_P8
  for (k = 0; k < 3; ++k) {
    wk = qweights[k] * dxyz[2];
#if 0
  }
#endif
#else
  k = 1;
  wk = 1.;
#endif
  for (j = 0; j < 3; ++j) {
    wkj = wk * qweights[j] * dxyz[1];
    for (i = 0; i < 3; ++i) {
      wkji = wkj * qweights[i] * dxyz[0];
      d += wkji * g->pidense (qpi[i], qpj[j], qpk[k], g->piddata);
    }
  }
#ifdef P4_TO_P8
#if 0
  {
#endif
  }
#endif
  return d;
}

static int
initrp_refine (p4est_t * p4est,
               p4est_topidx_t which_tree, p4est_quadrant_t * quadrant)
{
  qu_data_t          *qud = (qu_data_t *) quadrant->p.user_data;
  part_global_t      *g = (part_global_t *) p4est->user_pointer;
  p4est_locidx_t      ilem_particles;

  ilem_particles =
    (p4est_locidx_t) round (qud->u.d * g->num_particles / g->global_density);

  return (double) ilem_particles > g->elem_particles;
}

static void
initrp (part_global_t * g)
{
  int                 mpiret;
  int                 cycle, max_cycles;
  double              lxyz[3], hxyz[3], dxyz[3];
  double              d, ld;
  double              refine_maxd, refine_maxl;
  double              loclp[2], glolp[2];
  p4est_topidx_t      tt;
  p4est_locidx_t      lq;
  p4est_locidx_t      ilem_particles;
  p4est_gloidx_t      old_gnum, new_gnum;
  p4est_tree_t       *tree;
  p4est_quadrant_t   *quad;
  qu_data_t          *qud;

  glolp[0] = glolp[1] = 0.;
  max_cycles = g->maxlevel - g->minlevel;
  for (cycle = 0;; ++cycle) {
    /*** iterate through local cells to determine local particle density ***/
    ld = 0.;
    refine_maxd = refine_maxl = 0.;
    for (tt = g->p4est->first_local_tree; tt <= g->p4est->last_local_tree;
         ++tt) {
      tree = p4est_tree_array_index (g->p4est->trees, tt);
      for (lq = 0; lq < (p4est_locidx_t) tree->quadrants.elem_count; ++lq) {
        quad = p4est_quadrant_array_index (&tree->quadrants, lq);
        qud = (qu_data_t *) quad->p.user_data;
        loopquad (g, tt, quad, lxyz, hxyz, dxyz);

        /***  integrate density over quadrant ***/
        qud->u.d = d = integrate (g, lxyz, dxyz);
        ld += d;

        /*** maximum particle count and level ***/
        refine_maxd = SC_MAX (refine_maxd, d);
        refine_maxl = SC_MAX (refine_maxl, (double) quad->level);
      }
    }

    /*** get global integral over density ***/
    mpiret = sc_MPI_Allreduce (&ld, &g->global_density, 1, sc_MPI_DOUBLE,
                               sc_MPI_SUM, g->mpicomm);
    SC_CHECK_MPI (mpiret);
    P4EST_GLOBAL_STATISTICSF ("Global integral over density %g\n",
                              g->global_density);

    /*** get global maximum of particle count and level ***/
    loclp[0] = refine_maxd;
    loclp[1] = refine_maxl + g->bricklev;
    mpiret = sc_MPI_Allreduce (loclp, glolp, 2, sc_MPI_DOUBLE,
                               sc_MPI_MAX, g->mpicomm);
    SC_CHECK_MPI (mpiret);
    ilem_particles = (p4est_locidx_t) round
      (glolp[0] * g->num_particles / g->global_density);
    P4EST_GLOBAL_PRODUCTIONF
      ("Maximum particle number per quadrant %ld maxlevel %g\n",
       (long) ilem_particles, glolp[1]);

    /*** we have computed the density, this may be enough ***/
    if (cycle >= max_cycles || (double) ilem_particles <= g->elem_particles) {
      break;
    }

    /*** refine and balance ***/
    old_gnum = g->p4est->global_num_quadrants;
    p4est_refine_ext (g->p4est, 0, g->maxlevel - g->bricklev,
                      initrp_refine, NULL, NULL);
    new_gnum = g->p4est->global_num_quadrants;
    if (old_gnum == new_gnum) {
      /* done with refinement if no quadrants were added globally */
      /* cannot happen due to particle count above */
      SC_ABORT_NOT_REACHED ();
      break;
    }
#if 0                           /* we do not need balance for this application */
    if (cycle > 0) {
      p4est_balance (g->p4est, P4EST_CONNECT_FULL, NULL);
    }
#endif

    /*** unweighted partition ***/
    p4est_partition (g->p4est, 0, NULL);
  }

  P4EST_GLOBAL_ESSENTIALF ("Created %lld quadrants at maxlevel %g\n",
                           (long long) g->p4est->global_num_quadrants,
                           glolp[1]);
}

static void
srandquad (part_global_t * g, const double l[3])
{
  unsigned            u;

  P4EST_ASSERT (0 <= l[0] && l[0] < 1.);
  P4EST_ASSERT (0 <= l[1] && l[1] < 1.);
  P4EST_ASSERT (0 <= l[2] && l[2] < 1.);

  u = ((unsigned int) (l[2] * (1 << 10)) << 20) +
    ((unsigned int) (l[1] * (1 << 10)) << 10) +
    (unsigned int) (l[0] * (1 << 10));
  srand (u);
}

static void
create (part_global_t * g)
{
  int                 mpiret;
  int                 j;
  double              lxyz[3], hxyz[3], dxyz[3];
  double              r;
  p4est_topidx_t      tt;
  p4est_locidx_t      lpnum, lq;
  p4est_locidx_t      li, ilem_particles;
  p4est_gloidx_t      gpnum, gpoffset;
  p4est_tree_t       *tree;
  p4est_quadrant_t   *quad;
  qu_data_t          *qud;
  pa_data_t          *pad;

  P4EST_ASSERT (g->padata != NULL);

  /*** iterate through local cells and populate with particles ***/
  lpnum = 0;
  for (tt = g->p4est->first_local_tree; tt <= g->p4est->last_local_tree; ++tt) {
    tree = p4est_tree_array_index (g->p4est->trees, tt);
    for (lq = 0; lq < (p4est_locidx_t) tree->quadrants.elem_count; ++lq) {
      quad = p4est_quadrant_array_index (&tree->quadrants, lq);
      qud = (qu_data_t *) quad->p.user_data;

      /* calculate required number of particles */
      ilem_particles = (p4est_locidx_t) round
        (qud->u.d / g->global_density * g->num_particles);
      pad = (pa_data_t *) sc_array_push_count (g->padata, ilem_particles);

      /*** generate required number of particles ***/
      loopquad (g, tt, quad, lxyz, hxyz, dxyz);
      srandquad (g, lxyz);
      for (li = 0; li < ilem_particles; ++li) {
        for (j = 0; j < P4EST_DIM; ++j) {
          r = rand () / (double) RAND_MAX;
          pad->xv[j] = lxyz[j] + r * dxyz[j];

          /* begin with some velocity choice */
          pad->xv[3 + j] = 0.;
        }
#ifndef P4_TO_P8
        pad->xv[2] = pad->xv[5] = 0.;
#endif
        memset (pad->wo, 0, 6 * sizeof (double));
        memset (pad->up, 0, 6 * sizeof (double));
        pad->rm[0] = pad->rm[1] = pad->rm[2] = -1.;
        ++pad;
      }
      lpnum += ilem_particles;
      qud->u.lpend = lpnum;
      qud->premain = qud->preceive = 0;
    }
  }
  g->gplost = 0;
  gpnum = (p4est_gloidx_t) lpnum;
  mpiret = sc_MPI_Allreduce (&gpnum, &g->gpnum, 1, P4EST_MPI_GLOIDX,
                             sc_MPI_SUM, g->mpicomm);
  SC_CHECK_MPI (mpiret);
  P4EST_GLOBAL_ESSENTIALF ("Created %lld particles for %g\n",
                           (long long) g->gpnum, g->num_particles);

  /* create globally unique particle numbers */
  mpiret = sc_MPI_Exscan (&gpnum, &gpoffset, 1, P4EST_MPI_GLOIDX,
                          sc_MPI_SUM, g->mpicomm);
  SC_CHECK_MPI (mpiret);
  if (g->mpirank == 0) {
    gpoffset = 0;
  }
  pad = (pa_data_t *) sc_array_index_begin (g->padata);
  for (li = 0; li < lpnum; ++li) {
    (pad++)->id = gpoffset + li;
  }
  P4EST_ASSERT (pad == (pa_data_t *) sc_array_index_end (g->padata));
}

static void
rkrhs (part_global_t * g, const double xv[6], double rk[6])
{
  int                 i;
  int                 j;
  double              d;
  double              diff[3];

  for (i = 0; i < P4EST_DIM; ++i) {
    rk[i] = xv[3 + i];
    rk[3 + i] = 0.;
  }
#ifndef P4_TO_P8
  rk[2] = rk[5] = 0.;
#endif

  /* we use as many planets as we have dimensions */
  for (j = 0; j < P4EST_DIM; ++j) {
    d = 0.;
    /* distance is always computed in 3D space */
    for (i = 0; i < 3; ++i) {
      diff[i] = planet_xyz[j][i] - xv[i];
      d += SC_SQR (diff[i]);
    }
    d = planet_mass[j] * pow (d, -1.5);
    for (i = 0; i < P4EST_DIM; ++i) {
      rk[3 + i] += d * diff[i];
    }
  }
}

static void
rkstage (part_global_t * g, pa_data_t * pad, double h)
{
  const int           stage = g->stage;
  const int           order = g->order;
  int                 i;
  double              d;
  double              rk[6];

  /* evaluate right hand side */
  rkrhs (g, stage == 0 ? pad->xv : pad->wo, rk);

  /* compute new evaluation point if necessary */
  if (stage + 1 < order) {
    /* stage is not the last */
    d = h * prk[order - 1][0][stage];
    for (i = 0; i < 6; ++i) {
      pad->wo[i] = pad->xv[i] + d * rk[i];
    }
  }

  /* compute an update to the state */
  d = prk[order - 1][1][stage];
  if (stage == 0) {
    /* first stage */
    if (order > 1) {
      /* first stage is not the last */
      P4EST_ASSERT (stage + 1 < order);
      for (i = 0; i < 6; ++i) {
        pad->up[i] = d * rk[i];
      }
    }
    else {
      /* first stage is also the last */
      P4EST_ASSERT (stage + 1 == order);
      for (i = 0; i < 6; ++i) {
        pad->xv[i] += h * d * rk[i];
      }
    }
  }
  else {
    /* stage is not the first */
    if (stage + 1 < order) {
      /* stage is neither first nor last */
      P4EST_ASSERT (0 < stage);
      for (i = 0; i < 6; ++i) {
        pad->up[i] += d * rk[i];
      }
    }
    else {
      /* stage is last of several */
      P4EST_ASSERT (stage + 1 == order);
      for (i = 0; i < 6; ++i) {
        pad->xv[i] += h * (pad->up[i] + d * rk[i]);
      }
    }
  }
}

static int
psearch_quad (p4est_t * p4est, p4est_topidx_t which_tree,
              p4est_quadrant_t * quadrant, int pfirst, int plast,
              p4est_locidx_t local_num, void *point)
{
  part_global_t      *g = (part_global_t *) p4est->user_pointer;
#ifdef P4EST_ENABLE_DEBUG
  qu_data_t          *qud;

  if (local_num >= 0) {
    /* quadrant is a local leaf */
    qud = (qu_data_t *) quadrant->p.user_data;
    P4EST_ASSERT (qud->premain == 0);
    P4EST_ASSERT (qud->preceive == 0);
  }
#endif

  /* compute coordinate range of this quadrant */
  loopquad (g, which_tree, quadrant, g->lxyz, g->hxyz, g->dxyz);

  /* always return 1 to search particles individually */
  return 1;
}

static double      *
particle_lookfor (part_global_t * g, pa_data_t * pad)
{
  P4EST_ASSERT (0 <= g->stage && g->stage < g->order);
  P4EST_ASSERT (pad != NULL);

  return g->stage + 1 < g->order ? pad->wo : pad->xv;
}

static int
psearch_point (p4est_t * p4est, p4est_topidx_t which_tree,
               p4est_quadrant_t * quadrant, int pfirst, int plast,
               p4est_locidx_t local_num, void *point)
{
  int                 i;
  int                *pfn;
  size_t              zp;
  double             *x;
  part_global_t      *g = (part_global_t *) p4est->user_pointer;
  qu_data_t          *qud;
  pa_data_t          *pad = (pa_data_t *) point;

  /* access location of particle to be searched */
  x = particle_lookfor (g, pad);

  /* due to roundoff we call this even for a local leaf */
  for (i = 0; i < P4EST_DIM; ++i) {
    if (!(g->lxyz[i] <= x[i] && x[i] <= g->hxyz[i])) {
      /* the point is outside the search quadrant */
      return 0;
    }
  }

  /* convention for entries of pfound:
     -1              particle has not yet been found
     [0 .. mpisize)  particle found on that rank, me or other
   */

  /* find process/quadrant for this particle */
  if (local_num >= 0) {
    /* quadrant is a local leaf */
    P4EST_ASSERT (pfirst == g->mpirank && plast == g->mpirank);
    zp = sc_array_position (g->padata, point);
    pfn = (int *) sc_array_index (g->pfound, zp);
    /* first local match counts (due to roundoff there may be multiple) */
    if (*pfn != g->mpirank) {
      /* particle was either yet unfound, or found on another process */
      /* bump counter of particles in this local quadrant */

#if 0
      /* HACK */
      if (g->printn > 0 && !(pad->id % g->printn)) {
        pad = (pa_data_t *) point;
        P4EST_VERBOSEF ("Locremain particle %lld %d\n",
                        (long long) pad->id, (int) zp);
      }
#endif

      *pfn = g->mpirank;
      *(p4est_locidx_t *) sc_array_push (g->iremain) = (p4est_locidx_t) zp;
      qud = (qu_data_t *) quadrant->p.user_data;
      ++qud->premain;
    }
    /* return value will have no effect, but we must return */
    return 0;
  }
  if (pfirst == plast) {
    if (pfirst == g->mpirank) {
      /* continue recursion for local branch quadrant */
      P4EST_ASSERT (plast == g->mpirank);
      return 1;
    }
    /* found particle on a remote process */
    P4EST_ASSERT (plast != g->mpirank);
    zp = sc_array_position (g->padata, point);
    pfn = (int *) sc_array_index (g->pfound, zp);
    /* only count match if it has not been found locally or on lower rank */
    if (*pfn < 0 || (*pfn != g->mpirank && pfirst < *pfn)) {
      *pfn = pfirst;
    }
    /* return value will have no effect, but we must return */
    return 0;
  }

  /* the process for this particle has not yet been found */
  return 1;
}

static void
presearch (part_global_t * g)
{
  double              t0_searchp, t1;

  P4EST_ASSERT (g->padata != NULL);
  P4EST_ASSERT (g->pfound == NULL);
  P4EST_ASSERT (g->iremain == NULL);

  g->pfound = sc_array_new_count (sizeof (int), g->padata->elem_count);
  sc_array_memset (g->pfound, -1);

  g->iremain = sc_array_new (sizeof (p4est_locidx_t));

  /* STATS */
  t0_searchp = sc_MPI_Wtime ();

  /* search through partition for parallel branch boundaries */
  p4est_search_all (g->p4est, 0, psearch_quad, psearch_point, g->padata);

  /* STATS */
  t1 = sc_MPI_Wtime ();
  if (!g->scaling || g->verylast) {
    sc_stats_accumulate (g->si + PART_STATS_SEARCHP, t1 - t0_searchp);
  }
}

static unsigned
psend_hash (const void *v, const void *u)
{
  const comm_psend_t *ps = (const comm_psend_t *) v;

  P4EST_ASSERT (u == NULL);

  return ps->rank;
}

static int
psend_equal (const void *v1, const void *v2, const void *u)
{
  const comm_psend_t *ps1 = (const comm_psend_t *) v1;
  const comm_psend_t *ps2 = (const comm_psend_t *) v2;

  P4EST_ASSERT (u == NULL);

  return ps1->rank == ps2->rank;
}

static void
pack (part_global_t * g)
{
  int                 mpiret;
  int                 retval;
  int                *pfn;
  size_t              zz, numz;
  void              **hfound;
  p4est_locidx_t      lremain, lsend, llost;
  p4est_gloidx_t      loclrs[4], glolrs[4];
  comm_psend_t       *cps, *there;
  comm_prank_t       *trank;
#ifdef PART_SENDFULL
  pa_data_t          *pad;
#else
  double             *msg;
  double             *x;
#endif

  P4EST_ASSERT (g->psmem == NULL);
  g->psmem = sc_mempool_new (sizeof (comm_psend_t));

  P4EST_ASSERT (g->pfound != NULL);
  numz = g->pfound->elem_count;

  P4EST_ASSERT (g->padata != NULL);
  P4EST_ASSERT (g->padata->elem_count == numz);

  P4EST_ASSERT (g->psend == NULL);
  P4EST_ASSERT (g->recevs == NULL);

  g->psend = sc_hash_new (psend_hash, psend_equal, NULL, NULL);
  g->recevs = sc_array_new (sizeof (comm_prank_t));

  lremain = lsend = llost = 0;
  cps = (comm_psend_t *) sc_mempool_alloc (g->psmem);
  cps->rank = -1;
  for (zz = 0; zz < numz; ++zz) {
    pfn = (int *) sc_array_index (g->pfound, zz);

    /* ignore those that leave the domain or stay local */
    if (*pfn < 0) {
      P4EST_ASSERT (*pfn == -1);
      ++llost;
      continue;
    }
    if (*pfn == g->mpirank) {
      ++lremain;
      continue;
    }

    /* access message structure */
    P4EST_ASSERT (0 <= *pfn && *pfn < g->mpisize);
    cps->rank = *pfn;
    P4EST_ASSERT (cps->rank != g->mpirank);
    retval = sc_hash_insert_unique (g->psend, cps, &hfound);
    P4EST_ASSERT (hfound != NULL);
    there = *((comm_psend_t **) hfound);
    if (!retval) {
      /* message for this rank exists already */
      P4EST_ASSERT (there->message.elem_size == PART_MSGSIZE);
      P4EST_ASSERT (there->message.elem_count > 0);
    }
    else {
      /* message is added for this rank */
      P4EST_ASSERT (there == cps);
      trank = (comm_prank_t *) sc_array_push (g->recevs);
      trank->rank = there->rank;
      trank->psend = there;
      sc_array_init (&there->message, PART_MSGSIZE);
      cps = (comm_psend_t *) sc_mempool_alloc (g->psmem);
      cps->rank = -1;
    }

    /* add to message buffer */
#ifdef PART_SENDFULL
    pad = (pa_data_t *) sc_array_push (&there->message);
    memcpy (pad, sc_array_index (g->padata, zz), sizeof (pa_data_t));
#else
    msg = (double *) sc_array_push (&there->message);
    x = particle_lookfor (g, (pa_data_t *) sc_array_index (g->padata, zz));
    memcpy (msg, x, 3 * sizeof (double));
#endif

    /* this particle is to be sent to another process */
    ++lsend;
  }
  sc_mempool_free (g->psmem, cps);
  sc_array_sort (g->recevs, comm_prank_compare);
  /* TODO: hash table still needed?  Yes, unless we rearrange data sent */

  /* keep track of some global numbers */
  loclrs[0] = (p4est_gloidx_t) lremain;
  loclrs[1] = (p4est_gloidx_t) lsend;
  loclrs[2] = (p4est_gloidx_t) llost;
  loclrs[3] = (p4est_gloidx_t) g->recevs->elem_count;
  mpiret = sc_MPI_Allreduce (loclrs, glolrs, 4, P4EST_MPI_GLOIDX,
                             sc_MPI_SUM, g->mpicomm);
  SC_CHECK_MPI (mpiret);
  P4EST_GLOBAL_STATISTICSF
    ("Stage %d from %lld remain %lld sent %lld lost %lld avg peers %.3g\n",
     g->stage, (long long) g->gpnum, (long long) glolrs[0],
     (long long) glolrs[1], (long long) glolrs[2],
     glolrs[3] / (double) g->mpisize);
  P4EST_ASSERT (glolrs[0] + glolrs[1] + glolrs[2] == g->gpnum);
  g->gplost += glolrs[2];
  g->gpnum -= glolrs[2];

  /* count number of peers as statistics */
  if (!g->scaling || g->verylast) {
    sc_stats_accumulate (g->si + PART_STATS_PEERS, (double) loclrs[3]);
  }

  /* another array that is no longer needed */
  sc_array_destroy_null (&g->pfound);
}

static void
comm (part_global_t * g)
{
  int                 mpiret;
  int                 i;
  int                 num_receivers;
  int                 num_senders;
  int                 count, cucount;
  int                 msglen;
  double              t0_notify, t0_wait1, t0_comm, t1;
  sc_MPI_Request     *reqs;
  sc_array_t         *notif, *payl;
  sc_array_t         *arr;
  comm_psend_t       *cps;
  comm_prank_t       *trank;

  P4EST_ASSERT (g->psmem != NULL);
  P4EST_ASSERT (g->recevs != NULL);

  P4EST_ASSERT (g->prebuf == NULL);
  P4EST_ASSERT (g->recv_req == NULL);
  P4EST_ASSERT (g->send_req == NULL);

  /* STATS */
  t0_comm = sc_MPI_Wtime ();

  /* pass receiver ranks and message size to notify */
  num_receivers = (int) g->recevs->elem_count;
  P4EST_ASSERT (0 <= num_receivers && num_receivers < g->mpisize);
  notif = sc_array_new_count (sizeof (int), num_receivers);
  payl = sc_array_new_count (sizeof (int), num_receivers);
  if (g->olap_notify) {
    P4EST_ASSERT (g->send_req == NULL);
    g->send_req = sc_array_new_count (sizeof (sc_MPI_Request), num_receivers);
  }
  for (i = 0; i < num_receivers; ++i) {
    trank = (comm_prank_t *) sc_array_index_int (g->recevs, i);
    *(int *) sc_array_index_int (notif, i) = trank->rank;
    cps = trank->psend;
    P4EST_ASSERT (trank->rank == cps->rank);
    arr = &cps->message;
    P4EST_ASSERT (arr->elem_size == PART_MSGSIZE);
    P4EST_ASSERT (arr->elem_count > 0);
    *(int *) sc_array_index_int (payl, i) = (int) arr->elem_count;
    if (g->olap_notify) {
      msglen = (int) (arr->elem_count * arr->elem_size);
      mpiret = sc_MPI_Isend
        (arr->array, msglen, sc_MPI_BYTE, cps->rank, COMM_TAG_PART,
         g->mpicomm, (sc_MPI_Request *) sc_array_index_int (g->send_req, i));
      SC_CHECK_MPI (mpiret);
    }
  }

  /* STATS */
  t0_notify = sc_MPI_Wtime ();

  /* reverse communication pattern */
  sc_notify_ext (notif, NULL, payl, NULL, g->mpicomm);
  P4EST_ASSERT (payl->elem_count == notif->elem_count);
  num_senders = (int) notif->elem_count;
  P4EST_ASSERT (0 <= num_senders && num_senders < g->mpisize);

  /* STATS */
  t1 = sc_MPI_Wtime ();
  if (!g->scaling || g->verylast) {
    sc_stats_accumulate (g->si + PART_STATS_NOTIFY, t1 - t0_notify);
  }

  /* receive particles into a flat array over all processes */
  cucount = 0;
  for (i = 0; i < num_senders; ++i) {
    cucount += *(int *) sc_array_index_int (payl, i);
  }
  g->prebuf = sc_array_new_count (PART_MSGSIZE, cucount);

  /* post non-blocking receive */
  g->recv_req = sc_array_new_count (sizeof (sc_MPI_Request), num_senders);
  cucount = 0;
  for (i = 0; i < num_senders; ++i) {
    count = *(int *) sc_array_index_int (payl, i);
    msglen = count * (int) PART_MSGSIZE;
    mpiret = sc_MPI_Irecv
      (sc_array_index (g->prebuf, cucount), msglen, sc_MPI_BYTE,
       *(int *) sc_array_index_int (notif, i), COMM_TAG_PART, g->mpicomm,
       (sc_MPI_Request *) sc_array_index_int (g->recv_req, i));
    SC_CHECK_MPI (mpiret);
    cucount += count;
  }
  P4EST_ASSERT (cucount == (int) g->prebuf->elem_count);
  sc_array_destroy_null (&notif);
  sc_array_destroy_null (&payl);

  /* post non-blocking send if not done earlier */
  if (!g->olap_notify) {
    P4EST_ASSERT (g->send_req == NULL);
    g->send_req = sc_array_new_count (sizeof (sc_MPI_Request), num_receivers);
    for (i = 0; i < num_receivers; ++i) {
      trank = (comm_prank_t *) sc_array_index_int (g->recevs, i);
      cps = trank->psend;
      P4EST_ASSERT (trank->rank == cps->rank);
      arr = &cps->message;
      P4EST_ASSERT (arr->elem_size == PART_MSGSIZE);
      P4EST_ASSERT (arr->elem_count > 0);
      msglen = (int) (arr->elem_count * arr->elem_size);
      mpiret = sc_MPI_Isend
        (arr->array, msglen, sc_MPI_BYTE, cps->rank, COMM_TAG_PART,
         g->mpicomm, (sc_MPI_Request *) sc_array_index_int (g->send_req, i));
      SC_CHECK_MPI (mpiret);
    }
  }

  /* STATS */
  t0_wait1 = sc_MPI_Wtime ();

  /* wait for all incoming messages to complete */
  reqs = (sc_MPI_Request *) sc_array_index_begin (g->recv_req);
  mpiret = sc_MPI_Waitall (num_senders, reqs, sc_MPI_STATUSES_IGNORE);
  SC_CHECK_MPI (mpiret);
  sc_array_destroy_null (&g->recv_req);

  /* STATS */
  t1 = sc_MPI_Wtime ();
  if (!g->scaling || g->verylast) {
    sc_stats_accumulate (g->si + PART_STATS_WAIT1, t1 - t0_wait1);
    sc_stats_accumulate (g->si + PART_STATS_COMM, t1 - t0_comm);
  }
}

static int
slocal_quad (p4est_t * p4est, p4est_topidx_t which_tree,
             p4est_quadrant_t * quadrant, p4est_locidx_t local_num,
             void *point)
{
  part_global_t      *g = (part_global_t *) p4est->user_pointer;
#ifdef P4EST_ENABLE_DEBUG
  qu_data_t          *qud;

  if (local_num >= 0) {
    qud = (qu_data_t *) quadrant->p.user_data;
    P4EST_ASSERT (qud->preceive == 0);
  }
#endif

  /* compute coordinate range of this quadrant */
  loopquad (g, which_tree, quadrant, g->lxyz, g->hxyz, g->dxyz);

  /* always return 1 to search particles individually */
  return 1;
}

static int
slocal_point (p4est_t * p4est, p4est_topidx_t which_tree,
              p4est_quadrant_t * quadrant, p4est_locidx_t local_num,
              void *point)
{
  int                 i;
  char               *cf;
  size_t              zp;
  part_global_t      *g = (part_global_t *) p4est->user_pointer;
  qu_data_t          *qud;
#ifdef PART_SENDFULL
  double             *x;
  pa_data_t          *pad = (pa_data_t *) point;

  /* access location of particle to be searched */
  x = particle_lookfor (g, pad);
#else
  double             *x = (double *) point;
#endif

  /* due to roundoff we call this even for a local leaf */
  for (i = 0; i < P4EST_DIM; ++i) {
    if (!(g->lxyz[i] <= x[i] && x[i] <= g->hxyz[i])) {
      /* the point is outside the search quadrant */
      return 0;
    }
  }

  if (local_num >= 0) {
    /* quadrant is a local leaf */
    /* first local match counts (due to roundoff there may be multiple) */
    zp = sc_array_position (g->prebuf, point);
    cf = (char *) sc_array_index (g->cfound, zp);
    if (!*cf) {
      /* make sure this particle is not found twice */
      *cf = 1;

      /* count this particle in its target quadrant */
      *(p4est_locidx_t *) sc_array_push (g->ireceive) = (p4est_locidx_t) zp;
      qud = (qu_data_t *) quadrant->p.user_data;
      ++qud->preceive;
    }

    /* return value will have no effect */
    return 0;
  }

  /* the leaf for this particle has not yet been found */
  return 1;
}

static void
postsearch (part_global_t * g)
{
  double              t0_searchl, t1;

  P4EST_ASSERT (g->ireceive == NULL);
  P4EST_ASSERT (g->cfound == NULL);
  P4EST_ASSERT (g->prebuf != NULL);

  g->ireceive = sc_array_new (sizeof (p4est_locidx_t));
  g->cfound = sc_array_new_count (sizeof (char), g->prebuf->elem_count);
  sc_array_memset (g->cfound, 0);

  /* STATS */
  t0_searchl = sc_MPI_Wtime ();

  /* run local search to find particles sent to us */
  p4est_search_local (g->p4est, 0, slocal_quad, slocal_point, g->prebuf);
  sc_array_destroy_null (&g->cfound);

  /* STATS */
  t1 = sc_MPI_Wtime ();
  if (!g->scaling || g->verylast) {
    sc_stats_accumulate (g->si + PART_STATS_SEARCHL, t1 - t0_searchl);
  }
}

static int
adapt_coarsen (p4est_t * p4est, p4est_topidx_t which_tree,
               p4est_quadrant_t * quadrants[])
{
  int                 i;
  p4est_locidx_t      remain, receive;
  qu_data_t          *qud;
  part_global_t      *g = (part_global_t *) p4est->user_pointer;

  /* TODO: coarsen/refine on sum of still-there and to-be-there particles? */

  /* maybe this quadrant is just called for counting, or too big already */
  if (quadrants[1] == NULL ||
      quadrants[0]->level == g->minlevel - g->bricklev) {
    qud = (qu_data_t *) quadrants[0]->p.user_data;
#ifdef P4EST_ENABLE_DEBUG
    qud->u.lpend = -1;
#endif
    g->ireindex += qud->premain;
    g->irvindex += qud->preceive;
    return 0;
  }

  /* sum expected particle count over siblings */
  remain = receive = 0;
  for (i = 0; i < P4EST_CHILDREN; ++i) {
    qud = (qu_data_t *) quadrants[i]->p.user_data;
#ifdef P4EST_ENABLE_DEBUG
    qud->u.lpend = -1;
#endif
    remain += qud->premain;
    receive += qud->preceive;
  }
  if ((double) (remain + receive) < .5 * g->elem_particles) {
    /* we will coarsen and adjust ireindex, irvindex in adapt_replace */
    g->qremain = remain;
    g->qreceive = receive;
    return 1;
  }
  else {
    /* we will not coarsen and proceed with next quadrant */
    qud = (qu_data_t *) quadrants[0]->p.user_data;
    g->ireindex += qud->premain;
    g->irvindex += qud->preceive;
    return 0;
  }
}

static int
adapt_refine (p4est_t * p4est, p4est_topidx_t which_tree,
              p4est_quadrant_t * quadrant)
{
  qu_data_t          *qud = (qu_data_t *) quadrant->p.user_data;
  part_global_t      *g = (part_global_t *) p4est->user_pointer;

  /* we have set this to -1 in adapt_coarsen */
  P4EST_ASSERT (qud->u.lpend == -1);

  if ((double) (qud->premain + qud->preceive) > g->elem_particles) {
    /* we are trying to refine, we will possibly go into the replace function */
    g->ire2 = g->ireindex;
    g->ireindex += qud->premain;
    g->irv2 = g->irvindex;
    g->irvindex += qud->preceive;
    return 1;
  }
  else {
    /* maintain cumulative particle count for next quadrant */
    g->ireindex += qud->premain;
    g->irvindex += qud->preceive;
    return 0;
  }
}

typedef enum pa_mode
{
  PA_MODE_REMAIN,
  PA_MODE_RECEIVE,
  PA_MODE_LOCATE
}
pa_mode_t;

static void
split_by_coord (part_global_t * g, sc_array_t * in,
                sc_array_t * out[2], pa_mode_t mode, int component,
                const double lxyz[3], const double dxyz[3])
{
  p4est_locidx_t      ppos;
  const double       *x;
  size_t              zz, znum;
  pa_data_t          *pad;

  P4EST_ASSERT (in != NULL);
  P4EST_ASSERT (in->elem_size == sizeof (p4est_locidx_t));
  P4EST_ASSERT (out != NULL);
  P4EST_ASSERT (out[0] != NULL);
  P4EST_ASSERT (out[0]->elem_size == sizeof (p4est_locidx_t));
  sc_array_truncate (out[0]);
  P4EST_ASSERT (out[1] != NULL);
  P4EST_ASSERT (out[1]->elem_size == sizeof (p4est_locidx_t));
  sc_array_truncate (out[1]);

  znum = in->elem_count;
  for (zz = 0; zz < znum; ++zz) {
    ppos = *(p4est_locidx_t *) sc_array_index (in, zz);
    if (mode == PA_MODE_REMAIN) {
      pad = (pa_data_t *) sc_array_index (g->padata, ppos);
      x = particle_lookfor (g, pad);
    }
    else if (mode == PA_MODE_RECEIVE) {
#ifdef PART_SENDFULL
      pad = (pa_data_t *) sc_array_index (g->prebuf, ppos);
      x = particle_lookfor (g, pad);
#else
      x = (const double *) sc_array_index (g->prebuf, ppos);
#endif
    }
    else {
      P4EST_ASSERT (mode == PA_MODE_LOCATE);
      pad = (pa_data_t *) sc_array_index (g->padata, ppos);
      x = pad->xv;
    }
    if (x[component] <= lxyz[component] + .5 * dxyz[component]) {
      *(p4est_locidx_t *) sc_array_push (out[0]) = ppos;
    }
    else {
      *(p4est_locidx_t *) sc_array_push (out[1]) = ppos;
    }
  }
}

static void
adapt_replace (p4est_t * p4est, p4est_topidx_t which_tree,
               int num_outgoing, p4est_quadrant_t * outgoing[],
               int num_incoming, p4est_quadrant_t * incoming[])
{
#ifdef P4EST_ENABLE_DEBUG
  int                 i;
  p4est_locidx_t      remain, receive;
  qu_data_t          *qod;
#endif
  int                 wx, wy, wz;
  double              lxyz[3], hxyz[3], dxyz[3];
  sc_array_t          iview, *arr;
  p4est_locidx_t      irem, ibeg;
  p4est_quadrant_t  **pchild;
  qu_data_t          *qud;
  part_global_t      *g = (part_global_t *) p4est->user_pointer;

  if (num_outgoing == P4EST_CHILDREN) {
    P4EST_ASSERT (num_incoming == 1);
    /* we are coarsening */
    qud = (qu_data_t *) incoming[0]->p.user_data;
#ifdef P4EST_ENABLE_DEBUG
    qud->u.lpend = -1;

    /* sum counts over siblings */
    remain = receive = 0;
    for (i = 0; i < P4EST_CHILDREN; ++i) {
      qod = (qu_data_t *) outgoing[i]->p.user_data;
      P4EST_ASSERT (qod->u.lpend == -1);
      remain += qod->premain;
      receive += qod->preceive;
    }
    P4EST_ASSERT (remain == g->qremain);
    P4EST_ASSERT (receive == g->qreceive);
#endif
    g->ireindex += (qud->premain = g->qremain);
    g->irvindex += (qud->preceive = g->qreceive);
  }
  else {
    P4EST_ASSERT (num_outgoing == 1);
    P4EST_ASSERT (num_incoming == P4EST_CHILDREN);
    /* we are refining */

    /* access parent quadrant */
    loopquad (g, which_tree, outgoing[0], lxyz, hxyz, dxyz);
#ifdef P4EST_ENABLE_DEBUG
    qod = (qu_data_t *) outgoing[0]->p.user_data;
    P4EST_ASSERT (qod->u.lpend == -1);
#endif

    /* recover window onto remaining particles for the new family */
    ibeg = g->ire2;
    irem = g->ireindex - ibeg;
    P4EST_ASSERT (irem >= 0);
    sc_array_init_view (&iview, g->iremain, ibeg, irem);
    P4EST_ASSERT (qod->premain == irem);

    /* sort remaining particles into the children */
    pchild = incoming;
#ifdef P4_TO_P8
    split_by_coord (g, &iview, g->klh, PA_MODE_REMAIN, 2, lxyz, dxyz);
    for (wz = 0; wz < 2; ++wz) {
#if 0
    }
#endif
#else
    P4EST_ASSERT (g->klh[0] == NULL);
    P4EST_ASSERT (g->klh[1] == NULL);
    g->klh[0] = &iview;
    wz = 0;
#endif
    split_by_coord (g, g->klh[wz], g->jlh, PA_MODE_REMAIN, 1, lxyz, dxyz);
    for (wy = 0; wy < 2; ++wy) {
      split_by_coord (g, g->jlh[wy], g->ilh, PA_MODE_REMAIN, 0, lxyz, dxyz);
      for (wx = 0; wx < 2; ++wx) {
        /* we have a set of particles for child 4 * wz + 2 * wy + wx */
        arr = g->ilh[wx];
        sc_array_init_view (&iview, g->iremain, ibeg, arr->elem_count);
        sc_array_paste (&iview, arr);
        qud = (qu_data_t *) (*pchild++)->p.user_data;
#ifdef P4EST_ENABLE_DEBUG
        qud->u.lpend = -1;
#endif
        ibeg += (qud->premain = (p4est_locidx_t) arr->elem_count);
      }
    }
#ifdef P4_TO_P8
#if 0
    {
#endif
    }
#endif
    P4EST_ASSERT (ibeg == g->ireindex);
    P4EST_ASSERT (pchild == incoming + P4EST_CHILDREN);

    /* recover window onto received particles for the new family */
    ibeg = g->irv2;
    irem = g->irvindex - ibeg;
    P4EST_ASSERT (irem >= 0);
    sc_array_init_view (&iview, g->ireceive, ibeg, irem);
    P4EST_ASSERT (qod->preceive == irem);

    /* sort received particles into the children */
    pchild = incoming;
#ifdef P4_TO_P8
    split_by_coord (g, &iview, g->klh, PA_MODE_RECEIVE, 2, lxyz, dxyz);
    for (wz = 0; wz < 2; ++wz) {
#if 0
    }
#endif
#else
    P4EST_ASSERT (g->klh[0] == &iview);
    P4EST_ASSERT (g->klh[1] == NULL);
    wz = 0;
#endif
    split_by_coord (g, g->klh[wz], g->jlh, PA_MODE_RECEIVE, 1, lxyz, dxyz);
    for (wy = 0; wy < 2; ++wy) {
      split_by_coord (g, g->jlh[wy], g->ilh, PA_MODE_RECEIVE, 0, lxyz, dxyz);
      for (wx = 0; wx < 2; ++wx) {
        /* we have a set of particles for child 4 * wz + 2 * wy + wx */
        arr = g->ilh[wx];
        sc_array_init_view (&iview, g->ireceive, ibeg, arr->elem_count);
        sc_array_paste (&iview, arr);
        qud = (qu_data_t *) (*pchild++)->p.user_data;
        P4EST_ASSERT (qud->u.lpend == -1);
        ibeg += (qud->preceive = (p4est_locidx_t) arr->elem_count);
      }
    }
#ifdef P4_TO_P8
#if 0
    {
#endif
    }
#endif
    P4EST_ASSERT (ibeg == g->irvindex);
    P4EST_ASSERT (pchild == incoming + P4EST_CHILDREN);

#ifndef P4_TO_P8
    g->klh[0] = NULL;
    P4EST_ASSERT (g->klh[1] == NULL);
#endif
  }
}

static void
adapt (part_global_t * g)
{
  P4EST_ASSERT (g->padata != NULL);
  P4EST_ASSERT (g->prebuf != NULL);
  P4EST_ASSERT (g->iremain != NULL);
  P4EST_ASSERT (g->ireceive != NULL);

  /* coarsen the forest according to expected number of particles */
  g->ireindex = g->irvindex = 0;
  p4est_coarsen_ext (g->p4est, 0, 1, adapt_coarsen, NULL, adapt_replace);
  P4EST_ASSERT ((size_t) g->ireindex == g->iremain->elem_count);
  P4EST_ASSERT ((size_t) g->irvindex == g->ireceive->elem_count);

  /* refine the forest according to expected number of particles */
  g->ireindex = g->ire2 = 0;
  g->irvindex = g->irv2 = 0;
  p4est_refine_ext (g->p4est, 0, g->maxlevel - g->bricklev,
                    adapt_refine, NULL, adapt_replace);
  P4EST_ASSERT ((size_t) g->ireindex == g->iremain->elem_count);
  P4EST_ASSERT ((size_t) g->irvindex == g->ireceive->elem_count);

  /* TODO: coarsen and refine repeatedly if necessary */
}

static void
regroup (part_global_t * g)
{
  sc_array_t         *newpa;
  p4est_topidx_t      tt;
  p4est_locidx_t      newnum;
  p4est_locidx_t      ppos;
  p4est_locidx_t      lq, prev;
  p4est_locidx_t      qboth, li;
  p4est_locidx_t     *premain, *preceive;
  p4est_tree_t       *tree;
  p4est_quadrant_t   *quad;
  qu_data_t          *qud;
  pa_data_t          *pad;

  P4EST_ASSERT (g->padata != NULL);
  P4EST_ASSERT (g->prebuf != NULL);
  P4EST_ASSERT (g->iremain != NULL);
  P4EST_ASSERT (g->ireceive != NULL);

#ifndef PART_SENDFULL
#error "The following code is no longer compatible with SENDFULL"
#endif

  newnum =
    (p4est_locidx_t) (g->iremain->elem_count + g->ireceive->elem_count);
  P4EST_VERBOSEF ("New local particle number %lld\n", (long long) newnum);

  /* regroup remaining and received particles after adaptation */
  premain = (p4est_locidx_t *) sc_array_index_begin (g->iremain);
  preceive = (p4est_locidx_t *) sc_array_index_begin (g->ireceive);
  newpa = sc_array_new_count (sizeof (pa_data_t), newnum);
  pad = (pa_data_t *) sc_array_index_begin (newpa);
  prev = 0;
  for (tt = g->p4est->first_local_tree; tt <= g->p4est->last_local_tree; ++tt) {
    tree = p4est_tree_array_index (g->p4est->trees, tt);
    for (lq = 0; lq < (p4est_locidx_t) tree->quadrants.elem_count; ++lq) {
      quad = p4est_quadrant_array_index (&tree->quadrants, lq);
      qud = (qu_data_t *) quad->p.user_data;
      qboth = qud->premain + qud->preceive;
      if (qboth == 0) {
        qud->u.lpend = prev;
        qud->premain = qud->preceive = 0;
        continue;
      }
      prev += qboth;
      P4EST_ASSERT (prev <= newnum);
      for (li = 0; li < qud->premain; ++li) {
        P4EST_ASSERT (premain != NULL);
        P4EST_ASSERT
          (premain < (p4est_locidx_t *) sc_array_index_end (g->iremain));
        ppos = *premain++;
        memcpy (pad++, sc_array_index (g->padata, ppos), sizeof (pa_data_t));
#ifdef P4EST_ENABLE_DEBUG
        --qboth;
#endif
      }
      for (li = 0; li < qud->preceive; ++li) {
        P4EST_ASSERT (preceive != NULL);
        P4EST_ASSERT
          (preceive < (p4est_locidx_t *) sc_array_index_end (g->ireceive));
        ppos = *preceive++;
        memcpy (pad++, sc_array_index (g->prebuf, ppos), sizeof (pa_data_t));
#ifdef P4EST_ENABLE_DEBUG
        --qboth;
#endif
      }
      P4EST_ASSERT (qboth == 0);
      qud->u.lpend = prev;
      qud->premain = qud->preceive = 0;
    }
  }
  P4EST_ASSERT (prev == newnum);
  P4EST_ASSERT (pad == (pa_data_t *) sc_array_index_end (newpa));
  P4EST_ASSERT
    (premain == (p4est_locidx_t *) sc_array_index_end (g->iremain));
  sc_array_destroy_null (&g->iremain);
  P4EST_ASSERT
    (preceive == (p4est_locidx_t *) sc_array_index_end (g->ireceive));
  sc_array_destroy_null (&g->ireceive);
  sc_array_destroy_null (&g->prebuf);
  sc_array_destroy (g->padata);
  g->padata = newpa;
}

static void
pprint (part_global_t * g, double t)
{
  int                 k;
  double              d, ds;
  p4est_locidx_t      li, lpnum;
  pa_data_t          *pad;

  P4EST_ASSERT (g->padata != NULL);
  P4EST_ASSERT (g->stage == 0);

  /* only output when specified */
  if (g->printn <= 0) {
    return;
  }

  lpnum = (p4est_locidx_t) g->padata->elem_count;
  pad = (pa_data_t *) sc_array_index_begin (g->padata);
  for (li = 0; li < lpnum; ++li) {
    if (!(pad->id % g->printn)) {
      /* has the particle advanced far enough? */
      ds = 0.;
      for (k = 0; k < P4EST_DIM; ++k) {
        d = pad->xv[k] - pad->rm[k];
        ds += d * d;
      }
      if (ds >= SC_SQR (.005)) {
        memcpy (pad->rm, pad->xv, P4EST_DIM * sizeof (double));

        /* print current particle location */
        P4EST_ESSENTIALF ("T %g I %lld X %g %g %g V %g %g %g\n",
                          t, (long long) pad->id,
                          pad->xv[0], pad->xv[1], pad->xv[2],
                          pad->xv[3], pad->xv[4], pad->xv[5]);
      }
    }
    ++pad;
  }
  P4EST_ASSERT (pad == (pa_data_t *) sc_array_index_end (g->padata));
}

static void
waitmpi (part_global_t * g)
{
  int                 mpiret;
  int                 i;
  int                 num_receivers;
  double              t0_wait2, t1;
  sc_MPI_Request     *reqs;
  comm_psend_t       *cps;
  comm_prank_t       *trank;

  P4EST_ASSERT (g->recv_req == NULL);
  P4EST_ASSERT (g->send_req != NULL);
  P4EST_ASSERT (g->recevs != NULL);
  P4EST_ASSERT (g->psend != NULL);
  P4EST_ASSERT (g->psmem != NULL);

  /* STATS */
  t0_wait2 = sc_MPI_Wtime ();

  /* wait for sent messages to complete */
  num_receivers = (int) g->recevs->elem_count;
  reqs = (sc_MPI_Request *) sc_array_index_begin (g->send_req),
    mpiret = sc_MPI_Waitall (num_receivers, reqs, sc_MPI_STATUSES_IGNORE);
  SC_CHECK_MPI (mpiret);
  sc_array_destroy_null (&g->send_req);

  /* STATS */
  t1 = sc_MPI_Wtime ();
  if (!g->scaling || g->verylast) {
    sc_stats_accumulate (g->si + PART_STATS_WAIT2, t1 - t0_wait2);
  }

  /* free send buffer */
  for (i = 0; i < num_receivers; ++i) {
    trank = (comm_prank_t *) sc_array_index_int (g->recevs, i);
    cps = trank->psend;
    P4EST_ASSERT (cps->rank == trank->rank);
    P4EST_ASSERT (cps->message.elem_size == PART_MSGSIZE);
    P4EST_ASSERT (cps->message.elem_count > 0);
    sc_array_reset (&cps->message);
  }
  sc_array_destroy_null (&g->recevs);
  sc_hash_destroy (g->psend);
  g->psend = NULL;
  sc_mempool_destroy (g->psmem);
  g->psmem = NULL;
}

static int
part_weight (p4est_t * p4est,
             p4est_topidx_t which_tree, p4est_quadrant_t * quadrant)
{
  p4est_locidx_t      ilem_particles;
  part_global_t      *g = (part_global_t *) p4est->user_pointer;
  qu_data_t          *qud = (qu_data_t *) quadrant->p.user_data;

  ilem_particles = qud->u.lpend - g->prevlp;
  g->prevlp = qud->u.lpend;

  *(int *) sc_array_index (g->src_fixed, g->qcount++) =
    (int) (ilem_particles * sizeof (pa_data_t));

  return 1 + ilem_particles;
}

static void
part (part_global_t * g)
{
  double              t0_trans1, t1_trans1;
  double              t0_trans2, t1_trans2;
  sc_array_t         *dest_data;
  p4est_topidx_t      tt;
  p4est_locidx_t      ldatasiz, lcount;
  p4est_locidx_t      dest_quads, src_quads;
  p4est_locidx_t      dest_parts;
  p4est_locidx_t      lquad, lq;
  p4est_locidx_t      lpnum;
  p4est_gloidx_t      gshipped;
  p4est_gloidx_t     *src_gfq;
  p4est_tree_t       *tree;
  p4est_quadrant_t   *quad;
  qu_data_t          *qud;

  P4EST_ASSERT (g->src_fixed == NULL);
  P4EST_ASSERT (g->dest_fixed == NULL);

  if (g->mpisize == 1) {
    return;
  }

  /* remember current forest and its particle counts per quadrant */
  src_gfq = P4EST_ALLOC (p4est_gloidx_t, g->mpisize + 1);
  memcpy (src_gfq, g->p4est->global_first_quadrant,
          (g->mpisize + 1) * sizeof (p4est_gloidx_t));
  src_quads = g->p4est->local_num_quadrants;
  P4EST_ASSERT ((p4est_gloidx_t) src_quads ==
                src_gfq[g->mpirank + 1] - src_gfq[g->mpirank]);
  g->src_fixed = sc_array_new_count (sizeof (int), src_quads);

  /* count particles per quadrant in callback */
  g->qcount = 0;
  g->prevlp = 0;
  gshipped = p4est_partition_ext (g->p4est, 1, part_weight);
  P4EST_ASSERT (g->qcount == src_quads);
  dest_quads = g->p4est->local_num_quadrants;
  P4EST_ASSERT (g->prevlp == (int) g->padata->elem_count);

  /* if nothing happens, we're done */
  if (gshipped == 0) {
    sc_array_destroy_null (&g->src_fixed);
    P4EST_FREE (src_gfq);
    return;
  }

  /* STATS */
  t0_trans1 = sc_MPI_Wtime ();

  /* transfer particle counts per quadrant to new partition */
  g->dest_fixed = sc_array_new_count (sizeof (int), dest_quads);
  p4est_transfer_fixed (g->p4est->global_first_quadrant, src_gfq,
                        g->mpicomm, COMM_TAG_FIXED,
                        (int *) g->dest_fixed->array,
                        (const int *) g->src_fixed->array, sizeof (int));

  /* transfer particle data to new partition */
  ldatasiz = (p4est_locidx_t) sizeof (pa_data_t);
  dest_parts = 0;
  for (lq = 0; lq < dest_quads; ++lq) {
    dest_parts += *(int *) sc_array_index (g->dest_fixed, lq);
  }
  P4EST_ASSERT (dest_parts % ldatasiz == 0);

  /* STATS */
  t1_trans1 = t0_trans2 = sc_MPI_Wtime ();

  dest_parts /= ldatasiz;
  dest_data = sc_array_new_count (sizeof (pa_data_t), dest_parts);
  p4est_transfer_custom (g->p4est->global_first_quadrant, src_gfq,
                         g->mpicomm, COMM_TAG_CUSTOM,
                         (pa_data_t *) dest_data->array,
                         (const int *) g->dest_fixed->array,
                         (const pa_data_t *) g->padata->array,
                         (const int *) g->src_fixed->array);

  /* clean up and keep new particle data */
  sc_array_destroy_null (&g->src_fixed);
  P4EST_FREE (src_gfq);
  sc_array_destroy (g->padata);
  g->padata = dest_data;

  /* STATS */
  t1_trans2 = sc_MPI_Wtime ();
  if (!g->scaling || g->verylast) {
    sc_stats_accumulate (g->si + PART_STATS_TRANSF, t1_trans1 - t0_trans1);
    sc_stats_accumulate (g->si + PART_STATS_TRANSC, t1_trans2 - t0_trans2);
  }

  /* reassign cumulative particle counts */
  lpnum = 0;
  lquad = 0;
  for (tt = g->p4est->first_local_tree; tt <= g->p4est->last_local_tree; ++tt) {
    tree = p4est_tree_array_index (g->p4est->trees, tt);
    for (lq = 0; lq < (p4est_locidx_t) tree->quadrants.elem_count; ++lq) {
      /* access quadrant */
      quad = p4est_quadrant_array_index (&tree->quadrants, lq);
      qud = (qu_data_t *) quad->p.user_data;
      P4EST_ASSERT (qud->premain == 0);
      P4EST_ASSERT (qud->preceive == 0);

      /* back out particle count in quadrant from data size */
      lcount = *(int *) sc_array_index (g->dest_fixed, lquad);
      P4EST_ASSERT (lcount % ldatasiz == 0);
      lcount /= ldatasiz;
      lpnum += lcount;
      qud->u.lpend = lpnum;
      ++lquad;
    }
  }
  P4EST_ASSERT (lquad == dest_quads);
  P4EST_ASSERT (lpnum == dest_parts);
  sc_array_destroy_null (&g->dest_fixed);
}

static void
outp (part_global_t * g, int k)
{
  char                filename[BUFSIZ];
  sc_array_t         *pdata;
  p4est_topidx_t      tt;
  p4est_locidx_t      lpnum, lq;
  p4est_locidx_t      lall, ilem_particles;
  p4est_tree_t       *tree;
  p4est_quadrant_t   *quad;
  qu_data_t          *qud;
  p4est_vtk_context_t *cont;

  /* only output when specified */
  if (g->vtk <= 0 || k % g->vtk) {
    return;
  }

  /* run-once loop for clean return */
  cont = NULL;
  pdata = NULL;
  do {
    /* open files for output */
    snprintf (filename, BUFSIZ, "%s_%06d", g->prefix, k);
    if (!g->scaling) {
      cont = p4est_vtk_context_new (g->p4est, filename);
      if (NULL == p4est_vtk_write_header (cont)) {
        P4EST_LERRORF ("Failed to write header for %s\n", filename);
        break;
      }
    }

    /* prepare cell data for output */
    pdata = sc_array_new_count
      (sizeof (double), g->p4est->local_num_quadrants);
    for (lpnum = 0, lall = 0, tt = g->p4est->first_local_tree;
         tt <= g->p4est->last_local_tree; ++tt) {
      tree = p4est_tree_array_index (g->p4est->trees, tt);
      for (lq = 0; lq < (p4est_locidx_t) tree->quadrants.elem_count; ++lq) {

        /* fetch number of particles in quadrant */
        quad = p4est_quadrant_array_index (&tree->quadrants, lq);
        qud = (qu_data_t *) quad->p.user_data;
        ilem_particles = qud->u.lpend - lpnum;
        *(double *) sc_array_index (pdata, lall++) = (double) ilem_particles;

        /* move to next quadrant */
        lpnum = qud->u.lpend;
      }
    }

    /* write cell data to file */
    if (!g->scaling) {
      if (NULL == p4est_vtk_write_cell_dataf
          (cont, 1, 1, 1, g->mpiwrap, 1, 0, "particles", pdata, cont)) {
        P4EST_LERRORF ("Failed to write cell data for %s\n", filename);
        break;
      }
    }
    sc_array_destroy_null (&pdata);

    /* finish meta information and close files */
    if (!g->scaling) {
      if (p4est_vtk_write_footer (cont)) {
        P4EST_LERRORF ("Failed to write footer for %s\n", filename);
        break;
      }
    }
  }
  while (0);
  if (pdata != NULL) {
    sc_array_destroy_null (&pdata);
  }
}

typedef struct bu_data
{
  p4est_locidx_t      count;
}
bu_data_t;

static void
buildp_init_default (p4est_t * p4est,
                     p4est_topidx_t which_tree, p4est_quadrant_t * quadrant)
{
  bu_data_t          *bud = (bu_data_t *) quadrant->p.user_data;
  bud->count = 0;
}

static void
buildp_init_add (p4est_t * p4est,
                 p4est_topidx_t which_tree, p4est_quadrant_t * quadrant)
{
  part_global_t      *g = (part_global_t *) p4est->user_pointer;
  bu_data_t          *bud = (bu_data_t *) quadrant->p.user_data;

  /* tell new quadrant how many particles are in it */
  bud->count = g->add_count;
}

typedef struct pa_bitem
{
  p4est_quadrant_t    quad;
  sc_array_t          parr;
}
pa_bitem_t;

static void
buildp_add (part_global_t * g, p4est_build_t * bcon,
            p4est_topidx_t which_tree, pa_bitem_t * bit)
{
  int                 i;
  int                 cid, nhits;
  int                 wx, wy, wz;
  sc_array_t         *arr;
  pa_bitem_t          abita[P4EST_CHILDREN], *bita;

  P4EST_ASSERT (bcon != NULL);
  P4EST_ASSERT (bit != NULL);
  P4EST_ASSERT (p4est_quadrant_is_valid (&bit->quad));
  P4EST_ASSERT (bit->parr.elem_size == sizeof (p4est_locidx_t));

  if (bit->quad.level == g->maxlevel) {

    /*** we are at maximum level, all particles go into quadrant ***/

    g->add_count = (p4est_locidx_t) bit->parr.elem_count;
    p4est_build_add (bcon, which_tree, &bit->quad);
  }
  else {

    /*** search in children and call buildp_add recursively ***/

    /* access parent quadrant */
    loopquad (g, which_tree, &bit->quad, g->lxyz, g->hxyz, g->dxyz);

    /* number of child quadrant */
    cid = 0;

    /* number of child quadrants containing particles */
    nhits = 0;

    /* sort remaining particles into the children */
#ifdef P4_TO_P8
    split_by_coord
      (g, &bit->parr, g->klh, PA_MODE_LOCATE, 2, g->lxyz, g->dxyz);
    for (wz = 0; wz < 2; ++wz) {
#if 0
    }
#endif
#else
    P4EST_ASSERT (g->klh[0] == NULL);
    P4EST_ASSERT (g->klh[1] == NULL);
    g->klh[0] = &bit->parr;
    wz = 0;
#endif
    split_by_coord
      (g, g->klh[wz], g->jlh, PA_MODE_LOCATE, 1, g->lxyz, g->dxyz);
    for (wy = 0; wy < 2; ++wy) {
      split_by_coord
        (g, g->jlh[wy], g->ilh, PA_MODE_LOCATE, 0, g->lxyz, g->dxyz);
      for (wx = 0; wx < 2; ++wx) {
        /* we have a set of particles for child 4 * wz + 2 * wy + wx */
        P4EST_ASSERT (cid == 4 * wz + 2 * wy + wx);
        if ((arr = g->ilh[wx])->elem_count > 0) {

          /* grab next available slot in array */
          bita = &abita[nhits++];
          p4est_quadrant_child (&bit->quad, &bita->quad, cid);
          sc_array_init_count (&bita->parr, sizeof (p4est_locidx_t),
                               arr->elem_count);
          sc_array_paste (&bita->parr, arr);
        }
        cid++;
      }
    }
#ifdef P4_TO_P8
#if 0
    {
#endif
    }
#endif
    P4EST_ASSERT (cid == P4EST_CHILDREN);

#ifndef P4_TO_P8
    g->klh[0] = NULL;
    P4EST_ASSERT (g->klh[1] == NULL);
#endif

    /* go into recursion for non-empty children */
    for (i = 0; i < nhits; ++i) {
      buildp_add (g, bcon, which_tree, &abita[i]);
    }
  }

  /* this call is expected to clean up the array */
  sc_array_reset (&bit->parr);
}

static void
buildp (part_global_t * g, int k)
{
  double              t0_build, t0_pertreeF, t0_pertreeB, t1;
  char                filename[BUFSIZ];
  sc_array_t          inq;
  sc_array_t         *pdata;
  p4est_topidx_t      tt;
  p4est_locidx_t      lpnum, lq;
  p4est_locidx_t      lall, ilem_particles, li;
  p4est_gloidx_t     *pertree;
  p4est_tree_t       *tree;
  p4est_quadrant_t   *quad;
  p4est_build_t      *bcon;
  p4est_vtk_context_t *cont;
  p4est_t            *build;
  qu_data_t          *qud;
  bu_data_t          *bud;
  pa_data_t          *pad;
  pa_bitem_t          abit, *bit = &abit;

  P4EST_ASSERT (g->padata != NULL);
  P4EST_ASSERT (g->add_count == 0);

  /* only output when specified */
  if (g->scaling && g->verylast) {
    /* output on the last time step */
    if (g->build_part <= 0) {
      return;
    }
  }
  else {
    /* so this is the usual way of things: every so many steps */
    if (g->build_part <= 0 || g->build_step <= 0 || k % g->build_step) {
      return;
    }
  }

  /* STATS */
  t0_build = sc_MPI_Wtime ();

  /* iterate through particles to choose the relevant ones for new forest */
  bcon = p4est_build_new (g->p4est, sizeof (bu_data_t),
                          buildp_init_default, g);
  p4est_build_init_add (bcon, buildp_init_add);
  sc_array_init (&inq, sizeof (p4est_locidx_t));
  pad = (pa_data_t *) sc_array_index_begin (g->padata);
  for (lpnum = 0, lall = 0, tt = g->p4est->first_local_tree;
       tt <= g->p4est->last_local_tree; ++tt) {
    tree = p4est_tree_array_index (g->p4est->trees, tt);
    for (lq = 0; lq < (p4est_locidx_t) tree->quadrants.elem_count; ++lq) {

      /* fetch number of particles in quadrant */
      quad = p4est_quadrant_array_index (&tree->quadrants, lq);
      qud = (qu_data_t *) quad->p.user_data;
      if ((ilem_particles = qud->u.lpend - lpnum) == 0) {
        /* no particles in this quadrant at all */
        continue;
      }

      /* how many particles will we consider? */
      P4EST_ASSERT (inq.elem_count == 0);
      for (li = 0; li < ilem_particles; ++li, ++lall, ++pad) {
        if (!(pad->id % g->build_part)) {
          *(p4est_locidx_t *) sc_array_push (&inq) = lall;
        }
      }
      if (inq.elem_count > 0) {
        bit->quad = *quad;
        sc_array_swap_init (&bit->parr, &inq, sizeof (p4est_locidx_t));
        buildp_add (g, bcon, tt, bit);
      }

      /* move to next quadrant */
      lpnum = qud->u.lpend;
      P4EST_ASSERT (lpnum == lall);
    }
  }
  P4EST_ASSERT (lall == (p4est_locidx_t) g->padata->elem_count);
  P4EST_ASSERT (pad == (pa_data_t *) sc_array_index_end (g->padata));
  P4EST_ASSERT (inq.elem_count == 0 && inq.array == NULL);

  /* create a temporary sparse forest */
  cont = NULL;
  pdata = NULL;
  build = p4est_build_complete (bcon);
  if (!g->scaling || g->verylast) {
    P4EST_GLOBAL_ESSENTIALF ("Built forest with %lld quadrants from %lld\n",
                             (long long) build->global_num_quadrants,
                             (long long) g->p4est->global_num_quadrants);
  }

  /* STATS */
  t1 = sc_MPI_Wtime ();
  if (!g->scaling || g->verylast) {
    sc_stats_accumulate (g->si + PART_STATS_BUILD, t1 - t0_build);
  }

  /* count per-tree quadrants */
  pertree = P4EST_ALLOC (p4est_gloidx_t, g->conn->num_trees + 1);

  /* STATS */
  t0_pertreeF = sc_MPI_Wtime ();

  /* per-tree counts of full forest */
  p4est_comm_count_pertree (g->p4est, pertree);

  /* STATS */
  t1 = sc_MPI_Wtime ();
  if (!g->scaling || g->verylast) {
    sc_stats_accumulate (g->si + PART_STATS_PERTREEF, t1 - t0_pertreeF);
  }
  t0_pertreeB = sc_MPI_Wtime ();

  /* per-tree counts of full forest */
  p4est_comm_count_pertree (build, pertree);

  /* STATS */
  t1 = sc_MPI_Wtime ();
  if (!g->scaling || g->verylast) {
    sc_stats_accumulate (g->si + PART_STATS_PERTREEB, t1 - t0_pertreeB);
  }

  /* clean up per-tree counts */
  P4EST_FREE (pertree);
  pertree = NULL;

  /* write temporary sparse forest to disk */
  do {
    /* open files for output */
    snprintf (filename, BUFSIZ, "%s_W_%06d", g->prefix, k);
    if (!g->scaling) {
      cont = p4est_vtk_context_new (build, filename);
      if (NULL == p4est_vtk_write_header (cont)) {
        P4EST_LERRORF ("Failed to write header for %s\n", filename);
        break;
      }
    }

    /* prepare cell data for output */
    pdata = sc_array_new_count (sizeof (double), build->local_num_quadrants);
    for (lpnum = 0, lall = 0, tt = build->first_local_tree;
         tt <= build->last_local_tree; ++tt) {
      tree = p4est_tree_array_index (build->trees, tt);
      for (lq = 0; lq < (p4est_locidx_t) tree->quadrants.elem_count; ++lq) {

        /* fetch number of particles in quadrant */
        quad = p4est_quadrant_array_index (&tree->quadrants, lq);
        bud = (bu_data_t *) quad->p.user_data;
        *(double *) sc_array_index (pdata, lall++) = (double) bud->count;

        /* move to next quadrant */
        lpnum = qud->u.lpend;
      }
    }

    /* write cell data to file */
    if (!g->scaling) {
      if (NULL == p4est_vtk_write_cell_dataf
          (cont, 1, 1, 1, g->build_wrap, 1, 0, "particles", pdata, cont)) {
        P4EST_LERRORF ("Failed to write cell data for %s\n", filename);
        break;
      }
    }
    sc_array_destroy_null (&pdata);

    /* finish meta information and close files */
    if (!g->scaling) {
      if (p4est_vtk_write_footer (cont)) {
        P4EST_LERRORF ("Failed to write footer for %s\n", filename);
        break;
      }
    }
  }
  while (0);
  if (pdata != NULL) {
    sc_array_destroy_null (&pdata);
  }
  p4est_destroy (build);
  g->add_count = 0;
}

static void
sim (part_global_t * g)
{
  int                 k;
  double              t, h, f;
  double              t0_rk, t1;
  p4est_topidx_t      tt;
  p4est_locidx_t      lpnum, lq;
  p4est_locidx_t      li, ilem_particles;
  p4est_tree_t       *tree;
  p4est_quadrant_t   *quad;
  qu_data_t          *qud;
  pa_data_t          *pad;

  P4EST_ASSERT (g->padata != NULL);
  P4EST_ASSERT (g->stage == 0);

  /* output initial condition */
  t = 0.;
  k = 0;
  pprint (g, t);
  outp (g, k);
  buildp (g, k);
  g->verylast = 0;

  /*** loop over simulation time ***/
  while (t < g->finaltime) {
    h = g->deltat;
    f = t + h;
    if (f > g->finaltime - 1e-3 * g->deltat) {
      f = g->finaltime;
      h = f - t;
      P4EST_ASSERT (!g->verylast);
      g->verylast = 1;
      P4EST_GLOBAL_ESSENTIALF
        ("Last time step %d with %lld particles and %lld quadrants\n", k,
         (long long) g->gpnum, (long long) g->p4est->global_num_quadrants);
    }
    P4EST_GLOBAL_STATISTICSF ("Time %g into step %d with %g\n", t, k, h);

    /*** loop over Runge Kutta stages ***/
    for (g->stage = 0; g->stage < g->order; ++g->stage) {

      /* if a stage is not the last compute new evaluation location */
      /* for the last stage compute the new location of the particle */
      /* do parallel transfer at end of each stage */

      /*** time step for local particles ***/

      /* STATS */
      t0_rk = sc_MPI_Wtime ();

      pad = (pa_data_t *) sc_array_index_begin (g->padata);
      if (pad != NULL) {
        lpnum = 0;
        for (tt = g->p4est->first_local_tree; tt <= g->p4est->last_local_tree;
             ++tt) {
          tree = p4est_tree_array_index (g->p4est->trees, tt);
          for (lq = 0; lq < (p4est_locidx_t) tree->quadrants.elem_count; ++lq) {
            quad = p4est_quadrant_array_index (&tree->quadrants, lq);
            qud = (qu_data_t *) quad->p.user_data;
            ilem_particles = qud->u.lpend - lpnum;

            /*** loop through particles in this quadrant */
            for (li = 0; li < ilem_particles; ++li) {
              /* one Runge Kutta stage for this particle */
              rkstage (g, pad++, h);
            }

            /* move to next quadrant */
            lpnum = qud->u.lpend;
          }
        }
      }

      /* STATS */
      t1 = sc_MPI_Wtime ();
      if (!g->scaling || g->verylast) {
        sc_stats_accumulate (g->si + PART_STATS_RK, t1 - t0_rk);
      }

      /* begin loop -- currently there is no loop */

      /* find new local quadrant or process for each particle */
      presearch (g);

      /* send leaving particles to receiver processes */
      pack (g);
      comm (g);

      /* process remaining local and newly received particles */
      postsearch (g);
      adapt (g);
      regroup (g);

      /* wait for sent messages to complete */
      waitmpi (g);

      /* if no refinement occurred, store received particles and break loop */

      /* partition weighted by current + received count (?) */
      /* partition the forest and send particles along with the partition */
      /* transfer particles accordingly */
      /* think about partitioning and sending particles to future owner? */
      part (g);

      /* end loop -- currently there is no loop */

      if (g->gpnum == 0) {
        /* we have run out of particles */
        break;
      }
    }

    /*** finish up time step ***/
    ++k;
    t = f;

    /* maybe we have run out of particles and stop */
    if (g->gpnum == 0) {
      P4EST_GLOBAL_PRODUCTIONF
        ("We have lost all %lld particles\n", (long long) g->gplost);
      break;
    }
    P4EST_ASSERT (g->stage == g->order);
    g->stage = 0;

    /* write output files */
    pprint (g, t);
    outp (g, k);
    buildp (g, k);
  }

  P4EST_GLOBAL_ESSENTIALF
    ("Time %g is final after %d steps lost %lld remain %lld\n", t, k,
     (long long) g->gplost, (long long) g->gpnum);
}

static void
notif (part_global_t * g)
{
  int                 mpiret;
  int                 i, j;
  int                 q, r;
  double              t0_binary, t0_nary, t1;
  sc_array_t         *recv1, *send1, *payl1;
  sc_array_t         *recv2, *send2, *payl2;
  sc_notify_t        *notifyc;

  recv1 = sc_array_new (sizeof (int));
  send1 = sc_array_new (sizeof (int));
  payl1 = sc_array_new (sizeof (int));

  recv2 = sc_array_new (sizeof (int));
  send2 = sc_array_new (sizeof (int));
  payl2 = sc_array_new (sizeof (int));

  /* send to 7 ranks to the left, 11 apart */
  for (i = 0; i < 7; ++i) {
    q = g->mpirank - (7 - i) * 11;
    if (q >= 0) {
      *(int *) sc_array_push (recv1) = q;
      *(int *) sc_array_push (payl1) = g->mpirank % 19;
    }
  }
  /* send to 5 ranks to the right, 13 apart */
  for (i = 0; i < 5; ++i) {
    q = g->mpirank + (i + 1) * 13;
    if (q < g->mpisize) {
      *(int *) sc_array_push (recv1) = q;
      *(int *) sc_array_push (payl1) = g->mpirank % 19;
    }
  }
  sc_array_copy (recv2, recv1);
  sc_array_copy (payl2, payl1);

  /* allocate notify controller */
  notifyc = sc_notify_new (g->mpicomm);
  sc_notify_set_type (notifyc, SC_NOTIFY_NARY);

  mpiret = sc_MPI_Barrier (g->mpicomm);
  SC_CHECK_MPI (mpiret);

  /* STATS */
  t0_binary = sc_MPI_Wtime ();

  sc_notify_nary_set_widths (notifyc, 2, 2, 2);
  sc_notify_payload (recv1, send1, payl1, NULL, 1, notifyc);

  /* STATS */
  t1 = sc_MPI_Wtime ();
  sc_stats_accumulate (g->si + PART_STATS_BINARY, t1 - t0_binary);

  mpiret = sc_MPI_Barrier (g->mpicomm);
  SC_CHECK_MPI (mpiret);

  /* STATS */
  t0_nary = sc_MPI_Wtime ();

  sc_notify_nary_set_widths (notifyc, g->ntop, g->nint, g->nbot);
  sc_notify_payload (recv2, send2, payl2, NULL, 1, notifyc);

  /* STATS */
  t1 = sc_MPI_Wtime ();
  sc_stats_accumulate (g->si + PART_STATS_NARY, t1 - t0_nary);

  SC_CHECK_ABORT (send1->elem_count == payl1->elem_count, "Send Payl 1");
  SC_CHECK_ABORT (send2->elem_count == payl2->elem_count, "Send Payl 2");

  SC_CHECK_ABORT (sc_array_is_equal (send1, send2), "Send 1 2");
  SC_CHECK_ABORT (sc_array_is_equal (payl1, payl2), "Payl 1 2");

  j = 0;
  for (i = 0; i < 5; ++i) {
    q = g->mpirank - (5 - i) * 13;
    if (q >= 0) {
      SC_CHECK_ABORT (q == *(int *) sc_array_index_int (send1, j), "Left q");
      r = *(int *) sc_array_index_int (payl1, j);
      SC_CHECK_ABORT (q % 19 == r, "Left r");
      ++j;
    }
  }
  for (i = 0; i < 7; ++i) {
    q = g->mpirank + (i + 1) * 11;
    if (q < g->mpisize) {
      SC_CHECK_ABORT (q == *(int *) sc_array_index_int (send1, j), "Right q");
      r = *(int *) sc_array_index_int (payl1, j);
      SC_CHECK_ABORT (q % 19 == r, "Right r");
      ++j;
    }
  }
  SC_CHECK_ABORT (j == (int) send1->elem_count, "Count j");

  sc_notify_destroy (notifyc);
  sc_array_destroy (recv1);
  sc_array_destroy (send1);
  sc_array_destroy (payl1);
  sc_array_destroy (recv2);
  sc_array_destroy (send2);
  sc_array_destroy (payl2);
}

static void
run (part_global_t * g)
{
  pi_data_t           spiddata, *piddata = &spiddata;

  /*** initialize variables ***/
  run_pre (g, piddata);

  /*** initial mesh for domain ***/
  if (g->bricklev > 0) {
    g->conn = p4est_connectivity_new_brick (g->bricklength, g->bricklength
#ifdef P4_TO_P8
                                            , g->bricklength
#endif
                                            , 1, 1
#ifdef P4_TO_P8
                                            , 1
#endif
      );
  }
  else {
#ifndef P4_TO_P8
    g->conn = p4est_connectivity_new_unitsquare ();
#else
    g->conn = p8est_connectivity_new_unitcube ();
#endif
  }
  g->p4est = p4est_new_ext (g->mpicomm, g->conn, 0,
                            g->minlevel - g->bricklev, 1,
                            sizeof (qu_data_t), NULL, g);

  /*** initial refinement and partition ***/
  initrp (g);

  /*** create particles ***/
  g->padata = sc_array_new (sizeof (pa_data_t));
  create (g);

  /*** run simulation ***/
  sim (g);
  sc_array_destroy_null (&g->padata);

  /*** destroy mesh ***/
  p4est_destroy (g->p4est);
  g->p4est = NULL;
  p4est_connectivity_destroy (g->conn);
  g->conn = NULL;

  /*** extra timings for notify ***/
  notif (g);

  /*** clean up variables ***/
  run_post (g);
}

static int
usagerr (sc_options_t * opt, const char *msg)
{
  SC_GLOBAL_LERRORF ("Usage required: %s\n", msg);
  return 1;
}

int
main (int argc, char **argv)
{
  int                 mpiret;
  int                 first_argc;
  int                 ue;
  const char         *opt_notify, *opt_vtk, *opt_build;
  sc_options_t       *opt;
  part_global_t global, *g = &global;

  /*** setup mpi environment ***/

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  memset (g, 0, sizeof (*g));
  g->mpicomm = sc_MPI_COMM_WORLD;
  mpiret = sc_MPI_Comm_size (g->mpicomm, &g->mpisize);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (g->mpicomm, &g->mpirank);
  SC_CHECK_MPI (mpiret);
  sc_init (g->mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  /*** read command line parameters ***/

  opt = sc_options_new (argv[0]);
  sc_options_add_int (opt, 'l', "minlevel", &g->minlevel, 0, "Lowest level");
  sc_options_add_int (opt, 'L', "maxlevel", &g->maxlevel,
                      P4EST_QMAXLEVEL, "Highest level");
  sc_options_add_int (opt, 'b', "bricklev", &g->bricklev,
                      0, "Brick refinement level");
  sc_options_add_int (opt, 'r', "rkorder", &g->order,
                      1, "Order of Runge Kutta method");
  sc_options_add_bool (opt, 'p', "olap-notify", &g->olap_notify, 0,
                       "Overlap sending with notify");
  sc_options_add_string (opt, 'Y', "notify", &opt_notify, NULL,
                         "Notify ntop:nint:nbot");
  sc_options_add_double (opt, 'n', "particles", &g->num_particles,
                         1e3, "Global number of particles");
  sc_options_add_double (opt, 'e', "pperelem", &g->elem_particles,
                         P4EST_CHILDREN, "Number of particles per quadrant");
  sc_options_add_double (opt, 'h', "deltat", &g->deltat,
                         1e-1, "Time step size");
  sc_options_add_double (opt, 'T', "finaltime", &g->finaltime,
                         1., "Final time of simulation");
  sc_options_add_int (opt, 'R', "printn", &g->printn, 0,
                      "Print every nth particle");
  sc_options_add_string (opt, 'V', "vtk", &opt_vtk, NULL,
                         "VTK output everystep:wraprank");
  sc_options_add_string (opt, 'W', "build", &opt_build, NULL,
                         "Build output everystep:particle:wrap");
  sc_options_add_bool (opt, 'S', "scaling", &g->scaling, 0,
                       "Configure for scaling test");
  sc_options_add_bool (opt, 'C', "collapse", &g->collapse, 0,
                       "Collapse statistics over time");
#if 0
  sc_options_add_int (opt, 'C', "checkp", &g->checkp, 0,
                      "write checkpoint output");
#endif
  sc_options_add_string (opt, 'P', "prefix", &g->prefix,
                         "p" PARTICLES_48 ()"rticles",
                         "prefix for file output");

  /* proceed in run-once loop for clean abort */
  ue = 0;
  do {
    /* read command line and assign variables */
    first_argc = sc_options_parse (p4est_package_id, SC_LP_DEFAULT,
                                   opt, argc, argv);
    if (first_argc < 0 || first_argc != argc) {
      ue = usagerr (opt, "Invalid option format or non-option arguments");
      break;
    }
    g->ntop = g->nint = g->nbot = 2;
    part_string_to_int (opt_notify, 3, &g->ntop, &g->nint, &g->nbot);
    part_string_to_int (opt_vtk, 2, &g->vtk, &g->mpiwrap);
    part_string_to_int (opt_build,
                        3, &g->build_step, &g->build_part, &g->build_wrap);

    /* report option and parameter values */
    P4EST_GLOBAL_ESSENTIALF ("Dimension is %d\n", P4EST_DIM);
    sc_options_print_summary (p4est_package_id, SC_LP_ESSENTIAL, opt);
    P4EST_GLOBAL_ESSENTIALF ("Notify parameters: %d %d %d\n",
                             g->ntop, g->nint, g->nbot);
    P4EST_GLOBAL_ESSENTIALF ("VTK parameters: %d %d\n", g->vtk, g->mpiwrap);
    P4EST_GLOBAL_ESSENTIALF ("Build parameters: %d %d %d\n",
                             g->build_step, g->build_part, g->build_wrap);

    /* check options for consistency */
    if (g->minlevel < 0 || g->minlevel > P4EST_QMAXLEVEL) {
      ue = usagerr (opt, "Minlevel between 0 and P4EST_QMAXLEVEL");
    }
    if (g->maxlevel < g->minlevel || g->maxlevel > P4EST_QMAXLEVEL) {
      ue = usagerr (opt, "Maxlevel between minlevel and P4EST_QMAXLEVEL");
    }
    if (g->bricklev < 0 || g->bricklev > g->minlevel) {
      ue = usagerr (opt, "Brick level between 0 and minlevel");
    }
    if (g->order < 1 || g->order > 4) {
      ue = usagerr (opt, "Runge Kutta order between 1 and 4");
    }
    if (g->ntop < 2 || g->nint < 2 || g->nbot < 2) {
      ue = usagerr (opt, "Notify parameters greater equal 2");
    }
    if (g->num_particles <= 0.) {
      ue = usagerr (opt, "Global number of particles positive");
    }
    if (g->elem_particles <= 0.) {
      ue = usagerr (opt, "Number of particles per quadrant positive");
    }
    if (g->printn < 0) {
      ue = usagerr (opt, "Particle print interval non-negative");
    }
    if (g->vtk < 0 || g->mpiwrap < 0) {
      ue = usagerr (opt, "VTK output interval and wrap non-negative");
    }
    if (g->build_step < 0 || g->build_part < 0 || g->build_wrap < 0) {
      ue = usagerr (opt, "Build intervals and wrap non-negative");
    }
    if (ue) {
      break;
    }

    /*** run program ***/

    if (g->scaling) {
      sc_package_set_verbosity (sc_package_id, SC_LP_PRODUCTION);
      sc_package_set_verbosity (p4est_package_id, SC_LP_PRODUCTION);
    }
    run (g);
  }
  while (0);
  if (ue) {
    sc_options_print_usage (p4est_package_id, SC_LP_ERROR, opt, NULL);
  }

  /*** clean up and exit ***/

  sc_options_destroy (opt);
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return ue ? EXIT_FAILURE : EXIT_SUCCESS;
}
