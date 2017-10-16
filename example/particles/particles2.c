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
#include <p4est_extended.h>
#include <p4est_search.h>
#else
#include <p8est_bits.h>
#include <p8est_extended.h>
#include <p8est_search.h>
#endif /* P4_TO_P8 */
#include <sc_options.h>
#include "global.h"

#define PARTICLES_xstr(s) PARTICLES_str(s)
#define PARTICLES_str(s) #s
#define PARTICLES_48() PARTICLES_xstr(P4EST_CHILDREN)

#define PART_PLANETS (2)

typedef struct pi_data
{
  double              sigma;
  double              invs2;
  double              gnorm;
  double              center[3];
}
pi_data_t;

/** Data type for payload data inside each quadrant */
typedef union qu_data
{
  /** Offset into local array of all particles after this quadrant */
  long long           lpend;
  double              d;
}
qu_data_t;

/** Property data stored in a flat array over all particles */
typedef struct pa_data
{
  double              xv[6];
  double              wo[6];
  double              up[6];
}
pa_data_t;

/** Search metadata stored in a flat array over all particles */
typedef struct pa_found
{
  p4est_locidx_t      pori;
}
pa_found_t;

static const double simpson[3] = { 1. / 6, 2. / 3., 1. / 6. };

static const double planet_xyz[PART_PLANETS][3] = {
  {.48, .48, .56},
  {.58, .43, .59}
};
static const double planet_mass[PART_PLANETS] = { .1, .2 };

static const double rk1b[0] = { };
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

static double
integrate (part_global_t * g, const double lxyz[3], const double dxyz[3])
{
  int                 i, j, k;
  double              wk, wkj, wkji;
  double              d;

  /*** run Simpson's rule ***/
  d = 0.;
#ifdef P4_TO_P8
  for (k = 0; k < 3; ++k) {
    wk = simpson[k] * dxyz[2];
#if 0
  }
#endif
#else
  k = 0;
  wk = 1.;
#endif
  for (j = 0; j < 3; ++j) {
    wkj = wk * simpson[j] * dxyz[1];
    for (i = 0; i < 3; ++i) {
      wkji = wkj * simpson[i] * dxyz[0];
      d += wkji * g->pidense (lxyz[0] + .5 * i * dxyz[0],
                              lxyz[1] + .5 * j * dxyz[1],
                              lxyz[2] + .5 * k * dxyz[2], g->piddata);
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
  int                 ilem_particles;

  ilem_particles =
    (int) round (qud->d * g->num_particles / g->global_density);

  return (double) ilem_particles > g->elem_particles;
}

static void
initrp (part_global_t * g)
{
  int                 mpiret;
  int                 cycle, max_cycles;
  int                 ilem_particles;
  double              lxyz[3], hxyz[3], dxyz[3];
  double              d, ld;
  double              refine_maxd, refine_maxl;
  double              loclp[2], glolp[2];
  p4est_topidx_t      tt;
  p4est_locidx_t      lq;
  p4est_gloidx_t      old_gnum, new_gnum;
  p4est_tree_t       *tree;
  p4est_quadrant_t   *quad;
  qu_data_t          *qud;

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
        qud->d = d = integrate (g, lxyz, dxyz);
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
    P4EST_GLOBAL_INFOF ("Global integral over density %g\n",
                        g->global_density);

    /*** get global maximum of particle count and level ***/
    loclp[0] = refine_maxd;
    loclp[1] = refine_maxl + g->bricklev;
    mpiret = sc_MPI_Allreduce (loclp, glolp, 2, sc_MPI_DOUBLE,
                               sc_MPI_MAX, g->mpicomm);
    SC_CHECK_MPI (mpiret);
    ilem_particles =
      (int) round (glolp[0] * g->num_particles / g->global_density);
    P4EST_GLOBAL_INFOF ("Maximum particle number per quadrant %d"
                        " and level %g\n", ilem_particles, glolp[1]);

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

    /*** weighted partition ***/
    p4est_partition (g->p4est, 0, NULL);
  }
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
  int                 i, j;
  int                 ilem_particles;
  long long           lpnum;
  double              lxyz[3], hxyz[3], dxyz[3];
  double              r;
  p4est_topidx_t      tt;
  p4est_locidx_t      lq;
  p4est_tree_t       *tree;
  p4est_quadrant_t   *quad;
  qu_data_t          *qud;
  pa_data_t          *pad;

  /*** iterate through local cells and populate with particles ***/
  g->padata = sc_array_new (sizeof (pa_data_t));
  lpnum = 0;
  for (tt = g->p4est->first_local_tree; tt <= g->p4est->last_local_tree; ++tt) {
    tree = p4est_tree_array_index (g->p4est->trees, tt);
    for (lq = 0; lq < (p4est_locidx_t) tree->quadrants.elem_count; ++lq) {
      quad = p4est_quadrant_array_index (&tree->quadrants, lq);
      qud = (qu_data_t *) quad->p.user_data;
      ilem_particles =
        (int) round (qud->d / g->global_density * g->num_particles);
      pad = (pa_data_t *) sc_array_push_count (g->padata,
                                               (size_t) ilem_particles);

      /*** generate required number of particles ***/
      loopquad (g, tt, quad, lxyz, hxyz, dxyz);
      srandquad (g, lxyz);
      for (i = 0; i < ilem_particles; ++i) {
        for (j = 0; j < P4EST_DIM; ++j) {
          r = rand () / (double) RAND_MAX;
          pad->xv[j] = lxyz[j] + r * dxyz[j];
          pad->xv[3 + j] = 0.;
        }
#ifndef P4_TO_P8
        pad->xv[2] = pad->xv[5] = 0.;
#endif
#if 0
        P4EST_LDEBUGF ("Create particle <%g %g %g>\n",
                       pad->xv[0], pad->xv[1], pad->xv[2]);
#endif
        ++pad;
      }
      lpnum += (long long) ilem_particles;
      qud->lpend = lpnum;
    }
  }
  g->gplost = 0;
  mpiret = sc_MPI_Allreduce (&lpnum, &g->gpnum, 1, sc_MPI_LONG_LONG_INT,
                             sc_MPI_SUM, g->mpicomm);
  SC_CHECK_MPI (mpiret);
  P4EST_GLOBAL_INFOF ("Created %lld particles for %g\n",
                      g->gpnum, g->num_particles);
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

  for (j = 0; j < PART_PLANETS; ++j) {
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
  int                 stage = g->stage;
  int                 i;
  double              d;
  double              rk[6];

  /* evaluate right hand side */
  rkrhs (g, stage == 0 ? pad->xv : pad->wo, rk);

  /* compute new evaluation point if necessary */
  if (stage + 1 < g->order) {
    /* stage is not the last */
    d = h * prk[g->order - 1][0][stage];
    for (i = 0; i < 6; ++i) {
      pad->wo[i] = pad->xv[i] + d * rk[i];
    }
  }

  /* compute an update to the state */
  d = prk[g->order - 1][1][stage];
  if (stage == 0) {
    /* first stage */
    if (g->order > 1) {
      /* first stage is not the last */
      P4EST_ASSERT (stage + 1 < g->order);
      for (i = 0; i < 6; ++i) {
        pad->up[i] = d * rk[i];
      }
    }
    else {
      /* first stage is also the last */
      P4EST_ASSERT (stage + 1 == g->order);
      for (i = 0; i < 6; ++i) {
        pad->xv[i] += h * d * rk[i];
      }
    }
  }
  else {
    /* stage is not the first */
    if (stage + 1 < g->order) {
      /* stage is neither first nor last */
      P4EST_ASSERT (0 < stage);
      for (i = 0; i < 6; ++i) {
        pad->up[i] += d * rk[i];
      }
    }
    else {
      /* stage is last of several */
      P4EST_ASSERT (stage + 1 == g->order);
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
  return 1;
}

static int
psearch_point (p4est_t * p4est, p4est_topidx_t which_tree,
               p4est_quadrant_t * quadrant, int pfirst, int plast,
               p4est_locidx_t local_num, void *point)
{
  int                 i;
  size_t              zp;
  double              lxyz[3], hxyz[3], dxyz[3];
  const double       *x;
  part_global_t      *g = (part_global_t *) p4est->user_pointer;
  pa_data_t          *pad = (pa_data_t *) point;
  pa_found_t         *pfn;

  /* access location of particle to be searched */
  P4EST_ASSERT (pad != NULL);
  x = g->stage + 1 < g->order ? pad->wo : pad->xv;

  /* due to roundoff we call this even for a local leaf */
  loopquad (g, which_tree, quadrant, lxyz, hxyz, dxyz);
  for (i = 0; i < P4EST_DIM; ++i) {
    if (!(lxyz[i] <= x[i] && x[i] <= hxyz[i])) {
      /* the point is outside the search quadrant */
      return 0;
    }
  }

  /* find process/quadrant for this particle */
  if (local_num >= 0) {
    /* quadrant is a local leaf */
    P4EST_ASSERT (pfirst == g->mpirank && plast == g->mpirank);
    zp = sc_array_position (g->padata, pad);
    pfn = (pa_found_t *) sc_array_index (g->pfound, zp);
    /* first local match counts (due to roundoff there may be multiple) */
    if (pfn->pori < g->mpisize) {
      pfn->pori = (p4est_locidx_t) g->mpisize + local_num;
      /* TODO: bump counter of particles in this local quadrant */
    }
    /* return value will have no effect */
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
    zp = sc_array_position (g->padata, pad);
    pfn = (pa_found_t *) sc_array_index (g->pfound, zp);
    /* only count match if it has not been found locally or on lower rank */
    if (pfn->pori < 0 || (pfirst < pfn->pori && pfn->pori < g->mpisize)) {
      pfn->pori = (p4est_locidx_t) pfirst;
    }
    /* return value will have no effect */
    return 0;
  }

  /* the process for this particle has not yet been found */
  return 1;
}

static void
sim (part_global_t * g)
{
  int                 i, k;
  int                 ilem_particles;
  long long           lpnum;
  double              t, h, f;
  p4est_topidx_t      tt;
  p4est_locidx_t      lq;
  p4est_tree_t       *tree;
  p4est_quadrant_t   *quad;
  qu_data_t          *qud;
  pa_data_t          *pad;

  /*** loop over simulation time ***/
  k = 0;
  t = 0.;
  while (t < g->finaltime) {
    h = g->deltat;
    f = t + h;
    if (f > g->finaltime - 1e-3 * g->deltat) {
      f = g->finaltime;
      h = f - t;
    }
    P4EST_GLOBAL_INFOF ("Time %g into step %d with %g\n", t, k, h);

    /*** loop over Runge Kutta stages ***/
    for (g->stage = 0; g->stage < g->order; ++g->stage) {

      /* if a stage is not the last compute new evaluation location */
      /* for the last stage compute the new location of the particle */
      /* do parallel transfer at end of each stage */

      /*** time step for local particles ***/
      lpnum = 0;
      pad = (pa_data_t *) sc_array_index (g->padata, 0);
      for (tt = g->p4est->first_local_tree; tt <= g->p4est->last_local_tree;
           ++tt) {
        tree = p4est_tree_array_index (g->p4est->trees, tt);
        for (lq = 0; lq < (p4est_locidx_t) tree->quadrants.elem_count; ++lq) {
          quad = p4est_quadrant_array_index (&tree->quadrants, lq);
          qud = (qu_data_t *) quad->p.user_data;
          ilem_particles = (int) (qud->lpend - lpnum);

          /*** loop through particles in this element */
          for (i = 0; i < ilem_particles; ++i) {
            /* one Runge Kutta stage for this particle */
            rkstage (g, pad++, h);
          }

          /* move to next quadrant */
          lpnum = qud->lpend;
        }
      }

      /* TODO:
         reassign particle to other quadrant
         according to position in either wo (not last stage) or xv;
         notify quadrants of assignment to refine/coarsen
         parallel transfer and partition */

      /* begin loop */

      /* p4est_search_all to find new local element or process for each particle */
      g->pfound =
        sc_array_new_count (sizeof (pa_found_t), g->padata->elem_count);
      sc_array_memset (g->pfound, -1);
      p4est_search_all (g->p4est, 0, psearch_quad, psearch_point, g->padata);
      sc_array_destroy (g->pfound);

      /* send to-be-received particles to receiver processes */

      /* receive particles and run local search to count them per-quadrant */

      /* refine the mesh based on current + received count */

      /* if no refinement occurred, store received particles and break loop */

      /* partition weighted by current + received count (?) */

      /* transfer particles accordingly */

      /* end loop */
    }

    /*** finish up time step ***/
    ++k;
    t = f;
  }

  P4EST_GLOBAL_PRODUCTIONF ("Time %g is final after %d steps\n", t, k);
}

static void
run (part_global_t * g)
{
  int                 b;
  pi_data_t           spiddata, *piddata = &spiddata;

  /*** initial particle density ***/
  piddata->sigma = .1;
  piddata->invs2 = 1. / SC_SQR (piddata->sigma);
  piddata->gnorm = gaussnorm (piddata->sigma);
  piddata->center[0] = .3;
  piddata->center[1] = .4;
#ifndef P4_TO_P8
  piddata->center[2] = 0.;
#else
  piddata->center[2] = .5;
#endif
  g->pidense = pidense;
  g->piddata = piddata;

  /*** initial mesh for domain ***/
  b = g->bricklength = (1 << g->bricklev);
  if (g->bricklev > 0) {
    g->conn = p4est_connectivity_new_brick (b, b
#ifdef P4_TO_P8
                                            , b
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
  create (g);

  /*** run simulation ***/
  sim (g);

  /*** destroy data ***/
  sc_array_destroy (g->padata);

  /*** destroy mesh ***/
  p4est_destroy (g->p4est);
  p4est_connectivity_destroy (g->conn);
}

static int
usagerr (sc_options_t * opt, const char *msg)
{
  SC_GLOBAL_LERRORF ("Usage required: %s\n", msg);
  sc_options_print_usage (p4est_package_id, SC_LP_ERROR, opt, NULL);
  return 1;
}

int
main (int argc, char **argv)
{
  int                 mpiret;
  int                 first_argc;
  sc_options_t       *opt;
  part_global_t global, *g = &global;

  /*** setup mpi environment ***/

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

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
  sc_options_add_int (opt, 'r', "rk", &g->order,
                      1, "Order of Runge Kutta method");
  sc_options_add_double (opt, 'n', "particles", &g->num_particles,
                         1e3, "Global number of particles");
  sc_options_add_double (opt, 'e', "pperelem", &g->elem_particles,
                         3., "Number of particles per element");
  sc_options_add_double (opt, 'h', "deltat", &g->deltat,
                         1e-1, "Time step size");
  sc_options_add_double (opt, 'T', "finaltime", &g->finaltime,
                         1., "Final time of simulation");
  sc_options_add_switch (opt, 'V', "vtk", &g->vtk, "write VTK output");
  sc_options_add_switch (opt, 'C', "check", &g->check,
                         "write checkpoint output");
  sc_options_add_string (opt, 'P', "prefix", &g->prefix,
                         "p" PARTICLES_48 ()"rticles",
                         "prefix for file output");

  first_argc = sc_options_parse (p4est_package_id, SC_LP_DEFAULT,
                                 opt, argc, argv);
  if (first_argc < 0 || first_argc != argc) {
    return usagerr (opt, "No non-option arguments permitted");
  }
  if (g->minlevel < 0 || g->minlevel > P4EST_QMAXLEVEL) {
    return usagerr (opt, "Minlevel between 0 and P4EST_QMAXLEVEL");
  }
  if (g->maxlevel < g->minlevel || g->maxlevel > P4EST_QMAXLEVEL) {
    return usagerr (opt, "Maxlevel between minlevel and P4EST_QMAXLEVEL");
  }
  if (g->bricklev < 0 || g->bricklev > g->minlevel) {
    return usagerr (opt, "Brick level between 0 and minlevel");
  }
  if (g->order < 1 || g->order > 4) {
    return usagerr (opt, "Runge Kutta order between 1 and 4");
  }
  if (g->num_particles <= 0.) {
    return usagerr (opt, "Global number of particles positive");
  }
  if (g->elem_particles <= 0.) {
    return usagerr (opt, "Number of particles per element positive");
  }
  sc_options_print_summary (p4est_package_id, SC_LP_PRODUCTION, opt);
  sc_options_destroy (opt);

  /*** run program ***/

  run (g);

  /*** clean up and exit ***/

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
