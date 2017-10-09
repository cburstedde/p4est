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

typedef struct pi_data
{
  double              sigma;
  double              invs2;
  double              gnorm;
  double              center[3];
}
pi_data_t;

typedef struct qu_data
{
  long long           lpoffset;
  double              d;
}
qu_data_t;

typedef struct pa_data
{
  double              xyz[3];
}
pa_data_t;

static const double simpson[3] = { 1. / 6, 2. / 3., 1. / 6. };

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
    loclp[1] = refine_maxl;
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
srandquad (part_global_t * g, p4est_quadrant_t * quad)
{
  srand ((unsigned int) p4est_quadrant_linear_id (quad, P4EST_QMAXLEVEL));
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
      qud->lpoffset = lpnum;
      ilem_particles =
        (int) round (qud->d / g->global_density * g->num_particles);
      pad = (pa_data_t *) sc_array_push_count (g->padata,
                                               (size_t) ilem_particles);

      /*** generate required number of particles ***/
      loopquad (g, tt, quad, lxyz, hxyz, dxyz);
      srandquad (g, quad);
      for (i = 0; i < ilem_particles; ++i) {
        for (j = 0; j < P4EST_DIM; ++j) {
          r = rand () / (double) RAND_MAX;
          pad->xyz[j] = lxyz[j] + r * dxyz[j];
        }
        ++pad;
      }
      lpnum += (long long) ilem_particles;
    }
  }
  mpiret = sc_MPI_Allreduce (&lpnum, &g->gpnum, 1, sc_MPI_LONG_LONG_INT,
                             sc_MPI_SUM, g->mpicomm);
  SC_CHECK_MPI (mpiret);
  P4EST_GLOBAL_INFOF ("Created %lld particles for %g\n",
                      g->gpnum, g->num_particles);
}

static void
sim (part_global_t * g)
{
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
  sc_init (g->mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  /*** read command line parameters ***/

  opt = sc_options_new (argv[0]);
  sc_options_add_int (opt, 'l', "minlevel", &g->minlevel, 0, "Lowest level");
  sc_options_add_int (opt, 'L', "maxlevel", &g->maxlevel, 0, "Highest level");
  sc_options_add_int (opt, 'b', "bricklev", &g->bricklev,
                      0, "Brick refinement level");
  sc_options_add_double (opt, 'n', "particles", &g->num_particles,
                         1e3, "Global number of particles");
  sc_options_add_double (opt, 'e', "pperelem", &g->elem_particles,
                         3., "Number of particles per element");
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
