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

#include <sc_functions.h>
#include <sc_options.h>
#ifndef P4_TO_P8
#include <p4est_bits.h>
#include <p4est_extended.h>
#include <p4est_search.h>
#else
#include <p8est_bits.h>
#include <p8est_extended.h>
#include <p8est_search.h>
#endif /* P4_TO_P8 */
#include "spheres_global.h"

#if 0
#define SPHERES_CHATTY
#endif

#define SPHERES_xstr(s) SPHERES_str(s)
#define SPHERES_str(s) #s
#define SPHERES_48() SPHERES_xstr(P4EST_CHILDREN)

static const double irootlen = 1. / P4EST_ROOT_LEN;

static void
create_forest (spheres_global_t * g)
{
  int                 mpiret;
  double              rmin, rmax, rgeom2;
  double              coef, fact, vmult;
  double              Vexp, Nexp, r;
  double              vdensity, sumrd, gsrd;
#if 0
  int                 iscount;
  double              dnu, dnq;
  double              smin, smax, sexp, sdiv, s;
  size_t              zz, zno;
  sc_array_t         *points;
  sc_array_t         *notif;
  sc_array_t         *npayl;
  sc_array_t         *orecs;
  sc_MPI_Request     *sreqs, *rreqs;
#endif
  p4est_topidx_t      which_tree;
  p4est_locidx_t      ntrel, tin;
  p4est_locidx_t      qunsph;
  p4est_locidx_t      sph_excl, sph_incl;
  p4est_locidx_t      li;
  p4est_gloidx_t      lnsph, gnsph;
#if 0
  p4est_locidx_t      spoffs;
  p4est_gloidx_t      lml[5], lmg[5];
#endif
  p4est_qcoord_t      qh;
  p4est_tree_t       *tree;
  p4est_quadrant_t   *q;
  p4est_sphere_t     *sph;

  /* create empty initial forest */
  g->conn = p4est_connectivity_new_periodic ();
  g->p4est = p4est_new_ext (g->mpicomm, g->conn, 0, g->minlevel, 1,
                            sizeof (qu_data_t), NULL, g);

  /* minimum and maximum radius determined by target levels */
  rmax = g->rmax;
  rmin = .5 * g->spherelems * sc_intpowf (.5, g->maxlevel);
  rmin = SC_MIN (rmin, rmax);
  P4EST_ASSERT (rmin > 0.);
  P4EST_GLOBAL_PRODUCTIONF
    ("Generating spheres with radii %g to %g\n", rmin, rmax);

  /* expected volume */
  rgeom2 = rmin * rmax;
#ifndef P4_TO_P8
  Vexp = 4. * rgeom2;
  coef = 1. - rmin / rmax;
#else
  Vexp = 8. * rgeom2 * rgeom2 / (.5 * (rmin + rmax));
  coef = 1. - rmin * rmin / (rmax * rmax);
#endif

  /* the multiplication factor is rescaled to compare to a square/cube */
  vdensity = P4EST_DIM_POW (g->lfraction);
  vmult = vdensity / Vexp;

  /* generate the spheres on minimum refinement level */
  sumrd = 0.;
  sph_excl = sph_incl = 0;
  g->sphr = sc_array_new (sizeof (p4est_sphere_t));
  for (which_tree = g->p4est->first_local_tree;
       which_tree <= g->p4est->last_local_tree; ++which_tree) {
    tree = p4est_tree_array_index (g->p4est->trees, which_tree);
    ntrel = (p4est_locidx_t) tree->quadrants.elem_count;
    P4EST_ASSERT (ntrel > 0);
    for (tin = 0; tin < ntrel; ++tin) {

      /* initialize random number generator anew for each element */
      q = p4est_quadrant_array_index (&tree->quadrants, tin);
      p4est_quadrant_srand (q, &g->rstate);

      /* number of spheres relative to volume ratio of quadrant to tree */
      Nexp = vmult * sc_intpowf (.5, P4EST_DIM * q->level);
      qunsph = (p4est_locidx_t) sc_rand_poisson (&g->rstate, Nexp);

      /* generate spheres for this element */
      if (qunsph > 0) {
        qh = P4EST_QUADRANT_LEN (q->level);
        sph_incl = sph_excl + qunsph;
        sph = (p4est_sphere_t *) sc_array_push_count (g->sphr, qunsph);
        for (li = sph_excl; li < sph_incl; ++li, ++sph) {
          sph->center[0] = (q->x + qh * sc_rand (&g->rstate)) * irootlen;
          sph->center[1] = (q->y + qh * sc_rand (&g->rstate)) * irootlen;
#ifdef P4_TO_P8
          sph->center[2] = (q->z + qh * sc_rand (&g->rstate)) * irootlen;
#endif
          if (coef <= 0.) {
            /* this happens if rmin == rmax: no need to sample */
            r = rmin;
          }
          else {
            fact = 1. / (1. - coef * sc_rand (&g->rstate));
#ifndef P4_TO_P8
            r = rmin * fact;
#else
            r = rmin * sqrt (fact);
#endif
          }
          sph->radius = r;
          sumrd += P4EST_DIM_POW (2. * r);
        }
        sph_excl = sph_incl;
      }
    }
  }
  P4EST_ASSERT (sph_incl == sph_excl);
  P4EST_ASSERT (sph_incl == (p4est_locidx_t) g->sphr->elem_count);
  g->lsph = sph_incl;

  /* get globally unique numbers for the spheres */
  lnsph = (p4est_gloidx_t) g->lsph;
  mpiret = sc_MPI_Exscan (&lnsph, &gnsph, 1, P4EST_MPI_GLOIDX,
                          sc_MPI_SUM, g->mpicomm);
  SC_CHECK_MPI (mpiret);
  if (g->mpirank == 0) {
    gnsph = 0;                  /* value undefined so far */
  }
  if (g->mpirank == g->mpisize - 1) {
    P4EST_PRODUCTIONF
      ("Sphere expected volume %g ideal count %g generated %lld\n",
       Vexp, vmult * g->conn->num_trees, (long long) (gnsph + lnsph));
  }
  g->gsoff = gnsph;

  /* confirm expected volume */
  mpiret = sc_MPI_Allreduce (&sumrd, &gsrd, 1, sc_MPI_DOUBLE,
                             sc_MPI_SUM, g->mpicomm);
  SC_CHECK_MPI (mpiret);
  P4EST_GLOBAL_PRODUCTIONF ("Total volume ideal %g achieved %g\n",
                            vdensity * g->conn->num_trees, gsrd);
}

static int
spheres_partition_quadrant (p4est_t * p4est, p4est_topidx_t which_tree,
                            p4est_quadrant_t * quadrant, int pfirst,
                            int plast, void *point)
{
  spheres_global_t   *g = (spheres_global_t *) p4est->user_pointer;

  P4EST_ASSERT (g != NULL && g->p4est == p4est);
  P4EST_ASSERT (0 <= pfirst && pfirst <= plast && plast < g->mpisize);
  P4EST_ASSERT (point == NULL);

  /* we are not trying to find local spheres */
  if (pfirst == plast && pfirst == g->mpirank) {
    return 0;
  }

  /* record quadrant's position and size */
  p4est_quadrant_sphere_box (quadrant, &g->box);
  return 1;
}

static int
spheres_partition_point (p4est_t * p4est, p4est_topidx_t which_tree,
                         p4est_quadrant_t * quadrant, int pfirst, int plast,
                         void *point)
{
  spheres_global_t   *g = (spheres_global_t *) p4est->user_pointer;
  p4est_locidx_t      li;
  p4est_sphere_t     *sph;
  sph_item_t         *item;
  sr_buf_t           *to_proc;

  P4EST_ASSERT (g != NULL && g->p4est == p4est);
  P4EST_ASSERT (0 <= pfirst && pfirst <= plast && plast < g->mpisize);
  P4EST_ASSERT (point != NULL);
  P4EST_ASSERT (g->box.radius == .5 * P4EST_QUADRANT_LEN (quadrant->level));

  li = *(p4est_locidx_t *) point;
  sph = (p4est_sphere_t *) sc_array_index_int (g->sphr, li);

  /* we may be up in the tree branches */
  if (pfirst < plast) {
    if (!p4est_sphere_match_approx (&g->box, sph, g->thickness)) {
#ifdef SPHERES_CHATTY
      P4EST_INFOF ("No approx match in %d %d for %ld\n",
                   pfirst, plast, (long) li);
#endif
      return 0;
    }

    /* an optimistic match is good enough when we're walking the tree */
    return 1;
  }

  /* we have found the partition of one remote process */
  P4EST_ASSERT (pfirst == plast);
  P4EST_ASSERT (pfirst != g->mpirank);
  if (!p4est_sphere_match_exact (&g->box, sph, g->thickness)) {
#ifdef SPHERES_CHATTY
    P4EST_INFOF ("No exact match in %d %d for %ld\n",
                 pfirst, plast, (long) li);
#endif
    return 0;
  }

#ifdef SPHERES_CHATTY
  P4EST_INFOF ("Found in %d %d: %ld\n", pfirst, plast, (long) li);
#endif

  /* access send buffer for remote process */
  if (pfirst != g->last_to_rank) {
    P4EST_ASSERT (g->last_to_rank < pfirst);
    g->last_to_rank = pfirst;
    g->last_to_proc = to_proc = (sr_buf_t *) sc_array_push (g->to_procs);
    to_proc->rank = pfirst;
    to_proc->items = sc_array_new_count (sizeof (sph_item_t), 1);
    item = (sph_item_t *) sc_array_index (to_proc->items, 0);
  }
  else {
    to_proc = g->last_to_proc;
    P4EST_ASSERT (to_proc->rank == pfirst);
    P4EST_ASSERT (to_proc == (sr_buf_t *)
                  sc_array_index (g->to_procs, g->to_procs->elem_count - 1));
    P4EST_ASSERT (to_proc->items->elem_count > 0);
    item = (sph_item_t *) sc_array_push (to_proc->items);
  }

  /* pack sphere into send buffer */
  item->sph = *sph;
  item->gid = g->gsoff + (p4est_gloidx_t) li;

  /* this return value is ignored */
  return 0;
}

static void
refine_spheres (spheres_global_t * g)
{
  size_t              zz;
  p4est_locidx_t      li;
  sc_array_t         *points;
  sr_buf_t           *proc;

  /*---------------- search partition to find receivers --------------*/

  /* prepare data structures for sending spheres */
  g->last_to_proc = NULL;
  g->last_to_rank = -1;
  g->to_procs = sc_array_new (sizeof (sr_buf_t));

  /* search for remote quadrants that receive spheres */
  points = sc_array_new_count (sizeof (p4est_locidx_t), g->lsph);
  for (li = 0; li < g->lsph; ++li) {
    *(p4est_locidx_t *) sc_array_index_int (points, li) = li;
  }
#ifdef SPHERES_CHATTY
  P4EST_INFOF ("Searching partition for %ld local spheres\n", (long) g->lsph);
#endif
  p4est_search_partition (g->p4est, 0, spheres_partition_quadrant,
                          spheres_partition_point, points);
  sc_array_destroy_null (&points);

  /* send the spheres */

  /* reverse communication pattern */

  /* receive the spheres */

  for (zz = 0; zz < g->to_procs->elem_count; ++zz) {
    proc = (sr_buf_t *) sc_array_index (g->to_procs, zz);
    P4EST_ASSERT (proc->rank != g->mpirank);
    P4EST_ASSERT (proc->items != NULL);
    sc_array_destroy (proc->items);
  }
  sc_array_destroy_null (&g->to_procs);

  /* refine based on received spheres */

  /* partition and transfer owned spheres */
}

static void
destroy_forest (spheres_global_t * g)
{
  sc_array_destroy_null (&g->sphr);

  p4est_destroy (g->p4est);
  p4est_connectivity_destroy (g->conn);
}

static void
run (spheres_global_t * g)
{
  create_forest (g);
  refine_spheres (g);

  destroy_forest (g);
}

static int
usagerr (sc_options_t * opt, const char *msg)
{
  P4EST_GLOBAL_LERRORF ("Usage required: %s\n", msg);
  return 1;
}

int
main (int argc, char **argv)
{
  int                 mpiret;
  int                 ue;
  int                 first_argc;
#if 0
  const char         *opt_notify, *opt_vtk, *opt_build;
#endif
  sc_options_t       *opt;
  spheres_global_t global, *g = &global;

  /*** setup mpi environment ***/

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  /*** initialize global data ***/

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
  sc_options_add_int (opt, 'L', "maxlevel", &g->maxlevel, -1,
                      "Highest level");
  sc_options_add_double (opt, 'R', "rmax", &g->rmax, .5, "Max sphere radius");
  sc_options_add_double (opt, 't', "thickness", &g->thickness, .05,
                         "Relative sphere thickness");
  sc_options_add_double (opt, 'f', "lfraction", &g->lfraction, .3,
                         "Length density of spheres");
  sc_options_add_double (opt, 's', "spherelems", &g->spherelems, 1.,
                         "Min elements per sphere diameter");

  sc_options_add_bool (opt, 'S', "scaling", &g->scaling, 0,
                       "Configure for scaling test");

  sc_options_add_string (opt, 'P', "prefix", &g->prefix,
                         "sph" SPHERES_48 ()"res", "Prefix for file output");

  /* proceed in run-once loop for clean abort */
  ue = 0;
  do {
    /*** parse command line and assign configuration variables ***/

    first_argc = sc_options_parse (p4est_package_id, SC_LP_DEFAULT,
                                   opt, argc, argv);
    if (first_argc < 0 || first_argc != argc) {
      ue = usagerr (opt, "Invalid option format or non-option arguments");
      break;
    }
    P4EST_GLOBAL_ESSENTIALF ("Dimension is %d\n", P4EST_DIM);
    sc_options_print_summary (p4est_package_id, SC_LP_ESSENTIAL, opt);

    /*** check consistency of parameters ***/

    if (g->minlevel < 0 || g->minlevel > P4EST_QMAXLEVEL) {
      ue = usagerr (opt, "Minlevel between 0 and P4EST_QMAXLEVEL");
    }
    if (g->maxlevel == -1) {
      g->maxlevel = g->minlevel;
      P4EST_GLOBAL_ESSENTIALF ("Maxlevel set to %d\n", g->maxlevel);
    }
    if (g->maxlevel < g->minlevel || g->maxlevel > P4EST_QMAXLEVEL) {
      ue = usagerr (opt, "Maxlevel between minlevel and P4EST_QMAXLEVEL");
    }
    if (g->rmax <= 0.) {
      ue = usagerr (opt, "Maximum sphere radius positive");
    }
    if (g->lfraction < 0.) {
      ue = usagerr (opt, "Length density non-negative");
    }
    if (g->spherelems < 1.) {
      ue = usagerr (opt, "Elements per sphere no less than 1.");
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

  return ue;
}
