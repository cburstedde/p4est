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
#include <sc_notify.h>
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

typedef enum
{
  SPHERES_TAG_SPHERES = P4EST_COMM_TAG_LAST
}
spheres_tags;

static const double irootlen = 1. / P4EST_ROOT_LEN;

static void
create_forest (spheres_global_t * g)
{
  int                 mpiret;
  int                 p;
  double              rmin, rmax, rgeom2;
  double              coef, fact, vmult;
  double              Vexp, Nexp, r;
  double              vdensity, sumrd, gsrd;
  p4est_topidx_t      which_tree;
  p4est_locidx_t      ntrel, tin;
  p4est_locidx_t      qunsph;
  p4est_locidx_t      sph_excl, sph_incl;
  p4est_locidx_t      li;
  p4est_qcoord_t      qh;
  p4est_gloidx_t      lg, *gval;
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

  /* obtain globally unique numbers for the spheres */
  g->goffsets = sc_array_new_count (sizeof (p4est_gloidx_t), g->mpisize + 1);
  *(p4est_gloidx_t *) sc_array_index (g->goffsets, 0) = 0;
  lg = (p4est_gloidx_t) g->lsph;
  mpiret = sc_MPI_Allgather (&lg, 1, P4EST_MPI_GLOIDX,
                             sc_array_index (g->goffsets, 1),
                             1, P4EST_MPI_GLOIDX, g->mpicomm);
  SC_CHECK_MPI (mpiret);
  lg = 0;
  for (p = 1; p <= g->mpisize; ++p) {
    gval = (p4est_gloidx_t *) sc_array_index (g->goffsets, p);
    *gval = (lg += *gval);
  }
  g->gsoff = *(p4est_gloidx_t *) sc_array_index (g->goffsets, g->mpirank);
  lg = *(p4est_gloidx_t *) sc_array_index (g->goffsets, g->mpisize);
  P4EST_GLOBAL_PRODUCTIONF
    ("Sphere expected volume %g ideal count %g generated %lld\n",
     Vexp, vmult * g->conn->num_trees, (long long) lg);

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

    /* we have found a new receiver process */
    *(int *) sc_array_push (g->notify) = pfirst;
    *(g->last_payload = (int *) sc_array_push (g->payload)) = 1;

    /* establish new send buffer for this process */
    g->last_to_rank = pfirst;
    g->last_to_proc = to_proc = (sr_buf_t *) sc_array_push (g->to_procs);
    to_proc->rank = pfirst;
    to_proc->items = sc_array_new_count (sizeof (sph_item_t), 1);
    item = (sph_item_t *) sc_array_index (to_proc->items, 0);
  }
  else {
    /* pust to existing send buffer for current remote process */
    to_proc = g->last_to_proc;
    P4EST_ASSERT (to_proc->rank == pfirst);
    P4EST_ASSERT (to_proc == (sr_buf_t *)
                  sc_array_index (g->to_procs, g->to_procs->elem_count - 1));
    P4EST_ASSERT (to_proc->items->elem_count > 0);
    item = (sph_item_t *) sc_array_push (to_proc->items);
    ++*g->last_payload;
  }

  /* pack sphere into send buffer */
  item->sph = *sph;
  item->gid = g->gsoff + (p4est_gloidx_t) li;

  /* this return value is ignored */
  return 0;
}

static int
spheres_local_quadrant (p4est_t * p4est, p4est_topidx_t which_tree,
                        p4est_quadrant_t * quadrant, p4est_locidx_t local_num,
                        void *point)
{
  spheres_global_t   *g = (spheres_global_t *) p4est->user_pointer;

  P4EST_ASSERT (g != NULL && g->p4est == p4est);
  P4EST_ASSERT (point == NULL);

  /* record quadrant's position and size */
  p4est_quadrant_sphere_box (quadrant, &g->box);
  return 1;
}

static int
spheres_local_point (p4est_t * p4est, p4est_topidx_t which_tree,
                     p4est_quadrant_t * quadrant, p4est_locidx_t local_num,
                     void *point)
{
  spheres_global_t   *g = (spheres_global_t *) p4est->user_pointer;
  p4est_sphere_t     *sph;

  P4EST_ASSERT (g != NULL && g->p4est == p4est);
  P4EST_ASSERT (point != NULL);

  sph = *(p4est_sphere_t **) point;

  /* we may be up in the tree branches */
  if (local_num < 0) {
    if (!p4est_sphere_match_approx (&g->box, sph, g->thickness)) {
#ifdef SPHERES_CHATTY
      P4EST_INFOF ("No approx match in %ld\n", (long) local_num);
#endif
      return 0;
    }

    /* an optimistic match is good enough when we're walking the tree */
    return 1;
  }

  if (!p4est_sphere_match_exact (&g->box, sph, g->thickness)) {
#ifdef SPHERES_CHATTY
    P4EST_INFOF ("No approx match in %ld\n", (long) local_num);
#endif
    return 0;
  }

#ifdef SPHERES_CHATTY
  P4EST_INFOF ("Found in %ld\n", (long) local_num);
#endif

  /* this return value is ignored */
  return 0;
}

static void
refine_spheres (spheres_global_t * g)
{
  int                 mpiret;
  int                 q;
  int                 num_from_spheres;
  sc_array_t         *points;
  sc_array_t         *pi;
  p4est_locidx_t      li;
  p4est_locidx_t      sri, snum;
#ifdef P4EST_ENABLE_DEBUG
  p4est_gloidx_t      gcur, gnext;
#endif
  p4est_sphere_t    **psph;
  sph_item_t         *item;
  sr_buf_t           *proc;

  /*---------------- search partition to find receivers --------------*/

  /* prepare data structures for sending spheres */
  g->last_to_proc = NULL;
  g->last_to_rank = -1;
  g->to_procs = sc_array_new (sizeof (sr_buf_t));
  g->notify = sc_array_new (sizeof (int));
  g->payload = sc_array_new (sizeof (int));
  g->last_payload = NULL;

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

  /*------------------------ send the spheres ------------------------*/

  P4EST_ASSERT (g->notify->elem_count == g->to_procs->elem_count);
  P4EST_ASSERT (g->notify->elem_count == g->payload->elem_count);
  g->num_to_procs = (int) g->notify->elem_count;
  g->to_requests =
    sc_array_new_count (sizeof (sc_MPI_Request), g->num_to_procs);
  for (q = 0; q < g->num_to_procs; ++q) {
    proc = (sr_buf_t *) sc_array_index_int (g->to_procs, q);
    P4EST_ASSERT (proc->rank != g->mpirank);
    P4EST_ASSERT (proc->rank == *(int *) sc_array_index_int (g->notify, q));
    P4EST_ASSERT (proc->items != NULL);
    pi = proc->items;
    P4EST_ASSERT (pi->elem_size == sizeof (sph_item_t));
    mpiret = sc_MPI_Isend
      (pi->array, pi->elem_count * pi->elem_size, sc_MPI_BYTE,
       proc->rank, SPHERES_TAG_SPHERES, g->mpicomm,
       (sc_MPI_Request *) sc_array_index_int (g->to_requests, q));
    SC_CHECK_MPI (mpiret);
  }

  /*------------------ reverse communication pattern -----------------*/

  sc_notify_ext (g->notify, NULL, g->payload,
                 g->ntop, g->nint, g->nbot, g->mpicomm);

  /*---------------------- receive the spheres -----------------------*/

  P4EST_ASSERT (g->notify->elem_count == g->payload->elem_count);
  g->num_from_procs = (int) g->notify->elem_count;
  g->from_requests =
    sc_array_new_count (sizeof (sc_MPI_Request), g->num_from_procs);
  g->from_procs = sc_array_new_count (sizeof (sr_buf_t), g->num_from_procs);
  for (q = 0; q < g->num_from_procs; ++q) {
    proc = (sr_buf_t *) sc_array_index_int (g->from_procs, q);
    proc->rank = *(int *) sc_array_index_int (g->notify, q);
    P4EST_ASSERT (proc->rank != g->mpirank);
    num_from_spheres = *(int *) sc_array_index_int (g->payload, q);
    pi = proc->items =
      sc_array_new_count (sizeof (sph_item_t), num_from_spheres);
    mpiret = sc_MPI_Irecv
      (pi->array, pi->elem_count * pi->elem_size, sc_MPI_BYTE,
       proc->rank, SPHERES_TAG_SPHERES, g->mpicomm,
       (sc_MPI_Request *) sc_array_index_int (g->from_requests, q));
    SC_CHECK_MPI (mpiret);
  }

  /*--------------- complete receive and first cleanup ---------------*/

  sc_array_destroy_null (&g->payload);
  sc_array_destroy_null (&g->notify);

  mpiret = sc_MPI_Waitall
    (g->num_from_procs, (sc_MPI_Request *) g->from_requests->array,
     sc_MPI_STATUSES_IGNORE);
  SC_CHECK_MPI (mpiret);
  sc_array_destroy_null (&g->from_requests);

  /*---------------- refine based on received spheres ----------------*/

  points = sc_array_new (sizeof (p4est_sphere_t *));
  sri = 0;

  /* ranks less than ours */
  for (q = 0; q < g->num_from_procs; ++q) {
    proc = (sr_buf_t *) sc_array_index_int (g->from_procs, q);
    P4EST_ASSERT (proc->rank != g->mpirank);
    if (proc->rank >= g->mpirank) {
      break;
    }
#ifdef P4EST_ENABLE_DEBUG
    gcur = *(p4est_gloidx_t *) sc_array_index_int (g->goffsets, proc->rank);
    gnext =
      *(p4est_gloidx_t *) sc_array_index_int (g->goffsets, proc->rank + 1);
#endif
    snum = (p4est_locidx_t) proc->items->elem_count;
    psph = (p4est_sphere_t **) sc_array_push_count (points, snum);
    for (li = 0; li < snum; ++li) {
      item = (sph_item_t *) sc_array_index_int (proc->items, li);
      P4EST_ASSERT (gcur <= item->gid && item->gid < gnext);
      *psph++ = &item->sph;
    }
    sri += snum;
  }

  /* our local spheres */
  snum = g->lsph;
  P4EST_ASSERT (snum == (p4est_locidx_t) g->sphr->elem_count);
  psph = (p4est_sphere_t **) sc_array_push_count (points, snum);
  for (li = 0; li < snum; ++li) {
    *psph++ = (p4est_sphere_t *) sc_array_index_int (g->sphr, li);
  }
  sri += snum;

  /* ranks greater than ours */
  for (; q < g->num_from_procs; ++q) {
    proc = (sr_buf_t *) sc_array_index_int (g->from_procs, q);
    P4EST_ASSERT (proc->rank > g->mpirank);
#ifdef P4EST_ENABLE_DEBUG
    gcur = *(p4est_gloidx_t *) sc_array_index_int (g->goffsets, proc->rank);
    gnext =
      *(p4est_gloidx_t *) sc_array_index_int (g->goffsets, proc->rank + 1);
#endif
    snum = (p4est_locidx_t) proc->items->elem_count;
    psph = (p4est_sphere_t **) sc_array_push_count (points, snum);
    for (li = 0; li < snum; ++li) {
      item = (sph_item_t *) sc_array_index_int (proc->items, li);
      P4EST_ASSERT (gcur <= item->gid && item->gid < gnext);
      *psph++ = &item->sph;
    }
    sri += snum;
  }
  P4EST_ASSERT (sri == (p4est_locidx_t) points->elem_count);

  /* search through local elements and set refinement flag */
#ifdef SPHERES_CHATTY
  P4EST_INFOF ("Searching elements for %ld local spheres\n", (long) sri);
#endif
  p4est_search_local (g->p4est, 0, spheres_local_quadrant,
                      spheres_local_point, points);
  sc_array_destroy_null (&points);

  /*-------------- partition and transfer owned spheres --------------*/

  /*---------------- complete send and second cleanup ----------------*/

  for (q = 0; q < g->num_from_procs; ++q) {
    proc = (sr_buf_t *) sc_array_index_int (g->from_procs, q);
    P4EST_ASSERT (proc->rank != g->mpirank);
    P4EST_ASSERT (proc->items != NULL);
    sc_array_destroy (proc->items);
  }
  sc_array_destroy_null (&g->from_procs);

  mpiret = sc_MPI_Waitall
    (g->num_to_procs, (sc_MPI_Request *) g->to_requests->array,
     sc_MPI_STATUSES_IGNORE);
  SC_CHECK_MPI (mpiret);
  sc_array_destroy_null (&g->to_requests);

  for (q = 0; q < g->num_to_procs; ++q) {
    proc = (sr_buf_t *) sc_array_index_int (g->to_procs, q);
    P4EST_ASSERT (proc->rank != g->mpirank);
    P4EST_ASSERT (proc->items != NULL);
    sc_array_destroy (proc->items);
  }
  sc_array_destroy_null (&g->to_procs);
}

static void
destroy_forest (spheres_global_t * g)
{
  sc_array_destroy_null (&g->sphr);
  sc_array_destroy_null (&g->goffsets);

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

#ifdef SPHERES_CHATTY
  P4EST_GLOBAL_INFOF ("Sizeof sphere %u item %u\n",
                      (unsigned) sizeof (p4est_sphere_t),
                      (unsigned) sizeof (sph_item_t));
#endif

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

  /* set other parameters */
  g->ntop = 12;
  g->nint = g->nbot = 4;

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
