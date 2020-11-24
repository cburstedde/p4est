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
#include <p4est_communication.h>
#include <p4est_extended.h>
#include <p4est_search.h>
#include <p4est_vtk.h>
#define SPHERES_DIM_G           "%g %g"
#else
#include <p8est_bits.h>
#include <p8est_communication.h>
#include <p8est_extended.h>
#include <p8est_search.h>
#include <p8est_vtk.h>
#define SPHERES_DIM_G           "%g %g %g"
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
  SPHERES_TAG_SPHERES = P4EST_COMM_TAG_LAST,
  SPHERES_TAG_FIXED,
  SPHERES_TAG_CUSTOM
}
spheres_tags;

typedef enum
{
  STATS_ONCE_FIRST = 0,
  STATS_ONCE_WALL = STATS_ONCE_FIRST,
  STATS_ONCE_NEW,
  STATS_ONCE_GENERATE,
  STATS_ONCE_LAST
}
stats_once_t;

static const char  *once_names[] = { "Walltime", "New", "Generate" };

typedef enum
{
  STATS_LOOP_FIRST = 0,
  STATS_LOOP_PSEARCH = STATS_LOOP_FIRST,
  STATS_LOOP_SEND,
  STATS_LOOP_NOTIFY,
  STATS_LOOP_RECEIVE,
  STATS_LOOP_WAITSEND,
  STATS_LOOP_LSEARCH,
  STATS_LOOP_REFINE,
  STATS_LOOP_PARTITION,
  STATS_LOOP_TRANSFER,
  STATS_LOOP_LAST
}
stats_loop_t;

static const char  *loop_names[] = { "Psearch", "Send", "Notify",
  "Receive", "Waitsend", "Lsearch", "Refine", "Partition", "Transfer"
};

static sc_statinfo_t *
once_index (spheres_global_t * g, stats_once_t once)
{
  P4EST_ASSERT (0 <= once && once < STATS_ONCE_LAST);
  return (sc_statinfo_t *) sc_array_index_int (g->stats, once);
}

static sc_statinfo_t *
loop_index (spheres_global_t * g, int lev, stats_loop_t loop)
{
  P4EST_ASSERT (g->minlevel <= lev && lev < g->maxlevel);
  P4EST_ASSERT (0 <= loop && loop < STATS_LOOP_LAST);

  return (sc_statinfo_t *)
    sc_array_index_int (g->stats, STATS_ONCE_LAST +
                        (lev - g->minlevel) * STATS_LOOP_LAST + loop);
}

static void
accumulate_once (spheres_global_t * g, stats_once_t once, double t)
{
  sc_stats_accumulate (once_index (g, once), sc_MPI_Wtime () - t);
}

static void
accumulate_loop (spheres_global_t * g, int lev, stats_loop_t loop, double t)
{
  sc_stats_accumulate (loop_index (g, lev, loop), sc_MPI_Wtime () - t);
}

static const double irootlen = 1. / P4EST_ROOT_LEN;

static void
spheres_init_zero (p4est_t * p4est, p4est_topidx_t which_tree,
                   p4est_quadrant_t * quadrant)
{
  qu_data_t          *qud = (qu_data_t *) quadrant->p.user_data;
  qud->set_refine = 0;
}

static int
spheres_refine_callback (p4est_t * p4est, p4est_topidx_t which_tree,
                         p4est_quadrant_t * quadrant)
{
  spheres_global_t   *g = (spheres_global_t *) p4est->user_pointer;
  qu_data_t          *qud = (qu_data_t *) quadrant->p.user_data;
  int                 set_refine;

  /* update variables for the non-refining case */
  g->lsph_offset +=
    (*(int *) sc_array_index_int (g->lcounts_refined, g->lqindex_refined) =
     *(int *) sc_array_index_int (g->lcounts, g->lqindex));
  ++g->lqindex;
  ++g->lqindex_refined;

  /* the replace callback updates the variables if refined */
  set_refine = qud->set_refine;
  qud->set_refine = 0;
  return set_refine;
}

#ifdef P4EST_ENABLE_DEBUG

static int
center_in_box (const p4est_sphere_t * box, const p4est_sphere_t * sph)
{
  int                 i;

  for (i = 0; i < P4EST_DIM; ++i) {
    if (fabs (box->center[i] - sph->center[i]) > box->radius) {
      return 0;
    }
  }
  return 1;
}

#endif

static void
spheres_replace_callback (p4est_t * p4est, p4est_topidx_t which_tree,
                          int num_outgoing, p4est_quadrant_t * outgoing[],
                          int num_incoming, p4est_quadrant_t * incoming[])
{
  int                 c;
  double              discr[P4EST_DIM];
  sc_array_t          neworder[P4EST_CHILDREN];
  p4est_locidx_t      outg_spheres, prev_spheres, li;
  p4est_sphere_t     *sph;
#ifdef P4EST_ENABLE_DEBUG
  p4est_sphere_t      box;
#endif
  spheres_global_t   *g = (spheres_global_t *) p4est->user_pointer;

  P4EST_ASSERT (num_outgoing == 1);
  P4EST_ASSERT (num_incoming == P4EST_CHILDREN);

#ifdef P4EST_ENABLE_DEBUG
  p4est_quadrant_sphere_box (outgoing[0], &box);
#endif
  discr[0] = incoming[1]->x * irootlen;
  discr[1] = incoming[2]->y * irootlen;
#ifdef P4_TO_P8
  discr[2] = incoming[4]->z * irootlen;
#endif
  P4EST_ASSERT (!((qu_data_t *) outgoing[0]->p.user_data)->set_refine);
  for (c = 0; c < P4EST_CHILDREN; ++c) {
    sc_array_init (&neworder[c], sizeof (p4est_sphere_t));
    P4EST_ASSERT (!((qu_data_t *) incoming[c]->p.user_data)->set_refine);
  }

  /* we know that the refine callback for this quadrant returned true */
  P4EST_ASSERT (g->lqindex >= 1);
  P4EST_ASSERT (g->lqindex_refined >= 1);
  (void) sc_array_push_count (g->lcounts_refined, P4EST_CHILDREN - 1);
  outg_spheres = *(int *) sc_array_index (g->lcounts, g->lqindex - 1);
  prev_spheres = g->lsph_offset - outg_spheres;
  P4EST_ASSERT (0 <= prev_spheres && prev_spheres <= g->lsph);

  /* sort spheres in child quadrant order */
  for (li = 0; li < outg_spheres; ++li) {
    sph = (p4est_sphere_t *) sc_array_index_int (g->sphr, prev_spheres + li);
    P4EST_ASSERT (center_in_box (&box, sph));
    c = 0;
    c += (sph->center[0] > discr[0]);
    c += (sph->center[1] > discr[1]) << 1;
#ifdef P4_TO_P8
    c += (sph->center[2] > discr[2]) << 2;
#endif
    P4EST_ASSERT (0 <= c && c < P4EST_CHILDREN);
    *(p4est_sphere_t *) sc_array_push (&neworder[c]) = *sph;
  }

  /* assign sort results by child */
  li = prev_spheres;
  for (c = 0; c < P4EST_CHILDREN; ++c) {
    sc_array_copy_into (g->sphr, li, &neworder[c]);
    li += (*(int *) sc_array_index_int (g->lcounts_refined,
                                        g->lqindex_refined - 1 + c) =
           (int) neworder[c].elem_count);
    sc_array_reset (&neworder[c]);
  }
  P4EST_ASSERT (li == g->lsph_offset);

  /* finally update the running offset */
  g->lqindex_refined += P4EST_CHILDREN - 1;
}

static void
spheres_write_vtk (spheres_global_t * g, const char *str, int lev)
{
  char                filename[BUFSIZ];
  sc_array_t         *sdata;
  p4est_topidx_t      tt;
  p4est_locidx_t      lall, lq;
  p4est_tree_t       *tree;
  p4est_vtk_context_t *cont;

  P4EST_ASSERT (0 <= g->minlevel && g->maxlevel <= P4EST_QMAXLEVEL);
  P4EST_ASSERT (g->minlevel <= lev && lev <= g->maxlevel);
  P4EST_ASSERT (g->lcounts != NULL);

  /* write VTK output of sphere counts */
  if (!g->write_vtk) {
    return;
  }

  /* run-once loop for clean return */
  cont = NULL;
  sdata = NULL;
  do {
    /* open files for output */
    snprintf (filename, BUFSIZ, "%s_sph_%d_%d_%s_%d",
              g->prefix, g->minlevel, g->maxlevel, str, lev);
    cont = p4est_vtk_context_new (g->p4est, filename);
    if (NULL == p4est_vtk_write_header (cont)) {
      P4EST_LERRORF ("Failed to write header for %s\n", filename);
      break;
    }

    /* prepare cell data for output */
    sdata = sc_array_new_count
      (sizeof (double), g->p4est->local_num_quadrants);
    for (lall = 0, tt = g->p4est->first_local_tree;
         tt <= g->p4est->last_local_tree; ++tt) {
      tree = p4est_tree_array_index (g->p4est->trees, tt);
      for (lq = 0; lq < (p4est_locidx_t) tree->quadrants.elem_count; ++lq) {
        /* access number of spheres for this quadrant */
        *(double *) sc_array_index_int (sdata, lall) =
          *(int *) sc_array_index_int (g->lcounts, lall);
        ++lall;
      }
    }
    P4EST_ASSERT (lall == g->p4est->local_num_quadrants);

#if 1
    /* write cell data to file */
    if (NULL == p4est_vtk_write_cell_dataf
        (cont, 1, 1, 1, g->mpiwrap, 1, 0, "spheres", sdata, cont)) {
      P4EST_LERRORF ("Failed to write cell data for %s\n", filename);
      break;
    }
#endif
    sc_array_destroy_null (&sdata);

    /* finish meta information and close files */
    if (p4est_vtk_write_footer (cont)) {
      P4EST_LERRORF ("Failed to write footer for %s\n", filename);
      break;
    }
  }
  while (0);
  if (sdata != NULL) {
    sc_array_destroy_null (&sdata);
  }
}

static void
sphere_offsets (spheres_global_t * g)
{
  int                 p;
  int                 mpiret;
  p4est_gloidx_t      lg, *gval;

  P4EST_ASSERT (g != NULL && g->goffsets != NULL);
  P4EST_ASSERT (g->goffsets->elem_size == sizeof (p4est_gloidx_t));
  P4EST_ASSERT (g->goffsets->elem_count == (size_t) (g->mpisize + 1));
  P4EST_ASSERT (g->sphr->elem_count == (size_t) g->lsph);

  /* obtain globally unique numbers for the spheres */
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
  P4EST_ASSERT
    (lg == *(p4est_gloidx_t *) sc_array_index (g->goffsets, g->mpisize));
  g->gsoff = *(p4est_gloidx_t *) sc_array_index (g->goffsets, g->mpirank);
}

static void
create_forest (spheres_global_t * g)
{
  int                 mpiret;
  double              rmin, rmax, rgeom2;
  double              coef, fact, vmult;
  double              Vexp, Nexp, r;
  double              vdensity, sumrd, gsrd;
  double              tnew, tgenerate;
  char                filename[BUFSIZ];
  p4est_topidx_t      which_tree;
  p4est_locidx_t      ntrel, tin;
  p4est_locidx_t      qunsph;
  p4est_locidx_t      sph_excl, sph_incl;
  p4est_locidx_t      li;
  p4est_qcoord_t      qh;
  p4est_tree_t       *tree;
  p4est_quadrant_t   *q;
  p4est_sphere_t     *sph;

  /* create empty initial forest */
  tnew = sc_MPI_Wtime ();
  g->conn = p4est_connectivity_new_periodic ();
  g->p4est = p4est_new_ext (g->mpicomm, g->conn, 0, g->minlevel, 1,
                            sizeof (qu_data_t), spheres_init_zero, g);
  accumulate_once (g, STATS_ONCE_NEW, tnew);

  /* output initial mesh */
  if (g->write_vtk) {
    snprintf (filename, BUFSIZ, "%s_sph_%d_%d_%s_%d",
              g->prefix, g->minlevel, g->maxlevel, "new", g->minlevel);
    p4est_vtk_write_file (g->p4est, NULL, filename);
  }

  /* minimum and maximum radius determined by target levels */
  tgenerate = sc_MPI_Wtime ();
  rmax = .5 * g->spherelems * sc_intpowf (.5, g->minlevel);
  rmax = SC_MIN (rmax, g->rmax);
  P4EST_ASSERT (0. < rmax);
  rmin = .5 * g->spherelems * sc_intpowf (.5, g->maxlevel);
  rmin = SC_MIN (rmin, g->rmax);
  P4EST_ASSERT (0. < rmin && rmin <= rmax);
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
  g->lcounts = sc_array_new_count (sizeof (int),
                                   g->p4est->local_num_quadrants);
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
      *(int *) sc_array_index (g->lcounts, tree->quadrants_offset + tin) =
        (int) qunsph;

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
          P4EST_INFOF ("Created sphere at " SPHERES_DIM_G " radius %g\n",
                       sph->center[0], sph->center[1]
                       P4EST_ONLY_P8_COMMA (sph->center[2]), r);
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
  sphere_offsets (g);
  P4EST_GLOBAL_PRODUCTIONF
    ("Expected volume %g count %.1f generated %lld\n",
     Vexp, vmult * g->conn->num_trees, (long long)
     *(p4est_gloidx_t *) sc_array_index_int (g->goffsets, g->mpisize));

  /* confirm expected volume */
  mpiret = sc_MPI_Allreduce (&sumrd, &gsrd, 1, sc_MPI_DOUBLE,
                             sc_MPI_SUM, g->mpicomm);
  SC_CHECK_MPI (mpiret);
  P4EST_GLOBAL_PRODUCTIONF ("Total volume ideal %g achieved %g\n",
                            vdensity * g->conn->num_trees, gsrd);

  /* stop sphere generation timer */
  accumulate_once (g, STATS_ONCE_GENERATE, tgenerate);
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
  int                 last_sphere_proc;
  p4est_locidx_t      li;
  p4est_sphere_t     *sph;
  sph_item_t         *item;
  sr_buf_t           *to_proc;

  P4EST_ASSERT (g != NULL && g->p4est == p4est);
  P4EST_ASSERT (0 <= pfirst && pfirst <= plast && plast < g->mpisize);
  P4EST_ASSERT (point != NULL);
  P4EST_ASSERT (g->box.radius ==
                .5 * P4EST_QUADRANT_LEN (quadrant->level) * irootlen);

  /* access sphere under consideration */
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
  last_sphere_proc = *(int *) sc_array_index_int (g->sphere_procs, li);
  if (pfirst != g->last_to_rank) {
    P4EST_ASSERT (g->last_to_rank < pfirst);
    P4EST_ASSERT (last_sphere_proc <= g->last_to_rank);

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
    /* check whether the sphere is a duplicate */
    P4EST_ASSERT (last_sphere_proc <= pfirst);
    if (last_sphere_proc == pfirst) {
#ifdef SPHERES_CHATTY
      P4EST_INFOF ("Duplicate in %d %d: %ld\n", pfirst, plast, (long) li);
#endif
      return 0;
    }

    /* pust to existing send buffer for current remote process */
    to_proc = g->last_to_proc;
    P4EST_ASSERT (to_proc->rank == pfirst);
    P4EST_ASSERT (to_proc == (sr_buf_t *)
                  sc_array_index (g->to_procs, g->to_procs->elem_count - 1));
    P4EST_ASSERT (to_proc->items->elem_count > 0);
    item = (sph_item_t *) sc_array_push (to_proc->items);
    ++*g->last_payload;
  }
  *(int *) sc_array_index_int (g->sphere_procs, li) = pfirst;

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
  qu_data_t          *qud = (qu_data_t *) quadrant->p.user_data;
  p4est_sphere_t     *sph;

  P4EST_ASSERT (g != NULL && g->p4est == p4est);
  P4EST_ASSERT (point != NULL);

  /* access sphere under consideration */
  sph = *(p4est_sphere_t **) point;

  /* if quadrant is small enough already, skip this sphere */
  if (.5 * g->spherelems * P4EST_QUADRANT_LEN (quadrant->level) * irootlen <=
      sph->radius) {
#ifdef SPHERES_CHATTY
    P4EST_INFOF ("No radius match %g %g in %ld\n",
                 g->box.radius, sph->radius, (long) local_num);
#endif
    return 0;
  }

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
    P4EST_INFOF ("No exact match in %ld\n", (long) local_num);
#endif
    return 0;
  }

#ifdef SPHERES_CHATTY
  P4EST_INFOF ("Found in %ld\n", (long) local_num);
#endif
  qud->set_refine = 1;

  /* this return value is ignored */
  return 0;
}

static int
refine_spheres (spheres_global_t * g, int lev)
{
  int                 mpiret;
  int                 q;
  int                 num_from_spheres;
  int                 is_refined;
  double              tpsearch, tsend, tnotify, treceive, twaitsend;
  double              tlsearch, trefine, tpartition, ttransfer;
  sc_array_t         *points;
  sc_array_t         *pi;
  sc_array_t         *lcounts_partitioned;
  sc_array_t         *sphr_partitioned;
  sc_notify_t        *notifyc;
  p4est_locidx_t      li;
  p4est_locidx_t      sri, snum;
#ifdef P4EST_ENABLE_DEBUG
  p4est_gloidx_t      gcur, gnext;
#endif
  p4est_gloidx_t      gsloc, gsglo;
  p4est_gloidx_t      gnq;
  p4est_t            *post;
  p4est_sphere_t    **psph;
  sph_item_t         *item;
  sr_buf_t           *proc;

  P4EST_ASSERT (g->minlevel <= lev && lev < g->maxlevel);

  /*---------------- search partition to find receivers --------------*/

  /* prepare data structures for sending spheres */
  tpsearch = sc_MPI_Wtime ();
  g->last_to_proc = NULL;
  g->last_to_rank = -1;
  g->to_procs = sc_array_new (sizeof (sr_buf_t));
  g->notify = sc_array_new (sizeof (int));
  g->payload = sc_array_new (sizeof (int));
  g->last_payload = NULL;
  g->sphere_procs = sc_array_new_count (sizeof (int), g->lsph);
  sc_array_memset (g->sphere_procs, -1);

  /* search for remote quadrants that receive spheres */
  points = sc_array_new_count (sizeof (p4est_locidx_t), g->lsph);
  for (li = 0; li < g->lsph; ++li) {
    *(p4est_locidx_t *) sc_array_index_int (points, li) = li;
  }
  P4EST_INFOF ("Searching partition for %ld local spheres\n", (long) g->lsph);
  p4est_search_partition (g->p4est, 0, spheres_partition_quadrant,
                          spheres_partition_point, points);
  sc_array_destroy_null (&points);
  accumulate_loop (g, lev, STATS_LOOP_PSEARCH, tpsearch);

  /*------------------------ send the spheres ------------------------*/

  tsend = sc_MPI_Wtime ();
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
  sc_array_destroy_null (&g->sphere_procs);
  accumulate_loop (g, lev, STATS_LOOP_SEND, tsend);

  /*------------------ reverse communication pattern -----------------*/

  tnotify = sc_MPI_Wtime ();
  notifyc = sc_notify_new (g->mpicomm);
  if (!g->notify_alltoall) {
    sc_notify_set_type (notifyc, SC_NOTIFY_NARY);
    sc_notify_nary_set_widths (notifyc, g->ntop, g->nint, g->nbot);
  }
  else {
    sc_notify_set_type (notifyc, SC_NOTIFY_PEX);
  }
  sc_notify_payload (g->notify, NULL, g->payload, NULL, 1, notifyc);
  sc_notify_destroy (notifyc);
  accumulate_loop (g, lev, STATS_LOOP_NOTIFY, tnotify);

  /*---------------------- receive the spheres -----------------------*/

  treceive = sc_MPI_Wtime ();
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
  accumulate_loop (g, lev, STATS_LOOP_RECEIVE, treceive);

  /*---------------- refine based on received spheres ----------------*/

  tlsearch = sc_MPI_Wtime ();
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

  /* report (partially redundant) global sum of spheres */
  gsloc = (p4est_gloidx_t) sri;
  mpiret = sc_MPI_Allreduce (&gsloc, &gsglo, 1, P4EST_MPI_GLOIDX,
                             sc_MPI_SUM, g->mpicomm);
  SC_CHECK_MPI (mpiret);
  P4EST_GLOBAL_PRODUCTIONF ("Sphere redundant sum %lld\n", (long long) gsglo);

  /* search through local elements and set refinement flag */
  P4EST_INFOF ("Searching elements for %ld local spheres\n", (long) sri);
  p4est_search_local (g->p4est, 0, spheres_local_quadrant,
                      spheres_local_point, points);
  sc_array_destroy_null (&points);

  /* the information received is no longer necessary */
  for (q = 0; q < g->num_from_procs; ++q) {
    proc = (sr_buf_t *) sc_array_index_int (g->from_procs, q);
    P4EST_ASSERT (proc->rank != g->mpirank);
    P4EST_ASSERT (proc->items != NULL);
    sc_array_destroy (proc->items);
  }
  sc_array_destroy_null (&g->from_procs);
  accumulate_loop (g, lev, STATS_LOOP_LSEARCH, tlsearch);

  /* perform actual refinement */
  trefine = sc_MPI_Wtime ();
  g->lqindex = g->lqindex_refined = 0;
  g->lsph_offset = 0;
  P4EST_ASSERT ((p4est_locidx_t) g->lcounts->elem_count ==
                g->p4est->local_num_quadrants);
  g->lcounts_refined = sc_array_new_count (sizeof (int),
                                           g->p4est->local_num_quadrants);
  gnq = g->p4est->global_num_quadrants;
  p4est_refine_ext (g->p4est, 0, P4EST_QMAXLEVEL, spheres_refine_callback,
                    spheres_init_zero, spheres_replace_callback);
  P4EST_ASSERT (g->lqindex == (p4est_locidx_t) g->lcounts->elem_count);
  P4EST_ASSERT (g->lqindex_refined ==
                (p4est_locidx_t) g->lcounts_refined->elem_count);
  P4EST_ASSERT ((p4est_locidx_t) g->lcounts_refined->elem_count ==
                g->p4est->local_num_quadrants);
  P4EST_ASSERT (g->lsph_offset == g->lsph);
  P4EST_ASSERT (gnq <= g->p4est->global_num_quadrants);
  is_refined = (gnq < g->p4est->global_num_quadrants);
  sc_array_destroy (g->lcounts);
  g->lcounts = g->lcounts_refined;
  accumulate_loop (g, lev, STATS_LOOP_REFINE, trefine);

  /* output refined forest */
  if (is_refined) {
    spheres_write_vtk (g, "refined", lev + 1);
  }

  /*---------------- complete send and second cleanup ----------------*/

  twaitsend = sc_MPI_Wtime ();
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
  accumulate_loop (g, lev, STATS_LOOP_WAITSEND, twaitsend);

  /*-------------- partition and transfer owned spheres --------------*/

  if (is_refined) {
    /* copy the forest and partition it */
    tpartition = sc_MPI_Wtime ();
    post = p4est_copy (g->p4est, 1);
    (void) p4est_partition_ext (post, 0, NULL);
    accumulate_loop (g, lev, STATS_LOOP_PARTITION, tpartition);

    /* we go through this even if there was no change in partition. */
    ttransfer = sc_MPI_Wtime ();
    lcounts_partitioned =
      sc_array_new_count (sizeof (int), post->local_num_quadrants);
    p4est_transfer_fixed (post->global_first_quadrant,
                          g->p4est->global_first_quadrant,
                          g->mpicomm, SPHERES_TAG_FIXED,
                          lcounts_partitioned->array, g->lcounts->array,
                          sizeof (int));

    /* transfer spheres to their new owners */
    for (snum = 0, li = 0; li < post->local_num_quadrants; ++li) {
      snum += *(int *) sc_array_index_int (lcounts_partitioned, li);
    }
    sphr_partitioned = sc_array_new_count (sizeof (p4est_sphere_t), snum);
    p4est_transfer_items (post->global_first_quadrant,
                          g->p4est->global_first_quadrant,
                          g->mpicomm, SPHERES_TAG_CUSTOM,
                          sphr_partitioned->array,
                          (int *) lcounts_partitioned->array,
                          g->sphr->array, (int *) g->lcounts->array,
                          sizeof (p4est_sphere_t));

    /* reassign partitioned forest and data */
    p4est_destroy (g->p4est);
    g->p4est = post;
    sc_array_destroy (g->lcounts);
    g->lcounts = lcounts_partitioned;
    sc_array_destroy (g->sphr);
    g->sphr = sphr_partitioned;
    g->lsph = (p4est_locidx_t) g->sphr->elem_count;
    sphere_offsets (g);
    accumulate_loop (g, lev, STATS_LOOP_TRANSFER, ttransfer);

    /* output partitioned forest */
    spheres_write_vtk (g, "partitioned", lev + 1);
  }

  return is_refined;
}

static void
destroy_forest (spheres_global_t * g)
{
  sc_array_destroy_null (&g->sphr);
  sc_array_destroy_null (&g->lcounts);
  sc_array_destroy_null (&g->goffsets);

  p4est_destroy (g->p4est);
  p4est_connectivity_destroy (g->conn);
}

static void
run (spheres_global_t * g)
{
  int                 mpiret;
  int                 i;
  int                 lev;
  double              twall;
  char                name[BUFSIZ];
  sc_statinfo_t      *si;

  /* synchronize before timing begins */
  mpiret = sc_MPI_Barrier (g->mpicomm);
  SC_CHECK_MPI (mpiret);
  twall = sc_MPI_Wtime ();

  /* allocate statistics counters */
  g->num_stats =
    STATS_ONCE_LAST + (g->maxlevel - g->minlevel) * STATS_LOOP_LAST;
  g->stats = sc_array_new_count (sizeof (sc_statinfo_t), g->num_stats);
  for (i = 0; i < STATS_ONCE_LAST; ++i) {
    sc_stats_init_ext (once_index (g, (stats_once_t) i), once_names[i],
                       1, sc_stats_group_all, sc_stats_prio_all);
  }
  for (lev = g->minlevel; lev < g->maxlevel; ++lev) {
    for (i = 0; i < STATS_LOOP_LAST; ++i) {
      snprintf (name, BUFSIZ, "%10s_%02d", loop_names[i], lev);
      sc_stats_init_ext (loop_index (g, lev, (stats_loop_t) i), name,
                         1, sc_stats_group_all, sc_stats_prio_all);
    }
  }

  /* create forest, populate with spheres, loop refine and partition */
  create_forest (g);
  for (lev = g->minlevel; lev < g->maxlevel; ++lev) {
    P4EST_GLOBAL_PRODUCTIONF ("Trying refinement from level %d to %d\n",
                              lev, lev + 1);
    if (!refine_spheres (g, lev)) {
      P4EST_GLOBAL_PRODUCTIONF ("No refinement at level %d\n", lev);
      break;
    }
  }

  /* free all memory */
  destroy_forest (g);

  /* report performance statistics */
  accumulate_once (g, STATS_ONCE_WALL, twall);
  si = once_index (g, STATS_ONCE_FIRST);
  sc_stats_compute (sc_MPI_COMM_WORLD, g->num_stats, si);
  sc_stats_print (p4est_package_id, SC_LP_PRODUCTION, g->num_stats, si, 0, 1);
  for (i = 0; i < g->num_stats; ++i) {
    sc_stats_reset (si + i, 1);
  }
  sc_array_destroy (g->stats);
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
  int                 j;
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
  sc_options_add_double (opt, 'r', "rmax", &g->rmax, .5, "Max sphere radius");
  sc_options_add_double (opt, 't', "thickness", &g->thickness, .05,
                         "Relative sphere thickness");
  sc_options_add_double (opt, 'f', "lfraction", &g->lfraction, .3,
                         "Length density of spheres");
  sc_options_add_double (opt, 's', "spherelems", &g->spherelems, 1.,
                         "Min elements per sphere diameter");

  g->ntop = g->nint = P4EST_CHILDREN;
  sc_options_add_int (opt, 'N', "nbottom", &g->nbot, 24,
                      "Notify bottom multiplicator");
  sc_options_add_bool (opt, 'A', "alltoall", &g->notify_alltoall, 0,
                       "Notify alltoall implementation");
  sc_options_add_bool (opt, 'S', "scaling", &g->scaling, 0,
                       "Configure for scaling test");
  sc_options_add_int (opt, 'R', "repetitions", &g->repetitions, 1,
                      "Repeat run multiple times");

  sc_options_add_bool (opt, 'V', "write-vtk", &g->write_vtk, 0,
                       "Output VTK files");
  sc_options_add_string (opt, 'P', "prefix", &g->prefix,
                         "sph" SPHERES_48 ()"res", "Prefix for file output");

  /* set other parameters */
  g->mpiwrap = 16;

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
    P4EST_GLOBAL_PRODUCTIONF ("Dimension is %d\n", P4EST_DIM);
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
    for (j = 0; j < g->repetitions; ++j) {
      P4EST_GLOBAL_PRODUCTIONF ("Spheres run repetition %d of %d\n",
                                j, g->repetitions);
      run (g);
    }
    P4EST_GLOBAL_PRODUCTIONF ("Spheres run %d repetitions done\n",
                              g->repetitions);
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
