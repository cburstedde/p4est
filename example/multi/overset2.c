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
 * Usage: p4est_overset
 *
 * Split a communicator for a background and one or more overset meshes.
 * The background mesh is a standard p4est, each overset mesh is currently
 * represented by a partitioned set of triangles and nodes.
 * We implement overset algorithms required for e. g. windfarm simulation.
 */

#include <sc_notify.h>
#include <sc_options.h>
#ifndef P4_TO_P8
#include <p4est_extended.h>
#include <p4est_search.h>
#include <p4est_vtk.h>
#include <p4est_bits.h>
#include "p4est_multi_overset.h"
#else
#include <p8est_extended.h>
#include <p8est_search.h>
#include <p8est_vtk.h>
#include <p8est_bits.h>
#include "p8est_multi_overset.h"
#endif

typedef struct background
{
  int                 bgminl;
  p4est_connectivity_t *bgconn;
  p4est_t            *bgp4est;
}
background_t;

typedef struct overset
{
  int                 osi;      /* zero-based index of overset mesh */
}
overset_t;

typedef struct overset_global
{
  sc_MPI_Comm         glocomm;
  int                 glosize, glorank;
  int                 num_meshes;       /* background plus overset meshes */
  int                 num_overset;      /* number of overset meshes */
  int                 myrole;   /* 0 for background, index + 1 for overset */
  int                *roffsets; /* offset into sub-partition of meshes */
  int                *rcounts;  /* size of mesh partitions */
  union role {
    background_t        bg;
    overset_t           os;
  }
  r;
  sc_MPI_Comm         rolecomm;         /* mesh sub-communicator */
}
overset_global_t;

static void
overset_init_background (overset_global_t *g)
{
  P4EST_ASSERT (g->myrole == 0);
  background_t       *b = &g->r.bg;

#ifndef P4_TO_P8
  g->r.bg.bgconn = p4est_connectivity_new_unitsquare ();
#else
  g->r.bg.bgconn = p8est_connectivity_new_unitcube ();
#endif
  b->bgp4est = p4est_new_ext (g->rolecomm, b->bgconn, 0,
                              b->bgminl, 1, 0, NULL, g);
}

static void
overset_init_overset (overset_global_t *g, int osi)
{
  P4EST_ASSERT (g->myrole > 0 && osi + 1 == g->myrole);
}

static void
overset_apps_init (overset_global_t *g, sc_MPI_Comm mpicomm)
{
  int                mpiret;
  int                i, ro;

  /* determine parallel configuration */
  g->glocomm = mpicomm;
  mpiret = sc_MPI_Comm_size (g->glocomm, &g->glosize);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (g->glocomm, &g->glorank);
  SC_CHECK_MPI (mpiret);
  if (1 + g->num_overset > g->glosize) {
    g->num_overset = g->glosize - 1;
    P4EST_GLOBAL_PRODUCTIONF
      ("Processes provided %d: reducing num_overset to %d\n",
       g->glosize, g->num_overset);
  }
  g->num_meshes = 1 + g->num_overset;
  P4EST_ASSERT (g->num_meshes >= 1);

  /* compute simulation partition offsets */
  g->roffsets = P4EST_ALLOC (int, g->num_meshes + 1);
  g->rcounts = P4EST_ALLOC (int, g->num_meshes);
  g->myrole = -1;
  for (i = 0; i <= g->num_meshes; ++i) {
    ro = g->roffsets[i] =
      p4est_partition_cut_int (g->glosize, i, g->num_meshes);
    if (i > 0) {
      ro = g->rcounts[i - 1] = ro - g->roffsets[i - 1];
      P4EST_ASSERT (ro > 0);
      P4EST_GLOBAL_INFOF ("Mesh %d ranks %d\n", i - 1, ro);
      if (g->roffsets[i - 1] <= g->glorank &&
          g->glorank < g->roffsets[i]) {
        g->myrole = i - 1;
      }
    }
  }
  P4EST_ASSERT (0 <= g->myrole && g->myrole < g->num_meshes);
  P4EST_LDEBUGF ("My role %d\n", g->myrole);

  /* split communicator for the different meshes */
  mpiret = sc_MPI_Comm_split (g->glocomm, g->myrole, g->glorank,
                              &g->rolecomm);
  SC_CHECK_MPI (mpiret);
  if (g->myrole == 0) {
    overset_init_background (g);
  }
  else {
    overset_init_overset (g, g->myrole - 1);
  }
}

static void
overset_apps_reset (overset_global_t *g)
{
  if (g->myrole == 0) {
    p4est_destroy (g->r.bg.bgp4est);
    p4est_connectivity_destroy (g->r.bg.bgconn);
  }
  sc_MPI_Comm_free (&g->rolecomm);
  P4EST_FREE (g->roffsets);
  P4EST_FREE (g->rcounts);
}

static void
overset_overset (overset_global_t *g)
{
  p4est_multi_overset (g->glorank, g->myrole, g->num_meshes, g->roffsets);
}

int
main (int argc, char **argv)
{
  int                 mpiret;
  int                 first_argc;
  sc_MPI_Comm         mpicomm;
  sc_options_t       *opt;
  overset_global_t    global, *g = &global;

  memset (g, -1, sizeof (overset_global_t));

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  mpicomm = sc_MPI_COMM_WORLD;
  sc_init (mpicomm, 1, 1, NULL, SC_LP_DEFAULT);

  /* process command line arguments */
  opt = sc_options_new (argv[0]);
  sc_options_add_int (opt, 'b', "bg_minlevel", &g->r.bg.bgminl, 0,
                      "Lowest background level");
  sc_options_add_int (opt, 'o', "num_overset", &g->num_overset, 1,
                      "Number of overset meshes");

  first_argc = sc_options_parse (p4est_package_id, SC_LP_DEFAULT,
                                 opt, argc, argv);
  if (first_argc < 0 || first_argc != argc) {
    sc_options_print_usage (p4est_package_id, SC_LP_ERROR, opt, NULL);
    return EXIT_FAILURE;
  }
  sc_options_print_summary (p4est_package_id, SC_LP_ESSENTIAL, opt);

  overset_apps_init (g, mpicomm);

  overset_overset (g);

  overset_apps_reset (g);

  /* clean up application */
  sc_options_destroy (opt);
  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return EXIT_SUCCESS;
}
