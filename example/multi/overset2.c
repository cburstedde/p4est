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
#include <p4est_vtk.h>
#include <p4est_bits.h>
#include "p4est_multi_overset.h"
#else
#include <p8est_extended.h>
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

#ifdef P4_TO_P8
#define NUM_TRIELEM_CORNERS 4
#define NUM_BLOCK_TRIELEMS 6
static int          block_to_tri[24] = {
  0, 1, 3, 7,
  0, 1, 5, 7,
  0, 2, 3, 7,
  0, 2, 6, 7,
  0, 4, 5, 7,
  0, 4, 6, 7
};
#else
#define NUM_TRIELEM_CORNERS 3
#define NUM_BLOCK_TRIELEMS 2
static int          block_to_tri[6] = {
  0, 1, 3,
  0, 2, 3
};
#endif

typedef struct trielem
{
  /* triangle corners (2D) or tetrahedra corners (3D) */
  double              corners[NUM_TRIELEM_CORNERS * 3];
}
trielem_t;

typedef struct overset
{
  int                 osi;      /* zero-based index of overset mesh */
  sc_array_t         *elements;
}
overset_t;

typedef struct overset_global
{
  sc_MPI_Comm         glocomm;
  int                 glosize, glorank;
  int                 num_meshes;       /* background plus overset meshes */
  int                 num_overset;      /* number of overset meshes */
  int                 myrole;   /* 0 for background, index + 1 for overset */
  int                 ishead;   /* True if first rank of a mesh partition */
  int                *roffsets; /* offset into sub-partition of meshes */
  int                *rcounts;  /* size of mesh partitions */
  double              overset_scaling;  /* scale the non-shifting dimension(s) */
  int                 overset_gridconst;        /* block count of non-shifting dimension(s) */
  union role {
    background_t        bg;
    overset_t           os;
  }
  r;
  sc_MPI_Comm         rolecomm;         /* mesh sub-communicator */
  sc_MPI_Comm         headcomm;         /* sub-communicator for heads */
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
  P4EST_ASSERT (0 <= g->overset_scaling && g->overset_scaling <= 1);

  int                 i, j, tind, cind, btot;
  int                 yblocks, zblocks, num_blocks;
  double              xwidth;
  double              xoffset;
  double              yzblockwidth;
  double              yzboundary_distance;
  double              yoffset, zoffset;
  double              block_vertices[P4EST_CHILDREN][3];
  trielem_t          *tri;

  /* Create a simple prism mesh.
   * We create overset meshes that are evenly spaced in the first dimension.
   * They do not overlap each other, but they are contained in the background
   * mesh.
   * Every mesh is fit into a unifrom grid with overset_gridconst blocks in all
   * of the remaining dimension(s). Every block in the resulting grid is
   * divided into two triangle (in 2D) or six tetrahedra (in 3D) elements.
   * Furthermore, there is a scaling parameter for the remaining dimensions.
   * Hereby, the overset mesh is centered with respect to the remaining
   * dimension(s). */

  /* calculate shift in first dimension */
  xwidth = 1. / (2. * ((double) g->num_meshes - 1.) + 1.);
  xoffset = (2 * osi + 1) * xwidth;

  /* calculate offset in remaing dimensions */
  yzblockwidth = g->overset_scaling / g->overset_gridconst;
  yzboundary_distance = (1. - g->overset_scaling) / 2.;

  /* calculate vertices of a single block with lower left corner in (0,0,0) */
  for (i = 0; i < P4EST_CHILDREN; ++i) {
    block_vertices[i][0] = (i % 2 != 0) ? xwidth : 0.;
    block_vertices[i][1] = (i & 0x2) ? yzblockwidth : 0.;
    block_vertices[i][2] = (i & 0x4) ? yzblockwidth : 0.;
  }

  /* allocate space for the triangular elements in all blocks */
  yblocks = g->overset_gridconst;
#ifdef P4_TO_P8
  zblocks = g->overset_gridconst;
#else
  zblocks = 1;
#endif
  num_blocks = yblocks * zblocks;
  g->r.os.elements =
    sc_array_new_count (sizeof (trielem_t), NUM_BLOCK_TRIELEMS * num_blocks);

  /* iterate over all blocks and create the corresponding elements */
  for (i = 0; i < yblocks; i++) {
    yoffset = yzboundary_distance + i * yzblockwidth;
    for (j = 0; j < zblocks; j++) {
#ifdef P4_TO_P8
      zoffset = yzboundary_distance + j * yzblockwidth;
#else
      zoffset = 0.;
#endif

      for (tind = 0; tind < NUM_BLOCK_TRIELEMS; tind++) {
        tri =
          (trielem_t *) sc_array_index_int (g->r.os.elements,
                                            (i * zblocks +
                                             j) * NUM_BLOCK_TRIELEMS + tind);
        for (cind = 0; cind < NUM_TRIELEM_CORNERS; cind++) {
          /* assign element corners according to lookup table */
          btot = block_to_tri[tind * NUM_TRIELEM_CORNERS + cind];
          tri->corners[3 * cind + 0] = xoffset + block_vertices[btot][0];
          tri->corners[3 * cind + 1] = yoffset + block_vertices[btot][1];
          tri->corners[3 * cind + 2] = zoffset + block_vertices[btot][2];
        }
      }
    }
  }
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
  g->ishead = g->roffsets[g->myrole] == g->glorank;
  P4EST_LDEBUGF ("My role %d head %d\n", g->myrole, g->ishead);

  /* split communicator for the different meshes */
  mpiret = sc_MPI_Comm_split (g->glocomm, g->myrole, g->glorank,
                              &g->rolecomm);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_split (g->glocomm, g->ishead ? 0 : sc_MPI_UNDEFINED,
                              g->glorank, &g->headcomm);
  SC_CHECK_MPI (mpiret);
  P4EST_ASSERT (g->ishead || g->headcomm == sc_MPI_COMM_NULL);

  /* initialize the respective meshes */
  if (g->myrole == 0) {
    overset_init_background (g);
  }
  else {
    overset_init_overset (g, g->myrole - 1);
  }
}

static void
overset_create_query_points (overset_global_t *g, sc_array_t *query_points)
{
  int                 i, cind;
  int                 num_trielems;
  double             *xyzv;
  trielem_t          *tri;
  double              yboundary_distance;
  double              xwidth, yzwidth;
  double              trielem_volume;

  P4EST_ASSERT (g != NULL);
  P4EST_ASSERT (query_points != NULL);
  P4EST_ASSERT (query_points->elem_size == sizeof (double) * 4);
  P4EST_ASSERT (g->myrole != 0);

  /* allocate space for all query points */
  num_trielems = g->r.os.elements->elem_count;
  sc_array_resize (query_points, NUM_TRIELEM_CORNERS * num_trielems);

  /* compute relevant information about the triangle/tetrahedra mesh */
  yboundary_distance = (1. - g->overset_scaling) / 2.;
  xwidth = 1. / (2. * ((double) g->num_meshes - 1.) + 1.);
  yzwidth = g->overset_scaling / g->overset_gridconst;
#ifdef P4_TO_P8
  trielem_volume = xwidth * yzwidth * yzwidth / 6.;
#else
  trielem_volume = xwidth * yzwidth / 2.;
#endif
  for (i = 0; i < (int) num_trielems; ++i) {
    tri = (trielem_t *) sc_array_index_int (g->r.os.elements, i);
    for (cind = 0; cind < NUM_TRIELEM_CORNERS; cind++) {
      xyzv =
        (double *) sc_array_index_int (query_points,
                                       NUM_TRIELEM_CORNERS * i + cind);
      xyzv[0] = tri->corners[3 * cind + 0];
      xyzv[1] = tri->corners[3 * cind + 1];
      xyzv[2] = tri->corners[3 * cind + 2];

      if (xyzv[1] <= yboundary_distance + SC_1000_EPS) {
        /* set all points with minimal y-coordinate to wall boundary */
        xyzv[3] = -2.;
      }
      else if (xyzv[1] >= 1 - yboundary_distance - SC_1000_EPS) {
        /* set all points with maximal y-coordinate to overset boundary */
        xyzv[3] = -1.;
      }
      else {
        xyzv[3] = trielem_volume;
      }
    }
  }
}

static void
overset_apps_reset (overset_global_t *g)
{
  if (g->myrole == 0) {
    P4EST_ASSERT (g->r.bg.bgp4est != NULL);
    p4est_destroy (g->r.bg.bgp4est);
    p4est_connectivity_destroy (g->r.bg.bgconn);
  }
  else {
    P4EST_ASSERT (g->r.os.elements != NULL);
    sc_array_destroy (g->r.os.elements);
  }
  if (g->ishead) {
    sc_MPI_Comm_free (&g->headcomm);
  }
  sc_MPI_Comm_free (&g->rolecomm);
  P4EST_FREE (g->roffsets);
  P4EST_FREE (g->rcounts);
}

static void
overset_overset (overset_global_t *g)
{
  sc_array_t         *qpoints = NULL;
  p4est_t            *bgp4est = NULL;

  if (g->myrole == 0) {
    bgp4est = g->r.bg.bgp4est;
  }
  else {
    qpoints = sc_array_new_count (4 * sizeof (double), 0);
    overset_create_query_points (g, qpoints);
  }

  p4est_multi_overset (g->glocomm, g->headcomm, g->rolecomm,
                       g->myrole, g->num_meshes, g->roffsets,
                       bgp4est, qpoints, NULL, NULL, NULL, NULL, NULL);

  if (g->myrole > 0) {
    sc_array_destroy (qpoints);
  }
}

int
main (int argc, char **argv)
{
  int                 mpiret;
  int                 first_argc;
  sc_MPI_Comm         mpicomm;
  sc_options_t       *opt;
  overset_global_t global, *g = &global;

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
  sc_options_add_double (opt, 's', "scale_overset", &g->overset_scaling, 0.5,
                         "Scaling factor for overset meshes in [0,1]"
                         "for non-shifting dimension(s)");
  sc_options_add_int (opt, 'g', "gridconst_overset", &g->overset_gridconst, 1,
                      "Number of blocks of overset meshes"
                      "in the non-shifting dimension(s)");

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
