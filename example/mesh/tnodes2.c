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
 * Usage: p4est_tnodes <connectivity> <level>
 *        possible connectivities:
 *        o unit      Refinement on the unit square.
 *        o three     Refinement on a forest with three trees.
 *        o moebius   Refinement on a 5-tree Moebius band.
 *        o star      Refinement on a 6-tree star shaped domain.
 *        o periodic  Refinement on the unit square with all-periodic b.c.
 *        o rotwrap   Refinement on the unit square with weird periodic b.c.
 *        o cubed     Refinement on the 2D cubed sphere.
 *        o disk      Refinement on a 5-tree flat disk or square.
 *        o pdisk     Refinement on 5-tree flat disk or square, periodic b.c.
 */

#ifndef P4_TO_P8
#include <p4est_bits.h>
#include <p4est_extended.h>
#include <p4est_tnodes.h>
#include <p4est_vtk.h>
#else
#include <p8est_bits.h>
#include <p8est_extended.h>
#include <p8est_tnodes.h>
#include <p8est_vtk.h>
#endif

typedef enum
{
  P4EST_CONFIG_NULL,
#ifndef P4_TO_P8
  P4EST_CONFIG_UNIT,
  P4EST_CONFIG_THREE,
  P4EST_CONFIG_MOEBIUS,
  P4EST_CONFIG_STAR,
  P4EST_CONFIG_PERIODIC,
  P4EST_CONFIG_ROTWRAP,
  P4EST_CONFIG_CUBED,
  P4EST_CONFIG_DISK,
  P4EST_CONFIG_PDISK
#else
  P8EST_CONFIG_UNIT,
  P8EST_CONFIG_PERIODIC,
  P8EST_CONFIG_ROTWRAP,
  P8EST_CONFIG_TWOCUBES,
  P8EST_CONFIG_TWOWRAP,
  P8EST_CONFIG_ROTCUBES,
  P8EST_CONFIG_SHELL,
  P8EST_CONFIG_SPHERE
#endif
}
simple_config_t;

typedef struct
{
  int                 dummy;
}
user_data_t;

typedef struct
{
  sc_MPI_Comm         mpicomm;
  int                 mpisize;
  int                 mpirank;
}
mpi_context_t;

/* Global variable is not recommended, yes.
 * To clean this up in your own application, please put your context
 * into the p4est->user_pointer field and access it from the callbacks.
 */
static int          refine_level = 0;

static void
tmesh_meta (void)
{
  p4est_tnodes_context_t *econ;

  econ = p4est_tnodes_context_new ();
  p4est_tnodes_context_destroy (econ);
}

static void
init_fn (p4est_t * p4est, p4est_topidx_t which_tree,
         p4est_quadrant_t * quadrant)
{
  user_data_t        *data = (user_data_t *) quadrant->p.user_data;
  data->dummy = -1;
}

static int
refine_uniform (p4est_t * p4est, p4est_topidx_t which_tree,
                p4est_quadrant_t * quadrant)
{
  return (int) quadrant->level < refine_level;
}

static int
refine_normal (p4est_t * p4est, p4est_topidx_t which_tree,
               p4est_quadrant_t * quadrant)
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
#ifndef P4_TO_P8
  if (quadrant->x >= P4EST_QUADRANT_LEN (2)) {
    return 0;
  }
#else
  if (quadrant->z >= P8EST_QUADRANT_LEN (2)) {
    return 0;
  }
#endif

  return 1;
}

static void
tnodes_run (p4est_t * p4est, p4est_ghost_t * ghost,
            int full_style, int with_faces)
{
  p4est_lnodes_t     *ln;
  p4est_tnodes_t     *tm;
#if 0
#ifndef P4_TO_P8
  p4est_tnodes_iter_t *iter;
  p4est_locidx_t      lt;
#endif
#endif

  P4EST_GLOBAL_PRODUCTIONF ("tnodes run %d\n", with_faces);

  P4EST_ASSERT (p4est != NULL);
  P4EST_ASSERT (ghost != NULL);

  ln = p4est_lnodes_new (p4est, ghost, 2);
  tm = p4est_tnodes_new_Q2 (ln, 1, 0);

#if 0
  /* generate triangle mesh */
  tm = p4est_tnodes_new (p4est, ghost, full_style, with_faces
#ifdef P4_TO_P8
                         , 0
#endif
    );

#ifndef P4_TO_P8
  /* iterate through with triangle mesh */
  lt = 0;
  for (iter = p4est_tnodes_iter_new (p4est, tm);
       iter != NULL; p4est_tnodes_iter_next (&iter)) {
    P4EST_ASSERT (lt == iter->triangle);
    ++lt;
  }
  P4EST_ASSERT (lt == tm->global_tcount[p4est->mpirank]);
  P4EST_LDEBUGF ("Just iterated through %ld local triangles\n", (long) lt);
#endif
#endif

  /* free triangle mesh */
  p4est_tnodes_destroy (tm);
}

static void
forest_run (mpi_context_t * mpi,
            p4est_connectivity_t * connectivity, p4est_geometry_t * geom,
            int uniform)
{
  int                 l;
  unsigned            crc;
  char                msg[BUFSIZ];
  p4est_t            *p4est;
  p4est_ghost_t      *ghost;

  P4EST_GLOBAL_PRODUCTIONF ("Forest run uniform %d\n", uniform);

  /* create new coarse p4est from specified connectivity */
  p4est = p4est_new_ext (mpi->mpicomm, connectivity, 0, 0, 1,
                         sizeof (user_data_t), init_fn, NULL);
  p4est_vtk_write_file (p4est, geom, P4EST_STRING "_tnodes_new");

  /* non-recursive refinement loop */
  for (l = 1; l <= refine_level; ++l) {
    /* refine */
    p4est_refine (p4est, 0, uniform ? refine_uniform : refine_normal,
                  init_fn);
    snprintf (msg, BUFSIZ, P4EST_STRING "_tnodes_refined_%02d", l);
    p4est_vtk_write_file (p4est, geom, msg);

    if (!uniform) {
      /* balance */
      p4est_balance (p4est, P4EST_CONNECT_FULL, init_fn);
      snprintf (msg, BUFSIZ, P4EST_STRING "_tnodes_balanced_%02d", l);
      p4est_vtk_write_file (p4est, geom, msg);
    }

    /* partition */
    p4est_partition (p4est, 0, NULL);
    snprintf (msg, BUFSIZ, P4EST_STRING "_tnodes_partitioned_%02d", l);
    p4est_vtk_write_file (p4est, geom, msg);
  }
  crc = p4est_checksum (p4est);

  /* print and verify forest checksum */
  P4EST_GLOBAL_STATISTICSF ("Forest %s checksum 0x%08x\n",
                            uniform ? "uniform" : "adapted", crc);

  /* create ghost layer and triangle meshes */
  ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FULL);
  tnodes_run (p4est, ghost, 0, 0);
#if 0
  tnodes_run (p4est, ghost, 1, 0);
  tnodes_run (p4est, ghost, 0, 1);
#endif
  p4est_ghost_destroy (ghost);
#if 0
  tnodes_run (p4est, NULL, 1, 1);
#endif

  /* destroy the p4est structure */
  p4est_destroy (p4est);
}

int
main (int argc, char **argv)
{
  int                 mpiret;
  int                 wrongusage;
  const char         *usage;
  mpi_context_t       mpi_context, *mpi = &mpi_context;
  p4est_connectivity_t *connectivity;
  p4est_geometry_t   *geometry;
  simple_config_t     config;

  /* initialize MPI and p4est internals */
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  mpi->mpicomm = sc_MPI_COMM_WORLD;
  mpiret = sc_MPI_Comm_size (mpi->mpicomm, &mpi->mpisize);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (mpi->mpicomm, &mpi->mpirank);
  SC_CHECK_MPI (mpiret);

  sc_init (mpi->mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  /* process command line arguments */
  usage =
    "Arguments: <connectivity> <level>\n   The connectivity can be any of\n"
#ifndef P4_TO_P8
    "      unit|three|moebius|star|periodic|rotwrap|cubed|disk\n"
#else
    "      unit|periodic|rotwrap|twocubes|twowrap|rotcubes|shell|sphere\n"
#endif
    "   Level controls the maximum depth of refinement\n";
  wrongusage = 0;
  config = P4EST_CONFIG_NULL;
  if (!wrongusage && argc != 3) {
    wrongusage = 1;
  }
  if (!wrongusage) {
    if (!strcmp (argv[1], "unit")) {
#ifndef P4_TO_P8
      config = P4EST_CONFIG_UNIT;
#else
      config = P8EST_CONFIG_UNIT;
#endif
    }
#ifndef P4_TO_P8
    else if (!strcmp (argv[1], "three")) {
      config = P4EST_CONFIG_THREE;
    }
    else if (!strcmp (argv[1], "moebius")) {
      config = P4EST_CONFIG_MOEBIUS;
    }
    else if (!strcmp (argv[1], "star")) {
      config = P4EST_CONFIG_STAR;
    }
    else if (!strcmp (argv[1], "periodic")) {
      config = P4EST_CONFIG_PERIODIC;
    }
    else if (!strcmp (argv[1], "rotwrap")) {
      config = P4EST_CONFIG_ROTWRAP;
    }
    else if (!strcmp (argv[1], "cubed")) {
      config = P4EST_CONFIG_CUBED;
    }
    else if (!strcmp (argv[1], "disk")) {
      config = P4EST_CONFIG_DISK;
    }
    else if (!strcmp (argv[1], "pdisk")) {
      config = P4EST_CONFIG_PDISK;
    }
#else
    else if (!strcmp (argv[1], "periodic")) {
      config = P8EST_CONFIG_PERIODIC;
    }
    else if (!strcmp (argv[1], "rotwrap")) {
      config = P8EST_CONFIG_ROTWRAP;
    }
    else if (!strcmp (argv[1], "twocubes")) {
      config = P8EST_CONFIG_TWOCUBES;
    }
    else if (!strcmp (argv[1], "twowrap")) {
      config = P8EST_CONFIG_TWOWRAP;
    }
    else if (!strcmp (argv[1], "rotcubes")) {
      config = P8EST_CONFIG_ROTCUBES;
    }
    else if (!strcmp (argv[1], "shell")) {
      config = P8EST_CONFIG_SHELL;
    }
    else if (!strcmp (argv[1], "sphere")) {
      config = P8EST_CONFIG_SPHERE;
    }
#endif
    else {
      wrongusage = 1;
      P4EST_GLOBAL_LERROR ("Unknown connectivity\n");
    }
  }
  if (!wrongusage) {
    refine_level = atoi (argv[2]);
    if (refine_level < 0 || refine_level > P4EST_QMAXLEVEL) {
      wrongusage = 1;
      P4EST_GLOBAL_LERROR ("Refinement level out of range\n");
    }
  }
  if (wrongusage) {
    P4EST_GLOBAL_LERROR (usage);
    sc_abort_collective ("Usage error");
  }

  /* create connectivity and forest structures */
  connectivity = NULL;
  geometry = NULL;
  if (0) {
  }
#ifndef P4_TO_P8
  else if (config == P4EST_CONFIG_THREE) {
    connectivity = p4est_connectivity_new_corner ();
  }
  else if (config == P4EST_CONFIG_MOEBIUS) {
    connectivity = p4est_connectivity_new_moebius ();
  }
  else if (config == P4EST_CONFIG_STAR) {
    connectivity = p4est_connectivity_new_star ();
  }
  else if (config == P4EST_CONFIG_PERIODIC) {
    connectivity = p4est_connectivity_new_periodic ();
  }
  else if (config == P4EST_CONFIG_ROTWRAP) {
    connectivity = p4est_connectivity_new_rotwrap ();
  }
  else if (config == P4EST_CONFIG_CUBED) {
    connectivity = p4est_connectivity_new_cubed ();
  }
  else if (config == P4EST_CONFIG_DISK) {
    connectivity = p4est_connectivity_new_disk (0, 0);
  }
  else if (config == P4EST_CONFIG_PDISK) {
    connectivity = p4est_connectivity_new_disk (1, 1);
  }
#else
  else if (config == P8EST_CONFIG_PERIODIC) {
    connectivity = p8est_connectivity_new_periodic ();
  }
  else if (config == P8EST_CONFIG_ROTWRAP) {
    connectivity = p8est_connectivity_new_rotwrap ();
  }
  else if (config == P8EST_CONFIG_TWOCUBES) {
    connectivity = p8est_connectivity_new_twocubes ();
  }
  else if (config == P8EST_CONFIG_TWOWRAP) {
    connectivity = p8est_connectivity_new_twowrap ();
  }
  else if (config == P8EST_CONFIG_ROTCUBES) {
    connectivity = p8est_connectivity_new_rotcubes ();
  }
  else if (config == P8EST_CONFIG_SHELL) {
    connectivity = p8est_connectivity_new_shell ();
    geometry = p8est_geometry_new_shell (connectivity, 1., 1.5);
  }
  else if (config == P8EST_CONFIG_SPHERE) {
    connectivity = p8est_connectivity_new_sphere ();
    geometry = p8est_geometry_new_sphere (connectivity, 1., .6, .3);
  }
#endif
  else {
#ifndef P4_TO_P8
    connectivity = p4est_connectivity_new_unitsquare ();
#else
    connectivity = p8est_connectivity_new_unitcube ();
#endif
  }

  /* prepare simplex mesh metadata */
#if 0
  tmesh_meta ();
#endif

  /* run mesh tests */
  forest_run (mpi,              /* mpi context */
              connectivity,     /* p4est connectivity */
              geometry,         /* used for VTK output */
              1);               /* uniform refinement? */
#if 0
  forest_run (mpi, connectivity, geometry, 0);
#endif

  /* clean up and exit */
  if (geometry != NULL) {
    p4est_geometry_destroy (geometry);
  }
  p4est_connectivity_destroy (connectivity);
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
