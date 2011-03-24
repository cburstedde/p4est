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
 * Usage: p4est_mesh <configuration> <level>
 *        possible configurations:
 *        o unit      Refinement on the unit square.
 *        o three     Refinement on a forest with three trees.
 *        o moebius   Refinement on a 5-tree Moebius band.
 *        o star      Refinement on a 6-tree star shaped domain.
 *        o periodic  Refinement on the unit square with all-periodic b.c.
 *        o rotwrap   Refinement on the unit square with weird periodic b.c.
 */

#ifndef P4_TO_P8
#include <p4est_bits.h>
#include <p4est_mesh.h>
#include <p4est_vtk.h>
#else
#include <p8est_bits.h>
#include <p8est_mesh.h>
#include <p8est_vtk.h>
#endif

typedef enum
{
  P4EST_CONFIG_NULL,
  P4EST_CONFIG_UNIT,
  P4EST_CONFIG_THREE,
  P4EST_CONFIG_MOEBIUS,
  P4EST_CONFIG_STAR,
  P4EST_CONFIG_PERIODIC,
  P4EST_CONFIG_ROTWRAP,
}
simple_config_t;

typedef struct
{
  p4est_topidx_t      a;
}
user_data_t;

typedef struct
{
  MPI_Comm            mpicomm;
  int                 mpisize;
  int                 mpirank;
}
mpi_context_t;

static int          refine_level = 0;

static void
init_fn (p4est_t * p4est, p4est_topidx_t which_tree,
         p4est_quadrant_t * quadrant)
{
  user_data_t        *data = (user_data_t *) quadrant->p.user_data;

  data->a = which_tree;
}

static int
refine_uniform (p4est_t * p4est, p4est_topidx_t which_tree,
                p4est_quadrant_t * quadrant)
{
  return (int) quadrant->level < refine_level;
}

#if 0

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
  if (quadrant->x >= P4EST_QUADRANT_LEN (2)) {
    return 0;
  }

  return 1;
}

#endif

static void
test_mesh_uniform (p4est_t * p4est, p4est_mesh_t * mesh)
{
  int                 f, nf;
  p4est_locidx_t      K, kl;
  p4est_locidx_t      ql, QpG;

  K = mesh->local_num_quadrants;
  P4EST_ASSERT (K == p4est->local_num_quadrants);
  QpG = mesh->local_num_quadrants + mesh->ghost_num_quadrants;

  for (kl = 0; kl < K; ++kl) {
    for (f = 0; f < P4EST_FACES; ++f) {
      ql = mesh->quad_to_quad[P4EST_FACES * kl + f];
      SC_CHECK_ABORTF (0 <= ql && ql < QpG,
                       "quad %d face %d neighbor %d mismatch", kl, f, ql);
      nf = mesh->quad_to_face[P4EST_FACES * kl + f];
      SC_CHECK_ABORTF (0 <= nf && nf < P4EST_HALF * P4EST_FACES,
                       "quad %d face %d code %d mismatch", kl, f, nf);
    }
  }
}

static void
mesh_uniform (mpi_context_t * mpi, p4est_connectivity_t * connectivity)
{
  unsigned            crc;
  p4est_t            *p4est;
  p4est_ghost_t      *ghost;
  p4est_mesh_t       *mesh;

  p4est = p4est_new (mpi->mpicomm, connectivity,
                     sizeof (user_data_t), init_fn, NULL);
  p4est_vtk_write_file (p4est, NULL, "mesh2_uniform_new");

  /* refinement */
  p4est_refine (p4est, 1, refine_uniform, init_fn);
  p4est_vtk_write_file (p4est, NULL, "mesh2_uniform_refined");

  /* balance */
  p4est_balance (p4est, P4EST_BALANCE_FULL, init_fn);
  p4est_vtk_write_file (p4est, NULL, "mesh2_uniform_balanced");

  /* partition */
  p4est_partition (p4est, NULL);
  p4est_vtk_write_file (p4est, NULL, "mesh2_uniform_partition");
  crc = p4est_checksum (p4est);

  /* print and verify forest checksum */
  P4EST_GLOBAL_STATISTICSF ("Tree checksum 0x%08x\n", crc);

  /* create ghost layer and mesh */
  ghost = p4est_ghost_new (p4est, P4EST_BALANCE_FULL);
  mesh = p4est_mesh_new (p4est, ghost, P4EST_BALANCE_FULL);
  test_mesh_uniform (p4est, mesh);

  /* destroy ghost layer and mesh */
  p4est_mesh_destroy (mesh);
  p4est_ghost_destroy (ghost);

  /* destroy the p4est and its connectivity structure */
  p4est_destroy (p4est);
  p4est_connectivity_destroy (connectivity);
}

static void
mesh_adapted (mpi_context_t * mpi, p4est_connectivity_t * connectivity)
{
}

int
main (int argc, char **argv)
{
  int                 mpiret;
  int                 wrongusage;
  const char         *usage;
  mpi_context_t       mpi_context, *mpi = &mpi_context;
  p4est_connectivity_t *connectivity;
  simple_config_t     config;

  /* initialize MPI and p4est internals */
  mpiret = MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  mpi->mpicomm = MPI_COMM_WORLD;
  mpiret = MPI_Comm_size (mpi->mpicomm, &mpi->mpisize);
  SC_CHECK_MPI (mpiret);
  mpiret = MPI_Comm_rank (mpi->mpicomm, &mpi->mpirank);
  SC_CHECK_MPI (mpiret);

  sc_init (mpi->mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  /* process command line arguments */
  usage =
    "Arguments: <configuration> <level>\n"
    "   Configuration can be any of\n"
    "      unit|three|moebius|star|periodic|rotwrap\n"
    "   Level controls the maximum depth of refinement\n";
  wrongusage = 0;
  config = P4EST_CONFIG_NULL;
  if (!wrongusage && argc != 3) {
    wrongusage = 1;
  }
  if (!wrongusage) {
    if (!strcmp (argv[1], "unit")) {
      config = P4EST_CONFIG_UNIT;
    }
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
    else {
      wrongusage = 1;
    }
  }
  if (wrongusage) {
    P4EST_GLOBAL_LERROR (usage);
    sc_abort_collective ("Usage error");
  }

  /* assign variables based on configuration */
  refine_level = atoi (argv[2]);

  /* create connectivity and forest structures */
  if (config == P4EST_CONFIG_THREE) {
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
  else {
    connectivity = p4est_connectivity_new_unitsquare ();
  }

  /* run mesh tests */
  mesh_uniform (mpi, connectivity);
  mesh_adapted (mpi, connectivity);

  /* clean up and exit */
  sc_finalize ();

  mpiret = MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
