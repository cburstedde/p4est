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

/** \file circle.c
 * This example is used to generate the mesh of a circle consisting of six trees,
 * with edges connecting the trees to each other.
 * @image html /images/circle_nodes.png  width=30%
 * Usage: <level> \n
 * The images below show the mesh of the circle with a level of refinement 5 after each step. \n
 * Creation:
 * @image html /images/circle_new.png  width=30%
 * Refinement and Coarsening:
 * @image html /images/circle_refined.png  width=30%
 * Balance:
 * @image html /images/circle_balanced.png width=30%
 * Partition:
 * @image html /images/circle_partition.png  width=30%
 */

#include <p4est_bits.h>
#include <p4est_extended.h>
#include <p4est_vtk.h>

typedef enum
{
  P4EST_CONFIG_NULL,
  P4EST_CONFIG_CIRCLE,
}
simple_config_t;

typedef struct
{
  simple_config_t     config;
  int                 mpisize;
  int                 level;
  unsigned            checksum;
}
simple_regression_t;

typedef struct
{
  p4est_topidx_t      a;
}
user_data_t;

typedef struct
{
  sc_MPI_Comm         mpicomm;
  int                 mpisize;
  int                 mpirank;
}
mpi_context_t;

static int          refine_level = 0;

/* *INDENT-OFF* */
static const simple_regression_t regression[] =
  {{ P4EST_CONFIG_CIRCLE, 3, 6, 0x98ab6cb2U }};
/* *INDENT-ON* */

static void
init_fn (p4est_t * p4est, p4est_topidx_t which_tree,
         p4est_quadrant_t * quadrant)
{
  user_data_t        *data = (user_data_t *) quadrant->p.user_data;

  data->a = which_tree;
}

static int
refine_normal_fn (p4est_t * p4est, p4est_topidx_t which_tree,
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

p4est_connectivity_t *
p4est_connectivity_new_circle (void)
{
  const p4est_topidx_t num_vertices = 12;
  const p4est_topidx_t num_trees = 6;
  const p4est_topidx_t num_ctt = 0;
  const double        vertices[12 * 3] = {
    /* inner hexagon */
    0.0, 1.0, 0.0,
    0.866025404, 0.5, 0.0,
    0.866025404, -0.5, 0.0,
    1.2246468e-16, -1.0, 0.0,
    -0.866025404, -0.5, 0.0,
    -0.866025404, 0.5, 0.0,
    /* outer hexagon */
    0.0, 2.0, 0.0,
    1.73205081, 1.0, 0.0,
    1.73205081, -1.0, 0.0,
    2.4492936e-16, -2.0, 0.0,
    -1.73205081, -1.0, 0.0,
    -1.73205081, 1.0, 0.0,
  };
  const p4est_topidx_t tree_to_vertex[6 * 4] = {
    7, 6, 1, 0,
    11, 5, 6, 0,
    5, 11, 4, 10,
    9, 3, 10, 4,
    2, 3, 8, 9,
    8, 7, 2, 1,
  };
  const p4est_topidx_t tree_to_tree[6 * 4] = {
    5, 1, 0, 0,
    1, 1, 2, 0,
    2, 2, 1, 3,
    3, 3, 4, 2,
    5, 3, 4, 4,
    4, 0, 5, 5,
  };
  const int8_t        tree_to_face[6 * 4] = {
    1, 3, 2, 3,
    0, 1, 6, 1,
    0, 1, 6, 7,
    0, 1, 5, 7,
    4, 6, 2, 3,
    4, 0, 2, 3,
  };

  return p4est_connectivity_new_copy (num_vertices, num_trees, 0,
                                      vertices, tree_to_vertex,
                                      tree_to_tree, tree_to_face,
                                      NULL, &num_ctt, NULL, NULL);
}

int
main (int argc, char **argv)
{
  int                 mpiret;
  int                 wrongusage;
  unsigned            crc;
  const char         *usage;
  mpi_context_t       mpi_context, *mpi = &mpi_context;
  p4est_t            *p4est;
  p4est_connectivity_t *connectivity;
  p4est_geometry_t   *geom;
  p4est_refine_t      refine_fn;
  p4est_coarsen_t     coarsen_fn;
  simple_config_t     config;
  const simple_regression_t *r;

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
  usage = "Arguments: <level>\n";
  wrongusage = 0;
  if (!wrongusage && argc != 2) {
    wrongusage = 1;
  }
  config = P4EST_CONFIG_CIRCLE;
  if (wrongusage) {
    P4EST_GLOBAL_LERROR (usage);
    sc_abort_collective ("Usage error");
  }
  /* assign variables based on configuration */
  refine_level = atoi (argv[1]);
  refine_fn = refine_normal_fn;
  coarsen_fn = NULL;
  /* create connectivity and forest structures */
  geom = NULL;
  if (config == P4EST_CONFIG_CIRCLE) {
    connectivity = p4est_connectivity_new_circle ();
  }
  p4est = p4est_new_ext (mpi->mpicomm, connectivity, 15, 0, 0,
                         sizeof (user_data_t), init_fn, geom);
  p4est_vtk_write_file (p4est, geom, "circle_new");

  /* refinement and coarsening */
  p4est_refine (p4est, 1, refine_fn, init_fn);
  if (coarsen_fn != NULL) {
    p4est_coarsen (p4est, 1, coarsen_fn, init_fn);
  }
  p4est_vtk_write_file (p4est, geom, "circle_refined");

  /* balance */
  p4est_balance (p4est, P4EST_CONNECT_FULL, init_fn);
  p4est_vtk_write_file (p4est, geom, "circle_balanced");
  crc = p4est_checksum (p4est);

  /* partition */
  p4est_partition (p4est, 0, NULL);
  p4est_vtk_write_file (p4est, geom, "circle_partition");

#ifdef P4EST_ENABLE_DEBUG
  /* rebalance should not change checksum */
  p4est_balance (p4est, P4EST_CONNECT_FULL, init_fn);
  P4EST_ASSERT (p4est_checksum (p4est) == crc);
#endif

  /* print and verify forest checksum */
  P4EST_GLOBAL_STATISTICSF ("Tree checksum 0x%08x\n", crc);
  if (mpi->mpirank == 0) {
    for (r = regression; r->config != P4EST_CONFIG_NULL; ++r) {
      if (r->config != config || r->mpisize != mpi->mpisize
          || r->level != refine_level)
        continue;
      SC_CHECK_ABORT (crc == r->checksum, "Checksum mismatch");
      P4EST_GLOBAL_INFO ("Checksum regression OK\n");
      break;
    }
  }

  /* destroy the p4est and its connectivity structure */
  p4est_destroy (p4est);
  if (geom != NULL) {
    p4est_geometry_destroy (geom);
  }
  p4est_connectivity_destroy (connectivity);

  /* clean up and exit */
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
