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
 * Usage: p4est_simple <configuration> <level>
 *        possible configurations:
 *        o unit      Refinement on the unit square.
 *        o three     Refinement on a forest with three trees.
 *        o evil      Check second round of refinement with np=5 level=7
 *        o evil3     Check second round of refinement on three trees
 *        o pillow    Refinement on a 2-tree pillow-shaped domain.
 *        o moebius   Refinement on a 5-tree Moebius band.
 *        o star      Refinement on a 6-tree star shaped domain.
 *        o cubed     Refinement on a 6-tree cubed sphere surface.
 *        o disk      Refinement on a 5-tree spherical disk.
 *        o periodic  Refinement on the unit square with all-periodic b.c.
 *        o rotwrap   Refinement on the unit square with weird periodic b.c.
 */

#include <p4est_bits.h>
#include <p4est_extended.h>
#include <p4est_vtk.h>

typedef enum
{
  P4EST_CONFIG_NULL,
  P4EST_CONFIG_UNIT,
  P4EST_CONFIG_THREE,
  P4EST_CONFIG_EVIL,
  P4EST_CONFIG_EVIL3,
  P4EST_CONFIG_PILLOW,
  P4EST_CONFIG_MOEBIUS,
  P4EST_CONFIG_STAR,
  P4EST_CONFIG_CUBED,
  P4EST_CONFIG_DISK,
  P4EST_CONFIG_PERIODIC,
  P4EST_CONFIG_ROTWRAP
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
{{ P4EST_CONFIG_THREE, 1, 7, 0xa8d85863U },
 { P4EST_CONFIG_THREE, 2, 7, 0xa8d85863U },
 { P4EST_CONFIG_THREE, 3, 7, 0xa8d85863U },
 { P4EST_CONFIG_THREE, 4, 7, 0x20fb58edU },
 { P4EST_CONFIG_MOEBIUS, 1, 6, 0x98ab6cb2U },
 { P4EST_CONFIG_MOEBIUS, 3, 6, 0x98ab6cb2U },
 { P4EST_CONFIG_MOEBIUS, 5, 6, 0x98ab6cb2U },
 { P4EST_CONFIG_MOEBIUS, 6, 6, 0x6d2d6d6cU },
 { P4EST_CONFIG_STAR, 5, 6, 0x38d3736fU },
 { P4EST_CONFIG_STAR, 5, 7, 0xfb97aadfU },
 { P4EST_CONFIG_CUBED, 4, 3, 0x85581649U },
 { P4EST_CONFIG_CUBED, 5, 5, 0x64a1d105U },
 { P4EST_CONFIG_DISK, 5, 4, 0x4995411dU },
 { P4EST_CONFIG_DISK, 2, 6, 0x3f758706U },
 { P4EST_CONFIG_ROTWRAP, 1, 6, 0x9dd600c5U },
 { P4EST_CONFIG_ROTWRAP, 3, 6, 0x9dd600c5U },
 { P4EST_CONFIG_NULL, 0, 0, 0 }};
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

static int
refine_evil_fn (p4est_t * p4est, p4est_topidx_t which_tree,
                p4est_quadrant_t * quadrant)
{
  if ((int) quadrant->level >= refine_level) {
    return 0;
  }
  if (p4est->mpirank <= 1) {
    return 1;
  }

  return 0;
}

static int
refine_evil3_fn (p4est_t * p4est, p4est_topidx_t which_tree,
                 p4est_quadrant_t * quadrant)
{
  p4est_qcoord_t      u2;
  p4est_quadrant_t    ref;

  P4EST_QUADRANT_INIT (&ref);

  u2 = P4EST_QUADRANT_LEN (2);

  if (which_tree == 0) {
    ref.x = 3 * u2;
    ref.y = 2 * u2;
  }
  else if (which_tree == 1) {
    ref.x = 2 * u2;
    ref.y = 3 * u2;
  }
  ref.level = 2;

  if ((int) quadrant->level >= refine_level) {
    return 0;
  }
  if ((which_tree == 0 || which_tree == 1) &&
      (p4est_quadrant_is_equal (&ref, quadrant) ||
       p4est_quadrant_is_ancestor (&ref, quadrant))) {
    return 1;
  }

  return 0;
}

static int
coarsen_evil_fn (p4est_t * p4est, p4est_topidx_t which_tree,
                 p4est_quadrant_t * q[])
{
  if (p4est->mpirank >= 2) {
    return 1;
  }

  return 0;
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
  usage =
    "Arguments: <configuration> <level>\n"
    "   Configuration can be any of\n"
    "      unit|three|evil|evil3|pillow|moebius|\n"
    "         star|cubed|disk|periodic|rotwrap\n"
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
    else if (!strcmp (argv[1], "evil")) {
      config = P4EST_CONFIG_EVIL;
    }
    else if (!strcmp (argv[1], "evil3")) {
      config = P4EST_CONFIG_EVIL3;
    }
    else if (!strcmp (argv[1], "pillow")) {
      config = P4EST_CONFIG_PILLOW;
    }
    else if (!strcmp (argv[1], "moebius")) {
      config = P4EST_CONFIG_MOEBIUS;
    }
    else if (!strcmp (argv[1], "star")) {
      config = P4EST_CONFIG_STAR;
    }
    else if (!strcmp (argv[1], "cubed")) {
      config = P4EST_CONFIG_CUBED;
    }
    else if (!strcmp (argv[1], "disk")) {
      config = P4EST_CONFIG_DISK;
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
  if (config == P4EST_CONFIG_EVIL) {
    refine_fn = refine_evil_fn;
    coarsen_fn = coarsen_evil_fn;
  }
  else if (config == P4EST_CONFIG_EVIL3) {
    refine_fn = refine_evil3_fn;
    coarsen_fn = NULL;
  }
  else {
    refine_fn = refine_normal_fn;
    coarsen_fn = NULL;
  }

  /* create connectivity and forest structures */
  if (config == P4EST_CONFIG_THREE || config == P4EST_CONFIG_EVIL3) {
    connectivity = p4est_connectivity_new_corner ();
  }
  else if (config == P4EST_CONFIG_PILLOW) {
    connectivity = p4est_connectivity_new_pillow ();
  }
  else if (config == P4EST_CONFIG_MOEBIUS) {
    connectivity = p4est_connectivity_new_moebius ();
  }
  else if (config == P4EST_CONFIG_STAR) {
    connectivity = p4est_connectivity_new_star ();
  }
  else if (config == P4EST_CONFIG_CUBED) {
    connectivity = p4est_connectivity_new_cubed ();
  }
  else if (config == P4EST_CONFIG_DISK) {
    connectivity = p4est_connectivity_new_disk ();
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
  p4est = p4est_new_ext (mpi->mpicomm, connectivity, 15, 0, 0,
                         sizeof (user_data_t), init_fn, NULL);
  p4est_vtk_write_file (p4est, NULL, "simple2_new");

  /* refinement and coarsening */
  p4est_refine (p4est, 1, refine_fn, init_fn);
  if (coarsen_fn != NULL) {
    p4est_coarsen (p4est, 1, coarsen_fn, init_fn);
  }
  p4est_vtk_write_file (p4est, NULL, "simple2_refined");

  /* balance */
  p4est_balance (p4est, P4EST_CONNECT_FULL, init_fn);
  p4est_vtk_write_file (p4est, NULL, "simple2_balanced");
  crc = p4est_checksum (p4est);

  /* partition */
  p4est_partition (p4est, 0, NULL);
  p4est_vtk_write_file (p4est, NULL, "simple2_partition");

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
  p4est_connectivity_destroy (connectivity);

  /* clean up and exit */
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
