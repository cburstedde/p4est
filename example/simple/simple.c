/*
  This file is part of p4est.
  p4est is a C library to manage a parallel collection of quadtrees and/or
  octrees.

  Copyright (C) 2007,2008 Carsten Burstedde, Lucas Wilcox.

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/*
 * Usage: p4est_simple <configuration> <level>
 *        possible configurations:
 *        o unit      Refinement on the unit square.
 *        o three     Refinement on a forest with three trees.
 *        o evil      Check second round of refinement with np=5 level=7
 *        o evil3     Check second round of refinement on three trees
 *        o moebius   Refinement on a 5-tree Moebius band.
 *        o star      Refinement on a 6-tree star shaped domain.
 *        o periodic  Refinement on the unit square with periodic b.c.
 */

#include <p4est_algorithms.h>
#include <p4est_vtk.h>

enum
{
  P4EST_CONFIG_NULL,
  P4EST_CONFIG_UNIT,
  P4EST_CONFIG_THREE,
  P4EST_CONFIG_EVIL,
  P4EST_CONFIG_EVIL3,
  P4EST_CONFIG_MOEBIUS,
  P4EST_CONFIG_STAR,
  P4EST_CONFIG_PERIODIC,
};

typedef struct
{
  p4est_topidx_t      a;
}
user_data_t;

typedef struct
{
  MPI_Comm            mpicomm;
  int                 mpirank;
}
mpi_context_t;

static int          refine_level = 0;

static void
init_fn (p4est_t * p4est, p4est_topidx_t which_tree,
         p4est_quadrant_t * quadrant)
{
  user_data_t        *data = quadrant->p.user_data;

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
                 p4est_quadrant_t * q0, p4est_quadrant_t * q1,
                 p4est_quadrant_t * q2, p4est_quadrant_t * q3)
{
  if (p4est->mpirank >= 2) {
    return 1;
  }

  return 0;
}

static void
abort_fn (void *data)
{
  int                 mpiret;
  mpi_context_t      *mpi = data;

  fprintf (stderr, "[%d] p4est_simple abort handler\n", mpi->mpirank);

  mpiret = MPI_Abort (mpi->mpicomm, 1);
  SC_CHECK_MPI (mpiret);
}

int
main (int argc, char **argv)
{
  int                 mpiret;
  int                 wrongusage, config;
  unsigned            crc;
  const char         *usage, *errmsg;
  mpi_context_t       mpi_context, *mpi = &mpi_context;
  p4est_t            *p4est;
  p4est_connectivity_t *connectivity;
  p4est_refine_t      refine_fn;
  p4est_coarsen_t     coarsen_fn;

  /* initialize MPI and p4est internals */
  mpi->mpicomm = MPI_COMM_WORLD;
  mpiret = MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  mpi->mpicomm = MPI_COMM_WORLD;
  mpiret = MPI_Comm_rank (mpi->mpicomm, &mpi->mpirank);
  SC_CHECK_MPI (mpiret);

  sc_init (mpi->mpirank, abort_fn, mpi, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  /* process command line arguments */
  usage =
    "Arguments: <configuration> <level>\n"
    "   Configuration can be any of\n"
    "      unit|three|evil|evil3|moebius|star|periodic\n"
    "   Level controls the maximum depth of refinement\n";
  errmsg = NULL;
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
    else if (!strcmp (argv[1], "moebius")) {
      config = P4EST_CONFIG_MOEBIUS;
    }
    else if (!strcmp (argv[1], "star")) {
      config = P4EST_CONFIG_STAR;
    }
    else if (!strcmp (argv[1], "periodic")) {
      config = P4EST_CONFIG_PERIODIC;
    }
    else {
      wrongusage = 1;
    }
  }
  if (wrongusage) {
    if (mpi->mpirank == 0) {
      fputs ("Usage error\n", stderr);
      fputs (usage, stderr);
      if (errmsg != NULL) {
        fputs (errmsg, stderr);
      }
      sc_abort ();
    }
    mpiret = MPI_Barrier (mpi->mpicomm);
    SC_CHECK_MPI (mpiret);
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
  else if (config == P4EST_CONFIG_MOEBIUS) {
    connectivity = p4est_connectivity_new_moebius ();
  }
  else if (config == P4EST_CONFIG_STAR) {
    connectivity = p4est_connectivity_new_star ();
  }
  else if (config == P4EST_CONFIG_PERIODIC) {
    connectivity = p4est_connectivity_new_periodic ();
  }
  else {
    connectivity = p4est_connectivity_new_unitsquare ();
  }
  /* p4est_initial_quadrants_per_processor = 0; */
  p4est = p4est_new (mpi->mpicomm, connectivity,
                     sizeof (user_data_t), init_fn);
  p4est_vtk_write_file (p4est, "mesh_simple_new");

  /* refinement and coarsening */
  p4est_refine (p4est, refine_fn, init_fn);
  if (coarsen_fn != NULL) {
    p4est_coarsen (p4est, coarsen_fn, init_fn);
  }
  p4est_vtk_write_file (p4est, "mesh_simple_refined");

  /* balance */
  p4est_balance (p4est, init_fn);
  p4est_vtk_write_file (p4est, "mesh_simple_balanced");
  crc = p4est_checksum (p4est);

  /* partition */
  p4est_partition (p4est, NULL);
  p4est_vtk_write_file (p4est, "mesh_simple_partition");

  /* print forest checksum */
  P4EST_GLOBAL_STATISTICSF ("Tree checksum 0x%x\n", crc);

  /* destroy the p4est and its connectivity structure */
  p4est_destroy (p4est);
  p4est_connectivity_destroy (connectivity);

  /* clean up and exit */
  sc_finalize ();

  mpiret = MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}

/* EOF simple.c */
