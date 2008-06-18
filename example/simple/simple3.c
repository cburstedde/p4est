/*
  This file is part of p4est.
  p4est is a C library to manage a parallel collection of quadtrees and/or
  octrees.

  Copyright (C) 2008 Carsten Burstedde, Lucas Wilcox.

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
 * Usage: p8est_simple <configuration> <level>
 *        possible configurations:
 *        o unit      Refinement on the unit cube.
 *        o periodic  Refinement on the unit square with periodic b.c.
 */

#define VTK_OUTPUT 1

#include <p8est_bits.h>

#ifdef VTK_OUTPUT
#include <p8est_vtk.h>
#endif

enum
{
  P8EST_CONFIG_NULL,
  P8EST_CONFIG_UNIT,
  P8EST_CONFIG_PERIODIC,
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
init_fn (p8est_t * p8est, p4est_topidx_t which_tree,
         p8est_quadrant_t * quadrant)
{
  user_data_t        *data = quadrant->p.user_data;

  data->a = which_tree;
}

static int
refine_normal_fn (p8est_t * p8est, p4est_topidx_t which_tree,
                  p8est_quadrant_t * quadrant)
{
  if ((int) quadrant->level >= (refine_level - (int) (which_tree % 3))) {
    return 0;
  }
  if (quadrant->level == 1 && p8est_quadrant_child_id (quadrant) == 3) {
    return 1;
  }
  if (quadrant->x == P8EST_LAST_OFFSET (2) &&
      quadrant->y == P8EST_LAST_OFFSET (2)) {
    return 1;
  }
  if (quadrant->z >= P8EST_QUADRANT_LEN (2)) {
    return 0;
  }

  return 1;
}

static void
abort_fn (void *data)
{
  int                 mpiret;
  mpi_context_t      *mpi = data;

  fprintf (stderr, "[%d] p8est_simple abort handler\n", mpi->mpirank);

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
  p8est_t            *p8est;
  p8est_connectivity_t *connectivity;
  p8est_refine_t      refine_fn;
  p8est_coarsen_t     coarsen_fn;

  /* initialize MPI and p4est internals */
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
    "      unit|periodic\n"
    "   Level controls the maximum depth of refinement\n";
  errmsg = NULL;
  wrongusage = 0;
  config = P8EST_CONFIG_NULL;
  if (!wrongusage && argc != 3) {
    wrongusage = 1;
  }
  if (!wrongusage) {
    if (!strcmp (argv[1], "unit")) {
      config = P8EST_CONFIG_UNIT;
    }
    else if (!strcmp (argv[1], "periodic")) {
      config = P8EST_CONFIG_PERIODIC;
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
  refine_fn = refine_normal_fn;
  coarsen_fn = NULL;

  /* create connectivity and forest structures */
  if (config == P8EST_CONFIG_PERIODIC) {
    connectivity = p8est_connectivity_new_periodic ();
  }
  else {
    connectivity = p8est_connectivity_new_unitcube ();
  }
  p8est = p8est_new (mpi->mpicomm, connectivity,
                     sizeof (user_data_t), init_fn);
#ifdef VTK_OUTPUT
  p8est_vtk_write_file (p8est, "mesh_simple3_new");
#endif

  /* refinement and coarsening */
  p8est_refine (p8est, refine_fn, init_fn);
  if (coarsen_fn != NULL) {
    p8est_coarsen (p8est, coarsen_fn, init_fn);
  }
#ifdef VTK_OUTPUT
  p8est_vtk_write_file (p8est, "mesh_simple3_refined");
#endif

  /* balance */
  p8est_balance (p8est, init_fn);
#ifdef VTK_OUTPUT
  p8est_vtk_write_file (p8est, "mesh_simple3_balanced");
#endif

  crc = p8est_checksum (p8est);

  /* partition */
  p8est_partition (p8est, NULL);
#ifdef VTK_OUTPUT
  p8est_vtk_write_file (p8est, "mesh_simple3_partition");
#endif

  /* print forest checksum */
  P4EST_GLOBAL_STATISTICSF ("Tree checksum 0x%x\n", crc);

  /* destroy the p8est and its connectivity structure */
  p8est_destroy (p8est);
  p8est_connectivity_destroy (connectivity);

  /* clean up and exit */
  sc_finalize ();

  mpiret = MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}

/* EOF simple3.c */
