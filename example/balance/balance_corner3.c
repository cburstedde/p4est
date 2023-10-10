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

#include <p8est_bits.h>
#include <p8est_extended.h>
#include <p8est_vtk.h>

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

/* refinement level initialization */
static int          refine_level = 0;

/* connectivity initialization function */
static void
init_fn (p8est_t * p8est, p4est_topidx_t which_tree,
         p8est_quadrant_t * quadrant)
{
  user_data_t        *data = (user_data_t *) quadrant->p.user_data;

  data->a = which_tree;
}

/* refinement function */
static int
refine_fn (p8est_t * p8est, p4est_topidx_t which_tree,
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
  if (quadrant->x >= P8EST_QUADRANT_LEN (2)) {
    return 0;
  }

  return 1;
}

int
main (int argc, char **argv)
{
  int                 mpiret;
  int                 wrongusage;
  const char         *usage;
  p8est_connectivity_t *connectivity;
  mpi_context_t       mpi_context, *mpi = &mpi_context;
  p8est_geometry_t   *geom;
  p8est_t            *p8est;

  /* initialize MPI and p8est internals */
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  mpi->mpicomm = sc_MPI_COMM_WORLD;
  mpiret = sc_MPI_Comm_size (mpi->mpicomm, &mpi->mpisize);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (mpi->mpicomm, &mpi->mpirank);
  SC_CHECK_MPI (mpiret);

  sc_init (mpi->mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  /* usage error if the input is not in the correct format */
  usage =
    "Arguments: <level>\n"
    "   Level: controls the maximum depth of refinement\n";
  wrongusage = 0;
  if (!wrongusage && argc != 2) {
    wrongusage = 1;
  }
  if (wrongusage) {
    P4EST_GLOBAL_LERROR (usage);
    sc_abort_collective ("Usage error");
  }

  /* assign variables based on configuration */
  refine_level = atoi (argv[1]);

  /* create connectivity and forest structures */
  geom = NULL;
  connectivity = p8est_connectivity_new_drop ();
  p8est = p8est_new_ext (mpi->mpicomm, connectivity, 15, 0, 0,
                         sizeof (user_data_t), init_fn, geom);

  /* refinement */
  p8est_refine (p8est, 1, refine_fn, init_fn);

  /* balance */
  p8est_balance (p8est, P8EST_CONNECT_FULL, init_fn);

  /* partition */
  p8est_partition (p8est, 0, NULL);
  p8est_vtk_write_file (p8est, geom, "balance_corner3_partition");

  /* destroy p8est and its connectivity */
  p8est_destroy (p8est);
  p8est_connectivity_destroy (connectivity);
  /* clean up and exit */
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
