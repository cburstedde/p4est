/*
  This file is part of p4est.
  p4est is a C library to manage a parallel collection of quadtrees and/or
  octrees.

  Copyright (C) 2007 Carsten Burstedde, Lucas Wilcox.

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
 * Usage: p4est_simple <level>
 */

#include <p4est_algorithms.h>
#include <p4est_base.h>
#include <p4est_vtk.h>

typedef struct
{
  int32_t             a;
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
init_fn (p4est_t * p4est, int32_t which_tree, p4est_quadrant_t * quadrant)
{
  user_data_t        *data = quadrant->user_data;

  data->a = which_tree;
}

static int
refine_fn (p4est_t * p4est, int32_t which_tree, p4est_quadrant_t * quadrant)
{
  if (quadrant->level >= (refine_level - (which_tree % 3))) {
    return 0;
  }
  if (quadrant->x == (1 << (P4EST_MAXLEVEL)) - (1 << (P4EST_MAXLEVEL - 2)) &&
      quadrant->y == (1 << (P4EST_MAXLEVEL)) - (1 << (P4EST_MAXLEVEL - 2))) {
    return 1;
  }
  if (quadrant->x >= (1 << (P4EST_MAXLEVEL - 2))) {
    return 0;
  }

  return 1;
}

static void
abort_fn (void * data) 
{
#ifdef HAVE_MPI
  int                 mpiret;
#endif
  mpi_context_t      *mpi = data;

  fprintf (stderr, "[%d] p4est_simple abort handler\n", mpi->mpirank);

#ifdef HAVE_MPI
  mpiret = MPI_Abort (mpi->mpicomm, 1);
  P4EST_CHECK_MPI (mpiret);
#endif
}

int
main (int argc, char **argv)
{
#ifdef HAVE_MPI
  int                 use_mpi = 1;
  int                 mpiret;
#endif
  mpi_context_t       mpi_context, *mpi = &mpi_context;
  p4est_t            *p4est;
  p4est_connectivity_t *connectivity;

  /* initialize MPI */
  mpi->mpirank = 0;
  mpi->mpicomm = MPI_COMM_NULL;
#ifdef HAVE_MPI
  if (use_mpi) {
    mpiret = MPI_Init (&argc, &argv);
    P4EST_CHECK_MPI (mpiret);
    mpi->mpicomm = MPI_COMM_WORLD;
    mpiret = MPI_Comm_rank (mpi->mpicomm, &mpi->mpirank);
    P4EST_CHECK_MPI (mpiret);
  }
#endif

  /* register MPI abort handler */
  p4est_set_abort_handler (mpi->mpirank, abort_fn, mpi);

  /* get command line argument: maximum refinement level */
  if (mpi->mpirank == 0) {
    P4EST_CHECK_ABORT (argc == 2, "Give level");
  }
  refine_level = atoi (argv[1]);

  /* create connectivity and forest structures */
  /* connectivity = p4est_connectivity_new_unitsquare (); */
  connectivity = p4est_connectivity_new_corner ();
  p4est = p4est_new (mpi->mpicomm, stdout, connectivity,
                     sizeof (user_data_t), init_fn);
  p4est_tree_print (p4est_array_index (p4est->trees, 0),
                    mpi->mpirank, stdout);
  p4est_vtk_write_file (p4est, "mesh_simple_new");
  p4est_refine (p4est, refine_fn, init_fn);
  p4est_vtk_write_file (p4est, "mesh_simple_refined");
  p4est_balance (p4est, init_fn);
  p4est_vtk_write_file (p4est, "mesh_simple_balanced");

  /* destroy the p4est and its connectivity structure */
  p4est_destroy (p4est);
  p4est_connectivity_destroy (connectivity);

  /* clean up and exit */
  p4est_memory_check ();

#ifdef HAVE_MPI
  if (use_mpi) {
    mpiret = MPI_Finalize ();
    P4EST_CHECK_MPI (mpiret);
  }
#endif

  return 0;
}

/* EOF simple.c */
