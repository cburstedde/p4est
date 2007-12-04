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

#include <p4est_algorithms.h>
#include <p4est_base.h>

typedef struct
{
  int32_t             a;
}
user_data_t;

static void
init_fn (p4est_t * p4est, int32_t which_tree, p4est_quadrant_t * quadrant)
{
  user_data_t        *data = quadrant->user_data;

  data->a = which_tree;
}

int
main (int argc, char **argv)
{
#ifdef HAVE_MPI
  int                 use_mpi = 1;
  int                 mpiret;
#endif
  MPI_Comm            mpicomm;
  p4est_t            *p4est;
  p4est_connectivity_t *connectivity;

  mpicomm = MPI_COMM_NULL;
#ifdef HAVE_MPI
  if (use_mpi) {
    mpiret = MPI_Init (&argc, &argv);
    P4EST_CHECK_MPI (mpiret);
    mpicomm = MPI_COMM_WORLD;
  }
#endif

  /* create connectivity and forest structures */
  connectivity = p4est_connectivity_new_unitsquare ();
  p4est = p4est_new (mpicomm, stdout, connectivity,
                     sizeof (user_data_t), init_fn);
  p4est_tree_print (p4est_array_index (p4est->trees, 0),
                    p4est->mpirank, stdout);

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
