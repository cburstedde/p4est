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

#include <p4est_communication.h>
#include <p6est_communication.h>

void
p6est_comm_parallel_env_assign (p6est_t * p6est, sc_MPI_Comm mpicomm)
{
  /* set MPI communicator */
  p6est->mpicomm = mpicomm;
  p6est->mpicomm_owned = 0;

  /* retrieve MPI information */
  p6est_comm_parallel_env_get_info (p6est);
}

void
p6est_comm_parallel_env_duplicate (p6est_t * p6est)
{
  sc_MPI_Comm         mpicomm = p6est->mpicomm;
  int                 mpiret;

  /* duplicate MPI communicator */
  mpiret = sc_MPI_Comm_dup (mpicomm, &(p6est->mpicomm));
  SC_CHECK_MPI (mpiret);
  p6est->mpicomm_owned = 1;
}

void
p6est_comm_parallel_env_release (p6est_t * p6est)
{
  int                 mpiret;

  /* free MPI communicator if it's owned */
  if (p6est->mpicomm_owned) {
    mpiret = sc_MPI_Comm_free (&(p6est->mpicomm));
    SC_CHECK_MPI (mpiret);
  }
  p6est->mpicomm = sc_MPI_COMM_NULL;
  p6est->mpicomm_owned = 0;

  /* set MPI information */
  p6est->mpisize = 0;
  p6est->mpirank = sc_MPI_UNDEFINED;
}

void
p6est_comm_parallel_env_replace (p6est_t * p6est, sc_MPI_Comm mpicomm)
{
  /* check if input MPI communicator has same size and same rank order */
#ifdef P4EST_ENABLE_DEBUG
  {
    int                 mpiret, result;

    mpiret = sc_MPI_Comm_compare (p6est->mpicomm, mpicomm, &result);
    SC_CHECK_MPI (mpiret);

    P4EST_ASSERT (result == sc_MPI_IDENT || result == sc_MPI_CONGRUENT);
  }
#endif

  /* release the current parallel environment */
  p6est_comm_parallel_env_release (p6est);

  /* assign new MPI communicator */
  p6est_comm_parallel_env_assign (p6est, mpicomm);
}

void
p6est_comm_parallel_env_get_info (p6est_t * p6est)
{
  int                 mpiret;

  mpiret = sc_MPI_Comm_size (p6est->mpicomm, &(p6est->mpisize));
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (p6est->mpicomm, &(p6est->mpirank));
  SC_CHECK_MPI (mpiret);
}

int
p6est_comm_parallel_env_is_null (p6est_t * p6est)
{
  return (p6est->mpicomm == sc_MPI_COMM_NULL);
}

int
p6est_comm_parallel_env_reduce (p6est_t ** p6est_supercomm)
{
  return p6est_comm_parallel_env_reduce_ext (p6est_supercomm,
                                             sc_MPI_GROUP_NULL, 0, NULL);
}

int
p6est_comm_parallel_env_reduce_ext (p6est_t ** p6est_supercomm,
                                    sc_MPI_Group group_add,
                                    int add_to_beginning, int **ranks_subcomm)
{
  p6est_t            *p6est = *p6est_supercomm;
  int                 mpisize = p6est->mpisize;
  int                 mpiret;
  p4est_gloidx_t     *global_first_layer = p6est->global_first_layer;

  p4est_gloidx_t     *n_quadrants;
  int                 submpisize;
  sc_MPI_Comm         submpicomm;
  int                *ranks;
  int                 i;
  int                 is_nonempty;

  /* reduce MPI communicator of column layout */
  is_nonempty =
    p4est_comm_parallel_env_reduce_ext (&(p6est->columns), group_add,
                                        add_to_beginning, &ranks);

  /* destroy p4est and exit if this rank is empty */
  if (!is_nonempty) {
    p6est->columns = NULL;
    p6est_destroy (p6est);
    *p6est_supercomm = NULL;
    if (ranks_subcomm) {
      *ranks_subcomm = NULL;
    }
    P4EST_ASSERT (ranks == NULL);
    return 0;
  }

  /* get sub-communicator */
  submpicomm = p6est->columns->mpicomm;

  /* update size of new MPI communicator */
  mpiret = sc_MPI_Comm_size (submpicomm, &submpisize);
  SC_CHECK_MPI (mpiret);
  if (submpisize == p6est->mpisize) {
    P4EST_ASSERT (ranks == NULL);
    return 1;
  }

  /* set new parallel environment */
  p6est_comm_parallel_env_release (p6est);
  p6est_comm_parallel_env_assign (p6est, submpicomm);
  if (p6est->columns->mpicomm_owned) {
    p6est->columns->mpicomm_owned = 0;
    p6est->mpicomm_owned = 1;
  }
  P4EST_ASSERT (p6est->mpisize == submpisize);

  /* create array of non-empty processes that will be included to sub-comm */
  n_quadrants = P4EST_ALLOC (p4est_gloidx_t, mpisize);
  for (i = 0; i < mpisize; i++) {
    n_quadrants[i] = global_first_layer[i + 1] - global_first_layer[i];
  }

  /* allocate and set global layer count */
  P4EST_FREE (p6est->global_first_layer);
  p6est->global_first_layer = P4EST_ALLOC (p4est_gloidx_t, submpisize + 1);
  p6est->global_first_layer[0] = 0;
  for (i = 0; i < submpisize; i++) {
    P4EST_ASSERT (ranks[i] != sc_MPI_UNDEFINED);
    P4EST_ASSERT (group_add != sc_MPI_GROUP_NULL
                  || 0 < n_quadrants[ranks[i]]);
    p6est->global_first_layer[i + 1] =
      p6est->global_first_layer[i] + n_quadrants[ranks[i]];
  }
  P4EST_FREE (n_quadrants);
  if (ranks_subcomm) {
    *ranks_subcomm = ranks;
  }
  else {
    P4EST_FREE (ranks);
  }

  /* return that p6est exists on this rank */
  return 1;
}
