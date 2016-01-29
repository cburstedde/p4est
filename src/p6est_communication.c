/*
  This file is part of p4est.
  p4est is a C library to manage a collection (a forest) of multiple
  connected adaptive quadtrees or octrees in parallel.

  Copyright (C) 2015 The University of Texas System
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

#include <p6est_communication.h>

void
p6est_comm_parallel_env_create (p6est_t * p6est, sc_MPI_Comm mpicomm)
{
  int                 mpiret;

  /* duplicate MPI communicator */
  mpiret = sc_MPI_Comm_dup (mpicomm, &(p6est->mpicomm));
  SC_CHECK_MPI (mpiret);
  p6est->mpicomm_owned = 1;

  /* retrieve MPI information */
  mpiret = sc_MPI_Comm_size (mpicomm, &(p6est->mpisize));
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (mpicomm, &(p6est->mpirank));
  SC_CHECK_MPI (mpiret);
}

void
p6est_comm_parallel_env_free (p6est_t * p6est)
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

int
p6est_comm_parallel_env_is_null (p6est_t * p6est)
{
  return (p6est->mpicomm == sc_MPI_COMM_NULL);
}

void
p6est_comm_parallel_env_assign (p6est_t * p6est, sc_MPI_Comm mpicomm)
{
  int                 mpiret;
  int                 result;

  mpiret = sc_MPI_Comm_compare (p6est->mpicomm, mpicomm, &result);
  SC_CHECK_MPI (mpiret);
  /* check if input MPI communicator has same size and same rank order */
  if (result == sc_MPI_IDENT) {
    return;
  }

  P4EST_ASSERT (result == sc_MPI_CONGRUENT);

  /* free the current parallel environment */
  p6est_comm_parallel_env_free (p6est);

  /* assign MPI communicator of input, it is therefore not owned */
  p6est->mpicomm = mpicomm;
  p6est->mpicomm_owned = 0;

  /* retrieve MPI information */
  mpiret = sc_MPI_Comm_size (mpicomm, &(p6est->mpisize));
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (mpicomm, &(p6est->mpirank));
  SC_CHECK_MPI (mpiret);
}
