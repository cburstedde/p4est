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

/** \file p6est_communication.h
 *
 * MPI_Comm management.
 *
 * \ingroup p6est
 */

#include <p6est.h>

/** duplicate the mpi communicator in the p6est */
void                p6est_comm_parallel_env_create (p6est_t * p6est,
                                                    sc_MPI_Comm mpicomm);

/** free an owned mpi communicator */
void                p6est_comm_parallel_env_free (p6est_t * p6est);

/** test if a p6est has an mpi communicator */
int                 p6est_comm_parallel_env_is_null (p6est_t * p6est);

/** assign an mpi communicator */
void                p6est_comm_parallel_env_assign (p6est_t * p6est,
                                                    sc_MPI_Comm mpicomm);
