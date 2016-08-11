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

#ifndef P6EST_COMMUNICATION_H
#define P6EST_COMMUNICATION_H

/** \file p6est_communication.h
 *
 * MPI_Comm management.
 *
 * \ingroup p6est
 */

#include <p6est.h>

SC_EXTERN_C_BEGIN;

/** Assign an MPI communicator to p6est; retrieve parallel environment.
 *
 * \param [in] mpicomm    A valid MPI communicator.
 *
 * \note The provided MPI communicator is not owned by p6est.
 */
void                p6est_comm_parallel_env_assign (p6est_t * p6est,
                                                    sc_MPI_Comm mpicomm);

/** Duplicate MPI communicator and replace the current one by the duplicate.
 *
 * \note The duplicated MPI communicator is owned by p6est.
 */
void                p6est_comm_parallel_env_duplicate (p6est_t * p6est);

/** Release MPI communicator if it is owned by p6est.
 */
void                p6est_comm_parallel_env_release (p6est_t * p6est);

/** Replace the current MPI communicator by the one provided as input.
 *
 * \param [in] mpicomm    A valid MPI communicator.
 *
 * \note The provided MPI communicator is not owned by p6est.
 */
void                p6est_comm_parallel_env_replace (p6est_t * p6est,
                                                     sc_MPI_Comm mpicomm);

/** Retrieve parallel environment information.
 */
void                p6est_comm_parallel_env_get_info (p6est_t * p6est);

/** Check if the MPI communicator is valid.
 *
 * \return True if communicator is not NULL communicator, false otherwise.
 */
int                 p6est_comm_parallel_env_is_null (p6est_t * p6est);

/** Reduce MPI communicator to non-empty ranks (i.e., nonzero quadrant counts).
 *
 * \param [in/out] p6est_supercomm  Object which communicator is reduced.
 *                                  points to NULL if this p6est does not
 *                                  exists.
 *
 * \return True if p6est exists on this MPI rank after reduction.
 */
int                 p6est_comm_parallel_env_reduce (p6est_t **
                                                    p6est_supercomm);

/** Reduce MPI communicator to non-empty ranks and add a group of ranks that
 * will remain in the reduced communicator regardless whether they are empty
 * or not.
 *
 * \param [in/out] p6est_supercomm  Object which communicator is reduced.
 *                                  Points to NULL if this p6est does not
 *                                  exists.
 * \param [in] group_add         Group of ranks that will remain in
 *                               communicator.
 * \param [in] add_to_beginning  If true, ranks will be added to the beginning
 *                               of the reduced communicator, otherwise to the
 *                               end.
 * \param[out] ranks_subcomm     If not null, array of size 'subcommsize' with
 *                               subcommrank->supercommrank map.
 *
 * \return True if p6est exists on this MPI rank after reduction.
 */
int                 p6est_comm_parallel_env_reduce_ext (p6est_t **
                                                        p6est_supercomm,
                                                        sc_MPI_Group
                                                        group_add,
                                                        int add_to_beginning,
                                                        int **ranks_subcomm);

SC_EXTERN_C_END;

#endif /* !P6EST_COMMUNICATION_H */
