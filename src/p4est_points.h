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

#ifndef P4EST_POINTS_H
#define P4EST_POINTS_H

#include <p4est.h>

/** Create a new p4est based on a distributed set of points.
 *
 * \param [in] mpicomm       A valid MPI_Comm or MPI_COMM_NULL.
 * \param [in] connectivity  This is the connectivity information that
 *                           the forest is built with.  Note the p4est
 *                           does not take ownership of the memory.
 * \param [in] maxlevel      Level of the smallest possible quadrants.
 * \param [in] points        Unsorted collection of clamped quadrant nodes.
 *                           The tree id must be stored in p.which_tree.
 * \param [in] num_points    Number of local points provided in the array.
 * \param [in] data_size     This is the size of data for each quadrant which
 *                           can be zero.  Then user_data_pool is set to NULL.
 * \param [in] init_fn       Callback function to initialize the user_data
 *                           which is already allocated automatically.
 * \param [in] user_pointer  Assign to the user_pointer member of the p4est
 *                           before init_fn is called the first time.
 *
 * \return This returns a vaild forest.
 *
 * \note The connectivity structure must not be destroyed
 *       during the lifetime of this p4est.
 */
p4est_t            *p4est_new_points (MPI_Comm mpicomm,
                                      p4est_connectivity_t * connectivity,
                                      int maxlevel, p4est_quadrant_t * points,
                                      p4est_locidx_t num_points,
                                      size_t data_size, p4est_init_t init_fn,
                                      void *user_pointer);

#endif /* !P4EST_POINTS_H */
