/*
  This file is part of p4est.
  p4est is a C library to manage a collection (a forest) of multiple
  connected adaptive quadtrees or octrees in parallel.

  Copyright (C) 2010 The University of Texas System
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

/********************************************************************
 *                          IMPORTANT NOTE                          *
 *                                                                  *
 * These interfaces are intended for those who like finer control.  *
 * The API offers extended versions of some basic p4est functions.  *
 * The API may change without notice.                               *
 ********************************************************************/

#ifndef P4EST_EXTENDED_H
#define P4EST_EXTENDED_H

#include <p4est.h>

/** Extended callback function prototype to decide for refinement.
 * \param [in] quadrant     Quadrant under consideration for refinement.
 * \param [in] children     If true is returned these will be the children.
 * \return                  Nonzero if the quadrant shall be refined.
 */
typedef int         (*p4est_refine_ext_t) (p4est_t * p4est,
                                           p4est_topidx_t which_tree,
                                           p4est_quadrant_t * quadrant,
                                           p4est_quadrant_t * children[]);

/** Create a new forest.
 * This is a more general form of p4est_new.
 * See the documentation of p4est_new for basic usage.
 *
 * \param [in] min_quadrants    Minimum initial quadrants per processor.
 * \param [in] min_level        The forest is refined at least to this level.
 *                              May be negative or 0, then it has no effect.
 * \param [in] fill_uniform     If true, fill the forest with a uniform mesh
 *                              instead of the coarsest possible one.
 *                              The latter is partition-specific so that
 *                              is usually not a good idea.
 */
p4est_t            *p4est_new_ext (MPI_Comm mpicomm,
                                   p4est_connectivity_t * connectivity,
                                   p4est_locidx_t min_quadrants,
                                   int min_level, int fill_uniform,
                                   size_t data_size, p4est_init_t init_fn,
                                   void *user_pointer);

/** Refine a forest with a bounded maximum refinement level.
 * A quadrant is refined if either callback returns true.
 * \param [in] refine_fn Callback function that returns true
 *                       if a quadrant shall be refined, may be NULL.
 * \param [in] refine_fn Extended callback function that returns true
 *                       if a quadrant shall be refined, may be NULL.
 * \param [in] init_fn   Callback function to initialize the user_data
 *                       which is guaranteed to be allocated, may be NULL.
 * \param [in] maxlevel  Maximum allowed refinement level (inclusive).
 *                       If this is negative the level is unrestricted.
 */
void                p4est_refine_ext (p4est_t * p4est,
                                      int refine_recursive, int maxlevel,
                                      p4est_refine_t refine_fn,
                                      p4est_refine_ext_t refine_ext_fn,
                                      p4est_init_t init_fn);

/** Repartition the forest.
 *
 * The forest is partitioned between processors such that each processor
 * has an approximately equal number of quadrants (or weight).
 *
 * \param [in,out] p4est      The forest that will be partitioned.
 * \param [in]     partition_for_coarsening     If true, the partition
 *                            is modified to allow one level of coarsening.
 * \param [in]     weight_fn  A weighting function or NULL
 *                            for uniform partitioning.
 */
void                p4est_partition_ext (p4est_t * p4est,
                                         int partition_for_coarsening,
                                         p4est_weight_t weight_fn);

#endif /* !P4EST_EXTENDED_H */
