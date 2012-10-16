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

#ifndef P8EST_EXTENDED_H
#define P8EST_EXTENDED_H

#include <p8est.h>

SC_EXTERN_C_BEGIN;

/** Callback function prototype used by extended routines when the quadrants
 * of an existing, valid p8est are changed.  The callback allows the user to
 * make changes to newly initialized quadrants before the quadrants that they
 * replace are destroyed.
 *
 * \param [in] num_outgoing The number of outgoing quadrants.
 * \param [in] outgoing     The outgoing quadrants: after the callback, the
 *                          user_data, if \a p8est->data_size is nonzero,
 *                          will be destroyed.
 * \param [in] num_incoming The number of incoming quadrants.
 * \param [in/out] incoming The incoming quadrants: prior to the callback,
 *                          the user_data, if \a p8est->data_size is nonzero,
 *                          is allocated, and the p8est_init_t callback,
 *                          if it has been provided, will be called.
 *
 * Either num_outgoing or num_incoming will be 1, depending on whether the
 * change in the mesh is refinement or coarsening.  The incoming and outoing
 * quadrants always occupy the same region of the mesh.  See the description
 * of p8est_replace_type_t for more information.
 */
typedef void        (*p8est_replace_t) (p8est_t * p8est,
                                        p4est_topidx_t which_tree,
                                        p4est_locidx_t num_outgoing,
                                        p8est_quadrant_t * outgoing[],
                                        p4est_locidx_t num_incoming,
                                        p8est_quadrant_t * incoming[]);

/** A flag that is set by the user to define the capabilities of a particular
 * implementation of the p8est_replace_t prototype.
 *
 * Sometimes, the extended routines will change a mesh only by one level,
 * either replacing a single quadrant with its family of children
 * (p8est_refine_ext if refine_recursive if false), or replacing a family of
 * children with their parent (p8est_coarsen_ext if coarsen_recursive is
 * false).
 *
 * In other situations, the incoming or outgoing quadrants may contain
 * quadrants from multiple levels (when recursive refinement/coarsening
 * occurs, or when balancing).  Providing a p8est_replace_t callback that
 * handles multiple levels may be more efficient, but may also be more
 * difficult to implement.  If the user's callback can only handle one level
 * at a time, the user should set the flag to P8EST_REPLACE_FAMILY and the
 * callback will be called multiple times.
 *
 * For example, consider the following recursive refinement:
 *
 * o...............o     +-------+---+---+
 * :               :     |       |   |   |
 * :               :     |       +---+---+
 * :               :     |       |   |   |
 * :               : --> +-------+---+---+
 * :               :     |       |       |
 * :               :     |       |       |
 * :               :     |       |       |
 * o...............o     +-------+-------+
 *
 * If the flag is P8EST_REPLACE_BATCH, the p8est_replace_t callback will be
 * called once, with all of the descendants on the right in the incoming
 * array.
 *
 * If the flag is P8EST_REPLACE_FAMILY, the p8est_replace_t callback will be
 * called twice, first
 *
 * o...............o     +-------+-------+
 * :               :     |       |       |
 * :               :     |       |       |
 * :               :     |       |       |
 * :               : --> +-------+-------+ ,
 * :               :     |       |       |
 * :               :     |       |       |
 * :               :     |       |       |
 * o...............o     +-------+-------+
 *
 * and then
 *
 * o . . . o.......o     o . . . +---+---+
 * .       :       :     .       |   |   |
 * .       :       :     .       +---+---+
 * .       :       :     .       |   |   |
 * o . . . o.......o --> o . . . +---+---+ .
 * .       .       .     .       .       .
 * .       .       .     .       .       .
 * .       .       .     .       .       .
 * o . . . o . . . o     o . . . o . . . o
 */
typedef enum
{
  P8EST_REPLACE_FAMILY,
  P8EST_REPLACE_BATCH
}
p8est_replace_type_t;

/** Create a new forest.
 * This is a more general form of p8est_new.
 * See the documentation of p8est_new for basic usage.
 *
 * \param [in] min_quadrants    Minimum initial quadrants per processor.
 * \param [in] min_level        The forest is refined at least to this level.
 *                              May be negative or 0, then it has no effect.
 * \param [in] fill_uniform     If true, fill the forest with a uniform mesh
 *                              instead of the coarsest possible one.
 *                              The latter is partition-specific so that
 *                              is usually not a good idea.
 */
p8est_t            *p8est_new_ext (MPI_Comm mpicomm,
                                   p8est_connectivity_t * connectivity,
                                   p4est_locidx_t min_quadrants,
                                   int min_level, int fill_uniform,
                                   size_t data_size, p8est_init_t init_fn,
                                   void *user_pointer);

/** Refine a forest with a bounded maximum refinement level.
 * A quadrant is refined if either callback returns true.
 * \param [in] maxlevel   Maximum allowed refinement level (inclusive).
 *                        If this is negative the level is unrestricted.
 * \param [in] refine_fn  Callback function that returns true
 *                        if a quadrant shall be refined, may be NULL.
 * \param [in] init_fn    Callback function to initialize the user_data
 *                        which is guaranteed to be allocated, may be NULL.
 * \param [in] replace_fn Callback function that allows the user to change
 *                        incoming quadrants based on the quadrants they
 *                        replace.
 * \param [in] rtype      See description of p8est_replace_type_t
 */
void                p8est_refine_ext (p8est_t * p8est,
                                      int refine_recursive, int maxlevel,
                                      p8est_refine_t refine_fn,
                                      p8est_init_t init_fn,
                                      p8est_replace_t replace_fn,
                                      p8est_replace_type_t rtype);

/** Repartition the forest.
 *
 * The forest is partitioned between processors such that each processor
 * has an approximately equal number of quadrants (or weight).
 *
 * \param [in,out] p8est      The forest that will be partitioned.
 * \param [in]     partition_for_coarsening     If true, the partition
 *                            is modified to allow one level of coarsening.
 * \param [in]     weight_fn  A weighting function or NULL
 *                            for uniform partitioning.
 * \return         The global number of shipped quadrants
 */
p4est_gloidx_t      p8est_partition_ext (p8est_t * p8est,
                                         int partition_for_coarsening,
                                         p8est_weight_t weight_fn);

/** Save the complete connectivity/p8est data to disk.  This is a collective
 * operation that all MPI processes need to call.  All processes write
 * into the same file, so the filename given needs to be identical over
 * all parallel invocations.
 * See p8est_load_ext for information on the autopartition parameter.
 * \param [in] filename    Name of the file to write.
 * \param [in] p8est       Valid forest structure.
 * \param [in] save_data   If true, the element data is saved.
 *                         Otherwise, a data size of 0 is saved.
 * \param [in] save_partition   If false, save file as if 1 core was used.
 *                              If true, save core count and partition.
 *                         Advantage: Partition can be recovered on loading
 *                              with same mpisize and autopartition false.
 *                         Disadvantage: Makes the file depend on mpisize.
 *                  Either way the file can be loaded with autopartition true.
 * \note            Aborts on file errors.
 */
void                p8est_save_ext (const char *filename, p8est_t * p8est,
                                    int save_data, int save_partition);

/** Load the complete connectivity/p4est structure from disk.
 * It is possible to load the file with a different number of processors
 * than has been used to write it.  The partition will then be uniform.
 * \param [in] filename         Name of the file to read.
 * \param [in] mpicomm          A valid MPI communicator.
 * \param [in] data_size        Size of data for each quadrant which can be
 *                              zero.  Then user_data_pool is set to NULL.
 *                              If data_size is zero, load_data is ignored.
 * \param [in] load_data        If true, the element data is loaded.  This is
 *                              only permitted if the saved data size matches.
 *                              If false, the stored data size is ignored.
 * \param [in] autopartition    Ignore saved partition and make it uniform.
 * \param [in] broadcasthead    Have only rank 0 read headers and bcast them.
 * \param [in] user_pointer     Assign to the user_pointer member of the p4est
 *                              before init_fn is called the first time.
 * \param [out] connectivity    Connectivity must be destroyed separately.
 * \return          Returns a valid forest structure. A pointer to a valid
 *                  connectivity structure is returned through the last
 *                  argument.
 * \note            Aborts on file errors or invalid file contents.
 */
p8est_t            *p8est_load_ext (const char *filename, MPI_Comm mpicomm,
                                    size_t data_size, int load_data,
                                    int autopartition, int broadcasthead,
                                    void *user_pointer,
                                    p8est_connectivity_t ** connectivity);

SC_EXTERN_C_END;

#endif /* !P8EST_EXTENDED_H */
