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

#ifndef P6EST_EXTENDED_H
#define P6EST_EXTENDED_H

/********************************************************************
 *                          IMPORTANT NOTE                          *
 *                                                                  *
 * These interfaces are intended for those who like finer control.  *
 * The API offers extended versions of some basic p6est functions.  *
 * The API may change without notice.                               *
 ********************************************************************/

/** \file p6est_extended.h
 *
 * Interface routines with extended capabilities.
 *
 * \ingroup p6est
 */

#include <p6est.h>

SC_EXTERN_C_BEGIN;

/** Create a new forest.
 * This is a more general form of p6est_new().
 * See the documentation of p6est_new() for basic usage.
 *
 * \param [in] min_quadrants    Minimum initial quadrants per processor.
 *                              Makes the refinement pattern mpisize-specific.
 * \param [in] min_level        The forest is horizontally refined at least to
 *                              this level.  May be negative or 0, then it has
 *                              no effect.
 * \param [in] min_zlevel       The forest is vertically refined at least to
 *                              this level.  May be negative or 0, then it has
 *                              no effect.
 * \parem [in] num_zroot        The number of "root" vertical layers
 *                              (used when non-power-of-2 layers are desired)
 * \param [in] fill_uniform     If true, fill the forest with a uniform mesh
 *                              instead of the coarsest possible one.
 *                              The latter is partition-specific so that
 *                              is usually not a good idea.
 */
p6est_t            *p6est_new_ext (sc_MPI_Comm mpicomm,
                                   p6est_connectivity_t * connectivity,
                                   p4est_locidx_t min_quadrants,
                                   int min_level, int min_zlevel,
                                   int num_zroot,
                                   int fill_uniform, size_t data_size,
                                   p6est_init_t init_fn, void *user_pointer);

/** Make a deep copy of a p6est.
 * The connectivity is not duplicated.
 * Copying of quadrant user data is optional.
 * If old and new data sizes are 0, the user_data field is copied regardless.
 * The inspect member of the copy is set to NULL.
 *
 * \param [in]  copy_data  If true, data are copied.
 *                         If false, data_size is set to 0.
 * \param [in]  duplicate_mpicomm  If true, MPI communicator is copied.
 * \return  Returns a valid p6est that does not depend on the input.
 */
p6est_t            *p6est_copy_ext (p6est_t * input, int copy_data,
                                    int duplicate_mpicomm);

/** Save the complete connectivity/p6est data to disk.
 *
 * This is a collective operation that all MPI processes need to call.  All
 * processes write into the same file, so the filename given needs to be
 * identical over all parallel invocations.  See p6est_load_ext() for
 * information on the autopartition parameter.
 *
 * \param [in] filename    Name of the file to write.
 * \param [in] p6est       Valid forest structure.
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
void                p6est_save_ext (const char *filename, p6est_t * p6est,
                                    int save_data, int save_partition);

/** Load the complete connectivity/p6est structure from disk.
 *
 * It is possible to load the file with a different number of processors
 * than has been used to write it.  The partition will then be uniform.
 *
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
 * \param [in] user_pointer     Assign to the user_pointer member of the p6est
 *                              before init_fn is called the first time.
 * \param [out] connectivity    Connectivity must be destroyed separately.
 * \return          Returns a valid forest structure. A pointer to a valid
 *                  connectivity structure is returned through the last
 *                  argument.
 * \note            Aborts on file errors or invalid file contents.
 */
p6est_t            *p6est_load_ext (const char *filename, sc_MPI_Comm mpicomm,
                                    size_t data_size, int load_data,
                                    int autopartition, int broadcasthead,
                                    void *user_pointer,
                                    p6est_connectivity_t ** connectivity);

/** Horizontally refine a forest with a bounded refinement level and a replace option.
 *
 * \param [in,out] p6est The forest is changed in place.
 * \param [in] refine_recursive Boolean to decide on recursive refinement.
 * \param [in] maxlevel   Maximum allowed refinement level (inclusive).
 *                        If this is negative the level is restricted only
 *                        by the compile-time constant QMAXLEVEL in p4est.h.
 * \param [in] refine_fn  Callback function that must return true if a quadrant
 *                        shall be refined.  If refine_recursive is true,
 *                        refine_fn is called for every existing and newly
 *                        created quadrant.  Otherwise, it is called for every
 *                        existing quadrant.  It is possible that a refinement
 *                        request made by the callback is ignored.  To catch
 *                        this case, you can examine whether init_fn or
 *                        replace_fn gets called.
 * \param [in] init_fn    Callback function to initialize the user_data for
 *                        newly created quadrants, which is guaranteed to be
 *                        allocated.  This function pointer may be NULL.
 * \param [in] replace_fn Callback function that allows the user to change
 *                        incoming quadrants based on the quadrants they
 *                        replace; may be NULL.
 */
void                p6est_refine_columns_ext (p6est_t * p6est,
                                              int refine_recursive,
                                              int maxlevel,
                                              p6est_refine_column_t refine_fn,
                                              p6est_init_t init_fn,
                                              p6est_replace_t replace_fn);

/** Vertically refine a forest with a bounded refinement level and a replace option.
 *
 * \param [in,out] p6est The forest is changed in place.
 * \param [in] refine_recursive Boolean to decide on recursive refinement.
 * \param [in] maxlevel   Maximum allowed refinement level (inclusive).
 *                        If this is negative the level is restricted only
 *                        by the compile-time constant QMAXLEVEL in p4est.h.
 * \param [in] refine_fn  Callback function that must return true if a quadrant
 *                        shall be refined.  If refine_recursive is true,
 *                        refine_fn is called for every existing and newly
 *                        created quadrant.  Otherwise, it is called for every
 *                        existing quadrant.  It is possible that a refinement
 *                        request made by the callback is ignored.  To catch
 *                        this case, you can examine whether init_fn or
 *                        replace_fn gets called.
 * \param [in] init_fn    Callback function to initialize the user_data for
 *                        newly created quadrants, which is guaranteed to be
 *                        allocated.  This function pointer may be NULL.
 * \param [in] replace_fn Callback function that allows the user to change
 *                        incoming quadrants based on the quadrants they
 *                        replace; may be NULL.
 */
void                p6est_refine_layers_ext (p6est_t * p6est,
                                             int refine_recursive,
                                             int maxlevel,
                                             p6est_refine_layer_t refine_fn,
                                             p6est_init_t init_fn,
                                             p6est_replace_t replace_fn);

/** Horizontally coarsen a forest.
 * \param [in,out] p6est The forest is changed in place.
 * \param [in] coarsen_recursive Boolean to decide on recursive coarsening.
 * \param [in] callback_orphans Boolean to enable calling coarsen_fn even on
 *                        non-families.  In this case, the second quadrant
 *                        pointer in the argument list of the callback is NULL,
 *                        subsequent pointers are undefined, and the return
 *                        value is ignored.  If coarsen_recursive is true, it
 *                        is possible that a quadrant is called once or more as
 *                        an orphan and eventually becomes part of a family.
 * \param [in] coarsen_fn Callback function that returns true if a
 *                        family of quadrants shall be coarsened.
 * \param [in] init_fn    Callback function to initialize the user_data
 *                        which is already allocated automatically.
 * \param [in] replace_fn Callback function that allows the user to change
 *                        incoming quadrants based on the quadrants they
 *                        replace.
 */
void                p6est_coarsen_columns_ext (p6est_t * p6est,
                                               int coarsen_recursive,
                                               int callback_orphans,
                                               p6est_coarsen_column_t
                                               coarsen_fn,
                                               p6est_init_t init_fn,
                                               p6est_replace_t replace_fn);

/** Vertically coarsen a forest.
 * \param [in,out] p6est The forest is changed in place.
 * \param [in] coarsen_recursive Boolean to decide on recursive coarsening.
 * \param [in] callback_orphans Boolean to enable calling coarsen_fn even on
 *                        non-families.  In this case, the second quadrant
 *                        pointer in the argument list of the callback is NULL,
 *                        subsequent pointers are undefined, and the return
 *                        value is ignored.  If coarsen_recursive is true, it
 *                        is possible that a quadrant is called once or more as
 *                        an orphan and eventually becomes part of a family.
 * \param [in] coarsen_fn Callback function that returns true if a
 *                        family of quadrants shall be coarsened.
 * \param [in] init_fn    Callback function to initialize the user_data
 *                        which is already allocated automatically.
 * \param [in] replace_fn Callback function that allows the user to change
 *                        incoming quadrants based on the quadrants they
 *                        replace.
 */
void                p6est_coarsen_layers_ext (p6est_t * p6est,
                                              int coarsen_recursive,
                                              int callback_orphans,
                                              p6est_coarsen_layer_t
                                              coarsen_fn,
                                              p6est_init_t init_fn,
                                              p6est_replace_t replace_fn);

/** Repartition the forest.
 *
 * The forest is partitioned between processors such that each processor
 * has an approximately equal number of quadrants (or weight).
 *
 * \param [in,out] p6est      The forest that will be partitioned.
 * \param [in]     partition_for_coarsening     If true, the partition
 *                            is modified to allow one level of coarsening.
 * \param [in]     weight_fn  A weighting function or NULL
 *                            for uniform partitioning.
 * \return         The global number of shipped quadrants
 */
p4est_gloidx_t      p6est_partition_ext (p6est_t * p6est,
                                         int partition_for_coarsening,
                                         p6est_weight_t weight_fn);

/** 2:1 balance the size differences of neighboring elements in a forest.
 * \param [in,out] p6est  The p6est to be worked on.
 * \param [in] btype      Balance type (face or corner/full).
 *                        Corner balance is almost never required when
 *                        discretizing a PDE; just causes smoother mesh grading.
 * \param [in] max_diff   The maximum difference between the horizontal
 *                        refinement level and the vertical refinement level
 * \param [in] min_diff   The minimum difference between the horizontal
 *                        refinement level and the vertical refinement level
 * \param [in] init_fn    Callback function to initialize the user_data
 *                        which is already allocated automatically.
 * \param [in] replace_fn Callback function that allows the user to change
 *                        incoming quadrants based on the quadrants they
 *                        replace.
 */
void                p6est_balance_ext (p6est_t * p6est,
                                       p8est_connect_type_t btype,
                                       int max_diff, int min_diff,
                                       p6est_init_t init_fn,
                                       p6est_replace_t replace_fn);

SC_EXTERN_C_END;

#endif
