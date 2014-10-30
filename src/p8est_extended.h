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

/** \file p8est_extended.h
 *
 * Interface routines with extended capabilities.
 *
 * \ingroup p8est
 */

#ifndef P8EST_EXTENDED_H
#define P8EST_EXTENDED_H

#include <p8est.h>
#include <p8est_mesh.h>
#include <p8est_iterate.h>

SC_EXTERN_C_BEGIN;

/* Data pertaining to selecting, inspecting, and profiling algorithms.
 * A pointer to this structure is hooked into the p8est main structure.
 *
 * TODO: Describe the purpose of various switches, counters, and timings.
 *
 * The balance_ranges and balance_notify* times are collected
 * whenever an inspect structure is present in p8est.
 */
struct p8est_inspect
{
  /** Use sc_ranges to determine the asymmetric communication pattern.
   * If \a use_balance_ranges is false (the default), sc_notify is used. */
  int                 use_balance_ranges;
  /** If true, call both sc_ranges and sc_notify and verify consistency.
   * Which is actually used is still determined by \a use_balance_ranges. */
  int                 use_balance_ranges_notify;
  /** Verify sc_ranges and/or sc_notify as applicable. */
  int                 use_balance_verify;
  /** If positive and smaller than p8est_num ranges, overrides it */
  int                 balance_max_ranges;
  size_t              balance_A_count_in;
  size_t              balance_A_count_out;
  size_t              balance_comm_sent;
  size_t              balance_comm_nzpeers;
  size_t              balance_B_count_in;
  size_t              balance_B_count_out;
  size_t              balance_zero_sends[2], balance_zero_receives[2];
  double              balance_A;
  double              balance_comm;
  double              balance_B;
  double              balance_ranges;   /**< time spent in sc_ranges */
  double              balance_notify;   /**< time spent in sc_notify */
  /** time spent in sc_notify_allgather */
  double              balance_notify_allgather;
  int                 use_B;
};

/** Callback function prototype to replace one set of quadrants with another.
 *
 * This is used by extended routines when the quadrants of an existing, valid
 * p8est are changed.  The callback allows the user to make changes to newly
 * initialized quadrants before the quadrants that they replace are destroyed.
 *
 * \param [in] num_outgoing The number of outgoing quadrants.
 * \param [in] outgoing     The outgoing quadrants: after the callback, the
 *                          user_data, if \a p8est->data_size is nonzero,
 *                          will be destroyed.
 * \param [in] num_incoming The number of incoming quadrants.
 * \param [in,out] incoming The incoming quadrants: prior to the callback,
 *                          the user_data, if \a p8est->data_size is nonzero,
 *                          is allocated, and the p8est_init_t callback,
 *                          if it has been provided, will be called.
 *
 * If the mesh is being refined, num_outgoing will be 1 and num_incoming will
 * be 8, and vice versa if the mesh is being coarsened.
 */
typedef void        (*p8est_replace_t) (p8est_t * p8est,
                                        p4est_topidx_t which_tree,
                                        int num_outgoing,
                                        p8est_quadrant_t * outgoing[],
                                        int num_incoming,
                                        p8est_quadrant_t * incoming[]);

/** Create a new forest.
 * This is a more general form of p8est_new.
 * See the documentation of p8est_new for basic usage.
 *
 * \param [in] min_quadrants    Minimum initial quadrants per processor.
 *                              Makes the refinement pattern mpisize-specific.
 * \param [in] min_level        The forest is refined at least to this level.
 *                              May be negative or 0, then it has no effect.
 * \param [in] fill_uniform     If true, fill the forest with a uniform mesh
 *                              instead of the coarsest possible one.
 *                              The latter is partition-specific so that
 *                              is usually not a good idea.
 */
p8est_t            *p8est_new_ext (sc_MPI_Comm mpicomm,
                                   p8est_connectivity_t * connectivity,
                                   p4est_locidx_t min_quadrants,
                                   int min_level, int fill_uniform,
                                   size_t data_size, p8est_init_t init_fn,
                                   void *user_pointer);

/** Create a new mesh.
 * \param [in] p8est                A forest that is fully 2:1 balanced.
 * \param [in] ghost                The ghost layer created from the
 *                                  provided p4est.
 * \param [in] compute_tree_index   Boolean to decide whether to allocate and
 *                                  compute the quad_to_tree list.
 * \param [in] compute_level_lists  Boolean to decide whether to compute the
 *                                  level lists in quad_level.
 * \param [in] btype                Currently ignored, only face neighbors
 *                                  are stored.
 * \return                          A fully allocated mesh structure.
 */
p8est_mesh_t       *p8est_mesh_new_ext (p8est_t * p4est,
                                        p8est_ghost_t * ghost,
                                        int compute_tree_index,
                                        int compute_level_lists,
                                        p8est_connect_type_t btype);

/** Refine a forest with a bounded refinement level and a replace option.
 * \param [in,out] p8est The forest is changed in place.
 * \param [in] refine_recursive Boolean to decide on recursive refinement.
 * \param [in] maxlevel   Maximum allowed refinement level (inclusive).
 *                        If this is negative the level is restricted only
 *                        by the compile-time constant QMAXLEVEL in p8est.h.
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
void                p8est_refine_ext (p8est_t * p8est,
                                      int refine_recursive, int maxlevel,
                                      p8est_refine_t refine_fn,
                                      p8est_init_t init_fn,
                                      p8est_replace_t replace_fn);

/** Coarsen a forest.
 * \param [in,out] p8est The forest is changed in place.
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
void                p8est_coarsen_ext (p8est_t * p8est, int coarsen_recursive,
                                       int callback_orphans,
                                       p8est_coarsen_t coarsen_fn,
                                       p8est_init_t init_fn,
                                       p8est_replace_t replace_fn);

/** 2:1 balance the size differences of neighboring elements in a forest.
 * \param [in,out] p8est  The p8est to be worked on.
 * \param [in] btype      Balance type (face, edge, or corner/full).
 *                        Corner balance is almost never required when
 *                        discretizing a PDE; just causes smoother mesh grading.
 * \param [in] init_fn    Callback function to initialize the user_data
 *                        which is already allocated automatically.
 * \param [in] replace_fn Callback function that allows the user to change
 *                        incoming quadrants based on the quadrants they
 *                        replace.
 */
void                p8est_balance_ext (p8est_t * p8est,
                                       p8est_connect_type_t btype,
                                       p8est_init_t init_fn,
                                       p8est_replace_t replace_fn);

void                p8est_balance_subtree_ext (p8est_t * p8est,
                                               p8est_connect_type_t btype,
                                               p4est_topidx_t which_tree,
                                               p8est_init_t init_fn,
                                               p8est_replace_t replace_fn);

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

/** p8est_iterate_ext adds the option \a remote: if this is false, then it is
 * the same as p8est_iterate; if this is true, then corner/edge callbacks are
 * also called on corners/edges for hanging faces/edges touched by local
 * quadrants.
 */
void                p8est_iterate_ext (p8est_t * p8est,
                                       p8est_ghost_t * ghost_layer,
                                       void *user_data,
                                       p8est_iter_volume_t iter_volume,
                                       p8est_iter_face_t iter_face,
                                       p8est_iter_edge_t iter_edge,
                                       p8est_iter_corner_t iter_corner,
                                       int remote);

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
p8est_t            *p8est_load_ext (const char *filename, sc_MPI_Comm mpicomm,
                                    size_t data_size, int load_data,
                                    int autopartition, int broadcasthead,
                                    void *user_pointer,
                                    p8est_connectivity_t ** connectivity);

/** The same as p8est_load_ext, but reading the connectivity/p8est from an
 * open sc_io_source_t stream.
 */
p8est_t            *p8est_source_ext (sc_io_source_t * src,
                                      sc_MPI_Comm mpicomm, size_t data_size,
                                      int load_data, int autopartition,
                                      int broadcasthead, void *user_pointer,
                                      p8est_connectivity_t ** connectivity);

SC_EXTERN_C_END;

#endif /* !P8EST_EXTENDED_H */
