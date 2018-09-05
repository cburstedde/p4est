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

/********************************************************************
 *                          IMPORTANT NOTE                          *
 *                                                                  *
 * These interfaces are intended for those who like finer control.  *
 * The API offers extended versions of some basic p4est functions.  *
 * The API may change without notice.                               *
 ********************************************************************/

/** \file p4est_extended.h
 *
 * Interface routines with extended capabilities.
 *
 * \ingroup p4est
 */

#ifndef P4EST_EXTENDED_H
#define P4EST_EXTENDED_H

#include <p4est.h>
#include <p4est_mesh.h>
#include <p4est_iterate.h>
#include <p4est_lnodes.h>

SC_EXTERN_C_BEGIN;

/** Data pertaining to selecting, inspecting, and profiling algorithms.
 * A pointer to this structure is hooked into the p4est main structure.
 *
 *
 * The balance_ranges and balance_notify* times are collected
 * whenever an inspect structure is present in p4est.
 */
/* TODO: Describe the purpose of various switches, counters, and timings. */
struct p4est_inspect
{
  /** Use sc_ranges to determine the asymmetric communication pattern.
   * If \a use_balance_ranges is false (the default), sc_notify is used. */
  int                 use_balance_ranges;
  /** If true, call both sc_ranges and sc_notify and verify consistency.
   * Which is actually used is still determined by \a use_balance_ranges. */
  int                 use_balance_ranges_notify;
  /** Verify sc_ranges and/or sc_notify as applicable. */
  int                 use_balance_verify;
  /** If positive and smaller than p4est_num ranges, overrides it */
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
 * p4est are changed.  The callback allows the user to make changes to newly
 * initialized quadrants before the quadrants that they replace are destroyed.
 *
 * \param [in] num_outgoing The number of outgoing quadrants.
 * \param [in] outgoing     The outgoing quadrants: after the callback, the
 *                          user_data, if \a p4est->data_size is nonzero,
 *                          will be destroyed.
 * \param [in] num_incoming The number of incoming quadrants.
 * \param [in,out] incoming The incoming quadrants: prior to the callback,
 *                          the user_data, if \a p4est->data_size is nonzero,
 *                          is allocated, and the p4est_init_t callback,
 *                          if it has been provided, will be called.
 *
 * If the mesh is being refined, num_outgoing will be 1 and num_incoming will
 * be 4, and vice versa if the mesh is being coarsened.
 */
typedef void        (*p4est_replace_t) (p4est_t * p4est,
                                        p4est_topidx_t which_tree,
                                        int num_outgoing,
                                        p4est_quadrant_t * outgoing[],
                                        int num_incoming,
                                        p4est_quadrant_t * incoming[]);

/** Create a new forest.
 * This is a more general form of p4est_new.
 * See the documentation of p4est_new for basic usage.
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
p4est_t            *p4est_new_ext (sc_MPI_Comm mpicomm,
                                   p4est_connectivity_t * connectivity,
                                   p4est_locidx_t min_quadrants,
                                   int min_level, int fill_uniform,
                                   size_t data_size, p4est_init_t init_fn,
                                   void *user_pointer);

/** Create a new mesh.
 * \param [in] p4est                A forest that is fully 2:1 balanced.
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
p4est_mesh_t       *p4est_mesh_new_ext (p4est_t * p4est,
                                        p4est_ghost_t * ghost,
                                        int compute_tree_index,
                                        int compute_level_lists,
                                        p4est_connect_type_t btype);

/** Make a deep copy of a p4est.
 * The connectivity is not duplicated.
 * Copying of quadrant user data is optional.
 * If old and new data sizes are 0, the user_data field is copied regardless.
 * The inspect member of the copy is set to NULL.
 * The revision counter of the copy is set to zero.
 *
 * \param [in]  copy_data  If true, data are copied.
 *                         If false, data_size is set to 0.
 * \param [in]  duplicate_mpicomm  If true, MPI communicator is copied.
 * \return  Returns a valid p4est that does not depend on the input,
 *                         except for borrowing the same connectivity.
 *                         Its revision counter is 0.
 */
p4est_t            *p4est_copy_ext (p4est_t * input, int copy_data,
                                    int duplicate_mpicomm);

/** Refine a forest with a bounded refinement level and a replace option.
 * \param [in,out] p4est The forest is changed in place.
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
void                p4est_refine_ext (p4est_t * p4est,
                                      int refine_recursive, int maxlevel,
                                      p4est_refine_t refine_fn,
                                      p4est_init_t init_fn,
                                      p4est_replace_t replace_fn);

/** Coarsen a forest.
 * \param [in,out] p4est The forest is changed in place.
 * \param [in] coarsen_recursive Boolean to decide on recursive coarsening.
 * \param [in] callback_orphans Boolean to enable calling coarsen_fn even on
 *                        non-families.  In this case, the second quadrant
 *                        pointer in the argument list of the callback is NULL,
 *                        subsequent pointers are undefined, and the return
 *                        value is ignored.  If coarsen_recursive is true, it
 *                        is possible that a quadrant is called once or more as
 *                        an orphan and eventually becomes part of a family.
 *                        With coarsen_recursive false and callback_orphans true,
 *                        it is guaranteed that every quadrant is passed exactly
 *                        once into the coarsen_fn callback.
 * \param [in] coarsen_fn Callback function that returns true if a
 *                        family of quadrants shall be coarsened.
 * \param [in] init_fn    Callback function to initialize the user_data
 *                        which is already allocated automatically.
 * \param [in] replace_fn Callback function that allows the user to change
 *                        incoming quadrants based on the quadrants they
 *                        replace.
 */
void                p4est_coarsen_ext (p4est_t * p4est, int coarsen_recursive,
                                       int callback_orphans,
                                       p4est_coarsen_t coarsen_fn,
                                       p4est_init_t init_fn,
                                       p4est_replace_t replace_fn);

/** 2:1 balance the size differences of neighboring elements in a forest.
 * \param [in,out] p4est  The p4est to be worked on.
 * \param [in] btype      Balance type (face or corner/full).
 *                        Corner balance is almost never required when
 *                        discretizing a PDE; just causes smoother mesh grading.
 * \param [in] init_fn    Callback function to initialize the user_data
 *                        which is already allocated automatically.
 * \param [in] replace_fn Callback function that allows the user to change
 *                        incoming quadrants based on the quadrants they
 *                        replace.
 */
void                p4est_balance_ext (p4est_t * p4est,
                                       p4est_connect_type_t btype,
                                       p4est_init_t init_fn,
                                       p4est_replace_t replace_fn);

void                p4est_balance_subtree_ext (p4est_t * p4est,
                                               p4est_connect_type_t btype,
                                               p4est_topidx_t which_tree,
                                               p4est_init_t init_fn,
                                               p4est_replace_t replace_fn);

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
 * \return         The global number of shipped quadrants
 */
p4est_gloidx_t      p4est_partition_ext (p4est_t * p4est,
                                         int partition_for_coarsening,
                                         p4est_weight_t weight_fn);

/** Correct partition to allow one level of coarsening.
 *
 * \param [in] p4est                     forest whose partition is corrected
 * \param [in,out] num_quadrants_in_proc partition that will be corrected
 * \return                               absolute number of moved quadrants
 */
p4est_gloidx_t      p4est_partition_for_coarsening (p4est_t * p4est,
                                                    p4est_locidx_t *
                                                    num_quadrants_in_proc);

/** p4est_iterate_ext adds the option \a remote: if this is false, then it is
 * the same as p4est_iterate; if this is true, then corner callbacks are also
 * called on corners for hanging faces touched by local quadrants.
 */
void                p4est_iterate_ext (p4est_t * p4est,
                                       p4est_ghost_t * ghost_layer,
                                       void *user_data,
                                       p4est_iter_volume_t iter_volume,
                                       p4est_iter_face_t iter_face,
                                       p4est_iter_corner_t iter_corner,
                                       int remote);

/** Save the complete connectivity/p4est data to disk.  This is a collective
 * operation that all MPI processes need to call.  All processes write
 * into the same file, so the filename given needs to be identical over
 * all parallel invocations.
 * See p4est_load_ext for information on the autopartition parameter.
 * \param [in] filename    Name of the file to write.
 * \param [in] p4est       Valid forest structure.
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
void                p4est_save_ext (const char *filename, p4est_t * p4est,
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
p4est_t            *p4est_load_ext (const char *filename, sc_MPI_Comm mpicomm,
                                    size_t data_size, int load_data,
                                    int autopartition, int broadcasthead,
                                    void *user_pointer,
                                    p4est_connectivity_t ** connectivity);

/** The same as p4est_load_ext, but reading the connectivity/p4est from an
 * open sc_io_source_t stream.
 */
p4est_t            *p4est_source_ext (sc_io_source_t * src,
                                      sc_MPI_Comm mpicomm, size_t data_size,
                                      int load_data, int autopartition,
                                      int broadcasthead, void *user_pointer,
                                      p4est_connectivity_t ** connectivity);

/** Create the data necessary to create a PETsc DMPLEX representation of a
 * forest, as well as the accompanying lnodes and ghost layer.  The forest
 * must be at least face balanced (see p4est_balance()).  See
 * test/test_plex2.c for example usage.
 *
 * All arrays should be initialized to hold sizeof (p4est_locidx_t), except
 * for \a out_remotes, which should be initialized to hold
 * (2 * sizeof (p4est_locidx_t)).
 *
 * \param[in]     p4est                 the forest
 * \param[out]    ghost                 the ghost layer
 * \param[out]    lnodes                the lnodes
 * \param[in]     ctype                 the type of adjacency for the overlap
 * \param[in]     overlap               the number of layers of overlap (zero
 *                                      is acceptable)
 * \param[out]    first_local_quad      the local quadrants are assigned
 *                                      contiguous plex indices, starting with
 *                                      this index
 * \param[in,out] out_points_per_dim    filled with argument for
 *                                      DMPlexCreateFromDAG()
 * \param[in,out] out_cone_sizes        filled with argument for
 *                                      DMPlexCreateFromDAG()
 * \param[in,out] out_cones             filled with argument for
 *                                      DMPlexCreateFromDAG()
 * \param[in,out] out_cone_orientations filled with argument for
 *                                      DMPlexCreateFromDAG()
 * \param[in,out] out_vertex_coords     filled with argument for
 *                                      DMPlexCreateFromDAG()
 * \param[in,out] out_children          filled with argument for
 *                                      DMPlexSetTree()
 * \param[in,out] out_parents           filled with argument for
 *                                      DMPlexSetTree()
 * \param[in,out] out_childids          filled with argument for
 *                                      DMPlexSetTree()
 * \param[in,out] out_leaves            filled with argument for
 *                                      PetscSFSetGraph()
 * \param[in,out] out_remotes           filled with argument for
 *                                      PetscSFSetGraph()
 * \param[in]     custom_numbering      Whether or use the default numbering
 *                                      (0) of DMPlex child ids or the custom
 *                                      (1).
 */
void                p4est_get_plex_data_ext (p4est_t * p4est,
                                             p4est_ghost_t ** ghost,
                                             p4est_lnodes_t ** lnodes,
                                             p4est_connect_type_t ctype,
                                             int overlap,
                                             p4est_locidx_t *
                                             first_local_quad,
                                             sc_array_t * out_points_per_dim,
                                             sc_array_t * out_cone_sizes,
                                             sc_array_t * out_cones,
                                             sc_array_t *
                                             out_cone_orientations,
                                             sc_array_t * out_vertex_coords,
                                             sc_array_t * out_children,
                                             sc_array_t * out_parents,
                                             sc_array_t * out_childids,
                                             sc_array_t * out_leaves,
                                             sc_array_t * out_remotes,
                                             int custom_numbering);

SC_EXTERN_C_END;

#endif /* !P4EST_EXTENDED_H */
