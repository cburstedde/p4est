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

#ifndef P8EST_H
#define P8EST_H

/* p8est_connectivity.h includes p4est_base.h sc_containers.h */
#include <p8est_connectivity.h>

SC_EXTERN_C_BEGIN;

/* finest level of the octree for representing nodes */
#define P8EST_MAXLEVEL 19

/* finest level of the octree for representing octants */
#define P8EST_QMAXLEVEL 18

/* the length of a root quadrant */
#define P8EST_ROOT_LEN ((p4est_qcoord_t) 1 << P8EST_MAXLEVEL)

/* the length of a quadrant of level l */
#define P8EST_QUADRANT_LEN(l) ((p4est_qcoord_t) 1 << (P8EST_MAXLEVEL - (l)))

/* the offset of the highest quadrant at level l */
#define P8EST_LAST_OFFSET(l) (P8EST_ROOT_LEN - P8EST_QUADRANT_LEN (l))

typedef struct p8est_quadrant
{
  p4est_qcoord_t      x, y, z;
  int8_t              level, pad8;
  int16_t             pad16;
  union p8est_quadrant_data
  {
    void               *user_data;      /* never changed by p4est */
    long                user_long;      /* never changed by p4est */
    int                 user_int;       /* never changed by p4est */
    p4est_topidx_t      which_tree;
    struct
    {
      p4est_topidx_t      which_tree;
      int                 owner_rank;
    }
    piggy1;
    struct
    {
      p4est_topidx_t      which_tree;
      p4est_topidx_t      from_tree;
    }
    piggy2;
    struct
    {
      p4est_topidx_t      which_tree;
      p4est_locidx_t      local_num;
    }
    piggy3;
  }
  p;
}
p8est_quadrant_t;

typedef struct p8est_tree
{
  sc_array_t          quadrants;        /* locally stored quadrants */
  p8est_quadrant_t    first_desc, last_desc;    /* first and last descendant */
  p4est_locidx_t      quadrants_offset; /* cumulative sum over earlier trees */
  p4est_locidx_t      quadrants_per_level[P8EST_MAXLEVEL + 1];  /* locals only */
  int8_t              maxlevel; /* highest local quadrant level */
}
p8est_tree_t;

/* Data pertaining to selecting, inspecting, and profiling algorithms.
 * A pointer to this structure is hooked into the p8est main structure.
 *
 * TODO: Describe the purpose of various switches, counters, and timings.
 *
 * The balance_ranges and balance_notify* times are collected
 * whenever an inspect structure is present in p8est.
 */
typedef struct p8est_inspect
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
}
p8est_inspect_t;

typedef struct p8est
{
  sc_MPI_Comm         mpicomm;
  int                 mpisize, mpirank;

  size_t              data_size;        /* size of per-quadrant user_data */
  void               *user_pointer;     /* convenience pointer for users,
                                           will never be touched by p8est */

  p4est_topidx_t      first_local_tree; /* 0-based index of first local tree,
                                           must be -1 for an empty processor */
  p4est_topidx_t      last_local_tree;  /* 0-based index of last local tree,
                                           must be -2 for an empty processor */
  p4est_locidx_t      local_num_quadrants;      /* number of quadrants on all
                                                   trees on this processor */
  p4est_gloidx_t      global_num_quadrants;     /* number of quadrants on all
                                                   trees on all processors */
  p4est_gloidx_t     *global_first_quadrant;    /* first global quadrant index
                                                   for each proc and 1 beyond
                                                 */
  p8est_quadrant_t   *global_first_position;    /* first smallest possible quad
                                                   for each proc and 1 beyond
                                                 */
  p8est_connectivity_t *connectivity;   /* connectivity structure, not owned */
  sc_array_t         *trees;    /* list of all trees */

  sc_mempool_t       *user_data_pool;   /* memory allocator for user data
                                         * WARNING: This is NULL if data size
                                         *          equals zero.
                                         */
  sc_mempool_t       *quadrant_pool;    /* memory allocator
                                           for temporary quadrants */
  p8est_inspect_t    *inspect;  /* algorithmic switches */
}
p8est_t;

/** Calculate memory usage of a forest structure.
 * The connectivity structure is not counted since it is not owned;
 * use p8est_connectivity_memory_usage (p8est->connectivity).
 * \param [in] p8est    Forest structure.
 * \return              Memory used in bytes.
 */
size_t              p8est_memory_used (p8est_t * p8est);

/** Callback function prototype to initialize the quadrant's user data.
 */
typedef void        (*p8est_init_t) (p8est_t * p8est,
                                     p4est_topidx_t which_tree,
                                     p8est_quadrant_t * quadrant);

/** Callback function prototype to decide for refinement.
 * \return nonzero if the quadrant shall be refined.
 */
typedef int         (*p8est_refine_t) (p8est_t * p8est,
                                       p4est_topidx_t which_tree,
                                       p8est_quadrant_t * quadrant);

/** Callback function prototype to decide for coarsening.
 * \param [in] quadrants   Pointers to 8 siblings in Morton ordering.
 * \return nonzero if the quadrants shall be replaced with their parent.
 */
typedef int         (*p8est_coarsen_t) (p8est_t * p8est,
                                        p4est_topidx_t which_tree,
                                        p8est_quadrant_t * quadrants[]);

/** Callback function prototype to calculate weights for partitioning.
 * \return a 32bit integer >= 0 as the quadrant weight.
 * \note    Global sum of weights must fit into a 64bit integer.
 */
typedef int         (*p8est_weight_t) (p8est_t * p8est,
                                       p4est_topidx_t which_tree,
                                       p8est_quadrant_t * quadrant);

extern void        *P8EST_DATA_UNINITIALIZED;
extern const int    p8est_num_ranges;

/** set statically allocated quadrant to defined values */
#define P8EST_QUADRANT_INIT(q) \
  ((void) memset ((q), -1, sizeof (p8est_quadrant_t)))

/** Transform a quadrant coordinate into the space spanned by tree vertices.
 * \param [in] connectivity     Connectivity must provide the vertices.
 * \param [in] treeid           Identify the tree that contains x, y, z.
 * \param [in] x, y, z          Quadrant coordinates relative to treeid.
 * \param [out] vxyz            Transformed coordinates in vertex space.
 */
void                p8est_qcoord_to_vertex (p8est_connectivity_t *
                                            connectivity,
                                            p4est_topidx_t treeid,
                                            p4est_qcoord_t x,
                                            p4est_qcoord_t y,
                                            p4est_qcoord_t z, double vxyz[3]);

/** Create a new forest.
 * The new forest consists of equi-partitioned root quadrants.
 * When there are more processors than trees, some processors are empty.
 *
 * \param [in] mpicomm       A valid MPI communicator.
 * \param [in] connectivity  This is the connectivity information that
 *                           the forest is built with.  Note the p8est
 *                           does not take ownership of the memory.
 * \param [in] data_size     This is the size of data for each quadrant which
 *                           can be zero.  Then user_data_pool is set to NULL.
 * \param [in] init_fn       Callback function to initialize the user_data
 *                           which is already allocated automatically.
 * \param [in] user_pointer  Assign to the user_pointer member of the p8est
 *                           before init_fn is called the first time.
 *
 * \return This returns a valid forest.
 *
 * \note The connectivity structure must not be destroyed
 *       during the lifetime of this forest.
 */
p8est_t            *p8est_new (sc_MPI_Comm mpicomm,
                               p8est_connectivity_t * connectivity,
                               size_t data_size,
                               p8est_init_t init_fn, void *user_pointer);

/** Destroy a p8est.
 *
 * \note The connectivity structure is not destroyed with the p8est.
 */
void                p8est_destroy (p8est_t * p8est);

/** Make a deep copy of a p8est.
 * The connectivity is not duplicated.
 * Copying of quadrant user data is optional.
 * If old and new data sizes are 0, the user_data field is copied regardless.
 *
 * \param [in]  copy_data  If true, data are copied.
 *                         If false, data_size is set to 0.
 * \return  Returns a valid p8est that does not depend on the input.
 */
p8est_t            *p8est_copy (p8est_t * input, int copy_data);

/** Reset user pointer and element data.
 * When the data size is changed the quadrant data is freed and allocated.
 * The initialization callback is invoked on each quadrant.
 * Old user_data content is disregarded.
 *
 * \param [in] data_size     This is the size of data for each quadrant which
 *                           can be zero.  Then user_data_pool is set to NULL.
 * \param [in] init_fn       Callback function to initialize the user_data
 *                           which is already allocated automatically.
 *                           May be NULL.
 * \param [in] user_pointer  Assign to the user_pointer member of the p8est
 *                           before init_fn is called the first time.
 */
void                p8est_reset_data (p8est_t * p8est, size_t data_size,
                                      p8est_init_t init_fn,
                                      void *user_pointer);

/** Refine a forest.
 * \param [in,out] p8est The forest is changed in place.
 * \param [in] refine_recursive Boolean to decide on recursive refinement.
 * \param [in] refine_fn Callback function that must return true if a quadrant
 *                       shall be refined.  If refine_recursive is true,
 *                       refine_fn is called for every existing and newly
 *                       created quadrant.  Otherwise, it is called for every
 *                       existing quadrant.  It is possible that a refinement
 *                       request made by the callback is ignored.  To catch
 *                       this case, you can examine whether init_fn gets
 *                       called, or use p8est_refine_ext in p8est_extended.h
 *                       and examine whether replace_fn gets called.
 * \param [in] init_fn   Callback function to initialize the user_data of newly
 *                       created quadrants, which is already allocated.  This
 *                       function pointer may be NULL.
 */
void                p8est_refine (p8est_t * p8est,
                                  int refine_recursive,
                                  p8est_refine_t refine_fn,
                                  p8est_init_t init_fn);

/** Coarsen a forest.
 * \param [in,out] p8est  The forest is changed in place.
 * \param [in] coarsen_recursive Boolean to decide on recursive coarsening.
 * \param [in] coarsen_fn Callback function that returns true if a
 *                        family of quadrants shall be coarsened
 * \param [in] init_fn    Callback function to initialize the user_data
 *                        which is already allocated automatically.
 */
void                p8est_coarsen (p8est_t * p8est,
                                   int coarsen_recursive,
                                   p8est_coarsen_t coarsen_fn,
                                   p8est_init_t init_fn);

/** 2:1 balance the size differences of neighboring elements in a forest.
 * \param [in] p8est     The p8est to be worked on.
 * \param [in] btype     Balance type (face, edge, or corner/full).  Examples:
 *                       Finite volume or discontinous Galerkin methods only
 *                       require face balance.  Continuous finite element
 *                       methods usually require edge balance.  Corner balance
 *                       is almost never required mathematically; it just
 *                       produces a smoother mesh grading.
 * \param [in] init_fn   Callback function to initialize the user_data
 *                       which is already allocated automatically.
 */
void                p8est_balance (p8est_t * p8est,
                                   p8est_connect_type_t btype,
                                   p8est_init_t init_fn);

/** Equally partition the forest.
 * The partition can be by element count or by a user-defined weight.
 *
 * The forest will be partitioned between processors such that they
 * have an approximately equal number of quadrants (or sum of weights).
 *
 * \param [in,out] p8est      The forest that will be partitioned.
 * \param [in]     allow_for_coarsening Slightly modify partition such that
 *                            quadrant families are not split between ranks.
 * \param [in]     weight_fn  A weighting function or NULL
 *                            for uniform partitioning.
 */
void                p8est_partition (p8est_t * p8est,
                                     int allow_for_coarsening,
                                     p8est_weight_t weight_fn);

/** Compute the checksum for a forest.
 * Based on quadrant arrays only. It is independent of partition and mpisize.
 * \return  Returns the checksum on processor 0 only. 0 on other processors.
 */
unsigned            p8est_checksum (p8est_t * p8est);

/** Save the complete connectivity/p4est data to disk.  This is a collective
 * operation that all MPI processes need to call.  All processes write
 * into the same file, so the filename given needs to be identical over
 * all parallel invocations.
 *
 * \param [in] filename    Name of the file to write.
 * \param [in] p8est       Valid forest structure.
 * \param [in] save_data   If true, the element data is saved.
 *                         Otherwise, a data size of 0 is saved.
 * \note            Aborts on file errors.
 * \note            If p4est is not configured to use MPI-IO, some processes
 *                  return from this function before the file is complete, in
 *                  which case immediate read-access to the file may require a
 *                  call to sc_MPI_Barrier.
 */
void                p8est_save (const char *filename, p8est_t * p8est,
                                int save_data);

/** Load the complete connectivity/p4est structure from disk.
 * \param [in] filename         Name of the file to read.
 * \param [in] mpicomm          A valid MPI communicator.
 * \param [in] data_size        Size of data for each quadrant which can be
 *                              zero.  Then user_data_pool is set to NULL.
 *                              If data_size is zero, load_data is ignored.
 * \param [in] load_data        If true, the element data is loaded.  This is
 *                              only permitted if the saved data size matches.
 *                              If false, the stored data size is ignored.
 * \param [in] user_pointer     Assign to the user_pointer member of the p4est
 *                              before init_fn is called the first time.
 * \param [out] connectivity    Connectivity must be destroyed separately.
 * \return          Returns a valid forest structure. A pointer to a valid
 *                  connectivity structure is returned through the last
 *                  argument.
 * \note            Aborts on file errors or invalid file contents.
 */
p8est_t            *p8est_load (const char *filename, sc_MPI_Comm mpicomm,
                                size_t data_size, int load_data,
                                void *user_pointer,
                                p8est_connectivity_t ** connectivity);

/** Return a pointer to an array element indexed by a p4est_topidx_t.
 * \param [in] index needs to be in [0]..[elem_count-1].
 */
/*@unused@*/
static inline p8est_tree_t *
p8est_tree_array_index (sc_array_t * array, p4est_topidx_t it)
{
  P4EST_ASSERT (array->elem_size == sizeof (p8est_tree_t));
  P4EST_ASSERT (it >= 0 && (size_t) it < array->elem_count);

  return (p8est_tree_t *) (array->array +
                           sizeof (p8est_tree_t) * (size_t) it);
}

/** Return a pointer to a quadrant array element indexed by a size_t. */
/*@unused@*/
static inline p8est_quadrant_t *
p8est_quadrant_array_index (sc_array_t * array, size_t it)
{
  P4EST_ASSERT (array->elem_size == sizeof (p8est_quadrant_t));
  P4EST_ASSERT (it < array->elem_count);

  return (p8est_quadrant_t *) (array->array + sizeof (p8est_quadrant_t) * it);
}

/** Call sc_array_push for a quadrant array. */
/*@unused@*/
static inline p8est_quadrant_t *
p8est_quadrant_array_push (sc_array_t * array)
{
  P4EST_ASSERT (array->elem_size == sizeof (p8est_quadrant_t));

  return (p8est_quadrant_t *) sc_array_push (array);
}

/** Call sc_mempool_alloc for a mempool creating quadrants. */
/*@unused@*/
static inline p8est_quadrant_t *
p8est_quadrant_mempool_alloc (sc_mempool_t * mempool)
{
  P4EST_ASSERT (mempool->elem_size == sizeof (p8est_quadrant_t));

  return (p8est_quadrant_t *) sc_mempool_alloc (mempool);
}

/** Call sc_list pop for a quadrant array. */
/*@unused@*/
static inline p8est_quadrant_t *
p8est_quadrant_list_pop (sc_list_t * list)
{
  return (p8est_quadrant_t *) sc_list_pop (list);
}

SC_EXTERN_C_END;

#endif /* !P8EST_H */
