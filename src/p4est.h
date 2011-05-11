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

#ifndef P4EST_H
#define P4EST_H

#ifdef P4EST_TO_P8EST_H
#error "The include files p4est.h and p4est_to_p8est.h cannot be combined"
#endif

/* p4est_connectivity.h includes p4est_base.h sc_containers.h */
#include <p4est_connectivity.h>

SC_EXTERN_C_BEGIN;

/* size of insulation layer */
#define P4EST_INSUL 9

/* finest level of the quadtree for representing nodes */
#define P4EST_MAXLEVEL 30

/* finest level of the quadtree for representing quadrants */
#define P4EST_QMAXLEVEL 29

/* the length of a root quadrant */
#define P4EST_ROOT_LEN ((p4est_qcoord_t) 1 << P4EST_MAXLEVEL)

/* the length of a quadrant of level l */
#define P4EST_QUADRANT_LEN(l) ((p4est_qcoord_t) 1 << (P4EST_MAXLEVEL - (l)))

/* the offset of the highest quadrant at level l */
#define P4EST_LAST_OFFSET(l) (P4EST_ROOT_LEN - P4EST_QUADRANT_LEN (l))

typedef struct p4est_quadrant
{
  p4est_qcoord_t      x, y;
  int8_t              level, pad8;
  int16_t             pad16;
  union p4est_quadrant_data
  {
    void               *user_data;
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
p4est_quadrant_t;

typedef struct p4est_tree
{
  sc_array_t          quadrants;        /* locally stored quadrants */
  p4est_quadrant_t    first_desc, last_desc;    /* first and last descendant */
  p4est_locidx_t      quadrants_offset; /* cumulative sum over earlier trees */
  p4est_locidx_t      quadrants_per_level[P4EST_MAXLEVEL + 1];  /* locals only */
  int8_t              maxlevel; /* highest local quadrant level */
}
p4est_tree_t;

typedef struct p4est
{
  MPI_Comm            mpicomm;
  int                 mpisize, mpirank;

  size_t              data_size;        /* size of per-quadrant user_data */
  void               *user_pointer;     /* convenience pointer for users,
                                           will never be touched by p4est */

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
  p4est_quadrant_t   *global_first_position;    /* first smallest possible quad
                                                   for each proc and 1 beyond
                                                 */
  p4est_connectivity_t *connectivity;   /* connectivity structure, not owned */
  sc_array_t         *trees;    /* list of all trees */

  sc_mempool_t       *user_data_pool;   /* memory allocator for user data
                                         * WARNING: This is NULL if data size
                                         *          equals zero.
                                         */
  sc_mempool_t       *quadrant_pool;    /* memory allocator
                                           for temporary quadrants */
}
p4est_t;

/** Calculate memory usage of a forest structure.
 * The connectivity structure is not counted since it is not owned;
 * use p4est_connectivity_memory_usage (p4est->connectivity).
 * \param [in] p4est    Forest structure.
 * \return              Memory used in bytes.
 */
size_t              p4est_memory_used (p4est_t * p4est);

/** Callback function prototype to initialize the quadrant's user data.
 */
typedef void        (*p4est_init_t) (p4est_t * p4est,
                                     p4est_topidx_t which_tree,
                                     p4est_quadrant_t * quadrant);

/** Callback function prototype to decide for refinement.
 * \return nonzero if the quadrant shall be refined.
 */
typedef int         (*p4est_refine_t) (p4est_t * p4est,
                                       p4est_topidx_t which_tree,
                                       p4est_quadrant_t * quadrant);

/** Callback function prototype to decide for coarsening.
 * \param [in] quadrants   Pointers to 4 siblings in Morton ordering.
 * \return nonzero if the quadrants shall be replaced with their parent.
 */
typedef int         (*p4est_coarsen_t) (p4est_t * p4est,
                                        p4est_topidx_t which_tree,
                                        p4est_quadrant_t * quadrants[]);

/** Callback function prototype to calculate weights for partitioning.
 * \return a 32bit integer >= 0 as the quadrant weight.
 * \note    (global sum of weights * mpisize) must fit into a 64bit integer.
 */
typedef int         (*p4est_weight_t) (p4est_t * p4est,
                                       p4est_topidx_t which_tree,
                                       p4est_quadrant_t * quadrant);

extern void        *P4EST_DATA_UNINITIALIZED;
extern const int    p4est_num_ranges;

/** set statically allocated quadrant to defined values */
#define P4EST_QUADRANT_INIT(q) \
  ((void) memset ((q), -1, sizeof (p4est_quadrant_t)))

/** Transform a quadrant coordinate into the space spanned by tree vertices.
 * \param [in] connectivity     Connectivity must provide the vertices.
 * \param [in] treeid           Identify the tree that contains x, y.
 * \param [in] x, y             Quadrant coordinates relative to treeid.
 * \param [out] vxy             Transformed coordinates in vertex space.
 */
void                p4est_qcoord_to_vertex (p4est_connectivity_t *
                                            connectivity,
                                            p4est_topidx_t treeid,
                                            p4est_qcoord_t x,
                                            p4est_qcoord_t y, double vxyz[3]);

/** Create a new forest.
 * The new forest consists of equi-partitioned root quadrants.
 * When there are more processors than trees, some processors are empty.
 *
 * \param [in] mpicomm       A valid MPI communicator.
 * \param [in] connectivity  This is the connectivity information that
 *                           the forest is built with.  Note the p4est
 *                           does not take ownership of the memory.
 * \param [in] data_size     This is the size of data for each quadrant which
 *                           can be zero.  Then user_data_pool is set to NULL.
 * \param [in] init_fn       Callback function to initialize the user_data
 *                           which is already allocated automatically.
 * \param [in] user_pointer  Assign to the user_pointer member of the p4est
 *                           before init_fn is called the first time.
 *
 * \return This returns a valid forest.
 *
 * \note The connectivity structure must not be destroyed
 *       during the lifetime of this forest.
 */
p4est_t            *p4est_new (MPI_Comm mpicomm,
                               p4est_connectivity_t * connectivity,
                               size_t data_size,
                               p4est_init_t init_fn, void *user_pointer);

/** Destroy a p4est.
 *
 * \note The connectivity structure is not destroyed with the p4est.
 */
void                p4est_destroy (p4est_t * p4est);

/** Make a deep copy of a p4est.
 * The connectivity is not duplicated.
 * Copying of quadrant user data is optional.
 *
 * \param [in]  copy_data  If true, data are copied.
 *                         If false, data_size is set to 0.
 * \return  Returns a valid p4est that does not depend on the input.
 */
p4est_t            *p4est_copy (p4est_t * input, int copy_data);

/** Reset user pointer and element data.
 * When the data size is changed the quadrant data is freed and allocated.
 * The initialization callback is invoked on each quadrant.
 *
 * \param [in] data_size     This is the size of data for each quadrant which
 *                           can be zero.  Then user_data_pool is set to NULL.
 * \param [in] init_fn       Callback function to initialize the user_data
 *                           which is already allocated automatically.
 * \param [in] user_pointer  Assign to the user_pointer member of the p4est
 *                           before init_fn is called the first time.
 */
void                p4est_reset_data (p4est_t * p4est, size_t data_size,
                                      p4est_init_t init_fn,
                                      void *user_pointer);

/** Refine a forest.
 * \param [in] refine_fn Callback function that returns true
 *                       if a quadrant shall be refined
 * \param [in] init_fn   Callback function to initialize the user_data
 *                       which is already allocated automatically.
 */
void                p4est_refine (p4est_t * p4est,
                                  int refine_recursive,
                                  p4est_refine_t refine_fn,
                                  p4est_init_t init_fn);

/** Coarsen a forest.
 * \param [in] coarsen_fn Callback function that returns true if a
 *                        family of quadrants shall be coarsened
 * \param [in] init_fn    Callback function to initialize the user_data
 *                        which is already allocated automatically.
 */
void                p4est_coarsen (p4est_t * p4est,
                                   int coarsen_recursive,
                                   p4est_coarsen_t coarsen_fn,
                                   p4est_init_t init_fn);

/** Balance a forest.
 * \param [in] p4est     The p4est to be worked on.
 * \param [in] btype     Balance type (face, corner or default, full).
 * \param [in] init_fn   Callback function to initialize the user_data
 *                       which is already allocated automatically.
 */
void                p4est_balance (p4est_t * p4est,
                                   p4est_connect_type_t btype,
                                   p4est_init_t init_fn);

/** Equally partition the forest.
 *
 * The forest will be partitioned between processors where they each
 * have an approximately equal number of quadrants.
 *
 * \param [in,out] p4est      The forest that will be partitioned.
 * \param [in]     weight_fn  A weighting function or NULL
 *                            for uniform partitioning.
 */
void                p4est_partition (p4est_t * p4est,
                                     p4est_weight_t weight_fn);

/** Compute the checksum for a forest.
 * Based on quadrant arrays only. It is independent of partition and mpisize.
 * \return  Returns the checksum on processor 0 only. 0 on other processors.
 */
unsigned            p4est_checksum (p4est_t * p4est);

/** Save the complete connectivity/p4est data to disk.  This is a collective
 * operation that all MPI processes need to call.  All processes write
 * into the same file, so the filename given needs to be identical over
 * all parallel invocations.
 * \param [in] filename    Name of the file to write.
 * \param [in] p4est       Valid forest structure.
 * \param [in] save_data   If true, the element data is saved.
 *                         Otherwise, a data size of 0 is saved.
 * \note            Aborts on file errors.
 */
void                p4est_save (const char *filename, p4est_t * p4est,
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
p4est_t            *p4est_load (const char *filename, MPI_Comm mpicomm,
                                size_t data_size, int load_data,
                                void *user_pointer,
                                p4est_connectivity_t ** connectivity);

/** Return a pointer to an array element indexed by a p4est_topidx_t.
 * \param [in] index needs to be in [0]..[elem_count-1].
 */
/*@unused@*/
static inline p4est_tree_t *
p4est_tree_array_index (sc_array_t * array, p4est_topidx_t it)
{
  P4EST_ASSERT (array->elem_size == sizeof (p4est_tree_t));
  P4EST_ASSERT (it >= 0 && (size_t) it < array->elem_count);

  return (p4est_tree_t *) (array->array +
                           sizeof (p4est_tree_t) * (size_t) it);
}

/** Return a pointer to a quadrant array element indexed by a size_t. */
/*@unused@*/
static inline p4est_quadrant_t *
p4est_quadrant_array_index (sc_array_t * array, size_t it)
{
  P4EST_ASSERT (array->elem_size == sizeof (p4est_quadrant_t));
  P4EST_ASSERT (it < array->elem_count);

  return (p4est_quadrant_t *) (array->array + sizeof (p4est_quadrant_t) * it);
}

/** Call sc_array_push for a quadrant array. */
/*@unused@*/
static inline p4est_quadrant_t *
p4est_quadrant_array_push (sc_array_t * array)
{
  P4EST_ASSERT (array->elem_size == sizeof (p4est_quadrant_t));

  return (p4est_quadrant_t *) sc_array_push (array);
}

/** Call sc_mempool_alloc for a mempool creating quadrants. */
/*@unused@*/
static inline p4est_quadrant_t *
p4est_quadrant_mempool_alloc (sc_mempool_t * mempool)
{
  P4EST_ASSERT (mempool->elem_size == sizeof (p4est_quadrant_t));

  return (p4est_quadrant_t *) sc_mempool_alloc (mempool);
}

/** Call sc_list pop for a quadrant array. */
/*@unused@*/
static inline p4est_quadrant_t *
p4est_quadrant_list_pop (sc_list_t * list)
{
  return (p4est_quadrant_t *) sc_list_pop (list);
}

SC_EXTERN_C_END;

#endif /* !P4EST_H */
