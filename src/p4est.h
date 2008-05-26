/*
  This file is part of p4est.
  p4est is a C library to manage a parallel collection of quadtrees and/or
  octrees.

  Copyright (C) 2007,2008 Carsten Burstedde, Lucas Wilcox.

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

#ifndef P4EST_H
#define P4EST_H

/* p4est_connectivity.h includes p4est_base.h sc_containers.h */
#include <p4est_connectivity.h>

/* finest level of the quadtree */
#define P4EST_MAXLEVEL 30

/* the length of a root quadrant */
#define P4EST_ROOT_LEN ((p4est_qcoord_t) 1 << P4EST_MAXLEVEL)

/* the length of a quadrant of level l */
#define P4EST_QUADRANT_LEN(l) ((p4est_qcoord_t) 1 << (P4EST_MAXLEVEL - (l)))

/* the offset of the highest quadrant at level l */
#define P4EST_LAST_OFFSET(l) (P4EST_ROOT_LEN - P4EST_QUADRANT_LEN (l))

/* a negative magic number for consistency checks */
#define P4EST_NEG_MAGIC = -439623172;

typedef struct p4est_quadrant
{
  p4est_qcoord_t      x, y;
  int8_t              level;
  union
  {
    void               *user_data;
    struct
    {
      p4est_locidx_t      which_tree;
      int                 owner_rank;
    }
    piggy;
  }
  p;
}
p4est_quadrant_t;

typedef struct p4est_tree
{
  sc_array_t          quadrants;        /* locally stored quadrants */
  p4est_locidx_t      quadrants_per_level[P4EST_MAXLEVEL + 1];  /* locals only */
  int8_t              maxlevel; /* highest local quadrant level */
}
p4est_tree_t;

typedef struct p4est_position
{
  p4est_locidx_t      which_tree;
  p4est_qcoord_t      x, y;
}
p4est_position_t;

typedef struct p4est
{
  MPI_Comm            mpicomm;
  int                 mpisize, mpirank;

  size_t              data_size;        /* size of per-quadrant user_data */
  void               *user_global_pointer;      /* convenience pointer for users
                                                   will never be touched by p4est */

  p4est_locidx_t      first_local_tree; /* 0-based index of first local tree,
                                           must be -1 for an empty processor */
  p4est_locidx_t      last_local_tree;  /* 0-based index of last local tree,
                                           must be -2 for an empty processor */
  p4est_locidx_t      local_num_quadrants;      /* number of quadrants
                                                   on all trees on this processor */
  p4est_gloidx_t      global_num_quadrants;     /* number of quadrants
                                                   on all trees on all processors */
  p4est_gloidx_t     *global_last_quad_index;   /* Index in the total ordering
                                                   of all quadrants of the
                                                   last quadrant on each proc.
                                                 */
  p4est_position_t   *global_first_position;    /* first_tree, x, y of first quadrant
                                                   for each processor and 1 beyond */
  p4est_connectivity_t *connectivity;   /* connectivity structure */
  sc_array_t         *trees;    /* list of all trees */

  sc_mempool_t       *user_data_pool;   /* memory allocator for user data
                                         * WARNING: This is NULL if data size
                                         *          equals zero.
                                         */
  sc_mempool_t       *quadrant_pool;    /* memory allocator for temporary quadrants */
}
p4est_t;

/** Callback function prototype to initialize the quadrant's user data.
 */
typedef void        (*p4est_init_t) (p4est_t * p4est,
                                     p4est_locidx_t which_tree,
                                     p4est_quadrant_t * quadrant);

/** Callback function prototype to decide for refinement.
 * \return Returns 1 if the quadrant shall be refined.
 */
typedef int         (*p4est_refine_t) (p4est_t * p4est,
                                       p4est_locidx_t which_tree,
                                       p4est_quadrant_t * quadrant);

/** Callback function prototype to decide for coarsening.
 * The quadrants are siblings enumerated in Morton ordering.
 * \return Returns 1 if the quadrants shall be replaced with their parent.
 */
typedef int         (*p4est_coarsen_t) (p4est_t * p4est,
                                        p4est_locidx_t which_tree,
                                        p4est_quadrant_t * q0,
                                        p4est_quadrant_t * q1,
                                        p4est_quadrant_t * q2,
                                        p4est_quadrant_t * q3);

/** Callback function prototype to calculate weights for partitioning.
 * \return  Returns a 32bit integer >= 0 as the quadrant weight.
 * \note    (global sum of weights * mpisize) must fit into a 64bit integer.
 */
typedef int         (*p4est_weight_t) (p4est_t * p4est,
                                       p4est_locidx_t which_tree,
                                       p4est_quadrant_t * quadrant);

/** map right-hand corner rule to corner in z-ordering */
extern const int    p4est_corner_to_zorder[5];

/** set statically allocated quadrant to defined values */
#define P4EST_QUADRANT_INIT(q) \
  do { memset ((q), -1, sizeof (p4est_quadrant_t)); } while (0)

/** Create a new p4est.
 *
 * \param [in] mpicomm A valid MPI_Comm or MPI_COMM_NULL.
 * \param [in] connectivity This is the connectivity information that
 *                          the forest is build with.  Note the p4est
 *                          does not take ownership of the memory.
 * \param [in] data_size This is the size of data for each quadrant which
 *                       can be zero.  If zero the \c user_data_pool is
 *                       set to \c NULL.
 * \param [in] init_fn Callback function to initialize the user_data
 *                     which is already allocated automatically.
 *
 * \return This returns a vaild forest.
 *
 * \note The connectivity structure must not be destroyed
 *       during the lifetime of this p4est.
 */
p4est_t            *p4est_new (MPI_Comm mpicomm,
                               p4est_connectivity_t * connectivity,
                               size_t data_size, p4est_init_t init_fn);

/** Destroy a p4est.
 * \note The connectivity structure is not destroyed with the p4est.
 */
void                p4est_destroy (p4est_t * p4est);

/** Make a deep copy of a p4est. Copying of user data is optional.
 * \param [in]  copy_data  If true, data are copied.
 *                         If false, data_size is set to 0.
 * \return  Returns a valid p4est that does not depend on the input.
 */
p4est_t            *p4est_copy (p4est_t * input, int copy_data);

/** Refine a forest.
 * \param [in] refine_fn Callback function to decide
 *                       if a quadrant gets refined
 * \param [in] init_fn   Callback function to initialize the user_data
 *                       which is already allocated automatically.
 */
void                p4est_refine (p4est_t * p4est,
                                  p4est_refine_t refine_fn,
                                  p4est_init_t init_fn);

/** Coarsen a forest.
 * \param [in] coarsen_fn Callback function to decide
 *                        if quadrants get coarsened
 * \param [in] init_fn   Callback function to initialize the user_data
 *                       which is already allocated automatically.
 */
void                p4est_coarsen (p4est_t * p4est,
                                   p4est_coarsen_t coarsen_fn,
                                   p4est_init_t init_fn);

/** Balance a forest. Currently only doing local balance.
 * \param [in] init_fn   Callback function to initialize the user_data
 *                       which is already allocated automatically.
 * \note Balances edges and corners.
 *       Can be easily changed to edges only in p4est_algorithms.c.
 */
void                p4est_balance (p4est_t * p4est, p4est_init_t init_fn);

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

#endif /* !P4EST_H */
