
#ifndef __P4EST_H__
#define __P4EST_H__

/* finest level of the quadtree */
#define P4EST_MAXLEVEL 30

/* this will be changed to 1 by make install */
#define P4EST_CONFIG_INSTALLED 0

/* this will be changed to 1 by make install if mpi is configured in */
#define P4EST_CONFIG_MPI 0

/* do some magic to avoid using p4est_config.h in the installed header */
#if P4EST_CONFIG_INSTALLED
#if P4EST_CONFIG_MPI
#include <mpi.h>
#else
#include <p4est_mpi_dummy.h>
#endif
#else
#ifdef HAVE_CONFIG_H
#include <p4est_config.h>
#endif
#ifdef HAVE_MPI
#include <mpi.h>
#else
#include <p4est_mpi_dummy.h>
#endif
#endif

/* include necessary headers */
#include <p4est_connectivity.h>
#include <p4est_memory.h>
#include <stdio.h>

typedef struct p4est_quadrant
{
  int32_t             x, y;
  int8_t              level;
  void               *user_data;
}
p4est_quadrant_t;

typedef struct p4est_tree
{
  p4est_array_t      *quadrants;        /* locally stored quadrants */
  int32_t             quadrants_per_level[P4EST_MAXLEVEL + 1];  /* locals only */
  int8_t              maxlevel; /* highest local quadrant level */
}
p4est_tree_t;

typedef struct p4est
{
  MPI_Comm            mpicomm;
  int                 mpisize, mpirank;

  FILE               *nout;     /* log messages go here if not NULL */
  int                 data_size;        /* size of user_data */

  int32_t             first_local_tree; /* 0-based index of first local tree */
  int32_t             last_local_tree;  /* 0-based index of last local tree */
  int32_t             local_num_trees;  /* number of trees on this processor */
  int32_t             local_num_quadrants;      /* number of quadrants
                                                   on all trees on this processor */
  int64_t             global_num_quadrants;     /* number of quadrants
                                                   on all trees on all processors */
  p4est_connectivity_t *connectivity;   /* connectivity structure */
  p4est_array_t      *trees;    /* list of all trees */

  p4est_mempool_t    *user_data_pool;   /* memory allocator for user data
                                         * WARNING: This is NULL if data size
                                         *          equals zero.
                                         */
  p4est_mempool_t    *quadrant_pool;    /* memory allocator for temporary quadrants */
}
p4est_t;

/** Callback function prototype to initialize the quadrant's user data.
 */
typedef void        (*p4est_init_t) (p4est_t * p4est, int32_t which_tree,
                                     p4est_quadrant_t * quadrant);

/** Callback function prototype to decide for refinement.
 * \return Returns 1 if the quadrant shall be refined
 */
typedef int         (*p4est_refine_t) (p4est_t * p4est, int32_t which_tree,
                                       p4est_quadrant_t * quadrant);

/** Create a new p4est.
 *
 * \param [in] mpicomm A valid MPI_Comm or MPI_COMM_NULL.
 * \param [in] nout    Stream for log messages.  If NULL then no messages
 *                     are logged.
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
p4est_t            *p4est_new (MPI_Comm mpicomm, FILE * nout,
                               p4est_connectivity_t * connectivity,
                               int data_size, p4est_init_t init_fn);

/** Destroy a p4est.
 * \note The connectivity structure is not destroyed with the p4est.
 */
void                p4est_destroy (p4est_t * p4est);

/** Refine a forest.
 * \param [in] refine_fn Callback function to decide
 *                       if a quadrant gets refined
 * \param [in] init_fn   Callback function to initialize the user_data
 *                       which is already allocated automatically.
 */
void                p4est_refine (p4est_t * p4est,
                                  p4est_refine_t refine_fn,
                                  p4est_init_t init_fn);

#endif /* !__P4EST_H__ */
