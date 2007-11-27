
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
#include <p4est_memory.h>
#include <stdio.h>
#include <stdint.h>

typedef struct p4est_connectivity
{
  int32_t             num_trees;
  int32_t             num_vertices;
  int32_t            *tree_to_vertex;   /* allocated [0][0]..[0][3]
                                           .. [num_trees-1][0]..[num_trees-1][3] */
  int32_t            *tree_to_tree;     /* allocated [0][0]..[0][3]
                                           .. [num_trees-1][0]..[num_trees-1][3] */
  int8_t             *tree_to_face;     /* allocated [0][0]..[0][3]
                                           .. [num_trees-1][0]..[num_trees-1][3] */
}
p4est_connectivity_t;

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
  int32_t             local_num_trees;  /* number of trees on this processor */
  int32_t             local_num_quadrants;      /* number of quadrants
                                                   on all trees on this processor */
  int64_t             global_num_quadrants;     /* number of quadrants
                                                   on all trees on all processors */
  p4est_connectivity_t *connectivity;   /* connectivity structure */
  p4est_array_t      *trees;    /* list of all trees */

  p4est_mempool_t    *user_data_pool;   /* memory allocator for user data */
  p4est_mempool_t    *quadrant_pool;    /* memory allocator for temporary quadrants */
}
p4est_t;

/*
 * callback function to initialize the quadrant's user data
 */
typedef void        (*p4est_init_t) (p4est_t * p4est, int32_t which_tree,
                                     p4est_quadrant_t * quadrant);

p4est_connectivity_t *p4est_connectivity_new (int32_t num_trees,
                                              int32_t num_vertices);
void                p4est_connectivity_destroy (p4est_connectivity_t *
                                                connectivity);

/*
 * TODO document ownership of connectivity
 */
p4est_t            *p4est_new (MPI_Comm mpicomm, FILE * nout,
                               p4est_connectivity_t * connectivity,
                               int data_size, p4est_init_t init_fn);
void                p4est_destroy (p4est_t * p4est);

#endif /* !__P4EST_H__ */
