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

#include <sc_options.h>
#ifndef P4_TO_P8
#include <p4est_extended.h>
#include <p4est_search.h>
#include <p4est_vtk.h>
#include <p4est_bits.h>
#else
#include <p8est_extended.h>
#include <p8est_search.h>
#include <p8est_vtk.h>
#include <p8est_bits.h>
#endif

#define COORDINATE_IROOTLEN (1. / P4EST_ROOT_LEN)

/* define centers of refinement and point creation */
typedef struct search_partition_global
{
  /* p4est mesh */
  double              a[3], b[3], c[3]; /* refinement centers */
  int                 uniform_level;    /* level of initial uniform refinement */
  int                 max_level;        /* maximum level of adaptive refinement */
  p4est_connectivity_t *conn;   /* brick connectivity of the mesh */
  p4est_t            *p4est;    /* the resulting p4est */

  /* query points */
  p4est_locidx_t      num_global_queries;       /* global number of queries;
                                                 * of type p4est_locidx_t, since
                                                 * queries are replicated across
                                                 * all processes */
  int                 seed;     /* seed for random query creation */
  double              clustering_exponent;      /* affects the distribution of queries */
  sc_array_t         *queries;  /* array of query points */

  /* search statistics */
  int                 num_local_queries;        /* queries found in local search */
  sc_array_t         *num_queries_per_rank;     /* queries found in partition search */
  sc_array_t         *global_nlq;       /* num_local_queries gathered globally */

  /* vtk output */
  int                 write_vtk;        /* boolean to activate vtk output */
  sc_array_t         *num_queries_per_quad;     /* queries found in each local quad */
}
search_partition_global_t;

typedef struct query_point
{
  double              xyz[3];   /* 3D coordinates */
  int                 is_local; /* set to 1, if found in local search */
  int                 rank;     /* rank assigned during partition search */
}
query_point_t;

static void
map_coordinates (p4est_qcoord_t x, p4est_qcoord_t y, p4est_qcoord_t z,
                 p4est_topidx_t which_tree, double *xyz)
{
  P4EST_ASSERT (which_tree < P4EST_CHILDREN);   /* assert we have a 2x2(x2) brick */
  xyz[0] = 0.5 * (COORDINATE_IROOTLEN * x + which_tree % 2);
  xyz[1] = 0.5 * (COORDINATE_IROOTLEN * y + (which_tree / 2) % 2);
#ifndef P4_TO_P8
  xyz[2] = 0.;
#else
  xyz[2] = 0.5 * (COORDINATE_IROOTLEN * z + which_tree / 4);
#endif
}

static int
refine_fn (p4est_t *p4est, p4est_topidx_t which_tree,
           p4est_quadrant_t *quadrant)
{
  p4est_qcoord_t      coords[3];
  double              xyz[3];
  double              dist, min_dist;

  P4EST_ASSERT (p4est != NULL);
  P4EST_ASSERT (p4est->user_pointer != NULL);
  search_partition_global_t *g =
    (search_partition_global_t *) p4est->user_pointer;
  P4EST_ASSERT (quadrant != NULL);

  /* get quadrant center reference coordinates in the unit square */
  P4EST_ASSERT (which_tree < P4EST_CHILDREN);   /* assert we have a 2x2(x2) brick */
  p4est_quadrant_volume_coordinates (quadrant, coords); /* quadrant center */
  map_coordinates (coords[0], coords[1],
#ifndef P4_TO_P8
                   0,
#else
                   coords[2],
#endif
                   which_tree, xyz);

  /* compute distance to point a */
  dist = (g->a[0] - xyz[0]) * (g->a[0] - xyz[0]) +
    (g->a[1] - xyz[1]) * (g->a[1] - xyz[1]) +
    (g->a[2] - xyz[2]) * (g->a[2] - xyz[2]);
  min_dist = sqrt (dist);

  /* compute distance to point b */
  dist = (g->b[0] - xyz[0]) * (g->b[0] - xyz[0]) +
    (g->b[1] - xyz[1]) * (g->b[1] - xyz[1]) +
    (g->b[2] - xyz[2]) * (g->b[2] - xyz[2]);
  min_dist = SC_MIN (min_dist, sqrt (dist));

  /* refine if quadrant center is close enough to either point a or point b */
  return (quadrant->level <
          g->max_level -
          floor (min_dist * (g->max_level - g->uniform_level) / 0.2));
}

static void
create_p4est (search_partition_global_t *g)
{
  int                 il;
  p4est_gloidx_t      old_gnq;

  /* Create brick p4est. */
  g->conn = p4est_connectivity_new_brick (2, 2,
#ifdef P4_TO_P8
                                          2,
#endif
                                          0, 0
#ifdef P4_TO_P8
                                          , 0
#endif
    );
  g->p4est =
    p4est_new_ext (sc_MPI_COMM_WORLD, g->conn, 0, g->uniform_level, 1, 0,
                   NULL, g);

  /* We like the invariant that after partitioning, partition-independent
     coarsening is always possible since every family of siblings is placed
     on a single process; this is not guaranteed by p4est_new_ext. */
  p4est_partition (g->p4est, 1, NULL);

  /* refine the forest adaptively around two points g->a and g->b */
  for (il = g->uniform_level; il < g->max_level; il++) {
    old_gnq = g->p4est->global_num_quadrants;
    p4est_refine (g->p4est, 0, refine_fn, NULL);

    /* leave the loop as soon as no new refinement occurs */
    if (old_gnq == g->p4est->global_num_quadrants) {
      break;
    }

    p4est_partition (g->p4est, 0, NULL);
  }
}

static void
generate_queries (search_partition_global_t *g)
{
  p4est_locidx_t      iq, nqh;
  int                 id;
  query_point_t      *p;
  double              t;
  int                 mpiret;

  /* generate local queries */
  g->queries =
    sc_array_new_count (sizeof (query_point_t), g->num_global_queries);
  sc_array_memset (g->queries, 0);
  if (g->p4est->mpirank == 0) {
    srand (g->seed);
    nqh = g->num_global_queries / 2;
    for (iq = 0; iq < g->num_global_queries; iq++) {
      p = (query_point_t *) sc_array_index_int (g->queries, iq);
      p->is_local = 0;
      p->rank = -1;
      for (id = 0; id < P4EST_DIM; id++) {
        p->xyz[id] = rand () / (RAND_MAX + 1.);
      }

      /* move point closer to g->b or g->c depending on iq and random t */
      t = pow (rand () / (RAND_MAX + 1.), g->clustering_exponent);
      /* move the point to position sp->xyz * t + (1 - t) * {g->b,g->c} */
      if (iq < nqh) {
        p->xyz[0] = t * p->xyz[0] + (1 - t) * g->b[0];
        p->xyz[1] = t * p->xyz[1] + (1 - t) * g->b[1];
        p->xyz[2] = t * p->xyz[2] + (1 - t) * g->b[2];
      }
      else {
        p->xyz[0] = t * p->xyz[0] + (1 - t) * g->c[0];
        p->xyz[1] = t * p->xyz[1] + (1 - t) * g->c[1];
        p->xyz[2] = t * p->xyz[2] + (1 - t) * g->c[2];
      }
    }
  }

  /* broadcast queries to all processes */
  mpiret = sc_MPI_Bcast (g->queries->array,
                         g->num_global_queries * sizeof (query_point_t),
                         sc_MPI_BYTE, 0, g->p4est->mpicomm);
  SC_CHECK_MPI (mpiret);
  P4EST_GLOBAL_PRODUCTIONF ("Created %lld global queries.\n",
                            (unsigned long long) g->queries->elem_count);
}

static int
quadrant_contains_query (p4est_quadrant_t *quadrant,
                         p4est_topidx_t which_tree, query_point_t *p)
{
  double              qxyz[3];
  double              qlen;

  /* compute lower, left corners coords for quadrant bounds */
  map_coordinates (quadrant->x, quadrant->y,
#ifndef P4_TO_P8
                   0,
#else
                   quadrant->z,
#endif
                   which_tree, qxyz);
  qlen = 0.5 * P4EST_QUADRANT_LEN (quadrant->level) * COORDINATE_IROOTLEN;

  /* for the replicated points we can operate on half-open intervals without
   * adding a tolerance */
  if (p->xyz[0] < qxyz[0] || p->xyz[0] > qxyz[0] + qlen ||
      p->xyz[1] < qxyz[1] || p->xyz[1] > qxyz[1] + qlen ||
      p->xyz[2] < qxyz[2] || p->xyz[2] > qxyz[2] + qlen) {
    return 0;
  }
  return 1;
}

static int
local_callback (p4est_t *p4est, p4est_topidx_t which_tree,
                p4est_quadrant_t *quadrant, p4est_locidx_t local_num,
                void *point)
{
  P4EST_ASSERT (point != NULL);
  query_point_t      *p = (query_point_t *) point;
  P4EST_ASSERT (p4est != NULL);
  P4EST_ASSERT (p4est->user_pointer != NULL);
  search_partition_global_t *g =
    (search_partition_global_t *) p4est->user_pointer;

  if (!quadrant_contains_query (quadrant, which_tree, p)) {
    return 0;
  }

  if (local_num >= 0) {
    /* we are on a local leaf */
    p->is_local = 1;
    g->num_local_queries++;
    *(double *) sc_array_index (g->num_queries_per_quad, local_num) += 1.;
  }

  return 1;
}

static void
search_local (search_partition_global_t *g)
{
  int                 gnq;
  size_t              il;
  int                 mpiret;

  /* search queries locally */
  g->num_local_queries = 0;
  g->num_queries_per_quad =
    sc_array_new_count (sizeof (double), g->p4est->local_num_quadrants);
  sc_array_memset (g->num_queries_per_quad, 0);
  p4est_search_local (g->p4est, 0, NULL, local_callback, g->queries);
  P4EST_INFOF ("Queries found in local search = %d\n", g->num_local_queries);

  /* allgather local num queries for future comparison with partition search */
  g->global_nlq = sc_array_new_count (sizeof (int), g->p4est->mpisize);
  mpiret = sc_MPI_Allgather (&g->num_local_queries, 1, sc_MPI_INT,
                             g->global_nlq->array, 1, sc_MPI_INT,
                             g->p4est->mpicomm);
  SC_CHECK_MPI (mpiret);

  /* compute sum of local number of queries across all processes */
  gnq = 0;
  for (il = 0; il < g->global_nlq->elem_count; il++) {
    gnq += *(int *) sc_array_index (g->global_nlq, il);
  }
  P4EST_GLOBAL_PRODUCTIONF
    ("Queries found globally during local search: %d (expected %d)\n",
     gnq, (int) g->num_global_queries);
  P4EST_ASSERT (g->num_global_queries == (p4est_locidx_t) gnq);
}

static int
partition_callback (p4est_t *p4est, p4est_topidx_t which_tree,
                    p4est_quadrant_t *quadrant, int pfirst, int plast,
                    void *point)
{
  P4EST_ASSERT (point != NULL);
  query_point_t      *p = (query_point_t *) point;
  P4EST_ASSERT (p4est != NULL);
  P4EST_ASSERT (p4est->user_pointer != NULL);
  search_partition_global_t *g =
    (search_partition_global_t *) p4est->user_pointer;

  if (!quadrant_contains_query (quadrant, which_tree, p)) {
    return 0;
  }

  if (pfirst == plast) {
    /* we are on a local leaf */
    p->rank = pfirst;
    P4EST_ASSERT (((p->rank == p4est->mpirank) && p->is_local) ||
                  ((p->rank != p4est->mpirank) && !p->is_local));
    *(size_t *) sc_array_index (g->num_queries_per_rank, pfirst) += 1;
  }

  return 1;
}

static void
search_partition (search_partition_global_t *g)
{
  size_t              iz, lenz, buffer_size;
  int                 num_queries_found;
  int                 retb;
  char               *buffer;

  /* search queries in the partition */
  g->num_queries_per_rank =
    sc_array_new_count (sizeof (int), g->p4est->mpisize);
  sc_array_memset (g->num_queries_per_rank, 0);
  p4est_search_partition (g->p4est, 0, NULL, partition_callback, g->queries);

  /* output query points found per rank */
  buffer_size = 0;
  for (iz = 0; iz < g->num_queries_per_rank->elem_count; iz++) {
    if (iz % 10 == 0) {
      buffer_size += 1;
    }
    retb = snprintf (NULL, 0, "%7d ", *(int *)
                     sc_array_index (g->num_queries_per_rank, iz));
    SC_CHECK_ABORT (retb > 0, "Overflow in snprintf");
    buffer_size += (size_t) retb;
  }

  /* now we know the buffer size and allocate it */
  buffer = P4EST_ALLOC (char, buffer_size);

  lenz = 0;
  num_queries_found = 0;
  for (iz = 0; iz < g->num_queries_per_rank->elem_count; iz++) {
    /* add result to buffer */
    if (iz % 10 == 0) {
      P4EST_ASSERT (buffer_size >= lenz);
      retb = snprintf (buffer + lenz, buffer_size - lenz, "\n");
      SC_CHECK_ABORT (retb == 1, "Overflow in snprintf");
      lenz += retb;
    }
    retb = snprintf (buffer + lenz, buffer_size - lenz, "%7d ", *(int *)
                     sc_array_index (g->num_queries_per_rank, iz));
    SC_CHECK_ABORT (retb > 0, "Overflow in snprintf");
    lenz += retb;

    /* compute total number of queries found during partition search */
    num_queries_found +=
      *(int *) sc_array_index (g->num_queries_per_rank, iz);
  }
  P4EST_GLOBAL_PRODUCTIONF
    ("Queries found during partition search: %d (expected %d)\n",
     num_queries_found, (int) g->num_global_queries);
  P4EST_GLOBAL_STATISTICSF
    ("Partition search found the following query counts %s\n", buffer);

  /* the buffer is no longer accessed */
  P4EST_FREE (buffer);

  /* check results for consistency */
  for (iz = 0; iz < g->num_queries_per_rank->elem_count; iz++) {
    P4EST_ASSERT (*(int *) sc_array_index (g->num_queries_per_rank, iz) ==
                  *(int *) sc_array_index (g->global_nlq, iz));
  }
}

static void
write_vtk (search_partition_global_t *g)
{
  char                filename[BUFSIZ];
  p4est_vtk_context_t *cont;

  /* write VTK output of sphere counts */
  if (!g->write_vtk) {
    return;
  }

  /* run-once loop for clean return */
  cont = NULL;
  do {
    /* open files for output */
    snprintf (filename, BUFSIZ, "search_partition%d_%d_%d_%d_%d_%.2f",
              P4EST_DIM, g->uniform_level, g->max_level,
              (int) g->num_global_queries, g->seed, g->clustering_exponent);
    cont = p4est_vtk_context_new (g->p4est, filename);
    if (NULL == p4est_vtk_write_header (cont)) {
      P4EST_LERRORF ("Failed to write header for %s\n", filename);
      break;
    }

    /* write cell data to file */
    if (NULL == p4est_vtk_write_cell_dataf
        (cont, 1, 1, 1, 0, 1, 0, "num_queries", g->num_queries_per_quad,
         cont)) {
      P4EST_LERRORF ("Failed to write cell data for %s\n", filename);
      break;
    }

    /* finish meta information and close files */
    if (p4est_vtk_write_footer (cont)) {
      P4EST_LERRORF ("Failed to write footer for %s\n", filename);
      break;
    }
  }
  while (0);
}

static void
cleanup (search_partition_global_t *g)
{
  sc_array_destroy (g->queries);
  sc_array_destroy (g->global_nlq);
  sc_array_destroy (g->num_queries_per_quad);
  sc_array_destroy (g->num_queries_per_rank);
  p4est_destroy (g->p4est);
  p4est_connectivity_destroy (g->conn);
}

static void
run (search_partition_global_t *g)
{
  /* create a 2x2(x2) brick p4est mesh covering the unit square */
  create_p4est (g);

  /* generate search queries in the unit square */
  generate_queries (g);

  /* search queries in the local part of the mesh */
  search_local (g);

  /* search queries in the partition of the mesh and compare */
  search_partition (g);

  /* output forest to vtk */
  write_vtk (g);

  /* Free memory. */
  cleanup (g);
}

int
main (int argc, char **argv)
{
  int                 mpiret;
  sc_options_t       *opt;
  int                 first_argc, ue, help;
  search_partition_global_t global, *g = &global;
  int                 ngq;

  /* MPI initialization. */
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  /* Package init. */
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  opt = sc_options_new (argv[0]);
  sc_options_add_switch (opt, 'h', "help", &help,
                         "Print help message and exit cleanly");
  sc_options_add_int (opt, 'l', "minlevel", &g->uniform_level, 3,
                      "Level of uniform refinement");
  sc_options_add_int (opt, 'L', "maxlevel", &g->max_level, 5,
                      "Level of maximum refinement");
  sc_options_add_int (opt, 'q', "num-queries", &ngq, 100,
                      "Number of queries created per process");
  sc_options_add_int (opt, 's', "seed", &g->seed, 0,
                      "Seed for random queries");
  sc_options_add_bool (opt, 'v', "write-vtk", &g->write_vtk, 0,
                       "Activate vtk output");
  sc_options_add_double (opt, 'c', "clustering-exponent",
                         &g->clustering_exponent, 0.5,
                         "Clustering of queries");

  /* proceed in run-once loop for clean abort */
  ue = 0;
  do {
    first_argc = sc_options_parse (p4est_get_package_id (), SC_LP_DEFAULT,
                                   opt, argc, argv);
    /* num_global_queries is a p4est_locidx_t, as we replicate the queries on
     * all processes */
    g->num_global_queries = (p4est_locidx_t) ngq;

    if (first_argc < 0) {
      P4EST_GLOBAL_LERROR ("Invalid option format.\n");
      ue = 1;
      break;
    }

    /* check options for consistency */
    if (g->uniform_level < 0 || g->uniform_level > P4EST_QMAXLEVEL) {
      P4EST_GLOBAL_LERRORF ("Uniform level out of bounds 0..%d\n",
                            P4EST_QMAXLEVEL);
      ue = 1;
    }
    if (g->max_level < 0 || g->max_level > P4EST_QMAXLEVEL) {
      P4EST_GLOBAL_LERRORF ("Maximum level out of bounds 0..%d\n",
                            P4EST_QMAXLEVEL);
      ue = 1;
    }
    if (g->num_global_queries < 0) {
      P4EST_GLOBAL_LERROR ("Number of queries has to be non-negative.\n");
      ue = 1;
    }
    if ((long long) g->num_global_queries * sizeof (query_point_t) > INT_MAX) {
      P4EST_GLOBAL_LERROR ("Number of queries too large for MPI buffer.\n");
      ue = 1;
    }
    if (g->seed < 0) {
      P4EST_GLOBAL_LERROR ("Seed has to be non-negative.\n");
      ue = 1;
    }
    if (g->clustering_exponent < 0.) {
      P4EST_GLOBAL_LERROR ("Clustering exponent has to be non-negative.\n");
      ue = 1;
    }
    if (ue || help) {
      break;
    }
    sc_options_print_summary (p4est_get_package_id (), SC_LP_ESSENTIAL, opt);

    /* define centers of refinement and point creation */
    g->a[0] = 0.2;
    g->a[1] = 0.4;
    g->a[2] = 0.4;
    g->b[0] = 0.7;
    g->b[1] = 0.55;
    g->b[2] = 0.55;
    g->c[0] = 0.3;
    g->c[1] = 0.8;
    g->c[2] = 0.8;
#ifndef P4_TO_P8
    g->a[2] = g->b[2] = g->c[2] = 0.;   /* reset z-coordinate to 0 */
#endif

    /* run example */
    run (g);
  }
  while (0);
  if (ue || help) {
    sc_options_print_usage (p4est_get_package_id (), SC_LP_ERROR, opt, NULL);
  }

  /* Close MPI environment. */
  sc_options_destroy (opt);
  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return ue ? EXIT_FAILURE : EXIT_SUCCESS;
}
