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
#include <unistd.h>
#include <inttypes.h>
#ifndef P4_TO_P8
#include <p4est_extended.h>
#include <p4est_ghost.h>
#include <p4est_mesh.h>
#include <p4est_vtk.h>
#else /* !P4_TO_P8 */
#include <p8est_extended.h>
#include <p8est_ghost.h>
#include <p8est_mesh.h>
#include <p8est_vtk.h>
#endif /* !P4_TO_P8 */

/** Per-quadrant data
 *
 * Marker
 */
typedef struct test_mesh_data
{
  double              marker;
}
test_mesh_marker_t;

/** Init function for quadrants
 *
 * \param [in] p4est       The forest
 * \param [in] which_tree  Index of currently processed tree
 * \param [in] q           Quadrant which payload will be set
 */
static void
test_mesh_init (p4est_t * p4est, p4est_topidx_t which_tree,
                p4est_quadrant_t * q)
{
  test_mesh_marker_t *data = (test_mesh_marker_t *) q->p.user_data;
  data->marker = -1.;
}

/** Function to fetch payload of each quadrant on current process
 *
 * \param [in]      info        Volume callback info
 * \param [in][out] user_data   Array which will be set to value of markers
 */
static void
test_mesh_collect_markers (p4est_iter_volume_info_t * info, void *user_data)
{
  /* fetch payload */
  p4est_quadrant_t   *q = info->quad;
  test_mesh_marker_t *data = (test_mesh_marker_t *) q->p.user_data;

  /* fetch user_data */
  sc_array_t         *markers = (sc_array_t *) user_data;
  double             *marker_ptr;

  /* calculate position to write to and write data there */
  p4est_t            *p4est = info->p4est;
  p4est_topidx_t      which_tree = info->treeid;
  p4est_locidx_t      local_id = info->quadid;  /* this is the index of q
                                                   within its tree's numbering.
                                                   We want to convert the index
                                                   for all the quadrants on this
                                                   process, which we do below */
  p4est_tree_t       *tree;

  tree = p4est_tree_array_index (p4est->trees, which_tree);
  local_id += tree->quadrants_offset;   /* now the id is relative to the MPI
                                           process */

  marker_ptr = (double *) sc_array_index_int (markers, local_id);
  marker_ptr[0] = (double) data->marker;
}

/** Create VTK output
 *
 * \param [in] p4est     The forest
 * \param [in] filename  The base name of the vtk files
 */
static int
test_mesh_write_vtk (p4est_t * p4est, char *filename)
{
  sc_array_t         *marks;
  marks = sc_array_new_size (sizeof (double), p4est->local_num_quadrants);

  p4est_iterate (p4est, NULL, (void *) marks, test_mesh_collect_markers, NULL,
#ifdef P4_TO_P8
                 NULL,
#endif /* P4_TO_P8 */
                 NULL);

  /* create VTK output context and set its parameters */
  p4est_vtk_context_t *context = p4est_vtk_context_new (p4est, filename);
  p4est_vtk_context_set_scale (context, 1);     /* quadrant at almost full scale */

  /* begin writing the output files */
  context = p4est_vtk_write_header (context);
  SC_CHECK_ABORT (context != NULL,
                  P4EST_STRING "_vtk: Error writing vtk header");
  context = p4est_vtk_write_cell_dataf (context, 1,
                                        /* write tree indices */
                                        1,      /* write the refinement level */
                                        1,      /* write the mpi process id */
                                        0,      /* do not wrap the mpi rank */
                                        1,      /* write marks as scalar cell
                                                   data */
                                        0,      /* no custom cell vector data */
                                        "marks", marks, context);
  SC_CHECK_ABORT (context != NULL,
                  P4EST_STRING "_vtk: Error writing cell data");

  const int           retval = p4est_vtk_write_footer (context);
  SC_CHECK_ABORT (!retval, P4EST_STRING "_vtk: Error writing footer");

  sc_array_destroy (marks);

  return 0;
}

/* Function for refining a mesh exactly once, in the very first quadrant
 *
 * \param [in] p4est        The forest.
 * \param [in] which_tree   The tree index of the current quadrant \a q
 * \param [in] q            The currently considered quadrant
 */
static int
refineExactlyOnce (p4est_t * p4est, p4est_topidx_t which_tree,
                   p4est_quadrant_t * q)
{
  int                 dec;
#ifndef P4_TO_P8
  dec = q->x == 0 && q->y == 0 && which_tree == 0;
#else /* !P4_TO_P8 */
  dec = q->x == 0 && q->y == 0 && q->z == 0 && which_tree == 0;
#endif /* !P4_TO_P8 */
  if (dec) {
    return 1;
  }
  return 0;
}

/* Function prototype to check the created mesh
 *
 * \param [in] p4est    The forest.
 * \param [in] ghost    The process-local ghost-layer.
 * \param [in] mesh     The process-local mesh.
 */
int
check_mesh (p4est_t * p4est, p4est_ghost_t * ghost, p4est_mesh_t * mesh,
            char *scenario)
{
  int                 printFaces = 1;
#ifdef P4_TO_P8
  int                 printEdges = 1;
#endif /* P4_TO_P8 */
  int                 printCorners = 0;
  int                 direction = 0;
  int                 encoding;
  uint32_t            norm_quad;
  uint32_t            quad, i;
  int                 j;
#ifdef P4_TO_P8
  int                 dir_offset_edge, dir_offset_corner;
  dir_offset_edge = P4EST_FACES;
  dir_offset_corner = P4EST_FACES + P8EST_EDGES;
#else /* P4_TO_P8 */
  int                 dir_offset_corner;
  dir_offset_corner = P4EST_FACES;
#endif /* P4_TO_P8 */

  for (quad = 0; quad < p4est->global_num_quadrants; ++quad) {
    sc_MPI_Barrier (p4est->mpicomm);
    /* loop over all cells, set and unset only if cell is owned by processor */
    if (p4est->global_first_quadrant[p4est->mpirank] <= quad
        && quad < p4est->global_first_quadrant[p4est->mpirank + 1]) {

      /* norm global quad index to local index */
      norm_quad = quad - p4est->global_first_quadrant[p4est->mpirank];

      /* set */
      if (printFaces) {
        for (i = 0; i < P4EST_FACES; ++i) {
          direction = i;
          sc_array_t         *neighboring_quads, *neighboring_encs;
          sc_array_t         *neighboring_qids;
          neighboring_quads = sc_array_new (sizeof (p4est_quadrant_t *));
          neighboring_encs = sc_array_new (sizeof (int));
          neighboring_qids = sc_array_new (sizeof (int));
          p4est_mesh_get_neighbors (p4est, ghost, mesh, norm_quad,
                                    direction, neighboring_quads,
                                    neighboring_encs, neighboring_qids);

#ifdef P4EST_ENABLE_DEBUG
          /* print some debug info */
          printf ("rank %5i, local quad %5i, global quad %5i, direction: %2i,"
                  " number of neighboring cells: %zu\n",
                  p4est->mpirank, quad, norm_quad, direction,
                  neighboring_quads->elem_count);
#endif /* P4EST_ENABLE_DEBUG */

          for (j = 0; j < neighboring_quads->elem_count; ++j) {
            p4est_quadrant_t   *q =
              *(p4est_quadrant_t **) sc_array_index_int (neighboring_quads,
                                                         j);
            encoding = *(int *) sc_array_index_int (neighboring_encs, j);
            if (encoding > 0) {
              test_mesh_marker_t *marker = (test_mesh_marker_t *)
                q->p.user_data;
              marker->marker = (double) direction;
            }
          }
          sc_array_destroy (neighboring_quads);
          sc_array_destroy (neighboring_encs);
          sc_array_destroy (neighboring_qids);
        }
      }
#ifdef P4_TO_P8
      if (printEdges) {
        for (i = 0; i < P8EST_EDGES; ++i) {
          direction = dir_offset_edge + i;
          sc_array_t         *neighboring_quads, *neighboring_encs;
          sc_array_t         *neighboring_qids;
          neighboring_quads = sc_array_new (sizeof (p4est_quadrant_t *));
          neighboring_encs = sc_array_new (sizeof (int));
          neighboring_qids = sc_array_new (sizeof (int));
          p4est_mesh_get_neighbors (p4est, ghost, mesh, norm_quad,
                                    direction, neighboring_quads,
                                    neighboring_encs, neighboring_qids);

#ifdef P4EST_ENABLE_DEBUG
          /* print some debug info */
          printf ("rank %5i, local quad %5i, global quad %5i, direction: %2i,"
                  " number of neighboring cells: %zu\n",
                  p4est->mpirank, quad, norm_quad, direction,
                  neighboring_quads->elem_count);
#endif /* P4EST_ENABLE_DEBUG */

          for (j = 0; j < neighboring_quads->elem_count; ++j) {
            p4est_quadrant_t   *q =
              *(p4est_quadrant_t **) sc_array_index_int (neighboring_quads,
                                                         j);
            encoding = *(int *) sc_array_index_int (neighboring_encs, j);
            if (encoding > 0) {
              test_mesh_marker_t *marker = (test_mesh_marker_t *)
                q->p.user_data;
              marker->marker = (double) direction;
            }
          }
          sc_array_destroy (neighboring_quads);
          sc_array_destroy (neighboring_encs);
          sc_array_destroy (neighboring_qids);
        }
      }
#endif /* P4_TO_P8 */
      if (printCorners) {
        for (i = 0; i < P4EST_CHILDREN; ++i) {
          direction = dir_offset_corner + i;
          sc_array_t         *neighboring_quads, *neighboring_encs;
          sc_array_t         *neighboring_qids;
          neighboring_quads = sc_array_new (sizeof (p4est_quadrant_t *));
          neighboring_encs = sc_array_new (sizeof (int));
          neighboring_qids = sc_array_new (sizeof (int));
          p4est_mesh_get_neighbors (p4est, ghost, mesh, norm_quad,
                                    direction, neighboring_quads,
                                    neighboring_encs, neighboring_qids);

#ifdef P4EST_ENABLE_DEBUG
          /* print some debug info */
          printf ("rank %5i, local quad %5i, global quad %5i, direction: %2i,"
                  " number of neighboring cells: %zu\n",
                  p4est->mpirank, quad, norm_quad, direction,
                  neighboring_quads->elem_count);
#endif /* P4EST_ENABLE_DEBUG */

          for (j = 0; j < neighboring_quads->elem_count; ++j) {
            p4est_quadrant_t   *q =
              *(p4est_quadrant_t **) sc_array_index_int (neighboring_quads,
                                                         j);
            encoding = *(int *) sc_array_index_int (neighboring_encs, j);
            if (encoding > 0) {
              test_mesh_marker_t *marker = (test_mesh_marker_t *)
                q->p.user_data;
              marker->marker = (double) direction;
            }
          }
          sc_array_destroy (neighboring_quads);
          sc_array_destroy (neighboring_encs);
          sc_array_destroy (neighboring_qids);
        }
      }
    }

    /* write vtk output */
    char                filename[120];
    sprintf (filename, "%s_test_mesh_%s_neighbors_quad_%i", P4EST_STRING,
             scenario, quad);

    test_mesh_write_vtk (p4est, filename);

    /* unset */
    if (p4est->global_first_quadrant[p4est->mpirank] <= quad
        && p4est->global_first_quadrant[p4est->mpirank + 1]) {
      if (printFaces) {
        for (i = 0; i < P4EST_FACES; ++i) {
          direction = i;
          sc_array_t         *neighboring_quads, *neighboring_encs;
          sc_array_t         *neighboring_qids;
          neighboring_quads = sc_array_new (sizeof (p4est_quadrant_t *));
          neighboring_encs = sc_array_new (sizeof (int));
          neighboring_qids = sc_array_new (sizeof (int));
          p4est_mesh_get_neighbors (p4est, ghost, mesh, norm_quad,
                                    direction, neighboring_quads,
                                    neighboring_encs, neighboring_qids);
          for (j = 0; j < neighboring_quads->elem_count; ++j) {
            p4est_quadrant_t   *q =
              *(p4est_quadrant_t **) sc_array_index_int (neighboring_quads,
                                                         j);
            encoding = *(int *) sc_array_index_int (neighboring_encs, j);
            if (encoding > 0) {
              test_mesh_marker_t *marker =
                (test_mesh_marker_t *) q->p.user_data;
              marker->marker = -1.;
            }
          }
          sc_array_destroy (neighboring_quads);
          sc_array_destroy (neighboring_encs);
          sc_array_destroy (neighboring_qids);
        }
      }
#ifdef P4_TO_P8
      if (printEdges) {
        for (i = 0; i < P8EST_EDGES; ++i) {
          direction = dir_offset_edge + i;
          sc_array_t         *neighboring_quads, *neighboring_encs;
          sc_array_t         *neighboring_qids;
          neighboring_quads = sc_array_new (sizeof (p4est_quadrant_t *));
          neighboring_encs = sc_array_new (sizeof (int));
          neighboring_qids = sc_array_new (sizeof (int));
          p4est_mesh_get_neighbors (p4est, ghost, mesh, norm_quad,
                                    direction, neighboring_quads,
                                    neighboring_encs, neighboring_qids);

          for (j = 0; j < neighboring_quads->elem_count; ++j) {
            p4est_quadrant_t   *q =
              *(p4est_quadrant_t **) sc_array_index_int (neighboring_quads,
                                                         j);
            encoding = *(int *) sc_array_index_int (neighboring_encs, j);
            if (encoding > 0) {
              test_mesh_marker_t *marker = (test_mesh_marker_t *)
                q->p.user_data;
              marker->marker = -1.;
            }
          }
          sc_array_destroy (neighboring_quads);
          sc_array_destroy (neighboring_encs);
          sc_array_destroy (neighboring_qids);
        }
      }
#endif /* P4_TO_P8 */
      if (printCorners) {
        for (i = 0; i < P4EST_CHILDREN; ++i) {
          direction = dir_offset_corner + i;
          sc_array_t         *neighboring_quads, *neighboring_encs;
          sc_array_t         *neighboring_qids;
          neighboring_quads = sc_array_new (sizeof (p4est_quadrant_t *));
          neighboring_encs = sc_array_new (sizeof (int));
          neighboring_qids = sc_array_new (sizeof (int));
          p4est_mesh_get_neighbors (p4est, ghost, mesh, norm_quad,
                                    direction, neighboring_quads,
                                    neighboring_encs, neighboring_qids);

          for (j = 0; j < neighboring_quads->elem_count; ++j) {
            p4est_quadrant_t   *q =
              *(p4est_quadrant_t **) sc_array_index_int (neighboring_quads,
                                                         j);
            encoding = *(int *) sc_array_index_int (neighboring_encs, j);
            if (encoding > 0) {
              test_mesh_marker_t *marker = (test_mesh_marker_t *)
                q->p.user_data;
              marker->marker = -1.;
            }
          }
          sc_array_destroy (neighboring_quads);
          sc_array_destroy (neighboring_encs);
          sc_array_destroy (neighboring_qids);
        }
      }
    }
  }

  return 0;
}

/* Function for testing p4est-mesh for a single tree scenario
 *
 * \param [in] p4est     The forest.
 * \param [in] conn      The connectivity structure
 * \param [in] periodic  Flag for checking if we have periodic boundaries
 * \returns 0 for success, -1 for failure
 */
int
test_mesh_one_tree (p4est_t * p4est,
                    p4est_connectivity_t * conn,
                    int8_t periodic, sc_MPI_Comm mpicomm)
{
  /* ensure that we have null pointers at beginning and end of function */
  P4EST_ASSERT (p4est == NULL);
  P4EST_ASSERT (conn == NULL);
  /* create connectivity */
#ifndef P4_TO_P8
  conn =
    periodic == 1 ? p4est_connectivity_new_periodic () :
    p4est_connectivity_new_unitsquare ();
#else /* !P4_TO_P8 */
  conn =
    periodic == 1 ? p8est_connectivity_new_periodic () :
    p8est_connectivity_new_unitcube ();
#endif /* !P4_TO_P8 */
  /* setup p4est */
  int                 minLevel = 2;
  p4est = p4est_new_ext (mpicomm,
                         conn,
                         0,
                         minLevel,
                         1, sizeof (test_mesh_marker_t), test_mesh_init, 0);

#if 0
  p4est_refine (p4est, 0, refineExactlyOnce, test_mesh_init);
  p4est_partition (p4est, 0, 0);
  p4est_balance (p4est, P4EST_CONNECT_FULL, test_mesh_init);
#endif

  /* inspect setup of geometry and check if payload is set correctly */
  char                filename[35] = "test_mesh_setup_single_tree_";
  strcat (filename, P4EST_STRING);
  test_mesh_write_vtk (p4est, filename);

  /* create mesh */
#ifdef P4_TO_P8
  p4est_ghost_t      *ghost = p4est_ghost_new (p4est, P8EST_CONNECT_EDGE);
  p4est_mesh_t       *mesh =
    p4est_mesh_new_ext (p4est, ghost, 1, 1, P8EST_CONNECT_EDGE);
#else /* P4_TO_P8 */
  p4est_ghost_t      *ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FULL);
  p4est_mesh_t       *mesh =
    p4est_mesh_new_ext (p4est, ghost, 1, 1, P4EST_CONNECT_FULL);
#endif /* P4_TO_P8 */

  /* check mesh */
  char                scenario[30];
  snprintf (scenario, 30, (periodic ? "single_tree_p" : "single_tree_np"));
  check_mesh (p4est, ghost, mesh, scenario);

  /* cleanup */
  p4est_ghost_destroy (ghost);
  p4est_mesh_destroy (mesh);
  p4est_destroy (p4est);
  p4est_connectivity_destroy (conn);
  conn = 0;
  p4est = 0;
  P4EST_ASSERT (p4est == NULL);
  P4EST_ASSERT (conn == NULL);
  return 0;
}

/* Function for testing p4est-mesh for multiple trees in a brick scenario
 *
 * \param [in] p4est     The forest.
 * \param [in] conn      The connectivity structure
 * \param [in] periodic  Flag for checking if we have periodic boundaries
 * \returns 0 for success, -1 for failure
 */
int
test_mesh_multiple_trees_brick (p4est_t * p4est,
                                p4est_connectivity_t * conn,
                                int8_t periodic, sc_MPI_Comm mpicomm)
{
  /* ensure that we have null pointers at beginning and end of function */
  P4EST_ASSERT (p4est == NULL);
  P4EST_ASSERT (conn == NULL);
  /* create connectivity */
#ifndef P4_TO_P8
  conn = p4est_connectivity_new_brick (2, 1, periodic, periodic);
#else /* !P4_TO_P8 */
  conn = p8est_connectivity_new_brick (2, 1, 1, periodic, periodic, periodic);
#endif /* !P4_TO_P8 */
  /* setup p4est */
  int                 minLevel = 2;
  p4est = p4est_new_ext (mpicomm,
                         conn,
                         0,
                         minLevel,
                         1, sizeof (test_mesh_marker_t), test_mesh_init, 0);
  /*
     p4est_refine (p4est, 0, refineExactlyOnce, 0);
     p4est_partition (p4est, 0, 0);
     p4est_balance (p4est, P4EST_CONNECT_FULL, 0);
   */

  char                filename[29] = "test_mesh_setup_brick_";
  strcat (filename, P4EST_STRING);
  p4est_vtk_write_file (p4est, 0, filename);
  /* create mesh */

#ifdef P4_TO_P8
  p4est_ghost_t      *ghost = p4est_ghost_new (p4est, P8EST_CONNECT_EDGE);
  p4est_mesh_t       *mesh =
    p4est_mesh_new_ext (p4est, ghost, 1, 1, P8EST_CONNECT_EDGE);
#else /* P4_TO_P8 */
  p4est_ghost_t      *ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FULL);
  p4est_mesh_t       *mesh =
    p4est_mesh_new_ext (p4est, ghost, 1, 1, P4EST_CONNECT_FULL);
#endif /* P4_TO_P8 */

  /* check mesh */
  char                scenario[30];
  snprintf (scenario, 30,
            (periodic ? "multiple_tree_brick_p" : "multiple_tree_brick_np"));
  check_mesh (p4est, ghost, mesh, scenario);

  /* cleanup */
  p4est_ghost_destroy (ghost);
  p4est_mesh_destroy (mesh);
  p4est_destroy (p4est);
  p4est_connectivity_destroy (conn);
  conn = 0;
  p4est = 0;
  P4EST_ASSERT (p4est == NULL);
  P4EST_ASSERT (conn == NULL);
  return 0;
}

/* Function for testing p4est-mesh for multiple trees in a non-brick scenario
 *
 * \param [in] p4est     The forest.
 * \param [in] conn      The connectivity structure
 * \param [in] periodic  Flag for checking if we have periodic boundaries
 * \returns 0 for success, -1 for failure
 */
int
test_mesh_multiple_trees_nonbrick (p4est_t * p4est,
                                   p4est_connectivity_t *
                                   conn, int8_t periodic, sc_MPI_Comm mpicomm)
{
  return 0;
}

int
main (int argc, char **argv)
{
  sc_MPI_Comm         mpicomm;
  int                 mpiret;
  int                 mpisize, mpirank;
  p4est_t            *p4est;
  p4est_connectivity_t *conn;
  int8_t              periodic_boundaries;
  /* initialize MPI */
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  mpicomm = sc_MPI_COMM_WORLD;
  mpiret = sc_MPI_Comm_size (mpicomm, &mpisize);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (mpicomm, &mpirank);
  SC_CHECK_MPI (mpiret);
  sc_init (mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);
  p4est = 0;
  conn = 0;

  int                 test_single, test_multi_brick, test_multi_non_brick;
  test_single = 1;
  test_multi_brick = 0;
  test_multi_non_brick = 0;

  /* test both periodic and non-periodic boundaries */
  if (test_single) {
    /* test one tree */
    periodic_boundaries = 0;
    test_mesh_one_tree (p4est, conn, periodic_boundaries, mpicomm);

    periodic_boundaries = 1;
    test_mesh_one_tree (p4est, conn, periodic_boundaries, mpicomm);
  }

  if (test_multi_brick) {
    /* test multiple trees; brick */
    periodic_boundaries = 0;
    test_mesh_multiple_trees_brick (p4est, conn, periodic_boundaries,
                                    mpicomm);
    periodic_boundaries = 1;
    test_mesh_multiple_trees_brick (p4est, conn, periodic_boundaries,
                                    mpicomm);
  }

  if (test_multi_non_brick) {
    /* test multiple trees; non-brick */
    periodic_boundaries = 0;
    test_mesh_multiple_trees_nonbrick (p4est, conn,
                                       periodic_boundaries, mpicomm);
    periodic_boundaries = 1;
    test_mesh_multiple_trees_nonbrick (p4est, conn,
                                       periodic_boundaries, mpicomm);
  }
  /* exit */
  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;
}
