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

/*
 * Usage: p4est_mesh <configuration> <level>
 *        possible configurations:
 *        o unit      Refinement on the unit square.
 *        o three     Refinement on a forest with three trees.
 *        o moebius   Refinement on a 5-tree Moebius band.
 *        o star      Refinement on a 6-tree star shaped domain.
 *        o periodic  Refinement on the unit square with all-periodic b.c.
 *        o rotwrap   Refinement on the unit square with weird periodic b.c.
 */

#ifndef P4_TO_P8
#include <p4est_bits.h>
#include <p4est_extended.h>
#include <p4est_mesh.h>
#include <p4est_vtk.h>
#else
#include <p8est_bits.h>
#include <p8est_extended.h>
#include <p8est_mesh.h>
#include <p8est_vtk.h>
#endif

typedef enum
{
  P4EST_CONFIG_NULL,
#ifndef P4_TO_P8
  P4EST_CONFIG_UNIT,
  P4EST_CONFIG_THREE,
  P4EST_CONFIG_MOEBIUS,
  P4EST_CONFIG_STAR,
  P4EST_CONFIG_PERIODIC,
  P4EST_CONFIG_ROTWRAP
#else
  P8EST_CONFIG_UNIT,
  P8EST_CONFIG_PERIODIC,
  P8EST_CONFIG_ROTWRAP,
  P8EST_CONFIG_TWOCUBES,
  P8EST_CONFIG_TWOWRAP,
  P8EST_CONFIG_ROTCUBES,
  P8EST_CONFIG_SHELL,
  P8EST_CONFIG_SPHERE
#endif
}
simple_config_t;

typedef struct
{
  p4est_quadrant_t    quad;
}
user_data_t;

typedef struct
{
  sc_MPI_Comm         mpicomm;
  int                 mpisize;
  int                 mpirank;
}
mpi_context_t;

static int          refine_level = 0;

static void
init_fn (p4est_t * p4est, p4est_topidx_t which_tree,
         p4est_quadrant_t * quadrant)
{
  user_data_t        *data = (user_data_t *) quadrant->p.user_data;

  data->quad = *quadrant;
  data->quad.p.which_tree = which_tree;
}

static int
refine_uniform (p4est_t * p4est, p4est_topidx_t which_tree,
                p4est_quadrant_t * quadrant)
{
  return (int) quadrant->level < refine_level;
}

static int
refine_normal (p4est_t * p4est, p4est_topidx_t which_tree,
               p4est_quadrant_t * quadrant)
{
  if ((int) quadrant->level >= (refine_level - (int) (which_tree % 3))) {
    return 0;
  }
  if (quadrant->level == 1 && p4est_quadrant_child_id (quadrant) == 3) {
    return 1;
  }
  if (quadrant->x == P4EST_LAST_OFFSET (2) &&
      quadrant->y == P4EST_LAST_OFFSET (2)) {
    return 1;
  }
#ifndef P4_TO_P8
  if (quadrant->x >= P4EST_QUADRANT_LEN (2)) {
    return 0;
  }
#else
  if (quadrant->z >= P8EST_QUADRANT_LEN (2)) {
    return 0;
  }
#endif

  return 1;
}

static void
test_mesh (p4est_t * p4est, p4est_ghost_t * ghost, p4est_mesh_t * mesh,
           int compute_tree_index, int compute_level_lists,
           p4est_connect_type_t mesh_btype,
           user_data_t * ghost_data, int uniform)
{
  const int           HF = P4EST_HALF * P4EST_FACES;
  size_t              i;
  int                 level;
  int                 f, nf;
  int                 c;
  int                 nface;
  int                 nrank;
  p4est_topidx_t      which_tree;
  p4est_locidx_t      K, kl;
  p4est_locidx_t      ql, QpG, lnC;
  p4est_locidx_t      qlid, qumid, quadrant_id, which_quad;
  p4est_mesh_face_neighbor_t mfn, mfn2;
  p4est_quadrant_t   *q;
  p4est_tree_t       *tree;

  K = mesh->local_num_quadrants;
  P4EST_ASSERT (K == p4est->local_num_quadrants);
  QpG = mesh->local_num_quadrants + mesh->ghost_num_quadrants;
  lnC = mesh->local_num_corners;
  P4EST_ASSERT (lnC >= 0);

  P4EST_ASSERT (compute_tree_index == (mesh->quad_to_tree != NULL));
  P4EST_ASSERT (compute_level_lists == (mesh->quad_level != NULL));
  P4EST_ASSERT ((mesh_btype == P4EST_CONNECT_CORNER) ==
                (mesh->quad_to_corner != NULL));

  /* TODO: test the mesh relations in more depth */
  tree = NULL;
  for (kl = 0; kl < K; ++kl) {
    if (compute_tree_index) {
      tree = p4est_tree_array_index (p4est->trees, mesh->quad_to_tree[kl]);
      SC_CHECK_ABORTF (tree->quadrants_offset <= kl && kl <
                       tree->quadrants_offset +
                       (p4est_locidx_t) tree->quadrants.elem_count,
                       "Tree index mismatch %lld", (long long) kl);
    }

    if (mesh_btype == P4EST_CONNECT_CORNER) {
      for (c = 0; c < P4EST_CHILDREN; ++c) {
        qlid = mesh->quad_to_corner[P4EST_CHILDREN * kl + c];
        SC_CHECK_ABORTF (qlid >= -2
                         && qlid < QpG + lnC, "quad %lld corner %d mismatch",
                         (long long) kl, c);
      }
    }
    for (f = 0; f < P4EST_FACES; ++f) {
      ql = mesh->quad_to_quad[P4EST_FACES * kl + f];
      SC_CHECK_ABORTF (0 <= ql && ql < QpG,
                       "quad %d face %d neighbor %d mismatch", kl, f, ql);
      nf = mesh->quad_to_face[P4EST_FACES * kl + f];
      if (uniform) {
        SC_CHECK_ABORTF (0 <= nf && nf < HF,
                         "quad %d face %d code %d mismatch", kl, f, nf);
      }
      else {
        SC_CHECK_ABORTF (-HF <= nf && nf < (P4EST_HALF + 1) * HF,
                         "quad %d face %d code %d mismatch", kl, f, nf);
      }
    }
  }

  /* Test the level lists */
  if (compute_tree_index && compute_level_lists) {
    for (level = 0; level < P4EST_QMAXLEVEL; ++level) {
      for (i = 0; i < mesh->quad_level[level].elem_count; ++i) {
        /* get the local quadrant id */
        quadrant_id =
          *(p4est_locidx_t *) sc_array_index (&mesh->quad_level[level], i);

        /* get the tree it belongs to */
        kl = mesh->quad_to_tree[quadrant_id];
        tree = p4est_tree_array_index (p4est->trees, kl);

        /* and finally, get the actual quadrant from the tree quadrant list */
        quadrant_id -= tree->quadrants_offset;
        q =
          p4est_quadrant_array_index (&tree->quadrants, (size_t) quadrant_id);

        SC_CHECK_ABORTF (q->level == level,
                         "quad %d level %d mismatch", quadrant_id, level);
      }
    }
  }

  /* Test face neighbor iterator */
  for (qumid = 0; qumid < mesh->local_num_quadrants; ++qumid) {
    which_tree = -1;
    q = p4est_mesh_quadrant_cumulative (p4est, qumid,
                                        &which_tree, &quadrant_id);
    p4est_mesh_face_neighbor_init2 (&mfn, p4est, ghost, mesh,
                                    which_tree, quadrant_id);
    p4est_mesh_face_neighbor_init (&mfn2, p4est, ghost, mesh, which_tree, q);
    P4EST_ASSERT (mfn2.quadrant_id == quadrant_id);
    while ((q = p4est_mesh_face_neighbor_next (&mfn, &which_tree, &which_quad,
                                               &nface, &nrank)) != NULL) {
#ifdef P4EST_ENABLE_DEBUG
      user_data_t        *data;

      data = (user_data_t *) p4est_mesh_face_neighbor_data (&mfn, ghost_data);

      P4EST_ASSERT (p4est_quadrant_is_equal (q, &(data->quad)));
      P4EST_ASSERT (data->quad.p.which_tree == which_tree);
#endif
    }
  }
}

static void
mesh_run (mpi_context_t * mpi, p4est_connectivity_t * connectivity,
          int uniform, int compute_tree_index, int compute_level_lists,
          p4est_connect_type_t mesh_btype)
{
  int                 mpiret;
  unsigned            crc;
  long                local_used[4], global_used[4];
  p4est_t            *p4est;
  p4est_ghost_t      *ghost;
  p4est_mesh_t       *mesh;
  user_data_t        *ghost_data;

  p4est = p4est_new (mpi->mpicomm, connectivity,
                     sizeof (user_data_t), init_fn, NULL);
  if (!uniform)
    p4est_vtk_write_file (p4est, NULL, P4EST_STRING "_mesh_new");

  /* refinement */
  if (uniform) {
    p4est_refine (p4est, 1, refine_uniform, init_fn);
  }
  else {
    p4est_refine (p4est, 1, refine_normal, init_fn);
    p4est_vtk_write_file (p4est, NULL, P4EST_STRING "_mesh_refined");
  }

  /* balance */
  p4est_balance (p4est, P4EST_CONNECT_FULL, init_fn);
  if (!uniform)
    p4est_vtk_write_file (p4est, NULL, P4EST_STRING "_mesh_balanced");

  /* partition */
  p4est_partition (p4est, 0, NULL);
  if (!uniform) {
    p4est_vtk_write_file (p4est, NULL, P4EST_STRING "_mesh_partition");
  }
  crc = p4est_checksum (p4est);

  /* print and verify forest checksum */
  P4EST_GLOBAL_STATISTICSF ("Tree %s checksum 0x%08x\n",
                            uniform ? "uniform" : "adapted", crc);

  /* create ghost layer and mesh */
  ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FULL);
  ghost_data = P4EST_ALLOC (user_data_t, ghost->ghosts.elem_count);
  p4est_ghost_exchange_data (p4est, ghost, ghost_data);
  mesh = p4est_mesh_new_ext (p4est, ghost,
                             compute_tree_index, compute_level_lists,
                             mesh_btype);
  test_mesh (p4est, ghost, mesh,
             compute_tree_index, compute_level_lists, mesh_btype,
             ghost_data, uniform);

  /* compute memory used */
  local_used[0] = (long) p4est_connectivity_memory_used (p4est->connectivity);
  local_used[1] = (long) p4est_memory_used (p4est);
  local_used[2] = (long) p4est_ghost_memory_used (ghost);
  local_used[3] = (long) p4est_mesh_memory_used (mesh);
  mpiret = sc_MPI_Allreduce (local_used, global_used, 4, sc_MPI_LONG,
                             sc_MPI_SUM, mpi->mpicomm);
  SC_CHECK_MPI (mpiret);
  P4EST_GLOBAL_PRODUCTIONF ("Total %s memory used %ld %ld %ld %ld\n",
                            uniform ? "uniform" : "adapted",
                            global_used[0], global_used[1],
                            global_used[2], global_used[3]);

  /* destroy ghost layer and mesh */
  P4EST_FREE (ghost_data);
  p4est_mesh_destroy (mesh);
  p4est_ghost_destroy (ghost);

  /* destroy the p4est structure */
  p4est_destroy (p4est);
}

int
main (int argc, char **argv)
{
  int                 mpiret;
  int                 wrongusage;
  const char         *usage;
  mpi_context_t       mpi_context, *mpi = &mpi_context;
  p4est_connectivity_t *connectivity;
  simple_config_t     config;

  /* initialize MPI and p4est internals */
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  mpi->mpicomm = sc_MPI_COMM_WORLD;
  mpiret = sc_MPI_Comm_size (mpi->mpicomm, &mpi->mpisize);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (mpi->mpicomm, &mpi->mpirank);
  SC_CHECK_MPI (mpiret);

  sc_init (mpi->mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  /* process command line arguments */
  usage =
    "Arguments: <configuration> <level>\n   Configuration can be any of\n"
#ifndef P4_TO_P8
    "      unit|three|moebius|star|periodic|rotwrap\n"
#else
    "      unit|periodic|rotwrap|twocubes|twowrap|rotcubes|shell|sphere\n"
#endif
    "   Level controls the maximum depth of refinement\n";
  wrongusage = 0;
  config = P4EST_CONFIG_NULL;
  if (!wrongusage && argc != 3) {
    wrongusage = 1;
  }
  if (!wrongusage) {
    if (!strcmp (argv[1], "unit")) {
#ifndef P4_TO_P8
      config = P4EST_CONFIG_UNIT;
#else
      config = P8EST_CONFIG_UNIT;
#endif
    }
#ifndef P4_TO_P8
    else if (!strcmp (argv[1], "three")) {
      config = P4EST_CONFIG_THREE;
    }
    else if (!strcmp (argv[1], "moebius")) {
      config = P4EST_CONFIG_MOEBIUS;
    }
    else if (!strcmp (argv[1], "star")) {
      config = P4EST_CONFIG_STAR;
    }
    else if (!strcmp (argv[1], "periodic")) {
      config = P4EST_CONFIG_PERIODIC;
    }
    else if (!strcmp (argv[1], "rotwrap")) {
      config = P4EST_CONFIG_ROTWRAP;
    }
#else
    else if (!strcmp (argv[1], "periodic")) {
      config = P8EST_CONFIG_PERIODIC;
    }
    else if (!strcmp (argv[1], "rotwrap")) {
      config = P8EST_CONFIG_ROTWRAP;
    }
    else if (!strcmp (argv[1], "twocubes")) {
      config = P8EST_CONFIG_TWOCUBES;
    }
    else if (!strcmp (argv[1], "twowrap")) {
      config = P8EST_CONFIG_TWOWRAP;
    }
    else if (!strcmp (argv[1], "rotcubes")) {
      config = P8EST_CONFIG_ROTCUBES;
    }
    else if (!strcmp (argv[1], "shell")) {
      config = P8EST_CONFIG_SHELL;
    }
    else if (!strcmp (argv[1], "sphere")) {
      config = P8EST_CONFIG_SPHERE;
    }
#endif
    else {
      wrongusage = 1;
    }
  }
  if (wrongusage) {
    P4EST_GLOBAL_LERROR (usage);
    sc_abort_collective ("Usage error");
  }

  /* assign variables based on configuration */
  refine_level = atoi (argv[2]);

  /* create connectivity and forest structures */
  if (0) {
  }
#ifndef P4_TO_P8
  else if (config == P4EST_CONFIG_THREE) {
    connectivity = p4est_connectivity_new_corner ();
  }
  else if (config == P4EST_CONFIG_MOEBIUS) {
    connectivity = p4est_connectivity_new_moebius ();
  }
  else if (config == P4EST_CONFIG_STAR) {
    connectivity = p4est_connectivity_new_star ();
  }
  else if (config == P4EST_CONFIG_PERIODIC) {
    connectivity = p4est_connectivity_new_periodic ();
  }
  else if (config == P4EST_CONFIG_ROTWRAP) {
    connectivity = p4est_connectivity_new_rotwrap ();
  }
#else
  else if (config == P8EST_CONFIG_PERIODIC) {
    connectivity = p8est_connectivity_new_periodic ();
  }
  else if (config == P8EST_CONFIG_ROTWRAP) {
    connectivity = p8est_connectivity_new_rotwrap ();
  }
  else if (config == P8EST_CONFIG_TWOCUBES) {
    connectivity = p8est_connectivity_new_twocubes ();
  }
  else if (config == P8EST_CONFIG_TWOWRAP) {
    connectivity = p8est_connectivity_new_twowrap ();
  }
  else if (config == P8EST_CONFIG_ROTCUBES) {
    connectivity = p8est_connectivity_new_rotcubes ();
  }
  else if (config == P8EST_CONFIG_SHELL) {
    connectivity = p8est_connectivity_new_shell ();
  }
  else if (config == P8EST_CONFIG_SPHERE) {
    connectivity = p8est_connectivity_new_sphere ();
  }
#endif
  else {
#ifndef P4_TO_P8
    connectivity = p4est_connectivity_new_unitsquare ();
#else
    connectivity = p8est_connectivity_new_unitcube ();
#endif
  }

  /* run mesh tests */
  mesh_run (mpi, connectivity, 1, 0, 1, P4EST_CONNECT_FULL);
  mesh_run (mpi, connectivity, 0, 1, 0, P4EST_CONNECT_FULL);
  mesh_run (mpi, connectivity, 0, 0, 0, P4EST_CONNECT_FACE);
  mesh_run (mpi, connectivity, 1, 1, 1, P4EST_CONNECT_FACE);

  /* clean up and exit */
  p4est_connectivity_destroy (connectivity);
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
