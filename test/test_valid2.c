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

#ifndef P4_TO_P8
#include <p4est_bits.h>
#include <p4est_ghost.h>
#include <p4est_mesh.h>
#include <p4est_vtk.h>
#else
#include <p8est_bits.h>
#include <p8est_ghost.h>
#include <p8est_mesh.h>
#include <p8est_trilinear.h>
#include <p8est_vtk.h>
#endif

#ifndef P4_TO_P8
static const int    refine_level = 5;
#else
static const int    refine_level = 4;
#endif

static int
refine_fn (p4est_t * p4est, p4est_topidx_t which_tree,
           p4est_quadrant_t * quadrant)
{
  int                 cid;
  int                 endlevel = refine_level + (int) (which_tree % 3);

  if ((int) quadrant->level >= endlevel)
    return 0;

  cid = p4est_quadrant_child_id (quadrant);
  if (cid == 0 || cid == 3
#ifdef P4_TO_P8
      || cid == 5 || cid == 6
#endif
    ) {
    return 1;
  }

  return 0;
}

static int
coarsen_fn (p4est_t * p4est, p4est_topidx_t which_tree,
            p4est_quadrant_t * q0, p4est_quadrant_t * q1,
            p4est_quadrant_t * q2, p4est_quadrant_t * q3
#ifdef P4_TO_P8
            ,
            p4est_quadrant_t * q4, p4est_quadrant_t * q5,
            p4est_quadrant_t * q6, p4est_quadrant_t * q7
#endif
  )
{
  int                 pid;
  p4est_quadrant_t    p;

  if (q0->level <= 2)
    return 0;

  p4est_quadrant_parent (q0, &p);
  pid = p4est_quadrant_child_id (&p);

  return pid == 3;
}

static void
check_all (MPI_Comm mpicomm, p4est_connectivity_t * conn, const char *vtkname)
{
  p4est_t            *p4est;
  p4est_nodes_t      *nodes;
#ifdef P4_TO_P8
  trilinear_mesh_t   *mesh;
#endif
  sc_array_t          ghost_layer;

  p4est = p4est_new (mpicomm, conn, 0, NULL);
  p4est_refine (p4est, refine_fn, NULL);
  p4est_coarsen (p4est, coarsen_fn, NULL);
  p4est_balance (p4est, NULL);
  p4est_partition (p4est, NULL);
  p4est_vtk_write_file (p4est, vtkname);

  sc_array_init (&ghost_layer, sizeof (p4est_quadrant_t));
  p4est_build_ghost_layer (p4est, &ghost_layer);
  nodes = p4est_nodes_new (p4est, &ghost_layer);
#ifdef P4_TO_P8
  mesh = trilinear_mesh_extract (p4est, nodes);
  trilinear_mesh_delete (mesh);
#endif
  p4est_nodes_destroy (nodes);
  sc_array_reset (&ghost_layer);

  p4est_destroy (p4est);
  p4est_connectivity_destroy (conn);
}

int
main (int argc, char **argv)
{
  MPI_Comm            mpicomm;
  int                 mpiret;
  int                 size, rank;

  mpiret = MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  mpicomm = MPI_COMM_WORLD;
  mpiret = MPI_Comm_size (mpicomm, &size);
  SC_CHECK_MPI (mpiret);
  mpiret = MPI_Comm_rank (mpicomm, &rank);
  SC_CHECK_MPI (mpiret);

  sc_init (rank, sc_generic_abort_fn, &mpicomm, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  p4est_initial_quadrants_per_processor = 0;

#ifndef P4_TO_P8
  check_all (mpicomm, p4est_connectivity_new_unitsquare (),
             "test_unitsquare");
  check_all (mpicomm, p4est_connectivity_new_corner (), "test_corner");
  check_all (mpicomm, p4est_connectivity_new_moebius (), "test_moebius");
  check_all (mpicomm, p4est_connectivity_new_star (), "test_star");
  check_all (mpicomm, p4est_connectivity_new_periodic (), "test_periodic2");
#else
  check_all (mpicomm, p8est_connectivity_new_unitcube (), "test_unitcube");
  check_all (mpicomm, p8est_connectivity_new_periodic (), "test_periodic3");
  check_all (mpicomm, p8est_connectivity_new_twocubes (), "test_twocubes");
  check_all (mpicomm, p8est_connectivity_new_rotcubes (), "test_rotcubes");
#endif

  /* clean up and exit */
  sc_finalize ();

  mpiret = MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}

/* EOF test_valid2.c */
