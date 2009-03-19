/*
  This file is part of p4est.
  p4est is a C library to manage a parallel collection of quadtrees and/or
  octrees.

  Copyright (C) 2007-2009 Carsten Burstedde, Lucas Wilcox.

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

#include <p8est_vtk.h>

static int
refine_fn (p8est_t * p8est, p4est_topidx_t which_tree,
           p8est_quadrant_t * quadrant)
{
  const p4est_qcoord_t qh = P8EST_QUADRANT_LEN (quadrant->level);
  const p4est_qcoord_t rx = quadrant->x + qh;
  const p4est_qcoord_t ry = quadrant->y + qh;
  const p4est_qcoord_t rz = quadrant->z + qh;

  if (rx == P8EST_ROOT_LEN && quadrant->y == 0 && quadrant->z == 0
      && quadrant->level < 2)
    return 1;

  if (quadrant->x == 0 && ry == P8EST_ROOT_LEN && quadrant->z == 0
      && quadrant->level < 3)
    return 1;

  if (quadrant->x == 0 && quadrant->y == 0 && rz == P8EST_ROOT_LEN
      && quadrant->level < 4)
    return 1;

  return 0;
}

static void
write_vtk (p8est_t * p8est, p8est_geometry_t * geom, const char *name)
{
  const p4est_topidx_t *ttv = p8est->connectivity->tree_to_vertex;
  const p4est_locidx_t Ncells = p8est->local_num_quadrants;
  const p4est_locidx_t Ntotal = P8EST_CHILDREN * Ncells;        /* type ok */
  const double       *vertices = p8est->connectivity->vertices;
  int                 i, k;
  double             *xyz, center[3], transformed[3];
  size_t              zz;
  p4est_topidx_t      jt, vt;
  p4est_locidx_t      ql, tree_quads;
  p8est_tree_t       *tree;
  sc_array_t         *quadrants;

  P4EST_ASSERT (vertices != NULL);

  xyz = P4EST_ALLOC (double, 3 * Ntotal);

  ql = 0;
  for (jt = p8est->first_local_tree; jt <= p8est->last_local_tree; ++jt) {
    tree = p4est_array_index_topidx (p8est->trees, jt);
    for (i = 0; i < 3; ++i) {
      center[i] = 0.;
      for (k = 0; k < P8EST_CHILDREN; ++k) {
        vt = ttv[jt * P8EST_CHILDREN + k];
        center[i] += vertices[3 * vt + i] / P8EST_CHILDREN;
      }
    }
    if (geom != NULL) {
      geom->X (geom, jt, center, transformed);
    }
    else {
      memcpy (transformed, center, 3 * sizeof (double));
    }

    quadrants = &tree->quadrants;
    tree_quads = quadrants->elem_count;
    for (zz = 0; zz < tree_quads; ++ql, ++zz) {
      for (k = 0; k < 8; ++k) {
        xyz[8 * ql + k] = transformed[0];
        xyz[Ntotal + 8 * ql + k] = transformed[1];
        xyz[2 * Ntotal + 8 * ql + k] = transformed[2];
      }
    }
  }
  P4EST_ASSERT (ql == Ncells);

  p8est_vtk_write_all (p8est, geom, 3, 0, name,
                       "X", xyz, "Y", xyz + Ntotal, "Z", xyz + 2 * Ntotal);

  P4EST_FREE (xyz);
}

int
main (int argc, char **argv)
{
  int                 mpiret;
  MPI_Comm            mpicomm;
  p8est_connectivity_t *conn;
  p8est_geometry_t   *geye;
  p8est_geometry_t   *gshell;
  p8est_t            *p8est;

  mpiret = MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  mpicomm = MPI_COMM_WORLD;
  sc_init (mpicomm, true, true, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  p8est_vtk_default_wrap_rank = 16;

  geye = p8est_geometry_new_identity ();
  gshell = p8est_geometry_new_shell (1., 0.55);

  conn = p8est_connectivity_new_unitcube ();
  p8est = p8est_new (mpicomm, conn, 0, 0, NULL, NULL);
  p8est_refine (p8est, true, refine_fn, NULL);
  write_vtk (p8est, NULL, "unitcube_none");
  write_vtk (p8est, geye, "unitcube_identity");
  p8est_destroy (p8est);
  p8est_connectivity_destroy (conn);

  conn = p8est_connectivity_new_shell ();
  p8est = p8est_new (mpicomm, conn, 0, 0, NULL, NULL);
  p8est_refine (p8est, true, refine_fn, NULL);
  write_vtk (p8est, NULL, "shell_none");
  write_vtk (p8est, geye, "shell_identity");
  write_vtk (p8est, gshell, "shell_shell");
  p8est_destroy (p8est);
  p8est_connectivity_destroy (conn);

  P4EST_FREE (geye);
  P4EST_FREE (gshell);

  sc_finalize ();

  mpiret = MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
