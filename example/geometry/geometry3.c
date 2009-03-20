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
  const double       *v = p8est->connectivity->vertices;
  const double        intsize = 1.0 / P8EST_ROOT_LEN;
  int                 i, k;
  int                 zi, yi, xi;
  double              eta_x, eta_y, eta_z;
  double              xyz[3], XYZ[3];
  double              det1, det2, J[3][3];
  double             *double_data, *dptr[4];
  size_t              zz;
  p4est_topidx_t      jt, vt[P8EST_CHILDREN];
  p4est_locidx_t      ql, num_quads;
  p4est_qcoord_t      h;
  p8est_tree_t       *tree;
  p8est_quadrant_t   *quad;
  sc_array_t         *quadrants;

  P4EST_ASSERT (v != NULL);

  double_data = P4EST_ALLOC (double, 4 * Ntotal);
  dptr[0] = double_data;
  dptr[1] = double_data + Ntotal;
  dptr[2] = double_data + 2 * Ntotal;
  dptr[3] = double_data + 3 * Ntotal;

  ql = 0;
  for (jt = p8est->first_local_tree; jt <= p8est->last_local_tree; ++jt) {
    tree = p4est_array_index_topidx (p8est->trees, jt);
    quadrants = &tree->quadrants;
    num_quads = quadrants->elem_count;

    /* retrieve corners of the tree */
    for (k = 0; k < P8EST_CHILDREN; ++k)
      vt[k] = ttv[jt * P8EST_CHILDREN + k];

    /* loop over the elements in the tree and calculated vertex coordinates */
    for (zz = 0; zz < num_quads; ++ql, ++zz) {
      quad = sc_array_index (quadrants, zz);
      h = P8EST_QUADRANT_LEN (quad->level);
      k = 0;
      for (zi = 0; zi < 2; ++zi) {
        for (yi = 0; yi < 2; ++yi) {
          for (xi = 0; xi < 2; ++xi) {
            P4EST_ASSERT (0 <= k && k < P8EST_CHILDREN);
            eta_x = intsize * (quad->x + xi * h);
            eta_y = intsize * (quad->y + yi * h);
            eta_z = intsize * (quad->z + zi * h);

            for (i = 0; i < 3; ++i) {
              /* *INDENT-OFF* */
              xyz[i] =
          ((1. - eta_x) * ((1. - eta_y) * ((1. - eta_z) * v[3 * vt[0] + i] +
                                                 eta_z  * v[3 * vt[4] + i]) +
                                 eta_y  * ((1. - eta_z) * v[3 * vt[2] + i] +
                                                 eta_z  * v[3 * vt[6] + i])) +
                 eta_x  * ((1. - eta_y) * ((1. - eta_z) * v[3 * vt[1] + i] +
                                                 eta_z  * v[3 * vt[5] + i]) +
                                 eta_y  * ((1. - eta_z) * v[3 * vt[3] + i] +
                                                 eta_z  * v[3 * vt[7] + i])));
              /* *INDENT-ON* */
            }
            if (geom != NULL) {
              geom->X (geom, jt, xyz, XYZ);
              for (i = 0; i < 3; ++i) {
                *dptr[i]++ = XYZ[i];
              }
              det1 = geom->D (geom, jt, xyz);
              det2 = geom->J (geom, jt, xyz, J);
              SC_CHECK_ABORT (fabs ((det1 - det2) / SC_MAX (det1, det2)) <
                              1.e-8, "Determinant inconsistent");
              *dptr[3]++ = det1;
            }
            else {
              for (i = 0; i < 3; ++i) {
                *dptr[i]++ = xyz[i];
              }
              *dptr[3]++ = 1.;
            }
            ++k;
          }
        }
      }
    }
  }
  P4EST_ASSERT (ql == Ncells);
  P4EST_ASSERT (dptr[0] == double_data + Ntotal);
  P4EST_ASSERT (dptr[1] == double_data + 2 * Ntotal);
  P4EST_ASSERT (dptr[2] == double_data + 3 * Ntotal);
  P4EST_ASSERT (dptr[3] == double_data + 4 * Ntotal);

  p8est_vtk_write_all (p8est, geom, 4, 0, name,
                       "X", double_data, "Y", double_data + Ntotal,
                       "Z", double_data + 2 * Ntotal,
                       "D", double_data + 3 * Ntotal);

  P4EST_FREE (double_data);
}

int
main (int argc, char **argv)
{
  int                 mpiret;
  MPI_Comm            mpicomm;
  p8est_connectivity_t *conn;
  p8est_geometry_t   *geye, *gshell, *gsphere;
  p8est_t            *p8est;

  mpiret = MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  mpicomm = MPI_COMM_WORLD;
  sc_init (mpicomm, true, true, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  p8est_vtk_default_scale = 1.;
  p8est_vtk_default_wrap_rank = 16;

  geye = p8est_geometry_new_identity ();
  gshell = p8est_geometry_new_shell (1., 0.55);
  gsphere = p8est_geometry_new_sphere (1., 0.191728, 0.039856);

  P4EST_STATISTICS ("Unitcube\n");
  conn = p8est_connectivity_new_unitcube ();
  p8est = p8est_new (mpicomm, conn, 0, 0, NULL, NULL);
  p8est_refine (p8est, true, refine_fn, NULL);
  write_vtk (p8est, NULL, "unitcube_none");
  write_vtk (p8est, geye, "unitcube_identity");
  p8est_destroy (p8est);
  p8est_connectivity_destroy (conn);

  P4EST_STATISTICS ("Shell\n");
  conn = p8est_connectivity_new_shell ();
  p8est = p8est_new (mpicomm, conn, 0, 0, NULL, NULL);
  p8est_refine (p8est, true, refine_fn, NULL);
  write_vtk (p8est, NULL, "shell_none");
  write_vtk (p8est, geye, "shell_identity");
  write_vtk (p8est, gshell, "shell_shell");
  p8est_destroy (p8est);
  p8est_connectivity_destroy (conn);

  P4EST_STATISTICS ("Sphere\n");
  conn = p8est_connectivity_new_sphere ();
  p8est = p8est_new (mpicomm, conn, 0, 0, NULL, NULL);
  p8est_refine (p8est, true, refine_fn, NULL);
  write_vtk (p8est, NULL, "sphere_none");
  write_vtk (p8est, geye, "sphere_identity");
  write_vtk (p8est, gsphere, "sphere_sphere");
  p8est_destroy (p8est);
  p8est_connectivity_destroy (conn);

  P4EST_FREE (geye);
  P4EST_FREE (gshell);
  P4EST_FREE (gsphere);

  sc_finalize ();

  mpiret = MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
