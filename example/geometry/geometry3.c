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

static void
write_vtk (p8est_t * p8est, p8est_geometry_t * geom, const char *name)
{
  p8est_vtk_write_file (p8est, geom, name);
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

  geye = p8est_geometry_new_identity ();
  gshell = p8est_geometry_new_shell (1., 0.55);

  conn = p8est_connectivity_new_unitcube ();
  p8est = p8est_new (mpicomm, conn, 0, 0, NULL, NULL);
  write_vtk (p8est, NULL, "unitcube_none");
  write_vtk (p8est, geye, "unitcube_identity");
  p8est_destroy (p8est);
  p8est_connectivity_destroy (conn);

  conn = p8est_connectivity_new_shell ();
  p8est = p8est_new (mpicomm, conn, 0, 0, NULL, NULL);
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
