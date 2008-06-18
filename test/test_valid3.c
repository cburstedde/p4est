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

#include <p8est_vtk.h>

static void
check_all (MPI_Comm mpicomm, p8est_connectivity_t * conn, const char *vtkname)
{
  p8est_t            *p8est;

  p8est = p8est_new (mpicomm, conn, 0, NULL);
  p8est_vtk_write_file (p8est, vtkname);
  p8est_destroy (p8est);
  p8est_connectivity_destroy (conn);
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

  p8est_initial_quadrants_per_processor = 0;

  check_all (mpicomm, p8est_connectivity_new_unitcube (), "test_unitcube");
  check_all (mpicomm, p8est_connectivity_new_periodic (), "test_periodic");
  check_all (mpicomm, p8est_connectivity_new_rotcubes (), "test_rotcubes");

  /* clean up and exit */
  sc_finalize ();

  mpiret = MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}

/* EOF test_valid3.c */
