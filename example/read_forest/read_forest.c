/*
  This file is part of p4est.
  p4est is a C library to manage a parallel collection of quadtrees and/or octrees.

  Copyright (C) 2007 Carsten Burstedde, Lucas Wilcox.

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

#include <p4est.h>
#include <p4est_file.h>
#include <p4est_base.h>

int
main (void)
{
  int                 retval;

  p4est_connectivity_t *connectivity;

  retval = p4est_connectivity_read ("mesh.p4t", &connectivity);
  P4EST_CHECK_ABORT (!retval, "Unable to read the mesh file.");

  p4est_connectivity_print (connectivity);

  p4est_connectivity_destroy (connectivity);

  return 0;
}

/* EOF read_forest.c */
