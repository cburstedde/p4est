/*
  This file is part of p4est.
  p4est is a C library to manage a parallel collection of quadtrees and/or
  octrees.

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

#include <p4est_communication.h>
#include <p4est_base.h>

void
p4est_comm_count_quadrants (p4est_t * p4est)
{
#ifdef HAVE_MPI
  int                 mpiret;
#endif
  int64_t             qlocal;

  qlocal = p4est->local_num_quadrants;
  p4est->global_num_quadrants = qlocal;

#ifdef HAVE_MPI
  if (p4est->mpicomm != MPI_COMM_NULL) {
    mpiret = MPI_Allreduce (&qlocal, &p4est->global_num_quadrants,
                            1, MPI_LONG_LONG, MPI_SUM, p4est->mpicomm);
    P4EST_CHECK_MPI (mpiret);
  }
#endif
}

/* EOF p4est_communication.h */
