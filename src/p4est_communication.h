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

#ifndef __P4EST_COMMUNICATION_H__
#define __P4EST_COMMUNICATION_H__

#include <p4est.h>

/** Compute the global number of quadrants.
 * \parm [in,out] p4est  Adds all p4est->local_num_quadrant counters.
 */
void                p4est_comm_count_quadrants (p4est_t * p4est);

#endif /* !__P4EST_COMMUNICATION_H__ */
