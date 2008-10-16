/*
  This file is part of p4est.
  p4est is a C library to manage a parallel collection of quadtrees and/or
  octrees.

  Copyright (C) 2008 Carsten Burstedde, Lucas Wilcox.

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

#include <p8est.h>
#include <p4est_to_p8est.h>
#include "p4est.c"

static int          p8est_uninitialized_key;
void               *P8EST_DATA_UNINITIALIZED = &p8est_uninitialized_key;
const int           p8est_num_ranges = 25;

int
p8est_balance_type_int (p8est_balance_type_t btype)
{
  switch (btype) {
  case P8EST_BALANCE_FACE:
    return 1;
  case P8EST_BALANCE_EDGE:
    return 2;
  case P8EST_BALANCE_CORNER:
    return 3;
  default:
    SC_CHECK_NOT_REACHED ();
  }
}

const char         *
p8est_balance_type_string (p8est_balance_type_t btype)
{
  switch (btype) {
  case P8EST_BALANCE_FACE:
    return "FACE";
  case P8EST_BALANCE_EDGE:
    return "EDGE";
  case P8EST_BALANCE_CORNER:
    return "CORNER";
  default:
    SC_CHECK_NOT_REACHED ();
  }
}

/* EOF p8est.c */
