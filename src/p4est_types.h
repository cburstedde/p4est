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

#ifndef P4EST_TYPES_H
#define P4EST_TYPES_H

/* include p4est config header */

#include <p4est_config.h>

/* sc.h includes several headers with proper checks */

#include <sc.h>
#if 0                           /* if present then already included by sc.h */
#include <math.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <libgen.h>
#endif
/* also includes getopt, obstack, zlib and MPI */

/** Typedef for quadrant coordinates */
typedef int32_t     p4est_qcoord_t;
#define P4EST_MPI_QCOORD MPI_INT

/** Typedef for processor-local indexing */
typedef int32_t     p4est_locidx_t;
#define P4EST_MPI_LOCIDX MPI_INT

/** Typedef for globally unique indexing */
typedef int64_t     p4est_gloidx_t;
#define P4EST_MPI_GLOIDX MPI_LONG_LONG

#endif /* !P4EST_TYPES_H */
