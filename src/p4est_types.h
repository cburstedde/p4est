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

#ifndef __P4EST_TYPES_H__
#define __P4EST_TYPES_H__

/*
 * This header is included by p4est.h and p4est_base.h.
 * Do not include this header from any other .h file.
 *
 * This includes p4est_config.h if (!P4EST_CONFIG_INSTALLED).
 * For an installation of this library, the user may define
 *    P4EST_HAVE_DEBUG, P4EST_LOG_LEVEL if desired.
 *    Also, various HAVE_ macros need to be provided externally.
 */

/* this will be changed to 1 by make install */
#define P4EST_CONFIG_INSTALLED 0

/* this will be changed to 1 by make install if mpi is configured in */
#define P4EST_CONFIG_MPI 0

/* do some magic to avoid using p4est_config.h in the installed header */
#if (P4EST_CONFIG_INSTALLED)

#if (P4EST_CONFIG_MPI)
#include <mpi.h>
#else
#include <p4est_mpi_dummy.h>
#endif

#else /* !P4EST_CONFIG_INSTALLED */

#ifdef HAVE_CONFIG_H
#include <p4est_config.h>
#endif

#ifdef HAVE_MPI
#include <mpi.h>
#else
#include <p4est_mpi_dummy.h>
#endif

#ifdef P4EST_SPLINT
typedef void        (*sig_t) (int);
#define HAVE_ZLIB_H
#define HAVE_MPI
#define P4EST_HAVE_DEBUG
#endif /* P4EST_SPLINT */

#endif /* !P4EST_CONFIG_INSTALLED */

#if (defined HAVE_BACKTRACE && defined HAVE_BACKTRACE_SYMBOLS)
#define P4EST_BACKTRACE
#endif

#ifdef HAVE_CTYPE_H
#include <ctype.h>
#endif

#ifdef HAVE_STDARG_H
#include <stdarg.h>
#endif

#ifdef HAVE_STDDEF_H
#include <stddef.h>
#endif

#ifdef HAVE_STDINT_H
#include <stdint.h>
#endif

#ifdef HAVE_STDIO_H
#include <stdio.h>
#endif

#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif

#ifdef HAVE_STRING_H
#include <string.h>
#endif

/** Typedef for quadrant coordinates */
typedef int32_t     p4est_qcoord_t;
#define P4EST_MPI_QCOORD MPI_INT

/** Typedef for processor-local indexing */
typedef int32_t     p4est_locidx_t;
#define P4EST_MPI_LOCIDX MPI_INT

/** Typedef for globally unique indexing */
typedef int64_t     p4est_gloidx_t;
#define P4EST_MPI_GLOIDX MPI_LONG_LONG

#endif /* !__P4EST_TYPES_H__ */
