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

#ifndef __P4EST_BASE_H__
#define __P4EST_BASE_H__

/*
 * this header is the only one that includes p4est_config.h
 * it is not installed in the final include directory
 */

#ifdef HAVE_CONFIG_H
#include <p4est_config.h>
#endif

#ifdef HAVE_CTYPE_H
#include <ctype.h>
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

#if(defined HAVE_BACKTRACE && defined HAVE_BACKTRACE_SYMBOLS)
#define P4EST_BACKTRACE
#endif

#define P4EST_NOOP do { ; } while (0)
#define P4EST_CHECK_ABORT(c,s)                     \
  do {                                             \
    if (!(c)) {                                    \
      fprintf (stderr, "Abort: %s\n   in %s:%d\n", \
               (s), __FILE__, __LINE__);           \
      p4est_abort ();                              \
    }                                              \
  } while (0)
#define P4EST_CHECK_ALLOC(p) P4EST_CHECK_ABORT (((p) != NULL), "Allocation")
#define P4EST_CHECK_REALLOC(p,n) P4EST_CHECK_ABORT (((p) != NULL || (n) == 0), \
                                                    "Allocation")
#define P4EST_CHECK_MPI(r) P4EST_CHECK_ABORT ((r) == MPI_SUCCESS, "MPI operation")
#ifdef P4EST_HAVE_DEBUG
#define P4EST_ASSERT(c) P4EST_CHECK_ABORT ((c), "Assertion '" #c "'")
#else
#define P4EST_ASSERT(c) P4EST_NOOP
#endif
#define P4EST_ASSERT_NOT_REACHED() P4EST_CHECK_ABORT (0, "Unreachable code")

#define P4EST_ALLOC(t,n) (t *) p4est_malloc ((n) * sizeof(t))
#define P4EST_ALLOC_ZERO(t,n) (t *) p4est_calloc ((n), sizeof(t))
#define P4EST_REALLOC(p,t,n) (t *) p4est_realloc ((p), (n) * sizeof(t))

/* it is allowed to call P4EST_FREE (NULL) which does nothing. */
#define P4EST_FREE(p) p4est_free (p)

/* min and max helper macros */
#define P4EST_MIN(a,b) (((a) < (b)) ? (a) : (b))
#define P4EST_MAX(a,b) (((a) > (b)) ? (a) : (b))

/* hopefully fast binary logarithms and binary round up */
#define P4EST_LOG2_8(x) (p4est_log_lookup_table[(x)])
#define P4EST_LOG2_16(x) (((x) > 0xff) ? \
                          (P4EST_LOG2_8 ((x) >> 8) + 8) : P4EST_LOG2_8 (x))
#define P4EST_LOG2_32(x) (((x) > 0xffff) ? \
                          (P4EST_LOG2_16 ((x) >> 16)) + 16 : P4EST_LOG2_16 (x))
#define P4EST_ROUNDUP2_32(x) \
                    (((x) <= 0) ? 0 : (1 << (P4EST_LOG2_32 ((x) - 1) + 1)))

typedef void        (*p4est_handler_t) (void *data);

extern const int    p4est_log_lookup_table[256];

int                 p4est_int32_compare (const void *v1, const void *v2);

void               *p4est_malloc (size_t size);
void               *p4est_calloc (size_t nmemb, size_t size);
void               *p4est_realloc (void *ptr, size_t size);
void                p4est_free (void *ptr);
void                p4est_memory_check (void);

/** Sets stream to line buffered. Must be called before writing anything. */
void                p4est_set_linebuffered (FILE * stream);

/** Installs an abort handler and catches signals INT, SEGV, USR2. */
void                p4est_set_abort_handler (int identifier,
                                             p4est_handler_t handler,
                                             void *data);

/** Prints a stack trace, calls the abort handler and terminates. */
void                p4est_abort (void);

#endif /* !__P4EST_BASE_H__ */
