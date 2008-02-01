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
 *
 * do not include this header from any other .h file
 * do include this header topmost in any .c file
 */

#ifdef HAVE_CONFIG_H
#include <p4est_config.h>
#endif

#ifdef HAVE_CTYPE_H
#include <ctype.h>
#endif

#ifdef HAVE_STDARG_H
#include <stdarg.h>
#endif

#ifdef HAVE_STDINT_H
#include <stdint.h>
#endif

#ifdef HAVE_STDDEF_H
#include <stddef.h>
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

/* needs the include p4est_config.h directive from above */
#include <p4est_log.h>

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

/* log macros, for priorities see p4est_log.h */
#define P4EST_GLOBAL_LOG(p,s) \
  _LOG_PRE (p4est_log_category_global, (p), (s)) _LOG_POST
#define P4EST_GLOBAL_LOGF(p,f,...) \
  _LOG_PRE (p4est_log_category_global, (p), (f)) , __VA_ARGS__ _LOG_POST
#define P4EST_LOG(p,s) \
  _LOG_PRE (p4est_log_category_rank, (p), (s)) _LOG_POST
#define P4EST_LOGF(p,f,...) \
  _LOG_PRE (p4est_log_category_rank, (p), (f)) , __VA_ARGS__ _LOG_POST

/* convenience global log macros will only print if identifier <= 0 */
#define P4EST_GLOBAL_TRACE(s) P4EST_GLOBAL_LOG (P4EST_LP_TRACE, (s))
#define P4EST_GLOBAL_TRACEF(f,...) \
  P4EST_GLOBAL_LOGF (P4EST_LP_TRACE, (f), __VA_ARGS__)
#define P4EST_GLOBAL_DEBUG(s) P4EST_GLOBAL_LOG (P4EST_LP_DEBUG, (s))
#define P4EST_GLOBAL_DEBUGF(f,...) \
  P4EST_GLOBAL_LOGF (P4EST_LP_DEBUG, (f), __VA_ARGS__)
#define P4EST_GLOBAL_VERBOSE(s) P4EST_GLOBAL_LOG (P4EST_LP_VERBOSE, (s))
#define P4EST_GLOBAL_VERBOSEF(f,...) \
  P4EST_GLOBAL_LOGF (P4EST_LP_VERBOSE, (f), __VA_ARGS__)
#define P4EST_GLOBAL_INFO(s) P4EST_GLOBAL_LOG (P4EST_LP_INFO, (s))
#define P4EST_GLOBAL_INFOF(f,...) \
  P4EST_GLOBAL_LOGF (P4EST_LP_INFO, (f), __VA_ARGS__)
#define P4EST_GLOBAL_STATISTICS(s) P4EST_GLOBAL_LOG (P4EST_LP_STATISTICS, (s))
#define P4EST_GLOBAL_STATISTICSF(f,...) \
  P4EST_GLOBAL_LOGF (P4EST_LP_STATISTICS, (f), __VA_ARGS__)
#define P4EST_GLOBAL_PRODUCTION(s) P4EST_GLOBAL_LOG (P4EST_LP_PRODUCTION, (s))
#define P4EST_GLOBAL_PRODUCTIONF(f,...) \
  P4EST_GLOBAL_LOGF (P4EST_LP_PRODUCTION, (f), __VA_ARGS__)

/* convenience log macros that are active on every processor */
#define P4EST_TRACE(s) P4EST_LOG (P4EST_LP_TRACE, (s))
#define P4EST_TRACEF(f,...) \
  P4EST_LOGF (P4EST_LP_TRACE, (f), __VA_ARGS__)
#define P4EST_DEBUG(s) P4EST_LOG (P4EST_LP_DEBUG, (s))
#define P4EST_DEBUGF(f,...) \
  P4EST_LOGF (P4EST_LP_DEBUG, (f), __VA_ARGS__)
#define P4EST_VERBOSE(s) P4EST_LOG (P4EST_LP_VERBOSE, (s))
#define P4EST_VERBOSEF(f,...) \
  P4EST_LOGF (P4EST_LP_VERBOSE, (f), __VA_ARGS__)
#define P4EST_INFO(s) P4EST_LOG (P4EST_LP_INFO, (s))
#define P4EST_INFOF(f,...) \
  P4EST_LOGF (P4EST_LP_INFO, (f), __VA_ARGS__)
#define P4EST_STATISTICS(s) P4EST_LOG (P4EST_LP_STATISTICS, (s))
#define P4EST_STATISTICSF(f,...) \
  P4EST_LOGF (P4EST_LP_STATISTICS, (f), __VA_ARGS__)
#define P4EST_PRODUCTION(s) P4EST_LOG (P4EST_LP_PRODUCTION, (s))
#define P4EST_PRODUCTIONF(f,...) \
  P4EST_LOGF (P4EST_LP_PRODUCTION, (f), __VA_ARGS__)

typedef void        (*p4est_handler_t) (void *data);

extern const int    p4est_log_lookup_table[256];
extern struct LogCategory p4est_log_category_global;
extern struct LogCategory p4est_log_category_rank;

void                p4est_set_log_threshold (int log_priority);

#if(0)
int                 p4est_int32_compare (const void *v1, const void *v2);
int                 p4est_int64_compare (const void *v1, const void *v2);
#endif

/** Find the lowest position k in a sorted array such that array[k] >= target.
 * \param [in]  target  The target lower bound to binary search for.
 * \param [in]  array   The 64bit integer array to binary search in.
 * \param [in]  size    The number of int64_t's in the array.
 * \param [in]  guess   Initial array position to look at.
 * \return  Returns the matching position, or -1 if array[size-1] < target.
 */
int                 p4est_int64_lower_bound (int64_t target,
                                             const int64_t * array,
                                             int size, int guess);

void               *p4est_malloc (size_t size);
void               *p4est_calloc (size_t nmemb, size_t size);
void               *p4est_realloc (void *ptr, size_t size);
void                p4est_free (void *ptr);
void                p4est_memory_check (void);

/** Initializes p4est_log_category to print to stream.
 * \param [in]  identifier  Set this to mpirank, or -1 to avoid number prefix.
 */
void                p4est_init_logging (FILE * stream, int identifier);

/** Installs an abort handler and catches signals INT, SEGV, USR2. */
void                p4est_set_abort_handler (p4est_handler_t handler,
                                             void *data);

/** Combines calls to init_logging and set_abort_handler. */
void                p4est_init (FILE * stream, int identifier,
                                p4est_handler_t abort_handler,
                                void *abort_data);

/** Prints a stack trace, calls the abort handler and terminates. */
void                p4est_abort (void)
  __attribute__ ((noreturn));

#endif /* !__P4EST_BASE_H__ */
