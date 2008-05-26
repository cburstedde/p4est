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

#ifndef P4EST_BASE_H
#define P4EST_BASE_H

#include <p4est_types.h>

/* use error checking from libsc */
#define P4EST_CHECK_MPI(r)              SC_CHECK_MPI (r)
#define P4EST_CHECK_ABORT(c,s)          SC_CHECK_ABORT (c, s)
#define P4EST_CHECK_ABORTF(c,fmt,...)   SC_CHECK_ABORTF (c, fmt, __VA_ARGS__)
#ifdef P4EST_DEBUG
#define P4EST_ASSERT(c) SC_CHECK_ABORT ((c), "Assertion '" #c "'")
#else
#define P4EST_ASSERT(c) SC_NOOP ()
#endif
#define P4EST_ASSERT_NOT_REACHED() SC_ASSERT_NOT_REACHED ()

/* use memory allocation from libsc */
#define P4EST_ALLOC(t,n)        SC_ALLOC (t, n)
#define P4EST_ALLOC_ZERO(t,n)   SC_ALLOC_ZERO (t, n)
#define P4EST_REALLOC(p,t,n)    SC_REALLOC (p, t, n)
#define P4EST_STRDUP(s)         SC_STRDUP (s)
#define P4EST_FREE(p)           SC_FREE (p)

/* use libsc log macros for now */
#define P4EST_GLOBAL_LOG(p,s) SC_LOG ((p), SC_LC_GLOBAL, (s))
#define P4EST_GLOBAL_LOGF(p,f,...) SC_LOGF ((p), SC_LC_GLOBAL, (f), __VA_ARGS__)
#define P4EST_NORMAL_LOG(p,s) SC_LOG ((p), SC_LC_NORMAL, (s))
#define P4EST_NORMAL_LOGF(p,f,...) SC_LOGF ((p), SC_LC_NORMAL, (f), __VA_ARGS__)

/* convenience global log macros will only print if identifier <= 0 */
#define P4EST_GLOBAL_TRACE(s) P4EST_GLOBAL_LOG (SC_LP_TRACE, (s))
#define P4EST_GLOBAL_TRACEF(f,...) \
  P4EST_GLOBAL_LOGF (SC_LP_TRACE, (f), __VA_ARGS__)
#define P4EST_GLOBAL_LDEBUG(s) P4EST_GLOBAL_LOG (SC_LP_DEBUG, (s))
#define P4EST_GLOBAL_LDEBUGF(f,...) \
  P4EST_GLOBAL_LOGF (SC_LP_DEBUG, (f), __VA_ARGS__)
#define P4EST_GLOBAL_VERBOSE(s) P4EST_GLOBAL_LOG (SC_LP_VERBOSE, (s))
#define P4EST_GLOBAL_VERBOSEF(f,...) \
  P4EST_GLOBAL_LOGF (SC_LP_VERBOSE, (f), __VA_ARGS__)
#define P4EST_GLOBAL_INFO(s) P4EST_GLOBAL_LOG (SC_LP_INFO, (s))
#define P4EST_GLOBAL_INFOF(f,...) \
  P4EST_GLOBAL_LOGF (SC_LP_INFO, (f), __VA_ARGS__)
#define P4EST_GLOBAL_STATISTICS(s) P4EST_GLOBAL_LOG (SC_LP_STATISTICS, (s))
#define P4EST_GLOBAL_STATISTICSF(f,...) \
  P4EST_GLOBAL_LOGF (SC_LP_STATISTICS, (f), __VA_ARGS__)
#define P4EST_GLOBAL_PRODUCTION(s) P4EST_GLOBAL_LOG (SC_LP_PRODUCTION, (s))
#define P4EST_GLOBAL_PRODUCTIONF(f,...) \
  P4EST_GLOBAL_LOGF (SC_LP_PRODUCTION, (f), __VA_ARGS__)

/* convenience log macros that are active on every processor */
#define P4EST_TRACE(s) P4EST_NORMAL_LOG (SC_LP_TRACE, (s))
#define P4EST_TRACEF(f,...) \
  P4EST_NORMAL_LOGF (SC_LP_TRACE, (f), __VA_ARGS__)
#define P4EST_LDEBUG(s) P4EST_NORMAL_LOG (SC_LP_DEBUG, (s))
#define P4EST_LDEBUGF(f,...) \
  P4EST_NORMAL_LOGF (SC_LP_DEBUG, (f), __VA_ARGS__)
#define P4EST_VERBOSE(s) P4EST_NORMAL_LOG (SC_LP_VERBOSE, (s))
#define P4EST_VERBOSEF(f,...) \
  P4EST_NORMAL_LOGF (SC_LP_VERBOSE, (f), __VA_ARGS__)
#define P4EST_INFO(s) P4EST_NORMAL_LOG (SC_LP_INFO, (s))
#define P4EST_INFOF(f,...) \
  P4EST_NORMAL_LOGF (SC_LP_INFO, (f), __VA_ARGS__)
#define P4EST_STATISTICS(s) P4EST_NORMAL_LOG (SC_LP_STATISTICS, (s))
#define P4EST_STATISTICSF(f,...) \
  P4EST_NORMAL_LOGF (SC_LP_STATISTICS, (f), __VA_ARGS__)
#define P4EST_PRODUCTION(s) P4EST_NORMAL_LOG (SC_LP_PRODUCTION, (s))
#define P4EST_PRODUCTIONF(f,...) \
  P4EST_NORMAL_LOGF (SC_LP_PRODUCTION, (f), __VA_ARGS__)

typedef void        (*p4est_handler_t) (void *data);

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

#endif /* !P4EST_BASE_H */
