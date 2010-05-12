/*
  This file is part of p4est.
  p4est is a C library to manage a collection (a forest) of multiple
  connected adaptive quadtrees or octrees in parallel.

  Copyright (C) 2010 The University of Texas System
  Written by Carsten Burstedde, Lucas C. Wilcox, and Tobin Isaac

  p4est is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  p4est is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with p4est; if not, write to the Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

#ifndef P4EST_BASE_H
#define P4EST_BASE_H

/* include p4est config header */

#include <p4est_config.h>

/* indirectly also include sc.h and sc_config.h */

#include <sc_containers.h>

#if \
  (defined (P4EST_MPI) && !defined (SC_MPI)) || \
  (!defined (P4EST_MPI) && defined (SC_MPI))
#error "MPI configured differently in p4est and libsc"
#endif
#if \
  (defined (P4EST_MPIIO) && !defined (SC_MPIIO)) || \
  (!defined (P4EST_MPIIO) && defined (SC_MPIIO))
#error "MPI I/O configured differently in p4est and libsc"
#endif

SC_EXTERN_C_BEGIN;

/** Typedef for quadrant coordinates. */
typedef int32_t     p4est_qcoord_t;
#define P4EST_MPI_QCOORD MPI_INT
#define P4EST_VTK_QCOORD "Int32"
#define P4EST_QCOORD_MIN INT32_MIN
#define P4EST_QCOORD_MAX INT32_MAX

/** Typedef for counting topological entities (trees, vertices). */
typedef int32_t     p4est_topidx_t;
#define P4EST_MPI_TOPIDX MPI_INT
#define P4EST_VTK_TOPIDX "Int32"
#define P4EST_TOPIDX_MAX INT32_MAX

/** Typedef for processor-local indexing of quadrants and nodes. */
typedef int32_t     p4est_locidx_t;
#define p4est_locidx_compare sc_int32_compare
#define P4EST_MPI_LOCIDX MPI_INT
#define P4EST_VTK_LOCIDX "Int32"
#define P4EST_LOCIDX_MAX INT32_MAX

/** Typedef for globally unique indexing of quadrants. */
typedef int64_t     p4est_gloidx_t;
#define P4EST_MPI_GLOIDX MPI_LONG_LONG_INT
#define P4EST_VTK_GLOIDX "Int64"
#define P4EST_GLOIDX_MAX INT64_MAX

/* some error checking possibly specific to p4est */
#ifdef P4EST_DEBUG
#define P4EST_ASSERT(c) SC_CHECK_ABORT ((c), "Assertion '" #c "'")
#else
#define P4EST_ASSERT(c) SC_NOOP ()
#endif

/* macros for memory allocation, will abort if out of memory */
#define P4EST_ALLOC(t,n)          (t *) sc_malloc (p4est_package_id,    \
                                                   (n) * sizeof(t))
#define P4EST_ALLOC_ZERO(t,n)     (t *) sc_calloc (p4est_package_id,    \
                                                   (size_t) (n), sizeof(t))
#define P4EST_REALLOC(p,t,n)      (t *) sc_realloc (p4est_package_id,   \
                                                    (p), (n) * sizeof(t))
#define P4EST_STRDUP(s)                 sc_strdup (p4est_package_id, (s))
#define P4EST_FREE(p)                   sc_free (p4est_package_id, (p))

/* log helper macros */
#define P4EST_GLOBAL_LOG(p,s)                           \
  SC_GEN_LOG (p4est_package_id, SC_LC_GLOBAL, (p), (s))
#define P4EST_LOG(p,s)                                  \
  SC_GEN_LOG (p4est_package_id, SC_LC_NORMAL, (p), (s))
void                P4EST_GLOBAL_LOGF (int priority, const char *fmt, ...)
  __attribute__ ((format (printf, 2, 3)));
void                P4EST_LOGF (int priority, const char *fmt, ...)
  __attribute__ ((format (printf, 2, 3)));
#ifndef __cplusplus
#define P4EST_GLOBAL_LOGF(p,f,...)                                      \
  SC_GEN_LOGF (p4est_package_id, SC_LC_GLOBAL, (p), (f), __VA_ARGS__)
#define P4EST_LOGF(p,f,...)                                             \
  SC_GEN_LOGF (p4est_package_id, SC_LC_NORMAL, (p), (f), __VA_ARGS__)
#endif

/* convenience global log macros will only print if identifier <= 0 */
#define P4EST_GLOBAL_TRACE(s) P4EST_GLOBAL_LOG (SC_LP_TRACE, (s))
#define P4EST_GLOBAL_LDEBUG(s) P4EST_GLOBAL_LOG (SC_LP_DEBUG, (s))
#define P4EST_GLOBAL_VERBOSE(s) P4EST_GLOBAL_LOG (SC_LP_VERBOSE, (s))
#define P4EST_GLOBAL_INFO(s) P4EST_GLOBAL_LOG (SC_LP_INFO, (s))
#define P4EST_GLOBAL_STATISTICS(s) P4EST_GLOBAL_LOG (SC_LP_STATISTICS, (s))
#define P4EST_GLOBAL_PRODUCTION(s) P4EST_GLOBAL_LOG (SC_LP_PRODUCTION, (s))
#define P4EST_GLOBAL_ESSENTIAL(s) P4EST_GLOBAL_LOG (SC_LP_ESSENTIAL, (s))
#define P4EST_GLOBAL_LERROR(s) P4EST_GLOBAL_LOG (SC_LP_ERROR, (s))
void                P4EST_GLOBAL_TRACEF (const char *fmt, ...)
  __attribute__ ((format (printf, 1, 2)));
void                P4EST_GLOBAL_LDEBUGF (const char *fmt, ...)
  __attribute__ ((format (printf, 1, 2)));
void                P4EST_GLOBAL_VERBOSEF (const char *fmt, ...)
  __attribute__ ((format (printf, 1, 2)));
void                P4EST_GLOBAL_INFOF (const char *fmt, ...)
  __attribute__ ((format (printf, 1, 2)));
void                P4EST_GLOBAL_STATISTICSF (const char *fmt, ...)
  __attribute__ ((format (printf, 1, 2)));
void                P4EST_GLOBAL_PRODUCTIONF (const char *fmt, ...)
  __attribute__ ((format (printf, 1, 2)));
void                P4EST_GLOBAL_ESSENTIALF (const char *fmt, ...)
  __attribute__ ((format (printf, 1, 2)));
void                P4EST_GLOBAL_LERRORF (const char *fmt, ...)
  __attribute__ ((format (printf, 1, 2)));
#ifndef __cplusplus
#define P4EST_GLOBAL_TRACEF(f,...)                      \
  P4EST_GLOBAL_LOGF (SC_LP_TRACE, (f), __VA_ARGS__)
#define P4EST_GLOBAL_LDEBUGF(f,...)                     \
  P4EST_GLOBAL_LOGF (SC_LP_DEBUG, (f), __VA_ARGS__)
#define P4EST_GLOBAL_VERBOSEF(f,...)                    \
  P4EST_GLOBAL_LOGF (SC_LP_VERBOSE, (f), __VA_ARGS__)
#define P4EST_GLOBAL_INFOF(f,...)                       \
  P4EST_GLOBAL_LOGF (SC_LP_INFO, (f), __VA_ARGS__)
#define P4EST_GLOBAL_STATISTICSF(f,...)                         \
  P4EST_GLOBAL_LOGF (SC_LP_STATISTICS, (f), __VA_ARGS__)
#define P4EST_GLOBAL_PRODUCTIONF(f,...)                         \
  P4EST_GLOBAL_LOGF (SC_LP_PRODUCTION, (f), __VA_ARGS__)
#define P4EST_GLOBAL_ESSENTIALF(f,...)                          \
  P4EST_GLOBAL_LOGF (SC_LP_ESSENTIAL, (f), __VA_ARGS__)
#define P4EST_GLOBAL_LERRORF(f,...)                     \
  P4EST_GLOBAL_LOGF (SC_LP_ERROR, (f), __VA_ARGS__)
#endif
#define P4EST_GLOBAL_NOTICE     P4EST_GLOBAL_STATISTICS
#define P4EST_GLOBAL_NOTICEF    P4EST_GLOBAL_STATISTICSF

/* convenience log macros that are active on every processor */
#define P4EST_TRACE(s) P4EST_LOG (SC_LP_TRACE, (s))
#define P4EST_LDEBUG(s) P4EST_LOG (SC_LP_DEBUG, (s))
#define P4EST_VERBOSE(s) P4EST_LOG (SC_LP_VERBOSE, (s))
#define P4EST_INFO(s) P4EST_LOG (SC_LP_INFO, (s))
#define P4EST_STATISTICS(s) P4EST_LOG (SC_LP_STATISTICS, (s))
#define P4EST_PRODUCTION(s) P4EST_LOG (SC_LP_PRODUCTION, (s))
#define P4EST_ESSENTIAL(s) P4EST_LOG (SC_LP_ESSENTIAL, (s))
#define P4EST_LERROR(s) P4EST_LOG (SC_LP_ERROR, (s))
void                P4EST_TRACEF (const char *fmt, ...)
  __attribute__ ((format (printf, 1, 2)));
void                P4EST_LDEBUGF (const char *fmt, ...)
  __attribute__ ((format (printf, 1, 2)));
void                P4EST_VERBOSEF (const char *fmt, ...)
  __attribute__ ((format (printf, 1, 2)));
void                P4EST_INFOF (const char *fmt, ...)
  __attribute__ ((format (printf, 1, 2)));
void                P4EST_STATISTICSF (const char *fmt, ...)
  __attribute__ ((format (printf, 1, 2)));
void                P4EST_PRODUCTIONF (const char *fmt, ...)
  __attribute__ ((format (printf, 1, 2)));
void                P4EST_ESSENTIALF (const char *fmt, ...)
  __attribute__ ((format (printf, 1, 2)));
void                P4EST_LERRORF (const char *fmt, ...)
  __attribute__ ((format (printf, 1, 2)));
#ifndef __cplusplus
#define P4EST_TRACEF(f,...)                     \
  P4EST_LOGF (SC_LP_TRACE, (f), __VA_ARGS__)
#define P4EST_LDEBUGF(f,...)                    \
  P4EST_LOGF (SC_LP_DEBUG, (f), __VA_ARGS__)
#define P4EST_VERBOSEF(f,...)                   \
  P4EST_LOGF (SC_LP_VERBOSE, (f), __VA_ARGS__)
#define P4EST_INFOF(f,...)                      \
  P4EST_LOGF (SC_LP_INFO, (f), __VA_ARGS__)
#define P4EST_STATISTICSF(f,...)                        \
  P4EST_LOGF (SC_LP_STATISTICS, (f), __VA_ARGS__)
#define P4EST_PRODUCTIONF(f,...)                        \
  P4EST_LOGF (SC_LP_PRODUCTION, (f), __VA_ARGS__)
#define P4EST_ESSENTIALF(f,...)                         \
  P4EST_LOGF (SC_LP_ESSENTIAL, (f), __VA_ARGS__)
#define P4EST_LERRORF(f,...)                    \
  P4EST_LOGF (SC_LP_ERROR, (f), __VA_ARGS__)
#endif
#define P4EST_NOTICE            P4EST_STATISTICS
#define P4EST_NOTICEF           P4EST_STATISTICSF

/* extern declarations */
extern int          p4est_package_id;

/** Registers p4est with the SC Library and sets the logging behavior.
 * This function is optional.
 * If this function is not called or called with log_handler == NULL,
 * the default SC log handler will be used.
 * If this function is not called or called with log_threshold == SC_LP_DEFAULT,
 * the default SC log threshold will be used.
 * The default SC log settings can be changed with sc_set_log_defaults ().
 */
void                p4est_init (sc_log_handler_t log_handler,
                                int log_threshold);

SC_EXTERN_C_END;

#endif /* !P4EST_BASE_H */
