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

#include <p4est_base.h>

int                 p4est_package_id = -1;

void
p4est_init (sc_log_handler_t log_handler, int log_threshold)
{
  int                 w;

  p4est_package_id = sc_package_register (log_handler, log_threshold,
                                          "p4est", "A forest of octrees");

  w = 24;
  P4EST_GLOBAL_ESSENTIALF ("This is %s\n", P4EST_PACKAGE_STRING);
  P4EST_GLOBAL_PRODUCTIONF ("%-*s %s\n", w, "CPP", P4EST_CPP);
  P4EST_GLOBAL_PRODUCTIONF ("%-*s %s\n", w, "CPPFLAGS", P4EST_CPPFLAGS);
  P4EST_GLOBAL_PRODUCTIONF ("%-*s %s\n", w, "CC", P4EST_CC);
#if 0
  P4EST_GLOBAL_PRODUCTIONF ("%-*s %s\n", w, "C_VERSION", P4EST_C_VERSION);
#endif
  P4EST_GLOBAL_PRODUCTIONF ("%-*s %s\n", w, "CFLAGS", P4EST_CFLAGS);
  P4EST_GLOBAL_PRODUCTIONF ("%-*s %s\n", w, "LDFLAGS", P4EST_LDFLAGS);
  P4EST_GLOBAL_PRODUCTIONF ("%-*s %s\n", w, "LIBS", P4EST_LIBS);
}

#ifndef __cplusplus
#undef P4EST_GLOBAL_LOGF
#undef P4EST_LOGF
#undef P4EST_GLOBAL_TRACEF
#undef P4EST_GLOBAL_LDEBUGF
#undef P4EST_GLOBAL_VERBOSEF
#undef P4EST_GLOBAL_INFOF
#undef P4EST_GLOBAL_STATISTICSF
#undef P4EST_GLOBAL_PRODUCTIONF
#undef P4EST_GLOBAL_ESSENTIALF
#undef P4EST_GLOBAL_LERRORF
#undef P4EST_TRACEF
#undef P4EST_LDEBUGF
#undef P4EST_VERBOSEF
#undef P4EST_INFOF
#undef P4EST_STATISTICSF
#undef P4EST_PRODUCTIONF
#undef P4EST_ESSENTIALF
#undef P4EST_LERRORF
#endif

#ifndef SC_SPLINT

void
P4EST_GLOBAL_LOGF (int priority, const char *fmt, ...)
{
  va_list             ap;

  va_start (ap, fmt);
  sc_logv ("<unknown>", 0, p4est_package_id, SC_LC_GLOBAL, priority, fmt, ap);
  va_end (ap);
}

void
P4EST_LOGF (int priority, const char *fmt, ...)
{
  va_list             ap;

  va_start (ap, fmt);
  sc_logv ("<unknown>", 0, p4est_package_id, SC_LC_NORMAL, priority, fmt, ap);
  va_end (ap);
}

#define P4EST_LOG_IMP(n,p)                              \
  void                                                  \
  P4EST_GLOBAL_ ## n ## F (const char *fmt, ...)        \
  {                                                     \
    va_list             ap;                             \
    va_start (ap, fmt);                                 \
    sc_logv ("<unknown>", 0, p4est_package_id,          \
             SC_LC_GLOBAL, SC_LP_ ## p, fmt, ap);       \
    va_end (ap);                                        \
  }                                                     \
  void                                                  \
  P4EST_ ## n ## F (const char *fmt, ...)               \
  {                                                     \
    va_list             ap;                             \
    va_start (ap, fmt);                                 \
    sc_logv ("<unknown>", 0, p4est_package_id,          \
             SC_LC_NORMAL, SC_LP_ ## p, fmt, ap);       \
    va_end (ap);                                        \
  }

/* *INDENT-OFF* */
P4EST_LOG_IMP (TRACE, TRACE)
P4EST_LOG_IMP (LDEBUG, DEBUG)
P4EST_LOG_IMP (VERBOSE, VERBOSE)
P4EST_LOG_IMP (INFO, INFO)
P4EST_LOG_IMP (STATISTICS, STATISTICS)
P4EST_LOG_IMP (PRODUCTION, PRODUCTION)
P4EST_LOG_IMP (ESSENTIAL, ESSENTIAL)
P4EST_LOG_IMP (LERROR, ERROR)
/* *INDENT-ON* */

#endif
