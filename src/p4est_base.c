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

#include <p4est_base.h>

int                 p4est_package_id = -1;

ssize_t
p4est_int64_lower_bound (int64_t target, const int64_t * array,
                         size_t size, size_t guess)
{
  size_t              k_low, k_high;
  int64_t             cur;

  if (size == 0)
    return -1;

  k_low = 0;
  k_high = size - 1;
  for (;;) {
    P4EST_ASSERT (k_low <= k_high);
    P4EST_ASSERT (k_low < size && k_high < size);
    P4EST_ASSERT (k_low <= guess && guess <= k_high);

    /* compare two quadrants */
    cur = array[guess];

    /* check if guess is higher or equal target and there's room below it */
    if (target <= cur && (guess > 0 && target <= array[guess - 1])) {
      k_high = guess - 1;
      guess = (k_low + k_high + 1) / 2;
      continue;
    }

    /* check if guess is lower than target */
    if (target > cur) {
      k_low = guess + 1;
      if (k_low > k_high)
        return -1;

      guess = (k_low + k_high) / 2;
      continue;
    }

    /* otherwise guess is the correct position */
    break;
  }

  return (ssize_t) guess;
}

void
p4est_init (sc_log_handler_t log_handler, int log_threshold)
{
  int                 w;
  const char        **on, **ov;
  const char         *overrides[] = { P4EST_OVERRIDES NULL, NULL };

  p4est_package_id = sc_package_register (log_handler, log_threshold,
                                          "p4est", "A forest of octrees");

  w = 24;
  P4EST_GLOBAL_ESSENTIALF ("This is %s.\n", P4EST_PACKAGE_STRING);
  P4EST_GLOBAL_PRODUCTIONF ("%-*s %s\n", w, "CC", P4EST_CC);
  P4EST_GLOBAL_PRODUCTIONF ("%-*s %s\n", w, "C_VERSION", P4EST_C_VERSION);
  P4EST_GLOBAL_PRODUCTIONF ("%-*s %s\n", w, "CFLAGS", P4EST_CFLAGS);
  P4EST_GLOBAL_PRODUCTIONF ("%-*s %s\n", w, "CPP", P4EST_CPP);
  P4EST_GLOBAL_PRODUCTIONF ("%-*s %s\n", w, "CPPFLAGS", P4EST_CPPFLAGS);
  P4EST_GLOBAL_PRODUCTIONF ("%-*s %s\n", w, "LDFLAGS", P4EST_LDFLAGS);
  P4EST_GLOBAL_PRODUCTIONF ("%-*s %s\n", w, "LIBS", P4EST_LIBS);

  w = 32;
  for (on = overrides; *on != NULL; on = ov + 1) {
    ov = on + 1;
    SC_CHECK_ABORT (*ov != NULL, "P4EST_OVERRIDES should contain pairs");
    P4EST_GLOBAL_PRODUCTIONF ("%-*s %s\n", w, *on, *ov);
  }
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
