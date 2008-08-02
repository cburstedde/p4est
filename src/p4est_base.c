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

#if 0
struct p4est_log_appender
{
  struct LogAppender  appender;
  int                 identifier;
  FILE               *stream;
  FILE               *backup;
};

static struct p4est_log_appender p4est_log_appender_global;
static struct p4est_log_appender p4est_log_appender_rank;

struct LogCategory  p4est_log_category_global = {
  &_LOGV (LOG_ROOT_CAT), NULL, NULL,
  "P4EST_LOG_CATEGORY_GLOBAL", P4EST_LP_UNINITIALIZED, 1,
  NULL, 1
};

struct LogCategory  p4est_log_category_rank = {
  &_LOGV (LOG_ROOT_CAT), NULL, NULL,
  "P4EST_LOG_CATEGORY_RANK", P4EST_LP_UNINITIALIZED, 1,
  NULL, 1
};

static void
p4est_log_append_null (struct LogAppender *this0, struct LogEvent *ev)
{
  P4EST_ASSERT (ev->priority >= 0 && ev->priority <= P4EST_LP_SILENT);
}

static void
p4est_log_append (struct LogAppender *this0, struct LogEvent *ev)
{
  struct p4est_log_appender *this = (struct p4est_log_appender *) this0;
  int                 identifier = this->identifier;
  FILE               *stream = this->stream;
  char                prefix[BUFSIZ];
  char                basenm[BUFSIZ];
  char               *basept;
  va_list             vacopy;

  P4EST_ASSERT (ev->priority >= 0 && ev->priority <= P4EST_LP_SILENT);

  if (ev->priority == P4EST_LP_SILENT) {
    return;
  }

  if (identifier >= 0) {
    snprintf (prefix, BUFSIZ, "[%d] ", identifier);
  }
  else {
    prefix[0] = '\0';
  }

  if (this->backup != NULL) {
    va_copy (vacopy, ev->ap);
    snprintf (basenm, BUFSIZ, "%s", ev->fileName);
    basept = basename (basenm);
    fprintf (this->backup, "%s%s:%d: ", prefix, basept, ev->lineNum);
    vfprintf (this->backup, ev->fmt, vacopy);
    fflush (this->backup);
    va_end (vacopy);
  }

  if (ev->priority <= P4EST_LP_TRACE) {
    fprintf (stream, "%s%s:%d: ", prefix, ev->fileName, ev->lineNum);
  }
  else {
    fputs (prefix, stream);
  }
  vfprintf (stream, ev->fmt, ev->ap);
}

static void
p4est_set_linebuffered (FILE * stream)
{
  setvbuf (stream, NULL, _IOLBF, 0);
}

void
p4est_init_logging (FILE * stream, int identifier)
{
#ifdef P4EST_DEBUG
  char                filename[BUFSIZ];
  char               *job_id;
  char               *job_name;

  job_id = getenv ("JOB_ID");
  job_name = getenv ("JOB_NAME");
#endif

  if (stream == stdout) {
    p4est_set_linebuffered (stream);
  }
  p4est_base_identifier = identifier;

  p4est_log_appender_global.stream = stream;
  p4est_log_appender_global.backup = NULL;
  if (identifier > 0) {
    p4est_log_appender_global.appender.doAppend = p4est_log_append_null;
  }
  else {
#ifdef P4EST_DEBUG
    if (job_id != NULL && job_name != NULL) {
      snprintf (filename, BUFSIZ, "%s.%s_global", job_name, job_id);
    }
    else {
      snprintf (filename, BUFSIZ, "p4est.log_global");
    }
    p4est_log_appender_global.backup = fopen (filename, "wb");
#endif
    p4est_log_appender_global.appender.doAppend = p4est_log_append;
  }
  p4est_log_appender_global.identifier = -1;
  log_setAppender (&p4est_log_category_global,
                   (struct LogAppender *) &p4est_log_appender_global);

  p4est_log_appender_rank.stream = stream;
  p4est_log_appender_rank.backup = NULL;
#ifdef P4EST_DEBUG
  if (job_id != NULL && job_name != NULL) {
    snprintf (filename, BUFSIZ, "%s.%s_%d",
              job_name, job_id, SC_MAX (identifier, 0));
  }
  else {
    snprintf (filename, BUFSIZ, "p4est.log_%d", SC_MAX (identifier, 0));
  }
  p4est_log_appender_rank.backup = fopen (filename, "wb");
#endif
  p4est_log_appender_rank.appender.doAppend = p4est_log_append;
  p4est_log_appender_rank.identifier = identifier;
  log_setAppender (&p4est_log_category_rank,
                   (struct LogAppender *) &p4est_log_appender_rank);

#ifdef P4EST_DEBUG
  log_setThreshold (&p4est_log_category_global, P4EST_LP_DEBUG);
  log_setThreshold (&p4est_log_category_rank, P4EST_LP_DEBUG);
#else
  log_setThreshold (&p4est_log_category_global, P4EST_LP_INFO);
  log_setThreshold (&p4est_log_category_rank, P4EST_LP_INFO);
#endif /* !P4EST_DEBUG */
}
#endif /* 0 */

void
p4est_init (sc_log_handler_t log_handler, int log_threshold)
{
  p4est_package_id = sc_package_register (log_handler, log_threshold,
                                          "p4est", "A forest of octrees");
#ifdef P4EST_DEBUG
  P4EST_GLOBAL_PRODUCTION ("P4EST_DEBUG is ON\n");
#endif
}

/* EOF p4est_base.c */
