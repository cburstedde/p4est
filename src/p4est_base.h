
#ifndef __P4EST_BASE_H__
#define __P4EST_BASE_H__

#ifdef HAVE_CONFIG_H
#include <p4est_config.h>
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
#define P4EST_ASSERT(c) P4EST_CHECK_ABORT ((c), "Assertion '" #c "'")

#define P4EST_ALLOC(t,n) (t *) p4est_malloc ((n) * sizeof(t))
#define P4EST_ALLOC_ZERO(t,n) (t *) p4est_calloc ((n), sizeof(t))
#define P4EST_REALLOC(p,t,n) (t *) p4est_realloc ((p), (n) * sizeof(t))

/* it is allowed to call P4EST_FREE (NULL) which does nothing. */
#define P4EST_FREE(p) p4est_free (p)

#define P4EST_MIN(a,b) (((a) < (b)) ? (a) : (b))
#define P4EST_MAX(a,b) (((a) > (b)) ? (a) : (b))

void               *p4est_malloc (size_t size);
void               *p4est_calloc (size_t nmemb, size_t size);
void               *p4est_realloc (void * ptr, size_t size);
void                p4est_free (void * ptr);
void                p4est_memory_check (void);

void                p4est_abort (void);

#endif /* !__P4EST_BASE_H__ */
