
#ifndef __P4EST_BASE_H__
#define __P4EST_BASE_H__

#ifdef HAVE_CONFIG_H
#include <p4est_config.h>
#endif

#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif

#ifdef HAVE_STDIO_H
#include <stdio.h>
#endif

#define P4EST_ALLOC(t,n) (t *) malloc ((n) * sizeof(t))
#define P4EST_ALLOC_ZERO(t,n) (t *) calloc ((n), sizeof(t))
#define P4EST_REALLOC(p,t,n) (t *) realloc ((p), (n) * sizeof(t))
#define P4EST_FREE(p) free (p)

#endif /* !__P4EST_BASE_H__ */
