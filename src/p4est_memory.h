
#ifndef __P4EST_MEMORY_H__
#define __P4EST_MEMORY_H__

/* p4est memory management */

#include <obstack.h>

/*
 * Prototypes for memory allocation functions.
 * p4est_array and p4est_mempool print an error msg to stderr and call abort()
 *    if memory allocation functions fail. This can be prevented by providing
 *    functions that never return (i.e. call longjmp() or abort() internally).
 */
typedef void       *(*p4est_malloc_t) (size_t size);
typedef void       *(*p4est_realloc_t) (void *ptr, size_t size);
typedef void       *(*p4est_free_t) (void *ptr);

/*
 * The p4est_array object provides a large array of equal-size elements.
 * The array can be resized.
 * Elements are accessed by their 0-based index, their address may change.
 * Elements can be sorted.
 */
typedef struct p4est_array
{
  /* interface variables */
  int                 elem_size;        /* size of a single element */
  int                 elem_count;       /* number of valid elements */

  /* implementation variables */
  int                 elem_alloc;       /* number of allocated elements */
  void               *array;    /* linear array to store elements */

  /* allocation functions */
  p4est_malloc_t      my_malloc;
  p4est_realloc_t     my_realloc;
  p4est_free_t        my_free;
}
p4est_array_t;

p4est_array_t      *p4est_array_new (int elem_size, p4est_malloc_t * malloc,
                                     p4est_realloc_t * realloc,
                                     p4est_free_t * free);
void                p4est_array_destroy (p4est_array_t * array);

void                p4est_array_grow1 (p4est_array_t * array);
void                p4est_array_resize (p4est_array_t * array, int new_count);
void                p4est_array_sort (p4est_array_t * array,
                                      int (*compar) (const void *,
                                                     const void *));

void               *p4est_array_index (p4est_array_t * array, int index);

/*
 * The p4est_mempool object provides a large pool of equal-size elements.
 * The pool grows dynamically for element allocation.
 * Elements are referenced by their address which never changes.
 * Elements can be freed (that is, returned to the pool)
 *    and are transparently reused.
 */
typedef struct p4est_mempool
{
  /* interface variables */
  int                 elem_size;        /* size of a single element */
  int                 elem_count;       /* number of valid elements */

  /* implementation variables */
  struct obstack      stack;    /* holds the allocated elements */
  p4est_array_t       freed;    /* buffers the freed elements */

  /* allocation functions */
  p4est_malloc_t      my_malloc;
  p4est_realloc_t     my_realloc;
  p4est_free_t        my_free;
}
p4est_mempool_t;

/* create, destroy or reset to count=0 the memory pool */
p4est_mempool_t    *p4est_mempool_new (int elem_size, p4est_malloc_t * malloc,
                                       p4est_realloc_t * realloc,
                                       p4est_free_t * free);
void                p4est_mempool_destroy (p4est_mempool_t * mempool);
void                p4est_mempool_reset (p4est_mempool_t * mempool);

/* allocate or free a single element */
void               *p4est_mempool_alloc (p4est_mempool_t * mempool);
void               *p4est_mempool_free (p4est_mempool_t * mempool,
                                        void *elem);

#endif /* !__P4EST_MEMORY_H__ */
