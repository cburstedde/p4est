
#include <p4est_memory.h>
#include <p4est_base.h>

/* array routines */

p4est_array_t      *
p4est_array_new (int elem_size)
{
  p4est_array_t      *array;

  array = P4EST_ALLOC_ZERO (p4est_array_t, 1);
  P4EST_CHECK_ALLOC (array);

  array->elem_size = elem_size;
  array->elem_count = 0;
  array->elem_alloc = 0;
  array->array = NULL;

  return array;
}

void
p4est_array_destroy (p4est_array_t * array)
{
  P4EST_FREE (array->array);
  P4EST_FREE (array);
}

void
p4est_array_resize (p4est_array_t * array, int new_count)
{
  char               *ptr;
  int                 newsize;

  array->elem_count = new_count;

  if (new_count > array->elem_alloc) {
    array->elem_alloc = P4EST_MAX (8 + 2 * array->elem_alloc, new_count);
  }
  else if (new_count < (array->elem_alloc + 1) / 2) {
    array->elem_alloc = new_count;
  }
  else {
    return;
  }
  P4EST_ASSERT (array->elem_alloc >= 0 && array->elem_alloc >= new_count);

  newsize = array->elem_alloc * array->elem_size;
  ptr = P4EST_REALLOC (array->array, char, newsize);
  P4EST_CHECK_REALLOC (ptr, newsize);

  array->array = ptr;
}

void
p4est_array_sort (p4est_array_t * array,
                  int (*compar) (const void *, const void *))
{
  qsort (array->array, array->elem_count, array->elem_size, compar);
}

void               *
p4est_array_index (p4est_array_t * array, int index)
{
  P4EST_ASSERT (index >= 0 && index < array->elem_count);

  return (void *) (array->array + (array->elem_size * index));
}

/* mempool routines */

static void        *(*obstack_chunk_alloc) (size_t) = p4est_malloc;
static void         (*obstack_chunk_free) (void *) = p4est_free;

p4est_mempool_t    *
p4est_mempool_new (int elem_size)
{
  p4est_mempool_t    *mempool;

  mempool = P4EST_ALLOC_ZERO (p4est_mempool_t, 1);
  P4EST_CHECK_ALLOC (mempool);

  mempool->elem_size = elem_size;
  mempool->elem_count = 0;

  obstack_init (&mempool->obstack);
  obstack_chunk_size (&mempool->obstack) = 1000 * elem_size;
  mempool->freed = p4est_array_new (sizeof (void *));

  return mempool;
}

void
p4est_mempool_destroy (p4est_mempool_t * mempool)
{
  p4est_array_destroy (mempool->freed);
  obstack_free (&mempool->obstack, NULL);

  P4EST_FREE (mempool);
}

void
p4est_mempool_reset (p4est_mempool_t * mempool)
{
  p4est_array_resize (mempool->freed, 0);
  obstack_free (&mempool->obstack, NULL);
  mempool->elem_count = 0;
}

void               *
p4est_mempool_alloc (p4est_mempool_t * mempool)
{
  int                 new_count;
  void               *ret;

  ++mempool->elem_count;

  if (mempool->freed->elem_count > 0) {
    new_count = mempool->freed->elem_count - 1;
    ret = *(void **) p4est_array_index (mempool->freed, new_count);
    p4est_array_resize (mempool->freed, new_count);
    return ret;
  }

  return obstack_alloc (&mempool->obstack, mempool->elem_size);
}

void
p4est_mempool_free (p4est_mempool_t * mempool, void *elem)
{
  int                 old_count;

  --mempool->elem_count;

  old_count = mempool->freed->elem_count;
  p4est_array_resize (mempool->freed, old_count + 1);
  *(void **) p4est_array_index (mempool->freed, old_count) = elem;
}

/* EOF p4est_memory.c */
