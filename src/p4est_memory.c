
#include <p4est_memory.h>
#include <p4est_base.h>

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
p4est_array_grow1 (p4est_array_t * array)
{
  p4est_array_resize (array, array->elem_count + 1);
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

/* EOF p4est_memory.c */
