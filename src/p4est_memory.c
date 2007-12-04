/*
  This file is part of p4est.
  p4est is a C library to manage a parallel collection of quadtrees and/or
  octrees.

  Copyright (C) 2007 Carsten Burstedde, Lucas Wilcox.

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

#include <p4est_memory.h>
#include <p4est_base.h>

/* array routines */

static const int    elements_per_chunk = 1024;

p4est_array_t      *
p4est_array_new (int elem_size)
{
  p4est_array_t      *array;

  P4EST_ASSERT (elem_size > 0);

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

  P4EST_ASSERT (elem_size > 0);

  mempool = P4EST_ALLOC_ZERO (p4est_mempool_t, 1);
  P4EST_CHECK_ALLOC (mempool);

  mempool->elem_size = elem_size;
  mempool->elem_count = 0;

  obstack_init (&mempool->obstack);
  obstack_chunk_size (&mempool->obstack) = elements_per_chunk * elem_size;
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

/* list routines */

p4est_list_t       *
p4est_list_new (p4est_mempool_t * allocator)
{
  p4est_list_t       *list;

  list = P4EST_ALLOC_ZERO (p4est_list_t, 1);
  list->elem_count = 0;
  list->first = NULL;
  list->last = NULL;

  if (allocator != NULL) {
    list->allocator = allocator;
    list->allocator_owned = 0;
  }
  else {
    list->allocator = p4est_mempool_new (sizeof (p4est_link_t));
    list->allocator_owned = 1;
  }

  return list;
}

void
p4est_list_destroy (p4est_list_t * list)
{
  p4est_link_t       *link;
  p4est_link_t       *temp;

  link = list->first;
  while (link != NULL) {
    temp = link->next;
    p4est_mempool_free (list->allocator, link);
    link = temp;
    --list->elem_count;
  }
  P4EST_ASSERT (list->elem_count == 0);

  if (list->allocator_owned) {
    p4est_mempool_destroy (list->allocator);
  }
  P4EST_FREE (list);
}

void
p4est_list_prepend (p4est_list_t * list, void *data)
{
  p4est_link_t       *link;

  link = p4est_mempool_alloc (list->allocator);
  link->data = data;
  link->next = list->first;
  list->first = link;
  if (list->last == NULL) {
    list->last = link;
  }

  ++list->elem_count;
}

void
p4est_list_append (p4est_list_t * list, void *data)
{
  p4est_link_t       *link;

  link = p4est_mempool_alloc (list->allocator);
  link->data = data;
  link->next = NULL;
  if (list->last != NULL) {
    list->last->next = link;
  }
  else {
    list->first = link;
  }
  list->last = link;

  ++list->elem_count;
}

void
p4est_list_insert (p4est_list_t * list, p4est_link_t * after, void *data)
{
  p4est_link_t       *link;

  P4EST_ASSERT (after != NULL);

  link = p4est_mempool_alloc (list->allocator);
  link->data = data;
  link->next = after->next;
  after->next = link;
  if (after == list->last) {
    list->last = link;
  }

  ++list->elem_count;
}

void               *
p4est_list_pop (p4est_list_t * list)
{
  p4est_link_t       *link;
  void               *data;

  P4EST_ASSERT (list->first != NULL);

  link = list->first;
  list->first = link->next;
  data = link->data;
  p4est_mempool_free (list->allocator, link);
  if (list->first == NULL) {
    list->last = NULL;
  }

  --list->elem_count;
  return data;
}

/* EOF p4est_memory.c */
