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

#include <p4est_memory.h>
#include <p4est_base.h>

/* using sqrt to compute hash statistics */
#include <math.h>

/* require zlib header for adler32 checksums */
#ifdef HAVE_ZLIB_H
#include <zlib.h>
#endif

/* array routines */

static const int    elements_per_chunk = 1024;

p4est_array_t      *
p4est_array_new (int elem_size)
{
  p4est_array_t      *array;

  array = P4EST_ALLOC_ZERO (p4est_array_t, 1);
  P4EST_CHECK_ALLOC (array);

  p4est_array_init (array, elem_size);

  return array;
}

void
p4est_array_destroy (p4est_array_t * array)
{
  P4EST_FREE (array->array);
  P4EST_FREE (array);
}

void
p4est_array_init (p4est_array_t * array, int elem_size)
{
  P4EST_ASSERT (elem_size > 0);

  array->elem_size = elem_size;
  array->elem_count = 0;
  array->elem_alloc = 0;
  array->array = NULL;
}

void
p4est_array_reset (p4est_array_t * array)
{
  P4EST_FREE (array->array);
  array->array = NULL;

  array->elem_count = 0;
  array->elem_alloc = 0;
}

void
p4est_array_resize (p4est_array_t * array, int new_count)
{
  char               *ptr;
  int                 newsize, roundup;
#ifdef P4EST_HAVE_DEBUG
  int                 oldoffs, newoffs;
  int                 i, minoffs;
#endif

  P4EST_ASSERT (new_count >= 0);
#ifdef P4EST_HAVE_DEBUG
  oldoffs = array->elem_count * array->elem_size;
#endif
  array->elem_count = new_count;
#ifdef P4EST_HAVE_DEBUG
  newoffs = array->elem_count * array->elem_size;
#endif
  roundup = P4EST_ROUNDUP2_32 (new_count);
  P4EST_ASSERT (roundup >= new_count && roundup <= 2 * new_count);

  if (new_count > array->elem_alloc || roundup < array->elem_alloc) {
    array->elem_alloc = roundup;
  }
  else {
#ifdef P4EST_HAVE_DEBUG
    if (newoffs < oldoffs) {
      memset (array->array + newoffs, -1, oldoffs - newoffs);
    }
    for (i = oldoffs; i < newoffs; ++i) {
      P4EST_ASSERT (array->array[i] == (char) -1);
    }
#endif
    return;
  }
  P4EST_ASSERT (array->elem_alloc >= 0 &&
                array->elem_alloc >= new_count);

  newsize = array->elem_alloc * array->elem_size;
  ptr = P4EST_REALLOC (array->array, char, newsize);
  P4EST_CHECK_REALLOC (ptr, newsize);

  array->array = ptr;
#ifdef P4EST_HAVE_DEBUG
  minoffs = P4EST_MIN (oldoffs, newoffs);
  P4EST_ASSERT (minoffs <= newsize);
  memset (array->array + minoffs, -1, newsize - minoffs);
#endif
}

void
p4est_array_sort (p4est_array_t * array,
                  int (*compar) (const void *, const void *))
{
  qsort (array->array, array->elem_count, array->elem_size, compar);
}

void               *
p4est_array_bsearch (p4est_array_t * array, const void *key,
                     int (*compar) (const void *, const void *))
{
  return
    bsearch (key, array->array, array->elem_count, array->elem_size, compar);
}

unsigned
p4est_array_checksum (p4est_array_t * array, int first_elem)
{
  int                 first_byte;
  uInt                bytes;
  uLong               crc;

  P4EST_ASSERT (0 <= first_elem && first_elem <= array->elem_count);

  crc = adler32 (0L, Z_NULL, 0);
  if (array->elem_count == 0) {
    return (unsigned) crc;
  }

  first_byte = first_elem * array->elem_size;
  bytes = (array->elem_count - first_elem) * array->elem_size;
  crc = adler32 (crc, (const Bytef *) (array->array + first_byte), bytes);

  return (unsigned) crc;
}

int
p4est_array_pqueue_add (p4est_array_t * array, void *temp,
                        int (*compar) (const void *, const void *))
{
  int                 parent, child;
  int                 comp, swaps;
  const int           size = array->elem_size;
  void               *p, *c;

  swaps = 0;
  child = array->elem_count - 1;
  c = array->array + (size * child);
  while (child > 0) {
    parent = (child - 1) / 2;
    p = array->array + (size * parent);

    /* compare child to parent */
    comp = compar (p, c);
    if (comp <= 0) {
      break;
    }

    /* swap child and parent */
    memcpy (temp, c, size);
    memcpy (c, p, size);
    memcpy (p, temp, size);
    ++swaps;

    /* walk up the tree */
    child = parent;
    c = p;
  }

  return swaps;
}

int
p4est_array_pqueue_pop (p4est_array_t * array, void *result,
                        int (*compar) (const void *, const void *))
{
  int                 new_count;
  int                 parent, child, child1;
  int                 comp, swaps;
  const int           size = array->elem_size;
  void               *p, *c, *c1;
  void               *temp;

  swaps = 0;
  new_count = array->elem_count - 1;
  P4EST_ASSERT (new_count >= 0);

  /* extract root */
  parent = 0;
  p = array->array + (size * parent);
  memcpy (result, p, size);

  /* copy the last element to the top and reuse it as temp storage */
  temp = array->array + (size * new_count);
  if (new_count > 0) {
    memcpy (p, temp, size);
  }

  /* sift down the tree */
  while ((child = 2 * parent + 1) < new_count) {
    c = array->array + (size * child);

    /* check if child has a sibling and use that one if it is smaller */
    if ((child1 = 2 * parent + 2) < new_count) {
      c1 = array->array + (size * child1);
      comp = compar (c, c1);
      if (comp > 0) {
        child = child1;
        c = c1;
      }
    }

    /* sift down the parent if it is larger */
    comp = compar (p, c);
    if (comp <= 0) {
      break;
    }

    /* swap child and parent */
    memcpy (temp, c, size);
    memcpy (c, p, size);
    memcpy (p, temp, size);
    ++swaps;

    /* walk down the tree */
    parent = child;
    p = c;
  }

  /* we can resize down here only since we need the temp element above */
  p4est_array_resize (array, new_count);

  return swaps;
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
  p4est_array_reset (mempool->freed);
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
  }
  else {
    ret = obstack_alloc (&mempool->obstack, mempool->elem_size);
  }

#ifdef P4EST_HAVE_DEBUG
  memset (ret, -1, mempool->elem_size);
#endif

  return ret;
}

void
p4est_mempool_free (p4est_mempool_t * mempool, void *elem)
{
  int                 old_count;

#ifdef P4EST_HAVE_DEBUG
  memset (elem, -1, mempool->elem_size);
#endif

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
  P4EST_CHECK_ALLOC (list);

  list->elem_count = 0;
  list->first = NULL;
  list->last = NULL;

  if (allocator != NULL) {
    P4EST_ASSERT (allocator->elem_size == sizeof (p4est_link_t));
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
  p4est_list_reset (list);

  if (list->allocator_owned) {
    p4est_mempool_destroy (list->allocator);
  }
  P4EST_FREE (list);
}

void
p4est_list_init (p4est_list_t * list, p4est_mempool_t * allocator)
{
  list->elem_count = 0;
  list->first = NULL;
  list->last = NULL;

  P4EST_ASSERT (allocator != NULL);
  P4EST_ASSERT (allocator->elem_size == sizeof (p4est_link_t));

  list->allocator = allocator;
  list->allocator_owned = 0;
}

void
p4est_list_reset (p4est_list_t * list)
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

  list->first = list->last = NULL;
}

void
p4est_list_unlink (p4est_list_t * list)
{
  list->first = list->last = NULL;
  list->elem_count = 0;
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
p4est_list_insert (p4est_list_t * list, p4est_link_t * pred, void *data)
{
  p4est_link_t       *link;

  P4EST_ASSERT (pred != NULL);

  link = p4est_mempool_alloc (list->allocator);
  link->data = data;
  link->next = pred->next;
  pred->next = link;
  if (pred == list->last) {
    list->last = link;
  }

  ++list->elem_count;
}

void               *
p4est_list_remove (p4est_list_t * list, p4est_link_t * pred)
{
  p4est_link_t       *link;
  void               *data;

  if (pred == NULL) {
    return p4est_list_pop (list);
  }

  P4EST_ASSERT (pred->next != NULL);

  link = pred->next;
  pred->next = link->next;
  data = link->data;
  if (list->last == link) {
    list->last = pred;
  }
  p4est_mempool_free (list->allocator, link);

  --list->elem_count;
  return data;
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

/* hash table routines */

p4est_hash_t       *
p4est_hash_new (int table_size, p4est_hash_function_t hash_fn,
                p4est_equal_function_t equal_fn, p4est_mempool_t * allocator)
{
  int                 i;
  p4est_hash_t       *hash;
  p4est_list_t       *list;

  hash = P4EST_ALLOC_ZERO (p4est_hash_t, 1);
  P4EST_CHECK_ALLOC (hash);

  if (allocator != NULL) {
    P4EST_ASSERT (allocator->elem_size == sizeof (p4est_link_t));
    hash->allocator = allocator;
    hash->allocator_owned = 0;
  }
  else {
    hash->allocator = p4est_mempool_new (sizeof (p4est_link_t));
    hash->allocator_owned = 1;
  }

  hash->elem_count = 0;
  hash->table = p4est_array_new (sizeof (p4est_list_t));
  p4est_array_resize (hash->table, table_size);
  for (i = 0; i < table_size; ++i) {
    list = p4est_array_index (hash->table, i);
    p4est_list_init (list, hash->allocator);
  }
  hash->hash_fn = hash_fn;
  hash->equal_fn = equal_fn;

  return hash;
}

void
p4est_hash_destroy (p4est_hash_t * hash)
{
  p4est_hash_reset (hash);

  p4est_array_destroy (hash->table);
  if (hash->allocator_owned) {
    p4est_mempool_destroy (hash->allocator);
  }

  P4EST_FREE (hash);
}

void
p4est_hash_reset (p4est_hash_t * hash)
{
  int                 i, size, count;
  p4est_list_t       *list;

  if (hash->elem_count == 0) {
    return;
  }

  count = 0;
  size = hash->table->elem_count;

  for (i = 0; i < size; ++i) {
    list = p4est_array_index (hash->table, i);
    count += list->elem_count;
    p4est_list_reset (list);
  }
  P4EST_ASSERT (count == hash->elem_count);

  hash->elem_count = 0;
}

void
p4est_hash_unlink (p4est_hash_t * hash)
{
  int                 i, size, count;
  p4est_list_t       *list;

  count = 0;
  size = hash->table->elem_count;
  for (i = 0; i < size; ++i) {
    list = p4est_array_index (hash->table, i);
    count += list->elem_count;
    p4est_list_unlink (list);
  }
  P4EST_ASSERT (count == hash->elem_count);

  hash->elem_count = 0;
}

void
p4est_hash_unlink_destroy (p4est_hash_t * hash)
{
  p4est_array_destroy (hash->table);
  if (hash->allocator_owned) {
    p4est_mempool_destroy (hash->allocator);
  }

  P4EST_FREE (hash);
}

int
p4est_hash_lookup (p4est_hash_t * hash, void *v, void **found)
{
  int                 hval;
  p4est_list_t       *list;
  p4est_link_t       *link;

  hval = hash->hash_fn (v) % hash->table->elem_count;
  list = p4est_array_index (hash->table, hval);

  for (link = list->first; link != NULL; link = link->next) {
    /* check if an equal object is contained in the hash table */
    if (hash->equal_fn (link->data, v)) {
      if (found != NULL) {
        *found = link->data;
      }
      return 1;
    }
  }
  return 0;
}

int
p4est_hash_insert_unique (p4est_hash_t * hash, void *v, void **found)
{
  int                 hval;
  p4est_list_t       *list;
  p4est_link_t       *link;

  hval = hash->hash_fn (v) % hash->table->elem_count;
  list = p4est_array_index (hash->table, hval);

  for (link = list->first; link != NULL; link = link->next) {
    /* check if an equal object is already contained in the hash table */
    if (hash->equal_fn (link->data, v)) {
      if (found != NULL) {
        *found = link->data;
      }
      return 0;
    }
  }
  /* append new object to the list */
  p4est_list_append (list, v);
  ++hash->elem_count;
  return 1;
}

int
p4est_hash_remove (p4est_hash_t * hash, void *v, void **found)
{
  int                 hval;
  p4est_list_t       *list;
  p4est_link_t       *link, *prev;

  hval = hash->hash_fn (v) % hash->table->elem_count;
  list = p4est_array_index (hash->table, hval);

  prev = NULL;
  for (link = list->first; link != NULL; link = link->next) {
    /* check if an equal object is contained in the hash table */
    if (hash->equal_fn (link->data, v)) {
      if (found != NULL) {
        *found = link->data;
      }
      p4est_list_remove (list, prev);
      --hash->elem_count;
      return 1;
    }
    prev = link;
  }
  return 0;
}

void
p4est_hash_print_statistics (p4est_hash_t * hash, FILE * nout)
{
  int                 i;
  int64_t             a, sum, squaresum;
  double              divide, avg, sqr, std;
  p4est_list_t       *list;

  if (nout == NULL) {
    return;
  }

  sum = 0;
  squaresum = 0;
  for (i = 0; i < hash->table->elem_count; ++i) {
    list = p4est_array_index (hash->table, i);
    a = list->elem_count;
    sum += a;
    squaresum += a * a;
  }
  P4EST_ASSERT (sum == hash->elem_count);

  divide = hash->table->elem_count;
  avg = sum / divide;
  sqr = squaresum / divide - avg * avg;
  std = sqrt (sqr);
  fprintf (nout, "Hash size %d avg %.3g std %.3g\n",
           hash->table->elem_count, avg, std);
}

/* EOF p4est_memory.c */
