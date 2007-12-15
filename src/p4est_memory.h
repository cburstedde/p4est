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

#ifndef __P4EST_MEMORY_H__
#define __P4EST_MEMORY_H__

/* p4est memory management */

#include <p4est_obstack.h>
#include <stdio.h>

/** The p4est_array object provides a large array of equal-size elements.
 * The array can be resized.
 * Elements are accessed by their 0-based index, their address may change.
 * The size (== elem_count) of the array can be changed by array_resize.
 * Elements can be sorted with array_sort.
 * If the array is sorted elements can be binary searched with array_bsearch.
 * A priority queue is implemented with pqueue_add and pqueue_pop.
 * Use sort and search whenever possible, they are faster than the pqueue.
 */
typedef struct p4est_array
{
  /* interface variables */
  int                 elem_size;        /* size of a single element */
  int                 elem_count;       /* number of valid elements */

  /* implementation variables */
  int                 elem_alloc;       /* number of allocated elements */
  char               *array;    /* linear array to store elements */
}
p4est_array_t;

p4est_array_t      *p4est_array_new (int elem_size);
void                p4est_array_destroy (p4est_array_t * array);

/** Sets the element count to new_count.
 * Reallocation takes place only occasionally, so this function is usually fast.
 */
void                p4est_array_resize (p4est_array_t * array, int new_count);

/** Sorts the array in ascending order wrt. the comparison function.
 * \param [in] compar The comparison function to be used.
 */
void                p4est_array_sort (p4est_array_t * array,
                                      int (*compar) (const void *,
                                                     const void *));

/** Performs a binary search on an array. The array must be sorted.
 * \param [in] key     An element to be searched for.
 * \param [in] compar  The comparison function to be used.
 * \return Returns a pointer to the item found, or NULL.
 */
void               *p4est_array_bsearch (p4est_array_t * array,
                                         const void *key,
                                         int (*compar) (const void *,
                                                        const void *));

/** Adds an element to a priority queue.
 * The priority queue is implemented as a heap in ascending order.
 * A heap is a binary tree where the children are not less than their parent.
 * Assumes that elements [0]..[elem_count-2] form a valid heap.
 * Then propagates [elem_count-1] upward by swapping if necessary.
 * \param [in] temp    Pointer to unused allocated memory of elem_size.
 * \param [in] compar  The comparison function to be used.
 * \return Returns the number of swap operations.
 * \note  If the return value is zero for all elements in an array,
 *        the array is sorted linearly and unchanged.
 */
int                 p4est_array_pqueue_add (p4est_array_t * array,
                                            void *temp,
                                            int (*compar) (const void *,
                                                           const void *));

/** Pops the smallest element from a priority queue.
 * This function assumes that the array forms a valid heap in ascending order.
 * \param [out] result  Pointer to unused allocated memory of elem_size.
 * \param [in]  compar  The comparison function to be used.
 * \return Returns the number of swap operations.
 * \note This function resizes the array to elem_count-1.
 */
int                 p4est_array_pqueue_pop (p4est_array_t * array,
                                            void *result,
                                            int (*compar) (const void *,
                                                           const void *));

/** Returns a pointer to an array element.
 * \param [in] index needs to be in [0]..[elem_count-1].
 */
void               *p4est_array_index (p4est_array_t * array, int index);

/** The p4est_mempool object provides a large pool of equal-size elements.
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
  struct obstack      obstack;  /* holds the allocated elements */
  p4est_array_t      *freed;    /* buffers the freed elements */
}
p4est_mempool_t;

/* create, destroy or reset to count=0 the memory pool */
p4est_mempool_t    *p4est_mempool_new (int elem_size);
void                p4est_mempool_destroy (p4est_mempool_t * mempool);

/** Invalidates all previously returned pointers, resets count to 0. */
void                p4est_mempool_reset (p4est_mempool_t * mempool);

/* allocate or free a single element */
void               *p4est_mempool_alloc (p4est_mempool_t * mempool);
void                p4est_mempool_free (p4est_mempool_t * mempool,
                                        void *elem);

/** The p4est_link structure is one link of a linked list.
 */
typedef struct p4est_link
{
  void               *data;
  struct p4est_link  *next;
}
p4est_link_t;

/** The p4est_list object provides a linked list.
 */
typedef struct p4est_list
{
  /* interface variables */
  int                 elem_count;
  p4est_link_t       *first;
  p4est_link_t       *last;

  /* implementation variables */
  int                 allocator_owned;
  p4est_mempool_t    *allocator;        /* must allocate p4est_link_t */
}
p4est_list_t;

/** Allocate a linked list structure.
 * \param [in] allocator Memory allocator for p4est_link_t, can be NULL.
 */
p4est_list_t       *p4est_list_new (p4est_mempool_t * allocator);

/** Destroy a linked list structure in O(N).
 * \note If allocator was provided in p4est_list_new, it will not be destroyed.
 */
void                p4est_list_destroy (p4est_list_t * list);

/** Initializes an already allocated list structure.
 * \param [in,out]  list       List structure to be initialized.
 * \param [in]      allocator  External memory allocator for p4est_link_t.
 */
void                p4est_list_init (p4est_list_t * list,
                                     p4est_mempool_t * allocator);

/** Removes all elements from a list in O(N).
 * \param [in,out]  list       List structure to be resetted.
 * \note Calling p4est_list_init, then any list operations,
 *       then p4est_list_reset is memory neutral.
 */
void                p4est_list_reset (p4est_list_t * list);

/** Unliks all list elements without returning them to the mempool.
 * This runs in O(1) but is dangerous because of potential memory leaks.
 * \param [in,out]  list       List structure to be unlinked.
 */
void                p4est_list_unlink (p4est_list_t * list);

void                p4est_list_prepend (p4est_list_t * list, void *data);
void                p4est_list_append (p4est_list_t * list, void *data);

/** Insert an element after a given position.
 * \param [in] pred The predecessor of the element to be inserted.
 */
void                p4est_list_insert (p4est_list_t * list,
                                       p4est_link_t * pred, void *data);

/** Remove an element after a given position.
 * \param [in] pred  The predecessor of the element to be removed.
                     If \a pred == NULL, the first element is removed.
 * \return Returns the data of the removed element.
 */
void               *p4est_list_remove (p4est_list_t * list,
                                       p4est_link_t * pred);

/** Remove an element from the front of the list.
 * \return Returns the data of the removed first list element.
 */
void               *p4est_list_pop (p4est_list_t * list);

/** Function to compute a hash value of an object.
 * \return Returns a non-negative integer.
 */
typedef int         (*p4est_hash_function_t) (const void *v);

/** Function to check equality of two objects.
 * \return Returns 0 if *v1 is unequal *v2 and true otherwise.
 */
typedef int         (*p4est_equal_function_t) (const void *v1,
                                               const void *v2);

/** The p4est_hash implements a hash table.
 * It uses an array which has linked lists as elements.
 */
typedef struct p4est_hash
{
  /* interface variables */
  int                 elem_count;       /* total number of objects contained */

  /* implementation variables */
  p4est_array_t      *table;
  p4est_hash_function_t hash_fn;
  p4est_equal_function_t equal_fn;
  int                 allocator_owned;
  p4est_mempool_t    *allocator;        /* must allocate p4est_link_t */
}
p4est_hash_t;

/** Create a new hash table.
 * \param [in] table_size  The number of slots in the hash table.
 * \param [in] hash_fn     Function to compute the hash value.
 * \param [in] equal_fn    Function to test two objects for equality.
 * \param [in] allocator   Memory allocator for p4est_link_t, can be NULL.
 */
p4est_hash_t       *p4est_hash_new (int table_size,
                                    p4est_hash_function_t hash_fn,
                                    p4est_equal_function_t equal_fn,
                                    p4est_mempool_t * allocator);
/** Destroy a hash table in O(N).
 * \note If allocator was provided in p4est_hash_new, it will not be destroyed.
 */
void                p4est_hash_destroy (p4est_hash_t * hash);

/** Remove all entries from a hash table in O(N). */
void                p4est_hash_reset (p4est_hash_t * hash);

/** Unlinks all hash elements without returning them to the mempool.
 * This runs faster than p4est_hash_reset, still in O(N),
 *    but is dangerous because of potential memory leaks.
 * \param [in,out]  hash       Hash structure to be unlinked.
 */
void                p4est_hash_unlink (p4est_hash_t * hash);

/** Check if an object is contained in the hash table.
 * \param [in]  v      The object to be looked up.
 * \param [out] found  If found != NULL, *found is set to the object
 *                     if the object is found.
 * \return Returns 1 if object is found, 0 otherwise.
 */
int                 p4est_hash_lookup (p4est_hash_t * hash, void *v,
                                       void **found);

/** Insert an object into a hash table if it is not contained already.
 * \param [in]  v      The object to be inserted.
 * \param [out] found  If found != NULL, *found is set to the object
 *                     that is already contained if that exists.
 * \return Returns 1 if object is added, 0 if it is already contained.
 */
int                 p4est_hash_insert_unique (p4est_hash_t * hash, void *v,
                                              void **found);

/** Remove an object from a hash table.
 * \param [in]  v      The object to be removed.
 * \param [out] found  If found != NULL, *found is set to the object
                       that is removed if that exists.
 * \return Returns 1 if object is found, 0 if is not contained.
 */
int                 p4est_hash_remove (p4est_hash_t * hash, void *v,
                                       void **found);

/** Compute and print statistical information about the occupancy.
 * \param [in] nout  Stream to print to. May be NULL, then nothing happens.
 */
void                p4est_hash_print_statistics (p4est_hash_t * hash,
                                                 FILE * nout);

#endif /* !__P4EST_MEMORY_H__ */
