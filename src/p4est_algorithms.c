/*
  This file is part of p4est.
  p4est is a C library to manage a collection (a forest) of multiple
  connected adaptive quadtrees or octrees in parallel.

  Copyright (C) 2010 The University of Texas System
  Additional copyright (C) 2011 individual authors
  Written by Carsten Burstedde, Lucas C. Wilcox, and Tobin Isaac

  p4est is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  p4est is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with p4est; if not, write to the Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

#ifdef P4_TO_P8
#include <p8est_algorithms.h>
#include <p8est_bits.h>
#include <p8est_communication.h>
#include <p8est_search.h>
#include <p8est_balance.h>
#else
#include <p4est_algorithms.h>
#include <p4est_bits.h>
#include <p4est_communication.h>
#include <p4est_search.h>
#include <p4est_balance.h>
#endif /* !P4_TO_P8 */

/* htonl is in either of these three */
#ifdef P4EST_HAVE_ARPA_NET_H
#include <arpa/inet.h>
#endif
#ifdef P4EST_HAVE_NETINET_IN_H
#include <netinet/in.h>
#endif
#if defined P4EST_HAVE_WINSOCK2_H || defined _WIN32
#include <winsock2.h>
#endif

#ifndef P4_TO_P8

#if 0                           /* currently unused */

/* *INDENT-OFF* */

/** Store the number of quadrants to add for complete and balance stages. */
static const int    p4est_balance_count[P4EST_DIM + 1] =
{ 5, 7, 8 };

/** Store coordinates of quadrants to add for balancing. */
static const p4est_qcoord_t p4est_balance_coord[8][P4EST_DIM] =
{ /* faces */
  { -1,  1 },
  {  2,  0 },
  {  1, -1 },
  {  0,  2 },
  /* corners */
  { -1, -1 },
  {  2, -1 },
  { -1,  2 },
  {  2,  2 }};

/** Offset for corners into p4est_balance_coord */
static const int    pbco = P4EST_FACES;

/* *INDENT-ON* */

#endif /* currently unused */

#endif /* !P4_TO_P8 */

sc_mempool_t       *
p4est_quadrant_mempool_new (void)
{
  return sc_mempool_new_zero_and_persist (sizeof (p4est_quadrant_t));
}

void
p4est_quadrant_init_data (p4est_t * p4est, p4est_topidx_t which_tree,
                          p4est_quadrant_t * quad, p4est_init_t init_fn)
{
  P4EST_ASSERT (p4est_quadrant_is_extended (quad));

  if (p4est->data_size > 0) {
    quad->p.user_data = sc_mempool_alloc (p4est->user_data_pool);
  }
  else {
    quad->p.user_data = NULL;
  }
  if (init_fn != NULL && p4est_quadrant_is_inside_root (quad)) {
    init_fn (p4est, which_tree, quad);
  }
}

void
p4est_quadrant_free_data (p4est_t * p4est, p4est_quadrant_t * quad)
{
  P4EST_ASSERT (p4est_quadrant_is_extended (quad));

  if (p4est->data_size > 0) {
    sc_mempool_free (p4est->user_data_pool, quad->p.user_data);
  }
  quad->p.user_data = NULL;
}

unsigned
p4est_quadrant_checksum (sc_array_t * quadrants,
                         sc_array_t * checkarray, size_t first_quadrant)
{
  int                 own_check;
#ifdef P4_TO_P8
  int                 level_difference;
#endif
  size_t              kz, qcount;
  unsigned            crc;
  uint32_t           *check;
  p4est_quadrant_t   *q;

  qcount = quadrants->elem_count;

  P4EST_ASSERT (quadrants->elem_size == sizeof (p4est_quadrant_t));
  P4EST_ASSERT (first_quadrant <= qcount);

  if (checkarray == NULL) {
    checkarray = sc_array_new (4);
    own_check = 1;
  }
  else {
    P4EST_ASSERT (checkarray->elem_size == 4);
    own_check = 0;
  }

  sc_array_resize (checkarray, (qcount - first_quadrant) * (P4EST_DIM + 1));
  for (kz = first_quadrant; kz < qcount; ++kz) {
    q = p4est_quadrant_array_index (quadrants, kz);
    P4EST_ASSERT (p4est_quadrant_is_extended (q));
    check =
      (uint32_t *) sc_array_index (checkarray,
                                   (kz - first_quadrant) * (P4EST_DIM + 1));
#ifndef P4_TO_P8
    check[0] = htonl ((uint32_t) q->x);
    check[1] = htonl ((uint32_t) q->y);
#else
    if (q->level <= P4EST_OLD_QMAXLEVEL) {
      /* shift the quadrant coordinates to ensure backward compatibility */
      level_difference = P4EST_MAXLEVEL - P4EST_OLD_MAXLEVEL;
      /* *INDENT-OFF* */
      check[0] =
        htonl ((q->x < 0) ? -(((uint32_t) -q->x) >> level_difference) :
                              (((uint32_t) q->x) >> level_difference));
      check[1] =
        htonl ((q->y < 0) ? -(((uint32_t) -q->y) >> level_difference) :
                              (((uint32_t) q->y) >> level_difference));
      check[2] =
        htonl ((q->z < 0) ? -(((uint32_t) -q->z) >> level_difference) :
                              (((uint32_t) q->z) >> level_difference));
      /* *INDENT-ON* */
    }
    else {
      check[0] = htonl ((uint32_t) q->x);
      check[1] = htonl ((uint32_t) q->y);
      check[2] = htonl ((uint32_t) q->z);
    }
#endif
    check[P4EST_DIM] = htonl ((uint32_t) q->level);
  }
  crc = sc_array_checksum (checkarray);

  if (own_check) {
    sc_array_destroy (checkarray);
  }

  return crc;
}

int
p4est_quadrant_in_range (const p4est_quadrant_t * fd,
                         const p4est_quadrant_t * ld,
                         const p4est_quadrant_t * quadrant)
{
  p4est_quadrant_t    quad_last_desc;

  P4EST_ASSERT (p4est_quadrant_is_valid (fd));
  P4EST_ASSERT (fd->level == P4EST_QMAXLEVEL);
  P4EST_ASSERT (p4est_quadrant_is_valid (ld));
  P4EST_ASSERT (ld->level == P4EST_QMAXLEVEL);
  P4EST_ASSERT (p4est_quadrant_compare (fd, ld) <= 0);
  P4EST_ASSERT (p4est_quadrant_is_extended (quadrant));

  /* quadrants outside of the root tree cannot be in the range */
  if (!p4est_quadrant_is_valid (quadrant)) {
    return 0;
  }

  /* check that the quadrant's first descendant is not smaller than fd */
  if (p4est_quadrant_compare (fd, quadrant) > 0 &&
      (fd->x != quadrant->x || fd->y != quadrant->y
#ifdef P4_TO_P8
       || fd->z != quadrant->z
#endif
      )) {
    return 0;
  }

  /* check that the quadrant's last descendant is not bigger than ld */
  p4est_quadrant_last_descendant (quadrant, &quad_last_desc, P4EST_QMAXLEVEL);
  if (p4est_quadrant_compare (&quad_last_desc, ld) > 0) {
    return 0;
  }

  return 1;
}

int
p4est_tree_is_sorted (p4est_tree_t * tree)
{
  size_t              iz;
  p4est_quadrant_t   *q1, *q2;
  sc_array_t         *tquadrants = &tree->quadrants;

  if (tquadrants->elem_count <= 1) {
    return 1;
  }

  q1 = p4est_quadrant_array_index (tquadrants, 0);
  for (iz = 1; iz < tquadrants->elem_count; ++iz) {
    q2 = p4est_quadrant_array_index (tquadrants, iz);
    if (p4est_quadrant_compare (q1, q2) >= 0) {
      return 0;
    }
    q1 = q2;
  }

  return 1;
}

int
p4est_tree_is_linear (p4est_tree_t * tree)
{
  size_t              iz;
  p4est_quadrant_t   *q1, *q2;
  sc_array_t         *tquadrants = &tree->quadrants;

  if (tquadrants->elem_count <= 1) {
    return 1;
  }

  q1 = p4est_quadrant_array_index (tquadrants, 0);
  for (iz = 1; iz < tquadrants->elem_count; ++iz) {
    q2 = p4est_quadrant_array_index (tquadrants, iz);
    if (p4est_quadrant_compare (q1, q2) >= 0) {
      return 0;
    }
    if (p4est_quadrant_is_ancestor (q1, q2)) {
      return 0;
    }
    q1 = q2;
  }

  return 1;
}

int
p4est_tree_is_almost_sorted (p4est_tree_t * tree, int check_linearity)
{
  size_t              iz;
  int                 face_contact1;
  int                 face_contact2;
  int                 out_axis[P4EST_DIM];
  p4est_quadrant_t   *q1, *q2;
  sc_array_t         *tquadrants = &tree->quadrants;

  if (tquadrants->elem_count <= 1) {
    return 1;
  }

  q1 = p4est_quadrant_array_index (tquadrants, 0);
  face_contact1 = 0;
  face_contact1 |= ((q1->x < 0) ? 0x01 : 0);
  face_contact1 |= ((q1->x >= P4EST_ROOT_LEN) ? 0x02 : 0);
  face_contact1 |= ((q1->y < 0) ? 0x04 : 0);
  face_contact1 |= ((q1->y >= P4EST_ROOT_LEN) ? 0x08 : 0);
#ifdef P4_TO_P8
  face_contact1 |= ((q1->z < 0) ? 0x10 : 0);
  face_contact1 |= ((q1->z >= P4EST_ROOT_LEN) ? 0x20 : 0);
#endif
  for (iz = 1; iz < tquadrants->elem_count; ++iz) {
    q2 = p4est_quadrant_array_index (tquadrants, iz);
    face_contact2 = 0;
    face_contact2 |= ((q2->x < 0) ? 0x01 : 0);
    face_contact2 |= ((q2->x >= P4EST_ROOT_LEN) ? 0x02 : 0);
    face_contact2 |= ((q2->y < 0) ? 0x04 : 0);
    face_contact2 |= ((q2->y >= P4EST_ROOT_LEN) ? 0x08 : 0);
#ifdef P4_TO_P8
    face_contact2 |= ((q2->z < 0) ? 0x10 : 0);
    face_contact2 |= ((q2->z >= P4EST_ROOT_LEN) ? 0x20 : 0);
#endif
    out_axis[0] = face_contact2 & 0x03;
    out_axis[1] = face_contact2 & 0x0c;
#ifdef P4_TO_P8
    out_axis[2] = face_contact2 & 0x30;
#endif
    if (((out_axis[0] && out_axis[1])
#ifdef P4_TO_P8
         || (out_axis[0] && out_axis[2])
         || (out_axis[1] && out_axis[2])
#endif
        ) && face_contact1 == face_contact2) {
      /* both quadrants are outside the same edge/corner and can overlap */
    }
    else {
      if (p4est_quadrant_compare (q1, q2) >= 0) {
        return 0;
      }
      if (check_linearity && p4est_quadrant_is_ancestor (q1, q2)) {
        return 0;
      }
    }
    q1 = q2;
    face_contact1 = face_contact2;
  }

  return 1;
}

int
p4est_tree_is_complete (p4est_tree_t * tree)
{
  size_t              iz;
  p4est_quadrant_t   *q1, *q2;
  sc_array_t         *tquadrants = &tree->quadrants;

  if (tquadrants->elem_count <= 1) {
    return 1;
  }

  q1 = p4est_quadrant_array_index (tquadrants, 0);
  for (iz = 1; iz < tquadrants->elem_count; ++iz) {
    q2 = p4est_quadrant_array_index (tquadrants, iz);
    if (!p4est_quadrant_is_next (q1, q2)) {
      return 0;
    }
    q1 = q2;
  }

  return 1;
}

void
p4est_tree_print (int log_priority, p4est_tree_t * tree)
{
  size_t              jz;
  int                 l, childid, comp;
  char                buffer[BUFSIZ];
  p4est_quadrant_t   *q1, *q2;
  sc_array_t         *tquadrants = &tree->quadrants;

  q1 = NULL;
  for (jz = 0; jz < tquadrants->elem_count; ++jz) {
    q2 = p4est_quadrant_array_index (tquadrants, jz);
    childid = p4est_quadrant_child_id (q2);
#ifdef P4_TO_P8
    l = snprintf (buffer, BUFSIZ, "0x%llx 0x%llx 0x%llx %d",
                  (unsigned long long) q2->x, (unsigned long long) q2->y,
                  (unsigned long long) q2->z, (int) q2->level);
#else
    l = snprintf (buffer, BUFSIZ, "0x%llx 0x%llx %d",
                  (unsigned long long) q2->x, (unsigned long long) q2->y,
                  (int) q2->level);
#endif
    if (jz > 0) {
      comp = p4est_quadrant_compare (q1, q2);
      if (comp > 0) {
        l += snprintf (buffer + l, BUFSIZ - l, " R");
      }
      else if (comp == 0) {
        l += snprintf (buffer + l, BUFSIZ - l, " I");
      }
      else {
        if (p4est_quadrant_is_sibling (q1, q2)) {
          l += snprintf (buffer + l, BUFSIZ - l, " S%d", childid);
        }
        else if (p4est_quadrant_is_parent (q1, q2)) {
          l += snprintf (buffer + l, BUFSIZ - l, " C%d", childid);
        }
        else if (p4est_quadrant_is_ancestor (q1, q2)) {
          l += snprintf (buffer + l, BUFSIZ - l, " D");
        }
        else if (p4est_quadrant_is_next (q1, q2)) {
          l += snprintf (buffer + l, BUFSIZ - l, " N%d", childid);
        }
        else {
          l += snprintf (buffer + l, BUFSIZ - l, " q%d", childid);
        }
      }
    }
    else {
      l += snprintf (buffer + l, BUFSIZ - l, " F%d", childid);
    }
    l += snprintf (buffer + l, BUFSIZ - l, "\n");
    P4EST_LOG (log_priority, buffer);
    q1 = q2;
  }
}

int
p4est_is_equal (p4est_t * p4est1, p4est_t * p4est2, int compare_data)
{
  int                 i;
  size_t              zz;
  size_t              data_size;
  p4est_topidx_t      jt;
  p4est_tree_t       *tree1, *tree2;
  p4est_quadrant_t   *q1, *q2;
  sc_array_t         *tqs1, *tqs2;

  if (p4est1->mpisize != p4est2->mpisize)
    return 0;
  if (p4est1->mpirank != p4est2->mpirank)
    return 0;
  if (compare_data) {
    if (p4est1->data_size != p4est2->data_size)
      return 0;
    data_size = p4est1->data_size;
    if (data_size == 0) {
      compare_data = 0;
    }
  }
  else {
    data_size = 0;
  }

  if (p4est1->first_local_tree != p4est2->first_local_tree)
    return 0;
  if (p4est1->last_local_tree != p4est2->last_local_tree)
    return 0;
  if (p4est1->local_num_quadrants != p4est2->local_num_quadrants)
    return 0;
  if (p4est1->global_num_quadrants != p4est2->global_num_quadrants)
    return 0;

  if (memcmp (p4est1->global_first_quadrant, p4est2->global_first_quadrant,
              (p4est1->mpisize + 1) * sizeof (p4est_gloidx_t)))
    return 0;
  if (memcmp (p4est1->global_first_position, p4est2->global_first_position,
              (p4est1->mpisize + 1) * sizeof (p4est_quadrant_t)))
    return 0;

  for (jt = p4est1->first_local_tree; jt <= p4est1->last_local_tree; ++jt) {
    tree1 = p4est_tree_array_index (p4est1->trees, jt);
    tqs1 = &tree1->quadrants;
    tree2 = p4est_tree_array_index (p4est2->trees, jt);
    tqs2 = &tree2->quadrants;

    if (!p4est_quadrant_is_equal (&tree1->first_desc, &tree2->first_desc))
      return 0;
    if (!p4est_quadrant_is_equal (&tree1->last_desc, &tree2->last_desc))
      return 0;
    if (tree1->quadrants_offset != tree2->quadrants_offset)
      return 0;

    for (i = 0; i <= P4EST_MAXLEVEL; ++i) {
      if (tree1->quadrants_per_level[i] != tree2->quadrants_per_level[i])
        return 0;
    }
    if (tree1->maxlevel != tree2->maxlevel)
      return 0;

    if (tqs1->elem_count != tqs2->elem_count)
      return 0;
    for (zz = 0; zz < tqs1->elem_count; ++zz) {
      q1 = p4est_quadrant_array_index (tqs1, zz);
      q2 = p4est_quadrant_array_index (tqs2, zz);
      if (!p4est_quadrant_is_equal (q1, q2))
        return 0;
      if (compare_data
          && memcmp (q1->p.user_data, q2->p.user_data, data_size))
        return 0;
    }
  }

  return 1;
}

int
p4est_is_valid (p4est_t * p4est)
{
  const int           num_procs = p4est->mpisize;
  const int           rank = p4est->mpirank;
  const p4est_topidx_t first_tree = p4est->first_local_tree;
  const p4est_topidx_t last_tree = p4est->last_local_tree;
  int                 i, maxlevel;
  int                 failed;
  p4est_topidx_t      jt, next_tree;
  p4est_locidx_t      nquadrants, lquadrants, perlevel;
  p4est_qcoord_t      mh = P4EST_QUADRANT_LEN (P4EST_QMAXLEVEL);
  p4est_quadrant_t   *q;
  p4est_quadrant_t    mylow, nextlow, s;
  p4est_tree_t       *tree;

  failed = 0;
  P4EST_QUADRANT_INIT (&mylow);
  P4EST_QUADRANT_INIT (&nextlow);
  P4EST_QUADRANT_INIT (&s);

  /* we crash on NULL pointers */
  P4EST_ASSERT (p4est != NULL && p4est->connectivity != NULL);

  /* check parallel environment */
  if (p4est_comm_parallel_env_is_null (p4est)) {
    P4EST_NOTICE ("p4est invalid parallel environment");
    failed = 1;
    goto failtest;
  }

  /* make sure the revision counter is legitimate */
  if (p4est->revision < 0) {
    P4EST_NOTICE ("p4est invalid revision counter\n");
    failed = 1;
    goto failtest;
  }

  /* check last item of global partition */
  if (!(p4est->global_first_position[num_procs].p.which_tree ==
        p4est->connectivity->num_trees &&
        p4est->global_first_position[num_procs].x == 0 &&
        p4est->global_first_position[num_procs].y == 0
#ifdef P4_TO_P8
        && p4est->global_first_position[num_procs].z == 0
#endif
      )) {
    P4EST_NOTICE ("p4est invalid global first position");
    failed = 1;
    goto failtest;
  }

  /* tree count and quadrant first position level */
  if (p4est->connectivity->num_trees !=
      (p4est_topidx_t) p4est->trees->elem_count) {
    P4EST_NOTICE ("p4est invalid tree count");
    failed = 1;
    goto failtest;
  }
  for (i = 0; i <= num_procs; ++i) {
    if (p4est->global_first_position[i].level != P4EST_QMAXLEVEL) {
      failed = 1;
      break;
    }
  }
  if (failed) {
    P4EST_NOTICE ("p4est invalid global first position level");
    goto failtest;
  }

  /* check first tree in global partition */
  if (first_tree < 0) {
    if (!(first_tree == -1 && last_tree == -2)) {
      P4EST_NOTICE ("p4est invalid empty tree range A");
      failed = 1;
      goto failtest;
    }
  }
  else {
    if (p4est->global_first_position[rank].p.which_tree != first_tree) {
      P4EST_NOTICE ("p4est invalid first tree\n");
      failed = 1;
      goto failtest;
    }
    mylow.x = p4est->global_first_position[rank].x;
    mylow.y = p4est->global_first_position[rank].y;
#ifdef P4_TO_P8
    mylow.z = p4est->global_first_position[rank].z;
#endif
    mylow.level = P4EST_QMAXLEVEL;
    tree = p4est_tree_array_index (p4est->trees, first_tree);
    if (tree->quadrants.elem_count > 0) {
      q = p4est_quadrant_array_index (&tree->quadrants, 0);
      if (q->x != mylow.x || q->y != mylow.y ||
#ifdef P4_TO_P8
          q->z != mylow.z ||
#endif
          0) {
        P4EST_NOTICE ("p4est invalid low quadrant\n");
        failed = 1;
        goto failtest;
      }
      if (!p4est_quadrant_in_range (&tree->first_desc, &tree->last_desc, q)) {
        P4EST_NOTICE ("p4est invalid first quadrant range\n");
        failed = 1;
        goto failtest;
      }
    }
  }

  /* check last tree in global partition */
  if (last_tree < 0) {
    if (!(first_tree == -1 && last_tree == -2)) {
      P4EST_NOTICE ("p4est invalid empty tree range B");
      failed = 1;
      goto failtest;
    }
  }
  else {
    next_tree = p4est->global_first_position[rank + 1].p.which_tree;
    if (next_tree != last_tree && next_tree != last_tree + 1) {
      P4EST_NOTICE ("p4est invalid last tree\n");
      failed = 1;
      goto failtest;
    }
    nextlow.x = p4est->global_first_position[rank + 1].x;
    nextlow.y = p4est->global_first_position[rank + 1].y;
#ifdef P4_TO_P8
    nextlow.z = p4est->global_first_position[rank + 1].z;
#endif
    nextlow.level = P4EST_QMAXLEVEL;
    if (next_tree == last_tree + 1) {
      if (nextlow.x != 0 || nextlow.y != 0
#ifdef P4_TO_P8
          || nextlow.z != 0
#endif
        ) {
        P4EST_NOTICE ("p4est invalid next coordinates\n");
        failed = 1;
        goto failtest;
      }
    }
    tree = p4est_tree_array_index (p4est->trees, last_tree);
    if (tree->quadrants.elem_count > 0) {
      q =
        p4est_quadrant_array_index (&tree->quadrants,
                                    tree->quadrants.elem_count - 1);
      if (next_tree == last_tree) {
        if (!p4est_quadrant_is_next (q, &nextlow)) {
          P4EST_NOTICE ("p4est invalid next quadrant\n");
          failed = 1;
          goto failtest;
        }
      }
      else {
        p4est_quadrant_last_descendant (q, &s, P4EST_QMAXLEVEL);
        if (s.x + mh != P4EST_ROOT_LEN || s.y + mh != P4EST_ROOT_LEN ||
#ifdef P4_TO_P8
            s.z + mh != P4EST_ROOT_LEN ||
#endif
            0) {
          P4EST_NOTICE ("p4est invalid last quadrant\n");
          failed = 1;
          goto failtest;
        }
      }
      if (!p4est_quadrant_in_range (&tree->first_desc, &tree->last_desc, q)) {
        P4EST_NOTICE ("p4est invalid last quadrant range\n");
        failed = 1;
        goto failtest;
      }
    }
  }

  /* check individual trees */
  lquadrants = 0;
  for (jt = 0; jt < (p4est_topidx_t) p4est->trees->elem_count; ++jt) {
    tree = p4est_tree_array_index (p4est->trees, jt);
    if (tree->quadrants_offset != lquadrants) {
      P4EST_NOTICE ("p4est invalid quadrants offset\n");
      failed = 1;
      goto failtest;
    }
    if (!p4est_tree_is_complete (tree)) {
      P4EST_NOTICE ("p4est invalid not complete\n");
      failed = 1;
      goto failtest;
    }
    if (tree->quadrants.elem_count > 0) {
      if (jt < p4est->first_local_tree || jt > p4est->last_local_tree) {
        P4EST_NOTICE ("p4est invalid outside count\n");
        failed = 1;
        goto failtest;
      }
      q = p4est_quadrant_array_index (&tree->quadrants, 0);
      p4est_quadrant_first_descendant (q, &s, P4EST_QMAXLEVEL);
      if (!p4est_quadrant_is_equal (&s, &tree->first_desc)) {
        P4EST_NOTICE ("p4est invalid first tree descendant\n");
        failed = 1;
        goto failtest;
      }
      q =
        p4est_quadrant_array_index (&tree->quadrants,
                                    tree->quadrants.elem_count - 1);
      p4est_quadrant_last_descendant (q, &s, P4EST_QMAXLEVEL);
      if (!p4est_quadrant_is_equal (&s, &tree->last_desc)) {
        P4EST_NOTICE ("p4est invalid last tree descendant\n");
        failed = 1;
        goto failtest;
      }
    }
    else {
      P4EST_QUADRANT_INIT (&s);
      if (s.level != tree->first_desc.level ||
          s.level != tree->last_desc.level) {
        P4EST_NOTICE ("p4est invalid empty descendant\n");
        failed = 1;
        goto failtest;
      }
    }

    maxlevel = 0;
    nquadrants = 0;
    for (i = 0; i <= P4EST_QMAXLEVEL; ++i) {
      perlevel = tree->quadrants_per_level[i];
      if (perlevel < 0) {
        failed = 1;
        break;
      }
      nquadrants += perlevel;   /* same type */
      if (perlevel > 0) {
        maxlevel = i;
      }
    }
    if (!failed) {
      for (; i <= P4EST_MAXLEVEL; ++i) {
        if (tree->quadrants_per_level[i] != -1) {
          failed = 1;
          break;
        }
      }
      lquadrants += nquadrants; /* same type */
    }
    if (failed || maxlevel != (int) tree->maxlevel) {
      P4EST_NOTICE ("p4est invalid tree level\n");
      failed = 1;
      goto failtest;
    }

    if (nquadrants != (p4est_locidx_t) tree->quadrants.elem_count) {
      P4EST_NOTICE ("p4est invalid tree quadrant count\n");
      failed = 1;
      goto failtest;
    }
  }

  if (lquadrants != p4est->local_num_quadrants) {
    P4EST_NOTICE ("p4est invalid local quadrant count\n");
    failed = 1;
    goto failtest;
  }

  if (p4est->global_first_quadrant[0] != 0 ||
      p4est->global_first_quadrant[num_procs] !=
      p4est->global_num_quadrants) {
    P4EST_NOTICE ("p4est invalid global quadrant index\n");
    failed = 1;
    goto failtest;
  }

failtest:
  return !p4est_comm_sync_flag (p4est, failed, sc_MPI_BOR);
}

/* here come the heavyweight algorithms */
#ifndef P4_TO_P8
/* which face of the center quad touches this insul */
static const int    insul_to_f[9] = { -1, 2, -1, 0, -1, 1, -1, 3, -1 };

/* which corner of the center quad touches this insul */
static const int    insul_to_c[9] = { 0, -1, 1, -1, -1, -1, 2, -1, 3 };
#else
/* which face of the center quad touches this insul */
/* *INDENT-OFF* */
static const int insul_to_f[27] =
{-1, -1, -1, -1, 4, -1, -1, -1, -1,
 -1, 2, -1, 0, -1, 1, -1, 3, -1,
 -1, -1, -1, -1, 5, -1, -1, -1, -1};
/* which corner of the center quad touches this insul */
static const int insul_to_c[27] =
{0, -1, 1, -1, -1, -1, 2, -1, 3,
 -1, -1, -1, -1, -1, -1, -1, -1, -1,
 4, -1, 5, -1, -1, -1, 6, -1, 7};
/* which edge of the center quad touches this insul */
static const int insul_to_e[27] =
{-1, 0, -1, 4, -1, 5, -1, 1, -1,
  8, -1, 9, -1, -1, -1, 10, -1, 11,
  -1, 2, -1, 6, -1, 7, -1, 3, -1};
/* *INDENT-ON* */
#endif

static void
p4est_output_array_push_data (sc_array_t * out, const p4est_quadrant_t * src,
                              p4est_topidx_t which_tree)
{
  p4est_quadrant_t   *outq = p4est_quadrant_array_push (out);

  p4est_quadrant_sibling (src, outq, 0);
  outq->p.piggy2.which_tree = which_tree;
  /* *INDENT-OFF* HORRIBLE indent bug */
  outq->p.piggy2.from_tree = (p4est_topidx_t) -1;
  /* *INDENT-ON* */
}

void
p4est_tree_compute_overlap (p4est_t * p4est, sc_array_t * in,
                            sc_array_t * out, p4est_connect_type_t balance,
                            sc_array_t * borders, sc_array_t * inseeds)
{
  int                 k, l, m, which;
  int                 face, corner, level;
  int                 f = -1, c = -1;
  int                 ftransform[P4EST_FTRANSFORM];
  int                 face_axis[3];     /* 3 not P4EST_DIM */
  int                 contact_face_only, contact_edge_only;
  int                 inter_tree, outface[P4EST_FACES];
  size_t              iz, jz, kz, ctree;
  size_t              treecount, incount, seedcount;
  size_t              guess, split;
  ssize_t             first_index, last_index, js;
  p4est_topidx_t      qtree, ntree, first_tree, ftree = -1;
  p4est_qcoord_t      qh;
  p4est_quadrant_t    fd, ld, tempq, ins[P4EST_INSUL];
  p4est_quadrant_t   *treefd, *treeld;
  p4est_quadrant_t   *tq, *s, *u;
  p4est_quadrant_t   *inq;
  p4est_tree_t       *tree;
  p4est_connectivity_t *conn = p4est->connectivity;
#ifdef P4_TO_P8
  int                 edge;
  int                 e = -1;
  size_t              etree;
  p8est_edge_info_t   ei;
  p8est_edge_transform_t *et;
  sc_array_t         *eta;
  p4est_quadrant_t    tempq1, tempq2;
#endif
  p4est_corner_info_t ci;
  p4est_corner_transform_t *ct;
  sc_array_t         *cta;
  sc_array_t         *tquadrants;
  sc_array_t         *seeds = NULL;
  p4est_quadrant_t   *neigharray[P4EST_CHILDREN];
  size_t              nneigh = -1;

  P4EST_QUADRANT_INIT (&fd);
  P4EST_QUADRANT_INIT (&ld);
  P4EST_QUADRANT_INIT (&tempq);
#ifdef P4_TO_P8
  P4EST_QUADRANT_INIT (&tempq1);
  P4EST_QUADRANT_INIT (&tempq2);
#endif
  for (which = 0; which < P4EST_INSUL; ++which) {
    P4EST_QUADRANT_INIT (&ins[which]);
  }
#ifdef P4_TO_P8
  eta = &ei.edge_transforms;
  sc_array_init (eta, sizeof (p8est_edge_transform_t));
#endif
  cta = &ci.corner_transforms;
  sc_array_init (cta, sizeof (p4est_corner_transform_t));

  /* assign incoming quadrant count */
  incount = in->elem_count;

  /* initialize the tracking of trees */
  qtree = -1;
  tree = NULL;
  treefd = treeld = NULL;
  tquadrants = NULL;
  treecount = -1;

  seeds = sc_array_new (sizeof (p4est_quadrant_t));
  first_tree = p4est->first_local_tree;

  /* loop over input list of quadrants */
  for (iz = 0; iz < incount; ++iz) {
    inq = p4est_quadrant_array_index (in, iz);

    P4EST_ASSERT (inq->p.piggy2.from_tree >= 0 &&
                  inq->p.piggy2.from_tree < p4est->connectivity->num_trees);
    ftree = inq->p.piggy2.from_tree;
    nneigh = 0;

    /* potentially grab new tree */
    if (inq->p.piggy2.which_tree != qtree) {
      P4EST_ASSERT (qtree < inq->p.piggy2.which_tree);
      qtree = inq->p.piggy2.which_tree;

      tree = p4est_tree_array_index (p4est->trees, qtree);
      treefd = &tree->first_desc;
      treeld = &tree->last_desc;
      if (borders == NULL) {
        tquadrants = &tree->quadrants;
      }
      else {
        tquadrants = (sc_array_t *) sc_array_index_int (borders, (int)
                                                        (qtree - first_tree));
      }
      treecount = tquadrants->elem_count;
      P4EST_ASSERT (treecount > 0);
    }

    inter_tree = 0;
    ntree = -1;
    face = corner = -1;
#ifdef P4_TO_P8
    edge = -1;
    ei.iedge = -1;
    et = NULL;
#endif
    ci.icorner = -1;
    ct = NULL;
    contact_face_only = contact_edge_only = 0;
    if (!p4est_quadrant_is_inside_root (inq)) {
      /* this quadrant comes from a different tree */
      P4EST_ASSERT (p4est_quadrant_is_extended (inq));
      inter_tree = 1;
      outface[0] = (inq->x < 0);
      outface[1] = (inq->x >= P4EST_ROOT_LEN);
      face_axis[0] = outface[0] || outface[1];
      outface[2] = (inq->y < 0);
      outface[3] = (inq->y >= P4EST_ROOT_LEN);
      face_axis[1] = outface[2] || outface[3];
#ifndef P4_TO_P8
      face_axis[2] = 0;
#else
      outface[4] = (inq->z < 0);
      outface[5] = (inq->z >= P4EST_ROOT_LEN);
      face_axis[2] = outface[4] || outface[5];
#endif
      if (!face_axis[1] && !face_axis[2]) {
        contact_face_only = 1;
        face = 0 + outface[1];
      }
      else if (!face_axis[0] && !face_axis[2]) {
        contact_face_only = 1;
        face = 2 + outface[3];
      }
#ifdef P4_TO_P8
      else if (!face_axis[0] && !face_axis[1]) {
        contact_face_only = 1;
        face = 4 + outface[5];
      }
      else if (!face_axis[0]) {
        contact_edge_only = 1;
        edge = 0 + 2 * outface[5] + outface[3];
      }
      else if (!face_axis[1]) {
        contact_edge_only = 1;
        edge = 4 + 2 * outface[5] + outface[1];
      }
      else if (!face_axis[2]) {
        contact_edge_only = 1;
        edge = 8 + 2 * outface[3] + outface[1];
      }
#endif
      if (contact_face_only) {
        P4EST_ASSERT (!contact_edge_only && face >= 0 && face < P4EST_FACES);
        P4EST_ASSERT (outface[face]);
        ntree = p4est_find_face_transform (conn, qtree, face, ftransform);
        P4EST_ASSERT (ntree >= 0);
      }
#ifdef P4_TO_P8
      else if (contact_edge_only) {
        P4EST_ASSERT (!contact_face_only && edge >= 0 && edge < P8EST_EDGES);
        p8est_find_edge_transform (conn, qtree, edge, &ei);
        P4EST_ASSERT (ei.edge_transforms.elem_count > 0);
      }
#endif
      else {
        /* outside across a corner */
        P4EST_ASSERT (face_axis[0] && face_axis[1]);
        corner = outface[1] + 2 * outface[3];
#ifdef P4_TO_P8
        P4EST_ASSERT (face_axis[2]);
        corner += 4 * outface[5];
#endif
        P4EST_ASSERT (p4est_quadrant_touches_corner (inq, corner, 0));
        p4est_find_corner_transform (conn, qtree, corner, &ci);
        P4EST_ASSERT (ci.corner_transforms.elem_count > 0);
      }
    }
    qh = P4EST_QUADRANT_LEN (inq->level);

    /* loop over the insulation layer of inq */
#ifdef P4_TO_P8
    for (m = 0; m < 3; ++m) {
#if 0
    }
#endif
#else
    m = 0;
#endif
    for (k = 0; k < 3; ++k) {
      for (l = 0; l < 3; ++l) {
        which = m * 9 + k * 3 + l;      /* 2D: 0..8, 3D: 0..26 */

        /* exclude myself from the queries */
        if (which == P4EST_INSUL / 2) {
          continue;
        }
        s = &ins[which];
        *s = *inq;
        s->x += (l - 1) * qh;
        s->y += (k - 1) * qh;
#ifdef P4_TO_P8
        s->z += (m - 1) * qh;
#endif
        if ((s->x < 0 || s->x >= P4EST_ROOT_LEN) ||
            (s->y < 0 || s->y >= P4EST_ROOT_LEN) ||
#ifdef P4_TO_P8
            (s->z < 0 || s->z >= P4EST_ROOT_LEN) ||
#endif
            0) {
          /* this quadrant is outside this tree, no overlap */
          continue;
        }
        p4est_quadrant_first_descendant (s, &fd, P4EST_QMAXLEVEL);
        p4est_quadrant_last_descendant (s, &ld, P4EST_QMAXLEVEL);

        /* skip this insulation quadrant if there is no overlap */
        if (p4est_quadrant_compare (&ld, treefd) < 0 ||
            p4est_quadrant_compare (treeld, &fd) < 0) {
          continue;
        }

        /* Find last quadrant in tree <= ld */
        guess = treecount / 2;
        if (p4est_quadrant_compare (treeld, &ld) <= 0) {
          /* the last tree quadrant overlaps an insulation quadrant */
          last_index = (ssize_t) treecount - 1;
        }
        else {
          /* do a binary search for the highest tree quadrant <= ld */
          last_index = p4est_find_higher_bound (tquadrants, &ld, guess);
          if (last_index < 0) {
            SC_ABORT_NOT_REACHED ();
          }
          guess = (size_t) last_index;
        }

        if (p4est_quadrant_compare (&fd, treefd) < 0) {
          /* the first tree quadrant overlaps an insulation quadrant */
          first_index = 0;
        }
        else {
          /* Do a binary search for the lowest tree quadrant >= s.
             Does not accept a strict ancestor of s, which is on purpose. */
          first_index = p4est_find_lower_bound (tquadrants, s, guess);
        }

        if (first_index < 0 || first_index > last_index ||
            p4est_quadrant_compare (&fd, &ld) == 0) {
          /* The only possibility is that a quadrant tq larger than s
           * contains s, or that tq and s are of smallest possible size */
          tq = p4est_quadrant_array_index (tquadrants, last_index);
          P4EST_ASSERT (p4est_quadrant_is_ancestor (tq, s) ||
                        (p4est_quadrant_is_equal (tq, s) &&
                         tq->level == P4EST_QMAXLEVEL));
          if (tq->level < s->level - 1) {
            for (kz = 0; kz < nneigh; kz++) {
              if (neigharray[kz] == tq) {
                break;
              }
            }
            /* if this neighbor hasn't been calculated */
            if (kz == nneigh) {
              /* we should check to see if inq causes a split to tq */
              split = p4est_balance_seeds (inq, tq, balance, seeds);

              if (split) {
                seedcount = seeds->elem_count;
                for (jz = 0; jz < seedcount; jz++) {
                  u = p4est_quadrant_array_index (seeds, jz);
                  P4EST_ASSERT (p4est_quadrant_is_ancestor (tq, u));
                  p4est_output_array_push_data (inseeds, u, qtree);
                }
              }
              P4EST_ASSERT (nneigh < P4EST_CHILDREN - 1);
              neigharray[nneigh++] = tq;
            }
          }
          continue;
        }

        /* figure out the relationship of s to inq */
        f = insul_to_f[which];
#ifdef P4_TO_P8
        e = insul_to_e[which];
#endif
        c = insul_to_c[which];

        level = inq->level + 1;

        /* copy relevant quadrants into out */
        for (js = first_index; js <= last_index; ++js) {
          tq = p4est_quadrant_array_index (tquadrants, (size_t) js);
          if (tq->level <= level) {
            continue;
          }
          if (f >= 0) {
            p4est_quadrant_face_neighbor (tq, f ^ 1, &tempq);
            if (p4est_quadrant_is_ancestor (inq, &tempq)) {
              continue;
            }
            split = p4est_balance_seeds_face (tq, inq, f, balance, seeds);
          }
#ifdef P4_TO_P8
          else if (e >= 0) {
            p8est_quadrant_edge_neighbor (tq, e ^ 3, &tempq);
            if (p4est_quadrant_is_ancestor (inq, &tempq)) {
              continue;
            }
            split = p8est_balance_seeds_edge (tq, inq, e, balance, seeds);
          }
#endif
          else {
            P4EST_ASSERT (c >= 0);
            p4est_quadrant_corner_neighbor (tq, (P4EST_CHILDREN - 1) ^ c,
                                            &tempq);
            if (p4est_quadrant_is_ancestor (inq, &tempq)) {
              continue;
            }
            split = p4est_balance_seeds_corner (tq, inq, c, balance, seeds);
          }
          if (split) {
            seedcount = seeds->elem_count;
            for (jz = 0; jz < seedcount; jz++) {
              u = p4est_quadrant_array_index (seeds, jz);
              P4EST_ASSERT (p4est_quadrant_is_ancestor (inq, u));
              if (inter_tree) {
                if (contact_face_only) {
                  P4EST_ASSERT (!contact_edge_only);
                  P4EST_ASSERT (ntree == ftree);
                  p4est_quadrant_transform_face (u, &tempq, ftransform);
                  p4est_output_array_push_data (out, &tempq, ntree);
                }
#ifdef P4_TO_P8
                else if (contact_edge_only) {
                  P4EST_ASSERT (inq->pad16 >= 0 && inq->pad16 < P8EST_EDGES);
                  for (etree = 0; etree < eta->elem_count; ++etree) {
                    et = p8est_edge_array_index (eta, etree);
                    if (et->ntree == ftree && et->nedge == inq->pad16) {
                      p8est_quadrant_transform_edge (u, &tempq, &ei, et, 1);
                      p4est_output_array_push_data (out, &tempq, et->ntree);
                    }
                  }
                  et = NULL;
                }
#endif
                else {
                  P4EST_ASSERT (corner >= 0);
                  P4EST_ASSERT (inq->pad16 >= 0 &&
                                inq->pad16 < P4EST_CHILDREN);
                  for (ctree = 0; ctree < cta->elem_count; ++ctree) {
                    ct = p4est_corner_array_index (cta, ctree);
                    if (ct->ntree == ftree && ct->ncorner == inq->pad16) {
                      p4est_quadrant_transform_corner (u, (int) ct->ncorner,
                                                       1);
                      p4est_output_array_push_data (out, u, ct->ntree);
                    }
                  }
                  ct = NULL;
                }
              }
              else {
                p4est_output_array_push_data (out, u, qtree);
              }

              if (c >= 0) {
                level = SC_MAX (level, u->level);
              }
            }
          }
        }
      }
    }
#ifdef P4_TO_P8
#if 0
    {
#endif
    }
#endif
  }

#ifdef P4_TO_P8
  sc_array_reset (eta);
#endif
  sc_array_reset (cta);

  sc_array_destroy (seeds);
}

void
p4est_tree_uniqify_overlap (sc_array_t * out)
{
  size_t              iz, jz;
  size_t              outcount, dupcount, olcount;
  p4est_quadrant_t   *q, *p, tempq;

  outcount = out->elem_count;
  if (outcount == 0) {
    return;
  }

  /* sort array and remove duplicates */
  sc_array_sort (out, p4est_quadrant_compare_piggy);
  dupcount = olcount = 0;
  iz = 0;                       /* read counter */
  jz = 0;                       /* write counter */
  q = NULL;
  for (iz = 0; iz < outcount; iz++) {
    p = p4est_quadrant_array_index (out, iz);
    P4EST_ASSERT (p4est_quadrant_child_id (p) == 0);
    if (q != NULL && q->p.piggy2.which_tree == p->p.piggy2.which_tree) {
      p4est_nearest_common_ancestor (p, q, &tempq);
      if (tempq.level >= SC_MIN (q->level, p->level) - 1) {
        if (p->level > q->level) {
          olcount++;
          *q = *p;
        }
        else {
          P4EST_ASSERT (p->level == q->level);
          dupcount++;
        }
        continue;
      }
    }
    if (iz == jz) {
      q = p;
    }
    else {
      q = p4est_quadrant_array_index (out, jz);
      *q = *p;
    }
    jz++;
  }
  P4EST_ASSERT (jz + olcount + dupcount == outcount);
  sc_array_resize (out, jz);
}

size_t
p4est_tree_remove_nonowned (p4est_t * p4est, p4est_topidx_t which_tree)
{
  int                 full_tree[2];
  size_t              zz, incount, prev_good, removed;
#ifdef P4EST_ENABLE_DEBUG
  const p4est_topidx_t first_tree = p4est->first_local_tree;
  const p4est_topidx_t last_tree = p4est->last_local_tree;
#endif
  const p4est_quadrant_t *first_pos, *next_pos;
  p4est_quadrant_t   *q1, *q2;
  p4est_quadrant_t    ld;
  p4est_tree_t       *tree;
  sc_array_t         *quadrants;

  P4EST_ASSERT (first_tree <= which_tree && which_tree <= last_tree);
  tree = p4est_tree_array_index (p4est->trees, which_tree);
  P4EST_ASSERT (p4est_tree_is_almost_sorted (tree, 0));

  quadrants = &tree->quadrants;
  incount = quadrants->elem_count;
  if (incount == 0) {
    return 0;
  }

  P4EST_QUADRANT_INIT (&ld);
  p4est_comm_tree_info (p4est, which_tree, full_tree, NULL,
                        &first_pos, &next_pos);

  /* q1 is the last known good quadrant */
  q1 = NULL;
  prev_good = incount;
  removed = 0;
  for (zz = 0; zz < incount; ++zz) {
    q2 = p4est_quadrant_array_index (quadrants, zz);
    P4EST_ASSERT (p4est_quadrant_is_extended (q2));
    if (!p4est_quadrant_is_inside_root (q2) ||
        (!full_tree[0] &&
         (p4est_quadrant_compare (q2, first_pos) < 0 &&
          (q2->x != first_pos->x || q2->y != first_pos->y
#ifdef P4_TO_P8
           || q2->z != first_pos->z
#endif
          ))) ||
        (!full_tree[1] &&
         (p4est_quadrant_last_descendant (q2, &ld, P4EST_QMAXLEVEL),
          p4est_quadrant_compare (next_pos, &ld) <= 0))) {
      /* quadrant is outside of the unit square
         or at least partially outside of the tree bounds */
      --tree->quadrants_per_level[q2->level];
      p4est_quadrant_free_data (p4est, q2);
      ++removed;
#ifdef P4EST_ENABLE_DEBUG
      P4EST_QUADRANT_INIT (q2);
#endif
    }
    else {
      if (prev_good == incount) {
        /* this is the first good quadrant we find */
        prev_good = 0;
      }
      else {
        /* q1 at prev_good was the last known good */
        ++prev_good;
      }
      P4EST_ASSERT (prev_good <= zz);
      q1 = p4est_quadrant_array_index (quadrants, prev_good);
      if (zz > prev_good) {
        *q1 = *q2;
#ifdef P4EST_ENABLE_DEBUG
        P4EST_QUADRANT_INIT (q2);
#endif
      }
    }
  }

  if (prev_good == incount) {
    P4EST_ASSERT (removed == incount);
    incount = 0;
  }
  else {
    P4EST_ASSERT (prev_good + 1 + removed == incount);
    incount = prev_good + 1;
    q1 = p4est_quadrant_array_index (quadrants, 0);
  }
  sc_array_resize (quadrants, incount);

  tree->maxlevel = 0;
  for (zz = 0; zz < incount; ++zz) {
    q1 = p4est_quadrant_array_index (quadrants, zz);
    P4EST_ASSERT (p4est_quadrant_is_valid (q1));
    tree->maxlevel = (int8_t) SC_MAX (tree->maxlevel, q1->level);
  }

  P4EST_ASSERT (p4est_tree_is_sorted (tree));

  return removed;
}

void
p4est_complete_region (p4est_t * p4est,
                       const p4est_quadrant_t * q1,
                       int include_q1,
                       const p4est_quadrant_t * q2,
                       int include_q2,
                       p4est_tree_t * tree,
                       p4est_topidx_t which_tree, p4est_init_t init_fn)
{
#ifdef P4EST_ENABLE_DEBUG
  size_t              quadrant_pool_size, data_pool_size;
#endif

  p4est_tree_t       *R;
  sc_list_t          *W;

  p4est_quadrant_t    a = *q1;
  p4est_quadrant_t    b = *q2;

  p4est_quadrant_t    Afinest;
  p4est_quadrant_t   *c0, *c1, *c2, *c3;
#ifdef P4_TO_P8
  p4est_quadrant_t   *c4, *c5, *c6, *c7;
#endif

  sc_array_t         *quadrants;
  sc_mempool_t       *quadrant_pool = p4est->quadrant_pool;

  p4est_quadrant_t   *w;
  p4est_quadrant_t   *r;

  int                 comp;
  int                 maxlevel = 0;
  p4est_locidx_t     *quadrants_per_level;

  P4EST_QUADRANT_INIT (&Afinest);

  W = sc_list_new (NULL);
  R = tree;

  /* needed for sanity check */
#ifdef P4EST_ENABLE_DEBUG
  quadrant_pool_size = p4est->quadrant_pool->elem_count;
  data_pool_size = 0;
  if (p4est->user_data_pool != NULL) {
    data_pool_size = p4est->user_data_pool->elem_count;
  }
#endif

  quadrants = &R->quadrants;
  quadrants_per_level = R->quadrants_per_level;

  /* Assert that we have an empty tree */
  P4EST_ASSERT (quadrants->elem_count == 0);

  comp = p4est_quadrant_compare (&a, &b);
  /* Assert that a<b */
  P4EST_ASSERT (comp < 0);

  /* R <- R + a */
  if (include_q1) {
    r = p4est_quadrant_array_push_copy (quadrants, &a);
    p4est_quadrant_init_data (p4est, which_tree, r, init_fn);
    maxlevel = SC_MAX ((int) r->level, maxlevel);
    ++quadrants_per_level[r->level];
  }

  if (comp < 0) {
    /* W <- C(A_{finest}(a,b)) */
    p4est_nearest_common_ancestor (&a, &b, &Afinest);

    c0 = p4est_quadrant_mempool_alloc (quadrant_pool);
    c1 = p4est_quadrant_mempool_alloc (quadrant_pool);
    c2 = p4est_quadrant_mempool_alloc (quadrant_pool);
    c3 = p4est_quadrant_mempool_alloc (quadrant_pool);
#ifdef P4_TO_P8
    c4 = p4est_quadrant_mempool_alloc (quadrant_pool);
    c5 = p4est_quadrant_mempool_alloc (quadrant_pool);
    c6 = p4est_quadrant_mempool_alloc (quadrant_pool);
    c7 = p4est_quadrant_mempool_alloc (quadrant_pool);

    p8est_quadrant_children (&Afinest, c0, c1, c2, c3, c4, c5, c6, c7);
#else
    p4est_quadrant_children (&Afinest, c0, c1, c2, c3);
#endif

    (void) sc_list_append (W, c0);
    (void) sc_list_append (W, c1);
    (void) sc_list_append (W, c2);
    (void) sc_list_append (W, c3);
#ifdef P4_TO_P8
    (void) sc_list_append (W, c4);
    (void) sc_list_append (W, c5);
    (void) sc_list_append (W, c6);
    (void) sc_list_append (W, c7);
#endif

    /* for each w in W */
    while (W->elem_count > 0) {
      w = p4est_quadrant_list_pop (W);

      /* if (a < w < b) and (w not in {A(b)}) */
      if (((p4est_quadrant_compare (&a, w) < 0) &&
           (p4est_quadrant_compare (w, &b) < 0)
          ) && !p4est_quadrant_is_ancestor (w, &b)
        ) {
        /* R <- R + w */
        r = p4est_quadrant_array_push_copy (quadrants, w);
        p4est_quadrant_init_data (p4est, which_tree, r, init_fn);
        maxlevel = SC_MAX ((int) r->level, maxlevel);
        ++quadrants_per_level[r->level];
      }
      /* else if (w in {{A(a)}, {A(b)}}) */
      else if (p4est_quadrant_is_ancestor (w, &a)
               || p4est_quadrant_is_ancestor (w, &b)) {
        /* W <- W + C(w) */
        c0 = p4est_quadrant_mempool_alloc (quadrant_pool);
        c1 = p4est_quadrant_mempool_alloc (quadrant_pool);
        c2 = p4est_quadrant_mempool_alloc (quadrant_pool);
        c3 = p4est_quadrant_mempool_alloc (quadrant_pool);
#ifdef P4_TO_P8
        c4 = p4est_quadrant_mempool_alloc (quadrant_pool);
        c5 = p4est_quadrant_mempool_alloc (quadrant_pool);
        c6 = p4est_quadrant_mempool_alloc (quadrant_pool);
        c7 = p4est_quadrant_mempool_alloc (quadrant_pool);

        p8est_quadrant_children (w, c0, c1, c2, c3, c4, c5, c6, c7);
#else
        p4est_quadrant_children (w, c0, c1, c2, c3);
#endif

#ifdef P4_TO_P8
        (void) sc_list_prepend (W, c7);
        (void) sc_list_prepend (W, c6);
        (void) sc_list_prepend (W, c5);
        (void) sc_list_prepend (W, c4);
#endif
        (void) sc_list_prepend (W, c3);
        (void) sc_list_prepend (W, c2);
        (void) sc_list_prepend (W, c1);
        (void) sc_list_prepend (W, c0);
      }

      /* W <- W - w */
      sc_mempool_free (quadrant_pool, w);
    }                           /* end for */

    /* R <- R + b */
    if (include_q2) {
      r = p4est_quadrant_array_push_copy (quadrants, &b);
      p4est_quadrant_init_data (p4est, which_tree, r, init_fn);
      maxlevel = SC_MAX ((int) r->level, maxlevel);
      ++quadrants_per_level[r->level];
    }
  }

  R->maxlevel = (int8_t) maxlevel;

  P4EST_ASSERT (W->first == NULL && W->last == NULL);
  sc_list_destroy (W);

  P4EST_ASSERT (p4est_tree_is_complete (R));
  P4EST_ASSERT (quadrant_pool_size == p4est->quadrant_pool->elem_count);
  if (p4est->user_data_pool != NULL) {
    P4EST_ASSERT (data_pool_size + quadrants->elem_count ==
                  p4est->user_data_pool->elem_count);
  }
}

static int
p4est_quadrant_disjoint_parent (const void *a, const void *b)
{
  const p4est_quadrant_t *q = (p4est_quadrant_t *) a;
  const p4est_quadrant_t *r = (p4est_quadrant_t *) b;
  int8_t              level = SC_MIN (q->level - 1, r->level - 1);
  p4est_qcoord_t      mask = P4EST_QUADRANT_MASK (level);

  if (((q->x ^ r->x) & mask) || ((q->y ^ r->y) & mask)
#ifdef P4_TO_P8
      || ((q->z ^ r->z) & mask)
#endif
      || 0) {
    return p4est_quadrant_compare (a, b);
  }

  return 0;
}

/** Complete/balance a region of an tree.
 *
 * \param [in] inlist             List of quadrants to consider: should be
 *                                sorted and reduced, i.e., every quadrant
 *                                should have child_id == 0.
 * \param [in]     dom            Least common ancestor of all quadrants in
 *                                \a inlist.
 * \param [in]     bound          The number of quadrants in a neighborhood to
 *                                consider when balancing.
 *                                bound = 1 : just the quadrant itself, i.e.,
 *                                            completion.
 *                                bound = P4EST_DIM + 1 : face balance
 *                                bound = 2**P4EST_DIM  : full balance
 *                                bound = 2**P4EST_DIM - 1 : edge balance
 * \param [in/out] qpool          quadrant pool for temporary quadrants
 * \param [in/out] list_alloc     list mempool for hash tables
 * \param [in/out] out            the sorted, complete, balance quadrants in
 *                                the region will be appended to out
 * \param [in]     first_desc     the first quadrant defining the start of the
 *                                region.  if NULL, the region is understood
 *                                to start with the first descendant of \a
 *                                dom.
 * \param [in]     last_desc      the last quadrant defining the start of the
 *                                region.  if NULL, the region is understood
 *                                to end with the last descendant of \a
 *                                dom.
 * \param [in/out] count_in       If not NULL, points to an accumulator for
 *                                the number of times the balance algorithm
 *                                tries to insert a quadrant that already
 *                                exists
 * \param [in/out] count_out      If not NULL, points to an accumulator for
 *                                the number of times the balance algorithm
 *                                tries to duplicate the insertion of a new
 *                                quadrant
 * \param [in/out] count_an       If not NULL, points to an accumulator for
 *                                the number of times the balance algorithm
 *                                tries to insert the ancestor of an existing
 *                                quadrant
 */
static void
p4est_complete_or_balance_kernel (sc_array_t * inlist,
                                  p4est_quadrant_t * dom,
                                  int bound,
                                  sc_mempool_t * qpool,
                                  sc_mempool_t * list_alloc,
                                  sc_array_t * out,
                                  p4est_quadrant_t * first_desc,
                                  p4est_quadrant_t * last_desc,
                                  size_t *count_in, size_t *count_out,
                                  size_t *count_an)
{
  int                 inserted;
  size_t              iz, jz;
  size_t              incount, ocount;
#ifdef P4EST_ENABLE_DEBUG
  size_t              quadrant_pool_size;
  sc_array_t          outview;
  p4est_lid_t         lid, one;
  p4est_quadrant_t    ld_old;
#endif
  size_t              count_already_inlist, count_already_outlist;
  size_t              count_ancestor_inlist;
  p4est_quadrant_t   *q, *p, *r;
  int                 minlevel = dom->level + 1, maxlevel;
  int                 sid, pid;
  int                 duplicate = 1;
  int                 precluded = 2;
  int                 l;
  void              **vlookup;
  ssize_t             srindex, si;
  p4est_qcoord_t      ph;
  p4est_quadrant_t   *qalloc, *qlookup, **qpointer;
  p4est_quadrant_t    par, tempq, tempp, fd, ld;
  sc_array_t         *olist;
  sc_hash_t          *hash[P4EST_MAXLEVEL + 1];
  sc_array_t          outlist[P4EST_MAXLEVEL + 1];

  P4EST_QUADRANT_INIT (&par);
  par.p.user_int = 0;
  P4EST_QUADRANT_INIT (&tempq);
  P4EST_QUADRANT_INIT (&tempp);
  P4EST_QUADRANT_INIT (&fd);
  P4EST_QUADRANT_INIT (&ld);

#ifdef P4EST_ENABLE_DEBUG
  quadrant_pool_size = qpool->elem_count;
#endif

  count_already_inlist = count_already_outlist = 0;
  count_ancestor_inlist = 0;

#ifdef P4EST_ENABLE_DEBUG
  /* to increment linear id */
  p4est_lid_set_one (&one);
#endif

  incount = inlist->elem_count;
  maxlevel = minlevel;
  for (jz = 0; jz < incount; jz++) {
    q = p4est_quadrant_array_index (inlist, jz);
    q->p.user_int = 0;
    maxlevel = SC_MAX (maxlevel, q->level);
    P4EST_ASSERT (p4est_quadrant_is_ancestor (dom, q));
    P4EST_ASSERT (p4est_quadrant_child_id (q) == 0);
  }

  if (first_desc != NULL) {
    /* make sure that a quadrant at first_desc is represented in inlist */
    fd = *first_desc;
    while (fd.level > minlevel && p4est_quadrant_child_id (&fd) == 0) {
      p4est_quadrant_parent (&fd, &fd);
    }
    p4est_quadrant_sibling (&fd, &tempq, 0);
    si = p4est_find_lower_bound (inlist, &tempq, 0);
    if (si >= 0) {
      q = p4est_quadrant_array_index (inlist, si);
      p4est_nearest_common_ancestor (&tempq, q, &tempp);
      if (tempp.level < tempq.level - 1) {
        /* add tempq to inlist */
        sc_array_resize (inlist, inlist->elem_count + 1);
        memmove (sc_array_index (inlist, si + 1), sc_array_index (inlist, si),
                 (incount - si) * inlist->elem_size);
        q = p4est_quadrant_array_index (inlist, si);
        *q = tempq;
        q->p.user_int = 0;
        incount++;
      }
    }
    else {
      /* add tempq to inlist */
      q = p4est_quadrant_array_push_copy (inlist, &tempq);
      q->p.user_int = 0;
      incount++;
    }
  }
  else {
    p4est_quadrant_first_descendant (dom, &fd, minlevel);
  }

  if (last_desc != NULL) {
    /* make sure that a quadrant at last_desc is represented in inlist */
    tempq = *last_desc;
    while (tempq.level > minlevel &&
           p4est_quadrant_child_id (&tempq) == P4EST_CHILDREN - 1) {
      p4est_quadrant_parent (&tempq, &tempq);
    }
    p4est_quadrant_sibling (&tempq, &tempp, 0);
    si = p4est_find_higher_bound (inlist, last_desc, 0);
    if (si >= 0) {
      q = p4est_quadrant_array_index (inlist, si);
      p4est_nearest_common_ancestor (last_desc, q, &tempq);
      if (tempq.level < tempp.level - 1) {
        /* add tempp to inlist */
        sc_array_resize (inlist, inlist->elem_count + 1);
        if ((size_t) si < incount - 1) {
          memmove (sc_array_index (inlist, si + 2),
                   sc_array_index (inlist, si + 1),
                   (incount - (si + 1)) * inlist->elem_size);
        }
        q = p4est_quadrant_array_index (inlist, si + 1);
        *q = tempp;
        q->p.user_int = 0;
        incount++;
      }
    }
    else {
      /* add tempp to inlist */
      sc_array_resize (inlist, inlist->elem_count + 1);
      memmove (sc_array_index (inlist, 1), sc_array_index (inlist, 0),
               incount * inlist->elem_size);
      q = p4est_quadrant_array_index (inlist, 0);
      *q = tempp;
      q->p.user_int = 0;
      incount++;
    }
  }

  P4EST_ASSERT (sc_array_is_sorted (inlist, p4est_quadrant_compare));

  if (bound > 1) {
    /* initialize temporary storage */
    for (l = 0; l <= minlevel; ++l) {
      /* we don't need a hash table for minlevel, because all minlevel
       * quadrants will be created automatically by filling in gaps */
      hash[l] = NULL;
      memset (&outlist[l], -1, sizeof (sc_array_t));
    }
    for (; l < maxlevel; ++l) {
      hash[l] = sc_hash_new (p4est_quadrant_hash_fn, p4est_quadrant_equal_fn,
                             NULL, list_alloc);
      sc_array_init (&outlist[l], sizeof (p4est_quadrant_t *));
    }
    for (; l <= P4EST_MAXLEVEL; ++l) {
      /* we don't need a hash table for maxlevel because a quad only spawns
       * larger quads */
      hash[l] = NULL;
      memset (&outlist[l], -1, sizeof (sc_array_t));
    }
    outlist[maxlevel].elem_count = 0;

    /* walk through the input tree bottom-up */
    ph = 0;
    pid = -1;
    qalloc = p4est_quadrant_mempool_alloc (qpool);
    qalloc->p.user_int = 0;

    /* we don't need to run for minlevel + 1, because all of the quads that
     * would be created would be outside dom */
    for (l = maxlevel; l > minlevel + 1; l--) {
      ocount = outlist[l].elem_count;   /* ocount is not growing */
      olist = &outlist[l - 1];
      for (jz = 0; jz < incount + ocount; ++jz) {
        if (jz < incount) {
          q = p4est_quadrant_array_index (inlist, jz);
          if ((int) q->level != l || (q->p.user_int & duplicate)) {
            /* if a duplicate, don't run */
            continue;
          }
        }
        else {
          qpointer =
            (p4est_quadrant_t **) sc_array_index (&outlist[l], jz - incount);
          q = *qpointer;
          P4EST_ASSERT ((int) q->level == l);
        }
        P4EST_ASSERT (p4est_quadrant_is_ancestor (dom, q));
        P4EST_ASSERT (p4est_quadrant_child_id (q) == 0);

        p4est_quadrant_parent (q, &par);        /* get the parent */
        ph = P4EST_QUADRANT_LEN (par.level - 1);        /* twice its size */
        pid = p4est_quadrant_child_id (&par);   /* and position */
        p4est_quadrant_sibling (&par, &par, 0); /* now shift to 0 */

        for (sid = 0; sid < bound; sid++) {
          *qalloc = par;
          if (!sid) {
            qalloc->p.user_int = precluded;
            P4EST_ASSERT (p4est_quadrant_is_ancestor (dom, qalloc));
          }
          else if (sid <= P4EST_DIM) {
            /* include face neighbors */
            switch (sid - 1) {
            case 0:
              qalloc->x += ((pid & 1) ? ph : -ph);
              break;
            case 1:
              qalloc->y += ((pid & 2) ? ph : -ph);
              break;
#ifdef P4_TO_P8
            case 2:
              qalloc->z += ((pid & 4) ? ph : -ph);
              break;
#endif
            default:
              SC_ABORT_NOT_REACHED ();
            }
          }
#ifdef P4_TO_P8
          else if (sid < 7) {
            /* include edge neighbors */
            switch (sid - 4) {
            case 0:
              qalloc->y += ((pid & 2) ? ph : -ph);
              qalloc->z += ((pid & 4) ? ph : -ph);
              break;
            case 1:
              qalloc->x += ((pid & 1) ? ph : -ph);
              qalloc->z += ((pid & 4) ? ph : -ph);
              break;
            case 2:
              qalloc->x += ((pid & 1) ? ph : -ph);
              qalloc->y += ((pid & 2) ? ph : -ph);
              break;
            default:
              SC_ABORT_NOT_REACHED ();
            }
          }
#endif
          else {
            /* include corner neighbor */
            qalloc->x += ((pid & 1) ? ph : -ph);
            qalloc->y += ((pid & 2) ? ph : -ph);
#ifdef P4_TO_P8
            qalloc->z += ((pid & 4) ? ph : -ph);
#endif
          }

          P4EST_ASSERT (p4est_quadrant_is_extended (qalloc));
          P4EST_ASSERT (p4est_quadrant_child_id (qalloc) == 0);
          P4EST_ASSERT (!sid || qalloc->p.user_int == 0);
          P4EST_ASSERT (qalloc->level == l - 1);

          /* do not add quadrants outside of the domain */
          if (sid && !p4est_quadrant_is_ancestor (dom, qalloc)) {
            continue;
          }

          /* make sure that qalloc is not included more than once */
          inserted = sc_hash_insert_unique (hash[l - 1], qalloc, &vlookup);
          if (!inserted) {
            /* qalloc is already included in output list, this catches most */
            ++count_already_outlist;
            if (!sid) {
              /* we need to relay the fact that this octant is precluded */
              qlookup = (p4est_quadrant_t *) * vlookup;
              qlookup->p.user_int = precluded;
            }
            continue;
          }

          if (sid) {
            /* we do not need to search if we are adding the parent sibling: we
             * already know that the octant is precluded, and any other octant
             * we might find should already be marked duplicate */
            srindex = sc_array_bsearch (inlist, qalloc,
                                        p4est_quadrant_disjoint_parent);

            if (srindex != -1) {
              r = p4est_quadrant_array_index (inlist, srindex);

              if (r->level >= l - 1) {
                /* either qalloc duplicates r or is precluded by r: either way,
                 * we do not need to add qalloc to inlist in the final merge */
                qalloc->p.user_int = precluded;
                if (r->level > l - 1) {
                  ++count_ancestor_inlist;
                }
                else {
                  ++count_already_inlist;
                }
              }
              if (r->level <= l - 1) {
                /* either qalloc duplicates r, or an octant that can be traced to
                 * qalloc will duplicate r */
                r->p.user_int |= duplicate;
                if (r->level < l - 1) {
                  /* if qalloc precluded r, we can remove r before the final
                   * merge */
                  r->p.user_int |= precluded;
                }
              }
            }
          }

          qpointer = (p4est_quadrant_t **) sc_array_push (olist);
          *qpointer = qalloc;
          /* we need a new quadrant now, the old one is stored away */
          qalloc = p4est_quadrant_mempool_alloc (qpool);
          qalloc->p.user_int = 0;
        }
      }
    }
    sc_mempool_free (qpool, qalloc);

    /* remove unneeded octants */
    jz = 0;
    for (iz = 0; iz < incount; iz++) {
      q = p4est_quadrant_array_index (inlist, iz);
      if ((q->p.user_int & precluded) == 0) {
        if (jz != iz) {
          p = p4est_quadrant_array_index (inlist, jz++);
          *p = *q;
        }
        else {
          jz++;
        }
      }
    }
    sc_array_resize (inlist, jz);
    incount = jz;

    for (l = minlevel + 1; l < maxlevel; ++l) {
      /* print statistics and free hash tables */
#ifdef P4EST_ENABLE_DEBUG
      sc_hash_print_statistics (p4est_package_id, SC_LP_DEBUG, hash[l]);
#endif
      sc_hash_unlink_destroy (hash[l]);

      /* merge valid quadrants from outlist into inlist */
      ocount = outlist[l].elem_count;
      for (jz = 0; jz < ocount; ++jz) {
        /* go through output list */
        qpointer = (p4est_quadrant_t **) sc_array_index (&outlist[l], jz);
        qalloc = *qpointer;
        P4EST_ASSERT ((int) qalloc->level == l);
        P4EST_ASSERT (p4est_quadrant_is_ancestor (dom, qalloc));
        P4EST_ASSERT (p4est_quadrant_child_id (qalloc) == 0);
        /* copy temporary quadrant into inlist */
        if (first_desc != NULL && p4est_quadrant_compare (qalloc, &fd) < 0) {
          sc_mempool_free (qpool, qalloc);
          continue;
        }
        if (last_desc != NULL
            && p4est_quadrant_compare (qalloc, last_desc) > 0) {
          sc_mempool_free (qpool, qalloc);
          continue;
        }
        if (qalloc->p.user_int != precluded) {
          (void) p4est_quadrant_array_push_copy (inlist, qalloc);
        }
        sc_mempool_free (qpool, qalloc);
      }
      sc_array_reset (&outlist[l]);
    }
    P4EST_ASSERT (quadrant_pool_size == qpool->elem_count);
    sc_mempool_truncate (list_alloc);

    /* sort inlist */
    if (inlist->elem_count > incount) {
      sc_array_sort (inlist, p4est_quadrant_compare);
    }
  }

  /* step through inlist and fill in the gaps in out */
  /* note: we add to the end of out */
  ocount = out->elem_count;

  incount = inlist->elem_count;

  tempq = fd;
  if (first_desc == NULL) {
    pid = 0;
    jz = 0;
  }
  else {
    /* find the first quadrant after tempq */
    pid = p4est_quadrant_child_id (&tempq);
    si = p4est_find_lower_bound (inlist, &tempq, 0);
    if (si >= 0) {
      jz = si;
    }
    else {
      jz = incount;
    }
  }
  if (jz < incount) {
    q = p4est_quadrant_array_index (inlist, jz);
    P4EST_ASSERT (p4est_quadrant_child_id (q) == 0);
  }
  else if (last_desc != NULL) {
    p4est_quadrant_last_descendant (dom, &ld, P4EST_QMAXLEVEL);
    if (p4est_quadrant_is_equal (&ld, last_desc)) {
      q = NULL;
    }
    else {
      p4est_quadrant_successor (last_desc, &ld);
      P4EST_ASSERT (p4est_quadrant_is_ancestor (dom, &ld));
      q = &ld;
#ifdef P4EST_ENABLE_DEBUG
      P4EST_QUADRANT_INIT (&ld_old);
      p4est_quadrant_linear_id_ext128 (last_desc, P4EST_QMAXLEVEL, &lid);
      p4est_lid_add_inplace (&lid, &one);
      p4est_quadrant_set_morton_ext128 (&ld_old, P4EST_QMAXLEVEL, &lid);
      P4EST_ASSERT (p4est_quadrant_is_ancestor (dom, &ld_old));
      P4EST_ASSERT (p4est_quadrant_is_equal (&ld, &ld_old));
#endif
    }
  }
  else {
    q = NULL;
  }
  for (;;) {
    /* while tempq comes before q */
    while (q == NULL || (!p4est_quadrant_is_equal (&tempq, q) &&
                         !p4est_quadrant_is_ancestor (&tempq, q))) {
      P4EST_ASSERT (q == NULL || p4est_quadrant_compare (&tempq, q) < 0);

      /* stop once we're past last_desc */
      if (last_desc != NULL && p4est_quadrant_compare (&tempq, last_desc) > 0) {
        break;
      }
      /* add tempq to out */
      (void) p4est_quadrant_array_push_copy (out, &tempq);

      /* if tempq is a last sibling, go up a level */
      while (tempq.level >= minlevel && pid == P4EST_CHILDREN - 1) {
        p4est_quadrant_parent (&tempq, &tempp);
        tempq = tempp;
        pid = p4est_quadrant_child_id (&tempq);
      }

      /* if we've finished with all minlevel and smaller quadrants, we've
       * filled dom */
      if (tempq.level < minlevel) {
        break;
      }

      /* get the next sibling */
      p4est_quadrant_sibling (&tempq, &tempq, ++pid);
    }

    /* if we've finished with all minlevel and smaller quadrants, we've
     * filled dom. also stop if we're past last_desc */
    if (tempq.level < minlevel ||
        (last_desc != NULL &&
         p4est_quadrant_compare (&tempq, last_desc) > 0)) {
      break;
    }

    P4EST_ASSERT (q != NULL);
    P4EST_ASSERT (p4est_quadrant_is_equal (&tempq, q) ||
                  p4est_quadrant_is_ancestor (&tempq, q));

    if (q->x == tempq.x && q->y == tempq.y &&
#ifdef P4_TO_P8
        q->z == tempq.z &&
#endif
        1) {
      /* if q is the first descendant of tempq, set tempq = q and get the next
       * q */
      if (q->level > tempq.level) {
        pid = 0;
      }
      tempq.level = q->level;
      jz++;
      if (jz < incount) {
        q = p4est_quadrant_array_index (inlist, jz);
        P4EST_ASSERT (p4est_quadrant_child_id (q) == 0);
      }
      else if (last_desc != NULL) {
        p4est_quadrant_last_descendant (dom, &ld, P4EST_QMAXLEVEL);
        if (p4est_quadrant_is_equal (&ld, last_desc)) {
          q = NULL;
        }
        else {
          p4est_quadrant_successor (last_desc, &ld);
          P4EST_ASSERT (p4est_quadrant_is_ancestor (dom, &ld));
          q = &ld;
#ifdef P4EST_ENABLE_DEBUG
          P4EST_QUADRANT_INIT (&ld_old);
          p4est_quadrant_linear_id_ext128 (last_desc, P4EST_QMAXLEVEL, &lid);
          p4est_lid_add_inplace (&lid, &one);
          p4est_quadrant_set_morton_ext128 (&ld_old, P4EST_QMAXLEVEL, &lid);
          P4EST_ASSERT (p4est_quadrant_is_ancestor (dom, &ld_old));
          P4EST_ASSERT (p4est_quadrant_is_equal (&ld, &ld_old));
#endif
        }
      }
      else {
        q = NULL;
      }
    }
    else {
      /* get the largest first descendant of tempq that comes before
       * q */
      p4est_quadrant_first_descendant (&tempq, &tempp, P4EST_QMAXLEVEL);
      p4est_nearest_common_ancestor (&tempp, q, &tempq);
      tempq.level++;
      pid = 0;
      P4EST_ASSERT (p4est_quadrant_is_valid (&tempq));
      P4EST_ASSERT (p4est_quadrant_compare (&tempq, q) < 0);
      P4EST_ASSERT (!p4est_quadrant_is_ancestor (&tempq, q));
    }
  }

#ifdef P4EST_ENABLE_DEBUG
  sc_array_init_view (&outview, out, ocount, out->elem_count - ocount);
  P4EST_ASSERT (sc_array_is_sorted (&outview, p4est_quadrant_compare));
  P4EST_ASSERT (outview.elem_count > 1);

  for (jz = 0; jz < outview.elem_count - 1; jz++) {
    q = p4est_quadrant_array_index (&outview, jz);
    r = p4est_quadrant_array_index (&outview, jz + 1);
    P4EST_ASSERT (p4est_quadrant_is_next (q, r));
  }
#endif

  if (count_in != NULL) {
    *count_in += count_already_inlist;
  }
  if (count_out != NULL) {
    *count_out += count_already_outlist;
  }
  if (count_an != NULL) {
    *count_an += count_ancestor_inlist;
  }

}

static void
p4est_balance_replace_recursive (p4est_t * p4est, p4est_topidx_t nt,
                                 sc_array_t * array, size_t start, size_t end,
                                 p4est_quadrant_t * parent,
                                 p4est_init_t init_fn,
                                 p4est_replace_t replace_fn)
{
  p4est_quadrant_t    fam[P4EST_CHILDREN];
  p4est_quadrant_t   *famp[P4EST_CHILDREN];
  sc_array_t          view;
  size_t              jz;
  size_t              iz[P4EST_CHILDREN + 1];

  if (end - start == P4EST_CHILDREN) {
    for (jz = 0; jz < P4EST_CHILDREN; jz++) {
      famp[jz] = p4est_quadrant_array_index (array, start + jz);
    }
    P4EST_ASSERT (p4est_quadrant_is_familypv (famp));
    replace_fn (p4est, nt, 1, &parent, P4EST_CHILDREN, famp);
    p4est_quadrant_free_data (p4est, parent);
    return;
  }
  sc_array_init_view (&view, array, start, end - start);
  p4est_split_array (&view, parent->level, iz);

  for (jz = 0; jz < P4EST_CHILDREN; jz++) {
    if (iz[jz] + 1 == iz[jz + 1]) {
      famp[jz] = p4est_quadrant_array_index (array, start + iz[jz]);
      P4EST_ASSERT (p4est_quadrant_is_parent (parent, famp[jz]));
      P4EST_ASSERT (p4est_quadrant_child_id (famp[jz]) == (int) jz);
    }
    else {
      fam[jz] = *parent;
      famp[jz] = &fam[jz];
      famp[jz]->level++;
      p4est_quadrant_sibling (famp[jz], famp[jz], (int) jz);
      p4est_quadrant_init_data (p4est, nt, famp[jz], init_fn);
    }
  }
  replace_fn (p4est, nt, 1, &parent, P4EST_CHILDREN, famp);
  p4est_quadrant_free_data (p4est, parent);

  for (jz = 0; jz < P4EST_CHILDREN; jz++) {
    if (famp[jz] == &fam[jz]) {
      p4est_balance_replace_recursive (p4est, nt, array, start + iz[jz],
                                       start + iz[jz + 1], famp[jz],
                                       init_fn, replace_fn);
    }
  }
}

static void
p4est_complete_or_balance (p4est_t * p4est, p4est_topidx_t which_tree,
                           p4est_init_t init_fn, p4est_replace_t replace_fn,
                           int btype)
{
  p4est_tree_t       *tree;
  sc_array_t         *tquadrants;
  int                 bound;
  int8_t              maxlevel;
  sc_mempool_t       *qpool;
#ifdef P4EST_ENABLE_DEBUG
  size_t              data_pool_size;
#endif
  size_t              tcount;
  size_t              count_already_inlist, count_already_outlist;
  size_t              count_ancestor_inlist;
  p4est_quadrant_t   *q, *p;
  sc_mempool_t       *list_alloc;
  sc_array_t         *inlist, *outlist;
  size_t              iz, jz, jzstart = 0, jzend, ocount;
  p4est_quadrant_t    tempq, root;

  P4EST_ASSERT (which_tree >= p4est->first_local_tree);
  P4EST_ASSERT (which_tree <= p4est->last_local_tree);
  tree = p4est_tree_array_index (p4est->trees, which_tree);
  tquadrants = &(tree->quadrants);

  P4EST_ASSERT (0 <= btype && btype <= P4EST_DIM);
  P4EST_ASSERT (sc_array_is_sorted (tquadrants, p4est_quadrant_compare));

  P4EST_QUADRANT_INIT (&tempq);
  P4EST_QUADRANT_INIT (&root);

  switch (btype) {
  case 0:
    bound = 1;
    break;
  case 1:
    bound = P4EST_DIM + 1;
    break;
  case P4EST_DIM:
    bound = (1 << P4EST_DIM);
    break;
#ifdef P4_TO_P8
  case 2:
    bound = 7;
    break;
#endif
  default:
    SC_ABORT_NOT_REACHED ();
  }

  qpool = p4est->quadrant_pool;

#ifdef P4EST_ENABLE_DEBUG
  data_pool_size = 0;
  if (p4est->user_data_pool != NULL) {
    data_pool_size = p4est->user_data_pool->elem_count;
  }
#endif

  tcount = tquadrants->elem_count;
  /* if tree is empty, there is nothing to do */
  if (!tcount) {
    return;
  }

  /* initialize some counters */
  count_already_inlist = count_already_outlist = 0;
  count_ancestor_inlist = 0;

  /* get containing quadrant */
  P4EST_QUADRANT_INIT (&root);
  p4est_nearest_common_ancestor (&tree->first_desc, &tree->last_desc, &root);

  if (tcount == 1) {
    p = p4est_quadrant_array_index (tquadrants, 0);
    if (p4est_quadrant_is_equal (p, &root)) {
      /* nothing to be done, clean up and exit */
      return;
    }
  }

  /* initialize temporary storage */
  list_alloc = sc_mempool_new (sizeof (sc_link_t));

  inlist = sc_array_new (sizeof (p4est_quadrant_t));
  outlist = sc_array_new (sizeof (p4est_quadrant_t));

  /* get the reduced representation of the tree */
  q = p4est_quadrant_array_push (inlist);
  p = p4est_quadrant_array_index (tquadrants, 0);
  p4est_quadrant_sibling (p, q, 0);
  for (iz = 1; iz < tcount; iz++) {
    p = p4est_quadrant_array_index (tquadrants, iz);
    P4EST_ASSERT (p4est_quadrant_is_ancestor (&root, p));
    p4est_nearest_common_ancestor (p, q, &tempq);
    if (tempq.level >= SC_MIN (q->level, p->level) - 1) {
      if (p->level > q->level) {
        p4est_quadrant_sibling (p, q, 0);
      }
      continue;
    }
    q = p4est_quadrant_array_push (inlist);
    p4est_quadrant_sibling (p, q, 0);
  }

  /* balance */
  p4est_complete_or_balance_kernel (inlist, &root, bound, qpool,
                                    list_alloc, outlist,
                                    &(tree->first_desc),
                                    &(tree->last_desc),
                                    &count_already_inlist,
                                    &count_already_outlist,
                                    &count_ancestor_inlist);

  ocount = outlist->elem_count;

  iz = 0;                       /* tquadrants */
  jz = 0;                       /* outlist */
  maxlevel = tree->maxlevel;

  /* initialize quadrants in outlist */
  while (iz < tcount && jz < ocount) {
    q = p4est_quadrant_array_index (tquadrants, iz);
    p = p4est_quadrant_array_index (outlist, jz);

    /* watch out for gaps in tquadrants */
    while (p4est_quadrant_compare (p, q) < 0) {
      P4EST_ASSERT (!p4est_quadrant_is_ancestor (p, q));
      maxlevel = SC_MAX (maxlevel, p->level);
      ++tree->quadrants_per_level[p->level];
      p4est_quadrant_init_data (p4est, which_tree, p, init_fn);
      jz++;
      P4EST_ASSERT (jz < ocount);
      p = p4est_quadrant_array_index (outlist, jz);
    }

    /* watchout out for tquadrants that have been split */
    if (q->level < p->level) {
      P4EST_ASSERT (p4est_quadrant_is_ancestor (q, p));
      /* reset q */
      --tree->quadrants_per_level[q->level];
      if (replace_fn == NULL) {
        p4est_quadrant_free_data (p4est, q);
      }
      else {
        tempq = *q;
        jzstart = jz;
      }
      while (jz < ocount && p4est_quadrant_is_ancestor (q, p)) {
        maxlevel = SC_MAX (maxlevel, p->level);
        ++tree->quadrants_per_level[p->level];
        p4est_quadrant_init_data (p4est, which_tree, p, init_fn);
        if (++jz < ocount) {
          p = p4est_quadrant_array_index (outlist, jz);
        }
      }
      if (replace_fn != NULL) {
        jzend = jz;
        p4est_balance_replace_recursive (p4est, which_tree,
                                         outlist, jzstart, jzend, &tempq,
                                         init_fn, replace_fn);
      }
      iz++;
    }
    else {
      P4EST_ASSERT (p4est_quadrant_is_equal (q, p));
      p->p.user_data = q->p.user_data;
      iz++;
      jz++;
    }
  }

  P4EST_ASSERT (iz == tcount);

  /* initialize new quadrants after last tquadrant */
  for (; jz < ocount; jz++) {
    p = p4est_quadrant_array_index (outlist, jz);
    maxlevel = SC_MAX (maxlevel, p->level);
    ++tree->quadrants_per_level[p->level];
    p4est_quadrant_init_data (p4est, which_tree, p, init_fn);
  }

  /* resize tquadrants and copy */
  sc_array_resize (tquadrants, ocount);
  memcpy (tquadrants->array, outlist->array, outlist->elem_size * ocount);
  tree->maxlevel = maxlevel;

  /* sanity check */
  if (p4est->user_data_pool != NULL) {
    P4EST_ASSERT (data_pool_size + (ocount - tcount) ==
                  p4est->user_data_pool->elem_count);
  }

  P4EST_VERBOSEF
    ("Tree %lld inlist %llu outlist %llu ancestor %llu insert %llu\n",
     (long long) which_tree, (unsigned long long) count_already_inlist,
     (unsigned long long) count_already_outlist,
     (unsigned long long) count_ancestor_inlist,
     (unsigned long long) (ocount - tcount));

  sc_array_destroy (inlist);
  sc_array_destroy (outlist);
  sc_mempool_destroy (list_alloc);

  if (p4est->inspect) {
    if (!p4est->inspect->use_B) {
      p4est->inspect->balance_A_count_in += count_already_inlist;
      p4est->inspect->balance_A_count_in += count_ancestor_inlist;
      p4est->inspect->balance_A_count_out += count_already_outlist;
    }
    else {
      p4est->inspect->balance_B_count_in += count_already_inlist;
      p4est->inspect->balance_B_count_in += count_ancestor_inlist;
      p4est->inspect->balance_B_count_out += count_already_outlist;
    }
  }
}

void
p4est_balance_border (p4est_t * p4est, p4est_connect_type_t btype,
                      p4est_topidx_t which_tree, p4est_init_t init_fn,
                      p4est_replace_t replace_fn, sc_array_t * borders)
{
  size_t              iz, jz, kz;
  size_t              incount;
  size_t              count_already_inlist, count_already_outlist;
  size_t              count_ancestor_inlist;
  p4est_tree_t       *tree;
  p4est_quadrant_t   *q, *p, *r;
  p4est_quadrant_t    tempq, tempp;
  sc_array_t          qview;
  sc_array_t         *inlist, *flist, *tquadrants;
  sc_array_t          tqview;
  size_t              tqoffset, fcount;
  p4est_topidx_t      first_tree = p4est->first_local_tree;
  size_t              num_added, num_this_added;
  int                 bound;
  ssize_t             tqindex;
  size_t              tqorig;
  sc_mempool_t       *list_alloc, *qpool;
  /* get this tree's border */
  sc_array_t         *qarray = (sc_array_t *) sc_array_index (borders,
                                                              which_tree -
                                                              first_tree);
  size_t              qcount = qarray->elem_count;

  if (!qcount) {
    /* nothing to be done */
    return;
  }

  P4EST_QUADRANT_INIT (&tempq);
  P4EST_QUADRANT_INIT (&tempp);

  /* set up balance machinery */

  if (btype == P4EST_CONNECT_FULL) {
    bound = (1 << P4EST_DIM);
  }
#ifdef P4_TO_P8
  else if (btype == P8EST_CONNECT_EDGE) {
    bound = (1 << P4EST_DIM) - 1;
  }
#endif
  else {
    bound = P4EST_DIM + 1;
  }

  P4EST_ASSERT (which_tree >= p4est->first_local_tree);
  P4EST_ASSERT (which_tree <= p4est->last_local_tree);

  tree = p4est_tree_array_index (p4est->trees, which_tree);
  tquadrants = &(tree->quadrants);
  tqorig = tquadrants->elem_count;
  tqoffset = 0;
  sc_array_init_view (&tqview, tquadrants, tqoffset,
                      tquadrants->elem_count - tqoffset);

  qpool = p4est->quadrant_pool;

  count_already_inlist = count_already_outlist = 0;
  count_ancestor_inlist = 0;
  num_added = 0;

  /* initialize temporary storage */
  list_alloc = sc_mempool_new (sizeof (sc_link_t));

  inlist = sc_array_new (sizeof (p4est_quadrant_t));
  flist = sc_array_new (sizeof (p4est_quadrant_t));

  /* sort the border and remove duplicates */
  sc_array_sort (qarray, p4est_quadrant_compare);
  jz = 1;                       /* number included */
  kz = 0;                       /* number skipped */
  p = p4est_quadrant_array_index (qarray, 0);
  P4EST_ASSERT (p4est_quadrant_is_valid (p));
  for (iz = 1; iz < qcount; iz++) {
    q = p4est_quadrant_array_index (qarray, iz);
    P4EST_ASSERT (p4est_quadrant_is_extended (q));
    if (!p4est_quadrant_is_equal (q, p)) {
      p++;
      jz++;
      if (kz) {
        *p = *q;
      }
    }
    else {
      kz++;
    }
  }
  P4EST_ASSERT (kz + jz == qcount);
  sc_array_resize (qarray, jz);
  qcount = jz;

  /* step through border */
  for (iz = 0; iz < qcount; iz++) {
    p = p4est_quadrant_array_index (qarray, iz);

    if (p4est_quadrant_compare (p, &(tree->first_desc)) < 0 &&
        !p4est_quadrant_is_ancestor (p, &(tree->first_desc))) {
      continue;
    }
    if (p4est_quadrant_compare (p, &(tree->last_desc)) > 0) {
      continue;
    }

    P4EST_ASSERT (p4est_quadrant_is_valid (p));

    /* get a view of all of the quads that are descended from this quad */
    jz = iz + 1;
    kz = jz;

    if (kz < qcount) {
      q = p4est_quadrant_array_index (qarray, kz);
    }

    while (kz < qcount && p4est_quadrant_is_ancestor (p, q)) {
      kz++;

      P4EST_ASSERT (p4est_quadrant_child_id (q) == 0);

      if (kz < qcount) {
        q = p4est_quadrant_array_index (qarray, kz);
      }
    }

    incount = kz - jz;
    if (!incount) {
      continue;
    }

    /* find p in tquadrants */
    tqindex = sc_array_bsearch (&tqview, p, p4est_quadrant_compare);

    P4EST_ASSERT (tqindex >= 0);

    /* copy everything before p into flist */
    if (tqindex) {
      fcount = flist->elem_count;
      sc_array_resize (flist, fcount + tqindex);
      memcpy (sc_array_index (flist, fcount),
              tqview.array, tqindex * sizeof (p4est_quadrant_t));
    }

    /* update the view of tquadrants to be everything past p */
    tqindex += tqoffset;        /* tqindex is the index of p in tquadrants */
    tqoffset = tqindex + 1;
    sc_array_init_view (&tqview, tquadrants, tqoffset, tqorig - tqoffset);

    /* first, remove p */
    q = p4est_quadrant_array_index (tquadrants, tqindex);
    P4EST_ASSERT (p4est_quadrant_is_equal (q, p));
    /* reset the data, decrement level count */
    if (replace_fn == NULL) {
      p4est_quadrant_free_data (p4est, q);
    }
    else {
      tempp = *q;
    }
    --tree->quadrants_per_level[q->level];

    /* get all of the quadrants that descend from p into inlist */
    sc_array_init_view (&qview, qarray, jz, incount);
    sc_array_resize (inlist, 1);
    q = p4est_quadrant_array_index (inlist, 0);
    r = p4est_quadrant_array_index (&qview, 0);
    P4EST_ASSERT (p4est_quadrant_child_id (r) == 0);
    *q = *r;
    for (jz = 1; jz < incount; jz++) {
      r = p4est_quadrant_array_index (&qview, jz);
      P4EST_ASSERT (p4est_quadrant_child_id (r) == 0);
      p4est_nearest_common_ancestor (r, q, &tempq);
      if (tempq.level >= SC_MIN (r->level, q->level) - 1) {
        if (r->level > q->level) {
          *q = *r;
        }
        continue;
      }
      q = p4est_quadrant_array_push_copy (inlist, r);
    }

    fcount = flist->elem_count;

    /* balance them within the containing quad */
    p4est_complete_or_balance_kernel (inlist, p, bound, qpool, list_alloc,
                                      flist, NULL, NULL,
                                      &count_already_inlist,
                                      &count_already_outlist,
                                      &count_ancestor_inlist);

    /* count the amount we've added (-1 because we subtract p) */
    num_this_added = flist->elem_count - 1 - fcount;
    num_added += num_this_added;

    /* initialize */
    for (jz = fcount; jz < flist->elem_count; jz++) {
      q = p4est_quadrant_array_index (flist, jz);
      P4EST_ASSERT (p4est_quadrant_is_ancestor (p, q));
      ++tree->quadrants_per_level[q->level];
      tree->maxlevel = (int8_t) SC_MAX (tree->maxlevel, q->level);
      p4est_quadrant_init_data (p4est, which_tree, q, init_fn);
    }
    if (replace_fn != NULL) {
      p4est_balance_replace_recursive (p4est, which_tree,
                                       flist, fcount, flist->elem_count,
                                       &tempp, init_fn, replace_fn);
    }

    /* skip over the quadrants that we just operated on */
    iz = kz - 1;
  }

  /* copy the remaining tquadrants to flist */
  if (tqoffset < tqorig) {
    fcount = flist->elem_count;
    sc_array_resize (flist, fcount + tqorig - tqoffset);
    memcpy (sc_array_index (flist, fcount),
            tqview.array, (tqorig - tqoffset) * sizeof (p4est_quadrant_t));
  }

  /* copy flist into tquadrants */
  sc_array_resize (tquadrants, flist->elem_count);
  memcpy (tquadrants->array, flist->array,
          flist->elem_count * flist->elem_size);

  sc_mempool_destroy (list_alloc);
  P4EST_ASSERT (tqorig + num_added == tquadrants->elem_count);

  /* print more statistics */
  P4EST_VERBOSEF
    ("Tree border %lld inlist %llu outlist %llu ancestor %llu insert %llu\n",
     (long long) which_tree, (unsigned long long) count_already_inlist,
     (unsigned long long) count_already_outlist,
     (unsigned long long) count_ancestor_inlist,
     (unsigned long long) num_added);

  sc_array_destroy (inlist);
  sc_array_destroy (flist);

  P4EST_ASSERT (p4est_tree_is_complete (tree));

  if (p4est->inspect) {
    p4est->inspect->balance_B_count_in += count_already_inlist;
    p4est->inspect->balance_B_count_in += count_ancestor_inlist;
    p4est->inspect->balance_B_count_out += count_already_outlist;
  }
}

void
p4est_complete_subtree (p4est_t * p4est,
                        p4est_topidx_t which_tree, p4est_init_t init_fn)
{
  p4est_complete_or_balance (p4est, which_tree, init_fn, NULL, 0);
}

void
p4est_balance_subtree (p4est_t * p4est, p4est_connect_type_t btype,
                       p4est_topidx_t which_tree, p4est_init_t init_fn)
{
  p4est_complete_or_balance (p4est, which_tree, init_fn, NULL,
                             p4est_connect_type_int (btype));
}

void
p4est_balance_subtree_ext (p4est_t * p4est, p4est_connect_type_t btype,
                           p4est_topidx_t which_tree, p4est_init_t init_fn,
                           p4est_replace_t replace_fn)
{
  p4est_complete_or_balance (p4est, which_tree, init_fn, replace_fn,
                             p4est_connect_type_int (btype));
}

size_t
p4est_linearize_tree (p4est_t * p4est, p4est_tree_t * tree)
{
#ifdef P4EST_ENABLE_DEBUG
  size_t              data_pool_size;
#endif
  size_t              incount, removed;
  size_t              current, rest;
  p4est_locidx_t      num_quadrants;
  int                 i, maxlevel;
  p4est_quadrant_t   *q1, *q2;
  sc_array_t         *tquadrants = &tree->quadrants;

  P4EST_ASSERT (sc_array_is_sorted (tquadrants, p4est_quadrant_compare));

  incount = tquadrants->elem_count;
  if (incount <= 1) {
    return 0;
  }
#ifdef P4EST_ENABLE_DEBUG
  data_pool_size = 0;
  if (p4est->user_data_pool != NULL) {
    data_pool_size = p4est->user_data_pool->elem_count;
  }
#endif
  removed = 0;

  /* run through the array and remove ancestors */
  current = 0;
  rest = current + 1;
  q1 = p4est_quadrant_array_index (tquadrants, current);
  while (rest < incount) {
    q2 = p4est_quadrant_array_index (tquadrants, rest);
    if (p4est_quadrant_is_equal (q1, q2) ||
        p4est_quadrant_is_ancestor (q1, q2)) {
      --tree->quadrants_per_level[q1->level];
      p4est_quadrant_free_data (p4est, q1);
      *q1 = *q2;
      ++removed;
      ++rest;
    }
    else {
      ++current;
      if (current < rest) {
        q1 = p4est_quadrant_array_index (tquadrants, current);
        *q1 = *q2;
      }
      else {
        q1 = q2;
      }
      ++rest;
    }
  }

  /* resize array */
  sc_array_resize (tquadrants, current + 1);

  /* update level counters */
  maxlevel = 0;
  num_quadrants = 0;
  for (i = 0; i <= P4EST_QMAXLEVEL; ++i) {
    P4EST_ASSERT (tree->quadrants_per_level[i] >= 0);
    num_quadrants += tree->quadrants_per_level[i];      /* same type */
    if (tree->quadrants_per_level[i] > 0) {
      maxlevel = i;
    }
  }
  tree->maxlevel = (int8_t) maxlevel;

  /* sanity checks */
  P4EST_ASSERT (num_quadrants == (p4est_locidx_t) tquadrants->elem_count);
  P4EST_ASSERT (tquadrants->elem_count == incount - removed);
  if (p4est->user_data_pool != NULL) {
    P4EST_ASSERT (data_pool_size - removed ==
                  p4est->user_data_pool->elem_count);
  }
  P4EST_ASSERT (p4est_tree_is_sorted (tree));
  P4EST_ASSERT (p4est_tree_is_linear (tree));

  return removed;
}

p4est_locidx_t
p4est_partition_correction (p4est_gloidx_t * partition,
                            int num_procs, int rank,
                            p4est_gloidx_t min_quadrant_id,
                            p4est_gloidx_t max_quadrant_id)
{
  int                 i;
  int                 rank_with_max_quads = rank;
  p4est_gloidx_t      h;
  p4est_gloidx_t      max_num_quadrants =
    SC_MIN (max_quadrant_id, partition[rank + 1] - 1) - partition[rank] + 1;

  /* no correction if num quadrants not sufficient for family */
  if (max_quadrant_id - min_quadrant_id + 1 != P4EST_CHILDREN) {
    return 0;
  }

  /* decreasing search for process with highest amount of quadrants */
  i = rank_with_max_quads - 1;
  while (min_quadrant_id < partition[i + 1]) {
    h = partition[i + 1] - SC_MAX (min_quadrant_id, partition[i]);
    if (max_num_quadrants <= h) {
      max_num_quadrants = h;
      rank_with_max_quads = i;
    }
    i--;
  }

  /* increasing search for process with highest amount of quadrants */
  i = rank_with_max_quads + 1;
  while (partition[i] <= max_quadrant_id) {
    h = SC_MIN (max_quadrant_id, partition[i + 1] - 1) - partition[i] + 1;
    if (max_num_quadrants < h) {
      max_num_quadrants = h;
      rank_with_max_quads = i;
    }
    i++;
  }

  /* compute correction */
  if (rank_with_max_quads < rank) {
    return (p4est_locidx_t) (partition[rank] - max_quadrant_id - 1);
  }
  else {
    return (p4est_locidx_t) (partition[rank] - min_quadrant_id);
  }
}

int
p4est_next_nonempty_process (int rank, int num_procs,
                             p4est_locidx_t * num_quadrants_in_proc)
{
  if (rank >= num_procs) {      /* if `rank` is too high */
    /* return process id beyond scope */
    return num_procs;
  }

  /* search for next non empty process */
  while (rank < num_procs && num_quadrants_in_proc[rank] == 0) {
    rank++;
  }

  /* return non empty process id or last process id respectively */
  return rank;
}

p4est_gloidx_t
p4est_partition_given (p4est_t * p4est,
                       const p4est_locidx_t * new_num_quadrants_in_proc)
{
  const int           num_procs = p4est->mpisize;
  const int           rank = p4est->mpirank;
  const p4est_topidx_t first_local_tree = p4est->first_local_tree;
  const p4est_topidx_t last_local_tree = p4est->last_local_tree;
  const size_t        data_size = p4est->data_size;
  const size_t        quad_plus_data_size = sizeof (p4est_quadrant_t)
    + data_size;
  sc_array_t         *trees = p4est->trees;

  /* *INDENT-OFF* horrible indent bug */
  const p4est_topidx_t num_send_trees =
    p4est->global_first_position[rank + 1].p.which_tree - /* same type */
    p4est->global_first_position[rank].p.which_tree + 1;
  /* *INDENT-ON* */

  int                 i;
  int                 from_proc, to_proc;
  int                 num_proc_recv_from, num_proc_send_to;
  char               *user_data_send_buf;
  char               *user_data_recv_buf;
  char              **recv_buf, **send_buf;
  size_t              recv_size, send_size, zz, zoffset;
  p4est_topidx_t      it;
  p4est_topidx_t      which_tree;
  p4est_topidx_t      first_tree, last_tree;
  p4est_topidx_t      num_recv_trees;
  p4est_topidx_t      new_first_local_tree, new_last_local_tree;
  p4est_topidx_t      first_from_tree, last_from_tree, from_tree;
  p4est_locidx_t      il;
  p4est_locidx_t      num_copy;
  p4est_locidx_t      num_quadrants;
  p4est_locidx_t      new_local_num_quadrants;
  p4est_locidx_t     *num_recv_from, *num_send_to;
  p4est_locidx_t     *new_local_tree_elem_count;
  p4est_locidx_t     *new_local_tree_elem_count_before;
  p4est_locidx_t     *num_per_tree_local;
  p4est_locidx_t     *num_per_tree_send_buf;
  p4est_locidx_t     *num_per_tree_recv_buf;
  p4est_gloidx_t     *begin_send_to;
  p4est_gloidx_t      tree_from_begin, tree_from_end, num_copy_global;
  p4est_gloidx_t      from_begin, from_end, lower_bound,
    from_begin_global_quad, from_end_global_quad, to_begin_global_quad,
    to_end_global_quad;
  p4est_gloidx_t      to_begin, to_end;
  p4est_gloidx_t      my_base, my_begin, my_end;
  p4est_gloidx_t     *global_last_quad_index;
  p4est_gloidx_t     *new_global_last_quad_index;
  p4est_gloidx_t     *local_tree_last_quad_index;
  p4est_gloidx_t      diff64, total_quadrants_shipped;
  sc_array_t         *quadrants;
  p4est_quadrant_t   *quad_send_buf;
  p4est_quadrant_t   *quad_recv_buf;
  p4est_quadrant_t   *quad;
  p4est_tree_t       *tree;
#ifdef P4EST_ENABLE_MPI
  int                 sk;
  int                 mpiret;
  MPI_Comm            comm = p4est->mpicomm;
  MPI_Request        *recv_request, *send_request;
#endif
#ifdef P4EST_ENABLE_DEBUG
  int                 send_to_empty;
  unsigned            crc;
  p4est_gloidx_t      my_begin_comp, my_end_comp;
  p4est_gloidx_t      total_requested_quadrants = 0;
#endif

  P4EST_GLOBAL_INFOF
    ("Into " P4EST_STRING "_partition_given with %lld total quadrants\n",
     (long long) p4est->global_num_quadrants);

#ifdef P4EST_ENABLE_DEBUG
  /* Save a checksum of the original forest */
  crc = p4est_checksum (p4est);
#endif

  global_last_quad_index = P4EST_ALLOC (p4est_gloidx_t, num_procs);

  new_global_last_quad_index = P4EST_ALLOC (p4est_gloidx_t, num_procs);
  new_global_last_quad_index[0] = new_num_quadrants_in_proc[0] - 1;

  total_quadrants_shipped = 0;

  for (i = 0; i < num_procs; ++i) {
    /* Check for a valid requested partition and create last_quad_index */
    global_last_quad_index[i] = p4est->global_first_quadrant[i + 1] - 1;
#ifdef P4EST_ENABLE_DEBUG
    total_requested_quadrants += new_num_quadrants_in_proc[i];
    P4EST_ASSERT (new_num_quadrants_in_proc[i] >= 0);
#endif

    /* Print some diagnostics */
    if (rank == 0) {
      P4EST_GLOBAL_LDEBUGF ("partition global_last_quad_index[%d] = %lld\n",
                            i, (long long) global_last_quad_index[i]);
    }

    if (i >= 1) {
      /* Calculate the global_last_quad_index for the repartitioned forest */
      new_global_last_quad_index[i] = new_num_quadrants_in_proc[i] +
        new_global_last_quad_index[i - 1];

      /* Calculate the global number of shipped (= received) quadrants */
      diff64 =
        global_last_quad_index[i - 1] - new_global_last_quad_index[i - 1];
      if (diff64 >= 0) {
        total_quadrants_shipped +=
          SC_MIN (diff64, new_num_quadrants_in_proc[i]);
      }
      else {
        total_quadrants_shipped +=
          SC_MIN (-diff64, new_num_quadrants_in_proc[i - 1]);
      }
    }

    /* Print some diagnostics */
    if (rank == 0) {
      P4EST_GLOBAL_LDEBUGF
        ("partition new_global_last_quad_index[%d] = %lld\n",
         i, (long long) new_global_last_quad_index[i]);
    }
  }
  P4EST_ASSERT (total_requested_quadrants == p4est->global_num_quadrants);
  P4EST_ASSERT (global_last_quad_index[num_procs - 1] ==
                new_global_last_quad_index[num_procs - 1]);
  P4EST_ASSERT (0 <= total_quadrants_shipped &&
                total_quadrants_shipped <= p4est->global_num_quadrants);

  /* Calculate the local index of the end of each tree */
  local_tree_last_quad_index =
    P4EST_ALLOC_ZERO (p4est_gloidx_t, trees->elem_count);
  if (first_local_tree >= 0) {
    tree = p4est_tree_array_index (trees, first_local_tree);
    local_tree_last_quad_index[first_local_tree]
      = tree->quadrants.elem_count - 1;
  }
  else {
    /* empty processor */
    P4EST_ASSERT (first_local_tree == -1 && last_local_tree == -2);
  }
  for (which_tree = first_local_tree + 1;       /* same type */
       which_tree <= last_local_tree; ++which_tree) {
    tree = p4est_tree_array_index (trees, which_tree);
    local_tree_last_quad_index[which_tree] = tree->quadrants.elem_count
      + local_tree_last_quad_index[which_tree - 1];
  }

#ifdef P4EST_ENABLE_DEBUG
  for (which_tree = first_local_tree; which_tree <= last_local_tree;
       ++which_tree) {
    tree = p4est_tree_array_index (trees, which_tree);
    P4EST_LDEBUGF
      ("partition tree %lld local_tree_last_quad_index[%lld] = %lld\n",
       (long long) which_tree, (long long) which_tree,
       (long long) local_tree_last_quad_index[which_tree]);
  }
#endif

  /* Calculate where the new quadrants are coming from */
  num_recv_from = P4EST_ALLOC_ZERO (p4est_locidx_t, num_procs);

  my_begin = (rank == 0) ? 0 : (new_global_last_quad_index[rank - 1] + 1);
  my_end = new_global_last_quad_index[rank];

  num_proc_recv_from = 0;

  if (my_begin > my_end) {
    /* my_begin == my_end requires a search is legal for find_partition */
    from_begin = rank;
    from_end = rank;
  }
  else {
    p4est_find_partition (num_procs, global_last_quad_index,
                          my_begin, my_end, &from_begin, &from_end);
    for (from_proc = from_begin; from_proc <= from_end; ++from_proc) {
      lower_bound =
        (from_proc == 0) ? 0 : (global_last_quad_index[from_proc - 1] + 1);

      if ((lower_bound <= my_end)
          && (global_last_quad_index[from_proc] >= my_begin)) {
        num_recv_from[from_proc] =
          SC_MIN (my_end,
                  global_last_quad_index[from_proc]) - SC_MAX (my_begin,
                                                               lower_bound)
          + 1;
        P4EST_ASSERT (num_recv_from[from_proc] >= 0);
        if (from_proc != rank)
          ++num_proc_recv_from;
      }
    }
  }

  from_begin_global_quad = from_begin;
  from_end_global_quad = from_end;

#ifdef P4EST_ENABLE_DEBUG
  /* old calculation method */
  {
    int                 old_num_proc_recv_from;
    p4est_locidx_t     *old_num_recv_from;
    p4est_gloidx_t      old_from_begin, old_from_end;

    old_num_recv_from = P4EST_ALLOC_ZERO (p4est_locidx_t, num_procs);

    old_num_proc_recv_from = 0;
    for (from_proc = 0; from_proc < num_procs; ++from_proc) {
      old_from_begin = (from_proc == 0) ?
        0 : (global_last_quad_index[from_proc - 1] + 1);
      old_from_end = global_last_quad_index[from_proc];

      if (old_from_begin <= my_end && old_from_end >= my_begin) {
        /* from_proc sends to me but may be empty */
        old_num_recv_from[from_proc] = SC_MIN (my_end, old_from_end)
          - SC_MAX (my_begin, old_from_begin) + 1;
        P4EST_ASSERT (old_num_recv_from[from_proc] >= 0);
        if (from_proc != rank)
          ++old_num_proc_recv_from;
      }
    }

    P4EST_ASSERT ((my_begin <= my_end) ? num_proc_recv_from ==
                  old_num_proc_recv_from : 1);
    for (i = 0; i < num_procs; ++i) {
      if (num_recv_from[i] != 0) {
        P4EST_LDEBUGF ("partition num_recv_from[%d] = %lld\n", i,
                       (long long) num_recv_from[i]);
      }
      /* compare the results of the new and old calculation method */
      /* old_from_begin/end and from_begin/end are not the same */
      P4EST_ASSERT (num_recv_from[i] == old_num_recv_from[i]);
    }
    P4EST_FREE (old_num_recv_from);
  }
#endif

  /* Post receives for the quadrants and their data */
  recv_buf = P4EST_ALLOC (char *, num_procs);
#ifdef P4EST_ENABLE_MPI
  recv_request = P4EST_ALLOC (MPI_Request, num_proc_recv_from);
#endif

  /* Allocate space for receiving quadrants and user data */
  for (from_proc = from_begin_global_quad
#ifdef P4EST_ENABLE_MPI
       , sk = 0
#endif
       ; from_proc <= from_end_global_quad; ++from_proc) {
    if (from_proc != rank && num_recv_from[from_proc]) {
      num_recv_trees =          /* same type */
        p4est->global_first_position[from_proc + 1].p.which_tree
        - p4est->global_first_position[from_proc].p.which_tree + 1;
      recv_size = num_recv_trees * sizeof (p4est_locidx_t)
        + quad_plus_data_size * num_recv_from[from_proc];

      recv_buf[from_proc] = P4EST_ALLOC (char, recv_size);

      /* Post receives for the quadrants and their data */
#ifdef P4EST_ENABLE_MPI
      P4EST_LDEBUGF ("partition recv %lld quadrants from %d\n",
                     (long long) num_recv_from[from_proc], from_proc);
      mpiret = MPI_Irecv (recv_buf[from_proc], (int) recv_size, MPI_BYTE,
                          from_proc, P4EST_COMM_PARTITION_GIVEN,
                          comm, recv_request + sk);
      SC_CHECK_MPI (mpiret);
      ++sk;
#endif
    }
  }
#ifdef P4EST_ENABLE_MPI
  for (; sk < num_proc_recv_from; ++sk) {
    /* for empty processors in receiving range */
    recv_request[sk] = MPI_REQUEST_NULL;
  }
#endif

  /* For each processor calculate the number of quadrants sent */
  num_send_to = P4EST_ALLOC_ZERO (p4est_locidx_t, num_procs);
  begin_send_to = P4EST_ALLOC (p4est_gloidx_t, num_procs);
#ifdef P4EST_ENABLE_DEBUG
  /* For checking the window of relevant processes. */
  memset (begin_send_to, -1, num_procs * sizeof (p4est_gloidx_t));
#endif

  my_begin = (rank == 0) ? 0 : (global_last_quad_index[rank - 1] + 1);
  my_end = global_last_quad_index[rank];

  num_proc_send_to = 0;

  if (my_begin > my_end) {
    /* my_begin == my_end requires a search is legal for find_partition */
    to_begin = rank;
    to_end = rank;
    memset (begin_send_to, -1, num_procs * sizeof (p4est_gloidx_t));
#ifdef P4EST_ENABLE_DEBUG
    send_to_empty = 1;
#endif
  }
  else {
    p4est_find_partition (num_procs, new_global_last_quad_index,
                          my_begin, my_end, &to_begin, &to_end);
#ifdef P4EST_ENABLE_DEBUG
    send_to_empty = 0;
#endif
    for (to_proc = to_begin; to_proc <= to_end; ++to_proc) {
      /* I send to to_proc which may be empty */
      lower_bound =
        (to_proc == 0) ? 0 : (new_global_last_quad_index[to_proc - 1] + 1);

      if ((lower_bound <= my_end)
          && (new_global_last_quad_index[to_proc] >= my_begin)) {
        num_send_to[to_proc] =
          SC_MIN (my_end,
                  new_global_last_quad_index[to_proc]) - SC_MAX (my_begin,
                                                                 lower_bound)
          + 1;
        begin_send_to[to_proc] = SC_MAX (my_begin, lower_bound);
        P4EST_ASSERT (num_send_to[to_proc] >= 0);
        if (to_proc != rank)
          ++num_proc_send_to;
      }
    }
  }

  to_begin_global_quad = to_begin;
  to_end_global_quad = to_end;

#ifdef P4EST_ENABLE_DEBUG
  /* old calculation method */
  {
    int                 old_num_proc_send_to;
    p4est_locidx_t     *old_num_send_to;
    p4est_gloidx_t      old_to_end, old_to_begin, *old_begin_send_to;

    old_num_send_to = P4EST_ALLOC (p4est_locidx_t, num_procs);
    old_begin_send_to = P4EST_ALLOC (p4est_gloidx_t, num_procs);

    old_num_proc_send_to = 0;
    for (to_proc = 0; to_proc < num_procs; ++to_proc) {
      old_to_begin = (to_proc == 0)
        ? 0 : (new_global_last_quad_index[to_proc - 1] + 1);
      old_to_end = new_global_last_quad_index[to_proc];

      if (old_to_begin <= my_end && old_to_end >= my_begin) {
        /* I send to to_proc which may be empty */
        old_num_send_to[to_proc] = SC_MIN (my_end, old_to_end)
          - SC_MAX (my_begin, old_to_begin) + 1;
        old_begin_send_to[to_proc] = SC_MAX (my_begin, old_to_begin);
        P4EST_ASSERT (old_num_send_to[to_proc] >= 0);
        if (to_proc != rank)
          ++old_num_proc_send_to;
      }
      else {
        /* I don't send to to_proc */
        old_num_send_to[to_proc] = 0;
        old_begin_send_to[to_proc] = -1;
      }
    }

    P4EST_ASSERT ((my_begin <= my_end) ? num_proc_send_to ==
                  old_num_proc_send_to : 1);
    for (i = 0; i < num_procs; ++i) {
      if (num_send_to[i] != 0) {
        P4EST_LDEBUGF ("partition num_send_to[%d] = %lld\n",
                       i, (long long) num_send_to[i]);
      }
      P4EST_ASSERT (num_send_to[i] == old_num_send_to[i]);
      if (begin_send_to[i] != -1) {
        P4EST_LDEBUGF ("partition begin_send_to[%d] = %lld\n",
                       i, (long long) begin_send_to[i]);
      }
      P4EST_ASSERT ((my_begin <= my_end) ? begin_send_to[i] ==
                    old_begin_send_to[i] : 1);
    }

    P4EST_FREE (old_num_send_to);
    P4EST_FREE (old_begin_send_to);
  }
#endif

  /* Communicate the quadrants and their data */
  send_buf = P4EST_ALLOC (char *, num_procs);
#ifdef P4EST_ENABLE_MPI
  send_request = P4EST_ALLOC (MPI_Request, num_proc_send_to);
#endif

  /* Set the num_per_tree_local */
  num_per_tree_local = P4EST_ALLOC_ZERO (p4est_locidx_t, num_send_trees);
  if (num_send_to[rank] > 0) {
    to_proc = rank;
    my_base = (rank == 0) ? 0 : (global_last_quad_index[rank - 1] + 1);
    my_begin = begin_send_to[to_proc] - my_base;
    my_end = begin_send_to[to_proc] + num_send_to[to_proc] - 1 - my_base;
    for (which_tree = first_local_tree; which_tree <= last_local_tree;
         ++which_tree) {
      tree = p4est_tree_array_index (trees, which_tree);

      from_begin = (which_tree == first_local_tree) ? 0 :
        (local_tree_last_quad_index[which_tree - 1] + 1);
      from_end = local_tree_last_quad_index[which_tree];

      if (from_begin <= my_end && from_end >= my_begin) {
        /* Need to copy from tree which_tree */
        tree_from_begin = SC_MAX (my_begin, from_begin) - from_begin;
        tree_from_end = SC_MIN (my_end, from_end) - from_begin;
        num_copy_global = tree_from_end - tree_from_begin + 1;
        P4EST_ASSERT (num_copy_global >= 0);
        P4EST_ASSERT (num_copy_global <= (p4est_gloidx_t) P4EST_LOCIDX_MAX);
        num_copy = (p4est_locidx_t) num_copy_global;
        num_per_tree_local[which_tree - first_local_tree] = num_copy;
      }
    }
  }

  /* Allocate space for receiving quadrants and user data */
  for (to_proc = to_begin_global_quad
#ifdef P4EST_ENABLE_MPI
       , sk = 0
#endif
       ; to_proc <= to_end_global_quad; ++to_proc) {
    if (to_proc != rank && num_send_to[to_proc]) {
      send_size = num_send_trees * sizeof (p4est_locidx_t)
        + quad_plus_data_size * num_send_to[to_proc];

      send_buf[to_proc] = P4EST_ALLOC (char, send_size);

      num_per_tree_send_buf = (p4est_locidx_t *) send_buf[to_proc];
      memset (num_per_tree_send_buf, 0,
              num_send_trees * sizeof (p4est_locidx_t));
      quad_send_buf = (p4est_quadrant_t *) (send_buf[to_proc]
                                            +
                                            num_send_trees *
                                            sizeof (p4est_locidx_t));
      user_data_send_buf =
        send_buf[to_proc] + num_send_trees * sizeof (p4est_locidx_t)
        + num_send_to[to_proc] * sizeof (p4est_quadrant_t);

      /* Pack in the data to be sent */

      my_base = (rank == 0) ? 0 : (global_last_quad_index[rank - 1] + 1);
      my_begin = begin_send_to[to_proc] - my_base;
      my_end = begin_send_to[to_proc] + num_send_to[to_proc] - 1 - my_base;

      for (which_tree = first_local_tree; which_tree <= last_local_tree;
           ++which_tree) {
        tree = p4est_tree_array_index (trees, which_tree);

        from_begin = (which_tree == first_local_tree) ? 0 :
          (local_tree_last_quad_index[which_tree - 1] + 1);
        from_end = local_tree_last_quad_index[which_tree];

        if (from_begin <= my_end && from_end >= my_begin) {
          /* Need to copy from tree which_tree */
          tree_from_begin = SC_MAX (my_begin, from_begin) - from_begin;
          tree_from_end = SC_MIN (my_end, from_end) - from_begin;
          num_copy = tree_from_end - tree_from_begin + 1;

          num_per_tree_send_buf[which_tree - first_local_tree] = num_copy;

          /* copy quads to send buf */
          memcpy (quad_send_buf, tree->quadrants.array +
                  tree_from_begin * sizeof (p4est_quadrant_t),
                  num_copy * sizeof (p4est_quadrant_t));

          /* set tree in send buf and copy user data */
          P4EST_LDEBUGF ("partition send %lld [%lld,%lld] quadrants"
                         " from tree %lld to proc %d\n",
                         (long long) num_copy, (long long) tree_from_begin,
                         (long long) tree_from_end, (long long) which_tree,
                         to_proc);
          for (il = 0; il < num_copy; ++il) {
            memcpy (user_data_send_buf + il * data_size,
                    quad_send_buf[il].p.user_data, data_size);
            if (data_size) {
              quad_send_buf[il].p.user_data = NULL;
            }
          }

          /* move pointer to beginning of quads that need to be copied */
          my_begin += num_copy;
          quad_send_buf += num_copy;
          user_data_send_buf += num_copy * data_size;
        }
      }

      /* Post send operation for the quadrants and their data */
#ifdef P4EST_ENABLE_MPI
      P4EST_LDEBUGF ("partition send %lld quadrants to %d\n",
                     (long long) num_send_to[to_proc], to_proc);
      mpiret = MPI_Isend (send_buf[to_proc], (int) send_size, MPI_BYTE,
                          to_proc, P4EST_COMM_PARTITION_GIVEN,
                          comm, send_request + sk);
      SC_CHECK_MPI (mpiret);
      ++sk;
#endif
    }
  }
#ifdef P4EST_ENABLE_MPI
  for (; sk < num_proc_send_to; ++sk) {
    send_request[sk] = MPI_REQUEST_NULL;
  }

  /* Fill in forest */
  mpiret =
    sc_MPI_Waitall (num_proc_recv_from, recv_request, MPI_STATUSES_IGNORE);
  SC_CHECK_MPI (mpiret);
#endif

  /* Loop through and fill in */

  /* Calculate the local index of the end of each tree in the repartition */
  new_local_tree_elem_count =
    P4EST_ALLOC_ZERO (p4est_locidx_t, trees->elem_count);
  new_local_tree_elem_count_before =
    P4EST_ALLOC_ZERO (p4est_locidx_t, trees->elem_count);
  new_first_local_tree = (p4est_topidx_t) P4EST_TOPIDX_MAX;
  new_last_local_tree = 0;

  for (from_proc = from_begin_global_quad; from_proc <= from_end_global_quad;
       ++from_proc) {
    if (num_recv_from[from_proc] > 0) {
      first_from_tree = p4est->global_first_position[from_proc].p.which_tree;
      last_from_tree =
        p4est->global_first_position[from_proc + 1].p.which_tree;
      num_recv_trees =          /* same type */
        last_from_tree - first_from_tree + 1;

      P4EST_LDEBUGF
        ("partition from %d with trees [%lld,%lld] get %lld trees\n",
         from_proc, (long long) first_from_tree, (long long) last_from_tree,
         (long long) num_recv_trees);
      num_per_tree_recv_buf = (from_proc == rank) ?
        num_per_tree_local : (p4est_locidx_t *) recv_buf[from_proc];

      for (it = 0; it < num_recv_trees; ++it) {

        if (num_per_tree_recv_buf[it] > 0) {
          from_tree = first_from_tree + it;     /* same type */

          P4EST_ASSERT (from_tree >= 0
                        && from_tree < (p4est_topidx_t) trees->elem_count);
          P4EST_LDEBUGF ("partition recv %lld [%lld,%lld] quadrants"
                         " from tree %lld from proc %d\n",
                         (long long) num_per_tree_recv_buf[it],
                         (long long) new_local_tree_elem_count[from_tree],
                         (long long) new_local_tree_elem_count[from_tree]
                         + num_per_tree_recv_buf[it], (long long) from_tree,
                         from_proc);
          new_first_local_tree =        /* same type */
            SC_MIN (new_first_local_tree, from_tree);
          new_last_local_tree = /* same type */
            SC_MAX (new_last_local_tree, from_tree);
          new_local_tree_elem_count[from_tree] +=       /* same type */
            num_per_tree_recv_buf[it];
          new_local_tree_elem_count_before[from_tree] +=        /* same type */
            (from_proc < rank) ? num_per_tree_recv_buf[it] : 0;
        }
      }
    }
  }
  if (new_first_local_tree > new_last_local_tree) {
    new_first_local_tree = -1;
    new_last_local_tree = -2;
  }
  P4EST_VERBOSEF ("partition new forest [%lld,%lld]\n",
                  (long long) new_first_local_tree,
                  (long long) new_last_local_tree);

  /* Copy the local quadrants */
  if (first_local_tree >= 0 && new_first_local_tree >= 0) {
    P4EST_ASSERT (last_local_tree >= 0 && new_last_local_tree >= 0);
    first_tree =                /* same type */
      SC_MIN (first_local_tree, new_first_local_tree);
  }
  else {
    P4EST_ASSERT (last_local_tree == -2 || new_last_local_tree == -2);
    first_tree =                /* same type */
      SC_MAX (first_local_tree, new_first_local_tree);
  }
  last_tree =                   /* same type */
    SC_MAX (last_local_tree, new_last_local_tree);
  my_base = (rank == 0) ? 0 : (global_last_quad_index[rank - 1] + 1);
  /* Although begin_send_to[rank] may be -1 if my_begin > my_end was true above
   * it is valid to use just 0 in the conditional expressions below since
   * my_begin > my_end is also true for the new values and since it holds
   * 0 <= from_begin <= from_end and therefore the related if-statement is
   * false as in the old version of the code (cf. debug mode) and the exact
   * values of the new my_begin and my_end values are not needed for empty
   * processors.
   */
  my_begin = ((to_begin_global_quad <= rank && to_end_global_quad >= rank) ?
              begin_send_to[rank] : 0) - my_base;
  my_end = ((to_begin_global_quad <= rank && to_end_global_quad >= rank) ?
            begin_send_to[rank] : 0) + num_send_to[rank] - 1 - my_base;

  for (which_tree = first_tree; which_tree <= last_tree; ++which_tree) {
    tree = p4est_tree_array_index (trees, which_tree);
    quadrants = &tree->quadrants;

    if (new_local_tree_elem_count[which_tree] > 0) {
      if (which_tree >= first_local_tree && which_tree <= last_local_tree) {

        num_quadrants = new_local_tree_elem_count[which_tree];

        from_begin = (which_tree == first_local_tree) ? 0 :
          (local_tree_last_quad_index[which_tree - 1] + 1);
        from_end = local_tree_last_quad_index[which_tree];

#ifdef P4EST_ENABLE_DEBUG
        if (send_to_empty) {
          /* The corner case that begin_send_to[rank] == -1 holds. Therefore,
           * we need to ensure that from_begin <= my_end && from_end >= my_begin
           * is still evaluated to false even if begin_send_to[rank] == 0 was
           * assumed in the calculation of my_{begin,end}.
           * The reasoning of this is already presented above but
           * here we check this resoning again with an assertion.
           */
          P4EST_ASSERT (!(from_begin <= my_end && from_end >= my_begin));
          /* We also check if the evaluation of the expression mentioned
           * above coincides with the evaluation without the adjustment
           * in the calculation of my_{begin,end}. Both assertions combined
           * give us that the behaviour of the code is not affected.
           */
          my_begin_comp = -1 - my_base;
          my_end_comp = -1 + num_send_to[rank] - 1 - my_base;
          P4EST_ASSERT (!(from_begin <= my_end_comp &&
                          from_end >= my_begin_comp));
        }
#endif
        if (from_begin <= my_end && from_end >= my_begin) {
          /* Need to keep part of tree which_tree */
          tree_from_begin = SC_MAX (my_begin, from_begin) - from_begin;
          tree_from_end = SC_MIN (my_end, from_end) - from_begin;
          num_copy = tree_from_end - tree_from_begin + 1;
        }
        else {
          tree_from_begin = 0;
          tree_from_end = -1;
          num_copy = 0;
        }

        /* Free all userdata that no longer belongs to this process */
        zoffset = SC_MIN ((size_t) tree_from_begin, quadrants->elem_count);
        for (zz = 0; zz < zoffset; ++zz) {
          quad = p4est_quadrant_array_index (quadrants, zz);
          p4est_quadrant_free_data (p4est, quad);
        }
        zoffset = (size_t) tree_from_end + 1;
        for (zz = zoffset; zz < quadrants->elem_count; ++zz) {
          quad = p4est_quadrant_array_index (quadrants, zz);
          p4est_quadrant_free_data (p4est, quad);
        }

        if (num_quadrants > (p4est_locidx_t) quadrants->elem_count) {
          sc_array_resize (quadrants, num_quadrants);
        }

        P4EST_LDEBUGF ("copying %lld local quads to tree %lld\n",
                       (long long) num_copy, (long long) which_tree);
        P4EST_LDEBUGF
          ("   with %lld(%llu) quads from [%lld, %lld] to [%lld, %lld]\n",
           (long long) num_quadrants,
           (unsigned long long) quadrants->elem_count,
           (long long) tree_from_begin, (long long) tree_from_end,
           (long long) new_local_tree_elem_count_before[which_tree],
           (long long) new_local_tree_elem_count_before[which_tree] +
           num_copy - 1);
        memmove (quadrants->array +
                 new_local_tree_elem_count_before[which_tree] *
                 sizeof (p4est_quadrant_t),
                 quadrants->array +
                 tree_from_begin * sizeof (p4est_quadrant_t),
                 num_copy * sizeof (p4est_quadrant_t));

        if (num_quadrants < (p4est_locidx_t) quadrants->elem_count) {
          sc_array_resize (quadrants, num_quadrants);
        }
      }
    }
    else {
      /*
       * Check to see if we need to drop a tree because we no longer have
       * any quadrants in it.
       */
      if (which_tree >= first_local_tree && which_tree <= last_local_tree) {
        /* Free all userdata that no longer belongs to this process */
        for (zz = 0; zz < quadrants->elem_count; ++zz) {
          quad = p4est_quadrant_array_index (quadrants, zz);
          p4est_quadrant_free_data (p4est, quad);
        }

        /* The whole tree is dropped */
        P4EST_QUADRANT_INIT (&tree->first_desc);
        P4EST_QUADRANT_INIT (&tree->last_desc);
        sc_array_reset (quadrants);
        tree->quadrants_offset = 0;
        for (i = 0; i <= P4EST_QMAXLEVEL; ++i) {
          tree->quadrants_per_level[i] = 0;
        }
        tree->maxlevel = 0;
      }
    }
  }

  /* Copy in received quadrants */

  memset (new_local_tree_elem_count_before, 0,
          trees->elem_count * sizeof (p4est_locidx_t));
  for (from_proc = from_begin_global_quad; from_proc <= from_end_global_quad;
       ++from_proc) {
    if (num_recv_from[from_proc] > 0) {
      first_from_tree = p4est->global_first_position[from_proc].p.which_tree;
      last_from_tree =
        p4est->global_first_position[from_proc + 1].p.which_tree;
      num_recv_trees =          /* same type */
        last_from_tree - first_from_tree + 1;

      P4EST_LDEBUGF
        ("partition copy from %d with trees [%lld,%lld] get %lld trees\n",
         from_proc, (long long) first_from_tree,
         (long long) last_from_tree, (long long) num_recv_trees);
      num_per_tree_recv_buf =
        (from_proc == rank) ? num_per_tree_local :
        (p4est_locidx_t *) recv_buf[from_proc];

      quad_recv_buf = (p4est_quadrant_t *) (recv_buf[from_proc]
                                            + num_recv_trees *
                                            sizeof (p4est_locidx_t));
      user_data_recv_buf =
        recv_buf[from_proc] + num_recv_trees * sizeof (p4est_locidx_t)
        + num_recv_from[from_proc] * sizeof (p4est_quadrant_t);

      for (it = 0; it < num_recv_trees; ++it) {
        from_tree = first_from_tree + it;       /* same type */
        num_copy = num_per_tree_recv_buf[it];

        /* We might have sent trees that are not actual trees.  In
         * this case the num_copy should be zero
         */
        P4EST_ASSERT (num_copy == 0
                      || (num_copy > 0 && from_tree >= 0
                          && from_tree < (p4est_topidx_t) trees->elem_count));

        if (num_copy > 0 && rank != from_proc) {
          tree = p4est_tree_array_index (trees, from_tree);
          quadrants = &tree->quadrants;
          num_quadrants = new_local_tree_elem_count[from_tree];
          sc_array_resize (quadrants, num_quadrants);

          /* copy quadrants */
          P4EST_LDEBUGF ("copying %lld remote quads to tree %lld"
                         " with %lld quads from proc %d\n",
                         (long long) num_copy, (long long) from_tree,
                         (long long) num_quadrants, from_proc);
          memcpy (quadrants->array +
                  new_local_tree_elem_count_before[from_tree]
                  * sizeof (p4est_quadrant_t), quad_recv_buf,
                  num_copy * sizeof (p4est_quadrant_t));

          /* copy user data */
          zoffset = (size_t) new_local_tree_elem_count_before[from_tree];
          for (zz = 0; zz < (size_t) num_copy; ++zz) {
            quad = p4est_quadrant_array_index (quadrants, zz + zoffset);

            if (data_size > 0) {
              quad->p.user_data = sc_mempool_alloc (p4est->user_data_pool);
              memcpy (quad->p.user_data, user_data_recv_buf + zz * data_size,
                      data_size);
            }
          }
        }

        if (num_copy > 0) {
          P4EST_ASSERT (from_tree >= 0
                        && from_tree < (p4est_topidx_t) trees->elem_count);
          new_local_tree_elem_count_before[from_tree] +=        /* same type */
            num_copy;
        }

        /* increment the recv quadrant pointers */
        quad_recv_buf += num_copy;
        user_data_recv_buf += num_copy * data_size;
      }
    }
  }

  /* Set the global index and count of quadrants instead
   * of calling p4est_comm_count_quadrants
   */
  P4EST_FREE (global_last_quad_index);
  global_last_quad_index = new_global_last_quad_index;
  P4EST_ASSERT (p4est->global_num_quadrants ==
                new_global_last_quad_index[num_procs - 1] + 1);
  P4EST_ASSERT (p4est->global_first_quadrant[0] == 0);
  for (i = 0; i < num_procs; ++i) {
    p4est->global_first_quadrant[i + 1] = global_last_quad_index[i] + 1;
  }
  P4EST_FREE (new_global_last_quad_index);
  global_last_quad_index = new_global_last_quad_index = NULL;

  p4est->first_local_tree = new_first_local_tree;
  p4est->last_local_tree = new_last_local_tree;

  new_local_num_quadrants = 0;
  for (which_tree = 0; which_tree < new_first_local_tree; ++which_tree) {
    tree = p4est_tree_array_index (trees, which_tree);
    tree->quadrants_offset = 0;
    P4EST_QUADRANT_INIT (&tree->first_desc);
    P4EST_QUADRANT_INIT (&tree->last_desc);
  }
  for (; which_tree <= new_last_local_tree; ++which_tree) {
    tree = p4est_tree_array_index (trees, which_tree);
    tree->quadrants_offset = new_local_num_quadrants;
    quadrants = &tree->quadrants;
    P4EST_ASSERT (quadrants->elem_count > 0);

    new_local_num_quadrants +=  /* same type */
      (p4est_locidx_t) quadrants->elem_count;

    for (i = 0; i <= P4EST_QMAXLEVEL; ++i) {
      tree->quadrants_per_level[i] = 0;
    }
    tree->maxlevel = 0;
    for (zz = 0; zz < quadrants->elem_count; ++zz) {
      quad = p4est_quadrant_array_index (quadrants, zz);
      ++tree->quadrants_per_level[quad->level];
      tree->maxlevel = (int8_t) SC_MAX (quad->level, tree->maxlevel);
    }

    quad = p4est_quadrant_array_index (quadrants, 0);
    p4est_quadrant_first_descendant (quad, &tree->first_desc,
                                     P4EST_QMAXLEVEL);
    quad = p4est_quadrant_array_index (quadrants, quadrants->elem_count - 1);
    p4est_quadrant_last_descendant (quad, &tree->last_desc, P4EST_QMAXLEVEL);
  }
  for (; which_tree < p4est->connectivity->num_trees; ++which_tree) {
    tree = p4est_tree_array_index (trees, which_tree);
    tree->quadrants_offset = new_local_num_quadrants;
    P4EST_QUADRANT_INIT (&tree->first_desc);
    P4EST_QUADRANT_INIT (&tree->last_desc);
  }
  p4est->local_num_quadrants = new_local_num_quadrants;

  /* Clean up */

#ifdef P4EST_ENABLE_MPI
  mpiret =
    sc_MPI_Waitall (num_proc_send_to, send_request, MPI_STATUSES_IGNORE);
  SC_CHECK_MPI (mpiret);

#ifdef P4EST_ENABLE_DEBUG
  for (i = 0; i < num_proc_recv_from; ++i) {
    P4EST_ASSERT (recv_request[i] == MPI_REQUEST_NULL);
  }
  for (i = 0; i < num_proc_send_to; ++i) {
    P4EST_ASSERT (send_request[i] == MPI_REQUEST_NULL);
  }
#endif
  P4EST_FREE (recv_request);
  P4EST_FREE (send_request);
#endif

  for (i = from_begin_global_quad; i <= from_end_global_quad; ++i) {
    if (i != rank && num_recv_from[i])
      P4EST_FREE (recv_buf[i]);
  }

  for (i = to_begin_global_quad; i <= to_end_global_quad; ++i) {
    if (i != rank && num_send_to[i])
      P4EST_FREE (send_buf[i]);
  }

  P4EST_FREE (num_per_tree_local);
  P4EST_FREE (local_tree_last_quad_index);
  P4EST_FREE (new_local_tree_elem_count);
  P4EST_FREE (new_local_tree_elem_count_before);
  P4EST_FREE (recv_buf);
  P4EST_FREE (send_buf);
  P4EST_FREE (num_recv_from);
  P4EST_FREE (num_send_to);
  P4EST_FREE (begin_send_to);

  p4est_comm_global_partition (p4est, NULL);

  /* Assert that we have a valid partition */
  P4EST_ASSERT (crc == p4est_checksum (p4est));
  P4EST_GLOBAL_INFOF
    ("Done " P4EST_STRING
     "_partition_given shipped %lld quadrants %.3g%%\n",
     (long long) total_quadrants_shipped,
     total_quadrants_shipped * 100. / p4est->global_num_quadrants);

  return total_quadrants_shipped;
}

int
p4est_quadrant_on_face_boundary (p4est_t * p4est, p4est_topidx_t treeid,
                                 int face, const p4est_quadrant_t * q)
{
  p4est_qcoord_t      dh, xyz;
  p4est_connectivity_t *conn = p4est->connectivity;

  P4EST_ASSERT (0 <= face && face < P4EST_FACES);
  P4EST_ASSERT (p4est_quadrant_is_valid (q));

  if (conn->tree_to_tree[P4EST_FACES * treeid + face] != treeid ||
      (int) conn->tree_to_face[P4EST_FACES * treeid + face] != face) {
    return 0;
  }

  dh = P4EST_LAST_OFFSET (q->level);
  switch (face / 2) {
  case 0:
    xyz = q->x;
    break;
  case 1:
    xyz = q->y;
    break;
#ifdef P4_TO_P8
  case 2:
    xyz = q->z;
    break;
#endif
  default:
    SC_ABORT_NOT_REACHED ();
    break;
  }
  return xyz == ((face & 0x01) ? dh : 0);
}
