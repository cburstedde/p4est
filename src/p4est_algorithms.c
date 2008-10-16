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

#ifdef P4_TO_P8
#include <p8est_algorithms.h>
#include <p8est_bits.h>
#include <p8est_communication.h>
#else
#include <p4est_algorithms.h>
#include <p4est_bits.h>
#include <p4est_communication.h>
#endif /* !P4_TO_P8 */

/* htonl is in either of these two */
#ifdef P4EST_HAVE_ARPA_NET_H
#include <arpa/inet.h>
#endif
#ifdef P4EST_HAVE_NETINET_IN_H
#include <netinet/in.h>
#endif

/* *INDENT-OFF* */

#ifndef P4_TO_P8

/** The offsets of the 3 indirect neighbors in units of h.
 * Indexing [cid][neighbor][xy] where cid is the child id.
 * Neighbors are indexed in z-order.
 */
static const int    indirect_neighbors[4][3][2] =
{{{-1, -1}, { 1, -1}, {-1, 1}},
 {{ 0, -1}, { 2, -1}, { 1, 0}},
 {{-1,  0}, {-2,  1}, { 0, 1}},
 {{ 1, -1}, {-1,  1}, { 1, 1}}};

/** Indicate which neighbor to omit if edges are balanced, not corners
 * Indexing [cid] where cid is the child id.
 */
static const int    corners_omitted[4] =
{ 0, 1, 1, 2 };

#else /* P4_TO_P8 */

/** Store the number of quadrants to add for complete and balance stages. */
static const int    p8est_balance_count[P4EST_DIM + 1] =
{ 9, 12, 15, 16 };

/** Store coordinates of quadrants to add for balancing. */
static const p4est_qcoord_t p8est_balance_coord[26][P4EST_DIM] =
{ /* faces */
  { -1,  1,  1 },
  {  2,  0,  0 },
  {  1, -1,  1 },
  {  0,  2,  0 },
  {  1,  1, -1 },
  {  0,  0,  2 },
  /* edges */
  {  1, -1, -1 },
  {  1,  2, -1 },
  {  0, -1,  2 },
  {  0,  2,  2 },
  { -1,  1, -1 },
  {  2,  1, -1 },
  { -1,  0,  2 },
  {  2,  0,  2 },
  { -1, -1,  1 },
  {  2, -1,  1 },
  { -1,  2,  0 },
  {  2,  2,  0 },
  /* corners */
  { -1, -1, -1 },
  {  2, -1, -1 },
  { -1,  2, -1 },
  {  2,  2, -1 },
  { -1, -1,  2 },
  {  2, -1,  2 },
  { -1,  2,  2 },
  {  2,  2,  2 }};

#endif /* P4_TO_P8 */

/* *INDENT-ON* */

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
  bool                own_check;
  size_t              kz, qcount;
  unsigned            crc;
  uint32_t           *check;
  p4est_quadrant_t   *q;

  qcount = quadrants->elem_count;

  P4EST_ASSERT (quadrants->elem_size == sizeof (p4est_quadrant_t));
  P4EST_ASSERT (first_quadrant <= qcount);

  if (checkarray == NULL) {
    checkarray = sc_array_new (4);
    own_check = true;
  }
  else {
    P4EST_ASSERT (checkarray->elem_size == 4);
    own_check = false;
  }

  sc_array_resize (checkarray, (qcount - first_quadrant) * (P4EST_DIM + 1));
  for (kz = first_quadrant; kz < qcount; ++kz) {
    q = sc_array_index (quadrants, kz);
    P4EST_ASSERT (p4est_quadrant_is_extended (q));
    check =
      sc_array_index (checkarray, (kz - first_quadrant) * (P4EST_DIM + 1));
    check[0] = htonl ((uint32_t) q->x);
    check[1] = htonl ((uint32_t) q->y);
#ifdef P4_TO_P8
    check[2] = htonl ((uint32_t) q->z);
#endif
    check[P4EST_DIM] = htonl ((uint32_t) q->level);
  }
  crc = sc_array_checksum (checkarray, 0);

  if (own_check) {
    sc_array_destroy (checkarray);
  }

  return crc;
}

bool
p4est_tree_is_sorted (p4est_tree_t * tree)
{
  size_t              iz;
  p4est_quadrant_t   *q1, *q2;
  sc_array_t         *tquadrants = &tree->quadrants;

  if (tquadrants->elem_count <= 1) {
    return true;
  }

  q1 = sc_array_index (tquadrants, 0);
  for (iz = 1; iz < tquadrants->elem_count; ++iz) {
    q2 = sc_array_index (tquadrants, iz);
    if (p4est_quadrant_compare (q1, q2) >= 0) {
      return false;
    }
    q1 = q2;
  }

  return true;
}

bool
p4est_tree_is_linear (p4est_tree_t * tree)
{
  size_t              iz;
  p4est_quadrant_t   *q1, *q2;
  sc_array_t         *tquadrants = &tree->quadrants;

  if (tquadrants->elem_count <= 1) {
    return true;
  }

  q1 = sc_array_index (tquadrants, 0);
  for (iz = 1; iz < tquadrants->elem_count; ++iz) {
    q2 = sc_array_index (tquadrants, iz);
    if (p4est_quadrant_compare (q1, q2) >= 0) {
      return false;
    }
    if (p4est_quadrant_is_ancestor (q1, q2)) {
      return false;
    }
    q1 = q2;
  }

  return true;
}

bool
p4est_tree_is_almost_sorted (p4est_tree_t * tree, bool check_linearity)
{
  size_t              iz;
  int                 face_contact1;
  int                 face_contact2;
  int                 out_axis[P4EST_DIM];
  p4est_quadrant_t   *q1, *q2;
  sc_array_t         *tquadrants = &tree->quadrants;

  if (tquadrants->elem_count <= 1) {
    return true;
  }

  q1 = sc_array_index (tquadrants, 0);
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
    q2 = sc_array_index (tquadrants, iz);
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
        return false;
      }
      if (check_linearity && p4est_quadrant_is_ancestor (q1, q2)) {
        return false;
      }
    }
    q1 = q2;
    face_contact1 = face_contact2;
  }

  return true;
}

bool
p4est_tree_is_complete (p4est_tree_t * tree)
{
  size_t              iz;
  p4est_quadrant_t   *q1, *q2;
  sc_array_t         *tquadrants = &tree->quadrants;

  if (tquadrants->elem_count <= 1) {
    return true;
  }

  q1 = sc_array_index (tquadrants, 0);
  for (iz = 1; iz < tquadrants->elem_count; ++iz) {
    q2 = sc_array_index (tquadrants, iz);
    if (!p4est_quadrant_is_next (q1, q2)) {
      return false;
    }
    q1 = q2;
  }

  return true;
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
    q2 = sc_array_index (tquadrants, jz);
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
    P4EST_NORMAL_LOG (log_priority, buffer);
    q1 = q2;
  }
}

bool
p4est_is_valid (p4est_t * p4est)
{
#ifdef P4EST_DEBUG
  const int           num_procs = p4est->mpisize;
#endif
  const int           rank = p4est->mpirank;
  const p4est_topidx_t first_tree = p4est->first_local_tree;
  const p4est_topidx_t last_tree = p4est->last_local_tree;
  int                 i, maxlevel;
  bool                failed;
  size_t              jz;
  p4est_topidx_t      next_tree;
  p4est_locidx_t      nquadrants, lquadrants, perlevel;
  p4est_qcoord_t      mh = P4EST_QUADRANT_LEN (P4EST_QMAXLEVEL);
  p4est_quadrant_t   *q;
  p4est_quadrant_t    mylow, nextlow, s;
  p4est_tree_t       *tree;

  failed = false;
  P4EST_QUADRANT_INIT (&mylow);
  P4EST_QUADRANT_INIT (&nextlow);
  P4EST_QUADRANT_INIT (&s);

#ifdef P4EST_DEBUG
  /* check last item of global partition */
  P4EST_ASSERT (p4est->global_first_position[num_procs].p.which_tree ==
                p4est->connectivity->num_trees &&
                p4est->global_first_position[num_procs].x == 0 &&
                p4est->global_first_position[num_procs].y == 0);
#ifdef P4_TO_P8
  P4EST_ASSERT (p4est->global_first_position[num_procs].z == 0);
#endif
  P4EST_ASSERT (p4est->connectivity->num_trees ==
                (p4est_topidx_t) p4est->trees->elem_count);
  for (i = 0; i <= num_procs; ++i) {
    P4EST_ASSERT (p4est->global_first_position[i].level == P4EST_QMAXLEVEL);
  }
#endif /* P4EST_DEBUG */

  /* check first tree in global partition */
  if (first_tree < 0) {
    if (!(first_tree == -1 && last_tree == -2)) {
      P4EST_NOTICE ("p4est invalid empty tree range A");
      failed = true;
      goto failtest;
    }
  }
  else {
    if (p4est->global_first_position[rank].p.which_tree != first_tree) {
      P4EST_NOTICE ("p4est invalid first tree\n");
      failed = true;
      goto failtest;
    }
    mylow.x = p4est->global_first_position[rank].x;
    mylow.y = p4est->global_first_position[rank].y;
#ifdef P4_TO_P8
    mylow.z = p4est->global_first_position[rank].z;
#endif
    mylow.level = P4EST_QMAXLEVEL;
    tree = sc_array_index (p4est->trees, first_tree);
    if (tree->quadrants.elem_count > 0) {
      q = sc_array_index (&tree->quadrants, 0);
      if (q->x != mylow.x || q->y != mylow.y ||
#ifdef P4_TO_P8
          q->z != mylow.z ||
#endif
          false) {
        P4EST_NOTICE ("p4est invalid low quadrant\n");
        failed = true;
        goto failtest;
      }
    }
  }

  /* check last tree in global partition */
  if (last_tree < 0) {
    if (!(first_tree == -1 && last_tree == -2)) {
      P4EST_NOTICE ("p4est invalid empty tree range B");
      failed = true;
      goto failtest;
    }
  }
  else {
    next_tree = p4est->global_first_position[rank + 1].p.which_tree;
    if (next_tree != last_tree && next_tree != last_tree + 1) {
      P4EST_NOTICE ("p4est invalid last tree\n");
      failed = true;
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
        failed = true;
        goto failtest;
      }
    }
    tree = sc_array_index (p4est->trees, last_tree);
    if (tree->quadrants.elem_count > 0) {
      q = sc_array_index (&tree->quadrants, tree->quadrants.elem_count - 1);
      if (next_tree == last_tree) {
        if (!p4est_quadrant_is_next (q, &nextlow)) {
          P4EST_NOTICE ("p4est invalid next quadrant\n");
          failed = true;
          goto failtest;
        }
      }
      else {
        p4est_quadrant_last_descendent (q, &s, P4EST_QMAXLEVEL);
        if (s.x + mh != P4EST_ROOT_LEN || s.y + mh != P4EST_ROOT_LEN ||
#ifdef P4_TO_P8
            s.z + mh != P4EST_ROOT_LEN ||
#endif
            false) {
          P4EST_NOTICE ("p4est invalid last quadrant\n");
          failed = true;
          goto failtest;
        }
      }
    }
  }

  /* check individual trees */
  lquadrants = 0;
  for (jz = 0; jz < p4est->trees->elem_count; ++jz) {
    tree = sc_array_index (p4est->trees, jz);
    if (tree->quadrants_offset != lquadrants) {
      P4EST_NOTICE ("p4est invalid quadrants offset\n");
      failed = true;
      goto failtest;
    }
    if (!p4est_tree_is_complete (tree)) {
      P4EST_NOTICE ("p4est invalid not complete\n");
      failed = true;
      goto failtest;
    }
    if (tree->quadrants.elem_count > 0) {
      if (((p4est_topidx_t) jz < p4est->first_local_tree ||
           (p4est_topidx_t) jz > p4est->last_local_tree)) {
        P4EST_NOTICE ("p4est invalid outside count\n");
        failed = true;
        goto failtest;
      }
      q = sc_array_index (&tree->quadrants, 0);
      p4est_quadrant_first_descendent (q, &s, P4EST_QMAXLEVEL);
      if (!p4est_quadrant_is_equal (&s, &tree->first_desc)) {
        P4EST_NOTICE ("p4est invalid first tree descendent\n");
        failed = true;
        goto failtest;
      }
      q = sc_array_index (&tree->quadrants, tree->quadrants.elem_count - 1);
      p4est_quadrant_last_descendent (q, &s, P4EST_QMAXLEVEL);
      if (!p4est_quadrant_is_equal (&s, &tree->last_desc)) {
        P4EST_NOTICE ("p4est invalid last tree descendent\n");
        failed = true;
        goto failtest;
      }
    }
    else {
      P4EST_QUADRANT_INIT (&s);
      if (s.level != tree->first_desc.level ||
          s.level != tree->last_desc.level) {
        P4EST_NOTICE ("p4est invalid empty descendent\n");
        failed = true;
        goto failtest;
      }
    }

    maxlevel = 0;
    nquadrants = 0;
    for (i = 0; i <= P4EST_QMAXLEVEL; ++i) {
      perlevel = tree->quadrants_per_level[i];

      P4EST_ASSERT (perlevel >= 0);
      nquadrants += perlevel;   /* same type */
      if (perlevel > 0) {
        maxlevel = i;
      }
    }
    for (; i <= P4EST_MAXLEVEL; ++i) {
      P4EST_ASSERT (tree->quadrants_per_level[i] == -1);
    }
    lquadrants += nquadrants;   /* same type */

    if (maxlevel != (int) tree->maxlevel) {
      P4EST_NOTICE ("p4est invalid wrong maxlevel\n");
      failed = true;
      goto failtest;
    }
    if (nquadrants != (p4est_locidx_t) tree->quadrants.elem_count) {
      P4EST_NOTICE ("p4est invalid tree quadrant count\n");
      failed = true;
      goto failtest;
    }
  }

  if (lquadrants != p4est->local_num_quadrants) {
    P4EST_NOTICE ("p4est invalid local quadrant count\n");
    failed = true;
    goto failtest;
  }

failtest:
  return !p4est_comm_sync_flag (p4est, failed, MPI_BOR);
}

/* here come the heavyweight algorithms */

ssize_t
p4est_find_lower_bound (sc_array_t * array,
                        const p4est_quadrant_t * q, size_t guess)
{
  int                 comp;
  size_t              count;
  size_t              quad_low, quad_high;
  p4est_quadrant_t   *cur;

  count = array->elem_count;
  if (count == 0)
    return -1;

  quad_low = 0;
  quad_high = count - 1;

  for (;;) {
    P4EST_ASSERT (quad_low <= quad_high);
    P4EST_ASSERT (quad_low < count && quad_high < count);
    P4EST_ASSERT (quad_low <= guess && guess <= quad_high);

    /* compare two quadrants */
    cur = sc_array_index (array, guess);
    comp = p4est_quadrant_compare (q, cur);

    /* check if guess is higher or equal q and there's room below it */
    if (comp <= 0 && (guess > 0 && p4est_quadrant_compare (q, cur - 1) <= 0)) {
      quad_high = guess - 1;
      guess = (quad_low + quad_high + 1) / 2;
      continue;
    }

    /* check if guess is lower than q */
    if (comp > 0) {
      quad_low = guess + 1;
      if (quad_low > quad_high)
        return -1;

      guess = (quad_low + quad_high) / 2;
      continue;
    }

    /* otherwise guess is the correct quadrant */
    break;
  }

  return (ssize_t) guess;
}

ssize_t
p4est_find_higher_bound (sc_array_t * array,
                         const p4est_quadrant_t * q, size_t guess)
{
  int                 comp;
  size_t              count;
  size_t              quad_low, quad_high;
  p4est_quadrant_t   *cur;

  count = array->elem_count;
  if (count == 0)
    return -1;

  quad_low = 0;
  quad_high = count - 1;

  for (;;) {
    P4EST_ASSERT (quad_low <= quad_high);
    P4EST_ASSERT (quad_low < count && quad_high < count);
    P4EST_ASSERT (quad_low <= guess && guess <= quad_high);

    /* compare two quadrants */
    cur = sc_array_index (array, guess);
    comp = p4est_quadrant_compare (cur, q);

    /* check if guess is lower or equal q and there's room above it */
    if (comp <= 0 &&
        (guess < count - 1 && p4est_quadrant_compare (cur + 1, q) <= 0)) {
      quad_low = guess + 1;
      guess = (quad_low + quad_high) / 2;
      continue;
    }

    /* check if guess is higher than q */
    if (comp > 0) {
      if (guess == 0)
        return -1;

      quad_high = guess - 1;
      if (quad_high < quad_low)
        return -1;

      guess = (quad_low + quad_high + 1) / 2;
      continue;
    }

    /* otherwise guess is the correct quadrant */
    break;
  }

  return (ssize_t) guess;
}

void
p4est_tree_compute_overlap (p4est_t * p4est, p4est_topidx_t qtree,
                            sc_array_t * in, sc_array_t * out)
{
  int                 k, l, m, which;
  int                 face, corner, level;
  size_t              iz, ctree;
  size_t              treecount, incount;
  size_t              guess;
  ssize_t             first_index, last_index, js;
  bool                inter_tree;
  bool                outface[2 * P4EST_DIM];
  p4est_topidx_t      ntree;
  p4est_qcoord_t      qh;
  p4est_quadrant_t    fd, ld, tempq, ins[P4EST_INSUL];
  p4est_quadrant_t   *treefd, *treeld;
  p4est_quadrant_t   *tq, *s;
  p4est_quadrant_t   *inq, *outq;
  p4est_tree_t       *tree;
  p4est_connectivity_t *conn = p4est->connectivity;
#ifndef P4_TO_P8
  int                 transform;
  p4est_corner_transform_t *ct;
  sc_array_t          ctransforms, *cta;
#else
  bool                face_axis[3], contact_face_only, contact_edge_only;
  int                 ftransform[9], edge;
  size_t              etree;
  p8est_edge_info_t   ei;
  p8est_edge_transform_t *et;
  sc_array_t         *eta;
  p8est_corner_info_t ci;
  p8est_corner_transform_t *ct;
  sc_array_t         *cta;
#endif
  sc_array_t         *tquadrants;

  tree = sc_array_index (p4est->trees, qtree);
  P4EST_ASSERT (p4est_tree_is_complete (tree));
  tquadrants = &tree->quadrants;

  P4EST_QUADRANT_INIT (&fd);
  P4EST_QUADRANT_INIT (&ld);
  P4EST_QUADRANT_INIT (&tempq);
  for (which = 0; which < P4EST_INSUL; ++which) {
    P4EST_QUADRANT_INIT (&ins[which]);
  }
#ifndef P4_TO_P8
  cta = &ctransforms;
  sc_array_init (cta, sizeof (p4est_corner_transform_t));
#else
  eta = &ei.edge_transforms;
  sc_array_init (eta, sizeof (p8est_edge_transform_t));
  cta = &ci.corner_transforms;
  sc_array_init (cta, sizeof (p8est_corner_transform_t));
#endif

  /* assign some numbers */
  treecount = tquadrants->elem_count;
  P4EST_ASSERT (treecount > 0);
  incount = in->elem_count;

  /* return if there is nothing to do */
  if (treecount == 0 || incount == 0) {
    return;
  }

  /* retrieve first and last descendants in the tree */
  treefd = &tree->first_desc;
  treeld = &tree->last_desc;

  /* loop over input list of quadrants */
  for (iz = 0; iz < incount; ++iz) {
    inq = sc_array_index (in, iz);
    P4EST_ASSERT (0 <= inq->p.piggy2.which_tree);
    if (inq->p.piggy2.which_tree != qtree) {
      continue;
    }
    inter_tree = false;
    ntree = -1;
    face = corner = -1;
#ifndef P4_TO_P8
    transform = -1;
#else
    edge = -1;
    ftransform[0] = -1;
    ei.iedge = -1;
    et = NULL;
    ci.icorner = -1;
    ct = NULL;
    contact_face_only = contact_edge_only = false;
#endif
    if (!p4est_quadrant_is_inside_root (inq)) {
      /* this quadrant comes from a different tree */
      P4EST_ASSERT (p4est_quadrant_is_extended (inq));
      inter_tree = true;
#ifndef P4_TO_P8
      outface[0] = (inq->y < 0);
      outface[1] = (inq->x >= P4EST_ROOT_LEN);
      outface[2] = (inq->y >= P4EST_ROOT_LEN);
      outface[3] = (inq->x < 0);
      if ((outface[0] || outface[2]) && (outface[1] || outface[3])) {
        /* this quadrant is a corner neighbor */
        for (corner = 0; corner < 4; ++corner) {
          if (outface[(corner + 3) % 4] && outface[corner]) {
            break;
          }
        }
        p4est_find_corner_transform (conn, qtree, corner, cta);
        corner = p4est_corner_to_zorder[corner];
      }
      else {
        /* this quadrant is a face neighbor */
        for (face = 0; face < 4; ++face) {
          if (outface[face]) {
            break;
          }
        }
        P4EST_ASSERT (face < 4);
        ntree = conn->tree_to_tree[4 * qtree + face];
        transform = p4est_find_face_transform (conn, qtree, face);
      }
#else /* P4_TO_P8 */
      outface[0] = (inq->x < 0);
      outface[1] = (inq->x >= P4EST_ROOT_LEN);
      outface[2] = (inq->y < 0);
      outface[3] = (inq->y >= P4EST_ROOT_LEN);
      outface[4] = (inq->z < 0);
      outface[5] = (inq->z >= P4EST_ROOT_LEN);
      face_axis[0] = outface[0] || outface[1];
      face_axis[1] = outface[2] || outface[3];
      face_axis[2] = outface[4] || outface[5];
      if (!face_axis[1] && !face_axis[2]) {
        contact_face_only = true;
        face = 0 + outface[1];
      }
      else if (!face_axis[0] && !face_axis[2]) {
        contact_face_only = true;
        face = 2 + outface[3];
      }
      else if (!face_axis[0] && !face_axis[1]) {
        contact_face_only = true;
        face = 4 + outface[5];
      }
      else if (!face_axis[0]) {
        contact_edge_only = true;
        edge = 0 + 2 * outface[5] + outface[3];
      }
      else if (!face_axis[1]) {
        contact_edge_only = true;
        edge = 4 + 2 * outface[5] + outface[1];
      }
      else if (!face_axis[2]) {
        contact_edge_only = true;
        edge = 8 + 2 * outface[3] + outface[1];
      }
      if (contact_face_only) {
        P4EST_ASSERT (!contact_edge_only && face >= 0 && face < 6);
        P4EST_ASSERT (outface[face]);
        ntree = p8est_find_face_transform (conn, qtree, face, ftransform);
        P4EST_ASSERT (ntree >= 0);
      }
      else if (contact_edge_only) {
        P4EST_ASSERT (!contact_face_only && edge >= 0 && edge < 12);
        p8est_find_edge_transform (conn, qtree, edge, &ei);
        P4EST_ASSERT (ei.edge_transforms.elem_count > 0);
      }
      else {
        corner = 4 * outface[5] + 2 * outface[3] + outface[1];
        P4EST_ASSERT (p8est_quadrant_touches_corner (inq, corner, false));
        p8est_find_corner_transform (conn, qtree, corner, &ci);
        P4EST_ASSERT (ci.corner_transforms.elem_count > 0);
      }
#endif /* P4_TO_P8 */
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
            false) {
          /* this quadrant is outside this tree, no overlap */
          continue;
        }
        p4est_quadrant_first_descendent (s, &fd, P4EST_QMAXLEVEL);
        p4est_quadrant_last_descendent (s, &ld, P4EST_QMAXLEVEL);

        /* skip this insulation quadrant if there is no overlap */
        if (p4est_quadrant_compare (&ld, treefd) < 0 ||
            p4est_quadrant_compare (treeld, &fd) < 0) {
          continue;
        }

        /* find first quadrant in tree that fits between fd and ld */
        guess = treecount / 2;
        if (p4est_quadrant_compare (&fd, treefd) <= 0) {
          /* the first tree quadrant overlaps an insulation quadrant */
          first_index = 0;
        }
        else {
          /* do a binary search for the lowest tree quadrant >= s */
          first_index = p4est_find_lower_bound (tquadrants, s, guess);
          if (first_index < 0) {
            continue;
          }
          guess = (size_t) first_index;
        }

        /* find last quadrant in tree that fits between fd and ld */
        if (p4est_quadrant_compare (treeld, &ld) <= 0) {
          /* the last tree quadrant overlaps an insulation quadrant */
          last_index = (ssize_t) treecount - 1;
        }
        else {
          /* do a binary search for the highest tree quadrant <= ld */
          last_index = p4est_find_higher_bound (tquadrants, &ld, guess);
          if (last_index < 0) {
            SC_CHECK_NOT_REACHED ();
          }
        }

        /* skip if no overlap of sufficient level difference is found */
        if (first_index > last_index) {
          continue;
        }

        /* copy relevant quadrants into out */
        if (inter_tree && corner >= 0) {
          /* across a corner, find smallest quadrant to be sent */
          level = 0;
          for (js = first_index; js <= last_index; ++js) {
            tq = sc_array_index_ssize_t (tquadrants, js);
            if ((int) tq->level <= SC_MAX (level, (int) inq->level + 1)) {
              continue;
            }
            p4est_quadrant_shift_corner (tq, &tempq, corner);
            P4EST_ASSERT (p4est_quadrant_is_ancestor (s, &tempq));
            level = SC_MAX (level, (int) tempq.level);
          }
          if (level > 0) {
            /* send this small corner to all neighbor corner trees */
            for (ctree = 0; ctree < cta->elem_count; ++ctree) {
              outq = sc_array_push (out);
              outq->level = (int8_t) level;
              ct = sc_array_index (cta, ctree);
              p4est_quadrant_transform_corner (outq, (int) ct->ncorner,
                                               false);
              outq->p.piggy2.which_tree = ct->ntree;
            }
            ct = NULL;
          }
        }
        else {
          /* across face/edge or intra-tree, find quadrants that are small enough */
          P4EST_ASSERT (corner == -1);
          for (js = first_index; js <= last_index; ++js) {
            tq = sc_array_index_ssize_t (tquadrants, js);
            if (tq->level > inq->level + 1) {
              P4EST_ASSERT (p4est_quadrant_is_ancestor (s, tq));
              if (inter_tree) {
#ifndef P4_TO_P8
                outq = sc_array_push (out);
                tempq = *tq;
                p4est_quadrant_translate_face (&tempq, face);
                p4est_quadrant_transform_face (&tempq, outq, transform);
                outq->p.piggy2.which_tree = ntree;
#else
                if (contact_face_only) {
                  P4EST_ASSERT (!contact_edge_only);
                  outq = sc_array_push (out);
                  p8est_quadrant_transform_face (tq, outq, ftransform);
                  outq->p.piggy2.which_tree = ntree;
                }
                else {
                  P4EST_ASSERT (contact_edge_only);
                  p8est_quadrant_shift_edge (tq, &tempq, edge);
                  if (tempq.level > inq->level + 1) {
                    P4EST_ASSERT (p4est_quadrant_is_ancestor (s, &tempq));
                    for (etree = 0; etree < eta->elem_count; ++etree) {
                      outq = sc_array_push (out);
                      et = sc_array_index (eta, etree);
                      p8est_quadrant_transform_edge (&tempq, outq, &ei, et,
                                                     false);
                      outq->p.piggy2.which_tree = et->ntree;
                    }
                  }
                  et = NULL;
                }
#endif
              }
              else {
                outq = sc_array_push (out);
                *outq = *tq;
                outq->p.piggy2.which_tree = qtree;
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
}

void
p4est_tree_uniqify_overlap (sc_array_t * skip, sc_array_t * out)
{
  size_t              i, j;
  size_t              outcount, dupcount, skipcount;
  p4est_quadrant_t   *inq, *outq, *tq;

  outcount = out->elem_count;
  if (outcount == 0) {
    return;
  }

  /* sort array and remove duplicates */
  sc_array_sort (out, p4est_quadrant_compare);
  dupcount = skipcount = 0;
  i = 0;                        /* read counter */
  j = 0;                        /* write counter */
  inq = sc_array_index (out, i);
  while (i < outcount) {
    tq = ((i < outcount - 1) ? sc_array_index (out, i + 1) : NULL);
    if (i < outcount - 1 && p4est_quadrant_is_equal (inq, tq)) {
      ++dupcount;
      ++i;
    }
    else if (sc_array_bsearch (skip, inq, p4est_quadrant_compare_piggy) != -1) {
      ++skipcount;
      ++i;
    }
    else {
      if (i > j) {
        outq = sc_array_index (out, j);
        *outq = *inq;
      }
      ++i;
      ++j;
    }
    inq = tq;
  }
  P4EST_ASSERT (i == outcount);
  P4EST_ASSERT (j + dupcount + skipcount == outcount);
  sc_array_resize (out, j);
}

size_t
p4est_tree_remove_nonowned (p4est_t * p4est, p4est_topidx_t which_tree)
{
  bool                full_tree[2];
  size_t              zz, incount, prev_good, removed;
#ifdef P4EST_DEBUG
  const p4est_topidx_t first_tree = p4est->first_local_tree;
  const p4est_topidx_t last_tree = p4est->last_local_tree;
#endif
  const p4est_quadrant_t *first_pos, *next_pos;
  p4est_quadrant_t   *q1, *q2;
  p4est_quadrant_t    ld;
  p4est_tree_t       *tree;
  sc_array_t         *quadrants;

  P4EST_ASSERT (first_tree <= which_tree && which_tree <= last_tree);
  tree = p4est_array_index_topidx (p4est->trees, which_tree);
  P4EST_ASSERT (p4est_tree_is_almost_sorted (tree, false));

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
    q2 = sc_array_index (quadrants, zz);
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
         (p4est_quadrant_last_descendent (q2, &ld, P4EST_QMAXLEVEL),
          p4est_quadrant_compare (next_pos, &ld) <= 0))) {
      /* quadrant is outside of the unit square
         or at least partially outside of the tree bounds */
      --tree->quadrants_per_level[q2->level];
      p4est_quadrant_free_data (p4est, q2);
      ++removed;
#ifdef P4EST_DEBUG
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
      q1 = sc_array_index (quadrants, prev_good);
      if (zz > prev_good) {
        *q1 = *q2;
#ifdef P4EST_DEBUG
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
    q1 = sc_array_index (quadrants, 0);
  }
  sc_array_resize (quadrants, incount);

  tree->maxlevel = 0;
  for (zz = 0; zz < incount; ++zz) {
    q1 = sc_array_index (quadrants, zz);
    P4EST_ASSERT (p4est_quadrant_is_valid (q1));
    tree->maxlevel = (int8_t) SC_MAX (tree->maxlevel, q1->level);
  }

  P4EST_ASSERT (p4est_tree_is_sorted (tree));

  return removed;
}

void
p4est_complete_region (p4est_t * p4est,
                       const p4est_quadrant_t * q1,
                       bool include_q1,
                       const p4est_quadrant_t * q2,
                       bool include_q2,
                       p4est_tree_t * tree,
                       p4est_topidx_t which_tree, p4est_init_t init_fn)
{
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
  size_t              quadrant_pool_size;
  size_t              data_pool_size;
  int                 level, maxlevel = 0;
  p4est_locidx_t     *quadrants_per_level;
  p4est_locidx_t      num_quadrants = 0;

  P4EST_QUADRANT_INIT (&Afinest);

  W = sc_list_new (NULL);
  R = tree;

  /* needed for sanity check */
  quadrant_pool_size = p4est->quadrant_pool->elem_count;
  data_pool_size = 0;
  if (p4est->user_data_pool != NULL) {
    data_pool_size = p4est->user_data_pool->elem_count;
  }

  quadrants = &R->quadrants;
  quadrants_per_level = R->quadrants_per_level;

  /* Assert that we have an empty tree */
  P4EST_ASSERT (quadrants->elem_count == 0);

  comp = p4est_quadrant_compare (&a, &b);
  /* Assert that a<b */
  P4EST_ASSERT (comp < 0);

  /* R <- R + a */
  if (include_q1) {
    sc_array_resize (quadrants, 1);
    r = sc_array_index (quadrants, 0);
    *r = a;
    maxlevel = SC_MAX ((int) a.level, maxlevel);
    ++quadrants_per_level[a.level];
    ++num_quadrants;
  }

  if (comp < 0) {
    /* W <- C(A_{finest}(a,b)) */
    p4est_nearest_common_ancestor (&a, &b, &Afinest);

    c0 = sc_mempool_alloc (quadrant_pool);
    c1 = sc_mempool_alloc (quadrant_pool);
    c2 = sc_mempool_alloc (quadrant_pool);
    c3 = sc_mempool_alloc (quadrant_pool);
#ifdef P4_TO_P8
    c4 = sc_mempool_alloc (quadrant_pool);
    c5 = sc_mempool_alloc (quadrant_pool);
    c6 = sc_mempool_alloc (quadrant_pool);
    c7 = sc_mempool_alloc (quadrant_pool);

    p8est_quadrant_children (&Afinest, c0, c1, c2, c3, c4, c5, c6, c7);
#else
    p4est_quadrant_children (&Afinest, c0, c1, c2, c3);
#endif

    sc_list_append (W, c0);
    sc_list_append (W, c1);
    sc_list_append (W, c2);
    sc_list_append (W, c3);
#ifdef P4_TO_P8
    sc_list_append (W, c4);
    sc_list_append (W, c5);
    sc_list_append (W, c6);
    sc_list_append (W, c7);
#endif

    /* for each w in W */
    while (W->elem_count > 0) {
      w = sc_list_pop (W);
      level = (int) w->level;

      /* if (a < w < b) and (w not in {A(b)}) */
      if (((p4est_quadrant_compare (&a, w) < 0) &&
           (p4est_quadrant_compare (w, &b) < 0)
          ) && !p4est_quadrant_is_ancestor (w, &b)
        ) {
        /* R <- R + w */
        sc_array_resize (quadrants, num_quadrants + 1);
        r = sc_array_index (quadrants, num_quadrants);
        *r = *w;
        p4est_quadrant_init_data (p4est, which_tree, r, init_fn);
        maxlevel = SC_MAX (level, maxlevel);
        ++quadrants_per_level[level];
        ++num_quadrants;
      }
      /* else if (w in {{A(a)}, {A(b)}}) */
      else if (p4est_quadrant_is_ancestor (w, &a)
               || p4est_quadrant_is_ancestor (w, &b)) {
        /* W <- W + C(w) */
        c0 = sc_mempool_alloc (quadrant_pool);
        c1 = sc_mempool_alloc (quadrant_pool);
        c2 = sc_mempool_alloc (quadrant_pool);
        c3 = sc_mempool_alloc (quadrant_pool);
#ifdef P4_TO_P8
        c4 = sc_mempool_alloc (quadrant_pool);
        c5 = sc_mempool_alloc (quadrant_pool);
        c6 = sc_mempool_alloc (quadrant_pool);
        c7 = sc_mempool_alloc (quadrant_pool);

        p8est_quadrant_children (w, c0, c1, c2, c3, c4, c5, c6, c7);
#else
        p4est_quadrant_children (w, c0, c1, c2, c3);
#endif

#ifdef P4_TO_P8
        sc_list_prepend (W, c7);
        sc_list_prepend (W, c6);
        sc_list_prepend (W, c5);
        sc_list_prepend (W, c4);
#endif
        sc_list_prepend (W, c3);
        sc_list_prepend (W, c2);
        sc_list_prepend (W, c1);
        sc_list_prepend (W, c0);
      }

      /* W <- W - w */
      sc_mempool_free (quadrant_pool, w);
    }                           /* end for */

    /* R <- R + b */
    if (include_q2) {
      sc_array_resize (quadrants, num_quadrants + 1);
      r = sc_array_index (quadrants, num_quadrants);
      *r = b;
      maxlevel = SC_MAX ((int) b.level, maxlevel);
      ++quadrants_per_level[b.level];
      ++num_quadrants;
    }
  }

  R->maxlevel = (int8_t) maxlevel;

  P4EST_ASSERT (W->first == NULL && W->last == NULL);
  sc_list_destroy (W);

  P4EST_ASSERT (p4est_tree_is_complete (R));
  P4EST_ASSERT (quadrant_pool_size == p4est->quadrant_pool->elem_count);
  P4EST_ASSERT (num_quadrants == (p4est_locidx_t) quadrants->elem_count);
  if (p4est->user_data_pool != NULL) {
    P4EST_ASSERT (data_pool_size + quadrants->elem_count ==
                  p4est->user_data_pool->elem_count + (include_q1 ? 1 : 0)
                  + (include_q2 ? 1 : 0));
  }
}

/** Internal function to realize local completion / balancing.
 * \param [in] balance  can be 0: no balance only completion
 *                      and then in 2D:
 *                             1: balance across edges
 *                             2: balance across edges and corners
 *                      and in 3D:
 *                             1: balance across faces
 *                             2: balance across faces and edges
 *                             3: balance across faces, edges and corners
 */
static void
p4est_complete_or_balance (p4est_t * p4est, p4est_topidx_t which_tree,
                           p4est_init_t init_fn, int balance)
{
  bool                lookup, inserted;
  bool                isfamily, isoutroot;
#ifdef P4EST_BALANCE_OPTIMIZE
  bool                isintree;
#endif
  size_t              iz, jz;
  size_t              incount, ocount;
  size_t              quadrant_pool_size;
  size_t              data_pool_size;
  size_t              count_outside_root, count_outside_tree;
  size_t              count_already_inlist, count_already_outlist;
  size_t              count_moved1_outside, count_moved2_outside;
  size_t              num_added, num_nonowned, num_linearized;
  int                 qid, sid, pid;
  int                 bbound, fbound, rbound;
  int                 skey, *key = &skey;
  int                 pkey, *parent_key = &pkey;
  int                 l, inmaxl;
  void              **vlookup;
  ssize_t             srindex;
  p4est_qcoord_t      ph;
  p4est_quadrant_t   *family[P4EST_CHILDREN];
  p4est_quadrant_t   *q;
  p4est_quadrant_t   *qalloc, *qlookup, **qpointer;
  p4est_quadrant_t    ld, pshift;
  p4est_tree_t       *tree;
  sc_array_t         *inlist, *olist;
  sc_mempool_t       *list_alloc, *qpool;
  sc_hash_t          *hash[P4EST_MAXLEVEL + 1];
  sc_array_t          outlist[P4EST_MAXLEVEL + 1];
#ifdef P4_TO_P8
  int                 sindex;
#endif

  P4EST_ASSERT (which_tree >= p4est->first_local_tree);
  P4EST_ASSERT (which_tree <= p4est->last_local_tree);
  tree = p4est_array_index_topidx (p4est->trees, which_tree);

  P4EST_ASSERT (0 <= balance && balance <= P4EST_DIM);
  P4EST_ASSERT (p4est_tree_is_almost_sorted (tree, true));

  P4EST_QUADRANT_INIT (&ld);
  P4EST_QUADRANT_INIT (&pshift);
#ifdef P4_TO_P8
  P4EST_QUADRANT_INIT (&pshift);
#endif

  /*
   * Algorithm works with these data structures
   * inlist  --  sorted list of input quadrants
   * hash    --  hash table to hold additional quadrants not in inlist
   *             this is filled bottom-up to ensure balance condition
   * outlist --  filled simultaneously with hash, holding pointers
   *             don't rely on addresses of elements, it is resized
   * In the end, the elements of hash are appended to inlist
   * and inlist is sorted and linearized. This can be optimized later.
   */

  /* assign some shortcut variables */
#ifndef P4_TO_P8
  fbound = 8;
  bbound = ((balance == 0) ? 5 : fbound);
#else
  fbound = p8est_balance_count[P4EST_DIM];
  bbound = p8est_balance_count[balance];
#endif
  inlist = &tree->quadrants;
  incount = inlist->elem_count;
  inmaxl = (int) tree->maxlevel;
  qpool = p4est->quadrant_pool;

  /* needed for sanity check */
  quadrant_pool_size = qpool->elem_count;
  data_pool_size = 0;
  if (p4est->user_data_pool != NULL) {
    data_pool_size = p4est->user_data_pool->elem_count;
  }

  /* if tree is empty, there is nothing to do */
  if (incount == 0) {
    return;
  }

  /* initialize some counters */
  count_outside_root = count_outside_tree = 0;
  count_already_inlist = count_already_outlist = 0;
  count_moved1_outside = count_moved2_outside = 0;

  /* initialize temporary storage */
  list_alloc = sc_mempool_new (sizeof (sc_link_t));
  for (l = 0; l <= inmaxl; ++l) {
    hash[l] = sc_hash_new (p4est_quadrant_hash_fn, p4est_quadrant_equal_fn,
                           NULL, list_alloc);
    sc_array_init (&outlist[l], sizeof (p4est_quadrant_t *));
  }
  for (; l <= P4EST_MAXLEVEL; ++l) {
    hash[l] = NULL;
    memset (&outlist[l], -1, sizeof (sc_array_t));
  }

  /* walk through the input tree bottom-up */
  ph = 0;
  pid = -1;
  qalloc = sc_mempool_alloc (qpool);
  qalloc->p.user_data = key;
  for (l = inmaxl; l > 0; --l) {
    ocount = outlist[l].elem_count;     /* fix ocount here, it is growing */
    for (iz = 0; iz < incount + ocount; ++iz) {
      isfamily = false;
      if (iz < incount) {
        q = sc_array_index (inlist, iz);
        if ((int) q->level != l) {
          continue;
        }
        /* this is an optimization to catch adjacent siblings */
        if (iz + P4EST_CHILDREN <= incount) {
          family[0] = q;
          for (jz = 1; jz < P4EST_CHILDREN; ++jz) {
            family[jz] = sc_array_index (inlist, iz + jz);
          }
          if (p4est_quadrant_is_familypv (family)) {
            isfamily = true;
            iz += P4EST_CHILDREN - 1;   /* skip siblings */
          }
        }
      }
      else {
        qpointer = sc_array_index (&outlist[l], iz - incount);
        q = *qpointer;
        P4EST_ASSERT ((int) q->level == l);
      }
      P4EST_ASSERT (p4est_quadrant_is_extended (q));
      isoutroot = !p4est_quadrant_is_inside_root (q);
#ifdef P4EST_BALANCE_OPTIMIZE
      if (isoutroot) {
        isintree = false;
      }
      else {
        /* TODO: verify p4est_quadarant_is_inside_tree function */
        isintree = p4est_quadrant_is_inside_tree (tree, q);
        if (!isintree && p4est_quadrant_overlaps_tree (tree, q)) {
          ++count_moved1_outside;
          continue;
        }
      }
#endif
      /* TODO:
         For inter-tree quadrants always do full edge/corner balance
         May not be necessary and lead to too many quadrants */
      rbound = (isoutroot ? fbound : bbound);

      /*
       * check for q and its siblings,
       * then for q's parent and parent's indirect relevant neighbors
       * 2D
       * sid == 0..3    siblings including q
       *        4       parent of q
       *        5..7    relevant indirect neighbors of parent
       *                one of them is omitted if corner balance is off
       * 3D
       * sid == 0..7    siblings including q
       *        8       parent of q
       *        9..11   indirect face neighbors of parent
       *        12..14  indirect edge neighbors of parent
       *        15      indirect corner neighbor of parent
       *
       * if q is inside the tree, include all of the above.
       * if q is outside the tree, include only its parent and the neighbors.
       */
      qid = p4est_quadrant_child_id (q);        /* 0 <= qid < 4 resp. 6 */
      for (sid = 0; sid < rbound; ++sid) {
        /* stage 1: determine candidate qalloc */
        if (sid < P4EST_CHILDREN) {
          if (qid == sid || isfamily) {
            /* q (or its family) is included in inlist */
            continue;
          }
          if (isoutroot) {
            /* don't add siblings outside of the unit tree */
            continue;
          }
          p4est_quadrant_sibling (q, qalloc, sid);
        }
        else if (sid == P4EST_CHILDREN) {
          /* compute the parent */
          p4est_quadrant_parent (q, qalloc);
          if (balance > 0) {
            pshift = *qalloc;   /* copy parent for all balance cases */
            ph = P4EST_QUADRANT_LEN (pshift.level);     /* its size */
            pid = p4est_quadrant_child_id (&pshift);    /* and position */
#ifdef P4_TO_P8
            if (pid > 0 && pshift.level > 0)
              p4est_quadrant_sibling (&pshift, &pshift, 0);
#endif
          }
        }
        else {
          if (l == 1) {
            /* don't add tree-size quadrants as parent neighbors */
            break;
          }
#ifndef P4_TO_P8
          /* determine the 3 parent's relevant indirect neighbors */
          P4EST_ASSERT (sid >= 5 && sid < 8);
          if (balance < 2 && sid - 5 == corners_omitted[pid]) {
            /* this quadrant would only be needed for corner balance */
            continue;
          }
          qalloc->x = pshift.x + indirect_neighbors[pid][sid - 5][0] * ph;
          qalloc->y = pshift.y + indirect_neighbors[pid][sid - 5][1] * ph;
#else /* P4_TO_P8 */
          P4EST_ASSERT (sid >= p8est_balance_count[0]);
          if (sid < p8est_balance_count[1]) {
            /* face balance */
            sindex = p8est_corner_faces[pid][sid - p8est_balance_count[0]];
            P4EST_ASSERT (0 <= sindex && sindex < 6);
            qalloc->x = pshift.x + p8est_balance_coord[sindex][0] * ph;
            qalloc->y = pshift.y + p8est_balance_coord[sindex][1] * ph;
            qalloc->z = pshift.z + p8est_balance_coord[sindex][2] * ph;
          }
          else if (sid < p8est_balance_count[2]) {
            /* edge balance */
            sindex = p8est_corner_edges[pid][sid - p8est_balance_count[1]];
            P4EST_ASSERT (0 <= sindex && sindex < 12);
            qalloc->x = pshift.x + p8est_balance_coord[6 + sindex][0] * ph;
            qalloc->y = pshift.y + p8est_balance_coord[6 + sindex][1] * ph;
            qalloc->z = pshift.z + p8est_balance_coord[6 + sindex][2] * ph;
          }
          else {
            P4EST_ASSERT (sid == p8est_balance_count[3] - 1);
            /* corner balance */
            qalloc->x = pshift.x + p8est_balance_coord[18 + pid][0] * ph;
            qalloc->y = pshift.y + p8est_balance_coord[18 + pid][1] * ph;
            qalloc->z = pshift.z + p8est_balance_coord[18 + pid][2] * ph;
          }
#endif /* P4_TO_P8 */
          qalloc->level = pshift.level;

          /* TODO: Verify optimizations which may omit necessary quadrants */
          if (!isoutroot) {
            if (!p4est_quadrant_is_inside_root (qalloc)) {
              ++count_outside_root;
              continue;
            }
          }
          else {
            if (!p4est_quadrant_is_inside_3x3 (qalloc)) {
              ++count_outside_root;
              continue;
            }
#ifdef P4EST_BALANCE_OPTIMIZE
            if (!p4est_quadrant_is_inside_root (qalloc) &&
                (q->x / P4EST_ROOT_LEN != qalloc->x / P4EST_ROOT_LEN ||
                 q->y / P4EST_ROOT_LEN != qalloc->y / P4EST_ROOT_LEN ||
#ifdef P4_TO_P8
                 q->z / P4EST_ROOT_LEN != qalloc->z / P4EST_ROOT_LEN ||
#endif
                 false)) {
              ++count_outside_root;
              continue;
            }
#endif
          }
        }
        P4EST_ASSERT (p4est_quadrant_is_extended (qalloc));
        /*
           P4EST_DEBUGF ("Candidate level %d qxy 0x%x 0x%x at sid %d\n",
           qalloc->level, qalloc->x, qalloc->y, sid);
         */

        /* stage 2: include qalloc */
#if defined P4EST_BALANCE_WRONG && defined P4EST_BALANCE_OPTIMIZE
        if (isintree && p4est_quadrant_is_inside_root (qalloc) &&
            !p4est_quadrant_is_inside_tree (tree, qalloc)) {
          ++count_moved2_outside;
          continue;
        }
#endif
        /* make sure that qalloc is not included more than once */
        lookup = sc_hash_lookup (hash[qalloc->level], qalloc, &vlookup);
        if (lookup) {
          /* qalloc is already included in output list, this catches most */
          ++count_already_outlist;
          qlookup = *vlookup;
          if (sid == P4EST_CHILDREN && qlookup->p.user_data == parent_key) {
            break;              /* this parent has been triggered before */
          }
          continue;
        }
        srindex = sc_array_bsearch (inlist, qalloc, p4est_quadrant_compare);
        if (srindex != -1) {
          /* qalloc is included in inlist, this is more expensive to test */
          ++count_already_inlist;
          continue;
        }
        /* insert qalloc into the output list as well */
        if (sid == P4EST_CHILDREN) {
          qalloc->p.user_data = parent_key;
        }
        inserted = sc_hash_insert_unique (hash[qalloc->level], qalloc, NULL);
        P4EST_ASSERT (inserted);
        olist = &outlist[qalloc->level];
        qpointer = sc_array_push (olist);
        *qpointer = qalloc;
        /* we need a new quadrant now, the old one is stored away */
        qalloc = sc_mempool_alloc (qpool);
        qalloc->p.user_data = key;
      }
    }
  }
  sc_mempool_free (qpool, qalloc);

  /* merge outlist into input list and free temporary storage */
  P4EST_LDEBUGF ("Hash statistics for tree %lld\n", (long long) which_tree);
  num_added = 0;
  for (l = 0; l <= inmaxl; ++l) {
    /* print statistics and free hash tables */
#ifdef P4EST_DEBUG
    sc_hash_print_statistics (p4est_package_id, SC_LP_DEBUG, hash[l]);
#endif
    sc_hash_unlink_destroy (hash[l]);

    /* merge valid quadrants from outlist into inlist */
    ocount = outlist[l].elem_count;
    q = NULL;
    for (iz = 0; iz < ocount; ++iz) {
      /* go through output list */
      qpointer = sc_array_index (&outlist[l], iz);
      qalloc = *qpointer;
      P4EST_ASSERT ((int) qalloc->level == l);
      P4EST_ASSERT (qalloc->p.user_data == key ||
                    qalloc->p.user_data == parent_key);
      if (p4est_quadrant_is_inside_root (qalloc)) {
        /* copy temporary quadrant into final tree */
        q = sc_array_push (inlist);
        *q = *qalloc;
        ++num_added;
        ++tree->quadrants_per_level[l];

        /* complete quadrant initialization */
        p4est_quadrant_init_data (p4est, which_tree, q, init_fn);
      }
      else {
        P4EST_ASSERT (p4est_quadrant_is_extended (qalloc));
      }
      sc_mempool_free (qpool, qalloc);
    }
    if (q != NULL && l > (int) tree->maxlevel) {
      tree->maxlevel = (int8_t) l;
    }
    sc_array_reset (&outlist[l]);
  }
  sc_mempool_destroy (list_alloc);
  P4EST_ASSERT (incount + num_added == inlist->elem_count);

  /* print more statistics */
  P4EST_VERBOSEF ("Tree %lld Outside root %llu tree %llu\n",
                  (long long) which_tree,
                  (unsigned long long) count_outside_root,
                  (unsigned long long) count_outside_tree);
  P4EST_INFOF
    ("Tree %lld inlist %llu outlist %llu moved %llu %llu insert %llu\n",
     (long long) which_tree, (unsigned long long) count_already_inlist,
     (unsigned long long) count_moved1_outside,
     (unsigned long long) count_moved2_outside,
     (unsigned long long) count_already_outlist,
     (unsigned long long) num_added);

  /* sort and linearize tree */
  sc_array_sort (inlist, p4est_quadrant_compare);
  num_nonowned = p4est_tree_remove_nonowned (p4est, which_tree);
  num_linearized = p4est_linearize_tree (p4est, tree);

  /* run sanity checks */
  P4EST_ASSERT (quadrant_pool_size == qpool->elem_count);
  if (p4est->user_data_pool != NULL) {
    P4EST_ASSERT (data_pool_size + inlist->elem_count ==
                  p4est->user_data_pool->elem_count + incount);
  }
  P4EST_ASSERT (incount + num_added - num_nonowned - num_linearized ==
                tree->quadrants.elem_count);

  P4EST_ASSERT (p4est_tree_is_complete (tree));
}

void
p4est_complete_subtree (p4est_t * p4est,
                        p4est_topidx_t which_tree, p4est_init_t init_fn)
{
  p4est_complete_or_balance (p4est, which_tree, init_fn, 0);
}

void
p4est_balance_subtree (p4est_t * p4est, p4est_balance_type_t btype,
                       p4est_topidx_t which_tree, p4est_init_t init_fn)
{
  p4est_complete_or_balance (p4est, which_tree, init_fn,
                             p4est_balance_type_int (btype));
}

size_t
p4est_linearize_tree (p4est_t * p4est, p4est_tree_t * tree)
{
  size_t              data_pool_size;
  size_t              incount, removed;
  size_t              current, rest;
  p4est_locidx_t      num_quadrants;
  int                 i, maxlevel;
  p4est_quadrant_t   *q1, *q2;
  sc_array_t         *tquadrants = &tree->quadrants;

  P4EST_ASSERT (p4est_tree_is_sorted (tree));

  incount = tquadrants->elem_count;
  if (incount <= 1) {
    return 0;
  }
  data_pool_size = 0;
  if (p4est->user_data_pool != NULL) {
    data_pool_size = p4est->user_data_pool->elem_count;
  }
  removed = 0;

  /* run through the array and remove ancestors */
  current = 0;
  rest = current + 1;
  q1 = sc_array_index (tquadrants, current);
  while (rest < incount) {
    q2 = sc_array_index (tquadrants, rest);
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
        q1 = sc_array_index (tquadrants, current);
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

p4est_gloidx_t
p4est_partition_given (p4est_t * p4est,
                       const p4est_locidx_t * new_num_quadrants_in_proc)
{
  const int           num_procs = p4est->mpisize;
  const int           rank = p4est->mpirank;
  const p4est_gloidx_t *global_last_quad_index =
    p4est->global_last_quad_index;
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

  int                 i, sk;
  p4est_locidx_t      il, jl;
  p4est_topidx_t      it;
  p4est_topidx_t      which_tree;
  p4est_locidx_t      num_copy;
  p4est_topidx_t      first_tree, last_tree;
  int                 from_proc, to_proc;
  p4est_locidx_t      num_quadrants;
  int                 num_proc_recv_from, num_proc_send_to;
  p4est_topidx_t      num_recv_trees;
  p4est_topidx_t      new_first_local_tree, new_last_local_tree;
  p4est_locidx_t      new_local_num_quadrants;
  p4est_topidx_t      first_from_tree, last_from_tree, from_tree;
  p4est_locidx_t     *num_recv_from, *num_send_to;
  p4est_locidx_t     *new_local_tree_elem_count;
  p4est_locidx_t     *new_local_tree_elem_count_before;
  p4est_gloidx_t     *begin_send_to;
  p4est_gloidx_t      tree_from_begin, tree_from_end;
  p4est_gloidx_t      from_begin, from_end;
  p4est_gloidx_t      to_begin, to_end;
  p4est_gloidx_t      my_base, my_begin, my_end;
  p4est_gloidx_t     *new_global_last_quad_index;
  p4est_gloidx_t     *local_tree_last_quad_index;
  p4est_gloidx_t      diff64, total_quadrants_shipped;
  char              **recv_buf, **send_buf;
  sc_array_t         *quadrants;
  p4est_locidx_t     *num_per_tree_local;
  p4est_locidx_t     *num_per_tree_send_buf;
  p4est_locidx_t     *num_per_tree_recv_buf;
  p4est_quadrant_t   *quad_send_buf;
  p4est_quadrant_t   *quad_recv_buf;
  p4est_quadrant_t   *quad;
  p4est_tree_t       *tree;
  char               *user_data_send_buf;
  char               *user_data_recv_buf;
  size_t              recv_size, send_size;
#ifdef P4EST_MPI
  int                 mpiret;
  MPI_Comm            comm = p4est->mpicomm;
  MPI_Request        *recv_request, *send_request;
  MPI_Status         *recv_status, *send_status;
#endif
#ifdef P4EST_DEBUG
  unsigned            crc;
  p4est_gloidx_t      total_requested_quadrants = 0;
#endif

  P4EST_GLOBAL_INFOF
    ("Into " P4EST_STRING "_partition_given with %lld total quadrants\n",
     (long long) p4est->global_num_quadrants);

#ifdef P4EST_DEBUG
  /* Save a checksum of the original forest */
  crc = p4est_checksum (p4est);

  /* Check for a valid requested partition */
  for (i = 0; i < num_procs; ++i) {
    total_requested_quadrants += new_num_quadrants_in_proc[i];
    P4EST_ASSERT (new_num_quadrants_in_proc[i] >= 0);
  }
  P4EST_ASSERT (total_requested_quadrants == p4est->global_num_quadrants);
#endif

  /* Print some diagnostics */
  if (rank == 0) {
    for (i = 0; i < num_procs; ++i) {
      P4EST_GLOBAL_VERBOSEF ("partition global_last_quad_index[%d] = %lld\n",
                             i, (long long) global_last_quad_index[i]);
    }
  }

  /* Calculate the global_last_quad_index for the repartitioned forest */
  new_global_last_quad_index = P4EST_ALLOC (p4est_gloidx_t, num_procs);
  new_global_last_quad_index[0] = new_num_quadrants_in_proc[0] - 1;
  for (i = 1; i < num_procs; ++i) {
    new_global_last_quad_index[i] = new_num_quadrants_in_proc[i] +
      new_global_last_quad_index[i - 1];
  }
  P4EST_ASSERT (global_last_quad_index[num_procs - 1] ==
                new_global_last_quad_index[num_procs - 1]);

  /* Calculate the global number of shipped quadrants */
  total_quadrants_shipped = 0;
  for (i = 1; i < num_procs; ++i) {
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
  P4EST_ASSERT (0 <= total_quadrants_shipped &&
                total_quadrants_shipped <= p4est->global_num_quadrants);

  /* Print some diagnostics */
  if (rank == 0) {
    for (i = 0; i < num_procs; ++i) {
      P4EST_GLOBAL_VERBOSEF
        ("partition new_global_last_quad_index[%d] = %lld\n",
         i, (long long) new_global_last_quad_index[i]);
    }
  }

  /* Calculate the local index of the end of each tree */
  local_tree_last_quad_index =
    P4EST_ALLOC_ZERO (p4est_gloidx_t, trees->elem_count);
  if (first_local_tree >= 0) {
    tree = sc_array_index (p4est->trees, first_local_tree);
    local_tree_last_quad_index[first_local_tree]
      = tree->quadrants.elem_count - 1;
  }
  else {
    P4EST_ASSERT (first_local_tree == -1 && last_local_tree == -2);
  }
  for (which_tree = first_local_tree + 1;       /* same type */
       which_tree <= last_local_tree; ++which_tree) {
    tree = sc_array_index (p4est->trees, which_tree);
    local_tree_last_quad_index[which_tree] = tree->quadrants.elem_count
      + local_tree_last_quad_index[which_tree - 1];
  }

#ifdef P4EST_DEBUG
  for (which_tree = first_local_tree; which_tree <= last_local_tree;
       ++which_tree) {
    tree = sc_array_index (p4est->trees, which_tree);
    P4EST_LDEBUGF
      ("partition tree %lld local_tree_last_quad_index[%lld] = %lld\n",
       (long long) which_tree, (long long) which_tree,
       (long long) local_tree_last_quad_index[which_tree]);
  }
#endif

  /* Calculate where the quadrants are coming from */
  num_recv_from = P4EST_ALLOC (p4est_locidx_t, num_procs);

  my_begin = (rank == 0) ? 0 : (new_global_last_quad_index[rank - 1] + 1);
  my_end = new_global_last_quad_index[rank];

  num_proc_recv_from = 0;
  for (from_proc = 0; from_proc < num_procs; ++from_proc) {
    from_begin = (from_proc == 0) ?
      0 : (global_last_quad_index[from_proc - 1] + 1);
    from_end = global_last_quad_index[from_proc];

    if (from_begin <= my_end && from_end >= my_begin) {
      /* from_proc sends to me */
      num_recv_from[from_proc] = SC_MIN (my_end, from_end)
        - SC_MAX (my_begin, from_begin) + 1;
      if (from_proc != rank)
        ++num_proc_recv_from;
    }
    else {
      /* from_proc does not send to me */
      num_recv_from[from_proc] = 0;
    }
  }

#ifdef P4EST_DEBUG
  for (i = 0; i < num_procs; ++i) {
    if (num_recv_from[i] != 0) {
      P4EST_LDEBUGF ("partition num_recv_from[%d] = %lld\n", i,
                     (long long) num_recv_from[i]);
    }
  }
#endif

  /* Post receives for the quadrants and their data */
  recv_buf = P4EST_ALLOC (char *, num_procs);
#ifdef P4EST_MPI
  recv_request = P4EST_ALLOC (MPI_Request, num_proc_recv_from);
  recv_status = P4EST_ALLOC (MPI_Status, num_proc_recv_from);
#endif

  /* Allocate space for receiving quadrants and user data */
  for (from_proc = 0, sk = 0; from_proc < num_procs; ++from_proc) {
    if (from_proc != rank && num_recv_from[from_proc]) {
      num_recv_trees =          /* same type */
        p4est->global_first_position[from_proc + 1].p.which_tree
        - p4est->global_first_position[from_proc].p.which_tree + 1;
      recv_size = num_recv_trees * sizeof (p4est_locidx_t)
        + quad_plus_data_size * num_recv_from[from_proc];

      recv_buf[from_proc] = P4EST_ALLOC (char, recv_size);

      /* Post receives for the quadrants and their data */
#ifdef P4EST_MPI
      P4EST_LDEBUGF ("partition recv %lld quadrants from %d\n",
                     (long long) num_recv_from[from_proc], from_proc);
      mpiret = MPI_Irecv (recv_buf[from_proc], (int) recv_size, MPI_BYTE,
                          from_proc, P4EST_COMM_PARTITION_GIVEN,
                          comm, recv_request + sk);
      SC_CHECK_MPI (mpiret);
#endif
      ++sk;
    }
    else {
      recv_buf[from_proc] = NULL;
    }
  }
#ifdef P4EST_MPI
  for (; sk < num_proc_recv_from; ++sk) {
    recv_request[sk] = MPI_REQUEST_NULL;
  }
#endif

  /* For each processor calculate the number of quadrants sent */
  num_send_to = P4EST_ALLOC (p4est_locidx_t, num_procs);
  begin_send_to = P4EST_ALLOC (p4est_gloidx_t, num_procs);

  my_begin = (rank == 0) ? 0 : (global_last_quad_index[rank - 1] + 1);
  my_end = global_last_quad_index[rank];

  num_proc_send_to = 0;
  for (to_proc = 0; to_proc < num_procs; ++to_proc) {
    to_begin = (to_proc == 0)
      ? 0 : (new_global_last_quad_index[to_proc - 1] + 1);
    to_end = new_global_last_quad_index[to_proc];

    if (to_begin <= my_end && to_end >= my_begin) {
      /* I send to to_proc */
      num_send_to[to_proc] = SC_MIN (my_end, to_end)
        - SC_MAX (my_begin, to_begin) + 1;
      begin_send_to[to_proc] = SC_MAX (my_begin, to_begin);
      if (to_proc != rank)
        ++num_proc_send_to;
    }
    else {
      /* I don't send to to_proc */
      num_send_to[to_proc] = 0;
      begin_send_to[to_proc] = -1;
    }

  }

#ifdef P4EST_DEBUG
  for (i = 0; i < num_procs; ++i) {
    if (num_send_to[i] != 0) {
      P4EST_LDEBUGF ("partition num_send_to[%d] = %lld\n",
                     i, (long long) num_send_to[i]);
    }
  }
  for (i = 0; i < num_procs; ++i) {
    if (begin_send_to[i] >= 0) {
      P4EST_LDEBUGF ("partition begin_send_to[%d] = %lld\n",
                     i, (long long) begin_send_to[i]);
    }
  }
#endif

  /* Communicate the quadrants and their data */
  send_buf = P4EST_ALLOC (char *, num_procs);
#ifdef P4EST_MPI
  send_request = P4EST_ALLOC (MPI_Request, num_proc_send_to);
  send_status = P4EST_ALLOC (MPI_Status, num_proc_send_to);
#endif

  /* Set the num_per_tree_local */
  num_per_tree_local = P4EST_ALLOC_ZERO (p4est_locidx_t, num_send_trees);
  to_proc = rank;
  my_base = (rank == 0) ? 0 : (global_last_quad_index[rank - 1] + 1);
  my_begin = begin_send_to[to_proc] - my_base;
  my_end = begin_send_to[to_proc] + num_send_to[to_proc] - 1 - my_base;
  for (which_tree = first_local_tree; which_tree <= last_local_tree;
       ++which_tree) {
    tree = sc_array_index (p4est->trees, which_tree);

    from_begin = (which_tree == first_local_tree) ? 0 :
      (local_tree_last_quad_index[which_tree - 1] + 1);
    from_end = local_tree_last_quad_index[which_tree];

    if (from_begin <= my_end && from_end >= my_begin) {
      /* Need to copy from tree which_tree */
      tree_from_begin = SC_MAX (my_begin, from_begin) - from_begin;
      tree_from_end = SC_MIN (my_end, from_end) - from_begin;
      num_copy = tree_from_end - tree_from_begin + 1;

      num_per_tree_local[which_tree - first_local_tree] = num_copy;
    }
  }

  /* Allocate space for receiving quadrants and user data */
  for (to_proc = 0, sk = 0; to_proc < num_procs; ++to_proc) {
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
        tree = sc_array_index (p4est->trees, which_tree);

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
            quad_send_buf[il].p.user_data = NULL;

          }

          /* move the pointer to the begining of the quads that need copied */
          my_begin += num_copy;
          quad_send_buf += num_copy;
          user_data_send_buf += num_copy * data_size;
        }
      }

      /* Post receives for the quadrants and their data */
#ifdef P4EST_MPI
      P4EST_LDEBUGF ("partition send %lld quadrants to %d\n",
                     (long long) num_send_to[to_proc], to_proc);
      mpiret = MPI_Isend (send_buf[to_proc], (int) send_size, MPI_BYTE,
                          to_proc, P4EST_COMM_PARTITION_GIVEN,
                          comm, send_request + sk);
      SC_CHECK_MPI (mpiret);
      ++sk;
#endif
    }
    else {
      send_buf[to_proc] = NULL;
    }
  }
#ifdef P4EST_MPI
  for (; sk < num_proc_send_to; ++sk) {
    send_request[sk] = MPI_REQUEST_NULL;
  }

  /* Fill in forest */
  mpiret = MPI_Waitall (num_proc_recv_from, recv_request, recv_status);
  SC_CHECK_MPI (mpiret);
#endif

  /* Loop Through and fill in */

  /* Calculate the local index of the end of each tree in the repartition */
  new_local_tree_elem_count =
    P4EST_ALLOC_ZERO (p4est_locidx_t, trees->elem_count);
  new_local_tree_elem_count_before =
    P4EST_ALLOC_ZERO (p4est_locidx_t, trees->elem_count);
  new_first_local_tree = (p4est_topidx_t) P4EST_TOPIDX_MAX;
  new_last_local_tree = 0;

  for (from_proc = 0; from_proc < num_procs; ++from_proc) {
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
      num_per_tree_recv_buf =
        (from_proc ==
         rank) ? num_per_tree_local : (p4est_locidx_t *) recv_buf[from_proc];

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
  P4EST_INFOF ("partition new forest [%lld,%lld]\n",
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
  my_begin = begin_send_to[rank] - my_base;
  my_end = begin_send_to[rank] + num_send_to[rank] - 1 - my_base;

  for (which_tree = first_tree; which_tree <= last_tree; ++which_tree) {
    tree = sc_array_index (p4est->trees, which_tree);
    quadrants = &tree->quadrants;

    if (new_local_tree_elem_count[which_tree] > 0) {
      if (which_tree >= first_local_tree && which_tree <= last_local_tree) {

        num_quadrants = new_local_tree_elem_count[which_tree];

        from_begin = (which_tree == first_local_tree) ? 0 :
          (local_tree_last_quad_index[which_tree - 1] + 1);
        from_end = local_tree_last_quad_index[which_tree];

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
        for (il = 0; il < (p4est_locidx_t) quadrants->elem_count; ++il) {
          if (il < tree_from_begin || il > tree_from_end) {
            quad = sc_array_index (quadrants, il);
            p4est_quadrant_free_data (p4est, quad);
          }
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
        for (il = 0; il < (p4est_locidx_t) quadrants->elem_count; ++il) {
          quad = sc_array_index (quadrants, il);
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
  for (from_proc = 0; from_proc < num_procs; ++from_proc) {
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
          tree = sc_array_index (p4est->trees, from_tree);
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
          for (jl = 0; jl < num_copy; ++jl) {
            quad = sc_array_index (quadrants,
                                   jl +
                                   new_local_tree_elem_count_before
                                   [from_tree]);

            if (data_size > 0) {
              quad->p.user_data = sc_mempool_alloc (p4est->user_data_pool);
              memcpy (quad->p.user_data, user_data_recv_buf + jl * data_size,
                      data_size);
            }
            else {
              quad->p.user_data = NULL;
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
      if (recv_buf[from_proc] != NULL) {
        P4EST_FREE (recv_buf[from_proc]);
        recv_buf[from_proc] = NULL;
      }
    }
  }

  /* Set the global index and count of quadrants instead
   * of calling p4est_comm_count_quadrants
   */
  P4EST_FREE (p4est->global_last_quad_index);
  p4est->global_last_quad_index = new_global_last_quad_index;
  P4EST_ASSERT (p4est->global_num_quadrants ==
                new_global_last_quad_index[num_procs - 1] + 1);

  p4est->first_local_tree = new_first_local_tree;
  p4est->last_local_tree = new_last_local_tree;

  new_local_num_quadrants = 0;
  for (which_tree = 0; which_tree < new_first_local_tree; ++which_tree) {
    tree = sc_array_index (p4est->trees, which_tree);
    tree->quadrants_offset = 0;
    P4EST_QUADRANT_INIT (&tree->first_desc);
    P4EST_QUADRANT_INIT (&tree->last_desc);
  }
  for (; which_tree <= new_last_local_tree; ++which_tree) {
    tree = sc_array_index (p4est->trees, which_tree);
    tree->quadrants_offset = new_local_num_quadrants;
    quadrants = &tree->quadrants;
    P4EST_ASSERT (quadrants->elem_count > 0);

    new_local_num_quadrants +=  /* same type */
      (p4est_locidx_t) quadrants->elem_count;

    for (i = 0; i <= P4EST_QMAXLEVEL; ++i) {
      tree->quadrants_per_level[i] = 0;
    }
    tree->maxlevel = 0;
    for (il = 0; il < (p4est_locidx_t) quadrants->elem_count; ++il) {
      quad = sc_array_index (quadrants, il);
      ++tree->quadrants_per_level[quad->level];
      tree->maxlevel = (int8_t) SC_MAX (quad->level, tree->maxlevel);
    }

    quad = sc_array_index (quadrants, 0);
    p4est_quadrant_first_descendent (quad, &tree->first_desc,
                                     P4EST_QMAXLEVEL);
    quad = sc_array_index (quadrants, quadrants->elem_count - 1);
    p4est_quadrant_last_descendent (quad, &tree->last_desc, P4EST_QMAXLEVEL);
  }
  for (; which_tree < p4est->connectivity->num_trees; ++which_tree) {
    tree = sc_array_index (p4est->trees, which_tree);
    tree->quadrants_offset = new_local_num_quadrants;
    P4EST_QUADRANT_INIT (&tree->first_desc);
    P4EST_QUADRANT_INIT (&tree->last_desc);
  }
  p4est->local_num_quadrants = new_local_num_quadrants;

  /* Clean up */

#ifdef P4EST_MPI
  mpiret = MPI_Waitall (num_proc_send_to, send_request, send_status);
  SC_CHECK_MPI (mpiret);

#ifdef P4EST_DEBUG
  for (i = 0; i < num_proc_recv_from; ++i) {
    P4EST_ASSERT (recv_request[i] == MPI_REQUEST_NULL);
  }
  for (i = 0; i < num_proc_send_to; ++i) {
    P4EST_ASSERT (send_request[i] == MPI_REQUEST_NULL);
  }
#endif
  P4EST_FREE (recv_request);
  P4EST_FREE (recv_status);
  P4EST_FREE (send_request);
  P4EST_FREE (send_status);
#endif

  for (i = 0; i < num_procs; ++i) {
    if (recv_buf[i] != NULL)
      P4EST_FREE (recv_buf[i]);
    if (send_buf[i] != NULL)
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

  p4est_comm_global_partition (p4est);

  /* Assert that we have a valid partition */
  P4EST_ASSERT (crc == p4est_checksum (p4est));
  P4EST_GLOBAL_INFOF
    ("Done " P4EST_STRING
     "_partition_given shipped %lld quadrants %.3g%%\n",
     (long long) total_quadrants_shipped,
     total_quadrants_shipped * 100. / p4est->global_num_quadrants);

  return total_quadrants_shipped;
}

/* EOF p4est_algorithms.c */
