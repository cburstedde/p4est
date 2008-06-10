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

#ifndef P4_TO_P8

/* *INDENT-OFF* */

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
{0, 1, 1, 2};

/* *INDENT-ON* */

#endif /* !P4_TO_P8 */

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
  if (init_fn != NULL && p4est_quadrant_is_inside (quad)) {
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

void
p4est_quadrant_print (int log_priority, const p4est_quadrant_t * q)
{
#ifdef P4_TO_P8
  P4EST_NORMAL_LOGF (log_priority,
                     "x 0x%x y 0x%x z 0x%x level %d\n",
                     q->x, q->y, q->z, q->level);
#else
  P4EST_NORMAL_LOGF (log_priority,
                     "x 0x%x y 0x%x level %d\n", q->x, q->y, q->level);
#endif
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

#ifndef P4_TO_P8

bool
p4est_tree_is_almost_sorted (p4est_tree_t * tree, bool check_linearity)
{
  size_t              iz;
  int                 face_contact1;
  int                 face_contact2;
  p4est_quadrant_t   *q1, *q2;
  sc_array_t         *tquadrants = &tree->quadrants;

  if (tquadrants->elem_count <= 1) {
    return true;
  }

  q1 = sc_array_index (tquadrants, 0);
  face_contact1 = 0;
  face_contact1 |= ((q1->y < 0) ? 0x01 : 0);
  face_contact1 |= ((q1->x >= P4EST_ROOT_LEN) ? 0x02 : 0);
  face_contact1 |= ((q1->y >= P4EST_ROOT_LEN) ? 0x04 : 0);
  face_contact1 |= ((q1->x < 0) ? 0x08 : 0);
  for (iz = 1; iz < tquadrants->elem_count; ++iz) {
    q2 = sc_array_index (tquadrants, iz);
    face_contact2 = 0;
    face_contact2 |= ((q2->y < 0) ? 0x01 : 0);
    face_contact2 |= ((q2->x >= P4EST_ROOT_LEN) ? 0x02 : 0);
    face_contact2 |= ((q2->y >= P4EST_ROOT_LEN) ? 0x04 : 0);
    face_contact2 |= ((q2->x < 0) ? 0x08 : 0);

    if ((face_contact1 & 0x05) && (face_contact1 & 0x0a) &&
        face_contact1 == face_contact2) {
      /* both quadrants are outside the same corner and may overlap */
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

#endif /* !P4_TO_P8 */

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
  const int           num_procs = p4est->mpisize;
  const int           rank = p4est->mpirank;
  const p4est_topidx_t first_tree = p4est->first_local_tree;
  const p4est_topidx_t last_tree = p4est->last_local_tree;
  int                 i, maxlevel;
  size_t              jz;
  p4est_topidx_t      next_tree;
  p4est_locidx_t      nquadrants, lquadrants, perlevel;
  p4est_quadrant_t   *q;
  p4est_quadrant_t    mylow, nextlow, s;
  p4est_tree_t       *tree;

  P4EST_QUADRANT_INIT (&mylow);
  P4EST_QUADRANT_INIT (&nextlow);
  P4EST_QUADRANT_INIT (&s);

  /* check last item of global partition */
  P4EST_ASSERT (p4est->global_first_position[num_procs].which_tree ==
                p4est->connectivity->num_trees &&
                p4est->global_first_position[num_procs].x == 0 &&
                p4est->global_first_position[num_procs].y == 0 &&
#ifdef P4_TO_P8
                p4est->global_first_position[num_procs].z == 0 &&
#endif
                true);
  P4EST_ASSERT (p4est->connectivity->num_trees ==
                (p4est_topidx_t) p4est->trees->elem_count);

  /* check first tree in global partition */
  if (first_tree < 0) {
    if (!(first_tree == -1 && last_tree == -2)) {
      P4EST_INFO ("p4est invalid empty tree range A");
      return false;
    }
  }
  else {
    if (p4est->global_first_position[rank].which_tree != first_tree) {
      P4EST_INFO ("p4est invalid first tree\n");
      return false;
    }
    mylow.x = p4est->global_first_position[rank].x;
    mylow.y = p4est->global_first_position[rank].y;
#ifdef P4_TO_P8
    mylow.z = p4est->global_first_position[rank].z;
#endif
    mylow.level = P4EST_MAXLEVEL;
    tree = sc_array_index (p4est->trees, first_tree);
    if (tree->quadrants.elem_count > 0) {
      q = sc_array_index (&tree->quadrants, 0);
      if (q->x != mylow.x || q->y != mylow.y ||
#ifdef P4_TO_P8
          q->z != mylow.z ||
#endif
          false) {
        P4EST_INFO ("p4est invalid low quadrant\n");
        return false;
      }
    }
  }

  /* check last tree in global partition */
  if (last_tree < 0) {
    if (!(first_tree == -1 && last_tree == -2)) {
      P4EST_INFO ("p4est invalid empty tree range B");
      return false;
    }
  }
  else {
    next_tree = p4est->global_first_position[rank + 1].which_tree;
    if (next_tree != last_tree && next_tree != last_tree + 1) {
      P4EST_INFO ("p4est invalid last tree\n");
      return false;
    }
    nextlow.x = p4est->global_first_position[rank + 1].x;
    nextlow.y = p4est->global_first_position[rank + 1].y;
#ifdef P4_TO_P8
    nextlow.z = p4est->global_first_position[rank + 1].z;
#endif
    nextlow.level = P4EST_MAXLEVEL;
    tree = sc_array_index (p4est->trees, last_tree);
    if (tree->quadrants.elem_count > 0) {
      q = sc_array_index (&tree->quadrants, tree->quadrants.elem_count - 1);
      if (next_tree == last_tree) {
        if (!p4est_quadrant_is_next (q, &nextlow)) {
          P4EST_INFO ("p4est invalid next quadrant\n");
          return false;
        }
      }
      else {
        p4est_quadrant_last_descendent (q, &s, P4EST_MAXLEVEL);
        if (s.x + 1 != P4EST_ROOT_LEN || s.y + 1 != P4EST_ROOT_LEN ||
#ifdef P4_TO_P8
            s.z + 1 != P4EST_ROOT_LEN ||
#endif
            false) {
          P4EST_INFO ("p4est invalid last quadrant\n");
          return false;
        }
      }
    }
  }

  /* check individual trees */
  lquadrants = 0;
  for (jz = 0; jz < p4est->trees->elem_count; ++jz) {
    tree = sc_array_index (p4est->trees, jz);
    if (!p4est_tree_is_complete (tree)) {
      P4EST_INFO ("p4est invalid not complete\n");
      return false;
    }
    if (((p4est_topidx_t) jz < p4est->first_local_tree ||
         (p4est_topidx_t) jz > p4est->last_local_tree) &&
        tree->quadrants.elem_count > 0) {
      P4EST_INFO ("p4est invalid outside count\n");
      return false;
    }

    maxlevel = 0;
    nquadrants = 0;
    for (i = 0; i <= P4EST_MAXLEVEL; ++i) {
      perlevel = tree->quadrants_per_level[i];

      P4EST_ASSERT (perlevel >= 0);
      nquadrants += perlevel;   /* same type */
      if (perlevel > 0) {
        maxlevel = i;
      }
    }
    lquadrants += nquadrants;   /* same type */

    if (maxlevel != (int) tree->maxlevel) {
      P4EST_INFO ("p4est invalid wrong maxlevel\n");
      return false;
    }
    if (nquadrants != (p4est_locidx_t) tree->quadrants.elem_count) {
      P4EST_INFO ("p4est invalid tree quadrant count\n");
      return false;
    }
  }

  if (lquadrants != p4est->local_num_quadrants) {
    P4EST_INFO ("p4est invalid local quadrant count\n");
    return false;
  }

  return true;
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

#ifndef P4_TO_P8

void
p4est_tree_compute_overlap (p4est_t * p4est, p4est_topidx_t qtree,
                            sc_array_t * in, sc_array_t * out)
{
  int                 k, l, which;
  size_t              iz;
  size_t              ctree;
  size_t              treecount, incount, outcount;
  size_t              guess;
  ssize_t             first_index, last_index, js;
  bool                inter_tree;
  int                 transform, outface[4];
  int                 face, corner, zcorner = -1, level;
  p4est_topidx_t      ntree;
  p4est_qcoord_t      qh;
  sc_array_t          corner_info;
  p4est_quadrant_t    treefd, treeld;
  p4est_quadrant_t    fd, ld, tempq, ins[9], cq;
  p4est_quadrant_t   *tq, *s;
  p4est_quadrant_t   *inq, *outq;
  p4est_tree_t       *tree;
  p4est_connectivity_t *conn = p4est->connectivity;
  p4est_corner_info_t *ci;
  sc_array_t         *tquadrants;

  tree = sc_array_index (p4est->trees, qtree);
  P4EST_ASSERT (p4est_tree_is_complete (tree));
  tquadrants = &tree->quadrants;

  P4EST_QUADRANT_INIT (&treefd);
  P4EST_QUADRANT_INIT (&treeld);
  P4EST_QUADRANT_INIT (&fd);
  P4EST_QUADRANT_INIT (&ld);
  P4EST_QUADRANT_INIT (&tempq);
  P4EST_QUADRANT_INIT (&cq);
  for (which = 0; which < 9; ++which) {
    P4EST_QUADRANT_INIT (&ins[which]);
  }
  sc_array_init (&corner_info, sizeof (p4est_corner_info_t));

  /* assign some numbers */
  treecount = tquadrants->elem_count;
  P4EST_ASSERT (treecount > 0);
  incount = in->elem_count;
  outcount = out->elem_count;

  /* return if there is nothing to do */
  if (treecount == 0 || incount == 0) {
    return;
  }

  /* compute first and last descendants in the tree */
  tq = sc_array_index (tquadrants, 0);
  p4est_quadrant_first_descendent (tq, &treefd, P4EST_MAXLEVEL);
  tq = sc_array_index (tquadrants, treecount - 1);
  p4est_quadrant_last_descendent (tq, &treeld, P4EST_MAXLEVEL);

  /* loop over input list of quadrants */
  for (iz = 0; iz < incount; ++iz) {
    inq = sc_array_index (in, iz);
    if (inq->p.piggy.which_tree != qtree) {
      continue;
    }
    inter_tree = false;
    ntree = qtree;
    face = corner = -1;
    transform = -1;
    if (!p4est_quadrant_is_inside (inq)) {
      /* this quadrant comes from a different tree */
      P4EST_ASSERT (p4est_quadrant_is_extended (inq));
      inter_tree = true;
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
        p4est_find_corner_info (conn, qtree, corner, &corner_info);

        /* construct highest corner quadrant */
        cq.level = P4EST_MAXLEVEL;
        zcorner = p4est_corner_to_zorder[corner];
        p4est_quadrant_corner (&cq, zcorner, 1);
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
    }
    qh = P4EST_QUADRANT_LEN (inq->level);

    /* loop over the insulation layer of inq */
    for (k = 0; k < 3; ++k) {
      for (l = 0; l < 3; ++l) {
        which = k * 3 + l;      /* 0..8 */

        /* exclude myself from the queries */
        if (which == 4) {
          continue;
        }
        s = &ins[which];
        *s = *inq;
        s->x += (l - 1) * qh;
        s->y += (k - 1) * qh;
        if ((s->x < 0 || s->x >= P4EST_ROOT_LEN) ||
            (s->y < 0 || s->y >= P4EST_ROOT_LEN)) {
          /* this quadrant is outside this tree, no overlap */
          continue;
        }
        p4est_quadrant_first_descendent (s, &fd, P4EST_MAXLEVEL);
        p4est_quadrant_last_descendent (s, &ld, P4EST_MAXLEVEL);

        /* skip this insulation quadrant if there is no overlap */
        if (p4est_quadrant_compare (&ld, &treefd) < 0 ||
            p4est_quadrant_compare (&treeld, &fd) < 0) {
          continue;
        }

        /* find first quadrant in tree that fits between fd and ld */
        guess = treecount / 2;
        if (p4est_quadrant_compare (&fd, &treefd) <= 0) {
          /* the first tree quadrant is contained in insulation quadrant */
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
        if (p4est_quadrant_compare (&treeld, &ld) <= 0) {
          /* the last tree quadrant is contained in insulation layer */
          last_index = (ssize_t) treecount - 1;
        }
        else {
          /* do a binary search for the highest tree quadrant <= ld */
          last_index = p4est_find_higher_bound (tquadrants, &ld, guess);
          if (last_index < 0) {
            P4EST_VERBOSE ("Last index < 0\n");
            continue;
          }
        }

        /* skip if no overlap of sufficient level difference is found */
        if (first_index > last_index) {
          continue;
        }

        /* copy relevant quadrants into out */
        if (inter_tree && corner >= 0) {
          /* across the corner, find smallest corner quadrant to be sent */
          level = 0;
          for (js = first_index; js <= last_index; ++js) {
            tq = sc_array_index_ssize_t (tquadrants, js);
            if ((int) tq->level <= level) {
              continue;
            }
            level = p4est_quadrant_corner_level (tq, zcorner, level);
          }
          zcorner = -1;         /* will be recycled below */

          /* send this small corner to all neighbor corner trees */
          for (ctree = 0; ctree < corner_info.elem_count; ++ctree) {
            ci = sc_array_index (&corner_info, ctree);

            sc_array_resize (out, outcount + 1);
            outq = sc_array_index (out, outcount);
            outq->level = (int8_t) level;
            zcorner = p4est_corner_to_zorder[ci->ncorner];
            p4est_quadrant_corner (outq, zcorner, 0);
            outq->p.piggy.which_tree = ci->ntree;
            ++outcount;
          }
        }
        else {
          /* across face or intra-tree, find quadrants that are small enough */
          for (js = first_index; js <= last_index; ++js) {
            tq = sc_array_index_ssize_t (tquadrants, js);
            if (tq->level > inq->level + 1) {
              sc_array_resize (out, outcount + 1);
              outq = sc_array_index (out, outcount);
              if (inter_tree) {
                tempq = *tq;
                p4est_quadrant_translate (&tempq, face);
                p4est_quadrant_transform (&tempq, outq, transform);
              }
              else {
                *outq = *tq;
              }
              outq->p.piggy.which_tree = ntree;
              ++outcount;
            }
          }
        }
      }
    }
  }

  sc_array_reset (&corner_info);
}

#endif /* !P4_TO_P8 */

void
p4est_tree_uniqify_overlap (sc_array_t * not, sc_array_t * out)
{
  size_t              i, j;
  size_t              outcount, dupcount, notcount;
  p4est_quadrant_t   *inq, *outq, *tq;

  outcount = out->elem_count;
  if (outcount == 0) {
    return;
  }

  /* sort array and remove duplicates */
  sc_array_sort (out, p4est_quadrant_compare);
  dupcount = notcount = 0;
  i = 0;                        /* read counter */
  j = 0;                        /* write counter */
  inq = sc_array_index (out, i);
  while (i < outcount) {
    tq = ((i < outcount - 1) ? sc_array_index (out, i + 1) : NULL);
    if (i < outcount - 1 && p4est_quadrant_is_equal (inq, tq)) {
      ++dupcount;
      ++i;
    }
    else if (sc_array_bsearch (not, inq, p4est_quadrant_compare_piggy) != -1) {
      ++notcount;
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
  P4EST_ASSERT (j + dupcount + notcount == outcount);
  sc_array_resize (out, j);
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

#ifndef P4_TO_P8

/** Internal function to realize local completion / balancing.
 * \param [in] balance  can be 0: no balancing
 *                             1: balance across edges
 *                             2: balance across edges and corners
 */
static void
p4est_complete_or_balance (p4est_t * p4est, p4est_tree_t * tree, int balance,
                           p4est_topidx_t which_tree, p4est_init_t init_fn)
{
  size_t              iz, jz;
  size_t              incount, curcount, ocount;
  int                 comp;
  bool                lookup, inserted;
  bool                isfamily, isoutroot;
  size_t              quadrant_pool_size;
  size_t              data_pool_size;
  size_t              count_outside_root, count_outside_tree;
  size_t              count_already_inlist, count_already_outlist;
  size_t              first_inside, last_inside;
  int                 qid, sid, pid, bbound;
  int                 skey, *key = &skey;
  int                 pkey, *parent_key = &pkey;
  int                 outface[4];
  int                 l, inmaxl;
  void               *vlookup;
  ssize_t             srindex;
  p4est_qcoord_t      ph;
  p4est_quadrant_t   *family[4];
  p4est_quadrant_t   *q;
  p4est_quadrant_t   *qalloc, *qlookup, **qpointer;
  p4est_quadrant_t    ld, tree_first, tree_last, parent;
  sc_array_t         *inlist, *olist;
  sc_mempool_t       *list_alloc, *qpool;
  sc_hash_t          *hash[P4EST_MAXLEVEL + 1];
  sc_array_t          outlist[P4EST_MAXLEVEL + 1];

  P4EST_ASSERT (p4est_tree_is_almost_sorted (tree, 1));

  P4EST_QUADRANT_INIT (&ld);
  P4EST_QUADRANT_INIT (&tree_first);
  P4EST_QUADRANT_INIT (&tree_last);
  P4EST_QUADRANT_INIT (&parent);

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
  bbound = ((balance == 0) ? 5 : 8);
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

  /* if tree is empty or a single block, there is nothing to do */
  if (incount <= 1) {
    return;
  }

  /* determine the first and last small quadrants contained in the tree */
  first_inside = incount;
  q = NULL;
  for (iz = 0; iz < incount; ++iz) {
    q = sc_array_index (inlist, iz);
    if (p4est_quadrant_is_inside (q)) {
      first_inside = iz;
      p4est_quadrant_first_descendent (q, &tree_first, inmaxl);
      break;
    }
  }
  if (iz == incount) {
    /* only extended quadrants in the tree, there is nothing to do */
    return;
  }
  last_inside = incount - 1;
  p4est_quadrant_last_descendent (q, &tree_last, inmaxl);
  for (iz = first_inside + 1; iz < incount; ++iz) {
    q = sc_array_index (inlist, iz);
    if (!p4est_quadrant_is_inside (q)) {
      last_inside = iz - 1;
      break;
    }
    p4est_quadrant_last_descendent (q, &ld, inmaxl);
    comp = p4est_quadrant_compare (&tree_last, &ld);
    if (comp < 0) {
      tree_last = ld;
    }
  }
  P4EST_ASSERT (first_inside <= last_inside && last_inside < incount);
  P4EST_ASSERT (p4est_quadrant_is_valid (&tree_first));
  P4EST_ASSERT (p4est_quadrant_is_valid (&tree_last));

  /* initialize some counters */
  count_outside_root = count_outside_tree = 0;
  count_already_inlist = count_already_outlist = 0;

  /* initialize temporary storage */
  list_alloc = sc_mempool_new (sizeof (sc_link_t));
  for (l = 0; l <= inmaxl; ++l) {
    hash[l] = sc_hash_new (p4est_quadrant_hash, p4est_quadrant_is_equal,
                           list_alloc);
    sc_array_init (&outlist[l], sizeof (p4est_quadrant_t *));
  }
  for (l = inmaxl + 1; l <= P4EST_MAXLEVEL; ++l) {
    hash[l] = NULL;
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
        if (iz + 4 <= incount) {
          family[0] = q;
          for (jz = 1; jz < 4; ++jz) {
            family[jz] = sc_array_index (inlist, iz + jz);
          }
          if (p4est_quadrant_is_family (family[0], family[1],
                                        family[2], family[3])) {
            isfamily = true;
            iz += 3;            /* skip siblings */
          }
        }
      }
      else {
        qpointer = sc_array_index (&outlist[l], iz - incount);
        q = *qpointer;
        P4EST_ASSERT ((int) q->level == l);
      }
      P4EST_ASSERT (p4est_quadrant_is_extended (q));
      isoutroot = !p4est_quadrant_is_inside (q);

      /*
       * check for q and its siblings,
       * then for q's parent and parent's indirect relevant neighbors
       * sid == 0..3  siblings including q
       *        4     parent of q
       *        5..7  relevant indirect neighbors of parent
       *              one of them is omitted if corner balance is off
       *
       * if q is inside the tree, include all of the above.
       * if q is outside the tree, include only its parent and the neighbors.
       */
      qid = p4est_quadrant_child_id (q);        /* 0 <= qid < 4 */
      for (sid = 0; sid < bbound; ++sid) {
        /* stage 1: determine candidate qalloc */
        if (sid < 4) {
          if (qid == sid || isfamily) {
            /* q (or its family) is included in inlist */
            continue;
          }
          if (isoutroot) {
            /* q is outside the tree */
            continue;
          }
          p4est_quadrant_sibling (q, qalloc, sid);
        }
        else if (sid == 4) {
          /* compute the parent */
          p4est_quadrant_parent (q, qalloc);
          if (bbound > 5) {
            parent = *qalloc;   /* copy parent for cases 5..7 */
            ph = P4EST_QUADRANT_LEN (parent.level);     /* its size */
            pid = p4est_quadrant_child_id (&parent);    /* and position */
          }
        }
        else {
          /* determine the 3 parent's relevant indirect neighbors */
          P4EST_ASSERT (sid >= 5 && sid < 8);
          if (balance < 2 && sid - 5 == corners_omitted[pid]) {
            /* this quadrant would only be needed for corner balance */
            continue;
          }
          qalloc->x = parent.x + indirect_neighbors[pid][sid - 5][0] * ph;
          qalloc->y = parent.y + indirect_neighbors[pid][sid - 5][1] * ph;
          qalloc->level = parent.level;
          outface[0] = (qalloc->y < 0);
          outface[1] = (qalloc->x >= P4EST_ROOT_LEN);
          outface[2] = (qalloc->y >= P4EST_ROOT_LEN);
          outface[3] = (qalloc->x < 0);
          if (!isoutroot) {
            if (outface[0] || outface[1] || outface[2] || outface[3]) {
              /* q is inside and this quadrant is outside the root */
              ++count_outside_root;
              continue;
            }
          }
          else {
            if ((outface[0] || outface[2]) && (outface[1] || outface[3])) {
              /* quadrant is outside and across the corner */
              ++count_outside_root;
              continue;
            }
          }
        }
        /*
           P4EST_DEBUGF ("Candidate level %d qxy 0x%x 0x%x at sid %d\n",
           qalloc->level, qalloc->x, qalloc->y, sid);
         */

        /* stage 2: include qalloc if necessary */
        if (p4est_quadrant_is_inside (qalloc)) {
          p4est_quadrant_last_descendent (qalloc, &ld, inmaxl);
          if ((p4est_quadrant_compare (&tree_first, qalloc) > 0 &&
               (qalloc->x != tree_first.x || qalloc->y != tree_first.y)) ||
              p4est_quadrant_compare (&ld, &tree_last) > 0) {
            /* qalloc is outside the tree */
            ++count_outside_tree;
            continue;
          }
        }
        lookup = sc_hash_lookup (hash[qalloc->level], qalloc, &vlookup);
        if (lookup) {
          /* qalloc is already included in output list, this catches most */
          ++count_already_outlist;
          qlookup = vlookup;
          if (sid == 4 && qlookup->p.user_data == parent_key) {
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
        if (sid == 4) {
          qalloc->p.user_data = parent_key;
        }
        inserted = sc_hash_insert_unique (hash[qalloc->level], qalloc, NULL);
        P4EST_ASSERT (inserted);
        olist = &outlist[qalloc->level];
        sc_array_resize (olist, olist->elem_count + 1);
        qpointer = sc_array_index (olist, olist->elem_count - 1);
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
  curcount = inlist->elem_count;
  for (l = 0; l <= inmaxl; ++l) {
    /* print statistics and free hash tables */
#ifdef P4EST_DEBUG
    sc_hash_print_statistics (SC_LP_DEBUG, hash[l]);
#endif /* P4EST_DEBUG */
    sc_hash_unlink_destroy (hash[l]);   /* performance optimization */

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
      if (p4est_quadrant_is_inside (qalloc)) {
        /* copy temporary quadrant into final tree */
        sc_array_resize (inlist, curcount + 1);
        q = sc_array_index (inlist, curcount);
        *q = *qalloc;
        ++curcount;
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
#ifdef P4EST_DEBUG
  sc_mempool_reset (list_alloc);
#endif
  sc_mempool_destroy (list_alloc);
  P4EST_ASSERT (curcount >= incount);

  /* print more statistics */
  P4EST_VERBOSEF ("Tree %lld Outside root %llu tree %llu\n",
                  (long long) which_tree,
                  (unsigned long long) count_outside_root,
                  (unsigned long long) count_outside_tree);
  P4EST_INFOF
    ("Tree %lld Already in inlist %llu outlist %llu insertions %llu\n",
     (long long) which_tree, (unsigned long long) count_already_inlist,
     (unsigned long long) count_already_outlist,
     (unsigned long long) curcount - incount);

  /* sort and linearize tree */
  sc_array_sort (inlist, p4est_quadrant_compare);
  p4est_linearize_subtree (p4est, tree);

  /* run sanity checks */
  P4EST_ASSERT (quadrant_pool_size == qpool->elem_count);
  if (p4est->user_data_pool != NULL) {
    P4EST_ASSERT (data_pool_size + inlist->elem_count ==
                  p4est->user_data_pool->elem_count + incount);
  }
  P4EST_ASSERT (p4est_tree_is_linear (tree));
}

void
p4est_complete_subtree (p4est_t * p4est, p4est_tree_t * tree,
                        p4est_topidx_t which_tree, p4est_init_t init_fn)
{
  p4est_complete_or_balance (p4est, tree, 0, which_tree, init_fn);
}

void
p4est_balance_subtree (p4est_t * p4est, p4est_tree_t * tree,
                       p4est_topidx_t which_tree, p4est_init_t init_fn)
{
  p4est_complete_or_balance (p4est, tree, 2, which_tree, init_fn);
}

#endif /* !P4_TO_P8 */

void
p4est_linearize_subtree (p4est_t * p4est, p4est_tree_t * tree)
{
  size_t              data_pool_size;
  size_t              incount, removed;
  size_t              current, rest;
  p4est_locidx_t      num_quadrants;
  int                 i, maxlevel;
  p4est_quadrant_t   *q1, *q2;
  sc_array_t         *tquadrants = &tree->quadrants;

#ifdef P4_TO_P8
  /* TODO: need some variant of is_almost_sorted */
  P4EST_ASSERT (p4est_tree_is_sorted (tree));
#else
  P4EST_ASSERT (p4est_tree_is_almost_sorted (tree, 0));
#endif

  incount = tquadrants->elem_count;
  if (incount <= 1) {
    return;
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
  for (i = 0; i <= P4EST_MAXLEVEL; ++i) {
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
    p4est->global_first_position[rank + 1].which_tree - /* same type */
    p4est->global_first_position[rank].which_tree + 1;
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
    ("Into p4est_partition_given with %lld total quadrants\n",
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
        p4est->global_first_position[from_proc + 1].which_tree
        - p4est->global_first_position[from_proc].which_tree + 1;
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

  /* Fill in forest */
#ifdef P4EST_MPI
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
      first_from_tree = p4est->global_first_position[from_proc].which_tree;
      last_from_tree = p4est->global_first_position[from_proc + 1].which_tree;
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
        sc_array_reset (quadrants);
        for (i = 0; i <= P4EST_MAXLEVEL; ++i) {
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
      first_from_tree = p4est->global_first_position[from_proc].which_tree;
      last_from_tree = p4est->global_first_position[from_proc + 1].which_tree;
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
  for (which_tree = new_first_local_tree; which_tree <= new_last_local_tree;
       ++which_tree) {
    tree = sc_array_index (p4est->trees, which_tree);
    quadrants = &tree->quadrants;

    new_local_num_quadrants +=  /* same type */
      (p4est_locidx_t) quadrants->elem_count;

    for (i = 0; i <= P4EST_MAXLEVEL; ++i) {
      tree->quadrants_per_level[i] = 0;
    }
    tree->maxlevel = 0;
    for (il = 0; il < (p4est_locidx_t) quadrants->elem_count; ++il) {
      quad = sc_array_index (quadrants, il);
      ++tree->quadrants_per_level[quad->level];
      tree->maxlevel = (int8_t) SC_MAX (quad->level, tree->maxlevel);
    }
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
    ("Done p4est_partition_given shipped %lld quadrants %.3g%%\n",
     (long long) total_quadrants_shipped,
     total_quadrants_shipped * 100. / p4est->global_num_quadrants);

  return total_quadrants_shipped;
}

/* EOF p4est_algorithms.c */
