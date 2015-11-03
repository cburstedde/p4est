/*
  This file is part of p4est.
  p4est is a C library to manage a collection (a forest) of multiple
  connected adaptive quadtrees or octrees in parallel.

  Copyright (C) 2014 The University of Texas System
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

#include <p6est_profile.h>
#include <p4est_bits.h>

/* given two profiles (layers that have been reduced to just their levels),
 * take the union, i.e. combine them, taking the finer layers */
static void
p6est_profile_union (sc_array_t * a, sc_array_t * b, sc_array_t * c)
{
  size_t              az, bz, na;
  P4EST_ASSERT (SC_ARRAY_IS_OWNER (c));
  P4EST_ASSERT (a->elem_size == sizeof (int8_t));
  P4EST_ASSERT (b->elem_size == sizeof (int8_t));
  P4EST_ASSERT (c->elem_size == sizeof (int8_t));
  int8_t              al, bl, finel, *cc;
  p4est_qcoord_t      finesize, coarsesize;
  sc_array_t         *finer;
  size_t             *fineincr;

  sc_array_truncate (c);
  az = 0;
  bz = 0;
  na = a->elem_count;
  while (az < na) {
    P4EST_ASSERT (bz < b->elem_count);

    cc = (int8_t *) sc_array_push (c);

    al = *((int8_t *) sc_array_index (a, az++));
    bl = *((int8_t *) sc_array_index (b, bz++));
    if (al == bl) {
      *cc = al;
      continue;
    }
    else if (al > bl) {
      finer = a;
      finesize = P4EST_QUADRANT_LEN (al);
      fineincr = &az;
      finel = al;
      coarsesize = P4EST_QUADRANT_LEN (bl);
    }
    else {
      finer = b;
      finesize = P4EST_QUADRANT_LEN (bl);
      fineincr = &bz;
      finel = bl;
      coarsesize = P4EST_QUADRANT_LEN (al);
    }

    P4EST_ASSERT (finesize < coarsesize);

    do {
      *cc = finel;
      cc = (int8_t *) sc_array_push (c);
      finel = *((int8_t *) sc_array_index (finer, (*fineincr)++));
      finesize += P4EST_QUADRANT_LEN (finel);
    } while (finesize < coarsesize);
    P4EST_ASSERT (finesize == coarsesize);
    *cc = finel;
  }
}

/* given two profiles (layers that have been reduced to just their levels),
 * take the intersection, i.e. combine them, taking the coarser layers */
static void
p6est_profile_intersection (sc_array_t * a, sc_array_t * b, sc_array_t * c)
{
  size_t              az, bz, na;
  P4EST_ASSERT (SC_ARRAY_IS_OWNER (c));
  P4EST_ASSERT (a->elem_size == sizeof (int8_t));
  P4EST_ASSERT (b->elem_size == sizeof (int8_t));
  P4EST_ASSERT (c->elem_size == sizeof (int8_t));
  int8_t              al, bl, finel, *cc;
  p4est_qcoord_t      finesize, coarsesize;
  sc_array_t         *finer;
  size_t             *fineincr;

  sc_array_truncate (c);
  az = 0;
  bz = 0;
  na = a->elem_count;
  while (az < na) {
    P4EST_ASSERT (bz < b->elem_count);

    cc = (int8_t *) sc_array_push (c);

    al = *((int8_t *) sc_array_index (a, az++));
    bl = *((int8_t *) sc_array_index (b, bz++));
    if (al == bl) {
      *cc = al;
      continue;
    }
    else if (al > bl) {
      *cc = bl;
      finer = a;
      finesize = P4EST_QUADRANT_LEN (al);
      fineincr = &az;
      finel = al;
      coarsesize = P4EST_QUADRANT_LEN (bl);
    }
    else {
      *cc = al;
      finer = b;
      finesize = P4EST_QUADRANT_LEN (bl);
      fineincr = &bz;
      finel = bl;
      coarsesize = P4EST_QUADRANT_LEN (al);
    }

    P4EST_ASSERT (finesize < coarsesize);

    do {
      finel = *((int8_t *) sc_array_index (finer, (*fineincr)++));
      finesize += P4EST_QUADRANT_LEN (finel);
    } while (finesize < coarsesize);
    P4EST_ASSERT (finesize == coarsesize);
  }
}

static void
p6est_profile_balance_self_one_pass (sc_array_t * read, sc_array_t * write)
{
  int                 stackcount;
  int8_t              n, newn, p, l;
  int8_t             *wc;
  size_t              count = read->elem_count;
  size_t              zy;

  P4EST_ASSERT (SC_ARRAY_IS_OWNER (write));
  P4EST_ASSERT (read->elem_size == sizeof (int8_t));
  P4EST_ASSERT (write->elem_size == sizeof (int8_t));

  sc_array_truncate (write);
  wc = (int8_t *) sc_array_push (write);
  n = *((int8_t *) sc_array_index (read, count - 1));
  *wc = l = n;
  for (zy = 1; zy < count; zy++) {
    n = *((int8_t *) sc_array_index (read, count - 1 - zy));
    p = l - 1;
    newn = SC_MAX (p, n);
    stackcount = newn - n;
    wc = (int8_t *) sc_array_push_count (write, 1 + stackcount);
    *wc = l = newn;
    while (stackcount--) {
      *(++wc) = l = newn--;
    }
  }
}

static void
p6est_profile_balance_self (sc_array_t * a, sc_array_t * work)
{
  P4EST_ASSERT (SC_ARRAY_IS_OWNER (a));
  P4EST_ASSERT (SC_ARRAY_IS_OWNER (work));
  P4EST_ASSERT (a->elem_size == sizeof (int8_t));
  P4EST_ASSERT (work->elem_size == sizeof (int8_t));

  p6est_profile_balance_self_one_pass (a, work);
  p6est_profile_balance_self_one_pass (work, a);
}

static void
p6est_profile_balance_face_one_pass (sc_array_t * read, sc_array_t * write,
                                     p4est_qcoord_t readh)
{
  int8_t             *wc;
  size_t              count;
  int                 stackcount;
  int8_t              n, nn, newn, p, l;
  size_t              zy;

  P4EST_ASSERT (SC_ARRAY_IS_OWNER (write));
  P4EST_ASSERT (read->elem_size == sizeof (int8_t));
  P4EST_ASSERT (write->elem_size == sizeof (int8_t));

  count = read->elem_count;

  sc_array_truncate (write);
  l = 0;
  zy = 0;
  while (zy < count) {
    n = *((int8_t *) sc_array_index (read, count - 1 - zy++));
    if (n && !(readh & P4EST_QUADRANT_LEN (n))) {
      P4EST_ASSERT (zy < count);
      nn = *((int8_t *) sc_array_index (read, count - 1 - zy));
      if (n == nn) {
        zy++;
        n--;
      }
    }
    readh += P4EST_QUADRANT_LEN (n);
    p = l - 1;
    newn = SC_MAX (p, n);
    stackcount = newn - n;
    wc = (int8_t *) sc_array_push_count (write, 1 + stackcount);
    *wc = l = newn;
    while (stackcount--) {
      *(++wc) = l = newn--;
    }
  }
}

/* assumes a is already self balanced */
static void
p6est_profile_balance_face (sc_array_t * a, sc_array_t * b, sc_array_t * work,
                            p4est_qcoord_t diff)
{
  P4EST_ASSERT (SC_ARRAY_IS_OWNER (b));
  P4EST_ASSERT (SC_ARRAY_IS_OWNER (work));
  P4EST_ASSERT (a->elem_size == sizeof (int8_t));
  P4EST_ASSERT (b->elem_size == sizeof (int8_t));
  P4EST_ASSERT (work->elem_size == sizeof (int8_t));

  p6est_profile_balance_face_one_pass (a, work, diff);
  p6est_profile_balance_self_one_pass (work, b);
}

static void
p6est_profile_balance_full_one_pass (sc_array_t * read, sc_array_t * write,
                                     p4est_qcoord_t readh)
{
  int8_t             *wc;
  size_t              count;
  int                 stackcount;
  int8_t              n, nn, newn, p, l, prevl, nextl;
  size_t              zy;

  P4EST_ASSERT (SC_ARRAY_IS_OWNER (write));
  P4EST_ASSERT (read->elem_size == sizeof (int8_t));
  P4EST_ASSERT (write->elem_size == sizeof (int8_t));

  count = read->elem_count;

  sc_array_truncate (write);
  l = 0;
  zy = 0;
  while (zy < count) {
    n = *((int8_t *) sc_array_index (read, count - 1 - zy++));
    if (n && !(readh & P4EST_QUADRANT_LEN (n))) {
      P4EST_ASSERT (zy < count);
      nn = *((int8_t *) sc_array_index (read, count - 1 - zy));
      if (n == nn) {
        if (zy > 1) {
          prevl = *((int8_t *) sc_array_index (read, count - 1 - (zy - 2)));
        }
        else {
          prevl = -1;
        }
        if (zy < count - 1) {
          nextl = *((int8_t *) sc_array_index (read, count - 1 - (zy + 1)));
        }
        else {
          nextl = -1;
        }
        if (n >= SC_MAX (nextl, prevl) - 1) {
          zy++;
          n--;
        }
      }
    }
    readh += P4EST_QUADRANT_LEN (n);
    p = l - 1;
    newn = SC_MAX (p, n);
    stackcount = newn - n;
    wc = (int8_t *) sc_array_push_count (write, 1 + stackcount);
    *wc = l = newn;
    while (stackcount--) {
      *(++wc) = l = newn--;
    }
  }
}

/* assumes a is already self balanced */
static void
p6est_profile_balance_full (sc_array_t * a, sc_array_t * b, sc_array_t * work,
                            p4est_qcoord_t diff)
{
  P4EST_ASSERT (SC_ARRAY_IS_OWNER (b));
  P4EST_ASSERT (SC_ARRAY_IS_OWNER (work));
  P4EST_ASSERT (a->elem_size == sizeof (int8_t));
  P4EST_ASSERT (b->elem_size == sizeof (int8_t));
  P4EST_ASSERT (work->elem_size == sizeof (int8_t));

  p6est_profile_balance_full_one_pass (a, work, diff);
  p6est_profile_balance_self_one_pass (work, b);
}

static void
p6est_profile_element_to_node_single (sc_array_t * elem, sc_array_t * node,
                                      int degree, p4est_locidx_t offset,
                                      p4est_locidx_t ** elem_to_node,
                                      p6est_lnodes_code_t * fc, int fcoffset)
{
  size_t              nedge = node->elem_count;
  size_t              az, bz;
  int                 i;

  P4EST_ASSERT (degree > 1);

  az = 0;

  for (bz = 0; bz < nedge; bz++) {
    int8_t              a;
    int8_t              b = *((int8_t *) sc_array_index (node, bz));
    int                 loop = 0;

    do {
      a = *((int8_t *) sc_array_index (elem, az));
      P4EST_ASSERT (a == b || a == b + 1);
      loop = !loop && (a == b + 1);
      for (i = 0; i < degree + 1; i++) {
        elem_to_node[az][i] = offset + bz * degree + i;
      }
      if (fc && a == b + 1) {
        fc[az] |= (1 << (fcoffset + 5));
      }
      az++;
    } while (loop);
  }
  P4EST_ASSERT (az == elem->elem_count);
}

static void
p6est_profile_element_to_node_col (p6est_profile_t * profile,
                                   p4est_locidx_t cid,
                                   p4est_locidx_t * offsets,
                                   p4est_locidx_t * e_to_n,
                                   p6est_lnodes_code_t * fc)
{
  p4est_locidx_t (*lr)[2] = (p4est_locidx_t (*)[2]) profile->lnode_ranges;
  p4est_locidx_t      nelem;
  p4est_locidx_t    **elem_to_node;
  int                 i, j, k;
  p4est_locidx_t      ll;
  sc_array_t          elem, node;
  sc_array_t         *lc = profile->lnode_columns;
  p4est_locidx_t      ncid, nid;
  p4est_lnodes_code_t fc4 = profile->lnodes->face_code[cid];
  p4est_locidx_t     *en = profile->lnodes->element_nodes;
  int                 degree = profile->lnodes->degree;
  int                 Nrp = degree + 1;
  int                 Nfp = (degree + 1) * (degree + 1);

  P4EST_ASSERT (degree > 1);

  ncid = en[Nfp * cid + Nrp * (Nrp / 2) + (Nrp / 2)];
  nelem = lr[ncid][1];

  sc_array_init_view (&elem, lc, lr[ncid][0], nelem);

  elem_to_node = P4EST_ALLOC (p4est_locidx_t *, nelem);

  for (ll = 0; ll < nelem; ll++) {
    fc[ll] = (p6est_lnodes_code_t) fc4;
  }
  for (k = 0, j = 0; j < Nrp; j++) {
    for (i = 0; i < Nrp; i++, k++) {
      nid = en[Nfp * cid + k];
      sc_array_init_view (&node, lc, lr[nid][0], lr[nid][1]);
      for (ll = 0; ll < nelem; ll++) {
        elem_to_node[ll] = e_to_n +
          (degree + 1) * (degree + 1) * (degree + 1) * ll + (degree + 1) * k;
      }
      if (!(i % degree) && !(j % degree)) {
        int                 c = 2 * (! !j) + (! !i);

        p6est_profile_element_to_node_single (&elem, &node, degree,
                                              offsets[nid], elem_to_node, fc,
                                              4 + c);
      }
      else if ((i % degree) && (j % degree)) {
        p6est_profile_element_to_node_single (&elem, &elem, degree,
                                              offsets[nid], elem_to_node,
                                              NULL, -1);
      }
      else {
        int                 f = 2 * !(j % degree) + (i == degree
                                                     || j == degree);

        p6est_profile_element_to_node_single (&elem, &node, degree,
                                              offsets[nid], elem_to_node, fc,
                                              f);
      }
    }
  }
  P4EST_FREE (elem_to_node);
}

void
p6est_profile_element_to_node (p6est_t * p6est,
                               p6est_profile_t * profile,
                               p4est_locidx_t * offsets,
                               p4est_locidx_t * elem_to_node,
                               p6est_lnodes_code_t * fc)
{
  p4est_topidx_t      jt;
  p4est_t            *columns = p6est->columns;
  p4est_tree_t       *tree;
  p4est_quadrant_t   *col;
  sc_array_t         *tquadrants;
  p4est_locidx_t (*lr)[2] = (p4est_locidx_t (*)[2]) profile->lnode_ranges;
  p4est_locidx_t      cid;
  size_t              zz;
  p6est_lnodes_code_t mask = 0x1fe0;
  p6est_lnodes_code_t hbit = 0x0010;
  int                 degree = profile->lnodes->degree;
  int                 Nrp = (degree + 1);
  int                 Nfp = (degree + 1) * (degree + 1);
  sc_array_t         *layers = p6est->layers;

  for (cid = 0, jt = columns->first_local_tree;
       jt <= columns->last_local_tree; ++jt) {
    tree = p4est_tree_array_index (columns->trees, jt);
    tquadrants = &tree->quadrants;

    for (zz = 0; zz < tquadrants->elem_count; ++zz, cid++) {
      p4est_locidx_t      nlayers;
      p4est_locidx_t      nid =
        profile->lnodes->element_nodes[Nfp * cid + Nrp * (Nrp / 2) +
                                       (Nrp / 2)];
      size_t              first, last, zw, zy;

      col = p4est_quadrant_array_index (tquadrants, zz);
      P6EST_COLUMN_GET_RANGE (col, &first, &last);

      nlayers = lr[nid][1];
      p6est_profile_element_to_node_col (profile, cid, offsets,
                                         elem_to_node, fc);
      elem_to_node += nlayers * (degree + 1) * (degree + 1) * (degree + 1);

      for (zy = 0, zw = first; zw < last; zw++, zy++) {
        if (fc[zy] & mask) {
          /* this layer has vertical half faces, we need to set the bit that
           * says whether this is the upper half or the lower half */
          p2est_quadrant_t   *layer;

          layer = p2est_quadrant_array_index (layers, zw);

          if (layer->z & P4EST_QUADRANT_LEN (layer->level)) {
            /* upper half of a pair of layers */
            fc[zy] |= hbit;
          }
        }
      }
      fc += nlayers;
    }
  }
}

static void
p6est_profile_compress (p6est_profile_t * profile)
{
  p4est_locidx_t      nidx, il, old_off, nln =
    profile->lnodes->num_local_nodes;
  p4est_locidx_t (*lr)[2] = (p4est_locidx_t (*)[2]) profile->lnode_ranges;
  sc_array_t         *lc = profile->lnode_columns;
  size_t              old_count = lc->elem_count;
  size_t              new_count;
  sc_array_t         *perm;
  size_t             *newindex;
  size_t              zz, offset;

  if (!old_count) {
    return;
  }
  perm = sc_array_new_size (sizeof (size_t), old_count);
  newindex = (size_t *) sc_array_index (perm, 0);

  for (zz = 0; zz < old_count; zz++) {
    newindex[zz] = old_count;
  }

  offset = 0;

  for (nidx = 0; nidx < nln; nidx++) {
    old_off = lr[nidx][0];
    if (lr[nidx][1]) {
      lr[nidx][0] = offset;
    }
    else {
      P4EST_ASSERT (!lr[nidx][0]);
    }
    for (il = 0; il < lr[nidx][1]; il++) {
      newindex[il + old_off] = offset++;
    }
  }
  new_count = offset;

  for (zz = 0; zz < old_count; zz++) {
    if (newindex[zz] == old_count) {
      newindex[zz] = offset++;
    }
  }

  sc_array_permute (lc, perm, 0);
  sc_array_destroy (perm);
  sc_array_resize (lc, new_count);
}

p6est_profile_t    *
p6est_profile_new_local (p6est_t * p6est,
                         p6est_ghost_t * ghost,
                         p6est_profile_type_t ptype,
                         p8est_connect_type_t btype, int degree)
{
  p6est_profile_t    *profile = P4EST_ALLOC (p6est_profile_t, 1);
  p4est_lnodes_t     *lnodes;
  p4est_locidx_t      nln, nle;
  p4est_topidx_t      jt;
  p4est_t            *columns = p6est->columns;
  p4est_tree_t       *tree;
  sc_array_t         *tquadrants;
  p4est_quadrant_t   *col;
  p4est_qcoord_t      diff = P4EST_ROOT_LEN - p6est->root_len;
  size_t              first, last, count, zz, zy;
  p4est_locidx_t     *en, (*lr)[2];
  sc_array_t         *lc;
  int                 i, j;
  p2est_quadrant_t   *layer;
  sc_array_t         *layers = p6est->layers;
  p4est_locidx_t      nidx, enidx;
  p4est_connect_type_t hbtype;
  int8_t             *c;
  sc_array_t         *thisprof;
  sc_array_t         *selfprof;
  sc_array_t         *faceprof;
  sc_array_t         *cornerprof;
  sc_array_t         *work;
  sc_array_t          oldprof;
  const int           Nrp = degree + 1;

  P4EST_ASSERT (degree > 1);
  profile->ptype = ptype;
  profile->btype = btype;
  profile->lnode_changed[0] = NULL;
  profile->lnode_changed[1] = NULL;
  profile->enode_counts = NULL;
  profile->diff = diff;
  if (btype == P8EST_CONNECT_FACE) {
    hbtype = P4EST_CONNECT_FACE;
  }
  else {
    hbtype = P4EST_CONNECT_FULL;
  }
  if (ghost == NULL) {
    profile->cghost = p4est_ghost_new (p6est->columns, P4EST_CONNECT_FULL);
    profile->ghost_owned = 1;
  }
  else {
    P4EST_ASSERT (ghost->column_ghost->btype == P4EST_CONNECT_FULL);
    profile->cghost = ghost->column_ghost;
    profile->ghost_owned = 0;
  }
  if (ptype == P6EST_PROFILE_UNION) {
    P4EST_ASSERT (degree == 2);
  }
  profile->lnodes = lnodes = p4est_lnodes_new (p6est->columns,
                                               profile->cghost, degree);
  en = lnodes->element_nodes;
  nln = lnodes->num_local_nodes;
  nle = lnodes->num_local_elements;
  profile->lnode_ranges = P4EST_ALLOC_ZERO (p4est_locidx_t, 2 * nln);
  lr = (p4est_locidx_t (*)[2]) profile->lnode_ranges;
  profile->lnode_columns = lc = sc_array_new (sizeof (int8_t));
  selfprof = sc_array_new (sizeof (int8_t));
  work = sc_array_new (sizeof (int8_t));
  faceprof = sc_array_new (sizeof (int8_t));
  cornerprof = sc_array_new (sizeof (int8_t));
  if (ptype == P6EST_PROFILE_UNION) {
    profile->lnode_changed[0] = P4EST_ALLOC (p4est_locidx_t, nln);
    profile->lnode_changed[1] = P4EST_ALLOC (p4est_locidx_t, nln);
    profile->enode_counts = P4EST_ALLOC (p4est_locidx_t, P4EST_INSUL * nle);
    profile->evenodd = 0;
    memset (profile->lnode_changed[0], -1, nln * sizeof (int));
  }

  /* create the profiles for each node: layers are reduced to just their level
   * */
  for (enidx = 0, jt = columns->first_local_tree;
       jt <= columns->last_local_tree; ++jt) {
    tree = p4est_tree_array_index (columns->trees, jt);
    tquadrants = &tree->quadrants;

    for (zz = 0; zz < tquadrants->elem_count; ++zz) {
      col = p4est_quadrant_array_index (tquadrants, zz);
      P6EST_COLUMN_GET_RANGE (col, &first, &last);
      count = last - first;
      sc_array_truncate (selfprof);
      c = (int8_t *) sc_array_push_count (selfprof, count);
      for (zy = first; zy < last; zy++) {
        layer = p2est_quadrant_array_index (layers, zy);
        *(c++) = layer->level;
      }
      if (ptype == P6EST_PROFILE_UNION) {
        p6est_profile_balance_self (selfprof, work);
        if (btype == P8EST_CONNECT_FACE) {
          p6est_profile_balance_face (selfprof, faceprof, work, diff);
        }
        else {
          p6est_profile_balance_full (selfprof, faceprof, work, diff);
        }
        if (btype == P8EST_CONNECT_EDGE) {
          p6est_profile_balance_face (selfprof, cornerprof, work, diff);
        }
        else if (btype == P8EST_CONNECT_FULL) {
          p6est_profile_balance_full (selfprof, cornerprof, work, diff);
        }
      }
      for (j = 0; j < Nrp; j++) {
        for (i = 0; i < Nrp; i++, enidx++) {
          nidx = en[enidx];
          if (ptype == P6EST_PROFILE_UNION) {
            thisprof = NULL;
            if (!(i % degree) && !(j % degree)) {
              if (hbtype == P4EST_CONNECT_FACE) {
                /* skip corners if we don't need to balance them */
                P4EST_ASSERT (!lr[nidx][0]);
                P4EST_ASSERT (!lr[nidx][1]);
                continue;
              }
              else {
                thisprof = cornerprof;
              }
            }
            else if ((i % degree) && (j % degree)) {
              thisprof = selfprof;
            }
            else {
              thisprof = faceprof;
            }
            count = thisprof->elem_count;
            profile->enode_counts[enidx] = count;
            if (!lr[nidx][1]) {
              /* if this node has not yet been initialized, initialize it */
              lr[nidx][0] = lc->elem_count;
              lr[nidx][1] = count;
              c = (int8_t *) sc_array_push_count (lc, count);
              memcpy (c, thisprof->array, count * sizeof (int8_t));
            }
            else {
              /* if this node has been initialized, combine the two profiles,
               * taking the finer layers from each */
              sc_array_init_view (&oldprof, lc, lr[nidx][0], lr[nidx][1]);
              p6est_profile_union (thisprof, &oldprof, work);
              if (work->elem_count > oldprof.elem_count) {
                lr[nidx][0] = lc->elem_count;
                lr[nidx][1] = work->elem_count;
                c = (int8_t *) sc_array_push_count (lc, work->elem_count);
                memcpy (c, work->array, work->elem_count * work->elem_size);
              }
            }
          }
          else {
            count = selfprof->elem_count;
            if (!lr[nidx][1]) {
              /* if this node has not yet been initialized, initialize it */
              lr[nidx][0] = lc->elem_count;
              lr[nidx][1] = count;
              c = (int8_t *) sc_array_push_count (lc, count);
              memcpy (c, selfprof->array, count * sizeof (int8_t));
            }
            else {
              /* if this node has been initialized, combine the two profiles,
               * taking the coarser layers from each */
              sc_array_init_view (&oldprof, lc, lr[nidx][0], lr[nidx][1]);
              p6est_profile_intersection (selfprof, &oldprof, work);
              P4EST_ASSERT (work->elem_count <= oldprof.elem_count);
              if (work->elem_count < oldprof.elem_count) {
                lr[nidx][1] = work->elem_count;
                memcpy (oldprof.array, work->array,
                        work->elem_count * work->elem_size);
              }
            }
          }
        }
      }
    }
  }
  p6est_profile_compress (profile);

  sc_array_destroy (selfprof);
  sc_array_destroy (faceprof);
  sc_array_destroy (cornerprof);
  sc_array_destroy (work);

  return profile;
}

void
p6est_profile_balance_local (p6est_profile_t * profile)
{
  p4est_lnodes_t     *lnodes = profile->lnodes;
  p4est_locidx_t      nln, nle;
  p4est_locidx_t     *en, (*lr)[2];
  sc_array_t         *lc;
  int                 i, j;
  p4est_locidx_t      nidx, enidx, eidx;
  p8est_connect_type_t btype = profile->btype;
  p4est_connect_type_t hbtype;
  int8_t             *c;
  sc_array_t         *thisprof;
  sc_array_t         *selfprof;
  sc_array_t         *faceprof;
  sc_array_t         *cornerprof;
  sc_array_t         *work;
  sc_array_t          oldprof;
  sc_array_t          testprof;
  int                 any_prof_change;
  int                 any_local_change;
  int                 evenodd = profile->evenodd;
  p4est_qcoord_t      diff = profile->diff;

  P4EST_ASSERT (profile->lnodes->degree == 2);

  if (btype == P8EST_CONNECT_FACE) {
    hbtype = P4EST_CONNECT_FACE;
  }
  else {
    hbtype = P4EST_CONNECT_FULL;
  }
  en = lnodes->element_nodes;
  nln = lnodes->num_local_nodes;
  nle = lnodes->num_local_elements;
  lr = (p4est_locidx_t (*)[2]) profile->lnode_ranges;
  lc = profile->lnode_columns;
  selfprof = sc_array_new (sizeof (int8_t));
  work = sc_array_new (sizeof (int8_t));
  faceprof = sc_array_new (sizeof (int8_t));
  cornerprof = sc_array_new (sizeof (int8_t));

  do {
    /* We read from evenodd and write to evenodd ^ 1 */
    memset (&(profile->lnode_changed[evenodd ^ 1][0]), 0, sizeof (int) * nln);
    P4EST_GLOBAL_VERBOSE ("p6est_balance local loop\n");

    any_local_change = 0;
    for (eidx = 0, enidx = 0; eidx < nle; eidx++) {
      p4est_locidx_t      start_enidx = enidx;
      nidx = en[start_enidx + P4EST_INSUL / 2];
      P4EST_ASSERT (lr[nidx][1]);
      sc_array_init_view (&oldprof, lc, lr[nidx][0], lr[nidx][1]);
      thisprof = &oldprof;
      any_prof_change = 0;
      for (j = 0; j < 3; j++) {
        for (i = 0; i < 3; i++, enidx++) {
          nidx = en[enidx];
          if (!profile->lnode_changed[evenodd][nidx]) {
            /* if the profile hasn't changed since I wrote to it, there's no
             * need to balance against it */
            continue;
          }
          if (i != 1 && j != 1) {
            if (hbtype == P4EST_CONNECT_FACE) {
              /* skip corners if we don't need to balance them */
              P4EST_ASSERT (!lr[nidx][0]);
              P4EST_ASSERT (!lr[nidx][1]);
              continue;
            }
          }
          if (i == 1 && j == 1) {
            /* no need to further balance against oneself */
            continue;
          }
          P4EST_ASSERT (lr[nidx][1]);
          P4EST_ASSERT (profile->enode_counts[enidx] <= lr[nidx][1]);
          if (profile->enode_counts[enidx] == lr[nidx][1]) {
            /* if the profile hasn't changed since I wrote to it, there's no
             * need to balance against it */
            continue;
          }
          sc_array_init_view (&testprof, lc, lr[nidx][0], lr[nidx][1]);
          p6est_profile_union (thisprof, &testprof, work);
          if (work->elem_count > thisprof->elem_count) {
            P4EST_ASSERT (profile->lnode_changed[evenodd][nidx]);
            any_prof_change = 1;
            sc_array_copy (selfprof, work);
            thisprof = selfprof;
          }
        }
      }

      if (any_prof_change) {
        P4EST_ASSERT (thisprof == selfprof);
        P4EST_ASSERT (selfprof->elem_count > oldprof.elem_count);
        /* update */
        if (btype == P8EST_CONNECT_FACE) {
          p6est_profile_balance_face (selfprof, faceprof, work, diff);
        }
        else {
          p6est_profile_balance_full (selfprof, faceprof, work, diff);
        }
        if (btype == P8EST_CONNECT_EDGE) {
          p6est_profile_balance_face (selfprof, cornerprof, work, diff);
        }
        else if (btype == P8EST_CONNECT_FULL) {
          p6est_profile_balance_full (selfprof, cornerprof, work, diff);
        }
        enidx = start_enidx;
        for (j = 0; j < 3; j++) {
          for (i = 0; i < 3; i++, enidx++) {
            thisprof = NULL;
            nidx = en[enidx];
            if (i != 1 && j != 1) {
              if (hbtype == P4EST_CONNECT_FACE) {
                /* skip corners if we don't need to balance them */
                P4EST_ASSERT (!lr[nidx][0]);
                P4EST_ASSERT (!lr[nidx][1]);
                continue;
              }
              else {
                thisprof = cornerprof;
              }
            }
            else if (i == 1 && j == 1) {
              thisprof = selfprof;
            }
            else {
              thisprof = faceprof;
            }
            P4EST_ASSERT (lr[nidx][1]);
            /* if this node has been initialized, combine the two profiles,
             * taking the finer layers from each */
            sc_array_init_view (&oldprof, lc, lr[nidx][0], lr[nidx][1]);
            if (i == 1 && j == 1) {
              sc_array_copy (work, thisprof);
            }
            else {
              p6est_profile_union (thisprof, &oldprof, work);
            }
            if (work->elem_count > oldprof.elem_count) {
              if (!(i == 1 && j == 1)) {        /* we don't count changing self */
                profile->lnode_changed[evenodd ^ 1][nidx] = 1;
                any_local_change = 1;
              }
              lr[nidx][0] = lc->elem_count;
              lr[nidx][1] = work->elem_count;
              c = (int8_t *) sc_array_push_count (lc, work->elem_count);
              memcpy (c, work->array, work->elem_count * work->elem_size);
            }
            profile->enode_counts[enidx] = lr[nidx][1];
          }
        }
      }
    }
    p6est_profile_compress (profile);
    evenodd ^= 1;
  } while (any_local_change);

  profile->evenodd = evenodd;
  sc_array_destroy (selfprof);
  sc_array_destroy (faceprof);
  sc_array_destroy (cornerprof);
  sc_array_destroy (work);
}

int
p6est_profile_sync (p6est_profile_t * profile)
{
  p4est_lnodes_t     *lnodes = profile->lnodes;
  p4est_locidx_t      nln = lnodes->num_local_nodes;
  sc_array_t          lrview;
  p4est_lnodes_buffer_t *countbuf;
  sc_array_t         *sharers;
  size_t              zz, nsharers;
  int                 nleft;
  int8_t             *recv, *send;
  int                *array_of_indices;
  p4est_locidx_t      recv_total;
  p4est_locidx_t     *recv_offsets, recv_offset;
  p4est_locidx_t      send_total;
  p4est_locidx_t     *send_offsets, send_offset;
  p4est_locidx_t (*lr)[2];
  sc_array_t         *lc = profile->lnode_columns;
  sc_MPI_Request     *recv_request, *send_request;
  sc_array_t         *work;
  int                 any_change = 0;
  int                 any_global_change;
  int                 mpiret, mpirank;
  int                 evenodd = profile->evenodd;

  lr = (p4est_locidx_t (*)[2]) profile->lnode_ranges;
  sharers = lnodes->sharers;
  nsharers = sharers->elem_count;

  mpiret = sc_MPI_Comm_rank (lnodes->mpicomm, &mpirank);
  SC_CHECK_MPI (mpiret);

  sc_array_init_data (&lrview, lr, 2 * sizeof (p4est_locidx_t), nln);

  countbuf = p4est_lnodes_share_all_begin (&lrview, lnodes);
  send_offsets = P4EST_ALLOC (p4est_locidx_t, nsharers + 1);
  send_offset = 0;
  for (zz = 0; zz < nsharers; zz++) {
    p4est_lnodes_rank_t *sharer;
    sc_array_t         *send_buf;
    size_t              zy, nnodes;

    send_offsets[zz] = send_offset;
    sharer = p4est_lnodes_rank_array_index (sharers, zz);
    if (sharer->rank == mpirank) {
      continue;
    }
    send_buf = (sc_array_t *) sc_array_index (countbuf->send_buffers, zz);
    nnodes = sharer->shared_nodes.elem_count;

    P4EST_ASSERT (nnodes == send_buf->elem_count);

    P4EST_ASSERT (send_buf->elem_size == 2 * sizeof (p4est_locidx_t));
    for (zy = 0; zy < nnodes; zy++) {
      p4est_locidx_t     *lp =
        (p4est_locidx_t *) sc_array_index (send_buf, zy);
      P4EST_ASSERT (lp[0] >= 0);
      P4EST_ASSERT (lp[1] >= 0);
      send_offset += lp[1];
    }
  }
  send_total = send_offsets[nsharers] = send_offset;

  p4est_lnodes_share_all_end (countbuf);
  recv_offsets = P4EST_ALLOC (p4est_locidx_t, nsharers + 1);
  recv_offset = 0;
  for (zz = 0; zz < nsharers; zz++) {
    p4est_lnodes_rank_t *sharer;
    sc_array_t         *recv_buf;
    size_t              zy, nnodes;

    recv_offsets[zz] = recv_offset;
    sharer = p4est_lnodes_rank_array_index (sharers, zz);
    if (sharer->rank == mpirank) {
      continue;
    }
    recv_buf = (sc_array_t *) sc_array_index (countbuf->recv_buffers, zz);
    nnodes = sharer->shared_nodes.elem_count;

    P4EST_ASSERT (nnodes == recv_buf->elem_count);

    P4EST_ASSERT (recv_buf->elem_size == 2 * sizeof (p4est_locidx_t));
    for (zy = 0; zy < nnodes; zy++) {
      p4est_locidx_t     *lp =
        (p4est_locidx_t *) sc_array_index (recv_buf, zy);
      P4EST_ASSERT (lp[0] >= 0);
      P4EST_ASSERT (lp[1] >= 0);
      recv_offset += lp[1];
    }
  }
  recv_total = recv_offsets[nsharers] = recv_offset;

  recv = P4EST_ALLOC (int8_t, recv_total);
  recv_request = P4EST_ALLOC (sc_MPI_Request, nsharers);
  send = P4EST_ALLOC (int8_t, send_total);
  send_request = P4EST_ALLOC (sc_MPI_Request, nsharers);

  /* post receives */
  nleft = 0;
  for (zz = 0; zz < nsharers; zz++) {
    p4est_lnodes_rank_t *sharer;
    int                 icount = recv_offsets[zz + 1] - recv_offsets[zz];

    sharer = p4est_lnodes_rank_array_index (sharers, zz);
    if (sharer->rank == mpirank) {
      recv_request[zz] = sc_MPI_REQUEST_NULL;
      continue;
    }
    if (icount) {
      mpiret =
        sc_MPI_Irecv (recv + recv_offsets[zz], icount * sizeof (int8_t),
                      sc_MPI_BYTE, sharer->rank, P6EST_COMM_BALANCE,
                      lnodes->mpicomm, recv_request + zz);
      SC_CHECK_MPI (mpiret);
      nleft++;
    }
    else {
      recv_request[zz] = sc_MPI_REQUEST_NULL;
    }
  }

  /* post sends */
  for (zz = 0; zz < nsharers; zz++) {
    p4est_lnodes_rank_t *sharer;
    size_t              zy, nnodes;
    int                 icount;
    sc_array_t         *shared_nodes;

    sharer = p4est_lnodes_rank_array_index (sharers, zz);
    if (sharer->rank == mpirank) {
      send_request[zz] = sc_MPI_REQUEST_NULL;
      continue;
    }
    shared_nodes = &sharer->shared_nodes;
    nnodes = shared_nodes->elem_count;
    icount = 0;
    for (zy = 0; zy < nnodes; zy++) {
      p4est_locidx_t      nidx;
      int8_t             *c;

      nidx = *((p4est_locidx_t *) sc_array_index (shared_nodes, zy));

      if (lr[nidx][1]) {
        c = (int8_t *) sc_array_index (lc, lr[nidx][0]);
        memcpy (send + send_offsets[zz] + icount, c,
                lr[nidx][1] * sizeof (int8_t));
        icount += lr[nidx][1];
      }
      else {
        P4EST_ASSERT (!lr[nidx][0]);
      }
    }
    P4EST_ASSERT (icount == send_offsets[zz + 1] - send_offsets[zz]);
    if (icount) {
      mpiret =
        sc_MPI_Isend (send + send_offsets[zz], icount * sizeof (int8_t),
                      sc_MPI_BYTE, sharer->rank, P6EST_COMM_BALANCE,
                      lnodes->mpicomm, send_request + zz);
      SC_CHECK_MPI (mpiret);
    }
    else {
      send_request[zz] = sc_MPI_REQUEST_NULL;
    }
  }

  work = sc_array_new (sizeof (int8_t));
  array_of_indices = P4EST_ALLOC (int, nsharers);
  while (nleft) {
    int                 outcount;
    int                 i;

    mpiret = sc_MPI_Waitsome (nsharers, recv_request, &outcount,
                              array_of_indices, sc_MPI_STATUSES_IGNORE);
    SC_CHECK_MPI (mpiret);

    for (i = 0; i < outcount; i++) {
      p4est_lnodes_rank_t *sharer;
      size_t              zy, nnode;
      sc_array_t         *shared_nodes;
      sc_array_t         *recv_buf;

      zz = array_of_indices[i];
      sharer = p4est_lnodes_rank_array_index (sharers, zz);
      shared_nodes = &sharer->shared_nodes;
      recv_buf = (sc_array_t *) sc_array_index (countbuf->recv_buffers, zz);
      nnode = shared_nodes->elem_count;
      P4EST_ASSERT (nnode == recv_buf->elem_count);

      recv_offset = recv_offsets[zz];
      for (zy = 0; zy < nnode; zy++) {
        p4est_locidx_t     *lp;
        p4est_locidx_t      nidx;
        sc_array_t          oldview, newview;

        nidx = *((p4est_locidx_t *) sc_array_index (shared_nodes, zy));
        lp = (p4est_locidx_t *) sc_array_index (recv_buf, zy);

        sc_array_init_view (&oldview, lc, lr[nidx][0], lr[nidx][1]);
        sc_array_init_data (&newview, recv + recv_offset, sizeof (int8_t),
                            lp[1]);
        if (profile->ptype == P6EST_PROFILE_UNION) {
          p6est_profile_union (&oldview, &newview, work);

          if (work->elem_count > oldview.elem_count) {
            int8_t             *c;

            any_change = 1;
            lr[nidx][0] = lc->elem_count;
            lr[nidx][1] = work->elem_count;
            profile->lnode_changed[evenodd][nidx] = 1;

            c = (int8_t *) sc_array_push_count (lc, work->elem_count);
            memcpy (c, work->array, work->elem_count * work->elem_size);
          }
        }
        else {
          p6est_profile_intersection (&oldview, &newview, work);
          P4EST_ASSERT (work->elem_count <= oldview.elem_count);
          if (work->elem_count < oldview.elem_count) {
            lr[nidx][1] = work->elem_count;
            memcpy (oldview.array, work->array,
                    work->elem_count * work->elem_size);
          }
        }

        recv_offset += lp[1];
      }
      P4EST_ASSERT (recv_offset == recv_offsets[zz + 1]);
    }

    nleft -= outcount;
    P4EST_ASSERT (nleft >= 0);
  }
  P4EST_FREE (array_of_indices);
  sc_array_destroy (work);

  p6est_profile_compress (profile);
  p4est_lnodes_buffer_destroy (countbuf);

  P4EST_FREE (recv_request);
  P4EST_FREE (recv_offsets);
  P4EST_FREE (recv);

  {
    mpiret = sc_MPI_Waitall (nsharers, send_request, sc_MPI_STATUSES_IGNORE);

    SC_CHECK_MPI (mpiret);
    P4EST_FREE (send_request);
    P4EST_FREE (send_offsets);
    P4EST_FREE (send);

    any_global_change = any_change;
    mpiret = sc_MPI_Allreduce (&any_change, &any_global_change, 1, sc_MPI_INT,
                               sc_MPI_LOR, lnodes->mpicomm);

    SC_CHECK_MPI (mpiret);
  }

  return any_global_change;
}

void
p6est_profile_destroy (p6est_profile_t * profile)
{
  p4est_lnodes_destroy (profile->lnodes);
  if (profile->ghost_owned) {
    p4est_ghost_destroy (profile->cghost);
  }
  if (profile->lnode_changed[0]) {
    P4EST_ASSERT (profile->lnode_changed[1]);
    P4EST_FREE (profile->lnode_changed[0]);
    P4EST_FREE (profile->lnode_changed[1]);
    P4EST_ASSERT (profile->enode_counts);
    P4EST_FREE (profile->enode_counts);
  }
  P4EST_FREE (profile->lnode_ranges);
  sc_array_destroy (profile->lnode_columns);
  P4EST_FREE (profile);
}

void
p6est_refine_to_profile (p6est_t * p6est, p6est_profile_t * profile,
                         p6est_init_t init_fn, p6est_replace_t replace_fn)
{
  size_t              zz, zy, first, last;
  p4est_topidx_t      jt;
  p4est_quadrant_t   *col;
  p4est_tree_t       *tree;
  sc_array_t         *tquadrants;
  p4est_locidx_t      eidx;
  p4est_locidx_t     *en = profile->lnodes->element_nodes;
  p4est_locidx_t (*lr)[2];
  p4est_locidx_t      nidx, pidx, pfirst, plast;
  sc_array_t         *layers = p6est->layers;
  sc_array_t         *lc = profile->lnode_columns;
  sc_array_t         *work;

  P4EST_ASSERT (profile->lnodes->degree == 2);

  lr = (p4est_locidx_t (*)[2]) profile->lnode_ranges;
  work = sc_array_new (sizeof (p2est_quadrant_t));
  for (eidx = 0, jt = p6est->columns->first_local_tree;
       jt <= p6est->columns->last_local_tree; ++jt) {
    tree = p4est_tree_array_index (p6est->columns->trees, jt);
    tquadrants = &tree->quadrants;
    for (zz = 0; zz < tquadrants->elem_count; ++zz, eidx++) {

      col = p4est_quadrant_array_index (tquadrants, zz);
      P6EST_COLUMN_GET_RANGE (col, &first, &last);
      nidx = en[P4EST_INSUL * eidx + P4EST_INSUL / 2];
      P4EST_ASSERT ((size_t) lr[nidx][1] >= last - first);
      pfirst = lr[nidx][0];
      plast = pfirst + lr[nidx][1];
      if ((size_t) lr[nidx][1] > last - first) {
        p2est_quadrant_t    stack[P4EST_QMAXLEVEL];
        p2est_quadrant_t   *q, *r, s, t;
        int                 stackcount;

        sc_array_truncate (work);
        stackcount = 0;
        zy = first;
        for (pidx = pfirst; pidx < plast; pidx++) {
          int8_t              p;

          P4EST_ASSERT (stackcount || zy < last);

          p = *((int8_t *) sc_array_index (lc, pidx));

          if (stackcount) {
            q = &(stack[--stackcount]);
          }
          else {
            q = p2est_quadrant_array_index (layers, zy++);
          }

          P4EST_ASSERT (q->level <= p);
          while (q->level < p) {
            p2est_quadrant_t   *child[2];

            t = *q;
            s = *q;
            s.level++;
            stack[stackcount] = s;
            stack[stackcount].z += P4EST_QUADRANT_LEN (s.level);
            child[0] = &s;
            child[1] = &stack[stackcount++];
            p6est_layer_init_data (p6est, jt, col, child[0], init_fn);
            p6est_layer_init_data (p6est, jt, col, child[1], init_fn);
            q = &t;
            if (replace_fn) {
              replace_fn (p6est, jt, 1, 1, &col, &q, 1, 2, &col, child);
            }
            p6est_layer_free_data (p6est, &t);
            q = &s;
          }
          r = p2est_quadrant_array_push (work);
          *r = *q;
        }
        P4EST_ASSERT (work->elem_count == (size_t) lr[nidx][1]);
        first = layers->elem_count;
        last = first + work->elem_count;
        P6EST_COLUMN_SET_RANGE (col, first, last);
        q = (p2est_quadrant_t *) sc_array_push_count (layers,
                                                      work->elem_count);
        memcpy (q, work->array, work->elem_count * work->elem_size);
      }
    }
  }
  sc_array_destroy (work);
  p6est_compress_columns (p6est);
  p6est_update_offsets (p6est);
}
