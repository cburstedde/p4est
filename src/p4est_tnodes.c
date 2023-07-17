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

#include <p4est_bits.h>
#include <p4est_iterate.h>
#include <p4est_tnodes.h>

static const int n_ccorn[4] = { 0, 1, 2, 3 };
static const int n_center = 4;
static const int n_mface[4] = { 5, 6, 7, 8 };
static const int n_cface[4] = { 9, 10, 11, 12 };
static const int n_split[4] = { 14, 17, 20, 22 };
static const int n_hface[4][2] = {{ 13, 15 }, { 16, 18 }, { 19, 21 }, { 23, 24 }};
#ifdef P4EST_ENABLE_DEBUG
static const int alwaysowned[25] =
  { 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1,
    0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0 };
#endif

/** A single contributor process to a node under construction. */
typedef struct tnodes_contr
{
  int                 nodene;       /**< Relative to element. */
  int                 rank;         /**< The referring process. */
  p4est_locidx_t      le;           /**< Element/ghost number. */
}
tnodes_contr_t;

/** A node under construction may have several contributors. */
typedef struct tnodes_cnode
{
  p4est_locidx_t      runid;        /**< Running count of node. */
  p4est_connect_type_t bcon;        /**< Codimension of node. */
  tnodes_contr_t     *owner;        /**< Owning contributor. */
  sc_array_t          contr;        /**< Contributing processes. */
}
tnodes_cnode_t;

#ifdef P4EST_ENABLE_MPI

/** Record one communication partner and/or node sharer.
 *
 * A peer may be either a sender of queries, a replier to queries,
 * or a peer with which the current rank does not communicate, but
 * shares interest in a node with the current rank (both not owners).
 *
 * The three cases are distinguished as follows:
 *  - sender of queries adds to
 *     - \b localind the temporary local index (negative) of a node
 *     - \b querypos the position of the node within the recipiend
 *  - replier to queries increases \b bufcount
 *  - sharer keeps the above untouched
 *
 * A sharer may be promoted to a sender or replier at any time.
 */
typedef struct tnodes_peer
{
  int                 rank;
  int                 done;
  p4est_locidx_t      lastadd;
  p4est_locidx_t      bufcount;
  p4est_locidx_t      shacumul;
  sc_array_t          sharedno;
  sc_array_t          localind;
  sc_array_t          querypos;
}
tnodes_peer_t;

#endif

/** Global control structure for the tnodes algorithm. */
typedef struct tnodes_meta
{
  int                 full_style;
  int                 with_faces;
  int                 mpisize, mpirank;
  int                *ghost_rank;
  int                *proc_peer;
  int                 peer_with_self;
  uint8_t            *chilev;
  sc_MPI_Comm         mpicomm;
  sc_array_t          remotepos;
  sc_array_t          sortp;
  sc_array_t          peers;
  sc_array_t          pereq;
  sc_array_t          oldtolocal;
  sc_array_t          construct;
  p4est_locidx_t      lenum;
  p4est_locidx_t      num_owned;
  p4est_locidx_t      num_shared;
  p4est_locidx_t      szero[25];
  p4est_locidx_t      smone[25];
  p4est_gloidx_t     *goffset;
  p4est_t            *p4est;
  p4est_ghost_t      *ghost;
  p4est_tnodes_t     *tm;

  /* the above members will be cleaned up eventually */
  sc_array_t          ownsort;      /**< Sorted owned nodes. */
}
tnodes_meta_t;

#if 0

#if defined P4EST_ENABLE_MPI && defined P4EST_ENABLE_DEBUG

/* *INDENT-OFF* */
static const int    pos_is_boundary[25] =
{ 0, 1, 1, 1, 1, 1, 1, 1, 1,
  0, 0, 0, 0, 0, 0, 0, 0,
  1, 1, 1, 1, 1, 1, 1, 1
};

/* *INDENT_ON* */

#endif

static void
set_lnodes_face_full (tnodes_meta_t * me, p4est_locidx_t le,
                      int face, p4est_locidx_t lni)
{
  p4est_lnodes_t     *ln = me->tm->lnodes;
  p4est_locidx_t      lpos;

  P4EST_ASSERT (ln != NULL);
  P4EST_ASSERT (ln->vnodes == 25);
  P4EST_ASSERT (0 <= le && le < ln->num_local_elements);
  P4EST_ASSERT (0 <= face && face < P4EST_FACES);

  lpos = le * ln->vnodes + 123;
  P4EST_ASSERT (ln->element_nodes[lpos] == 0);
  P4EST_ASSERT (ln->element_nodes[lpos + 1] == 0);
  ln->element_nodes[lpos] = lni;

  if (lni < 0) {
    /* save every element node position with remotely owned node */
    *(p4est_locidx_t *) sc_array_push (&me->remotepos) = lpos;
  }
}

#ifdef P4EST_ENABLE_MPI

static tnodes_peer_t *
peer_access (tnodes_meta_t * me, int q)
{
  int                 pi;
  tnodes_peer_t      *peer;

  P4EST_ASSERT (me != NULL);
  P4EST_ASSERT (me->ghost_rank != NULL);
  P4EST_ASSERT (me->proc_peer != NULL);
  P4EST_ASSERT (0 <= q && q < me->mpisize);

  if ((pi = me->proc_peer[q]) == 0) {
    peer = (tnodes_peer_t *) sc_array_push (&me->peers);
    me->proc_peer[q] = (int) me->peers.elem_count;
    peer->rank = q;
    peer->done = 0;
    peer->lastadd = 0;
    peer->bufcount = 0;
    sc_array_init (&peer->sharedno, sizeof (p4est_locidx_t));
    sc_array_init (&peer->localind, sizeof (p4est_locidx_t));
    sc_array_init (&peer->querypos, sizeof (p4est_locidx_t));
  }
  else {
    P4EST_ASSERT (0 < pi && pi <= me->mpisize);
    peer = (tnodes_peer_t *) sc_array_index_int (&me->peers, pi - 1);
    P4EST_ASSERT (peer->rank == q);
  }
  return peer;
}

static void
peer_add_share (tnodes_peer_t * peer, p4est_locidx_t lni)
{
  P4EST_ASSERT (peer != NULL);
  P4EST_ASSERT (lni != 0);      /*< owned node 0 is always non-shared */

  if (peer->lastadd != lni) {
    *(p4est_locidx_t *) sc_array_push (&peer->sharedno) = peer->lastadd = lni;
  }
}

static void
peer_add_reply (tnodes_peer_t * peer, p4est_locidx_t lni)
{
  P4EST_ASSERT (peer != NULL);
  P4EST_ASSERT (lni > 0);

  P4EST_ASSERT (peer->lastadd <= lni);
  if (peer->lastadd != lni) {
    ++peer->bufcount;
    *(p4est_locidx_t *) sc_array_push (&peer->sharedno) = peer->lastadd = lni;
  }
}

static void
peer_add_query (tnodes_peer_t * peer, p4est_locidx_t gpos, p4est_locidx_t lni)
{
  P4EST_ASSERT (peer != NULL);
  P4EST_ASSERT (gpos >= 0);
  P4EST_ASSERT (lni < 0);
  P4EST_ASSERT (peer->localind.elem_count == peer->querypos.elem_count);

  P4EST_ASSERT (peer->lastadd >= lni);
  if (peer->lastadd != lni) {
    *(p4est_locidx_t *) sc_array_push (&peer->localind) = lni;
    *(p4est_locidx_t *) sc_array_push (&peer->querypos) = gpos;
    *(p4est_locidx_t *) sc_array_push (&peer->sharedno) = peer->lastadd = lni;
  }
}

#endif /* P4EST_ENABLE_MPI */

#endif /* 0 */

static p4est_locidx_t
tree_quad_to_le (p4est_t *p4est, p4est_topidx_t which_tree, p4est_locidx_t quadid)
{
  p4est_tree_t      *tree;

  P4EST_ASSERT (p4est != NULL);
  P4EST_ASSERT (p4est->trees != NULL);
  tree = p4est_tree_array_index (p4est->trees, which_tree);

  P4EST_ASSERT (0 <= quadid && quadid < (p4est_locidx_t) tree->quadrants.elem_count);
  return tree->quadrants_offset + quadid;
}

/** Register a node position relative to an element.
 * The element is either process local or a ghost.
 * Multiple positions may reference the same local node.
 * We store only the lowes referrer for each process.
 */
static void
node_register (tnodes_meta_t *me, p4est_locidx_t *lni,
               int rank, p4est_locidx_t le, int nodene, p4est_connect_type_t bcon)
{
  tnodes_cnode_t    *cnode;
  tnodes_contr_t    *contr;
  p4est_lnodes_t    *ln;
  p4est_locidx_t     lnis;
  int                owner_rank;
  size_t             zz, siz;

  /* basic checks */
  P4EST_ASSERT (me != NULL);
  P4EST_ASSERT (me->tm != NULL);
  ln = me->tm->lnodes;

  /* a new node is to be created or an existing one is passed in */
  if (lni == NULL) {
    lnis = -1;
    lni = &lnis;
  }
  P4EST_ASSERT (-1 <= *lni && *lni < (p4est_locidx_t) me->construct.elem_count);

  /* abbreviate local rank */
  P4EST_ASSERT (rank >= -1 && rank != me->mpirank);
  if (rank == -1) {
    rank = me->mpirank;
  }
  P4EST_ASSERT (0 <= rank && rank < me->mpisize);

  /* check remaining arguments */
  P4EST_ASSERT (0 <= le && le < (rank == me->mpirank ?
                ln->num_local_elements :
                (p4est_locidx_t) me->ghost->ghosts.elem_count));
  P4EST_ASSERT (0 <= nodene && nodene < ln->vnodes);
  P4EST_ASSERT (bcon == P4EST_CONNECT_FACE || bcon == P4EST_CONNECT_CORNER);
  P4EST_ASSERT (bcon != P4EST_CONNECT_FACE || (nodene >= 4));
  P4EST_ASSERT (bcon != P4EST_CONNECT_CORNER || (0 <= nodene && nodene < 9));

  if (*lni == -1) {
    /* create a new node with one instance */
    *lni = (p4est_locidx_t) me->construct.elem_count;
    cnode = (tnodes_cnode_t *) sc_array_push (&me->construct);
    cnode->runid = *lni;
    cnode->bcon = bcon;
    cnode->owner = NULL;
    sc_array_init (&cnode->contr, sizeof (tnodes_contr_t));
  }
  else {
    /* create a new instance of an existing node */
    cnode = (tnodes_cnode_t *) sc_array_index (&me->construct, (size_t) *lni);
    P4EST_ASSERT (cnode->runid == *lni);
    P4EST_ASSERT (cnode->bcon == bcon);
    P4EST_ASSERT (cnode->contr.elem_size == sizeof (tnodes_contr_t));
    P4EST_ASSERT (cnode->contr.elem_count > 0);
    P4EST_ASSERT (cnode->owner != NULL);
  }

  /* assign node to the local element position */
  if (rank == me->mpirank) {
    P4EST_ASSERT (ln->element_nodes[le * ln->vnodes + nodene] == -1);
    ln->element_nodes[le * ln->vnodes + nodene] = *lni;
  }

  /* iterate through instances to find matching process */
  siz = cnode->contr.elem_count;
  P4EST_ASSERT (siz == 0 || cnode->owner != NULL);
  for (zz = 0; zz < siz; ++zz) {
    contr = (tnodes_contr_t *) sc_array_index (&cnode->contr, zz);
    if (contr->rank == rank) {
      /* rank is found and we remember the smallest node position */
      if (le < contr->le || (le == contr->le && nodene < contr->nodene)) {
        contr->nodene = nodene;
        contr->le = le;
      }
      return;
    }
  }

  /* add new node process to the list for this node */
  owner_rank = cnode->owner == NULL ? -1 : cnode->owner->rank;
  contr = (tnodes_contr_t *) sc_array_push (&cnode->contr);
  contr->nodene = nodene;
  contr->rank = rank;
  contr->le = le;
  if (cnode->owner == NULL || rank <= owner_rank) {
    /* the array was empty before or we know to set a new owner */
    P4EST_ASSERT (rank != owner_rank);
    cnode->owner = contr;
  }
  else {
    /* pushing to the array has invalidated the previous owner pointer */
    cnode->owner = (tnodes_contr_t *) sc_array_index (&cnode->contr, 0);
    siz = cnode->contr.elem_count;
    for (zz = 1; zz < siz - 1; ++zz) {
      contr = (tnodes_contr_t *) sc_array_index (&cnode->contr, zz);
      if (contr->rank < cnode->owner->rank) {
        cnode->owner = contr;
      }
    }
  }
}

static void
node_lregister (tnodes_meta_t *me, p4est_locidx_t *lni,
                p4est_locidx_t le, int nodene, p4est_connect_type_t bcon)
{
  node_register (me, lni, -1, le, nodene, bcon);
}

static void
node_lfacetocorner (tnodes_meta_t *me, p4est_locidx_t le, int nodene)
{
  tnodes_cnode_t    *cnode;
  p4est_lnodes_t    *ln;
  p4est_locidx_t     lni;

  P4EST_ASSERT (me != NULL);
  P4EST_ASSERT (me->tm != NULL);
  ln = me->tm->lnodes;
  P4EST_ASSERT (0 <= le && le < ln->num_local_elements);

  P4EST_ASSERT (nodene == n_center);

  /* access node in local element */
  lni = ln->element_nodes[le * ln->vnodes + nodene];
  P4EST_ASSERT (lni >= 0);
  cnode = (tnodes_cnode_t *) sc_array_index (&me->construct, lni);
  P4EST_ASSERT (cnode->owner != NULL);
  P4EST_ASSERT (cnode->runid == lni);
  P4EST_ASSERT (cnode->bcon == P4EST_CONNECT_FACE);
  P4EST_ASSERT (cnode->contr.elem_size == sizeof (tnodes_contr_t));
  P4EST_ASSERT (cnode->contr.elem_count == 1);

  /* change instance of an existing node */
  cnode->bcon = P4EST_CONNECT_CORNER;
}

static void
node_gregister (tnodes_meta_t *me, p4est_locidx_t *lni,
                p4est_locidx_t ghostid, int nodene, p4est_connect_type_t bcon)
{
  P4EST_ASSERT (me != NULL);
  P4EST_ASSERT (me->tm != NULL);
  P4EST_ASSERT (0 <= nodene && nodene < me->tm->lnodes->vnodes);
  P4EST_ASSERT (!alwaysowned[nodene]);

  if (me->ghost != NULL) {
    P4EST_ASSERT (me->ghost_rank != NULL);
    P4EST_ASSERT (0 <= ghostid &&
                  ghostid < (p4est_locidx_t) me->ghost->ghosts.elem_count);

    node_register (me, lni, me->ghost_rank[ghostid], ghostid, nodene, bcon);
  }
}

static void
iter_volume1 (p4est_iter_volume_info_t * vi, void *user_data)
{
  tnodes_meta_t      *me = (tnodes_meta_t *) user_data;
  p4est_tnodes_t     *tm = me->tm;
  p4est_locidx_t      le;
  int                 j;
  int                 childid;
  int8_t              level;
#ifdef P4EST_ENABLE_DEBUG
  p4est_lnodes_t     *ln = tm->lnodes;
  p4est_tree_t       *tree;

  /* initial checks  */
  P4EST_ASSERT (vi->p4est == me->p4est);
  tree = p4est_tree_array_index (vi->p4est->trees, vi->treeid);
  P4EST_ASSERT (tree->quadrants_offset + vi->quadid == me->lenum);
#endif

  /* store quadrant level and child id */
  le = me->lenum++;
  level = vi->quad->level;
  childid = p4est_quadrant_child_id (vi->quad);
  me->chilev[le] = (((uint8_t) level) << 3) | ((uint8_t) childid);
  P4EST_ASSERT (tm->configuration[le] == 0);
  P4EST_ASSERT (ln->face_code[le] == 0);
  P4EST_ASSERT (!memcmp (ln->element_nodes + le * ln->vnodes,
                         me->smone, sizeof (p4est_locidx_t) * ln->vnodes));

  /* add nodes as required */
  if (me->full_style || level == 0) {
    tm->configuration[le] = (((uint8_t) 1) << 5);
    node_lregister (me, NULL, le, n_center, P4EST_CONNECT_CORNER);
    if (me->with_faces) {
      for (j = 0; j < 4; ++j) {
        node_lregister (me, NULL, le, n_cface[j], P4EST_CONNECT_FACE);
      }
    }
  }
  else {
    if (childid == 1 || childid == 2) {
      tm->configuration[le] = (((uint8_t) 1) << 4);
    }
    if (me->with_faces) {
      node_lregister (me, NULL, le, n_center, P4EST_CONNECT_FACE);
    }
  }
}

static void
iter_face1 (p4est_iter_face_info_t * fi, void *user_data)
{
  tnodes_meta_t      *me = (tnodes_meta_t *) user_data;
  int                 i, j;
  int                 face;
  int                 nodene;
  int                 childid;
  int                 swapi;
#if 0
  /* each face connection produces at most 3 nodes: 1 corner, 2 face */
  int                 nunodes;          /**< nodes on interface */
  int                 codim[3];         /**< codimension of a node */
  int                 is_owned[3];      /**< is that node locally owned */
  int                 is_shared[3];     /**< does the node have sharers */
  int                 sharers[3][3];    /**< sharer processes for each node */
  int                 owner[3];         /**< owner process for each node */
#endif
  p4est_locidx_t      le;               /**< local element number */
  p4est_locidx_t      lni;              /**< local node number */
  p4est_locidx_t      lnh[2];           /**< hanging face mids */
  p4est_iter_face_side_t *fs;
  p4est_iter_face_side_t *fss[2];
  p4est_iter_face_side_full_t *fu;
  p4est_iter_face_side_hanging_t *fh;
  p4est_lnodes_t     *ln = me->tm->lnodes;
#if 0
#ifdef P4EST_ENABLE_MPI
  p4est_locidx_t      gpos[3][3];       /**< position within ghost */
  p4est_locidx_t      igi;              /**< iterator ghost index */
  p4est_quadrant_t   *gquad;
  tnodes_peer_t      *peer;
#endif
#endif

  /* initial checks  */
  P4EST_ASSERT (fi->p4est == me->p4est);

  /* a boundary face is the easiest case */
  fs = fss[0] = fss[1] = NULL;
  if (fi->sides.elem_count == 1) {
    P4EST_ASSERT (fi->orientation == 0);
    P4EST_ASSERT (fi->tree_boundary == P4EST_CONNECT_FACE);
    fs = (p4est_iter_face_side_t *) sc_array_index_int (&fi->sides, 0);
    P4EST_ASSERT (!fs->is_hanging);
    fu = &fs->is.full;
    P4EST_ASSERT (!fu->is_ghost);
    /* a boundary face does not contribute to the configuration */
    if (me->with_faces) {
      le = tree_quad_to_le (fi->p4est, fs->treeid, fu->quadid);
      node_lregister (me, NULL, le, n_mface[fs->face], P4EST_CONNECT_FACE);
    }
    return;
  }

  /* we have two sides to the face connection */
  P4EST_ASSERT (fi->sides.elem_count == 2);
  fss[0] = (p4est_iter_face_side_t *) sc_array_index_int (&fi->sides, 0);
  fss[1] = (p4est_iter_face_side_t *) sc_array_index_int (&fi->sides, 1);
  P4EST_ASSERT (!fss[0]->is_hanging || !fss[1]->is_hanging);
  if (!fss[0]->is_hanging && !fss[1]->is_hanging) {
    /* same size face connection does not contribute to configuration */
    if (me->with_faces) {
      lni = -1;
      for (i = 0; i < 2; ++i) {
        fu = &fss[i]->is.full;
        nodene = n_mface[fss[i]->face];
        if (!fu->is_ghost) {
          le = tree_quad_to_le (fi->p4est, fss[i]->treeid, fu->quadid);
          node_lregister (me, &lni, le, nodene, P4EST_CONNECT_FACE);
        }
        else {
          node_gregister (me, &lni, fu->quadid, nodene, P4EST_CONNECT_FACE);
        }
      }
    }
    return;
  }

  /* one of the two sides is hanging */
  lni = lnh[0] = lnh[1] = -1;
  for (i = 0; i < 2; ++i) {
    swapi = (i == 0 || fi->orientation == 0) ? 0 : 1;
    if (!fss[i]->is_hanging) {
      fu = &fss[i]->is.full;
      face = fss[i]->face;
      nodene = n_mface[face];

      if (!fu->is_ghost) {
        /* this is a large local quadrant which must insert the face midpoint */
        le = tree_quad_to_le (fi->p4est, fss[i]->treeid, fu->quadid);
        if ((me->tm->configuration[le] & ~(1 << 4)) == 0) {
          /* this is a half refinement, which must be promoted to full */
          P4EST_ASSERT (!me->full_style && fu->quad->level > 0);
          if (!me->with_faces) {
            node_lregister (me, NULL, le, n_center, P4EST_CONNECT_CORNER);
          }
          else {
            node_lfacetocorner (me, le, n_center);
            for (j = 0; j < 4; ++j) {
              node_lregister (me, NULL, le, n_cface[j], P4EST_CONNECT_FACE);
            }
          }
        }
        me->tm->configuration[le] &= ~((1 << 4) | (1 << 5));
        me->tm->configuration[le] |= (1 << face);
        node_lregister (me, &lni, le, nodene, P4EST_CONNECT_CORNER);
        if (me->with_faces) {
          node_lregister (me, NULL, le, n_split[face], P4EST_CONNECT_FACE);
          for (j = 0; j < 2; ++j) {
            node_lregister (me, &lnh[j ^ swapi], le, n_hface[face][j], P4EST_CONNECT_FACE);
          }
        }
      }
      else {
        node_gregister (me, &lni, fu->quadid, nodene, P4EST_CONNECT_CORNER);
        if (me->with_faces) {
          for (j = 0; j < 2; ++j) {
            node_gregister (me, &lnh[j ^ swapi], fu->quadid, n_hface[face][j], P4EST_CONNECT_FACE);
          }
        }
      }
    }
    else {
      P4EST_ASSERT (fss[i]->is_hanging);

      /* for each small local quadrant contribute to its face code */
      fh = &fss[i]->is.hanging;
      face = fss[i]->face;
      for (j = 0; j < P4EST_HALF; ++j) {
        nodene = n_ccorn[p4est_face_corners[face][j ^ (P4EST_HALF - 1)]];
        if (!fh->is_ghost[j]) {
          le = tree_quad_to_le (fi->p4est, fss[i]->treeid, fh->quadid[j]);
          node_lregister (me, &lni, le, nodene, P4EST_CONNECT_CORNER);
          if (me->with_faces) {
            node_lregister (me, &lnh[j ^ swapi], le, n_mface[face], P4EST_CONNECT_FACE);
          }

          /* update the face code */
          childid = p4est_quadrant_child_id (fh->quad[j]);
          P4EST_ASSERT (p4est_face_corners[face][j] == childid);
          P4EST_ASSERT (face == p4est_corner_faces[childid][face >> 1]);
          P4EST_ASSERT ((ln->face_code[le] &
                         (1 << (P4EST_DIM + (face >> 1)))) == 0);
          ln->face_code[le] |= (1 << (P4EST_DIM + (face >> 1))) | childid;
        }
        else {
          node_gregister (me, &lni, fh->quadid[j], nodene, P4EST_CONNECT_CORNER);
          if (me->with_faces) {
            node_gregister (me, &lnh[j ^ swapi], fh->quadid[j], n_mface[face], P4EST_CONNECT_FACE);
          }
        }
      }
    }
  }

#if 0
  {
    if (me->with_faces) {
      /* one face node on same-size connection */
      nunodes = 1;
      codim[0] = 1;
      is_owned[0] = 1;
      for (i = 0; i < 2; ++i) {
        fu = &fss[i]->is.full;

        /* examine ownership situation */
        q = -1;
        if (!fu->is_ghost) {
          q = sharers[0][i] = me->mpirank;
        }
#ifdef P4EST_ENABLE_MPI
        else if ((igi = fu->quadid) >= 0) {
          P4EST_ASSERT (me->ghost != NULL);
          q = sharers[0][i] = me->ghost_rank[igi];
          gquad =
            (p4est_quadrant_t *) sc_array_index (&me->ghost->ghosts, igi);
          P4EST_ASSERT (gquad->p.piggy3.which_tree == fss[i]->treeid);
          gpos[0][i] = gquad->p.piggy3.local_num * ln->vnodes +
            pos_lnodes_face_full (fss[i]->face);
          is_shared[0] = 1;
        }
        if (q >= 0) {
          /* this side face is local or found in ghost layer */
          if (q < owner[0]) {
            is_owned[0] = 0;
            owner[0] = q;
          }
        }
#endif
      }
      if (is_owned[0]) {
        P4EST_ASSERT (owner[0] == me->mpirank);
        lni = me->num_owned++;
      }
      else {
        P4EST_ASSERT (owner[0] < me->mpirank);
        lni = -1 - me->num_shared++;
      }
      for (i = 0; i < 2; ++i) {
        if ((q = sharers[0][i]) == me->mpirank) {
          /* this is a local element */
          tree = p4est_tree_array_index (fi->p4est->trees, fss[i]->treeid);
          le = tree->quadrants_offset + fss[i]->is.full.quadid;
          set_lnodes_face_full (me, le, fss[i]->face, lni);
#ifdef P4EST_ENABLE_MPI
          if (is_shared[0]) {
            /* the current rank partakes in sharing the node */
            peer = peer_access (me, q);
            peer_add_share (peer, lni);
          }
#endif
        }
#ifdef P4EST_ENABLE_MPI
        else if (q >= 0) {
          /* this is a remote element */
          peer = peer_access (me, q);
          if (is_owned[0]) {
            P4EST_ASSERT (me->mpirank < q);
            /* add space for query from q to receive buffer and add sharing */
            peer_add_reply (peer, lni);
          }
          else if (q == owner[0]) {
            P4EST_ASSERT (q < me->mpirank);
            /* add query to send buffer to q and add sharing */
            peer_add_query (peer, gpos[0][i], lni);
          }
          else {
            P4EST_ASSERT (owner[0] < me->mpirank && owner[0] < q);
            /* no message but add shared node entry */
            peer_add_share (peer, lni);
          }
        }
#else
        else {
          SC_ABORT_NOT_REACHED ();
        }
#endif
      }
    }
    return;
  }
#endif

#if 0
  /* this is a hanging face connection */
  nunodes = 1 + (me->with_faces ? 2 : 0);
  codim[0] = P4EST_DIM;
  if (me->with_faces) {
    codim[1] = codim[2] = 1;
  }
  for (j = 0, i = 0; i < 2; ++i) {
    fs = (p4est_iter_face_side_t *) sc_array_index_int (&fi->sides, i);
    if (!fs->is_hanging) {
      /* add midface corner and possibly two half face nodes */
      fu = &fs->is.full;
      q = -1;
      if (!fu->is_ghost) {
        q = sharers[0][j] = me->mpirank;
      }
      else if (fu->quadid >= 0) {
        P4EST_ASSERT (me->ghost != NULL);
        q = sharers[0][j] = me->ghost_rank[fu->quadid];
        is_shared[0] = 1;
      }
      if (q >= 0) {
        if (me->with_faces) {
          sharers[1][j] = sharers[2][j] = q;
        }
        if (me->with_faces) {

        }
        for (j = 0; j < nunodes; ++j) {
          if (q < me->mpirank) {
            is_owned[j] = 0;
          }
        }

      }

    }

  }
#endif
}

static void
iter_corner1 (p4est_iter_corner_info_t * ci, void *user_data)
{
  tnodes_meta_t      *me = (tnodes_meta_t *) user_data;
  p4est_iter_corner_side_t *cs;

  p4est_locidx_t      le;               /**< local element number */
  p4est_locidx_t      lni;              /**< local node number */
  int                 nodene;
  size_t              zz;

  /* initial checks  */
  P4EST_ASSERT (ci->p4est == me->p4est);

  lni = -1;
  for (zz = 0; zz < ci->sides.elem_count; ++zz) {
    cs = (p4est_iter_corner_side_t *) sc_array_index (&ci->sides, zz);
    nodene = n_ccorn[cs->corner];
    if (!cs->is_ghost) {
      le = tree_quad_to_le (ci->p4est, cs->treeid, cs->quadid);
      node_lregister (me, &lni, le, nodene, P4EST_CONNECT_CORNER);
    }
    else {
      node_gregister (me, &lni, cs->quadid, nodene, P4EST_CONNECT_CORNER);
    }
  }
}

static int
cnode_compare (const void *v1, const void *v2)
{
  const tnodes_cnode_t **cc1 = (const tnodes_cnode_t **) v1;
  const tnodes_cnode_t **cc2 = (const tnodes_cnode_t **) v2;
  const tnodes_contr_t *o1, *o2;
  p4est_locidx_t ldiff;

  P4EST_ASSERT (cc1 != NULL && *cc1 != NULL);
  P4EST_ASSERT (cc2 != NULL && *cc2 != NULL);

  P4EST_ASSERT ((*cc1)->contr.elem_size == sizeof (tnodes_contr_t));
  P4EST_ASSERT ((*cc2)->contr.elem_size == sizeof (tnodes_contr_t));

  /* we sort within the same owner process */
  o1 = (*cc1)->owner;
  o2 = (*cc2)->owner;
  P4EST_ASSERT (o1->rank == o2->rank);

  /* if the elements are the same, compare node position */
  if (o1->le == o2->le) {
    return o1->nodene - o2->nodene;
  }

  /* nodes are sorted according to their element number */
  ldiff = o1->le - o2->le;
  return ldiff == 0 ? 0 : ldiff < 0 ? -1 : 1;
}

static void
sort_owned (tnodes_meta_t *me)
{
  tnodes_cnode_t     *cnode, **ccn;
  p4est_lnodes_t     *ln = me->tm->lnodes;
  p4est_gloidx_t      gc;
  const int           s = me->mpisize;
  int                 mpiret;
  int                 q;
  size_t              zz, siz;



  /* lookup nodes separately per process */
  P4EST_ASSERT (me->num_owned == 0);
  P4EST_ASSERT (me->num_shared == 0);
  siz = me->construct.elem_count;
  for (zz = 0; zz < siz; ++zz) {
    cnode = (tnodes_cnode_t *) sc_array_index (&me->construct, zz);
    P4EST_ASSERT (cnode->runid == (p4est_locidx_t) zz);
    if (cnode->owner->rank == me->mpirank) {
      ccn = (tnodes_cnode_t **) sc_array_push (&me->ownsort);
      *ccn = cnode;
      ++me->num_owned;
    }
    else {
      ++me->num_shared;
    }
  }

  /* send node number queries to remote owners */



  /* sort local node list */
  sc_array_sort (&me->ownsort, cnode_compare);

  /* share owned count */
  ln->num_local_nodes = (ln->owned_count = me->num_owned) + me->num_shared;
  ln->nonlocal_nodes = P4EST_ALLOC (p4est_gloidx_t, me->num_shared);
  ln->global_owned_count = P4EST_ALLOC (p4est_locidx_t, s);
  mpiret = sc_MPI_Allgather (&ln->owned_count, 1, P4EST_MPI_LOCIDX,
                             ln->global_owned_count, 1, P4EST_MPI_LOCIDX,
                             me->p4est->mpicomm);
  SC_CHECK_MPI (mpiret);
  me->goffset = P4EST_ALLOC (p4est_gloidx_t, s + 1);
  gc = me->goffset[0] = 0;
  for (q = 0; q < s; ++q) {
    gc = me->goffset[q + 1] = gc + ln->global_owned_count[q];
  }
  ln->global_offset = me->goffset[me->mpirank];
}

#if 0

#ifdef P4EST_ENABLE_MPI

static int
peer_compare (const void *v1, const void *v2)
{
  const tnodes_peer_t **p1 = (const tnodes_peer_t **) v1;
  const tnodes_peer_t **p2 = (const tnodes_peer_t **) v2;
  return (*p1)->rank - (*p2)->rank;
}

static void
sort_peers (tnodes_meta_t * me)
{
  int                 i;
  int                 num_peers = (int) me->peers.elem_count;
  p4est_locidx_t      nonlofs;
  tnodes_peer_t      *tp;

  /* make it possible to iterate through peers in rank order */
  sc_array_resize (&me->sortp, num_peers);
  for (i = 0; i < num_peers; ++i) {
    *(tnodes_peer_t **) sc_array_index_int (&me->sortp, i) =
      (tnodes_peer_t *) sc_array_index_int (&me->peers, i);
  }
  sc_array_sort (&me->sortp, peer_compare);
  nonlofs = 0;
  for (i = 0; i < num_peers; ++i) {
    tp = *(tnodes_peer_t **) sc_array_index_int (&me->sortp, i);
    tp->shacumul = nonlofs;
    if (tp->rank < me->mpirank) {
      nonlofs += tp->bufcount;
    }
#if 0
    P4EST_LDEBUGF ("Peer in order %d: %d count %ld offset %ld\n",
                   i, tp->rank, (long) tp->bufcount, (long) tp->shacumul);
#endif
  }
  P4EST_ASSERT (nonlofs == me->num_shared);

  /* create lookup list to find local node index of shared node */
  sc_array_resize (&me->oldtolocal, me->num_shared);
#ifdef P4EST_ENABLE_DEBUG
  sc_array_memset (&me->oldtolocal, -1);
#endif
}

static int
twop_compare (const void *v1, const void *v2)
{
  const p4est_locidx_t *p1 = (const p4est_locidx_t *) v1;
  const p4est_locidx_t *p2 = (const p4est_locidx_t *) v2;
  return p1[1] - p2[1];
}

static void
post_query_reply (tnodes_meta_t * me)
{
  int                 mpiret;
  size_t              zp, iz;
  sc_MPI_Request     *preq;
  tnodes_peer_t      *peer;

  P4EST_ASSERT (!me->peer_with_self);
  zp = me->peers.elem_count;
  sc_array_resize (&me->pereq, zp);
  for (iz = 0; iz < zp; ++iz) {
    peer = (tnodes_peer_t *) sc_array_index (&me->peers, iz);
    preq = (sc_MPI_Request *) sc_array_index (&me->pereq, iz);
    if (peer->rank == me->mpirank) {
      P4EST_ASSERT (!me->peer_with_self);
      P4EST_ASSERT (peer->bufcount == 0);
      P4EST_ASSERT (peer->querypos.elem_count == 0);
      *preq = sc_MPI_REQUEST_NULL;
      me->peer_with_self = 1;
    }
    else if (peer->rank > me->mpirank) {
      /* expecting query from higher rank */
      P4EST_ASSERT (peer->bufcount > 0);
      P4EST_ASSERT (peer->querypos.elem_count == 0);
      sc_array_resize (&peer->querypos, peer->bufcount);
      mpiret = sc_MPI_Irecv (sc_array_index (&peer->querypos, 0),
                             peer->bufcount, P4EST_MPI_LOCIDX, peer->rank,
                             P4EST_COMM_TNODES_QUERY, me->mpicomm, preq);
      SC_CHECK_MPI (mpiret);
      peer->done = 1;
    }
    else {
      /* address query to lower rank */
      P4EST_ASSERT (peer->bufcount == 0);
      P4EST_ASSERT (peer->querypos.elem_count > 0);
      peer->bufcount = (p4est_locidx_t) peer->querypos.elem_count;
      mpiret = sc_MPI_Isend (sc_array_index (&peer->querypos, 0),
                             peer->bufcount, P4EST_MPI_LOCIDX, peer->rank,
                             P4EST_COMM_TNODES_QUERY, me->mpicomm, preq);
      SC_CHECK_MPI (mpiret);
      peer->done = 3;
    }
  }
}

static void
wait_query_reply (tnodes_meta_t * me)
{
  int                 i, j;
  int                 mpiret;
  int                 nwalloc;
  int                 nwtotal;
  int                 nwaited;
  int                *waitind;
  sc_MPI_Request     *preq;
  sc_array_t          rindloc;
  p4est_locidx_t      lbc, lni, nonloc;
  p4est_locidx_t      gpos, oind, rind;
  p4est_locidx_t     *twop;
  p4est_gloidx_t      gof;
  p4est_lnodes_t     *ln = me->tm->lnodes;
  tnodes_peer_t      *peer;

  nwalloc = (int) me->peers.elem_count;
  P4EST_ASSERT (nwalloc > 0 || !me->peer_with_self);
  nwtotal = nwalloc - (me->peer_with_self ? 1 : 0);
  waitind = P4EST_ALLOC (int, nwalloc);
  while (nwtotal > 0) {
    mpiret = sc_MPI_Waitsome
      (nwalloc, (sc_MPI_Request *) sc_array_index (&me->pereq, 0),
       &nwaited, waitind, sc_MPI_STATUSES_IGNORE);
    SC_CHECK_MPI (mpiret);
    SC_CHECK_ABORT (nwaited > 0, "Invalid count after MPI_Waitsome");
    for (i = 0; i < nwaited; ++i) {
      j = waitind[i];
      peer = (tnodes_peer_t *) sc_array_index (&me->peers, j);
      P4EST_ASSERT (peer->rank != me->mpirank);
      preq = (sc_MPI_Request *) sc_array_index (&me->pereq, j);
      P4EST_ASSERT (*preq == sc_MPI_REQUEST_NULL);
      if (peer->rank > me->mpirank) {
        P4EST_ASSERT (peer->shacumul == me->num_shared);
        if (peer->done == 1) {
          /* we have received a request and shall send a reply */
          lbc = peer->bufcount;
          for (lni = 0; lni < lbc; ++lni) {
            gpos = *(p4est_locidx_t *) sc_array_index (&peer->querypos, lni);

            P4EST_LDEBUGF ("Got %d gquad %d pos %d\n from %d\n", lni,
                           gpos / ln->vnodes, gpos % ln->vnodes, peer->rank);

            P4EST_ASSERT (0 <= gpos && gpos < ln->vnodes * ln->owned_count);
            P4EST_ASSERT (pos_is_boundary[gpos % ln->vnodes]);
            oind = ln->element_nodes[gpos];
            P4EST_ASSERT (0 <= oind && oind < ln->owned_count);
            *(p4est_locidx_t *) sc_array_index (&peer->querypos, lni) = oind;
          }
          mpiret = sc_MPI_Isend (sc_array_index (&peer->querypos, 0),
                                 peer->bufcount, P4EST_MPI_LOCIDX, peer->rank,
                                 P4EST_COMM_TNODES_REPLY, me->mpicomm, preq);
          SC_CHECK_MPI (mpiret);
          peer->done = 2;
        }
        else {
          /* our reply has been received */
          P4EST_ASSERT (peer->done == 2);
          peer->done = 0;
          --nwtotal;
        }
      }
      else {
        if (peer->done == 3) {
          /* our request has been sent and we await the reply */
          mpiret = sc_MPI_Irecv (sc_array_index (&peer->querypos, 0),
                                 peer->bufcount, P4EST_MPI_LOCIDX, peer->rank,
                                 P4EST_COMM_TNODES_REPLY, me->mpicomm, preq);
          SC_CHECK_MPI (mpiret);
          peer->done = 4;
        }
        else {
          /* process information in reply received */
          P4EST_ASSERT (peer->done == 4);
          gof = me->goffset[peer->rank];
          lbc = peer->bufcount;

          /* sort old local nodes by owner's node index */
          sc_array_init_count (&rindloc, 2 * sizeof (p4est_locidx_t), lbc);
          for (lni = 0; lni < lbc; ++lni) {
            rind = *(p4est_locidx_t *) sc_array_index (&peer->localind, lni);
            P4EST_ASSERT (rind <= -1);
            rind = -rind - 1;
            P4EST_ASSERT (0 <= rind && rind < me->num_shared);
            oind = *(p4est_locidx_t *) sc_array_index (&peer->querypos, lni);
            twop = (p4est_locidx_t *) sc_array_index (&rindloc, lni);
            twop[0] = rind;
            twop[1] = oind;
          }
          sc_array_sort (&rindloc, twop_compare);

          /* store local node's global position and old local index */
          for (lni = 0; lni < lbc; ++lni) {
            twop = (p4est_locidx_t *) sc_array_index (&rindloc, lni);
            nonloc = peer->shacumul + lni;
            ln->nonlocal_nodes[nonloc] = gof + twop[1];
            *(p4est_locidx_t *) sc_array_index (&me->oldtolocal, twop[0])
              = nonloc;
          }
          sc_array_reset (&rindloc);
          peer->done = 0;
          --nwtotal;
        }
      }
    }
  }
  P4EST_FREE (waitind);
}

static void
finalize_nodes (tnodes_meta_t * me)
{
  int                 i;
  int                 num_peers = (int) me->peers.elem_count;
  p4est_locidx_t      pnum, pind, lpos;
  p4est_locidx_t      rind, nonloc;
  p4est_locidx_t      lbc, lni;
  p4est_lnodes_t     *ln = me->tm->lnodes;
  tnodes_peer_t      *peer;

  /* retrieve final numbers of nonlocal element nodes */
  pnum = (p4est_locidx_t) me->remotepos.elem_count;
  for (pind = 0; pind < pnum; ++pind) {
    lpos = *(p4est_locidx_t *) sc_array_index (&me->remotepos, pind);
    P4EST_ASSERT (0 <= lpos && lpos < ln->num_local_elements * ln->vnodes);
    rind = ln->element_nodes[lpos];
    P4EST_ASSERT (rind <= -1);
    rind = -1 - rind;
    P4EST_ASSERT (0 <= rind && rind < me->num_shared);
    nonloc = *(p4est_locidx_t *) sc_array_index (&me->oldtolocal, rind);
    P4EST_ASSERT (0 <= nonloc && nonloc < me->num_shared);
    ln->element_nodes[lpos] = ln->owned_count + nonloc;
  }

  /* sort sharer node arrays */
  for (i = 0; i < num_peers; ++i) {
    peer = *(tnodes_peer_t **) sc_array_index_int (&me->sortp, i);
    lbc = peer->sharedno.elem_count;
    for (lni = 0; lni < lbc; ++lni) {
      rind = *(p4est_locidx_t *) sc_array_index (&peer->sharedno, lni);
      if (rind > 0) {
        P4EST_ASSERT (rind < ln->owned_count);
        continue;
      }
      P4EST_ASSERT (rind <= -1);
      rind = -rind - 1;
      P4EST_ASSERT (0 <= rind && rind < me->num_shared);
      nonloc = *(p4est_locidx_t *) sc_array_index (&me->oldtolocal, rind);
      P4EST_ASSERT (0 <= nonloc && nonloc < me->num_shared);
      *(p4est_locidx_t *) sc_array_index (&peer->sharedno, lni)
        = ln->owned_count + nonloc;
    }
    sc_array_sort (&peer->sharedno, p4est_locidx_compare);
  }
}

#endif /* P4EST_ENABLE_MPI */

#endif /* 0 */

static void
clean_construct (tnodes_meta_t *me)
{
  tnodes_cnode_t    *cnode;
  size_t             zz;

  for (zz = 0; zz < me->construct.elem_count; ++zz) {
    cnode = (tnodes_cnode_t *) sc_array_index (&me->construct, zz);
    sc_array_reset (&cnode->contr);
  }
}

p4est_tnodes_t     *
p4est_tnodes_new (p4est_t * p4est, p4est_ghost_t * ghost,
                  int full_style, int with_faces)
{
  int                 q, s;
  int                 vn;
  p4est_locidx_t      lel, lg, ng;
  p4est_tnodes_t     *tm;
  p4est_lnodes_t     *ln;
  tnodes_meta_t       tmeta, *me = &tmeta;
#ifdef P4EST_ENABLE_DEBUG
  uint8_t             configure;
#if 0
  p4est_lnodes_t     *testnodes;
#endif
  p4est_locidx_t      le;
#endif
#ifdef P4EST_ENABLE_MPI
  size_t              nz, zi;
  tnodes_peer_t      *peer;
#endif

  P4EST_ASSERT (p4est_is_balanced (p4est, P4EST_CONNECT_FACE));

  /* basic assignment of members */
  memset (me, 0, sizeof (tnodes_meta_t));
  memset (me->smone, -1, 25 * sizeof (p4est_locidx_t));
  me->p4est = p4est;
  me->mpicomm = p4est->mpicomm;
  s = me->mpisize = p4est->mpisize;
  me->mpirank = p4est->mpirank;
  tm = me->tm = P4EST_ALLOC_ZERO (p4est_tnodes_t, 1);
  tm->full_style = me->full_style = full_style;
  tm->with_faces = me->with_faces = with_faces;
  ln = tm->lnodes = P4EST_ALLOC_ZERO (p4est_lnodes_t, 1);

  /* lookup structure for ghost owner rank */
  if ((me->ghost = ghost) != NULL) {
    P4EST_ASSERT (ghost->proc_offsets[0] == 0);
    P4EST_ASSERT (ghost->proc_offsets[s] ==
                  (p4est_locidx_t) ghost->ghosts.elem_count);
    me->ghost_rank = P4EST_ALLOC (int, ghost->ghosts.elem_count);
    lg = 0;
    for (q = 0; q < s; ++q) {
      ng = ghost->proc_offsets[q + 1];
      for (; lg < ng; ++lg) {
        me->ghost_rank[lg] = q;
      }
    }
    P4EST_ASSERT (lg == (p4est_locidx_t) ghost->ghosts.elem_count);
#ifdef P4EST_ENABLE_MPI
    me->proc_peer = P4EST_ALLOC_ZERO (int, s);
    sc_array_init (&me->remotepos, sizeof (p4est_locidx_t));
    sc_array_init (&me->sortp, sizeof (tnodes_peer_t *));
    sc_array_init (&me->peers, sizeof (tnodes_peer_t));
    sc_array_init (&me->pereq, sizeof (sc_MPI_Request));
    sc_array_init (&me->oldtolocal, sizeof (p4est_locidx_t));
#endif
  }
  sc_array_init (&me->construct, sizeof (tnodes_cnode_t));
  sc_array_init (&me->ownsort, sizeof (tnodes_cnode_t *));

  /* prepare node information */
  ln->mpicomm = p4est->mpicomm;
  ln->sharers = sc_array_new (sizeof (p4est_lnodes_rank_t));
  ln->degree = 0;
  vn = ln->vnodes = 9 + (with_faces ? 16 : 0);
  lel = ln->num_local_elements = p4est->local_num_quadrants;
  P4EST_ASSERT ((size_t) lel * (size_t) vn <= (size_t) P4EST_LOCIDX_MAX);
  me->chilev = P4EST_ALLOC_ZERO (uint8_t, lel);
  ln->face_code = P4EST_ALLOC_ZERO (p4est_lnodes_code_t, lel);
  ln->element_nodes = P4EST_ALLOC (p4est_locidx_t, lel * vn);
  memset (ln->element_nodes, -1, lel * vn * sizeof (p4est_locidx_t));

  /* allocate arrays for node encoding */
  tm->configuration = P4EST_ALLOC_ZERO (uint8_t, lel);
  tm->local_toffset = P4EST_ALLOC (p4est_locidx_t, lel + 1);
  tm->global_toffset = P4EST_ALLOC (p4est_gloidx_t, s + 1);

  /* determine triangle configuration of each element */
  me->lenum = 0;
  p4est_iterate (p4est, ghost, me, iter_volume1, iter_face1, iter_corner1);
  P4EST_ASSERT (me->lenum == lel);
  sort_owned (me);
  P4EST_INFOF ("p4est_tnodes_new: owned %ld shared %ld\n",
               (long) me->num_owned, (long) me->num_shared);
  P4EST_GLOBAL_PRODUCTIONF ("p4est_tnodes_new: global owned %lld\n",
                            (long long) me->goffset[s]);

#if 0
#ifdef P4EST_ENABLE_MPI
  /* post messages */
  post_query_reply (me);

  /* sort peers by process */
  sort_peers (me);

  /* receive messages */
  wait_query_reply (me);

  /* finalize element nodes */
  finalize_nodes (me);
#endif
#endif

  /* free memory */
  P4EST_FREE (me->goffset);
  if (me->ghost != NULL) {
#ifdef P4EST_ENABLE_MPI
    nz = me->peers.elem_count;
    for (zi = 0; zi < nz; ++zi) {
      peer = (tnodes_peer_t *) sc_array_index (&me->peers, zi);
      P4EST_ASSERT (!peer->done);
      sc_array_reset (&peer->sharedno);
      sc_array_reset (&peer->localind);
      sc_array_reset (&peer->querypos);
    }
    sc_array_reset (&me->remotepos);
    sc_array_reset (&me->sortp);
    sc_array_reset (&me->peers);
    sc_array_reset (&me->pereq);
    sc_array_reset (&me->oldtolocal);
    P4EST_FREE (me->proc_peer);
#endif
    P4EST_FREE (me->ghost_rank);
  }
  clean_construct (me);
  sc_array_reset (&me->construct);
  sc_array_reset (&me->ownsort);
  P4EST_FREE (me->chilev);

#ifdef P4EST_ENABLE_DEBUG
  for (le = 0; le < lel; ++le) {
    configure = tm->configuration[le];
#if 0
    P4EST_LDEBUGF ("Configuration %ld %ld is %u\n", (long) le, (long) lel, (unsigned) configure);
#endif
    P4EST_ASSERT (configure <= 16 || configure == 32);
  }
  if (me->ghost != NULL) {
#if 0
    testnodes = p4est_lnodes_new (p4est, me->ghost, 2);
    P4EST_ASSERT (testnodes->num_local_elements == lel);
    for (le = 0; le < lel; ++le) {
      P4EST_ASSERT (testnodes->face_code[le] == ln->face_code[le]);
    }
    p4est_lnodes_destroy (testnodes);
#endif
  }
#endif

  return tm;
}

void
p4est_tnodes_destroy (p4est_tnodes_t * tm)
{
  P4EST_ASSERT (tm != NULL);
  P4EST_ASSERT (tm->lnodes != NULL);

  p4est_lnodes_destroy (tm->lnodes);
  P4EST_FREE (tm->global_toffset);
  P4EST_FREE (tm->local_toffset);
  P4EST_FREE (tm->configuration);
  P4EST_FREE (tm);
}