/*
  This file is part of p4est.
  p4est is a C library to manage a collection (a forest) of multiple
  connected adaptive quadtrees or octrees in parallel.

  Copyright (C) 2011 The University of Texas System
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
#include <p8est_bits.h>
#include <p8est_communication.h>
#include <p8est_search.h>
#else
#include <p4est_bits.h>
#include <p4est_communication.h>
#include <p4est_search.h>
#endif /* !P4_TO_P8 */

/* a quadrant q is outside of quadrant p, as a descendent of p's neighbor n
 * outside of \a corner.  We know the \a child of n that contains q.
 * Depending on the \a balance type, q is either relevant (return true) or
 * irrelevant (return false) to the balance of p. */
#ifndef P4_TO_P8
static inline int
p4est_balance_corner_relevant (int balance, int corner, int child)
{
  /* if balance == 1 (i.e. face and corner balance), then every descendent of
   * n is relevant; if balance == 0 (face balance), the child of n farthest
   * from p's corner (i.e. that has the same number as p's corner) is the only
   * child that is not relevant) */

  P4EST_ASSERT (0 <= corner && corner < P4EST_CHILDREN);
  P4EST_ASSERT (0 <= child && child < P4EST_CHILDREN);
  return (balance || corner != child);
}
#else

#define p4est_balance_corner_relevant p8est_balance_corner_relevant
/* assuming face only balance \a corner == 0, only the closest child of n (child 7) and its face neighbors (3, 5, 6) are relevant to p's balance */
static const int    p8est_balance_corner_face_relevant[8] =
  { 0, 0, 0, 1, 0, 1, 1, 1 };

static inline int
p8est_balance_corner_relevant (int balance, int corner, int child)
{
  /* if balance == 2 (i.e. face, edge and corner balance), then every
   * descendent of n is relevant; if balance == 1 (face and edge balance),
   * the child of n farthest from p's corner (i.e. that has the same number as
   * p's corner) is the only child that is not relevant); if balance == 0
   * (face only balance), we transform \a child as though \a corner == 0 and
   * use the table above */

  P4EST_ASSERT (0 <= corner && corner < P4EST_CHILDREN);
  P4EST_ASSERT (0 <= child && child < P4EST_CHILDREN);
  return (balance == 2 ||
          (balance == 1 && child != corner) ||
          p8est_balance_corner_face_relevant[child ^ corner]);
}

/* a quadrant q is outside of quadrant p, as a descendent of p's neighbor n
 * outside of \a edge.  We know the \a child of n that contains q.
 * Depending on the \a balance type, q is either relevant (return true) or
 * irrelevant (return false) to the balance of p. */

static inline int
p8est_balance_edge_relevant (int balance, int edge, int child)
{
  /* if there is edge balancing (balance > 0), then every descendent of n is
   * relevant; if there is only face balancing, \a child is not relevant if it
   * touches the farthest edge of n away from p */

  P4EST_ASSERT (0 <= edge && edge < P8EST_EDGES);
  P4EST_ASSERT (0 <= child && child < P4EST_CHILDREN);
  return (balance || p8est_corner_edges[child][edge >> 2] != edge);
}
#endif

static void         p4est_bal_corner_con_internal (p4est_quadrant_t const
                                                   *restrict q,
                                                   p4est_quadrant_t *
                                                   restrict p, int corner,
                                                   int balance,
                                                   int *consistent,
                                                   int recurse);

static void         p4est_bal_face_con_internal (p4est_quadrant_t const
                                                 *restrict q,
                                                 p4est_quadrant_t *
                                                 restrict p, int face,
                                                 int balance, int *consistent,
                                                 int *conextra, int recurse);

#ifdef P4_TO_P8
static void         p8est_bal_edge_con_internal (p4est_quadrant_t const
                                                 *restrict q,
                                                 p4est_quadrant_t *
                                                 restrict p, int edge,
                                                 int balance, int *consistent,
                                                 int *conextra, int recurse);
#endif

/* the recursive call of p4est_bal_corner_con_internal: we have determined
 * that q causes p to split: this means that they are not consistent, so if \a
 * consistent is given, it is set to false.  If we are to \a recurse, we set p
 * to be its own child in \a corner and recur */
#ifndef P4_TO_P8
#define p4est_bal_corner_recurse(q, p, c, b, cn, r)                    \
do {                                                                   \
  if ((cn) != NULL) {                                                  \
    *(cn) = 0;                                                         \
  }                                                                    \
  if (r) {                                                             \
    (p)->level++;                                                      \
    (p)->x += (((c) & 1) ? P4EST_QUADRANT_LEN ((p)->level) : 0);       \
    (p)->y += (((c) & 2) ? P4EST_QUADRANT_LEN ((p)->level) : 0);       \
    P4EST_ASSERT (p4est_quadrant_is_extended (p));                     \
    p4est_bal_corner_con_internal (q, p, c, b, NULL, r);               \
  }                                                                    \
} while (0)
#else
#define p4est_bal_corner_recurse(q, p, c, b, cn, r)                    \
do {                                                                   \
  if ((cn) != NULL) {                                                  \
    *(cn) = 0;                                                         \
  }                                                                    \
  if (r) {                                                             \
    (p)->level++;                                                      \
    (p)->x += (((c) & 1) ? P4EST_QUADRANT_LEN ((p)->level) : 0);       \
    (p)->y += (((c) & 2) ? P4EST_QUADRANT_LEN ((p)->level) : 0);       \
    (p)->z += (((c) & 4) ? P4EST_QUADRANT_LEN ((p)->level) : 0);       \
    P4EST_ASSERT (p4est_quadrant_is_extended (p));                     \
    p4est_bal_corner_con_internal (q, p, c, b, NULL, r);               \
  }                                                                    \
} while (0)
#endif

#ifndef P4_TO_P8
#define p4est_bal_corner_no_recurse(q, p, c, b, cn, r)                 \
do {                                                                   \
  if ((cn) != NULL) {                                                  \
    *(cn) = 0;                                                         \
  }                                                                    \
  if (r) {                                                             \
    (p)->level++;                                                      \
    (p)->x += (((c) & 1) ? P4EST_QUADRANT_LEN ((p)->level) : 0);       \
    (p)->y += (((c) & 2) ? P4EST_QUADRANT_LEN ((p)->level) : 0);       \
    P4EST_ASSERT (p4est_quadrant_is_extended (p));                     \
  }                                                                    \
} while (0)
#else
#define p4est_bal_corner_no_recurse(q, p, c, b, cn, r)                 \
do {                                                                   \
  if ((cn) != NULL) {                                                  \
    *(cn) = 0;                                                         \
  }                                                                    \
  if (r) {                                                             \
    (p)->level++;                                                      \
    (p)->x += (((c) & 1) ? P4EST_QUADRANT_LEN ((p)->level) : 0);       \
    (p)->y += (((c) & 2) ? P4EST_QUADRANT_LEN ((p)->level) : 0);       \
    (p)->z += (((c) & 4) ? P4EST_QUADRANT_LEN ((p)->level) : 0);       \
    P4EST_ASSERT (p4est_quadrant_is_extended (p));                     \
  }                                                                    \
} while (0)
#endif

/* \a corner is the corner of \a p closest to \a q: the corner of \a q closest
 * to \a p is the dual of \a corner */
static void
p4est_bal_corner_con_internal (p4est_quadrant_t const *restrict q,
                               p4est_quadrant_t * restrict p,
                               int corner, int balance,
                               int *consistent, int recurse)
{
  int                 level = p->level + 1;     /* level is now the size of a
                                                   child */
  int                 child;
  int                 dual = corner ^ (P4EST_CHILDREN - 1);
  int                 relative;
  int                 face;
#ifdef P4_TO_P8
  int                 edge;
  int                 grandchild;
#endif
  int                 recon;
  p4est_quadrant_t    temp;
  p4est_quadrant_t    a;

#ifdef P4EST_DEBUG
  P4EST_ASSERT (p4est_quadrant_is_valid (q));
  P4EST_ASSERT (p4est_quadrant_is_extended (p));
  P4EST_ASSERT (q->level >= p->level);
  /* get the ancestor of q at the same level as p */
  a = *q;
  a.x &= (-1 << (P4EST_MAXLEVEL - p->level));
  a.y &= (-1 << (P4EST_MAXLEVEL - p->level));
#ifdef P4_TO_P8
  a.z &= (-1 << (P4EST_MAXLEVEL - p->level));
#endif
  a.level = p->level;
  /* check that q really is a descendent of the correct neighbor of p */
  p4est_quadrant_corner_neighbor (&a, dual, &temp);
  P4EST_ASSERT (p4est_quadrant_is_equal (&temp, p));
#endif

  /* start with the assumption that the two are consistent */
  if (consistent != NULL) {
    *consistent = 1;
  }

  /* if q is larger than or equal to half of p,
   * or p is a smallest quadrant, then of
   * course q and p are balance consistent */
  if (q->level <= level || level > P4EST_QMAXLEVEL) {
    return;
  }

  /* early exit if q is in an irrelevant descendent of n */
  child = p4est_quadrant_ancestor_id (q, level);
  if (!p4est_balance_corner_relevant (balance, corner, child)) {
    return;
  }

  /* if q is in the child of n closest to p, i.e. directly outside \a corner
   * */
  if (child == dual) {
    if (balance == (P4EST_DIM - 1)) {
      /* q is smaller than a child of n, which means that a quadrant smaller
       * than a child of p is directly outside of \a corner, which means
       * that q and p are not balance compatible */
      p4est_bal_corner_recurse (q, p, corner, balance, consistent, recurse);
      return;
    }
    if (++level == q->level) {  /* level is now the size of a grandchild */
      /* in all of the remaining cases, a grandchild is consistent */
      return;
    }
#ifdef P4_TO_P8
    if (balance == 1) {         /* edge balance */
#endif
      do {
        child = p4est_quadrant_ancestor_id (q, level);
        if (child != corner) {
          p4est_bal_corner_recurse (q, p, corner, balance, consistent,
                                    recurse);
          return;
        }
      } while (++level < q->level);
      return;
#ifdef P4_TO_P8
    }

    P4EST_ASSERT (balance == 0);
    P4EST_ASSERT (q->level > level);
    /* this is the trickiest case: outside a corner, 3D, face only balance */
    child = p4est_quadrant_ancestor_id (q, level);
    relative = child ^ dual;
    switch (relative) {
    case 7:
      /* if we are in the (grand)child farthest from p, we are consistent */
      return;
    case 0:
      /* if we are in the (grand)child closest to p, it looks like the edge
       * balance case */
      while (++level < q->level) {
        child = p4est_quadrant_ancestor_id (q, level);
        if (child != corner) {
          p4est_bal_corner_recurse (q, p, corner, balance, consistent,
                                    recurse);
          return;
        }
      }
      return;
    case 1:
    case 2:
    case 4:
      /* we are in a (grand)child that is face adjacent to the closest
       * (grand)child to p, which looks like trying to balance
       * outside an edge */
      edge = p8est_corner_edges[corner][relative >> 1];
      while (++level < q->level) {
        child = p4est_quadrant_ancestor_id (q, level);
        if (p8est_balance_edge_relevant (balance, edge, child)) {
          p4est_bal_corner_recurse (q, p, corner, balance, consistent,
                                    recurse);
          return;
        }
      }
      return;
    case 3:
    case 5:
    case 6:
      /* we are in a (grand)child that is edge adjacent to the closest
       * (grand)child to p.  this is the bizarre fractal case */
      while (++level < q->level) {
        grandchild = p4est_quadrant_ancestor_id (q, level);
        if (grandchild == child || grandchild == (child ^ 7)) {
          /* this is the self-similar case */
          continue;
        }
        relative = grandchild ^ dual;
        switch (relative) {
        case 7:
          /* consistent */
          return;
        case 0:
          while (++level < q->level) {
            /* this looks like the outside corner, balance edge case */
            grandchild = p4est_quadrant_ancestor_id (q, level);
            if (grandchild != corner) {
              p4est_bal_corner_recurse (q, p, corner, balance, consistent,
                                        recurse);
              return;
            }
          }
          return;
        case 1:
        case 2:
        case 4:
          /* this case looks like the outside edge, balance face case */
          edge = p8est_corner_edges[corner][relative >> 1];
          while (++level < q->level) {
            grandchild = p4est_quadrant_ancestor_id (q, level);
            if (p8est_balance_edge_relevant (balance, edge, grandchild)) {
              p4est_bal_corner_recurse (q, p, corner, balance, consistent,
                                        recurse);
              return;
            }
          }
          return;
        case 3:
        case 5:
        case 6:
          /* this case looks sort of like the outside edge,
           * balance face case */
          edge = p8est_corner_edges[corner][((child ^ grandchild) ^ 7) >> 1];
          while (++level < q->level) {
            grandchild = p4est_quadrant_ancestor_id (q, level);
            if (!p8est_balance_edge_relevant (balance, edge, grandchild)) {
              /* consistent */
              return;
            }
            if (!p8est_balance_edge_relevant (balance, edge ^ 3, grandchild)) {
              while (++level < q->level) {
                grandchild = p4est_quadrant_ancestor_id (q, level);
                if (p8est_balance_edge_relevant (balance, edge, grandchild)) {
                  p4est_bal_corner_recurse (q, p, corner, balance, consistent,
                                            recurse);
                  return;
                }
              }
              return;
            }
          }
          return;
        default:
          SC_ABORT_NOT_REACHED ();
        }
      }
      return;
    default:
      SC_ABORT_NOT_REACHED ();
    }
#endif
  }

  /* set a to be the child of p's neighbor that contains q */
  a.x = (q->x & (-1 << (P4EST_MAXLEVEL - level)));
  a.y = (q->y & (-1 << (P4EST_MAXLEVEL - level)));
#ifdef P4_TO_P8
  a.z = (q->z & (-1 << (P4EST_MAXLEVEL - level)));
#endif
  a.level = level;
  P4EST_ASSERT (p4est_quadrant_is_ancestor (&a, q));

  if (balance == (P4EST_DIM - 1)) {

    relative = child ^ dual;
    switch (relative) {
    case (P4EST_CHILDREN - 1):
      p4est_quadrant_corner_neighbor (&a, dual, &temp);
      p4est_bal_corner_con_internal (q, &temp, corner, balance, &recon, 0);
      if (!recon) {
        p4est_bal_corner_no_recurse (q, p, corner, balance, consistent,
                                     recurse);
        return;
      }
      return;
    case 1:
    case 2:
#ifdef P4_TO_P8
    case 4:
#endif
      face = p4est_corner_faces[corner][relative >> 1];
      p4est_quadrant_face_neighbor (&a, face ^ 1, &temp);
      p4est_bal_face_con_internal (q, &temp, face, balance, &recon, NULL, 0);
      if (!recon) {
        p4est_bal_corner_no_recurse (q, p, corner, balance, consistent,
                                     recurse);
        return;
      }
      return;
#ifdef P4_TO_P8
    case 3:
    case 5:
    case 6:
      edge = p8est_corner_edges[corner][(relative ^ 7) >> 1];
      p8est_quadrant_edge_neighbor (&a, edge ^ 3, &temp);
      p8est_bal_edge_con_internal (q, &temp, edge, balance, &recon, NULL, 0);
      if (!recon) {
        p4est_bal_corner_no_recurse (q, p, corner, balance, consistent,
                                     recurse);
        return;
      }
      return;
#endif
    default:
      SC_ABORT_NOT_REACHED ();
    }
  }

#ifndef P4_TO_P8
  P4EST_ASSERT (!balance);
#else
  if (!balance) {
#endif
    p4est_quadrant_corner_neighbor (&a, dual, &temp);
    p4est_bal_corner_con_internal (q, &temp, corner, balance, &recon, 0);
    if (!recon) {
      p4est_bal_corner_no_recurse (q, p, corner, balance, consistent,
                                   recurse);
      return;
    }
    return;
#ifdef P4_TO_P8
  }

  P4EST_ASSERT (balance == 1);

  relative = child ^ dual;
  switch (relative) {
  case 1:
  case 2:
  case 4:
    relative >>= 1;
    edge = p8est_corner_edges[corner][(relative + 1) % 3];
    p8est_quadrant_edge_neighbor (&a, edge ^ 3, &temp);
    p8est_bal_edge_con_internal (q, &temp, edge, balance, &recon, NULL, 0);
    if (!recon) {
      p4est_bal_corner_no_recurse (q, p, corner, balance, consistent,
                                   recurse);
      return;
    }
    edge = p8est_corner_edges[corner][(relative + 2) % 3];
    p8est_quadrant_edge_neighbor (&a, edge ^ 3, &temp);
    p8est_bal_edge_con_internal (q, &temp, edge, balance, &recon, NULL, 0);
    if (!recon) {
      p4est_bal_corner_no_recurse (q, p, corner, balance, consistent,
                                   recurse);
      return;
    }
    return;
  case 3:
  case 5:
  case 6:
    p4est_quadrant_corner_neighbor (&a, dual, &temp);
    p4est_bal_corner_con_internal (q, &temp, corner, balance, &recon, 0);
    if (!recon) {
      p4est_bal_corner_no_recurse (q, p, corner, balance, consistent,
                                   recurse);
      return;
    }
    return;
  default:
    SC_ABORT_NOT_REACHED ();
  }
#endif
}

static void
p4est_bal_face_con_internal (p4est_quadrant_t const *restrict q,
                             p4est_quadrant_t * restrict p, int face,
                             int balance, int *consistent, int *conextra,
                             int recurse)
{
  int                 level = p->level + 1;     /* level is now the size of a
                                                   child */
  int                 child;
  int                 dual = face ^ 1;
#ifdef P4_TO_P8
  int                 edge;
#endif
  int                 recon;
  p4est_quadrant_t    temp;
  p4est_quadrant_t    a;
  int                 cfc;
  int                 cfcextra;
#ifndef P4_TO_P8
  int                 nconextra = 3;
#else
  int                 nconextra = 9;
  int                 reconextra[3];
  int                 i;
#endif
  int                 extraid;

#ifdef P4EST_DEBUG
  P4EST_ASSERT (p4est_quadrant_is_valid (q));
  P4EST_ASSERT (p4est_quadrant_is_extended (p));
  P4EST_ASSERT (q->level >= p->level);
  /* get the ancestor of q at the same level as p */
  a = *q;
  a.x &= (-1 << (P4EST_MAXLEVEL - p->level));
  a.y &= (-1 << (P4EST_MAXLEVEL - p->level));
#ifdef P4_TO_P8
  a.z &= (-1 << (P4EST_MAXLEVEL - p->level));
#endif
  a.level = p->level;

  /* check that q really is a descendent of the correct neighbor of p */
  p4est_quadrant_face_neighbor (&a, dual, &temp);
  P4EST_ASSERT (p4est_quadrant_is_equal (&temp, p));
#endif

  /* start with the assumption that the two are consistent */
  if (consistent != NULL) {
    *consistent = 1;
    if (conextra != NULL) {
      conextra[nconextra / 2] = p->level;
    }
  }

  /* if q is larger than half of p, or p is a smallest quadrant, then of
   * course q and p are balance consistent */
  if (q->level <= level || level > P4EST_QMAXLEVEL) {
    return;
  }

  child = p4est_quadrant_ancestor_id (q, level);
  cfc = p4est_corner_face_corners[child][dual];
  if (cfc >= 0) {
    if (consistent != NULL) {
      *consistent = 0;
    }
    if (recurse) {
      child = p4est_face_corners[face][cfc];
      p->level++;
      p->x += ((child & 1) ? P4EST_QUADRANT_LEN (p->level) : 0);
      p->y += ((child & 2) ? P4EST_QUADRANT_LEN (p->level) : 0);
#ifdef P4_TO_P8
      p->z += ((child & 4) ? P4EST_QUADRANT_LEN (p->level) : 0);
#endif
      P4EST_ASSERT (p4est_quadrant_is_extended (p));

      if (conextra != NULL) {
        conextra[nconextra / 2] = p->level;
#ifndef P4_TO_P8
        cfcextra = cfc ^ 1;
        extraid = (cfcextra < cfc) ? 0 : 2;
#else
        cfcextra = cfc ^ 3;
        extraid = ((cfcextra & 1) < (cfc & 1)) ? 0 : 2;
        extraid += ((cfcextra & 2) < (cfc & 2)) ? 0 : 6;
#endif

        p4est_quadrant_sibling (p, &temp, p4est_face_corners[face][cfcextra]);
        p4est_bal_corner_con_internal (q, &temp, child, balance, &recon,
                                       recurse);
        if (temp.level > p->level) {
          conextra[extraid] = temp.level;
        }

#ifdef P4_TO_P8
        cfcextra = cfc ^ 1;
        extraid = 3 + ((cfcextra < cfc) ? 0 : 2);

        p4est_quadrant_sibling (p, &temp, p4est_face_corners[face][cfcextra]);
        edge = p8est_face_edges[face][2 + (cfc & 1)];
        reconextra[0] = reconextra[1] = reconextra[2] = p->level;
        p8est_bal_edge_con_internal (q, &temp, edge, balance, &recon,
                                     reconextra, recurse);

        for (i = 0; i < 3; i++) {
          if (reconextra[i] > p->level) {
            conextra[extraid + (i - 1) * 3] = reconextra[i];
          }
        }

        cfcextra = cfc ^ 2;
        extraid = 1 + ((cfcextra < cfc) ? 0 : 6);

        p4est_quadrant_sibling (p, &temp, p4est_face_corners[face][cfcextra]);
        edge = p8est_face_edges[face][((cfc & 2) >> 1)];
        reconextra[0] = reconextra[1] = reconextra[2] = p->level;
        p8est_bal_edge_con_internal (q, &temp, edge, balance, &recon,
                                     reconextra, recurse);

        for (i = 0; i < 3; i++) {
          if (reconextra[i] > p->level) {
            conextra[extraid + (i - 1)] = reconextra[i];
          }
        }
#endif
      }

      p4est_bal_face_con_internal (q, p, face, balance, NULL, conextra,
                                   recurse);
    }
    return;
  }

  /* set a to be the child of p's neighbor that contains q */
  a.x = (q->x & (-1 << (P4EST_MAXLEVEL - level)));
  a.y = (q->y & (-1 << (P4EST_MAXLEVEL - level)));
#ifdef P4_TO_P8
  a.z = (q->z & (-1 << (P4EST_MAXLEVEL - level)));
#endif
  a.level = level;
  P4EST_ASSERT (p4est_quadrant_is_ancestor (&a, q));

  p4est_quadrant_face_neighbor (&a, dual, &temp);
  p4est_bal_face_con_internal (q, &temp, face, balance, &recon, NULL, 0);
  if (!recon) {
    if (consistent != NULL) {
      *consistent = 0;
    }
    if (recurse) {
      p->level++;
      p->x += ((child & 1) ? P4EST_QUADRANT_LEN (p->level) : 0);
      p->y += ((child & 2) ? P4EST_QUADRANT_LEN (p->level) : 0);
#ifdef P4_TO_P8
      p->z += ((child & 4) ? P4EST_QUADRANT_LEN (p->level) : 0);
#endif
      P4EST_ASSERT (p4est_quadrant_is_extended (p));

      if (conextra != NULL) {
        conextra[nconextra / 2] = p->level;
      }
    }
    return;
  }
}

#ifdef P4_TO_P8
static void
p8est_bal_edge_con_internal (p4est_quadrant_t const *restrict q,
                             p4est_quadrant_t * restrict p, int edge,
                             int balance, int *consistent, int *conextra,
                             int recurse)
{
  int                 level = p->level + 1;     /* level is now the size of a
                                                   child */
  int                 child;
  int                 child2;
  int                 cec;
  int                 cecextra;
  int                 dual = edge ^ 3;
  int                 recon;
  p4est_quadrant_t    temp;
  p4est_quadrant_t    a;
  int                 extraid;
  int                 face;
  int                 relative;

#ifdef P4EST_DEBUG
  P4EST_ASSERT (p4est_quadrant_is_valid (q));
  P4EST_ASSERT (p4est_quadrant_is_extended (p));
  P4EST_ASSERT (q->level >= p->level);
  /* get the ancestor of q at the same level as p */
  a = *q;
  a.x &= (-1 << (P4EST_MAXLEVEL - p->level));
  a.y &= (-1 << (P4EST_MAXLEVEL - p->level));
#ifdef P4_TO_P8
  a.z &= (-1 << (P4EST_MAXLEVEL - p->level));
#endif
  a.level = p->level;

  /* check that q really is a descendent of the correct neighbor of p */
  p8est_quadrant_edge_neighbor (&a, dual, &temp);
  P4EST_ASSERT (p4est_quadrant_is_equal (&temp, p));
#endif

  /* start with the assumption that the two are consistent */
  if (consistent != NULL) {
    *consistent = 1;
    if (conextra != NULL) {
      conextra[1] = p->level;
    }
  }

  /* if q is larger than half of p, or p is a smallest quadrant, then of
   * course q and p are balance consistent */
  if (q->level <= level || level > P4EST_QMAXLEVEL) {
    return;
  }

  child = p4est_quadrant_ancestor_id (q, level);
  if (!p8est_balance_edge_relevant (balance, edge, child)) {
    return;
  }

  cec = -1;
  if (p8est_edge_corners[dual][0] == child) {
    cec = 0;
  }
  else if (p8est_edge_corners[dual][1] == child) {
    cec = 1;
  }

  if (cec >= 0) {
    if (balance) {
      if (consistent != NULL) {
        *consistent = 0;
      }
      if (recurse) {
        child = p8est_edge_corners[edge][cec];
        p->level++;
        p->x += ((child & 1) ? P4EST_QUADRANT_LEN ((p)->level) : 0);
        p->y += ((child & 2) ? P4EST_QUADRANT_LEN ((p)->level) : 0);
        p->z += ((child & 4) ? P4EST_QUADRANT_LEN ((p)->level) : 0);
        P4EST_ASSERT (p4est_quadrant_is_extended (p));

        if (conextra != NULL) {
          conextra[1] = p->level;

          cecextra = cec ^ 1;
          extraid = (cecextra < cec) ? 0 : 2;

          p4est_quadrant_sibling (p, &temp,
                                  p8est_edge_corners[edge][cecextra]);
          p4est_bal_corner_con_internal (q, &temp, child, balance, &recon,
                                         recurse);

          if (temp.level > p->level) {
            conextra[extraid] = temp.level;
          }
        }

        p8est_bal_edge_con_internal (q, p, edge, balance, NULL, conextra,
                                     recurse);
      }
      return;
    }

    while (++level < q->level) {
      child = p4est_quadrant_ancestor_id (q, level);
      if (p8est_balance_edge_relevant (balance, edge, child)) {
        if (consistent != NULL) {
          *consistent = 0;
        }
        if (recurse) {
          child = p8est_edge_corners[edge][cec];
          p->level++;
          p->x += ((child & 1) ? P4EST_QUADRANT_LEN ((p)->level) : 0);
          p->y += ((child & 2) ? P4EST_QUADRANT_LEN ((p)->level) : 0);
          p->z += ((child & 4) ? P4EST_QUADRANT_LEN ((p)->level) : 0);
          P4EST_ASSERT (p4est_quadrant_is_extended (p));

          if (conextra != NULL) {
            conextra[1] = p->level;

            cecextra = cec ^ 1;
            extraid = (cecextra < cec) ? 0 : 2;

            p4est_quadrant_sibling (p, &temp,
                                    p8est_edge_corners[edge][cecextra]);
            p4est_bal_corner_con_internal (q, &temp, child, balance, &recon,
                                           recurse);

            if (temp.level > p->level) {
              conextra[extraid] = temp.level;
            }
          }

          p8est_bal_edge_con_internal (q, p, edge, balance, NULL, conextra,
                                       recurse);
        }
        return;
      }
    }
    return;
  }

  /* set a to be the child of p's neighbor that contains q */
  a.x = (q->x & (-1 << (P4EST_MAXLEVEL - level)));
  a.y = (q->y & (-1 << (P4EST_MAXLEVEL - level)));
  a.z = (q->z & (-1 << (P4EST_MAXLEVEL - level)));
  a.level = level;
  P4EST_ASSERT (p4est_quadrant_is_ancestor (&a, q));

  if (!balance) {
    p8est_quadrant_edge_neighbor (&a, dual, &temp);
    p8est_bal_edge_con_internal (q, &temp, edge, balance, &recon, NULL, 0);
    if (!recon) {
      if (consistent != NULL) {
        *consistent = 0;
      }
      if (recurse) {
        p->level++;
        p->x += ((child & 1) ? P4EST_QUADRANT_LEN (p->level) : 0);
        p->y += ((child & 2) ? P4EST_QUADRANT_LEN (p->level) : 0);
        p->z += ((child & 4) ? P4EST_QUADRANT_LEN (p->level) : 0);
        P4EST_ASSERT (p4est_quadrant_is_extended (p));

        if (conextra != NULL) {
          conextra[1] = p->level;
        }
      }
      return;
    }
    return;
  }

  if (p8est_edge_corners[p8est_corner_edges[child][edge / 4]][0] == child) {
    cec = 0;
  }
  else {
    cec = 1;
  }
  child2 = p8est_edge_corners[dual][cec];
  face = p8est_corner_faces[child][edge / 4];
  P4EST_ASSERT (face == p8est_corner_faces[child2][edge / 4]);

  relative = (p8est_corner_face_corners[child][face] ^
              p8est_corner_face_corners[child2][face]);

  switch (relative) {
  case 3:
    p8est_quadrant_edge_neighbor (&a, dual, &temp);
    p8est_bal_edge_con_internal (q, &temp, edge, balance, &recon, NULL, 0);
    if (!recon) {
      if (consistent != NULL) {
        *consistent = 0;
      }
      if (recurse) {
        p->level++;
        p->x += ((child & 1) ? P4EST_QUADRANT_LEN (p->level) : 0);
        p->y += ((child & 2) ? P4EST_QUADRANT_LEN (p->level) : 0);
        p->z += ((child & 4) ? P4EST_QUADRANT_LEN (p->level) : 0);
        P4EST_ASSERT (p4est_quadrant_is_extended (p));

        if (conextra != NULL) {
          conextra[1] = p->level;
        }
      }
      return;
    }
    return;
  case 1:
  case 2:
    face = p8est_edge_faces[edge][relative >> 1];
    p8est_quadrant_face_neighbor (&a, face ^ 1, &temp);
    p4est_bal_face_con_internal (q, &temp, face, balance, &recon, NULL, 0);
    if (!recon) {
      if (consistent != NULL) {
        *consistent = 0;
      }
      if (recurse) {
        p->level++;
        p->x += ((child & 1) ? P4EST_QUADRANT_LEN (p->level) : 0);
        p->y += ((child & 2) ? P4EST_QUADRANT_LEN (p->level) : 0);
        p->z += ((child & 4) ? P4EST_QUADRANT_LEN (p->level) : 0);
        P4EST_ASSERT (p4est_quadrant_is_extended (p));

        if (conextra != NULL) {
          conextra[1] = p->level;
        }
      }
      return;
    }
    return;
  default:
    SC_ABORT_NOT_REACHED ();
  }
}
#endif

int
p4est_balance_face_test (p4est_quadrant_t * restrict q,
                         p4est_quadrant_t * restrict p,
                         int face, p4est_balance_type_t balance,
                         sc_array_t * seeds)
{
  p4est_quadrant_t    temp = *p;
  p4est_quadrant_t   *s;
  int                 ibalance;
  int                 consistent;
#ifndef P4_TO_P8
  int                 nextra = 3;
  int                 conextra[3];
#else
  int                 nextra = 9;
  int                 conextra[9];
#endif
  int                 i;

  P4EST_ASSERT (seeds == NULL ||
                seeds->elem_size == sizeof (p4est_quadrant_t));

  if (balance == P4EST_BALANCE_FULL) {
    ibalance = P4EST_DIM - 1;
  }
#ifdef P4_TO_P8
  else if (balance == P8EST_BALANCE_EDGE) {
    ibalance = 1;
  }
#endif
  else {
    ibalance = 0;
  }

  if (seeds == NULL) {
    p4est_bal_face_con_internal (q, &temp, face, ibalance, &consistent, NULL,
                                 0);
    return consistent;
  }
  else {
    for (i = 0; i < nextra; i++) {
      conextra[i] = p->level;
    }
    p4est_bal_face_con_internal (q, &temp, face, ibalance, &consistent,
                                 conextra, 1);

    sc_array_resize (seeds, 0);
    if (!consistent) {
      for (i = 0; i < nextra; i++) {
        if (conextra[i] == temp.level) {
          sc_array_resize (seeds, seeds->elem_count + 1);
          s = p4est_quadrant_array_index (seeds, seeds->elem_count - 1);

          *s = temp;

          switch (i % 3) {
          case 0:
            switch (face / 2) {
            case 0:
              s->y += -P4EST_QUADRANT_LEN (temp.level - 1);
              break;
            case 1:
#ifdef P4_TO_P8
            case 2:
#endif
              s->x += -P4EST_QUADRANT_LEN (temp.level - 1);
              break;
            default:
              SC_ABORT_NOT_REACHED ();
            }
            break;
          case 1:
            break;
          case 2:
            switch (face / 2) {
            case 0:
              s->y += P4EST_QUADRANT_LEN (temp.level - 1);
              break;
            case 1:
#ifdef P4_TO_P8
            case 2:
#endif
              s->x += P4EST_QUADRANT_LEN (temp.level - 1);
              break;
            default:
              SC_ABORT_NOT_REACHED ();
            }
            break;
          default:
            SC_ABORT_NOT_REACHED ();
          }

#ifdef P4_TO_P8
          switch ((i / 3) % 3) {
          case 0:
            switch (face / 2) {
            case 0:
            case 1:
              s->z += -P4EST_QUADRANT_LEN (temp.level - 1);
              break;
            case 2:
              s->y += -P4EST_QUADRANT_LEN (temp.level - 1);
              break;
            default:
              SC_ABORT_NOT_REACHED ();
            }
            break;
          case 1:
            break;
          case 2:
            switch (face / 2) {
            case 0:
            case 1:
              s->z += P4EST_QUADRANT_LEN (temp.level - 1);
              break;
            case 2:
              s->y += P4EST_QUADRANT_LEN (temp.level - 1);
              break;
            default:
              SC_ABORT_NOT_REACHED ();
            }
            break;
          default:
            SC_ABORT_NOT_REACHED ();
          }
#endif
        }
      }
    }

    return consistent;
  }
}

int
p4est_balance_corner_test (p4est_quadrant_t * restrict q,
                           p4est_quadrant_t * restrict p,
                           int corner, p4est_balance_type_t balance,
                           sc_array_t * seeds)
{
  p4est_quadrant_t    temp = *p;
  p4est_quadrant_t   *s;
  int                 ibalance;
  int                 consistent;

  P4EST_ASSERT (seeds == NULL ||
                seeds->elem_size == sizeof (p4est_quadrant_t));

  if (balance == P4EST_BALANCE_FULL) {
    ibalance = P4EST_DIM - 1;
  }
#ifdef P4_TO_P8
  else if (balance == P8EST_BALANCE_EDGE) {
    ibalance = 1;
  }
#endif
  else {
    ibalance = 0;
  }

  if (seeds == NULL) {
    p4est_bal_corner_con_internal (q, &temp, corner, ibalance, &consistent,
                                   0);
    return consistent;
  }
  else {
    p4est_bal_corner_con_internal (q, &temp, corner, ibalance, &consistent,
                                   1);

    sc_array_resize (seeds, 0);
    if (!consistent) {
      sc_array_resize (seeds, seeds->elem_count + 1);
      s = p4est_quadrant_array_index (seeds, seeds->elem_count - 1);
      *s = temp;
    }

    return consistent;
  }
}

#ifdef P4_TO_P8
int
p8est_balance_edge_test (p4est_quadrant_t * restrict q,
                         p4est_quadrant_t * restrict p,
                         int edge, p4est_balance_type_t balance,
                         sc_array_t * seeds)
{
  p4est_quadrant_t    temp = *p;
  p4est_quadrant_t   *s;
  int                 ibalance;
  int                 consistent;
  int                 nextra = 3;
  int                 conextra[3];
  int                 i;

  P4EST_ASSERT (seeds == NULL ||
                seeds->elem_size == sizeof (p4est_quadrant_t));

  if (balance == P4EST_BALANCE_FULL) {
    ibalance = P4EST_DIM - 1;
  }
#ifdef P4_TO_P8
  else if (balance == P8EST_BALANCE_EDGE) {
    ibalance = 1;
  }
#endif
  else {
    ibalance = 0;
  }

  if (seeds == NULL) {
    p8est_bal_edge_con_internal (q, &temp, edge, ibalance, &consistent, NULL,
                                 0);
    return consistent;
  }
  else {
    for (i = 0; i < nextra; i++) {
      conextra[i] = p->level;
    }
    p8est_bal_edge_con_internal (q, &temp, edge, ibalance, &consistent,
                                 conextra, 1);

    sc_array_resize (seeds, 0);
    if (!consistent) {
      for (i = 0; i < nextra; i++) {
        if (conextra[i] == temp.level) {
          sc_array_resize (seeds, seeds->elem_count + 1);
          s = p4est_quadrant_array_index (seeds, seeds->elem_count - 1);

          *s = temp;

          switch (i % 3) {
          case 0:
            switch (edge / 4) {
            case 0:
              s->x += -P4EST_QUADRANT_LEN (temp.level - 1);
              break;
            case 1:
              s->y += -P4EST_QUADRANT_LEN (temp.level - 1);
              break;
            case 2:
              s->z += -P4EST_QUADRANT_LEN (temp.level - 1);
              break;
            default:
              SC_ABORT_NOT_REACHED ();
            }
            break;
          case 1:
            break;
          case 2:
            switch (edge / 4) {
            case 0:
              s->x += P4EST_QUADRANT_LEN (temp.level - 1);
              break;
            case 1:
              s->y += P4EST_QUADRANT_LEN (temp.level - 1);
              break;
            case 2:
              s->z += P4EST_QUADRANT_LEN (temp.level - 1);
              break;
            default:
              SC_ABORT_NOT_REACHED ();
            }
            break;
          default:
            SC_ABORT_NOT_REACHED ();
          }
        }
      }
    }

    return consistent;
  }
}
#endif
