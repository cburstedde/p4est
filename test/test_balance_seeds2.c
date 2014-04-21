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

#ifndef P4_TO_P8
#include <p4est_balance.h>
#include <p4est_bits.h>
#else
#include <p8est_balance.h>
#include <p8est_bits.h>
#endif

static int
is_farther (p4est_quadrant_t * orig, p4est_quadrant_t * targ,
            p4est_quadrant_t * newq)
{
  p4est_qcoord_t      ox1, ox2, nx1, nx2;
  p4est_qcoord_t      oy1, oy2, ny1, ny2;
#ifdef P4_TO_P8
  p4est_qcoord_t      oz1, oz2, nz1, nz2;
#endif
  p4est_qcoord_t      nl = P4EST_QUADRANT_LEN (newq->level);
  p4est_qcoord_t      tl = P4EST_QUADRANT_LEN (targ->level);

  P4EST_ASSERT (newq->level == orig->level);

  ox1 = targ->x - orig->x;
  ox2 = (orig->x + nl) - (targ->x + tl);
  nx1 = targ->x - newq->x;
  nx2 = (newq->x + nl) - (targ->x + tl);
  if (ox1 > 0 && nx1 > ox1) {
    return 1;
  }
  if (ox2 > 0 && nx2 > ox2) {
    return 1;
  }

  oy1 = targ->y - orig->y;
  oy2 = (orig->y + nl) - (targ->y + tl);
  ny1 = targ->y - newq->y;
  ny2 = (newq->y + nl) - (targ->y + tl);
  if (oy1 > 0 && ny1 > oy1) {
    return 1;
  }
  if (oy2 > 0 && ny2 > oy2) {
    return 1;
  }

#ifdef P4_TO_P8
  oz1 = targ->z - orig->z;
  oz2 = (orig->z + nl) - (targ->z + tl);
  nz1 = targ->z - newq->z;
  nz2 = (newq->z + nl) - (targ->z + tl);
  if (oz1 > 0 && nz1 > oz1) {
    return 1;
  }
  if (oz2 > 0 && nz2 > oz2) {
    return 1;
  }
#endif

  return 0;

}

int
check_balance_seeds (p4est_quadrant_t * q, p4est_quadrant_t * p,
                     p4est_connect_type_t b, sc_array_t * seeds)
{
  int                 ib;
  int                 level = q->level;
  p4est_quadrant_t   *s, *t;
  sc_array_t         *thislevel = sc_array_new (sizeof (p4est_quadrant_t));
  sc_array_t         *nextlevel = sc_array_new (sizeof (p4est_quadrant_t));
  sc_array_t         *temparray;
  p4est_quadrant_t    temp1, temp2;
  int                 f, c;
#ifdef P4_TO_P8
  int                 e;
#endif
  int                 stop = 0;

  sc_array_resize (seeds, 0);

  s = (p4est_quadrant_t *) sc_array_push (thislevel);
  p4est_quadrant_sibling (q, s, 0);

#ifndef P4_TO_P8
  if (b == P4EST_CONNECT_FACE) {
    ib = 0;
  }
  else {
    ib = 1;
  }
#else
  if (b == P8EST_CONNECT_FACE) {
    ib = 0;
  }
  else if (b == P8EST_CONNECT_EDGE) {
    ib = 1;
  }
  else {
    ib = 2;
  }
#endif

  while (level > p->level + 1) {
    size_t              nlast = thislevel->elem_count;
    size_t              zz;

    stop = 0;

    for (zz = 0; zz < nlast; zz++) {
      s = p4est_quadrant_array_index (thislevel, zz);
      P4EST_ASSERT (p4est_quadrant_child_id (s) == 0);
      p4est_quadrant_parent (s, &temp1);
      for (f = 0; f < P4EST_FACES; f++) {
        p4est_quadrant_face_neighbor (&temp1, f, &temp2);
        if (is_farther (&temp1, p, &temp2)) {
          continue;
        }
        if (p4est_quadrant_is_ancestor (p, &temp2)) {
          stop = 1;
          sc_array_resize (seeds, seeds->elem_count + 1);
          t = p4est_quadrant_array_index (seeds, seeds->elem_count - 1);
          p4est_quadrant_sibling (&temp2, t, 0);
        }
        else if (p4est_quadrant_is_inside_root (&temp2)) {
          t = (p4est_quadrant_t *) sc_array_push (nextlevel);
          p4est_quadrant_sibling (&temp2, t, 0);
        }
      }

      if (ib == 0) {
        continue;
      }

#ifdef P4_TO_P8
      for (e = 0; e < P8EST_EDGES; e++) {
        p8est_quadrant_edge_neighbor (&temp1, e, &temp2);
        if (is_farther (&temp1, p, &temp2)) {
          continue;
        }
        if (p4est_quadrant_is_ancestor (p, &temp2)) {
          stop = 1;
          sc_array_resize (seeds, seeds->elem_count + 1);
          t = p4est_quadrant_array_index (seeds, seeds->elem_count - 1);
          p4est_quadrant_sibling (&temp2, t, 0);
        }
        else if (p4est_quadrant_is_inside_root (&temp2)) {
          t = (p4est_quadrant_t *) sc_array_push (nextlevel);
          p4est_quadrant_sibling (&temp2, t, 0);
        }
      }

      if (ib == 1) {
        continue;
      }
#endif

      for (c = 0; c < P4EST_CHILDREN; c++) {
        p4est_quadrant_corner_neighbor (&temp1, c, &temp2);
        if (is_farther (&temp1, p, &temp2)) {
          continue;
        }
        if (p4est_quadrant_is_ancestor (p, &temp2)) {
          stop = 1;
          sc_array_resize (seeds, seeds->elem_count + 1);
          t = p4est_quadrant_array_index (seeds, seeds->elem_count - 1);
          p4est_quadrant_sibling (&temp2, t, 0);
        }
        else if (p4est_quadrant_is_inside_root (&temp2)) {
          t = (p4est_quadrant_t *) sc_array_push (nextlevel);
          p4est_quadrant_sibling (&temp2, t, 0);
        }
      }
    }

    if (stop) {
      sc_array_sort (seeds, p4est_quadrant_compare);
      sc_array_uniq (seeds, p4est_quadrant_compare);

#ifdef P4_TO_P8
      if (!ib && seeds->elem_count == 1) {
        sc_array_sort (nextlevel, p4est_quadrant_compare);
        sc_array_uniq (nextlevel, p4est_quadrant_compare);
        temparray = thislevel;
        thislevel = nextlevel;
        nextlevel = temparray;
        sc_array_reset (nextlevel);
        level--;

        nlast = thislevel->elem_count;
        for (zz = 0; zz < nlast; zz++) {
          s = p4est_quadrant_array_index (thislevel, zz);
          P4EST_ASSERT (p4est_quadrant_child_id (s) == 0);
          p4est_quadrant_parent (s, &temp1);
          for (f = 0; f < P4EST_FACES; f++) {
            p4est_quadrant_face_neighbor (&temp1, f, &temp2);
            if (p4est_quadrant_is_ancestor (p, &temp2)) {
              int                 f2;
              p4est_quadrant_t    a;
              p4est_quadrant_t    u;

              t = p4est_quadrant_array_index (seeds, 0);

              p8est_quadrant_parent (t, &a);

              for (f2 = 0; f2 < P8EST_FACES; f2++) {
                if (f2 / 2 == f / 2) {
                  continue;
                }
                p8est_quadrant_face_neighbor (&a, f2, &u);

                if (p8est_quadrant_is_equal (&temp2, &u) ||
                    p8est_quadrant_is_sibling (&temp2, &u)) {
                  break;
                }
              }

              if (f2 == P8EST_FACES) {
                sc_array_resize (seeds, seeds->elem_count + 1);
                t = p4est_quadrant_array_index (seeds, seeds->elem_count - 1);
                p4est_quadrant_sibling (&temp2, t, 0);
              }
            }
          }
        }
      }
#endif
      sc_array_sort (seeds, p4est_quadrant_compare);
      sc_array_uniq (seeds, p4est_quadrant_compare);

      break;
    }
    sc_array_sort (nextlevel, p4est_quadrant_compare);
    sc_array_uniq (nextlevel, p4est_quadrant_compare);
    temparray = thislevel;
    thislevel = nextlevel;
    nextlevel = temparray;
    sc_array_reset (nextlevel);
    level--;
  }

  sc_array_destroy (thislevel);
  sc_array_destroy (nextlevel);

  return stop;
}

void
standard_seeds (sc_array_t * seeds)
{
  size_t              count = seeds->elem_count;
  size_t              zz;
  p4est_quadrant_t   *q, temp;

  for (zz = 0; zz < count; zz++) {
    q = p4est_quadrant_array_index (seeds, zz);
    p4est_quadrant_sibling (q, &temp, 0);
    *q = temp;
  }
  sc_array_sort (seeds, p4est_quadrant_compare);
  sc_array_uniq (seeds, p4est_quadrant_compare);
}

void
compare_seeds (sc_array_t * seeds, sc_array_t * seeds_check)
{
  size_t              count = seeds->elem_count;
  size_t              zz;
  p4est_quadrant_t   *s, *t;

  SC_CHECK_ABORT (seeds_check->elem_count == count, "seed count");

  for (zz = 0; zz < count; zz++) {
    s = p4est_quadrant_array_index (seeds, zz);
    t = p4est_quadrant_array_index (seeds_check, zz);
    SC_CHECK_ABORT (p4est_quadrant_is_equal (s, t), "seed equality");
  }
}

int
main (int argc, char **argv)
{
  p4est_quadrant_t    root;
  p4est_quadrant_t    p;
  p4est_quadrant_t    q;
  p4est_quadrant_t    desc;
  int                 face, corner;
#ifndef P4_TO_P8
  int                 maxlevel = 9;
#else
  int                 edge;
  int                 maxlevel = 6;
#endif
  int                 mpiret, mpisize, mpirank;
  sc_MPI_Comm         mpicomm;
  uint64_t            i, ifirst, ilast;
  int                 level;
  sc_array_t         *seeds, *seeds_check;
  int                 testval;
  int                 checkval;
  int                 j, nrand = 1000;

  /* initialize MPI */
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  mpicomm = sc_MPI_COMM_WORLD;
  mpiret = sc_MPI_Comm_size (mpicomm, &mpisize);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (mpicomm, &mpirank);
  SC_CHECK_MPI (mpiret);

  srandom (9212007);
  sc_init (mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  seeds = sc_array_new (sizeof (p4est_quadrant_t));
  seeds_check = sc_array_new (sizeof (p4est_quadrant_t));

  memset (&root, 0, sizeof (p4est_quadrant_t));
  root.level = 2;
  root.x = P4EST_QUADRANT_LEN (2);
  root.y = P4EST_QUADRANT_LEN (2);
#ifdef P4_TO_P8
  root.z = P4EST_QUADRANT_LEN (2);
#endif
  P4EST_QUADRANT_INIT (&p);
  P4EST_QUADRANT_INIT (&q);

#if 1
  for (face = 0; face < P4EST_FACES; face++) {
    p4est_quadrant_face_neighbor (&root, face ^ 1, &p);
    P4EST_GLOBAL_VERBOSEF ("Testing face %d\n", face);
    for (level = 4; level <= maxlevel; level++) {
      P4EST_GLOBAL_VERBOSEF (" level %d\n", level);
      p4est_quadrant_first_descendant (&root, &desc, level);
      ifirst = p4est_quadrant_linear_id (&desc, level);
      p4est_quadrant_last_descendant (&root, &desc, level);
      ilast = p4est_quadrant_linear_id (&desc, level);
      for (i = ifirst; i <= ilast; i += P4EST_CHILDREN) {
        p4est_quadrant_set_morton (&q, level, i);
#ifndef P4_TO_P8
        testval = p4est_balance_seeds_face (&q, &p, face, P4EST_CONNECT_FACE,
                                            seeds);
        standard_seeds (seeds);
        checkval = check_balance_seeds (&q, &p, P4EST_CONNECT_FACE,
                                        seeds_check);
        SC_CHECK_ABORT (testval == checkval,
                        "p4est_balance_seeds_face error");
        compare_seeds (seeds, seeds_check);
#else
        testval = p4est_balance_seeds_face (&q, &p, face, P8EST_CONNECT_FACE,
                                            seeds);
        standard_seeds (seeds);
        checkval = check_balance_seeds (&q, &p, P8EST_CONNECT_FACE,
                                        seeds_check);
        SC_CHECK_ABORT (testval == checkval,
                        "p8est_balance_seeds_face error");
        compare_seeds (seeds, seeds_check);
        testval = p4est_balance_seeds_face (&q, &p, face, P8EST_CONNECT_EDGE,
                                            seeds);
        standard_seeds (seeds);
        checkval = check_balance_seeds (&q, &p, P8EST_CONNECT_EDGE,
                                        seeds_check);
        SC_CHECK_ABORT (testval == checkval,
                        "p8est_balance_seeds_face error");
        compare_seeds (seeds, seeds_check);
#endif
        testval = p4est_balance_seeds_face (&q, &p, face, P4EST_CONNECT_FULL,
                                            seeds);
        standard_seeds (seeds);
        checkval = check_balance_seeds (&q, &p, P4EST_CONNECT_FULL,
                                        seeds_check);
        SC_CHECK_ABORT (testval == checkval,
                        "p4est_balance_seeds_face error");
        compare_seeds (seeds, seeds_check);
      }
    }
    if (!face) {
      P4EST_GLOBAL_VERBOSE (" random levels\n");
      for (j = 0; j < (int) nrand; j++) {
        level = ((random ()) % (P4EST_QMAXLEVEL - maxlevel)) + maxlevel + 1;
        p4est_quadrant_first_descendant (&root, &desc, level);
        ifirst = p4est_quadrant_linear_id (&desc, level);
        p4est_quadrant_last_descendant (&root, &desc, level);
        ilast = p4est_quadrant_linear_id (&desc, level);
        i = ((random ()) % (ilast + 1 - ifirst)) + ifirst;
        p4est_quadrant_set_morton (&q, level, i);
#ifndef P4_TO_P8
        testval = p4est_balance_seeds_face (&q, &p, face, P4EST_CONNECT_FACE,
                                            seeds);
        standard_seeds (seeds);
        checkval = check_balance_seeds (&q, &p, P4EST_CONNECT_FACE,
                                        seeds_check);
        SC_CHECK_ABORT (testval == checkval,
                        "p4est_balance_seeds_face error");
        compare_seeds (seeds, seeds_check);
#else
        testval = p4est_balance_seeds_face (&q, &p, face, P8EST_CONNECT_FACE,
                                            seeds);
        standard_seeds (seeds);
        checkval = check_balance_seeds (&q, &p, P8EST_CONNECT_FACE,
                                        seeds_check);
        SC_CHECK_ABORT (testval == checkval,
                        "p8est_balance_seeds_face error");
        compare_seeds (seeds, seeds_check);
        testval = p4est_balance_seeds_face (&q, &p, face, P8EST_CONNECT_EDGE,
                                            seeds);
        standard_seeds (seeds);
        checkval = check_balance_seeds (&q, &p, P8EST_CONNECT_EDGE,
                                        seeds_check);
        SC_CHECK_ABORT (testval == checkval,
                        "p8est_balance_seeds_face error");
        compare_seeds (seeds, seeds_check);
#endif
        testval = p4est_balance_seeds_face (&q, &p, face, P4EST_CONNECT_FULL,
                                            seeds);
        standard_seeds (seeds);
        checkval = check_balance_seeds (&q, &p, P4EST_CONNECT_FULL,
                                        seeds_check);
        SC_CHECK_ABORT (testval == checkval,
                        "p4est_balance_seeds_face error");
        compare_seeds (seeds, seeds_check);
      }
    }
  }

#ifdef P4_TO_P8
  for (edge = 0; edge < P8EST_EDGES; edge++) {
    p8est_quadrant_edge_neighbor (&root, edge ^ 3, &p);
    P4EST_GLOBAL_VERBOSEF ("Testing edge %d\n", edge);
    for (level = 4; level <= maxlevel; level++) {
      P4EST_GLOBAL_VERBOSEF (" level %d\n", level);
      p4est_quadrant_first_descendant (&root, &desc, level);
      ifirst = p4est_quadrant_linear_id (&desc, level);
      p4est_quadrant_last_descendant (&root, &desc, level);
      ilast = p4est_quadrant_linear_id (&desc, level);
      for (i = ifirst; i <= ilast; i += P4EST_CHILDREN) {
        p4est_quadrant_set_morton (&q, level, i);
        testval = p8est_balance_seeds_edge (&q, &p, edge, P8EST_CONNECT_FACE,
                                            seeds);
        standard_seeds (seeds);
        checkval = check_balance_seeds (&q, &p, P8EST_CONNECT_FACE,
                                        seeds_check);
        SC_CHECK_ABORT (testval == checkval,
                        "p8est_balance_seeds_edge error");
        compare_seeds (seeds, seeds_check);
        testval = p8est_balance_seeds_edge (&q, &p, edge, P8EST_CONNECT_EDGE,
                                            seeds);
        standard_seeds (seeds);
        checkval = check_balance_seeds (&q, &p, P8EST_CONNECT_EDGE,
                                        seeds_check);
        SC_CHECK_ABORT (testval == checkval,
                        "p8est_balance_seeds_edge error");
        compare_seeds (seeds, seeds_check);
        testval = p8est_balance_seeds_edge (&q, &p, edge, P8EST_CONNECT_FULL,
                                            seeds);
        standard_seeds (seeds);
        checkval = check_balance_seeds (&q, &p, P8EST_CONNECT_FULL,
                                        seeds_check);
        SC_CHECK_ABORT (testval == checkval,
                        "p8est_balance_seeds_edge error");
        compare_seeds (seeds, seeds_check);
      }
    }
    if (!edge) {
      P4EST_GLOBAL_VERBOSE (" random levels\n");
      for (j = 0; j < (int) nrand; j++) {
        level = ((random ()) % (P4EST_QMAXLEVEL - maxlevel)) + maxlevel + 1;
        p4est_quadrant_first_descendant (&root, &desc, level);
        ifirst = p4est_quadrant_linear_id (&desc, level);
        p4est_quadrant_last_descendant (&root, &desc, level);
        ilast = p4est_quadrant_linear_id (&desc, level);
        i = ((random ()) % (ilast + 1 - ifirst)) + ifirst;
        p4est_quadrant_set_morton (&q, level, i);
        testval = p8est_balance_seeds_edge (&q, &p, edge, P8EST_CONNECT_FACE,
                                            seeds);
        standard_seeds (seeds);
        checkval = check_balance_seeds (&q, &p, P8EST_CONNECT_FACE,
                                        seeds_check);
        SC_CHECK_ABORT (testval == checkval,
                        "p8est_balance_seeds_edge error");
        compare_seeds (seeds, seeds_check);
        testval = p8est_balance_seeds_edge (&q, &p, edge, P8EST_CONNECT_EDGE,
                                            seeds);
        standard_seeds (seeds);
        checkval = check_balance_seeds (&q, &p, P8EST_CONNECT_EDGE,
                                        seeds_check);
        SC_CHECK_ABORT (testval == checkval,
                        "p8est_balance_seeds_edge error");
        compare_seeds (seeds, seeds_check);
        testval = p8est_balance_seeds_edge (&q, &p, edge, P8EST_CONNECT_FULL,
                                            seeds);
        standard_seeds (seeds);
        checkval = check_balance_seeds (&q, &p, P8EST_CONNECT_FULL,
                                        seeds_check);
        SC_CHECK_ABORT (testval == checkval,
                        "p8est_balance_seeds_edge error");
        compare_seeds (seeds, seeds_check);
      }
    }
  }
#endif
#endif

  for (corner = 0; corner < P4EST_FACES; corner++) {
    p4est_quadrant_corner_neighbor (&root, corner ^ (P4EST_CHILDREN - 1), &p);
    P4EST_GLOBAL_VERBOSEF ("Testing corner %d\n", corner);
    for (level = 4; level <= maxlevel; level++) {
      P4EST_GLOBAL_VERBOSEF (" level %d\n", level);
      p4est_quadrant_first_descendant (&root, &desc, level);
      ifirst = p4est_quadrant_linear_id (&desc, level);
      p4est_quadrant_last_descendant (&root, &desc, level);
      ilast = p4est_quadrant_linear_id (&desc, level);
      for (i = ifirst; i <= ilast; i += P4EST_CHILDREN) {
        p4est_quadrant_set_morton (&q, level, i);
#ifndef P4_TO_P8
        testval =
          p4est_balance_seeds_corner (&q, &p, corner, P4EST_CONNECT_FACE,
                                      seeds);
        standard_seeds (seeds);
        checkval = check_balance_seeds (&q, &p, P4EST_CONNECT_FACE,
                                        seeds_check);
        SC_CHECK_ABORT (testval == checkval,
                        "p4est_balance_seeds_corner error");
        compare_seeds (seeds, seeds_check);
#else
        testval = p4est_balance_seeds_corner (&q, &p, corner,
                                              P8EST_CONNECT_FACE, seeds);
        standard_seeds (seeds);
        checkval = check_balance_seeds (&q, &p, P8EST_CONNECT_FACE,
                                        seeds_check);
        SC_CHECK_ABORT (testval == checkval,
                        "p8est_balance_seeds_corner error");
        compare_seeds (seeds, seeds_check);
        testval =
          p4est_balance_seeds_corner (&q, &p, corner, P8EST_CONNECT_EDGE,
                                      seeds);
        standard_seeds (seeds);
        checkval = check_balance_seeds (&q, &p, P8EST_CONNECT_EDGE,
                                        seeds_check);
        SC_CHECK_ABORT (testval == checkval,
                        "p8est_balance_seeds_corner error");
        compare_seeds (seeds, seeds_check);
#endif
        testval =
          p4est_balance_seeds_corner (&q, &p, corner, P4EST_CONNECT_FULL,
                                      seeds);
        standard_seeds (seeds);
        checkval = check_balance_seeds (&q, &p, P4EST_CONNECT_FULL,
                                        seeds_check);
        SC_CHECK_ABORT (testval == checkval,
                        "p4est_balance_seeds_corner error");
        compare_seeds (seeds, seeds_check);
      }
    }
    if (!corner) {
      P4EST_GLOBAL_VERBOSE (" random levels\n");
      for (j = 0; j < (int) nrand; j++) {
        level = ((random ()) % (P4EST_QMAXLEVEL - maxlevel)) + maxlevel + 1;
        p4est_quadrant_first_descendant (&root, &desc, level);
        ifirst = p4est_quadrant_linear_id (&desc, level);
        p4est_quadrant_last_descendant (&root, &desc, level);
        ilast = p4est_quadrant_linear_id (&desc, level);
        i = ((random ()) % (ilast + 1 - ifirst)) + ifirst;
        p4est_quadrant_set_morton (&q, level, i);
#ifndef P4_TO_P8
        testval =
          p4est_balance_seeds_corner (&q, &p, corner, P4EST_CONNECT_FACE,
                                      seeds);
        standard_seeds (seeds);
        checkval = check_balance_seeds (&q, &p, P4EST_CONNECT_FACE,
                                        seeds_check);
        SC_CHECK_ABORT (testval == checkval,
                        "p4est_balance_seeds_corner error");
        compare_seeds (seeds, seeds_check);
#else
        testval = p4est_balance_seeds_corner (&q, &p, corner,
                                              P8EST_CONNECT_FACE, seeds);
        standard_seeds (seeds);
        checkval = check_balance_seeds (&q, &p, P8EST_CONNECT_FACE,
                                        seeds_check);
        SC_CHECK_ABORT (testval == checkval,
                        "p8est_balance_seeds_corner error");
        compare_seeds (seeds, seeds_check);
        testval =
          p4est_balance_seeds_corner (&q, &p, corner, P8EST_CONNECT_EDGE,
                                      seeds);
        standard_seeds (seeds);
        checkval = check_balance_seeds (&q, &p, P8EST_CONNECT_EDGE,
                                        seeds_check);
        SC_CHECK_ABORT (testval == checkval,
                        "p8est_balance_seeds_corner error");
        compare_seeds (seeds, seeds_check);
#endif
        testval =
          p4est_balance_seeds_corner (&q, &p, corner, P4EST_CONNECT_FULL,
                                      seeds);
        standard_seeds (seeds);
        checkval = check_balance_seeds (&q, &p, P4EST_CONNECT_FULL,
                                        seeds_check);
        SC_CHECK_ABORT (testval == checkval,
                        "p4est_balance_seeds_corner error");
        compare_seeds (seeds, seeds_check);
      }
    }
  }

  sc_array_destroy (seeds);
  sc_array_destroy (seeds_check);

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
