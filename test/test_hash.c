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

#include <p4est_algorithms.h>
#include <p4est_base.h>

static int
int_hash_fn (const void *v)
{
  return (int) v;
}

static int
int_equal_fn (const void *v1, const void *v2)
{
  return (int) v1 == (int) v2;
}

int
main (int argc, char **argv)
{
  int                 i, k, inserted;
  int                 i1, i2, i3;
  void               *v1, *v2, *v3;
  p4est_quadrant_t    q1, q2, q3;
  p4est_quadrant_t   *f1, *f2, *f3;
  p4est_hash_t       *ihash;
  p4est_hash_t       *qhash;

  p4est_init (stdout, 0, NULL, NULL);

  for (k = 0; k < 3; ++k) {
    ihash = p4est_hash_new (int_hash_fn, int_equal_fn, NULL);

    inserted = 0;
    for (i = 0; i < 347; ++i) {
      inserted += p4est_hash_insert_unique (ihash, (void *) (i % 91), NULL);
    }
    printf ("Integers inserted %d total %d\n", inserted, ihash->elem_count);
    P4EST_CHECK_ABORT (inserted == ihash->elem_count, "Integer hash");

    p4est_hash_destroy (ihash);
  }

  qhash = p4est_hash_new (p4est_quadrant_hash, p4est_quadrant_is_equal, NULL);

  p4est_quadrant_set_morton (&q1, 3, 15);
  p4est_quadrant_set_morton (&q2, 3, 18);
  p4est_quadrant_set_morton (&q3, 3, 18);
  q1.user_data = NULL;
  q2.user_data = (void *) 5;
  q3.user_data = NULL;

  f1 = f2 = f3 = NULL;
  i1 = p4est_hash_insert_unique (qhash, &q1, &v1);
  f1 = v1;
  i2 = p4est_hash_insert_unique (qhash, &q2, &v2);
  f2 = v2;
  i3 = p4est_hash_insert_unique (qhash, &q3, &v3);
  f3 = v3;
  printf ("Quadrants inserted %d %d %d total %d\n", i1, i2, i3,
          qhash->elem_count);

  P4EST_CHECK_ABORT (i1 + i2 + i3 == qhash->elem_count, "Quadrant hash");
  P4EST_CHECK_ABORT (f3 == &q2
                     && f3->user_data == (void *) 5, "Insert return");

  f1 = f2 = f3 = NULL;
  p4est_quadrant_set_morton (&q1, 3, 19);
  i1 = p4est_hash_lookup (qhash, &q1, NULL);
  i2 = p4est_hash_lookup (qhash, &q2, NULL);
  i3 = p4est_hash_lookup (qhash, &q3, &v3);
  f3 = v3;
  printf ("Quadrants lookup %d %d %d total %d\n", i1, i2, i3,
          qhash->elem_count);
  P4EST_CHECK_ABORT (i1 == 0 && i2 == 1 && i3 == 1, "Quadrant lookup");
  P4EST_CHECK_ABORT (f3 == &q2
                     && f3->user_data == (void *) 5, "Lookup return");

  f1 = f2 = f3 = NULL;
  i1 = p4est_hash_remove (qhash, &q1, &v1);
  f1 = v1;
  i2 = p4est_hash_remove (qhash, &q2, &v2);
  f2 = v2;
  i3 = p4est_hash_remove (qhash, &q3, NULL);
  P4EST_CHECK_ABORT (i1 == 0 && i2 == 1 && i3 == 0, "Quadrant remove");
  P4EST_CHECK_ABORT (f2 == &q2
                     && f2->user_data == (void *) 5, "Remove return");
  f2 = f1;

  p4est_hash_destroy (qhash);

  p4est_memory_check ();

  return 0;
}

/* EOF test_hash.c */
