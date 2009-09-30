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

#include <p4est_bits.h>

static unsigned
int_hash_fn (const void *v, const void *u)
{
  return (unsigned) (unsigned long) v;
}

static int
int_equal_fn (const void *v1, const void *v2, const void *u)
{
  return (long) v1 == (long) v2;
}

static void
check_hash_array (void)
{
  sc_hash_array_t    *ha;

  ha = sc_hash_array_new (sizeof (int), NULL, NULL, NULL);

  sc_hash_array_destroy (ha);
}

int
main (int argc, char **argv)
{
  int                 i, k, inserted;
  int                 i1, i2, i3;
  void              **vv1, **vv2, **vv3;
  void               *v1, *v2;
  p4est_quadrant_t    q1, q2, q3;
  p4est_quadrant_t   *f1, *f2, *f3;
  sc_hash_t          *ihash;
  sc_hash_t          *qhash;

  for (k = 0; k < 3; ++k) {
    ihash = sc_hash_new (int_hash_fn, int_equal_fn, NULL, NULL);

    inserted = 0;
    for (i = 0; i < 347; ++i) {
      inserted +=
        sc_hash_insert_unique (ihash, (void *) ((long) i % 91), NULL);
    }
    printf ("Integers inserted %d total %llu\n",
            inserted, (unsigned long long) ihash->elem_count);
    SC_CHECK_ABORT (inserted == (int) ihash->elem_count, "Integer hash");

    sc_hash_destroy (ihash);
  }

  qhash = sc_hash_new (p4est_quadrant_hash_fn, p4est_quadrant_equal_fn,
                       NULL, NULL);

  p4est_quadrant_set_morton (&q1, 3, 15);
  p4est_quadrant_set_morton (&q2, 3, 18);
  p4est_quadrant_set_morton (&q3, 3, 18);
  q1.p.piggy1.owner_rank = 0;
  q2.p.piggy1.owner_rank = 5;
  q3.p.piggy1.owner_rank = 0;

  f1 = f2 = f3 = NULL;
  i1 = sc_hash_insert_unique (qhash, &q1, &vv1);
  f1 = *vv1;
  i2 = sc_hash_insert_unique (qhash, &q2, &vv2);
  f2 = *vv2;
  i3 = sc_hash_insert_unique (qhash, &q3, &vv3);
  f3 = *vv3;
  printf ("Quadrants inserted %d %d %d total %lu\n",
          i1, i2, i3, (unsigned long) qhash->elem_count);

  SC_CHECK_ABORT (i1 + i2 + i3 == (int) qhash->elem_count, "Quadrant hash");
  SC_CHECK_ABORT (f3 == &q2 && f3->p.piggy1.owner_rank == 5, "Insert return");

  f1 = f2 = f3 = NULL;
  p4est_quadrant_set_morton (&q1, 3, 19);
  i1 = sc_hash_lookup (qhash, &q1, NULL);
  i2 = sc_hash_lookup (qhash, &q2, NULL);
  i3 = sc_hash_lookup (qhash, &q3, &vv3);
  f3 = *vv3;
  printf ("Quadrants lookup %d %d %d total %lu\n",
          i1, i2, i3, (unsigned long) qhash->elem_count);
  SC_CHECK_ABORT (i1 == 0 && i2 == 1 && i3 == 1, "Quadrant lookup");
  SC_CHECK_ABORT (f3 == &q2 && f3->p.piggy1.owner_rank == 5, "Lookup return");

  f1 = f2 = f3 = NULL;
  i1 = sc_hash_remove (qhash, &q1, &v1);
  f1 = v1;
  i2 = sc_hash_remove (qhash, &q2, &v2);
  f2 = v2;
  i3 = sc_hash_remove (qhash, &q3, NULL);
  SC_CHECK_ABORT (i1 == 0 && i2 == 1 && i3 == 0, "Quadrant remove");
  SC_CHECK_ABORT (f2 == &q2 && f2->p.piggy1.owner_rank == 5, "Remove return");
  f2 = f1;

  sc_hash_destroy (qhash);

  check_hash_array ();

  sc_finalize ();

  return 0;
}
