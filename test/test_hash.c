
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
  p4est_quadrant_t    q1, q2, q3;
  p4est_quadrant_t   *f1, *f2, *f3;
  p4est_hash_t       *ihash;
  p4est_hash_t       *qhash;

  for (k = 0; k < 3; ++k) {
    ihash = p4est_hash_new (21 + k * 70, int_hash_fn, int_equal_fn, NULL);

    inserted = 0;
    for (i = 0; i < 347; ++i) {
      inserted += p4est_hash_insert_unique (ihash, (void *) (i % 91), NULL);
    }
    printf ("Integers inserted %d total %d\n", inserted, ihash->elem_count);
    P4EST_CHECK_ABORT (inserted == ihash->elem_count, "Integer hash");

    p4est_hash_destroy (ihash);
  }

  qhash = p4est_hash_new (17, p4est_quadrant_hash_fn,
                          p4est_quadrant_is_equal, NULL);

  p4est_quadrant_set_morton (&q1, 3, 15);
  p4est_quadrant_set_morton (&q2, 3, 18);
  p4est_quadrant_set_morton (&q3, 3, 18);
  q1.user_data = NULL;
  q2.user_data = (void *) 5;
  q3.user_data = NULL;

  f1 = f2 = f3 = NULL;
  i1 = p4est_hash_insert_unique (qhash, &q1, (void **) &f1);
  i2 = p4est_hash_insert_unique (qhash, &q2, (void **) &f2);
  i3 = p4est_hash_insert_unique (qhash, &q3, (void **) &f3);
  printf ("Quadrants inserted %d %d %d total %d\n", i1, i2, i3,
          qhash->elem_count);

  P4EST_CHECK_ABORT (i1 + i2 + i3 == qhash->elem_count, "Quadrant hash");
  P4EST_CHECK_ABORT (f3 == &q2
                     && f3->user_data == (void *) 5, "Insert return");

  f1 = f2 = f3 = NULL;
  p4est_quadrant_set_morton (&q1, 3, 19);
  i1 = p4est_hash_lookup (qhash, &q1, NULL);
  i2 = p4est_hash_lookup (qhash, &q2, NULL);
  i3 = p4est_hash_lookup (qhash, &q3, (void **) &f3);
  printf ("Quadrants lookup %d %d %d total %d\n", i1, i2, i3,
          qhash->elem_count);
  P4EST_CHECK_ABORT (i1 == 0 && i2 == 1 && i3 == 1, "Quadrant lookup");
  P4EST_CHECK_ABORT (f3 == &q2
                     && f3->user_data == (void *) 5, "Lookup return");

  p4est_hash_destroy (qhash);

  p4est_memory_check ();

  return 0;
}

/* EOF test_hash.c */
