
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

static int
quadrant_hash_fn (const void *v)
{
  const p4est_quadrant_t *q = v;

  return p4est_quadrant_linear_id (q, q->level) % (1LL << 30);
}

static int
quadrant_equal_fn (const void *v1, const void *v2)
{
  return p4est_quadrant_is_equal (v1, v2);
}

int
main (int argc, char **argv)
{
  int                 i, k, inserted;
  int                 i1, i2, i3;
  p4est_quadrant_t    q1, q2, q3;
  p4est_hash_t       *ihash;
  p4est_hash_t       *qhash;

  for (k = 0; k < 3; ++k) {
    ihash = p4est_hash_new (21 + k * 70, int_hash_fn, int_equal_fn, NULL);

    inserted = 0;
    for (i = 0; i < 347; ++i) {
      inserted += p4est_hash_insert_unique (ihash, (void *) (i % 91));
    }
    printf ("Integers inserted %d total %d\n", inserted, ihash->elem_count);
    P4EST_CHECK_ABORT (inserted == ihash->elem_count, "Integer hash");

    p4est_hash_destroy (ihash);
  }

  qhash = p4est_hash_new (17, quadrant_hash_fn, quadrant_equal_fn, NULL);

  p4est_quadrant_set_morton (&q1, 3, 15);
  p4est_quadrant_set_morton (&q2, 3, 18);
  p4est_quadrant_set_morton (&q3, 3, 18);

  i1 = p4est_hash_insert_unique (qhash, &q1);
  i2 = p4est_hash_insert_unique (qhash, &q2);
  i3 = p4est_hash_insert_unique (qhash, &q3);
  printf ("Quadrants inserted %d %d %d total %d\n", i1, i2, i3,
          qhash->elem_count);
  P4EST_CHECK_ABORT (i1 + i2 + i3 == qhash->elem_count, "Quadrant hash");

  p4est_hash_destroy (qhash);

  p4est_memory_check ();

  return 0;
}

/* EOF test_hash.c */
