
#include <p4est_algorithms.h>
#include <p4est_base.h>

static const int8_t log_lookup_table[256] =
  { -1, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3,
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
  5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
  6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
  6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
  6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
  6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
  7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
  7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
  7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
  7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
  7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
  7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
  7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
  7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
};

#define P4EST_LOG2_8(x) (log_lookup_table[(x)])
#define P4EST_LOG2_16(x) (((x) > 0xff) ? \
                          P4EST_LOG2_8 ((x) >> 8) : P4EST_LOG2_8 (x))
#define P4EST_LOG2_32(x) (((x) > 0xffff) ? \
                          P4EST_LOG2_16 ((x) >> 16) : P4EST_LOG2_16 (x))

int
p4est_quadrant_compare (const void *v1, const void *v2)
{
  const p4est_quadrant_t *q1 = v1;
  const p4est_quadrant_t *q2 = v2;

  int32_t             exclorx, exclory;

  P4EST_ASSERT (p4est_quadrant_is_valid (q1));
  P4EST_ASSERT (p4est_quadrant_is_valid (q2));

  exclorx = q1->x ^ q2->x;
  exclory = q1->y ^ q2->y;

  if (exclory == 0 && exclorx == 0) {
    return q1->level - q2->level;
  }
  else if (exclory >= exclorx) {
    return q1->y - q2->y;
  }
  else {
    return q1->x - q2->x;
  }
}

int
p4est_quadrant_is_valid (const p4est_quadrant_t * q)
{
  return
    (q->level >= 0 && q->level <= P4EST_MAXLEVEL) &&
    (q->x >= 0 && q->x < (1 << P4EST_MAXLEVEL)) &&
    (q->y >= 0 && q->y < (1 << P4EST_MAXLEVEL)) &&
    ((q->x & ((1 << (P4EST_MAXLEVEL - q->level)) - 1)) == 0) &&
    ((q->y & ((1 << (P4EST_MAXLEVEL - q->level)) - 1)) == 0);
}

int
p4est_quadrant_is_equal (const p4est_quadrant_t * q1,
                         const p4est_quadrant_t * q2)
{
  P4EST_ASSERT (p4est_quadrant_is_valid (q1));
  P4EST_ASSERT (p4est_quadrant_is_valid (q2));

  return (q1->level == q2->level && q1->x == q2->x && q1->y == q2->y);
}

int
p4est_quadrant_is_ancestor (const p4est_quadrant_t * q,
                            const p4est_quadrant_t * r)
{
  int32_t             exclorx;
  int32_t             exclory;

  P4EST_ASSERT (p4est_quadrant_is_valid (q));
  P4EST_ASSERT (p4est_quadrant_is_valid (r));

  if (q->level > r->level) {
    return 0;
  }

  exclorx = (q->x ^ r->x) >> (P4EST_MAXLEVEL - q->level);
  exclory = (q->y ^ r->y) >> (P4EST_MAXLEVEL - q->level);

  return (exclorx == 0 && exclory == 0);
}

int
p4est_quadrant_is_ancestor_D (const p4est_quadrant_t * q,
                              const p4est_quadrant_t * r)
{
  p4est_quadrant_t    s;

  p4est_nearest_common_ancestor_D (q, r, &s);

  return p4est_quadrant_is_equal (q, &s);
}

void
p4est_quadrant_parent (const p4est_quadrant_t * q, p4est_quadrant_t * r)
{
  P4EST_ASSERT (p4est_quadrant_is_valid (q));
  P4EST_ASSERT (q->level > 0);

  r->x = q->x & ~(1 << (P4EST_MAXLEVEL - q->level));
  r->y = q->y & ~(1 << (P4EST_MAXLEVEL - q->level));
  r->level = (int8_t) (q->level - 1);

  P4EST_ASSERT (p4est_quadrant_is_valid (r));
}

void
p4est_quadrant_children (const p4est_quadrant_t * q,
                         p4est_quadrant_t * c0, p4est_quadrant_t * c1,
                         p4est_quadrant_t * c2, p4est_quadrant_t * c3)
{
  P4EST_ASSERT (p4est_quadrant_is_valid (q));
  P4EST_ASSERT (q->level < P4EST_MAXLEVEL);

  c0->x = q->x;
  c0->y = q->y;
  c0->level = (int8_t) (q->level + 1);

  c1->x = c0->x | (1 << (P4EST_MAXLEVEL - c0->level));
  c1->y = c0->y;
  c1->level = c0->level;

  c2->x = c0->x;
  c2->y = c0->y | (1 << (P4EST_MAXLEVEL - c0->level));
  c2->level = c0->level;

  c3->x = c1->x;
  c3->y = c2->y;
  c3->level = c0->level;

  P4EST_ASSERT (p4est_quadrant_is_valid (c0));
  P4EST_ASSERT (p4est_quadrant_is_valid (c1));
  P4EST_ASSERT (p4est_quadrant_is_valid (c2));
  P4EST_ASSERT (p4est_quadrant_is_valid (c3));
}

void
p4est_nearest_common_ancestor (const p4est_quadrant_t * q1,
                               const p4est_quadrant_t * q2,
                               p4est_quadrant_t * r)
{
  int32_t             exclorx, exclory;
  int32_t             maxclor, maxlevel;

  P4EST_ASSERT (p4est_quadrant_is_valid (q1));
  P4EST_ASSERT (p4est_quadrant_is_valid (q2));

  exclorx = q1->x ^ q2->x;
  exclory = q1->y ^ q2->y;
  maxclor = exclorx | exclory;
  maxlevel = P4EST_LOG2_32 (maxclor) + 1;

  r->x = q1->x & ~((1 << maxlevel) - 1);
  r->y = q1->y & ~((1 << maxlevel) - 1);
  r->level = (int8_t) P4EST_MIN (P4EST_MAXLEVEL - maxlevel,
                                 P4EST_MIN (q1->level, q2->level));

  P4EST_ASSERT (p4est_quadrant_is_valid (r));
}

void
p4est_nearest_common_ancestor_D (const p4est_quadrant_t * q1,
                                 const p4est_quadrant_t * q2,
                                 p4est_quadrant_t * r)
{
  p4est_quadrant_t    s1 = *q1;
  p4est_quadrant_t    s2 = *q2;

  P4EST_ASSERT (p4est_quadrant_is_valid (q1));
  P4EST_ASSERT (p4est_quadrant_is_valid (q2));

  /* first stage: promote the deepest one to the same level */
  while (s1.level > s2.level) {
    p4est_quadrant_parent (&s1, &s1);
  }
  while (s1.level < s2.level) {
    p4est_quadrant_parent (&s2, &s2);
  }

  /* second stage: simultaneously go through their parents */
  while (!p4est_quadrant_is_equal (&s1, &s2)) {
    p4est_quadrant_parent (&s1, &s1);
    p4est_quadrant_parent (&s2, &s2);
  }

  /* don't overwrite r's user_data */
  r->x = s1.x;
  r->y = s1.y;
  r->level = s1.level;

  P4EST_ASSERT (p4est_quadrant_is_valid (r));
}

void
p4est_quadrant_set_morton (p4est_quadrant_t * quadrant,
                           int8_t level, int32_t index)
{
  int8_t              i;

  P4EST_ASSERT (level < 16 && index < (1 << (2 * level)));

  quadrant->level = level;
  quadrant->x = 0;
  quadrant->y = 0;

  for (i = 0; i < level; ++i) {
    quadrant->x |= ((index & (1 << (2 * i))) >> i);
    quadrant->y |= ((index & (1 << (2 * i + 1))) >> (i + 1));
  }

  quadrant->x <<= (P4EST_MAXLEVEL - level);
  quadrant->y <<= (P4EST_MAXLEVEL - level);

  P4EST_ASSERT (p4est_quadrant_is_valid (quadrant));
}

void
p4est_quadrant_init_data (p4est_t * p4est, int32_t which_tree,
                          p4est_quadrant_t * quad, p4est_init_t init_fn)
{
  P4EST_ASSERT (p4est_quadrant_is_valid (quad));

  if (p4est->data_size > 0) {
    quad->user_data = p4est_mempool_alloc (p4est->user_data_pool);
  }
  else {
    quad->user_data = NULL;
  }
  if (init_fn != NULL) {
    init_fn (p4est, which_tree, quad);
  }
}

void
p4est_quadrant_free_data (p4est_t * p4est, p4est_quadrant_t * quad)
{
  P4EST_ASSERT (p4est_quadrant_is_valid (quad));

  if (p4est->data_size > 0) {
    p4est_mempool_free (p4est->user_data_pool, quad->user_data);
  }
  quad->user_data = NULL;
}

/* EOF p4est_algorithms.c */
