
#include <p4est_algorithms.h>
#include <p4est_base.h>

int
p4est_quadrant_compare (const void *v1, const void *v2)
{
  const p4est_quadrant_t *q1 = v1;
  const p4est_quadrant_t *q2 = v2;

  int32_t             exclorx, exclory;

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
p4est_quadrant_is_ancestor (const p4est_quadrant_t *q,
                            const p4est_quadrant_t *r)
{
  int32_t             exclorx;
  int32_t             exclory;

  if (q->level > r->level) {
    return 0;
  }

  exclorx = (q->x ^ r->x) >> (P4EST_MAXLEVEL - q->level);
  exclory = (q->y ^ r->y) >> (P4EST_MAXLEVEL - q->level);

  return (exclorx == 0 && exclory == 0);
}

void
p4est_quadrant_parent (const p4est_quadrant_t * q, p4est_quadrant_t * r)
{
  P4EST_ASSERT (q->level > 0);

  r->x = q->x & ~(1 << (P4EST_MAXLEVEL - q->level));
  r->y = q->y & ~(1 << (P4EST_MAXLEVEL - q->level));
  r->level = (int8_t) (q->level - 1);
}

void
p4est_nearest_common_ancestor (const p4est_quadrant_t * q1,
                               const p4est_quadrant_t * q2,
                               p4est_quadrant_t * r)
{
  p4est_quadrant_t    s1 = *q1;
  p4est_quadrant_t    s2 = *q2;

  /* first stage: promote the deepest one to the same level */
  while (s1.level > s2.level) {
    p4est_quadrant_parent (&s1, &s1);
  }
  while (s1.level < s2.level) {
    p4est_quadrant_parent (&s2, &s2);
  }

  /* second stage: simultaneously go through their parents */
  while (p4est_quadrant_compare (&s1, &s2) != 0) {
    p4est_quadrant_parent (&s1, &s1);
    p4est_quadrant_parent (&s2, &s2);
  }

  *r = s1;
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
}

void
p4est_quadrant_init_data (p4est_t * p4est, int32_t which_tree,
                          p4est_quadrant_t * quad, p4est_init_t init_fn)
{
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
  if (p4est->data_size > 0) {
    p4est_mempool_free (p4est->user_data_pool, quad->user_data);
  }
  quad->user_data = NULL;
}

/* EOF p4est_algorithms.c */
