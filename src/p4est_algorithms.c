
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
}

/* EOF p4est_algorithms.c */
