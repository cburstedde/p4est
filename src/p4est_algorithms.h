
#ifndef __P4EST_ALGORITHMS_H__
#define __P4EST_ALGORITHMS_H__

#include <p4est.h>

/*
 * comparison function in Morton ordering
 */
int                 p4est_quadrant_compare (const void *v1, const void *v2);

/*
 * set quadrant Morton indices based on linear index in uniform grid
 * uniform grid implies level < 16 and thus morton_xy < INT32_MAX
 */
void                p4est_quadrant_set_morton (p4est_quadrant_t * quadrant,
                                               int8_t level, int32_t index);

#endif /* !__P4EST_ALGORITHMS_H__ */
