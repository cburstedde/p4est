
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

/*
 * initialize the user data of a quadrant which already has Morton indices
 */
void                p4est_quadrant_init_data (p4est_t * p4est,
                                              int32_t which_tree,
                                              p4est_quadrant_t * quad,
                                              p4est_init_t init_fn);

/*
 * free the user data of a quadrant
 */
void                p4est_quadrant_free_data (p4est_t * p4est,
                                              p4est_quadrant_t * quad);

#endif /* !__P4EST_ALGORITHMS_H__ */
