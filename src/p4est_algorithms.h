
#ifndef __P4EST_ALGORITHMS_H__
#define __P4EST_ALGORITHMS_H__

#include <p4est.h>

/** Compare two quadrants in their Morton ordering
 * \return Returns < 0 if \a v1 < \a v2,
 *                   0 if \a v1 == \a v2,
 *                 > 0 if \a v1 > \a v2
 */
int                 p4est_quadrant_compare (const void *v1, const void *v2);

/** Test if a quadrant is an ancestor of another quadrant
 * \param [in] q Quadrant to be tested.
 * \param [in] r Descendent quadrant.
 * \return True if \a q is an ancestor of \a r or if \a q == \a r.
 */
int                 p4est_quadrant_is_ancestor (const p4est_quadrant_t *q,
                                                const p4est_quadrant_t *r);

/** Compute the parent of a quadrant
 * \param [in]  q Input quadrant.
 * \param [out] r Existing quadrant whose Morton index will be filled
 *                with the Morton index of the parent of \a q.
 *                Its user_data will be untouched.
 * \note \a q == \a r is permitted.
 */
void                p4est_quadrant_parent (const p4est_quadrant_t * q,
                                           p4est_quadrant_t * r);

/** Computes the nearest common ancestor of two quadrants
 * \param [in]  q1 First input quadrant.
 * \param [in]  q2 Second input quadrant.
 * \param [out] r Existing quadrant whose Morton index will be filled.
 *                Its user_data will be untouched.
 * \note \a q1, \a q2, \a r may point to the same quadrant.
 */
void                p4est_nearest_common_ancestor (const p4est_quadrant_t *
                                                   q1,
                                                   const p4est_quadrant_t *
                                                   q2, p4est_quadrant_t * r);

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
