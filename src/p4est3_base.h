/*
  This file is part of p4est, version 3.
  p4est is a C library to manage a collection (a forest) of multiple
  connected adaptive quadtrees or octrees in parallel.

  Copyright (C) 2019 individual authors
  Originally written by Carsten Burstedde, Lucas C. Wilcox, and Tobin Isaac

  p4est is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  p4est is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with p4est; if not, write to the Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/** \file p4est3_base.h
 *
 * The first include file to any p4est version 3 code.
 * This file pulls in important sc3 and system headers.
 * It is included indirectly in all p4est3 headers and source files.
 *
 * \ingroup p4est3
 */

/** \defgroup p4est3 p4est version 3
 *
 * Version 3 is an alternative, rewritten interface to p4est.
 * p4est version 3 depends on the version 3 rewrite of libsc.
 * While the standard versions up to 2 of p4est remain available,
 * This version 3 adds interfaces that are independent of the older ones.
 * Both versions may safely coexist in the same source file and program.
 *
 * For version 3, we are using new conventions on error returns, for
 * the construction and destruction of objects, and parallel shared memory.
 * p4est objects of version 2 may be wrapped in version 3 objects and thus
 * used from all version 3 code.  When the version 3 wrapper is destroyed,
 * the original version 2 object may be used by version 2 code as before.
 *
 * The interface of version 3 is shared between 2D and 3D.
 *
 * The reference counting mechanism of p4est is different from libsc.
 * We use \c p4est3_<object>_ref and \c p4est3_<object>_unref purely to
 * count up and down, respectively.  It is forbidden to count below one.
 * Destruction of an object is exclusively performed with \c
 * p4est3_<object>_destroy.  This, in turn, is only legal if the object's
 * reference count is one.  The intention is that the user is responsible to
 * know where in the program the last instance of an object lives.
 * We see reference counting as an additional safety feature, not to provide
 * additional flexibility by magical auto-destruction.
 */

#ifndef P4EST3_BASE_H
#define P4EST3_BASE_H

#include <sc3_alloc.h>
#include <sc3_error.h>
#include <p4est_config.h>

/** Integer type for topology counts, such as the number of trees. */
typedef int         p4est3_topidx;
#define P4EST3_TOPIDX_MAX INT_MAX       /**< Maximum value for \ref p4est3_topidx. */

/** Integer type for process-local object counts, such as local quadrants. */
typedef int         p4est3_locidx;
#define P4EST3_LOCIDX_MAX INT_MAX       /**< Maximum value for \ref p4est3_locidx. */
#define p4est3_loccut sc3_intcut        /**< Suitable partition cut function. */

/** Integer type for global sums of object counts,
    such as the total quadrants in a mesh. */
typedef long        p4est3_gloidx;
#define P4EST3_GLOIDX_MAX LONG_MAX      /**< Maximum value for \ref p4est3_gloidx. */
#define p4est3_glopow sc3_longpow       /**< Suitable integer power function. */
#define p4est3_glocut sc3_longcut       /**< Suitable partition cut function. */

#define p4est3_uint64cut sc3_uint64cut  /**< Suitable partition cut function. */

#ifdef __cplusplus
extern              "C"
{
#if 0
}
#endif
#endif

/* no function prototypes yet */

#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif /* !P4EST3_BASE_H */
