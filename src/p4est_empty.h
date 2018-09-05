/*
  This file is part of p4est.
  p4est is a C library to manage a collection (a forest) of multiple
  connected adaptive quadtrees or octrees in parallel.

  Copyright (C) 2010 The University of Texas System
  Additional copyright (C) 2011 individual authors
  Written by Carsten Burstedde, Lucas C. Wilcox, and Tobin Isaac

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

/******************************************************************************
 * The files:                                                                 *
 *                                                                            *
 *   p4est_empty.h p4est_empty.c                                              *
 *   p8est_empty.h p8est_empty.c p4est_to_p8est_empty.h                       *
 *   p6est_empty.h p6est_empty.c                                              *
 *                                                                            *
 * are intended to allow users to add functionality to the p4est library for  *
 * specific applications.                                                     *
 *                                                                            *
 * A possible use case would be for an application to distribute its version  *
 * of these files and copy them into the p4est source directory prior to      *
 * configuration and build of p4est.  By doing this, applications can use     *
 * distributed tarballs of p4est without having to re-bootstrap p4est.        *
 ******************************************************************************/

#ifndef P4EST_EMPTY_H
#define P4EST_EMPTY_H

#include <p4est_base.h>

SC_EXTERN_C_BEGIN;

/* This is a dummy .h file that the user can replace as needed */

/** This function does nothing. */
void                p4est_empty_noop (void);

SC_EXTERN_C_END;

#endif /* !P4EST_EMPTY_H */
