/*
  This file is part of p4est.
  p4est is a C library to manage a parallel collection of quadtrees and/or
  octrees.

  Copyright (C) 2007,2008 Carsten Burstedde, Lucas Wilcox.

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef __P4EST_VTK_H__
#define __P4EST_VTK_H__

#include <p4est.h>

/** This will write out one field in VTK format.
 *
 * This is a convenience function for the special
 * case of writing out one field.  Note this
 * function will abort if there is a file error.
 *
 *  \param p4est     The p4est to be outputted.
 *  \param baseName  The first part of the name which will have
 *                   the proc number appended to it (i.e., the
 *                   output file will be baseName_procNum.vtu).
 */
void                p4est_vtk_write_file (p4est_t * p4est,
                                          const char *baseName);

/** This will write the header of the vtu file.
 *
 * Writing a VTK file is split into a couple of routines.
 * The allows there to be an arbitrary number of
 * fields.  To write out two fields the
 * calling sequence would be something like
 *
 * \begincode
 * p4est_vtk_write_header(p4est, "output");
 * p4est_vtk_write_footer(p4est, "output");
 * \endcode
 *
 *  \param p4est    The p4est to be outputted.
 *  \param baseName The first part of the name which will have
 *                  the proc number appended to it (i.e., the
 *                  output file will be baseName_procNum.vtu).
 *
 *  \return         This returns 0 if no error and -1 if there is an error.
 */
int                 p4est_vtk_write_header (p4est_t * p4est,
                                            const char *baseName);

/** This will write the footer of the vtu file.
 *
 * Writing a VTK file is split into a couple of routines.
 * The allows there to be an arbitrary number of
 * fields.  To write out two fields the
 * calling sequence would be something like
 *
 * \begincode
 * p4est_vtk_write_header(p4est, "output");
 * p4est_vtk_write_footer(p4est, "output");
 * \endcode
 *
 *  \param p4est    The p4est to be outputted.
 *  \param baseName The first part of the name which will have
 *                  the proc number appended to it (i.e., the
 *                  output file will be baseName_procNum.vtu).
 *
 *  \return         This returns 0 if no error and -1 if there is an error.
 */
int                 p4est_vtk_write_footer (p4est_t * p4est,
                                            const char *baseName);

#endif /* !__P4EST_VTK_H__ */
