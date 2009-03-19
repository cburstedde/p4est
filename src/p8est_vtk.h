/*
  This file is part of p4est.
  p4est is a C library to manage a parallel collection of quadtrees and/or
  octrees.

  Copyright (C) 2007-2009 Carsten Burstedde, Lucas Wilcox.

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

#ifndef P8EST_VTK_H
#define P8EST_VTK_H

#include <p8est_geometry.h>

SC_EXTERN_C_BEGIN;

extern double       p8est_vtk_default_scale;
extern bool         p8est_vtk_default_write_tree;
extern bool         p8est_vtk_default_write_rank;
extern int          p8est_vtk_default_wrap_rank;

/** This will write out the p8est in VTK format.
 * The p8est_vtk_default_* variables will be honored.
 *
 * This is a convenience function for the special
 * case of writing out the MPI rank only.  Note this
 * function will abort if there is a file error.
 *
 * \param p8est     The p8est to be written.
 * \param geom      A p8est_geometry_t structure or NULL for identity.
 * \param baseName  The first part of the name which will have
 *                  the proc number appended to it (i.e., the
 *                  output file will be baseName_procNum.vtu).
 */
void                p8est_vtk_write_file (p8est_t * p8est,
                                          p8est_geometry_t * geom,
                                          const char *baseName);

/** This will write out the p8est and any numner of point fields.
 * The p8est_vtk_default_* variables will be honored.
 *
 * This is a convenience function that will abort if there is a file error.
 *
 * \param p8est     The p8est to be written.
 * \param geom      A p8est_geometry_t structure or NULL for identity.
 * \param baseName  The first part of the name which will have
 *                  the proc number appended to it (i.e., the
 *                  output file will be baseName_procNum.vtu).
 * \param num_scalars   Number of scalar fields to write.
 * \param num_vectors   Number of vector fields to write.
 *
 * The variable arguments need to be pairs of (fieldname, fieldvalues)
 * where the scalars come first, then the vectors.
 */
void                p8est_vtk_write_all (p8est_t * p8est,
                                         p8est_geometry_t * geom,
                                         int num_scalars, int num_vectors,
                                         const char *baseName, ...);

/** This will write the header of the vtu file.
 *
 * Writing a VTK file is split into a couple of routines.
 * The allows there to be an arbitrary number of
 * fields.  The calling sequence would be something like
 *
 * \begincode
 * p8est_vtk_write_header(p8est, 1., true, "output");
 * p8est_vtk_write_point_scalar (...);
 * ...
 * p8est_vtk_write_footer(p8est, "output");
 * \endcode
 *
 * \param p8est     The p8est to be written.
 * \param geom      A p8est_geometry_t structure or NULL for identity.
 * \param scale     The relative length factor of the quadrants.
 *                  Use 1.0 to fit quadrants exactly, less to create gaps.
 * \param write_tree    Boolean to determine if the tree id should be output.
 * \param write_rank    Boolean to determine if the MPI rank should be output.
 * \param wrap_rank Number to wrap around the rank with a modulo operation.
 *                  Can be 0 for no wrapping.
 * \param pointScalars  Comma-separated list of point scalar fields, or NULL.
 * \param pointVectors  Comma-separated list of point vector fields, or NULL.
 * \param baseName  The first part of the name which will have
 *                  the proc number appended to it (i.e., the
 *                  output file will be baseName_procNum.vtu).
 *
 * \return          This returns 0 if no error and -1 if there is an error.
 */
int                 p8est_vtk_write_header (p8est_t * p8est,
                                            p8est_geometry_t * geom,
                                            double scale, bool write_tree,
                                            bool write_rank, int wrap_rank,
                                            const char *pointScalars,
                                            const char *pointVectors,
                                            const char *baseName);

/** This will write a scalar field to the vtu file.
 *
 * It is good practice to make sure that the scalar field also
 * exists in the comma separated string \a pointscalars passed
 * to \c p8est_vtk_write_header.
 *
 * Writing a VTK file is split into a couple of routines.
 * The allows there to be an arbitrary number of fields.
 *
 * \param p8est     The p8est to be written.
 * \param geom      A p8est_geometry_t structure or NULL for identity.
 * \param baseName  The first part of the name which will have
 *                  the proc number appended to it (i.e., the
 *                  output file will be baseName_procNum.vtu).
 * \param scalarName The name of the scalar field.
 * \param values    The point values that will be written.
 *
 * \return          This returns 0 if no error and -1 if there is an error.
 */
int                 p8est_vtk_write_point_scalar (p8est_t * p8est,
                                                  p8est_geometry_t * geom,
                                                  const char *baseName,
                                                  const char *scalarName,
                                                  const double *values);

/** This will write a 3-vector field to the vtu file.
 *
 * It is good practice to make sure that the vector field also
 * exists in the comma separated string \a pointvectors passed
 * to \c p8est_vtk_write_header.
 *
 * Writing a VTK file is split into a couple of routines.
 * The allows there to be an arbitrary number of fields.
 *
 * \param p8est     The p8est to be written.
 * \param geom      A p8est_geometry_t structure or NULL for identity.
 * \param baseName  The first part of the name which will have
 *                  the proc number appended to it (i.e., the
 *                  output file will be baseName_procNum.vtu).
 * \param vectorName The name of the vector field.
 * \param values    The point values that will be written.
 *
 * \return          This returns 0 if no error and -1 if there is an error.
 */
int                 p8est_vtk_write_point_vector (p8est_t * p8est,
                                                  p8est_geometry_t * geom,
                                                  const char *baseName,
                                                  const char *vectorName,
                                                  const double *values);

/** This will write the footer of the vtu file.
 *
 * Writing a VTK file is split into a couple of routines.
 * The allows there to be an arbitrary number of
 * fields.  To write out two fields the
 * calling sequence would be something like
 *
 * \begincode
 * p8est_vtk_write_header(p8est, "output");
 * p8est_vtk_write_footer(p8est, "output");
 * \endcode
 *
 * \param p8est     The p8est to be written.
 * \param baseName  The first part of the name which will have
 *                  the proc number appended to it (i.e., the
 *                  output file will be baseName_procNum.vtu).
 *
 * \return          This returns 0 if no error and -1 if there is an error.
 */
int                 p8est_vtk_write_footer (p8est_t * p8est,
                                            const char *baseName);

SC_EXTERN_C_END;

#endif /* !P8EST_VTK_H */
