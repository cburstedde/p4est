/*
  This file is part of p4est.
  p4est is a C library to manage a collection (a forest) of multiple
  connected adaptive quadtrees or octrees in parallel.

  Copyright (C) 2010 The University of Texas System
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

/** \file p8est_vtk.h
 *
 * Routines for printing a forest and associated fields to vtk format.
 *
 * \ingroup p8est
 */

#ifndef P8EST_VTK_H
#define P8EST_VTK_H

#include <p8est_geometry.h>
#include <p8est.h>

SC_EXTERN_C_BEGIN;

/** This writes out the p8est in VTK format.
 *
 * This is a convenience function for the special
 * case of writing out the tree id, quadrant level, and MPI rank only.
 * One file is written per MPI rank, and one meta file on rank 0.
 * This function will abort if there is a file error.
 *
 * \param [in] p8est    The p8est to be written.
 * \param [in] geom     A p8est_geometry_t structure or NULL for vertex space.
 * \param [in] filename The first part of the file name which will have the
 *                      MPI rank appended to it: The output file will be
 *                      filename_procNum.vtu, and the meta file filename.pvtu.
 */
void                p8est_vtk_write_file (p8est_t * p8est,
                                          p8est_geometry_t * geom,
                                          const char *filename);

/** This writes out the p8est and any number of point fields in VTK format.
 *
 * This is a convenience function that will abort if there is a file error.
 *
 * \param [in] p8est    The p8est to be written.
 * \param [in] geom     A p8est_geometry_t structure or NULL for vertex space.
 * \param [in] scale    Double value between 0 and 1 to scale each quadrant.
 * \param [in] write_tree   Include the tree id as output field.
 * \param [in] write_level  Include the tree levels as output field.
 * \param [in] write_rank   Include the MPI rank as output field.
 * \param [in] wrap_tree    The MPI rank is written module wrap_tree, or 0.
 * \param filename      First part of the name, see p8est_vtk_write_file.
 * \param num_scalars   Number of scalar fields to write.
 * \param num_vectors   Number of vector fields to write.
 *
 * The variable arguments need to be pairs of (fieldname, fieldvalues)
 * where the scalars come first, then the vectors.
 */
void                p8est_vtk_write_all (p8est_t * p8est,
                                         p8est_geometry_t * geom,
                                         double scale,
                                         int write_tree, int write_level,
                                         int write_rank, int wrap_rank,
                                         int num_scalars, int num_vectors,
                                         const char *filename, ...);

void                p8est_vtk_writeAll (p8est_t * p8est,
                                        p8est_geometry_t * geom,
                                        double scale,
                                        int write_tree, int write_level,
                                        int write_rank, int wrap_rank,
                                        int num_cell_scalars,
                                        int num_cell_vectors,
                                        int num_point_scalars,
                                        int num_point_vectors,
                                        const char * filename);

/** This will write the header of the vtu file.
 *
 * Writing a VTK file is split into a couple of routines.
 * The allows there to be an arbitrary number of
 * fields.  The calling sequence would be something like
 *
 * \begincode
 * p8est_vtk_write_header(p8est, geom, 1., 1, 1, 1, 0, "output");
 * p8est_vtk_write_point_scalar (...);
 * ...
 * p8est_vtk_write_footer(p8est, "output");
 * \endcode
 *
 * \param p8est     The p8est to be written.
 * \param geom      A p8est_geometry_t structure or NULL for vertex space.
 * \param scale     The relative length factor of the quadrants.
 *                  Use 1.0 to fit quadrants exactly, less to create gaps.
 * \param write_tree    Boolean to determine if the tree id should be output.
 * \param write_level   Boolean to determine if the tree levels should be output.
 * \param write_rank    Boolean to determine if the MPI rank should be output.
 * \param wrap_rank Number to wrap around the rank with a modulo operation.
 *                  Can be 0 for no wrapping.
 * \param point_scalars  Comma-separated list of point scalar fields, or NULL.
 * \param point_vectors  Comma-separated list of point vector fields, or NULL.
 * \param filename  The first part of the name which will have
 *                  the proc number appended to it (i.e., the
 *                  output file will be filename_procNum.vtu).
 *
 * \return          This returns 0 if no error and -1 if there is an error.
 */
int                 p8est_vtk_write_header (p8est_t * p8est,
                                            p8est_geometry_t * geom,
                                            double scale,
                                            int write_tree, int write_level,
                                            int write_rank, int wrap_rank,
                                            const char *point_scalars,
                                            const char *point_vectors,
                                            const char *filename);

int                 p8est_vtk_writeHeader (p8est_t * p8est,
                                           p8est_geometry_t * geom,
                                           double scale,
                                           int write_tree, int write_level,
                                           int write_rank, int wrap_rank,
                                           const char * cell_scalars,
                                           const char * cell_vectors,
                                           const char * filename);

int                 p8est_vtk_write_cell_scalar (p8est_t * p8est,
                                                 p8est_geometry_t * geom,
                                                 const char * filename,
                                                 const char * scalar_name,
                                                 const double * values);

int                 p8est_vtk_write_cell_vector (p8est_t * p8est,
                                                 p8est_geometry_t * geom,
                                                 const char * filename,
                                                 const char * vector_name,
                                                 const double * values);

int                 p8est_vtk_switch_to_point_data (p8est_t * p8est,
                                                    p8est_geometry_t * geom,
                                                    int write_tree,
                                                    int write_level,
                                                    int write_rank,
                                                    int num_cell_data,
                                                    const char * point_scalars,
                                                    const char * point_vectors,
                                                    const char * filename);

/** This will write a scalar field to the vtu file.
 *
 * It is good practice to make sure that the scalar field also
 * exists in the comma separated string \a point_scalars passed
 * to \c p8est_vtk_write_header.
 *
 * Writing a VTK file is split into a couple of routines.
 * The allows there to be an arbitrary number of fields.
 *
 * \param p8est     The p8est to be written.
 * \param geom      A p8est_geometry_t structure or NULL for vertex space.
 * \param filename  The first part of the name which will have
 *                  the proc number appended to it (i.e., the
 *                  output file will be filename_procNum.vtu).
 * \param scalar_name The name of the scalar field.
 * \param values    The point values that will be written.
 *
 * \return          This returns 0 if no error and -1 if there is an error.
 */
int                 p8est_vtk_write_point_scalar (p8est_t * p8est,
                                                  p8est_geometry_t * geom,
                                                  const char *filename,
                                                  const char *scalar_name,
                                                  const double *values);

/** This will write a 3-vector field to the vtu file.
 *
 * It is good practice to make sure that the vector field also
 * exists in the comma separated string \a point_vectors passed
 * to \c p8est_vtk_write_header.
 *
 * Writing a VTK file is split into a couple of routines.
 * The allows there to be an arbitrary number of fields.
 *
 * \param p8est     The p8est to be written.
 * \param geom      A p8est_geometry_t structure or NULL for vertex space.
 * \param filename  The first part of the name which will have
 *                  the proc number appended to it (i.e., the
 *                  output file will be filename_procNum.vtu).
 * \param vector_name The name of the vector field.
 * \param values    The point values that will be written.
 *
 * \return          This returns 0 if no error and -1 if there is an error.
 */
int                 p8est_vtk_write_point_vector (p8est_t * p8est,
                                                  p8est_geometry_t * geom,
                                                  const char *filename,
                                                  const char *vector_name,
                                                  const double *values);

/** This will write the footer of the vtu file.
 *
 * Writing a VTK file is split into a couple of routines.
 * The allows there to be an arbitrary number of
 * fields.  To write out two fields the
 * calling sequence would be something like
 *
 * \begincode
 * p8est_vtk_write_header(p8est, ..., "output");
 * p8est_vtk_write_footer(p8est, "output");
 * \endcode
 *
 * \param p8est     The p8est to be written.
 * \param filename  The first part of the name which will have
 *                  the proc number appended to it (i.e., the
 *                  output file will be filename_procNum.vtu).
 *
 * \return          This returns 0 if no error and -1 if there is an error.
 */
int                 p8est_vtk_write_footer (p8est_t * p8est,
                                            const char *filename);

SC_EXTERN_C_END;

#endif /* !P8EST_VTK_H */
