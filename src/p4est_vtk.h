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

/** \file p4est_vtk.h
 *
 * Routines for printing a forest and associated fields to vtk format.
 *
 * \ingroup p4est
 */

#ifndef P4EST_VTK_H
#define P4EST_VTK_H

#include <p4est_geometry.h>
#include <p4est.h>

SC_EXTERN_C_BEGIN;

/** This writes out the p4est in VTK format.
 *
 * This is a convenience function for the special case of writing out
 * the tree id, quadrant level, and MPI rank only.
 * One file is written per MPI rank, and one meta file on rank 0.
 * This function will abort if there is a file error.
 *
 * \param [in] p4est    The p4est to be written.
 * \param [in] geom     A p4est_geometry_t structure or NULL for vertex space
 *                      as defined by p4est->connectivity.
 * \param [in] filename The first part of the file name which will have the
 *                      MPI rank appended to it: The output file will be
 *                      filename_rank.vtu, and the meta file filename.pvtu.
 */
void                p4est_vtk_write_file (p4est_t * p4est,
                                          p4est_geometry_t * geom,
                                          const char *filename);

/** This writes out the p4est and any number of point fields in VTK format.
 *
 * This is a convenience function that will abort if there is a file error.
 *
 * \param [in] p4est    The p4est to be written.
 * \param [in] geom     A p4est_geometry_t structure or NULL for vertex space.
 * \param [in] scale    Double value between 0 and 1 to scale each quadrant.
 * \param [in] write_tree  Include the tree id as cell output field.
 * \param [in] write_level Include the tree levels as cell output field.
 * \param [in] write_rank  Include the MPI rank as cell output field.
 * \param [in] wrap_rank   The MPI rank is written module wrap_tree, or 0.
 * \param [in] num_point_scalars  Number of point scalar fields to write.
 * \param [in] num_point_vectors  Number of point vector fields to write.
 * \param [in] filename           First part of the name;
 *                                see p4est_vtk_write_file.
 *
 * The variable arguments need to be pairs of (fieldname, fieldvalues)
 * where the point scalars come first, followed by the point vectors.
 *
 * \note this function only supports point data; if cell data is required see
 * p4est_vtk_write_cell_data.
 */
void                p4est_vtk_write_all (p4est_t * p4est,
                                         p4est_geometry_t * geom,
                                         double scale,
                                         int write_tree, int write_level,
                                         int write_rank, int wrap_rank,
                                         int num_point_scalars,
                                         int num_point_vectors,
                                         const char *filename, ...);

/** This will write the header of the vtu file.
 *
 * Writing a VTK file is split into a few routines.
 * This allows there to be an arbitrary number of
 * fields.  The calling sequence would be something like
 *
 * \begincode
 * p4est_vtk_write_header (p4est, geom, 1., "output");
 * p4est_vtk_write_point_data (...);
 * p4est_vtk_write_cell_data (...);
 * ...
 * p4est_vtk_write_footer (p4est, "output");
 * \endcode
 *
 * \param p4est     The p4est to be written.
 * \param geom      A p4est_geometry_t structure or NULL for vertex space.
 * \param scale     The relative length factor of the quadrants.
 *                  Use 1.0 to fit quadrants exactly, less to create gaps.
 * \param filename  The first part of the name which will have
 *                  the proc number appended to it (i.e., the
 *                  output file will be filename_rank.vtu).
 *
 * \return          This returns 0 if no error and -1 if there is an error.
 */
int                 p4est_vtk_write_header (p4est_t * p4est,
                                            p4est_geometry_t * geom,
                                            double scale,
                                            const char *filename);

/** This will write custom cell data to the vtu file.
 *
 * There are options to have this function write
 * the tree id, quadrant level, or MPI rank without explicit input data.
 *
 * Writing a VTK file is split into a few routines.
 * This allows there to be an arbitrary number of
 * fields.
 *
 * \param p4est     The p4est to be written.
 * \param geom      A p4est_geometry_t structure or NULL for vertex space.
 * \param write_tree    Boolean to determine if the tree id should be output.
 * \param write_level   Boolean to determine if the tree levels should be output.
 * \param write_rank    Boolean to determine if the MPI rank should be output.
 * \param wrap_rank Number to wrap around the rank with a modulo operation.
 *                  Can be 0 for no wrapping.
 * \param num_cell_scalars Number of cell scalar datasets to output.
 * \param num_cell_vectors Number of cell vector datasets to output.
 * \param filename  The first part of the name which will have
 *                  the proc number appended to it (i.e., the
 *                  output file will be filename_rank.vtu).
 *
 * The variable arguments need to be pairs of (fieldname, fieldvalues)
 * where the cell scalar pairs come first, followed by the cell vector pairs.
 *
 * \return          This returns 0 if no error and -1 if there is an error.
 */
int                 p4est_vtk_write_cell_data (p4est_t * p4est,
                                               p4est_geometry_t * geom,
                                               const int write_tree,
                                               const int write_level,
                                               const int write_rank,
                                               const int wrap_rank,
                                               const int num_cell_scalars,
                                               const int num_cell_vectors,
                                               const char *filename, ...);

/** This will write custom point data to the vtu file.
 *
 * Writing a VTK file is split into a few routines.
 * This allows there to be an arbitrary number of
 * fields.
 *
 * \param p4est     The p4est to be written.
 * \param geom      A p4est_geometry_t structure or NULL for vertex space.
 * \param num_point_scalars Number of point scalar datasets to output.
 * \param num_point_vectors Number of point vector datasets to output.
 * \param filename  The first part of the name which will have
 *                  the proc number appended to it (i.e., the
 *                  output file will be filename_rank.vtu).
 *
 * The variable arguments need to be pairs of (fieldname, fieldvalues)
 * where the point scalar pairs come first, followed by the point vector pairs.
 *
 * \return          This returns 0 if no error and -1 if there is an error.
 */
int                 p4est_vtk_write_point_data (p4est_t * p4est,
                                                p4est_geometry_t * geom,
                                                const int num_point_scalars,
                                                const int num_point_vectors,
                                                const char *filename, ...);

/** TODO: Please document this function and
 *        add analogous function for cell data.
 */
int                 p4est_vtk_write_point_datav (p4est_t * p4est,
                                                 p4est_geometry_t * geom,
                                                 const int num_point_scalars,
                                                 const int num_point_vectors,
                                                 const char *filename,
                                                 va_list ap);

/** This will write a point scalar field to the vtu file.
 *
 * Writing a VTK file is split into a few routines.
 * This allows there to be an arbitrary number of fields.
 *
 * \param p4est     The p4est to be written.
 * \param geom      A p4est_geometry_t structure or NULL for vertex space.
 * \param filename  The first part of the name which will have
 *                  the proc number appended to it (i.e., the
 *                  output file will be filename_rank.vtu).
 * \param scalar_name The name of the scalar field.
 * \param values    The point values that will be written.
 *
 * \return          This returns 0 if no error and -1 if there is an error.
 */
int                 p4est_vtk_write_point_scalar (p4est_t * p4est,
                                                  p4est_geometry_t * geom,
                                                  const char *filename,
                                                  const char *scalar_name,
                                                  const double *values);

/** This will write a cell scalar field to the vtu file.
 *
 * Writing a VTK file is split into a few routines.
 * This allows there to be an arbitrary number of fields.
 *
 * \param p4est     The p4est to be written.
 * \param geom      A p4est_geometry_t structure or NULL for vertex space.
 * \param filename  The first part of the name which will have
 *                  the proc number appended to it (i.e., the
 *                  output file will be filename_rank.vtu).
 * \param scalar_name The name of the scalar field.
 * \param values    The cell values that will be written.
 *
 * \return          This returns 0 if no error and -1 if there is an error.
 */
int                 p4est_vtk_write_cell_scalar (p4est_t * p4est,
                                                 p4est_geometry_t * geom,
                                                 const char *filename,
                                                 const char *scalar_name,
                                                 const double *values);

/** This will write a 3-vector point field to the vtu file.
 *
 * Writing a VTK file is split into a few routines.
 * This allows there to be an arbitrary number of fields.
 *
 * \param p4est     The p4est to be written.
 * \param geom      A p4est_geometry_t structure or NULL for vertex space.
 * \param filename  The first part of the name which will have
 *                  the proc number appended to it (i.e., the
 *                  output file will be filename_rank.vtu).
 * \param vector_name The name of the vector field.
 * \param values    The point values that will be written.
 *
 * \return          This returns 0 if no error and -1 if there is an error.
 */
int                 p4est_vtk_write_point_vector (p4est_t * p4est,
                                                  p4est_geometry_t * geom,
                                                  const char *filename,
                                                  const char *vector_name,
                                                  const double *values);

/** This will write a 3-vector cell field to the vtu file.
 *
 * Writing a VTK file is split into a few routines.
 * This allows there to be an arbitrary number of fields.
 *
 * \param p4est     The p4est to be written.
 * \param geom      A p4est_geometry_t structure or NULL for vertex space.
 * \param filename  The first part of the name which will have
 *                  the proc number appended to it (i.e., the
 *                  output file will be filename_rank.vtu).
 * \param vector_name The name of the vector field.
 * \param values    The cell values that will be written.
 *
 * \return          This returns 0 if no error and -1 if there is an error.
 */
int                 p4est_vtk_write_cell_vector (p4est_t * p4est,
                                                 p4est_geometry_t * geom,
                                                 const char *filename,
                                                 const char *vector_name,
                                                 const double *values);

/** This will write the footer of the vtu file.
 *
 * Writing a VTK file is split into a few routines.
 * This allows there to be an arbitrary number of
 * fields.  To write out two fields the
 * calling sequence would be something like
 *
 * \begincode
 * p4est_vtk_write_header (p4est, ..., "output");
 * p4est_vtk_write_footer (p4est, "output");
 * \endcode
 *
 * \param p4est     The p4est to be written.
 * \param filename  The first part of the name which will have
 *                  the proc number appended to it (i.e., the
 *                  output file will be filename_rank.vtu).
 *
 * \return          This returns 0 if no error and -1 if there is an error.
 */
int                 p4est_vtk_write_footer (p4est_t * p4est,
                                            const char *filename);

SC_EXTERN_C_END;

#endif /* !P4EST_VTK_H */
