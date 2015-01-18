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

typedef struct p4est_vtk_context p4est_vtk_context_t;

/** Cleanly destroy a \a p4est_vtk_context_t structure.
 *
 * This function closes all the file pointers and free's the \a
 * p4est_vtk_context_t structure.
 *
 * \param[in] context the \a p4est_vtk_context_t to be destroyed.
 */
void                p4est_vtk_context_destroy (p4est_vtk_context_t * context);

/** Write the p4est in VTK format.
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

/** Write the VTK header.
 *
 * Writing a VTK file is split into a few routines.
 * This allows there to be an arbitrary number of
 * fields.  The calling sequence would be something like
 *
 * \begincode
 * vtk_context = p4est_vtk_write_header (p4est, geom, 1., "output");
 * vtk_context = p4est_vtk_write_point_data (vtk_context, ...);
 * vtk_context = p4est_vtk_write_cell_data (vtk_context, ...);
 * ...
 * p4est_vtk_write_footer (vtk_context);
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
 * \return          On success, an opaque context (p4est_vtk_context_t) pointer
 *                  that must be passed to subsequent p4est_vtk calls.  It is
 *                  required to call p4est_vtk_write_footer eventually with
 *                  this value.  Returns NULL on error.
 */
p4est_vtk_context_t *p4est_vtk_write_header (p4est_t * p4est,
                                             p4est_geometry_t * geom,
                                             double scale,
                                             const char *filename);

/** Write VTK cell data.
 *
 * There are options to have this function write
 * the tree id, quadrant level, or MPI rank without explicit input data.
 *
 * Writing a VTK file is split into a few routines.
 * This allows there to be an arbitrary number of
 * fields.
 *
 * \param [in,out] cont    A vtk context created by p4est_vtk_write_header.
 * \param [in] write_tree  Boolean to determine if the tree id should be output.
 * \param [in] write_level Boolean to determine if the tree levels should be output.
 * \param [in] write_rank  Boolean to determine if the MPI rank should be output.
 * \param [in] wrap_rank   Number to wrap around the rank with a modulo operation.
 *                         Can be 0 for no wrapping.
 * \param [in] num_cell_scalars Number of cell scalar datasets to output.
 * \param [in] num_cell_vectors Number of cell vector datasets to output.
 *
 * The variable arguments need to be pairs of (fieldname, fieldvalues) where
 * the cell scalar pairs come first, followed by the cell vector pairs.  Each
 * 'fieldname' argument shall be a char string containing the name of the data
 * contained in the following 'fieldvalues'. Each of the 'fieldvalues'
 * arguments shall be an sc_array_t * holding double variables.  The number of
 * doubles in each sc_array must be exactly \a p4est->local_num_quadrants for
 * scalar data and \a 3*p4est->local_num_quadrants for vector data.
 *
 * TODO: For safety reasons, there must be a final parameter holding the context
 *       after the list of pairs.
 *
 * \return          On success, the context that has been passed in.
 *                  On failure, returns NULL and deallocates the context.
 */
p4est_vtk_context_t *p4est_vtk_write_cell_data (p4est_vtk_context_t * cont,
                                                int write_tree,
                                                int write_level,
                                                int write_rank,
                                                int wrap_rank,
                                                int num_cell_scalars,
                                                int num_cell_vectors, ...);

/** TODO: This is an alternate version of the varargs function above.
 * Works exactly the same otherwise.
 * TODO: not sure about this one yet, let's wait and think it over again soon.
 */
p4est_vtk_context_t *p4est_vtk_write_cell_data2 (p4est_vtk_context_t * cont,
                                                 int foo_bar_etc,
                                                 const char *filenames[],
                                                 const sc_array_t *
                                                 scalars[]);

/** Write VTK point data.
 *
 * Writing a VTK file is split into a few routines.
 * This allows there to be an arbitrary number of
 * fields.
 *
 * \param [in,out] cont          A vtk context created by p4est_vtk_write_header.
 * \param [in] num_point_scalars Number of point scalar datasets to output.
 * \param [in] num_point_vectors Number of point vector datasets to output.
 *
 * The variable arguments need to be pairs of (fieldname, fieldvalues) where
 * the point scalar pairs come first, followed by the point vector pairs.  Each
 * 'fieldname' argument shall be a char string containing the name of the data
 * contained in the following 'fieldvalues'. Each of the 'fieldvalues'
 * arguments shall be an sc_array_t * holding double variables. The number of
 * doubles in each sc_array must be exactly \a cont->num_nodes for scalar data
 * and \a 3*cont->num_nodes for vector data.
 *
 * \note \a cont->num_nodes is set in \b p4est_vtk_write_header based on the \a
 * scale parameter.
 *
 * \return          On success, the context that has been passed in.
 *                  On failure, returns NULL and deallocates the context.
 */
p4est_vtk_context_t *p4est_vtk_write_point_data (p4est_vtk_context_t * cont,
                                                 int num_point_scalars,
                                                 int num_point_vectors, ...);

/** Write VTK point data.
 *
 * This function exports custom point data to the vtk file; it is functionally
 * the same as \b p4est_vtk_write_point_data with the only difference being
 * that instead of a variable argument list, an initialized \a va_list is
 * passed as the last argument. The \a va_list is initialized from the variable
 * argument list of the calling function.
 *
 * \note This function is actually called from \b p4est_vtk_write_point_data
 * and does all of the work.
 *
 * \param [in,out] cont    A vtk context created by p4est_vtk_write_header.
 * \param [in] num_point_scalars Number of point scalar datasets to output.
 * \param [in] num_point_vectors Number of point vector datasets to output.
 * \param [in,out] ap      An initialized va_list used to access the
 *                         scalar/vector data.
 *
 * \return          On success, the context that has been passed in.
 *                  On failure, returns NULL and deallocates the context.
 */
p4est_vtk_context_t *p4est_vtk_write_point_datav (p4est_vtk_context_t * cont,
                                                  int num_point_scalars,
                                                  int num_point_vectors,
                                                  va_list ap);

/** Write VTK cell data.
 *
 * This function exports custom cell data to the vtk file; it is functionally
 * the same as \b p4est_vtk_write_cell_data with the only difference being
 * that instead of a variable argument list, an initialized \a va_list is
 * passed as the last argument. The \a va_list is initialized from the variable
 * argument list of the calling function.
 *
 * \note This function is actually called from \b p4est_vtk_write_cell_data
 * and does all of the work.
 *
 * \param [in,out] cont    A vtk context created by p4est_vtk_write_header.
 * \param [in] num_point_scalars Number of point scalar datasets to output.
 * \param [in] num_point_vectors Number of point vector datasets to output.
 * \param [in,out] ap      An initialized va_list used to access the
 *                         scalar/vector data.
 *
 * \return          On success, the context that has been passed in.
 *                  On failure, returns NULL and deallocates the context.
 */
p4est_vtk_context_t *p4est_vtk_write_cell_datav (p4est_vtk_context_t * cont,
                                                 int write_tree,
                                                 int write_level,
                                                 int write_rank,
                                                 int wrap_rank,
                                                 int num_cell_scalars,
                                                 int num_cell_vectors,
                                                 va_list ap);

/** Write a point scalar field to the VTU file.
 *
 * Writing a VTK file is split into a few routines.
 * This allows there to be an arbitrary number of fields.
 *
 * \param [in,out] cont    A vtk context created by p4est_vtk_write_header.
 * \param [in] scalar_name The name of the scalar field.
 * \param [in] values      The point values that will be written.
 *
 * \return          On success, the context that has been passed in.
 *                  On failure, returns NULL and deallocates the context.
 */
p4est_vtk_context_t *p4est_vtk_write_point_scalar (p4est_vtk_context_t * cont,
                                                   const char *scalar_name,
                                                   const sc_array_t * values);

/** Write a cell scalar field to the VTU file.
 *
 * Writing a VTK file is split into a few routines.
 * This allows there to be an arbitrary number of fields.
 *
 * \param [in,out] cont    A vtk context created by p4est_vtk_write_header.
 * \param [in] scalar_name The name of the scalar field.
 * \param [in] values      The cell values that will be written.
 *
 * \return          On success, the context that has been passed in.
 *                  On failure, returns NULL and deallocates the context.
 */
p4est_vtk_context_t *p4est_vtk_write_cell_scalar (p4est_vtk_context_t * cont,
                                                  const char *scalar_name,
                                                  const sc_array_t * values);

/** Write a 3-vector point field to the VTU file.
 *
 * Writing a VTK file is split into a few routines.
 * This allows there to be an arbitrary number of fields.
 *
 * \param [in,out] cont    A vtk context created by p4est_vtk_write_header.
 * \param [in] vector_name The name of the vector field.
 * \param [in] values      The point values that will be written.
 *
 * \return          On success, the context that has been passed in.
 *                  On failure, returns NULL and deallocates the context.
 */
p4est_vtk_context_t *p4est_vtk_write_point_vector (p4est_vtk_context_t * cont,
                                                   const char *vector_name,
                                                   const sc_array_t * values);

/** Write a 3-vector cell field to the VTU file.
 *
 * Writing a VTK file is split into a few routines.
 * This allows there to be an arbitrary number of fields.
 *
 * \param [in,out] cont    A vtk context created by p4est_vtk_write_header.
 * \param [in] vector_name The name of the vector field.
 * \param [in] values      The cell values that will be written.
 *
 * \return          On success, the context that has been passed in.
 *                  On failure, returns NULL and deallocates the context.
 */
p4est_vtk_context_t *p4est_vtk_write_cell_vector (p4est_vtk_context_t * cont,
                                                  const char *vector_name,
                                                  const sc_array_t * values);

/** Write the VTU footer and clean up.
 *
 * Writing a VTK file is split into a few routines.
 * This allows there to be an arbitrary number of
 * fields.  To write out two fields the
 * calling sequence would be something like
 *
 * \begincode
 * vtk_context = p4est_vtk_write_header (p4est, ..., "output");
 * vtk_context = p4est_vtk_write_point_data (vtk_context, ...);
 * vtk_context = p4est_vtk_write_cell_data (vtk_context, ...);
 * ...
 * p4est_vtk_write_footer (vtk_context);
 * \endcode
 *
 * This function writes the footer information to the vtk file and cleanly
 * destroys the vtk context.
 *
 * \param [in] cont Context is deallocated before the function returns.
 *
 * \return          This returns 0 if no error and -1 if there is an error.
 */
int                 p4est_vtk_write_footer (p4est_vtk_context_t * cont);

SC_EXTERN_C_END;

#endif /* !P4EST_VTK_H */
