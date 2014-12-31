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

//typedef struct p8est_vtk_context p8est_vtk_context_t;
/** \struct p8est_vtk_context_t
 * Opaque context type for writing VTK output with multiple function calls.
 *
 * This structure holds all the information needed for the p8est vtk context.
 * It is used to relay necessary vtk information to the \b p8est_vtk_write_*
 * functions. This structure is initialized by \b p8est_vtk_write_header and
 * destroyed by \b p8est_vtk_write_footer; it can also be destroyed manually
 * using the \b p8est_vtk_context_destroy function if necessary.
 *
 * The \a p4est member is a pointer to the local p8est.
 * The \a geom member is a pointer to the geometry used to create the p8est.
 * The \a num_nodes member holds the number of nodes present in the vtk output;
 * this is determined in \b p8est_vtk_write_header using the \a scale parameter
 * and is used to assure the proper number of point variables are provided.
 * The \a filename member holds the vtk file basename: for error reporting.
 * The \a vtufilename, \a pvtufilename, and \a visitfilename members are the
 * vtk file names.
 * The \a vtufile, \a pvtufile, and \a visitfile members are the vtk file
 * pointers; opened by \b p8est_vtk_write_header and closed by \b
 * p8est_vtk_write_footer.
 *
 */
typedef struct p8est_vtk_context
{
  /* TODO: Add members as needed */
  p8est_t            *p4est;
  p8est_geometry_t   *geom;
  p4est_locidx_t      num_nodes;
  char               *filename;
  char                vtufilename[BUFSIZ], pvtufilename[BUFSIZ],
    visitfilename[BUFSIZ];
  FILE               *vtufile, *pvtufile, *visitfile;
}
p8est_vtk_context_t;

/** This writes out the p8est in VTK format.
 *
 * This is a convenience function for the special case of writing out
 * the tree id, quadrant level, and MPI rank only.
 * One file is written per MPI rank, and one meta file on rank 0.
 * This function will abort if there is a file error.
 *
 * \param [in] p8est    The p8est to be written.
 * \param [in] geom     A p8est_geometry_t structure or NULL for vertex space
 *                      as defined by p8est->connectivity.
 * \param [in] filename The first part of the file name which will have the
 *                      MPI rank appended to it: The output file will be
 *                      filename_rank.vtu, and the meta file filename.pvtu.
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
 * \param [in] write_tree  Include the tree id as cell output field.
 * \param [in] write_level Include the tree levels as cell output field.
 * \param [in] write_rank  Include the MPI rank as cell output field.
 * \param [in] wrap_rank   The MPI rank is written module wrap_tree, or 0.
 * \param [in] num_point_scalars  Number of point scalar fields to write.
 * \param [in] num_point_vectors  Number of point vector fields to write.
 * \param [in] filename           First part of the name;
 *                                see p8est_vtk_write_file.
 *
 * The variable arguments need to be pairs of (fieldname, fieldvalues)
 * where the point scalars come first, followed by the point vectors.
 *
 * \note this function only supports point data; if cell data is required see
 * p8est_vtk_write_cell_data.
 */
void                p8est_vtk_write_all (p8est_t * p8est,
                                         p8est_geometry_t * geom,
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
 * p8est_vtk_write_header (p8est, geom, 1., "output");
 * p8est_vtk_write_point_data (...);
 * p8est_vtk_write_cell_data (...);
 * ...
 * p8est_vtk_write_footer (p8est, "output");
 * \endcode
 *
 * \param p8est     The p8est to be written.
 * \param geom      A p8est_geometry_t structure or NULL for vertex space.
 * \param scale     The relative length factor of the quadrants.
 *                  Use 1.0 to fit quadrants exactly, less to create gaps.
 * \param filename  The first part of the name which will have
 *                  the proc number appended to it (i.e., the
 *                  output file will be filename_rank.vtu).
 *
 * \return          This returns 0 if no error and -1 if there is an error.
 */
p8est_vtk_context_t *p8est_vtk_write_header (p8est_t * p8est,
                                             p8est_geometry_t * geom,
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
 * \param [in,out] cont    A vtk context created by p8est_vtk_write_header.
 * \param [in] write_tree  Boolean to determine if the tree id should be output.
 * \param [in] write_level Boolean to determine if the tree levels should be output.
 * \param [in] write_rank  Boolean to determine if the MPI rank should be output.
 * \param [in] wrap_rank   Number to wrap around the rank with a modulo operation.
 *                         Can be 0 for no wrapping.
 * \param [in] num_cell_scalars Number of cell scalar datasets to output.
 * \param [in] num_cell_vectors Number of cell vector datasets to output.
 *
 * The variable arguments need to be pairs of (fieldname, fieldvalues)
 * where the cell scalar pairs come first, followed by the cell vector pairs.
 *
 * \return          On success, the context that has been passed in.
 *                  On failure, returns NULL and deallocates the context.
 */
p8est_vtk_context_t *p8est_vtk_write_cell_data (p8est_vtk_context_t * cont,
                                                int write_tree,
                                                int write_level,
                                                int write_rank,
                                                int wrap_rank,
                                                int num_cell_scalars,
                                                int num_cell_vectors, ...);

/** This will write custom point data to the vtu file.
 *
 * Writing a VTK file is split into a few routines.
 * This allows there to be an arbitrary number of
 * fields.
 *
 * \param [in,out] cont          A vtk context created by p8est_vtk_write_header.
 * \param [in] num_point_scalars Number of point scalar datasets to output.
 * \param [in] num_point_vectors Number of point vector datasets to output.
 *
 * The variable arguments need to be pairs of (fieldname, fieldvalues)
 * where the point scalar pairs come first, followed by the point vector pairs.
 *
 * \return          On success, the context that has been passed in.
 *                  On failure, returns NULL and deallocates the context.
 */
p8est_vtk_context_t *p8est_vtk_write_point_data (p8est_vtk_context_t * cont,
                                                 int num_point_scalars,
                                                 int num_point_vectors, ...);

/** TODO: Please document this function and
 *        add analogous function for cell data.
 */
p8est_vtk_context_t *p8est_vtk_write_point_datav (p8est_vtk_context_t * cont,
                                                  int num_point_scalars,
                                                  int num_point_vectors,
                                                  va_list ap);

p8est_vtk_context_t *p8est_vtk_write_cell_datav (p8est_vtk_context_t * cont,
                                                 int write_tree,
                                                 int write_level,
                                                 int write_rank,
                                                 int wrap_rank,
                                                 int num_cell_scalars,
                                                 int num_cell_vectors,
                                                 va_list ap);

/** This will write a point scalar field to the vtu file.
 *
 * Writing a VTK file is split into a few routines.
 * This allows there to be an arbitrary number of fields.
 *
 * \param [in,out] cont    A vtk context created by p8est_vtk_write_header.
 * \param [in] scalar_name The name of the scalar field.
 * \param [in] values      The point values that will be written.
 *
 * \return          On success, the context that has been passed in.
 *                  On failure, returns NULL and deallocates the context.
 */
p8est_vtk_context_t *p8est_vtk_write_point_scalar (p8est_vtk_context_t * cont,
                                                   const char *scalar_name,
                                                   const double *values);

/** This will write a cell scalar field to the vtu file.
 *
 * Writing a VTK file is split into a few routines.
 * This allows there to be an arbitrary number of fields.
 *
 * \param [in,out] cont    A vtk context created by p8est_vtk_write_header.
 * \param [in] scalar_name The name of the scalar field.
 * \param [in] values      The cell values that will be written.
 *
 * \return          On success, the context that has been passed in.
 *                  On failure, returns NULL and deallocates the context.
 */
p8est_vtk_context_t *p8est_vtk_write_cell_scalar (p8est_vtk_context_t * cont,
                                                  const char *scalar_name,
                                                  const double *values);

/** This will write a 3-vector point field to the vtu file.
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
p8est_vtk_context_t *p8est_vtk_write_point_vector (p8est_vtk_context_t * cont,
                                                   const char *vector_name,
                                                   const double *values);

/** This will write a 3-vector cell field to the vtu file.
 *
 * Writing a VTK file is split into a few routines.
 * This allows there to be an arbitrary number of fields.
 *
 * \param [in,out] cont    A vtk context created by p8est_vtk_write_header.
 * \param [in] vector_name The name of the vector field.
 * \param [in] values      The cell values that will be written.
 *
 * \return          On success, the context that has been passed in.
 *                  On failure, returns NULL and deallocates the context.
 */
p8est_vtk_context_t *p8est_vtk_write_cell_vector (p8est_vtk_context_t * cont,
                                                  const char *vector_name,
                                                  const double *values);

/** Write the footer of the vtu file and clean up.
 *
 * Writing a VTK file is split into a few routines.
 * This allows there to be an arbitrary number of
 * fields.  To write out two fields the
 * calling sequence would be something like
 *
 * \begincode
 * p8est_vtk_write_header (p8est, ..., "output");
 * p8est_vtk_write_point_data (...);
 * p8est_vtk_write_cell_data (...);
 * ...
 * p8est_vtk_write_footer (vtk_context);
 * \endcode
 *
 * This function writes the footer information to the vtk file and closes all
 * destroys the vtk context.
 *
 * \param [in] cont Context is deallocated before the function returns.
 *
 * \return          This returns 0 if no error and -1 if there is an error.
 */
int                 p8est_vtk_write_footer (p8est_vtk_context_t * cont);

SC_EXTERN_C_END;

#endif /* !P8EST_VTK_H */
