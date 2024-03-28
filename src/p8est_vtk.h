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

/** \file p8est_vtk.h
 *
 * Routines for printing a forest and associated fields to VTK format.
 *
 * \ingroup p8est
 */

#ifndef P8EST_VTK_H
#define P8EST_VTK_H

#include <p8est_geometry.h>
#include <p8est.h>

SC_EXTERN_C_BEGIN;

/** Opaque context type for writing VTK output with multiple function calls.
 */
typedef struct p8est_vtk_context p8est_vtk_context_t;

/** Write the p8est in VTK format.
 *
 * This is a convenience function for the special case of writing out
 * the tree id, quadrant level, and MPI rank of each quadrant as cell data.
 * One file is written per MPI rank, and one meta file on rank 0.
 * The quadrants are scaled to length .95; see \ref p8est_vtk_write_header.
 * This function will abort if there is a file error.
 *
 * \param [in] p8est    The p8est to be written.
 * \param [in] geom     A \ref p8est_geometry_t structure or NULL for vertex space
 *                      as defined by the \a p8est's \ref p8est_connectivity_t member.
 * \param [in] filename The first part of the file name which will have the
 *                      MPI rank appended to it: The output file will be
 *                      filename_rank.vtu, and the meta file filename.pvtu.
 */
void                p8est_vtk_write_file (p8est_t * p8est,
                                          p8est_geometry_t * geom,
                                          const char *filename);

/** The first call to write a VTK file using individual functions.
 *
 * Writing a VTK file is split into multiple functions that keep a context.
 * This is the first function that allocates the opaque context structure.
 * After allocation, further parameters can be set for the context.
 * Then, the header, possible data fields, and the footer must be written.
 * The process can be aborted any time by destroying the context.  In this
 * case, open files are closed cleanly with only partially written content.
 *
 * \param p4est     The p8est to be written.
 *                  If no geometry is specified in
 *                  \ref p8est_vtk_context_set_geom, we require
 *                  \b p8est->connectivity to have valid vertex arrays.
 * \param filename  The first part of the name which will have the processor
 *                  number appended to it (i.e., the output file will be
 *                  filename_rank.vtu).  The parallel meta-files for Paraview
 *                  and Visit use this basename too.
 *                  We copy this filename to internal storage, so it is not
 *                  needed to remain alive after calling this function.
 * \return          A VTK context fur further use.
 */
p8est_vtk_context_t *p8est_vtk_context_new (p8est_t * p4est,
                                            const char *filename);

/** Modify the geometry transformation registered in the context.
 * After \ref p8est_vtk_context_new, it is at the default NULL.
 * \param [in,out] cont         The context is modified.
 *                              It must not yet have been used to start writing
 *                              in \ref p8est_vtk_write_header.
 * \param geom      A \ref p8est_geometry_t structure, or NULL for vertex space.
 *                  If NULL, \b p8est->connectivity->vertices and
 *                  \b tree_to_vertex must be non-NULL.
 */
void                p8est_vtk_context_set_geom (p8est_vtk_context_t * cont,
                                                p8est_geometry_t * geom);

/** Modify the context parameter for scaling the quadrants.
 * A scale < 1 places a visual gap between adjacent quadrants.
 * After \ref p8est_vtk_context_new, it is at the default 0.95.
 * \param [in,out] cont         The context is modified.
 *                              It must not yet have been used to start writing
 *                              in \ref p8est_vtk_write_header.
 * \param [in] scale            Scale parameter must be in (0, 1].
 */
void                p8est_vtk_context_set_scale (p8est_vtk_context_t * cont,
                                                 double scale);

/** Modify the context parameter for expecting continuous point data.
 * If set to true, the point data is understood as a continuous field.
 * In this case, we can significantly reduce the file size when scale == 1.
 * For discontinuous point data, it should be set to false.
 * After \ref p8est_vtk_context_new, it is at the default false.
 * \param [in,out] cont         The context is modified.
 *                              It must not yet have been used to start writing
 *                              in \ref p8est_vtk_write_header.
 * \param [in] continuous       Boolean parameter.
 */
void                p8est_vtk_context_set_continuous (p8est_vtk_context_t *
                                                      cont, int continuous);
/** Cleanly destroy a \ref p8est_vtk_context_t structure.
 *
 * This function closes all the file pointers and frees the context.
 * Tt can be called even if the VTK output
 * has only been partially written, the files' content will be incomplete.
 *
 * \param[in] context     The VTK file context to be destroyed.
 */
void                p8est_vtk_context_destroy (p8est_vtk_context_t * context);

/** Write the VTK header.
 *
 * Writing a VTK file is split into a few routines.
 * This allows there to be an arbitrary number of
 * fields.  The calling sequence would be something like
 *
 *     vtk_context = p8est_vtk_context_new (p8est, "output");
 *     p8est_vtk_context_set_* (vtk_context, parameter);
 *     vtk_context = p8est_vtk_write_header (vtk_context, ...);
 *     if (vtk_context == NULL) { error; }
 *     vtk_context = p8est_vtk_write_cell_data (vtk_context, ...);
 *     if (vtk_context == NULL) { error; }
 *     vtk_context = p8est_vtk_write_point_data (vtk_context, ...);
 *     if (vtk_context == NULL) { error; }
 *     retval = p8est_vtk_write_footer (vtk_context);
 *     if (retval) { error; }
 *
 * Each of these functions opens and closes files as necessary.
 * Generally, each output file is written in subsequent chunks.
 *
 * This function writes point positions for the quadrants' vertices
 * and the maps of elements to types and vertices as required by VTK.
 *
 * \param [in,out] cont    A VTK context created by \ref p8est_vtk_context_new.
 *                         None of the vtk_write functions must have been called.
 *                         This context is the return value if no error occurs.
 *
 * \return          On success, an opaque context (p8est_vtk_context_t) pointer
 *                  that must be passed to subsequent p8est_vtk calls.  It is
 *                  required to call \ref p8est_vtk_write_footer eventually with
 *                  this value.  Returns NULL on error.
 */
p8est_vtk_context_t *p8est_vtk_write_header (p8est_vtk_context_t * cont);

/** Write the VTK header for higher order visualization.
 *
 * This function follows the same routines as p8est_vtk_write_header.
 * In addition, the caller must pass in an array containing coordinates for
 * each point, as well as an integer representing the number of points in
 * one direction each element has (for example, in an 8x8x8 cell, pass in 8).
 *
 * \param [in,out] cont    A VTK context created by \ref p8est_vtk_context_new.
 *                         None of the vtk_write functions must have been called.
 *                         This context is the return value if no error occurs.
 * \param [in] positions   An sc_array_t of doubles containing the coordinates
 *                         of all points to be written. Ordering of data is
 *                         [ x_0, y_0, (z_0) ... x_n, y_n, (z_n) ]
 * \param [in] Nnodes1D    Integer number of points in each element in 1D.
 *
 * \return          On success, an opaque context (p8est_vtk_context_t) pointer
 *                  that must be passed to subsequent p8est_vtk calls.  It is
 *                  required to call \ref p8est_vtk_write_footer eventually with
 *                  this value.  Returns NULL on error.
 */
p8est_vtk_context_t *p8est_vtk_write_header_ho (p8est_vtk_context_t * cont,
                                                sc_array_t * positions,
                                                int Nnodes1D);

 /** Write VTK cell data.
 *
 * There are options to have this function write
 * the tree id, quadrant level, or MPI rank without explicit input data.
 *
 * Writing a VTK file is split into a few routines.
 * This allows there to be an arbitrary number of fields.
 *
 * \param [in,out] cont    A VTK context created by \ref p8est_vtk_context_new.
 * \param [in] write_tree  Boolean to determine if the tree id should be output.
 * \param [in] write_level Boolean to determine if the tree levels should be output.
 * \param [in] write_rank  Boolean to determine if the MPI rank should be output.
 * \param [in] wrap_rank   Number to wrap around the rank with a modulo operation.
 *                         Can be 0 for no wrapping.
 * \param [in] num_cell_scalars Number of cell scalar datasets to output.
 * \param [in] num_cell_vectors Number of cell vector datasets to output.
 *
 * \param [in] fieldnames  Array of char strings containing the name
 *                         of each data field.
 * \param [in] values      Array of pointers to sc_array_t holding variables
 *                         of type double, one value (scalar data) or three
 *                         values (vector data) for each local quadrant.
 *
 * The number of doubles in each sc_array must be exactly \a
 * p8est->local_num_quadrants for scalar data and \a
 * 3*p8est->local_num_quadrants for vector data.
 * The cell scalar data come first, followed by the cell vector data.
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
                                                int num_cell_vectors,
                                                const char *fieldnames[],
                                                sc_array_t * values[]);

/** Write VTK cell data.
 *
 * There are options to have this function write
 * the tree id, quadrant level, or MPI rank without explicit input data.
 *
 * Writing a VTK file is split into a few routines.
 * This allows there to be an arbitrary number of fields.
 *
 * \param [in,out] cont    A VTK context created by \ref p8est_vtk_context_new.
 * \param [in] write_tree  Boolean to determine if the tree id should be output.
 * \param [in] write_level Boolean to determine if the tree levels should be output.
 * \param [in] write_rank  Boolean to determine if the MPI rank should be output.
 * \param [in] wrap_rank   Number to wrap around the rank with a modulo operation.
 *                         Can be 0 for no wrapping.
 * \param [in] num_cell_scalars Number of cell scalar datasets to output.
 * \param [in] num_cell_vectors Number of cell vector datasets to output.
 *
 * The variable arguments need to be pairs of (fieldname, fieldvalues), followed
 * by a final argument of the VTK context cont (same as the first argument).
 * The cell scalar pairs come first, followed by the cell vector pairs, then cont.
 * Each 'fieldname' argument shall be a char string containing the name of the data
 * contained in the following 'fieldvalues'.  Each of the 'fieldvalues'
 * arguments shall be an sc_array_t * holding double variables.  The number of
 * doubles in each sc_array must be exactly \a p8est->local_num_quadrants for
 * scalar data and \a 3*p8est->local_num_quadrants for vector data.
 *
 * \note The current p8est_vtk_context_t structure, \a cont, must be the first
 * and the last argument of any call to this function; this argument is used to
 * verify that the correct number of variable arguments have been provided.
 *
 * \return          On success, the context that has been passed in.
 *                  On failure, returns NULL and deallocates the context.
 */
p8est_vtk_context_t *p8est_vtk_write_cell_dataf (p8est_vtk_context_t * cont,
                                                 int write_tree,
                                                 int write_level,
                                                 int write_rank,
                                                 int wrap_rank,
                                                 int num_cell_scalars,
                                                 int num_cell_vectors, ...);

/** Write VTK cell data.
 *
 * This function exports custom cell data to the vtk file; it is functionally
 * the same as \b p8est_vtk_write_cell_dataf with the only difference being
 * that instead of a variable argument list, an initialized \a va_list is
 * passed as the last argument. That means \a va_start has already been called.
 * The \a va_list is initialized from the variable argument list of the calling
 * function. Elements of va_list are processed as "pairs" of (fieldname, fieldvalues).
 * That means <va_list[0], va_list[1]> represents one pair, <va_list[2],
 * va_list[3]> next one and so on.
 * Each 'fieldname' shall be a char string containing the name of the data
 * contained in the following 'fieldvalues'. Each of the 'fieldvalues'
 * shall be an sc_array_t * holding double variables.
 * The cell scalar pairs come first, followed by the cell vector pairs, followed
 * by VTK context \a cont (same as the first argument).
 * The number of * doubles in each sc_array must be exactly \a
 * p8est->local_num_quadrants for scalar data and \a
 * 3*p8est->local_num_quadrants for vector data.
 *
 * \param [in,out] cont    A VTK context created by \ref p8est_vtk_context_new.
 * \param [in] write_tree  Boolean to determine if the tree id should be output.
 * \param [in] write_level Boolean to determine if the tree levels should be output.
 * \param [in] write_rank  Boolean to determine if the MPI rank should be output.
 * \param [in] wrap_rank   Number to wrap around the rank with a modulo operation.
 *                         Can be 0 for no wrapping.
 * \param [in] num_cell_scalars Number of cell scalar datasets to output.
 * \param [in] num_cell_vectors Number of cell vector datasets to output.
 * \param [in,out] ap      An initialized va_list used to access the
 *                         scalar/vector data.
 *
 * \return          On success, the context that has been passed in.
 *                  On failure, returns NULL and deallocates the context.
 */

p8est_vtk_context_t *p8est_vtk_write_cell_datav (p8est_vtk_context_t * cont,
                                                 int write_tree,
                                                 int write_level,
                                                 int write_rank,
                                                 int wrap_rank,
                                                 int num_cell_scalars,
                                                 int num_cell_vectors,
                                                 va_list ap);

/** Write VTK point data.
 *
 * Writing a VTK file is split into a few routines.
 * This allows there to be an arbitrary number of
 * fields.
 *
 * \param [in,out] cont    A VTK context created by \ref p8est_vtk_context_new.
 * \param [in] num_point_scalars Number of point scalar datasets to output.
 * \param [in] num_point_vectors Number of point vector datasets to output.
 *
 * The variable arguments need to be pairs of (fieldname, fieldvalues) where
 * the point scalar pairs come first, followed by the point vector pairs.  Each
 * 'fieldname' argument shall be a char string containing the name of the data
 * contained in the following 'fieldvalues'. Each of the 'fieldvalues'
 * arguments shall be an sc_array_t * holding double variables. The number of
 * doubles in each sc_array must be exactly the number of components (1 for
 * scalar and 3 for vector) times 8 times number of elements.
 *
 * \note The current p8est_vtk_context_t structure, cont, must be the last
 * argument of any call to this function; this argument is used to verify
 * that the correct number of variable arguments have been provided.
 *
 * \note The number of point scalar data in each
 * sc_array must be exactly \a P8EST_CHILDREN*local_num_quadrants, and the
 * number of point vector data must be exactly \a
 * 3*P8EST_CHILDREN*local_num_quadrants. I.e. there must be data for every
 * corner of every quadrant in the \a p8est, even if the corner is shared by
 * multiple quadrants.
 *
 * \return          On success, the context that has been passed in.
 *                  On failure, returns NULL and deallocates the context.
 */
p8est_vtk_context_t *p8est_vtk_write_point_dataf (p8est_vtk_context_t * cont,
                                                  int num_point_scalars,
                                                  int num_point_vectors, ...);

/** Write the VTU footer and clean up.
 *
 * Writing a VTK file is split into a few routines.
 * This function writes the footer information to the VTK file and cleanly
 * destroys the VTK context.
 *
 * \param [in] cont Context is deallocated before the function returns.
 *
 * \return          This returns 0 if no error and -1 if there is an error.
 */
int                 p8est_vtk_write_footer (p8est_vtk_context_t * cont);

SC_EXTERN_C_END;

#endif /* !P8EST_VTK_H */
