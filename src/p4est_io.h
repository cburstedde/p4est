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

/** \file p4est_io.h
 *
 * Provide functions to serialize/deserialize a forest.
 * Some are used as building blocks for \ref p4est_load and \ref p4est_save.
 * Others allow for saving and loading user-defined data to a parallel file.
 */

#ifndef P4EST_IO_H
#define P4EST_IO_H

#include <p4est.h>

SC_EXTERN_C_BEGIN;

/** Extract processor local quadrants' x y level data.
 * Optionally extracts the quadrant data as well into a separate array.
 * \param [in] p4est    The forest is not modified.
 * \param [in,out] data If not NULL, pointer to a pointer that will be set
 *                      to a newly allocated array with per-quadrant data.
 *                      Must be NULL if p4est->data_size == 0.
 * \return              An array of type p4est_qcoord_t that contains
 *                      x y level for each quadrant on this processor.
 *                      The tree information is not extracted.
 */
sc_array_t         *p4est_deflate_quadrants (p4est_t * p4est,
                                             sc_array_t ** data);

/** Create a new p4est based on serialized data.
 * Its revision counter is set to zero.
 * See p4est.h and p4est_communication.h for more information on parameters.
 * \param [in] mpicomm       A valid MPI communicator.
 * \param [in] connectivity  This is the connectivity information that
 *                           the forest is built with.  Note that p4est
 *                           does not take ownership of the memory.
 * \param [in] global_first_quadrant First global quadrant on each proc and
 *                           one beyond.  Copied into global_first_quadrant.
 *                           Local count on rank is gfq[rank + 1] - gfq[rank].
 * \param [in] pertree       The cumulative quadrant counts per tree.
 * \param [in] quadrants     Array as returned by p4est_deflate_quadrants.
 * \param [in] data          Array as from p4est_deflate_quadrants or NULL.
 *                           The elem_size of this array informs data_size.
 *                           Its elem_count equals the number of local quads.
 * \param [in] user_pointer  Assign to the user_pointer member of the p4est.
 * \return              The newly created p4est with a zero revision counter.
 */
p4est_t            *p4est_inflate (sc_MPI_Comm mpicomm,
                                   p4est_connectivity_t * connectivity,
                                   const p4est_gloidx_t *
                                   global_first_quadrant,
                                   const p4est_gloidx_t * pertree,
                                   sc_array_t * quadrants, sc_array_t * data,
                                   void *user_pointer);

/** Opaque context used for writing a p4est data file. */
typedef struct p4est_file_context p4est_file_context_t;

/** Context for the quadrant-wise callback function. */
typedef struct callback_context /* TODO: Write a get callback_ct function */
{
  p4est_quadrant_t   *quad;
  void               *user;
} callback_context_t;

/** Callback to write a user-defined output data into an internal buffer.
 * We expect exactly \a data_size bytes to be copied into \a buffer.
 */
typedef void        (*p4est_file_write_data_t)
                    (size_t data_size, char *buffer, void *user);

/** Callback to read from a file into a user-defined input buffer.
 * The buffer must exist and be at least of length \a data_size.
 */
typedef void        (*p4est_file_read_data_t)
                    (size_t data_size, char *buffer, void *user);

/** Begin saving forest header and per-quadrant data into a parallel file.
 *
 * This function creates a new file or overwrites an existing one.
 * It is collective and creates the file on a parallel file system.
 * It takes an (optional) callback to write a header of given size.
 * This function leaves the file open.  It is necessary to call \ref
 * p4est_file_close (possibly after writing one or more data sets).
 * The file is opened in a write-only mode.
 *
 * We add some basic metadata to the file.
 * The file written contains the header data and data sets
 * exactly as specified by the open/write functions called.
 * The header consists of the metadata header specified by p4est
 * followed by a user-defined header. 
 *
 * It is the application's responsibility to write sufficient header
 * information to determine the number and size of the data sets
 * if such information is not recorded and maintained externally.
 *
 * This function aborts on I/O and MPI errors.
 *
 * \param [in] p4est          Valid forest.
 * \param [in] filename       Path to parallel file that is to be created.
 * \param [in] header_size    This number of bytes is written at the start
 *                            of the file on rank zero.  May be 0.
 * \param [in] quadrant_data TODO A pointer to a array of header_size many
 *                            bytes. The data is written to the file as a
 *                            header. 
 *                            For header_size == 0
 *                            the function does not write a user-header. 
 *                            May be NULL if header_size == 0. 
 *                            Must not be NULL on rank zero when
 *                            \a header_size is greater zero.
 * \return                    Newly allocated context to continue writing
 *                            and eventually closing the file.
 */
p4est_file_context_t *p4est_file_open_create
  (p4est_t * p4est, const char *filename,
   size_t header_size, void *header_data);

/** Similar to \ref p4est_file_open_create except the header exists.
 * The file specified must exist and is opened.  Its header is preserved.
 */
p4est_file_context_t *p4est_file_open_append
  (p4est_t * p4est, const char *filename, size_t header_size);

/** Open a file for reading and read its header on rank zero.
 * The header data is broadcast to all ranks after reading.
 * The file must exist and be at least of the size of the header.
 * In practice, the header size should match the one writing the file.
 *
 * \parma [in] p4est        The forest must be of the same refinement
 *                          pattern as the one used for writing the file.
 *                          Its global number of quadrants must match.
 *                          It is possible, however, to use a different
 *                          partition or number of ranks from writing it.
 * \param [in] hcall        Callback executed on all ranks.
 *                          It supplies the header data just read.
 */
p4est_file_context_t *p4est_file_open_read (p4est_t * p4est,
                                            const char *filename,
                                            size_t header_size,
                                            void *header_data);

/** Write one (more) per-quadrant data set to a parallel output file.
 *
 * This function requires an opened file context.
 * The data set is appended to the header/previously written data sets.
 * This function writes a block of the size number of quadrants * data_size.
 *
 * This function aborts on I/O and MPI errors.
 *
 * \param [in,out] fc         Context previously created by \ref
 *                            p4est_file_open_create or \ref
 *                            p4est_file_open_append.
 * \param [in] quadrant_data  An array of the length number of local quadrants
 *                            with the element size equal to number of bytes
 *                            written per quadrant. The quadrant data is expected
 *                            to be stored according to the Morton order of
 *                            the quadrants. For quadrant_data->elem_size == 0
 *                            the function does nothing and returns the unchanged
 *                            file context.
 * \return                    Return the input context to continue writing
 *                            and eventually closing the file.
 */
void                p4est_file_write (p4est_file_context_t * fc,
                                      sc_array_t * quadrant_data);

/** Read one (more) per-quadrant data set from a parallel input file.
 * This function requires the appropriate number of readable bytes.
 * In practice, the data size to read should match the size written.
 * This function aborts if the number of bytes to read is bigger than the
 * datatset that corresponds to the processor.
 * The data size to read is encoded by the element size of quadrant_data
 * It is possible to skip over a data set to read by a NULL \ref sc_array.
 * It is legal to close a file before all data sets have been read.
 *
 * \param [in,out] fc         Context previously created by \ref
 *                            p4est_file_open_read.  It keeps track
 *                            of the data sets read one after another.
 * \param [in,out] quadrant_data  An array of the length number of local quadrants
 *                            with the element size equal to number of bytes
 *                            read per quadrant. The quadrant data is read
 *                            according to the Morton order of the quadrants.
 *                            For quadrant_data->elem_size == 0
 *                            the function does nothing and returns the unchanged
 *                            file context. The same holds for
 *                            quadrant_data == NULL.
 */
void                p4est_file_read (p4est_file_context_t * fc,
                                     sc_array_t * quadrant_data);

/** Read the metadata information.
 * The file must be opened by \ref p4est_file_open_read.
 * \param [in]  fc                  File context previously created by \ref
 *                                  p4est_file_open_read.
 * \param [in,out] global_num_quad  After the function call this variable will
 *                                  hold the number of global quadrants of the
 *                                  mesh that was used to create the given file.
 * \param [in,out] p4est_version    The p4est version that was used to create
 *                                  this file.
 * \param [in,out] file_io_rev      After the function call this variable will
 *                                  hold the revision number of file io data
 *                                  structure of p4est that was used to create
 *                                  the given file.
 * \param [in,out] magic_num        After the function call this variable will
 *                                  hold the magic number of the file.
 * \param [in,out] elem_size        After the function call this variable will
 *                                  hold an array of the lentgh number of
 *                                  arrays in the given file and the associated
 *                                  values of the array are the number of bytes
 *                                  per quadrant-wise data. TODO: nicht expotieren
 */
void                p4est_file_info (p4est_file_context_t * fc,
                                     p4est_gloidx_t * global_num_quads,
                                     char p4est_version[16], int *file_io_rev,
                                     int *magic_num, sc_array_t * elem_size);

/** Close a file opened for parallel write/read and free the context.
 * \param [in,out] fc       Context previously created by \ref
 *                          p4est_file_open_create, \ref
 *                          p4est_file_open_append, or \ref
 *                          p4est_file_open_read.  Is freed.
 */
void                p4est_file_close (p4est_file_context_t * fc);

SC_EXTERN_C_END;

#endif /* !P4EST_IO_H */
