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
 *
 * Furthermore, this module provides functions to write and read general data
 * files associated with a p4est.
 */

#ifndef P4EST_IO_H
#define P4EST_IO_H

#include <p4est.h>

SC_EXTERN_C_BEGIN;

#define P4EST_DATA_FILE_EXT "p4d" /**< file extension of p4est data files */
#define P4EST_MAGIC_NUMBER "p4data0" /**< magic string for p4est data files */
#define P4EST_NUM_METADATA_BYTES 96 /**< number of file metadata bytes */
#define P4EST_NUM_MAGIC_BYTES 8 /**< number of bytes of the magic number */
#define P4EST_NUM_VERSION_STR_BYTES 24 /**< number of bytes of the version string*/
#define P4EST_NUM_ARRAY_METADATA_BYTES 14 /**< number of array metadata bytes */
/* subtract 2 for '\n' at the beginning and end of the array metadata */
#define P4EST_NUM_ARRAY_METADATA_CHARS (P4EST_NUM_ARRAY_METADATA_BYTES - 2) /**< number of array metadata chars */
#define P4EST_BYTE_DIV 16 /**< All data blocks are padded to be divisible by this. */
#define P4EST_MAX_NUM_PAD_BYTES (P4EST_BYTE_DIV + 1) /**< We enforce to pad in any
                                                          case and the padding string
                                                          needs to contain two
                                                          newline characters and
                                                          therefore this is the
                                                          maximal number of pad
                                                          bytes. */
#define P4EST_NUM_USER_STRING_BYTES 48 /**< number of user string bytes */
#define P4EST_NUM_FIELD_HEADER_BYTES (2 + P4EST_NUM_ARRAY_METADATA_BYTES + P4EST_NUM_USER_STRING_BYTES)
                                     /**< number of bytes of one field header */
#define P4EST_FILE_MAX_GLOBAL_QUAD 9999999999999999 /**< maximal number of global quadrants */
#define P4EST_FILE_MAX_BLOCK_SIZE 9999999999999 /**< maximal number of block bytes */
#define P4EST_FILE_MAX_FIELD_ENTRY_SIZE 9999999999999 /**< maximal numeber of bytes per field entry */

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

/** p4est data file format
 * All p4est data files have 64 bytes file header section at the beginning of the file.
 * The file header section is written to the file as string without NUL-termination
 * (called string*) and is therefore readable in a text editor.
 *
 * File Header (96 bytes):
 * 7 bytes magic number (p4data0) and 1 byte new line char.
 * 23 bytes p4est version string* and 1 byte new line char.
 * 47 bytes user string*  and 1 byte new line char.
 * 16 bytes number of global quadrants.
 *
 * The file header section is padded by 16 bytes consisting of 1 byte
 * new line char succeeded by 14 bytes of spaces and 1 trailing byte
 * new line char.
 *
 * The actual data is stored in arrays corresponding to a mesh of a p4est
 * or in header sections that have a fixed user-defined size. The header
 * sections are written and read on rank 0.
 * One data field stores a fixed number of bytes of user-
 * defined data per quadrant of a certain p4est. Therefore, one user-defined
 * data field is of the size p4est->global_num_quadrants * data_size, where
 * data_size is set by the user. The file format is partition independent.
 * The data fields are padded such that the number of bytes for
 * an array is divisible by 16. The padding also enforced for data blocks
 * that have a size that is divisble by 16.
 * The p4est data file consists of a variable number (including 0) of
 * these two types of data sections.
 * Every data section includes 64 bytes of section header written at the beginning
 * by p4est. These 64 bytes are again written to the file as string* and can
 * be read using a text editor.
 *
 * Data section Header (64 bytes):
 * One byte data section type specific character (H for a header section and F for
 * a data field), 1 byte space and 13 bytes size in number of bytes for a
 * header section and data size per element in byte for a field section
 * and one trailing byte new line char.
 * 47 bytes user-defined string* and 1 byte new line char.
 *
 * The structure of p4est and p8est data files differs only by the magic number.
 *
 * The p4est metadata of a p4est data file can be accessed by \ref p4est_file_info().
 */

/** Opaque context used for writing a p4est data file. */
typedef struct p4est_file_context p4est_file_context_t;

/** Error values for p4est_file functions.
 */
typedef enum p4est_file_error
{
  P4EST_FILE_ERR_SUCCESS = sc_MPI_ERR_LASTCODE, /**< file function completed with success */
  /* MPI I/O error classes or its replacement without MPI I/O */
  P4EST_FILE_ERR_FILE, /**< invalid file handle */
  P4EST_FILE_ERR_NOT_SAME, /**< collective arg not identical */
  P4EST_FILE_ERR_AMODE, /**< access mode error */
  P4EST_FILE_ERR_NO_SUCH_FILE, /**< file does not exist */
  P4EST_FILE_ERR_FILE_EXIST, /**< file exists already */
  P4EST_FILE_ERR_BAD_FILE, /**< invaild file name */
  P4EST_FILE_ERR_ACCESS, /**< permission denied */
  P4EST_FILE_ERR_NO_SPACE, /**< not enough space */
  P4EST_FILE_ERR_QUOTA, /**< quota exceeded */
  P4EST_FILE_ERR_READ_ONLY, /**< read only file (system) */
  P4EST_FILE_ERR_IN_USE, /**< file currently open by other process */
  P4EST_FILE_ERR_IO, /**< other I/O error */
  /* the following error codes are only defined in p4est */
  P4EST_FILE_ERR_FORMAT,  /**< read file has a wrong format */
  P4EST_FILE_ERR_IN_DATA, /**< input data of file function is invalid */
  P4EST_FILE_ERR_COUNT,   /**< read or write count error that was not
                                 classified as a format error */
  P4EST_FILE_ERR_UNKNOWN, /**< unknown error */
  P4EST_FILE_ERR_LASTCODE /**< to define own error codes for
                                  a higher level application
                                  that is using p4est_file
                                  functions */
}
p4est_file_error_t;

/** Begin writing file header and saving data blocks into a parallel file.
 *
 * This function creates a new file or overwrites an existing one.
 * It is collective and creates the file on a parallel file system.
 * It takes an (optional) pointer to write a header of given size.
 * This function leaves the file open if MPI I/O is available.
 * It is necessary to call \ref
 * p4est_file_close (possibly after writing one or more data sets).
 * The file is opened in a write-only mode.
 *
 * We add some basic metadata to the file.
 * The file written contains the file header and data sets
 * as specified by the open/write functions called.
 * The file header consists of the metadata specified by p4est.
 *
 * The number of global quadrants must be less or equal
 * \ref P4EST_FILE_MAX_GLOBAL_QUAD.
 *
 * It is the application's responsibility to write sufficient header
 * information (cf. \ref p4est_file_write_header) to determine the number and
 * size of the data sets if such information is not recorded and maintained
 * externally.
 * However, p4est makes some metadata accessible via
 * \ref p4est_file_info.
 *
 * This function does not abort on MPI I/O errors but returns NULL.
 *
 * \param [in] p4est          Valid forest.
 * \param [in] filename       Path to parallel file that is to be created.
 * \param [in] user_string    A user string that is written to the file header.
 *                            Only \ref P4EST_NUM_USER_STRING_BYTES
 *                            bytes without NUL-termination are
 *                            written to the file. If the user gives less
 *                            bytes the user_string in the file header is padded
 *                            by spaces.
 * \param [out] errcode       An errcode that can be interpreted by
 *                            \ref p4est_file_error_string.
 * \return                    Newly allocated context to continue writing
 *                            and eventually closing the file. NULL in
 *                            case of error.
 */
p4est_file_context_t *p4est_file_open_create
  (p4est_t * p4est, const char *filename,
   const char user_string[P4EST_NUM_USER_STRING_BYTES], int *errcode);

/** Open a file for reading and read its user string on rank zero.
 * The user string is broadcasted to all ranks after reading.
 * The file must exist and be at least of the size of the file header.
 *
 * If the file has wrong metadata the function reports the error using
 * /ref P4EST_LERRORF, collectively close the file and deallocate
 * the file context. In this case the function returns NULL on all ranks.
 * The wrong file format or a wrong file header causes \ref P4EST_FILE_ERR_IO
 * as errcode.
 *
 * This function does not abort on MPI I/O errors but returns NULL.
 *
 * \param [in] p4est            The forest must be of the same refinement
 *                              pattern as the one used for writing the file.
 *                              Its global number of quadrants must match.
 *                              It is possible, however, to use a different
 *                              partition or number of ranks from writing it.
 * \param [in] filename         The path to the file that is opened.
 * \param [in,out] user_string  At least \ref P4EST_NUM_USER_STRING_BYTES
 *                              bytes. The user string is written
 *                              to the passed array including padding spaces
 *                              and a trailing NUL-termination.
 * \param [out] errcode         An errcode that can be interpreted by \ref
 *                              p4est_file_error_string.
 * \return                      Newly allocated context to continue reading
 *                              and eventually closing the file. NULL in
 *                              case of error.
 */
p4est_file_context_t *p4est_file_open_read (p4est_t * p4est,
                                            const char *filename,
                                            char
                                            user_string
                                            [P4EST_NUM_USER_STRING_BYTES],
                                            int *errcode);

/** Write a header block to an opened file.
 * This function requires an opened file context.
 * The header data and its metadata are written on rank 0.
 * The number of header bytes must be less or equal
 * \ref P4EST_FILE_MAX_BLOCK_SIZE.
 *
 * \param [out] fc            Context previously created by \ref
 *                            p4est_file_open_create.
 * \param [in]  header_size   The size of header_data in bytes.
 *                            This function returns the passed fc
 *                            parameter and sets errcode to
 *                            sc_MPI_SUCCESS if it is called
 *                            for header_size == 0. Nothing is
 *                            written to the file and fc stays
 *                            untouched.
 * \param [in]  header_data   A pointer to the header data. The user is
 *                            responsible for the validality of the header
 *                            data.
 * \param [in]  user_string   Maximal \ref P4EST_NUM_USER_STRING_BYTES bytes.
 *                            These chars are written to the block
 *                            header and padded to 
 *                            \ref P4EST_NUM_USER_STRING_BYTES - 1 chars
 *                            by adding spaces. The '\0' is not written
 *                            to the file.
 * \param [out] errcode       An errcode that can be interpreted by \ref
 *                            p4est_file_error_string.
 * \return                    Return the input context to continue writing
 *                            and eventually closing the file. The return
 *                            value is NULL in case of error, then
 *                            it also holds errcode != 0 and the file is
 *                            tried to close and fc is freed.
 */
p4est_file_context_t *p4est_file_write_header (p4est_file_context_t * fc,
                                               size_t header_size,
                                               const void *header_data,
                                               const char
                                               user_string
                                               [P4EST_NUM_USER_STRING_BYTES],
                                               int *errcode);

/** Read a header block from an opened file.
 * This function requires an opened file context.
 * The header data is read on rank 0.
 *
 * If the user does not have the header_size to call this function, the user
 * can user \ref p4est_file_info to obtain the required information.
 *
 * The passed header_size is compared to the header_size stored in the file.
 * If the values do not equal each other, the function reports details via
 * /ref P4EST_LERRORF and closes and deallocate the file context. The return
 * value in this case is NULL.
 * If the block header information is not matching the passed parameters
 * the function sets \ref P4EST_FILE_ERR_IO for errcode.
 *
 * \param [out] fc              Context previously created by \ref
 *                              p4est_file_open_create.
 * \param [in]  header_size     The size of the header that is read.
 * \param [in, out] header_data header_size allocated bytes. This data will be
 *                              filled with the header data from file. If this
 *                              is NULL it means that the current header block
 *                              is skipped and the internal file pointer of the
 *                              file context is set to the next data block. If
 *                              current data block is not a header block, the
 *                              file is closed and the file context is
 *                              deallocated. Furthermore, in this case the
 *                              function returns NULL and sets errcode to
 *                              \ref P4EST_FILE_ERR_IO.
 * \param [in,out] user_string  At least \ref P4EST_NUM_USER_STRING_BYTES bytes.
 *                              Filled by the padded user string and
 *                              a trailing NUL-termination char.
 * \param [out] errcode         An errcode that can be interpreted by \ref
 *                              p4est_file_error_string.
 * \return                      Return the input context to continue reading
 *                              and eventually closing the file. The return value
 *                              is NULL if the function was called for
 *                              header_size == 0. The return
 *                              value is also NULL in case of error but then
 *                              it also holds errcode != 0 and the file is
 *                              tried to close and fc is freed.
 */
p4est_file_context_t *p4est_file_read_header (p4est_file_context_t * fc,
                                              size_t header_size,
                                              void *header_data,
                                              char
                                              user_string
                                              [P4EST_NUM_USER_STRING_BYTES],
                                              int *errcode);

/** Write one (more) per-quadrant data set to a parallel output file.
 *
 * This function requires an opened file context.
 * The data set is appended to the header/previously written data sets.
 * This function writes a block of the size number of quadrants * data_size.
 *
 * This function does not abort on MPI I/O errors but returns NULL.
 *
 * The number of bytes per field entry must be less or equal
 * \ref P4EST_FILE_MAX_FIELD_ENTRY_SIZE.
 *
 * \param [out] fc            Context previously created by \ref
 *                            p4est_file_open_create.
 * \param [in] quadrant_data  An array of the length number of local quadrants
 *                            with the element size equal to number of bytes
 *                            written per quadrant. The quadrant data is expected
 *                            to be stored according to the Morton order of
 *                            the quadrants. For quadrant_data->elem_size == 0
 *                            the function does nothing and returns the unchanged
 *                            file context. In this case errcode is set
 *                            to sc_MPI_SUCCESS.
 * \param [in] user_string    An array of maximal \ref
 *                            P4EST_NUM_USER_STRING_BYTES bytes that
 *                            is written without the NUL-termination
 *                            after the array-dependent metadata and before
 *                            the actual data. If the array is shorter the
 *                            written char array will be padded to the
 *                            right by spaces. The user_string is
 *                            written on rank 0 and therefore also only
 *                            required on rank 0. Can be NULL for other
 *                            ranks.
 * \param [out] errcode       An errcode that can be interpreted by \ref
 *                            p4est_file_error_string.
 * \return                    Return the input context to continue writing
 *                            and eventually closing the file. The return value
 *                            is NULL if the function was called for
 *                            quadrant_data->elem_size == 0. The return
 *                            value is also NULL in case of error but then
 *                            it also holds errcode != 0 and the file is
 *                            tried to close and fc is freed.
 */
p4est_file_context_t *p4est_file_write_field (p4est_file_context_t * fc,
                                              sc_array_t * quadrant_data,
                                              const char
                                              user_string
                                              [P4EST_NUM_USER_STRING_BYTES],
                                              int *errcode);

/** Read one (more) per-quadrant data set from a parallel input file.
 * This function requires an opened file context.
 * This function requires the appropriate number of readable bytes.
 * In practice, the data size to read should match the size written.
 * This function reports an error if the number of bytes to read is
 * bigger than the dataset that corresponds to the processor.
 * The data size to read is encoded by the element size of quadrant_data
 * It is possible to skip over a data set to read by a NULL \ref sc_array.
 * It is legal to close a file before all data sets have been read.
 *
 * The function closes and deallocates the file context and returns NULL
 * if the bytes the user wants to read exceed the given file and/or
 * the element size of the array given by quadrant_data->elem_size does not
 * coincide with the element size according to the array metadata given in
 * the file.
 *
 * If the block header information is not matching the passed parameters
 * the function sets \ref P4EST_FILE_ERR_IO for errcode.
 *
 * This function does not abort on MPI I/O errors but returns NULL.
 *
 * \param [in,out] fc         Context previously created by \ref
 *                            p4est_file_open_read (_ext).  It keeps track
 *                            of the data sets read one after another.
 * \param [in,out] quadrant_data  An array of the length number of local quadrants
 *                            with the element size equal to number of bytes
 *                            read per quadrant. The quadrant data is read
 *                            according to the Morton order of the quadrants.
 *                            For quadrant_data->elem_size == 0
 *                            the function does nothing and returns the unchanged
 *                            file context. For quadrant_data == NULL the
 *                            function skips one data array in the file.
 *                            If fc was opened by \ref p4est_file_open_read_ext
 *                            and fc->global_first_quadrant was not set by the
 *                            user, the function uses a uniform partition to read
 *                            the data field in parallel.
 *                            quadrant_data is resized by \ref sc_array_resize.
 * \param [in,out]  user_string At least \ref P4EST_NUM_USER_STRING_BYTES bytes.
 *                            The user string is read on rank 0 and internally
 *                            broadcasted to all ranks.
 * \param [out] errcode       An errcode that can be interpreted by \ref
 *                            p4est_file_error_string.
 * \return                    Return a pointer to input context or NULL in case
 *                            of errors that does not abort the program or if
 *                            the function was called with quadrant_data == NULL.
 *                            In case of error the file is tried to close
 *                            and fc is freed.
 */
p4est_file_context_t *p4est_file_read_field (p4est_file_context_t * fc,
                                             sc_array_t * quadrant_data,
                                             char
                                             user_string
                                             [P4EST_NUM_USER_STRING_BYTES],
                                             int *errcode);

/** A data type that encodes the metadata of one data block in a p4est data file.
 */
typedef struct p4est_file_section_metadata
{
  char                block_type; /**< 'H' (header) or 'F' (data file) */
  size_t              data_size;  /**< data size in bytes per array element ('F')
                                       or of the header section ('H') */
  char                user_string[P4EST_NUM_USER_STRING_BYTES]; /**< user string of the data section */
}
p4est_file_section_metadata_t;

/** Read metadata information of a file written by a matching forest.
 * Matching refers to the global count of quadrants; partition is irrelevant.
 *
 * This function parses the given file on rank 0 and broadcasts the information
 * on the number of data fields contained to all other ranks.  Collective call.
 *
 * This function catches all I/O and file format errors and returns a valid MPI
 * error class related to file handling.  Errors are collectively synchronized.
 *
 * If the number of bytes that the user intend to read is larger than the number
 * bytes left in the file, the function prints out an information about this
 * situation using P4EST_LERROR. In this case the function reads the bytes
 * that are possible to read but returns NULL to indicate an error.
 * If the file or block header information is not matching the passed parameters
 * the function sets \ref P4EST_FILE_ERR_IO for errcode.
 *
 * \param [in]  p4est               A p4est that is only required for the
 *                                  MPI communicator, and to verify the
 *                                  global quadrant count found in the file.
 * \param [in]  filename            Path to parallel file.
 * \param [in,out] user_string      At least \ref P4EST_NUM_USER_STRING_BYTES
 *                                  bytes. This array will be filled with the
 *                                  user string of the file after a successful
 *                                  call of this function.
 * \param [in,out] data_sections    After a successful function call this
 *                                  variable holds an array with a length
 *                                  corresponding to the number of data section
 *                                  in the file that are successfully found
 *                                  and seeked. The values in the array are the
 *                                  number of bytes of stored data per quadrant.
 *                                  Require elem_size->elem_size
 *                                  == sizeof (p4est_file_block_metadata_t)
 *                                  on input and preserve it on output.
 *                                  See p4est_file_block_metadata_t to obtain
 *                                  detailed information about the data blocks
 *                                  of the file.
 * \param [out] errcode             An errcode that can be interpreted by \ref
 *                                  p4est_file_error_string.
 * \return                          0 for a successful call and -1 in case of
 *                                  an error. See also errcode argument.
 */
int                 p4est_file_info (p4est_t * p4est, const char *filename,
                                     char
                                     user_string[P4EST_NUM_USER_STRING_BYTES],
                                     sc_array_t * data_sections,
                                     int *errcode);

/** Turn p4est_file errcode into a string.
 * 
 * \param [in] errclass     An errcode that is output by a
 *                          p4est_file function.
 * \param [in,out] string   At least sc_MPI_MAX_ERROR_STRING bytes.
 * \param [out] resultlen   Length of string on return.
 * \return                  P4EST_FILE_ERR_SUCCESS on success or
 *                          something else on invalid arguments.
 */
int                 p4est_file_error_string (int errclass,
                                             char
                                             string[sc_MPI_MAX_ERROR_STRING],
                                             int *resultlen);

/** Close a file opened for parallel write/read and free the context.
 * \param [in,out] fc       Context previously created by \ref
 *                          p4est_file_open_create or \ref
 *                          p4est_file_open_read (_ext).  Is freed.
 * \param [out] errcode     An errcode that can be interpreted by \ref
 *                          p4est_file_error_string.
 * \return                  0 for a successful call and -1 in case of
 *                          an error. See also errcode argument.
 */
int                 p4est_file_close (p4est_file_context_t * fc,
                                      int *errcode);

SC_EXTERN_C_END;

#endif /* !P4EST_IO_H */
