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

/*
 * There are several configure options that affect the VTK output.
 * --enable-vtk-binary (default) uses binary/base64 encoding.
 * --enable-vtk-compression (default) uses zlib compression.
 * --enable-vtk-doubles (default) uses 64-bit binary reals.
 * These may be disabled using the --disable-vtk-* forms.
 *
 * These options translate into the preprocessor #defines
 * P4EST_ENABLE_VTK_BINARY,
 * P4EST_ENABLE_VTK_COMPRESSION,
 * P4EST_ENABLE_VTK_DOUBLES.
 * The older version of the macros without _ENABLE_ is obsolete.
 *
 * We define further macros of the form P4EST_VTK_ in this file.
 * These are local to the code in this file and not exported.
 */

#ifdef P4_TO_P8
#include <p8est_vtk.h>
#include <p8est_nodes.h>
#define P4EST_VTK_CELL_TYPE     11      /* VTK_VOXEL */
#define P4EST_VTK_CELL_TYPE_HO  72      /* VTK_LAGRANGE_HEXAHEDRON */
#else
#include <p4est_vtk.h>
#include <p4est_nodes.h>
#define P4EST_VTK_CELL_TYPE      8      /* VTK_PIXEL */
#define P4EST_VTK_CELL_TYPE_HO  70      /* VTK_LAGRANGE_QUADRILATERAL */
#endif /* !P4_TO_P8 */

/* default parameters for the vtk context */
static const double p4est_vtk_scale = 0.95;
static const int    p4est_vtk_continuous = 0;

/* default parameters for p4est_vtk_write_file */
static const int    p4est_vtk_write_tree = 1;
static const int    p4est_vtk_write_level = 1;
static const int    p4est_vtk_write_rank = 1;
static const int    p4est_vtk_wrap_rank = 0;

/** Write a cell scalar field to the VTU file.
 *
 * Writing a VTK file is split into a few routines.
 * This allows there to be an arbitrary number of fields.
 * When in doubt, please use \ref p4est_vtk_write_cell_data instead.
 *
 * \param [in,out] cont    A VTK context created by \ref p4est_vtk_context_new.
 * \param [in] scalar_name The name of the scalar field.
 * \param [in] values      The cell values that will be written.
 *
 * \return          On success, the context that has been passed in.
 *                  On failure, returns NULL and deallocates the context.
 */
static
  p4est_vtk_context_t *p4est_vtk_write_cell_scalar (p4est_vtk_context_t *
                                                    cont,
                                                    const char *scalar_name,
                                                    sc_array_t * values);

/** Write a 3-vector cell field to the VTU file.
 *
 * Writing a VTK file is split into a few routines.
 * This allows there to be an arbitrary number of fields.
 * When in doubt, please use \ref p4est_vtk_write_cell_data instead.
 *
 * \param [in,out] cont    A VTK context created by \ref p4est_vtk_context_new.
 * \param [in] vector_name The name of the vector field.
 * \param [in] values      The cell values that will be written.
 *
 * \return          On success, the context that has been passed in.
 *                  On failure, returns NULL and deallocates the context.
 */
static
  p4est_vtk_context_t *p4est_vtk_write_cell_vector (p4est_vtk_context_t *
                                                    cont,
                                                    const char *vector_name,
                                                    sc_array_t * values);

/** Write a point scalar field to the VTU file.
 *
 * Writing a VTK file is split into a few routines.
 * This allows there to be an arbitrary number of fields.
 * When in doubt, please use \ref p4est_vtk_write_point_data instead.
 *
 * \param [in,out] cont    A VTK context created by \ref p4est_vtk_context_new.
 * \param [in] scalar_name The name of the scalar field.
 * \param [in] values      The point values that will be written.
 *
 * \return          On success, the context that has been passed in.
 *                  On failure, returns NULL and deallocates the context.
 */
static
  p4est_vtk_context_t *p4est_vtk_write_point_scalar (p4est_vtk_context_t *
                                                     cont,
                                                     const char *scalar_name,
                                                     sc_array_t * values);

/** Write a 3-vector point field to the VTU file.
 *
 * Writing a VTK file is split into a few routines.
 * This allows there to be an arbitrary number of fields.
 * When in doubt, please use \ref p4est_vtk_write_point_data instead.
 *
 * \param [in,out] cont    A VTK context created by \ref p4est_vtk_context_new.
 * \param [in] vector_name The name of the vector field.
 * \param [in] values      The point values that will be written.
 *
 * \return          On success, the context that has been passed in.
 *                  On failure, returns NULL and deallocates the context.
 */
static
  p4est_vtk_context_t *p4est_vtk_write_point_vector (p4est_vtk_context_t *
                                                     cont,
                                                     const char *vector_name,
                                                     sc_array_t * values);

#ifdef P4_TO_P8
#define p4est_vtk_context               p8est_vtk_context
#endif

#ifndef P4EST_ENABLE_VTK_DOUBLES
#define P4EST_VTK_FLOAT_NAME "Float32"
#define P4EST_VTK_FLOAT_TYPE float
#else
#define P4EST_VTK_FLOAT_NAME "Float64"
#define P4EST_VTK_FLOAT_TYPE double
#endif

/** Return an array of the proper P4EST_VTK_FLOAT_TYPE.
 * \param [in] input  Array of 3-tuples of type double.
 * \return            Array of 3-tuples of the proper
 *                    P4EST_VTK_FLOAT_TYPE.  This may
 *                    be a view if the types coincide.
 */
static sc_array_t   *
p4est_vtk_vector_array (sc_array_t *input)
{
#ifndef P4EST_ENABLE_VTK_DOUBLES
  P4EST_VTK_FLOAT_TYPE *floats;
  double             *doubles;
  size_t              zz, zn;
#endif
  sc_array_t         *output;

  P4EST_ASSERT (input != NULL);
  P4EST_ASSERT (input->elem_size == 3 * sizeof (double));

#ifndef P4EST_ENABLE_VTK_DOUBLES
  output = sc_array_new_count (3 * sizeof (P4EST_VTK_FLOAT_TYPE),
                               input->elem_count);
  floats = (P4EST_VTK_FLOAT_TYPE *) sc_array_index (output, 0);
  doubles = (double *) sc_array_index (input, 0);
  zn = 3 * input->elem_count;
  for (zz = 0; zz < zn; ++zz) {
    *floats++ = (P4EST_VTK_FLOAT_TYPE) *doubles++;
  }
#else
  output = sc_array_new_view (input, 0, input->elem_count);
#endif
  return output;
}

#ifndef P4EST_ENABLE_VTK_BINARY
#define P4EST_VTK_ASCII 1
#define P4EST_VTK_FORMAT_STRING "ascii"
#else
#define P4EST_VTK_FORMAT_STRING "binary"

static int
p4est_vtk_write_binary (FILE * vtkfile, char *numeric_data,
                        size_t byte_length)
{
  int                 retval;

  /* leading spaces */
  if (fprintf (vtkfile, "          ") < 0) {
    P4EST_LERROR ("Error writing pre-padding to binary VTK\n");
    return -1;
  }

  /* write binary data */
#ifndef P4EST_ENABLE_VTK_COMPRESSION
  retval = sc_vtk_write_binary (vtkfile, numeric_data, byte_length);
#else
  retval = sc_vtk_write_compressed (vtkfile, numeric_data, byte_length);
#endif /* P4EST_ENABLE_VTK_COMPRESSION */
  if (retval) {
    P4EST_LERROR ("Error encoding/writing data to binary VTK\n");
    return -1;
  }

  /* trailing line break */
  if (fprintf (vtkfile, "\n") < 0) {
    P4EST_LERROR ("Error writing post-padding to binary VTK\n");
    return -1;
  }

  /* success */
  return 0;
}

#endif /* P4EST_ENABLE_VTK_BINARY */

/** Opaque context type for writing VTK output with multiple function calls.
 *
 * This structure holds all the information needed for the p4est vtk context.
 * It is used to relay necessary vtk information to the \b p4est_vtk_write_*
 * functions. This structure is initialized by \ref p4est_vtk_write_header and
 * destroyed by \b p4est_vtk_write_footer; it can also be destroyed manually
 * using the \b p4est_vtk_context_destroy function if necessary.
 *
 * The \a p4est member is a pointer to the local p4est.
 * The \a geom member is a pointer to the geometry used to create the p4est.
 * The \a Ncells member holds the number of VTK elements displayed.
 * The \a Npoints member holds the number of nodes present in the vtk output;
 * this is determined in \ref p4est_vtk_write_header using the \a scale parameter
 * and is used to assure the proper number of point variables are provided.
 * The \a filename member holds the vtk file basename: for error reporting.
 * The \a vtufilename, \a pvtufilename, and \a visitfilename members are the
 * vtk file names.
 * The \a vtufile, \a pvtufile, and \a visitfile members are the vtk file
 * pointers; opened by \ref p4est_vtk_write_header and closed by \b
 * p4est_vtk_write_footer.
 *
 */
struct p4est_vtk_context
{
  /* data passed initially */
  p4est_t            *p4est;       /**< The p4est structure must be alive. */
  char               *filename;    /**< Original filename provided is copied. */

  /* parameters that can optionally be set in a context */
  p4est_geometry_t   *geom;        /**< The geometry may be NULL. */
  double              scale;       /**< Parameter to shrink quadrants. */
  int                 continuous;  /**< Assume continuous point data? */

  /* cell meta data to override defaults */
  p4est_tnodes_t     *tnodes;      /**< If we're visualizing simplices. */
  p4est_locidx_t      Npoints;     /**< Total number of point locations. */
  p4est_locidx_t      Ncells;      /**< Remembered from writing header.
                                        May differ from the given
                                        \ref p4est local element count. */

  /* internal context data */
  int                 writing;     /**< True after p4est_vtk_write_header. */
  p4est_locidx_t      num_corners; /**< Number of local element corners. */
  p4est_locidx_t     *node_to_corner;     /**< Map a node to an element corner. */
  p4est_nodes_t      *nodes;       /**< NULL? depending on scale/continuous. */
  char                vtufilename[BUFSIZ];   /**< Each process writes one. */
  char                pvtufilename[BUFSIZ];  /**< Only root writes this one. */
  char                visitfilename[BUFSIZ]; /**< Only root writes this one. */
  FILE               *vtufile;     /**< File pointer for the VTU file. */
  FILE               *pvtufile;    /**< Paraview meta file. */
  FILE               *visitfile;   /**< Visit meta file. */
};

p4est_vtk_context_t *
p4est_vtk_context_new (p4est_t * p4est, const char *filename)
{
  p4est_vtk_context_t *cont;

  P4EST_ASSERT (p4est != NULL);
  P4EST_ASSERT (filename != NULL);

  /* Allocate, initialize the vtk context.  Important to zero all fields. */
  cont = P4EST_ALLOC_ZERO (p4est_vtk_context_t, 1);

  cont->p4est = p4est;
  cont->filename = P4EST_STRDUP (filename);

  cont->scale = p4est_vtk_scale;
  cont->continuous = p4est_vtk_continuous;

  return cont;
}

void
p4est_vtk_context_set_geom (p4est_vtk_context_t * cont,
                            p4est_geometry_t * geom)
{
  P4EST_ASSERT (cont != NULL);
  P4EST_ASSERT (!cont->writing);

  cont->geom = geom;
}

void
p4est_vtk_context_set_scale (p4est_vtk_context_t * cont, double scale)
{
  P4EST_ASSERT (cont != NULL);
  P4EST_ASSERT (!cont->writing);
  P4EST_ASSERT (0. < scale && scale <= 1.);

  cont->scale = scale;
}

void
p4est_vtk_context_set_continuous (p4est_vtk_context_t * cont, int continuous)
{
  P4EST_ASSERT (cont != NULL);
  P4EST_ASSERT (!cont->writing);

  cont->continuous = continuous;
}

void
p4est_vtk_context_destroy (p4est_vtk_context_t * context)
{
  P4EST_ASSERT (context != NULL);
  P4EST_ASSERT (context->p4est != NULL);

  /* since this function is called inside write_header and write_footer,
   * we cannot assume a consistent state of all member variables */

  P4EST_ASSERT (context->filename != NULL);
  P4EST_FREE (context->filename);

  /* deallocate node storage */
  if (context->nodes != NULL) {
    p4est_nodes_destroy (context->nodes);
  }
  P4EST_FREE (context->node_to_corner);

  /* Close all file pointers. */
  if (context->vtufile != NULL) {
    if (fclose (context->vtufile)) {
      P4EST_LERRORF (P4EST_STRING "_vtk: Error closing <%s>.\n",
                     context->vtufilename);
    }
    context->vtufile = NULL;
  }

  /* Close paraview master file */
  if (context->pvtufile != NULL) {
    /* Only the root process opens/closes these files. */
    P4EST_ASSERT (context->p4est->mpirank == 0);
    if (fclose (context->pvtufile)) {
      P4EST_LERRORF (P4EST_STRING "_vtk: Error closing <%s>.\n",
                     context->pvtufilename);
    }
    context->pvtufile = NULL;
  }

  /* Close visit master file */
  if (context->visitfile != NULL) {
    /* Only the root process opens/closes these files. */
    P4EST_ASSERT (context->p4est->mpirank == 0);
    if (fclose (context->visitfile)) {
      P4EST_LERRORF (P4EST_STRING "_vtk: Error closing <%s>.\n",
                     context->visitfilename);
    }
    context->visitfile = NULL;
  }

  /* Free context structure. */
  P4EST_FREE (context);
}

void
p4est_vtk_write_file (p4est_t * p4est, p4est_geometry_t * geom,
                      const char *filename)
{
  int                 retval;
  p4est_vtk_context_t *cont;

  /* allocate context and set parameters */
  cont = p4est_vtk_context_new (p4est, filename);
  p4est_vtk_context_set_geom (cont, geom);

  /* We do not write point data, so it is safe to set continuous to true.
   * This will not save any space though since the default scale is < 1. */
  p4est_vtk_context_set_continuous (cont, 1);

  /* write header, that is, vertex positions and quadrant-to-vertex map */
  cont = p4est_vtk_write_header (cont);
  SC_CHECK_ABORT (cont != NULL, P4EST_STRING "_vtk: Error writing header");

  /* write the tree/level/rank data */
  cont =
    p4est_vtk_write_cell_dataf (cont, p4est_vtk_write_tree,
                                p4est_vtk_write_level, p4est_vtk_write_rank,
                                p4est_vtk_wrap_rank, 0, 0, cont);
  SC_CHECK_ABORT (cont != NULL, P4EST_STRING "_vtk: Error writing cell data");

  /* properly write rest of the files' contents */
  retval = p4est_vtk_write_footer (cont);
  SC_CHECK_ABORT (!retval, P4EST_STRING "_vtk: Error writing footer");
}

/** Open a VTK file on each rank and write part of the topmatter.
 *
 * \param [in,out] cont         Valid context not yet used for writing.
 *                              It is constructed and options have been set,
 *                              but there must not have been a write call.
 *                              In this function is set to writing.
 * \param [in] points           An array of 3-tuples of the datatype
 *                              P4EST_VTK_FLOAT_TYPE, one tuple
 *                              for each point written to the VTK file.
 * \param [in] Ncells           Number of cells anchored at respective points.
 * \return                      The context on success, NULL otherwise.
 *                              On NULL, we close the file and VTK state,
 *                              the caller just deallocates their inputs.
 */
static p4est_vtk_context_t *
p4est_vtk_write_header_points (p4est_vtk_context_t * cont,
                               sc_array_t *points, p4est_locidx_t Ncells)
{
  int                 mpirank;
  const char         *filename;
#ifdef P4EST_VTK_ASCII
  p4est_locidx_t      il;
  double              wx, wy, wz;
  P4EST_VTK_FLOAT_TYPE *float_data;
#endif
  p4est_locidx_t      Npoints;

  /* basic checks */
  P4EST_ASSERT (cont != NULL && cont->p4est != NULL);
  P4EST_ASSERT (cont->filename != NULL);
  P4EST_ASSERT (!cont->writing);
  P4EST_ASSERT (points != NULL);
  P4EST_ASSERT (points->elem_size == 3 * sizeof (P4EST_VTK_FLOAT_TYPE));
  P4EST_ASSERT (Ncells >= 0);

  /* variables from context */
  mpirank = cont->p4est->mpirank;
  filename = cont->filename;
  Npoints = (p4est_locidx_t) points->elem_count;

  /* we enter writing mode */
  cont->Npoints = Npoints;
  cont->Ncells = Ncells;
  cont->writing = 1;

  /* have each proc write to its own file */
  snprintf (cont->vtufilename, BUFSIZ, "%s_%04d.vtu", filename, mpirank);
  /*
   * Use "w" for writing the initial part of the file.
   * For further parts, use "r+" and fseek so write_compressed succeeds.
   */
  cont->vtufile = fopen (cont->vtufilename, "wb");
  if (cont->vtufile == NULL) {
    P4EST_LERRORF ("Fail opening %s for output\n", cont->vtufilename);
    p4est_vtk_context_destroy (cont);
    return NULL;
  }

  /* begin with VTK header */
  fprintf (cont->vtufile, "<?xml version=\"1.0\"?>\n");
  fprintf (cont->vtufile,
           "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\"");
#if defined P4EST_ENABLE_VTK_BINARY && defined P4EST_ENABLE_VTK_COMPRESSION
  fprintf (cont->vtufile, " compressor=\"vtkZLibDataCompressor\"");
#endif
#ifdef SC_IS_BIGENDIAN
  fprintf (cont->vtufile, " byte_order=\"BigEndian\">\n");
#else
  fprintf (cont->vtufile, " byte_order=\"LittleEndian\">\n");
#endif
  fprintf (cont->vtufile, "  <UnstructuredGrid>\n");
  fprintf (cont->vtufile,
           "    <Piece NumberOfPoints=\"%lld\" NumberOfCells=\"%lld\">\n",
           (long long) Npoints, (long long) Ncells);

  /* begin point context */
  fprintf (cont->vtufile, "      <Points>\n");
  fprintf (cont->vtufile, "        <DataArray type=\"%s\" Name=\"Position\""
           " NumberOfComponents=\"3\" format=\"%s\">\n",
           P4EST_VTK_FLOAT_NAME, P4EST_VTK_FORMAT_STRING);

#ifdef P4EST_VTK_ASCII
  for (il = 0; il < Npoints; ++il) {
    float_data =
      (P4EST_VTK_FLOAT_TYPE *) sc_array_index (points, (size_t) il);
    wx = float_data[0];
    wy = float_data[1];
    wz = float_data[2];

    fprintf (cont->vtufile,
#ifdef P4EST_ENABLE_VTK_DOUBLES
             "     %24.16e %24.16e %24.16e\n",
#else
             "          %16.8e %16.8e %16.8e\n",
#endif
             wx, wy, wz);
  }
#else
  /* TO DO: Don't allocate the full size of the array, only allocate
   * the chunk that will be passed to zlib and do this a chunk
   * at a time.
   */
  if (p4est_vtk_write_binary
      (cont->vtufile, points->array,
       sizeof (P4EST_VTK_FLOAT_TYPE) * Npoints * 3)) {
    P4EST_LERROR (P4EST_STRING "_vtk: Error encoding points\n");
    p4est_vtk_context_destroy (cont);
    return NULL;
  }
#endif

  /* close point context */
  fprintf (cont->vtufile, "        </DataArray>\n");
  fprintf (cont->vtufile, "      </Points>\n");

  /* catch all file errors */
  if (ferror (cont->vtufile)) {
    P4EST_LERROR (P4EST_STRING "_vtk: Error writing point header\n");
    p4est_vtk_context_destroy (cont);
    return NULL;
  }

  /* successful return */
  return cont;
}

/** Write the cell data to the header of a VTK file.
 *
 * \param [in,out] cont         Valid context whose header has been started.
 * \param [in] vtk_cell_type    The cell type according to the VTK format.
 * \param [in] num_cell_corners Number of corners to each data cell.
 * \param [in] cells            An array of p4est_locidx n-tuples, where
 *                              each tuples has \c num_cell_corners entries.
 *                              The entries index into the header point data.
 * \return                      The context on success, NULL otherwise.
 *                              On NULL, we close the file and VTK state,
 *                              the caller just deallocates their inputs.
 */
static p4est_vtk_context_t *
p4est_vtk_write_header_cells (p4est_vtk_context_t *cont, int vtk_cell_type,
                              int num_cell_corners, sc_array_t *cells)
{
#ifdef P4EST_VTK_ASCII
  int                 k;
  p4est_locidx_t      sk;
#else
  uint8_t            *uint8_data;
  p4est_locidx_t     *locidx_data;
#endif
  p4est_locidx_t      il;
  p4est_locidx_t      Ncells;

  /* basic checks */
  P4EST_ASSERT (cont != NULL && cont->p4est != NULL);
  P4EST_ASSERT (cont->filename != NULL);
  P4EST_ASSERT (cont->writing);
  P4EST_ASSERT (num_cell_corners > 0);
  P4EST_ASSERT (cells != NULL);
  P4EST_ASSERT (cells->elem_size ==
                num_cell_corners * sizeof (p4est_locidx_t));
  P4EST_ASSERT (cells->elem_count == (size_t) cont->Ncells);

  /* variables from context */
  Ncells = (p4est_locidx_t) cells->elem_count;
  P4EST_ASSERT (Ncells == cont->Ncells);
  cont->num_corners = Ncells * num_cell_corners;

  /* begin cell context */
  fprintf (cont->vtufile, "      <Cells>\n");

  /* write connectivity data */
  fprintf (cont->vtufile,
           "        <DataArray type=\"%s\" Name=\"connectivity\""
           " format=\"%s\">\n", P4EST_VTK_LOCIDX, P4EST_VTK_FORMAT_STRING);
#ifdef P4EST_VTK_ASCII
  for (il = 0; il < Ncells; ++il) {
    p4est_locidx_t *cell_data = (p4est_locidx_t *) sc_array_index (cells, il);
    fprintf (cont->vtufile, "         ");
    for (k = 0; k < num_cell_corners; ++k) {
      fprintf (cont->vtufile, " %lld", (long long) cell_data[k]);
    }
    fprintf (cont->vtufile, "\n");
  }
#else
  if (p4est_vtk_write_binary
      (cont->vtufile, cells->array,
       sizeof (p4est_locidx_t) * Ncells * num_cell_corners)) {
    P4EST_LERROR (P4EST_STRING "_vtk: Error encoding connectivity\n");
    p4est_vtk_context_destroy (cont);
    return NULL;
  }
#endif
  fprintf (cont->vtufile, "        </DataArray>\n");

  /* write offset data */
  fprintf (cont->vtufile, "        <DataArray type=\"%s\" Name=\"offsets\""
           " format=\"%s\">\n", P4EST_VTK_LOCIDX, P4EST_VTK_FORMAT_STRING);
#ifdef P4EST_VTK_ASCII
  fprintf (cont->vtufile, "         ");
  for (il = 1, sk = 1; il <= Ncells; ++il, ++sk) {
    fprintf (cont->vtufile, " %lld", (long long) (num_cell_corners * il));
    if (!(sk % 8) && il != Ncells) {
      fprintf (cont->vtufile, "\n         ");
    }
  }
  fprintf (cont->vtufile, "\n");
#else
  locidx_data = P4EST_ALLOC (p4est_locidx_t, Ncells);
  for (il = 1; il <= Ncells; ++il) {
    locidx_data[il - 1] = num_cell_corners * il;
  }
  if (p4est_vtk_write_binary (cont->vtufile, (char *) locidx_data,
                              sizeof (p4est_locidx_t) * Ncells)) {
    P4EST_LERROR (P4EST_STRING "_vtk: Error encoding offsets\n");
    P4EST_FREE (locidx_data);
    p4est_vtk_context_destroy (cont);
    return NULL;
  }
  P4EST_FREE (locidx_data);
#endif
  fprintf (cont->vtufile, "        </DataArray>\n");

  /* write type data */
  fprintf (cont->vtufile, "        <DataArray type=\"UInt8\" Name=\"types\""
           " format=\"%s\">\n", P4EST_VTK_FORMAT_STRING);
#ifdef P4EST_VTK_ASCII
  fprintf (cont->vtufile, "         ");
  for (il = 0, sk = 1; il < Ncells; ++il, ++sk) {
    fprintf (cont->vtufile, " %d", vtk_cell_type);
    if (!(sk % 20) && il != (Ncells - 1))
      fprintf (cont->vtufile, "\n         ");
  }
  fprintf (cont->vtufile, "\n");
#else
  uint8_data = P4EST_ALLOC (uint8_t, Ncells);
  for (il = 0; il < Ncells; ++il) {
    uint8_data[il] = vtk_cell_type;
  }
  if (p4est_vtk_write_binary (cont->vtufile, (char *) uint8_data,
                              sizeof (*uint8_data) * Ncells)) {
    P4EST_LERROR (P4EST_STRING "_vtk: Error encoding types\n");
    P4EST_FREE (uint8_data);
    p4est_vtk_context_destroy (cont);
    return NULL;
  }
  P4EST_FREE (uint8_data);
#endif

  /* close cell context */
  fprintf (cont->vtufile, "        </DataArray>\n");
  fprintf (cont->vtufile, "      </Cells>\n");

  /* catch all file errors */
  if (ferror (cont->vtufile)) {
    P4EST_LERROR (P4EST_STRING "_vtk: Error writing cell header\n");
    p4est_vtk_context_destroy (cont);
    return NULL;
  }

  /* successful return */
  return cont;
}

static p4est_vtk_context_t *
p4est_vtk_write_header_root (p4est_vtk_context_t *cont)
{
  const char         *filename;

  /* basic checks */
  P4EST_ASSERT (cont != NULL && cont->p4est != NULL);
  P4EST_ASSERT (cont->p4est->mpirank == 0);
  P4EST_ASSERT (cont->filename != NULL);
  P4EST_ASSERT (cont->writing);

  /* access variables from context */
  filename = cont->filename;

  /* open parallel VTU file */
  snprintf (cont->pvtufilename, BUFSIZ, "%s.pvtu", filename);
  cont->pvtufile = fopen (cont->pvtufilename, "wb");
  if (!cont->pvtufile) {
    P4EST_LERRORF ("Fail opening %s for output\n", cont->pvtufilename);
    p4est_vtk_context_destroy (cont);
    return NULL;
  }

  /* write parallel VTU header */
  fprintf (cont->pvtufile, "<?xml version=\"1.0\"?>\n");
  fprintf (cont->pvtufile,
           "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\"");
#if defined P4EST_ENABLE_VTK_BINARY && defined P4EST_ENABLE_VTK_COMPRESSION
  fprintf (cont->pvtufile, " compressor=\"vtkZLibDataCompressor\"");
#endif
#ifdef SC_IS_BIGENDIAN
  fprintf (cont->pvtufile, " byte_order=\"BigEndian\">\n");
#else
  fprintf (cont->pvtufile, " byte_order=\"LittleEndian\">\n");
#endif

  fprintf (cont->pvtufile, "  <PUnstructuredGrid GhostLevel=\"0\">\n");
  fprintf (cont->pvtufile, "    <PPoints>\n");
  fprintf (cont->pvtufile, "      <PDataArray type=\"%s\" Name=\"Position\""
           " NumberOfComponents=\"3\" format=\"%s\"/>\n",
           P4EST_VTK_FLOAT_NAME, P4EST_VTK_FORMAT_STRING);
  fprintf (cont->pvtufile, "    </PPoints>\n");

  if (ferror (cont->pvtufile)) {
    P4EST_LERROR (P4EST_STRING "_vtk: Error writing parallel header\n");
    p4est_vtk_context_destroy (cont);
    return NULL;
  }

  /* Create a master file for visualization in Visit; this will be used
   * only in p4est_vtk_write_footer ().
   */
  snprintf (cont->visitfilename, BUFSIZ, "%s.visit", filename);
  cont->visitfile = fopen (cont->visitfilename, "wb");
  if (!cont->visitfile) {
    P4EST_LERRORF ("Fail opening %s for output\n", cont->visitfilename);
    p4est_vtk_context_destroy (cont);
    return NULL;
  }

  /* successful return */
  return cont;
}

#ifndef P4_TO_P8
#define P4EST_VTK_SIMPLEX_TYPE 5
#else
#define P4EST_VTK_SIMPLEX_TYPE 10
#endif

p4est_vtk_context_t *
p4est_vtk_write_header_simplices (p4est_vtk_context_t * cont, sc_array_t *simplices, sc_array_t *vertices)
{
  int                 mpirank;
  const char         *filename;
  p4est_locidx_t      Ncells;
  p4est_t            *p4est;
#ifdef P4EST_VTK_ASCII
  int                 k;
  p4est_locidx_t      sk;
  double              wx, wy, wz;
#else
  uint8_t            *uint8_data;
  p4est_locidx_t     *locidx_data;
#endif
  p4est_locidx_t      il, Npoints;
  P4EST_VTK_FLOAT_TYPE *float_data;

  /* check a whole bunch of assertions, here and below */
  P4EST_ASSERT (cont != NULL);
  P4EST_ASSERT (!cont->writing);

  /* from now on this context is officially in use for writing */
  cont->writing = 1;

  /* grab context variables */
  p4est = cont->p4est;
  filename = cont->filename;
  P4EST_ASSERT (filename != NULL);

  /* grab details from the forest */
  P4EST_ASSERT (p4est != NULL);
  mpirank = p4est->mpirank;
  Ncells = cont->Ncells = simplices->elem_count;
  Npoints = cont->Npoints = vertices->elem_count;
  cont->num_corners = (P4EST_DIM + 1) * Ncells;

  /* Have each proc write to its own file */
  snprintf (cont->vtufilename, BUFSIZ, "%s_%04d.vtu", filename, mpirank);
  /* Use "w" for writing the initial part of the file.
   * For further parts, use "r+" and fseek so write_compressed succeeds.
   */
  cont->vtufile = fopen (cont->vtufilename, "wb");
  if (cont->vtufile == NULL) {
    P4EST_LERRORF ("Could not open %s for output\n", cont->vtufilename);
    p4est_vtk_context_destroy (cont);
    return NULL;
  }

  fprintf (cont->vtufile, "<?xml version=\"1.0\"?>\n");
  fprintf (cont->vtufile,
           "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\"");
#if defined P4EST_ENABLE_VTK_BINARY && defined P4EST_ENABLE_VTK_COMPRESSION
  fprintf (cont->vtufile, " compressor=\"vtkZLibDataCompressor\"");
#endif
#ifdef SC_IS_BIGENDIAN
  fprintf (cont->vtufile, " byte_order=\"BigEndian\">\n");
#else
  fprintf (cont->vtufile, " byte_order=\"LittleEndian\">\n");
#endif
  fprintf (cont->vtufile, "  <UnstructuredGrid>\n");
  fprintf (cont->vtufile,
           "    <Piece NumberOfPoints=\"%lld\" NumberOfCells=\"%lld\">\n",
           (long long) Npoints, (long long) Ncells);
  fprintf (cont->vtufile, "      <Points>\n");

  float_data = P4EST_ALLOC (P4EST_VTK_FLOAT_TYPE, 3 * Npoints);

  /* write point position data */
  fprintf (cont->vtufile, "        <DataArray type=\"%s\" Name=\"Position\""
           " NumberOfComponents=\"3\" format=\"%s\">\n",
           P4EST_VTK_FLOAT_NAME, P4EST_VTK_FORMAT_STRING);

  {
    const double *v = (const double *) vertices->array;
    for (size_t i = 0; i < vertices->elem_count * 3; i++) {
      float_data[i] = v[i];
    }
  }

#ifdef P4EST_VTK_ASCII
  for (il = 0; il < Npoints; ++il) {
    wx = float_data[3 * il + 0];
    wy = float_data[3 * il + 1];
    wz = float_data[3 * il + 2];

    fprintf (cont->vtufile,
#ifdef P4EST_VTK_DOUBLES
             "     %24.16e %24.16e %24.16e\n",
#else
             "          %16.8e %16.8e %16.8e\n",
#endif
             wx, wy, wz);
  }
#else
  /* TODO: Don't allocate the full size of the array, only allocate
   * the chunk that will be passed to zlib and do this a chunk
   * at a time.
   */
  if (p4est_vtk_write_binary (cont->vtufile, (char *) float_data,
                              sizeof (*float_data) * 3 * Npoints)) {
    P4EST_LERROR (P4EST_STRING "_vtk: Error encoding points\n");
    P4EST_FREE (float_data);
    p4est_vtk_context_destroy (cont);
    return NULL;
  }
#endif
  P4EST_FREE (float_data);

  fprintf (cont->vtufile, "        </DataArray>\n");
  fprintf (cont->vtufile, "      </Points>\n");
  fprintf (cont->vtufile, "      <Cells>\n");

  /* write connectivity data */
  fprintf (cont->vtufile,
           "        <DataArray type=\"%s\" Name=\"connectivity\""
           " format=\"%s\">\n", P4EST_VTK_LOCIDX, P4EST_VTK_FORMAT_STRING);
#ifdef P4EST_VTK_ASCII
  for (sk = 0, il = 0; il < Ncells; ++il) {
    p4est_locidx_t *simplex = (p4est_locidx_t *) sc_array_index (simplices, il);
    fprintf (cont->vtufile, "         ");
    for (k = 0; k < P4EST_DIM+1; ++sk, ++k) {
      fprintf (cont->vtufile, " %lld", (long long) simplex[k]);
    }
    fprintf (cont->vtufile, "\n");
  }
#else
  if (p4est_vtk_write_binary (cont->vtufile, (char *) simplices->array,
                              sizeof (p4est_locidx_t) * Ncells * (P4EST_DIM + 1))) {
    P4EST_LERROR (P4EST_STRING "_vtk: Error encoding connectivity\n");
    p4est_vtk_context_destroy (cont);
    return NULL;
  }
#endif
  fprintf (cont->vtufile, "        </DataArray>\n");

  /* write offset data */
  fprintf (cont->vtufile, "        <DataArray type=\"%s\" Name=\"offsets\""
           " format=\"%s\">\n", P4EST_VTK_LOCIDX, P4EST_VTK_FORMAT_STRING);
#ifdef P4EST_VTK_ASCII
  fprintf (cont->vtufile, "         ");
  for (il = 1, sk = 1; il <= Ncells; ++il, ++sk) {
    fprintf (cont->vtufile, " %lld", (long long) ((P4EST_DIM + 1) * il));
    if (!(sk % 8) && il != Ncells)
      fprintf (cont->vtufile, "\n         ");
  }
  fprintf (cont->vtufile, "\n");
#else
  locidx_data = P4EST_ALLOC (p4est_locidx_t, Ncells);
  for (il = 1; il <= Ncells; ++il) {
    locidx_data[il - 1] = (P4EST_DIM + 1) * il;  /* same type */
  }
  if (p4est_vtk_write_binary (cont->vtufile, (char *) locidx_data,
                              sizeof (p4est_locidx_t) * Ncells)) {
    P4EST_LERROR (P4EST_STRING "_vtk: Error encoding offsets\n");
    P4EST_FREE (locidx_data);
    p4est_vtk_context_destroy (cont);
    return NULL;
  }
  P4EST_FREE (locidx_data);
#endif
  fprintf (cont->vtufile, "        </DataArray>\n");

  /* write type data */
  fprintf (cont->vtufile, "        <DataArray type=\"UInt8\" Name=\"types\""
           " format=\"%s\">\n", P4EST_VTK_FORMAT_STRING);
#ifdef P4EST_VTK_ASCII
  fprintf (cont->vtufile, "         ");
  for (il = 0, sk = 1; il < Ncells; ++il, ++sk) {
    fprintf (cont->vtufile, " %d", P4EST_VTK_SIMPLEX_TYPE);
    if (!(sk % 20) && il != (Ncells - 1))
      fprintf (cont->vtufile, "\n         ");
  }
  fprintf (cont->vtufile, "\n");
#else
  uint8_data = P4EST_ALLOC (uint8_t, Ncells);
  for (il = 0; il < Ncells; ++il) {
    uint8_data[il] = P4EST_VTK_SIMPLEX_TYPE;
  }
  if (p4est_vtk_write_binary (cont->vtufile, (char *) uint8_data,
                              sizeof (*uint8_data) * Ncells)) {
    P4EST_LERROR (P4EST_STRING "_vtk: Error encoding types\n");
    P4EST_FREE (uint8_data);
    p4est_vtk_context_destroy (cont);
    return NULL;
  }
  P4EST_FREE (uint8_data);
#endif
  fprintf (cont->vtufile, "        </DataArray>\n");
  fprintf (cont->vtufile, "      </Cells>\n");

  if (ferror (cont->vtufile)) {
    P4EST_LERROR (P4EST_STRING "_vtk: Error writing header\n");
    p4est_vtk_context_destroy (cont);
    return NULL;
  }

  /* Only have the root write to the parallel vtk file */
  if (mpirank == 0) {
    snprintf (cont->pvtufilename, BUFSIZ, "%s.pvtu", filename);

    cont->pvtufile = fopen (cont->pvtufilename, "wb");
    if (!cont->pvtufile) {
      P4EST_LERRORF ("Could not open %s for output\n", cont->pvtufilename);
      p4est_vtk_context_destroy (cont);
      return NULL;
    }

    fprintf (cont->pvtufile, "<?xml version=\"1.0\"?>\n");
    fprintf (cont->pvtufile,
             "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\"");
#if defined P4EST_ENABLE_VTK_BINARY && defined P4EST_ENABLE_VTK_COMPRESSION
    fprintf (cont->pvtufile, " compressor=\"vtkZLibDataCompressor\"");
#endif
#ifdef SC_IS_BIGENDIAN
    fprintf (cont->pvtufile, " byte_order=\"BigEndian\">\n");
#else
    fprintf (cont->pvtufile, " byte_order=\"LittleEndian\">\n");
#endif

    fprintf (cont->pvtufile, "  <PUnstructuredGrid GhostLevel=\"0\">\n");
    fprintf (cont->pvtufile, "    <PPoints>\n");
    fprintf (cont->pvtufile, "      <PDataArray type=\"%s\" Name=\"Position\""
             " NumberOfComponents=\"3\" format=\"%s\"/>\n",
             P4EST_VTK_FLOAT_NAME, P4EST_VTK_FORMAT_STRING);
    fprintf (cont->pvtufile, "    </PPoints>\n");

    if (ferror (cont->pvtufile)) {
      P4EST_LERROR (P4EST_STRING "_vtk: Error writing parallel header\n");
      p4est_vtk_context_destroy (cont);
      return NULL;
    }

    /* Create a master file for visualization in Visit; this will be used
     * only in p4est_vtk_write_footer().
     */
    snprintf (cont->visitfilename, BUFSIZ, "%s.visit", filename);
    cont->visitfile = fopen (cont->visitfilename, "wb");
    if (!cont->visitfile) {
      P4EST_LERRORF ("Could not open %s for output\n", cont->visitfilename);
      p4est_vtk_context_destroy (cont);
      return NULL;
    }
  }

  return cont;
}

p4est_vtk_context_t *
p4est_vtk_write_header_tnodes (p4est_vtk_context_t *cont,
                               p4est_tnodes_t *tnodes)
{
  int                 mpirank;
  p4est_locidx_t      Ncells;
  p4est_vtk_context_t *retcont;
  p4est_geometry_t   *geom;
  sc_array_t         *simplices;
  sc_array_t         *points;
  sc_array_t         *cells;

  /* check a whole bunch of assertions, here and below */
  P4EST_ASSERT (cont != NULL);
  P4EST_ASSERT (cont->p4est != NULL);

  /* check tnodes input and remember it */
  P4EST_ASSERT (tnodes != NULL);
  P4EST_ASSERT (tnodes->simplices != NULL);
  P4EST_ASSERT (tnodes->coordinates != NULL);
  cont->tnodes = tnodes;

  /* grab details from the forest */
  mpirank = cont->p4est->mpirank;
  Ncells = (p4est_locidx_t)
    (simplices = tnodes->simplices)->elem_count;

  /* access input point locations */
  if ((geom = cont->geom) == NULL) {
    /*
     * When visualizing a p4est, a NULL geometry means we default
     * to the geometry defined by the p4est->connectivity.  We do this
     * as a convenience, and since there would be no way to input a
     * p4est with coordinates already mapped to geometry space.
     *
     * With a simplex mesh, this is different: the coordinate input
     * may full well be already transformed.  Thus, with a NULL geometry
     * in the VTK context we do not touch the coordinates at all.
     *
     * If the user wishes to transform by the p4est connectivity, they
     * may use p4est_geometry_new_connectivity and pass us the result.
     */
    points = p4est_vtk_vector_array (tnodes->coordinates);
  }
  else {
    int                 k;
    int8_t             *done;
#ifndef P4EST_ENABLE_VTK_DOUBLES
    double              xyz[3];
#else
    double             *xyz;
#endif
    P4EST_VTK_FLOAT_TYPE *pco;
    p4est_topidx_t      tt, lftm;
    p4est_locidx_t      is, ci, numc;
    p4est_locidx_t     *simc, stoff;

    /* we loop by tree since it is a parameter to the geometry transform */
    numc = tnodes->coordinates->elem_count;
    done = P4EST_ALLOC_ZERO (int8_t, numc);
    points = sc_array_new_count (3 * sizeof (P4EST_VTK_FLOAT_TYPE), numc);

    /*
     * Loop through the trees and their simplices in order to
     * reference all coordinates eventually, which we transform.
     */
    P4EST_ASSERT (tnodes->local_tree_offset[0] == 0);
    for (is = 0, lftm = (tt = tnodes->local_first_tree) - 1;
         tt <= tnodes->local_last_tree; ++tt) {
      /* tree offsets point into the range of simplices */
      stoff = tnodes->local_tree_offset[tt - lftm];
      while (is < stoff) {
        simc = (p4est_locidx_t *) sc_array_index (simplices, is);
        for (k = 0; k < P4EST_DIM + 1; ++k) {
          P4EST_ASSERT (0 <= simc[k] && simc[k] < numc);
          if (done[ci = simc[k]]) {
            /* this coordinate is computed already */
            continue;
          }
          pco = (P4EST_VTK_FLOAT_TYPE *) sc_array_index (points, ci);
#ifdef P4EST_ENABLE_VTK_DOUBLES
          xyz = pco;
#endif
          geom->X (geom, tt, (double *)
                   sc_array_index (tnodes->coordinates, ci), xyz);
#ifndef P4EST_ENABLE_VTK_DOUBLES
          pco[0] = (P4EST_VTK_FLOAT_TYPE) xyz[0];
          pco[1] = (P4EST_VTK_FLOAT_TYPE) xyz[1];
          pco[2] = (P4EST_VTK_FLOAT_TYPE) xyz[2];
#endif
          done[ci] = 1;
        }
        ++is;
      }
    }
    P4EST_ASSERT (is == Ncells);
    P4EST_FREE (done);
  }

  /* possibly expand point locations to simplex vertices */
  if (cont->scale < 1. || !cont->continuous) {
    int                 k, i;
    p4est_locidx_t      is, *simc, *cell;
    P4EST_VTK_FLOAT_TYPE *vinput[P4EST_DIM + 1], vbary[3];
    P4EST_VTK_FLOAT_TYPE *voutput[P4EST_DIM + 1];
    P4EST_VTK_FLOAT_TYPE fscale, fbary;
    sc_array_t          *vpoints;

    /* basic settings */
    fbary = (1. - (fscale = cont->scale)) / (P4EST_DIM + 1);

    /* when we scale the quadrants we need each corner separately */
    cells = sc_array_new_count ((P4EST_DIM + 1) * sizeof (p4est_locidx_t), Ncells);
    vpoints = sc_array_new_count
      (3 * sizeof (P4EST_VTK_FLOAT_TYPE), (P4EST_DIM + 1) * Ncells);
    for (is = 0; is < Ncells; ++is) {
      /* simplex indexes into coordinates array mapped to points array */
      simc = (p4est_locidx_t *) sc_array_index (simplices, is);
      cell = (p4est_locidx_t *) sc_array_index (cells, is);

      /* compute simplex barycenter */
      vbary[0] = vbary[1] = vbary[2] = 0.;
      for (k = 0; k < P4EST_DIM + 1; ++k) {
        vinput[k] = (P4EST_VTK_FLOAT_TYPE *) sc_array_index (points, simc[k]);
        for (i = 0; i < 3; ++i) {
          vbary[i] += vinput[k][i];
        }
      }
      for (i = 0; i < 3; ++i) {
        vbary[i] *= fbary;
      }

      /* move corners towards barycenter */
      for (k = 0; k < P4EST_DIM + 1; ++k) {
        voutput[k] = (P4EST_VTK_FLOAT_TYPE *)
          sc_array_index (vpoints, cell[k] = (P4EST_DIM + 1) * is + k);
        for (i = 0; i < 3; ++i) {
          voutput[k][i] = fscale * vinput[k][i] + vbary[i];
        }
      }
    }

    /* assign new array to points */
    sc_array_destroy (points);
    points = vpoints;
  }
  else {
    cells = sc_array_new_view (simplices, 0, simplices->elem_count);
  }

  /* write final point locations to file */
  retcont = p4est_vtk_write_header_points (cont, points, Ncells);
  sc_array_destroy_null (&points);
  if (retcont != cont) {
    P4EST_ASSERT (retcont == NULL);
    P4EST_LERRORF ("Fail writing header points to %s\n", cont->vtufilename);
    return NULL;
  }

  /* write cell data to file */
  retcont = p4est_vtk_write_header_cells (cont, P4EST_VTK_SIMPLEX_TYPE,
                                          P4EST_DIM + 1, cells);
  sc_array_destroy (cells);
  if (retcont != cont) {
    P4EST_ASSERT (retcont == NULL);
    P4EST_LERRORF ("Fail writing header cells to %s\n", cont->vtufilename);
    return NULL;
  }

  /* Only have the root write to the parallel vtk file */
  if (mpirank == 0) {

    /* write parallel data on root rank */
    retcont = p4est_vtk_write_header_root (cont);
    if (retcont != cont) {
      P4EST_ASSERT (retcont == NULL);
      P4EST_LERROR ("Fail writing root headers\n");
      return NULL;
    }
  }

  /* successful return */
  return cont;
}

p4est_vtk_context_t *
p4est_vtk_write_header (p4est_vtk_context_t * cont)
{
  const double        intsize = 1.0 / P4EST_ROOT_LEN;
  int                 mpirank;
  int                 conti;
  double              scale;
  const char         *filename;
  const double       *v;
  const p4est_topidx_t *tree_to_vertex;
  p4est_topidx_t      first_local_tree, last_local_tree;
  p4est_locidx_t      Ncells, Ncorners;
  p4est_t            *p4est;
  p4est_connectivity_t *connectivity;
  p4est_geometry_t   *geom;
#ifdef P4EST_VTK_ASCII
  double              wx, wy, wz;
#else
  int                 retval;
  uint8_t            *uint8_data;
  p4est_locidx_t     *locidx_data;
#endif
  int                 xi, yi, j, k;
#ifdef P4_TO_P8
  int                 zi;
#endif
  double              h2, eta_x, eta_y, eta_z = 0.;
  double              xyz[3], XYZ[3];   /* 3 not P4EST_DIM */
  size_t              num_quads, zz;
  p4est_topidx_t      jt;
  p4est_topidx_t      vt[P4EST_CHILDREN];
  p4est_locidx_t      quad_count, Npoints;
  p4est_locidx_t      sk, il, ntcid, *ntc;
  P4EST_VTK_FLOAT_TYPE *float_data;
  sc_array_t         *quadrants, *indeps;
  sc_array_t         *trees;
  p4est_tree_t       *tree;
  p4est_quadrant_t   *quad;
  p4est_nodes_t      *nodes;
  p4est_indep_t      *in;

  /* check a whole bunch of assertions, here and below */
  P4EST_ASSERT (cont != NULL);
  P4EST_ASSERT (!cont->writing);

  /* avoid uninitialized warning */
  for (k = 0; k < P4EST_CHILDREN; ++k) {
    vt[k] = -(k + 1);
  }

  /* from now on this context is officially in use for writing */
  cont->writing = 1;

  /* grab context variables */
  p4est = cont->p4est;
  filename = cont->filename;
  geom = cont->geom;
  scale = cont->scale;
  conti = cont->continuous;
  P4EST_ASSERT (filename != NULL);

  /* grab details from the forest */
  P4EST_ASSERT (p4est != NULL);
  mpirank = p4est->mpirank;
  connectivity = p4est->connectivity;
  P4EST_ASSERT (connectivity != NULL);
  v = connectivity->vertices;
  tree_to_vertex = connectivity->tree_to_vertex;
  if (geom == NULL) {
    SC_CHECK_ABORT (connectivity->num_vertices > 0,
                    "Must provide connectivity with vertex information");
    P4EST_ASSERT (v != NULL && tree_to_vertex != NULL);
  }
  trees = p4est->trees;
  first_local_tree = p4est->first_local_tree;
  last_local_tree = p4est->last_local_tree;
  Ncells = cont->Ncells = p4est->local_num_quadrants;

  cont->num_corners = Ncorners = P4EST_CHILDREN * Ncells;
  if (scale < 1. || !conti) {
    /* when we scale the quadrants we need each corner separately */
    cont->nodes = nodes = NULL;
    cont->Npoints = Npoints = Ncorners;
    cont->node_to_corner = ntc = NULL;
    indeps = NULL;
  }
  else {
    /* if scale == 1. and the point data is continuous,
     * we can reuse shared quadrant corners */
    cont->nodes = nodes = p4est_nodes_new (p4est, NULL);
    indeps = &nodes->indep_nodes;
    cont->Npoints = Npoints = nodes->num_owned_indeps;
    P4EST_ASSERT ((size_t) Npoints == indeps->elem_count);

    /* Establish a reverse lookup table from a node to its first reference.
     * It is slow to run twice through memory like this.  However, we also know
     * that writing data to disk is slower still, so we do not optimize.
     */
    cont->node_to_corner = ntc = P4EST_ALLOC (p4est_locidx_t, Npoints);
    memset (ntc, -1, Npoints * sizeof (p4est_locidx_t));
    for (sk = 0, il = 0; il < Ncells; ++il) {
      for (k = 0; k < P4EST_CHILDREN; ++sk, ++k) {
        ntcid = nodes->local_nodes[sk];
        P4EST_ASSERT (0 <= ntcid && ntcid < Npoints);
        if (ntc[ntcid] < 0) {
          ntc[ntcid] = sk;
        }
      }
    }
#ifdef P4EST_ENABLE_DEBUG
    /* the particular version of nodes we call makes sure they are tight */
    for (ntcid = 0; ntcid < Npoints; ++ntcid) {
      P4EST_ASSERT (0 <= ntc[ntcid] && ntc[ntcid] < Ncorners);
    }
#endif
  }

  /* Have each proc write to its own file */
  snprintf (cont->vtufilename, BUFSIZ, "%s_%04d.vtu", filename, mpirank);
  /* Use "w" for writing the initial part of the file.
   * For further parts, use "r+" and fseek so write_compressed succeeds.
   */
  cont->vtufile = fopen (cont->vtufilename, "wb");
  if (cont->vtufile == NULL) {
    P4EST_LERRORF ("Could not open %s for output\n", cont->vtufilename);
    p4est_vtk_context_destroy (cont);
    return NULL;
  }

  fprintf (cont->vtufile, "<?xml version=\"1.0\"?>\n");
  fprintf (cont->vtufile,
           "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\"");
#if defined P4EST_ENABLE_VTK_BINARY && defined P4EST_ENABLE_VTK_COMPRESSION
  fprintf (cont->vtufile, " compressor=\"vtkZLibDataCompressor\"");
#endif
#ifdef SC_IS_BIGENDIAN
  fprintf (cont->vtufile, " byte_order=\"BigEndian\">\n");
#else
  fprintf (cont->vtufile, " byte_order=\"LittleEndian\">\n");
#endif
  fprintf (cont->vtufile, "  <UnstructuredGrid>\n");
  fprintf (cont->vtufile,
           "    <Piece NumberOfPoints=\"%lld\" NumberOfCells=\"%lld\">\n",
           (long long) Npoints, (long long) Ncells);
  fprintf (cont->vtufile, "      <Points>\n");

  float_data = P4EST_ALLOC (P4EST_VTK_FLOAT_TYPE, 3 * Npoints);

  /* write point position data */
  fprintf (cont->vtufile, "        <DataArray type=\"%s\" Name=\"Position\""
           " NumberOfComponents=\"3\" format=\"%s\">\n",
           P4EST_VTK_FLOAT_NAME, P4EST_VTK_FORMAT_STRING);

  if (nodes == NULL) {
    /* loop over the trees */
    for (jt = first_local_tree, quad_count = 0; jt <= last_local_tree; ++jt) {
      tree = p4est_tree_array_index (trees, jt);
      quadrants = &tree->quadrants;
      num_quads = quadrants->elem_count;

      /* retrieve corners of the tree */
      if (geom == NULL) {
        for (k = 0; k < P4EST_CHILDREN; ++k) {
          vt[k] = tree_to_vertex[jt * P4EST_CHILDREN + k];
        }
      }
      else {
        /* provoke crash on logic bug */
        P4EST_ASSERT (vt[0] == -1);
        v = NULL;
      }

      /* loop over the elements in tree and calculate vertex coordinates */
      for (zz = 0; zz < num_quads; ++zz, ++quad_count) {
        quad = p4est_quadrant_array_index (quadrants, zz);
        h2 = .5 * intsize * P4EST_QUADRANT_LEN (quad->level);
        k = 0;
#ifdef P4_TO_P8
        for (zi = 0; zi < 2; ++zi) {
          eta_z = intsize * quad->z + h2 * (1. + (zi * 2 - 1) * scale);
#endif
          for (yi = 0; yi < 2; ++yi) {
            eta_y = intsize * quad->y + h2 * (1. + (yi * 2 - 1) * scale);
            for (xi = 0; xi < 2; ++xi) {
              P4EST_ASSERT (0 <= k && k < P4EST_CHILDREN);
              eta_x = intsize * quad->x + h2 * (1. + (xi * 2 - 1) * scale);
              if (geom != NULL) {
                xyz[0] = eta_x;
                xyz[1] = eta_y;
                xyz[2] = eta_z;
                geom->X (geom, jt, xyz, XYZ);
                for (j = 0; j < 3; ++j) {
                  float_data[3 * (P4EST_CHILDREN * quad_count + k) + j] =
                    (P4EST_VTK_FLOAT_TYPE) XYZ[j];
                }
              }
              else {
                for (j = 0; j < 3; ++j) {
                  /* *INDENT-OFF* */
                xyz[j] =
            ((1. - eta_z) * ((1. - eta_y) * ((1. - eta_x) * v[3 * vt[0] + j] +
                                                   eta_x  * v[3 * vt[1] + j]) +
                                   eta_y  * ((1. - eta_x) * v[3 * vt[2] + j] +
                                                   eta_x  * v[3 * vt[3] + j]))
#ifdef P4_TO_P8
             +     eta_z  * ((1. - eta_y) * ((1. - eta_x) * v[3 * vt[4] + j] +
                                                   eta_x  * v[3 * vt[5] + j]) +
                                   eta_y  * ((1. - eta_x) * v[3 * vt[6] + j] +
                                                   eta_x  * v[3 * vt[7] + j]))
#endif
            );
                  /* *INDENT-ON* */
                  float_data[3 * (P4EST_CHILDREN * quad_count + k) + j] =
                    (P4EST_VTK_FLOAT_TYPE) xyz[j];
                }
              }
              ++k;
            }
          }
#ifdef P4_TO_P8
        }
#endif
        P4EST_ASSERT (k == P4EST_CHILDREN);
      }
    }
    P4EST_ASSERT (P4EST_CHILDREN * quad_count == Npoints);
  }
  else {
    for (zz = 0; zz < indeps->elem_count; ++zz) {
      in = (p4est_indep_t *) sc_array_index (indeps, zz);

      /* retrieve corners of the tree */
      jt = in->p.which_tree;
      if (geom == NULL) {
        for (k = 0; k < P4EST_CHILDREN; ++k) {
          vt[k] = tree_to_vertex[jt * P4EST_CHILDREN + k];
        }
      }
      else {
        /* provoke crash on logic bug */
        P4EST_ASSERT (vt[0] == -1);
        v = NULL;
      }

      /* calculate vertex coordinates */
      eta_x = intsize * in->x;
      eta_y = intsize * in->y;
#ifdef P4_TO_P8
      eta_z = intsize * in->z;
#endif
      if (geom != NULL) {
        xyz[0] = eta_x;
        xyz[1] = eta_y;
        xyz[2] = eta_z;
        geom->X (geom, jt, xyz, XYZ);
        for (j = 0; j < 3; ++j) {
          float_data[3 * zz + j] = (P4EST_VTK_FLOAT_TYPE) XYZ[j];
        }
      }
      else {
        for (j = 0; j < 3; ++j) {
          /* *INDENT-OFF* */
          xyz[j] =
            ((1. - eta_z) * ((1. - eta_y) * ((1. - eta_x) * v[3 * vt[0] + j] +
                                                   eta_x  * v[3 * vt[1] + j]) +
                                   eta_y  * ((1. - eta_x) * v[3 * vt[2] + j] +
                                                   eta_x  * v[3 * vt[3] + j]))
  #ifdef P4_TO_P8
             +     eta_z  * ((1. - eta_y) * ((1. - eta_x) * v[3 * vt[4] + j] +
                                                   eta_x  * v[3 * vt[5] + j]) +
                                   eta_y  * ((1. - eta_x) * v[3 * vt[6] + j] +
                                                   eta_x  * v[3 * vt[7] + j]))
  #endif
            );
          /* *INDENT-ON* */

          float_data[3 * zz + j] = (P4EST_VTK_FLOAT_TYPE) xyz[j];
        }
      }
    }
  }

#ifdef P4EST_VTK_ASCII
  for (il = 0; il < Npoints; ++il) {
    wx = float_data[3 * il + 0];
    wy = float_data[3 * il + 1];
    wz = float_data[3 * il + 2];

    fprintf (cont->vtufile,
#ifdef P4EST_ENABLE_VTK_DOUBLES
             "     %24.16e %24.16e %24.16e\n",
#else
             "          %16.8e %16.8e %16.8e\n",
#endif
             wx, wy, wz);
  }
#else
  /* TODO: Don't allocate the full size of the array, only allocate
   * the chunk that will be passed to zlib and do this a chunk
   * at a time.
   */
  if (p4est_vtk_write_binary (cont->vtufile, (char *) float_data,
                              sizeof (*float_data) * 3 * Npoints)) {
    P4EST_LERROR (P4EST_STRING "_vtk: Error encoding points\n");
    P4EST_FREE (float_data);
    p4est_vtk_context_destroy (cont);
    return NULL;
  }
#endif
  P4EST_FREE (float_data);

  fprintf (cont->vtufile, "        </DataArray>\n");
  fprintf (cont->vtufile, "      </Points>\n");
  fprintf (cont->vtufile, "      <Cells>\n");

  /* write connectivity data */
  fprintf (cont->vtufile,
           "        <DataArray type=\"%s\" Name=\"connectivity\""
           " format=\"%s\">\n", P4EST_VTK_LOCIDX, P4EST_VTK_FORMAT_STRING);
#ifdef P4EST_VTK_ASCII
  for (sk = 0, il = 0; il < Ncells; ++il) {
    fprintf (cont->vtufile, "         ");
    for (k = 0; k < P4EST_CHILDREN; ++sk, ++k) {
      fprintf (cont->vtufile, " %lld", nodes == NULL ?
               (long long) sk : (long long) nodes->local_nodes[sk]);
    }
    fprintf (cont->vtufile, "\n");
  }
#else
  if (nodes == NULL) {
    locidx_data = P4EST_ALLOC (p4est_locidx_t, Ncorners);
    for (il = 0; il < Ncorners; ++il) {
      locidx_data[il] = il;
    }
    retval =
      p4est_vtk_write_binary (cont->vtufile, (char *) locidx_data,
                              sizeof (p4est_locidx_t) * Ncorners);
    P4EST_FREE (locidx_data);
  }
  else {
    retval =
      p4est_vtk_write_binary (cont->vtufile, (char *) nodes->local_nodes,
                              sizeof (p4est_locidx_t) * Ncorners);
  }
  if (retval) {
    P4EST_LERROR (P4EST_STRING "_vtk: Error encoding connectivity\n");
    p4est_vtk_context_destroy (cont);
    return NULL;
  }
#endif
  fprintf (cont->vtufile, "        </DataArray>\n");

  /* write offset data */
  fprintf (cont->vtufile, "        <DataArray type=\"%s\" Name=\"offsets\""
           " format=\"%s\">\n", P4EST_VTK_LOCIDX, P4EST_VTK_FORMAT_STRING);
#ifdef P4EST_VTK_ASCII
  fprintf (cont->vtufile, "         ");
  for (il = 1, sk = 1; il <= Ncells; ++il, ++sk) {
    fprintf (cont->vtufile, " %lld", (long long) (P4EST_CHILDREN * il));
    if (!(sk % 8) && il != Ncells)
      fprintf (cont->vtufile, "\n         ");
  }
  fprintf (cont->vtufile, "\n");
#else
  locidx_data = P4EST_ALLOC (p4est_locidx_t, Ncells);
  for (il = 1; il <= Ncells; ++il) {
    locidx_data[il - 1] = P4EST_CHILDREN * il;  /* same type */
  }
  if (p4est_vtk_write_binary (cont->vtufile, (char *) locidx_data,
                              sizeof (p4est_locidx_t) * Ncells)) {
    P4EST_LERROR (P4EST_STRING "_vtk: Error encoding offsets\n");
    P4EST_FREE (locidx_data);
    p4est_vtk_context_destroy (cont);
    return NULL;
  }
  P4EST_FREE (locidx_data);
#endif
  fprintf (cont->vtufile, "        </DataArray>\n");

  /* write type data */
  fprintf (cont->vtufile, "        <DataArray type=\"UInt8\" Name=\"types\""
           " format=\"%s\">\n", P4EST_VTK_FORMAT_STRING);
#ifdef P4EST_VTK_ASCII
  fprintf (cont->vtufile, "         ");
  for (il = 0, sk = 1; il < Ncells; ++il, ++sk) {
    fprintf (cont->vtufile, " %d", P4EST_VTK_CELL_TYPE);
    if (!(sk % 20) && il != (Ncells - 1))
      fprintf (cont->vtufile, "\n         ");
  }
  fprintf (cont->vtufile, "\n");
#else
  uint8_data = P4EST_ALLOC (uint8_t, Ncells);
  for (il = 0; il < Ncells; ++il) {
    uint8_data[il] = P4EST_VTK_CELL_TYPE;
  }
  if (p4est_vtk_write_binary (cont->vtufile, (char *) uint8_data,
                              sizeof (*uint8_data) * Ncells)) {
    P4EST_LERROR (P4EST_STRING "_vtk: Error encoding types\n");
    P4EST_FREE (uint8_data);
    p4est_vtk_context_destroy (cont);
    return NULL;
  }
  P4EST_FREE (uint8_data);
#endif
  fprintf (cont->vtufile, "        </DataArray>\n");
  fprintf (cont->vtufile, "      </Cells>\n");

  if (ferror (cont->vtufile)) {
    P4EST_LERROR (P4EST_STRING "_vtk: Error writing header\n");
    p4est_vtk_context_destroy (cont);
    return NULL;
  }

  /* Only have the root write to the parallel vtk file */
  if (mpirank == 0) {
    snprintf (cont->pvtufilename, BUFSIZ, "%s.pvtu", filename);

    cont->pvtufile = fopen (cont->pvtufilename, "wb");
    if (!cont->pvtufile) {
      P4EST_LERRORF ("Could not open %s for output\n", cont->pvtufilename);
      p4est_vtk_context_destroy (cont);
      return NULL;
    }

    fprintf (cont->pvtufile, "<?xml version=\"1.0\"?>\n");
    fprintf (cont->pvtufile,
             "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\"");
#if defined P4EST_ENABLE_VTK_BINARY && defined P4EST_ENABLE_VTK_COMPRESSION
    fprintf (cont->pvtufile, " compressor=\"vtkZLibDataCompressor\"");
#endif
#ifdef SC_IS_BIGENDIAN
    fprintf (cont->pvtufile, " byte_order=\"BigEndian\">\n");
#else
    fprintf (cont->pvtufile, " byte_order=\"LittleEndian\">\n");
#endif

    fprintf (cont->pvtufile, "  <PUnstructuredGrid GhostLevel=\"0\">\n");
    fprintf (cont->pvtufile, "    <PPoints>\n");
    fprintf (cont->pvtufile, "      <PDataArray type=\"%s\" Name=\"Position\""
             " NumberOfComponents=\"3\" format=\"%s\"/>\n",
             P4EST_VTK_FLOAT_NAME, P4EST_VTK_FORMAT_STRING);
    fprintf (cont->pvtufile, "    </PPoints>\n");

    if (ferror (cont->pvtufile)) {
      P4EST_LERROR (P4EST_STRING "_vtk: Error writing parallel header\n");
      p4est_vtk_context_destroy (cont);
      return NULL;
    }

    /* Create a master file for visualization in Visit; this will be used
     * only in p4est_vtk_write_footer().
     */
    snprintf (cont->visitfilename, BUFSIZ, "%s.visit", filename);
    cont->visitfile = fopen (cont->visitfilename, "wb");
    if (!cont->visitfile) {
      P4EST_LERRORF ("Could not open %s for output\n", cont->visitfilename);
      p4est_vtk_context_destroy (cont);
      return NULL;
    }
  }

  /* the nodes object is no longer needed */
  if (nodes != NULL) {
    p4est_nodes_destroy (cont->nodes);
    cont->nodes = NULL;
  }

  return cont;
}

#ifdef P4_TO_P8
/* Based on
 * https://github.com/Kitware/VTK/blob/99770c75c2df471c456323d66a4a0bd154cf3a82
 * /Common/DataModel/vtkHigherOrderHexahedron.cxx#L611
 */
static int
point_index_from_ijk (int i, int j, int k, const int *order)
{
  int                 ibdy = (i == 0 || i == order[0]);
  int                 jbdy = (j == 0 || j == order[1]);
  int                 kbdy = (k == 0 || k == order[2]);
  /* How many boundaries do we lie on at once? */
  int                 nbdy = (ibdy ? 1 : 0) + (jbdy ? 1 : 0) + (kbdy ? 1 : 0);

  if (nbdy == 3) {              /* Vertex DOF */
    /* ijk is a corner node. Return the proper index (somewhere in [0,7]): */
    return (i ? (j ? 2 : 1) : (j ? 3 : 0)) + (k ? 4 : 0);
  }

  int                 offset = 8;
  if (nbdy == 2) {              /* Edge DOF */
    if (!ibdy) {                /* On i axis */
      return (i - 1) + (j ? order[0] + order[1] - 2 : 0) +
        (k ? 2 * (order[0] + order[1] - 2) : 0) + offset;
    }
    if (!jbdy) {                /* On j axis */
      return (j - 1) + (i ? order[0] - 1 : 2 * (order[0] - 1) + order[1] -
                        1) + (k ? 2 * (order[0] + order[1] - 2) : 0) + offset;
    }
    /* !kbdy, On k axis */
    offset += 4 * (order[0] - 1) + 4 * (order[1] - 1);
    return (k - 1) + (order[2] - 1) * (i ? (j ? 3 : 1) : (j ? 2 : 0)) +
      offset;
  }

  offset += 4 * (order[0] + order[1] + order[2] - 3);
  if (nbdy == 1) {              /* Face DOF */
    if (ibdy) {                 /* On i-normal face */
      return (j - 1) + ((order[1] - 1) * (k - 1)) +
        (i ? (order[1] - 1) * (order[2] - 1) : 0) + offset;
    }
    offset += 2 * (order[1] - 1) * (order[2] - 1);
    if (jbdy) {                 /* On j-normal face */
      return (i - 1) + ((order[0] - 1) * (k - 1)) +
        (j ? (order[2] - 1) * (order[0] - 1) : 0) + offset;
    }
    offset += 2 * (order[2] - 1) * (order[0] - 1);
    /* kbdy, On k-normal face */
    return (i - 1) + ((order[0] - 1) * (j - 1)) +
      (k ? (order[0] - 1) * (order[1] - 1) : 0) + offset;
  }

  /* nbdy == 0: Body DOF */
  offset += 2 *
    ((order[1] - 1) * (order[2] - 1) + (order[2] - 1) * (order[0] - 1) +
     (order[0] - 1) * (order[1] - 1));
  return offset + (i - 1) + (order[0] - 1) *
    ((j - 1) + (order[1] - 1) * ((k - 1)));
}
#else
/* Based on
 * https://github.com/Kitware/VTK/blob/99770c75c2df471c456323d66a4a0bd154cf3a82
 * /Common/DataModel/vtkHigherOrderQuadrilateral.cxx#L446
 */
static int
point_index_from_ijk (int i, int j, const int *order)
{
  int                 ibdy = (i == 0 || i == order[0]);
  int                 jbdy = (j == 0 || j == order[1]);
  /* How many boundaries do we lie on at once? */
  int                 nbdy = (ibdy ? 1 : 0) + (jbdy ? 1 : 0);

  if (nbdy == 2) {              /* Vertex DOF */
    /* ijk is a corner node. Return the proper index (somewhere in [0,7]): */
    return (i ? (j ? 2 : 1) : (j ? 3 : 0));
  }

  int                 offset = 4;
  if (nbdy == 1) {              /* Edge DOF */
    if (!ibdy) {                /* On i axis */
      return (i - 1) + (j ? order[0] - 1 + order[1] - 1 : 0) + offset;
    }
    if (!jbdy) {                /* On j axis */
      return (j - 1) + (i ? order[0] - 1 : 2 * (order[0] - 1) + order[1] -
                        1) + offset;
    }
  }

  offset += 2 * (order[0] - 1 + order[1] - 1);
  /* nbdy == 0: Face DOF */
  return offset + (i - 1) + (order[0] - 1) * ((j - 1));
}
#endif

p4est_vtk_context_t *
p4est_vtk_write_header_ho (p4est_vtk_context_t * cont, sc_array_t * positions,
                           int Nnodes1D)
{
  /* positions in order of: x,y,z,x,y,z */
  int                 mpirank;
  const char         *filename;
  p4est_locidx_t      Ncells;
  p4est_t            *p4est;
  p4est_locidx_t     *locidx_data;
#ifdef P4EST_VTK_ASCII
  double              wx, wy, wz = 0.;
#else
  P4EST_VTK_FLOAT_TYPE *float_data;
  uint8_t            *uint8_data;
#endif
  int                 i, j;
#if defined P4_TO_P8 || defined P4EST_VTK_ASCII
  int                 k;
#endif
  int                 order[P4EST_DIM];
  p4est_locidx_t      Npoints, Npointscell;
  p4est_locidx_t      sk, il;

  /* check a whole bunch of assertions, here and below */
  P4EST_ASSERT (cont != NULL);
  P4EST_ASSERT (!cont->writing);

  /* from now on this context is officially in use for writing */
  cont->writing = 1;

  /* grab context variables */
  p4est = cont->p4est;
  filename = cont->filename;
  P4EST_ASSERT (filename != NULL);

  /* grab details from the forest */
  P4EST_ASSERT (p4est != NULL);
  mpirank = p4est->mpirank;
  Ncells = cont->Ncells = p4est->local_num_quadrants;

  cont->num_corners = P4EST_CHILDREN * Ncells;
  cont->nodes = NULL;
#ifdef P4_TO_P8
  Npointscell = Nnodes1D * Nnodes1D * Nnodes1D;
#else
  Npointscell = Nnodes1D * Nnodes1D;
#endif
  cont->Npoints = Npoints = Npointscell * Ncells;
  cont->node_to_corner = NULL;

  /* Have each proc write to its own file */
  snprintf (cont->vtufilename, BUFSIZ, "%s_%04d.vtu", filename, mpirank);
  /* Use "w" for writing the initial part of the file.
   * For further parts, use "r+" and fseek so write_compressed succeeds.
   */
  cont->vtufile = fopen (cont->vtufilename, "wb");
  if (cont->vtufile == NULL) {
    P4EST_LERRORF ("Could not open %s for output\n", cont->vtufilename);
    p4est_vtk_context_destroy (cont);
    return NULL;
  }

  fprintf (cont->vtufile, "<?xml version=\"1.0\"?>\n");
  fprintf (cont->vtufile,
           "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\"");
#if defined P4EST_ENABLE_VTK_BINARY && defined P4EST_ENABLE_VTK_COMPRESSION
  fprintf (cont->vtufile, " compressor=\"vtkZLibDataCompressor\"");
#endif
#ifdef SC_IS_BIGENDIAN
  fprintf (cont->vtufile, " byte_order=\"BigEndian\">\n");
#else
  fprintf (cont->vtufile, " byte_order=\"LittleEndian\">\n");
#endif
  fprintf (cont->vtufile, "  <UnstructuredGrid>\n");
  fprintf (cont->vtufile,
           "    <Piece NumberOfPoints=\"%lld\" NumberOfCells=\"%lld\">\n",
           (long long) Npoints, (long long) Ncells);
  fprintf (cont->vtufile, "      <Points>\n");

  /* write point position data */
  fprintf (cont->vtufile, "        <DataArray type=\"%s\" Name=\"Position\""
           " NumberOfComponents=\"3\" format=\"%s\">\n",
           P4EST_VTK_FLOAT_NAME, P4EST_VTK_FORMAT_STRING);

#ifdef P4EST_VTK_ASCII
  for (il = 0; il < Npoints; ++il) {
    wx = *(double *) sc_array_index (positions, (il * P4EST_DIM));
    wy = *(double *) sc_array_index (positions, (il * P4EST_DIM) + 1);
#ifdef P4_TO_P8
    wz = *(double *) sc_array_index (positions, (il * P4EST_DIM) + 2);
#endif

    fprintf (cont->vtufile,
#ifdef P4EST_ENABLE_VTK_DOUBLES
             "     %24.16e %24.16e %24.16e\n",
#else
             "          %16.8e %16.8e %16.8e\n",
#endif
             wx, wy, wz);
  }
#else
  float_data = P4EST_ALLOC (P4EST_VTK_FLOAT_TYPE, 3 * Npoints);
  /* TODO: Don't allocate the full size of the array, only allocate
   * the chunk that will be passed to zlib and do this a chunk
   * at a time.
   */
  for (il = 0; il < Npoints; ++il) {
    for (j = 0; j < P4EST_DIM; ++j) {
      float_data[il * 3 + j] = (P4EST_VTK_FLOAT_TYPE) *
        ((double *) sc_array_index (positions, (il * P4EST_DIM) + j));
    }
#ifndef P4_TO_P8
    float_data[il * 3 + 2] = 0.;
#endif
  }
  if (p4est_vtk_write_binary (cont->vtufile, (char *) float_data,
                              sizeof (*float_data) * 3 * Npoints)) {
    P4EST_LERROR (P4EST_STRING "_vtk: Error encoding points\n");
    P4EST_FREE (float_data);
    p4est_vtk_context_destroy (cont);
    return NULL;
  }
  P4EST_FREE (float_data);
#endif

  fprintf (cont->vtufile, "        </DataArray>\n");
  fprintf (cont->vtufile, "      </Points>\n");
  fprintf (cont->vtufile, "      <Cells>\n");

  /* write connectivity data */
  fprintf (cont->vtufile,
           "        <DataArray type=\"%s\" Name=\"connectivity\""
           " format=\"%s\">\n", P4EST_VTK_LOCIDX, P4EST_VTK_FORMAT_STRING);
  locidx_data = P4EST_ALLOC (p4est_locidx_t, Npoints);
  order[0] = order[1] = Nnodes1D - 1;
#ifdef P4_TO_P8
  order[2] = order[0];
#endif
  for (sk = 0, il = 0; il < Ncells; ++il) {
#ifdef P4_TO_P8
    for (k = 0; k < Nnodes1D; ++k) {
#endif
      for (j = 0; j < Nnodes1D; ++j) {
        for (i = 0; i < Nnodes1D; ++sk, ++i) {
          locidx_data[il * Npointscell +
#ifdef P4_TO_P8
                      point_index_from_ijk (i, j, k, order)
#else
                      point_index_from_ijk (i, j, order)
#endif
            ] = sk;
        }
      }
#ifdef P4_TO_P8
    }
#endif
  }
#ifdef P4EST_VTK_ASCII
  for (sk = 0, il = 0; il < Ncells; ++il) {
    fprintf (cont->vtufile, "         ");
    for (k = 0; k < Npointscell; ++sk, ++k) {
      fprintf (cont->vtufile, " %lld", (long long) locidx_data[sk]);
    }
    fprintf (cont->vtufile, "\n");
  }
#else
  if (p4est_vtk_write_binary (cont->vtufile, (char *) locidx_data,
                              sizeof (p4est_locidx_t) * Npoints)) {
    P4EST_LERROR (P4EST_STRING "_vtk: Error encoding connectivity\n");
    P4EST_FREE (locidx_data);
    p4est_vtk_context_destroy (cont);
    return NULL;
  }
#endif
  P4EST_FREE (locidx_data);
  fprintf (cont->vtufile, "        </DataArray>\n");

  /* write offset data */
  fprintf (cont->vtufile, "        <DataArray type=\"%s\" Name=\"offsets\""
           " format=\"%s\">\n", P4EST_VTK_LOCIDX, P4EST_VTK_FORMAT_STRING);
#ifdef P4EST_VTK_ASCII
  fprintf (cont->vtufile, "         ");
  for (il = 1, sk = 1; il <= Ncells; ++il, ++sk) {
    fprintf (cont->vtufile, " %lld", (long long) (Npointscell * il));
    if (!(sk % 8) && il != Ncells)
      fprintf (cont->vtufile, "\n         ");
  }
  fprintf (cont->vtufile, "\n");
#else
  locidx_data = P4EST_ALLOC (p4est_locidx_t, Ncells);
  for (il = 1; il <= Ncells; ++il) {
    locidx_data[il - 1] = Npointscell * il;     /* same type */
  }
  if (p4est_vtk_write_binary (cont->vtufile, (char *) locidx_data,
                              sizeof (p4est_locidx_t) * Ncells)) {
    P4EST_LERROR (P4EST_STRING "_vtk: Error encoding offsets\n");
    P4EST_FREE (locidx_data);
    p4est_vtk_context_destroy (cont);
    return NULL;
  }
  P4EST_FREE (locidx_data);
#endif
  fprintf (cont->vtufile, "        </DataArray>\n");

  /* write type data */
  fprintf (cont->vtufile, "        <DataArray type=\"UInt8\" Name=\"types\""
           " format=\"%s\">\n", P4EST_VTK_FORMAT_STRING);
#ifdef P4EST_VTK_ASCII
  fprintf (cont->vtufile, "         ");
  for (il = 0, sk = 1; il < Ncells; ++il, ++sk) {
    fprintf (cont->vtufile, " %d", P4EST_VTK_CELL_TYPE_HO);
    if (!(sk % 20) && il != (Ncells - 1))
      fprintf (cont->vtufile, "\n         ");
  }
  fprintf (cont->vtufile, "\n");
#else
  uint8_data = P4EST_ALLOC (uint8_t, Ncells);
  for (il = 0; il < Ncells; ++il) {
    uint8_data[il] = P4EST_VTK_CELL_TYPE_HO;
  }
  if (p4est_vtk_write_binary (cont->vtufile, (char *) uint8_data,
                              sizeof (*uint8_data) * Ncells)) {
    P4EST_LERROR (P4EST_STRING "_vtk: Error encoding types\n");
    P4EST_FREE (uint8_data);
    p4est_vtk_context_destroy (cont);
    return NULL;
  }
  P4EST_FREE (uint8_data);
#endif
  fprintf (cont->vtufile, "        </DataArray>\n");
  fprintf (cont->vtufile, "      </Cells>\n");

  if (ferror (cont->vtufile)) {
    P4EST_LERROR (P4EST_STRING "_vtk: Error writing header\n");
    p4est_vtk_context_destroy (cont);
    return NULL;
  }

  /* Only have the root write to the parallel vtk file */
  if (mpirank == 0) {
    snprintf (cont->pvtufilename, BUFSIZ, "%s.pvtu", filename);

    cont->pvtufile = fopen (cont->pvtufilename, "wb");
    if (!cont->pvtufile) {
      P4EST_LERRORF ("Could not open %s for output\n", cont->pvtufilename);
      p4est_vtk_context_destroy (cont);
      return NULL;
    }

    fprintf (cont->pvtufile, "<?xml version=\"1.0\"?>\n");
    fprintf (cont->pvtufile,
             "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\"");
#if defined P4EST_ENABLE_VTK_BINARY && defined P4EST_ENABLE_VTK_COMPRESSION
    fprintf (cont->pvtufile, " compressor=\"vtkZLibDataCompressor\"");
#endif
#ifdef SC_IS_BIGENDIAN
    fprintf (cont->pvtufile, " byte_order=\"BigEndian\">\n");
#else
    fprintf (cont->pvtufile, " byte_order=\"LittleEndian\">\n");
#endif

    fprintf (cont->pvtufile, "  <PUnstructuredGrid GhostLevel=\"0\">\n");
    fprintf (cont->pvtufile, "    <PPoints>\n");
    fprintf (cont->pvtufile, "      <PDataArray type=\"%s\" Name=\"Position\""
             " NumberOfComponents=\"3\" format=\"%s\"/>\n",
             P4EST_VTK_FLOAT_NAME, P4EST_VTK_FORMAT_STRING);
    fprintf (cont->pvtufile, "    </PPoints>\n");

    if (ferror (cont->pvtufile)) {
      P4EST_LERROR (P4EST_STRING "_vtk: Error writing parallel header\n");
      p4est_vtk_context_destroy (cont);
      return NULL;
    }

    /* Create a master file for visualization in Visit; this will be used
     * only in p4est_vtk_write_footer().
     */
    snprintf (cont->visitfilename, BUFSIZ, "%s.visit", filename);
    cont->visitfile = fopen (cont->visitfilename, "wb");
    if (!cont->visitfile) {
      P4EST_LERRORF ("Could not open %s for output\n", cont->visitfilename);
      p4est_vtk_context_destroy (cont);
      return NULL;
    }
  }

  return cont;
}

/** Write VTK point data.
 *
 * This function exports custom point data to the vtk file; it is functionally
 * the same as \b p4est_vtk_write_point_dataf with the only difference being
 * that instead of a variable argument list, an initialized \a va_list is
 * passed as the last argument. The \a va_list is initialized from the variable
 * argument list of the calling function.
 *
 * \note This function is actually called from \b p4est_vtk_write_point_dataf
 *
 * \param [in,out] cont    A vtk context created by \ref p4est_vtk_context_new.
 * \param [in] num_point_scalars Number of point scalar datasets to output.
 * \param [in] num_point_vectors Number of point vector datasets to output.
 * \param [in,out] ap      An initialized va_list used to access the
 *                         scalar/vector data.
 *
 * \return          On success, the context that has been passed in.
 *                  On failure, returns NULL and deallocates the context.
 */
static p4est_vtk_context_t *
p4est_vtk_write_point_datav (p4est_vtk_context_t * cont,
                             int num_point_scalars,
                             int num_point_vectors, va_list ap)
{
  const int           num_point_all = num_point_scalars + num_point_vectors;
  int                 mpirank;
  int                 retval;
  int                 i, all;
  int                 scalar_strlen, vector_strlen;
  char                point_scalars[BUFSIZ], point_vectors[BUFSIZ];
  const char         *name, **names;
  p4est_vtk_context_t *list_end;
  sc_array_t        **values;
  p4est_locidx_t      ecount;

  P4EST_ASSERT (cont != NULL && cont->writing);
  P4EST_ASSERT (cont->p4est != NULL);

  /* This function needs to do nothing if there is no data. */
  if (!(num_point_scalars || num_point_vectors)) {
    return cont;
  }
  mpirank = cont->p4est->mpirank;

  /* Allocate storage to manage the data fields. */
  values = P4EST_ALLOC (sc_array_t *, num_point_all);
  names = P4EST_ALLOC (const char *, num_point_all);

  ecount =
    cont->node_to_corner == NULL ? cont->Npoints : cont->num_corners;

  /* Gather point data. */
  all = 0;
  scalar_strlen = 0;
  point_scalars[0] = '\0';
  for (i = 0; i < num_point_scalars; ++all, ++i) {
    name = names[all] = va_arg (ap, const char *);
    retval = snprintf (point_scalars + scalar_strlen, BUFSIZ - scalar_strlen,
                       "%s%s", i == 0 ? "" : ",", name);
    SC_CHECK_ABORT (retval > 0,
                    P4EST_STRING "_vtk: Error collecting point scalars");
    scalar_strlen += retval;
    values[all] = va_arg (ap, sc_array_t *);

    /* Validate input. */
    SC_CHECK_ABORT (values[all]->elem_size == sizeof (double),
                    P4EST_STRING
                    "_vtk: Error: incorrect point scalar data type;"
                    " scalar data must contain doubles.");
    SC_CHECK_ABORT (values[all]->elem_count == (size_t) ecount,
                    P4EST_STRING
                    "_vtk: Error: incorrect point scalar data count; see "
                    P4EST_STRING "_vtk.h for more details.");
  }

  /* keep variable all at current value */
  vector_strlen = 0;
  point_vectors[0] = '\0';
  for (i = 0; i < num_point_vectors; ++all, ++i) {
    name = names[all] = va_arg (ap, const char *);
    retval = snprintf (point_vectors + vector_strlen, BUFSIZ - vector_strlen,
                       "%s%s", i == 0 ? "" : ",", name);
    SC_CHECK_ABORT (retval > 0,
                    P4EST_STRING "_vtk: Error collecting point vectors");
    vector_strlen += retval;
    values[all] = va_arg (ap, sc_array_t *);

    /* Validate input. */
    SC_CHECK_ABORT (values[all]->elem_size == sizeof (double),
                    P4EST_STRING
                    "_vtk: Error: incorrect point vector data type;"
                    " vector data must contain doubles.");
    SC_CHECK_ABORT (values[all]->elem_count == 3 * (size_t) ecount,
                    P4EST_STRING
                    "_vtk: Error: incorrect point vector data count; see "
                    P4EST_STRING "_vtk.h for more details.");
  }

  /* Check for pointer variable marking the end of variable data input. */
  list_end = va_arg (ap, p4est_vtk_context_t *);
  SC_CHECK_ABORT (list_end == cont,
                  P4EST_STRING "_vtk Error: the end of variable data must be"
                  " specified by passing, as the last argument, the current "
                  P4EST_STRING "_vtk_context_t pointer.  See " P4EST_STRING
                  "_vtk.h for more information.");

  fprintf (cont->vtufile, "      <PointData");
  fprintf (cont->vtufile, " Scalars=\"%s\"", point_scalars);
  fprintf (cont->vtufile, " Vectors=\"%s\"", point_vectors);
  fprintf (cont->vtufile, ">\n");

  if (ferror (cont->vtufile)) {
    P4EST_LERRORF (P4EST_STRING "_vtk: Error writing %s\n",
                   cont->vtufilename);
    p4est_vtk_context_destroy (cont);

    P4EST_FREE (values);
    P4EST_FREE (names);

    return NULL;
  }

  /* now we are done counting and checking we write the fields */
  all = 0;
  for (i = 0; i < num_point_scalars; ++all, ++i) {
    cont = p4est_vtk_write_point_scalar (cont, names[all], values[all]);
    SC_CHECK_ABORT (cont != NULL,
                    P4EST_STRING "_vtk: Error writing point scalars");
  }

  for (i = 0; i < num_point_vectors; ++all, ++i) {
    cont = p4est_vtk_write_point_vector (cont, names[all], values[all]);
    SC_CHECK_ABORT (cont != NULL,
                    P4EST_STRING "_vtk: Error writing point vectors");
  }

  fprintf (cont->vtufile, "      </PointData>\n");

  P4EST_FREE (values);

  if (ferror (cont->vtufile)) {
    P4EST_LERRORF (P4EST_STRING "_vtk: Error writing %s\n",
                   cont->vtufilename);
    p4est_vtk_context_destroy (cont);

    P4EST_FREE (names);

    return NULL;
  }

  /* Only have the root write to the parallel vtk file */
  if (mpirank == 0) {
    fprintf (cont->pvtufile, "    <PPointData>\n");

    all = 0;
    for (i = 0; i < num_point_scalars; ++all, i++)
      fprintf (cont->pvtufile, "      "
               "<PDataArray type=\"%s\" Name=\"%s\" format=\"%s\"/>\n",
               P4EST_VTK_FLOAT_NAME, names[all], P4EST_VTK_FORMAT_STRING);

    for (i = 0; i < num_point_vectors; ++all, i++)
      fprintf (cont->pvtufile, "      "
               "<PDataArray type=\"%s\" Name=\"%s\" NumberOfComponents=\"3\" "
               "format=\"%s\"/>\n",
               P4EST_VTK_FLOAT_NAME, names[all], P4EST_VTK_FORMAT_STRING);

    fprintf (cont->pvtufile, "    </PPointData>\n");

    if (ferror (cont->pvtufile)) {
      P4EST_LERROR (P4EST_STRING "_vtk: Error writing parallel header\n");
      p4est_vtk_context_destroy (cont);

      P4EST_FREE (names);

      return NULL;
    }
  }

  P4EST_FREE (names);

  return cont;
}

p4est_vtk_context_t *
p4est_vtk_write_point_dataf (p4est_vtk_context_t * cont,
                             int num_point_scalars, int num_point_vectors,
                             ...)
{
  va_list             ap;

  P4EST_ASSERT (cont != NULL && cont->writing);
  P4EST_ASSERT (num_point_scalars >= 0 && num_point_vectors >= 0);

  va_start (ap, num_point_vectors);
  cont = p4est_vtk_write_point_datav (cont, num_point_scalars,
                                      num_point_vectors, ap);
  va_end (ap);

  return cont;
}

p4est_vtk_context_t *
p4est_vtk_write_cell_datav (p4est_vtk_context_t * cont,
                            int write_tree, int write_level,
                            int write_rank, int wrap_rank,
                            int num_cell_scalars,
                            int num_cell_vectors, va_list ap)
{
  int                 i;
  const char        **names;
  sc_array_t        **values;
  p4est_vtk_context_t *list_end = NULL, *retContext = NULL;

  P4EST_ASSERT (cont != NULL && cont->writing);
  P4EST_ASSERT (num_cell_scalars >= 0 && num_cell_vectors >= 0);
  P4EST_ASSERT (wrap_rank >= 0);

  names = P4EST_ALLOC (const char *, num_cell_scalars + num_cell_vectors);
  values = P4EST_ALLOC (sc_array_t *, num_cell_scalars + num_cell_vectors);

  for (i = 0; i < num_cell_scalars + num_cell_vectors; ++i) {
    names[i] = va_arg (ap, const char *);
    values[i] = va_arg (ap, sc_array_t *);
  }

  /* Check for pointer variable marking the end of variable data input. */
  list_end = va_arg (ap, p4est_vtk_context_t *);
  SC_CHECK_ABORT (list_end == cont,
                  P4EST_STRING "_vtk Error: the end of variable data must be"
                  " specified by passing, as the last argument, the current "
                  P4EST_STRING "_vtk_context_t pointer.  See " P4EST_STRING
                  "_vtk.h for more information.");

  retContext = p4est_vtk_write_cell_data (cont, write_tree, write_level,
                                          write_rank, wrap_rank,
                                          num_cell_scalars, num_cell_vectors,
                                          names, values);

  P4EST_FREE (values);
  P4EST_FREE (names);

  return retContext;
}

p4est_vtk_context_t *
p4est_vtk_write_cell_dataf (p4est_vtk_context_t * cont,
                            int write_tree, int write_level,
                            int write_rank, int wrap_rank,
                            int num_cell_scalars, int num_cell_vectors, ...)
{
  va_list             ap;

  P4EST_ASSERT (cont != NULL && cont->writing);
  P4EST_ASSERT (num_cell_scalars >= 0 && num_cell_vectors >= 0);
  P4EST_ASSERT (wrap_rank >= 0);

  va_start (ap, num_cell_vectors);
  cont = p4est_vtk_write_cell_datav (cont,
                                     write_tree, write_level,
                                     write_rank, wrap_rank,
                                     num_cell_scalars, num_cell_vectors, ap);
  va_end (ap);

  return cont;
}

p4est_vtk_context_t *
p4est_vtk_write_cell_data (p4est_vtk_context_t * cont,
                           int write_tree,
                           int write_level,
                           int write_rank,
                           int wrap_rank,
                           int num_cell_scalars,
                           int num_cell_vectors,
                           const char *fieldnames[], sc_array_t * values[])
{
  const int           mpirank = cont->p4est->mpirank;
  int                 retval;
  int                 i, all = 0;
  int                 scalar_strlen, vector_strlen;
  sc_array_t         *trees = cont->p4est->trees;
  p4est_tree_t       *tree;
  const p4est_topidx_t first_local_tree = cont->p4est->first_local_tree;
  const p4est_topidx_t last_local_tree = cont->p4est->last_local_tree;
  const p4est_locidx_t Ncells = cont->Ncells;
  char                cell_scalars[BUFSIZ], cell_vectors[BUFSIZ];
  const char         *name, **names;
  size_t              num_quads, zz;
  sc_array_t         *quadrants;
  p4est_quadrant_t   *quad;
#ifdef P4EST_VTK_ASCII
  p4est_locidx_t      sk;
#endif
  p4est_topidx_t      jt;
  p4est_locidx_t      il;

  char                vtkCellDataString[BUFSIZ] = "";
  int                 printed = 0;

  P4EST_ASSERT (cont != NULL && cont->writing);
  P4EST_ASSERT (num_cell_scalars >= 0 && num_cell_vectors >= 0);
  P4EST_ASSERT (wrap_rank >= 0);

  /* This function needs to do nothing if there is no data. */
  if (!(write_tree || write_level || write_rank ||
        num_cell_scalars || num_cell_vectors)) {
    return cont;
  }
  names = P4EST_ALLOC (const char *, num_cell_scalars + num_cell_vectors);

  /* Gather cell data. */
  scalar_strlen = 0;
  cell_scalars[0] = '\0';
  for (i = 0; i < num_cell_scalars; ++all, ++i) {
    name = names[all] = fieldnames[all];
    retval = snprintf (cell_scalars + scalar_strlen, BUFSIZ - scalar_strlen,
                       "%s%s", i == 0 ? "" : ",", name);
    SC_CHECK_ABORT (retval > 0,
                    P4EST_STRING "_vtk: Error collecting cell scalars");
    scalar_strlen += retval;

    /* Validate input. */
    SC_CHECK_ABORT (values[all]->elem_size == sizeof (double),
                    P4EST_STRING
                    "_vtk: Error: incorrect cell scalar data type;"
                    " scalar data must contain doubles.");
    SC_CHECK_ABORT (values[all]->elem_count ==
                    (size_t) cont->p4est->local_num_quadrants,
                    P4EST_STRING
                    "_vtk: Error: incorrect cell scalar data count;"
                    " scalar data must contain exactly"
                    " p4est->local_num_quadrants doubles.");
  }

  vector_strlen = 0;
  cell_vectors[0] = '\0';
  for (i = 0; i < num_cell_vectors; ++all, ++i) {
    name = names[all] = fieldnames[all];
    retval = snprintf (cell_vectors + vector_strlen, BUFSIZ - vector_strlen,
                       "%s%s", i == 0 ? "" : ",", name);
    SC_CHECK_ABORT (retval > 0,
                    P4EST_STRING "_vtk: Error collecting cell vectors");
    vector_strlen += retval;

    /* Validate input. */
    SC_CHECK_ABORT (values[all]->elem_size == sizeof (double),
                    P4EST_STRING
                    "_vtk: Error: incorrect cell vector data type;"
                    " vector data must contain doubles.");
    SC_CHECK_ABORT (values[all]->elem_count ==
                    3 * (size_t) cont->p4est->local_num_quadrants,
                    P4EST_STRING
                    "_vtk: Error: incorrect cell vector data count;"
                    " vector data must contain exactly"
                    " 3 * p4est->local_num_quadrants doubles.");
  }

  /* assemble cell scalar label string */
  if (write_tree)
    printed +=
      snprintf (vtkCellDataString + printed, BUFSIZ - printed, "treeid");
  if (write_level)
    printed +=
      snprintf (vtkCellDataString + printed, BUFSIZ - printed,
                printed > 0 ? ",level" : "level");
  if (write_rank)
    printed +=
      snprintf (vtkCellDataString + printed, BUFSIZ - printed,
                printed > 0 ? ",mpirank" : "mpirank");
  if (num_cell_scalars)
    printed +=
      snprintf (vtkCellDataString + printed, BUFSIZ - printed,
                printed > 0 ? ",%s" : "%s", cell_scalars);
  if (num_cell_vectors)
    printed +=
      snprintf (vtkCellDataString + printed, BUFSIZ - printed,
                printed > 0 ? ",%s" : "%s", cell_vectors);

  /* enter cell data context */
  fprintf (cont->vtufile, "      <CellData Scalars=\"%s\">\n",
           vtkCellDataString);

  if (write_tree) {
    p4est_locidx_t     *locidx_data = P4EST_ALLOC (p4est_locidx_t, Ncells);

    /* fill tree buffer */
    if (cont->tnodes == NULL) {
      P4EST_ASSERT (Ncells == cont->p4est->local_num_quadrants);
      for (il = 0, jt = first_local_tree; jt <= last_local_tree; ++jt) {
        tree = p4est_tree_array_index (trees, jt);
        num_quads = tree->quadrants.elem_count;
        for (zz = 0; zz < num_quads; ++zz) {
          locidx_data[il++] = (p4est_locidx_t) jt;
        }
      }
    }
    else {
      p4est_topidx_t     lftm;
      p4est_locidx_t     stoff;

      P4EST_ASSERT (cont->tnodes->local_tree_offset[0] == 0);
      for (il = 0, lftm = (jt = cont->tnodes->local_first_tree) - 1;
           jt <= cont->tnodes->local_last_tree; ++jt) {
        /* local simplices are stored in order of tree, then element */
        stoff = cont->tnodes->local_tree_offset[jt - lftm];
        while (il < stoff) {
          locidx_data[il++] = (p4est_locidx_t) jt;
        }
      }
    }
    P4EST_ASSERT (il == Ncells);

    /* write tree data */
    fprintf (cont->vtufile, "        <DataArray type=\"%s\" Name=\"treeid\""
             " format=\"%s\">\n", P4EST_VTK_LOCIDX, P4EST_VTK_FORMAT_STRING);
#ifdef P4EST_VTK_ASCII
    fprintf (cont->vtufile, "         ");
    for (il = 0, sk = 1; il < Ncells; ++il, ++sk) {
      fprintf (cont->vtufile, " %lld", (long long) locidx_data[il]);
      if (!(sk % 20) && il != (Ncells - 1)) {
        fprintf (cont->vtufile, "\n         ");
      }
    }
    fprintf (cont->vtufile, "\n");
#else
    if (p4est_vtk_write_binary (cont->vtufile, (char *) locidx_data,
                                sizeof (*locidx_data) * Ncells)) {
      P4EST_LERROR (P4EST_STRING "_vtk: Error encoding trees\n");
      P4EST_FREE (names);
      P4EST_FREE (locidx_data);
      p4est_vtk_context_destroy (cont);
      return NULL;
    }
#endif
    P4EST_FREE (locidx_data);
    fprintf (cont->vtufile, "        </DataArray>\n");
  }

  if (write_level) {
    int8_t             *int8_data = NULL;
    uint8_t            *uint8_data;

    /* fill level buffer */
    if (cont->tnodes == NULL) {
      P4EST_ASSERT (Ncells == cont->p4est->local_num_quadrants);
      int8_data = P4EST_ALLOC (int8_t, Ncells);
      for (il = 0, jt = first_local_tree; jt <= last_local_tree; ++jt) {
        tree = p4est_tree_array_index (trees, jt);
        quadrants = &tree->quadrants;
        num_quads = quadrants->elem_count;
        for (zz = 0; zz < num_quads; ++zz, ++il) {
          quad = p4est_quadrant_array_index (quadrants, zz);
          int8_data[il] = quad->level;
        }
      }
      P4EST_ASSERT (il == Ncells);

      /* casting int8_t to uint8_t is fine */
      uint8_data = (uint8_t *) int8_data;
    }
    else {
      if (cont->tnodes->simplex_level == NULL) {
        P4EST_LERROR (P4EST_STRING "_vtk: Simplices without levels\n");
        p4est_vtk_context_destroy (cont);
        P4EST_FREE (names);
        return NULL;
      }
      P4EST_ASSERT (cont->tnodes->simplex_level->elem_count ==
                    (size_t) Ncells);

      /* casting int8_t to uint8_t is fine */
      uint8_data = (uint8_t *) cont->tnodes->simplex_level->array;
    }

    /* write level data */
    fprintf (cont->vtufile, "        <DataArray type=\"%s\" Name=\"level\""
             " format=\"%s\">\n", "UInt8", P4EST_VTK_FORMAT_STRING);
#ifdef P4EST_VTK_ASCII
    fprintf (cont->vtufile, "         ");
    for (il = 0, sk = 1; il < Ncells; ++il, ++sk) {
      fprintf (cont->vtufile, " %d", (int) uint8_data[il]);
      if (!(sk % 20) && il != (Ncells - 1)) {
        fprintf (cont->vtufile, "\n         ");
      }
    }
    fprintf (cont->vtufile, "\n");
#else
    if (p4est_vtk_write_binary (cont->vtufile, (char *) uint8_data,
                                sizeof (*uint8_data) * Ncells)) {
      P4EST_LERROR (P4EST_STRING "_vtk: Error encoding levels\n");
      P4EST_FREE (names);
      P4EST_FREE (int8_data);
      p4est_vtk_context_destroy (cont);
      return NULL;
    }
#endif
    P4EST_FREE (int8_data);
    fprintf (cont->vtufile, "        </DataArray>\n");
  }

  if (write_rank) {
#ifndef P4EST_VTK_ASCII
    p4est_locidx_t     *locidx_data = P4EST_ALLOC (p4est_locidx_t, Ncells);
#endif
    const int           wrapped_rank =
      wrap_rank > 0 ? mpirank % wrap_rank : mpirank;

    fprintf (cont->vtufile, "        <DataArray type=\"%s\" Name=\"mpirank\""
             " format=\"%s\">\n", P4EST_VTK_LOCIDX, P4EST_VTK_FORMAT_STRING);
#ifdef P4EST_VTK_ASCII
    fprintf (cont->vtufile, "         ");
    for (il = 0, sk = 1; il < Ncells; ++il, ++sk) {
      fprintf (cont->vtufile, " %d", wrapped_rank);
      if (!(sk % 20) && il != (Ncells - 1))
        fprintf (cont->vtufile, "\n         ");
    }
    fprintf (cont->vtufile, "\n");
#else
    for (il = 0; il < Ncells; ++il) {
      locidx_data[il] = (p4est_locidx_t) wrapped_rank;
    }
    if (p4est_vtk_write_binary (cont->vtufile, (char *) locidx_data,
                                sizeof (*locidx_data) * Ncells)) {
      P4EST_LERROR (P4EST_STRING "_vtk: Error encoding types\n");
      P4EST_FREE (names);
      P4EST_FREE (locidx_data);
      p4est_vtk_context_destroy (cont);
      return NULL;
    }
    P4EST_FREE (locidx_data);
#endif
    fprintf (cont->vtufile, "        </DataArray>\n");
  }

  if (ferror (cont->vtufile)) {
    P4EST_LERRORF (P4EST_STRING "_vtk: Error writing %s\n",
                   cont->vtufilename);
    p4est_vtk_context_destroy (cont);

    P4EST_FREE (names);
    return NULL;
  }

  all = 0;
  for (i = 0; i < num_cell_scalars; ++all, ++i) {
    cont = p4est_vtk_write_cell_scalar (cont, names[all], values[all]);
    SC_CHECK_ABORT (cont != NULL,
                    P4EST_STRING "_vtk: Error writing cell scalars");
  }

  for (i = 0; i < num_cell_vectors; ++all, ++i) {
    cont = p4est_vtk_write_cell_vector (cont, names[all], values[all]);
    SC_CHECK_ABORT (cont != NULL,
                    P4EST_STRING "_vtk: Error writing cell vectors");
  }

  /* leave cell data context */
  fprintf (cont->vtufile, "      </CellData>\n");

  if (ferror (cont->vtufile)) {
    P4EST_LERRORF (P4EST_STRING "_vtk: Error writing %s\n",
                   cont->vtufilename);
    p4est_vtk_context_destroy (cont);
    P4EST_FREE (names);
    return NULL;
  }

  /* Only have the root write to the parallel vtk file */
  if (mpirank == 0) {
    fprintf (cont->pvtufile, "    <PCellData Scalars=\"%s\">\n",
             vtkCellDataString);

    if (write_tree)
      fprintf (cont->pvtufile, "      "
               "<PDataArray type=\"%s\" Name=\"treeid\" format=\"%s\"/>\n",
               P4EST_VTK_LOCIDX, P4EST_VTK_FORMAT_STRING);

    if (write_level)
      fprintf (cont->pvtufile, "      "
               "<PDataArray type=\"%s\" Name=\"level\" format=\"%s\"/>\n",
               "UInt8", P4EST_VTK_FORMAT_STRING);

    if (write_rank)
      fprintf (cont->pvtufile, "      "
               "<PDataArray type=\"%s\" Name=\"mpirank\" format=\"%s\"/>\n",
               P4EST_VTK_LOCIDX, P4EST_VTK_FORMAT_STRING);

    all = 0;
    for (i = 0; i < num_cell_scalars; ++all, i++)
      fprintf (cont->pvtufile, "      "
               "<PDataArray type=\"%s\" Name=\"%s\" format=\"%s\"/>\n",
               P4EST_VTK_FLOAT_NAME, names[all], P4EST_VTK_FORMAT_STRING);

    for (i = 0; i < num_cell_vectors; ++all, i++)
      fprintf (cont->pvtufile, "      "
               "<PDataArray type=\"%s\" NumberOfComponents=\"3\" Name=\"%s\" "
               "format=\"%s\"/>\n",
               P4EST_VTK_FLOAT_NAME, names[all], P4EST_VTK_FORMAT_STRING);

    fprintf (cont->pvtufile, "    </PCellData>\n");

    if (ferror (cont->pvtufile)) {
      P4EST_LERROR (P4EST_STRING "_vtk: Error writing parallel header\n");
      p4est_vtk_context_destroy (cont);
      P4EST_FREE (names);
      return NULL;
    }
  }

  P4EST_FREE (names);
  return cont;
}

/* Write either a scalar or a vector field.
 * If is_vector is true then a vector field with 3 values per point
 * is written. */
static p4est_vtk_context_t *
p4est_vtk_write_point (p4est_vtk_context_t * cont,
                       const char *field_name, sc_array_t * values,
                       int is_vector)
{
  p4est_locidx_t      il, ddl;
  int                 use_nodes;
#ifdef P4EST_ENABLE_DEBUG
  int                 ecount;
  int                 Nentries;
#endif
  int                 Npoints;
#ifndef P4EST_VTK_ASCII
  P4EST_VTK_FLOAT_TYPE *float_data;
#endif
  p4est_locidx_t     *ntc;

  P4EST_ASSERT (cont != NULL && cont->writing);
  Npoints = cont->Npoints;
  ntc = cont->node_to_corner;
#ifdef P4EST_ENABLE_DEBUG
  ecount = ntc == NULL ? Npoints : cont->num_corners;
  Nentries = ecount * (is_vector ? 3 : 1);
#endif

  P4EST_ASSERT (values != NULL && values->elem_count == (size_t) Nentries);

  if (ntc == NULL) {
    /* may be writing a discontinuous field, possibly due to vertex scaling */
    /* or we are writing higher order cells, relying on user coordinates */
    use_nodes = 0;
  }
  else {
    /* we are definitely writing a continuous field, reusing corner values */
    P4EST_ASSERT (cont->scale == 1. && cont->continuous);
    use_nodes = 1;
  }

  /* write point data */
  fprintf (cont->vtufile, "        <DataArray type=\"%s\" %s Name=\"%s\""
           " format=\"%s\">\n",
           P4EST_VTK_FLOAT_NAME,
           is_vector ? "NumberOfComponents=\"3\"" : "", field_name,
           P4EST_VTK_FORMAT_STRING);

#ifdef P4EST_VTK_ASCII
  if (!is_vector) {
    for (il = 0; il < Npoints; ++il) {
      ddl = use_nodes ? ntc[il] : il;
      P4EST_ASSERT (0 <= ddl && ddl < ecount);

      fprintf (cont->vtufile,
#ifdef P4EST_ENABLE_VTK_DOUBLES
               "     %24.16e\n",
#else
               "          %16.8e\n",
#endif
               *(double *) sc_array_index (values, ddl));
    }
  }
  else {
    for (il = 0; il < Npoints; ++il) {
      ddl = use_nodes ? ntc[il] : il;
      P4EST_ASSERT (0 <= ddl && ddl < ecount);

      fprintf (cont->vtufile,
#ifdef P4EST_ENABLE_VTK_DOUBLES
               "     %24.16e %24.16e %24.16e\n",
#else
               "          %16.8e %16.8e %16.8e\n",
#endif
               *(double *) sc_array_index (values, 3 * ddl),
               *(double *) sc_array_index (values, 3 * ddl + 1),
               *(double *) sc_array_index (values, 3 * ddl + 2));
    }
  }
#else
  if (!is_vector) {
    float_data = P4EST_ALLOC (P4EST_VTK_FLOAT_TYPE, Npoints);
    for (il = 0; il < Npoints; ++il) {
      ddl = use_nodes ? ntc[il] : il;
      P4EST_ASSERT (0 <= ddl && ddl < ecount);
      float_data[il] =
        (P4EST_VTK_FLOAT_TYPE) * ((double *) sc_array_index (values, ddl));
    }
  }
  else {
    float_data = P4EST_ALLOC (P4EST_VTK_FLOAT_TYPE, Npoints * 3);
    for (il = 0; il < Npoints; ++il) {
      ddl = use_nodes ? ntc[il] : il;
      P4EST_ASSERT (0 <= ddl && ddl < ecount);
      float_data[3 * il] =
        (P4EST_VTK_FLOAT_TYPE) *
        ((double *) sc_array_index (values, 3 * ddl));
      float_data[3 * il + 1] =
        (P4EST_VTK_FLOAT_TYPE) *
        ((double *) sc_array_index (values, 3 * ddl + 1));
      float_data[3 * il + 2] =
        (P4EST_VTK_FLOAT_TYPE) *
        ((double *) sc_array_index (values, 3 * ddl + 2));
    }
  }

  /* TODO: Don't allocate the full size of the array, only allocate
   * the chunk that will be passed to zlib and do this a chunk
   * at a time.
   */
  if (p4est_vtk_write_binary
      (cont->vtufile, (char *) float_data,
       sizeof (*float_data) * Npoints * (is_vector ? 3 : 1))) {
    P4EST_LERROR (P4EST_STRING "_vtk: Error encoding points\n");
    P4EST_FREE (float_data);
    p4est_vtk_context_destroy (cont);
    return NULL;
  }
  P4EST_FREE (float_data);
#endif
  fprintf (cont->vtufile, "        </DataArray>\n");

  if (ferror (cont->vtufile)) {
    P4EST_LERROR (P4EST_STRING "_vtk: Error writing point scalar\n");
    p4est_vtk_context_destroy (cont);
    return NULL;
  }

  return cont;
}

p4est_vtk_context_t *
p4est_vtk_write_point_scalar (p4est_vtk_context_t * cont,
                              const char *scalar_name, sc_array_t * values)
{
  return p4est_vtk_write_point (cont, scalar_name, values, 0);
}

p4est_vtk_context_t *
p4est_vtk_write_point_vector (p4est_vtk_context_t * cont,
                              const char *vector_name, sc_array_t * values)
{
  return p4est_vtk_write_point (cont, vector_name, values, 1);
}

/* Write either a scalar or a vector field.
 * If is_vector is true then a vector field with 3 values per cell
 * is written. */
static p4est_vtk_context_t *
p4est_vtk_write_cell (p4est_vtk_context_t * cont,
                      const char *field_name, sc_array_t * values,
                      int is_vector)
{
  const p4est_locidx_t Ncells = cont->p4est->local_num_quadrants;
  p4est_locidx_t      il;
#ifndef P4EST_VTK_ASCII
  P4EST_VTK_FLOAT_TYPE *float_data;
#endif

  P4EST_ASSERT (cont != NULL && cont->writing);

  /* Write cell data. */
  fprintf (cont->vtufile, "        <DataArray type=\"%s\" %s Name=\"%s\""
           " format=\"%s\">\n",
           P4EST_VTK_FLOAT_NAME, is_vector ? "NumberOfComponents=\"3\"" : "",
           field_name, P4EST_VTK_FORMAT_STRING);
#ifdef P4EST_VTK_ASCII
  if (!is_vector) {
    for (il = 0; il < Ncells; ++il) {
      fprintf (cont->vtufile,
#ifdef P4EST_ENABLE_VTK_DOUBLES
               "     %24.16e\n",
#else
               "          %16.8e\n",
#endif
               *(double *) sc_array_index (values, il));
    }
  }
  else {
    /* Write vector data */
    for (il = 0; il < Ncells; ++il) {
      fprintf (cont->vtufile,
#ifdef P4EST_ENABLE_VTK_DOUBLES
               "     %24.16e  %24.16e  %24.16e\n",
#else
               "          %16.8e  %16.8e  %16.8e\n",
#endif
               *(double *) sc_array_index (values, 3 * il),
               *(double *) sc_array_index (values, 3 * il + 1),
               *(double *) sc_array_index (values, 3 * il + 2));
    }
  }
#else
  if (!is_vector) {
    float_data = P4EST_ALLOC (P4EST_VTK_FLOAT_TYPE, Ncells);
    for (il = 0; il < Ncells; ++il) {
      float_data[il] =
        (P4EST_VTK_FLOAT_TYPE) * ((double *) sc_array_index (values, il));
    }
  }
  else {
    /* vector data */
    float_data = P4EST_ALLOC (P4EST_VTK_FLOAT_TYPE, 3 * Ncells);
    for (il = 0; il < Ncells; ++il) {
      float_data[3 * il] =
        (P4EST_VTK_FLOAT_TYPE) * ((double *) sc_array_index (values, 3 * il));
      float_data[3 * il + 1] =
        (P4EST_VTK_FLOAT_TYPE) *
        ((double *) sc_array_index (values, 3 * il + 1));
      float_data[3 * il + 2] =
        (P4EST_VTK_FLOAT_TYPE) *
        ((double *) sc_array_index (values, 3 * il + 2));
    }
  }

  /* TODO: Don't allocate the full size of the array, only allocate
   * the chunk that will be passed to zlib and do this a chunk
   * at a time.
   */
  if (p4est_vtk_write_binary
      (cont->vtufile, (char *) float_data,
       sizeof (*float_data) * Ncells * (is_vector ? 3 : 1))) {
    P4EST_LERROR (P4EST_STRING "_vtk: Error encoding scalar cell data\n");
    P4EST_FREE (float_data);
    p4est_vtk_context_destroy (cont);
    return NULL;
  }
  P4EST_FREE (float_data);
#endif
  fprintf (cont->vtufile, "        </DataArray>\n");

  if (ferror (cont->vtufile)) {
    P4EST_LERROR (P4EST_STRING "_vtk: Error writing cell scalar file\n");
    p4est_vtk_context_destroy (cont);
    return NULL;
  }

  return cont;
}

p4est_vtk_context_t *
p4est_vtk_write_cell_scalar (p4est_vtk_context_t * cont,
                             const char *scalar_name, sc_array_t * values)
{
  return p4est_vtk_write_cell (cont, scalar_name, values, 0);
}

p4est_vtk_context_t *
p4est_vtk_write_cell_vector (p4est_vtk_context_t * cont,
                             const char *vector_name, sc_array_t * values)
{
  return p4est_vtk_write_cell (cont, vector_name, values, 1);
}

int
p4est_vtk_write_footer (p4est_vtk_context_t * cont)
{
  int                 p;
  int                 procRank = cont->p4est->mpirank;
  int                 numProcs = cont->p4est->mpisize;
  char               *filename_basename, filename_cpy[BUFSIZ];

  P4EST_ASSERT (cont != NULL && cont->writing);

  fprintf (cont->vtufile, "    </Piece>\n");
  fprintf (cont->vtufile, "  </UnstructuredGrid>\n");
  fprintf (cont->vtufile, "</VTKFile>\n");

  if (ferror (cont->vtufile)) {
    P4EST_LERROR (P4EST_STRING "_vtk: Error writing footer\n");
    p4est_vtk_context_destroy (cont);
    return -1;
  }

  /* Only have the root write to the parallel vtk file */
  if (procRank == 0) {
    fprintf (cont->visitfile, "!NBLOCKS %d\n", numProcs);

    /* Write data about the parallel pieces into both files */
    for (p = 0; p < numProcs; ++p) {
      /* We want to write the basename of each processes file, since
       * the basename function could modify its argument, we create a
       * temporary copy. */
      snprintf (filename_cpy, BUFSIZ, "%s", cont->filename);
#ifdef _MSC_VER
      _splitpath (filename_cpy, NULL, NULL, NULL, NULL);
      filename_basename = filename_cpy;
#else
      filename_basename = basename (filename_cpy);
#endif
      fprintf (cont->pvtufile,
               "    <Piece Source=\"%s_%04d.vtu\"/>\n", filename_basename, p);
      fprintf (cont->visitfile, "%s_%04d.vtu\n", filename_basename, p);
    }
    fprintf (cont->pvtufile, "  </PUnstructuredGrid>\n");
    fprintf (cont->pvtufile, "</VTKFile>\n");

    if (ferror (cont->pvtufile)) {
      P4EST_LERROR (P4EST_STRING "_vtk: Error writing parallel footer\n");
      p4est_vtk_context_destroy (cont);
      return -1;
    }

    if (ferror (cont->visitfile)) {
      P4EST_LERROR (P4EST_STRING "_vtk: Error writing parallel footer\n");
      p4est_vtk_context_destroy (cont);
      return -1;
    }
  }

  /* Destroy context structure. */
  p4est_vtk_context_destroy (cont);

  return 0;
}
