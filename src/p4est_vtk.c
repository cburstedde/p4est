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

#ifdef P4_TO_P8
#include <p8est_vtk.h>
#include <p8est_nodes.h>
#define P4EST_VTK_CELL_TYPE     11      /* VTK_VOXEL */
#else
#include <p4est_vtk.h>
#include <p4est_nodes.h>
#define P4EST_VTK_CELL_TYPE      8      /* VTK_PIXEL */
#endif /* !P4_TO_P8 */

#include <sc_io.h>

static const double p4est_vtk_scale = 0.95;
static const int    p4est_vtk_write_tree = 1;
static const int    p4est_vtk_write_level = 1;
static const int    p4est_vtk_write_rank = 1;
static const int    p4est_vtk_wrap_rank = 0;

#ifndef P4EST_VTK_DOUBLES
#define P4EST_VTK_FLOAT_NAME "Float32"
#define P4EST_VTK_FLOAT_TYPE float
#else
#define P4EST_VTK_FLOAT_NAME "Float64"
#define P4EST_VTK_FLOAT_TYPE double
#endif

#ifndef P4EST_VTK_BINARY
#define P4EST_VTK_ASCII 1
#define P4EST_VTK_FORMAT_STRING "ascii"
#else
#define P4EST_VTK_FORMAT_STRING "binary"

static int
p4est_vtk_write_binary (FILE * vtkfile, char *numeric_data,
                        size_t byte_length)
{
#ifndef P4EST_VTK_COMPRESSION
  return sc_vtk_write_binary (vtkfile, numeric_data, byte_length);
#else
  return sc_vtk_write_compressed (vtkfile, numeric_data, byte_length);
#endif /* P4EST_VTK_COMPRESSION */
}

#endif /* P4EST_VTK_BINARY */

void
p4est_vtk_write_file (p4est_t * p4est, p4est_geometry_t * geom,
                      const char *filename)
{
  p4est_vtk_write_all (p4est, geom,
                       p4est_vtk_scale,
                       p4est_vtk_write_tree, p4est_vtk_write_level,
                       p4est_vtk_write_rank, p4est_vtk_wrap_rank,
                       0, 0, filename);
}

void
p4est_vtk_write_all (p4est_t * p4est, p4est_geometry_t * geom,
                     double scale,
                     int write_tree, int write_level,
                     int write_rank, int wrap_rank,
                     int num_scalars, int num_vectors,
                     const char *filename, ...)
{
  int                 retval;
  int                 i, all;
  int                 scalar_strlen, vector_strlen;
  char                point_scalars[BUFSIZ], point_vectors[BUFSIZ];
  const char         *name, **names;
  double            **values;
  va_list             ap;

  P4EST_ASSERT (num_scalars >= 0 && num_vectors >= 0);

  values = P4EST_ALLOC (double *, num_scalars + num_vectors);
  names = P4EST_ALLOC (const char *, num_scalars + num_vectors);

  va_start (ap, filename);
  all = 0;
  scalar_strlen = 0;
  point_scalars[0] = '\0';
  for (i = 0; i < num_scalars; ++all, ++i) {
    name = names[all] = va_arg (ap, const char *);
    retval = snprintf (point_scalars + scalar_strlen, BUFSIZ - scalar_strlen,
                       "%s%s", i == 0 ? "" : ",", name);
    SC_CHECK_ABORT (retval > 0,
                    P4EST_STRING "_vtk: Error collecting point scalars");
    scalar_strlen += retval;
    values[all] = va_arg (ap, double *);
  }
  vector_strlen = 0;
  point_vectors[0] = '\0';
  for (i = 0; i < num_vectors; ++all, ++i) {
    name = names[all] = va_arg (ap, const char *);
    retval = snprintf (point_vectors + vector_strlen, BUFSIZ - vector_strlen,
                       "%s%s", i == 0 ? "" : ",", name);
    SC_CHECK_ABORT (retval > 0,
                    P4EST_STRING "_vtk: Error collecting point vectors");
    vector_strlen += retval;
    values[all] = va_arg (ap, double *);
  }
  va_end (ap);

  retval = p4est_vtk_write_header (p4est, geom, scale,
                                   write_tree, write_level,
                                   write_rank, wrap_rank,
                                   num_scalars > 0 ? point_scalars : NULL,
                                   num_vectors > 0 ? point_vectors : NULL,
                                   filename);
  SC_CHECK_ABORT (!retval, P4EST_STRING "_vtk: Error writing header");

  all = 0;
  for (i = 0; i < num_scalars; ++all, ++i) {
    retval = p4est_vtk_write_point_scalar (p4est, geom, filename,
                                           names[all], values[all]);
    SC_CHECK_ABORT (!retval,
                    P4EST_STRING "_vtk: Error writing point scalars");
  }
  for (i = 0; i < num_vectors; ++all, ++i) {
    retval = p4est_vtk_write_point_vector (p4est, geom, filename,
                                           names[all], values[all]);
    SC_CHECK_ABORT (!retval,
                    P4EST_STRING "_vtk: Error writing point vectors");
  }

  retval = p4est_vtk_write_footer (p4est, filename);
  SC_CHECK_ABORT (!retval, P4EST_STRING "_vtk: Error writing footer");

  P4EST_FREE (values);
  P4EST_FREE (names);
}

int
p4est_vtk_write_header (p4est_t * p4est, p4est_geometry_t * geom,
                        double scale,
                        int write_tree, int write_level,
                        int write_rank, int wrap_rank,
                        const char *point_scalars, const char *point_vectors,
                        const char *filename)
{
  p4est_connectivity_t *connectivity = p4est->connectivity;
  sc_array_t         *trees = p4est->trees;
  const int           mpirank = p4est->mpirank;
  const double        intsize = 1.0 / P4EST_ROOT_LEN;
  const double       *v = connectivity->vertices;
  const p4est_topidx_t first_local_tree = p4est->first_local_tree;
  const p4est_topidx_t last_local_tree = p4est->last_local_tree;
  const p4est_topidx_t *tree_to_vertex = connectivity->tree_to_vertex;
  const p4est_locidx_t Ncells = p4est->local_num_quadrants;
  const p4est_locidx_t Ncorners = P4EST_CHILDREN * Ncells;
#ifdef P4EST_VTK_ASCII
  double              wx, wy, wz;
  p4est_locidx_t      sk;
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
  p4est_locidx_t      quad_count, Ntotal;
  p4est_locidx_t      il;
  P4EST_VTK_FLOAT_TYPE *float_data;
  sc_array_t         *quadrants, *indeps;
  p4est_tree_t       *tree;
  p4est_quadrant_t   *quad;
  p4est_nodes_t      *nodes;
  p4est_indep_t      *in;
  char                vtufilename[BUFSIZ];
  FILE               *vtufile;

  SC_CHECK_ABORT (p4est->connectivity->num_vertices > 0,
                  "Must provide connectivity with vertex information");

  P4EST_ASSERT (0. <= scale && scale <= 1. && wrap_rank >= 0);
  P4EST_ASSERT (v != NULL && tree_to_vertex != NULL);

  if (scale < 1.) {
    /* when we scale the quadrants we need each corner separately */
    nodes = NULL;
    indeps = NULL;
    Ntotal = Ncorners;
  }
  else {
    /* when scale == 1. we can reuse shared quadrant corners */
    nodes = p4est_nodes_new (p4est, NULL);
    indeps = &nodes->indep_nodes;
    Ntotal = nodes->num_owned_indeps;
    P4EST_ASSERT ((size_t) Ntotal == indeps->elem_count);
  }

  /* Have each proc write to its own file */
  snprintf (vtufilename, BUFSIZ, "%s_%04d.vtu", filename, mpirank);
  /* Use "w" for writing the initial part of the file.
   * For further parts, use "r+" and fseek so write_compressed succeeds.
   */
  vtufile = fopen (vtufilename, "wb");
  if (vtufile == NULL) {
    P4EST_LERRORF ("Could not open %s for output\n", vtufilename);
    return -1;
  }

  fprintf (vtufile, "<?xml version=\"1.0\"?>\n");
  fprintf (vtufile, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\"");
#if defined P4EST_VTK_BINARY && defined P4EST_VTK_COMPRESSION
  fprintf (vtufile, " compressor=\"vtkZLibDataCompressor\"");
#endif
#ifdef SC_IS_BIGENDIAN
  fprintf (vtufile, " byte_order=\"BigEndian\">\n");
#else
  fprintf (vtufile, " byte_order=\"LittleEndian\">\n");
#endif
  fprintf (vtufile, "  <UnstructuredGrid>\n");
  fprintf (vtufile,
           "    <Piece NumberOfPoints=\"%lld\" NumberOfCells=\"%lld\">\n",
           (long long) Ntotal, (long long) Ncells);
  fprintf (vtufile, "      <Points>\n");

  float_data = P4EST_ALLOC (P4EST_VTK_FLOAT_TYPE, 3 * Ntotal);

  /* write point position data */
  fprintf (vtufile, "        <DataArray type=\"%s\" Name=\"Position\""
           " NumberOfComponents=\"3\" format=\"%s\">\n",
           P4EST_VTK_FLOAT_NAME, P4EST_VTK_FORMAT_STRING);

  if (nodes == NULL) {
    /* loop over the trees */
    for (jt = first_local_tree, quad_count = 0; jt <= last_local_tree; ++jt) {
      tree = p4est_tree_array_index (trees, jt);
      quadrants = &tree->quadrants;
      num_quads = quadrants->elem_count;

      /* retrieve corners of the tree */
      for (k = 0; k < P4EST_CHILDREN; ++k)
        vt[k] = tree_to_vertex[jt * P4EST_CHILDREN + k];

      /* loop over the elements in tree and calculated vertex coordinates */
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
    P4EST_ASSERT (P4EST_CHILDREN * quad_count == Ntotal);
  }
  else {
    for (zz = 0; zz < indeps->elem_count; ++zz) {
      in = (p4est_indep_t *) sc_array_index (indeps, zz);

      /* retrieve corners of the tree */
      jt = in->p.which_tree;
      for (k = 0; k < P4EST_CHILDREN; ++k)
        vt[k] = tree_to_vertex[jt * P4EST_CHILDREN + k];

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
  for (il = 0; il < Ntotal; ++il) {
    wx = float_data[3 * il + 0];
    wy = float_data[3 * il + 1];
    wz = float_data[3 * il + 2];

#ifdef P4EST_VTK_DOUBLES
    fprintf (vtufile, "     %24.16e %24.16e %24.16e\n", wx, wy, wz);
#else
    fprintf (vtufile, "          %16.8e %16.8e %16.8e\n", wx, wy, wz);
#endif
  }
#else
  fprintf (vtufile, "          ");
  /* TODO: Don't allocate the full size of the array, only allocate
   * the chunk that will be passed to zlib and do this a chunk
   * at a time.
   */
  retval = p4est_vtk_write_binary (vtufile, (char *) float_data,
                                   sizeof (*float_data) * 3 * Ntotal);
  fprintf (vtufile, "\n");
  if (retval) {
    P4EST_LERROR (P4EST_STRING "_vtk: Error encoding points\n");
    fclose (vtufile);
    return -1;
  }
#endif
  P4EST_FREE (float_data);
  fprintf (vtufile, "        </DataArray>\n");
  fprintf (vtufile, "      </Points>\n");
  fprintf (vtufile, "      <Cells>\n");

  /* write connectivity data */
  fprintf (vtufile, "        <DataArray type=\"%s\" Name=\"connectivity\""
           " format=\"%s\">\n", P4EST_VTK_LOCIDX, P4EST_VTK_FORMAT_STRING);
#ifdef P4EST_VTK_ASCII
  for (sk = 0, il = 0; il < Ncells; ++il) {
    fprintf (vtufile, "         ");
    for (k = 0; k < P4EST_CHILDREN; ++sk, ++k) {
      fprintf (vtufile, " %lld", nodes == NULL ?
               (long long) sk : (long long) nodes->local_nodes[sk]);
    }
    fprintf (vtufile, "\n");
  }
#else
  locidx_data = P4EST_ALLOC (p4est_locidx_t, Ncorners);
  fprintf (vtufile, "          ");
  if (nodes == NULL) {
    for (il = 0; il < Ncorners; ++il) {
      locidx_data[il] = il;
    }
    retval =
      p4est_vtk_write_binary (vtufile, (char *) locidx_data,
                              sizeof (*locidx_data) * Ncorners);
  }
  else {
    retval =
      p4est_vtk_write_binary (vtufile, (char *) nodes->local_nodes,
                              sizeof (*locidx_data) * Ncorners);
  }
  fprintf (vtufile, "\n");
  if (retval) {
    P4EST_LERROR (P4EST_STRING "_vtk: Error encoding connectivity\n");
    fclose (vtufile);
    return -1;
  }
#endif
  fprintf (vtufile, "        </DataArray>\n");

  /* write offset data */
  fprintf (vtufile, "        <DataArray type=\"%s\" Name=\"offsets\""
           " format=\"%s\">\n", P4EST_VTK_LOCIDX, P4EST_VTK_FORMAT_STRING);
#ifdef P4EST_VTK_ASCII
  fprintf (vtufile, "         ");
  for (il = 1, sk = 1; il <= Ncells; ++il, ++sk) {
    fprintf (vtufile, " %lld", (long long) (P4EST_CHILDREN * il));
    if (!(sk % 8) && il != Ncells)
      fprintf (vtufile, "\n         ");
  }
  fprintf (vtufile, "\n");
#else
  for (il = 1; il <= Ncells; ++il)
    locidx_data[il - 1] = P4EST_CHILDREN * il;  /* same type */

  fprintf (vtufile, "          ");
  retval = p4est_vtk_write_binary (vtufile, (char *) locidx_data,
                                   sizeof (*locidx_data) * Ncells);
  fprintf (vtufile, "\n");
  if (retval) {
    P4EST_LERROR (P4EST_STRING "_vtk: Error encoding offsets\n");
    fclose (vtufile);
    return -1;
  }
#endif
  fprintf (vtufile, "        </DataArray>\n");

  /* write type data */
  fprintf (vtufile, "        <DataArray type=\"UInt8\" Name=\"types\""
           " format=\"%s\">\n", P4EST_VTK_FORMAT_STRING);
#ifdef P4EST_VTK_ASCII
  fprintf (vtufile, "         ");
  for (il = 0, sk = 1; il < Ncells; ++il, ++sk) {
    fprintf (vtufile, " %d", P4EST_VTK_CELL_TYPE);
    if (!(sk % 20) && il != (Ncells - 1))
      fprintf (vtufile, "\n         ");
  }
  fprintf (vtufile, "\n");
#else
  uint8_data = P4EST_ALLOC (uint8_t, Ncells);
  for (il = 0; il < Ncells; ++il)
    uint8_data[il] = P4EST_VTK_CELL_TYPE;

  fprintf (vtufile, "          ");
  retval = p4est_vtk_write_binary (vtufile, (char *) uint8_data,
                                   sizeof (*uint8_data) * Ncells);
  fprintf (vtufile, "\n");
  if (retval) {
    P4EST_LERROR (P4EST_STRING "_vtk: Error encoding types\n");
    fclose (vtufile);
    return -1;
  }
#endif
  fprintf (vtufile, "        </DataArray>\n");
  fprintf (vtufile, "      </Cells>\n");

  if (write_tree || write_level || write_rank) {
    char                vtkCellDataString[BUFSIZ] = "";
    int                 printed = 0;

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

    fprintf (vtufile, "      <CellData Scalars=\"%s\">\n", vtkCellDataString);
  }

  if (write_tree) {
    fprintf (vtufile, "        <DataArray type=\"%s\" Name=\"treeid\""
             " format=\"%s\">\n", P4EST_VTK_LOCIDX, P4EST_VTK_FORMAT_STRING);
#ifdef P4EST_VTK_ASCII
    fprintf (vtufile, "         ");
    for (il = 0, sk = 1, jt = first_local_tree; jt <= last_local_tree; ++jt) {
      tree = p4est_tree_array_index (trees, jt);
      num_quads = tree->quadrants.elem_count;
      for (zz = 0; zz < num_quads; ++zz, ++sk, ++il) {
        fprintf (vtufile, " %lld", (long long) jt);
        if (!(sk % 20) && il != (Ncells - 1))
          fprintf (vtufile, "\n         ");
      }
    }
    fprintf (vtufile, "\n");
#else
    for (il = 0, jt = first_local_tree; jt <= last_local_tree; ++jt) {
      tree = p4est_tree_array_index (trees, jt);
      num_quads = tree->quadrants.elem_count;
      for (zz = 0; zz < num_quads; ++zz, ++il) {
        locidx_data[il] = (p4est_locidx_t) jt;
      }
    }
    fprintf (vtufile, "          ");
    retval = p4est_vtk_write_binary (vtufile, (char *) locidx_data,
                                     sizeof (*locidx_data) * Ncells);
    fprintf (vtufile, "\n");
    if (retval) {
      P4EST_LERROR (P4EST_STRING "_vtk: Error encoding types\n");
      fclose (vtufile);
      return -1;
    }
#endif
    fprintf (vtufile, "        </DataArray>\n");
    P4EST_ASSERT (il == Ncells);
  }

  if (write_level) {
    fprintf (vtufile, "        <DataArray type=\"%s\" Name=\"level\""
             " format=\"%s\">\n", "UInt8", P4EST_VTK_FORMAT_STRING);
#ifdef P4EST_VTK_ASCII
    fprintf (vtufile, "         ");
    for (il = 0, sk = 1, jt = first_local_tree; jt <= last_local_tree; ++jt) {
      tree = p4est_tree_array_index (trees, jt);
      quadrants = &tree->quadrants;
      num_quads = quadrants->elem_count;
      for (zz = 0; zz < num_quads; ++zz, ++sk, ++il) {
        quad = p4est_quadrant_array_index (quadrants, zz);
        fprintf (vtufile, " %d", (int) quad->level);
        if (!(sk % 20) && il != (Ncells - 1))
          fprintf (vtufile, "\n         ");
      }
    }
    fprintf (vtufile, "\n");
#else
    for (il = 0, jt = first_local_tree; jt <= last_local_tree; ++jt) {
      tree = p4est_tree_array_index (trees, jt);
      quadrants = &tree->quadrants;
      num_quads = quadrants->elem_count;
      for (zz = 0; zz < num_quads; ++zz, ++il) {
        quad = p4est_quadrant_array_index (quadrants, zz);
        uint8_data[il] = (uint8_t) quad->level;
      }
    }

    fprintf (vtufile, "          ");
    retval = p4est_vtk_write_binary (vtufile, (char *) uint8_data,
                                     sizeof (*uint8_data) * Ncells);
    fprintf (vtufile, "\n");
    if (retval) {
      P4EST_LERROR (P4EST_STRING "_vtk: Error encoding types\n");
      fclose (vtufile);
      return -1;
    }
#endif
    fprintf (vtufile, "        </DataArray>\n");
  }

  if (write_rank) {
    const int           wrapped_rank =
      wrap_rank > 0 ? mpirank % wrap_rank : mpirank;

    fprintf (vtufile, "        <DataArray type=\"%s\" Name=\"mpirank\""
             " format=\"%s\">\n", P4EST_VTK_LOCIDX, P4EST_VTK_FORMAT_STRING);
#ifdef P4EST_VTK_ASCII
    fprintf (vtufile, "         ");
    for (il = 0, sk = 1; il < Ncells; ++il, ++sk) {
      fprintf (vtufile, " %d", wrapped_rank);
      if (!(sk % 20) && il != (Ncells - 1))
        fprintf (vtufile, "\n         ");
    }
    fprintf (vtufile, "\n");
#else
    for (il = 0; il < Ncells; ++il)
      locidx_data[il] = (p4est_locidx_t) wrapped_rank;

    fprintf (vtufile, "          ");
    retval = p4est_vtk_write_binary (vtufile, (char *) locidx_data,
                                     sizeof (*locidx_data) * Ncells);
    fprintf (vtufile, "\n");
    if (retval) {
      P4EST_LERROR (P4EST_STRING "_vtk: Error encoding types\n");
      fclose (vtufile);
      return -1;
    }
#endif
    fprintf (vtufile, "        </DataArray>\n");
  }

  if (write_tree || write_level || write_rank) {
    fprintf (vtufile, "      </CellData>\n");
  }

#ifndef P4EST_VTK_ASCII
  P4EST_FREE (locidx_data);
  P4EST_FREE (uint8_data);
#endif

  if (nodes != NULL) {
    p4est_nodes_destroy (nodes);
  }

  fprintf (vtufile, "      <PointData");
  if (point_scalars != NULL)
    fprintf (vtufile, " Scalars=\"%s\"", point_scalars);
  if (point_vectors != NULL)
    fprintf (vtufile, " Vectors=\"%s\"", point_vectors);
  fprintf (vtufile, ">\n");

  if (ferror (vtufile)) {
    P4EST_LERROR (P4EST_STRING "_vtk: Error writing header\n");
    fclose (vtufile);
    return -1;
  }
  if (fclose (vtufile)) {
    P4EST_LERROR (P4EST_STRING "_vtk: Error closing header\n");
    return -1;
  }
  vtufile = NULL;

  /* Only have the root write to the parallel vtk file */
  if (mpirank == 0) {
    char                pvtufilename[BUFSIZ];
    FILE               *pvtufile;

    snprintf (pvtufilename, BUFSIZ, "%s.pvtu", filename);

    pvtufile = fopen (pvtufilename, "wb");
    if (!pvtufile) {
      P4EST_LERRORF ("Could not open %s for output\n", vtufilename);
      return -1;
    }

    fprintf (pvtufile, "<?xml version=\"1.0\"?>\n");
    fprintf (pvtufile, "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\"");
#if defined P4EST_VTK_BINARY && defined P4EST_VTK_COMPRESSION
    fprintf (pvtufile, " compressor=\"vtkZLibDataCompressor\"");
#endif
#ifdef SC_IS_BIGENDIAN
    fprintf (pvtufile, " byte_order=\"BigEndian\">\n");
#else
    fprintf (pvtufile, " byte_order=\"LittleEndian\">\n");
#endif

    fprintf (pvtufile, "  <PUnstructuredGrid GhostLevel=\"0\">\n");
    fprintf (pvtufile, "    <PPoints>\n");
    fprintf (pvtufile, "      <PDataArray type=\"%s\" Name=\"Position\""
             " NumberOfComponents=\"3\" format=\"%s\"/>\n",
             P4EST_VTK_FLOAT_NAME, P4EST_VTK_FORMAT_STRING);
    fprintf (pvtufile, "    </PPoints>\n");
    if (write_tree || write_level || write_rank) {
      char                vtkCellDataString[BUFSIZ] = "";
      int                 printed = 0;

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

      fprintf (pvtufile, "    <PCellData Scalars=\"%s\">\n",
               vtkCellDataString);
    }
    if (write_tree) {
      fprintf (pvtufile, "      "
               "<PDataArray type=\"%s\" Name=\"treeid\" format=\"%s\"/>\n",
               P4EST_VTK_LOCIDX, P4EST_VTK_FORMAT_STRING);
    }
    if (write_level) {
      fprintf (pvtufile, "      "
               "<PDataArray type=\"%s\" Name=\"level\" format=\"%s\"/>\n",
               "UInt8", P4EST_VTK_FORMAT_STRING);
    }
    if (write_rank) {
      fprintf (pvtufile, "      "
               "<PDataArray type=\"%s\" Name=\"mpirank\" format=\"%s\"/>\n",
               P4EST_VTK_LOCIDX, P4EST_VTK_FORMAT_STRING);
    }
    if (write_tree || write_level || write_rank) {
      fprintf (pvtufile, "    </PCellData>\n");
    }
    fprintf (pvtufile, "    <PPointData>\n");

    if (ferror (pvtufile)) {
      P4EST_LERROR (P4EST_STRING "_vtk: Error writing parallel header\n");
      fclose (pvtufile);
      return -1;
    }
    if (fclose (pvtufile)) {
      P4EST_LERROR (P4EST_STRING "_vtk: Error closing parallel header\n");
      return -1;
    }
  }

  return 0;
}

int
p4est_vtk_write_point_scalar (p4est_t * p4est, p4est_geometry_t * geom,
                              const char *filename,
                              const char *scalar_name, const double *values)
{
  const int           mpirank = p4est->mpirank;
  const p4est_locidx_t Ncells = p4est->local_num_quadrants;
  const p4est_locidx_t Ncorners = P4EST_CHILDREN * Ncells;      /* type ok */
  int                 retval;
  p4est_locidx_t      il;
#ifndef P4EST_VTK_ASCII
  P4EST_VTK_FLOAT_TYPE *float_data;
#endif
  char                vtufilename[BUFSIZ];
  FILE               *vtufile;

  /* Have each proc write to its own file */
  snprintf (vtufilename, BUFSIZ, "%s_%04d.vtu", filename, mpirank);
  /* To be able to fseek in a file you cannot open in append mode.
   * so you need to open with "r+" and fseek to SEEK_END.
   */
  vtufile = fopen (vtufilename, "rb+");
  if (vtufile == NULL) {
    P4EST_LERRORF ("Could not open %s for output\n", vtufilename);
    return -1;
  }
  retval = fseek (vtufile, 0L, SEEK_END);
  if (retval) {
    P4EST_LERRORF ("Could not fseek %s for output\n", vtufilename);
    fclose (vtufile);
    return -1;
  }

  /* write point position data */
  fprintf (vtufile, "        <DataArray type=\"%s\" Name=\"%s\""
           " format=\"%s\">\n",
           P4EST_VTK_FLOAT_NAME, scalar_name, P4EST_VTK_FORMAT_STRING);

#ifdef P4EST_VTK_ASCII
  for (il = 0; il < Ncorners; ++il) {
#ifdef P4EST_VTK_DOUBLES
    fprintf (vtufile, "     %24.16e\n", values[il]);
#else
    fprintf (vtufile, "          %16.8e\n", values[il]);
#endif
  }
#else
  float_data = P4EST_ALLOC (P4EST_VTK_FLOAT_TYPE, Ncorners);
  for (il = 0; il < Ncorners; ++il) {
    float_data[il] = (P4EST_VTK_FLOAT_TYPE) values[il];
  }

  fprintf (vtufile, "          ");
  /* TODO: Don't allocate the full size of the array, only allocate
   * the chunk that will be passed to zlib and do this a chunk
   * at a time.
   */
  retval = p4est_vtk_write_binary (vtufile, (char *) float_data,
                                   sizeof (*float_data) * Ncorners);
  fprintf (vtufile, "\n");
  if (retval) {
    P4EST_LERROR (P4EST_STRING "_vtk: Error encoding points\n");
    fclose (vtufile);
    return -1;
  }
  P4EST_FREE (float_data);
#endif
  fprintf (vtufile, "        </DataArray>\n");

  if (ferror (vtufile)) {
    P4EST_LERROR (P4EST_STRING "_vtk: Error writing point scalar\n");
    fclose (vtufile);
    return -1;
  }
  if (fclose (vtufile)) {
    P4EST_LERROR (P4EST_STRING "_vtk: Error closing point scalar\n");
    return -1;
  }
  vtufile = NULL;

  /* Only have the root write to the parallel vtk file */
  if (mpirank == 0) {
    char                pvtufilename[BUFSIZ];
    FILE               *pvtufile;
    snprintf (pvtufilename, BUFSIZ, "%s.pvtu", filename);

    pvtufile = fopen (pvtufilename, "ab");
    if (!pvtufile) {
      P4EST_LERRORF ("Could not open %s for output\n", vtufilename);
      return -1;
    }

    fprintf (pvtufile, "      <PDataArray type=\"%s\" Name=\"%s\""
             " format=\"%s\"/>\n",
             P4EST_VTK_FLOAT_NAME, scalar_name, P4EST_VTK_FORMAT_STRING);

    if (ferror (pvtufile)) {
      P4EST_LERROR (P4EST_STRING
                    "_vtk: Error writing parallel point scalar\n");
      fclose (pvtufile);
      return -1;
    }
    if (fclose (pvtufile)) {
      P4EST_LERROR (P4EST_STRING
                    "_vtk: Error closing parallel point scalar\n");
      return -1;
    }
  }

  return 0;
}

int
p4est_vtk_write_point_vector (p4est_t * p4est, p4est_geometry_t * geom,
                              const char *filename,
                              const char *vector_name, const double *values)
{
  SC_ABORT (P4EST_STRING "_vtk_write_point_vector not implemented");
}

int
p4est_vtk_write_footer (p4est_t * p4est, const char *filename)
{
  char                vtufilename[BUFSIZ];
  int                 p;
  int                 procRank = p4est->mpirank;
  int                 numProcs = p4est->mpisize;
  FILE               *vtufile;

  /* Have each proc write to its own file */
  snprintf (vtufilename, BUFSIZ, "%s_%04d.vtu", filename, procRank);
  vtufile = fopen (vtufilename, "ab");
  if (vtufile == NULL) {
    P4EST_LERRORF ("Could not open %s for output!\n", vtufilename);
    return -1;
  }

  fprintf (vtufile, "      </PointData>\n");
  fprintf (vtufile, "    </Piece>\n");
  fprintf (vtufile, "  </UnstructuredGrid>\n");
  fprintf (vtufile, "</VTKFile>\n");

  if (ferror (vtufile)) {
    P4EST_LERROR ("p4est_vtk: Error writing footer\n");
    fclose (vtufile);
    return -1;
  }
  if (fclose (vtufile)) {
    P4EST_LERROR ("p4est_vtk: Error closing footer\n");
    return -1;
  }
  vtufile = NULL;

  /* Only have the root write to the parallel vtk file */
  if (procRank == 0) {
    char                visitfilename[BUFSIZ];
    char                pvtufilename[BUFSIZ];
    FILE               *pvtufile, *visitfile;

    /* Reopen paraview master file for writing bottom half */
    snprintf (pvtufilename, BUFSIZ, "%s.pvtu", filename);
    pvtufile = fopen (pvtufilename, "ab");
    if (!pvtufile) {
      P4EST_LERRORF ("Could not open %s for output!\n", vtufilename);
      return -1;
    }

    /* Create a master file for visualization in Visit */
    snprintf (visitfilename, BUFSIZ, "%s.visit", filename);
    visitfile = fopen (visitfilename, "wb");
    if (!visitfile) {
      P4EST_LERRORF ("Could not open %s for output\n", visitfilename);
      fclose (pvtufile);
      return -1;
    }
    fprintf (visitfile, "!NBLOCKS %d\n", numProcs);

    /* Write data about the parallel pieces into both files */
    fprintf (pvtufile, "    </PPointData>\n");
    for (p = 0; p < numProcs; ++p) {
      fprintf (pvtufile,
               "    <Piece Source=\"%s_%04d.vtu\"/>\n", filename, p);
      fprintf (visitfile, "%s_%04d.vtu\n", filename, p);
    }
    fprintf (pvtufile, "  </PUnstructuredGrid>\n");
    fprintf (pvtufile, "</VTKFile>\n");

    /* Close paraview master file */
    if (ferror (pvtufile)) {
      P4EST_LERROR ("p4est_vtk: Error writing parallel footer\n");
      fclose (visitfile);
      fclose (pvtufile);
      return -1;
    }
    if (fclose (pvtufile)) {
      fclose (visitfile);
      P4EST_LERROR ("p4est_vtk: Error closing parallel footer\n");
      return -1;
    }

    /* Close visit master file */
    if (ferror (visitfile)) {
      P4EST_LERROR ("p4est_vtk: Error writing parallel footer\n");
      fclose (visitfile);
      return -1;
    }
    if (fclose (visitfile)) {
      P4EST_LERROR ("p4est_vtk: Error closing parallel footer\n");
      return -1;
    }
  }

  return 0;
}
