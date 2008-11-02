/*
  This file is part of p4est.
  p4est is a C library to manage a parallel collection of quadtrees and/or
  octrees.

  Copyright (C) 2007,2008 Carsten Burstedde, Lucas Wilcox.

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

#ifdef P4_TO_P8
#include <p8est_vtk.h>
#else
#include <p4est_vtk.h>
#include <p4est_mesh.h>

bool                p4est_vtk_default_write_rank = true;
#endif /* !P4_TO_P8 */

#include <sc_io.h>

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

#ifndef P4_TO_P8

void
p4est_vtk_write_file (p4est_t * p4est, const char *baseName)
{
  int                 retval;

  SC_CHECK_ABORT (p4est->connectivity->vertices != NULL,
                  "Must provide connectivity with vertex information");

  retval =
    p4est_vtk_write_header (p4est, p4est_vtk_default_write_rank, baseName);
  SC_CHECK_ABORT (!retval, "VTK: write header");
  retval = p4est_vtk_write_footer (p4est, baseName);
  SC_CHECK_ABORT (!retval, "VTK: write footer");
}

int
p4est_vtk_write_header (p4est_t * p4est, bool write_rank,
                        const char *baseName)
{
  p4est_locidx_t      il, jl;
#ifdef P4EST_VTK_ASCII
  p4est_locidx_t      sk;
#else
  int                 retval;
  uint8_t            *uint8_data;
#endif
  p4est_locidx_t     *locidx_data;
  P4EST_VTK_FLOAT_TYPE *float_data;
  char                vtufilename[BUFSIZ];
  int                 rootRank = 0;
  int                 procRank = p4est->mpirank;
  FILE               *vtufile;
  p4est_connectivity_t *connectivity = p4est->connectivity;
  p4est_locidx_t      Ncells = p4est->local_num_quadrants;
  p4est_locidx_t      Ntotal = 0;
  p4est_qcoord_t      inth;
  double              h, eta1, eta2;
  double              v0x, v0y, v0z, v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y,
    v3z;
  double              w0x, w0y, w0z, w1x, w1y, w1z, w2x, w2y, w2z, w3x, w3y,
    w3z;
  p4est_topidx_t      v0, v1, v2, v3;
  p4est_locidx_t      lv0, lv1, lv2, lv3;
  p4est_topidx_t      first_local_tree = p4est->first_local_tree;
  p4est_topidx_t      last_local_tree = p4est->last_local_tree;
  sc_array_t         *trees = p4est->trees;
  sc_array_t         *quadrants;
  p4est_tree_t       *tree;
  p4est_locidx_t      num_quads;
  p4est_locidx_t      quad_count;
  p4est_topidx_t     *tree_to_vertex = connectivity->tree_to_vertex;
  double             *vertices = connectivity->vertices;
  p4est_quadrant_t   *quad;
  double              intsize = 1.0 / (1 << P4EST_MAXLEVEL);

  P4EST_ASSERT (connectivity->vertices != NULL);

  /* Have each proc write to its own file */
  snprintf (vtufilename, BUFSIZ, "%s_%04d.vtu", baseName, procRank);
  vtufile = fopen (vtufilename, "wb");
  if (vtufile == NULL) {
    P4EST_LERRORF ("Could not open %s for output!\n", vtufilename);
    return -1;
  }

  locidx_data = P4EST_ALLOC (p4est_locidx_t, 4 * Ncells);

  p4est_order_local_vertices (p4est, 0, &Ntotal, locidx_data);

  fprintf (vtufile, "<?xml version=\"1.0\"?>\n");
  fprintf (vtufile, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\"");
#if defined(P4EST_VTK_BINARY) && defined(P4EST_VTK_COMPRESSION)
  fprintf (vtufile, " compressor=\"vtkZLibDataCompressor\"");
#endif
#ifdef SC_WORDS_BIGENDIAN
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
#ifdef P4EST_VTK_ASCII
  fprintf (vtufile, "        <DataArray type=\"%s\" Name=\"Position\""
           " NumberOfComponents=\"3\" format=\"ascii\">\n",
           P4EST_VTK_FLOAT_NAME);
#else
  uint8_data = NULL;

  fprintf (vtufile, "        <DataArray type=\"%s\" Name=\"Position\""
           " NumberOfComponents=\"3\" format=\"binary\">\n",
           P4EST_VTK_FLOAT_NAME);
#endif
  /* loop over the trees */
  for (jl = first_local_tree, quad_count = 0; jl <= last_local_tree; ++jl) {
    tree = sc_array_index (trees, jl);
    /*
     * Note that we switch from right-hand-rule order for tree_to_vertex
     * to pixel order for v
     */
    v0 = tree_to_vertex[jl * 4 + 0];
    v1 = tree_to_vertex[jl * 4 + 1];
    v2 = tree_to_vertex[jl * 4 + 3];
    v3 = tree_to_vertex[jl * 4 + 2];

    v0x = vertices[v0 * 3 + 0];
    v0y = vertices[v0 * 3 + 1];
    v0z = vertices[v0 * 3 + 2];
    v1x = vertices[v1 * 3 + 0];
    v1y = vertices[v1 * 3 + 1];
    v1z = vertices[v1 * 3 + 2];
    v2x = vertices[v2 * 3 + 0];
    v2y = vertices[v2 * 3 + 1];
    v2z = vertices[v2 * 3 + 2];
    v3x = vertices[v3 * 3 + 0];
    v3y = vertices[v3 * 3 + 1];
    v3z = vertices[v3 * 3 + 2];

    quadrants = &tree->quadrants;
    num_quads = quadrants->elem_count;

    /* loop over the elements in the tree */
    for (il = 0; il < num_quads; ++il, ++quad_count) {
      quad = sc_array_index (quadrants, il);
      inth = (p4est_qcoord_t) (1 << (P4EST_MAXLEVEL - quad->level));
      h = intsize * inth;
      eta1 = intsize * quad->x;
      eta2 = intsize * quad->y;

      w0x = v0x * (1.0 - eta1) * (1.0 - eta2)
        + v1x * (eta1) * (1.0 - eta2)
        + v2x * (1.0 - eta1) * (eta2)
        + v3x * (eta1) * (eta2);

      w0y = v0y * (1.0 - eta1) * (1.0 - eta2)
        + v1y * (eta1) * (1.0 - eta2)
        + v2y * (1.0 - eta1) * (eta2)
        + v3y * (eta1) * (eta2);

      w0z = v0z * (1.0 - eta1) * (1.0 - eta2)
        + v1z * (eta1) * (1.0 - eta2)
        + v2z * (1.0 - eta1) * (eta2)
        + v3z * (eta1) * (eta2);

      w1x = v0x * (1.0 - eta1 - h) * (1.0 - eta2)
        + v1x * (eta1 + h) * (1.0 - eta2)
        + v2x * (1.0 - eta1 - h) * (eta2)
        + v3x * (eta1 + h) * (eta2);

      w1y = v0y * (1.0 - eta1 - h) * (1.0 - eta2)
        + v1y * (eta1 + h) * (1.0 - eta2)
        + v2y * (1.0 - eta1 - h) * (eta2)
        + v3y * (eta1 + h) * (eta2);

      w1z = v0z * (1.0 - eta1 - h) * (1.0 - eta2)
        + v1z * (eta1 + h) * (1.0 - eta2)
        + v2z * (1.0 - eta1 - h) * (eta2)
        + v3z * (eta1 + h) * (eta2);

      w2x = v0x * (1.0 - eta1) * (1.0 - eta2 - h)
        + v1x * (eta1) * (1.0 - eta2 - h)
        + v2x * (1.0 - eta1) * (eta2 + h)
        + v3x * (eta1) * (eta2 + h);

      w2y = v0y * (1.0 - eta1) * (1.0 - eta2 - h)
        + v1y * (eta1) * (1.0 - eta2 - h)
        + v2y * (1.0 - eta1) * (eta2 + h)
        + v3y * (eta1) * (eta2 + h);

      w2z = v0z * (1.0 - eta1) * (1.0 - eta2 - h)
        + v1z * (eta1) * (1.0 - eta2 - h)
        + v2z * (1.0 - eta1) * (eta2 + h)
        + v3z * (eta1) * (eta2 + h);

      w3x = v0x * (1.0 - eta1 - h) * (1.0 - eta2 - h)
        + v1x * (eta1 + h) * (1.0 - eta2 - h)
        + v2x * (1.0 - eta1 - h) * (eta2 + h)
        + v3x * (eta1 + h) * (eta2 + h);

      w3y = v0y * (1.0 - eta1 - h) * (1.0 - eta2 - h)
        + v1y * (eta1 + h) * (1.0 - eta2 - h)
        + v2y * (1.0 - eta1 - h) * (eta2 + h)
        + v3y * (eta1 + h) * (eta2 + h);

      w3z = v0z * (1.0 - eta1 - h) * (1.0 - eta2 - h)
        + v1z * (eta1 + h) * (1.0 - eta2 - h)
        + v2z * (1.0 - eta1 - h) * (eta2 + h)
        + v3z * (eta1 + h) * (eta2 + h);

      lv0 = locidx_data[4 * quad_count + 0];
      lv1 = locidx_data[4 * quad_count + 1];
      lv2 = locidx_data[4 * quad_count + 2];
      lv3 = locidx_data[4 * quad_count + 3];

      float_data[3 * lv0 + 0] = (P4EST_VTK_FLOAT_TYPE) w0x;
      float_data[3 * lv0 + 1] = (P4EST_VTK_FLOAT_TYPE) w0y;
      float_data[3 * lv0 + 2] = (P4EST_VTK_FLOAT_TYPE) w0z;

      float_data[3 * lv1 + 0] = (P4EST_VTK_FLOAT_TYPE) w1x;
      float_data[3 * lv1 + 1] = (P4EST_VTK_FLOAT_TYPE) w1y;
      float_data[3 * lv1 + 2] = (P4EST_VTK_FLOAT_TYPE) w1z;

      float_data[3 * lv2 + 0] = (P4EST_VTK_FLOAT_TYPE) w2x;
      float_data[3 * lv2 + 1] = (P4EST_VTK_FLOAT_TYPE) w2y;
      float_data[3 * lv2 + 2] = (P4EST_VTK_FLOAT_TYPE) w2z;

      float_data[3 * lv3 + 0] = (P4EST_VTK_FLOAT_TYPE) w3x;
      float_data[3 * lv3 + 1] = (P4EST_VTK_FLOAT_TYPE) w3y;
      float_data[3 * lv3 + 2] = (P4EST_VTK_FLOAT_TYPE) w3z;
    }
  }

#ifdef P4EST_VTK_ASCII
  for (il = 0; il < Ntotal; ++il) {
    w0x = float_data[3 * il + 0];
    w0y = float_data[3 * il + 1];
    w0z = float_data[3 * il + 2];

#ifdef P4EST_VTK_DOUBLES
    fprintf (vtufile, "     %24.16e %24.16e %24.16e\n", w0x, w0y, w0z);
#else
    fprintf (vtufile, "          %16.8e %16.8e %16.8e\n", w0x, w0y, w0z);
#endif
  }
#else
  fprintf (vtufile, "          ");
  retval = p4est_vtk_write_binary (vtufile, (char *) float_data,
                                   sizeof (*float_data) * 3 * Ntotal);
  fprintf (vtufile, "\n");
  if (retval) {
    P4EST_LERROR ("p4est_vtk: Error encoding points\n");
    fclose (vtufile);
    return -1;
  }
#endif
  P4EST_FREE (float_data);
  fprintf (vtufile, "        </DataArray>\n");
  fprintf (vtufile, "      </Points>\n");
  fprintf (vtufile, "      <Cells>\n");

  /* write connectivity data */
#ifdef P4EST_VTK_ASCII
  fprintf (vtufile, "        <DataArray type=\"%s\" Name=\"connectivity\""
           " format=\"ascii\">\n", P4EST_VTK_LOCIDX);
  for (il = 0; il < Ncells; ++il) {
    fprintf (vtufile, "          %lld %lld %lld %lld\n",
             (long long) locidx_data[4 * il + 0],
             (long long) locidx_data[4 * il + 1],
             (long long) locidx_data[4 * il + 2],
             (long long) locidx_data[4 * il + 3]);
  }
#else
  fprintf (vtufile, "        <DataArray type=\"%s\" Name=\"connectivity\""
           " format=\"binary\">\n", P4EST_VTK_LOCIDX);
  fprintf (vtufile, "          ");
  retval = p4est_vtk_write_binary (vtufile, (char *) locidx_data,
                                   sizeof (*locidx_data) * 4 * Ncells);
  fprintf (vtufile, "\n");
  if (retval) {
    P4EST_LERROR ("p4est_vtk: Error encoding connectivity\n");
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
    P4EST_LERROR ("p4est_vtk: Error encoding offsets\n");
    fclose (vtufile);
    return -1;
  }
#endif
  fprintf (vtufile, "        </DataArray>\n");

  /* write type data */
#ifdef P4EST_VTK_ASCII
  fprintf (vtufile, "        <DataArray type=\"UInt8\" Name=\"types\""
           " format=\"ascii\">\n");
  fprintf (vtufile, "         ");
  for (il = 0, sk = 1; il < Ncells; ++il, ++sk) {
    fprintf (vtufile, " 8");    /* 0:VTK_PIXEL */
    if (!(sk % 20) && il != (Ncells - 1))
      fprintf (vtufile, "\n         ");
  }
  fprintf (vtufile, "\n");
#else
  fprintf (vtufile, "        <DataArray type=\"UInt8\" Name=\"types\""
           " format=\"binary\">\n");
  uint8_data = P4EST_ALLOC (uint8_t, Ncells);
  for (il = 0; il < Ncells; ++il) {
    uint8_data[il] = 8;
  }
  fprintf (vtufile, "          ");
  retval = p4est_vtk_write_binary (vtufile, (char *) uint8_data,
                                   sizeof (*uint8_data) * Ncells);
  P4EST_FREE (uint8_data);
  fprintf (vtufile, "\n");
  if (retval) {
    P4EST_LERROR ("p4est_vtk: Error encoding types\n");
    fclose (vtufile);
    return -1;
  }
#endif
  fprintf (vtufile, "        </DataArray>\n");
  fprintf (vtufile, "      </Cells>\n");

  if (write_rank) {
    fprintf (vtufile, "      <CellData Scalars=\"mpirank\">\n");
#ifdef P4EST_VTK_ASCII
    fprintf (vtufile, "        <DataArray type=\"%s\" Name=\"mpirank\""
             " format=\"ascii\">\n", P4EST_VTK_LOCIDX);
    fprintf (vtufile, "         ");
    for (il = 0, sk = 1; il < Ncells; ++il, ++sk) {
      fprintf (vtufile, " %d", procRank);
      if (!(sk % 20) && il != (Ncells - 1))
        fprintf (vtufile, "\n         ");
    }
    fprintf (vtufile, "\n");
#else
    fprintf (vtufile, "        <DataArray type=\"%s\" Name=\"mpirank\""
             " format=\"binary\">\n", P4EST_VTK_LOCIDX);
    for (il = 0; il < Ncells; ++il) {
      locidx_data[il] = (p4est_locidx_t) procRank;
    }
    fprintf (vtufile, "          ");
    retval = p4est_vtk_write_binary (vtufile, (char *) locidx_data,
                                     sizeof (*locidx_data) * Ncells);
    fprintf (vtufile, "\n");
    if (retval) {
      P4EST_LERROR ("p4est_vtk: Error encoding types\n");
      fclose (vtufile);
      return -1;
    }
#endif
    fprintf (vtufile, "        </DataArray>\n");
    fprintf (vtufile, "      </CellData>\n");
  }

  fprintf (vtufile, "      <PointData>\n");
  P4EST_FREE (locidx_data);

  if (ferror (vtufile)) {
    P4EST_LERROR ("p4est_vtk: Error writing header\n");
    fclose (vtufile);
    return -1;
  }
  if (fclose (vtufile)) {
    P4EST_LERROR ("p4est_vtk: Error closing header\n");
    return -1;
  }
  vtufile = NULL;

  /* Only have the root write to the parallel vtk file */
  if (procRank == rootRank) {
    char                pvtufilename[BUFSIZ];
    FILE               *pvtufile;
    snprintf (pvtufilename, BUFSIZ, "%s.pvtu", baseName);

    pvtufile = fopen (pvtufilename, "wb");
    if (!pvtufile) {
      P4EST_LERRORF ("Could not open %s for output!\n", vtufilename);
      return -1;
    }

    fprintf (pvtufile, "<?xml version=\"1.0\"?>\n");
    fprintf (pvtufile, "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\"");
#if defined(P4EST_VTK_BINARY) && defined(P4EST_VTK_COMPRESSION)
    fprintf (pvtufile, " compressor=\"vtkZLibDataCompressor\"");
#endif
#ifdef SC_WORDS_BIGENDIAN
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
    if (write_rank) {
      fprintf (pvtufile, "    <PCellData Scalars=\"mpirank\">\n");
      fprintf (pvtufile,
               "      <PDataArray type=\"%s\" Name=\"mpirank\" format=\"%s\"/>\n",
               P4EST_VTK_LOCIDX, P4EST_VTK_FORMAT_STRING);
      fprintf (pvtufile, "    </PCellData>\n");
    }
    fprintf (pvtufile, "    <PPointData>\n");

    if (ferror (pvtufile)) {
      P4EST_LERROR ("p4est_vtk: Error writing parallel header\n");
      fclose (pvtufile);
      return -1;
    }
    if (fclose (pvtufile)) {
      P4EST_LERROR ("p4est_vtk: Error closing parallel header\n");
      return -1;
    }
  }

  return 0;
}

#endif /* !P4_TO_P8 */

int
p4est_vtk_write_footer (p4est_t * p4est, const char *baseName)
{
  char                vtufilename[BUFSIZ];
  int                 p;
  int                 procRank = p4est->mpirank;
  int                 numProcs = p4est->mpisize;
  FILE               *vtufile;

  P4EST_ASSERT (p4est->connectivity->vertices != NULL);

  /* Have each proc write to its own file */
  snprintf (vtufilename, BUFSIZ, "%s_%04d.vtu", baseName, procRank);
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
    char                pvtufilename[BUFSIZ];
    FILE               *pvtufile;
    snprintf (pvtufilename, BUFSIZ, "%s.pvtu", baseName);

    pvtufile = fopen (pvtufilename, "ab");
    if (!pvtufile) {
      P4EST_LERRORF ("Could not open %s for output!\n", vtufilename);
      return -1;
    }

    fprintf (pvtufile, "    </PPointData>\n");
    for (p = 0; p < numProcs; ++p) {
      fprintf (pvtufile, "    <Piece Source=\"%s_%04d.vtu\"/>\n", baseName,
               p);
    }
    fprintf (pvtufile, "  </PUnstructuredGrid>\n");
    fprintf (pvtufile, "</VTKFile>\n");

    if (ferror (pvtufile)) {
      P4EST_LERROR ("p4est_vtk: Error writing parallel footer\n");
      fclose (pvtufile);
      return -1;
    }
    if (fclose (pvtufile)) {
      P4EST_LERROR ("p4est_vtk: Error closing parallel footer\n");
      return -1;
    }
  }

  return 0;
}

/* EOF p4est_vtk.c */
