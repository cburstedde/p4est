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

#include <p4est_to_p8est.h>
#include "p4est_vtk.c"

double              p8est_vtk_default_scale = 0.97;
bool                p8est_vtk_default_write_rank = true;
bool                p8est_vtk_default_write_tree = true;

void
p8est_vtk_write_file (p4est_t * p4est, const char *baseName)
{
  int                 retval;

  SC_CHECK_ABORT (p4est->connectivity->vertices != NULL,
                  "Must provide connectivity with vertex information");

  retval = p8est_vtk_write_header (p4est, p8est_vtk_default_scale,
                                   p8est_vtk_default_write_rank,
                                   p8est_vtk_default_write_tree, baseName);
  SC_CHECK_ABORT (!retval, "VTK: write header");
  retval = p4est_vtk_write_footer (p4est, baseName);
  SC_CHECK_ABORT (!retval, "VTK: write footer");
}

int
p8est_vtk_write_header (p4est_t * p4est, double scale, bool write_rank,
                        bool write_tree, const char *baseName)
{
  p4est_connectivity_t *connectivity = p4est->connectivity;
  sc_array_t         *trees = p4est->trees;
  const int           mpirank = p4est->mpirank;
  const double        intsize = 1.0 / P8EST_ROOT_LEN;
  const double       *v = connectivity->vertices;
  const p4est_locidx_t Ncells = p4est->local_num_quadrants;
  const p4est_locidx_t Ntotal = P4EST_CHILDREN * Ncells;        /* type ok */
  const p4est_topidx_t first_local_tree = p4est->first_local_tree;
  const p4est_topidx_t last_local_tree = p4est->last_local_tree;
  const p4est_topidx_t *tree_to_vertex = connectivity->tree_to_vertex;
#ifdef P4EST_VTK_ASCII
  double              wx, wy, wz;
  p4est_locidx_t      sk;
#else
  int                 retval;
  uint8_t            *uint8_data;
  p4est_locidx_t     *locidx_data;
#endif
  int                 xi, yi, zi, j, k;
  double              h2, eta_x, eta_y, eta_z;
  size_t              num_quads, zz;
  p4est_topidx_t      jt;
  p4est_topidx_t      vt[P4EST_CHILDREN];
  p4est_locidx_t      quad_count;
  p4est_locidx_t      il;
  P4EST_VTK_FLOAT_TYPE *float_data;
  sc_array_t         *quadrants;
  p4est_tree_t       *tree;
  p4est_quadrant_t   *quad;
  char                vtufilename[BUFSIZ];
  FILE               *vtufile;

  P4EST_ASSERT (scale >= 0.0);
  P4EST_ASSERT (connectivity->vertices != NULL);

  /* Have each proc write to its own file */
  snprintf (vtufilename, BUFSIZ, "%s_%04d.vtu", baseName, mpirank);
  vtufile = fopen (vtufilename, "w");
  if (vtufile == NULL) {
    P4EST_LERRORF ("Could not open %s for output!\n", vtufilename);
    return -1;
  }

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

  float_data = P4EST_ALLOC (P4EST_VTK_FLOAT_TYPE, P4EST_DIM * Ntotal);

  /* write point position data */
  fprintf (vtufile, "        <DataArray type=\"%s\" Name=\"Position\""
           " NumberOfComponents=\"3\" format=\"%s\">\n",
           P4EST_VTK_FLOAT_NAME, P4EST_VTK_FORMAT_STRING);

  /* loop over the trees */
  for (jt = first_local_tree, quad_count = 0; jt <= last_local_tree; ++jt) {
    tree = p4est_array_index_topidx (trees, jt);
    quadrants = &tree->quadrants;
    num_quads = quadrants->elem_count;

    /* retrieve corners of the tree */
    for (k = 0; k < P4EST_CHILDREN; ++k)
      vt[k] = tree_to_vertex[jt * P4EST_CHILDREN + k];

    /* loop over the elements in the tree and calculated vertex coordinates */
    for (zz = 0; zz < num_quads; ++zz, ++quad_count) {
      quad = sc_array_index (quadrants, zz);
      h2 = .5 * intsize * P8EST_QUADRANT_LEN (quad->level);
      k = 0;
      for (zi = 0; zi < 2; ++zi)
        for (yi = 0; yi < 2; ++yi)
          for (xi = 0; xi < 2; ++xi) {
            P4EST_ASSERT (0 <= k && k < P4EST_CHILDREN);
            eta_x = intsize * quad->x + h2 * (1. + (xi * 2 - 1) * scale);
            eta_y = intsize * quad->y + h2 * (1. + (yi * 2 - 1) * scale);
            eta_z = intsize * quad->z + h2 * (1. + (zi * 2 - 1) * scale);

            for (j = 0; j < P4EST_DIM; ++j) {
              /* *INDENT-OFF* */
              float_data[P4EST_DIM * (P4EST_CHILDREN * quad_count + k) + j] =
                (P4EST_VTK_FLOAT_TYPE)
                (eta_x * (eta_y * (eta_z * v[3 * vt[0] + j] +
                                   (1. - eta_z) * v[3 * vt[4] + j]) +
                          (1. - eta_y) * (eta_z * v[3 * vt[2] + j] +
                                          (1. - eta_z) * v[3 * vt[6] + j])) +
                 (1. - eta_x) * (eta_y * (eta_z * v[3 * vt[1] + j] +
                                          (1. - eta_z) * v[3 * vt[5] + j]) +
                                 (1. - eta_y) * (eta_z * v[3 * vt[3] + j] +
                                                 (1. - eta_z) * v[3 * vt[7] + j])));
              /* *INDENT-ON* */
            }
            ++k;
          }
      P4EST_ASSERT (k == P4EST_CHILDREN);
    }
  }
  P4EST_ASSERT (P4EST_CHILDREN * quad_count == Ntotal);

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
  retval = p4est_vtk_write_binary (vtufile, (char *) float_data,
                                   sizeof (*float_data) * P4EST_DIM * Ntotal);
  fprintf (vtufile, "\n");
  if (retval) {
    P4EST_LERROR ("p8est_vtk: Error encoding points\n");
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
  for (il = 0; il < Ncells; ++il) {
    sk = P4EST_CHILDREN * il;   /* type ok */
    fprintf (vtufile, "          %lld %lld %lld %lld %lld %lld %lld %lld\n",
             (long long) sk + 0, (long long) sk + 1,
             (long long) sk + 2, (long long) sk + 3,
             (long long) sk + 4, (long long) sk + 5,
             (long long) sk + 6, (long long) sk + 7);
  }
#else
  locidx_data = P4EST_ALLOC (p4est_locidx_t, Ntotal);
  for (il = 0; il < Ntotal; ++il)
    locidx_data[il] = il;

  fprintf (vtufile, "          ");
  retval = p4est_vtk_write_binary (vtufile, (char *) locidx_data,
                                   sizeof (*locidx_data) * Ntotal);
  fprintf (vtufile, "\n");
  if (retval) {
    P4EST_LERROR ("p8est_vtk: Error encoding connectivity\n");
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
    P4EST_LERROR ("p8est_vtk: Error encoding offsets\n");
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
    fprintf (vtufile, " 11");
    if (!(sk % 20) && il != (Ncells - 1))
      fprintf (vtufile, "\n         ");
  }
  fprintf (vtufile, "\n");
#else
  uint8_data = P4EST_ALLOC (uint8_t, Ncells);
  for (il = 0; il < Ncells; ++il)
    uint8_data[il] = 11;

  fprintf (vtufile, "          ");
  retval = p4est_vtk_write_binary (vtufile, (char *) uint8_data,
                                   sizeof (*uint8_data) * Ncells);
  P4EST_FREE (uint8_data);
  fprintf (vtufile, "\n");
  if (retval) {
    P4EST_LERROR ("p8est_vtk: Error encoding types\n");
    fclose (vtufile);
    return -1;
  }
#endif
  fprintf (vtufile, "        </DataArray>\n");
  fprintf (vtufile, "      </Cells>\n");

  if (write_rank || write_tree) {
    fprintf (vtufile, "      <CellData Scalars=\"%s\">\n",
             !write_tree ? "mpirank" : !write_rank ? "treeid" :
             "mpirank,treeid");
  }
  if (write_rank) {
    fprintf (vtufile, "        <DataArray type=\"%s\" Name=\"mpirank\""
             " format=\"%s\">\n", P4EST_VTK_LOCIDX, P4EST_VTK_FORMAT_STRING);
#ifdef P4EST_VTK_ASCII
    fprintf (vtufile, "         ");
    for (il = 0, sk = 1; il < Ncells; ++il, ++sk) {
      fprintf (vtufile, " %d", mpirank);
      if (!(sk % 20) && il != (Ncells - 1))
        fprintf (vtufile, "\n         ");
    }
    fprintf (vtufile, "\n");
#else
    for (il = 0; il < Ncells; ++il)
      locidx_data[il] = (p4est_locidx_t) mpirank;

    fprintf (vtufile, "          ");
    retval = p4est_vtk_write_binary (vtufile, (char *) locidx_data,
                                     sizeof (*locidx_data) * Ncells);
    fprintf (vtufile, "\n");
    if (retval) {
      P4EST_LERROR ("p8est_vtk: Error encoding types\n");
      fclose (vtufile);
      return -1;
    }
#endif
    fprintf (vtufile, "        </DataArray>\n");
  }
  if (write_tree) {
    il = 0;
    fprintf (vtufile, "        <DataArray type=\"%s\" Name=\"treeid\""
             " format=\"%s\">\n", P4EST_VTK_LOCIDX, P4EST_VTK_FORMAT_STRING);
#ifdef P4EST_VTK_ASCII
    fprintf (vtufile, "         ");

    sk = 1;
    for (jt = first_local_tree; jt <= last_local_tree; ++jt) {
      tree = p4est_array_index_topidx (trees, jt);
      num_quads = tree->quadrants.elem_count;
      for (zz = 0; zz < num_quads; ++zz, ++sk, ++il) {
        fprintf (vtufile, " %lld", (long long) jt);
        if (!(sk % 20) && il != (Ncells - 1))
          fprintf (vtufile, "\n         ");
      }
    }
    fprintf (vtufile, "\n");
#else
    for (jt = first_local_tree; jt <= last_local_tree; ++jt) {
      tree = p4est_array_index_topidx (trees, jt);
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
      P4EST_LERROR ("p8est_vtk: Error encoding types\n");
      fclose (vtufile);
      return -1;
    }
#endif
    fprintf (vtufile, "        </DataArray>\n");
    P4EST_ASSERT (il == Ncells);
  }
  if (write_rank || write_tree) {
    fprintf (vtufile, "      </CellData>\n");
  }
#ifndef P4EST_VTK_ASCII
  P4EST_FREE (locidx_data);
#endif

  fprintf (vtufile, "      <PointData>\n");

  if (ferror (vtufile)) {
    P4EST_LERROR ("p8est_vtk: Error writing header\n");
    fclose (vtufile);
    return -1;
  }
  if (fclose (vtufile)) {
    P4EST_LERROR ("p8est_vtk: Error closing header\n");
    return -1;
  }
  vtufile = NULL;

  /* Only have the root write to the parallel vtk file */
  if (mpirank == 0) {
    char                pvtufilename[BUFSIZ];
    FILE               *pvtufile;
    snprintf (pvtufilename, BUFSIZ, "%s.pvtu", baseName);

    pvtufile = fopen (pvtufilename, "w");
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
    if (write_rank || write_tree) {
      fprintf (pvtufile, "    <PCellData Scalars=\"%s\">\n",
               !write_tree ? "mpirank" : !write_rank ? "treeid" :
               "mpirank,treeid");
    }
    if (write_rank) {
      fprintf (pvtufile,
               "      <PDataArray type=\"%s\" Name=\"mpirank\" format=\"%s\"/>\n",
               P4EST_VTK_LOCIDX, P4EST_VTK_FORMAT_STRING);
    }
    if (write_tree) {
      fprintf (pvtufile,
               "      <PDataArray type=\"%s\" Name=\"treeid\" format=\"%s\"/>\n",
               P4EST_VTK_LOCIDX, P4EST_VTK_FORMAT_STRING);
    }
    if (write_rank || write_tree) {
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

/* EOF p8est_vtk.c */
