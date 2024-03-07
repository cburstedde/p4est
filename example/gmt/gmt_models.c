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

#include "gmt_models.h"
#include "gmt_global.h"
#include <sc_notify.h>
#include <p4est_search.h>

/************************ generic model code *********************/

static void
model_set_geom (p4est_gmt_model_t *model,
                const char *name, p4est_geometry_X_t X)
{
  model->sgeom.name = name;
  model->sgeom.user = model;
  model->sgeom.X = X;
  model->sgeom.destroy = NULL;
  model->model_geom = &model->sgeom;
}

/************** demonstration: synthetic model setup ***********************/

typedef struct p4est_gmt_model_synth
{
  int                 synthno;
  int                 resolution;
  size_t              num_points;
  double             *points;
}
p4est_gmt_model_synth_t;

static void
model_synth_destroy_data (void *vmodel_data)
{
  p4est_gmt_model_synth_t *sdata = (p4est_gmt_model_synth_t *) vmodel_data;
  P4EST_FREE (sdata->points);
  P4EST_FREE (sdata);
}

static int
model_synth_intersect (p4est_topidx_t which_tree, const double coord[4],
                       size_t m, void *vmodel)
{
  p4est_gmt_model_t  *model = (p4est_gmt_model_t *) vmodel;
  p4est_gmt_model_synth_t *sdata;
  const double       *pco;
  double              hx, hy;

  P4EST_ASSERT (model != NULL);
  P4EST_ASSERT (m < model->M);
  sdata = (p4est_gmt_model_synth_t *) model->model_data;
  P4EST_ASSERT (sdata != NULL && sdata->points != NULL);
  pco = sdata->points + 2 * m;
  P4EST_ASSERT (sdata->resolution >= 0);

  /* In this model we have only one tree, the unit square. */
  P4EST_ASSERT (which_tree == 0);

  /* Rectangle coordinates are in [0, 1] for the numbered reference tree and
   * stored as { lower left x, lower left y, upper right x, upper right y }. */

  /* We do not refine if target resolution is reached. */
  hx = coord[2] - coord[0];
  hy = coord[3] - coord[1];
  if (SC_MAX (hx, hy) <= pow (.5, sdata->resolution)) {
    return 0;
  }

  /* In this synthetic example the point IS the object.  There are no lines. */
  if ((coord[0] <= pco[0] && pco[0] <= coord[2]) &&
      (coord[1] <= pco[1] && pco[1] <= coord[3])) {
    return 1;
  }

  /* We have exhausted the refinement criteria. */
  return 0;
}

static void
model_synth_geom_X (p4est_geometry_t *geom, p4est_topidx_t which_tree,
                    const double abc[3], double xyz[3])
{
  /* In this model we have only one tree, the unit square. */
  P4EST_ASSERT (which_tree == 0);

  /* We work with the unit square as physical space. */
  memcpy (xyz, abc, 3 * sizeof (double));
}

p4est_gmt_model_t  *
p4est_gmt_model_synth_new (int synthno, int resolution)
{
  p4est_gmt_model_t  *model = P4EST_ALLOC_ZERO (p4est_gmt_model_t, 1);
  p4est_gmt_model_synth_t *sdata = NULL;
  double             *p;

  /* initalize model */
  switch (synthno) {
  case 0:
    model->output_prefix = "triangle";
    model->conn = p4est_connectivity_new_unitsquare ();
    model->model_data = sdata = P4EST_ALLOC (p4est_gmt_model_synth_t, 1);
    sdata->synthno = synthno;
    sdata->resolution = resolution;
    sdata->num_points = model->M = 3;
    p = sdata->points = P4EST_ALLOC (double, 6);
    p[0] = 0.2;
    p[1] = 0.1;
    p[2] = 0.7;
    p[3] = 0.4;
    p[4] = 0.5;
    p[5] = 0.8;
    model->destroy_data = model_synth_destroy_data;
    model->intersect = model_synth_intersect;
    model_set_geom (model, model->output_prefix, model_synth_geom_X);
    break;
    /* possibly add more cases that work with polygon segments */
  default:
    SC_ABORT_NOT_REACHED ();
  }

  /* return initialized model */
  P4EST_ASSERT (sdata != NULL);
  sdata->synthno = synthno;
  return model;
}

/************** demonstration: latitude/longitude rectangle ***********************/

/** reads the binary GSHHG data file (*.b) */
/** polygons for which their bounding box does not intersect with the bounding box lon[2] = {lon_min lon_max}, */
/** lat[2] = {lat_min lat_max} are discarded.  */
/** NOTE: only the bounding box is tested, not the polygon (there might by false positves)! */
static coastline_polygon_list_t *read_land_polygons_bin (const char *filename,
                                                         double lon[2],
                                                         double lat[2]);

static int
model_latlong_intersect (p4est_topidx_t which_tree, const double coord[4],
                         size_t m, void *vmodel)
{
#ifdef P4EST_ENABLE_DEBUG
  p4est_gmt_model_t  *model = (p4est_gmt_model_t *) vmodel;
#endif

  P4EST_ASSERT (model != NULL);
  P4EST_ASSERT (m < model->M);

  /* Rectangle coordinates are in [0, 1] for the numbered reference tree and
   * stored as { lower left x, lower left y, upper right x, upper right y }. */

  return 0;
}

static void
model_latlong_geom_X (p4est_geometry_t *geom, p4est_topidx_t which_tree,
                      const double abc[3], double xyz[3])
{
#if LATLONG_DATA_HAS_BEEN_PROGRAMMED
  p4est_gmt_model_t  *model = (p4est_gmt_model_t *) geom->user;

  /* put the parameters latitude, longitude into the model data */
  longitude =
    ((typecast into gmt lanlong model data *) model->model_data)->longitude;
  latitude = ...;

  xyz[0] = longitude[0] + (longitude->[1] - longitude[0]) * abc[0];
  xyz[1] = latitude[0] + (latitude->[1] - latitude[0]) * abc[1];
#else
  xyz[0] = abc[0];
  xyz[1] = abc[1];
#endif
  xyz[2] = 0.;
}

static void
p4est_gmt_model_latlong_destroy (void *model_data)
{
  int                 i;
  coastline_polygon_list_t *coast_poly_data =
    (coastline_polygon_list_t *) model_data;

  for (i = 0; i < coast_poly_data->num_polygons; i++) {
    P4EST_FREE (coast_poly_data->polygon_headers[i].pts);
  }
  P4EST_FREE (coast_poly_data->polygon_headers);
  P4EST_FREE (coast_poly_data);
}

p4est_gmt_model_t  *
p4est_gmt_model_latlong_new (p4est_gmt_model_latlong_params_t *params)
{
  p4est_gmt_model_t  *model = P4EST_ALLOC_ZERO (p4est_gmt_model_t, 1);

  /* the latlong models live on the unit square as reference domain */
  model->conn = p4est_connectivity_new_unitsquare ();

  /* load model properties */
  coastline_polygon_list_t *coast_poly =
    read_land_polygons_bin (params->load_filename, params->longitude,
                            params->latitude);
  model->model_data = coast_poly;

  /* set virtual functions */
  model->intersect = model_latlong_intersect;
  model->destroy_data = p4est_gmt_model_latlong_destroy;        /* <- needs to free whatever is in model_data */

  /* setup input/output parameters */
  model->output_prefix = params->output_prefix; /*< Prefix defaults to NULL */
  model_set_geom (model, params->output_prefix, model_latlong_geom_X);

  /* the model is ready */
  model->M = coast_poly->num_line_segments;
  return model;
}

/** are two bounding boxes overlapping ? */
static int
is_overlapping (double x1min, double x1max, double y1min, double y1max,
                double x2min, double x2max, double y2min, double y2max)
{
  return ((x1min < x2max) && (x2min < x1max) && (y1min < y2max)
          && (y2min < y1max));
}

/* convert endianess from big to little */
static int
to_little_end (int i)
{
  unsigned char      *data = (unsigned char *) &(i);
  int                 j =
    (data[3] << 0) | (data[2] << 8) | (data[1] << 16) | ((unsigned) data[0] <<
                                                         24);
  return j;
}

static coastline_polygon_list_t *
read_land_polygons_bin (const char *filename, double lon[2], double lat[2])
{
  printf ("Reading land poygons in BIN format from %s\n", filename);

  FILE               *infile = fopen (filename, "r");
  int                 num_polygons = 0;
  int                 num_line_segments = 0;
  int                 global_line_segment_index = 0;

  gshhg_header_t     *all_used = P4EST_ALLOC (gshhg_header_t, 500000);

  while (!feof (infile)) {
    gshhg_header_t      poly_header;
    int                 i;
    int                 h[11];
    int                 s = fread (h, 11, sizeof (int), infile);
    if (s > 0) {
      poly_header.id = to_little_end (h[0]);
      poly_header.n = to_little_end (h[1]);
      poly_header.flag = to_little_end (h[2]);
      poly_header.west = to_little_end (h[3]) / 1.0e6;
      poly_header.east = to_little_end (h[4]) / 1.0e6;
      poly_header.south = to_little_end (h[5]) / 1.0e6;
      poly_header.north = to_little_end (h[6]) / 1.0e6;
      poly_header.area = to_little_end (h[7]);
      poly_header.area_full = to_little_end (h[8]);
      poly_header.container = to_little_end (h[9]);
      poly_header.ancestor = to_little_end (h[10]);
      poly_header.global_line_segment_index = -1;

#if 0
      printf ("Id %d with %d pts\n", poly_header.id, poly_header.n);
#endif

      int                *pts = P4EST_ALLOC (int, 2 * poly_header.n);
      double             *coord_list =
        P4EST_ALLOC (double, 2 * poly_header.n);

      fread (pts, 2 * poly_header.n, sizeof (int), infile);

      for (i = 0; i < poly_header.n; i++) {
        coord_list[2 * i] = to_little_end (pts[2 * i]) / 1.0e6;
        if (coord_list[2 * i] > 180.0) {
          coord_list[2 * i] -= 360.0;
        }
        coord_list[2 * i + 1] = to_little_end (pts[2 * i + 1]) / 1.0e6;
      }
      poly_header.pts = coord_list;
      int                 level = poly_header.flag & 255;
      if ((level == 1) && (poly_header.container == -1)) {
        /* ceck if bbox of polygon overlaps with region of intrest */
        if (is_overlapping
            (poly_header.west, poly_header.east, poly_header.south,
             poly_header.north, lon[0], lon[1], lat[0], lat[1])) {
          poly_header.global_line_segment_index = global_line_segment_index;
          all_used[num_polygons] = poly_header;
          num_polygons++;
          /* polygons are closed, i.e. line segemnts are number of points - 1 */
          num_line_segments += poly_header.n - 1;
          /* start index of global indexed segements */
          global_line_segment_index += poly_header.n - 1;
        }
        else {
          P4EST_FREE (coord_list);
        }
      }
      else {
        /* printf("Level: %d und cont %d\n", level, poly_header.container); */
        P4EST_FREE (coord_list);
      }

      P4EST_FREE (pts);
    }
  }
  printf ("We have %d polygons meeting the requests\n", num_polygons);
  coastline_polygon_list_t *pl_ptr =
    P4EST_ALLOC (coastline_polygon_list_t, 1);

  pl_ptr->polygon_headers = all_used;
  pl_ptr->num_polygons = num_polygons;
  pl_ptr->num_line_segments = num_line_segments;
  pl_ptr->west = lon[0];
  pl_ptr->east = lon[1];
  pl_ptr->south = lat[0];
  pl_ptr->north = lat[1];
  fclose (infile);
  return pl_ptr;
}

/************** demonstration: geodesics on a spherical surface ***********************/

typedef struct p4est_gmt_model_sphere
{
  int                 resolution;
}
p4est_gmt_model_sphere_t;

static void
model_sphere_destroy_data (void *vmodel_data)
{
  p4est_gmt_model_sphere_t *sdata = (p4est_gmt_model_sphere_t *) vmodel_data;
  P4EST_FREE (sdata);
}

/** Returns 1 if the line segments (p0 to p1) and (p2 to p3) intersect, otherwise 0 */
/* would this function be general enough/of interest across models? */
static int
lines_intersect (double p0_x, double p0_y, double p1_x, double p1_y,
                 double p2_x, double p2_y, double p3_x, double p3_y)
{
  /* We solve the matrix equation (p1-p0, p2-p3) (s, t)^T = (p2-p0),
   * by inverting the matrix (p1-p0, p2-p3). */

  /* Precompute reused values for efficiency */
  double              s1_x, s1_y, s2_x, s2_y;
  s1_x = p1_x - p0_x;
  s1_y = p1_y - p0_y;
  s2_x = p3_x - p2_x;
  s2_y = p3_y - p2_y;

  /* Compute line intersection */
  double              s, t;
  s =
    (-s1_y * (p0_x - p2_x) + s1_x * (p0_y - p2_y)) / (-s2_x * s1_y +
                                                      s1_x * s2_y);
  t =
    (s2_x * (p0_y - p2_y) - s2_y * (p0_x - p2_x)) / (-s2_x * s1_y +
                                                     s1_x * s2_y);

  /* Check intersection lies on relevant segment */
  if (s >= 0 && s <= 1 && t >= 0 && t <= 1) {
    return 1;
  }

  return 0;
}

/** Returns 1 if the given geodesic intersects the given rectangle and 0 otherwise.
 *
 * \param[in] which_tree  tree id inside forest
 * \param[in] coord       rectangle for intersection checking. Rectangle coordinates
 *                        are in [0, 1] for the numbered reference tree and stored as
 *                        { lower left x, lower left y, upper right x, upper right y }.
 * \param[in] m           index of the geodesic we are checking
 * \param[in] vmodel      spherical model
 *
 */
static int
model_sphere_intersect (p4est_topidx_t which_tree, const double coord[4],
                        size_t m, void *vmodel)
{
  p4est_gmt_model_t  *model = (p4est_gmt_model_t *) vmodel;
  p4est_gmt_model_sphere_t *sdata;
  const p4est_gmt_sphere_geoseg_t *pco; /* mth geodesic segment */
  double              hx, hy;   /* width, height */

  P4EST_ASSERT (model != NULL);
  P4EST_ASSERT (m < model->M);
  sdata = (p4est_gmt_model_sphere_t *) model->model_data;
  P4EST_ASSERT (sdata != NULL && model->points != NULL);
  pco = ((p4est_gmt_sphere_geoseg_t *) model->points) + m;
  P4EST_ASSERT (sdata->resolution >= 0);

  /* In this model we have 6 trees */
  P4EST_ASSERT (which_tree >= 0 && which_tree <= 5);

  /* Check the segment is on the relevant tree */
  if (pco->which_tree != which_tree) {
    return 0;
  }

  /* We do not refine if target resolution is reached. */
  hx = coord[2] - coord[0];
  hy = coord[3] - coord[1];
  if (SC_MAX (hx, hy) <= pow (.5, sdata->resolution)) {
    return 0;
  }

  /* Check if the line segment L between p1 and p2 intersects the edges of
   * the rectangle.
   */

  /* Check if L intersects the bottom edge of rectangle */
  if (lines_intersect
      (pco->p1x, pco->p1y, pco->p2x, pco->p2y, coord[0], coord[1], coord[2],
       coord[1])) {
    return 1;
  }
  /* Check if L intersects the top edge of rectangle */
  if (lines_intersect
      (pco->p1x, pco->p1y, pco->p2x, pco->p2y, coord[0], coord[3], coord[2],
       coord[3])) {
    return 1;
  }
  /* Check if L intersects the left edge of rectangle */
  if (lines_intersect
      (pco->p1x, pco->p1y, pco->p2x, pco->p2y, coord[0], coord[1], coord[0],
       coord[3])) {
    return 1;
  }
  /* Check if L intersects the right edge of rectangle */
  if (lines_intersect
      (pco->p1x, pco->p1y, pco->p2x, pco->p2y, coord[2], coord[1], coord[2],
       coord[3])) {
    return 1;
  }

  /* Check if L is contained in the interior of rectangle.
   * Since we have already ruled out intersections it suffices
   * to check if one of the endpoints of L is in the interior.
   */
  if (pco->p1x >= coord[0] && pco->p1x <= coord[2] && pco->p1y >= coord[1]
      && pco->p1y <= coord[3]) {
    return 1;
  }

  /* We have exhausted the refinement criteria. */
  return 0;
}

p4est_gmt_model_t  *
p4est_gmt_model_sphere_new (int resolution, const char *input,
                            const char *output_prefix,
                            int dist, sc_MPI_Comm mpicomm)
{
  sc_MPI_File         file_handle;
  p4est_gmt_model_t  *model;
  p4est_gmt_model_sphere_t *sdata = NULL;
  p4est_gloidx_t      offset_mine, offset_next;
  p4est_gloidx_t      global_num_points = 0;
  p4est_locidx_t      local_num_points = 0;
  int                 local_int_bytes;
  int                 rank, num_procs;
  int                 mpiret;
  int                 mpival;
  int                 mpiall;
  int                 ocount;
  int                 mpireslen;
  char                mpierrstr[sc_MPI_MAX_ERROR_STRING];
  sc_MPI_Offset       mpi_offset;
  const char         *count_mismatch_message =
    "This should only occur when attempting to read "
    "beyond the bounds of the input file. "
    "If you correctly specified your input as the "
    "output of the preprocessing script then we "
    "expect that this error should never occur.\n";

  /* Get rank and number of processes */
  mpiret = sc_MPI_Comm_size (mpicomm, &num_procs);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (mpicomm, &rank);
  SC_CHECK_MPI (mpiret);

  /* clean initialization */
  mpival = sc_MPI_SUCCESS;
  mpiall = sc_MPI_SUCCESS;
  file_handle = sc_MPI_FILE_NULL;
  model = NULL;

  /* check for required parameters */
  if (input == NULL) {
    P4EST_GLOBAL_LERROR ("Sphere model expects non-NULL input filename.\n");
    P4EST_GLOBAL_LERROR ("Use the -F flag to set a filename.\n");
    return NULL;
  }

  /* collectively open file of precomputed geodesic segments */
  mpiret = sc_io_open (mpicomm, input, SC_IO_READ, sc_MPI_INFO_NULL,
                       &file_handle);

  /* check file open errors */
  if (mpiret != sc_MPI_SUCCESS) {
    mpiret = sc_MPI_Error_string (mpiret, mpierrstr, &mpireslen);
    SC_CHECK_MPI (mpiret);
    P4EST_GLOBAL_LERRORF ("Could not open input file: %s\n", input);
    P4EST_GLOBAL_LERRORF ("Error Code: %s\n", mpierrstr);
    P4EST_GLOBAL_LERROR ("Check you have run the preprocessing script.\n");
    P4EST_GLOBAL_LERROR ("Check you specified the input path correctly\n");
    return NULL;
  }

  if (rank == 0) {
    size_t              gnp;

    /* read the global number of points from file */
    mpiall = sc_io_read_at (file_handle, 0, &gnp,
                            sizeof (size_t), sc_MPI_BYTE, &ocount);
    global_num_points = (p4est_gloidx_t) gnp;

    /* check we read the expected number of bytes */
    if (mpiall == sc_MPI_SUCCESS && ocount != (int) sizeof (size_t)) {
      P4EST_GLOBAL_LERROR ("Count mismatch: reading number of points\n");
      P4EST_GLOBAL_LERROR (count_mismatch_message);
      mpiall = sc_MPI_ERR_OTHER;
    }
  }

  /* broadcast possible read errors */
  mpiret = sc_MPI_Bcast (&mpiall, 1, sc_MPI_INT, 0, mpicomm);
  SC_CHECK_MPI (mpiret);

  /* check read errors */
  if (mpiall != sc_MPI_SUCCESS) {
    mpiret = sc_MPI_Error_string (mpiall, mpierrstr, &mpireslen);
    SC_CHECK_MPI (mpiret);
    P4EST_GLOBAL_LERROR ("Error reading number of global points\n");
    P4EST_GLOBAL_LERRORF ("Error Code: %s\n", mpierrstr);

    /* cleanup on error */
    (void) sc_io_close (&file_handle);
    return NULL;
  }

  /* broadcast the global number of points */
  mpiret = sc_MPI_Bcast (&global_num_points, 1, P4EST_MPI_GLOIDX, 0, mpicomm);
  SC_CHECK_MPI (mpiret);

  /* set read offsets */
  mpi_offset = 0;

  /* set read offsets depending on whether we are running distributed */
  if (!dist) {
    mpi_offset = 0;
    local_num_points = (p4est_locidx_t) global_num_points;
  }
  else {
    /* offset to first point of current MPI process */
    offset_mine = p4est_partition_cut_gloidx (global_num_points,
                                              rank, num_procs);

    /* offset to first point of successor MPI process */
    offset_next = p4est_partition_cut_gloidx (global_num_points,
                                              rank + 1, num_procs);

    /* set file offset (in bytes) for this calling process */
    mpi_offset =
      (sc_MPI_Offset) (offset_mine * sizeof (p4est_gmt_sphere_geoseg_t));
    local_num_points = (p4est_locidx_t) (offset_next - offset_mine);
  }

  /* Check that the number of bytes being read does not overflow int.
   * We do this because by convention we record data size with a size_t,
   * whereas MPIIO uses int.
   * Note: we could define a custom MPI datatype rather than sending
   * bytes to increase this maximum
   */
  if (local_num_points * sizeof (p4est_gmt_sphere_geoseg_t) >
      (size_t) INT_MAX) {
    P4EST_GLOBAL_LERRORF ("Local number of points %lld is too big.\n",
                          (long long) local_num_points);

    /* cleanup on error */
    (void) sc_io_close (&file_handle);
    return NULL;
  }

  local_int_bytes =
    (int) (local_num_points * sizeof (p4est_gmt_sphere_geoseg_t));
  P4EST_ASSERT (local_int_bytes >= 0);

  /* allocate model */
  model = P4EST_ALLOC_ZERO (p4est_gmt_model_t, 1);
  model->model_data = sdata = P4EST_ALLOC (p4est_gmt_model_sphere_t, 1);
  model->points = P4EST_ALLOC (p4est_gmt_sphere_geoseg_t, local_num_points);

  /* Set point size */
  model->point_size = sizeof (p4est_gmt_sphere_geoseg_t);

  /* Set final geodesic count */
  model->M = local_num_points;

  /* Assign resolution, intersector and destructor */
  sdata->resolution = resolution;
  model->destroy_data = model_sphere_destroy_data;
  model->intersect = model_sphere_intersect;

  /* Assign connectivity */
  model->conn = p4est_connectivity_new_cubed ();

  /* Assign geometry */
  /* Note: the following allocates memory externally, rather than in sgeom. */
  model->model_geom = p4est_geometry_new_sphere2d (model->conn, 1.0);
  model->geom_allocated = 1;

  /* set default output prefix */
  if (output_prefix == NULL) {
    model->output_prefix = "sphere";
  }
  else {
    model->output_prefix = output_prefix;
  }

  /* each mpi process reads its data for its own offset */
  mpival = sc_io_read_at_all (file_handle, mpi_offset + sizeof (size_t),
                              model->points,
                              local_int_bytes, sc_MPI_BYTE, &ocount);

  /* check for read errors */
  if (mpival != sc_MPI_SUCCESS) {
    mpiret = sc_MPI_Error_string (mpiall, mpierrstr, &mpireslen);
    SC_CHECK_MPI (mpiret);
    P4EST_GLOBAL_LERROR ("Error reading geodesics from file\n");
    P4EST_GLOBAL_LERRORF ("Error Code: %s\n", mpierrstr);

    /* cleanup on error */
    (void) sc_io_close (&file_handle);
    p4est_gmt_model_destroy (model);
    return NULL;
  }

  /* check we read the expected number of bytes */
  if (ocount != local_int_bytes) {
    mpival = sc_MPI_ERR_OTHER;
  }

  /* communicate read count errors */
  /* note: LOR is correct as the standard mandates that MPI_SUCCESS == 0 */
  mpiret =
    sc_MPI_Allreduce (&mpival, &mpiall, 1, sc_MPI_INT, sc_MPI_LOR, mpicomm);
  SC_CHECK_MPI (mpiret);

  /* check read count errors */
  if (mpiall != sc_MPI_SUCCESS) {
    P4EST_GLOBAL_LERROR ("Count mismatch: reading geodesics\n");
    P4EST_GLOBAL_LERROR (count_mismatch_message);
    P4EST_GLOBAL_LERROR ("Error reading geodesics from file\n");

    /* cleanup on error */
    (void) sc_io_close (&file_handle);
    p4est_gmt_model_destroy (model);
    return NULL;
  }

  /* close the file collectively */
  mpival = sc_io_close (&file_handle);

  /* check file close error */
  if (mpival != sc_MPI_SUCCESS) {
    mpiret = sc_MPI_Error_string (mpival, mpierrstr, &mpireslen);
    SC_CHECK_MPI (mpiret);
    P4EST_GLOBAL_LERROR ("Error closing file\n");
    P4EST_GLOBAL_LERRORF ("Error Code: %s\n", mpierrstr);

    /* cleanup on error */
    p4est_gmt_model_destroy (model);
    return NULL;
  }

  /* the model is ready */
  return model;
}

/************************ generic model code *********************/

void
p4est_gmt_model_destroy (p4est_gmt_model_t *model)
{
  /* free geometry when dynamically allocated */
  if (model->geom_allocated) {
    p4est_geometry_destroy (model->model_geom);
  }

  /* free model-specific data */
  if (model->destroy_data != NULL) {
    /* only needed for non-trivial free code */
    model->destroy_data (model->model_data);
  }
  else {
    /* the default clears a standard allocation or respects NULL */
    P4EST_FREE (model->model_data);
  }

  /* free point data */
  P4EST_FREE (model->points);

  /* destroy connectivity outside of the specific model */
  p4est_connectivity_destroy (model->conn);
  P4EST_FREE (model);
}

/************************ generic communication code *********************/

static const double irootlen = 1. / (double) P4EST_ROOT_LEN;

/** quadrant callback for search_partition in compute_outgoing_points */
static int
partition_quadrant (p4est_t *p4est, p4est_topidx_t which_tree,
                    p4est_quadrant_t *quadrant, int pfirst,
                    int plast, void *point)
{

  P4EST_ASSERT (0 <= pfirst && pfirst <= plast);
  P4EST_ASSERT (point == NULL);

  /* always continue, as recursion is controlled by point callbacks */
  return 1;
}

/** point callback for search_partition in compute_outgoing_points */
static int
partition_point (p4est_t *p4est, p4est_topidx_t which_tree,
                 p4est_quadrant_t *quadrant, int pfirst, int plast,
                 void *point)
{
  /* global context */
  global_t           *g = (global_t *) p4est->user_pointer;
  p4est_gmt_model_t  *model = g->model;

  /* last process's send buffer this point was added to */
  int                 last_proc;
  /* quadrant coords */
  double              coord[4];
  p4est_qcoord_t      qh;
  /* point index */
  size_t              pi = *(size_t *) point;

  /* sanity checks */
  P4EST_ASSERT (g != NULL && g->p4est == p4est);
  P4EST_ASSERT (0 <= pfirst && pfirst <= plast);
  P4EST_ASSERT (point != NULL);
  P4EST_ASSERT (pi < model->num_resp);

  /* quadrant coordinates */
  qh = P4EST_QUADRANT_LEN (quadrant->level);
  coord[0] = irootlen * quadrant->x;
  coord[1] = irootlen * quadrant->y;
  coord[2] = irootlen * (quadrant->x + qh);
  coord[3] = irootlen * (quadrant->y + qh);

  /* if current quadrant has multiple owners */
  if (pfirst < plast) {
    /* point follows recursion when it intersects the quadrant */
    return model->intersect (quadrant->p.which_tree, coord, pi, model);
  }

  /* current quadrant has a single owner */
  P4EST_ASSERT (pfirst == plast);

  if (!model->intersect (quadrant->p.which_tree, coord, pi, model)) {
    /* point does not intersect this quadrant */
    return 0;
  }

  /* get last process whose domain we have already recorded as intersecting 
   * this point
   */
  last_proc = model->last_procs[pi];

  /* since we traverse in order we expect not to have seen this point in
   * in higher process domains yet
   */
  P4EST_ASSERT (last_proc <= pfirst);

  if (last_proc == pfirst) {
    /* we have found an already recorded process */
    return 0;
  }
  /* otherwise we have found a new process intersecting the point */

  /* record this new process */
  model->last_procs[pi] = pfirst;

  /* add point to corresponding send buffer */
  if (last_proc == -1) {
    /* process should own point and be responsible for its propagation */

    /* initialise (resp) message to pfirst if it not already initialised */
    if (model->resp.to_send[pfirst] == NULL) {
      model->resp.to_send[pfirst] = sc_array_new (model->point_size);
    }

    /* add point to message */
    memcpy (sc_array_push (model->resp.to_send[pfirst]),
            ((char *) model->points) + pi * model->point_size,
            model->point_size);
  }
  else {
    /* process should own point but not be responsible for its propagation */

    /* initialise (own) message to pfirst if it not already initialised */
    if (model->own.to_send[pfirst] == NULL) {
      model->own.to_send[pfirst] = sc_array_new (model->point_size);
    }

    /* add point to message */
    memcpy (sc_array_push (model->own.to_send[pfirst]),
            ((char *) model->points) + pi * model->point_size,
            model->point_size);
  }

  /* end recursion */
  return 0;
}

/** Prepare outgoing buffers of points to propagate.
 * 
 * \param[out] resp Comm data for points whose receiver will be responsible
 *                  for propagating them in the next iteration
 * \param[out] own  Comm data for points whose receiver will *not* be
 *                  responsible for propagating in the next iteration
 * \param[in] p4est The forest
 * \param[in] model Sphere model
 * \param[in] num_procs Number of MPI processes
 */
static void
compute_outgoing_points (p4est_gmt_comm_t *resp,
                         p4est_gmt_comm_t *own,
                         p4est_t *p4est,
                         p4est_gmt_model_t *model, int num_procs)
{
  sc_array_t         *points;

  /* initialise to -1 to signify no points have been added to send buffers */
  /* Here we are relying on the fact that the char -1 is 11111111 in bits,
   * and so the resulting int array will be filled with -1.
   */
  model->last_procs = P4EST_ALLOC (int, model->num_resp);
  memset (model->last_procs, -1, model->num_resp * sizeof (int));

  /* initialise index of outgoing message buffers */
  own->to_send = P4EST_ALLOC (sc_array_t *, num_procs);
  resp->to_send = P4EST_ALLOC (sc_array_t *, num_procs);

  /* initialise outgoing message buffers to NULL */
  for (int q = 0; q < num_procs; q++) {
    own->to_send[q] = NULL;
    resp->to_send[q] = NULL;
  }

  /* set up search objects for partition search */
  points = sc_array_new_count (sizeof (size_t), model->num_resp);
  for (size_t zz = 0; zz < model->num_resp; ++zz) {
    *(size_t *) sc_array_index (points, zz) = zz;
  }

  /* add points to the relevant send buffers (by partition search) */
  p4est_search_partition (p4est, 0, partition_quadrant,
                          partition_point, points);

  /* clean up */
  sc_array_destroy_null (&points);
  P4EST_FREE (model->last_procs);
}

/** Update communication data with which processes p is sending points to and
 *  how many points are being sent to each of these. 
 * 
 *  The output is stored in the fields comm->receivers and comm->recvs_counts.
 *  We assume comm->to_send is already populated.
 * 
 * \param[in,out] comm communication data.
 * \param[in] num_procs number of mpi processes
 * \param[in] rank rank of the local process
 * \param[in] point_size byte size of points in this model 
 */
static int
compute_receivers (p4est_gmt_comm_t *comm, int num_procs, int rank,
                   size_t point_size)
{
  int                 err = 0;

  /* initialize receivers and receiver counts */
  comm->receivers = sc_array_new (sizeof (int));
  comm->recvs_counts = sc_array_new (sizeof (size_t));

  /* compute receivers and counts */
  for (int q = 0; q < num_procs; q++) {
    if (comm->to_send[q] != NULL) {
      /* check that the number of points communicated will not overflow int */
      if (comm->to_send[q]->elem_count * point_size > (size_t) INT_MAX) {
        P4EST_LERRORF ("Message of %lld points from rank %d to rank %d is "
                       "too large.\n",
                       (long long) comm->to_send[q]->elem_count, rank, q);
        err = 1;
        break;
      }

      /* add q to receivers */
      *(int *) sc_array_push (comm->receivers) = q;

      /* record how many points p is sending to q */
      *(size_t *) sc_array_push (comm->recvs_counts) =
        comm->to_send[q]->elem_count;
    }
  }
  P4EST_ASSERT (comm->receivers->elem_count ==
                comm->recvs_counts->elem_count);

  return err;
}

/** Update communication data with total number of incoming points, and 
 *  offsets to receive incoming points at. 
 * 
 *  The outputs are stored in the fields comm->num_incoming and comm->offsets.
 *  We assume that comm->senders and comm->senders_counts are already
 *  populated.
 * 
 * \param[in,out] comm communication data.
 */
static void
compute_offsets (p4est_gmt_comm_t *comm, size_t point_size)
{
  /* initialize offset array */
  comm->num_incoming = 0;
  comm->offsets = P4EST_ALLOC (size_t, comm->senders->elem_count);

  /* compute offsets */
  for (int i = 0; i < (int) comm->senders->elem_count; i++) {
    comm->offsets[i] = comm->num_incoming * point_size;
    comm->num_incoming +=
      *(size_t *) sc_array_index_int (comm->senders_counts, i);
  }
}

/** Post non-blocking sends for points in the given communication data.
 * 
 * To each rank q in comm->receivers we send the points stored at 
 * comm->to_send[q]
 * 
 * \param[in]   comm        communication data
 * \param[in]   mpicomm     MPI communicator
 * \param[out]  req         request storage of same length as comm->receivers
 * \param[in]   point_size  size of a single point
 */
static void
post_sends (p4est_gmt_comm_t *comm,
            sc_MPI_Comm mpicomm, sc_MPI_Request *req, size_t point_size)
{
  int                 mpiret;
  int                 q;
  int                 rank;

  /* get rank */
  mpiret = sc_MPI_Comm_rank (mpicomm, &rank);
  SC_CHECK_MPI (mpiret);

  /* for each receiver q */
  for (int i = 0; i < (int) comm->receivers->elem_count; i++) {
    q = *(int *) sc_array_index_int (comm->receivers, i);

    if (q != rank) {
      /* post non-blocking send of points to q */
      mpiret = sc_MPI_Isend (sc_array_index (comm->to_send[q], 0),
                             comm->to_send[q]->elem_count * point_size,
                             sc_MPI_BYTE, q, 0, mpicomm, req + i);
      SC_CHECK_MPI (mpiret);
    }
    else {
      /* we do not send a message to ourself with MPI; copy is faster */
      req[i] = sc_MPI_REQUEST_NULL;
    }
  }
}

/** Post non-blocking receives for senders in the given communication data.
 *  If there is a message for ourself then we copy it manually here rather
 *  than receiving it through MPI.
 * 
 *  We expect to receive points from each sender in comm->senders. The number
 *  of points each sender is sending is stored in comm->senders_counts (with
 *  corresponding indexing). We receive each message at the offset stored in
 *  comm->offsets (again with corresponding indexing).
 * 
 *  \param[in] comm communication data
 *  \param[in,out] recv_buffer points to array where received points are stored
 *  \param[out] req request storage of same length as comm->senders
 *  \param[in] mpicomm MPI communicator
 *  \param[in] point_size size of a single point
*/
static void
post_receives (p4est_gmt_comm_t *comm,
               void *recv_buffer,
               sc_MPI_Request *req, sc_MPI_Comm mpicomm, size_t point_size)
{
  int                 mpiret;
  int                 q;
  int                 rank;
  void               *self_dest = NULL;

  /* get rank */
  mpiret = sc_MPI_Comm_rank (mpicomm, &rank);
  SC_CHECK_MPI (mpiret);

  /* for each sender q */
  for (int i = 0; i < (int) comm->senders->elem_count; i++) {
    q = *(int *) sc_array_index_int (comm->senders, i);

    if (q != rank) {
      /* post non-blocking receive for points from q */
      mpiret = sc_MPI_Irecv (((char *) recv_buffer) + comm->offsets[i],
                             (*(size_t *)
                              sc_array_index_int (comm->senders_counts,
                                                  i)) * point_size,
                             sc_MPI_BYTE, q, 0, mpicomm, req + i);
      SC_CHECK_MPI (mpiret);
    }
    else {
      /* we do not receive a message from ourself with MPI; copy is faster */
      req[i] = sc_MPI_REQUEST_NULL;
      /* set receive buffer */
      self_dest = ((char *) recv_buffer) + comm->offsets[i];
    }
  }
  /* receive requests have been posted */

  /* if the message to ourself was non-empty */
  if (self_dest != NULL) {
    /* copy message to self manually rather than post receive request */
    memcpy (self_dest, sc_array_index (comm->to_send[rank], 0),
            comm->to_send[rank]->elem_count * point_size);
  }
}

/** Clean up data used to send after sends have completed */
static void
destroy_send_data (p4est_gmt_comm_t *comm, int num_procs)
{
  sc_array_destroy_null (&comm->receivers);
  sc_array_destroy_null (&comm->recvs_counts);
  for (int q = 0; q < num_procs; q++) {
    if (comm->to_send[q] != NULL) {
      sc_array_destroy_null (&comm->to_send[q]);
    }
  }
  P4EST_FREE (comm->to_send);
}

/** Clean up data used to receive after receives have completed */
static void
destroy_recv_data (p4est_gmt_comm_t *comm)
{
  sc_array_destroy (comm->senders);
  sc_array_destroy (comm->senders_counts);
  P4EST_FREE (comm->offsets);
}

int
p4est_gmt_communicate_points (p4est_t *p4est, p4est_gmt_model_t *model)
{
  int                 mpiret;
  int                 num_procs, rank;
  sc_MPI_Comm         mpicomm = p4est->mpicomm;
  int                 errsend = 0;
  int                 err = 0;

  /* Communication data */
  p4est_gmt_comm_t   *own = &model->own;        /**< not responsible */
  p4est_gmt_comm_t   *resp = &model->resp;      /**< responsible */
  sc_MPI_Request     *send_req = NULL; /**< requests for sending to receivers */
  sc_MPI_Request     *recv_req = NULL; /**< requests for receiving from senders */

  /* get rank and total process count */
  mpiret = sc_MPI_Comm_rank (mpicomm, &rank);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_size (mpicomm, &num_procs);
  SC_CHECK_MPI (mpiret);

  /* prepare outgoing buffers of points to send */
  compute_outgoing_points (resp, own, p4est, model, num_procs);

  /* free the points we received last iteration */
  P4EST_FREE (model->points);
  model->points = NULL;

  /* record which processes p is sending points to and how many points each
     process receives */
  /* note: an error is recorded here if p is attempting to send more than 
     INT_MAX bytes in an own or resp message to another process. We defer
     synchronising these errors until just before calling sc_notify_ext
     to avoid creating an unnecessary synchronisation point */
  errsend = compute_receivers (own, num_procs, rank, model->point_size);
  errsend = errsend
    || compute_receivers (resp, num_procs, rank, model->point_size);

  /* if messages from this process are valid then continue optimistically */
  if (!errsend) {
    /* initialise outgoing request arrays */
    send_req =
      P4EST_ALLOC (sc_MPI_Request,
                   own->receivers->elem_count + resp->receivers->elem_count);

    /* post non-blocking sends */
    post_sends (resp, mpicomm, send_req, model->point_size);
    post_sends (own, mpicomm, send_req + resp->receivers->elem_count,
                model->point_size);
  }

  /* synchronise possible message errors */
  mpiret =
    sc_MPI_Allreduce (&errsend, &err, 1, sc_MPI_INT, sc_MPI_LOR, mpicomm);
  SC_CHECK_MPI (mpiret);

  /* clean up and exit on message error */
  if (err) {
    /* clean up send data */
    destroy_send_data (own, num_procs);
    destroy_send_data (resp, num_procs);
    P4EST_FREE (send_req);

    /* return failure */
    return 1;
  }

  /* initialize buffers for receiving communication data with sc_notify_ext */
  own->senders = sc_array_new (sizeof (int));
  resp->senders = sc_array_new (sizeof (int));
  own->senders_counts = sc_array_new (sizeof (size_t));
  resp->senders_counts = sc_array_new (sizeof (size_t));

  /* notify processes receiving points from p and determine processes sending
     to p. Also exchange counts of points being sent. */
  sc_notify_ext (own->receivers, own->senders, own->recvs_counts,
                 own->senders_counts, mpicomm);
  sc_notify_ext (resp->receivers, resp->senders, resp->recvs_counts,
                 resp->senders_counts, mpicomm);

  /* sanity checks */
  P4EST_ASSERT (own->senders->elem_count == own->senders_counts->elem_count);
  P4EST_ASSERT (resp->senders->elem_count ==
                resp->senders_counts->elem_count);

  /* compute offsets for storing incoming points */
  compute_offsets (own, model->point_size);
  compute_offsets (resp, model->point_size);

  /* update counts of known points */
  model->M = resp->num_incoming + own->num_incoming;
  model->num_resp = resp->num_incoming;
  model->num_own = own->num_incoming;

  /* allocate memory for incoming points */
  model->points =
    sc_malloc (p4est_package_id, (model->M) * model->point_size);

  /* initialise incoming request arrays */
  recv_req =
    P4EST_ALLOC (sc_MPI_Request,
                 own->senders->elem_count + resp->senders->elem_count);

  /* post non-blocking receives */
  post_receives (resp, model->points, recv_req, mpicomm, model->point_size);
  post_receives (own,
                 ((char *) model->points) +
                 resp->num_incoming * model->point_size,
                 recv_req + resp->senders->elem_count, mpicomm,
                 model->point_size);

  /* wait for messages to send */
  mpiret =
    sc_MPI_Waitall (own->receivers->elem_count + resp->receivers->elem_count,
                    send_req, sc_MPI_STATUSES_IGNORE);
  SC_CHECK_MPI (mpiret);

  /* clean up data used for sending */
  destroy_send_data (own, num_procs);
  destroy_send_data (resp, num_procs);
  P4EST_FREE (send_req);

  /* Wait to receive messages */
  mpiret =
    sc_MPI_Waitall (own->senders->elem_count + resp->senders->elem_count,
                    recv_req, sc_MPI_STATUSES_IGNORE);
  SC_CHECK_MPI (mpiret);

  /* clean up data used for receiving */
  destroy_recv_data (own);
  destroy_recv_data (resp);
  P4EST_FREE (recv_req);

  /* return success */
  return 0;
}
