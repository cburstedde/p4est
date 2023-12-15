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

static void
model_set_geom (p4est_gmt_model_t * model,
                const char *name, p4est_geometry_X_t X)
{
  model->sgeom.name = name;
  model->sgeom.user = model;
  model->sgeom.X = X;
  model->sgeom.destroy = NULL;
  model->model_geom = &model->sgeom;
}

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
model_synth_geom_X (p4est_geometry_t * geom, p4est_topidx_t which_tree,
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

static int
model_latlong_intersect (p4est_topidx_t which_tree, const double coord[4],
                         size_t m, void *vmodel)
{
  p4est_gmt_model_t  *model = (p4est_gmt_model_t *) vmodel;

  P4EST_ASSERT (model != NULL);
  P4EST_ASSERT (m < model->M);

  /* Rectangle coordinates are in [0, 1] for the numbered reference tree and
   * stored as { lower left x, lower left y, upper right x, upper right y }. */

  return 0;
}

static void
model_latlong_geom_X (p4est_geometry_t * geom, p4est_topidx_t which_tree,
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
  coastline_polygon_list_t * coast_poly_data = (coastline_polygon_list_t *) model_data;
  for( int i = 0; i < coast_poly_data->num_polygons; i++ ) {
    free (coast_poly_data->polygon_headers[i].pts);
  }
  free (coast_poly_data->polygon_headers);
  free (coast_poly_data);
}

p4est_gmt_model_t  *
p4est_gmt_model_latlong_new (p4est_gmt_model_latlong_params_t * params)
{
  p4est_gmt_model_t  *model = P4EST_ALLOC_ZERO (p4est_gmt_model_t, 1);

  /* the latlong models live on the unit square as reference domain */
  model->conn = p4est_connectivity_new_unitsquare ();

  /* load model properties */
  coastline_polygon_list_t * coast_poly = read_land_polygons_bin(params->load_filename, params->longitude , params->latitude); 
  model->model_data = coast_poly; 

  /* set virtual functions */
  model->intersect = model_latlong_intersect;
  model->destroy_data = p4est_gmt_model_latlong_destroy;   /* <- needs to free whatever is in model_data */

  /* setup input/output parameters */
  model->output_prefix = params->output_prefix; /*< Prefix defaults to NULL */
  model_set_geom (model, params->output_prefix, model_latlong_geom_X);

  /* the model is ready */
  model->M = coast_poly->num_line_segments; 
  return model;
}

/** are two bounding boxes overlapping ? */
int is_overlapping(double x1min,
                   double x1max,
                   double y1min,
                   double y1max,
                   double x2min,
                   double x2max,
                   double y2min,
                   double y2max)
{
  return ((x1min < x2max) && (x2min < x1max) && (y1min < y2max) && (y2min < y1max));
}

/* convert endianess from big to little */
int to_little_end(int i)
{
  unsigned char* data = (unsigned char*)&(i);
  int j = (data[3] << 0) | (data[2] << 8) | (data[1] << 16) | ((unsigned)data[0] << 24);
  return j;
}

coastline_polygon_list_t * read_land_polygons_bin(const char* filename,
                                                  double lon[2],
                                                  double lat[2])
{
  printf("Reading land poygons in BIN format from %s\n", filename);

  FILE* infile                  = fopen(filename, "r");
  int num_polygons              = 0;
  int num_line_segments         = 0;
  int global_line_segment_index = 0;
  
  gshhg_header_t* all_used       = (gshhg_header_t*)malloc(500000 * sizeof(gshhg_header_t));
  while (!feof(infile)) {
    gshhg_header_t poly_header;
    int h[11];
    int s = fread(h, 11, sizeof(int), infile);
    if (s > 0) {
      poly_header.id                        = to_little_end(h[0]);
      poly_header.n                         = to_little_end(h[1]);
      poly_header.flag                      = to_little_end(h[2]);
      poly_header.west                      = to_little_end(h[3]) / 1.0e6;
      poly_header.east                      = to_little_end(h[4]) / 1.0e6;
      poly_header.south                     = to_little_end(h[5]) / 1.0e6;
      poly_header.north                     = to_little_end(h[6]) / 1.0e6;
      poly_header.area                      = to_little_end(h[7]);
      poly_header.area_full                 = to_little_end(h[8]);
      poly_header.container                 = to_little_end(h[9]);
      poly_header.ancestor                  = to_little_end(h[10]);
      poly_header.global_line_segment_index = -1;

      // printf("Id %d with %d pts\n", poly_header.id, poly_header.n);

      int* pts           = (int*)malloc(2 * poly_header.n * sizeof(int));
      double* coord_list = (double*)malloc(2 * poly_header.n * sizeof(double));
      fread(pts, 2 * poly_header.n, sizeof(int), infile);

      for (int i = 0; i < poly_header.n; i++) {
        coord_list[2 * i] = to_little_end(pts[2 * i]) / 1.0e6;
        if (coord_list[2 * i] > 180.0) {
          coord_list[2 * i] -= 360.0;
        }
        coord_list[2 * i + 1] = to_little_end(pts[2 * i + 1]) / 1.0e6;
      }
      poly_header.pts = coord_list;
      int level       = poly_header.flag & 255;
      if ((level == 1) && (poly_header.container == -1)) {
        // ceck if bbox of polygon overlaps with region of intrest
        if (is_overlapping(poly_header.west, poly_header.east, poly_header.south, poly_header.north, lon[0], lon[1], lat[0], lat[1]))
        {
          poly_header.global_line_segment_index = global_line_segment_index;
          all_used[num_polygons]                = poly_header;
          num_polygons++;
          // polygons are closed, i.e. line segemnts are number of points - 1
          num_line_segments += poly_header.n - 1;
          // start index of global indexed segements
          global_line_segment_index += poly_header.n - 1;
        } else {
          free(coord_list);
        }
      } else {
        // printf("Level: %d und cont %d\n", level, poly_header.container);
        free(coord_list);
      }

      free(pts);
    }
  }
  printf("We have %d polygons meeting the requests\n", num_polygons);
  coastline_polygon_list_t* pl_ptr = (coastline_polygon_list_t*)malloc(1*sizeof(coastline_polygon_list_t));
  pl_ptr->polygon_headers   = all_used;
  pl_ptr->num_polygons      = num_polygons;
  pl_ptr->num_line_segments = num_line_segments;
  pl_ptr->west = lon[0];
  pl_ptr->east = lon[1];
  pl_ptr->south = lat[0];
  pl_ptr->north = lat[1];
  fclose(infile);
  return pl_ptr;
}

void
p4est_gmt_model_destroy (p4est_gmt_model_t * model)
{
  if (model->destroy_data != NULL) {
    /* only needed for non-trivial free code */
    model->destroy_data (model->model_data);
  }
  else {
    /* the default clears a standard allocation or respects NULL */
    P4EST_FREE (model->model_data);
  }
  P4EST_FREE (model);
}
