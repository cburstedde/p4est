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
  //TODO: in p4est_geometry_connectivity_X user is assumed to point to
  //the relevant connectivity. The following line will cause geometries
  // to violate that assumption.
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

p4est_gmt_model_t  *
p4est_gmt_model_latlong_new (p4est_gmt_model_latlong_params_t * params)
{
  p4est_gmt_model_t  *model = P4EST_ALLOC_ZERO (p4est_gmt_model_t, 1);

  /* the latlong models live on the unit square as reference domain */
  model->conn = p4est_connectivity_new_unitsquare ();

  /* load model properties */
  model->model_data = NULL;     /* <- Load something from params->load_filename,
                                   also deep copy the parameters into it. */

  /* set virtual functions */
  model->intersect = model_latlong_intersect;
  model->destroy_data = NULL;   /* <- needs to free whatever is in model_data */

  /* setup input/output parameters */
  model->output_prefix = params->output_prefix;
  model_set_geom (model, params->output_prefix, model_latlong_geom_X);

  /* the model is ready */
  model->M = 17;                /* <- update to actual value */
  return model;
}

/** Represents the intersection of a geodesic with one of the cubed
 * connectivity faces.
 * 
 * The endpoints p1 and p2 are given in tree-local coordinates and
 * in local coordinates the geodesic is just the line segment between
 * them.
*/
typedef struct sphere_geodesic_segment
{
  int                 which_tree;
  double              p1x, p1y, p2x, p2y;
} sphere_geodesic_segment_t;

typedef struct p4est_gmt_model_sphere
{
  int                         resolution;
  size_t                      num_geodesics;
  sphere_geodesic_segment_t  *geodesics;
}
p4est_gmt_model_sphere_t;

static void
model_sphere_destroy_data (void *vmodel_data)
{
  p4est_gmt_model_sphere_t *sdata = (p4est_gmt_model_sphere_t *) vmodel_data;
  P4EST_FREE (sdata->geodesics);
  P4EST_FREE (sdata);
}

/* Convert from angular coordinates to corresponding point on cube face*/
static void angular_to_cube(const double angular[2], double xyz[3]) {
    double inf_norm;
    double phi, theta;

    /* Convert to radians*/
    phi = angular[0]*M_PI/180.0;
    theta = angular[1]*M_PI/180.0;

    xyz[0] = sin(theta) * cos(phi);
    xyz[1] = sin(theta) * sin(phi);
    xyz[2] = cos(theta);

    inf_norm = fmax(fabs(xyz[0]), fmax(fabs(xyz[1]), fabs(xyz[2])))*2.0;
    xyz[0] /= inf_norm;
    xyz[1] /= inf_norm;
    xyz[2] /= inf_norm;
    return;
}

/**
 * Which face does a given point belong to
 * 
 * Note: ties are decided arbitrarily for the time being
 * 
 * \param[in] xyz  cartesian coordinates on surface of cube [-0.5,0.5]x[-0.5,0.5]
 */
static int point_to_tree(const double xyz[3]) {
    if (fabs(xyz[0]+0.5) < SC_EPS) return 2;
    if (fabs(xyz[0]-0.5) < SC_EPS) return 5;
    if (fabs(xyz[1]+0.5) < SC_EPS) return 4;
    if (fabs(xyz[1]-0.5) < SC_EPS) return 1;
    if (fabs(xyz[2]+0.5) < SC_EPS) return 0;
    if (fabs(xyz[2]-0.5) < SC_EPS) return 3;
    return -1; // This should not happen
}

/**
 * coordinate transformation from the surface of cube [-0.5,0.5]x[-0.5,0.5]
 * to tree-local coordinates. This is the inverse of p4est_geometry_cubed_X.
 *
 * \param[in]  xyz  cartesian coordinates in physical space
 * \param[out] rst  coordinates in AMR space : [0,1]^3
 * 
 */
static void p4est_geometry_cubed_Y(const double xyz[3], double rst[3]) {
    rst[2] = 0.0;

    int tree = point_to_tree(xyz);

    /* align center with origin */
    switch (tree)
    {
    case 0:
        rst[0] = xyz[1]+0.5;
        rst[1] = xyz[0]+0.5;
        break;
    case 1:
        rst[0] = xyz[2]+0.5;
        rst[1] = xyz[0]+0.5;
        break;
    case 2:
        rst[0] = xyz[2]+0.5;
        rst[1] = xyz[1]+0.5;
        break;
    case 3:
        rst[0] = xyz[0]+0.5;
        rst[1] = xyz[1]+0.5;
        break;
    case 4:
        rst[0] = xyz[0]+0.5;
        rst[1] = xyz[2]+0.5;
        break;
    case 5:
        rst[0] = xyz[1]+0.5;
        rst[1] = xyz[2]+0.5;
        break;
    default:
        break;
    }
}

/** Returns 1 if the line segments intersect, otherwise 0 */
static int lines_intersect(double p0_x, double p0_y, double p1_x, double p1_y, 
    double p2_x, double p2_y, double p3_x, double p3_y)
{
    double s1_x, s1_y, s2_x, s2_y;
    s1_x = p1_x - p0_x;     s1_y = p1_y - p0_y;
    s2_x = p3_x - p2_x;     s2_y = p3_y - p2_y;

    double s, t;
    s = (-s1_y * (p0_x - p2_x) + s1_x * (p0_y - p2_y)) / (-s2_x * s1_y + s1_x * s2_y);
    t = ( s2_x * (p0_y - p2_y) - s2_y * (p0_x - p2_x)) / (-s2_x * s1_y + s1_x * s2_y);

    if (s >= 0 && s <= 1 && t >= 0 && t <= 1)
    {
        // Collision detected
        return 1;
    }

    return 0; // No collision
}

static int
model_sphere_intersect (p4est_topidx_t which_tree, const double coord[4],
                         size_t m, void *vmodel)
{
  p4est_gmt_model_t  *model = (p4est_gmt_model_t *) vmodel;
  p4est_gmt_model_sphere_t *sdata;
  const sphere_geodesic_segment_t  *pco; /* mth geodesic segment */
  double                            hx, hy; /* width, height */

  P4EST_ASSERT (model != NULL);
  P4EST_ASSERT (m < model->M);
  sdata = (p4est_gmt_model_sphere_t *) model->model_data;
  P4EST_ASSERT (sdata != NULL && sdata->geodesics != NULL);
  pco = sdata->geodesics + m;
  P4EST_ASSERT (sdata->resolution >= 0);

  /* In this model we have 6 trees */
  P4EST_ASSERT (which_tree >= 0 && which_tree <= 5);

  /* Check the segment is on the relevant tree */
  if (pco->which_tree != which_tree) {
    return 0;
  }

  /* Rectangle coordinates are in [0, 1] for the numbered reference tree and
   * stored as { lower left x, lower left y, upper right x, upper right y }. */
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
  if (lines_intersect(pco->p1x, pco->p1y, pco->p2x, pco->p2y, coord[0], coord[1], coord[2], coord[1])) {
    return 1;
  }
  /* Check if L intersects the top edge of rectangle */
  if (lines_intersect(pco->p1x, pco->p1y, pco->p2x, pco->p2y, coord[0], coord[3], coord[2], coord[3])) {
    return 1;
  }
  /* Check if L intersects the left edge of rectangle */
  if (lines_intersect(pco->p1x, pco->p1y, pco->p2x, pco->p2y, coord[0], coord[1], coord[0], coord[3])) {
    return 1;
  }
  /* Check if L intersects the right edge of rectangle */
  if (lines_intersect(pco->p1x, pco->p1y, pco->p2x, pco->p2y, coord[2], coord[1], coord[2], coord[3])) {
    return 1;
  }
  
  /* Check if L is contained in the interior of rectangle.
   * Since we have already ruled out intersections it suffices
   * to check if one of the endpoints of L is in the interior.
  */
  if (pco->p1x >= coord[0] && pco->p1x <= coord[2] && pco->p1y >= coord[1] && pco->p1y <= coord[3]) {
    return 1;
  }

  /* We have exhausted the refinement criteria. */
  return 0;
}

p4est_gmt_model_t  *
p4est_gmt_model_sphere_new (int resolution)
{
  p4est_gmt_model_t  *model = P4EST_ALLOC_ZERO (p4est_gmt_model_t, 1);
  p4est_gmt_model_sphere_t *sdata = NULL;
  sphere_geodesic_segment_t *p;
  FILE               *fp;
  //double phi1, theta1, phi2, theta2;
  double angular1[2], xyz1[3], rst1[3];
  double angular2[2], xyz2[3], rst2[3];
  int which_tree_1, which_tree_2;


  /* the sphere model lives on the unit square as reference domain */
  model->conn = p4est_connectivity_new_cubed ();
  model->output_prefix = "sphere";
  model->model_data = sdata = P4EST_ALLOC (p4est_gmt_model_sphere_t, 1);

  /* Load geodesics */
  int n_geodesics = 406708;
  sdata->num_geodesics = model->M = n_geodesics;
  printf("Allocating\n");
  //p = sdata->geodesics = P4EST_ALLOC (double, n_geodesics*4);
  p = sdata->geodesics = P4EST_ALLOC (sphere_geodesic_segment_t, n_geodesics);
  printf("Allocated\n");

  fp = fopen("coastlines.csv","r");
  printf("File opened\n");

  int r = 0;
  while (!feof(fp) && r < n_geodesics)
  {
      /* Read geodesic in angular coordinates */
      fscanf(fp, "%lf,%lf,%lf,%lf", &angular1[0], &angular1[1], &angular2[0], &angular2[1]);
      // printf("Scanned %f, %f, %f, %f\n", angular1[0], angular1[1], angular2[0], angular2[1]);

      /* Convert to cartesian coordinates */
      angular_to_cube(angular1, xyz1);
      angular_to_cube(angular2, xyz2);
      // printf("Cube coords 1: %f, %f, %f\n", xyz1[0], xyz1[1], xyz1[2]);
      // printf("Cube coords 2: %f, %f, %f\n", xyz2[0], xyz2[1], xyz2[2]);


      /* Find which face the geodesic endpoints belong to */
      which_tree_1 = point_to_tree(xyz1);
      which_tree_2 = point_to_tree(xyz2);
      P4EST_ASSERT(which_tree_1 == which_tree_2);

      /* TODO eventually this should support inter-face geodesics.
       * Currently we require that geodesics do not cross faces.
       */

      /* Convert to tree-local coordinates*/
      p4est_geometry_cubed_Y(xyz1, rst1);
      p4est_geometry_cubed_Y(xyz2, rst2);
      // printf("Local coords 1: %f, %f\n", rst1[0], rst1[1]);
      // printf("Local coords 2: %f, %f\n", rst2[0], rst2[1]);
      // printf("\n");

      /* Store geodesic as a sphere_geodesic_segment_t */
      p[r].which_tree = which_tree_1;
      p[r].p1x = rst1[0];
      p[r].p1y = rst1[1];
      p[r].p2x = rst2[0];
      p[r].p2y = rst2[1];
      
      r += 1;
  }

  fclose(fp); /* Finished loading geodesics */
  
  /* Assign resolution, intersector and destructor */
  sdata->resolution = resolution;
  model->destroy_data = model_sphere_destroy_data;
  model->intersect = model_sphere_intersect;

  /* Assign geometry */
  //p4est_geometry_t *temp = p4est_geometry_new_sphere2d(model->conn, 1.0);
  //model->sgeom = *temp; /* We cannot directly create the geometry in sgeom */
  //model->model_geom = &model->sgeom;
  //p4est_geometry_destroy(temp); /* Delete temporary  */

  
  /* Note: the problem with the following is that it allocates memory externally,
   * rather than in sgeom. The destructor then doesn't know about this
  */
  model->model_geom = p4est_geometry_new_sphere2d(model->conn, 1.0);

  /* the model is ready */
  return model;
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
  p4est_geometry_destroy(model->model_geom);
  P4EST_FREE (model);
}