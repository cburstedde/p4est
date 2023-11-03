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
 * \param[in]  which_tree face of cube
 * \param[out] rst  coordinates in AMR space : [0,1]^3
 * 
 */
static void p4est_geometry_cubed_Y(const double xyz[3], double rst[3], int which_tree) {
    rst[2] = 0.0;

    //int tree = point_to_tree(xyz); //TODO this is resulting in the wrong tree for edge values
    //printf("Tree in local transform %d\n", tree);

    /* align center with origin */
    switch (which_tree)
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

/** Returns 1 if the line segments (p0 to p1) and (p2 to p3) intersect, otherwise 0 */
static int lines_intersect(double p0_x, double p0_y, double p1_x, double p1_y, 
    double p2_x, double p2_y, double p3_x, double p3_y)
{
    /* We solve the matrix equation (p1-p0, p2-p3) (s, t)^T = (p2-p0),
     * by inverting the matrix (p1-p0, p2-p3). */

    /* Precompute reused values for efficiency */
    double s1_x, s1_y, s2_x, s2_y;
    s1_x = p1_x - p0_x;     s1_y = p1_y - p0_y;
    s2_x = p3_x - p2_x;     s2_y = p3_y - p2_y;

    /* Compute line intersection */
    double s, t;
    s = (-s1_y * (p0_x - p2_x) + s1_x * (p0_y - p2_y)) / (-s2_x * s1_y + s1_x * s2_y);
    t = ( s2_x * (p0_y - p2_y) - s2_y * (p0_x - p2_x)) / (-s2_x * s1_y + s1_x * s2_y);

    /* Check intersection lies on relevant segment */
    if (s >= 0 && s <= 1 && t >= 0 && t <= 1)
    {
        return 1;
    }

    return 0;
}

/** Returns true if the cone spanned by v1 and v2 intersects the line segment
 *  between p1 and p2. If an intersection is detected then p_intersect is
 *  set to the computed intersection point.
 * 
*/
static int cone_line_intersection(const double v1[3], const double v2[3], const double p1[3], 
                              const double p2[3], double p_intersect[3]) 
{
  /* We solve the matrix equation (v1, v2, p1-p2) x = p1 by inverting (v1, v2, p1-p2) */
  double A[3][3]; /* Matrix we are inverting */
  double cofactor[3][3]; /* Cofactor matrix */
  double det_A, det_A_inv;
  double x[3];
  
  A[0][0] = v1[0];
  A[1][0] = v1[1];
  A[2][0] = v1[2];
  A[0][1] = v2[0];
  A[1][1] = v2[1];
  A[2][1] = v2[2];
  A[0][2] = p1[0]-p2[0];
  A[1][2] = p1[1]-p2[1];
  A[2][2] = p1[2]-p2[2];

  /* Compute minors */
  for (int r = 0; r < 3; r++) {
    for (int c = 0; c < 3; c++) {
      cofactor[r][c] = A[(r+1)%3][(c+1)%3] * A[(r+2)%3][(c+2)%3]
                        - A[(r+1)%3][(c+2)%3] * A[(r+2)%3][(c+1)%3];
    }
  }

  /* Compute the determinant by Laplace expansion along the first column */
  det_A = 0;
  for (int r=0; r<3; r++) {
    det_A += A[r][0] * cofactor[r][0];
  }
  det_A_inv = 1/det_A; /* TODO: What if det_A = 0?*/
  

  /* Multiply p1 with inverse of A */
  for (int i=0; i<3; i++) {
    /* Compute x[i] */
    x[i] = 0;
    for (int j=0; j<3; j++) {
      x[i] += cofactor[j][i] * p1[j];
    }
    x[i] *= det_A_inv;
  }
  
  if (x[0] < 0 || x[1] < 0 || x[2] < 0 || x[2] > 1) {
    return 0; /* Invalid intersection */
  }

  p_intersect[0] = x[0] * v1[0] + x[1] * v2[0];
  p_intersect[1] = x[0] * v1[1] + x[1] * v2[1];
  p_intersect[2] = x[0] * v1[2] + x[1] * v2[2];
  printf("Intersection at coords %f, %f, %f\n", p_intersect[0], p_intersect[1], p_intersect[2]);
  return 1; /* Valid intersection */
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

/** If the geodesic between xyz1 and xyz2 intersects the given edge then add this
 *  intersection point to endpoints and increment the corresponding entry in
 *  endpoints_count. To deal with edge cases coming from corners we should only update
 *  if the new intersection point is distinct to previously seen intersection points.
 */
static void update_endpoints(const double xyz1[3], const double xyz2[3], int edge,
                              double endpoints[6][2][3], int endpoints_count[6])
{
  double p_intersect[3];
  int detected;

  /* Which cube faces are adjacent to the given edge */
  const int edge_to_face[12][2] = {
    {0,1},
    {0,2},
    {0,4},
    {0,5},
    {1,2},
    {1,3},
    {1,5},
    {2,3},
    {2,4},
    {3,4},
    {3,5},
    {4,5}
  }; 

  /* Cube edge endpoint coordinates */
  const double edge_endpoints[12][2][3] = {
    {{-0.5, 0.5, -0.5}, {0.5, 0.5, -0.5}}, /* 0,1 edge */
    {{-0.5, -0.5, -0.5}, {-0.5, 0.5, -0.5}}, /* 0,2 edge */
    {{-0.5, -0.5, -0.5}, {0.5, -0.5, -0.5}}, /* 0,4 edge*/
    {{0.5, -0.5, -0.5}, {0.5, 0.5, -0.5}}, /* 0,5 edge */
    {{-0.5, 0.5, -0.5}, {-0.5, 0.5, 0.5}}, /* 1,2 edge */
    {{-0.5, 0.5, 0.5}, {0.5, 0.5, 0.5}}, /* 1,3 edge */
    {{0.5, 0.5, -0.5}, {0.5, 0.5, 0.5}}, /* 1,5 edge */
    {{-0.5, -0.5, 0.5}, {-0.5, 0.5, 0.5}}, /* 2,3 edge */
    {{-0.5, -0.5, -0.5}, {-0.5, -0.5, 0.5}}, /* 2,4 edge */
    {{-0.5, -0.5, 0.5}, {0.5, -0.5, 0.5}}, /* 3,4 edge */
    {{0.5, -0.5, 0.5}, {0.5, 0.5, 0.5}}, /* 3,5 edge */
    {{0.5, -0.5, -0.5}, 0.5, -0.5, 0.5} /* 4,5 edge*/
  }; 

  detected = cone_line_intersection(xyz1, xyz2, edge_endpoints[edge][0], 
                            edge_endpoints[edge][1], p_intersect);
  
  if (detected == 1) {
    printf("Detection was on edge %d\n", edge);
    for (int i = 0; i < 2; i++) {
      if (endpoints_count[edge_to_face[edge][i]] == 0) { 
        /* Record first endpoint */
        endpoints[edge_to_face[edge][i]][0][0] = p_intersect[0];
        endpoints[edge_to_face[edge][i]][0][1] = p_intersect[1];
        endpoints[edge_to_face[edge][i]][0][2] = p_intersect[2];
        endpoints_count[edge_to_face[edge][i]] += 1; /* update count */
      }
      if (endpoints_count[edge_to_face[edge][i]] == 1) { 
        /* Check if distinct from first endpoint */
        if ( fabs(endpoints[edge_to_face[edge][i]][0][0] - p_intersect[0]) > SC_EPS
            || fabs(endpoints[edge_to_face[edge][i]][0][1] - p_intersect[1]) > SC_EPS
            || fabs(endpoints[edge_to_face[edge][i]][0][2] - p_intersect[2]) > SC_EPS ) 
        {
          /* Record second endpoint */
          endpoints[edge_to_face[edge][i]][1][0] = p_intersect[0];
          endpoints[edge_to_face[edge][i]][1][1] = p_intersect[1];
          endpoints[edge_to_face[edge][i]][1][2] = p_intersect[2];
          endpoints_count[edge_to_face[edge][i]] += 1; /* update count */
        }
      }
    }
  }
}

/** The sphere model refines a spherical mesh based on geodesics. More specifically,
 * squares in the mesh are recursively refined as long as they intersect a geodesic and
 * have refinement level less than the desired resolution. An example application is
 * refining a map of the globe based on coastlines.
 * 
 * A geodesic is represented by its endpoints given in spherical coordinates. We take
 * the convention described here https://en.wikipedia.org/wiki/Spherical_coordinate_system 
 * so that a spherical coordinate is a pair (phi, theta) where:
 *  0 <= theta <= 180     is the polar angle
 *  0 <= phi <= 360   is the azimuth
 * The input is a CSV file where each line
 *    phi1,theta1,phi2,theta2
 * represents a geodesic between endpoints (phi1, theta1) and (phi2, theta2).
 * 
 * \warning Currently geodesics are assumed not to cross faces. Full geodesic support will
 * be implemented soon.
 * 
 * \param[in] resolution maximum refinement level
 * 
*/
p4est_gmt_model_t  *
p4est_gmt_model_sphere_new (int resolution)
{
  p4est_gmt_model_t  *model = P4EST_ALLOC_ZERO (p4est_gmt_model_t, 1);
  p4est_gmt_model_sphere_t *sdata = NULL;
  sphere_geodesic_segment_t *p;
  FILE               *fp;
  /* The following variables get reused for each geodesic we read in. */
  double angular1[2], xyz1[3], rst1[3]; /* angular, cartesian, and tree-local coords respectively*/
  double angular2[2], xyz2[3], rst2[3];
  int which_tree_1, which_tree_2;
  int n_geodesics, capacity; /* capacity is the size of our dynamic array */
  double endpoints[6][2][3]; /* stores endpoints of split geodesics in cartesian coords */
  int endpoints_count[6]; /* counts endpoints of split geodesics assigned to each face */

  /* the sphere model lives on the cube surface reference */
  model->conn = p4est_connectivity_new_cubed ();
  model->output_prefix = "sphere";
  model->model_data = sdata = P4EST_ALLOC (p4est_gmt_model_sphere_t, 1);

  /* Load geodesics */
  n_geodesics = 0; /* Start with 0 and allocate memory as needed */
  capacity = 1;
  p = sdata->geodesics = P4EST_ALLOC (sphere_geodesic_segment_t, capacity);

  fp = fopen("coastlines.csv","r");

  /* Each iteration reads a single geodesic in angular coordinates */
  while (fscanf(fp, "%lf,%lf,%lf,%lf", &angular1[0], &angular1[1], &angular2[0], &angular2[1]) != EOF)
  {
      printf("Loop\n");
      /* Expand capacity when necessary */
      // if (n_geodesics == capacity) {
      //       capacity = capacity * 2;
      //       p = sdata->geodesics = P4EST_REALLOC(sdata->geodesics, sphere_geodesic_segment_t, capacity);
      // }
      while (n_geodesics + 5 >= capacity) { /* We split a geodesic into less than 5 faces */
          capacity = capacity * 2;
          p = sdata->geodesics = P4EST_REALLOC(sdata->geodesics, sphere_geodesic_segment_t, capacity);
      }

      /* Convert to cartesian coordinates */
      angular_to_cube(angular1, xyz1);
      angular_to_cube(angular2, xyz2);
      printf("Scanned 1: %f, %f, %f\n", xyz1[0], xyz1[1], xyz1[2]);
      printf("Scanned 2: %f, %f, %f\n", xyz2[0], xyz2[1], xyz2[2]);
      
      /* Find which face the geodesic endpoints belong to */
      which_tree_1 = point_to_tree(xyz1);
      which_tree_2 = point_to_tree(xyz2);

      if (which_tree_1 == which_tree_2) { /* Geodesic is contained on one face*/
        /* Convert to tree-local coordinates*/
        p4est_geometry_cubed_Y(xyz1, rst1, which_tree_1);
        p4est_geometry_cubed_Y(xyz2, rst2, which_tree_2);

        /* Store geodesic as a sphere_geodesic_segment_t */
        p[n_geodesics].which_tree = which_tree_1;
        p[n_geodesics].p1x = rst1[0];
        p[n_geodesics].p1y = rst1[1];
        p[n_geodesics].p2x = rst2[0];
        p[n_geodesics].p2y = rst2[1];

        n_geodesics++;
      }
      else { /* Geodesic spans multiple faces, so we must split geodesic into segments*/
        /* Reset the mapping {trees : {endpoints}} */
        memset(endpoints_count, 0, sizeof(endpoints_count[0]) * 6);
        memset(endpoints, 0.0, sizeof(endpoints[0][0][0]) * 6 * 2 * 3);

        /* Add endpoints */
        endpoints_count[which_tree_1] += 1;
        endpoints[which_tree_1][0][0] = xyz1[0];
        endpoints[which_tree_1][0][1] = xyz1[1];
        endpoints[which_tree_1][0][2] = xyz1[2];

        endpoints_count[which_tree_2] += 1;
        endpoints[which_tree_2][0][0] = xyz2[0];
        endpoints[which_tree_2][0][1] = xyz2[1];
        endpoints[which_tree_2][0][2] = xyz2[2];

        /* For the 12 edges of the cube compute intersection points and add them to mapping */
        for (int edge = 0; edge < 12; edge++) {
            /* Compute the intersection of geodesic with the given edge*/
            update_endpoints(xyz1, xyz2, edge, endpoints, endpoints_count);
        }

        for (int tree = 0; tree < 6; tree++) {
          if (endpoints_count[tree] == 0) {
            printf("Continuing on tree %d\n", tree);
            continue; /* The geodesic does not cross this cube face */
          }
          if (endpoints_count[tree] == 1) {
            printf("One point on tree %d\n", tree);
            /* The geodesic has an endpoint on the edge of this face (edge case) */
            xyz1[0] = endpoints[tree][0][0];
            xyz1[1] = endpoints[tree][0][1];
            xyz1[2] = endpoints[tree][0][2];
            xyz2[0] = endpoints[tree][0][0];
            xyz2[1] = endpoints[tree][0][1];
            xyz2[2] = endpoints[tree][0][2];
          }
          if (endpoints_count[tree] == 2) {
            printf("Two poins on tree %d\n", tree);
            printf("Point one: %f, %f, %f\n", endpoints[tree][0][0], endpoints[tree][0][1], endpoints[tree][0][2]);
            printf("Point two: %f, %f, %f\n", endpoints[tree][1][0], endpoints[tree][1][1], endpoints[tree][1][2]);
            /* The geodesic has a generic segment on this face */
            /* Set xyz1 and xyz2 to computed endpoints*/
            xyz1[0] = endpoints[tree][0][0];
            xyz1[1] = endpoints[tree][0][1];
            xyz1[2] = endpoints[tree][0][2];
            xyz2[0] = endpoints[tree][1][0];
            xyz2[1] = endpoints[tree][1][1];
            xyz2[2] = endpoints[tree][1][2];
          }

          /* Convert to tree-local coordinates*/
          p4est_geometry_cubed_Y(xyz1, rst1, tree);
          p4est_geometry_cubed_Y(xyz2, rst2, tree);

          /* Store geodesic as a sphere_geodesic_segment_t */
          p[n_geodesics].which_tree = tree;
          p[n_geodesics].p1x = rst1[0];
          p[n_geodesics].p1y = rst1[1];
          p[n_geodesics].p2x = rst2[0];
          p[n_geodesics].p2y = rst2[1];

          n_geodesics++;
        }
      }
  }

  /* Free extra capacity */
  p = sdata->geodesics = P4EST_REALLOC(sdata->geodesics, sphere_geodesic_segment_t, n_geodesics);
  
  fclose(fp); /* Finished loading geodesics */
  sdata->num_geodesics = model->M = n_geodesics; /* Set final geodesic count */
  printf("n_geodesics %d\n", n_geodesics);

  /* Assign resolution, intersector and destructor */
  sdata->resolution = resolution;
  model->destroy_data = model_sphere_destroy_data;
  model->intersect = model_sphere_intersect;

  /* Assign geometry */
  /* Note: the problem with the following is that it allocates memory externally,
   * rather than in sgeom. 
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