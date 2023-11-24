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

#include "gmt_models.c"

/* Convert from angular coordinates to corresponding point on cube face*/
static void angular_to_cube(const double angular[2], double xyz[3])
{
  double inf_norm;
  double phi, theta;

  /* Convert to radians*/
  phi = angular[0] * M_PI / 180.0;
  theta = angular[1] * M_PI / 180.0;

  xyz[0] = sin(theta) * cos(phi);
  xyz[1] = sin(theta) * sin(phi);
  xyz[2] = cos(theta);

  inf_norm = fmax(fabs(xyz[0]), fmax(fabs(xyz[1]), fabs(xyz[2]))) * 2.0;
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
static int point_to_tree(const double xyz[3])
{
  if (fabs(xyz[0] + 0.5) < SC_EPS)
    return 2;
  if (fabs(xyz[0] - 0.5) < SC_EPS)
    return 5;
  if (fabs(xyz[1] + 0.5) < SC_EPS)
    return 4;
  if (fabs(xyz[1] - 0.5) < SC_EPS)
    return 1;
  if (fabs(xyz[2] + 0.5) < SC_EPS)
    return 0;
  if (fabs(xyz[2] - 0.5) < SC_EPS)
    return 3;
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
static void p4est_geometry_cubed_Y(const double xyz[3], double rst[3], int which_tree)
{
  rst[2] = 0.0;

  /* align center with origin */
  switch (which_tree)
  {
  case 0:
    rst[0] = xyz[1] + 0.5;
    rst[1] = xyz[0] + 0.5;
    break;
  case 1:
    rst[0] = xyz[2] + 0.5;
    rst[1] = xyz[0] + 0.5;
    break;
  case 2:
    rst[0] = xyz[2] + 0.5;
    rst[1] = xyz[1] + 0.5;
    break;
  case 3:
    rst[0] = xyz[0] + 0.5;
    rst[1] = xyz[1] + 0.5;
    break;
  case 4:
    rst[0] = xyz[0] + 0.5;
    rst[1] = xyz[2] + 0.5;
    break;
  case 5:
    rst[0] = xyz[1] + 0.5;
    rst[1] = xyz[2] + 0.5;
    break;
  default:
    break;
  }
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
  double A[3][3];        /* Matrix we are inverting */
  double cofactor[3][3]; /* Cofactor matrix */
  double det_A, det_A_inv;
  double x[3];

  A[0][0] = v1[0];
  A[1][0] = v1[1];
  A[2][0] = v1[2];
  A[0][1] = v2[0];
  A[1][1] = v2[1];
  A[2][1] = v2[2];
  A[0][2] = p1[0] - p2[0];
  A[1][2] = p1[1] - p2[1];
  A[2][2] = p1[2] - p2[2];

  /* Compute minors */
  for (int r = 0; r < 3; r++)
  {
    for (int c = 0; c < 3; c++)
    {
      cofactor[r][c] = A[(r + 1) % 3][(c + 1) % 3] * A[(r + 2) % 3][(c + 2) % 3] - A[(r + 1) % 3][(c + 2) % 3] * A[(r + 2) % 3][(c + 1) % 3];
    }
  }

  /* Compute the determinant by Laplace expansion along the first column */
  det_A = 0;
  for (int r = 0; r < 3; r++)
  {
    det_A += A[r][0] * cofactor[r][0];
  }
  det_A_inv = 1 / det_A; /* TODO: What if det_A = 0?*/

  /* Multiply p1 with inverse of A */
  for (int i = 0; i < 3; i++)
  {
    /* Compute x[i] */
    x[i] = 0;
    for (int j = 0; j < 3; j++)
    {
      x[i] += cofactor[j][i] * p1[j];
    }
    x[i] *= det_A_inv;
  }

  if (x[0] < 0 || x[1] < 0 || x[2] < 0 || x[2] > 1)
  {
    return 0; /* Invalid intersection */
  }

  p_intersect[0] = x[0] * v1[0] + x[1] * v2[0];
  p_intersect[1] = x[0] * v1[1] + x[1] * v2[1];
  p_intersect[2] = x[0] * v1[2] + x[1] * v2[2];
  return 1; /* Valid intersection */
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
      {0, 1},
      {0, 2},
      {0, 4},
      {0, 5},
      {1, 2},
      {1, 3},
      {1, 5},
      {2, 3},
      {2, 4},
      {3, 4},
      {3, 5},
      {4, 5}};

  /* Cube edge endpoint coordinates */
  const double edge_endpoints[12][2][3] = {
      {{-0.5, 0.5, -0.5}, {0.5, 0.5, -0.5}},   /* 0,1 edge */
      {{-0.5, -0.5, -0.5}, {-0.5, 0.5, -0.5}}, /* 0,2 edge */
      {{-0.5, -0.5, -0.5}, {0.5, -0.5, -0.5}}, /* 0,4 edge*/
      {{0.5, -0.5, -0.5}, {0.5, 0.5, -0.5}},   /* 0,5 edge */
      {{-0.5, 0.5, -0.5}, {-0.5, 0.5, 0.5}},   /* 1,2 edge */
      {{-0.5, 0.5, 0.5}, {0.5, 0.5, 0.5}},     /* 1,3 edge */
      {{0.5, 0.5, -0.5}, {0.5, 0.5, 0.5}},     /* 1,5 edge */
      {{-0.5, -0.5, 0.5}, {-0.5, 0.5, 0.5}},   /* 2,3 edge */
      {{-0.5, -0.5, -0.5}, {-0.5, -0.5, 0.5}}, /* 2,4 edge */
      {{-0.5, -0.5, 0.5}, {0.5, -0.5, 0.5}},   /* 3,4 edge */
      {{0.5, -0.5, 0.5}, {0.5, 0.5, 0.5}},     /* 3,5 edge */
      {{0.5, -0.5, -0.5}, 0.5, -0.5, 0.5}      /* 4,5 edge*/
  };

  detected = cone_line_intersection(xyz1, xyz2, edge_endpoints[edge][0],
                                    edge_endpoints[edge][1], p_intersect);

  if (detected == 1)
  {
    for (int i = 0; i < 2; i++)
    {
      if (endpoints_count[edge_to_face[edge][i]] == 0)
      {
        /* Record first endpoint */
        endpoints[edge_to_face[edge][i]][0][0] = p_intersect[0];
        endpoints[edge_to_face[edge][i]][0][1] = p_intersect[1];
        endpoints[edge_to_face[edge][i]][0][2] = p_intersect[2];
        endpoints_count[edge_to_face[edge][i]] += 1; /* update count */
      }
      if (endpoints_count[edge_to_face[edge][i]] == 1)
      {
        /* Check if distinct from first endpoint */
        if (fabs(endpoints[edge_to_face[edge][i]][0][0] - p_intersect[0]) > SC_EPS || fabs(endpoints[edge_to_face[edge][i]][0][1] - p_intersect[1]) > SC_EPS || fabs(endpoints[edge_to_face[edge][i]][0][2] - p_intersect[2]) > SC_EPS)
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

/** Load geodesics from coastlines.csv, convert to Cartesian coordinates, split into
 *  segments corresponding to the trees in the cubed connectivity, then write to array
 *  of type geodesic_segment_t.
 * 
 * The input is a CSV file where each line
 *    phi1,theta1,phi2,theta2
 * represents a geodesic between endpoints (phi1, theta1) and (phi2, theta2).
 *
*/
int main(int argc, char **argv)
{
  const char *usage;
  FILE *geodesic_file;
  FILE *fp;
  sphere_geodesic_segment_t *geodesics;
  /* The following variables get reused for each geodesic we read in. */
  double angular1[2], xyz1[3], rst1[3]; /* angular, cartesian, and tree-local coords respectively*/
  double angular2[2], xyz2[3], rst2[3];
  int which_tree_1, which_tree_2;
  int n_geodesics, capacity; /* capacity is the size of our dynamic array */
  double endpoints[6][2][3]; /* stores endpoints of split geodesics in cartesian coords */
  int endpoints_count[6];    /* counts endpoints of split geodesics assigned to each face */

  usage = "Arguments: <input.csv>\n";

  if (argc != 2) {
    P4EST_GLOBAL_LERROR (usage);
    sc_abort_collective ("incorrect number of arguments");
  } 

  /* Load geodesics */
  n_geodesics = 0; /* Start with 0 and allocate memory as needed */
  capacity = 1;
  geodesics = P4EST_ALLOC(sphere_geodesic_segment_t, capacity);

  fp = fopen(argv[1], "r");

  if (fp == NULL) {
    P4EST_GLOBAL_LERROR (usage);
    sc_abort_collective ("<input.csv> not found");
  }

  /* Each iteration reads a single geodesic in angular coordinates */
  while (fscanf(fp, "%lf,%lf,%lf,%lf", &angular1[0], &angular1[1], &angular2[0], &angular2[1]) != EOF)
  {
    while (n_geodesics + 5 >= capacity)
    { /* We split a geodesic into less than 5 faces */
      capacity = capacity * 2;
      geodesics = P4EST_REALLOC(geodesics, sphere_geodesic_segment_t, capacity);
    }

    /* Convert to cartesian coordinates */
    angular_to_cube(angular1, xyz1);
    angular_to_cube(angular2, xyz2);

    /* Find which face the geodesic endpoints belong to */
    which_tree_1 = point_to_tree(xyz1);
    which_tree_2 = point_to_tree(xyz2);

    if (which_tree_1 == which_tree_2)
    { /* Geodesic is contained on one face*/
      /* Convert to tree-local coordinates*/
      p4est_geometry_cubed_Y(xyz1, rst1, which_tree_1);
      p4est_geometry_cubed_Y(xyz2, rst2, which_tree_2);

      /* Store geodesic as a sphere_geodesic_segment_t */
      geodesics[n_geodesics].which_tree = which_tree_1;
      geodesics[n_geodesics].p1x = rst1[0];
      geodesics[n_geodesics].p1y = rst1[1];
      geodesics[n_geodesics].p2x = rst2[0];
      geodesics[n_geodesics].p2y = rst2[1];

      n_geodesics++;
    }
    else
    { /* Geodesic spans multiple faces, so we must split geodesic into segments*/
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
      for (int edge = 0; edge < 12; edge++)
      {
        /* Compute the intersection of geodesic with the given edge*/
        update_endpoints(xyz1, xyz2, edge, endpoints, endpoints_count);
      }
      for (int tree = 0; tree < 6; tree++)
      {
        if (endpoints_count[tree] == 0)
        {
          continue; /* The geodesic does not cross this cube face */
        }
        if (endpoints_count[tree] == 1)
        {
          /* The geodesic has an endpoint on the edge of this face (edge case) */
          xyz1[0] = endpoints[tree][0][0];
          xyz1[1] = endpoints[tree][0][1];
          xyz1[2] = endpoints[tree][0][2];
          xyz2[0] = endpoints[tree][0][0];
          xyz2[1] = endpoints[tree][0][1];
          xyz2[2] = endpoints[tree][0][2];
        }
        if (endpoints_count[tree] == 2)
        {
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
        geodesics[n_geodesics].which_tree = tree;
        geodesics[n_geodesics].p1x = rst1[0];
        geodesics[n_geodesics].p1y = rst1[1];
        geodesics[n_geodesics].p2x = rst2[0];
        geodesics[n_geodesics].p2y = rst2[1];

        n_geodesics++;
      }
    }
  }

  /* Free extra capacity */
  geodesics = P4EST_REALLOC(geodesics, sphere_geodesic_segment_t, n_geodesics);

  fclose(fp);                                    /* Finished loading geodesics */
  printf("n_geodesics %d\n", n_geodesics);

  /* Write geodesics to disk */
  geodesic_file = sc_fopen("geodesics", "w", "opening geodesics file");
  /* Write n_geodesics */
  sc_fwrite(&n_geodesics, sizeof(int), 1, geodesic_file, "writing n_geodesics");
  /* Write geodesics */
  sc_fwrite(geodesics, sizeof(sphere_geodesic_segment_t), n_geodesics, 
              geodesic_file, "writing geodesics");
  fclose(geodesic_file); /* Finished writing geodesics*/

  return 0;
}