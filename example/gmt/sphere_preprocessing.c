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

/** \file sphere_preprocessing.c
 *
 * Preprocessing for the sphere model. The main function reads in a list of
 * geodesics (straight paths on the sphere) and splits them into segments, such
 * that each segment is contained entirely in a single cube face.
 * 
 * Usage: p4est_sphere_preprocessing <input.csv> <output_file_name>
 * 
 * Here <input.csv> is a CSV file where each line
 *    phi1,theta1,phi2,theta2
 * represents the unique geodesic between endpoints (phi1, theta1) and
 * (phi2, theta2). We require that the endpoints are not antipodal, as in this
 * case there is not a unique geodesic between them.
 * 
 * Endpoints are given in spherical coordinates. We take the convention
 * described here:
 * https://en.wikipedia.org/wiki/Spherical_coordinate_system.
 * That is, a spherical coordinate is a pair (phi, theta) where:
 *  0 <= theta <= 180 is the polar angle
 *  0 <= phi <= 360   is the azimuth
 * 
 * An example input file sphere_hello_world.csv is included.
 * 
 * Geodesic segments are defined in the following way: 
 * 
 * Let v1 and v2 be two distinct points on the sphere that are not antipodal,
 * given in cartesian coordinates. The (linear) cone spanned by v1 and v2 is
 * the set of points of the form a*v1 + b*v2 with a and b non-negative. The 
 * geodesic between v1 and v2 is exactly the intersection of the sphere with
 * the cone spanned by v1 and v2.
 * 
 * We map the sphere to the surface of a cube by scaling each point by the
 * inverse of its uniform norm. Under this transformation the image of a
 * the geodesic between v1 and v2 is the intersection of the cube surface
 * with the cone spanned by v1 and v2. On a single face of a cube this is
 * simply the intersection of a square with the cone, which is either a line
 * segment or empty. Thus we can split the geodesic into segments for each
 * face of the cube that it intersects, and each segment is just a line 
 * segment. These line segments are represents by the face they belong to and
 * their endpoints in the local coordinate system of that face.
 * 
 * An endpoint of a geodesic segment is either an endpoint of the entire
 * geodesic, or it lies on the edge of its face. Thus to segment geodesics it
 * suffices to compute their intersections with the 12 edges of the
 * connectivity \ref p4est_connectivity_new_cubed. This is performed by
 * \ref cone_line_intersection which solves a system of linear equations to
 * find the intersection of a cone and a line segment.
 */

#include "gmt_models.h"

/** Convert from angular coordinates to corresponding point on cube face */
static void
angular_to_cube (const double angular[2], double xyz[3])
{
  double              inf_norm_inv;
  double              phi, theta;

  /* Convert to radians */
  phi = angular[0] * M_PI / 180.0;
  theta = angular[1] * M_PI / 180.0;

  xyz[0] = sin (theta) * cos (phi);
  xyz[1] = sin (theta) * sin (phi);
  xyz[2] = cos (theta);

  inf_norm_inv =
    0.5 / SC_MAX (fabs (xyz[0]), SC_MAX (fabs (xyz[1]), fabs (xyz[2])));
  xyz[0] *= inf_norm_inv;
  xyz[1] *= inf_norm_inv;
  xyz[2] *= inf_norm_inv;

  return;
}

/** Which tree in the connectivity \ref p4est_connectivity_new_cubed does a
 *  given point belong to.
 * 
 * See \ref p4est_connectivity_new_cubed for a description of cube face 
 * numbering. The cube is embedded in R^3 via the vertex coordinates
 * coordinates specified in this connectivity, and then translated so that it
 * is centred at the origin. 
 *
 * \note Ties are decided arbitrarily. 
 *
 * \param[in] xyz  cartesian coordinates on surface of the cube [-0.5,0.5]^3
 */
static int
point_to_tree (const double xyz[3])
{
  if (fabs (xyz[0] + 0.5) < SC_EPS)
    return 2;
  if (fabs (xyz[0] - 0.5) < SC_EPS)
    return 5;
  if (fabs (xyz[1] + 0.5) < SC_EPS)
    return 4;
  if (fabs (xyz[1] - 0.5) < SC_EPS)
    return 1;
  if (fabs (xyz[2] + 0.5) < SC_EPS)
    return 0;
  if (fabs (xyz[2] - 0.5) < SC_EPS)
    return 3;
  SC_ABORT_NOT_REACHED ();
}

/** Coordinate transformation from the surface of cube [-0.5,0.5]^3 to 
 * tree-local coordinates. After translation this is the inverse of the 
 * coordinate transform created by \ref p4est_geometry_new_connectivity
 * applied to the connectivity \ref p4est_connectivity_new_cubed.
 *
 * \param[in]  xyz  cartesian coordinates in physical space
 * \param[in]  which_tree face of cube
 * \param[out] rst  tree-local reference coordinates : [0,1]^3
 *
 */
static void
p4est_geometry_cubed_Y (const double xyz[3], double rst[3],
                        p4est_topidx_t which_tree)
{
  P4EST_ASSERT (0 <= which_tree && which_tree <= 5);
  P4EST_ASSERT (-0.5 <= xyz[0] && xyz[0] <= 0.5);
  P4EST_ASSERT (-0.5 <= xyz[1] && xyz[1] <= 0.5);
  P4EST_ASSERT (-0.5 <= xyz[2] && xyz[2] <= 0.5);
  P4EST_ASSERT (fabs (fabs (xyz[0]) - 0.5) < SC_1000_EPS
                || fabs (fabs (xyz[1]) - 0.5) < SC_1000_EPS
                || fabs (fabs (xyz[2]) - 0.5) < SC_1000_EPS);

  rst[2] = 0.0;                 /* third coordinate is never used by p4est */

  /* align center with origin */
  switch (which_tree) {
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
    SC_ABORT_NOT_REACHED ();
  }
}

/** Solves a system of linear equations to find the intersection of a cone
 * and a line segment in R3. 
 * 
 * Returns 1 if the cone spanned by v1 and v2 intersects the line segment
 * between p1 and p2. If an intersection is detected then p_intersect is
 * set to the computed intersection point.
 * 
 * We assume that v1 and v2 are not colinear.
 * 
 * This is used to calculate where a geodesic intersects the edges of the 
 * connectivity \ref p4est_connectivity_new_cubed so that we can determine the
 * segment of the geodesic lying on a particular face. In this case v1 and v2
 * are the endpoints of the geodesic, and p1 and p2 are the vertices of the
 * particular edge we are interested in calculating the intersection with.
 * 
 * \param[in] v1 generator of the cone
 * \param[in] v2 generator of the cone
 * \param[in] p1 line segment endpoint
 * \param[in] p2 line segment endpoint
 * \param[out] p_intersect computed intersection point
 */
static int
cone_line_intersection (const double v1[3], const double v2[3],
                        const double p1[3], const double p2[3],
                        double p_intersect[3])
{
  /** We solve the matrix equation (v1, v2, p1-p2) x = p1 by inverting 
   * (v1, v2, p1-p2) */
  double              A[3][3];  /* Matrix we are inverting */
  double              cofactor[3][3];   /* Cofactor matrix */
  double              det_A, det_A_inv;
  double              x[3];     /* Solving for x */

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
  for (int r = 0; r < 3; r++) {
    for (int c = 0; c < 3; c++) {
      cofactor[r][c] =
        A[(r + 1) % 3][(c + 1) % 3] * A[(r + 2) % 3][(c + 2) % 3] -
        A[(r + 1) % 3][(c + 2) % 3] * A[(r + 2) % 3][(c + 1) % 3];
    }
  }

  /* Compute the determinant by Laplace expansion along the first column */
  det_A = 0;
  for (int r = 0; r < 3; r++) {
    det_A += A[r][0] * cofactor[r][0];
  }

  /** If the determinant is zero then the line segment is parallel to the 
   * cone. We count this case as not intersecting.
  */
  if (fabs (det_A) < SC_1000_EPS) {
    return 0;
  }

  det_A_inv = 1 / det_A;

  /* Multiply p1 with inverse of A */
  for (int i = 0; i < 3; i++) {
    /* Compute x[i] */
    x[i] = 0;
    for (int j = 0; j < 3; j++) {
      x[i] += cofactor[j][i] * p1[j];
    }
    x[i] *= det_A_inv;
  }

  if (x[0] < 0 || x[1] < 0 || x[2] < 0 || x[2] > 1) {
    return 0;                   /* No intersection */
  }

  p_intersect[0] = x[0] * v1[0] + x[1] * v2[0];
  p_intersect[1] = x[0] * v1[1] + x[1] * v2[1];
  p_intersect[2] = x[0] * v1[2] + x[1] * v2[2];
  return 1;                     /* Valid intersection */
}

/** Used to clamp coordinates so that they lie in [-0.5,0.5]^3 */
static double
clamp (double x)
{
  if (x < -0.5) {
    return -0.5;
  }
  else if (x > 0.5) {
    return 0.5;
  }
  return x;
}

/** If the geodesic between xyz1 and xyz2 intersects the given edge then add 
 * this intersection point to endpoints and increment the corresponding entry
 * in endpoints_count. To deal with edge cases coming from corners we should
 * only update if the new intersection point is distinct to previously seen 
 * intersection points.
 */
static void
update_endpoints (const double xyz1[3], const double xyz2[3], int edge,
                  double endpoints[6][2][3], int endpoints_count[6])
{
  double              p_intersect[3];
  int                 detected;

  /* Which cube faces are adjacent to the given edge */
  const int           edge_to_face[12][2] = {
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
    {4, 5}
  };

  /* Cube edge endpoint coordinates */
  const double        edge_endpoints[12][2][3] = {
    {{-0.5, 0.5, -0.5}, {0.5, 0.5, -0.5}},      /* 0,1 edge */
    {{-0.5, -0.5, -0.5}, {-0.5, 0.5, -0.5}},    /* 0,2 edge */
    {{-0.5, -0.5, -0.5}, {0.5, -0.5, -0.5}},    /* 0,4 edge */
    {{0.5, -0.5, -0.5}, {0.5, 0.5, -0.5}},      /* 0,5 edge */
    {{-0.5, 0.5, -0.5}, {-0.5, 0.5, 0.5}},      /* 1,2 edge */
    {{-0.5, 0.5, 0.5}, {0.5, 0.5, 0.5}},        /* 1,3 edge */
    {{0.5, 0.5, -0.5}, {0.5, 0.5, 0.5}},        /* 1,5 edge */
    {{-0.5, -0.5, 0.5}, {-0.5, 0.5, 0.5}},      /* 2,3 edge */
    {{-0.5, -0.5, -0.5}, {-0.5, -0.5, 0.5}},    /* 2,4 edge */
    {{-0.5, -0.5, 0.5}, {0.5, -0.5, 0.5}},      /* 3,4 edge */
    {{0.5, -0.5, 0.5}, {0.5, 0.5, 0.5}},        /* 3,5 edge */
    {{0.5, -0.5, -0.5}, {0.5, -0.5, 0.5}}       /* 4,5 edge */
  };

  /* Solve for intersection */
  detected = cone_line_intersection (xyz1, xyz2, edge_endpoints[edge][0],
                                     edge_endpoints[edge][1], p_intersect);

  if (detected) {

    /* Correct for numerical instabilities */
    for (int i = 0; i < 3; i++) {
      if (edge_endpoints[edge][0][i] == edge_endpoints[edge][1][i]) {
        /* two of the three coordinates should be exactly 0.5 or -0.5 */
        p_intersect[i] = edge_endpoints[edge][0][i];
      }
      else {
        /* the final coordinate is free to range between -0.5 and 0.5 */
        p_intersect[i] = clamp (p_intersect[i]);
      }
    }

    /* Record endpoints */
    for (int i = 0; i < 2; i++) {
      if (endpoints_count[edge_to_face[edge][i]] == 0) {
        /* Record first endpoint */
        endpoints[edge_to_face[edge][i]][0][0] = p_intersect[0];
        endpoints[edge_to_face[edge][i]][0][1] = p_intersect[1];
        endpoints[edge_to_face[edge][i]][0][2] = p_intersect[2];
        endpoints_count[edge_to_face[edge][i]] += 1;    /* update count */
      }
      if (endpoints_count[edge_to_face[edge][i]] == 1) {
        /* Check if distinct from first endpoint */
        if (fabs (endpoints[edge_to_face[edge][i]][0][0] - p_intersect[0]) >
            SC_EPS
            || fabs (endpoints[edge_to_face[edge][i]][0][1] -
                     p_intersect[1]) > SC_EPS
            || fabs (endpoints[edge_to_face[edge][i]][0][2] -
                     p_intersect[2]) > SC_EPS) {
          /* Record second endpoint */
          endpoints[edge_to_face[edge][i]][1][0] = p_intersect[0];
          endpoints[edge_to_face[edge][i]][1][1] = p_intersect[1];
          endpoints[edge_to_face[edge][i]][1][2] = p_intersect[2];
          endpoints_count[edge_to_face[edge][i]] += 1;  /* update count */
        }
      }
    }
  }
}

/** Split geodesics in input into segments and store these in an array
 * 
 * \param[in] input csv file containing geodesics
 * \param[out] geodesics_out array storing geodesic segments
 * \param[out] n_geodesics_out number of geodesic segments 
*/
static int
compute_geodesic_splits (size_t *n_geodesics_out,
                         p4est_gmt_sphere_geoseg_t **geodesics_out,
                         FILE *input)
{
  int                 scanned;
  /* capacity is the size of our dynamic array */
  size_t              n_geodesics, capacity;
  p4est_gmt_sphere_geoseg_t *geodesics;
  /* The following variables get reused for each geodesic we read in. */
  /* angular, cartesian, and tree-local coords respectively */
  double              angular1[2], xyz1[3], rst1[3];
  double              angular2[2], xyz2[3], rst2[3];
  p4est_topidx_t      which_tree_1, which_tree_2;
  /* stores endpoints of split geodesics in cartesian coords */
  double              endpoints[6][2][3];
  /* counts endpoints of split geodesics assigned to each face */
  int                 endpoints_count[6];

  /* Start with 0 and allocate memory as needed */
  n_geodesics = 0;
  capacity = 1;
  geodesics = P4EST_ALLOC (p4est_gmt_sphere_geoseg_t, capacity);

  /* Each iteration reads a single geodesic in angular coordinates */
  while ((scanned = fscanf
          (input, "%lf,%lf,%lf,%lf", &angular1[0], &angular1[1], &angular2[0],
           &angular2[1])) != EOF) {
    if (scanned != 4) {
      P4EST_GLOBAL_LERROR ("Badly formatted input\n"
                           "Expected csv with 4 doubles per line\n");
      return 1;
    }

    /* Ensure we will have enough capacity to store our split geodesic */
    while (n_geodesics + 5 >= capacity) {
      capacity = capacity * 2;
      geodesics =
        P4EST_REALLOC (geodesics, p4est_gmt_sphere_geoseg_t, capacity);
    }

    /* Convert to cartesian coordinates */
    angular_to_cube (angular1, xyz1);
    angular_to_cube (angular2, xyz2);

    /* Check that endpoints are not antipodal */
    if (fabs (xyz1[0] + xyz2[0]) < SC_EPS
        && fabs (xyz1[1] + xyz2[1]) < SC_EPS
        && fabs (xyz1[2] + xyz2[2]) < SC_EPS) {
      P4EST_GLOBAL_LERROR
        ("Endpoints of a geodesic should not be antipodal\n");
      return 1;
    }

    /* Find which face the geodesic endpoints belong to */
    which_tree_1 = point_to_tree (xyz1);
    which_tree_2 = point_to_tree (xyz2);

    if (which_tree_1 == which_tree_2) {
      /* Geodesic is contained on one face, no splitting required */

      /* Convert to tree-local coordinates */
      p4est_geometry_cubed_Y (xyz1, rst1, which_tree_1);
      p4est_geometry_cubed_Y (xyz2, rst2, which_tree_2);

      /* Store geodesic as a p4est_gmt_sphere_geoseg_t */
      geodesics[n_geodesics].which_tree = which_tree_1;
      geodesics[n_geodesics].p1x = rst1[0];
      geodesics[n_geodesics].p1y = rst1[1];
      geodesics[n_geodesics].p2x = rst2[0];
      geodesics[n_geodesics].p2y = rst2[1];
      /* initialise padding so valgrind is happy */
      geodesics[n_geodesics].pad4 = 0;

      n_geodesics++;
    }
    else {
      /* Geodesic spans multiple faces, so we must split it into segments */

      /* Reset the mapping [tree -> {endpoints}] */
      memset (endpoints_count, 0, sizeof (endpoints_count[0]) * 6);
      memset (endpoints, 0.0, sizeof (endpoints[0][0][0]) * 6 * 2 * 3);

      /* Add initial endpoints */
      endpoints_count[which_tree_1] += 1;
      endpoints[which_tree_1][0][0] = xyz1[0];
      endpoints[which_tree_1][0][1] = xyz1[1];
      endpoints[which_tree_1][0][2] = xyz1[2];

      endpoints_count[which_tree_2] += 1;
      endpoints[which_tree_2][0][0] = xyz2[0];
      endpoints[which_tree_2][0][1] = xyz2[1];
      endpoints[which_tree_2][0][2] = xyz2[2];

      /* Add endpoints for intersections of the geodesic with cube edges */
      for (int edge = 0; edge < 12; edge++) {
        update_endpoints (xyz1, xyz2, edge, endpoints, endpoints_count);
      }

      /* Record segments given by the computed endpoints */
      for (p4est_topidx_t tree = 0; tree < 6; tree++) {
        if (endpoints_count[tree] == 0) {
          /* The geodesic does not cross this cube face */
          continue;
        }
        if (endpoints_count[tree] == 1) {
          /* The geodesic has an endpoint on the edge of this face */
          /* This is an edge case and probably never happens */
          /* Set xyz1 = xyz2 to this point */
          xyz1[0] = xyz2[0] = endpoints[tree][0][0];
          xyz1[1] = xyz2[1] = endpoints[tree][0][1];
          xyz1[2] = xyz2[2] = endpoints[tree][0][2];
        }
        if (endpoints_count[tree] == 2) {
          /* The geodesic has a generic segment on this face */
          /* Set xyz1 and xyz2 to computed endpoints */
          xyz1[0] = endpoints[tree][0][0];
          xyz1[1] = endpoints[tree][0][1];
          xyz1[2] = endpoints[tree][0][2];
          xyz2[0] = endpoints[tree][1][0];
          xyz2[1] = endpoints[tree][1][1];
          xyz2[2] = endpoints[tree][1][2];
        }

        /* Convert to tree-local coordinates */
        p4est_geometry_cubed_Y (xyz1, rst1, tree);
        p4est_geometry_cubed_Y (xyz2, rst2, tree);

        /* Store geodesic segment */
        geodesics[n_geodesics].which_tree = tree;
        geodesics[n_geodesics].p1x = rst1[0];
        geodesics[n_geodesics].p1y = rst1[1];
        geodesics[n_geodesics].p2x = rst2[0];
        geodesics[n_geodesics].p2y = rst2[1];
        /* initialise padding so valgrind is happy */
        geodesics[n_geodesics].pad4 = 0;

        n_geodesics++;
      }
    }
  }

  /* Free extra capacity */
  geodesics =
    P4EST_REALLOC (geodesics, p4est_gmt_sphere_geoseg_t, n_geodesics);

  *n_geodesics_out = n_geodesics;
  *geodesics_out = geodesics;

  return 0;
}

/** Load geodesics from input csv, convert to Cartesian coordinates, split into
 *  segments corresponding to the trees in the cubed connectivity, then write
 *  to array of type p4est_gmt_sphere_geoseg_t.
*/
int
main (int argc, char **argv)
{
  int                 progerr;
  int                 close_err;
  const char         *usage;
  FILE               *output = NULL;
  FILE               *input = NULL;
  p4est_gmt_sphere_geoseg_t *geodesics = NULL;
  size_t              n_geodesics;
  size_t              nwritten;
  sc_MPI_Comm         mpicomm;
  int                 mpiret;
  int                 rank, num_procs;

  /* initialize MPI */
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  mpicomm = sc_MPI_COMM_WORLD;

  /* Get rank and number of processes */
  mpiret = sc_MPI_Comm_size (mpicomm, &num_procs);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (mpicomm, &rank);
  SC_CHECK_MPI (mpiret);

  if (num_procs > 1 && rank == 0) {
    P4EST_GLOBAL_INFOF ("Warning: you ran this script with %d processes. "
                      "Currently no distributed version is offered. "
                      "All work is being performed on rank 0.\n",
                      num_procs);
  }

  /* all work is performed on rank 0 */
  if (rank == 0) {
    usage = "Arguments: <input.csv> <output>\n";
    progerr = 0;

    if (argc != 3) {
      P4EST_GLOBAL_LERROR ("Incorrect number of arguments\n");
      P4EST_GLOBAL_LERROR (usage);
      progerr = 1;
    }

    /* open input file for reading */
    if (!progerr) {
      input = fopen (argv[1], "r");
      if (input == NULL) {
        P4EST_GLOBAL_LERRORF ("Could not find input file: %s\n", argv[1]);
        P4EST_GLOBAL_LERROR (usage);
        progerr = 1;
      }
    }

    /* split geodesics into segments */
    if (!progerr) {
      progerr = compute_geodesic_splits (&n_geodesics, &geodesics, input);
      if (progerr) {
        P4EST_GLOBAL_LERRORF ("Failed parsing input file: %s\n", argv[1]);
      }
    }

    /* close input file */
    if (input != NULL) {
      close_err = fclose (input);
      if (close_err) {
        P4EST_GLOBAL_LERRORF ("Error closing input file: %s\n", argv[1]);
        progerr = 1;
      }
    }

    /* open output file for writing */
    if (!progerr) {
      output = fopen (argv[2], "w");
      if (output == NULL) {
        progerr = 1;
        P4EST_GLOBAL_LERRORF ("File open fail: %s\n", argv[2]);
      }
    }

    /* write number of geodesic segments */
    if (!progerr) {
      nwritten = fwrite (&n_geodesics, sizeof (size_t), 1, output);
      if (nwritten != 1) {
        progerr = 1;
        P4EST_GLOBAL_LERRORF ("File write fail: "
                              "writing n_geodesics to %s\n", argv[2]);
      }
    }

    /* write geodesic segments */
    if (!progerr) {
      nwritten = fwrite (geodesics, sizeof (p4est_gmt_sphere_geoseg_t),
                        n_geodesics, output);
      if (nwritten != n_geodesics) {
        progerr = 1;
        P4EST_GLOBAL_LERRORF ("File write fail: "
                              "writing geodesics to %s\n", argv[2]);
      }
    }

    /* free written geodesics */
    P4EST_FREE (geodesics);

    /* close output file */
    if (output != NULL) {
      close_err = fclose (output);
      if (close_err) {
        P4EST_GLOBAL_LERRORF ("Error closing output file: %s\n", argv[2]);
        progerr = 1;
      }
    }
  }

  /* broadcast rank 0 error state */
  mpiret = sc_MPI_Bcast(&progerr, sizeof (int), sc_MPI_BYTE, 0, mpicomm);
  SC_CHECK_MPI (mpiret);

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return progerr ? EXIT_FAILURE : EXIT_SUCCESS;
}
