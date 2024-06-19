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

/**
 * \file p4est_geometry.c
 * We provide several transformations for reference.
 * Please implement p4est_geometry_t as you see fit.
 */

#ifndef P4_TO_P8
#include <p4est_bits.h>
#include <p4est_geometry.h>
#else
#include <p8est_bits.h>
#include <p8est_geometry.h>
#endif

#ifndef P4_TO_P8

typedef enum
{
  P4EST_GEOMETRY_BUILTIN_MAGIC = 0x20F2F8DE,
  P4EST_GEOMETRY_BUILTIN_ICOSAHEDRON,
  P4EST_GEOMETRY_BUILTIN_SHELL2D,
  P4EST_GEOMETRY_BUILTIN_DISK2D,
  P4EST_GEOMETRY_BUILTIN_SPHERE2D,
  P4EST_GEOMETRY_LAST
}
p4est_geometry_builtin_type_t;

typedef struct p4est_geometry_builtin_icosahedron
{
  p4est_geometry_builtin_type_t type;
  double              R;        /* sphere radius */
}
p4est_geometry_builtin_icosahedron_t;

typedef struct p4est_geometry_builtin_shell2d
{
  p4est_geometry_builtin_type_t type;
  double              R2, R1;
  double              R2byR1, R1sqrbyR2, Rlog;
}
p4est_geometry_builtin_shell2d_t;

typedef struct p4est_geometry_builtin_disk2d
{
  p4est_geometry_builtin_type_t type;
  double              R0, R1;
  double              R1byR0, R0sqrbyR1, R0log;
  double              Clength;
}
p4est_geometry_builtin_disk2d_t;

typedef struct p4est_geometry_builtin_sphere2d
{
  p4est_geometry_builtin_type_t type;
  double              R;
}
p4est_geometry_builtin_sphere2d_t;

typedef struct p4est_geometry_builtin
{
  /** The geom member needs to come first; we cast to p4est_geometry_t * */
  p4est_geometry_t    geom;
  union
  {
    p4est_geometry_builtin_type_t type;
    p4est_geometry_builtin_icosahedron_t icosahedron;
    p4est_geometry_builtin_shell2d_t shell2d;
    p4est_geometry_builtin_disk2d_t disk2d;
    p4est_geometry_builtin_sphere2d_t sphere2d;
  }
  p;
}
p4est_geometry_builtin_t;

#endif /* !P4_TO_P8 */

void
p4est_geometry_destroy (p4est_geometry_t * geom)
{
  if (geom->destroy != NULL) {
    geom->destroy (geom);
  }
  else {
    P4EST_FREE (geom);
  }
}

void
p4est_geometry_connectivity_X (p4est_geometry_t * geom,
                               p4est_topidx_t which_tree,
                               const double abc[3], double xyz[3])
{
  P4EST_ASSERT (geom->user != NULL);
  p4est_connectivity_t *connectivity = (p4est_connectivity_t *) geom->user;
  P4EST_ASSERT (connectivity->tree_to_vertex != NULL);
  const p4est_topidx_t *tree_to_vertex = connectivity->tree_to_vertex;
  const double       *v = connectivity->vertices;
  double              eta_x, eta_y, eta_z = 0.;
  int                 j, k;
  p4est_topidx_t      vt[P4EST_CHILDREN];

  /* retrieve corners of the tree */
  for (k = 0; k < P4EST_CHILDREN; ++k) {
    vt[k] = tree_to_vertex[which_tree * P4EST_CHILDREN + k];
  }

  /* these are reference coordinates in [0, 1]**d */
  eta_x = abc[0];
  eta_y = abc[1];
#ifdef P4_TO_P8
  eta_z = abc[2];
#endif

  /* bi/trilinear transformation */
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
  }
}

p4est_geometry_t   *
p4est_geometry_new_connectivity (p4est_connectivity_t * conn)
{
  p4est_geometry_t   *geom;

  P4EST_ASSERT (conn->vertices != NULL);

  geom = P4EST_ALLOC_ZERO (p4est_geometry_t, 1);

  geom->name = P4EST_STRING "_connectivity";
  geom->user = conn;
  geom->X = p4est_geometry_connectivity_X;

  return geom;
}

#ifndef P4_TO_P8

/**
 * Geometric coordinate transformation for icosahedron geometry.
 *
 * Define the geometric transformation from tree-local reference coordinates
 * to physical space
 *
 * \param[in]  geom       associated geometry
 * \param[in]  which_tree tree id inside forest
 * \param[in]  rst        tree-local reference coordinates : [0,1]^2.
 *                        Note: rst[2] is never accessed
 * \param[out] xyz        Cartesian coordinates in physical space after geometry
 *
 */
static void
p4est_geometry_icosahedron_X (p4est_geometry_t * geom,
                              p4est_topidx_t which_tree,
                              const double rst[3], double xyz[3])
{
  const struct p4est_geometry_builtin_icosahedron *icosahedron
    = &((p4est_geometry_builtin_t *) geom)->p.icosahedron;
  double              a = 1.0;
  double              g = (1.0 + sqrt (5.0)) * 0.5;     /* golden ratio */
  double              ga = a / g;
  double              r = sqrt ((5.0 - sqrt (5.0)) * 0.5);      /* sqrt(a*a+ga*ga) -> current radius */

  double              R = icosahedron->R;       /* target sphere radius */
  double              radius_ratio = R / r;

  /* these are reference coordinates in [0, 1]**d */
  double              eta_x, eta_y;
  eta_x = rst[0];
  eta_y = rst[1];

  /*
   * icosahedron node Cartesian coordinates
   * used for mapping connectivity vertices to 3D nodes.
   */
  const double        N[12 * 3] = {
    0, -ga, a,                  /*  N0 */
    ga, -a, 0,                  /*  N1 */
    a, 0, ga,                   /*  N2 */
    0, ga, a,                   /*  N3 */
    -a, 0, ga,                  /*  N4 */
    -ga, -a, 0,                 /*  N5 */
    a, 0, -ga,                  /*  N6 */
    ga, a, 0,                   /*  N7 */
    -ga, a, 0,                  /*  N8 */
    -a, 0, -ga,                 /*  N9 */
    0, -ga, -a,                 /* N10 */
    0, ga, -a,                  /* N11 */
  };

  /*
   * tree to nodes:
   *
   * tree 0:  1  6  0  2
   * tree 1:  2  7  0  3
   * tree 2:  3  8  0  4
   * tree 3:  4  9  0  5
   * tree 4:  5 10  0  1
   * tree 5:  6 11  2  7
   * tree 6:  7 11  3  8
   * tree 7:  8 11  4  9
   * tree 8:  9 11  5 10
   * tree 9: 10 11  1  6
   *
   */
  const int           tree_to_nodes[10 * 4] = {
    1, 6, 0, 2,
    6, 11, 2, 7,
    2, 7, 0, 3,
    7, 11, 3, 8,
    3, 8, 0, 4,
    8, 11, 4, 9,
    4, 9, 0, 5,
    9, 11, 5, 10,
    5, 10, 0, 1,
    10, 11, 1, 6,
  };

  /* assert that input points are in the expected range */
  P4EST_ASSERT (icosahedron->type == P4EST_GEOMETRY_BUILTIN_ICOSAHEDRON);
  P4EST_ASSERT (0 <= which_tree && which_tree < 10);

  /* transform from the reference cube into vertex space */

  /* assign correct coordinates based on patch id */
  /* use bilinear SLERP :  spherical bilinear interpolation */
  {
    int                 j;

    /* use tree to nodes mapping to get nodes index of current tree */
    const int           i0 = tree_to_nodes[which_tree * 4 + 0];
    const int           i1 = tree_to_nodes[which_tree * 4 + 1];
    const int           i2 = tree_to_nodes[which_tree * 4 + 2];
    const int           i3 = tree_to_nodes[which_tree * 4 + 3];

    /* get 3D Cartesian coordinates of our face */
    const double        n0[3] =
      { N[i0 * 3 + 0], N[i0 * 3 + 1], N[i0 * 3 + 2] };
    const double        n1[3] =
      { N[i1 * 3 + 0], N[i1 * 3 + 1], N[i1 * 3 + 2] };
    const double        n2[3] =
      { N[i2 * 3 + 0], N[i2 * 3 + 1], N[i2 * 3 + 2] };
    const double        n3[3] =
      { N[i3 * 3 + 0], N[i3 * 3 + 1], N[i3 * 3 + 2] };
    double              norme2 =
      n0[0] * n0[0] + n0[1] * n0[1] + n0[2] * n0[2];

    /* 1. apply slerp
     * - between n0 and n1
     * - between n2 and n3
     */
    double              xyz01[3];       /* slerp along n0 -> n1 */
    double              xyz23[3];       /* slerp along n2 -> n3 */
    double              dot1 = n0[0] * n1[0] + n0[1] * n1[1] + n0[2] * n1[2];
    double              theta1 = acos (dot1 / norme2);

    for (j = 0; j < 3; ++j) {
      xyz01[j] =
        sin ((1.0 - eta_x) * theta1) / sin (theta1) * n0[j] +
        sin ((eta_x) * theta1) / sin (theta1) * n1[j];
      xyz23[j] =
        sin ((1.0 - eta_x) * theta1) / sin (theta1) * n2[j] +
        sin ((eta_x) * theta1) / sin (theta1) * n3[j];
    }

    /* apply slerp between xyz01 and xyz23 */
    double              dot2 =
      xyz01[0] * xyz23[0] + xyz01[1] * xyz23[1] + xyz01[2] * xyz23[2];
    double              theta2 = acos (dot2 / norme2);
    for (j = 0; j < 3; ++j) {
      xyz[j] =
        sin ((1.0 - eta_y) * theta2) / sin (theta2) * xyz01[j] +
        sin ((eta_y) * theta2) / sin (theta2) * xyz23[j];
      xyz[j] *= radius_ratio;   /* rescale coordinates to target radius */
    }

    /* printf("DEBUG : %g | %g %g %g | \n", */
    /*     xyz[0]*xyz[0] + xyz[1]*xyz[1] + xyz[2]*xyz[2], */
    /*     xyz[0],xyz[1],xyz[2] */
    /*     ); */

  }                             /* end of bilinear slerp */

}                               /* p4est_geometry_icosahedron_X */

p4est_geometry_t   *
p4est_geometry_new_icosahedron (p4est_connectivity_t * conn, double R)
{
  p4est_geometry_builtin_t *builtin;
  struct p4est_geometry_builtin_icosahedron *icosahedron;

  builtin = P4EST_ALLOC_ZERO (p4est_geometry_builtin_t, 1);

  icosahedron = &builtin->p.icosahedron;
  icosahedron->type = P4EST_GEOMETRY_BUILTIN_ICOSAHEDRON;
  icosahedron->R = R;

  builtin->geom.name = "p4est_icosahedron";
  builtin->geom.user = conn;
  builtin->geom.X = p4est_geometry_icosahedron_X;

  return (p4est_geometry_t *) builtin;

}                               /* p4est_geometry_new_icosahedron */

/**
 * Geometric coordinate transformation for shell2d geometry.
 *
 * Define the geometric transformation from tree-local reference coordinates
 * to physical space.
 *
 * \param[in]  geom       associated geometry
 * \param[in]  which_tree tree id inside forest
 * \param[in]  rst        tree-local reference coordinates : [0,1]^2.
 *                        Note: rst[2] is never accessed
 * \param[out] xyz        Cartesian coordinates in physical space after geometry
 *
 */
static void
p4est_geometry_shell2d_X (p4est_geometry_t * geom,
                          p4est_topidx_t which_tree,
                          const double rst[3], double xyz[3])
{
  const struct p4est_geometry_builtin_shell2d *shell2d
    = &((p4est_geometry_builtin_t *) geom)->p.shell2d;
  double              x, R, q;
  double              abc[3];

  xyz[2] = 0.0;

  /* transform from the reference cube into vertex space */
  p4est_geometry_connectivity_X (geom, which_tree, rst, abc);

  /* assert that input points are in the expected range */
  P4EST_ASSERT (shell2d->type == P4EST_GEOMETRY_BUILTIN_SHELL2D);
  P4EST_ASSERT (0 <= which_tree && which_tree < 8);
  P4EST_ASSERT (abc[0] < 1.0 + SC_1000_EPS && abc[0] > -1.0 - SC_1000_EPS);
  P4EST_ASSERT (abc[1] < 2.0 + SC_1000_EPS && abc[1] > 1.0 - SC_1000_EPS);

  /* abc[2] is always 0 here ... */

  /* transform abc[0] in-place for nicer grading */
  x = tan (abc[0] * M_PI_4);

  /* compute transformation ingredients */
  R = shell2d->R1sqrbyR2 * pow (shell2d->R2byR1, abc[1]);
  q = R / sqrt (x * x + 1.);

  /* assign correct coordinates based on patch id */
  switch (which_tree / 2) {
  case 0:                      /* bottom */
    xyz[0] = +q;                /*  R*cos(theta) */
    xyz[1] = +q * x;            /*  R*sin(theta) */
    break;
  case 1:                      /* right */
    xyz[0] = -q * x;            /* -R*sin(theta) = R*cos(theta+PI/2) */
    xyz[1] = +q;                /*  R*cos(theta) = R*sin(theta+PI/2) */
    break;
  case 2:                      /* top */
    xyz[0] = -q;                /* - R*cos(theta) = R*cos(theta+PI) */
    xyz[1] = -q * x;            /* - R*sin(theta) = R*sin(theta+PI) */
    break;
  case 3:                      /* left */
    xyz[0] = +q * x;            /*  R*sin(theta) = R*cos(theta+3*PI/2) */
    xyz[1] = -q;                /* -R*cos(theta) = R*sin(theta+3*PI/2) */
    break;
  default:
    SC_ABORT_NOT_REACHED ();
  }
}                               /* p4est_geometry_shell2d_X */

p4est_geometry_t   *
p4est_geometry_new_shell2d (p4est_connectivity_t * conn, double R2, double R1)
{
  p4est_geometry_builtin_t *builtin;
  struct p4est_geometry_builtin_shell2d *shell2d;

  builtin = P4EST_ALLOC_ZERO (p4est_geometry_builtin_t, 1);

  shell2d = &builtin->p.shell2d;
  shell2d->type = P4EST_GEOMETRY_BUILTIN_SHELL2D;
  shell2d->R2 = R2;
  shell2d->R1 = R1;
  shell2d->R2byR1 = R2 / R1;
  shell2d->R1sqrbyR2 = R1 * R1 / R2;
  shell2d->Rlog = log (R2 / R1);

  builtin->geom.name = "p4est_shell2d";
  builtin->geom.user = conn;
  builtin->geom.X = p4est_geometry_shell2d_X;

  return (p4est_geometry_t *) builtin;

}                               /* p4est_geometry_new_shell2d */

/**
 * geometric coordinate transformation for disk2d geometry.
 *
 * Define the geometric transformation from tree-local reference coordinates
 * to physical space.
 *
 * \param[in]  geom       associated geometry
 * \param[in]  which_tree tree id inside forest
 * \param[in]  rst        tree-local reference coordinates : [0,1]^2.
 *                        Note: rst[2] is never accessed.
 * \param[out] xyz        Cartesian coordinates in physical space after geometry
 *
 * Note abc[3] contains Cartesian coordinates in logical
 * vertex space (before geometry).
 */
static void
p4est_geometry_disk2d_X (p4est_geometry_t * geom,
                         p4est_topidx_t which_tree,
                         const double rst[3], double xyz[3])
{
  const p4est_geometry_builtin_disk2d_t *disk2d
    = &((p4est_geometry_builtin_t *) geom)->p.disk2d;
  double              x, y, R, q;
  double              abc[3];

  (void) y;

  /* in 2D z is ZERO ! */
  xyz[2] = 0.0;

  /* transform from the reference cube [0,1]^3 into logical vertex space
     using bi/trilinear transformation */
  p4est_geometry_connectivity_X (geom, which_tree, rst, abc);

  /*
   * assert that input points are in the expected range
   * Note: maybe we should remove these assert, this would allow
   * ghost quadrant at external boundary to call this routine ?
   */
  P4EST_ASSERT (disk2d->type == P4EST_GEOMETRY_BUILTIN_DISK2D);
  P4EST_ASSERT (0 <= which_tree && which_tree < 5);
  P4EST_ASSERT (abc[0] < 1.0 + SC_1000_EPS && abc[0] > -1.0 - SC_1000_EPS);
  if (which_tree < 4)
    P4EST_ASSERT (abc[1] < 2.0 + SC_1000_EPS && abc[1] > 1.0 - SC_1000_EPS);
  else
    P4EST_ASSERT (abc[1] < 1.0 + SC_1000_EPS && abc[1] > -1.0 - SC_1000_EPS);

  /* abc[2] is always 0 here and so unused in 2D ... */

  if (which_tree < 4) {
    double              p, tanx;

    p = 2.0 - abc[1];
    tanx = -tan (abc[0] * M_PI_4);      /* x = tan (theta) */

    x = p * (-abc[0]) + (1. - p) * tanx;

    /* compute transformation ingredients */
    R = disk2d->R0sqrbyR1 * pow (disk2d->R1byR0, abc[1]);

    /* R*cos(theta) */
    /* q = R / sqrt (x * x + 1.); */
    q = R / sqrt (1. + (1. - p) * (tanx * tanx) + 1. * p);

    /* assign correct coordinates based on patch id */
    switch (which_tree) {
    case 0:                    /* bottom */
      xyz[0] = +q;              /*   R*cos(theta) */
      xyz[1] = +q * x;          /*   R*sin(theta) */
      break;
    case 1:                    /* right */
      xyz[0] = +q * x;          /*   R*sin(theta) = R*cos(theta-PI/2) */
      xyz[1] = -q;              /* - R*cos(theta) = R*sin(theta-PI/2) */
      break;
    case 2:                    /* top */
      xyz[0] = -q;              /* - R*cos(theta) = R*cos(theta-PI) */
      xyz[1] = -q * x;          /* - R*sin(theta) = R*sin(theta-PI) */
      break;
    case 3:                    /* left */
      xyz[0] = -q * x;          /* -R*sin(theta) = R*cos(theta-3*PI/2) */
      xyz[1] = +q;              /*  R*cos(theta) = R*sin(theta-3*PI/2) */
      break;
    default:
      SC_ABORT_NOT_REACHED ();
    }

  }
  else {

    /* center square */
    xyz[0] = abc[0] * disk2d->Clength;
    xyz[1] = abc[1] * disk2d->Clength;
    xyz[2] = 0.0;

  }

}                               /* p4est_geometry_disk2d_X */

p4est_geometry_t   *
p4est_geometry_new_disk2d (p4est_connectivity_t * conn, double R0, double R1)
{
  p4est_geometry_builtin_t *builtin;
  struct p4est_geometry_builtin_disk2d *disk2d;

  builtin = P4EST_ALLOC_ZERO (p4est_geometry_builtin_t, 1);

  disk2d = &builtin->p.disk2d;
  disk2d->type = P4EST_GEOMETRY_BUILTIN_DISK2D;
  disk2d->R0 = R0;
  disk2d->R1 = R1;

  /* variables useful for the outer shell */
  disk2d->R1byR0 = R1 / R0;
  disk2d->R0sqrbyR1 = R0 * R0 / R1;
  disk2d->R0log = log (R1 / R0);

  /* variables useful for the center square */
  disk2d->Clength = R0 / sqrt (2.);
  /* disk2d->CdetJ = pow (R0 / sqrt (3.), 3.); */

  builtin->geom.name = "p4est_disk2d";
  builtin->geom.user = conn;
  builtin->geom.X = p4est_geometry_disk2d_X;

  return (p4est_geometry_t *) builtin;

}                               /* p4est_geometry_new_disk2d */

/**
 * geometric coordinate transformation for sphere2d geometry.
 *
 * Define the geometric transformation from tree-local reference coordinates to the
 * physical space.
 *
 * \param[in]  geom       associated geometry
 * \param[in]  which_tree tree id inside forest
 * \param[in]  rst        tree-local reference coordinates : [0,1]^2.
 *                        Note: rst[2] is never accessed
 * \param[out] xyz        Cartesian coordinates in physical space after geometry
 *
 */
static void
p4est_geometry_sphere2d_X (p4est_geometry_t * geom,
                           p4est_topidx_t which_tree,
                           const double rst[3], double xyz[3])
{
  const struct p4est_geometry_builtin_sphere2d *sphere2d
    = &((p4est_geometry_builtin_t *) geom)->p.sphere2d;
  double              R;

  /* transform from the tree-local reference coordinates into the cube-surface
   * in physical space using vertex bi/trilinear transformation.
   */
  p4est_geometry_connectivity_X (geom, which_tree, rst, xyz);

  /* align cube center with origin */
  xyz[0] -= 0.5;
  xyz[1] -= 0.5;
  xyz[2] -= 0.5;

  /* normalise to radius R sphere */
  R = sphere2d->R;
  double              R_on_norm =
    R / sqrt (xyz[0] * xyz[0] + xyz[1] * xyz[1] + xyz[2] * xyz[2]);
  xyz[0] *= R_on_norm;
  xyz[1] *= R_on_norm;
  xyz[2] *= R_on_norm;
}                               /* p4est_geometry_sphere2d_X */

p4est_geometry_t   *
p4est_geometry_new_sphere2d (p4est_connectivity_t * conn, double R)
{
  p4est_geometry_builtin_t *builtin;
  struct p4est_geometry_builtin_sphere2d *sphere2d;

  builtin = P4EST_ALLOC_ZERO (p4est_geometry_builtin_t, 1);

  sphere2d = &builtin->p.sphere2d;
  sphere2d->type = P4EST_GEOMETRY_BUILTIN_SPHERE2D;
  sphere2d->R = R;

  builtin->geom.name = "p4est_sphere2d";
  builtin->geom.user = conn;
  builtin->geom.X = p4est_geometry_sphere2d_X;

  return (p4est_geometry_t *) builtin;
}                               /* p4est_geometry_new_sphere2d */

#endif /* !P4_TO_P8 */

static void
p4est_geometry_Q2_node (const p4est_quadrant_t *quad,
#ifdef P4_TO_P8
                        int k,
#endif
                        int j, int i, p4est_quadrant_t *cnode)
{
  const p4est_qcoord_t h2 = P4EST_QUADRANT_LEN (quad->level + 1);

  P4EST_ASSERT (p4est_quadrant_is_valid (quad));
  P4EST_ASSERT (0 <= i && i <= 2);
  P4EST_ASSERT (0 <= j && j <= 2);
#ifdef P4_TO_P8
  P4EST_ASSERT (0 <= k && k <= 2);
#endif

  cnode->x = quad->x + i * h2;
  cnode->y = quad->y + j * h2;
#ifdef P4_TO_P8
  cnode->z = quad->z + k * h2;
#endif
  cnode->level = P4EST_MAXLEVEL;
}

typedef struct p4est_geometry_bounds
{
  p4est_topidx_t      end_which_tree;
  p4est_locidx_t      end_local_node;
}
p4est_geometry_bounds_t;

/** A geometry coordinate tuple with tree and node information. */
typedef struct p4est_geometry_node_coordinate
{
  p4est_locidx_t      coord_index;      /**< Index of node coordinate. */
  p4est_locidx_t      local_node;       /**< Local node index of point. */
  p4est_topidx_t      which_tree;       /**< Tree number for this point. */

  /** Tree boundary index in [0, \ref P4EST_INSUL) of a tree node.
   * To ensure correct node coordinates even in the case of periodic meshes,
   * we hash the tree number, the local node, and this tree boundary index.
   */
  int8_t              tree_bound;
}
p4est_geometry_node_coordinate_t;

static unsigned int
p4est_geometry_node_hash (const void *v, const void *u)
{
  uint32_t            utt, uln, c;
  const p4est_geometry_node_coordinate_t *nt =
    (const p4est_geometry_node_coordinate_t *) v;
#ifdef P4EST_ENABLE_DEBUG
  const p4est_geometry_bounds_t *tb = (const p4est_geometry_bounds_t *) u;

  P4EST_ASSERT (v != NULL);
  if (tb != NULL) {
    /* may be NULL due to use outside of constructor */
    P4EST_ASSERT (0 <= nt->which_tree && nt->which_tree < tb->end_which_tree);
    P4EST_ASSERT (0 <= nt->local_node && nt->local_node < tb->end_local_node);
  }
  P4EST_ASSERT (0 <= nt->tree_bound && nt->tree_bound < P4EST_INSUL);
#endif

  /* execute primitive hash function */
  utt = (uint32_t) nt->which_tree;
  uln = (uint32_t) nt->local_node;
  c = (uint32_t) nt->tree_bound;
  sc_hash_mix (utt, uln, c);
  sc_hash_final (utt, uln, c);

  return (unsigned int) c;
}

static int
p4est_geometry_node_equal (const void *v1, const void *v2, const void *u)
{
  const p4est_geometry_node_coordinate_t *nt1 =
    (const p4est_geometry_node_coordinate_t *) v1;
  const p4est_geometry_node_coordinate_t *nt2 =
    (const p4est_geometry_node_coordinate_t *) v2;
#ifdef P4EST_ENABLE_DEBUG
  const p4est_geometry_bounds_t *tb = (const p4est_geometry_bounds_t *) u;

  P4EST_ASSERT (v1 != NULL);
  P4EST_ASSERT (v2 != NULL);
  if (tb != NULL) {
    /* may be NULL due to use outside of constructor */
    P4EST_ASSERT (0 <= nt1->which_tree
                  && nt1->which_tree < tb->end_which_tree);
    P4EST_ASSERT (0 <= nt1->local_node
                  && nt1->local_node < tb->end_local_node);
    P4EST_ASSERT (0 <= nt2->which_tree
                  && nt2->which_tree < tb->end_which_tree);
    P4EST_ASSERT (0 <= nt2->local_node
                  && nt2->local_node < tb->end_local_node);
  }
  P4EST_ASSERT (0 <= nt1->tree_bound && nt1->tree_bound < P4EST_INSUL);
  P4EST_ASSERT (0 <= nt2->tree_bound && nt2->tree_bound < P4EST_INSUL);
#endif

  return nt1->which_tree == nt2->which_tree &&
    nt1->local_node == nt2->local_node && nt1->tree_bound == nt2->tree_bound;
}

void
p4est_geometry_coordinates_new_lnodes
  (p4est_t *p4est, p4est_geometry_t *geom,
   p4est_lnodes_t *lnodes, const double *refloc,
   sc_array_t *coordinates, sc_array_t *element_coordinates)
{
  static const double irlen = 1. / P4EST_ROOT_LEN;
  int                 vno, vd, deg;
  int                 i, kjixo, j, kjxo, kji, kxo;
  int                 cid, fcd;
  int                 n;
#ifdef P4_TO_P8
  int                 k;
  int                 l;
  int                 e;
#endif
  int                 dtb[P4EST_DIM], dth[P4EST_DIM], dts;
  double              abc[3], *xyz;
  p4est_topidx_t      tt;
  p4est_locidx_t      el, ne;
  p4est_locidx_t      collected, duplicates;
  p4est_locidx_t     *enodes, *ecoords;
  p4est_lnodes_code_t *fcodes, fc;
  p4est_quadrant_t    sparent, *parent = &sparent, *quad;
  p4est_quadrant_t    scnode, *cnode = &scnode, *thequad;
  p4est_tree_t       *tree;
  sc_mempool_t       *pool;
  sc_hash_t          *hash;
  void              **found;
#ifndef P4EST_ENABLE_DEBUG
  void               *hash_user = NULL;
#else
  p4est_geometry_bounds_t stb, *tb = &stb;
  void               *hash_user = (void *) tb;
#endif
  p4est_geometry_node_coordinate_t stn, *tn = &stn, *inserted;

  /* basic checks */
  P4EST_ASSERT (p4est != NULL);
  P4EST_ASSERT (p4est->connectivity != NULL);
  P4EST_ASSERT (geom != NULL || p4est->connectivity->vertices != NULL);
  P4EST_ASSERT (geom == NULL || geom->X != NULL);
  P4EST_ASSERT (lnodes != NULL);
  P4EST_ASSERT (lnodes->degree == 1 || lnodes->degree == 2);
  P4EST_ASSERT ((lnodes->degree == 1 && lnodes->vnodes == P4EST_CHILDREN) ||
                (lnodes->degree == 2 && lnodes->vnodes == P4EST_INSUL));
  P4EST_ASSERT (p4est->local_num_quadrants == lnodes->num_local_elements);
  P4EST_ASSERT (coordinates != NULL &&
                coordinates->elem_size == 3 * sizeof (double));
  P4EST_ASSERT (element_coordinates != NULL &&
                element_coordinates->elem_size == sizeof (p4est_locidx_t));

  /* prepare lookup information */
  pool = sc_mempool_new (sizeof (p4est_geometry_node_coordinate_t));
#ifdef P4EST_ENABLE_DEBUG
  memset (tb, 0, sizeof (stb));
  tb->end_which_tree = p4est->connectivity->num_trees;
  tb->end_local_node = lnodes->num_local_nodes;
#endif
  hash = sc_hash_new
    (p4est_geometry_node_hash, p4est_geometry_node_equal, hash_user, NULL);

  /* loop through p4est qaudrants in natural order */
  abc[2] = 0.;
  memset (tn, 0, sizeof (stn));
  vno = lnodes->vnodes;
  vd = (deg = lnodes->degree) + 1;
  fcodes = lnodes->face_code;
  enodes = lnodes->element_nodes;
  sc_array_resize (coordinates, 0);
  sc_array_resize (element_coordinates, lnodes->num_local_elements * vno);
  ecoords = (p4est_locidx_t *) sc_array_index (element_coordinates, 0);
  collected = duplicates = 0;
  for (tt = p4est->first_local_tree; tt <= p4est->last_local_tree; ++tt) {
    tn->which_tree = tt;
    tree = p4est_tree_array_index (p4est->trees, tt);
    ne = (p4est_locidx_t) tree->quadrants.elem_count;
    for (el = 0; el < ne; ++el, enodes += vno, ecoords += vno) {

      /* access quadrant in forest */
      quad = p4est_quadrant_array_index (&tree->quadrants, (size_t) el);
      fc = *(fcodes++);
      if (fc) {
        /* this element is hanging */
        cid = fc & (P4EST_CHILDREN - 1);
        P4EST_ASSERT (cid == p4est_quadrant_child_id (quad));
        fc >>= P4EST_DIM;
        P4EST_ASSERT (fc);
        p4est_quadrant_parent (quad, parent);
      }
      else {
        /* non-hanging element */
        cid = 0;
        P4EST_QUADRANT_INIT (parent);
      }

      /* iterate through the local nodes referenced by the element */
      kji = 0;
#ifndef P4_TO_P8
      kxo = 0;
#else
      for (k = 0; k < vd; ++k) {
        kxo = (((cid & 4) ? deg - k : k) != 0) << 2;
#if 0
      }
#endif
#endif
      for (j = 0; j < vd; ++j) {
        kjxo = kxo + ((((cid & 2) ? deg - j : j) != 0) << 1);
        for (i = 0; i < vd; ++i, ++kji) {
          kjixo = kjxo + (((cid & 1) ? deg - i : i) != 0);
          P4EST_ASSERT (0 <= kjixo && kjixo < P4EST_CHILDREN);
          P4EST_ASSERT (deg != 1 || (kjixo ^ cid) == kji);
          tn->local_node = enodes[kji];

          /* compute relevant quadrant to determine node coordinates */
          thequad = quad;
          if (fc) {
            fcd = p4est_lnodes_corner_hanging[kjixo];
            if (fcd >= 0 && (fc & (1 << fcd))) {
              thequad = parent;
            }
          }

          /* compute coordinates of quadrant node */
          if (deg == 1) {
            p4est_quadrant_corner_node (thequad, kji, cnode);
          }
          else {
            P4EST_ASSERT (deg == 2);
            p4est_geometry_Q2_node (thequad,
#ifdef P4_TO_P8
                                    k,
#endif
                                    j, i, cnode);
          }

          /* compute tree boundary status of quadrant node */
          dts =
#ifdef P4_TO_P8
            (dtb[2] =
             (dth[2] = (cnode->z == P4EST_ROOT_LEN)) || cnode->z == 0) +
#endif
            (dtb[1] =
             (dth[1] = (cnode->y == P4EST_ROOT_LEN)) || cnode->y == 0) +
            (dtb[0] =
             (dth[0] = (cnode->x == P4EST_ROOT_LEN)) || cnode->x == 0);
          switch (dts) {
          case 0:
            /* the most frequent case comes first */
            tn->tree_bound = P4EST_INSUL / 2;
            break;
          case 1:
            /* the node sits inside a tree face */
            for (n = 0; n < P4EST_DIM; ++n) {
              if (dtb[n]) {
                tn->tree_bound = p4est_face_points[(n << 1) + dth[n]];
                break;
              }
            }
            P4EST_ASSERT (n < P4EST_DIM);
            break;
#ifdef P4_TO_P8
          case 2:
            /* the node sits inside a tree edge */
            e = 0;
            l = 0;
            for (n = 0; n < P4EST_DIM; ++n) {
              if (!dtb[n]) {
                e += n << 2;
              }
              else {
                e += dth[n] ? (1 << l) : 0;
                ++l;
              }
            }
            P4EST_ASSERT (l == 2);
            tn->tree_bound = p8est_edge_points[e];
            break;
#endif
          case P4EST_DIM:
            /* the node sits on a tree corner equal a quadrant corner */
            tn->tree_bound = (deg == 1) ? p4est_corner_points[kji] : kji;
            break;
          default:
            SC_ABORT_NOT_REACHED ();
          }
          P4EST_ASSERT (0 <= tn->tree_bound && tn->tree_bound < P4EST_INSUL);

          /* determine whether this local node has already been processed */
          if (!sc_hash_insert_unique (hash, tn, &found)) {

            /* this tree node is already computed */
            inserted = *(p4est_geometry_node_coordinate_t **) found;
            P4EST_ASSERT (p4est_geometry_node_equal (inserted, tn, hash_user));
            P4EST_ASSERT (0 <= inserted->coord_index && (size_t)
                          inserted->coord_index < coordinates->elem_count);
            ecoords[kji] = inserted->coord_index;

            /* no more to be done for this element node */
            ++duplicates;
#if 0
            P4EST_LDEBUGF
              ("Duplicate tree %ld element %ld boundary %d cid %d index %d node %ld coord %ld\n",
               (long) tn->which_tree, (long) tn->local_node,
               (int) tn->tree_bound, cid, kji, (long) enodes[kji], (long) ecoords[kji]);
#endif
            continue;
          }

          /* install a new element node coordinate in output array */
          inserted = *(p4est_geometry_node_coordinate_t **) found =
            (p4est_geometry_node_coordinate_t *) sc_mempool_alloc (pool);
          ecoords[kji] = inserted->coord_index =
            (p4est_locidx_t) coordinates->elem_count;

          /* remember node coordinate key */
          inserted->local_node = tn->local_node;
          inserted->which_tree = tn->which_tree;
          inserted->tree_bound = tn->tree_bound;

          /* we write to this element node coordinate below */
          xyz = (double *) sc_array_push (coordinates);
          ++collected;
#if 0
          P4EST_LDEBUGF
            ("Added tree %ld element %ld boundary %d cid %d index %d node %ld coord %ld\n",
             (long) tn->which_tree, (long) tn->local_node,
             (int) tn->tree_bound, cid, kji, (long) enodes[kji], (long) ecoords[kji]);
#endif

          /* apply geometry transformation */
          if (geom == NULL) {
            p4est_qcoord_to_vertex (p4est->connectivity, tt,
                                    cnode->x, cnode->y,
#ifdef P4_TO_P8
                                    cnode->z,
#endif
                                    xyz);
          }
          else {
            abc[0] = cnode->x * irlen;
            abc[1] = cnode->y * irlen;
#ifdef P4_TO_p8
            abc[2] = cnode->z * irlen;
#endif
            geom->X (geom, tt, abc, xyz);
          }

        }                       /* i loop */
      }                         /* j loop */
#ifdef P4_TO_P8
#if 0
      {
#endif
      }                         /* k loop */
#endif
    }                           /* element loop */
  }                             /* tree loop */
  P4EST_ASSERT (fcodes - lnodes->face_code ==
                (ptrdiff_t) lnodes->num_local_elements);
  P4EST_ASSERT (enodes - lnodes->element_nodes ==
                (ptrdiff_t) lnodes->num_local_elements * vno);
  P4EST_ASSERT (collected + duplicates ==
                lnodes->num_local_elements * vno);
  P4EST_ASSERT (collected >= lnodes->num_local_nodes);

  sc_hash_destroy (hash);
  sc_mempool_destroy (pool);
}
