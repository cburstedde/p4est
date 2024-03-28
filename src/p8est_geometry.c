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
 * \file p8est_geometry.c
 * We provide some geometry transformations for the builtin connectivities.
 * They are not meant as blueprints for future user code.
 * Please implement p8est_geometry_t as you see fit.
 */

#include <p4est_to_p8est.h>
#include "p4est_geometry.c"

typedef enum
{
  P8EST_GEOMETRY_BUILTIN_MAGIC = 0x30F3F8DF,
  P8EST_GEOMETRY_BUILTIN_SHELL,
  P8EST_GEOMETRY_BUILTIN_SPHERE,
  P8EST_GEOMETRY_BUILTIN_TORUS,
  P8EST_GEOMETRY_LAST
}
p8est_geometry_builtin_type_t;

typedef struct p8est_geometry_builtin_shell
{
  p8est_geometry_builtin_type_t type;
  double              R2, R1;
  double              R2byR1, R1sqrbyR2, Rlog;
}
p8est_geometry_builtin_shell_t;

typedef struct p8est_geometry_builtin_sphere
{
  p8est_geometry_builtin_type_t type;
  double              R2, R1, R0;
  double              R2byR1, R1sqrbyR2, R1log;
  double              R1byR0, R0sqrbyR1, R0log;
  double              Clength, CdetJ;
}
p8est_geometry_builtin_sphere_t;

typedef struct p8est_geometry_builtin_torus
{
  p8est_geometry_builtin_type_t type;
  int                 nSegments;
  double              R2, R1, R0;
  double              R1byR0, R0sqrbyR1, R0log;
  double              Clength;
}
p8est_geometry_builtin_torus_t;

typedef struct p8est_geometry_builtin
{
  /** The geom member needs to come first; we cast to p8est_geometry_t * */
  p8est_geometry_t    geom;
  union
  {
    p8est_geometry_builtin_type_t type;
    p8est_geometry_builtin_shell_t shell;
    p8est_geometry_builtin_sphere_t sphere;
    p8est_geometry_builtin_torus_t torus;
  }
  p;
}
p8est_geometry_builtin_t;

static void
p8est_geometry_shell_X (p8est_geometry_t * geom,
                        p4est_topidx_t which_tree,
                        const double rst[3], double xyz[3])
{
  const struct p8est_geometry_builtin_shell *shell
    = &((p8est_geometry_builtin_t *) geom)->p.shell;
  double              x, y, R, q;
  double              abc[3];

  /* transform from the reference cube into vertex space */
  p8est_geometry_connectivity_X (geom, which_tree, rst, abc);

  /* assert that input points are in the expected range */
  P4EST_ASSERT (shell->type == P8EST_GEOMETRY_BUILTIN_SHELL);
  P4EST_ASSERT (0 <= which_tree && which_tree < 24);
  P4EST_ASSERT (abc[0] < 1.0 + SC_1000_EPS && abc[0] > -1.0 - SC_1000_EPS);
  P4EST_ASSERT (abc[1] < 1.0 + SC_1000_EPS && abc[1] > -1.0 - SC_1000_EPS);
  P4EST_ASSERT (abc[2] < 2.0 + SC_1000_EPS && abc[2] > 1.0 - SC_1000_EPS);

  /* transform abc[0] and y in-place for nicer grading */
  x = tan (abc[0] * M_PI_4);
  y = tan (abc[1] * M_PI_4);

  /* compute transformation ingredients */
  R = shell->R1sqrbyR2 * pow (shell->R2byR1, abc[2]);
  q = R / sqrt (x * x + y * y + 1.);

  /* assign correct coordinates based on patch id */
  switch (which_tree / 4) {
  case 3:                      /* top */
    xyz[0] = +q * y;
    xyz[1] = -q * x;
    xyz[2] = +q;
    break;
  case 2:                      /* left */
    xyz[0] = -q;
    xyz[1] = -q * x;
    xyz[2] = +q * y;
    break;
  case 1:                      /* bottom */
    xyz[0] = -q * y;
    xyz[1] = -q * x;
    xyz[2] = -q;
    break;
  case 0:                      /* right */
    xyz[0] = +q;
    xyz[1] = -q * x;
    xyz[2] = -q * y;
    break;
  case 4:                      /* back */
    xyz[0] = -q * x;
    xyz[1] = +q;
    xyz[2] = +q * y;
    break;
  case 5:                      /* front */
    xyz[0] = +q * x;
    xyz[1] = -q;
    xyz[2] = +q * y;
    break;
  default:
    SC_ABORT_NOT_REACHED ();
  }
}

p8est_geometry_t   *
p8est_geometry_new_shell (p8est_connectivity_t * conn, double R2, double R1)
{
  p8est_geometry_builtin_t *builtin;
  struct p8est_geometry_builtin_shell *shell;

  builtin = P4EST_ALLOC_ZERO (p8est_geometry_builtin_t, 1);

  shell = &builtin->p.shell;
  shell->type = P8EST_GEOMETRY_BUILTIN_SHELL;
  shell->R2 = R2;
  shell->R1 = R1;
  shell->R2byR1 = R2 / R1;
  shell->R1sqrbyR2 = R1 * R1 / R2;
  shell->Rlog = log (R2 / R1);

  builtin->geom.name = "p8est_shell";
  builtin->geom.user = conn;
  builtin->geom.X = p8est_geometry_shell_X;

  return (p8est_geometry_t *) builtin;
}

static void
p8est_geometry_sphere_X (p8est_geometry_t * geom,
                         p4est_topidx_t which_tree,
                         const double rst[3], double xyz[3])
{
  const struct p8est_geometry_builtin_sphere *sphere
    = &((p8est_geometry_builtin_t *) geom)->p.sphere;
  double              x, y, R, q;
  double              abc[3];

  /* transform from the reference cube into vertex space */
  p8est_geometry_connectivity_X (geom, which_tree, rst, abc);

  /* assert that input points are in the expected range */
  P4EST_ASSERT (sphere->type == P8EST_GEOMETRY_BUILTIN_SPHERE);
  P4EST_ASSERT (0 <= which_tree && which_tree < 13);
  P4EST_ASSERT (abc[0] < 1.0 + SC_1000_EPS && abc[0] > -1.0 - SC_1000_EPS);
  P4EST_ASSERT (abc[1] < 1.0 + SC_1000_EPS && abc[1] > -1.0 - SC_1000_EPS);
#ifdef P4EST_ENABLE_DEBUG
  if (which_tree < 12) {
    P4EST_ASSERT (abc[2] < 2.0 + SC_1000_EPS && abc[2] > 1.0 - SC_1000_EPS);
  }
  else {
    P4EST_ASSERT (abc[2] < 1.0 + SC_1000_EPS && abc[2] > -1.0 - SC_1000_EPS);
  }
#endif /* P4EST_ENABLE_DEBUG */

  if (which_tree < 6) {         /* outer shell */
    const double        z_cmb = abc[2] - (1. + 5. / 8.);
    const double        dist = 1. / 8.; /* keep it inside the tree */

    x = tan (abc[0] * M_PI_4);
    y = tan (abc[1] * M_PI_4);
    if (fabs (z_cmb) < dist) {
      /* correct z grading for the PREM model */
      const double        correction = 0.008873;

      R = sphere->R1sqrbyR2 * pow (sphere->R2byR1,
                                   abc[2] + correction *
                                   exp (1. / (dist * dist) -
                                        1. / ((z_cmb + dist) *
                                              (dist - z_cmb))));
    }
    else {
      R = sphere->R1sqrbyR2 * pow (sphere->R2byR1, abc[2]);
    }
    q = R / sqrt (x * x + y * y + 1.);
  }
  else if (which_tree < 12) {   /* inner shell */
    double              p, tanx, tany;

    p = 2. - abc[2];
    tanx = tan (abc[0] * M_PI_4);
    tany = tan (abc[1] * M_PI_4);
    x = p * abc[0] + (1. - p) * tanx;
    y = p * abc[1] + (1. - p) * tany;
    R = sphere->R0sqrbyR1 * pow (sphere->R1byR0, abc[2]);
    q = R / sqrt (1. + (1. - p) * (tanx * tanx + tany * tany) + 2. * p);
  }
  else {                        /* center cube */
    xyz[0] = abc[0] * sphere->Clength;
    xyz[1] = abc[1] * sphere->Clength;
    xyz[2] = abc[2] * sphere->Clength;

    return;
  }

  /* assign correct coordinates based on direction */
  switch (which_tree % 6) {
  case 0:                      /* front */
    xyz[0] = +q * x;
    xyz[1] = -q;
    xyz[2] = +q * y;
    break;
  case 1:                      /* top */
    xyz[0] = +q * x;
    xyz[1] = +q * y;
    xyz[2] = +q;
    break;
  case 2:                      /* back */
    xyz[0] = +q * x;
    xyz[1] = +q;
    xyz[2] = -q * y;
    break;
  case 3:                      /* right */
    xyz[0] = +q;
    xyz[1] = -q * x;
    xyz[2] = -q * y;
    break;
  case 4:                      /* bottom */
    xyz[0] = -q * y;
    xyz[1] = -q * x;
    xyz[2] = -q;
    break;
  case 5:                      /* left */
    xyz[0] = -q;
    xyz[1] = -q * x;
    xyz[2] = +q * y;
    break;
  default:
    SC_ABORT_NOT_REACHED ();
  }
}

p8est_geometry_t   *
p8est_geometry_new_sphere (p8est_connectivity_t * conn,
                           double R2, double R1, double R0)
{
  p8est_geometry_builtin_t *builtin;
  struct p8est_geometry_builtin_sphere *sphere;

  builtin = P4EST_ALLOC_ZERO (p8est_geometry_builtin_t, 1);

  sphere = &builtin->p.sphere;
  sphere->type = P8EST_GEOMETRY_BUILTIN_SPHERE;
  sphere->R2 = R2;
  sphere->R1 = R1;
  sphere->R0 = R0;

  /* variables useful for the outer shell */
  sphere->R2byR1 = R2 / R1;
  sphere->R1sqrbyR2 = R1 * R1 / R2;
  sphere->R1log = log (R2 / R1);

  /* variables useful for the inner shell */
  sphere->R1byR0 = R1 / R0;
  sphere->R0sqrbyR1 = R0 * R0 / R1;
  sphere->R0log = log (R1 / R0);

  /* variables useful for the center cube */
  sphere->Clength = R0 / sqrt (3.);
  sphere->CdetJ = pow (R0 / sqrt (3.), 3.);

  builtin->geom.name = "p8est_sphere";
  builtin->geom.user = conn;
  builtin->geom.X = p8est_geometry_sphere_X;

  return (p8est_geometry_t *) builtin;
}

/**
 * geometric coordinate transformation for torus of revolution geometry.
 *
 * Define the geometric transformation from logical space (where AMR
 * is performed) to the physical space.
 *
 * \param[in]  p4est      the forest
 * \param[in]  which_tree tree id inside forest
 * \param[in]  rst        coordinates in AMR space : [0,1]^3
 * \param[out] xyz        cartesian coordinates in physical space after geometry
 *
 * Note abc[3] contains cartesian coordinates in logical
 * vertex space (before geometry).
 */
static void
p8est_geometry_torus_X (p8est_geometry_t * geom,
                        p4est_topidx_t which_tree,
                        const double rst[3], double xyz[3])
{
  const struct p8est_geometry_builtin_torus *torus
    = &((p8est_geometry_builtin_t *) geom)->p.torus;
  double              x, y, R, q;
  double              abc[3];

  (void) y;

  /* transform from the reference cube [0,1]^3 into logical vertex space
     using bi/trilinear transformation */
  p8est_geometry_connectivity_X (geom, which_tree, rst, abc);

  /*
   * assert that input points are in the expected range
   * Note: maybe we should remove these assert, this would allow
   * ghost quadrant at external boundary to call this routine ?
   */

  p4est_topidx_t      which_tree_local = which_tree % 5;

  P4EST_ASSERT (torus->type == P8EST_GEOMETRY_BUILTIN_TORUS);
  P4EST_ASSERT (abc[0] < 1.0 + SC_1000_EPS && abc[0] > -1.0 - SC_1000_EPS);
  if (which_tree_local < 4)
    P4EST_ASSERT (abc[1] < 2.0 + SC_1000_EPS && abc[1] > 1.0 - SC_1000_EPS);
  else
    P4EST_ASSERT (abc[1] < 1.0 + SC_1000_EPS && abc[1] > -1.0 - SC_1000_EPS);

  /* abc[2] is always 0 here and so unused in 2D ... */

  if (which_tree_local < 4) {
    double              p, tanx;

    p = 2.0 - abc[1];
    tanx = -tan (abc[0] * M_PI_4);      /* x = tan (theta) */

    x = p * (-abc[0]) + (1. - p) * tanx;

    /* compute transformation ingredients */
    R = torus->R0sqrbyR1 * pow (torus->R1byR0, abc[1]);

    /* R*cos(theta) */
    /* q = R / sqrt (x * x + 1.); */
    q = R / sqrt (1. + (1. - p) * (tanx * tanx) + 1. * p);

    /* assign correct coordinates based on patch id */

    switch (which_tree_local) {
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
    xyz[0] = abc[0] * torus->Clength;
    xyz[1] = abc[1] * torus->Clength;

  }

  /* translation */
  xyz[0] += torus->R2;

  /* rotate around Y-axis */
  {

    int                 iSegment = which_tree / 5;

    double              tmp = xyz[0];
    double              phi =
      2 * M_PI / torus->nSegments * (iSegment + abc[2]);
    xyz[0] = tmp * cos (phi);
    xyz[2] = tmp * sin (phi);
  }

}                               /* p8est_geometry_torus_X */

p8est_geometry_t   *
p8est_geometry_new_torus (p8est_connectivity_t * conn,
                          double R0, double R1, double R2)
{
  p8est_geometry_builtin_t *builtin;
  struct p8est_geometry_builtin_torus *torus;

  builtin = P4EST_ALLOC_ZERO (p8est_geometry_builtin_t, 1);

  torus = &builtin->p.torus;
  torus->type = P8EST_GEOMETRY_BUILTIN_TORUS;
  torus->R0 = R0;
  torus->R1 = R1;
  torus->R2 = R2;
  torus->nSegments = conn->num_trees / 5;

  /* variables useful for the outer shell */
  torus->R1byR0 = R1 / R0;
  torus->R0sqrbyR1 = R0 * R0 / R1;
  torus->R0log = log (R1 / R0);

  /* variables useful for the center square */
  torus->Clength = R0 / sqrt (2.);
  /* torus->CdetJ = pow (R0 / sqrt (3.), 3.); */

  builtin->geom.name = "p8est_torus";
  builtin->geom.user = conn;
  builtin->geom.X = p8est_geometry_torus_X;

  return (p8est_geometry_t *) builtin;

}                               /* p8est_geometry_new_torus */
