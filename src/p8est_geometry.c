/*
  This file is part of p4est.
  p4est is a C library to manage a parallel collection of quadtrees and/or
  octrees.

  Copyright (C) 2007-2009 Carsten Burstedde, Lucas Wilcox.

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

#include <p8est_geometry.h>

typedef enum
{
  P8EST_GEOMETRY_BUILTIN_MAGIC = 0x65F2F8DF,
  P8EST_GEOMETRY_BUILTIN_SHELL,
  P8EST_GEOMETRY_BUILTIN_SPHERE,
}
p8est_geometry_builtin_type_t;

typedef struct p8est_geometry_builtin
{
  p8est_geometry_t    geom;
  union
  {
    p8est_geometry_builtin_type_t type;
    struct p8est_geometry_builtin_shell
    {
      p8est_geometry_builtin_type_t type;
      double              R2, R1;
      double              R2byR1, R1sqrbyR2, Rlog;
    }
    shell;
    struct p8est_geometry_builtin_sphere
    {
      p8est_geometry_builtin_type_t type;
      double              R2, R1, R0;
      double              R2byR1, R1sqrbyR2, R1log;
      double              R1byR0, R0sqrbyR1, R0log;
      double              Clength;
    }
    sphere;
  }
  p;
}
p8est_geometry_builtin_t;

double
p8est_geometry_Jit (p8est_geometry_t * geom,
                    p4est_topidx_t which_tree,
                    const double xyz[3], double Jit[3][3])
{
  double              J[3][3];
  double              detJ, idetJ;

  idetJ = 1. / (detJ = geom->J (geom, which_tree, xyz, J));

  Jit[0][0] = (J[1][1] * J[2][2] - J[1][2] * J[2][1]) * idetJ;
  Jit[0][1] = (J[1][2] * J[2][0] - J[1][0] * J[2][2]) * idetJ;
  Jit[0][2] = (J[1][0] * J[2][1] - J[1][1] * J[2][0]) * idetJ;

  Jit[1][0] = (J[0][2] * J[2][1] - J[0][1] * J[2][2]) * idetJ;
  Jit[1][1] = (J[0][0] * J[2][2] - J[0][2] * J[2][0]) * idetJ;
  Jit[1][2] = (J[0][1] * J[2][0] - J[0][0] * J[2][1]) * idetJ;

  Jit[2][0] = (J[0][1] * J[1][2] - J[1][1] * J[0][2]) * idetJ;
  Jit[2][1] = (J[0][2] * J[1][0] - J[1][2] * J[0][0]) * idetJ;
  Jit[2][2] = (J[0][0] * J[1][1] - J[1][0] * J[0][1]) * idetJ;

  return detJ;
}

void
p8est_geometry_identity_X (p8est_geometry_t * geom,
                           p4est_topidx_t which_tree,
                           const double xyz[3], double XYZ[3])
{
  memcpy (XYZ, xyz, 3 * sizeof (double));
}

double
p8est_geometry_identity_D (p8est_geometry_t * geom,
                           p4est_topidx_t which_tree, const double xyz[3])
{
  return 1.;
}

double
p8est_geometry_identity_J (p8est_geometry_t * geom,
                           p4est_topidx_t which_tree,
                           const double xyz[3], double J[3][3])
{
  J[0][0] = J[1][1] = J[2][2] = 1.;
  J[0][1] = J[1][2] = J[2][0] = 0.;
  J[1][0] = J[2][1] = J[0][2] = 0.;

  return 1.;
}

p8est_geometry_t   *
p8est_geometry_new_identity (void)
{
  p8est_geometry_t   *geom;

  geom = P4EST_ALLOC (p8est_geometry_t, 1);

  geom->X = p8est_geometry_identity_X;
  geom->D = p8est_geometry_identity_D;
  geom->J = geom->Jit = p8est_geometry_identity_J;      /* identical here */

  return geom;
}

static void
p8est_geometry_shell_X (p8est_geometry_t * geom,
                        p4est_topidx_t which_tree,
                        const double xyz[3], double XYZ[3])
{
  const struct p8est_geometry_builtin_shell *shell
    = &((p8est_geometry_builtin_t *) geom)->p.shell;
  double              x, y, R, q;

  /* assert that input points are in the expected range */
  P4EST_ASSERT (shell->type == P8EST_GEOMETRY_BUILTIN_SHELL);
  P4EST_ASSERT (0 <= which_tree && which_tree < 24);
  P4EST_ASSERT (xyz[0] < 1.0 + 1e-12 && xyz[0] > -1.0 - 1e-12);
  P4EST_ASSERT (xyz[1] < 1.0 + 1e-12 && xyz[1] > -1.0 - 1e-12);
  P4EST_ASSERT (xyz[2] < 2.0 + 1e-12 && xyz[2] > 1.0 - 1e-12);

  /* transform xyz[0] and y in-place for nicer grading */
  x = tan (xyz[0] * M_PI_4);
  y = tan (xyz[1] * M_PI_4);

  /* compute transformation ingredients */
  R = shell->R1sqrbyR2 * pow (shell->R2byR1, xyz[2]);
  q = R / sqrt (x * x + y * y + 1.);

  /* assign correct coordinates based on patch id */
  switch (which_tree / 4) {
  case 3:                      /* top */
    XYZ[0] = +q * y;
    XYZ[1] = -q * x;
    XYZ[2] = +q;
    break;
  case 2:                      /* left */
    XYZ[0] = -q;
    XYZ[1] = -q * x;
    XYZ[2] = +q * y;
    break;
  case 1:                      /* bottom */
    XYZ[0] = -q * y;
    XYZ[1] = -q * x;
    XYZ[2] = -q;
    break;
  case 0:                      /* right */
    XYZ[0] = +q;
    XYZ[1] = -q * x;
    XYZ[2] = -q * y;
    break;
  case 4:                      /* back */
    XYZ[0] = -q * x;
    XYZ[1] = +q;
    XYZ[2] = +q * y;
    break;
  case 5:                      /* front */
    XYZ[0] = +q * x;
    XYZ[1] = -q;
    XYZ[2] = +q * y;
    break;
  default:
    SC_CHECK_NOT_REACHED ();
  }
}

static double
p8est_geometry_shell_D (p8est_geometry_t * geom,
                        p4est_topidx_t which_tree, const double xyz[3])
{
  const struct p8est_geometry_builtin_shell *shell
    = &((p8est_geometry_builtin_t *) geom)->p.shell;
  double              cx, cy, x, y, R, t, q;
  double              derx, dery;
  double              detJ;
  double              J[3][3];

  /* assert that input points are in the expected range */
  P4EST_ASSERT (shell->type == P8EST_GEOMETRY_BUILTIN_SHELL);
  P4EST_ASSERT (0 <= which_tree && which_tree < 24);
  P4EST_ASSERT (xyz[0] < 1.0 + 1e-12 && xyz[0] > -1.0 - 1e-12);
  P4EST_ASSERT (xyz[1] < 1.0 + 1e-12 && xyz[1] > -1.0 - 1e-12);
  P4EST_ASSERT (xyz[2] < 2.0 + 1e-12 && xyz[2] > 1.0 - 1e-12);

  /* transform x and y in-place for nicer grading */
  cx = cos (xyz[0] * M_PI_4);
  derx = M_PI_4 / (cx * cx);
  x = tan (xyz[0] * M_PI_4);
  cy = cos (xyz[1] * M_PI_4);
  dery = M_PI_4 / (cy * cy);
  y = tan (xyz[1] * M_PI_4);

  /* compute transformation ingredients */
  R = shell->R1sqrbyR2 * pow (shell->R2byR1, xyz[2]);
  t = 1. / (x * x + y * y + 1.);
  q = R * sqrt (t);

  /* compute Jacobian in xyz space aligned to the octree modulo scaling */
  J[0][0] = (1. - x * x * t);
  J[0][1] = -x * y * t;
  J[0][2] = x;
  J[1][0] = -x * y * t;
  J[1][1] = (1. - y * y * t);
  J[1][2] = y;
  J[2][0] = -x * t;
  J[2][1] = -y * t;
  J[2][2] = 1.;

  /* compute the determinant */
  detJ = (J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1])
          + J[0][1] * (J[1][2] * J[2][0] - J[1][0] * J[2][2])
          + J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][0]))
    * q * q * q * derx * dery * Rlog;
  P4EST_ASSERT (detJ > 0.);

  return detJ;
}

static double
p8est_geometry_shell_J (p8est_geometry_t * geom,
                        p4est_topidx_t which_tree,
                        const double xyz[3], double J[3][3])
{
  const struct p8est_geometry_builtin_shell *shell
    = &((p8est_geometry_builtin_t *) geom)->p.shell;
  const double        Rlog = shell->Rlog;
  double              cx, cy, x, y, R, t, q;
  double              derx, dery;
  double              detJ;

  /* assert that input points are in the expected range */
  P4EST_ASSERT (shell->type == P8EST_GEOMETRY_BUILTIN_SHELL);
  P4EST_ASSERT (0 <= which_tree && which_tree < 24);
  P4EST_ASSERT (xyz[0] < 1.0 + 1e-12 && xyz[0] > -1.0 - 1e-12);
  P4EST_ASSERT (xyz[1] < 1.0 + 1e-12 && xyz[1] > -1.0 - 1e-12);
  P4EST_ASSERT (xyz[2] < 2.0 + 1e-12 && xyz[2] > 1.0 - 1e-12);

  /* transform x and y in-place for nicer grading */
  cx = cos (xyz[0] * M_PI_4);
  derx = M_PI_4 / (cx * cx);
  x = tan (xyz[0] * M_PI_4);
  cy = cos (xyz[1] * M_PI_4);
  dery = M_PI_4 / (cy * cy);
  y = tan (xyz[1] * M_PI_4);

  /* compute transformation ingredients */
  R = shell->R1sqrbyR2 * pow (shell->R2byR1, xyz[2]);
  t = 1. / (x * x + y * y + 1.);
  q = R * sqrt (t);

  /* compute Jacobian in xyz space aligned to the octree */
  /* assign correct coordinates based on patch id */
  switch (which_tree / 4) {
  case 3:                      /* top */
    J[0][0] = -q * x * y * t * derx;
    J[0][1] = q * (1. - y * y * t) * dery;
    J[0][2] = q * y * Rlog;
    J[1][0] = -q * (1. - x * x * t) * derx;
    J[1][1] = q * x * y * t * dery;
    J[1][2] = -q * x * Rlog;
    J[2][0] = -q * x * t * derx;
    J[2][1] = -q * y * t * dery;
    J[2][2] = q * Rlog;
    break;
  case 2:                      /* left */
    J[0][0] = q * x * t * derx;
    J[0][1] = q * y * t * dery;
    J[0][2] = -q * Rlog;
    J[1][0] = -q * (1. - x * x * t) * derx;
    J[1][1] = q * x * y * t * dery;
    J[1][2] = -q * x * Rlog;
    J[2][0] = -q * x * y * t * derx;
    J[2][1] = q * (1. - y * y * t) * dery;
    J[2][2] = q * y * Rlog;
    break;
  case 1:                      /* bottom */
    J[0][0] = q * x * y * t * derx;
    J[0][1] = -q * (1. - y * y * t) * dery;
    J[0][2] = -q * y * Rlog;
    J[1][0] = -q * (1. - x * x * t) * derx;
    J[1][1] = q * x * y * t * dery;
    J[1][2] = -q * x * Rlog;
    J[2][0] = q * x * t * derx;
    J[2][1] = q * y * t * dery;
    J[2][2] = -q * Rlog;
    break;
  case 0:                      /* right */
    J[0][0] = -q * x * t * derx;
    J[0][1] = -q * y * t * dery;
    J[0][2] = q * Rlog;
    J[1][0] = -q * (1. - x * x * t) * derx;
    J[1][1] = q * x * y * t * dery;
    J[1][2] = -q * x * Rlog;
    J[2][0] = q * x * y * t * derx;
    J[2][1] = -q * (1. - y * y * t) * dery;
    J[2][2] = -q * y * Rlog;
    break;
  case 4:                      /* back */
    J[0][0] = -q * (1. - x * x * t) * derx;
    J[0][1] = q * x * y * t * dery;
    J[0][2] = -q * x * Rlog;
    J[1][0] = -q * x * t * derx;
    J[1][1] = -q * y * t * dery;
    J[1][2] = q * Rlog;
    J[2][0] = -q * x * y * t * derx;
    J[2][1] = q * (1. - y * y * t) * dery;
    J[2][2] = q * y * Rlog;
    break;
  case 5:                      /* front */
    J[0][0] = q * (1. - x * x * t) * derx;
    J[0][1] = -q * x * y * t * dery;
    J[0][2] = q * x * Rlog;
    J[1][0] = q * x * t * derx;
    J[1][1] = q * y * t * dery;
    J[1][2] = -q * Rlog;
    J[2][0] = -q * x * y * t * derx;
    J[2][1] = q * (1. - y * y * t) * dery;
    J[2][2] = q * y * Rlog;
    break;
  default:
    SC_CHECK_NOT_REACHED ();
  }

  /* compute the determinant */
  detJ = J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1])
    + J[0][1] * (J[1][2] * J[2][0] - J[1][0] * J[2][2])
    + J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][0]);
  P4EST_ASSERT (detJ > 0.);

  return detJ;
}

p8est_geometry_t   *
p8est_geometry_new_shell (double R2, double R1)
{
  p8est_geometry_builtin_t *builtin;
  struct p8est_geometry_builtin_shell *shell;

  builtin = P4EST_ALLOC (p8est_geometry_builtin_t, 1);

  shell = &builtin->p.shell;
  shell->type = P8EST_GEOMETRY_BUILTIN_SHELL;
  shell->R2 = R2;
  shell->R1 = R1;
  shell->R2byR1 = R2 / R1;
  shell->R1sqrbyR2 = R1 * R1 / R2;
  shell->Rlog = log (R2 / R1);

  builtin->geom.X = p8est_geometry_shell_X;
  builtin->geom.D = p8est_geometry_shell_D;
  builtin->geom.J = p8est_geometry_shell_J;
  builtin->geom.Jit = p8est_geometry_Jit;

  return (p8est_geometry_t *) builtin;
}

static void
p8est_geometry_sphere_X (p8est_geometry_t * geom,
                         p4est_topidx_t which_tree,
                         const double xyz[3], double XYZ[3])
{
  const struct p8est_geometry_builtin_sphere *sphere
    = &((p8est_geometry_builtin_t *) geom)->p.sphere;
  double              x, y, R, q;
  double              p, tpx, tpy;

  /* assert that input points are in the expected range */
  P4EST_ASSERT (sphere->type == P8EST_GEOMETRY_BUILTIN_SPHERE);
  P4EST_ASSERT (0 <= which_tree && which_tree < 13);
  P4EST_ASSERT (xyz[0] < 1.0 + 1e-12 && xyz[0] > -1.0 - 1e-12);
  P4EST_ASSERT (xyz[1] < 1.0 + 1e-12 && xyz[1] > -1.0 - 1e-12);
#ifdef P4EST_DEBUG
  if (which_tree < 12) {
    P4EST_ASSERT (xyz[2] < 2.0 + 1e-12 && xyz[2] > 1.0 - 1e-12);
  }
  else {
    P4EST_ASSERT (xyz[2] < 1.0 + 1e-12 && xyz[2] > -1.0 - 1e-12);
  }
#endif /* P4EST_DEBUG */

  if (which_tree < 6) {         /* outer shell */
    x = tan (xyz[0] * M_PI_4);
    y = tan (xyz[1] * M_PI_4);
    R = sphere->R1sqrbyR2 * pow (sphere->R2byR1, xyz[2]);
    q = R / sqrt (x * x + y * y + 1.);
  }
  else if (which_tree < 12) {   /* inner shell */
    p = 2. - xyz[2];
    x = p * xyz[0] + (tpx = (1. - p) * tan (xyz[0] * M_PI_4));
    y = p * xyz[1] + (tpy = (1. - p) * tan (xyz[1] * M_PI_4));
    R = sphere->R0sqrbyR1 * pow (sphere->R1byR0, xyz[2]);
    q = R / sqrt ((p + tpx) * (p + tpx) + (p + tpy) * (p + tpy) + 1.);
  }
  else {                        /* center cube */
    XYZ[0] = xyz[0] * sphere->Clength;
    XYZ[1] = xyz[1] * sphere->Clength;
    XYZ[2] = xyz[2] * sphere->Clength;
    return;
  }

  /* assign correct coordinates based on direction */
  switch (which_tree % 6) {
  case 0:                      /* front */
    XYZ[0] = +q * x;
    XYZ[1] = -q;
    XYZ[2] = +q * y;
    break;
  case 1:                      /* top */
    XYZ[0] = +q * x;
    XYZ[1] = +q * y;
    XYZ[2] = +q;
    break;
  case 2:                      /* back */
    XYZ[0] = +q * x;
    XYZ[1] = +q;
    XYZ[2] = -q * y;
    break;
  case 3:                      /* right */
    XYZ[0] = +q;
    XYZ[1] = -q * x;
    XYZ[2] = -q * y;
    break;
  case 4:                      /* bottom */
    XYZ[0] = -q * y;
    XYZ[1] = -q * x;
    XYZ[2] = -q;
    break;
  case 5:                      /* left */
    XYZ[0] = -q;
    XYZ[1] = -q * x;
    XYZ[2] = +q * y;
    break;
  default:
    SC_CHECK_NOT_REACHED ();
  }
}

p8est_geometry_t   *
p8est_geometry_new_sphere (double R2, double R1, double R0)
{
  p8est_geometry_builtin_t *builtin;
  struct p8est_geometry_builtin_sphere *sphere;

  builtin = P4EST_ALLOC (p8est_geometry_builtin_t, 1);

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

  builtin->geom.X = p8est_geometry_sphere_X;
  builtin->geom.D = p8est_geometry_identity_D;
  builtin->geom.J = p8est_geometry_identity_J;
  builtin->geom.Jit = p8est_geometry_Jit;

  return (p8est_geometry_t *) builtin;
}
