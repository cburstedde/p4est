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
      double              Clength, CdetJ;
    }
    sphere;
  }
  p;
}
p8est_geometry_builtin_t;

int                 p8est_geometry_max_newton = 20;

double
p8est_geometry_Jit (p8est_geometry_t * geom,
                    p4est_topidx_t which_tree,
                    const double abc[3], double Jit[3][3])
{
  double              J[3][3];
  double              detJ, idetJ;

  idetJ = 1. / (detJ = geom->J (geom, which_tree, abc, J));

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

int
p8est_geometry_I (p8est_geometry_t * geom, p4est_topidx_t which_tree,
                  const double txyz[3], double cabc[8][3],
                  double abc[3], double rst[3])
{
  int                 i, j, k;
  int                 ri, si, ti;
  double              w[3], factor, xyz[3], rhs[3], residual;
  double              Jit[3][3], habc[3], d[3], dfactor[3], AR[3][3];
  double              idetAR, ARit[3][3], step[3], steplength;
  const double       *vkd;

  P4EST_VERBOSEF ("Target location XYZ %g %g %g tree %lld\n",
                  txyz[0], txyz[1], txyz[2], (long long) which_tree);

  /* use the center of this octree as initial guess */
  rst[0] = rst[1] = rst[2] = 0.;

  /* run Newton's method on the geometry transformation */
  residual = 1.;
  steplength = 1.;
  for (k = 0;; ++k) {
    P4EST_LDEBUGF ("Guess %d is RST %g %g %g\n", k, rst[0], rst[1], rst[2]);

    /* transform reference coordinate into physical space */
    abc[0] = abc[1] = abc[2] = 0.;
    for (ti = 0; ti < 2; ++ti) {
      w[2] = 1. + (2 * ti - 1.) * rst[2];
      for (si = 0; si < 2; ++si) {
        w[1] = 1. + (2 * si - 1.) * rst[1];
        for (ri = 0; ri < 2; ++ri) {
          w[0] = 1. + (2 * ri - 1.) * rst[0];
          vkd = cabc[4 * ti + 2 * si + ri];
          factor = .125 * w[0] * w[1] * w[2];
          abc[0] += factor * vkd[0];
          abc[1] += factor * vkd[1];
          abc[2] += factor * vkd[2];
        }
      }
    }
    geom->X (geom, which_tree, abc, xyz);

    /* compute residual in physical space */
    rhs[0] = xyz[0] - txyz[0];
    rhs[1] = xyz[1] - txyz[1];
    rhs[2] = xyz[2] - txyz[2];
    residual = sqrt (rhs[0] * rhs[0] + rhs[1] * rhs[1] + rhs[2] * rhs[2]);
    P4EST_LDEBUGF ("Guess %d is XYZ %g %g %g residual %g\n",
                   k, xyz[0], xyz[1], xyz[2], residual);

    /* delayed examination of stopping criterion */
    if (steplength < 1. / P8EST_ROOT_LEN) {
      break;
    }
    if (k > p8est_geometry_max_newton) {
      P4EST_NOTICE ("Geometry inverse transformation failed\n");
      return -1;
    }

    /* apply inverse Jacobian to residual */
    (void) geom->Jit (geom, which_tree, abc, Jit);      /* physical to tree */
    habc[0] = Jit[0][0] * rhs[0] + Jit[1][0] * rhs[1] + Jit[2][0] * rhs[2];
    habc[1] = Jit[0][1] * rhs[0] + Jit[1][1] * rhs[1] + Jit[2][1] * rhs[2];
    habc[2] = Jit[0][2] * rhs[0] + Jit[1][2] * rhs[1] + Jit[2][2] * rhs[2];

    memset (AR, 0, 3 * 3 * sizeof (double));    /* tree to reference */
    for (ti = 0; ti < 2; ++ti) {
      d[2] = 2 * ti - 1.;
      w[2] = 1. + d[2] * rst[2];
      for (si = 0; si < 2; ++si) {
        d[1] = 2 * si - 1.;
        w[1] = 1. + d[1] * rst[1];
        for (ri = 0; ri < 2; ++ri) {
          d[0] = 2 * ri - 1.;
          w[0] = 1. + d[0] * rst[0];
          vkd = cabc[4 * ti + 2 * si + ri];
          dfactor[0] = .125 * d[0] * w[1] * w[2];
          dfactor[1] = .125 * w[0] * d[1] * w[2];
          dfactor[2] = .125 * w[0] * w[1] * d[2];
          for (i = 0; i < 3; ++i) {
            for (j = 0; j < 3; ++j) {
              AR[i][j] += dfactor[j] * vkd[i];
            }
          }
        }
      }
    }
    idetAR = 1. / (AR[0][0] * (AR[1][1] * AR[2][2] - AR[1][2] * AR[2][1]) +
                   AR[0][1] * (AR[1][2] * AR[2][0] - AR[1][0] * AR[2][2]) +
                   AR[0][2] * (AR[1][0] * AR[2][1] - AR[1][1] * AR[2][0]));

    ARit[0][0] = (AR[1][1] * AR[2][2] - AR[1][2] * AR[2][1]) * idetAR;
    ARit[0][1] = (AR[1][2] * AR[2][0] - AR[1][0] * AR[2][2]) * idetAR;
    ARit[0][2] = (AR[1][0] * AR[2][1] - AR[1][1] * AR[2][0]) * idetAR;

    ARit[1][0] = (AR[0][2] * AR[2][1] - AR[0][1] * AR[2][2]) * idetAR;
    ARit[1][1] = (AR[0][0] * AR[2][2] - AR[0][2] * AR[2][0]) * idetAR;
    ARit[1][2] = (AR[0][1] * AR[2][0] - AR[0][0] * AR[2][1]) * idetAR;

    ARit[2][0] = (AR[0][1] * AR[1][2] - AR[1][1] * AR[0][2]) * idetAR;
    ARit[2][1] = (AR[0][2] * AR[1][0] - AR[1][2] * AR[0][0]) * idetAR;
    ARit[2][2] = (AR[0][0] * AR[1][1] - AR[1][0] * AR[0][1]) * idetAR;

    /* compute step in reference tree */
    step[0] =
      ARit[0][0] * habc[0] + ARit[1][0] * habc[1] + ARit[2][0] * habc[2];
    step[1] =
      ARit[0][1] * habc[0] + ARit[1][1] * habc[1] + ARit[2][1] * habc[2];
    step[2] =
      ARit[0][2] * habc[0] + ARit[1][2] * habc[1] + ARit[2][2] * habc[2];
    steplength =
      sqrt (step[0] * step[0] + step[1] * step[1] + step[2] * step[2]);
    P4EST_LDEBUGF ("Step %d is RST %g %g %g length %g\n", k,
                   -step[0], -step[1], -step[2], steplength);

    rst[0] -= step[0];
    rst[1] -= step[1];
    rst[2] -= step[2];
    P4EST_LDEBUGF ("Attempt %d is RST %g %g %g\n", k, rst[0], rst[1], rst[2]);

    /* make sure we don't leave the reference tree */
    rst[0] = SC_MIN (rst[0], +1.);
    rst[0] = SC_MAX (rst[0], -1.);
    rst[1] = SC_MIN (rst[1], +1.);
    rst[1] = SC_MAX (rst[1], -1.);
    rst[2] = SC_MIN (rst[2], +1.);
    rst[2] = SC_MAX (rst[2], -1.);
  }

  return k;
}

void
p8est_geometry_identity_X (p8est_geometry_t * geom,
                           p4est_topidx_t which_tree,
                           const double abc[3], double xyz[3])
{
  memcpy (xyz, abc, 3 * sizeof (double));
}

double
p8est_geometry_identity_D (p8est_geometry_t * geom,
                           p4est_topidx_t which_tree, const double abc[3])
{
  return 1.;
}

double
p8est_geometry_identity_J (p8est_geometry_t * geom,
                           p4est_topidx_t which_tree,
                           const double abc[3], double J[3][3])
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

  geom->name = "p8est:identity";
  geom->X = p8est_geometry_identity_X;
  geom->D = p8est_geometry_identity_D;
  geom->J = geom->Jit = p8est_geometry_identity_J;      /* identical here */
  geom->I = p8est_geometry_I;

  return geom;
}

static void
p8est_geometry_shell_X (p8est_geometry_t * geom,
                        p4est_topidx_t which_tree,
                        const double abc[3], double xyz[3])
{
  const struct p8est_geometry_builtin_shell *shell
    = &((p8est_geometry_builtin_t *) geom)->p.shell;
  double              x, y, R, q;

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
    SC_CHECK_NOT_REACHED ();
  }
}

static double
p8est_geometry_shell_D (p8est_geometry_t * geom,
                        p4est_topidx_t which_tree, const double abc[3])
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
  P4EST_ASSERT (abc[0] < 1.0 + SC_1000_EPS && abc[0] > -1.0 - SC_1000_EPS);
  P4EST_ASSERT (abc[1] < 1.0 + SC_1000_EPS && abc[1] > -1.0 - SC_1000_EPS);
  P4EST_ASSERT (abc[2] < 2.0 + SC_1000_EPS && abc[2] > 1.0 - SC_1000_EPS);

  /* transform x and y in-place for nicer grading */
  cx = cos (abc[0] * M_PI_4);
  derx = M_PI_4 / (cx * cx);
  x = tan (abc[0] * M_PI_4);
  cy = cos (abc[1] * M_PI_4);
  dery = M_PI_4 / (cy * cy);
  y = tan (abc[1] * M_PI_4);

  /* compute transformation ingredients */
  R = shell->R1sqrbyR2 * pow (shell->R2byR1, abc[2]);
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
    * q * q * q * derx * dery * shell->Rlog;
  P4EST_ASSERT (detJ > 0.);

  return detJ;
}

static double
p8est_geometry_shell_J (p8est_geometry_t * geom,
                        p4est_topidx_t which_tree,
                        const double abc[3], double J[3][3])
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
  P4EST_ASSERT (abc[0] < 1.0 + SC_1000_EPS && abc[0] > -1.0 - SC_1000_EPS);
  P4EST_ASSERT (abc[1] < 1.0 + SC_1000_EPS && abc[1] > -1.0 - SC_1000_EPS);
  P4EST_ASSERT (abc[2] < 2.0 + SC_1000_EPS && abc[2] > 1.0 - SC_1000_EPS);

  /* transform x and y in-place for nicer grading */
  cx = cos (abc[0] * M_PI_4);
  derx = M_PI_4 / (cx * cx);
  x = tan (abc[0] * M_PI_4);
  cy = cos (abc[1] * M_PI_4);
  dery = M_PI_4 / (cy * cy);
  y = tan (abc[1] * M_PI_4);

  /* compute transformation ingredients */
  R = shell->R1sqrbyR2 * pow (shell->R2byR1, abc[2]);
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

  builtin->geom.name = "p8est:shell";
  builtin->geom.X = p8est_geometry_shell_X;
  builtin->geom.D = p8est_geometry_shell_D;
  builtin->geom.J = p8est_geometry_shell_J;
  builtin->geom.Jit = p8est_geometry_Jit;
  builtin->geom.I = p8est_geometry_I;

  return (p8est_geometry_t *) builtin;
}

static void
p8est_geometry_sphere_X (p8est_geometry_t * geom,
                         p4est_topidx_t which_tree,
                         const double abc[3], double xyz[3])
{
  const struct p8est_geometry_builtin_sphere *sphere
    = &((p8est_geometry_builtin_t *) geom)->p.sphere;
  double              x, y, R, q;

  /* assert that input points are in the expected range */
  P4EST_ASSERT (sphere->type == P8EST_GEOMETRY_BUILTIN_SPHERE);
  P4EST_ASSERT (0 <= which_tree && which_tree < 13);
  P4EST_ASSERT (abc[0] < 1.0 + SC_1000_EPS && abc[0] > -1.0 - SC_1000_EPS);
  P4EST_ASSERT (abc[1] < 1.0 + SC_1000_EPS && abc[1] > -1.0 - SC_1000_EPS);
#ifdef P4EST_DEBUG
  if (which_tree < 12) {
    P4EST_ASSERT (abc[2] < 2.0 + SC_1000_EPS && abc[2] > 1.0 - SC_1000_EPS);
  }
  else {
    P4EST_ASSERT (abc[2] < 1.0 + SC_1000_EPS && abc[2] > -1.0 - SC_1000_EPS);
  }
#endif /* P4EST_DEBUG */

  if (which_tree < 6) {         /* outer shell */
    x = tan (abc[0] * M_PI_4);
    y = tan (abc[1] * M_PI_4);
    R = sphere->R1sqrbyR2 * pow (sphere->R2byR1, abc[2]);
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
    SC_CHECK_NOT_REACHED ();
  }
}

static double
p8est_geometry_sphere_D (p8est_geometry_t * geom,
                         p4est_topidx_t which_tree, const double abc[3])
{
  const struct p8est_geometry_builtin_sphere *sphere
    = &((p8est_geometry_builtin_t *) geom)->p.sphere;
  double              Rlog;
  double              cx, cy, x, y, R, t, q;
  double              derx, dery;
  double              factor, detJ;
  double              J[3][3];

  /* assert that input points are in the expected range */
  P4EST_ASSERT (sphere->type == P8EST_GEOMETRY_BUILTIN_SPHERE);
  P4EST_ASSERT (0 <= which_tree && which_tree < 13);
  P4EST_ASSERT (abc[0] < 1.0 + SC_1000_EPS && abc[0] > -1.0 - SC_1000_EPS);
  P4EST_ASSERT (abc[1] < 1.0 + SC_1000_EPS && abc[1] > -1.0 - SC_1000_EPS);
#ifdef P4EST_DEBUG
  if (which_tree < 12) {
    P4EST_ASSERT (abc[2] < 2.0 + SC_1000_EPS && abc[2] > 1.0 - SC_1000_EPS);
  }
  else {
    P4EST_ASSERT (abc[2] < 1.0 + SC_1000_EPS && abc[2] > -1.0 - SC_1000_EPS);
  }
#endif /* P4EST_DEBUG */

  if (which_tree < 6) {         /* outer shell */
    cx = cos (abc[0] * M_PI_4);
    derx = M_PI_4 / (cx * cx);
    x = tan (abc[0] * M_PI_4);

    cy = cos (abc[1] * M_PI_4);
    dery = M_PI_4 / (cy * cy);
    y = tan (abc[1] * M_PI_4);

    R = sphere->R1sqrbyR2 * pow (sphere->R2byR1, abc[2]);
    t = 1. / (x * x + y * y + 1.);
    q = R * sqrt (t);
    Rlog = sphere->R1log;

    J[0][0] = (1. - x * x * t);
    J[0][1] = -x * y * t;
    J[0][2] = x;
    J[1][0] = -x * y * t;
    J[1][1] = (1. - y * y * t);
    J[1][2] = y;
    J[2][0] = -x * t;
    J[2][1] = -y * t;
    J[2][2] = 1.;
    factor = q * q * q * derx * dery * Rlog;
  }
  else if (which_tree < 12) {   /* inner shell */
    double              p, tanx, tany, tsqr;

    p = 2. - abc[2];

    cx = cos (abc[0] * M_PI_4);
    derx = (1. - p) * M_PI_4 / (cx * cx);
    tanx = tan (abc[0] * M_PI_4);
    x = p * abc[0] + (1. - p) * tanx;

    cy = cos (abc[1] * M_PI_4);
    dery = (1. - p) * M_PI_4 / (cy * cy);
    tany = tan (abc[1] * M_PI_4);
    y = p * abc[1] + (1. - p) * tany;

    R = sphere->R0sqrbyR1 * pow (sphere->R1byR0, abc[2]);
    tsqr = tanx * tanx + tany * tany;
    t = 1. / (1. + (1. - p) * tsqr + 2. * p);
    q = R * sqrt (t);
    Rlog = sphere->R0log + t * (1. - .5 * tsqr);

    J[0][0] = p + (1. - x * tanx * t) * derx;
    J[0][1] = -x * tany * t * dery;
    J[0][2] = x * Rlog - abc[0] + tanx;
    J[1][0] = -y * tanx * t * derx;
    J[1][1] = p + (1. - y * tany * t) * dery;
    J[1][2] = y * Rlog - abc[1] + tany;
    J[2][0] = -tanx * t * derx;
    J[2][1] = -tany * t * dery;
    J[2][2] = Rlog;
    factor = q * q * q;
  }
  else {                        /* center cube */
    return sphere->CdetJ;
  }

  /* compute the determinant */
  detJ = (J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1])
          + J[0][1] * (J[1][2] * J[2][0] - J[1][0] * J[2][2])
          + J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][0]))
    * factor;
  P4EST_ASSERT (detJ > 0.);

  return detJ;
}

static double
p8est_geometry_sphere_J (p8est_geometry_t * geom,
                         p4est_topidx_t which_tree,
                         const double abc[3], double J[3][3])
{
  /* *INDENT-OFF* HORRIBLE indent bug */
  const struct p8est_geometry_builtin_sphere *sphere
    = &((p8est_geometry_builtin_t *) geom)->p.sphere;
  const int           mapJ[6][3] =
    {{  0,  2,  1, },
     {  0,  1,  2, },
     {  0,  2,  1, },
     {  1,  2,  0, },
     {  1,  0,  2, },
     {  1,  2,  0, }};
  const double        mapM[6][3] =
    {{  1,  1, -1, },
     {  1,  1,  1, },
     {  1, -1,  1, },
     { -1, -1,  1, },
     { -1, -1, -1, },
     { -1,  1, -1, }};
  /* *INDENT-ON* */
  int                 pid, j0, j1, j2;
  double              Rlog;
  double              cx, cy, x, y, R, t, q;
  double              q0, q1, q2;
  double              derx, dery;
  double              detJ;

  /* assert that input points are in the expected range */
  P4EST_ASSERT (sphere->type == P8EST_GEOMETRY_BUILTIN_SPHERE);
  P4EST_ASSERT (0 <= which_tree && which_tree < 13);
  P4EST_ASSERT (abc[0] < 1.0 + SC_1000_EPS && abc[0] > -1.0 - SC_1000_EPS);
  P4EST_ASSERT (abc[1] < 1.0 + SC_1000_EPS && abc[1] > -1.0 - SC_1000_EPS);
#ifdef P4EST_DEBUG
  if (which_tree < 12) {
    P4EST_ASSERT (abc[2] < 2.0 + SC_1000_EPS && abc[2] > 1.0 - SC_1000_EPS);
  }
  else {
    P4EST_ASSERT (abc[2] < 1.0 + SC_1000_EPS && abc[2] > -1.0 - SC_1000_EPS);
  }
#endif /* P4EST_DEBUG */

  if (which_tree < 6) {         /* outer shell */
    cx = cos (abc[0] * M_PI_4);
    derx = M_PI_4 / (cx * cx);
    x = tan (abc[0] * M_PI_4);

    cy = cos (abc[1] * M_PI_4);
    dery = M_PI_4 / (cy * cy);
    y = tan (abc[1] * M_PI_4);

    R = sphere->R1sqrbyR2 * pow (sphere->R2byR1, abc[2]);
    t = 1. / (x * x + y * y + 1.);
    q = R * sqrt (t);
    Rlog = sphere->R1log;

    pid = (int) which_tree;
    j0 = mapJ[pid][0];
    j1 = mapJ[pid][1];
    j2 = mapJ[pid][2];
    q0 = mapM[pid][0] * q;
    q1 = mapM[pid][1] * q;
    q2 = mapM[pid][2] * q;
    J[j0][0] = q0 * (1. - x * x * t) * derx;
    J[j0][1] = -q0 * x * y * t * dery;
    J[j0][2] = q0 * x * Rlog;
    J[j1][0] = -q1 * x * y * t * derx;
    J[j1][1] = q1 * (1. - y * y * t) * dery;
    J[j1][2] = q1 * y * Rlog;
    J[j2][0] = -q2 * x * t * derx;
    J[j2][1] = -q2 * y * t * dery;
    J[j2][2] = q2 * Rlog;
  }
  else if (which_tree < 12) {   /* inner shell */
    double              p, tanx, tany, tsqr;

    p = 2. - abc[2];

    cx = cos (abc[0] * M_PI_4);
    derx = (1. - p) * M_PI_4 / (cx * cx);
    tanx = tan (abc[0] * M_PI_4);
    x = p * abc[0] + (1. - p) * tanx;

    cy = cos (abc[1] * M_PI_4);
    dery = (1. - p) * M_PI_4 / (cy * cy);
    tany = tan (abc[1] * M_PI_4);
    y = p * abc[1] + (1. - p) * tany;

    R = sphere->R0sqrbyR1 * pow (sphere->R1byR0, abc[2]);
    tsqr = tanx * tanx + tany * tany;
    t = 1. / (1. + (1. - p) * tsqr + 2. * p);
    q = R * sqrt (t);
    Rlog = sphere->R0log + t * (1. - .5 * tsqr);

    pid = (int) which_tree - 6;
    j0 = mapJ[pid][0];
    j1 = mapJ[pid][1];
    j2 = mapJ[pid][2];
    q0 = mapM[pid][0] * q;
    q1 = mapM[pid][1] * q;
    q2 = mapM[pid][2] * q;
    J[j0][0] = q0 * (p + (1. - x * tanx * t) * derx);
    J[j0][1] = -q0 * x * tany * t * dery;
    J[j0][2] = q0 * (x * Rlog - abc[0] + tanx);
    J[j1][0] = -q1 * y * tanx * t * derx;
    J[j1][1] = q1 * (p + (1. - y * tany * t) * dery);
    J[j1][2] = q1 * (y * Rlog - abc[1] + tany);
    J[j2][0] = -q2 * tanx * t * derx;
    J[j2][1] = -q2 * tany * t * dery;
    J[j2][2] = q2 * Rlog;
  }
  else {                        /* center cube */
    J[0][0] = J[1][1] = J[2][2] = sphere->Clength;
    J[0][1] = J[1][2] = J[2][0] = 0.;
    J[1][0] = J[2][1] = J[0][2] = 0.;

    return sphere->CdetJ;
  }

  /* compute the determinant */
  detJ = J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1])
    + J[0][1] * (J[1][2] * J[2][0] - J[1][0] * J[2][2])
    + J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][0]);
  P4EST_ASSERT (detJ > 0.);

  return detJ;
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
  sphere->CdetJ = pow (R0 / sqrt (3.), 3.);

  builtin->geom.name = "p8est:sphere";
  builtin->geom.X = p8est_geometry_sphere_X;
  builtin->geom.D = p8est_geometry_sphere_D;
  builtin->geom.J = p8est_geometry_sphere_J;
  builtin->geom.Jit = p8est_geometry_Jit;
  builtin->geom.I = p8est_geometry_I;

  return (p8est_geometry_t *) builtin;
}
