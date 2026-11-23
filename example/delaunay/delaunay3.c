/*
 * Usage: p8est_delaunay [OPTIONS] [-c <configuration>] [-v <VTK basename>]
 *        possible options:
 *        -l <minlevel>
 *        -L <maxlevel>
 *        possible configurations:
 *        o unit      Refinement on the unit cube.
 *        o twocubes  Refinement on a forest with two trees.
 *        o brick235  Refinement on a 2x3x5 brick shape.
 *        o torus8    Refinement on a 40-tree volumetric torus.
 *        o periodic  Refinement on the unit cube with all-periodic b.c.
 *        o rotcubes  Refinement on six weirdly connected trees.
 *        o shell     Refinement on a 24-tree spherical shell.
 *        o sphere    Refinement on a 13-tree volumetric sphere.
 *        without the -v option, we do not write any VTK data at all.
 */
#include <p4est_to_p8est.h>
#include "delaunay2.c"
