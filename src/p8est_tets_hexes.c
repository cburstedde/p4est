/*
  This file is part of p4est.
  p4est is a C library to manage a collection (a forest) of multiple
  connected adaptive quadtrees or octrees in parallel.

  Copyright (C) 2010 The University of Texas System
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

#include <p8est_tets_hexes.h>

sc_array_t         *
p8est_tetgen_read_node (const char *nodefilename)
{
  int                 retval;
  int                 i, j;
  int                 num_nodes, dims, num_attributes, boundary_marker;
  double             *pc;
  sc_array_t         *nodes;
  FILE               *nodefile;

  /* prepare cleanup on error */
  nodefile = NULL;
  nodes = NULL;

  /* open node file */
  nodefile = fopen (nodefilename, "rb");
  if (nodefile == NULL) {
    P4EST_LERRORF ("Failed to open %s\n", nodefilename);
    goto dead;
  }

  /* read header information */
  retval = fscanf (nodefile, "%d %d %d %d",
                   &num_nodes, &dims, &num_attributes, &boundary_marker);
  if (retval != 4 || num_nodes < 0 || dims != 3 || num_attributes < 0) {
    P4EST_LERROR ("Failed to read node header\n");
    goto dead;
  }

  /* read node coordinates */
  nodes = sc_array_new_size (sizeof (double), 3 * num_nodes);
  for (i = 0; i < num_nodes; ++i) {
    pc = (double *) sc_array_index_int (nodes, 3 * i);
    retval = fscanf (nodefile, "%d %lf %lf %lf", &j, pc, pc + 1, pc + 2);
    if (retval != 4 || i != j) {
      P4EST_LERRORF ("Failed to read node %d coordinates\n", i);
      goto dead;
    }
    for (j = 0; j < num_attributes; ++j) {
      retval = fscanf (nodefile, "%*f");
      if (retval != 0) {
        P4EST_LERRORF ("Failed to read node %d attribute %d\n", i, j);
        goto dead;
      }
    }
    if (boundary_marker) {
      retval = fscanf (nodefile, "%*d");
      if (retval != 0) {
        P4EST_LERRORF ("Failed to read node %d boundary marker\n", i);
        goto dead;
      }
    }
  }

  /* close node file and return */
  retval = fclose (nodefile);
  nodefile = NULL;
  if (retval) {
    P4EST_LERRORF ("Failed to close %s\n", nodefilename);
    goto dead;
  }
  return nodes;

dead:
  /* clean up on error */
  if (nodefile != NULL) {
    fclose (nodefile);
  }
  if (nodes != NULL) {
    sc_array_destroy (nodes);
  }
  return NULL;
}

sc_array_t         *
p8est_tetgen_read_ele (const char *elefilename, int num_nodes,
                       sc_array_t ** attributes)
{
  int                 retval;
  int                 i, j;
  int                 num_tets, nodespertet, region;
  int                *pi;
  sc_array_t         *tets, *attr;
  FILE               *elefile;

  /* prepare cleanup on error */
  elefile = NULL;
  tets = NULL;
  attr = NULL;
  if (attributes != NULL) {
    *attributes = NULL;
  }

  /* open ele file */
  elefile = fopen (elefilename, "rb");
  if (elefile == NULL) {
    P4EST_LERRORF ("Failed to open %s\n", elefilename);
    goto dead;
  }

  /* read header information */
  retval = fscanf (elefile, "%d %d %d", &num_tets, &nodespertet, &region);
  if (retval != 3 || num_tets < 0 || nodespertet != 4) {
    P4EST_LERROR ("Failed to read tet header\n");
    goto dead;
  }

  /* read tet coordinates */
  tets = sc_array_new_size (sizeof (int), 4 * num_tets);
  if (region && attributes != NULL) {
    attr = *attributes = sc_array_new_size (sizeof (int), num_tets);
  }
  for (i = 0; i < num_tets; ++i) {
    pi = (int *) sc_array_index_int (tets, 4 * i);
    retval = fscanf (elefile, "%d %d %d %d %d",
                     &j, pi, pi + 1, pi + 2, pi + 3);
    if (retval != 5 || i != j) {
      P4EST_LERRORF ("Failed to read tet %d node numbers\n", i);
      goto dead;
    }
    for (j = 0; j < 4; ++j) {
      if (pi[j] < 0 || (num_nodes >= 0 && pi[j] >= num_nodes)) {
        P4EST_LERRORF ("Tet %d has invalid node number %d\n", i, j);
        goto dead;
      }
    }
    if (region) {
      if (attr != NULL) {
        pi = (int *) sc_array_index_int (attr, i);
        retval = fscanf (elefile, "%d", pi) - 1;
      }
      else {
        retval = fscanf (elefile, "%*d");
      }
      if (retval != 0) {
        P4EST_LERRORF ("Failed to read tet %d region attribute\n", i);
        goto dead;
      }
    }
  }

  /* close ele file */
  retval = fclose (elefile);
  elefile = NULL;
  if (retval) {
    P4EST_LERRORF ("Failed to close %s\n", elefilename);
    goto dead;
  }
  return tets;

dead:
  /* clean up on error */
  if (elefile != NULL) {
    fclose (elefile);
  }
  if (tets != NULL) {
    sc_array_destroy (tets);
  }
  if (attr != NULL) {
    sc_array_destroy (attr);
    *attributes = NULL;
  }
  return NULL;
}

p8est_tetgen_t     *
p8est_tetgen_read (const char *tetgenbasename)
{
  int                 num_nodes;
  char                nodefilename[BUFSIZ];
  char                elefilename[BUFSIZ];
  sc_array_t         *nodes, *tets, *attr;
  p8est_tetgen_t     *ptg;

  /* prepare cleanup */
  nodes = tets = attr = NULL;
  ptg = P4EST_ALLOC (p8est_tetgen_t, 1);

  /* read nodes */
  snprintf (nodefilename, BUFSIZ, "%s.node", tetgenbasename);
  nodes = ptg->nodes = p8est_tetgen_read_node (nodefilename);
  if (nodes == NULL) {
    P4EST_LERRORF ("Failed to read nodes for %s\n", tetgenbasename);
    goto dead;
  }
  num_nodes = (int) nodes->elem_count;

  /* read tetrahedra */
  snprintf (elefilename, BUFSIZ, "%s.ele", tetgenbasename);
  tets = ptg->tets = p8est_tetgen_read_ele (elefilename, num_nodes, &attr);
  if (tets == NULL) {
    P4EST_ASSERT (attr == NULL);
    P4EST_LERRORF ("Failed to read tetrahedra for %s\n", tetgenbasename);
    goto dead;
  }
  ptg->tet_attributes = attr;

  /* we are successful */
  return ptg;

dead:
  /* clean up on error */
  if (nodes != NULL) {
    sc_array_destroy (nodes);
  }
  if (tets != NULL) {
    sc_array_destroy (tets);
  }
  if (attr != NULL) {
    sc_array_destroy (attr);
  }
  P4EST_FREE (ptg);
  return NULL;
}

void
p8est_tetgen_destroy (p8est_tetgen_t * ptg)
{
  sc_array_destroy (ptg->nodes);
  sc_array_destroy (ptg->tets);
  if (ptg->tet_attributes != NULL) {
    sc_array_destroy (ptg->tet_attributes);
  }

  P4EST_FREE (ptg);
}

static int
p8est_tet_is_righthanded (sc_array_t * nodes, int *tet)
{
  int                 i, j;
  double             *nc[4];
  double              v0[3], v1[3], v2[3], cross01[3];
  double              vol;

  /* compute tet volume */
  for (i = 0; i < 4; ++i) {
    nc[i] = (double *) sc_array_index_int (nodes, 3 * tet[i]);
  }
  for (j = 0; j < 3; ++j) {
    v0[j] = nc[1][j] - nc[0][j];
    v1[j] = nc[2][j] - nc[0][j];
    v2[j] = nc[3][j] - nc[0][j];
  }
  cross01[0] = v0[1] * v1[2] - v0[2] * v1[1];
  cross01[1] = v0[2] * v1[0] - v0[0] * v1[2];
  cross01[2] = v0[0] * v1[1] - v0[1] * v1[0];
  vol = 0.;
  for (j = 0; j < 3; ++j) {
    vol += cross01[j] * v2[j];
  }
  vol *= 1. / 3.;

  return vol >= 0.;
}

static void
p8est_tet_make_righthanded (int *tet)
{
  int                 tempn;

  tempn = tet[3];
  tet[3] = tet[2];
  tet[2] = tempn;
}

int
p8est_tetgen_make_righthanded (p8est_tetgen_t * ptg)
{
  int                 i, num_tets;
  int                 num_flips;
  int                *tet;

  num_flips = 0;
  num_tets = (int) ptg->tets->elem_count / 4;
  for (i = 0; i < num_tets; ++i) {
    tet = (int *) sc_array_index_int (ptg->tets, 4 * i);
    if (!p8est_tet_is_righthanded (ptg->nodes, tet)) {
      p8est_tet_make_righthanded (tet);
      P4EST_ASSERT (p8est_tet_is_righthanded (ptg->nodes, tet));
      ++num_flips;
    }
  }

  return num_flips;
}
