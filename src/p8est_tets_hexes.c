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
p8est_tets_read_node (const char *nodefilename)
{
  int                 retval;
  int                 k;
  int                 dims, num_attributes, boundary_marker;
  long                jl, lnum_nodes;
  double             *pc;
  size_t              iz, znum_nodes;
  FILE               *nodefile;
  sc_array_t         *nodes;

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
  retval = fscanf (nodefile, "%ld %d %d %d",
                   &lnum_nodes, &dims, &num_attributes, &boundary_marker) - 4;
  if (retval || lnum_nodes < 0 || lnum_nodes > P4EST_TOPIDX_MAX
      || dims != 3 || num_attributes < 0) {
    P4EST_LERROR ("Failed to read node header\n");
    goto dead;
  }
  znum_nodes = (size_t) lnum_nodes;

  /* read node coordinates */
  nodes = sc_array_new_size (sizeof (double), 3 * znum_nodes);
  for (iz = 0; iz < znum_nodes; ++iz) {
    pc = (double *) sc_array_index (nodes, 3 * iz);
    retval = fscanf (nodefile, "%ld %lf %lf %lf",
                     &jl, pc, pc + 1, pc + 2) - 4;
    if (retval || (long) iz != jl) {
      P4EST_LERRORF ("Failed to read node %ld coordinates\n", (long) iz);
      goto dead;
    }
    for (k = 0; k < num_attributes; ++k) {
      retval = fscanf (nodefile, "%*f");
      if (retval != 0) {
        P4EST_LERRORF ("Failed to read node %ld attribute %d\n", (long) iz,
                       k);
        goto dead;
      }
    }
    if (boundary_marker) {
      retval = fscanf (nodefile, "%*d");
      if (retval != 0) {
        P4EST_LERRORF ("Failed to read node %ld boundary marker\n",
                       (long) iz);
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
p8est_tets_read_ele (const char *elefilename, p4est_topidx_t num_nodes,
                     sc_array_t ** attributes)
{
  int                 retval;
  int                 k;
  int                 nodespertet, region;
  long                jl, lnum_tets, lmax_nodes;
  long                nl[4];
  size_t              iz, znum_tets;
  int                *pi;
  p4est_topidx_t     *pt;
  FILE               *elefile;
  sc_array_t         *tets, *attr;

  /* prepare cleanup on error */
  elefile = NULL;
  tets = NULL;
  attr = NULL;
  if (attributes != NULL) {
    *attributes = NULL;
  }
  lmax_nodes = (num_nodes >= 0 ? (long) num_nodes : P4EST_TOPIDX_MAX);

  /* open ele file */
  elefile = fopen (elefilename, "rb");
  if (elefile == NULL) {
    P4EST_LERRORF ("Failed to open %s\n", elefilename);
    goto dead;
  }

  /* read header information */
  retval = fscanf (elefile, "%ld %d %d",
                   &lnum_tets, &nodespertet, &region) - 3;
  if (retval || lnum_tets < 0 || lnum_tets > P4EST_TOPIDX_MAX
      || nodespertet != 4) {
    P4EST_LERROR ("Failed to read tet header\n");
    goto dead;
  }
  znum_tets = (size_t) lnum_tets;

  /* read tet coordinates */
  tets = sc_array_new_size (sizeof (p4est_topidx_t), 4 * znum_tets);
  if (region && attributes != NULL) {
    attr = *attributes = sc_array_new_size (sizeof (int), znum_tets);
  }
  for (iz = 0; iz < znum_tets; ++iz) {
    pt = (p4est_topidx_t *) sc_array_index (tets, 4 * iz);
    retval = fscanf (elefile, "%ld %ld %ld %ld %ld",
                     &jl, nl, nl + 1, nl + 2, nl + 3) - 5;
    if (retval || (long) iz != jl) {
      P4EST_LERRORF ("Failed to read tet %ld node numbers\n", (long) iz);
      goto dead;
    }
    for (k = 0; k < 4; ++k) {
      if (nl[k] < 0 || nl[k] > lmax_nodes) {
        P4EST_LERRORF ("Tet %ld has invalid node number %d\n", (long) iz, k);
        goto dead;
      }
      pt[k] = (p4est_topidx_t) nl[k];
    }
    if (region) {
      if (attr != NULL) {
        pi = (int *) sc_array_index (attr, iz);
        retval = fscanf (elefile, "%d", pi) - 1;
      }
      else {
        retval = fscanf (elefile, "%*d");
      }
      if (retval != 0) {
        P4EST_LERRORF ("Failed to read tet %ld region attribute\n",
                       (long) iz);
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

p8est_tets_t       *
p8est_tets_read (const char *tetgenbasename)
{
  p4est_topidx_t      num_nodes;
  char                nodefilename[BUFSIZ];
  char                elefilename[BUFSIZ];
  sc_array_t         *nodes, *tets, *attr;
  p8est_tets_t       *ptg;

  /* prepare cleanup */
  nodes = tets = attr = NULL;
  ptg = P4EST_ALLOC (p8est_tets_t, 1);

  /* read nodes */
  snprintf (nodefilename, BUFSIZ, "%s.node", tetgenbasename);
  nodes = ptg->nodes = p8est_tets_read_node (nodefilename);
  if (nodes == NULL) {
    P4EST_LERRORF ("Failed to read nodes for %s\n", tetgenbasename);
    goto dead;
  }
  num_nodes = (p4est_topidx_t) (nodes->elem_count / 3);

  /* read tetrahedra */
  snprintf (elefilename, BUFSIZ, "%s.ele", tetgenbasename);
  tets = ptg->tets = p8est_tets_read_ele (elefilename, num_nodes, &attr);
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
p8est_tets_destroy (p8est_tets_t * ptg)
{
  sc_array_destroy (ptg->nodes);
  sc_array_destroy (ptg->tets);
  if (ptg->tet_attributes != NULL) {
    sc_array_destroy (ptg->tet_attributes);
  }

  P4EST_FREE (ptg);
}

static int
p8est_tet_is_righthanded (sc_array_t * nodes, p4est_topidx_t * tet)
{
  int                 i, j;
  double             *nc[4];
  double              v0[3], v1[3], v2[3], cross01[3];
  double              vol;

  /* compute tet volume */
  for (i = 0; i < 4; ++i) {
    nc[i] = (double *) sc_array_index (nodes, (size_t) (3 * tet[i]));
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
p8est_tet_flip (p4est_topidx_t * tet)
{
  p4est_topidx_t      temp;

  temp = tet[3];
  tet[3] = tet[2];
  tet[2] = temp;
}

p4est_topidx_t
p8est_tets_make_righthanded (p8est_tets_t * ptg)
{
  size_t              iz, znum_tets;
  p4est_topidx_t      tnum_flips;
  p4est_topidx_t     *tet;

  tnum_flips = 0;
  znum_tets = ptg->tets->elem_count / 4;
  for (iz = 0; iz < znum_tets; ++iz) {
    tet = (p4est_topidx_t *) sc_array_index (ptg->tets, 4 * iz);
    if (!p8est_tet_is_righthanded (ptg->nodes, tet)) {
      p8est_tet_flip (tet);
      P4EST_ASSERT (p8est_tet_is_righthanded (ptg->nodes, tet));
      ++tnum_flips;
    }
  }

  return tnum_flips;
}
