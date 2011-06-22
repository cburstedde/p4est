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
  int                 i;
  int                 num_nodes, dims, num_attributes, boundary_marker;
  sc_array_t         *nodes;
  FILE               *nodefile;

  /* open node file */
  nodefile = fopen (nodefilename, "rb");
  if (nodefile == NULL) {
    P4EST_LERRORF ("Failed to open %s\n", nodefilename);
    return NULL;
  }

  /* read header information */
  retval = fscanf (nodefile, "%d %d %d %d",
                   &num_nodes, &dims, &num_attributes, &boundary_marker);
  if (retval != 4 || num_nodes < 0 || dims != 3 || num_attributes < 0) {
    P4EST_LERROR ("Failed to read node header\n");
    fclose (nodefile);
    return NULL;
  }

  /* read node coordinates */
  nodes = sc_array_new_size (sizeof (double), 3 * num_nodes);
  for (i = 0; i < num_nodes; ++i) {

  }

  /* close node file */
  retval = fclose (nodefile);
  if (retval) {
    P4EST_LERRORF ("Failed to close %s\n", nodefilename);
    sc_array_destroy (nodes);
    return NULL;
  }

  return nodes;
}

sc_array_t         *
p8est_tetgen_read_ele (const char *elefilename, sc_array_t ** attributes)
{
  int                 retval;
  int                 i;
  int                 num_tets, nodespertet, region;
  sc_array_t         *tets, *attr;
  FILE               *elefile;

  attr = NULL;
  if (attributes != NULL) {
    *attributes = NULL;
  }

  /* open ele file */
  elefile = fopen (elefilename, "rb");
  if (elefile == NULL) {
    P4EST_LERRORF ("Failed to open %s\n", elefilename);
    return NULL;
  }

  /* read header information */
  retval = fscanf (elefile, "%d %d %d", &num_tets, &nodespertet, &region);
  if (retval != 3 || num_tets < 0 || nodespertet != 4) {
    P4EST_LERROR ("Failed to read tet header\n");
    fclose (elefile);
    return NULL;
  }

  /* read tet coordinates */
  tets = sc_array_new_size (sizeof (int), 4 * num_tets);
  if (region && attributes != NULL) {
    attr = *attributes = sc_array_new_size (sizeof (int), num_tets);
  }
  for (i = 0; i < num_tets; ++i) {

  }

  /* close ele file */
  retval = fclose (elefile);
  if (retval) {
    P4EST_LERRORF ("Failed to close %s\n", elefilename);
    sc_array_destroy (tets);
    return NULL;
  }

  return tets;
}

p8est_tetgen_t     *
p8est_tetgen_read (const char *tetgenbasename)
{
  char                nodefilename[BUFSIZ];
  char                elefilename[BUFSIZ];
  p8est_tetgen_t     *ptg;

  ptg = P4EST_ALLOC_ZERO (p8est_tetgen_t, 1);

  /* read nodes */
  snprintf (nodefilename, BUFSIZ, "%s.node", tetgenbasename);
  ptg->nodes = p8est_tetgen_read_node (nodefilename);
  if (ptg->nodes == NULL) {
    P4EST_LERRORF ("Failed to read nodes for %s\n", tetgenbasename);
    return NULL;
  }

  /* read tetrahedra */
  snprintf (elefilename, BUFSIZ, "%s.ele", tetgenbasename);
  ptg->tets = p8est_tetgen_read_ele (elefilename, &ptg->tet_attributes);
  if (ptg->tets == NULL) {
    P4EST_ASSERT (ptg->tet_attributes == NULL);
    P4EST_LERRORF ("Failed to read tetrahedra for %s\n", tetgenbasename);
    return NULL;
  }

  return ptg;
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
