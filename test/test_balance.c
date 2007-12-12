/*
  This file is part of p4est.
  p4est is a C library to manage a parallel collection of quadtrees and/or
  octrees.

  Copyright (C) 2007 Carsten Burstedde, Lucas Wilcox.

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

#include <p4est_algorithms.h>
#include <p4est_base.h>

int
main (int argc, char ** argv) 
{
  int8_t              l;
  p4est_t            *p4est;
  p4est_tree_t        stree, *tree = &stree;
  p4est_quadrant_t   *q;
  p4est_connectivity_t *connectivity;

  connectivity = p4est_connectivity_new_unitsquare ();
  p4est = p4est_new (MPI_COMM_NULL, stdout, connectivity, 4, NULL);

  /* build empty tree */
  tree->quadrants = p4est_array_new (sizeof (p4est_quadrant_t));
  for (l = 0; l <= P4EST_MAXLEVEL; ++l) {
    tree->quadrants_per_level[l] = 0;
  }
  tree->maxlevel = 0;

  /* insert two quadrants */
  p4est_array_resize (tree->quadrants, 2);
  q = p4est_array_index (tree->quadrants, 0);
  p4est_quadrant_set_morton (q, 0, 0);
  ++tree->quadrants_per_level[0];
  q = p4est_array_index (tree->quadrants, 1);
  p4est_quadrant_set_morton (q, 3, 11);
  ++tree->quadrants_per_level[3];
  tree->maxlevel = 3;

  /* balance the tree */   
  p4est_balance_subtree (p4est, tree, 0, NULL);

  /* clean up memory */
  p4est_array_destroy (tree->quadrants);
  p4est_destroy (p4est);
  p4est_connectivity_destroy (connectivity);

  p4est_memory_check ();

  return 0;
}

/* EOF test_balance.c */
