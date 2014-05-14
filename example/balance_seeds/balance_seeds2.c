/*
  This file is part of p4est.
  p4est is a C library to manage a collection (a forest) of multiple
  connected adaptive quadtrees or octrees in parallel.

  Copyright (C) 2011 The University of Texas System
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

#include <sc_dmatrix.h>
#ifndef P4_TO_P8
#include <p4est_balance.h>
#include <p4est_bits.h>
#include <p4est_vtk.h>
#include <p4est_extended.h>
#else
#include <p8est_balance.h>
#include <p8est_bits.h>
#include <p8est_vtk.h>
#include <p8est_extended.h>
#endif

typedef struct
{
  int                 flag;
}
balance_seeds_elem_t;

#ifndef P4_TO_P8
static const int    refine_level = 8;
static p4est_quadrant_t center = { 0x10000000, 0x10000000, 2, 0, 0, {NULL} };
#else
static const int    refine_level = 8;
static p4est_quadrant_t center =
  { 0x20000, 0x20000, 0x20000, 2, 0, 0, {NULL} };
#endif

static void
init_fn (p4est_t * p4est, p4est_topidx_t which_tree,
         p4est_quadrant_t * quadrant)
{
  ((balance_seeds_elem_t *) (quadrant->p.user_data))->flag = -2;
}

static int
refine_fn (p4est_t * p4est, p4est_topidx_t which_tree,
           p4est_quadrant_t * quadrant)
{
#ifndef P4_TO_P8
  p4est_connect_type_t balance = P4EST_CONNECT_FACE;
#else
  p4est_connect_type_t balance = P8EST_CONNECT_EDGE;
#endif
  p4est_quadrant_t    desc;
  int                 i;

  if (((balance_seeds_elem_t *) (quadrant->p.user_data))->flag > -2) {
    return 0;
  }

  if (p4est_quadrant_is_ancestor (quadrant, &center)) {
    return 1;
  }
  if (p4est_quadrant_is_equal (quadrant, &center)) {
    ((balance_seeds_elem_t *) (quadrant->p.user_data))->flag = center.level;
    return 0;
  }
#ifndef P4_TO_P8
  if (quadrant->x >= 0x30000000 || quadrant->y >= 0x30000000) {
    ((balance_seeds_elem_t *) (quadrant->p.user_data))->flag = -1;
    return 0;
  }
#else
  if (quadrant->x >= 0x60000 || quadrant->y >= 0x60000 ||
      quadrant->z >= 0x60000) {
    ((balance_seeds_elem_t *) (quadrant->p.user_data))->flag = -1;
    return 0;
  }
#endif

  for (i = 0; i < P4EST_CHILDREN; i++) {
    p4est_quadrant_corner_descendant (quadrant, &desc, i, P4EST_QMAXLEVEL);
    if (p4est_balance_seeds (&desc, &center, balance, NULL)) {
      break;
    }
  }
  if (i == P4EST_CHILDREN) {
    P4EST_ASSERT (!p4est_balance_seeds (quadrant, &center, balance, NULL));
    ((balance_seeds_elem_t *) (quadrant->p.user_data))->flag = -1;
    return 0;
  }
  p4est_quadrant_corner_descendant (quadrant, &desc, i, quadrant->level + 1);
  if (!p4est_balance_seeds (&desc, &center, balance, NULL)) {
    if (quadrant->level < refine_level) {
      return 1;
    }
  }
  ((balance_seeds_elem_t *) (quadrant->p.user_data))->flag = quadrant->level;
  return 0;
}

int
main (int argc, char **argv)
{
  sc_MPI_Comm         mpicomm;
  int                 mpiret;
  int                 mpisize, mpirank;
  p4est_t            *p4est;
  p4est_connectivity_t *connectivity;
  sc_dmatrix_t       *vtkvec;
  p4est_tree_t       *tree;
  sc_array_t         *quadrants;
  size_t              zz, count;
  p4est_quadrant_t   *q;
  int                 i;
#ifndef P4_TO_P8
  const char          filename[] = "p4est_balance_face";
#else
  const char          filename[] = "p8est_balance_edge";
#endif

  /* initialize MPI */
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  mpicomm = sc_MPI_COMM_WORLD;
  mpiret = sc_MPI_Comm_size (mpicomm, &mpisize);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (mpicomm, &mpirank);
  SC_CHECK_MPI (mpiret);

  sc_init (mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

#ifndef P4_TO_P8
  connectivity = p4est_connectivity_new_unitsquare ();
#else
  connectivity = p8est_connectivity_new_unitcube ();
#endif

  p4est = p4est_new_ext (mpicomm, connectivity, 0, 2, 1,
                         sizeof (balance_seeds_elem_t), init_fn, NULL);

  p4est_refine (p4est, 1, refine_fn, init_fn);

  p4est_vtk_write_header (p4est, NULL, 1. - 2. * SC_EPS,
                          0, 0, 0, 0, "level", NULL, filename);
  vtkvec = sc_dmatrix_new (p4est->local_num_quadrants, P4EST_CHILDREN);
  tree = p4est_tree_array_index (p4est->trees, 0);
  quadrants = &(tree->quadrants);
  count = quadrants->elem_count;
  for (zz = 0; zz < count; zz++) {
    q = p4est_quadrant_array_index (quadrants, zz);
    for (i = 0; i < P4EST_CHILDREN; i++) {
      vtkvec->e[zz][i] = (double)
        ((balance_seeds_elem_t *) (q->p.user_data))->flag;
    }
  }
  p4est_vtk_write_point_scalar (p4est, NULL, filename, "level", vtkvec->e[0]);
  p4est_vtk_write_footer (p4est, filename);

  sc_dmatrix_destroy (vtkvec);
  p4est_destroy (p4est);
  p4est_connectivity_destroy (connectivity);

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
