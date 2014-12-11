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

#ifndef P4_TO_P8
#include <p4est_lnodes.h>
#include <p4est_bits.h>
#else
#include <p8est_lnodes.h>
#include <p8est_bits.h>
#endif

static void
trilinear_interp (double (*v)[3], double eta[3], double xyz[3])
{
  int                 i;

  for (i = 0; i < 3; i++) {
    xyz[i] = (1. - eta[2]) * ((1. - eta[1]) * ((1. - eta[0]) * v[0][i]
                                               + eta[0] * v[1][i]
                              )
                              + eta[1] * ((1. - eta[0]) * v[2][i]
                                          + eta[0] * v[3][i]
                              )
      )
#ifdef P4_TO_P8
      + eta[2] * ((1. - eta[1]) * ((1. - eta[0]) * v[4][i]
                                   + eta[0] * v[5][i]
                  )
                  + eta[1] * ((1. - eta[0]) * v[6][i]
                              + eta[0] * v[7][i]
                  )
      )
#endif
      ;
  }
}

p4est_connectivity_t *
p4est_connectivity_refine (p4est_connectivity_t * conn_in, int num_per_edge)
{
  p4est_t            *dummy_forest;
  p4est_ghost_t      *dummy_ghost;
  p4est_lnodes_t     *dummy_lnodes;
  p4est_connectivity_t *conn_out;
  p4est_topidx_t      num_old_trees = conn_in->num_trees;
  int                 ceillog = SC_LOG2_32 (num_per_edge - 1) + 1;
  int                 ceil = (1 << ceillog);
#ifndef P4_TO_P8
  int                 N = num_per_edge * num_per_edge;
  int                 M = ceil * ceil;
#else
  int                 N = num_per_edge * num_per_edge * num_per_edge;
  int                 M = ceil * ceil * ceil;
#endif
  p4est_topidx_t      num_new_trees = num_old_trees * N;
  p4est_topidx_t      num_new_vertices, ti, count;
  int                 j;

  P4EST_ASSERT (num_per_edge >= 1);

  /* each processor redundantly creates the new connectivity */
  dummy_forest = p4est_new (MPI_COMM_SELF, conn_in, 0, 0, NULL);
  dummy_ghost = p4est_ghost_new (dummy_forest, P4EST_CONNECT_FULL);
  dummy_lnodes = p4est_lnodes_new (dummy_forest, dummy_ghost, num_per_edge);

  num_new_vertices = (p4est_topidx_t) dummy_lnodes->num_local_nodes;

  conn_out = p4est_connectivity_new (num_new_vertices, num_new_trees,
#ifdef P4_TO_P8
                                     0, 0,
#endif
                                     0, 0);

  for (ti = 0; ti < num_new_trees; ti++) {
    for (j = 0; j < P4EST_FACES; j++) {
      conn_out->tree_to_tree[P4EST_FACES * ti + j] = ti;
      conn_out->tree_to_face[P4EST_FACES * ti + j] = j;
    }
  }
  for (count = 0, ti = 0; ti < num_old_trees; ti++) {
    double              v[P4EST_CHILDREN][3];

    for (j = 0; j < P4EST_CHILDREN; j++) {
      int                 k;

      for (k = 0; k < 3; k++) {
        v[j][k] =
          conn_in->vertices[3 *
                            conn_in->tree_to_vertex[P4EST_CHILDREN * ti + j] +
                            k];
      }
    }
    for (j = 0; j < M; j++) {
      p4est_quadrant_t    dummy;
      uint64_t            R = j;
      int                 x[P4EST_DIM], k;
      int                 id, pow;
      double              xyz[3];
      p4est_topidx_t      thisvert;

      p4est_quadrant_set_morton (&dummy, ceillog, R);

      x[0] = (dummy.x >> (P4EST_MAXLEVEL - ceillog));
      x[1] = (dummy.y >> (P4EST_MAXLEVEL - ceillog));
#ifdef P4_TO_P8
      x[2] = (dummy.z >> (P4EST_MAXLEVEL - ceillog));
#endif
      for (k = 0; k < P4EST_DIM; k++) {
        if (x[k] >= num_per_edge) {
          break;
        }
      }
      if (k < P4EST_DIM) {
        continue;
      }

      id = 0;
      pow = 1;
      for (k = 0; k < P4EST_DIM; k++) {
        id += x[k] * pow;
        pow *= (num_per_edge + 1);
      }

      for (k = 0; k < P4EST_CHILDREN; k++) {
        int                 thisid = id, l;
        double              eta[3] = { 0. };

        pow = 1;
        for (l = 0; l < P4EST_DIM; l++) {
          int                 thisx = x[l];
          int                 thisincr = (! !(k & 1 << l));

          thisid += pow * thisincr;
          pow *= (num_per_edge + 1);
          eta[l] = ((double) (thisx + thisincr)) / ((double) num_per_edge);
        }
        P4EST_ASSERT (thisid < dummy_lnodes->vnodes);
        trilinear_interp (v, eta, xyz);
        conn_out->tree_to_vertex[P4EST_CHILDREN * count + k] = thisvert =
          dummy_lnodes->element_nodes[dummy_lnodes->vnodes * ti + thisid];
        for (l = 0; l < 3; l++) {
          conn_out->vertices[3 * thisvert + l] = xyz[l];
        }
      }

      count++;
    }
  }
  P4EST_ASSERT (count == num_new_trees);

  p4est_lnodes_destroy (dummy_lnodes);
  p4est_ghost_destroy (dummy_ghost);
  p4est_destroy (dummy_forest);

  p4est_connectivity_complete (conn_out);

  return conn_out;
}
