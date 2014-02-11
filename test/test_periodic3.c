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

#include <p8est.h>

static void
test_periodic (p8est_connectivity_t * conn)
{
  int                 i;
  int                 iface, nface;
  int                 iedge, icorner;
  int                 ft[9];
  size_t              zz;
  p4est_topidx_t      itree, ntree;
  p8est_edge_info_t   ei;
  p8est_edge_transform_t *et;
  p8est_corner_info_t ci;
  p8est_corner_transform_t *ct;
  sc_array_t         *eta, *cta;

  itree = 0;
  for (iface = 0; iface < 6; ++iface) {
    ntree = p8est_find_face_transform (conn, itree, iface, ft);
    nface = p8est_face_dual[iface];
    SC_CHECK_ABORT (ntree == itree, "PF tree");
    for (i = 0; i < 3; ++i) {
      SC_CHECK_ABORTF (ft[i] == ft[i + 3], "PF axis %d", i);
    }
    SC_CHECK_ABORT (ft[6] == 0 && ft[7] == 0, "PF reverse");
    SC_CHECK_ABORT (ft[8] == 2 * (iface % 2) + nface % 2, "PF code");
  }

  eta = &ei.edge_transforms;
  sc_array_init (eta, sizeof (p8est_edge_transform_t));
  for (iedge = 0; iedge < 12; ++iedge) {
    p8est_find_edge_transform (conn, itree, iedge, &ei);
    SC_CHECK_ABORT ((int) ei.iedge == iedge, "PE ei");
    SC_CHECK_ABORT (eta->elem_count == 1, "PE count");
    for (zz = 0; zz < eta->elem_count; ++zz) {
      et = p8est_edge_array_index (eta, zz);
      SC_CHECK_ABORT (et->ntree == itree, "PE tree");
      SC_CHECK_ABORT ((int) et->nedge + iedge == 8 * (iedge / 4) + 3,
                      "PE edge");
      SC_CHECK_ABORT (et->nflip == 0, "PE flip");
      SC_CHECK_ABORT (et->corners == et->nedge % 4, "PE corners");
      SC_CHECK_ABORT ((int) et->naxis[0] == iedge / 4 &&
                      et->naxis[1] < et->naxis[2] &&
                      et->naxis[0] + et->naxis[1] + et->naxis[2] == 3,
                      "PE axis");
    }
  }
  sc_array_reset (eta);

  cta = &ci.corner_transforms;
  sc_array_init (cta, sizeof (p8est_corner_transform_t));
  for (icorner = 0; icorner < 8; ++icorner) {
    p8est_find_corner_transform (conn, itree, icorner, &ci);
    SC_CHECK_ABORT ((int) ci.icorner == icorner, "PC ci");
    SC_CHECK_ABORT (cta->elem_count == 1, "PC count");
    for (zz = 0; zz < cta->elem_count; ++zz) {
      ct = p8est_corner_array_index (cta, zz);
      SC_CHECK_ABORT (ct->ntree == itree, "PC tree");
      SC_CHECK_ABORT (ct->ncorner + icorner == 7, "PC corner");
    }
  }
  sc_array_reset (cta);
}

/* *INDENT-OFF* */
static const int
rotwrap_edges[8][2] =
{{ 1,  6 }, { 0,  7 }, { 3,  5 }, { 2,  4 },
 { 3, -1 }, { 2, -1 }, { 0, -1 }, { 1, -1 }};

static const int
rotwrap_axes[8][2] =
{{ 0,  1 }, { 0,  1 }, { 0,  1 }, { 0,  1 },
 { 0, -1 }, { 0, -1 }, { 0, -1 }, { 0, -1 }};

static const int
rotwrap_flip[8][2] =
{{ 0,  0 }, { 0,  0 }, { 0,  1 }, { 0,  1 },
 { 1, -1 }, { 1, -1 }, { 0, -1 }, { 0, -1 }};

static const int
rotwrap_corners[8][2] =
{{ 3, 6 }, { 2, 4 }, { 1, 7 }, { 0, 5 },
 { 1, 7 }, { 3, 6 }, { 0, 5 }, { 2, 4 }};
/* *INDENT-ON* */

static void
test_rotwrap (p8est_connectivity_t * conn)
{
  int                 i;
  int                 iface, nface;
  int                 iedge, icorner;
  int                 ft[9];
  size_t              zz;
  p4est_topidx_t      itree, ntree;
  p8est_edge_info_t   ei;
  p8est_edge_transform_t *et;
  p8est_corner_info_t ci;
  p8est_corner_transform_t *ct;
  sc_array_t         *eta, *cta;

  itree = 0;
  for (iface = 0; iface < 6; ++iface) {
    ntree = p8est_find_face_transform (conn, itree, iface, ft);
    if (iface == 2 || iface == 3) {
      SC_CHECK_ABORT (ntree == -1, "RF tree");
      continue;
    }
    nface = p8est_face_dual[iface];
    SC_CHECK_ABORT (ntree == itree, "RF tree");
    if (iface == 0 || iface == 1) {
      for (i = 0; i < 3; ++i) {
        SC_CHECK_ABORTF (ft[i] == (i + 1) % 3, "RFA axis A%d", i);
        SC_CHECK_ABORTF (ft[i] == ft[i + 3], "RFA axis B%d", i);
      }
      SC_CHECK_ABORT (ft[6] == 0 && ft[7] == 0, "RFA reverse");
    }
    else {
      for (i = 0; i < 3; ++i) {
        SC_CHECK_ABORTF (ft[i] == i, "RFB axis A%d", i);
      }
      SC_CHECK_ABORT (ft[0] == ft[4], "RFB axis B0");
      SC_CHECK_ABORT (ft[1] == ft[3], "RFB axis B1");
      SC_CHECK_ABORT (ft[2] == ft[5], "RFB axis B2");
      SC_CHECK_ABORT (ft[6] == (iface != 4) &&
                      ft[7] == (iface != 5), "RFB reverse");
    }
    SC_CHECK_ABORT (ft[8] == 2 * (iface % 2) + nface % 2, "RF code");
  }

  eta = &ei.edge_transforms;
  sc_array_init (eta, sizeof (p8est_edge_transform_t));
  for (iedge = 0; iedge < 12; ++iedge) {
    p8est_find_edge_transform (conn, itree, iedge, &ei);
    SC_CHECK_ABORT ((int) ei.iedge == iedge, "RE ei");
    SC_CHECK_ABORT ((int) eta->elem_count == 2 - (iedge / 4), "RE count AB");
    for (zz = 0; zz < eta->elem_count; ++zz) {
      et = p8est_edge_array_index (eta, zz);
      SC_CHECK_ABORT (et->ntree == itree, "RE tree");
      SC_CHECK_ABORT ((int) et->nedge == rotwrap_edges[iedge][zz], "RE edge");
      SC_CHECK_ABORT ((int) et->nflip == rotwrap_flip[iedge][zz], "RE flip");
      SC_CHECK_ABORT (et->corners == et->nedge % 4, "RE corners");
      SC_CHECK_ABORT ((int) et->naxis[0] == rotwrap_axes[iedge][zz] &&
                      et->naxis[1] < et->naxis[2] &&
                      et->naxis[0] + et->naxis[1] + et->naxis[2] == 3,
                      "RE axis");
    }
  }
  sc_array_reset (eta);

  cta = &ci.corner_transforms;
  sc_array_init (cta, sizeof (p8est_corner_transform_t));
  for (icorner = 0; icorner < 8; ++icorner) {
    p8est_find_corner_transform (conn, itree, icorner, &ci);
    SC_CHECK_ABORT ((int) ci.icorner == icorner, "RC ci");
    SC_CHECK_ABORT (cta->elem_count == 2, "RC count");
    for (zz = 0; zz < cta->elem_count; ++zz) {
      ct = p8est_corner_array_index (cta, zz);
      SC_CHECK_ABORT (ct->ntree == itree, "RC tree");
      SC_CHECK_ABORT ((int) ct->ncorner == rotwrap_corners[icorner][zz],
                      "RC corner");
    }
  }
  sc_array_reset (cta);
}

/* *INDENT-OFF* */
static const int
weird_edges[2][2] = {{ 5, 7 }, { 7, 5 }};
/* *INDENT-ON* */

static void
test_weird (void)
{
  const p4est_topidx_t num_edges = 1, num_ett = 2;
  const p4est_topidx_t num_corners = 1, num_ctt = 4;
  int                 i;
  size_t              zz;
  p8est_edge_info_t   ei;
  p8est_edge_transform_t *et;
  p8est_corner_info_t ci;
  p8est_corner_transform_t *ct;
  sc_array_t         *eta, *cta;
  p8est_connectivity_t *conn;

  conn = p8est_connectivity_new (0, 1,
                                 num_edges, num_ett, num_corners, num_ctt);
  for (i = 0; i < 6; ++i) {
    conn->tree_to_tree[i] = 0;
    conn->tree_to_face[i] = (int8_t) i;
  }
  conn->tree_to_face[4] = 5;
  conn->tree_to_face[5] = 4;

  for (i = 0; i < 12; ++i) {
    conn->tree_to_edge[i] = -1;
  }
  conn->tree_to_edge[5] = 0;
  conn->tree_to_edge[7] = 0;
  conn->edge_to_tree[0] = 0;
  conn->edge_to_tree[1] = 0;
  conn->edge_to_edge[0] = 5;
  conn->edge_to_edge[1] = 19;
  conn->ett_offset[0] = 0;

  for (i = 0; i < 8; ++i) {
    conn->tree_to_corner[i] = -1;
  }
  conn->tree_to_corner[0] = 0;
  conn->tree_to_corner[1] = 0;
  conn->tree_to_corner[4] = 0;
  conn->tree_to_corner[5] = 0;
  conn->corner_to_tree[0] = 0;
  conn->corner_to_tree[1] = 0;
  conn->corner_to_tree[2] = 0;
  conn->corner_to_tree[3] = 0;
  conn->corner_to_corner[0] = 0;
  conn->corner_to_corner[1] = 1;
  conn->corner_to_corner[2] = 4;
  conn->corner_to_corner[3] = 5;
  conn->ctt_offset[0] = 0;

  P4EST_ASSERT (p8est_connectivity_is_valid (conn));

  eta = &ei.edge_transforms;
  sc_array_init (eta, sizeof (p8est_edge_transform_t));
  for (i = 0; i < 2; ++i) {
    p8est_find_edge_transform (conn, 0, weird_edges[i][0], &ei);
    SC_CHECK_ABORT ((int) ei.iedge == weird_edges[i][0], "WE ei");
    SC_CHECK_ABORT (eta->elem_count == 1, "WE count A");
    for (zz = 0; zz < eta->elem_count; ++zz) {
      et = p8est_edge_array_index (eta, zz);
      SC_CHECK_ABORT (et->ntree == 0, "WE tree");
      SC_CHECK_ABORT ((int) et->nedge == weird_edges[i][1], "WE edge");
      SC_CHECK_ABORT (et->nflip == 1, "WE flip");
      SC_CHECK_ABORT (et->corners == et->nedge % 4, "WE corners");
      SC_CHECK_ABORT (et->naxis[0] == 1 && et->naxis[1] == 0 &&
                      et->naxis[2] == 2, "WE axis");
    }
  }
  sc_array_reset (eta);

  cta = &ci.corner_transforms;
  sc_array_init (cta, sizeof (p8est_corner_transform_t));
  for (i = 0; i < 8; ++i) {
    p8est_find_corner_transform (conn, 0, i, &ci);
    SC_CHECK_ABORT ((int) ci.icorner == i, "WC ci");
    SC_CHECK_ABORT ((int) cta->elem_count == 2 - (i & 0x02), "WC count");
    for (zz = 0; zz < cta->elem_count; ++zz) {
      ct = p8est_corner_array_index (cta, zz);
      SC_CHECK_ABORT (ct->ntree == 0, "WC tree");
      SC_CHECK_ABORT ((size_t) ct->ncorner == 4 * zz + !(i % 2), "WC corner");
    }
  }
  sc_array_reset (cta);

  p8est_connectivity_destroy (conn);
}

/*
 * Purpose of this program is to verify that
 * p8est_find_edge_transform and p8est_find_corner_transform
 * work as expected for several periodic connectivities.
 */
int
main (int argc, char **argv)
{
  p8est_connectivity_t *conn;

  sc_init (sc_MPI_COMM_NULL, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  conn = p8est_connectivity_new_periodic ();
  test_periodic (conn);
  p8est_connectivity_destroy (conn);

  conn = p8est_connectivity_new_rotwrap ();
  test_rotwrap (conn);
  p8est_connectivity_destroy (conn);

  test_weird ();

  sc_finalize ();

  return 0;
}
