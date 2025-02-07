/*
  This file is part of p4est.
  p4est is a C library to manage a collection (a forest) of multiple
  connected adaptive quadtrees or octrees in parallel.

  Copyright (C) 2010 The University of Texas System
  Additional copyright (C) 2011 individual authors
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

/* Generate coordinate tuples for quadrants and hash them. */

#ifndef P4_TO_P8
#include <p4est_bits.h>
#include <p4est_connectivity.h>
#else
#include <p8est_bits.h>
#include <p8est_connectivity.h>
#endif

/* A coordinate point is uniquefied by adding the lowest touching tree */
typedef struct coordinates_hash_key
{
  p4est_qcoord_t      coords[P4EST_DIM];        /* in/on the unit tree */
  p4est_topidx_t      which_tree;       /* legal wrt. connectivity */
}
coordinates_hash_key_t;

/* Provide an allocator for the key data as well as the hash map itself */
typedef struct coordinates_hash_data
{
  sc_mempool_t       *ckeys;    /* memory pool for allocating the hash keys */
  sc_hash_t          *chash;    /* the hash map links keys without copying */
  p4est_locidx_t      added;    /* count each coordinate point just once */
  p4est_locidx_t      duped;    /* count attempts to add more than once */
}
coordinates_hash_data_t;

/* Calculate a hash function for a coordinate point */
static unsigned
coordinates_hash_fn (const void *v, const void *u)
{
  const coordinates_hash_key_t *k = (coordinates_hash_key_t *) v;
  uint32_t            a, b, c;

  P4EST_ASSERT (k != NULL);

  a = (uint32_t) k->coords[0];
  b = (uint32_t) k->coords[1];
#ifndef P4_TO_P8
  c = (uint32_t) k->which_tree;
#else
  c = (uint32_t) k->coords[2];
  sc_hash_mix (a, b, c);
  a += (uint32_t) k->which_tree;
#endif
  sc_hash_final (a, b, c);

  return (unsigned) c;
}

/* Determine whether two coordinate points are equal */
static int
coordinates_equal_fn (const void *v1, const void *v2, const void *u)
{
  const coordinates_hash_key_t *k1 = (coordinates_hash_key_t *) v1;
  const coordinates_hash_key_t *k2 = (coordinates_hash_key_t *) v2;

  P4EST_ASSERT (k1 != NULL);
  P4EST_ASSERT (k2 != NULL);

  return k1->which_tree == k2->which_tree &&
    !memcmp (k1->coords, k2->coords, P4EST_DIM * sizeof (p4est_qcoord_t));
}

/* Check whether a p4est node quadrant sits at a given coordinate tuple */
static int
test_node_coordinates (const p4est_quadrant_t *r,
                       const p4est_qcoord_t coords[])
{
  P4EST_ASSERT (r != NULL);
  P4EST_ASSERT (p4est_quadrant_is_node (r, 0));
  P4EST_ASSERT (coords != NULL);
  return (r->x == coords[0] && r->y == coords[1] &&
#ifdef P4_TO_P8
          r->z == coords[2] &&
#endif
          1);
}

/* Allocate and initialize hash key and try to add coordinates to map */
static void
test_hash (coordinates_hash_data_t *hdata,
           p4est_topidx_t tt, const p4est_qcoord_t coords[])
{
  void              **found;
  coordinates_hash_key_t *k;

  k = (coordinates_hash_key_t *) sc_mempool_alloc (hdata->ckeys);
  k->which_tree = tt;
  k->coords[0] = coords[0];
  k->coords[1] = coords[1];
#ifdef P4_TO_P8
  k->coords[2] = coords[2];
#endif
  if (sc_hash_insert_unique (hdata->chash, k, &found)) {
    /* The key is newly linked into the hash table: count it */
    P4EST_ASSERT (*found == k);
#ifndef P4_TO_P8
    P4EST_INFOF ("First time adding tree %ld, coordinates %lx %lx\n",
                 (long) tt, (long) coords[0], (long) coords[1]);
#else
    P4EST_INFOF ("First time adding tree %ld, coordinates %lx %lx %lx\n",
                 (long) tt, (long) coords[0], (long) coords[1],
                 (long) coords[2]);
#endif
    ++hdata->added;
  }
  else {
    /* The key for this coordinate had already been stored earlier */
    P4EST_ASSERT (*found != k);
    sc_mempool_free (hdata->ckeys, k);
    ++hdata->duped;
  }
}

/* The main routine to test our coordinate functions */
static void
test_connectivity (sc_MPI_Comm mpicomm, p4est_connectivity_t *conn)
{
  int                 mpiret;
  int                 size, rank;
  int                 face, corner;
#ifdef P4_TO_P8
  int                 edge;
#endif
  p4est_topidx_t      num_trees, tt, nt;
  p4est_qcoord_t      coords[P4EST_DIM];
  p4est_qcoord_t      coords_out[P4EST_DIM];
  p4est_quadrant_t    root, *q, node, *r;
  coordinates_hash_data_t shdata, *hdata = &shdata;

  mpiret = sc_MPI_Comm_size (mpicomm, &size);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (mpicomm, &rank);
  SC_CHECK_MPI (mpiret);

  P4EST_ASSERT (p4est_connectivity_is_valid (conn));
  num_trees = conn->num_trees;

  if (rank != 0) {
    /* the connectivity is the same on every rank */
    return;
  }

  /* hash table for coordinate counting and checking */
  hdata->ckeys = sc_mempool_new (sizeof (coordinates_hash_key_t));
  hdata->chash = sc_hash_new (coordinates_hash_fn, coordinates_equal_fn,
                              hdata, NULL);
  hdata->added = hdata->duped = 0;

  /* initialize quadrant storage */
  p4est_quadrant_root (q = &root);
  r = &node;

  /* loop over all trees in the connectivity */
  for (tt = 0; tt < num_trees; ++tt) {
    P4EST_INFOF ("Going through tree %ld\n", (long) tt);

    /* verify volume midpoint */
    p4est_quadrant_volume_coordinates (q, coords);
    p4est_connectivity_coordinates_canonicalize
      (conn, tt, coords, &nt, coords_out);
    SC_CHECK_ABORT (nt == tt, "Mysterious volume tree");
    SC_CHECK_ABORT (!p4est_coordinates_compare (coords, coords_out),
                    "Mysterious volume coordinates");
    test_hash (hdata, nt, coords_out);

    for (face = 0; face < P4EST_FACES; ++face) {
      /* verify face midpoints */
      p4est_quadrant_face_coordinates (q, face, coords);
      p4est_connectivity_coordinates_canonicalize
        (conn, tt, coords, &nt, coords_out);
      SC_CHECK_ABORT (nt <= tt, "Mysterious face tree");
      SC_CHECK_ABORT (nt < tt ||
                      p4est_coordinates_compare (coords, coords_out) >= 0,
                      "Mysterious face coordinates");
      test_hash (hdata, nt, coords_out);
    }

#ifdef P4_TO_P8
    for (edge = 0; edge < P8EST_EDGES; ++edge) {
      /* verify edge midpoints */
      p8est_quadrant_edge_coordinates (q, edge, coords);
      p4est_connectivity_coordinates_canonicalize
        (conn, tt, coords, &nt, coords_out);
      SC_CHECK_ABORT (nt <= tt, "Mysterious edge tree");
      SC_CHECK_ABORT (nt < tt ||
                      p4est_coordinates_compare (coords, coords_out) >= 0,
                      "Mysterious edge coordinates");
      test_hash (hdata, nt, coords_out);
    }
#endif

    for (corner = 0; corner < P4EST_CHILDREN; ++corner) {
      /* verify tree corners */
      p4est_quadrant_corner_node (q, corner, r);
      p4est_quadrant_corner_coordinates (q, corner, coords);
      SC_CHECK_ABORT (test_node_coordinates (r, coords), "Node coordinates");
      p4est_connectivity_coordinates_canonicalize
        (conn, tt, coords, &nt, coords_out);
      SC_CHECK_ABORT (nt <= tt, "Mysterious corner tree");
      SC_CHECK_ABORT (nt < tt ||
                      p4est_coordinates_compare (coords, coords_out) >= 0,
                      "Mysterious corner coordinates");
      test_hash (hdata, nt, coords_out);
    }
  }

  /* clean up memory */
  sc_hash_destroy (hdata->chash);
  sc_mempool_destroy (hdata->ckeys);
  P4EST_PRODUCTIONF ("Added %ld coordinates, duplicates %ld\n",
                     (long) hdata->added, (long) hdata->duped);
}

/* wrapper around test routine */
static void
test_coordinates (sc_MPI_Comm mpicomm)
{
  /* a connectivity stores connections between the trees in the mesh */
  p4est_connectivity_t *conn;

#ifndef P4_TO_P8
  /* this connectivity has three 2D trees all connected at one corner */
  conn = p4est_connectivity_new_corner ();
#else
  /* this connectivity has six 3D trees connected in some blocky way */
  conn = p8est_connectivity_new_rotcubes ();
#endif

  /* run the actual test routine */
  test_connectivity (mpicomm, conn);

  /* the connectivity is no longer used */
  p4est_connectivity_destroy (conn);
}

/* boilerplate code and calling the test routine */
int
main (int argc, char **argv)
{
  sc_MPI_Comm         mpicomm;
  int                 mpiret;

  /* initialize MPI environment */
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  mpicomm = sc_MPI_COMM_WORLD;

  /* establish parallel logging */
  sc_init (mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  /* proceed with tests */
  test_coordinates (mpicomm);

  /* clean up and exit */
  sc_finalize ();

  /* finalize MPI environment */
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  /* by default exit cleanly */
  return EXIT_SUCCESS;
}
