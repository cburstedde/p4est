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
#include <p4est_bits.h>
#include <p4est_extended.h>
#include <p4est_iterate.h>
#include <p4est_mesh.h>
#include <p4est_search.h>
#else /* P4_TO_P8 */
#include <p8est_bits.h>
#include <p8est_extended.h>
#include <p8est_iterate.h>
#include <p8est_mesh.h>
#include <p8est_search.h>
#endif /* P4_TO_P8 */

/*********************** constructor functions ***********************/

/** For a quadrant that touches a tree face with a corner inside the face,
 * get the number of the touching face.
 *
 * \param [in] q        currently considered quadrant
 * \param [in] corner   corner index
 * \return              face index
 */
static int
tree_face_quadrant_corner_face (const p4est_quadrant_t * q, int corner)
{
  int                 which;
  p4est_qcoord_t      end = P4EST_LAST_OFFSET (q->level);

  P4EST_ASSERT (p4est_quadrant_is_valid (q));
  P4EST_ASSERT (0 <= corner && corner < P4EST_CHILDREN);

  which = corner & 1;
  if (q->x == (which ? end : 0)) {
    return which;
  }
  which = corner & 2;
  if (q->y == (which ? end : 0)) {
    return 2 + (which >> 1);
  }
#ifdef P4_TO_P8
  which = corner & 4;
  if (q->z == (which ? end : 0)) {
    return 4 + (which >> 2);
  }
#endif /* P4_TO_P8 */
  SC_ABORT_NOT_REACHED ();
}

/** Populate mesh information for corners across tree boundaries, i.e. every
 *  neighborhood scenario where we need more information (like orientation) than
 *  a single index. Note that this function only allocates the memory, correct
 *  values have to be set separately
 *
 * \param [in][out] mesh     The mesh structure to which we will add corner
 *                           information
 * \param [in]      clen     Number of quadrants to be added
 * \param [in]      pcquad   List of quadrant indices
 * \param [in]      pccorner List of quadrant encodings
 */
static              p4est_locidx_t
mesh_corner_allocate (p4est_mesh_t * mesh, p4est_locidx_t clen,
                      p4est_locidx_t ** pcquad, int8_t ** pccorner)
{
  p4est_locidx_t      cornerid, cstart, cend;

  P4EST_ASSERT (clen > 0);
  P4EST_ASSERT (mesh->corner_offset->elem_count ==
                (size_t) (mesh->local_num_corners + 1));

  cornerid = mesh->local_num_corners++;
  cstart = *(p4est_locidx_t *) sc_array_index (mesh->corner_offset, cornerid);
  cend = cstart + clen;
  *(p4est_locidx_t *) sc_array_push (mesh->corner_offset) = cend;

  P4EST_ASSERT (mesh->corner_offset->elem_count ==
                (size_t) (mesh->local_num_corners + 1));

  P4EST_ASSERT (mesh->corner_quad->elem_count == (size_t) cstart);
  *pcquad = (p4est_locidx_t *) sc_array_push_count (mesh->corner_quad, clen);
  P4EST_ASSERT (mesh->corner_quad->elem_count == (size_t) cend);

  P4EST_ASSERT (mesh->corner_corner->elem_count == (size_t) cstart);
  *pccorner = (int8_t *) sc_array_push_count (mesh->corner_corner, clen);
  P4EST_ASSERT (mesh->corner_corner->elem_count == (size_t) cend);

  return cornerid;
}

#ifdef P4_TO_P8
/** Find face neighbors of an edge neighborship
 *
 * \param side    Contains the quadrant to search face neighbors for
 * \param conn    p4est connectivity information encoding the tree relations
 * \param nftree  Result vector containing the tree index of the face neighbors
 * \param nface   Result vector containing the face index of the face neighbors
 * \param nedgef  Result vector containing the edge index of the face neighbors
 */
static int
mesh_edge_find_face_neighbors (p8est_iter_edge_side_t * side,
                               p4est_connectivity_t * conn,
                               p4est_locidx_t * nftree,
                               p4est_locidx_t * nface,
                               p4est_locidx_t * nedgef)
{
  int                 i;
  p4est_locidx_t      t1 = side->treeid;
  int                 e1 = (int) side->edge;
  p4est_locidx_t      f1, faceOrientation;

  /* Get all local quadrant faces touching this edge */
  for (i = 0; i < 2; ++i) {
    f1 = p8est_edge_faces[e1][i];
    nftree[i] = conn->tree_to_tree[P4EST_FACES * t1 + f1];
    nface[i] = conn->tree_to_face[P4EST_FACES * t1 + f1];

    if (nftree[i] == t1 && nface[i] == f1) {
      /* If the quadrant sees itself, we are at a physical face
       * boundary, i.e. there is no face neighbor */
      nedgef[i] = -1;
    }
    else {
      /* calculate orientation and index of the adjacent face and derive
       * the currently processed edge's edge index w.r.t. the adjacent
       * quadrant */
      faceOrientation = nface[i] / P4EST_FACES;
      nface[i] %= P4EST_FACES;
      nedgef[i] =
        p8est_connectivity_face_neighbor_edge_orientation (e1, f1,
                                                           nface[i],
                                                           faceOrientation);
    }
  }
  return 0;
}

/** Populate mesh information for hanging edges and edges across tree
 *  boundaries, i.e. every neighborhood scenario where we need more information
 *  (like orientation) than a single index. Note that this function only pushes
 *  an address whose data has to be set separately.
 *
 * \param [in][out] mesh     The mesh structure to which we will add edge
 *                           information
 * \param [in]      clen     Number of quadrants to be added
 * \param [in]      pequad   List of quadrant indices
 * \param [in]      peedge   List of quadrant encodings
 */
static              p4est_locidx_t
mesh_edge_allocate (p4est_mesh_t * mesh, p4est_locidx_t elen,
                    p4est_locidx_t ** pequad, int8_t ** peedge)
{
  p4est_locidx_t      edgeid, estart, eend;

  P4EST_ASSERT (elen > 0);
  P4EST_ASSERT (mesh->edge_offset->elem_count ==
                (size_t) (mesh->local_num_edges + 1));
  edgeid = mesh->local_num_edges++;
  estart = *(p4est_locidx_t *) sc_array_index (mesh->edge_offset, edgeid);
  eend = estart + elen;
  *(p4est_locidx_t *) sc_array_push (mesh->edge_offset) = eend;

  P4EST_ASSERT (mesh->edge_offset->elem_count ==
                (size_t) (mesh->local_num_edges + 1));

  P4EST_ASSERT (mesh->edge_quad->elem_count == (size_t) estart);
  *pequad = (p4est_locidx_t *) sc_array_push_count (mesh->edge_quad, elen);
  P4EST_ASSERT (mesh->edge_quad->elem_count == (size_t) eend);

  P4EST_ASSERT (mesh->edge_edge->elem_count == (size_t) estart);
  *peedge = (int8_t *) sc_array_push_count (mesh->edge_edge, elen);
  P4EST_ASSERT (mesh->edge_edge->elem_count == (size_t) eend);

  return edgeid;
}

/** Process inter-tree edge neighbors
 * \param[in]      info       General iterator information
 * \param[in]      side1      Currently processed edge
 * \param[in]      subedge_id Sub edge index in 0..1 for hanging edges, -1
 *                            for full edges
 * \param[in][out] mesh       The mesh structure that will be filled
 *                            with the edge neighbors along current edge
 * \param[in]      nftree     Tree indices of face neighbors wrt. the
 *                            current edge
 * \param[in]      nedgef     Edge indices of face neighbors wrt. the
 *                            current edge
 * \param[in]      cz         Number of adjacent trees
 * \param[in]      zz         internal number of currently processed edge
 */
static int
mesh_edge_process_inter_tree_edges (p8est_iter_edge_info_t * info,
                                    p8est_iter_edge_side_t * side1,
                                    p4est_topidx_t subedge_id,
                                    p4est_mesh_t * mesh,
                                    p4est_locidx_t * nftree,
                                    p4est_locidx_t * nedgef, int cz, int zz)
{
  int                 ignore, j, z2;
  int                 nAdjacentQuads, qid1, qid2, edgeid;
  p8est_iter_edge_side_t *side2;
  p4est_tree_t       *tree1, *tree2;
  int8_t             *eedges;
  p4est_locidx_t     *equads;
  int8_t             *peedge;
  p4est_locidx_t     *pequad;
  int                 edgeid_offset =
    mesh->local_num_quadrants + mesh->ghost_num_quadrants;

  /* sanity check for subedge_id */
  P4EST_ASSERT (-1 <= subedge_id && subedge_id < 2);

  /* overestimate number of adjacent quads:
   * Due to hanging edges we cannot avoid iterating over all
   * quadrants, if we want to know the exact number. so assume
   * that each edge is hanging (which cannot be the case; it would
   * be 2 separate edges if they were) and allocate space for
   * twice the number of adjacent quadrants.
   * However, we can consider that the quadrant itself is not
   * stored as its own neighbor. */
  nAdjacentQuads = 2 * cz - 1;

  equads = P4EST_ALLOC (p4est_locidx_t, nAdjacentQuads);
  eedges = P4EST_ALLOC (int8_t, nAdjacentQuads);

  tree1 = p4est_tree_array_index (info->p4est->trees, side1->treeid);
  if (0 <= subedge_id) {
    qid1 = side1->is.hanging.quadid[subedge_id] + tree1->quadrants_offset;
  }
  else {
    qid1 = side1->is.full.quadid + tree1->quadrants_offset;
  }
  /* Go through edge neighbors and collect the true edges */
  int                 goodones = 0;
  for (z2 = 0; z2 < cz; ++z2) {
    if (z2 == zz) {
      /* We do not count ourselves as a neighbor */
      continue;
    }
    ignore = 0;
    side2 = (p8est_iter_edge_side_t *) sc_array_index (&info->sides, z2);
    P4EST_ASSERT (0 <= side2->treeid &&
                  side2->treeid < info->p4est->connectivity->num_trees);
    P4EST_ASSERT (side2->edge >= 0 && side2->edge < P8EST_EDGES);

    /* check if current side2 is among the face neighbors:
     * We have either found it via connectivity or it is in the same
     * tree as side1.
     */
    for (j = 0; j < 2; ++j) {
      if (info->tree_boundary <= P8EST_CONNECT_EDGE) {
        if ((nedgef[j] == (int) side2->edge && nftree[j] == side2->treeid) ||
            (side1->treeid == side2->treeid)) {
          ignore = 1;
          break;
        }
      }
    }
    if (!ignore) {
      /* Record this edge neighbor if we don't ignore it */
      int8_t              localOri =
        (side1->orientation + side2->orientation) % 2;
      tree2 = p4est_tree_array_index (info->p4est->trees, side2->treeid);

      if (side1->is_hanging && side2->is_hanging) {
        /* both sides are hanging; only store the one that shares
         * the same part of the the half edge */
        int8_t              subEdgeIdx = (subedge_id + localOri) % 2;
        qid2 =
          side2->is.hanging.quadid[subEdgeIdx] +
          (side2->is.hanging.is_ghost[subEdgeIdx] ?
           mesh->local_num_quadrants : tree2->quadrants_offset);
        equads[goodones] = qid2;
        eedges[goodones] = P8EST_EDGES * localOri + (int) side2->edge;
        ++goodones;
      }
      else if (!side1->is_hanging && side2->is_hanging) {
        /* check if we have to swap hanging quads in order to store
         * them wrt. edge corners, i.e. first subquad that is adjacent
         * to corner of edge with lower corner index
         */
        for (j = 0; j < 2; ++j) {
          int                 pos = (localOri + j) % 2;
          qid2 =
            side2->is.hanging.quadid[pos] +
            (side2->is.hanging.is_ghost[pos] ?
             mesh->local_num_quadrants : tree2->quadrants_offset);
          equads[goodones] = qid2;
          eedges[goodones] = -24 + P8EST_EDGES * localOri + side2->edge;
          ++goodones;
        }
      }
      else if (side1->is_hanging && !side2->is_hanging) {
        qid2 = side2->is.full.quadid +
          (side2->is.full.is_ghost ?
           mesh->local_num_quadrants : tree2->quadrants_offset);
        int                 pos = (localOri + subedge_id) % 2;
        equads[goodones] = qid2;
        eedges[goodones] =
          24 + 24 * pos + P8EST_EDGES * localOri + side2->edge;
        ++goodones;
      }
      else {
        qid2 = side2->is.full.quadid +
          (side2->is.full.is_ghost ?
           mesh->local_num_quadrants : tree2->quadrants_offset);
        equads[goodones] = qid2;
        eedges[goodones] = P8EST_EDGES * localOri + side2->edge;
        ++goodones;
      }
    }
  }

  /* we have excluded between 0 and all quadrants. */
  P4EST_ASSERT (0 <= (size_t) goodones && (size_t) goodones < nAdjacentQuads);

  if (goodones > 0) {
    /* Allocate and fill edge information in the mesh structure */
    edgeid = mesh_edge_allocate (mesh, goodones, &pequad, &peedge);
    /* "link" to arrays encoding inter-tree edge-neighborhood */
    P4EST_ASSERT (mesh->quad_to_edge[P8EST_EDGES * qid1 + side1->edge] == -1);
    mesh->quad_to_edge[P8EST_EDGES * qid1 + side1->edge] =
      edgeid_offset + edgeid;
    /* populate allocated memory */
    memcpy (pequad, equads, goodones * sizeof (p4est_locidx_t));
    memcpy (peedge, eedges, goodones * sizeof (int8_t));
  }
  P4EST_FREE (equads);
  P4EST_FREE (eedges);
}
#endif /* P4_TO_P8 */

static void
mesh_iter_corner (p4est_iter_corner_info_t * info, void *user_data)
{
  int                 i, j;
  int                 f1, f2, code, faceOrientation;
  int                 fc1, fc2, diagonal;
#ifdef P4_TO_P8
  int                 e1, edgeOrientation;
  int                 pref, pset;
#endif /* P4_TO_P8 */
  int                 visited[P4EST_CHILDREN];
  int8_t             *pccorner;
  size_t              cz, zz;
  p4est_locidx_t      qoffset, qid1, qid2;
  p4est_locidx_t      cornerid_offset, cornerid;
  p4est_locidx_t     *pcquad;
  p4est_mesh_t       *mesh = (p4est_mesh_t *) user_data;
  p4est_iter_corner_side_t *side1, *side2;
  p4est_tree_t       *tree1, *tree2;
  p4est_connectivity_t *conn = info->p4est->connectivity;
  sc_array_t         *trees = info->p4est->trees;

  cornerid_offset = mesh->local_num_quadrants + mesh->ghost_num_quadrants;

  /* Check the case when the corner does not involve neighbors */
  cz = info->sides.elem_count;
  P4EST_ASSERT (cz > 0);
  P4EST_ASSERT (info->tree_boundary || cz == P4EST_CHILDREN);

  if (cz == 1) {
    side1 = (p4est_iter_corner_side_t *) sc_array_index_int (&info->sides, 0);
    P4EST_ASSERT (!side1->is_ghost);
    tree1 = p4est_tree_array_index (trees, side1->treeid);
    qid1 = side1->quadid + tree1->quadrants_offset;
    P4EST_ASSERT (0 <= qid1 && qid1 < mesh->local_num_quadrants);

    P4EST_ASSERT (mesh->quad_to_corner[P4EST_CHILDREN * qid1 +
                                       side1->corner] == -1);
    mesh->quad_to_corner[P4EST_CHILDREN * qid1 + side1->corner] = -3;

    return;
  }

  if (info->tree_boundary == P4EST_CONNECT_FACE) {
    /* This corner is inside an inter-tree face */
    if (cz == P4EST_HALF) {
      /* This is a tree face boundary, no corner neighbors exist */
      for (i = 0; i < P4EST_HALF; ++i) {
        side1 =
          (p4est_iter_corner_side_t *) sc_array_index_int (&info->sides, i);
        if (!side1->is_ghost) {
          tree1 = p4est_tree_array_index (trees, side1->treeid);
          qid1 = side1->quadid + tree1->quadrants_offset;
          P4EST_ASSERT (0 <= qid1 && qid1 < mesh->local_num_quadrants);

          P4EST_ASSERT (mesh->quad_to_corner[P4EST_CHILDREN * qid1 +
                                             side1->corner] == -1);
          mesh->quad_to_corner[P4EST_CHILDREN * qid1 + side1->corner] = -3;
        }
      }
      return;
    }
    P4EST_ASSERT (cz == P4EST_CHILDREN);

    /* Process a corner in pairs of diagonal inter-tree neighbors */
    memset (visited, 0, P4EST_CHILDREN * sizeof (int));
    for (i = 0; i < P4EST_HALF; ++i) {
      side1 = side2 = NULL;
      f1 = -1;
      fc1 = -1;
      qid1 = -3;
      for (j = 0; j < P4EST_CHILDREN; ++j) {
        if (visited[j]) {
          continue;
        }

        /* Remember the first side we want to pair up */
        if (side1 == NULL) {
          side1 =
            (p4est_iter_corner_side_t *) sc_array_index_int (&info->sides, j);
          P4EST_ASSERT (side1->quad != NULL);
          f1 = tree_face_quadrant_corner_face (side1->quad, side1->corner);
          fc1 = p4est_corner_face_corners[side1->corner][f1];
          P4EST_ASSERT (0 <= fc1 && fc1 < P4EST_HALF);
          tree1 = p4est_tree_array_index (trees, side1->treeid);
          qid1 = side1->quadid + (side1->is_ghost ?
                                  mesh->local_num_quadrants :
                                  tree1->quadrants_offset);
          visited[j] = 1;
          continue;
        }

        /* Examine a potential second side */
        P4EST_ASSERT (side2 == NULL);
        side2 =
          (p4est_iter_corner_side_t *) sc_array_index_int (&info->sides, j);
        P4EST_ASSERT (side2->quad != NULL);
        f2 = tree_face_quadrant_corner_face (side2->quad, side2->corner);
        if (side1->treeid == side2->treeid && f1 == f2) {
          /* Periodicity allows for equal trees and unequal faces */
          side2 = NULL;
          continue;
        }

        /* This side as in the opposite tree */
        fc2 = p4est_corner_face_corners[side2->corner][f2];
        P4EST_ASSERT (0 <= fc2 && fc2 < P4EST_HALF);
        code = conn->tree_to_face[P4EST_FACES * side1->treeid + f1];
        faceOrientation = code / P4EST_FACES;
        P4EST_ASSERT (f2 == code % P4EST_FACES);
#ifdef P4_TO_P8
        pref = p8est_face_permutation_refs[f1][f2];
        pset = p8est_face_permutation_sets[pref][faceOrientation];
        diagonal = (p8est_face_permutations[pset][fc1] ^ fc2) == 3;
#else /* P4_TO_P8 */
        diagonal = (fc1 ^ fc2) != faceOrientation;
#endif /* P4_TO_P8 */
        if (!diagonal) {
          side2 = NULL;
          continue;
        }

        /* We have found a diagonally opposite second side */
        tree2 = p4est_tree_array_index (trees, side2->treeid);
        qid2 = side2->quadid + (side2->is_ghost ? mesh->local_num_quadrants
                                : tree2->quadrants_offset);
        if (!side1->is_ghost) {
          P4EST_ASSERT (0 <= qid1 && qid1 < mesh->local_num_quadrants);
          P4EST_ASSERT (mesh->quad_to_corner[P4EST_CHILDREN * qid1 +
                                             side1->corner] == -1);
          cornerid = mesh_corner_allocate (mesh, 1, &pcquad, &pccorner);
          mesh->quad_to_corner[P4EST_CHILDREN * qid1 + side1->corner] =
            cornerid_offset + cornerid;
          *pcquad = qid2;
          *pccorner = side2->corner;
        }
        if (!side2->is_ghost) {
          P4EST_ASSERT (0 <= qid2 && qid2 < mesh->local_num_quadrants);
          P4EST_ASSERT (mesh->quad_to_corner[P4EST_CHILDREN * qid2 +
                                             side2->corner] == -1);
          cornerid = mesh_corner_allocate (mesh, 1, &pcquad, &pccorner);
          mesh->quad_to_corner[P4EST_CHILDREN * qid2 + side2->corner] =
            cornerid_offset + cornerid;
          *pcquad = qid1;
          *pccorner = side1->corner;
        }
        visited[j] = 1;
        break;
      }
      P4EST_ASSERT (side1 != NULL && side2 != NULL);
    }
    return;
  }

#ifdef P4_TO_P8
  if (info->tree_boundary == P8EST_CONNECT_EDGE) {
    if (cz == 2) {
      /* This is a tree edge boundary, no corner neighbors exist */
      for (i = 0; i < 2; ++i) {
        side1 =
          (p4est_iter_corner_side_t *) sc_array_index_int (&info->sides, i);
        if (!side1->is_ghost) {
          tree1 = p4est_tree_array_index (trees, side1->treeid);
          qid1 = side1->quadid + tree1->quadrants_offset;
          P4EST_ASSERT (0 <= qid1 && qid1 < mesh->local_num_quadrants);

          P4EST_ASSERT (mesh->quad_to_corner[P4EST_CHILDREN * qid1 +
                                             side1->corner] == -1);
          mesh->quad_to_corner[P4EST_CHILDREN * qid1 + side1->corner] = -3;
        }
      }
      return;
    }
    /* Tree corner neighbors across an edge are not implemented: set to -2 */
    for (zz = 0; zz < cz; ++zz) {
      side1 = (p4est_iter_corner_side_t *) sc_array_index (&info->sides, zz);
      if (!side1->is_ghost) {
        tree1 = p4est_tree_array_index (trees, side1->treeid);
        qid1 = side1->quadid + tree1->quadrants_offset;
        P4EST_ASSERT (0 <= qid1 && qid1 < mesh->local_num_quadrants);
        P4EST_ASSERT (mesh->quad_to_corner[P4EST_CHILDREN * qid1 +
                                           side1->corner] == -1);
        mesh->quad_to_corner[P4EST_CHILDREN * qid1 + side1->corner] = -2;
      }
    }
    return;
  }
#endif /* P4_TO_P8 */

  if (info->tree_boundary == P4EST_CONNECT_CORNER) {
    int                 c1, ncornerf[P4EST_DIM];
    int                 nface[P4EST_DIM];
#ifdef P4_TO_P8
    int8_t              which_corner;
    int                 ncornere[P4EST_DIM], nedge[P4EST_DIM];
    p4est_locidx_t      netree[P4EST_DIM];
#endif /* P4_TO_P8 */
    int                 ignore;
    size_t              z2;
    int8_t             *ccorners;
    p4est_topidx_t      t1, nftree[P4EST_DIM];
    p4est_locidx_t      goodones;
    p4est_locidx_t     *cquads;

    /* initialize ncornerf, ncornere to zero */
    SC_BZERO (ncornerf, P4EST_DIM);
#ifdef P4_TO_P8
    SC_BZERO (ncornere, P4EST_DIM);
#endif /* P4_TO_P8 */

    /* Loop through all corner sides, that is the quadrants touching it.  For
     * each of these quadrants, determine the corner sides that can potentially
     * occur by being a face neighbor as well.  Exclude these face neighbors
     * and the quadrant itself, record all others as corner neighbors.
     */
    cquads = P4EST_ALLOC (p4est_locidx_t, cz - 1);
    ccorners = P4EST_ALLOC (int8_t, cz - 1);
    for (zz = 0; zz < cz; ++zz) {
      side1 = (p4est_iter_corner_side_t *) sc_array_index (&info->sides, zz);
      if (!side1->is_ghost) {
        /* We only create corner information for processor-local quadrants */
        t1 = side1->treeid;
        c1 = side1->corner;
        tree1 = p4est_tree_array_index (trees, t1);
        qid1 = side1->quadid + tree1->quadrants_offset;
        P4EST_ASSERT (0 <= qid1 && qid1 < mesh->local_num_quadrants);
        P4EST_ASSERT (mesh->quad_to_corner[P4EST_CHILDREN * qid1 + c1] == -1);

        /* Get all quadrant faces and edges touching this corner */
        for (i = 0; i < P4EST_DIM; ++i) {
          /* begin with faces */
          f1 = p4est_corner_faces[c1][i];
          nftree[i] = conn->tree_to_tree[P4EST_FACES * t1 + f1];
          nface[i] = conn->tree_to_face[P4EST_FACES * t1 + f1];

          if (nftree[i] == t1 && nface[i] == f1) {
            /* This is a physical face boundary, no face neighbor present */
            ncornerf[i] = -1;
          }
          else {
            /* We have a face neighbor */
            faceOrientation = nface[i] / P4EST_FACES;
            nface[i] %= P4EST_FACES;
            ncornerf[i] =
              p4est_connectivity_face_neighbor_corner_orientation (c1, f1,
                                                                   nface[i],
                                                                   faceOrientation);
          }

#ifdef P4_TO_P8
          /* Check all quadrant edges that touch this corner */
          e1 = p8est_corner_edges[c1][i];
          if (conn->edge_to_tree != 0) {
            netree[i] = conn->edge_to_tree[P8EST_EDGES * t1 + e1];
            nedge[i] = conn->edge_to_edge[P8EST_EDGES * t1 + e1];
            if (netree[i] == t1 && nedge[i] == e1) {
              /* This is a physical edge boundary, no edge neighbor present */
              ncornere[i] = -1;
            }

            if (ncornere[i] == 0) {
              /* We have an edge neighbor */
              edgeOrientation = nedge[i] / P8EST_EDGES;
              nedge[i] %= P8EST_EDGES;
              which_corner = (e1 == p8est_edge_corners[c1][0] ? 0 : 1);
              which_corner = (which_corner + edgeOrientation) % 2;
              ncornere[i] = p8est_edge_corners[nedge[i]][which_corner];
            }
          }
          else {
            for (j = 0; j < cz; ++j) {
              side1 =
                (p4est_iter_corner_side_t *) sc_array_index_int (&info->sides,
                                                                 i);
              if (!side1->is_ghost) {
                tree1 = p4est_tree_array_index (trees, side1->treeid);
                qid1 = side1->quadid + tree1->quadrants_offset;
                P4EST_ASSERT (0 <= qid1 && qid1 < mesh->local_num_quadrants);

                P4EST_ASSERT (mesh->quad_to_corner[P4EST_CHILDREN * qid1 +
                                                   side1->corner] == -1);
                mesh->quad_to_corner[P4EST_CHILDREN * qid1 + side1->corner] =
                  -3;
              }
              P4EST_FREE (cquads);
              P4EST_FREE (ccorners);
              return;
            }
          }
#endif /* P4_TO_P8 */
        }

        /* Go through corner neighbors and collect the true corners */
        goodones = 0;
        for (z2 = 0; z2 < cz; ++z2) {
          if (z2 == zz) {
            /* We do not count ourselves as a neighbor */
            continue;
          }
          ignore = 0;
          side2 =
            (p4est_iter_corner_side_t *) sc_array_index (&info->sides, z2);
          P4EST_ASSERT (side2->corner >= 0);
          for (i = 0; i < P4EST_DIM; ++i) {
            /* 2D: Ignore if this is one of the face neighbors' corners */
            /* 3D: Ignore if this is either one of the face neighbors' or one of
             *     the edge neighbors' corners */
#ifdef P4_TO_P8
            if ((ncornerf[i] == (int) side2->corner &&
                 nftree[i] == side2->treeid) ||
                (ncornere[i] == (int) side2->corner &&
                 netree[i] == side2->treeid))
#else /* P4_TO_P8 */
            if (ncornerf[i] == (int) side2->corner &&
                nftree[i] == side2->treeid)
#endif /* P4_TO_P8 */
            {
              ignore = 1;
              break;
            }
          }
          if (!ignore) {
            /* Record this corner neighbor if we don't ignore it */
            tree2 = p4est_tree_array_index (trees, side2->treeid);
            qid2 =
              side2->quadid +
              (side2->is_ghost ?
               mesh->local_num_quadrants : tree2->quadrants_offset);
            cquads[goodones] = qid2;
            ccorners[goodones] = (int) side2->corner;
            ++goodones;
          }
        }

        /* we have excluded between 0 and all quadrants. */
        P4EST_ASSERT (0 <= (size_t) goodones && (size_t) goodones < cz);

        if (goodones > 0) {
          /* Allocate and fill corner information in the mesh structure */
          cornerid =
            mesh_corner_allocate (mesh, goodones, &pcquad, &pccorner);
          /* "link" to arrays encoding inter-tree corner-neighborhood */
          P4EST_ASSERT (mesh->quad_to_corner[P4EST_CHILDREN * qid1 + c1] ==
                        -1);
          mesh->quad_to_corner[P4EST_CHILDREN * qid1 + c1] =
            cornerid_offset + cornerid;
          /* populate allocated memory */
          memcpy (pcquad, cquads, goodones * sizeof (p4est_locidx_t));
          memcpy (pccorner, ccorners, goodones * sizeof (int8_t));
        }
        else {
          for (i = 0; i < cz; ++i) {
            side1 =
              (p4est_iter_corner_side_t *) sc_array_index_int (&info->sides,
                                                               i);
            if (!side1->is_ghost) {
              tree1 = p4est_tree_array_index (trees, side1->treeid);
              qid1 = side1->quadid + tree1->quadrants_offset;
              P4EST_ASSERT (0 <= qid1 && qid1 < mesh->local_num_quadrants);

              P4EST_ASSERT (mesh->quad_to_corner[P4EST_CHILDREN * qid1 +
                                                 side1->corner] == -1);
              mesh->quad_to_corner[P4EST_CHILDREN * qid1 + side1->corner] =
                -3;
            }
          }
          P4EST_FREE (cquads);
          P4EST_FREE (ccorners);
          return;
        }
      }
    }
    P4EST_FREE (cquads);
    P4EST_FREE (ccorners);
    return;
  }

  /* Process a corner inside the tree in pairs of diagonal neighbors */
  P4EST_ASSERT (!info->tree_boundary);
  side1 = (p4est_iter_corner_side_t *) sc_array_index (&info->sides, 0);
  tree1 = p4est_tree_array_index (trees, side1->treeid);
  qoffset = tree1->quadrants_offset;
  memset (visited, 0, P4EST_CHILDREN * sizeof (int));
  for (i = 0; i < P4EST_HALF; ++i) {
    side1 = side2 = NULL;
    qid1 = -3;
    for (j = 0; j < P4EST_CHILDREN; ++j) {
      if (visited[j]) {
        continue;
      }

      /* Remember the first side we want to pair up */
      if (side1 == NULL) {
        side1 =
          (p4est_iter_corner_side_t *) sc_array_index_int (&info->sides, j);
        qid1 = side1->quadid +
          (side1->is_ghost ? mesh->local_num_quadrants : qoffset);
        visited[j] = 1;
        continue;
      }

      /* Examine a potential second side */
      /* The diagonally opposite corner is found by adding the corner indices
       * of z-ordering convention. diagonally opposite corner indeces add up
       * to 3 in 2D and to 7 in 3D, i.e. P4EST_CHILDREN - 1. */
      P4EST_ASSERT (side2 == NULL);
      side2 =
        (p4est_iter_corner_side_t *) sc_array_index_int (&info->sides, j);
      P4EST_ASSERT (side1->treeid == side2->treeid);
      if (side1->corner + side2->corner != P4EST_CHILDREN - 1) {
        side2 = NULL;
        continue;
      }

      /* We have found a diagonally opposite second side */
      qid2 = side2->quadid +
        (side2->is_ghost ? mesh->local_num_quadrants : qoffset);
      if (!side1->is_ghost) {
        P4EST_ASSERT (0 <= qid1 && qid1 < mesh->local_num_quadrants);
        P4EST_ASSERT (mesh->quad_to_corner[P4EST_CHILDREN * qid1 +
                                           side1->corner] == -1);
        mesh->quad_to_corner[P4EST_CHILDREN * qid1 + side1->corner] = qid2;
      }
      if (!side2->is_ghost) {
        P4EST_ASSERT (0 <= qid2 && qid2 < mesh->local_num_quadrants);
        P4EST_ASSERT (mesh->quad_to_corner[P4EST_CHILDREN * qid2 +
                                           side2->corner] == -1);
        mesh->quad_to_corner[P4EST_CHILDREN * qid2 + side2->corner] = qid1;
      }
      visited[j] = 1;
      break;
    }
    P4EST_ASSERT (side1 != NULL && side2 != NULL);
  }
}

#ifdef P4_TO_P8
static void
mesh_iter_edge (p8est_iter_edge_info_t * info, void *user_data)
{
  int8_t              visited[P4EST_HALF];
  size_t              i, j, k, cz, zz;
  int                 swapsides;
  p4est_locidx_t      qid1, qid2, qls1[2], qoffset;
  p4est_locidx_t      eid1, eid2;
  p4est_locidx_t      edgeid;
  p4est_locidx_t      in_qtoe, edgeid_offset;
  p4est_mesh_t       *mesh = (p4est_mesh_t *) user_data;
  p8est_iter_edge_side_t *side1, *side2, *tempside;
  p4est_tree_t       *tree1, *tree2;

  edgeid_offset = mesh->local_num_quadrants + mesh->ghost_num_quadrants;

  /* general sanity checks */
  cz = info->sides.elem_count;
  P4EST_ASSERT (cz > 0);
  P4EST_ASSERT (info->tree_boundary || cz == P4EST_HALF);

  /* edge limits domain or is located on a face limitting the domain */
  if (cz <= 2) {
    /* sanity checks */
    if (cz == 1) {
      P4EST_ASSERT (info->tree_boundary);
    }
    for (i = 0; i < cz; ++i) {
      side1 = (p8est_iter_edge_side_t *) sc_array_index_int (&info->sides, i);
      P4EST_ASSERT (0 <= side1->treeid
                    && side1->treeid < info->p4est->connectivity->num_trees);
      P4EST_ASSERT (0 <= side1->edge && side1->edge < P8EST_EDGES);
      if (!side1->is_hanging) {
        if (!side1->is.full.is_ghost) {
          tree1 = p4est_tree_array_index (info->p4est->trees, side1->treeid);
          qid1 = side1->is.full.quadid + tree1->quadrants_offset;

          P4EST_ASSERT (0 <= qid1 && qid1 < mesh->local_num_quadrants);
          P4EST_ASSERT
            (mesh->quad_to_edge[P8EST_EDGES * qid1 + side1->edge] == -1);

          mesh->quad_to_edge[P8EST_EDGES * qid1 + side1->edge] = -3;
        }
      }
      else {
        for (j = 0; j < 2; ++j) {
          if (!side1->is.hanging.is_ghost[j]) {
            tree1 =
              p4est_tree_array_index (info->p4est->trees, side1->treeid);
            qid1 = side1->is.hanging.quadid[j] + tree1->quadrants_offset;

            P4EST_ASSERT (0 <= qid1 && qid1 < mesh->local_num_quadrants);

            P4EST_ASSERT
              (mesh->quad_to_edge[P8EST_EDGES * qid1 + side1->edge] == -1);
            mesh->quad_to_edge[P8EST_EDGES * qid1 + side1->edge] = -3;
          }
        }
      }
    }
    return;
  }
  else {
    /* edges inside of domain have at least 4 adjacent quadrants */
    P4EST_ASSERT (4 <= cz);

    /* edges on tree boundaries */
    if (info->tree_boundary) {
      p4est_locidx_t      nedgef[2];
      p4est_locidx_t      nface[2];
      int                 ignore;
      size_t              z2;
      p4est_topidx_t      nftree[2];
      p4est_locidx_t      goodones;
      int8_t             *eedges;
      p4est_locidx_t     *equads;

      /* initialize nedgef to zero */
      SC_BZERO (nedgef, 2);

      /* Loop through all edge sides, that is the quadrants touching it.  For
       * each of these quadrants, determine the edge sides that can potentially
       * occur by being a face neighbor as well.  Exclude these face neighbors
       * and the quadrant itself, record all others as edge neighbors.
       */
      for (zz = 0; zz < cz; ++zz) {
        side1 = (p8est_iter_edge_side_t *) sc_array_index (&info->sides, zz);
        P4EST_ASSERT (0 <= side1->treeid &&
                      side1->treeid < info->p4est->connectivity->num_trees);
        P4EST_ASSERT (0 <= side1->edge && side1->edge < P8EST_EDGES);

        /* look for face neighbors of current quadrant */
        mesh_edge_find_face_neighbors (side1,
                                       info->p4est->connectivity,
                                       nftree, nface, nedgef);

        /* We only create edge information for processor-local quadrants */
        if (side1->is_hanging) {
          for (i = 0; i < 2; ++i) {
            if (!side1->is.hanging.is_ghost[i]) {
              mesh_edge_process_inter_tree_edges (info, side1, i, mesh,
                                                  nftree, nedgef, cz, zz);
            }
          }
        }
        else {
          /* side1 is full */
          if (!side1->is.full.is_ghost) {
            mesh_edge_process_inter_tree_edges (info, side1, -1, mesh, nftree,
                                                nedgef, cz, zz);

          }
        }
      }
      return;
    }

    /* intra-tree */
    else {
      int8_t             *peedge;
      p4est_locidx_t     *pequad;

      P4EST_ASSERT (!info->tree_boundary);

      side1 = (p8est_iter_edge_side_t *) sc_array_index (&info->sides, 0);
      tree1 = p4est_tree_array_index (info->p4est->trees, side1->treeid);
      qoffset = tree1->quadrants_offset;
      memset (visited, 0, P4EST_HALF * sizeof (int8_t));

      for (i = 0; i < 0.5 * cz; ++i) {
        side1 = side2 = NULL;
        qid1 = -1;

        for (j = 0; j < P4EST_HALF; ++j) {
          if (visited[j]) {
            continue;
          }

          /* remember first side */
          if (side1 == NULL) {
            side1 =
              (p8est_iter_edge_side_t *) sc_array_index_int (&info->sides, j);
            eid1 = side1->edge;
            visited[j] = 1;
            continue;
          }

          /* Examine second side */
          P4EST_ASSERT (side2 == NULL);
          side2 =
            (p8est_iter_edge_side_t *) sc_array_index_int (&info->sides, j);
          eid2 = side2->edge;
          P4EST_ASSERT (side2->treeid == side1->treeid);

          /* edge is diagonal opposite if e1 XOR 3 = e2 */
          if ((eid1 ^ 3) != eid2) {
            side2 = NULL;
            continue;
          }

          /* we found a pair of diagonal edges */
          if (!side1->is_hanging && !side2->is_hanging) {
            /* case 1: equally sized cells */
            /* this is the standard case. we do not store any encoding for
             * this case */
            if (side1->is.full.is_ghost && side2->is.full.is_ghost) {
              continue;
            }

            /* get correct id depending on whether cell is a ghost or not */
            if (!side1->is.full.is_ghost) {
              tree1 =
                p4est_tree_array_index (info->p4est->trees, side1->treeid);
              qid1 = side1->is.full.quadid + tree1->quadrants_offset;
              P4EST_ASSERT (0 <= qid1 && qid1 < mesh->local_num_quadrants);
            }
            else {
              P4EST_ASSERT (side1->is.full.quad != NULL);
              P4EST_ASSERT (side1->is.full.quadid >= 0);
              qid1 = mesh->local_num_quadrants + side1->is.full.quadid;
              P4EST_ASSERT (mesh->local_num_quadrants <= qid1 &&
                            qid1 < mesh->local_num_quadrants +
                            mesh->ghost_num_quadrants);
            }
            if (!side2->is.full.is_ghost) {
              tree2 =
                p4est_tree_array_index (info->p4est->trees, side2->treeid);
              qid2 = side2->is.full.quadid + tree2->quadrants_offset;
              P4EST_ASSERT (0 <= qid2 && qid2 < mesh->local_num_quadrants);
            }
            else {
              P4EST_ASSERT (side2->is.full.quad != NULL);
              P4EST_ASSERT (side2->is.full.quadid >= 0);
              qid2 = mesh->local_num_quadrants + side2->is.full.quadid;
              P4EST_ASSERT (mesh->local_num_quadrants <= qid2 &&
                            qid2 < mesh->local_num_quadrants +
                            mesh->ghost_num_quadrants);
            }

            /* write info to correct position in quad_to_edge array if cell
             * is not part of the ghost layer */
            if (!side1->is.full.is_ghost) {
              P4EST_ASSERT (mesh->quad_to_edge[P8EST_EDGES * qid1 +
                                               side1->edge] == -1);
              mesh->quad_to_edge[P8EST_EDGES * qid1 + side1->edge] = qid2;
            }
            if (!side2->is.full.is_ghost) {
              P4EST_ASSERT (mesh->quad_to_edge[P8EST_EDGES * qid2 +
                                               side2->edge] == -1);
              mesh->quad_to_edge[P8EST_EDGES * qid2 + side2->edge] = qid1;
            }
          }
          else if ((!side1->is_hanging && side2->is_hanging) ||
                   (side1->is_hanging && !side2->is_hanging)) {
            /* one side is hanging. assert it is always side2. */
            swapsides = side1->is_hanging;
            if (swapsides) {
              tempside = side1;
              side1 = side2;
              side2 = tempside;
            }
            P4EST_ASSERT (!side1->is_hanging && side2->is_hanging);

            /* determine quadrant number for non-hanging large quadrant */
            if (!side1->is.full.is_ghost) {
              tree1 =
                p4est_tree_array_index (info->p4est->trees, side1->treeid);
              qid1 = side1->is.full.quadid + tree1->quadrants_offset;
              P4EST_ASSERT (0 <= qid1 && qid1 < mesh->local_num_quadrants);
            }
            else {
              P4EST_ASSERT (side1->is.full.quad != NULL);
              P4EST_ASSERT (side1->is.full.quadid >= 0);
              qid1 = mesh->local_num_quadrants + side1->is.full.quadid;
              P4EST_ASSERT (mesh->local_num_quadrants <= qid1
                            && qid1 < mesh->local_num_quadrants +
                            mesh->ghost_num_quadrants);
            }

            /* determine quadrant numbers for both hanging quadrants */
            for (k = 0; k < 2; ++k) {
              if (!side2->is.hanging.is_ghost[k]) {
                tree2 =
                  p4est_tree_array_index (info->p4est->trees, side2->treeid);
                qls1[k] =
                  side2->is.hanging.quadid[k] + tree2->quadrants_offset;
                P4EST_ASSERT (0 <= qls1[k]
                              && qls1[k] < mesh->local_num_quadrants);
              }
              else {
                P4EST_ASSERT (side2->is.hanging.quad[k] != NULL);
                P4EST_ASSERT (side2->is.hanging.quadid[k] >= 0);
                qls1[k] =
                  mesh->local_num_quadrants + side2->is.hanging.quadid[k];
                P4EST_ASSERT (qls1[k] >= mesh->local_num_quadrants
                              && qls1[k] <
                              mesh->local_num_quadrants +
                              mesh->ghost_num_quadrants);
              }
            }

            /* encode quadrant neighborhood:
             * inside an octree the orientation is always 0. All we need to
             * encode is the information, if it is the first or second quadrant
             * along the edge */
            /* again: check before writing that values are untouched */
            if (!side1->is.full.is_ghost) {
              /* determine position for quad_to_edge */
              in_qtoe = P8EST_EDGES * qid1 + side1->edge;

              /* allocate space for encodings */
              int8_t             *eedges;
              eedges = P4EST_ALLOC (int8_t, 2);

              /* assert there is no entry for this cell yet */
              P4EST_ASSERT (mesh->quad_to_edge[in_qtoe] == -1);

              /* calculate encodings */
              for (k = 0; k < 2; ++k) {
                eedges[k] = -24 + side2->edge;
              }

              edgeid = mesh_edge_allocate (mesh, 2, &pequad, &peedge);
              memcpy (pequad, qls1, 2 * sizeof (p4est_locidx_t));
              memcpy (peedge, eedges, 2 * sizeof (int8_t));

              P4EST_FREE (eedges);

              /* refer to "special" arrays */
              mesh->quad_to_edge[in_qtoe] = edgeid_offset + edgeid;
            }

            for (k = 0; k < 2; ++k) {
              if (!side2->is.hanging.is_ghost[k]) {
                in_qtoe = P8EST_EDGES * qls1[k] + side2->edge;

                P4EST_ASSERT (mesh->quad_to_edge[in_qtoe] == -1);

                edgeid = mesh_edge_allocate (mesh, 1, &pequad, &peedge);
                *pequad = qid1;

                /* orientation is 0 as we are in the same tree */
                *peedge = 24 + 24 * k + side1->edge;

                mesh->quad_to_edge[in_qtoe] = edgeid_offset + edgeid;
              }
            }
          }
          else {
            /* both sides are hanging with respect to a bigger edge that is not
             * diagonally opposite */
            /* determine quadrant numbers for both "hanging" edges and write
             * them directly to the corresponding positions */
            for (k = 0; k < 2; ++k) {
              if (!side1->is.hanging.is_ghost[k]) {
                tree1 =
                  p4est_tree_array_index (info->p4est->trees, side1->treeid);
                qid1 = side1->is.hanging.quadid[k] + tree1->quadrants_offset;
                P4EST_ASSERT (0 <= qid1 && qid1 < mesh->local_num_quadrants);
              }
              else {
                P4EST_ASSERT (side1->is.hanging.quad[k] != NULL);
                P4EST_ASSERT (side1->is.hanging.quadid[k] >= 0);
                qid1 =
                  mesh->local_num_quadrants + side1->is.hanging.quadid[k];
                P4EST_ASSERT (mesh->local_num_quadrants <= qid1
                              && qid1 <
                              mesh->local_num_quadrants +
                              mesh->ghost_num_quadrants);
              }

              if (!side2->is.hanging.is_ghost[k]) {
                tree2 =
                  p4est_tree_array_index (info->p4est->trees, side2->treeid);
                qid2 = side2->is.hanging.quadid[k] + tree2->quadrants_offset;
                P4EST_ASSERT (0 <= qid2 && qid2 < mesh->local_num_quadrants);
              }
              else {
                P4EST_ASSERT (side2->is.hanging.quad[k] != NULL);
                P4EST_ASSERT (side2->is.hanging.quadid[k] >= 0);
                qid2 =
                  mesh->local_num_quadrants + side2->is.hanging.quadid[k];
                P4EST_ASSERT (mesh->local_num_quadrants <= qid2
                              && qid2 <
                              mesh->local_num_quadrants +
                              mesh->ghost_num_quadrants);
              }

              if (!side1->is.hanging.is_ghost[k]) {
                in_qtoe = P8EST_EDGES * qid1 + side1->edge;
                P4EST_ASSERT (mesh->quad_to_edge[in_qtoe] == -1);
                mesh->quad_to_edge[in_qtoe] = qid2;
              }
              if (!side2->is.hanging.is_ghost[k]) {
                in_qtoe = P8EST_EDGES * qid2 + side2->edge;
                P4EST_ASSERT (mesh->quad_to_edge[in_qtoe] == -1);
                mesh->quad_to_edge[in_qtoe] = qid1;
              }
            }
          }
          visited[j] = 1;
          side1 = side2 = NULL;
        }
      }
    }
  }
}
#endif /* P4_TO_P8 */

static void
mesh_iter_face (p4est_iter_face_info_t * info, void *user_data)
{
  int                 h;
  int                 swapsides;
  p4est_mesh_t       *mesh = (p4est_mesh_t *) user_data;
  p4est_locidx_t      jl, jl2, jls[P4EST_HALF];
  p4est_locidx_t      in_qtoq, halfindex;
  p4est_locidx_t     *halfentries;
  p4est_tree_t       *tree;
  p4est_iter_face_side_t *side, *side2, *tempside;

  if (info->sides.elem_count == 1) {
    /* this face is on an outside boundary of the forest */
    P4EST_ASSERT (info->orientation == 0);
    P4EST_ASSERT (info->tree_boundary);

    side = (p4est_iter_face_side_t *) sc_array_index (&info->sides, 0);
    P4EST_ASSERT (0 <= side->treeid &&
                  side->treeid < info->p4est->connectivity->num_trees);
    P4EST_ASSERT (0 <= side->face && side->face < P4EST_FACES);
    P4EST_ASSERT (!side->is_hanging && !side->is.full.is_ghost);

    tree = p4est_tree_array_index (info->p4est->trees, side->treeid);
    jl = side->is.full.quadid + tree->quadrants_offset;
    P4EST_ASSERT (0 <= jl && jl < mesh->local_num_quadrants);
    in_qtoq = P4EST_FACES * jl + side->face;
    mesh->quad_to_quad[in_qtoq] = jl;   /* put in myself and my own face */
    mesh->quad_to_face[in_qtoq] = side->face;
  }
  else {
    /* this face is between two quadrants */
    P4EST_ASSERT (info->orientation == 0 || info->tree_boundary);
    P4EST_ASSERT (info->sides.elem_count == 2);
    side = (p4est_iter_face_side_t *) sc_array_index (&info->sides, 0);
    side2 = (p4est_iter_face_side_t *) sc_array_index (&info->sides, 1);
    P4EST_ASSERT (info->tree_boundary || side->treeid == side2->treeid);
    P4EST_ASSERT (!side->is_hanging || !side2->is_hanging);
    if (!side->is_hanging && !side2->is_hanging) {
      /* same-size face neighbors */
      P4EST_ASSERT (!side->is.full.is_ghost || !side2->is.full.is_ghost);

      /* determine both quadrant numbers */
      if (!side->is.full.is_ghost) {
        tree = p4est_tree_array_index (info->p4est->trees, side->treeid);
        jl = side->is.full.quadid + tree->quadrants_offset;
        P4EST_ASSERT (0 <= jl && jl < mesh->local_num_quadrants);
      }
      else {
        P4EST_ASSERT (side->is.full.quad != NULL);
        P4EST_ASSERT (side->is.full.quadid >= 0);
        jl = mesh->local_num_quadrants + side->is.full.quadid;
      }
      if (!side2->is.full.is_ghost) {
        tree = p4est_tree_array_index (info->p4est->trees, side2->treeid);
        jl2 = side2->is.full.quadid + tree->quadrants_offset;
        P4EST_ASSERT (0 <= jl2 && jl2 < mesh->local_num_quadrants);
      }
      else {
        P4EST_ASSERT (side2->is.full.quad != NULL);
        P4EST_ASSERT (side2->is.full.quadid >= 0);
        jl2 = mesh->local_num_quadrants + side2->is.full.quadid;
      }

      /* encode quadrant neighborhood */
      /* check that nothing has been written to that spot, i.e. nothing is
       * overwritten and at the position we want to write to we still find the
       * initially set values. */
      if (!side->is.full.is_ghost) {
        in_qtoq = P4EST_FACES * jl + side->face;
        P4EST_ASSERT (mesh->quad_to_quad[in_qtoq] == -1);
        P4EST_ASSERT (mesh->quad_to_face[in_qtoq] == -25);
        mesh->quad_to_quad[in_qtoq] = jl2;
        mesh->quad_to_face[in_qtoq] =
          P4EST_FACES * info->orientation + side2->face;
      }
      if (!side2->is.full.is_ghost) {
        in_qtoq = P4EST_FACES * jl2 + side2->face;
        P4EST_ASSERT (mesh->quad_to_quad[in_qtoq] == -1);
        P4EST_ASSERT (mesh->quad_to_face[in_qtoq] == -25);
        mesh->quad_to_quad[in_qtoq] = jl;
        mesh->quad_to_face[in_qtoq] =
          P4EST_FACES * info->orientation + side->face;
      }
    }
    else {
      /* one of the faces is hanging, rename so it's always side2 */
      swapsides = side->is_hanging;
      if (swapsides) {
        tempside = side;
        side = side2;
        side2 = tempside;
      }
      P4EST_ASSERT (!side->is_hanging && side2->is_hanging);

      /* determine quadrant number for non-hanging large face */
      if (!side->is.full.is_ghost) {
        tree = p4est_tree_array_index (info->p4est->trees, side->treeid);
        jl = side->is.full.quadid + tree->quadrants_offset;
        P4EST_ASSERT (0 <= jl && jl < mesh->local_num_quadrants);
      }
      else {
        P4EST_ASSERT (side->is.full.quad != NULL);
        P4EST_ASSERT (side->is.full.quadid >= 0);
        jl = mesh->local_num_quadrants + side->is.full.quadid;
      }

      /* determine quadrant numbers for all hanging faces */
      for (h = 0; h < P4EST_HALF; ++h) {
        int                 pos =
          p4est_connectivity_face_neighbor_face_corner_orientation (h,
                                                                    side->
                                                                    face,
                                                                    side2->
                                                                    face,
                                                                    info->
                                                                    orientation);
        if (!side2->is.hanging.is_ghost[h]) {
          tree = p4est_tree_array_index (info->p4est->trees, side2->treeid);
          jls[h] = side2->is.hanging.quadid[pos] + tree->quadrants_offset;
          P4EST_ASSERT (0 <= jls[h] && jls[h] < mesh->local_num_quadrants);
        }
        else {
          P4EST_ASSERT (side2->is.hanging.quad[pos] != NULL);
          P4EST_ASSERT (side2->is.hanging.quadid[pos] >= 0);
          jls[h] = mesh->local_num_quadrants + side2->is.hanging.quadid[pos];
        }
      }

      /* encode quadrant neighborhood */
      /* again: check before writing that corresponding values are untouched */
      if (!side->is.full.is_ghost) {
        in_qtoq = P4EST_FACES * jl + side->face;
        P4EST_ASSERT (mesh->quad_to_quad[in_qtoq] == -1);
        P4EST_ASSERT (mesh->quad_to_face[in_qtoq] == -25);
        halfindex = (p4est_locidx_t) mesh->quad_to_half->elem_count;
        mesh->quad_to_quad[in_qtoq] = halfindex;
        mesh->quad_to_face[in_qtoq] =
          P4EST_FACES * (info->orientation - P4EST_HALF) + side2->face;
        halfentries = (p4est_locidx_t *) sc_array_push (mesh->quad_to_half);
        for (h = 0; h < P4EST_HALF; ++h) {
          halfentries[h] = jls[h];
        }
      }
      for (h = 0; h < P4EST_HALF; ++h) {
        int                 pos =
          p4est_connectivity_face_neighbor_face_corner_orientation (h,
                                                                    side->
                                                                    face,
                                                                    side2->
                                                                    face,
                                                                    info->
                                                                    orientation);
        if (!side2->is.hanging.is_ghost[pos]) {
          in_qtoq = P4EST_FACES * jls[h] + side2->face;
          P4EST_ASSERT (mesh->quad_to_quad[in_qtoq] == -1);
          P4EST_ASSERT (mesh->quad_to_face[in_qtoq] == -25);
          mesh->quad_to_quad[in_qtoq] = jl;
          mesh->quad_to_face[in_qtoq] =
            P4EST_FACES * (info->orientation + (h + 1) * P4EST_HALF) +
            side->face;
        }
      }
    }
  }
}

static void
mesh_iter_volume (p4est_iter_volume_info_t * info, void *user_data)
{
  p4est_mesh_t       *mesh = (p4est_mesh_t *) user_data;
  p4est_tree_t       *tree;
  p4est_locidx_t     *quadid, qid;
  int                 level = info->quad->level, i;

  /* We could use a static quadrant counter, but that gets uglier */
  tree = p4est_tree_array_index (info->p4est->trees, info->treeid);
  P4EST_ASSERT (0 <= info->quadid &&
                info->quadid < (p4est_locidx_t) tree->quadrants.elem_count);

  if (mesh->quad_to_tree != NULL) {
    mesh->quad_to_tree[tree->quadrants_offset + info->quadid] = info->treeid;
  }

  if (mesh->quad_level != NULL) {
    quadid = (p4est_locidx_t *) sc_array_push (mesh->quad_level + level);
    qid = tree->quadrants_offset + info->quadid;
    *quadid = qid;
  }
}

size_t
p4est_mesh_memory_used (p4est_mesh_t * mesh)
{
  size_t              lqz, ngz;
  int                 level;
  size_t              qtt_memory = 0;
  size_t              ql_memory = 0;
  size_t              all_memory;

  lqz = (size_t) mesh->local_num_quadrants;
  ngz = (size_t) mesh->ghost_num_quadrants;

  if (mesh->quad_to_tree != NULL) {
    qtt_memory = sizeof (p4est_locidx_t) * lqz;
  }

  if (mesh->quad_level != NULL) {
    ql_memory = sizeof (sc_array_t) * (P4EST_QMAXLEVEL + 1);
    for (level = 0; level <= P4EST_QMAXLEVEL; ++level) {
      ql_memory += sc_array_memory_used (mesh->quad_level + level, 0);
    }
  }

  /* basic memory plus face information */
  all_memory =
    sizeof (p4est_mesh_t) + qtt_memory + ql_memory +
    P4EST_FACES * lqz * (sizeof (p4est_locidx_t) + sizeof (int8_t)) +
    ngz * sizeof (int) + sc_array_memory_used (mesh->quad_to_half, 1);

  /* add edge information */

  /* add corner information */
  if (mesh->quad_to_corner != NULL) {
    all_memory +=
      P4EST_CHILDREN * lqz * sizeof (p4est_locidx_t) +
      sc_array_memory_used (mesh->corner_offset, 1) +
      sc_array_memory_used (mesh->corner_quad, 1) +
      sc_array_memory_used (mesh->corner_corner, 1);
  }

  return all_memory;
}

p4est_mesh_t       *
p4est_mesh_new (p4est_t * p4est, p4est_ghost_t * ghost,
                p4est_connect_type_t btype)
{
  return p4est_mesh_new_ext (p4est, ghost, 0, 0, btype);
}

p4est_mesh_t       *
p4est_mesh_new_ext (p4est_t * p4est, p4est_ghost_t * ghost,
                    int compute_tree_index, int compute_level_lists,
                    p4est_connect_type_t btype)
{
  int                 do_corner = 0;
#ifdef P4_TO_P8
  int                 do_edge = 0;
#endif /* P4_TO_P8 */
  int                 do_volume = 0;
  int                 rank;
  p4est_locidx_t      lq, ng;
  p4est_locidx_t      jl;
  p4est_mesh_t       *mesh;

  /* check whether input condition for p4est is met */
#ifdef P4EST_DEBUG
  P4EST_ASSERT (p4est_is_balanced (p4est, btype));
#endif /* P4EST_DEBUG */

  mesh = P4EST_ALLOC_ZERO (p4est_mesh_t, 1);

  /* number of local quadrants and number of local ghost cells */
  lq = mesh->local_num_quadrants = p4est->local_num_quadrants;
  ng = mesh->ghost_num_quadrants = (p4est_locidx_t) ghost->ghosts.elem_count;

  /* decide which callback function have to be activated */
#ifdef P4_TO_P8
  if (btype >= P8EST_CONNECT_EDGE) {
    do_edge = 1;
  }
#endif
  if (btype >= P4EST_CONNECT_FULL) {
    do_corner = 1;
  }
  do_volume = compute_tree_index || compute_level_lists;

  /* Optional map of tree index for each quadrant */
  if (compute_tree_index) {
    mesh->quad_to_tree = P4EST_ALLOC (p4est_topidx_t, lq);
  }

  /* allocate data structures */
  mesh->ghost_to_proc = P4EST_ALLOC (int, ng);
  mesh->quad_to_quad = P4EST_ALLOC (p4est_locidx_t, P4EST_FACES * lq);
  mesh->quad_to_face = P4EST_ALLOC (int8_t, P4EST_FACES * lq);
  mesh->quad_to_half = sc_array_new (P4EST_HALF * sizeof (p4est_locidx_t));

  /* Allocate optional per-level lists of quadrants */
  if (compute_level_lists) {
    mesh->quad_level = P4EST_ALLOC (sc_array_t, P4EST_QMAXLEVEL + 1);

    for (jl = 0; jl <= P4EST_QMAXLEVEL; ++jl) {
      sc_array_init (mesh->quad_level + jl, sizeof (p4est_locidx_t));
    }
  }

  /* Populate ghost information */
  rank = 0;
  for (jl = 0; jl < ng; ++jl) {
    while (ghost->proc_offsets[rank + 1] <= jl) {
      ++rank;
      P4EST_ASSERT (rank < p4est->mpisize);
    }
    mesh->ghost_to_proc[jl] = rank;
  }

  /* Fill face arrays with default values */
  memset (mesh->quad_to_quad, -1, P4EST_FACES * lq * sizeof (p4est_locidx_t));
  memset (mesh->quad_to_face, -25, P4EST_FACES * lq * sizeof (int8_t));

#ifdef P4_TO_P8
  if (do_edge) {
    /* Allocate optional lists for edge information */
    mesh->quad_to_edge = P4EST_ALLOC (p4est_locidx_t, P8EST_EDGES * lq);
    mesh->edge_offset = sc_array_new (sizeof (p4est_locidx_t));
    mesh->edge_quad = sc_array_new (sizeof (p4est_locidx_t));
    mesh->edge_edge = sc_array_new (sizeof (int8_t));

    /* Initialize lists with default values */
    memset (mesh->quad_to_edge,
            -1, P8EST_EDGES * lq * sizeof (p4est_locidx_t));
    *(p4est_locidx_t *) sc_array_push (mesh->edge_offset) = 0;
  }
#endif /* P4_TO_P8 */

  /* Allocate optional lists of corner information */
  if (do_corner) {
    /* Initialize corner information to a consistent state */
    mesh->quad_to_corner = P4EST_ALLOC (p4est_locidx_t, P4EST_CHILDREN * lq);
    memset (mesh->quad_to_corner, -1,
            P4EST_CHILDREN * lq * sizeof (p4est_locidx_t));

    mesh->corner_offset = sc_array_new (sizeof (p4est_locidx_t));
    *(p4est_locidx_t *) sc_array_push (mesh->corner_offset) = 0;

    mesh->corner_quad = sc_array_new (sizeof (p4est_locidx_t));
    mesh->corner_corner = sc_array_new (sizeof (int8_t));
  }

  /* Call the forest iterator to collect face connectivity */
  p4est_iterate (p4est,         /* p4est */
                 ghost,         /* ghost layer */
                 mesh,          /* user_data */
                 (do_volume ? mesh_iter_volume : NULL), mesh_iter_face,
#ifdef P4_TO_P8
                 (do_edge ? mesh_iter_edge : NULL),
#endif /* P4_TO_P8 */
                 (do_corner ? mesh_iter_corner : NULL));

  return mesh;
}

void
p4est_mesh_destroy (p4est_mesh_t * mesh)
{
  int                 level = 0;

  if (mesh->quad_to_tree != NULL) {
    P4EST_FREE (mesh->quad_to_tree);
  }

  if (mesh->quad_level != NULL) {
    for (level = 0; level <= P4EST_QMAXLEVEL; ++level) {
      sc_array_reset (mesh->quad_level + level);
    }
    P4EST_FREE (mesh->quad_level);
  }

  P4EST_FREE (mesh->ghost_to_proc);
  P4EST_FREE (mesh->quad_to_quad);
  P4EST_FREE (mesh->quad_to_face);
  sc_array_destroy (mesh->quad_to_half);

#ifdef P4_TO_P8
  if (mesh->quad_to_edge != NULL) {
    P4EST_FREE (mesh->quad_to_edge);
    sc_array_destroy (mesh->edge_offset);
    sc_array_destroy (mesh->edge_quad);
    sc_array_destroy (mesh->edge_edge);
  }
#endif /* P4_TO_P8 */

  if (mesh->quad_to_corner != NULL) {
    P4EST_FREE (mesh->quad_to_corner);
    sc_array_destroy (mesh->corner_offset);
    sc_array_destroy (mesh->corner_quad);
    sc_array_destroy (mesh->corner_corner);
  }

  P4EST_FREE (mesh);
}

/************************* accessor functions ************************/

p4est_quadrant_t   *
p4est_mesh_get_quadrant (p4est_t * p4est, p4est_mesh_t * mesh,
                         p4est_locidx_t qid)
{
  p4est_tree_t       *tree;
  p4est_quadrant_t   *quad;
  tree =
    (p4est_tree_t *) sc_array_index_int (p4est->trees,
                                         mesh->quad_to_tree[qid]);
  quad =
    (p4est_quadrant_t *) sc_array_index_int (&tree->quadrants,
                                             qid - tree->quadrants_offset);
  return quad;
}

p4est_locidx_t
p4est_mesh_get_neighbors (p4est_t * p4est,
                          p4est_ghost_t * ghost,
                          p4est_mesh_t * mesh,
                          p4est_locidx_t curr_quad_id,
                          p4est_locidx_t direction,
                          sc_array_t * neighboring_quads,
                          sc_array_t * neighboring_encs)
{
  int                 i;
  p4est_locidx_t      lq = mesh->local_num_quadrants;
  p4est_locidx_t      gq = mesh->ghost_num_quadrants;

#ifdef P4EST_DEBUG
  /* Integrity checks: */
  /*  result arrays should be empty, */
  P4EST_ASSERT (neighboring_quads->elem_count == 0);
  P4EST_ASSERT (neighboring_encs->elem_count == 0);
  /*  mesh has to be created, i.e. not NULL, */
  P4EST_ASSERT (mesh != NULL);
  /*  curr_quad_id must be part of the processors quadrants, */
  P4EST_ASSERT (curr_quad_id < lq);

  /*  and direction must be within the allowed range */
  p4est_locidx_t      limit;
  switch (ghost->btype) {
  case P4EST_CONNECT_FACE:
    P4EST_ASSERT (direction >= 0 && direction < P4EST_FACES);
    break;

#ifdef P4_TO_P8
  case P8EST_CONNECT_EDGE:
    limit = P4EST_FACES + P8EST_EDGES;
    P4EST_ASSERT (direction >= 0 && direction < limit);
    break;
#endif /* P4_TO_P8 */

  case P4EST_CONNECT_CORNER:
#ifdef P4_TO_P8
    limit = P4EST_FACES + P8EST_EDGES + P4EST_CHILDREN;
#else /* P4_TO_P8 */
    limit = P4EST_FACES + P4EST_CHILDREN;
#endif /* P4_TO_P8 */
    P4EST_ASSERT (direction >= 0 && direction < limit);
    break;

  default:
    SC_ABORT_NOT_REACHED ();
  }
#endif /* P4EST_DEBUG */

  /* tools for decoding direction */
#ifndef P4_TO_P8
  p4est_locidx_t      lFace, uFace, lCorner, uCorner;
  p4est_locidx_t      convFace = 33;

  lFace = 0;
  uFace = lCorner = P4EST_FACES;
  uCorner = P4EST_FACES + P4EST_CHILDREN;
#else /* !P4_TO_P8 */
  p4est_locidx_t      lFace, uFace, lEdge, uEdge, lCorner, uCorner;
  p4est_locidx_t      convFace, convEdge;
  convFace = 145;
  convEdge = 97;

  lFace = 0;
  uFace = lEdge = P4EST_FACES;
  uEdge = lCorner = lEdge + P8EST_EDGES;
  uCorner = lCorner + P4EST_CHILDREN;
#endif /* P4_TO_P8 */

  p4est_locidx_t      neighbor_idx, neighbor_encoding;
  p4est_quadrant_t  **quad_ins;
  p4est_quadrant_t   *quad;
  int                *enc_ptr;

  /* obtain quadrants and write them into result array */
  /* faces */
  if (lFace <= direction && direction < uFace) {
    neighbor_idx = mesh->quad_to_quad[P4EST_FACES * curr_quad_id + direction];
    neighbor_encoding =
      mesh->quad_to_face[P4EST_FACES * curr_quad_id + direction];

    /* no neighbor present */
    if ((neighbor_idx < 0 || neighbor_idx == curr_quad_id)
        && neighbor_encoding < P4EST_FACES) {
      return 0;
    }

    if (neighbor_encoding < 0) {
      /* half size */
      p4est_locidx_t      quad_idx;
      p4est_locidx_t     *quad_ptr;
      quad_ptr =
        (p4est_locidx_t *) sc_array_index_int (mesh->quad_to_half,
                                               neighbor_idx);
      for (i = 0; i < P4EST_HALF; ++i) {
        quad_idx = quad_ptr[i];
        if (quad_idx < lq) {
          quad_ins = (p4est_quadrant_t **) sc_array_push (neighboring_quads);
          quad = p4est_mesh_get_quadrant (p4est, mesh, quad_idx);
          *quad_ins = quad;

          enc_ptr = (int *) sc_array_push (neighboring_encs);

          /* convert encoding */
          *enc_ptr = neighbor_encoding + convFace;
        }
        else {
          /* neighbor is part of ghost layer */
          quad_ins = (p4est_quadrant_t **) sc_array_push (neighboring_quads);
          quad =
            (p4est_quadrant_t *) sc_array_index_int (&ghost->ghosts,
                                                     quad_idx - lq);
          *quad_ins = quad;
          enc_ptr = (int *) sc_array_push (neighboring_encs);

          /* convert encoding */
          *enc_ptr = -neighbor_encoding - convFace;
        }
      }
    }
    else {
      /* same or double size */
      if (neighbor_idx < lq) {
        /* neighbor is part of quadrants owned by processor */
        quad_ins = (p4est_quadrant_t **) sc_array_push (neighboring_quads);
        quad = p4est_mesh_get_quadrant (p4est, mesh, neighbor_idx);
        *quad_ins = quad;

        enc_ptr = (int *) sc_array_push (neighboring_encs);

        /* convert encoding */
        neighbor_encoding++;
        *enc_ptr = neighbor_encoding;
      }
      /* neighbor is part of ghost layer */
      else {
        quad_ins = (p4est_quadrant_t **) sc_array_push (neighboring_quads);
        quad =
          (p4est_quadrant_t *) sc_array_index_int (&ghost->ghosts,
                                                   neighbor_idx - lq);
        *quad_ins = quad;

        enc_ptr = (int *) sc_array_push (neighboring_encs);

        /* convert encoding */
        neighbor_encoding++;
        *enc_ptr = -neighbor_encoding;
      }
    }
    return 0;
  }

#ifdef P4_TO_P8
  /* edges (3D only) */
  else if (lEdge <= direction && direction < uEdge) {
    neighbor_idx =
      mesh->quad_to_edge[P8EST_EDGES * curr_quad_id + (direction - lEdge)];

    /* no neighbor present */
    if (neighbor_idx < 0 || neighbor_idx == curr_quad_id) {
      return 0;
    }

    if (neighbor_idx < lq) {
      /* same size neighbor, same proc */
      quad_ins = (p4est_quadrant_t **) sc_array_push (neighboring_quads);
      quad = p4est_mesh_get_quadrant (p4est, mesh, neighbor_idx);
      *quad_ins = quad;
      enc_ptr = (int *) sc_array_push (neighboring_encs);

      /* create implicitly saved encoding */
      neighbor_encoding = (direction - lEdge) ^ 3;

      /* convert encoding */
      neighbor_encoding++;
      *enc_ptr = neighbor_encoding;
    }

    if (lq < neighbor_idx && neighbor_idx < lq + gq) {
      /* same size neighbor, ghost layer */
      quad_ins = (p4est_quadrant_t **) sc_array_push (neighboring_quads);
      quad =
        (p4est_quadrant_t *) sc_array_index_int (&ghost->ghosts,
                                                 neighbor_idx - lq);
      *quad_ins = quad;

      enc_ptr = (int *) sc_array_push (neighboring_encs);
      /* create implicitly saved encoding */
      neighbor_encoding = (direction - lEdge) ^ 3;
      /* convert encoding */
      neighbor_encoding++;
      *enc_ptr = -neighbor_encoding;
    }

    if ((lq + gq) <= neighbor_idx) {
      /* anything else */
      /* normalize neighbor index */
      neighbor_idx -= (lq + gq);
      P4EST_ASSERT (0 <= neighbor_idx
                    && neighbor_idx < mesh->local_num_edges);
      p4est_locidx_t      n_adj_quads, offset;

      /* get offset and number of adjacent quads */
      offset = *((p4est_locidx_t *)
                 sc_array_index_int (mesh->edge_offset, neighbor_idx));
      n_adj_quads = *((p4est_locidx_t *)
                      sc_array_index_int (mesh->edge_offset,
                                          neighbor_idx + 1)) - offset;

      p4est_locidx_t      quad_idx;
      for (i = 0; i < n_adj_quads; ++i) {
        quad_idx = *((p4est_locidx_t *)
                     sc_array_index_int (mesh->edge_quad, offset + i));
        neighbor_encoding = *((int8_t *)
                              sc_array_index_int (mesh->edge_edge,
                                                  offset + i));

        if (quad_idx < lq) {
          /* neighbor is part of quadrants owned by processor */
          quad_ins = (p4est_quadrant_t **) sc_array_push (neighboring_quads);
          quad = p4est_mesh_get_quadrant (p4est, mesh, quad_idx);
          *quad_ins = quad;

          enc_ptr = (int *) sc_array_push (neighboring_encs);
          /* convert encoding */
          neighbor_encoding += (neighbor_encoding < 0 ? convEdge : 1);
          *enc_ptr = neighbor_encoding;
        }
        else {
          /* neighbor is part of ghost layer */
          quad_ins = (p4est_quadrant_t **) sc_array_push (neighboring_quads);
          quad =
            (p4est_quadrant_t *) sc_array_index_int (&ghost->ghosts,
                                                     quad_idx - lq);
          *quad_ins = quad;

          enc_ptr = (int *) sc_array_push (neighboring_encs);
          /* convert encoding */
          neighbor_encoding += (neighbor_encoding < 0 ? convEdge : 1);
          *enc_ptr = -neighbor_encoding;
        }
      }
    }
    return 0;
  }
#endif /* P4_TO_P8 */
  /* corners */
  else if (lCorner <= direction && direction < uCorner) {
    p4est_locidx_t      encHelper;
#ifdef P4_TO_P8
    encHelper = 7;
#else /* P4_TO_P8 */
    encHelper = 3;
#endif /* P4_TO_P8 */

    neighbor_idx =
      mesh->quad_to_corner[P4EST_CHILDREN * curr_quad_id +
                           (direction - lCorner)];

    /* no neighbor present */
    if (neighbor_idx < 0 || neighbor_idx == curr_quad_id) {
      return 0;
    }

    if (neighbor_idx < lq) {
      /* same size neighbor, same proc */
      quad_ins = (p4est_quadrant_t **) sc_array_push (neighboring_quads);
      quad = p4est_mesh_get_quadrant (p4est, mesh, neighbor_idx);
      *quad_ins = quad;

      enc_ptr = (int *) sc_array_push (neighboring_encs);
      /* create implicitly saved encoding */
      neighbor_encoding = (direction - lCorner) ^ encHelper;
      /* convert encoding */
      neighbor_encoding++;
      *enc_ptr = neighbor_encoding;
    }

    else if (lq < neighbor_idx && neighbor_idx < lq + gq) {
      /* same size neighbor, ghost layer */
      quad_ins = (p4est_quadrant_t **) sc_array_push (neighboring_quads);
      quad =
        (p4est_quadrant_t *) sc_array_index_int (&ghost->ghosts,
                                                 neighbor_idx - lq);
      *quad_ins = quad;

      enc_ptr = (int *) sc_array_push (neighboring_encs);
      /* create implicitly saved encoding */
      neighbor_encoding = (direction - lCorner) ^ encHelper;
      /* convert encoding */
      neighbor_encoding++;
      *enc_ptr = -neighbor_encoding;
    }

    else if ((lq + gq) <= neighbor_idx) {
      /* anything else */
      /* normalize neighbor index */
      neighbor_idx -= (lq + gq);
      P4EST_ASSERT (0 <= neighbor_idx
                    && neighbor_idx < mesh->local_num_corners);
      p4est_locidx_t      n_adj_quads, offset;

      /* get offset and number of adjacent quads */
      offset = *((p4est_locidx_t *)
                 sc_array_index_int (mesh->corner_offset, neighbor_idx));
      n_adj_quads = *((p4est_locidx_t *)
                      sc_array_index_int (mesh->corner_offset,
                                          neighbor_idx + 1)) - offset;

      p4est_locidx_t      quad_idx;
      for (i = 0; i < n_adj_quads; ++i) {
        quad_idx = *((p4est_locidx_t *)
                     sc_array_index_int (mesh->corner_quad, offset + i));
        neighbor_encoding = *((int8_t *)
                              sc_array_index_int (mesh->corner_corner,
                                                  offset + i));

        if (quad_idx < lq) {
          /* neighbor is part of quadrants owned by processor */
          quad_ins = (p4est_quadrant_t **) sc_array_push (neighboring_quads);
          quad = p4est_mesh_get_quadrant (p4est, mesh, quad_idx);
          *quad_ins = quad;

          enc_ptr = (int *) sc_array_push (neighboring_encs);
          /* convert encoding */
          neighbor_encoding++;
          *enc_ptr = neighbor_encoding;
        }
        else {
          /* neighbor is part of ghost layer */
          quad_ins = (p4est_quadrant_t **) sc_array_push (neighboring_quads);
          quad =
            (p4est_quadrant_t *) sc_array_index_int (&ghost->ghosts,
                                                     quad_idx - lq);
          *quad_ins = quad;

          enc_ptr = (int *) sc_array_push (neighboring_encs);
          /* convert encoding */
          neighbor_encoding++;
          *enc_ptr = -neighbor_encoding;
        }
      }
    }
    return 0;
  }

  /* must not occur */
  else {
    SC_ABORT_NOT_REACHED ();
  }

  /* this is not reached anyways */
  return -420;
}

p4est_quadrant_t   *
p4est_mesh_quadrant_cumulative (p4est_t * p4est, p4est_mesh_t * mesh,
                                p4est_locidx_t cumulative_id,
                                p4est_topidx_t * pwhich_tree,
                                p4est_locidx_t * pquadrant_id)
{
  p4est_topidx_t      which_tree;
  p4est_locidx_t      quadrant_id;
  p4est_quadrant_t   *quadrant;
#ifdef P4EST_ENABLE_DEBUG
  p4est_topidx_t      dwhich_tree;
  p4est_locidx_t      dquadrant_id;
  p4est_quadrant_t   *dquadrant;
#endif
  p4est_tree_t       *tree;

  P4EST_ASSERT (p4est != NULL);
  P4EST_ASSERT (mesh != NULL);
  P4EST_ASSERT (p4est->local_num_quadrants == mesh->local_num_quadrants);

  P4EST_ASSERT (0 <= cumulative_id &&
                cumulative_id < mesh->local_num_quadrants);
  P4EST_ASSERT (NULL == pwhich_tree || *pwhich_tree == -1 ||
                (0 <= *pwhich_tree &&
                 *pwhich_tree < p4est->connectivity->num_trees));

  if (mesh->quad_to_tree != NULL) {
    /* in this case we can do an O(1) lookup */
    which_tree = mesh->quad_to_tree[cumulative_id];
    if (pwhich_tree != NULL) {
#ifdef P4EST_ENABLE_DEBUG
      dwhich_tree = *pwhich_tree;
#endif
      *pwhich_tree = which_tree;
    }
    tree = p4est_tree_array_index (p4est->trees, which_tree);
    P4EST_ASSERT (tree->quadrants_offset <= cumulative_id);
    quadrant_id = cumulative_id - tree->quadrants_offset;
    P4EST_ASSERT (quadrant_id < (p4est_locidx_t) tree->quadrants.elem_count);
    if (pquadrant_id != NULL) {
      *pquadrant_id = quadrant_id;
    }
    quadrant = p4est_quadrant_array_index (&tree->quadrants, quadrant_id);

#ifdef P4EST_ENABLE_DEBUG
    /* we use the more expensive binary search for debugging */
    dquadrant = p4est_find_quadrant_cumulative (p4est, cumulative_id,
                                                &dwhich_tree, &dquadrant_id);
    P4EST_ASSERT (dwhich_tree == which_tree);
    P4EST_ASSERT (dquadrant_id == quadrant_id);
    P4EST_ASSERT (dquadrant == quadrant);
#endif
  }
  else {
    /* we do not have the O(1) lookup table and need to binary search */
    quadrant = p4est_find_quadrant_cumulative (p4est, cumulative_id,
                                               pwhich_tree, pquadrant_id);
  }

  return quadrant;
}

void
p4est_mesh_face_neighbor_init2 (p4est_mesh_face_neighbor_t * mfn,
                                p4est_t * p4est, p4est_ghost_t * ghost,
                                p4est_mesh_t * mesh,
                                p4est_topidx_t which_tree,
                                p4est_locidx_t quadrant_id)
{
  p4est_tree_t       *tree;

  mfn->p4est = p4est;
  mfn->ghost = ghost;
  mfn->mesh = mesh;

  P4EST_ASSERT (0 <= which_tree &&
                which_tree < p4est->connectivity->num_trees);
  mfn->which_tree = which_tree;
  tree = p4est_tree_array_index (p4est->trees, which_tree);

  P4EST_ASSERT (0 <= quadrant_id &&
                (size_t) quadrant_id < tree->quadrants.elem_count);
  mfn->quadrant_id = quadrant_id;
  mfn->quadrant_code = P4EST_FACES * (tree->quadrants_offset + quadrant_id);

  mfn->face = 0;
  mfn->subface = 0;
  mfn->current_qtq = -1;
}

void
p4est_mesh_face_neighbor_init (p4est_mesh_face_neighbor_t * mfn,
                               p4est_t * p4est, p4est_ghost_t * ghost,
                               p4est_mesh_t * mesh, p4est_topidx_t which_tree,
                               p4est_quadrant_t * quadrant)
{
  p4est_locidx_t      quadrant_id;
  p4est_tree_t       *tree;

  mfn->p4est = p4est;
  mfn->ghost = ghost;
  mfn->mesh = mesh;

  P4EST_ASSERT (0 <= which_tree &&
                which_tree < p4est->connectivity->num_trees);
  mfn->which_tree = which_tree;
  tree = p4est_tree_array_index (p4est->trees, which_tree);

  quadrant_id =
    (p4est_locidx_t) sc_array_position (&tree->quadrants, quadrant);
  mfn->quadrant_id = quadrant_id;
  mfn->quadrant_code = P4EST_FACES * (tree->quadrants_offset + quadrant_id);

  mfn->face = 0;
  mfn->subface = 0;
  mfn->current_qtq = -1;
}

p4est_quadrant_t   *
p4est_mesh_face_neighbor_next (p4est_mesh_face_neighbor_t * mfn,
                               p4est_topidx_t * ntree, p4est_locidx_t * nquad,
                               int *nface, int *nrank)
{
  int                 qtf;
  p4est_topidx_t      which_tree;
  p4est_locidx_t      qtq, quadfacecode;
  p4est_locidx_t      lnq, *halfs;
#ifdef P4EST_ENABLE_DEBUG
  p4est_locidx_t      ngh;
#endif
  p4est_quadrant_t   *q;

  /* We have already processed the last quadrant */
  if (mfn->face == P4EST_FACES) {
    mfn->current_qtq = -1;
    P4EST_ASSERT (mfn->subface == 0);
    return NULL;
  }

  /* Make sure we have a valid quadrant face and iterator */
  lnq = mfn->mesh->local_num_quadrants;
#ifdef P4EST_ENABLE_DEBUG
  ngh = mfn->mesh->ghost_num_quadrants;
#endif
  P4EST_ASSERT (mfn->face >= 0 && mfn->face < P4EST_FACES);
  P4EST_ASSERT (mfn->subface >= 0 && mfn->subface < P4EST_HALF);
  P4EST_ASSERT (mfn->p4est->local_num_quadrants == lnq);
  P4EST_ASSERT (mfn->ghost->ghosts.elem_count == (size_t) ngh);

  /* Retrieve face and quadrant codes */
  quadfacecode = mfn->quadrant_code + (p4est_locidx_t) mfn->face;
  qtq = mfn->mesh->quad_to_quad[quadfacecode];
  qtf = (int) mfn->mesh->quad_to_face[quadfacecode];
  if (qtf >= 0) {
    /* Neighbor is same or double size */
    ;

    /* Advance to next quadrant */
    ++mfn->face;
  }
  else {
    /* Neighbors across this face are half size */
    P4EST_ASSERT (qtq >= 0);
    halfs = (p4est_locidx_t *) sc_array_index (mfn->mesh->quad_to_half,
                                               (size_t) qtq);
    qtq = halfs[mfn->subface];

    /* Advance to next quadrant */
    if (++mfn->subface == P4EST_HALF) {
      mfn->subface = 0;
      ++mfn->face;
    }
  }

  mfn->current_qtq = qtq;
  /* From here on face and subface have advanced and can no longer be used */
  P4EST_ASSERT (qtq >= 0);
  if (qtq < lnq) {
    /* Local quadrant */
    which_tree = mfn->which_tree;
    q = p4est_mesh_quadrant_cumulative (mfn->p4est, mfn->mesh,
                                        qtq, &which_tree, nquad);
    if (ntree != NULL) {
      *ntree = which_tree;
    }
    if (nrank != NULL) {
      *nrank = mfn->p4est->mpirank;
    }
  }
  else {
    /* Ghost quadrant */
    qtq -= lnq;
    P4EST_ASSERT (qtq < ngh);
    q = p4est_quadrant_array_index (&mfn->ghost->ghosts, (size_t) qtq);
    if (ntree != NULL) {
      *ntree = q->p.piggy3.which_tree;
    }
    if (nquad != NULL) {
      *nquad = qtq;             /* number of ghost in the ghost layer */
    }
    if (nrank != NULL) {
      *nrank = mfn->mesh->ghost_to_proc[qtq];
    }
  }
  if (nface != NULL) {
    *nface = qtf;
  }

  return q;
}

void               *
p4est_mesh_face_neighbor_data (p4est_mesh_face_neighbor_t * mfn,
                               void *ghost_data)
{
  p4est_locidx_t      qtq = mfn->current_qtq;
  p4est_locidx_t      lnq = mfn->mesh->local_num_quadrants;
  size_t              data_size = mfn->p4est->data_size;

  P4EST_ASSERT (qtq >= 0);

  if (qtq < lnq) {
    p4est_topidx_t      which_tree;
    p4est_quadrant_t   *q;
    /* Local quadrant */
    which_tree = mfn->which_tree;
    q = p4est_mesh_quadrant_cumulative (mfn->p4est, mfn->mesh,
                                        qtq, &which_tree, NULL);
    return q->p.user_data;
  }
  else {
    qtq -= lnq;
    return (void *) ((char *) ghost_data + data_size * qtq);
  }
}
