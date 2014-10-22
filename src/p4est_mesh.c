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
#else
#include <p8est_bits.h>
#include <p8est_extended.h>
#include <p8est_iterate.h>
#include <p8est_mesh.h>
#endif

/** For a quadrant that touches a tree face with a corner inside the face,
 * get the number of the touching face.
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
#endif
  SC_ABORT_NOT_REACHED ();
}

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

static void
mesh_iter_corner (p4est_iter_corner_info_t * info, void *user_data)
{
  int                 i, j;
  int                 f1, f2, code, orientation;
  int                 fc1, fc2, diagonal;
#ifdef P4_TO_P8
  int                 pref, pset;
#endif
  int                 visited[P4EST_CHILDREN];
  int8_t             *pccorner;
  size_t              cz, zz;
  p4est_locidx_t      qoffset, qid1, qid2;
  p4est_locidx_t      cornerid_offset, cornerid;
  p4est_locidx_t     *pcquad;
  p4est_mesh_t       *mesh = (p4est_mesh_t *) user_data;
  p4est_iter_corner_side_t *side1, *side2;
  p4est_tree_t       *tree1, *tree2;

  /* Check the case when the corner does not involve neighbors */
  cz = info->sides.elem_count;
  P4EST_ASSERT (cz > 0);
  P4EST_ASSERT (info->tree_boundary || cz == P4EST_CHILDREN);
  if (cz == 1) {
    return;
  }

  if (info->tree_boundary == P4EST_CONNECT_FACE) {
    /* This corner is inside an inter-tree face */
    if (cz == P4EST_HALF) {
      /* This is a tree face boundary, no corner neighbors exist */
      return;
    }
    P4EST_ASSERT (cz == P4EST_CHILDREN);
    cornerid_offset = mesh->local_num_quadrants + mesh->ghost_num_quadrants;

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
          tree1 = p4est_tree_array_index (info->p4est->trees, side1->treeid);
          qid1 = side1->quadid + (side1->is_ghost ? mesh->local_num_quadrants
                                  : tree1->quadrants_offset);
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
        code = info->p4est->connectivity->tree_to_face[P4EST_FACES *
                                                       side1->treeid + f1];
        orientation = code / P4EST_FACES;
        P4EST_ASSERT (f2 == code % P4EST_FACES);
#ifdef P4_TO_P8
        pref = p8est_face_permutation_refs[f1][f2];
        pset = p8est_face_permutation_sets[pref][orientation];
        diagonal = (p8est_face_permutations[pset][fc1] ^ fc2) == 3;
#else
        diagonal = (fc1 ^ fc2) != orientation;
#endif
        if (!diagonal) {
          side2 = NULL;
          continue;
        }

        /* We have found a diagonally opposite second side */
        tree2 = p4est_tree_array_index (info->p4est->trees, side2->treeid);
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
    /* Tree corner neighbors across an edge are not implemented: set to -2 */
    for (zz = 0; zz < cz; ++zz) {
      side1 = (p4est_iter_corner_side_t *) sc_array_index (&info->sides, zz);
      if (!side1->is_ghost) {
        tree1 = p4est_tree_array_index (info->p4est->trees, side1->treeid);
        qid1 = side1->quadid + tree1->quadrants_offset;
        P4EST_ASSERT (0 <= qid1 && qid1 < mesh->local_num_quadrants);
        P4EST_ASSERT (mesh->quad_to_corner[P4EST_CHILDREN * qid1 +
                                           side1->corner] == -1);
        mesh->quad_to_corner[P4EST_CHILDREN * qid1 + side1->corner] = -2;
      }
    }
    return;
  }
#endif

  if (info->tree_boundary == P4EST_CONNECT_CORNER) {
    /* True tree corner neighbors are not implemented yet: set to -2 */
    for (zz = 0; zz < cz; ++zz) {
      side1 = (p4est_iter_corner_side_t *) sc_array_index (&info->sides, zz);
      if (!side1->is_ghost) {
        tree1 = p4est_tree_array_index (info->p4est->trees, side1->treeid);
        qid1 = side1->quadid + tree1->quadrants_offset;
        P4EST_ASSERT (0 <= qid1 && qid1 < mesh->local_num_quadrants);
        P4EST_ASSERT (mesh->quad_to_corner[P4EST_CHILDREN * qid1 +
                                           side1->corner] == -1);
        mesh->quad_to_corner[P4EST_CHILDREN * qid1 + side1->corner] = -2;
      }
    }
    return;
  }

  /* Process a corner inside the tree in pairs of diagonal neighbors */
  P4EST_ASSERT (!info->tree_boundary);
  side1 = (p4est_iter_corner_side_t *) sc_array_index (&info->sides, 0);
  tree1 = p4est_tree_array_index (info->p4est->trees, side1->treeid);
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
        if (!side2->is.hanging.is_ghost[h]) {
          tree = p4est_tree_array_index (info->p4est->trees, side2->treeid);
          jls[h] = side2->is.hanging.quadid[h] + tree->quadrants_offset;
          P4EST_ASSERT (0 <= jls[h] && jls[h] < mesh->local_num_quadrants);
        }
        else {
          P4EST_ASSERT (side2->is.hanging.quad[h] != NULL);
          P4EST_ASSERT (side2->is.hanging.quadid[h] >= 0);
          jls[h] = mesh->local_num_quadrants + side2->is.hanging.quadid[h];
        }
      }

      /* encode quadrant neighborhood */
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
        if (!side2->is.hanging.is_ghost[h]) {
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
  p4est_locidx_t     *quadid;
  int                 level = info->quad->level;

  /* We could use a static quadrant counter, but that gets uglier */
  tree = p4est_tree_array_index (info->p4est->trees, info->treeid);
  P4EST_ASSERT (0 <= info->quadid &&
                info->quadid < (p4est_locidx_t) tree->quadrants.elem_count);

  if (mesh->quad_to_tree != NULL) {
    mesh->quad_to_tree[tree->quadrants_offset + info->quadid] = info->treeid;
  }

  if (mesh->quad_level != NULL) {
    quadid = (p4est_locidx_t *) sc_array_push (mesh->quad_level + level);
    *quadid = tree->quadrants_offset + info->quadid;
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
  int                 do_volume = 0;
  int                 rank;
  p4est_locidx_t      lq, ng;
  p4est_locidx_t      jl;
  p4est_mesh_t       *mesh;

  P4EST_ASSERT (p4est_is_balanced (p4est, P4EST_CONNECT_FULL));

  mesh = P4EST_ALLOC_ZERO (p4est_mesh_t, 1);

  lq = mesh->local_num_quadrants = p4est->local_num_quadrants;
  ng = mesh->ghost_num_quadrants = (p4est_locidx_t) ghost->ghosts.elem_count;

  if (btype == P4EST_CONNECT_FULL) {
    do_corner = 1;
  }
  do_volume = (compute_tree_index || compute_level_lists ? 1 : 0);

  /* Optional map of tree index for each quadrant */
  if (compute_tree_index) {
    mesh->quad_to_tree = P4EST_ALLOC (p4est_topidx_t, lq);
  }

  mesh->ghost_to_proc = P4EST_ALLOC (int, ng);
  mesh->quad_to_quad = P4EST_ALLOC (p4est_locidx_t, P4EST_FACES * lq);
  mesh->quad_to_face = P4EST_ALLOC (int8_t, P4EST_FACES * lq);
  mesh->quad_to_half = sc_array_new (P4EST_HALF * sizeof (p4est_locidx_t));

  /* Optional per-level lists of quadrants */
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

  if (do_corner) {
    /* Initialize corner information to a consistent state */
    mesh->quad_to_corner = P4EST_ALLOC (p4est_locidx_t, P4EST_CHILDREN * lq);
    memset (mesh->quad_to_corner, -1, P4EST_CHILDREN * lq *
            sizeof (p4est_locidx_t));

    mesh->corner_offset = sc_array_new (sizeof (p4est_locidx_t));
    *(p4est_locidx_t *) sc_array_push (mesh->corner_offset) = 0;

    mesh->corner_quad = sc_array_new (sizeof (p4est_locidx_t));
    mesh->corner_corner = sc_array_new (sizeof (int8_t));
  }

  /* Call the forest iterator to collect face connectivity */
  p4est_iterate (p4est, ghost, mesh,
                 (do_volume ? mesh_iter_volume : NULL), mesh_iter_face,
#ifdef P4_TO_P8
                 NULL,
#endif
                 do_corner ? mesh_iter_corner : NULL);

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

  if (mesh->quad_to_corner != NULL) {
    P4EST_FREE (mesh->quad_to_corner);
    sc_array_destroy (mesh->corner_offset);
    sc_array_destroy (mesh->corner_quad);
    sc_array_destroy (mesh->corner_corner);
  }

  P4EST_FREE (mesh);
}

p4est_quadrant_t   *
p4est_mesh_quadrant_cumulative (p4est_t * p4est, p4est_locidx_t cumulative_id,
                                p4est_topidx_t * which_tree,
                                p4est_locidx_t * quadrant_id)
{
  int                 the_quadrant_id;
  p4est_topidx_t      low_tree, high_tree, guess_tree;
  p4est_tree_t       *tree;

  P4EST_ASSERT (0 <= cumulative_id &&
                cumulative_id < p4est->local_num_quadrants);

  low_tree = p4est->first_local_tree;
  high_tree = p4est->last_local_tree;
  if (which_tree != NULL && *which_tree != -1) {
    guess_tree = *which_tree;
  }
  else {
    guess_tree = (low_tree + high_tree) / 2;
  }
  for (;;) {
    P4EST_ASSERT (p4est->first_local_tree <= low_tree);
    P4EST_ASSERT (high_tree <= p4est->last_local_tree);
    P4EST_ASSERT (low_tree <= guess_tree && guess_tree <= high_tree);

    tree = p4est_tree_array_index (p4est->trees, guess_tree);
    if (cumulative_id < tree->quadrants_offset) {
      high_tree = guess_tree - 1;
    }
    else if (cumulative_id >= tree->quadrants_offset +
             (p4est_locidx_t) tree->quadrants.elem_count) {
      low_tree = guess_tree + 1;
    }
    else {
      the_quadrant_id = cumulative_id - tree->quadrants_offset;
      P4EST_ASSERT (0 <= the_quadrant_id);

      if (which_tree != NULL) {
        *which_tree = guess_tree;
      }
      if (quadrant_id != NULL) {
        *quadrant_id = the_quadrant_id;
      }
      return p4est_quadrant_array_index (&tree->quadrants,
                                         (size_t) the_quadrant_id);
    }
    guess_tree = (low_tree + high_tree) / 2;
  }
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
    q = p4est_mesh_quadrant_cumulative (mfn->p4est, qtq, &which_tree, nquad);
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
    q = p4est_mesh_quadrant_cumulative (mfn->p4est, qtq, &which_tree, NULL);
    return q->p.user_data;
  }
  else {
    qtq -= lnq;
    return (void *) ((char *) ghost_data + data_size * qtq);
  }
}
