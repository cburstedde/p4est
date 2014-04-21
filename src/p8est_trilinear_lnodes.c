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

#include <p4est_to_p8est.h>
#include <p8est_bits.h>
#include <p8est_trilinear.h>

static void
p8est_mesh_fcode_to_ntype (p8est_lnodes_code_t face_code, int8_t node_type[])
{
  int                 i, c;
  int                 edge[12];
  int                 face[6];

  memset (node_type, 0, 8);

  if (p8est_lnodes_decode (face_code, face, edge)) {
    for (i = 0; i < 6; i++) {
      if (face[i] >= 0) {
        c = p8est_face_corners[i][3 - face[i]];
        node_type[c] = 1;
      }
    }
    for (i = 0; i < 12; i++) {
      if (edge[i] >= 0 && edge[i] < 4) {
        c = p8est_edge_corners[i][1 - (edge[i] % 2)];
        node_type[c] = 2;
      }
    }
  }
}

static void
p8est_mesh_points (point_t * points, trilinear_mesh_pid_t * pids,
                   p4est_locidx_t * local_nodes,
                   p4est_locidx_t num_local_nodes, p4est_lnodes_t * nodes,
                   p4est_t * p4est)
{
  p4est_locidx_t      nin = nodes->num_local_nodes;
  p4est_locidx_t      owned_offset = 0;
  p4est_locidx_t      owned_count = nodes->owned_count;
  p4est_locidx_t      elid, nid;
  p4est_locidx_t      lz;
  int                 qid;
  int                 i, j, k, f, e;
  int                 nf, nf2, ne, nc, cid;
  int                 o, ref, set;
  int8_t              node_type[8];
  p8est_lnodes_code_t *face_code = nodes->face_code;
  int                 ownerc;
  p4est_topidx_t      t, nt, ownert;
  p4est_topidx_t      ft = p4est->first_local_tree;
  p4est_topidx_t      lt = p4est->last_local_tree;
  sc_array_t         *trees = p4est->trees;
  p4est_tree_t       *tree;
  size_t              count, zz;
  p4est_quadrant_t   *q, p, r, s, ownerq;
  p4est_qcoord_t      shift;
  p8est_connectivity_t *conn = p4est->connectivity;
  p4est_topidx_t     *ttt = conn->tree_to_tree;
  int8_t             *ttf = conn->tree_to_face;
  int                 ftr[P8EST_FTRANSFORM];
  p8est_edge_info_t   ei;
  sc_array_t         *eta = &ei.edge_transforms;
  p8est_edge_transform_t *et;
  size_t              etree;
  p8est_corner_info_t ci;
  sc_array_t         *cta = &ci.corner_transforms;
  p8est_corner_transform_t *ct;
  size_t              ctree;

  P4EST_ASSERT (nodes->degree == 1 && nodes->vnodes == 8);

  memset (pids, -1, sizeof (trilinear_mesh_pid_t) * num_local_nodes);
  for (nid = 0, elid = 0, t = ft; t <= lt; t++) {
    tree = p4est_tree_array_index (trees, t);
    count = tree->quadrants.elem_count;
    q = p4est_quadrant_array_index (&tree->quadrants, 0);
    for (zz = 0; zz < count; zz++, q++, elid++) {
      p8est_mesh_fcode_to_ntype (face_code[elid], node_type);
      /* for every corner of every quadrant */
      for (i = 0; i < 8; i++, nid++) {
        for (k = 0; k < 2; k++) {
          if (k == 0) {
            lz = local_nodes[nid];
            /* either a node is independent or hanging */
            P4EST_ASSERT ((node_type[i] == 0) ^ (lz >= nin));
          }
          /* if a corner is hanging, the point values of both the hanging node
           * and the independent node must be computed */
          else if (local_nodes[nid] >= nin) {
            lz = nodes->element_nodes[nid];
            P4EST_ASSERT (lz < nin);
            /* if the independent node is local, it must (by the ownership
             * rules) be touches by a local quadrant, so its point values will
             * be set by that quadrant: the quadrant that owns that hanging node
             * needn't bother */
            if (lz >= owned_offset && lz - owned_offset < owned_count) {
              continue;
            }
          }
          else {
            continue;
          }
          /* pids are set whenever points are: if pid is set, nothing to do */
          if (pids[lz] >= 0) {
            continue;
          }
          /* if the node is local or hanging, its location is the current
           * corner of the current quadrant */
          if ((lz >= owned_offset && lz - owned_offset < owned_count) ||
              lz >= nin) {
            shift = P8EST_QUADRANT_LEN (q->level);
            points[lz].x = (tick_t) (q->x + ((i & 1) ? shift : 0));
            points[lz].y = (tick_t) (q->y + (((i >> 1) & 1) ? shift : 0));
            points[lz].z = (tick_t) (q->z + ((i >> 2) ? shift : 0));
            pids[lz] = (trilinear_mesh_pid_t) t;
          }
          else {
            if (k == 0) {
              p8est_quadrant_corner_descendant (q, &p, i, P8EST_QMAXLEVEL);
            }
            /* if k == 1, then that means that the current node is indpendent,
             * but is not touched by the current quadrant: rather, it is
             * touched by the parent of the current quadrant */
            else {
              p8est_quadrant_parent (q, &r);
              p8est_quadrant_corner_descendant (&r, &p, i, P8EST_QMAXLEVEL);
            }
            p4est_quadrant_corner_neighbor (&p, i, &r);
            /* if the corner in question is in the interior of the tree, its
             * point values are unambiguous */
            if (p8est_quadrant_is_inside_root (&r)) {
              shift = P8EST_QUADRANT_LEN (P8EST_QMAXLEVEL);
              points[lz].x = (tick_t) (p.x + ((i & 1) ? shift : 0));
              points[lz].y = (tick_t) (p.y + (((i >> 1) & 1) ? shift : 0));
              points[lz].z = (tick_t) (p.z + ((i >> 2) ? shift : 0));
              pids[lz] = (trilinear_mesh_pid_t) t;
            }
            else {
              ownert = t;
              ownerq = p;
              ownerc = i;
              /* if q->level == P8EST_QMAXLEVEL,
               * qid and i may not be the same */
              qid = p8est_quadrant_child_id (&p);
              for (j = 0; j < 3; j++) {
                f = p8est_corner_faces[i][j];
                p8est_quadrant_face_neighbor (&p, f, &r);
                if (p8est_quadrant_is_outside_face (&r)) {
                  nt = ttt[t * P4EST_FACES + f];
                  if (ownert < nt) {
                    continue;
                  }
                  nf = (int) ttf[t * P4EST_FACES + f];
                  o = nf / P8EST_FACES;
                  nf %= P8EST_FACES;
                  if (nt == t && nf == f) {
                    continue;
                  }
                  cid = p8est_corner_face_corners[i][f];
                  ref = p8est_face_permutation_refs[f][nf];
                  set = p8est_face_permutation_sets[ref][o];
                  cid = p8est_face_permutations[set][cid];
                  nc = p8est_face_corners[nf][cid];
                  nf2 = (int) p8est_find_face_transform (conn, t, f, ftr);
                  P4EST_ASSERT (nf == nf2);
                  p8est_quadrant_transform_face (&r, &s, ftr);
                  if (nt < ownert || ((nt == ownert) &&
                                      p8est_quadrant_compare (&s,
                                                              &ownerq) < 0)) {
                    ownerq = s;
                    ownert = nt;
                    ownerc = nc;
                  }
                }
              }
              for (j = 0; j < 3; j++) {
                e = p8est_corner_edges[i][j];
                p8est_quadrant_edge_neighbor (&p, e, &r);
                if (p8est_quadrant_is_outside_edge (&r)) {
                  sc_array_init (eta, sizeof (p8est_edge_transform_t));
                  p8est_find_edge_transform (conn, t, e, &ei);
                  for (etree = 0; etree < eta->elem_count; etree++) {
                    et = (p8est_edge_transform_t *) sc_array_index (eta,
                                                                    etree);
                    nt = et->ntree;
                    if (ownert < nt) {
                      continue;
                    }
                    p8est_quadrant_transform_edge (&r, &s, &ei, et, 1);
                    ne = (int) et->nedge;
                    if (nt < ownert || ((nt == ownert) &&
                                        p8est_quadrant_compare (&s,
                                                                &ownerq) <
                                        0)) {
                      ownerq = s;
                      ownert = nt;
                      ownerc = p8est_quadrant_child_id (&s);
                      if (qid != i) {
                        P4EST_ASSERT ((p8est_edge_corners[e][0] == i
                                       && p8est_edge_corners[e][1] == qid) ||
                                      (p8est_edge_corners[e][0] == qid
                                       && p8est_edge_corners[e][1] == i));
                        P4EST_ASSERT (p8est_edge_corners[ne][0] == ownerc ||
                                      p8est_edge_corners[ne][1] == ownerc);
                        ownerc = (p8est_edge_corners[ne][0] == ownerc) ?
                          p8est_edge_corners[ne][1] :
                          p8est_edge_corners[ne][0];
                      }
                    }
                  }
                  sc_array_reset (eta);
                }
              }
              p8est_quadrant_corner_neighbor (&p, i, &r);
              if (p8est_quadrant_is_outside_corner (&r)) {
                sc_array_init (cta, sizeof (p8est_corner_transform_t));
                p8est_find_corner_transform (conn, t, i, &ci);
                for (ctree = 0; ctree < cta->elem_count; ctree++) {
                  ct = (p8est_corner_transform_t *) sc_array_index (cta,
                                                                    ctree);
                  nt = ct->ntree;
                  if (ownert < nt) {
                    continue;
                  }
                  p8est_quadrant_transform_corner (&r, (int) ct->ncorner, 1);
                  if (nt < ownert || ((nt == ownert) &&
                                      p8est_quadrant_compare (&r,
                                                              &ownerq) < 0)) {
                    ownerq = r;
                    ownert = nt;
                    ownerc = (int) ct->ncorner;
                  }
                }
                sc_array_reset (cta);
              }
              P4EST_ASSERT (ownerq.level == P8EST_QMAXLEVEL);
              shift = P8EST_QUADRANT_LEN (P8EST_QMAXLEVEL);
              points[lz].x = (tick_t) (ownerq.x + ((ownerc & 1) ? shift : 0));
              points[lz].y = (tick_t) (ownerq.y +
                                       (((ownerc >> 1) & 1) ? shift : 0));
              points[lz].z =
                (tick_t) (ownerq.z + ((ownerc >> 2) ? shift : 0));
              pids[lz] = (trilinear_mesh_pid_t) ownert;
            }
          }
        }
      }
    }
  }
}

static unsigned
p8est_mesh_indep_hash_fn (const void *v, const void *u)
{
  const p4est_locidx_t *indep = (p4est_locidx_t *) v;
  uint32_t            a, b, c;

  a = (uint32_t) indep[0];
  b = (uint32_t) indep[1];
  c = (uint32_t) indep[2];
  sc_hash_mix (a, b, c);
  a += (uint32_t) indep[3];
  sc_hash_final (a, b, c);

  return (unsigned) c;
}

static int
p8est_mesh_indep_equal_fn (const void *v1, const void *v2, const void *u)
{
  const p4est_locidx_t *a = (p4est_locidx_t *) v1;
  const p4est_locidx_t *b = (p4est_locidx_t *) v2;
  return a[0] == b[0] && a[1] == b[1] && a[2] == b[2] && a[3] == b[3];
}

static              p4est_locidx_t
p8est_mesh_dcount (p4est_lnodes_t * nodes, p4est_locidx_t ** Local_nodes,
                   sc_hash_array_t ** Hash)
{
  p4est_locidx_t      nel = nodes->num_local_elements;
  p4est_locidx_t      nin = nodes->num_local_nodes;
  p4est_locidx_t      nln = nel * 8;
  sc_array_t          local_nodes;
  p4est_locidx_t      elid, nid;
  int                 i, j, f, e;
  p8est_lnodes_code_t *face_code = nodes->face_code;
  int8_t              ntype[8];
  int                 faces[6], edges[12];
  p4est_locidx_t      indep[4], *r;
  sc_hash_array_t    *hanging;
  size_t              position;
  p4est_locidx_t      num_nodes = nin;
  p4est_locidx_t     *lp;

  sc_array_init (&local_nodes, sizeof (p4est_locidx_t));
  sc_array_resize (&local_nodes, nln);
  memcpy (local_nodes.array, nodes->element_nodes,
          nln * sizeof (p4est_locidx_t));

  *Hash = hanging = sc_hash_array_new (4 * sizeof (p4est_locidx_t),
                                       p8est_mesh_indep_hash_fn,
                                       p8est_mesh_indep_equal_fn, NULL);

  for (nid = 0, elid = 0; elid < nel; elid++) {
    if (p8est_lnodes_decode (face_code[elid], faces, edges)) {
      p8est_mesh_fcode_to_ntype (face_code[elid], ntype);
      for (i = 0; i < 8; i++, nid++) {
        if (ntype[i] == 1) {
          for (j = 0; j < 3; j++) {
            f = p8est_corner_faces[i][j];
            if (faces[f] >= 0) {
              P4EST_ASSERT ((faces[f] ^ p8est_corner_face_corners[i][f]) ==
                            3);
              break;
            }
          }
          P4EST_ASSERT (j < 3);
          for (j = 0; j < 4; j++) {
            indep[j] =
              nodes->element_nodes[elid * 8 + p8est_face_corners[f][j]];
          }
          qsort (indep, 4, sizeof (p4est_locidx_t), p4est_locidx_compare);
          lp = (p4est_locidx_t *) sc_array_index (&local_nodes, nid);
          r = (p4est_locidx_t *) sc_hash_array_insert_unique (hanging, indep,
                                                              &position);
          if (r != NULL) {
            r[0] = indep[0], r[1] = indep[1], r[2] = indep[2], r[3] =
              indep[3];
            *lp = num_nodes++;
          }
          else {
            *lp = nin + position;
          }
        }
        else if (ntype[i] == 2) {
          for (j = 0; j < 3; j++) {
            e = p8est_corner_edges[i][j];
            if (edges[e] >= 0 && edges[e] < 4) {
              indep[0] =
                nodes->element_nodes[elid * 8 +
                                     p8est_edge_corners[e][(edges[e] % 2)]];
              indep[1] = nodes->element_nodes[elid * 8 + i];
              break;
            }
          }
          P4EST_ASSERT (j < 3);
          indep[2] = -1;
          indep[3] = -1;
          qsort (indep, 4, sizeof (p4est_locidx_t), p4est_locidx_compare);
          lp = (p4est_locidx_t *) sc_array_index (&local_nodes, nid);
          r = (p4est_locidx_t *) sc_hash_array_insert_unique (hanging, indep,
                                                              &position);
          if (r != NULL) {
            P4EST_ASSERT (position == (size_t) (num_nodes - nin));
            r[0] = indep[0], r[1] = indep[1], r[2] = indep[2], r[3] =
              indep[3];
            *lp = num_nodes++;
          }
          else {
            *lp = nin + position;
          }
        }
      }
    }
    else {
      nid += 8;
    }
  }

  *Local_nodes = (p4est_locidx_t *) local_nodes.array;
  return num_nodes;
}

trilinear_mesh_t   *
p8est_trilinear_mesh_new_from_lnodes (p4est_t * p4est, p4est_lnodes_t * nodes)
{
  const int           num_procs = p4est->mpisize;
  const int           rank = p4est->mpirank;
  int                 mpiret;
  int                 k;
  size_t              current, zz, zy;
  int32_t             e, n;
  int64_t             global_borrowed, global_shared;
  int64_t             local_counts[5], global_counts[5];
  int32link_t        *lynk, **tail;
  p4est_topidx_t      which_tree;
  p4est_locidx_t     *local_nodes;
  p4est_locidx_t      num_local_nodes;
  p4est_tree_t       *tree;
  p4est_quadrant_t   *q;
  trilinear_elem_t   *elem;
  trilinear_anode_t  *anode;
  trilinear_dnode_t  *dnode;
  trilinear_mesh_t   *mesh;
  sc_array_t         *sharers = nodes->sharers;
  size_t              scount = sharers->elem_count, count;
  p4est_locidx_t      num_owned_shared = 0;
  p8est_lnodes_rank_t *lrank;
  point_t            *points;
  trilinear_mesh_pid_t *pids;
  p4est_locidx_t      nid;
  sc_hash_array_t    *hanging;
  p4est_locidx_t     *indep, indep_count;

  P4EST_GLOBAL_PRODUCTIONF
    ("Into trilinear_mesh_extract with %lld total elements\n",
     (long long) p4est->global_num_quadrants);

  /* Allocate output data structure. */
  mesh = P4EST_ALLOC_ZERO (trilinear_mesh_t, 1);
  memset (mesh, -1, sizeof (*mesh));

  /* Count dnodes */
  num_local_nodes = p8est_mesh_dcount (nodes, &local_nodes, &hanging);
  indep = (p4est_locidx_t *) hanging->a.array;
  points = P4EST_ALLOC (point_t, num_local_nodes);
  pids = P4EST_ALLOC (trilinear_mesh_pid_t, num_local_nodes);
  p8est_mesh_points (points, pids, local_nodes, num_local_nodes, nodes,
                     p4est);

  /* Get number of owned shared. */
  for (zz = 0; zz < scount; zz++) {
    lrank = (p8est_lnodes_rank_t *) sc_array_index (sharers, zz);
    if (lrank->rank == rank) {
      num_owned_shared = lrank->shared_mine_count;
      break;
    }
  }

  /* Assign local counts. */
  P4EST_ASSERT (nodes->num_local_elements == p4est->local_num_quadrants);
  mesh->local_elem_num = p4est->local_num_quadrants;
  mesh->local_anode_num = nodes->num_local_nodes;
  mesh->local_dnode_num = num_local_nodes - nodes->num_local_nodes;
  mesh->local_onode_num = nodes->owned_count;
  mesh->local_owned_offset = 0;
  mesh->local_node_num = num_local_nodes;

  /* Communicate global counts. */
  local_counts[0] = mesh->local_elem_num;
  local_counts[1] = mesh->local_anode_num;
  local_counts[2] = mesh->local_onode_num;
  local_counts[3] = mesh->local_dnode_num;
  local_counts[4] = num_owned_shared;
  mpiret = sc_MPI_Allreduce (local_counts, global_counts, 5,
                             sc_MPI_LONG_LONG_INT, sc_MPI_SUM,
                             p4est->mpicomm);
  SC_CHECK_MPI (mpiret);
  P4EST_ASSERT (global_counts[0] == p4est->global_num_quadrants);
  mesh->total_elem_num = global_counts[0];
  global_borrowed = global_counts[1] - global_counts[2];
  mesh->total_anode_num = global_counts[2];
  mesh->total_dnode_num = global_counts[3];
  global_shared = global_counts[4];
  mesh->total_node_num = mesh->total_anode_num + mesh->total_dnode_num;

  /* Allocate the mesh memory. */
  mesh->elem_table = P4EST_ALLOC (trilinear_elem_t, mesh->local_elem_num);
  mesh->node_table = P4EST_ALLOC (trilinear_node_t, mesh->local_node_num);
  mesh->fvnid_count_table = P4EST_ALLOC (int64_t, num_procs + 1);
  mesh->fvnid_interval_table = P4EST_ALLOC (int64_t, num_procs + 1);
  mesh->all_fvnid_start = mesh->fvnid_interval_table;
  mesh->sharer_pool = sc_mempool_new (sizeof (int32link_t));
  mesh->elem_pids = P4EST_ALLOC (trilinear_mesh_pid_t, mesh->local_elem_num);
  mesh->node_pids = P4EST_ALLOC (trilinear_mesh_pid_t, mesh->local_node_num);

  /* Assign global free variable information. */
  mesh->fvnid_interval_table[0] = 0;
  for (k = 0; k < num_procs; ++k) {
    mesh->fvnid_count_table[k] = nodes->global_owned_count[k];
    mesh->fvnid_interval_table[k + 1] = mesh->fvnid_interval_table[k] +
      mesh->fvnid_count_table[k];
  }
  mesh->fvnid_count_table[num_procs] = -1;
  mesh->global_fvnid_num = mesh->fvnid_interval_table[num_procs];
  mesh->global_fvnid_start = 0;
  mesh->global_fvnid_end = mesh->global_fvnid_num - 1;
  P4EST_ASSERT (mesh->global_fvnid_num == mesh->total_anode_num);

  /* Assign element information. */
  which_tree = p4est->first_local_tree;
  nid = 0;
  if (which_tree >= 0) {
    tree = p4est_tree_array_index (p4est->trees, which_tree);
    current = 0;
    for (e = 0; e < mesh->local_elem_num; ++e) {
      if (current == tree->quadrants.elem_count) {
        ++which_tree;
        tree = p4est_tree_array_index (p4est->trees, which_tree);
        current = 0;
      }
      q = p8est_quadrant_array_index (&tree->quadrants, current);
      elem = mesh->elem_table + e;
      for (k = 0; k < P4EST_CHILDREN; ++k) {
        elem->local_node_id[k] = local_nodes[nid++];
      }
      elem->lx = (tick_t) q->x;
      elem->ly = (tick_t) q->y;
      elem->lz = (tick_t) q->z;
      elem->size = P4EST_QUADRANT_LEN (q->level);
      elem->data = q->p.user_data;
      mesh->elem_pids[e] = (trilinear_mesh_pid_t) which_tree;
      ++current;
    }
    P4EST_ASSERT (which_tree == p4est->last_local_tree);
    P4EST_ASSERT (current == tree->quadrants.elem_count);
  }

  /* Assign anchored node information. */
  mesh->anode_table = mesh->node_table;
  mesh->onode_table = mesh->node_table + mesh->local_owned_offset;
  mesh->dnode_table = mesh->node_table + mesh->local_anode_num;
  for (n = 0; n < mesh->local_anode_num; ++n) {
    anode = &mesh->node_table[n].anchored;
    anode->point.x = points[n].x;
    anode->point.y = points[n].y;
    anode->point.z = points[n].z;
    mesh->node_pids[n] = pids[n];
    P4EST_ASSERT (pids[n] >= 0
                  && (size_t) pids[n] < p4est->trees->elem_count);
    anode->fvnid = p4est_lnodes_global_index (nodes, n);
    anode->share = NULL;
  }

  for (zz = 0; zz < scount; zz++) {
    lrank = (p8est_lnodes_rank_t *) sc_array_index (sharers, zz);
    count = lrank->shared_nodes.elem_count;
    if (lrank->rank == rank) {
      continue;
    }
    for (zy = 0; zy < count; zy++) {
      nid = *((p4est_locidx_t *) sc_array_index (&lrank->shared_nodes, zy));
      anode = &mesh->node_table[nid].anchored;
      tail = &anode->share;
      while (*tail != NULL) {
        tail = &((*tail)->next);
      }
      *tail = lynk = (int32link_t *) sc_mempool_alloc (mesh->sharer_pool);
      lynk->id = (int32_t) lrank->rank;
      lynk->next = NULL;
    }
  }

  /* Assign hanging node information. */
  indep_count = 0;
  for (; n < num_local_nodes; n++) {
    dnode = &mesh->node_table[n].dangling;
    dnode->point.x = points[n].x;
    dnode->point.y = points[n].y;
    dnode->point.z = points[n].z;
    mesh->node_pids[n] = pids[n];
    P4EST_ASSERT (pids[n] >= 0
                  && (size_t) pids[n] < p4est->trees->elem_count);
    dnode->type = 0;            /* Not used in Rhea. */
    if (indep[indep_count] == -1) {
      dnode->local_anode_id[2] = indep[indep_count++];
      P4EST_ASSERT (indep[indep_count] == -1);
      dnode->local_anode_id[3] = indep[indep_count++];
      dnode->local_anode_id[0] = indep[indep_count++];
      dnode->local_anode_id[1] = indep[indep_count++];
    }
    else {
      dnode->local_anode_id[0] = indep[indep_count++];
      dnode->local_anode_id[1] = indep[indep_count++];
      dnode->local_anode_id[2] = indep[indep_count++];
      dnode->local_anode_id[3] = indep[indep_count++];
    }
  }

  /* Assign the remaining variables. */
  mesh->mpicomm = p4est->mpicomm;
  mesh->mpisize = (int32_t) num_procs;
  mesh->mpirank = (int32_t) rank;
  mesh->recsize = (int32_t) p4est->data_size;
  mesh->destructor = p8est_trilinear_mesh_destroy;

  /* These members are incomplete and need to be filled later. */
  memset (mesh->bounds, 0, 6 * sizeof (int));
  memset (mesh->sizes, 0, 3 * sizeof (int));
  mesh->minsize = mesh->maxsize = 0;
  mesh->ticksize = 0.;
  mesh->extra_info = NULL;
  mesh->gid = -1;

  /* We are done */
  P4EST_GLOBAL_PRODUCTIONF ("Done trilinear_mesh_extract"
                            " with %lld anodes %lld %lld\n",
                            (long long) mesh->total_anode_num,
                            (long long) global_borrowed,
                            (long long) global_shared);

  P4EST_FREE (points);
  P4EST_FREE (pids);
  SC_FREE (local_nodes);
  sc_hash_array_destroy (hanging);
  return mesh;
}
