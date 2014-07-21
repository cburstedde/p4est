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
#include <p8est_trilinear.h>

trilinear_mesh_t   *
p8est_trilinear_mesh_new_from_nodes (p4est_t * p4est, p4est_nodes_t * nodes)
{
  const int           num_procs = p4est->mpisize;
  const int           rank = p4est->mpirank;
  int                 mpiret;
  int                 k, owner;
#ifdef P4EST_DEBUG
  int                 prev_owner = 0;
  int64_t             prev_fvnid = -1;
#endif
  int                *sharers;
  size_t              current, zz, num_sharers;
  int32_t             e, n, local_owned_end;
  int64_t             global_borrowed, global_shared;
  int64_t             local_counts[5], global_counts[5];
  int32link_t        *lynk, **tail;
  p4est_topidx_t      which_tree;
  p4est_locidx_t     *local_node, *shared_offsets;
  p4est_tree_t       *tree;
  p4est_quadrant_t   *q;
  p4est_indep_t      *in;
  p8est_hang4_t      *fh;
  p8est_hang2_t      *eh;
  trilinear_elem_t   *elem;
  trilinear_anode_t  *anode;
  trilinear_dnode_t  *dnode;
  trilinear_mesh_t   *mesh;
  trilinear_mesh_pid_t *elem_pids;
  trilinear_mesh_pid_t *node_pids;
  sc_recycle_array_t *rarr;

  P4EST_GLOBAL_PRODUCTIONF
    ("Into trilinear_mesh_extract with %lld total elements\n",
     (long long) p4est->global_num_quadrants);
  p4est_log_indent_push ();

  /* Allocate output data structure. */
  mesh = P4EST_ALLOC_ZERO (trilinear_mesh_t, 1);
  memset (mesh, -1, sizeof (*mesh));
  shared_offsets = nodes->shared_offsets;

  /* Assign local counts. */
  P4EST_ASSERT (nodes->num_local_quadrants == p4est->local_num_quadrants);
  mesh->local_elem_num = p4est->local_num_quadrants;
  mesh->local_anode_num = nodes->indep_nodes.elem_count;
  mesh->local_dnode_num =
    nodes->face_hangings.elem_count + nodes->edge_hangings.elem_count;
  mesh->local_onode_num = nodes->num_owned_indeps;
  mesh->local_owned_offset = nodes->offset_owned_indeps;
  mesh->local_node_num = mesh->local_anode_num + mesh->local_dnode_num;
  local_owned_end = mesh->local_owned_offset + mesh->local_onode_num;

  /* Communicate global counts. */
  local_counts[0] = mesh->local_elem_num;
  local_counts[1] = mesh->local_anode_num;
  local_counts[2] = mesh->local_onode_num;
  local_counts[3] = mesh->local_dnode_num;
  local_counts[4] = nodes->num_owned_shared;
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
  mesh->elem_pids = P4EST_ALLOC (trilinear_mesh_pid_t, mesh->local_node_num);
  mesh->node_pids = P4EST_ALLOC (trilinear_mesh_pid_t, mesh->local_node_num);

  /* Assign global free variable information. */
  mesh->fvnid_interval_table[0] = 0;
  for (k = 0; k < num_procs; ++k) {
    mesh->fvnid_interval_table[k + 1] = mesh->fvnid_interval_table[k] +
      (mesh->fvnid_count_table[k] = nodes->global_owned_indeps[k]);
  }
  mesh->fvnid_count_table[num_procs] = -1;
  mesh->global_fvnid_num = mesh->fvnid_interval_table[num_procs];
  mesh->global_fvnid_start = 0;
  mesh->global_fvnid_end = mesh->global_fvnid_num - 1;
  P4EST_ASSERT (mesh->global_fvnid_num == mesh->total_anode_num);

  /* Assign element information. */
  local_node = nodes->local_nodes;
  which_tree = p4est->first_local_tree;
  elem_pids = mesh->elem_pids;
  if (which_tree >= 0) {
    tree = p4est_tree_array_index (p4est->trees, which_tree);
    current = 0;
    for (e = 0; e < mesh->local_elem_num; ++e) {
      if (current == tree->quadrants.elem_count) {
        ++which_tree;
        tree = p4est_tree_array_index (p4est->trees, which_tree);
        current = 0;
      }
      q = p4est_quadrant_array_index (&tree->quadrants, current);
      elem = mesh->elem_table + e;
      for (k = 0; k < P4EST_CHILDREN; ++k) {
        elem->local_node_id[k] = *local_node++;
      }
      elem->lx = (tick_t) q->x;
      elem->ly = (tick_t) q->y;
      elem->lz = (tick_t) q->z;
      elem->size = P4EST_QUADRANT_LEN (q->level);
      elem->data = q->p.user_data;
      elem_pids[e] = (trilinear_mesh_pid_t) which_tree;
      ++current;
    }
    P4EST_ASSERT (which_tree == p4est->last_local_tree);
    P4EST_ASSERT (current == tree->quadrants.elem_count);
  }

  /* Assign anchored node information. */
  mesh->anode_table = mesh->node_table;
  mesh->onode_table = mesh->node_table + mesh->local_owned_offset;
  mesh->dnode_table = mesh->node_table + mesh->local_anode_num;
  node_pids = mesh->node_pids;
  for (n = 0; n < mesh->local_anode_num; ++n) {
    anode = &mesh->node_table[n].anchored;
    in = (p4est_indep_t *) sc_array_index (&nodes->indep_nodes, (size_t) n);
    anode->point.x = in->x;
    anode->point.y = in->y;
    anode->point.z = in->z;
    node_pids[n] = (trilinear_mesh_pid_t) in->p.piggy3.which_tree;
    if (n < mesh->local_owned_offset) {
      owner = nodes->nonlocal_ranks[n];
      P4EST_ASSERT (owner < rank && owner >= prev_owner);
    }
    else if (n >= local_owned_end) {
      owner = nodes->nonlocal_ranks[n - mesh->local_onode_num];
      P4EST_ASSERT (owner > rank && owner >= prev_owner);
    }
    else {
      owner = rank;
    }
    anode->fvnid = mesh->all_fvnid_start[owner] + in->p.piggy3.local_num;
    P4EST_ASSERT (anode->fvnid > prev_fvnid);
    if (in->pad8 == 0) {
      P4EST_ASSERT (in->pad16 == -1);
      P4EST_ASSERT (shared_offsets == NULL || shared_offsets[n] == -1);
      anode->share = NULL;
    }
    else {
      P4EST_ASSERT (in->pad8 > 0);
      num_sharers = (size_t) in->pad8;
      rarr =
        (sc_recycle_array_t *) sc_array_index (&nodes->shared_indeps,
                                               num_sharers - 1);
      if (nodes->shared_offsets == NULL) {
        P4EST_ASSERT (in->pad16 >= 0);
        zz = (size_t) in->pad16;
      }
      else {
        P4EST_ASSERT (in->pad16 == -1);
        zz = (size_t) shared_offsets[n];
      }
      sharers = (int *) sc_array_index (&rarr->a, zz);
      tail = &anode->share;
      for (zz = 0; zz < num_sharers; ++zz) {
        *tail = lynk = (int32link_t *) sc_mempool_alloc (mesh->sharer_pool);
        lynk->id = (int32_t) sharers[zz];
        tail = &lynk->next;
      }
      *tail = NULL;
    }
#ifdef P4EST_DEBUG
    prev_owner = owner;
    prev_fvnid = anode->fvnid;
#endif
  }

  /* Assign face hanging node information. */
  for (zz = 0; zz < nodes->face_hangings.elem_count; ++n, ++zz) {
    dnode = &mesh->node_table[n].dangling;
    fh = (p8est_hang4_t *) sc_array_index (&nodes->face_hangings, zz);
    dnode->point.x = fh->x;
    dnode->point.y = fh->y;
    dnode->point.z = fh->z;
    dnode->type = 0;            /* Not used in Rhea. */
    dnode->local_anode_id[0] = fh->p.piggy.depends[0];
    dnode->local_anode_id[1] = fh->p.piggy.depends[1];
    dnode->local_anode_id[2] = fh->p.piggy.depends[2];
    dnode->local_anode_id[3] = fh->p.piggy.depends[3];
    node_pids[n] = (trilinear_mesh_pid_t) fh->p.piggy.which_tree;
  }

  /* Assign edge hanging node information. */
  for (zz = 0; zz < nodes->edge_hangings.elem_count; ++n, ++zz) {
    dnode = &mesh->node_table[n].dangling;
    eh = (p8est_hang2_t *) sc_array_index (&nodes->edge_hangings, zz);
    dnode->point.x = eh->x;
    dnode->point.y = eh->y;
    dnode->point.z = eh->z;
    dnode->type = 0;            /* Not used in Rhea. */
    dnode->local_anode_id[0] = eh->p.piggy.depends[0];
    dnode->local_anode_id[1] = eh->p.piggy.depends[1];
    dnode->local_anode_id[2] = dnode->local_anode_id[3] = -1;
    node_pids[n] = (trilinear_mesh_pid_t) eh->p.piggy.which_tree;
  }
  P4EST_ASSERT (n == mesh->local_node_num);

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
  p4est_log_indent_pop ();
  P4EST_GLOBAL_PRODUCTIONF ("Done trilinear_mesh_extract"
                            " with %lld anodes %lld %lld\n",
                            (long long) mesh->total_anode_num,
                            (long long) global_borrowed,
                            (long long) global_shared);

  return mesh;
}
