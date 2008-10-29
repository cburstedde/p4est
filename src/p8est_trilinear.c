/*
  This file is part of p4est.
  p4est is a C library to manage a parallel collection of quadtrees and/or
  octrees.

  Copyright (C) 2007,2008 Carsten Burstedde, Lucas Wilcox.

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

#include <p4est_to_p8est.h>
#include <p8est_trilinear.h>

trilinear_mesh_t   *
p8est_trilinear_mesh_new (p4est_t * p4est, p4est_nodes_t * nodes)
{
  const int           num_procs = p4est->mpisize;
  const int           rank = p4est->mpirank;
#ifdef P4EST_MPI
  int                 mpiret;
#endif
  int                 k, prev_owner, owner;
  int                *sharers;
  size_t              current, zz, num_sharers;
  int32_t             e, n, local_owned_end;
  int64_t             global_borrowed, global_shared, prev_fvnid;
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
  sc_recycle_array_t *rarr;

  P4EST_GLOBAL_PRODUCTIONF
    ("Into trilinear_mesh_extract with %lld total elements\n",
     (long long) p4est->global_num_quadrants);

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
  memcpy (global_counts, local_counts, 5 * sizeof (int64_t));
#ifdef P4EST_MPI
  if (p4est->mpicomm != MPI_COMM_NULL) {
    mpiret = MPI_Allreduce (local_counts, global_counts, 5, MPI_LONG_LONG_INT,
                            MPI_SUM, p4est->mpicomm);
    SC_CHECK_MPI (mpiret);
  }
#endif
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
  if (which_tree >= 0) {
    tree = p4est_array_index_topidx (p4est->trees, which_tree);
    current = 0;
    for (e = 0; e < mesh->local_elem_num; ++e) {
      if (current == tree->quadrants.elem_count) {
        ++which_tree;
        tree = p4est_array_index_topidx (p4est->trees, which_tree);
        current = 0;
      }
      q = sc_array_index (&tree->quadrants, current);
      elem = mesh->elem_table + e;
      for (k = 0; k < P4EST_CHILDREN; ++k) {
        elem->local_node_id[k] = *local_node++;
      }
      elem->size = P4EST_QUADRANT_LEN (q->level);
      elem->data = q->p.user_data;
      ++current;
    }
    P4EST_ASSERT (which_tree == p4est->last_local_tree);
    P4EST_ASSERT (current == tree->quadrants.elem_count);
  }

  /* Assign anchored node information. */
  prev_owner = 0;
  prev_fvnid = -1;
  mesh->anode_table = mesh->node_table;
  mesh->onode_table = mesh->node_table + mesh->local_owned_offset;
  mesh->dnode_table = mesh->node_table + mesh->local_anode_num;
  for (n = 0; n < mesh->local_anode_num; ++n) {
    anode = &mesh->node_table[n].anchored;
    in = sc_array_index (&nodes->indep_nodes, n);
    anode->point.x = in->x;
    anode->point.y = in->y;
    anode->point.z = in->z;
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
      rarr = sc_array_index (&nodes->shared_indeps, num_sharers - 1);
      if (nodes->shared_offsets == NULL) {
        P4EST_ASSERT (in->pad16 >= 0);
        zz = (size_t) in->pad16;
      }
      else {
        P4EST_ASSERT (in->pad16 == -1);
        zz = (size_t) shared_offsets[n];
      }
      sharers = sc_array_index (&rarr->a, zz);
      tail = &anode->share;
      for (zz = 0; zz < num_sharers; ++zz) {
        *tail = lynk = sc_mempool_alloc (mesh->sharer_pool);
        lynk->id = (int32_t) sharers[zz];
        tail = &lynk->next;
      }
      *tail = NULL;
    }
    prev_owner = owner;
    prev_fvnid = anode->fvnid;
  }

  /* Assign face hanging node information. */
  for (zz = 0; zz < nodes->face_hangings.elem_count; ++n, ++zz) {
    dnode = &mesh->node_table[n].dangling;
    fh = sc_array_index (&nodes->face_hangings, zz);
    dnode->point.x = fh->x;
    dnode->point.y = fh->y;
    dnode->point.z = fh->z;
    dnode->type = 0;            /* Not used in Rhea. */
    dnode->local_anode_id[0] = fh->p.piggy.depends[0];
    dnode->local_anode_id[1] = fh->p.piggy.depends[1];
    dnode->local_anode_id[2] = fh->p.piggy.depends[2];
    dnode->local_anode_id[3] = fh->p.piggy.depends[3];
  }

  /* Assign edge hanging node information. */
  for (zz = 0; zz < nodes->edge_hangings.elem_count; ++n, ++zz) {
    dnode = &mesh->node_table[n].dangling;
    eh = sc_array_index (&nodes->edge_hangings, zz);
    dnode->point.x = eh->x;
    dnode->point.y = eh->y;
    dnode->point.z = eh->z;
    dnode->type = 0;            /* Not used in Rhea. */
    dnode->local_anode_id[0] = eh->p.piggy.depends[0];
    dnode->local_anode_id[1] = eh->p.piggy.depends[1];
    dnode->local_anode_id[2] = dnode->local_anode_id[3] = -1;
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
  mesh->ticksize = mesh->volume = 0.;
  mesh->extra_info = NULL;
  mesh->patch_ids = NULL;

  /* We are done */
  P4EST_GLOBAL_PRODUCTIONF ("Done trilinear_mesh_extract"
                            " with %lld anodes %lld %lld\n",
                            (long long) mesh->total_anode_num,
                            (long long) global_borrowed,
                            (long long) global_shared);

  return mesh;
}

void
p8est_trilinear_mesh_destroy (trilinear_mesh_t * mesh)
{
  P4EST_ASSERT (mesh->destructor == p8est_trilinear_mesh_destroy);

  P4EST_FREE (mesh->elem_table);
  P4EST_FREE (mesh->node_table);
  P4EST_FREE (mesh->fvnid_count_table);
  P4EST_FREE (mesh->fvnid_interval_table);
  sc_mempool_destroy (mesh->sharer_pool);

  P4EST_FREE (mesh->patch_ids);

  P4EST_FREE (mesh);
}

/* EOF p8est_trilinear.c */
