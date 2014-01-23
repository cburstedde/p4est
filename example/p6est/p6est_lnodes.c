/*
  This file is part of p4est.
  p4est is a C library to manage a collection (a forest) of multiple
  connected adaptive quadtrees or octrees in parallel.

  Copyright (C) 2014 The University of Texas System
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

#include <p6est_lnodes.h>
#include <p6est_profile.h>

p6est_lnodes_t     *
p6est_lnodes_new (p6est_t * p6est, p6est_ghost_t * ghost, int degree)
{
  p6est_lnodes_t     *lnodes;
  p6est_profile_t    *profile;
  p4est_lnodes_t     *clnodes;
  int                 nperelem = (degree + 1) * (degree + 1) * (degree + 1);
  //int nperface = (degree - 1) * (degree - 1);
  //int nperedge = (degree - 1);
  p4est_locidx_t      ncid, cid, enid, *en;
  p4est_locidx_t      nnodecols;
  p4est_locidx_t      nelemcols;
  p4est_locidx_t      nll;
  p4est_locidx_t      nlayers;
  p4est_locidx_t     *layernodecount;
  p4est_locidx_t     *layernodeoffsets;
  p4est_locidx_t (*lr)[2];
  p4est_locidx_t      ncolnodes;
  p4est_locidx_t     *global_owned_count;
  p4est_locidx_t      num_owned, num_local;
  p4est_gloidx_t      gnum_owned, offset;
  p4est_gloidx_t     *owned_offsets;
  int                 i, j;
  int                 is_owned;
  int                 mpisize = p6est->mpisize;
  int                 mpiret;
  sc_array_t          lnoview;
  size_t              zz, nsharers;

  P4EST_ASSERT (degree >= 1);

  lnodes = P4EST_ALLOC (p6est_lnodes_t, 1);

  /* first get the profile */
  profile = p6est_profile_new_local (p6est, ghost, P6EST_PROFILE_INTERSECTION,
                                     P8EST_CONNECT_DEFAULT);
  p6est_profile_sync (profile);

  lr = (p4est_locidx_t (*)[2]) profile->lnode_ranges;

  clnodes = profile->lnodes;

  nnodecols = clnodes->num_local_nodes;
  nelemcols = clnodes->num_local_elements;
  en = clnodes->element_nodes;
  layernodecount = P4EST_ALLOC_ZERO (p4est_locidx_t, nnodecols);
  layernodeoffsets = P4EST_ALLOC_ZERO (p4est_locidx_t, nnodecols);
  num_owned = 0;
  num_local = 0;
  for (cid = 0, enid = 0; cid < nelemcols; cid++) {
    for (j = 0; j < 3; j++) {
      for (i = 0; i < 3; i++, enid++) {
        ncid = en[enid];
        is_owned = (ncid < clnodes->owned_count);
        nlayers = lr[ncid][1];
        ncolnodes = nlayers * degree + 1;
        if (i != 1 && j != 1) {
          /* this is a corner column */
          layernodecount[ncid] = ncolnodes;
        }
        else if (i == 1 && j == 1) {
          /* this is a quad column */
          layernodecount[ncid] = ncolnodes * (degree - 1) * (degree - 1);
        }
        else {
          /* this is a face column */
          layernodecount[ncid] = ncolnodes * (degree - 1);
        }
        num_local += layernodecount[ncid];
        if (is_owned) {
          num_owned += layernodecount[ncid];
        }
      }
    }
  }

  if (nnodecols) {
    layernodeoffsets[0] = 0;
    for (ncid = 0; ncid < nnodecols; ncid++) {
      layernodeoffsets[ncid + 1] = layernodeoffsets[ncid] +
        layernodecount[ncid];
    }
  }

  gnum_owned = num_owned;

  owned_offsets = P4EST_ALLOC (p4est_gloidx_t, mpisize + 1);
  global_owned_count = P4EST_ALLOC (p4est_locidx_t, mpisize);

  mpiret = MPI_Allgather (&gnum_owned, 1, P4EST_MPI_GLOIDX,
                          owned_offsets, 1, P4EST_MPI_GLOIDX, p6est->mpicomm);
  SC_CHECK_MPI (mpiret);

  offset = 0;
  for (i = 0; i < mpisize; i++) {
    global_owned_count[i] = (p4est_locidx_t) owned_offsets[i];
    gnum_owned = owned_offsets[i];
    owned_offsets[i] = offset;
    offset += gnum_owned;
  }
  owned_offsets[mpisize] = offset;

  nll = p6est->layers->elem_count;
  nsharers = clnodes->sharers->elem_count;

  lnodes->mpicomm = p6est->mpicomm;
  lnodes->num_local_nodes = num_local;
  lnodes->owned_count = num_owned;
  lnodes->global_offset = owned_offsets[p6est->mpirank];
  lnodes->nonlocal_nodes = P4EST_ALLOC (p4est_gloidx_t, num_local - num_owned); /* TODO */
  lnodes->sharers =
    sc_array_new_size (sizeof (p6est_lnodes_rank_t), nsharers);
  lnodes->global_owned_count = global_owned_count;

  lnodes->degree = degree;
  lnodes->vnodes = nperelem;
  lnodes->num_local_elements = nll;
  lnodes->face_code = P4EST_ALLOC (p6est_lnodes_code_t, nll);
  lnodes->element_nodes = P4EST_ALLOC (p4est_locidx_t, nperelem * nll);

  p6est_profile_element_to_node (p6est, profile, degree, layernodeoffsets,
                                 lnodes->element_nodes, lnodes->face_code);

  for (zz = 0; zz < nsharers; zz++) {
    p4est_lnodes_rank_t *crank = p4est_lnodes_rank_array_index
      (clnodes->sharers, zz);
    p6est_lnodes_rank_t *rank = p6est_lnodes_rank_array_index
      (lnodes->sharers, zz);
    size_t              zy;
    size_t              nshared;

    rank->rank = crank->rank;
    sc_array_init (&rank->shared_nodes, sizeof (p4est_locidx_t));
    nshared = crank->shared_nodes.elem_count;

    rank->owned_offset = -1;
    rank->owned_count = 0;
    rank->shared_mine_count = 0;
    rank->shared_mine_offset = -1;
    for (zy = 0; zy < nshared; zy++) {
      p4est_locidx_t      cnid =
        *((p4est_locidx_t *) sc_array_index (&crank->shared_nodes, zy));
      p4est_locidx_t     *lp;
      p4est_locidx_t      nthis, il;
      p4est_locidx_t      old_count = rank->shared_nodes.elem_count;

      nthis = layernodecount[cnid];
      lp =
        (p4est_locidx_t *) sc_array_push_count (&rank->shared_nodes, nthis);

      for (il = 0; il < nthis; il++) {
        lp[il] = layernodeoffsets[cnid] + il;
        if (zy >= crank->shared_mine_offset
            && (p4est_locidx_t) zy - crank->shared_mine_offset <
            crank->shared_mine_count) {
          rank->shared_mine_count++;
          if (rank->shared_mine_offset == -1) {
            rank->shared_mine_offset = old_count + il;
          }
        }
        if (cnid >= crank->owned_offset
            && cnid - crank->owned_offset < crank->owned_count) {
          rank->owned_count++;
          if (rank->owned_offset == -1) {
            rank->owned_offset = lp[il];
          }
        }
      }
    }
  }

  memcpy (layernodecount, layernodeoffsets,
          nnodecols * sizeof (p4est_locidx_t));
  sc_array_init_data (&lnoview, layernodecount, sizeof (p4est_locidx_t),
                      (size_t) nnodecols);

  p4est_lnodes_share_owned (&lnoview, clnodes);

  for (zz = 0; zz < nsharers; zz++) {
    p4est_lnodes_rank_t *crank = p4est_lnodes_rank_array_index
      (clnodes->sharers, zz);

    if (crank->rank == p6est->mpirank) {
      continue;
    }

    for (ncid = crank->owned_offset;
         ncid < crank->owned_offset + crank->owned_count; ncid++) {
      p4est_gloidx_t      owners_offset;
      p4est_locidx_t      nid;

      owners_offset = owned_offsets[crank->rank] + layernodecount[ncid];
      for (nid = layernodeoffsets[ncid]; nid < layernodeoffsets[ncid + 1];
           nid++) {
        P4EST_ASSERT (nid >= num_owned);
        P4EST_ASSERT (nid < num_local);
        lnodes->nonlocal_nodes[nid - num_owned] = owners_offset++;
      }
    }
  }

  p6est_profile_destroy (profile);

  P4EST_FREE (owned_offsets);
  P4EST_FREE (layernodecount);
  P4EST_FREE (layernodeoffsets);

  return lnodes;
}
