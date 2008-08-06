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

/*
 * This file provides an interface that is compatible
 * with the trilinear finite elements used in the Rhea code.
 */

#ifndef P8EST_TRILINEAR_H
#define P8EST_TRILINEAR_H

#include <p8est_mesh.h>

/** tick_t: The unit of the embeded integer domain. */
typedef int32_t     tick_t;

/** point_t: 3D coordinates in tick units. */
typedef struct point
{
  tick_t              x, y, z;
}
point_t;

/** int32link_t: Link list for int32_t's. */
typedef struct int32link
{
  int32_t             id;
  struct int32link   *next;
}
int32link_t;

/* *INDENT-OFF* */
/** Encode the boundary status of a node. */
typedef enum trilinear_boundary_enum
{
  TRILINEAR_BOUNDARY_NONE      =      0,
  TRILINEAR_BOUNDARY_IS_LEFT   = 0x0001,
  TRILINEAR_BOUNDARY_IS_RIGHT  = 0x0002,
  TRILINEAR_BOUNDARY_IS_FRONT  = 0x0004,
  TRILINEAR_BOUNDARY_IS_BACK   = 0x0008,
  TRILINEAR_BOUNDARY_IS_BOTTOM = 0x0010,
  TRILINEAR_BOUNDARY_IS_TOP    = 0x0020,
  TRILINEAR_BOUNDARY_IS_EDGE   = 0x0040,
  TRILINEAR_BOUNDARY_IS_CORNER = 0x0080,
  TRILINEAR_BOUNDARY_IS_ORIGIN = 0x0100,
  TRILINEAR_BOUNDARY_IS_XBC    = (TRILINEAR_BOUNDARY_IS_LEFT |
                                  TRILINEAR_BOUNDARY_IS_RIGHT),
  TRILINEAR_BOUNDARY_IS_YBC    = (TRILINEAR_BOUNDARY_IS_FRONT |
                                  TRILINEAR_BOUNDARY_IS_BACK),
  TRILINEAR_BOUNDARY_IS_ZBC    = (TRILINEAR_BOUNDARY_IS_BOTTOM |
                                  TRILINEAR_BOUNDARY_IS_TOP),
  TRILINEAR_BOUNDARY_IS_FACE   = (TRILINEAR_BOUNDARY_IS_XBC |
                                  TRILINEAR_BOUNDARY_IS_YBC |
                                  TRILINEAR_BOUNDARY_IS_ZBC)
}
trilinear_boundary_enum_t;
/* *INDENT-ON* */

/* This integer is big enough to hold above flags. */
typedef uint16_t    trilinear_boundary_flag_t;

/** */
typedef struct trilinear_elem
{
  int32_t             local_node_id[8]; /* indices into local node table */
  tick_t              size;     /* size of the element in ticks */
  void               *data;     /* pointer to its data record */
}
trilinear_elem_t;

/** trilinear_anode_t: anchored node (associated with free variables). */
typedef struct trilinear_anode
{
  point_t             point;
  int64_t             fvnid;
  int32link_t        *share;    /* processors that share this anchored node */
}
trilinear_anode_t;

/** trilinear_dnode_t: dangling node (associcated with dependent variables). */
typedef struct trilinear_dnode
{
  point_t             point;
  int32_t             type;     /* not used in Rhea */
  int32_t             local_anode_id[4];        /* [2] is -1 for edge nodes */
}
trilinear_dnode_t;

/** trilinear_node_t: union of both node types. */
typedef union trilinear_node
{
  point_t             point;
  trilinear_anode_t   anchored;
  trilinear_dnode_t   dangling;
}
trilinear_node_t;

/*
 * This is the structure precomputed after mesh extraction.
 * It is ordered differently than trilinear_element_into_t.
 * First come the anchored nodes, then the dangling nodes.
 * Variables beginning with Q are over direct and indirect anchored nodes (iqc).
 * Variable dQcolumn is indexed by dangling nodes only as
 *   dQcolumn[nd][l] = iqc of depended anchored node l.
 * If interior_anchors_only is nonzero, then the element has only
 *   anchored nodes as corners that are not on the domain boundary.
 */
typedef struct trilinear_element_info2
{
  int8_t              nanchored, ndangling;
  int8_t              interior_only, interior_anchors_only;
  int8_t              corner[8];
  int8_t              Qisdirect[8];
  int32_t             Qindices[8];
  int64_t             Qfvnids[8];
  trilinear_boundary_flag_t Qboundary[8];
  int8_t              dQcolumn[8][4];
}
trilinear_element_info2_t;

typedef struct trilinear_mesh_extra
{
  int32_t             shared_elem_num;
  int32_t            *shared_elem_ids;
  trilinear_element_info2_t *info2;
}
trilinear_mesh_extra_t;

typedef struct trilinear_mesh
{
  /* Global mesh statistics */
  int64_t             total_elem_num;
  int64_t             total_node_num;
  int64_t             total_anode_num;
  int64_t             total_dnode_num;

  /* Local mesh parameters */
  int32_t             local_elem_num;   /* number of element on this processor */
  int32_t             local_node_num;   /* number of anchored and dangling nodes */

  /* The first part of the node table contains anchored nodes. The free
     variables (anchored nodes) owned by a processor is clustered together
     with no holes within the first part of the table.

     The second part of the node table contains dangling nodes. */

  int32_t             local_anode_num;  /* number of anchored nodes */
  int32_t             local_onode_num;  /* number of owned anchored nodes */
  int32_t             local_dnode_num;  /* number of dangling nodes */

  int32_t             local_owned_offset;       /* offset to the first
                                                   owned anchored node */

  /* Memory allocated to hold the trilinear elements and nodes */
  trilinear_elem_t   *elem_table;
  trilinear_node_t   *node_table;

  /* Memory allocated for free variable interval table.
   * The first two are allocated with (num_procs + 1) entries each.
   * The third (all_fvnid_start) is a convenience pointer.
   */
  int64_t            *fvnid_count_table;
  int64_t            *fvnid_interval_table;
  int64_t            *all_fvnid_start;

  /* Convenience pointers into the node table. */
  trilinear_node_t   *anode_table;
  trilinear_node_t   *onode_table;
  trilinear_node_t   *dnode_table;

  /* Convenience variables recording the total number of free variables,
     the starting id and the ending id. identical on all processors */
  int64_t             global_fvnid_num; /* total number of fvnids */
  int64_t             global_fvnid_start;       /* first global fvnid */
  int64_t             global_fvnid_end; /* last global fvnid */

  tick_t              bounds[3][2];
  tick_t              sizes[3], minsize, maxsize;
  double              ticksize, volume;

  MPI_Comm            mpicomm;
  int32_t             mpisize, mpirank;
  int32_t             recsize;

  void                (*destructor) (struct trilinear_mesh *);
  trilinear_mesh_extra_t *extra_info;

  /* Everything below here is not present in Rhea. */
  sc_mempool_t       *sharer_pool;      /* allocator for node sharers */
}
trilinear_mesh_t;

/** Creates a trilinear mesh structure from a p8est and its node data.
 */
trilinear_mesh_t   *p8est_trilinear_mesh_new (p8est_t * p8est,
                                              p8est_nodes_t * nodes);

/** Frees a trilinear mesh structure.
 */
void                p8est_trilinear_mesh_destroy (trilinear_mesh_t * mesh);

/** p8est_trilinear_neighbor_t
 *
 * The face neighbors are ordered in -x, +x, -y, +y, -z, and +z directions.
 *
 * For each direction, the neighbor(s) may be:
 *
 * 1. Out of the domain.
 *
 *      face_neighbor_eid[direction][0] = -1.
 *      face_neighbor_eid[direction][1--3] are undefined.
 *
 * 2. As large or twice as large as the current element:
 *
 *      face_neighbor_eid[direction][0] = index
 *
 *    where ((index >= 0) && (index < local_elem_num)) if the neighbor is
 *    LOCAL, or
 *    ((index >= local_elem_num) &&
 *    (index < (local_elem_num + ghost_elem_num))) if the neighor is
 *    REMOTE.
 *
 *      face_neighbor_eid[direction][1--3] are undefined.
 *
 * 3. Half as large as the current element:
 *
 *      face_neighbor_eid[direction][i] = index
 *
 *    where index is defined as above. Note that in this case all four
 *    neighbors (half as large) must exist.
 */

typedef struct local_neighbor
{
  int32_t             face_neighbor_eid[6][4];
}
local_neighbor_t;

typedef struct ghost_elem
{
  tick_t              lx, ly, lz;
  tick_t              size;
  int32_t             owner_procid;     /* remote processor id */
  int32_t             reid;     /* remote processor element index */
}
ghost_elem_t;

typedef struct trilinear_neighborhood
{
  int32_t             ghost_elem_num;
  ghost_elem_t       *ghost_elem_table;
  local_neighbor_t   *local_elem_neighbor_table;
}
trilinear_neighborhood_t;

trilinear_neighborhood_t *trilinear_neighborhood_new (p8est_t * p8est);

/* *INDENT-OFF* */
void
trilinear_neighborhood_destroy (trilinear_neighborhood_t * neighborhood);
/* *INDENT-ON* */

#endif /* !P8EST_TRILINEAR_H */
