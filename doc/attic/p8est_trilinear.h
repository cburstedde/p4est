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

/*
 * This file provides an interface that is compatible
 * with the trilinear finite elements used in the Rhea code.
 */

#ifndef P8EST_TRILINEAR_H
#define P8EST_TRILINEAR_H

#include <p8est_nodes.h>
#include <p8est_lnodes.h>

/*
 * BEGIN verbatim copy of trilinear_mesh_types.h from the Rhea code
 */

#ifndef TRILINEAR_MESH_TYPES_H
#define TRILINEAR_MESH_TYPES_H

#include <sc_containers.h>

/*
 * This file contains typedefs and macros related to the trilinear mesh.
 * It should not contain specific dependencies to either octor of p4est.
 */

#ifndef OCTOR_TICKT_TYPES
#define OCTOR_TICKT_TYPES

typedef int32_t     tick_t;

typedef struct point_t
{
  tick_t              x, y, z;
}
point_t;

#endif /* !OCTOR_TICKT_TYPES */

typedef int16_t     trilinear_mesh_pid_t;

typedef enum
{
  ANCHORED = 0,
  DANGLING_ON_XEDGE = -1,
  DANGLING_ON_YEDGE = -2,
  DANGLING_ON_ZEDGE = -3,
  DANGLING_ON_XFACE = -4,
  DANGLING_ON_YFACE = -5,
  DANGLING_ON_ZFACE = -6
}
trilinear_node_type_t;

typedef struct trilinear_elem_t
{
  int32_t             local_node_id[8]; /* indices into local node table */
  tick_t              lx, ly, lz;       /* lower-left element coordinates */
  tick_t              size;     /* size of the element in ticks */
  void               *data;     /* pointer to its record, managed by octor */
}
trilinear_elem_t;

/**
 * int32link_t: linked list entry recording processor ids
 */
typedef struct int32link_t
{
  int32_t             id;
  struct int32link_t *next;
}
int32link_t;

/**
 * trilinear_anode_t: anchored node (associated with free variables)
 */
typedef struct trilinear_anode_t
{
  point_t             point;
  int64_t             fvnid;
  int32link_t        *share;    /* processors that share this anchored node */
}
trilinear_anode_t;

typedef struct trilinear_dnode_t
{
  point_t             point;
  int32_t             type;
  int32_t             local_anode_id[4];        /* _id[2] is -1 for an edge node */
}
trilinear_dnode_t;

typedef union trilinear_node_t
{
  point_t             point;
  trilinear_anode_t   anchored;
  trilinear_dnode_t   dangling;
}
trilinear_node_t;

/* *INDENT-OFF* */
/* define below the smallest integer type that holds these flags */
typedef enum trilinear_boundary_enum
{
  TRILINEAR_BOUNDARY_NONE       =      0,
  TRILINEAR_BOUNDARY_IS_LEFT    = 0x0001,
  TRILINEAR_BOUNDARY_IS_RIGHT   = 0x0002,
  TRILINEAR_BOUNDARY_IS_FRONT   = 0x0004,
  TRILINEAR_BOUNDARY_IS_BACK    = 0x0008,
  TRILINEAR_BOUNDARY_IS_BOTTOM  = 0x0010,
  TRILINEAR_BOUNDARY_IS_TOP     = 0x0020,
  TRILINEAR_BOUNDARY_IS_EDGE    = 0x0040,
  TRILINEAR_BOUNDARY_IS_CORNER  = 0x0080,
  TRILINEAR_BOUNDARY_IS_ORIGIN  = 0x0100,
  TRILINEAR_BOUNDARY_IS_3EDGE   = 0x0200,
  TRILINEAR_BOUNDARY_IS_3CORNER = 0x0400,
  TRILINEAR_BOUNDARY_IS_PRDCX   = 0x0800,
  TRILINEAR_BOUNDARY_IS_PRDCY   = 0x1000,
  TRILINEAR_BOUNDARY_IS_PRDCZ   = 0x2000,

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

/* this integer is big enough to hold a trilinear_boundary_enum_t */
typedef uint16_t    trilinear_boundary_flag_t;

/**
 * trilinear_element_info2: a structure precomputed after mesh extraction.
 * First come the anchored nodes, then the dangling nodes.
 * Variables beginning with Q are over direct and indirect anchored nodes (iqc).
 * Variable dQcolumn is indexed by dangling nodes only as
 * dQcolumn[nd][l] = iqc of depended anchored node l.
 * If interior_anchors_only is nonzero, then the element has only
 * anchored nodes as corners that are not on the domain boundary.
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

/**
 * trilinear_mesh_t: main mesh structure used in Rhea
 */
typedef struct trilinear_mesh_t
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

  /* Memory allocated by Octor to hold the trilinear elements and nodes */
  trilinear_elem_t   *elem_table;
  trilinear_node_t   *node_table;

  /* Memory allocated for free variable interval table. Has
     (groupsize + 1) entries. The last entry records the total
     number of free variables plus 1 */
  int64_t            *fvnid_count_table;
  int64_t            *fvnid_interval_table;
  int64_t            *all_fvnid_start;

  /* convenience variables pointing to the node table.
     Don't try to free the memory pointed to */
  trilinear_node_t   *anode_table;
  trilinear_node_t   *onode_table;
  trilinear_node_t   *dnode_table;

  /* convenience variables recording the total number of free variables,
     the starting id and the ending id. identical on all processors */
  int64_t             global_fvnid_num; /* total number of fvnids */
  int64_t             global_fvnid_start;       /* first global fvnid */
  int64_t             global_fvnid_end; /* last global fvnid */

  tick_t              bounds[3][2];
  tick_t              sizes[3], minsize, maxsize;
  double              ticksize;

  sc_MPI_Comm         mpicomm;
  int32_t             mpisize, mpirank;
  int32_t             recsize;

  /* geometry type and element and node patch ids */
  int8_t              gid;
  trilinear_mesh_pid_t *elem_pids;
  trilinear_mesh_pid_t *node_pids;

  /* placeholder pointers */
  void                (*destructor) (struct trilinear_mesh_t *);
  trilinear_mesh_extra_t *extra_info;

  /* this is used in p4est only and must not be touched */
  sc_mempool_t       *sharer_pool;      /* allocator for node sharers */
}
trilinear_mesh_t;

/**
 * octor_neighbor_t:
 *
 * The face neighbors are ordered in -x, +x, -y, +y, -z, and +z directions.
 *
 * For each direction, the neighbor(s) may be:
 *
 * 1. Out of the domain.
 *
 *      face_neighbor_eid[direction][0] = -1.
 *      face_neighbor_eid[direction][1--3] = -1.
 *
 * 2. As large or twice as large as the current element:
 *
 *      face_neighbor_eid[direction][0] = index
 *
 *    where ((index >= 0) && (index < local_elem_num)) if the neighbor is
 *    LOCAL, or
 *    ((index >= local_elem_num) &&
 *    (index < (local_elem_num + phantom_elem_num))) if the neighor is
 *    REMOTE.
 *
 *
 *      face_neighbor_eid[direction][1--3] = -1.
 *
 * 3. Half as large as the current element:
 *
 *      face_neighbor_eid[direction][i] = index
 *
 *    where index is defined as above. Note that in this case all four
 *    neighbors (half as large) must exist.
 *
 */

typedef struct octor_neighor_t
{
  int32_t             face_neighbor_eid[6][4];
}
octor_neighbor_t;

typedef struct phantom_elem_t
{
  tick_t              lx, ly, lz;
  tick_t              size;
  int32_t             owner_procid;     /* remote processor id */
  int32_t             reid;     /* element index on the remote processor  */
}
phantom_elem_t;

typedef struct octor_neighborhood_t
{
  int32_t             phantom_elem_num;
  phantom_elem_t     *phantom_elem_table;

  octor_neighbor_t   *local_elem_neighbor_table;
}
octor_neighborhood_t;

#endif /* !TRILINEAR_MESH_TYPES_H */

/*
 * END verbatim copy of trilinear_mesh_types.h from the Rhea code
 */

SC_EXTERN_C_BEGIN;

/** Creates a trilinear mesh structure from a p8est and its node data.
 */
trilinear_mesh_t   *p8est_trilinear_mesh_new_from_nodes (p8est_t * p8est,
                                                         p8est_nodes_t *
                                                         nodes);

/** Creates a trilinear mesh structure from a p8est and its lnode data.
 */
trilinear_mesh_t   *p8est_trilinear_mesh_new_from_lnodes (p8est_t * p8est,
                                                          p8est_lnodes_t *
                                                          lnodes);
/** Frees a trilinear mesh structure.
 */
void                p8est_trilinear_mesh_destroy (trilinear_mesh_t * mesh);

SC_EXTERN_C_END;

#endif /* !P8EST_TRILINEAR_H */
