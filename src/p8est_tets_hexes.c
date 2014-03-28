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

#include <p8est_tets_hexes.h>

/* *INDENT-OFF* */
static const int
p8est_tet_edge_corners[6][2] =
{{0, 1},
 {0, 2},
 {0, 3},
 {1, 2},
 {1, 3},
 {2, 3}};

static const int
p8est_tet_face_corners[4][3] =
{{0, 1, 2},
 {0, 1, 3},
 {0, 2, 3},
 {1, 2, 3}};

static const int
p8est_tet_tree_nodes[4][8] =
{{ 0, 4, 5, 10, 6, 11, 12, 14 },
 { 4, 1, 10, 7, 11, 8, 14, 13 },
 { 5, 10, 2, 7, 12, 14, 9, 13 },
 { 6, 11, 12, 14, 3, 8, 9, 13 }};
/* *INDENT-ON* */

static inline double *
p8est_tets_node_index (p8est_tets_t * ptg, size_t nodeno)
{
  return (double *) sc_array_index (ptg->nodes, 3 * nodeno);
}

static inline p4est_topidx_t *
p8est_tets_tet_index (p8est_tets_t * ptg, size_t tetno)
{
  return (p4est_topidx_t *) sc_array_index (ptg->tets, 4 * tetno);
}

sc_array_t         *
p8est_tets_read_node (const char *nodefilename)
{
  int                 retval;
  int                 k;
  int                 dims, num_attributes, boundary_marker;
  long                jl, lnum_nodes;
  double             *pc;
  size_t              iz, znum_nodes;
  FILE               *nodefile;
  sc_array_t         *nodes;

  /* prepare cleanup on error */
  nodefile = NULL;
  nodes = NULL;

  /* open node file */
  nodefile = fopen (nodefilename, "rb");
  if (nodefile == NULL) {
    P4EST_LERRORF ("Failed to open %s\n", nodefilename);
    goto dead;
  }

  /* read header information */
  retval = fscanf (nodefile, "%ld %d %d %d",
                   &lnum_nodes, &dims, &num_attributes, &boundary_marker) - 4;
  if (retval || lnum_nodes < 0 || lnum_nodes > P4EST_TOPIDX_MAX
      || dims != 3 || num_attributes < 0) {
    P4EST_LERROR ("Failed to read node header\n");
    goto dead;
  }
  znum_nodes = (size_t) lnum_nodes;

  /* read node coordinates */
  nodes = sc_array_new_size (sizeof (double), 3 * znum_nodes);
  for (iz = 0; iz < znum_nodes; ++iz) {
    pc = (double *) sc_array_index (nodes, 3 * iz);
    retval = fscanf (nodefile, "%ld %lf %lf %lf",
                     &jl, pc, pc + 1, pc + 2) - 4;
    if (retval || (long) iz != jl) {
      P4EST_LERRORF ("Failed to read node %ld coordinates\n", (long) iz);
      goto dead;
    }
    for (k = 0; k < num_attributes; ++k) {
      retval = fscanf (nodefile, "%*f");
      if (retval != 0) {
        P4EST_LERRORF ("Failed to read node %ld attribute %d\n", (long) iz,
                       k);
        goto dead;
      }
    }
    if (boundary_marker) {
      retval = fscanf (nodefile, "%*d");
      if (retval != 0) {
        P4EST_LERRORF ("Failed to read node %ld boundary marker\n",
                       (long) iz);
        goto dead;
      }
    }
  }

  /* close node file and return */
  retval = fclose (nodefile);
  nodefile = NULL;
  if (retval) {
    P4EST_LERRORF ("Failed to close %s\n", nodefilename);
    goto dead;
  }
  return nodes;

dead:
  /* clean up on error */
  if (nodefile != NULL) {
    fclose (nodefile);
  }
  if (nodes != NULL) {
    sc_array_destroy (nodes);
  }
  return NULL;
}

sc_array_t         *
p8est_tets_read_ele (const char *elefilename, p4est_topidx_t num_nodes,
                     sc_array_t ** attributes)
{
  int                 retval;
  int                 k;
  int                 nodespertet, region;
  long                jl, lnum_tets, lmax_nodes;
  long                nl[4];
  size_t              iz, znum_tets;
  int                *pi;
  p4est_topidx_t     *pt;
  FILE               *elefile;
  sc_array_t         *tets, *attr;

  /* prepare cleanup on error */
  elefile = NULL;
  tets = NULL;
  attr = NULL;
  if (attributes != NULL) {
    *attributes = NULL;
  }
  lmax_nodes = (num_nodes >= 0 ? (long) num_nodes : P4EST_TOPIDX_MAX);

  /* open ele file */
  elefile = fopen (elefilename, "rb");
  if (elefile == NULL) {
    P4EST_LERRORF ("Failed to open %s\n", elefilename);
    goto dead;
  }

  /* read header information */
  retval = fscanf (elefile, "%ld %d %d",
                   &lnum_tets, &nodespertet, &region) - 3;
  if (retval || lnum_tets < 0 || lnum_tets > P4EST_TOPIDX_MAX
      || nodespertet != 4) {
    P4EST_LERROR ("Failed to read tet header\n");
    goto dead;
  }
  znum_tets = (size_t) lnum_tets;

  /* read tet coordinates */
  tets = sc_array_new_size (sizeof (p4est_topidx_t), 4 * znum_tets);
  if (region && attributes != NULL) {
    attr = *attributes = sc_array_new_size (sizeof (int), znum_tets);
  }
  for (iz = 0; iz < znum_tets; ++iz) {
    pt = (p4est_topidx_t *) sc_array_index (tets, 4 * iz);
    retval = fscanf (elefile, "%ld %ld %ld %ld %ld",
                     &jl, nl, nl + 1, nl + 2, nl + 3) - 5;
    if (retval || (long) iz != jl) {
      P4EST_LERRORF ("Failed to read tet %ld node numbers\n", (long) iz);
      goto dead;
    }
    for (k = 0; k < 4; ++k) {
      if (nl[k] < 0 || nl[k] > lmax_nodes) {
        P4EST_LERRORF ("Tet %ld has invalid node number %d\n", (long) iz, k);
        goto dead;
      }
      pt[k] = (p4est_topidx_t) nl[k];
    }
    if (region) {
      if (attr != NULL) {
        pi = (int *) sc_array_index (attr, iz);
        retval = fscanf (elefile, "%d", pi) - 1;
      }
      else {
        retval = fscanf (elefile, "%*d");
      }
      if (retval != 0) {
        P4EST_LERRORF ("Failed to read tet %ld region attribute\n",
                       (long) iz);
        goto dead;
      }
    }
  }

  /* close ele file */
  retval = fclose (elefile);
  elefile = NULL;
  if (retval) {
    P4EST_LERRORF ("Failed to close %s\n", elefilename);
    goto dead;
  }
  return tets;

dead:
  /* clean up on error */
  if (elefile != NULL) {
    fclose (elefile);
  }
  if (tets != NULL) {
    sc_array_destroy (tets);
  }
  if (attr != NULL) {
    sc_array_destroy (attr);
    *attributes = NULL;
  }
  return NULL;
}

p8est_tets_t       *
p8est_tets_read (const char *tetgenbasename)
{
  p4est_topidx_t      num_nodes;
  char                nodefilename[BUFSIZ];
  char                elefilename[BUFSIZ];
  sc_array_t         *nodes, *tets, *attr;
  p8est_tets_t       *ptg;

  /* prepare cleanup */
  nodes = tets = attr = NULL;
  ptg = P4EST_ALLOC (p8est_tets_t, 1);

  /* read nodes */
  snprintf (nodefilename, BUFSIZ, "%s.node", tetgenbasename);
  nodes = ptg->nodes = p8est_tets_read_node (nodefilename);
  if (nodes == NULL) {
    P4EST_LERRORF ("Failed to read nodes for %s\n", tetgenbasename);
    goto dead;
  }
  num_nodes = (p4est_topidx_t) (nodes->elem_count / 3);

  /* read tetrahedra */
  snprintf (elefilename, BUFSIZ, "%s.ele", tetgenbasename);
  tets = ptg->tets = p8est_tets_read_ele (elefilename, num_nodes, &attr);
  if (tets == NULL) {
    P4EST_ASSERT (attr == NULL);
    P4EST_LERRORF ("Failed to read tetrahedra for %s\n", tetgenbasename);
    goto dead;
  }
  ptg->tet_attributes = attr;

  /* we are successful */
  return ptg;

dead:
  /* clean up on error */
  if (nodes != NULL) {
    sc_array_destroy (nodes);
  }
  if (tets != NULL) {
    sc_array_destroy (tets);
  }
  if (attr != NULL) {
    sc_array_destroy (attr);
  }
  P4EST_FREE (ptg);
  return NULL;
}

void
p8est_tets_destroy (p8est_tets_t * ptg)
{
  sc_array_destroy (ptg->nodes);
  sc_array_destroy (ptg->tets);
  if (ptg->tet_attributes != NULL) {
    sc_array_destroy (ptg->tet_attributes);
  }

  P4EST_FREE (ptg);
}

static int
p8est_tet_is_righthanded (sc_array_t * nodes, p4est_topidx_t * tet)
{
  int                 i, j;
  double             *nc[4];
  double              v0[3], v1[3], v2[3], cross01[3];
  double              vol;

  /* compute tet volume */
  for (i = 0; i < 4; ++i) {
    nc[i] = (double *) sc_array_index (nodes, (size_t) (3 * tet[i]));
  }
  for (j = 0; j < 3; ++j) {
    v0[j] = nc[1][j] - nc[0][j];
    v1[j] = nc[2][j] - nc[0][j];
    v2[j] = nc[3][j] - nc[0][j];
  }
  cross01[0] = v0[1] * v1[2] - v0[2] * v1[1];
  cross01[1] = v0[2] * v1[0] - v0[0] * v1[2];
  cross01[2] = v0[0] * v1[1] - v0[1] * v1[0];
  vol = 0.;
  for (j = 0; j < 3; ++j) {
    vol += cross01[j] * v2[j];
  }
  vol *= 1. / 3.;

  return vol >= 0.;
}

static void
p8est_tet_flip (p4est_topidx_t * tet)
{
  p4est_topidx_t      temp;

  temp = tet[3];
  tet[3] = tet[2];
  tet[2] = temp;
}

p4est_topidx_t
p8est_tets_make_righthanded (p8est_tets_t * ptg)
{
  size_t              iz, znum_tets;
  p4est_topidx_t      tnum_flips;
  p4est_topidx_t     *tet;

  tnum_flips = 0;
  znum_tets = ptg->tets->elem_count / 4;
  for (iz = 0; iz < znum_tets; ++iz) {
    tet = p8est_tets_tet_index (ptg, iz);
    if (!p8est_tet_is_righthanded (ptg->nodes, tet)) {
      p8est_tet_flip (tet);
      P4EST_ASSERT (p8est_tet_is_righthanded (ptg->nodes, tet));
      ++tnum_flips;
    }
  }

  return tnum_flips;
}

typedef struct p8est_tet_edge_info
{
  p4est_topidx_t      ek[2];
  sc_array_t          tets;
  sc_array_t          tet_edges;
}
p8est_tet_edge_info_t;

/** Create unique edge key for a given edge of a tetrahedron.
 * \param [out] ek      The edge key consists of two node numbers.
 * \param [in] tet      A tetrahedron referring to node indices.
 * \param [in] edge     Tetrahedron edge number in [ 0 .. 5 ].
 */
static void
p8est_tet_edge_key (p4est_topidx_t * ek, p4est_topidx_t * tet, int edge)
{
  P4EST_ASSERT (0 <= edge && edge < 6);

  ek[0] = tet[p8est_tet_edge_corners[edge][0]];
  ek[1] = tet[p8est_tet_edge_corners[edge][1]];

  P4EST_ASSERT (ek[0] != ek[1]);
  p4est_topidx_bsort (ek, 2);
}

static unsigned
p8est_tet_edge_hash (const void *v, const void *u)
{
  const p8est_tet_edge_info_t *ei = (p8est_tet_edge_info_t *) v;

  return p4est_topidx_hash2 (ei->ek);
}

static int
p8est_tet_edge_equal (const void *v1, const void *v2, const void *u)
{
  const p8est_tet_edge_info_t *ei1 = (p8est_tet_edge_info_t *) v1;
  const p8est_tet_edge_info_t *ei2 = (p8est_tet_edge_info_t *) v2;

  return !memcmp (ei1->ek, ei2->ek, 2 * sizeof (p4est_topidx_t));
}

static sc_hash_array_t *
p8est_tets_identify_edges (p8est_tets_t * ptg)
{
  int                 edge, *pi;
  size_t              iz, znum_tets, pz;
  sc_hash_array_t    *edge_ha;
  p4est_topidx_t     *tet, *pt;
  p8est_tet_edge_info_t eikey, *ei;

  /* create hash array for shared edges */
  edge_ha = sc_hash_array_new (sizeof (p8est_tet_edge_info_t),
                               p8est_tet_edge_hash, p8est_tet_edge_equal,
                               NULL);

  /* loop through all edges and identify edge-neighbor tet groups */
  znum_tets = ptg->tets->elem_count / 4;
  for (iz = 0; iz < znum_tets; ++iz) {
    tet = p8est_tets_tet_index (ptg, iz);
    for (edge = 0; edge < 6; ++edge) {
      p8est_tet_edge_key (eikey.ek, tet, edge);
      ei = (p8est_tet_edge_info_t *)
        sc_hash_array_insert_unique (edge_ha, &eikey, &pz);
      if (ei != NULL) {
        /* added new edge group ei to hash array */
        P4EST_ASSERT (sc_array_position (&edge_ha->a, ei) == pz);
        memcpy (ei->ek, eikey.ek, 2 * sizeof (p4est_topidx_t));
        sc_array_init (&ei->tets, sizeof (p4est_topidx_t));
        pt = (p4est_topidx_t *) sc_array_push (&ei->tets);
        *pt = (p4est_topidx_t) iz;
        sc_array_init (&ei->tet_edges, sizeof (int));
        pi = (int *) sc_array_push (&ei->tet_edges);
        *pi = edge;
      }
      else {
        /* found existing entry from earlier edge */
        ei = (p8est_tet_edge_info_t *) sc_array_index (&edge_ha->a, pz);
        P4EST_ASSERT (p8est_tet_edge_equal (ei->ek, eikey.ek, NULL));
        P4EST_ASSERT (ei->tets.elem_count > 0);
        pt = (p4est_topidx_t *) sc_array_push (&ei->tets);
        *pt = (p4est_topidx_t) iz;
        P4EST_ASSERT (ei->tet_edges.elem_count > 0);
        pi = (int *) sc_array_push (&ei->tet_edges);
        *pi = edge;
      }
    }
  }

  return edge_ha;
}

typedef struct p8est_tet_face_info
{
  p4est_topidx_t      fk[3];
  p4est_topidx_t      tets[2];
  int                 tet_faces[2];
}
p8est_tet_face_info_t;

/** Create unique face key for a given face of a tetrahedron.
 * \param [out] fk      The face key consists of three node numbers.
 * \param [in] tet      A tetrahedron referring to node indices.
 * \param [in] face     Tetrahedron face number in [ 0 .. 3 ].
 */
static void
p8est_tet_face_key (p4est_topidx_t * fk, p4est_topidx_t * tet, int face)
{
  P4EST_ASSERT (0 <= face && face < 4);

  fk[0] = tet[p8est_tet_face_corners[face][0]];
  fk[1] = tet[p8est_tet_face_corners[face][1]];
  fk[2] = tet[p8est_tet_face_corners[face][2]];

  P4EST_ASSERT (fk[0] != fk[1] && fk[0] != fk[2] && fk[1] != fk[2]);
  p4est_topidx_bsort (fk, 3);
}

static unsigned
p8est_tet_face_hash (const void *v, const void *u)
{
  const p8est_tet_face_info_t *fi = (p8est_tet_face_info_t *) v;

  return p4est_topidx_hash3 (fi->fk);
}

static int
p8est_tet_face_equal (const void *v1, const void *v2, const void *u)
{
  const p8est_tet_face_info_t *fi1 = (p8est_tet_face_info_t *) v1;
  const p8est_tet_face_info_t *fi2 = (p8est_tet_face_info_t *) v2;

  return !memcmp (fi1->fk, fi2->fk, 3 * sizeof (p4est_topidx_t));
}

static sc_hash_array_t *
p8est_tets_identify_faces (p8est_tets_t * ptg)
{
  int                 face;
  size_t              iz, znum_tets, pz;
  sc_hash_array_t    *face_ha;
  p4est_topidx_t     *tet;
  p8est_tet_face_info_t fikey, *fi;

  /* create hash array for shared faces */
  face_ha = sc_hash_array_new (sizeof (p8est_tet_face_info_t),
                               p8est_tet_face_hash, p8est_tet_face_equal,
                               NULL);

  /* loop through all faces and identify face-neighbor tet pairs */
  znum_tets = ptg->tets->elem_count / 4;
  for (iz = 0; iz < znum_tets; ++iz) {
    tet = p8est_tets_tet_index (ptg, iz);
    for (face = 0; face < 4; ++face) {
      p8est_tet_face_key (fikey.fk, tet, face);
      fi = (p8est_tet_face_info_t *)
        sc_hash_array_insert_unique (face_ha, &fikey, &pz);
      if (fi != NULL) {
        /* added fi to hash array as the first of two tets */
        P4EST_ASSERT (sc_array_position (&face_ha->a, fi) == pz);
        memcpy (fi->fk, fikey.fk, 3 * sizeof (p4est_topidx_t));
        fi->tets[0] = (p4est_topidx_t) iz;
        fi->tets[1] = -1;
        fi->tet_faces[0] = face;
        fi->tet_faces[1] = -1;
      }
      else {
        /* found existing entry from the first face */
        fi = (p8est_tet_face_info_t *) sc_array_index (&face_ha->a, pz);
        P4EST_ASSERT (p8est_tet_face_equal (fi->fk, fikey.fk, NULL));
        P4EST_ASSERT (fi->tets[0] >= 0 && fi->tet_faces[0] >= 0);
        P4EST_ASSERT (fi->tets[1] == -1 && fi->tet_faces[1] == -1);
        fi->tets[1] = (p4est_topidx_t) iz;
        fi->tet_faces[1] = face;
      }
    }
  }

  return face_ha;
}

/** Create a connectivity where the trees are not connected to each other. */
static p8est_connectivity_t *
p8est_tets_connectivity_new (p8est_tets_t * ptg,
                             sc_hash_array_t * edge_ha,
                             sc_hash_array_t * face_ha)
{
  int                 j, k;
  int                 edge, face;
  size_t              nvz, evzoffset, fvzoffset, vvzoffset;
  size_t              iz, pz;
  double             *vp, *n[4];
  int8_t             *ttf;
  p4est_topidx_t      tt, *tet, node;
  p4est_topidx_t     *ttv, *ttt;
  p4est_topidx_t      nid[15];
  p8est_connectivity_t *conn;
  p8est_tet_edge_info_t *ei, eikey;
  p8est_tet_face_info_t *fi, fikey;

  /* arrange vertices by tet corners, edges, faces, and volumes */
  evzoffset = ptg->nodes->elem_count / 3;
  fvzoffset = evzoffset + edge_ha->a.elem_count;
  vvzoffset = fvzoffset + face_ha->a.elem_count;
  nvz = vvzoffset + ptg->tets->elem_count / 4;

  /* allocate connectivity */
  conn = p8est_connectivity_new (nvz, ptg->tets->elem_count, 0, 0, 0, 0);

  /* populate vertices */
  memcpy (conn->vertices, ptg->nodes->array, 3 * evzoffset * sizeof (double));
  vp = conn->vertices + 3 * evzoffset;
  for (iz = 0; iz < edge_ha->a.elem_count; ++iz) {
    ei = (p8est_tet_edge_info_t *) sc_array_index (&edge_ha->a, iz);
    tt = *((p4est_topidx_t *) sc_array_index (&ei->tets, 0));
    edge = *((int *) sc_array_index (&ei->tet_edges, 0));
    tet = p8est_tets_tet_index (ptg, tt);
    for (j = 0; j < 2; ++j) {
      node = tet[p8est_tet_edge_corners[edge][j]];
      n[j] = p8est_tets_node_index (ptg, node);
    }
    vp[0] = .5 * (n[0][0] + n[1][0]);
    vp[1] = .5 * (n[0][1] + n[1][1]);
    vp[2] = .5 * (n[0][2] + n[1][2]);
    vp += 3;
  }
  for (iz = 0; iz < face_ha->a.elem_count; ++iz) {
    fi = (p8est_tet_face_info_t *) sc_array_index (&face_ha->a, iz);
    tt = fi->tets[0];
    face = fi->tet_faces[0];
    tet = p8est_tets_tet_index (ptg, tt);
    for (j = 0; j < 3; ++j) {
      node = tet[p8est_tet_face_corners[face][j]];
      n[j] = p8est_tets_node_index (ptg, node);
    }
    vp[0] = (1. / 3.) * (n[0][0] + n[1][0] + n[2][0]);
    vp[1] = (1. / 3.) * (n[0][1] + n[1][1] + n[2][1]);
    vp[2] = (1. / 3.) * (n[0][2] + n[1][2] + n[2][2]);
    vp += 3;
  }
  for (iz = 0; iz < ptg->tets->elem_count / 4; ++iz) {
    tet = p8est_tets_tet_index (ptg, iz);
    for (j = 0; j < 4; ++j) {
      n[j] = p8est_tets_node_index (ptg, tet[j]);
    }
    vp[0] = .25 * (n[0][0] + n[1][0] + n[2][0] + n[3][0]);
    vp[1] = .25 * (n[0][1] + n[1][1] + n[2][1] + n[3][1]);
    vp[2] = .25 * (n[0][2] + n[1][2] + n[2][2] + n[3][2]);
    vp += 3;
  }

  /* associate forest trees with vertices */
  ttv = conn->tree_to_vertex;
  for (iz = 0; iz < ptg->tets->elem_count / 4; ++iz) {
    tet = p8est_tets_tet_index (ptg, iz);

    /* look up node numbers for all vertices in this tetrahedron */
    for (j = 0; j < 4; ++j) {
      nid[j] = tet[j];
    }
    for (edge = 0; edge < 6; ++edge) {
      p8est_tet_edge_key (eikey.ek, tet, edge);
      P4EST_EXECUTE_ASSERT_TRUE (sc_hash_array_lookup (edge_ha, &eikey, &pz));
      nid[4 + edge] = (p4est_topidx_t) (evzoffset + pz);
    }
    for (face = 0; face < 4; ++face) {
      p8est_tet_face_key (fikey.fk, tet, face);
      P4EST_EXECUTE_ASSERT_TRUE (sc_hash_array_lookup (face_ha, &fikey, &pz));
      nid[10 + face] = (p4est_topidx_t) (fvzoffset + pz);
    }
    nid[14] = (p4est_topidx_t) (vvzoffset + iz);

    /* create four trees from this tetrahedron */
    for (j = 0; j < 4; ++j) {
      for (k = 0; k < P8EST_CHILDREN; ++k) {
        *ttv++ = nid[p8est_tet_tree_nodes[j][k]];
      }
    }
  }

  /* create neighborhood information for isolated trees */
  ttt = conn->tree_to_tree;
  ttf = conn->tree_to_face;
  for (tt = 0; tt < conn->num_trees; ++tt) {
    for (face = 0; face < P8EST_FACES; ++face) {
      ttt[face] = tt;
      ttf[face] = (int8_t) face;
    }
    ttt += P8EST_FACES;
    ttf += P8EST_FACES;
  }

  return conn;
}

p8est_connectivity_t *
p8est_connectivity_new_tets (p8est_tets_t * ptg)
{
  int                *pint, i;
  int8_t              attr;
  size_t              ez, znum_edges;
  size_t              tz, znum_tets;
  sc_hash_array_t    *edge_ha, *face_ha;
  sc_array_t          edge_array;
  p8est_tet_edge_info_t *ei;
  p8est_connectivity_t *conn;

  /* identify unique edges and faces */
  edge_ha = p8est_tets_identify_edges (ptg);
  P4EST_GLOBAL_LDEBUGF ("Added %ld unique tetrahedron edges\n",
                        (long) edge_ha->a.elem_count);

  face_ha = p8est_tets_identify_faces (ptg);
  P4EST_GLOBAL_LDEBUGF ("Added %ld unique tetrahedron faces\n",
                        (long) face_ha->a.elem_count);

  /* add vertex information to connectivity */
  conn = p8est_tets_connectivity_new (ptg, edge_ha, face_ha);
  P4EST_GLOBAL_LDEBUGF ("Connectivity has %ld vertices and %ld trees\n",
                        (long) conn->num_vertices, (long) conn->num_trees);

  /* clean unique edges and faces */
  sc_hash_array_rip (edge_ha, &edge_array);
  znum_edges = edge_array.elem_count;
  for (ez = 0; ez < znum_edges; ++ez) {
    ei = (p8est_tet_edge_info_t *) sc_array_index (&edge_array, ez);
    sc_array_reset (&ei->tets);
    sc_array_reset (&ei->tet_edges);
  }
  sc_array_reset (&edge_array);
  sc_hash_array_destroy (face_ha);

  /* transfer tree tags */
  if (ptg->tet_attributes != NULL) {
    znum_tets = ptg->tet_attributes->elem_count;
    P4EST_ASSERT (4 * znum_tets == (size_t) conn->num_trees);
    p8est_connectivity_set_attr (conn, 1);
    for (tz = 0; tz < znum_tets; ++tz) {
      pint = (int *) sc_array_index (ptg->tet_attributes, tz);
      attr = (int8_t) pint[0];
      for (i = 0; i < 4; ++i) {
        conn->tree_to_attr[4 * tz + i] = attr;
      }
    }
  }

  /* connect p4est tree through faces, edges, and corners */
  p8est_connectivity_complete (conn);
  P4EST_GLOBAL_LDEBUGF ("Connectivity has %ld edges and %ld corners\n",
                        (long) conn->num_edges, (long) conn->num_corners);

  return conn;
}
