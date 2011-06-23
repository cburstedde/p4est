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
p8est_tet_edge_to_corner[6][2] =
{{0, 1},
 {0, 2},
 {0, 3},
 {1, 2},
 {1, 3},
 {2, 3}};

static const int
p8est_tet_face_to_corner[4][3] =
{{0, 1, 2},
 {0, 1, 3},
 {0, 2, 3},
 {1, 2, 3}};
/* *INDENT-ON* */

static int
p8est_topidx_is_sorted (p4est_topidx_t * t, int length)
{
  int                 i;

  for (i = 1; i < length; ++i) {
    if (t[i - 1] > t[i]) {
      return 0;
    }
  }
  return 1;
}

static void
p8est_topidx_bsort (p4est_topidx_t * t, int length)
{
  int                 i, j;
  p4est_topidx_t      tswap;

  /* go through all elements except the last */
  for (i = length - 1; i > 0; --i) {
    /* bubble up the first element until before position i */
    for (j = 0; j < i; ++j) {
      if (t[j] > t[j + 1]) {
        tswap = t[j + 1];
        t[j + 1] = t[j];
        t[j] = tswap;
      }
    }
  }
  P4EST_ASSERT (p8est_topidx_is_sorted (t, length));
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

/** Create unique edge key for a given edge of a tetrahedron.
 * \param [out] ek      The edge key consists of two node numbers.
 * \param [in] tet      A tetrahedron referring to node indices.
 * \param [in] edge     Tetrahedron edge number in [ 0 .. 5 ].
 */
static void
p8est_tet_edge_key (p4est_topidx_t * ek, p4est_topidx_t * tet, int edge)
{
  P4EST_ASSERT (0 <= edge && edge < 6);

  ek[0] = tet[p8est_tet_edge_to_corner[edge][0]];
  ek[1] = tet[p8est_tet_edge_to_corner[edge][1]];

  P4EST_ASSERT (ek[0] != ek[1]);
  p8est_topidx_bsort (ek, 2);
}

/** Create unique face key for a given face of a tetrahedron.
 * \param [out] fk      The edge key consists of three node numbers.
 * \param [in] tet      A tetrahedron referring to node indices.
 * \param [in] face     Tetrahedron face number in [ 0 .. 3 ].
 */
static void
p8est_tet_face_key (p4est_topidx_t * fk, p4est_topidx_t * tet, int face)
{
  P4EST_ASSERT (0 <= face && face < 4);

  fk[0] = tet[p8est_tet_face_to_corner[face][0]];
  fk[1] = tet[p8est_tet_face_to_corner[face][1]];
  fk[2] = tet[p8est_tet_face_to_corner[face][2]];

  P4EST_ASSERT (fk[0] != fk[1] && fk[0] != fk[2] && fk[1] != fk[2]);
  p8est_topidx_bsort (fk, 3);
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
    tet = (p4est_topidx_t *) sc_array_index (ptg->tets, 4 * iz);
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

static unsigned
p8est_tet_edge_hash (const void *v, const void *u)
{
  const p8est_tet_edge_info_t *ei = (p8est_tet_edge_info_t *) v;
  uint32_t            a, b, c;

#if (P4EST_TOPIDX_FITS_32)
  a = (uint32_t) ei->ek[0];
  b = (uint32_t) ei->ek[1];
  c = 0;
#else
  a = (uint32_t) (ei->ek[0] && 0xFFFFFFFF);
  b = (uint32_t) (ei->ek[0] >> 32);
  c = (uint32_t) (ei->ek[1] && 0xFFFFFFFF);
  sc_hash_mix (a, b, c);
  a += (uint32_t) (ei->ek[1] >> 32);
#endif
  sc_hash_final (a, b, c);

  return (unsigned) c;
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
    tet = (p4est_topidx_t *) sc_array_index (ptg->tets, 4 * iz);
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

static unsigned
p8est_tet_face_hash (const void *v, const void *u)
{
  const p8est_tet_face_info_t *fi = (p8est_tet_face_info_t *) v;
  uint32_t            a, b, c;

#if (P4EST_TOPIDX_FITS_32)
  a = (uint32_t) fi->fk[0];
  b = (uint32_t) fi->fk[1];
  c = (uint32_t) fi->fk[2];
#else
  a = (uint32_t) (fi->fk[0] && 0xFFFFFFFF);
  b = (uint32_t) (fi->fk[0] >> 32);
  c = (uint32_t) (fi->fk[1] && 0xFFFFFFFF);
  sc_hash_mix (a, b, c);
  a += (uint32_t) (fi->fk[1] >> 32);
  b += (uint32_t) (fi->fk[2] && 0xFFFFFFFF);
  c += (uint32_t) (fi->fk[2] >> 32);
#endif
  sc_hash_final (a, b, c);

  return (unsigned) c;
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
    tet = (p4est_topidx_t *) sc_array_index (ptg->tets, 4 * iz);
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

p8est_connectivity_t *
p8est_connectivity_new_tets (p8est_tets_t * ptg)
{
  size_t              ez, znum_edges;
  sc_hash_array_t    *edge_ha, *face_ha;
  sc_array_t          edge_array;
  p8est_tet_edge_info_t *ei;

  /* identify unique edges and faces */
  edge_ha = p8est_tets_identify_edges (ptg);
  P4EST_GLOBAL_LDEBUGF ("Added %ld unique tetrahedron edges\n",
                        (long) edge_ha->a.elem_count);

  face_ha = p8est_tets_identify_faces (ptg);
  P4EST_GLOBAL_LDEBUGF ("Added %ld unique tetrahedron faces\n",
                        (long) face_ha->a.elem_count);

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

  return NULL;
}
