/*
  This file is part of p4est.
  p4est is a C library to manage a collection (a forest) of multiple
  connected adaptive quadtrees or octrees in parallel.

  Copyright (C) 2010 The University of Texas System
  Additional copyright (C) 2011 individual authors
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
#include <p4est_iterate.h>
#include <p4est_tnodes.h>
#else
#include <p8est_bits.h>
#include <p8est_iterate.h>
#include <p8est_tnodes.h>
#define p4est_tnodes_private            p8est_tnodes_private
#define p4est_tnodes_iter_private       p8est_tnodes_iter_private
#endif

#ifndef P4_TO_P8

#define P4EST_TNODES_MAXNE 25

/* *INDENT-OFF* */

/* cube corners */
static const int    n_cornr[4] = {  0,  1,  2,  3 };

/* cube center */
static const int    n_center =      4;

/* face midpoints */
static const int    n_mface[4] = {  5,  6,  7,  8 };

/* faces between center and corners */
static const int    n_cface[4] = {  9, 10, 11, 12 };

/* faces between center and face midpoints */
static const int    n_split[4] = { 14, 17, 20, 22 };

/* quarter cube faces */
static const int    n_hface[4][2] =
  {{ 13, 15 }, { 16, 18 }, { 19, 21 }, { 23, 24 }};

#ifdef P4EST_ENABLE_DEBUG
static const int    alwaysowned[25] =
  { 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1,
    0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0 };
#endif

const int p4est_tnodes_triangle_nodes[16][6] =
  {{ 0, 1, 3,  7,  6,  4 },     /* lower half of subconfig 0 */
   { 0, 3, 2,  4,  8,  5 },     /* upper half of subconfig 0 */
   { 0, 1, 2,  7,  4,  5 },     /* near half of subconfig 1 */
   { 1, 3, 2,  6,  8,  4 },     /* far half of subconfig 1 */
   { 0, 4, 2,  9, 11,  5 },     /* left triangle of subconfig 2 */
   { 1, 3, 4,  6, 12, 10 },     /* right triangle of subconfig 2 */
   { 0, 1, 4,  7, 10,  9 },     /* bottom triangle of subconfig 2 */
   { 2, 4, 3, 11, 12,  8 },     /* top triangle of subconfig 2 */
   { 0, 4, 5,  9, 14, 13 },     /* lower left center triangle */
   { 2, 5, 4, 15, 14, 11 },     /* upper left center triangle */
   { 1, 6, 4, 16, 17, 10 },     /* lower right center triangle */
   { 3, 4, 6, 12, 17, 18 },     /* upper right center triangle */
   { 0, 7, 4, 19, 20,  9 },     /* left bottom center triangle */
   { 1, 4, 7, 10, 20, 21 },     /* right bottom center triangle */
   { 2, 4, 8, 11, 22, 23 },     /* left top center triangle */
   { 3, 8, 4, 24, 22, 12 }};    /* right top center triangle */

const int p4est_tnodes_config_triangles[18][8] =
  {{ 0,  1, -1, -1, -1, -1, -1, -1, },
   { 8,  9,  5,  6,  7, -1, -1, -1, },
   { 4, 10, 11,  6,  7, -1, -1, -1, },  /*  2 */
   { 8,  9, 10, 11,  6,  7, -1, -1, },
   { 4,  5, 12, 13,  7, -1, -1, -1, },  /*  4 */
   { 8,  9,  5, 12, 13,  7, -1, -1, },
   { 4, 10, 11, 12, 13,  7, -1, -1, },
   { 8,  9, 10, 11, 12, 13,  7, -1, },
   { 4,  5,  6, 14, 15, -1, -1, -1, },  /*  8 */
   { 8,  9,  5,  6, 14, 15, -1, -1, },
   { 4, 10, 11,  6, 14, 15, -1, -1, },
   { 8,  9, 10, 11,  6, 14, 15, -1, },
   { 4,  5, 12, 13, 14, 15, -1, -1, },  /* 12 */
   { 8,  9,  5, 12, 13, 14, 15, -1, },
   { 4, 10, 11, 12, 13, 14, 15, -1, },
   { 8,  9, 10, 11, 12, 13, 14, 15, },  /* 15 */
   { 2,  3, -1, -1, -1, -1, -1, -1, },
   { 4,  5,  6,  7, -1, -1, -1, -1, }};

const int p4est_tnodes_lookup_counts[6][3] =
  {{ 4,  5, 2 },                        /* 0, subconfig 0, 1 */
   { 5,  8, 4 },                        /* 0, subconfig 2 */
   { 6, 10, 5 },                        /* 1, 2, 4, 8) */
   { 7, 12, 6 },                        /* 3, 5, 6, 9, 10, 12 */
   { 8, 14, 7 },                        /* 7, 11, 13, 14 */
   { 9, 16, 8 }};                       /* 15 */

const int p4est_tnodes_config_lookup[18] =
  { 0,                                 /* 0, subconfig 0 */
    2, 2,                              /* 1, 2 (rotated: 4, 8) */
    3,                                 /* 3 (rotated: 12) */
    2,                                 /* 4 (see 0, 1, 8) */
    3, 3,                              /* 5, 6 (rotated: 9, 10) */
    4,                                 /* 7 (rotated: 11, 13, 14) */
    2,                                 /* 8 (see 1, 2, 4) */
    3, 3,                              /* 9, 10 (see 5, 6) */
    4,                                 /* 11 (see: 7, 13, 14) */
    3,                                 /* 12 (see: 3) */
    4, 4,                              /* 13, 14 (see: 7, 11) */
    5,                                 /* 15 */
    0,                                 /* 0, subconfig 1 */
    1 };                               /* 0, subconfig 2 */

const int p4est_tnodes_config_corners[18][9] =
  {{ 0, 1, 2, 3, -1, -1, -1, -1, -1 },
   { 0, 1, 2, 3,  4,  5, -1, -1, -1 },  /*  1 */
   { 0, 1, 2, 3,  4,  6, -1, -1, -1 },  /*  2 */
   { 0, 1, 2, 3,  4,  5,  6, -1, -1 },
   { 0, 1, 2, 3,  4,  7, -1, -1, -1 },  /*  4 */
   { 0, 1, 2, 3,  4,  5,  7, -1, -1 },
   { 0, 1, 2, 3,  4,  6,  7, -1, -1 },
   { 0, 1, 2, 3,  4,  5,  6,  7, -1 },
   { 0, 1, 2, 3,  4,  8, -1, -1, -1 },  /*  8 */
   { 0, 1, 2, 3,  4,  5,  8, -1, -1 },
   { 0, 1, 2, 3,  4,  6,  8, -1, -1 },  /* 10 */
   { 0, 1, 2, 3,  4,  5,  6,  8, -1 },
   { 0, 1, 2, 3,  4,  7,  8, -1, -1 },  /* 12 */
   { 0, 1, 2, 3,  4,  5,  7,  8, -1 },
   { 0, 1, 2, 3,  4,  6,  7,  8, -1 },
   { 0, 1, 2, 3,  4,  5,  6,  7,  8 },  /* 15 */
   { 0, 1, 2, 3, -1, -1, -1, -1, -1 },
   { 0, 1, 2, 3,  4, -1, -1, -1, -1 }};

/** For each configuration the list of face nodes padded with -1. */
const int p4est_tnodes_config_faces[18][16] = {
  { 4,  5,  6,  7,   8, -1, -1, -1,   -1, -1, -1, -1,  -1, -1, -1, -1 },
  { 6,  7,  8,  9,  10, 11, 12, 13,   14, 15, -1, -1,  -1, -1, -1, -1 }, /*  1 */
  { 5,  7,  8,  9,  10, 11, 12, 16,   17, 18, -1, -1,  -1, -1, -1, -1 }, /*  2 */
  { 7,  8,  9, 10,  11, 12, 13, 14,   15, 16, 17, 18,  -1, -1, -1, -1 },
  { 5,  6,  8,  9,  10, 11, 12, 19,   20, 21, -1, -1,  -1, -1, -1, -1 }, /*  4 */
  { 6,  8,  9, 10,  11, 12, 13, 14,   15, 19, 20, 21,  -1, -1, -1, -1 },
  { 5,  8,  9, 10,  11, 12, 16, 17,   18, 19, 20, 21,  -1, -1, -1, -1 },
  { 8,  9, 10, 11,  12, 13, 14, 15,   16, 17, 18, 19,  20, 21, -1, -1 },
  { 5,  6,  7,  9,  10, 11, 12, 22,   23, 24, -1, -1,  -1, -1, -1, -1 }, /*  8 */
  { 6,  7,  9, 10,  11, 12, 13, 14,   15, 22, 23, 24,  -1, -1, -1, -1 },
  { 5,  7,  9, 10,  11, 12, 16, 17,   18, 22, 23, 24,  -1, -1, -1, -1 }, /* 10 */
  { 7,  9, 10, 11,  12, 13, 14, 15,   16, 17, 18, 22,  23, 24, -1, -1 },
  { 5,  6,  9, 10,  11, 12, 19, 20,   21, 22, 23, 24,  -1, -1, -1, -1 }, /* 12 */
  { 6,  9, 10, 11,  12, 13, 14, 15,   19, 20, 21, 22,  23, 24, -1, -1 },
  { 5,  9, 10, 11,  12, 16, 17, 18,   19, 20, 21, 22,  23, 24, -1, -1 },
  { 9, 10, 11, 12,  13, 14, 15, 16,   17, 18, 19, 20,  21, 22, 23, 24 }, /* 15 */
  { 4,  5,  6,  7,   8, -1, -1, -1,   -1, -1, -1, -1,  -1, -1, -1, -1 },
  { 5,  6,  7,  8,   9, 10, 11, 12,   -1, -1, -1, -1,  -1, -1, -1, -1 }};
/* *INDENT-ON* */

typedef struct p4est_tnodes_iter_private
{
  p4est_tree_t       *tree;         /**< Pointer to current tree. */
  p4est_locidx_t      numtreeq;     /**< Number of quadrants in tree. */
  p4est_locidx_t      treequad;     /**< Quadrant within local tree. */
  p4est_locidx_t      numtris;      /**< Nunmber of local triangles. */
  int                 numqtri;      /**< Number triangles within quadrant. */
  int                 quadtri;      /**< Triangle within quadrant. */
  int                 cind;         /**< Current configuration index. */
}
p4est_tnodes_iter_private_t;

#endif /* !P4_TO_P8 */

/** Normalized unit length for element node coordinate */
#define P4EST_TNODES_ESHIFT (1 << 4)

/** Number of node coordinates (high end included) */
#define P4EST_TNODES_FACTOR (P4EST_TNODES_ESHIFT + 1)

/** Number of elementary forest roots. */
#define P4EST_TNODES_NUM_SROOTS ((P4EST_DIM - 1) * P4EST_DIM)

#ifndef P4_TO_P8

/** Number of corners in a simplex. */
#define P4EST_TNODES_NUM_SCORNERS 3

/** Number of edges of a simplex. */
#define P4EST_TNODES_NUM_SEDGES 3

/** Number of leaves in deepest simplex refinement */
#define P4EST_TNODES_NUM_LEAVES 8

/** Number of nodes in simplex refinement forest */
#define P4EST_TNODES_NUM_SIMPLICES 14

/** High end (exclusive) of element node index in 2D */
#define P4EST_TNODES_ERANGE \
  (P4EST_TNODES_FACTOR * P4EST_TNODES_FACTOR)

#else

/** Number of corners in a simplex. */
#define P4EST_TNODES_NUM_SCORNERS 4

/** Number of edges of a simplex. */
#define P4EST_TNODES_NUM_SEDGES 6

/** Number of leaves in deepest simplex refinement */
#define P4EST_TNODES_NUM_LEAVES 48

/** Number of nodes in simplex refinement forest */
#define P4EST_TNODES_NUM_SIMPLICES 90

/** High end (exclusive) of element node index in 3D */
#define P4EST_TNODES_ERANGE \
  (P4EST_TNODES_FACTOR * P4EST_TNODES_FACTOR * P4EST_TNODES_FACTOR)

#endif

#define P4EST_TNODES_IS_ECO(eco)                                \
  (0 <= (eco) && (eco) < P4EST_TNODES_FACTOR)

#define P4EST_TNODES_IS_EIN(ein)                                \
  (0 <= (ein) && (ein) < P4EST_TNODES_ERANGE)

#ifndef P4_TO_P8

/* Transform cube corner number into element node index */
#define P4EST_TNODES_CTOEIN(c)                                  \
  (P4EST_TNODES_FACTOR * (((c) & 2) ? P4EST_TNODES_ESHIFT : 0) +        \
                         (((c) & 1) ? P4EST_TNODES_ESHIFT : 0))

#else

/* Transform cube corner number into element node index */
#define P4EST_TNODES_CTOEIN(c) (                                        \
  (((c) & 1) ? P4EST_TNODES_ESHIFT : 0) +                               \
  (                                                                     \
    (((c) & 2) ? P4EST_TNODES_ESHIFT : 0) +                             \
    (((c) & 4) ? P4EST_TNODES_ESHIFT : 0) * P4EST_TNODES_FACTOR         \
  ) * P4EST_TNODES_FACTOR                                               \
)

#endif /* P4_TO_P8 */

#ifndef P4_TO_P8

/** All edges of a simplex by their edge corners */
static const int    p4est_tnodes_sedge[3][2] = {
  {0, 1},
  {0, 2},
  {1, 2}
};

/** Cube corner numbers for every root simplex */
static const int    p4est_tnodes_rsim[2][3] = {
  {0, 1, 3},
  {0, 3, 2}
};

/** Corners of cube diagonal for depth-0 triangle subdivision */
static const int    p4est_tnodes_cdiag[2] = {
  0, 3
};

#ifdef P4EST_ENABLE_DEBUG

/** Sequence of cube faces for depth-1 triangle subdivision */
static const int    p4est_tnodes_cface[2][2] = {
  { 2, 1 },
  { 0, 3 }
};

#endif /* P4EST_ENABLE_DEBUG */

static const int    p4est_tnodes_codim_bits[3] = {
  0, 2, 4
};

static const int    p4est_tnodes_corner_index[4] = {
  0, 2, 6, 8
};

static const int    p4est_tnodes_face_index[4] = {
  3, 5, 1, 7
};

static const int    p4est_tnodes_volume_index = 4;

#else

/** All edges of a simplex by their edge corners */
static const int    p4est_tnodes_sedge[6][2] = {
  {0, 1},
  {0, 2},
  {0, 3},
  {1, 2},
  {1, 3},
  {2, 3}
};

/** Cube corner numbers for every root simplex */
static const int    p4est_tnodes_rsim[6][4] = {
  {0, 1, 3, 7},
  {0, 5, 1, 7},
  {0, 3, 2, 7},
  {0, 2, 6, 7},
  {0, 4, 5, 7},
  {0, 6, 4, 7}
};

/** Corners of cube diagonal for depth-0 tetrahedron subdivision */
static const int    p4est_tnodes_cdiag[2] = {
  0, 7
};

#ifdef P4EST_ENABLE_DEBUG

/** Sequence of cube faces for depth-1 tetrahedron subdivision */
static const int    p4est_tnodes_cface[6][2] = {
  { 4, 1 },
  { 2, 1 },
  { 4, 3 },
  { 0, 3 },
  { 2, 5 },
  { 0, 5 }
};

/** Sequence of cube edges for depth-2 tetrahedron subdivision */
static const int    p4est_tnodes_cedge[6][2][2] = {
  {{  0,  5 },
   {  5, 11 }},
  {{  0,  9 },
   {  9,  7 }},
  {{  4,  1 },
   {  1, 11 }},
  {{  4, 10 },
   { 10,  3 }},
  {{  8,  2 },
   {  2,  7 }},
  {{  8,  6 },
   {  6,  3 }}
};

#endif /* P4EST_ENABLE_DEBUG */

static const int    p4est_tnodes_codim_bits[4] = {
  0, 3, 7, 10
};

static const int    p4est_tnodes_corner_index[8] = {
  0, 2, 6, 8, 18, 20, 24, 26
};

static const int    p4est_tnodes_edge_index[12] = {
  1,  7, 19, 25,
  3,  5, 21, 23,
  9, 11, 15, 17
};

static const int    p4est_tnodes_face_index[6] = {
  12, 14,
  10, 16,
   4, 22
};

static const int    p4est_tnodes_volume_index = 13;

#endif /* P4_TO_P8 */

#define P4EST_TNODES_SIMPLEX_ENDKEY (1 << p4est_tnodes_codim_bits[P4EST_DIM])

/** Transform element node index into integer coordinates */
static void
p4est_tnodes_ein_to_eco (p4est_tnodes_eindex_t ein,
                         p4est_tnodes_eindex_t enodes[P4EST_DIM])
{
  P4EST_ASSERT (P4EST_TNODES_IS_EIN (ein));
  P4EST_ASSERT (enodes != NULL);

#ifdef P4_TO_P8
  enodes[2] = ein / (P4EST_TNODES_FACTOR * P4EST_TNODES_FACTOR);
  enodes[1] = (ein / P4EST_TNODES_FACTOR) % P4EST_TNODES_FACTOR;
#else
  enodes[1] = ein / P4EST_TNODES_FACTOR;
#endif
  enodes[0] = ein % P4EST_TNODES_FACTOR;
}

static p4est_tnodes_eind_code_t *
p4est_tnodes_eind_code_new (void)
{
  int                 c, f, fc;
#ifdef P4_TO_P8
  int                 e, ec;
#endif
  p4est_tnodes_eindex_t ein, zwei;
  p4est_tnodes_eind_code_t *eic;

  /* initialize array to invalid entries */
  eic = P4EST_ALLOC (p4est_tnodes_eind_code_t, P4EST_TNODES_ERANGE);
  memset (eic, -1, P4EST_TNODES_ERANGE * sizeof (p4est_tnodes_eind_code_t));

  /* set corner nodes */
  for (c = 0; c < P4EST_CHILDREN; ++c) {
    ein = P4EST_TNODES_CTOEIN (c);
    P4EST_ASSERT (P4EST_TNODES_IS_EIN (ein));

    P4EST_LDEBUGF ("Corner %d ein %d eic %d\n", c, ein, eic[ein]);

    P4EST_ASSERT (eic[ein] == -1);
    eic[ein] = (P4EST_DIM << 4) + c;
  }

#ifdef P4_TO_P8
  /* set edge nodes */
  for (e = 0; e < P8EST_EDGES; ++e) {
    /* average corner indices over edge endpoints */
    ein = 0;
    for (ec = 0; ec < 2; ++ec) {
      zwei = P4EST_TNODES_CTOEIN (p8est_edge_corners[e][ec]);
      P4EST_ASSERT (P4EST_TNODES_IS_EIN (zwei));
      P4EST_ASSERT ((eic[zwei] >> 4) == P4EST_DIM);

      P4EST_LDEBUGF ("Edge %d edge corner %d zwei %d eic %d\n", e, ec, zwei,
                     eic[zwei]);

      ein += zwei;
    }
    P4EST_ASSERT ((ein & 1) == 0);
    ein >>= 1;
    P4EST_ASSERT (P4EST_TNODES_IS_EIN (ein));

    P4EST_LDEBUGF ("Edge %d point %d eic %d\n", e, ein, eic[ein]);

    P4EST_ASSERT (eic[ein] == -1);
    eic[ein] = (2 << 4) + e;
  }
#endif

  /* set face nodes */
  for (f = 0; f < P4EST_FACES; ++f) {
    /* average corner indices over face corners */
    ein = 0;
    for (fc = 0; fc < P4EST_HALF; ++fc) {
      zwei = P4EST_TNODES_CTOEIN (p4est_face_corners[f][fc]);
      P4EST_ASSERT (P4EST_TNODES_IS_EIN (zwei));
      P4EST_ASSERT ((eic[zwei] >> 4) == P4EST_DIM);

      P4EST_LDEBUGF ("Face %d face corner %d zwei %d eic %d\n", f, fc, zwei,
                     eic[zwei]);

      ein += zwei;
    }
    P4EST_ASSERT ((ein & (P4EST_HALF - 1)) == 0);
    ein >>= (P4EST_DIM - 1);
    P4EST_ASSERT (P4EST_TNODES_IS_EIN (ein));

    P4EST_LDEBUGF ("Face %d point %d eic %d\n", f, ein, eic[ein]);

    P4EST_ASSERT (eic[ein] == -1);
    eic[ein] = (1 << 4) + f;
  }

  /* set volume node */
  ein = (P4EST_TNODES_CTOEIN (p4est_tnodes_cdiag[0]) +
         P4EST_TNODES_CTOEIN (p4est_tnodes_cdiag[1])) >> 1;
  P4EST_ASSERT (P4EST_TNODES_IS_EIN (ein));

  P4EST_LDEBUGF ("Volume point %d eic %d\n", ein, eic[ein]);

  P4EST_ASSERT (eic[ein] == -1);
  eic[ein] = (0 << 4) + 0;

  /* all required entries are computed */
  return eic;
}

#ifdef P4EST_ENABLE_DEBUG

static int
p4est_tnodes_simplex_is_valid (p4est_tnodes_simplex_t *sim)
{
  int                 i;

  P4EST_ASSERT (sim != NULL);
  P4EST_ASSERT (sim->level == 0 ||
                p4est_tnodes_simplex_is_valid (sim->parent));

  /* check hierarchical structure */
  if (!(0 <= sim->index && sim->index < P4EST_TNODES_NUM_SIMPLICES) ||
      !(0 <= sim->level && sim->level <= P4EST_DIM) ||
      !(0 <= sim->key && sim->key <= P4EST_TNODES_SIMPLEX_ENDKEY)) {
    return 0;
  }
  if (!(0 <= sim->sibid && sim->sibid <=
        (sim->level == 0 ? P4EST_TNODES_NUM_SROOTS : 2))) {
    return 0;
  }

  /* check simplex corner nodes */
  for (i = 0; i < P4EST_TNODES_NUM_SCORNERS; ++i) {
    if (!P4EST_TNODES_IS_EIN (sim->nodes[i])) {
      return 0;
    }
  }
  if (sim->nodes[0] == sim->nodes[1] ||
      sim->nodes[0] == sim->nodes[2] || sim->nodes[1] == sim->nodes[2]) {
    return 0;
  }
#ifdef P4_TO_P8
  if (sim->nodes[0] == sim->nodes[3] ||
      sim->nodes[1] == sim->nodes[3] || sim->nodes[2] == sim->nodes[3]) {
    return 0;
  }
#endif
  return 1;
}

#endif /* P4EST_ENABLE_DEBUG */

/** Compute the two simplex corners of its longest edge */
static void
p4est_tnodes_longest_edge (p4est_tnodes_simplex_t * sim)
{
  int                 i, j;
  int                 mind;
  int8_t             *loedge, lo;
  int32_t             esum, msqr;
  p4est_tnodes_eindex_t enode[P4EST_TNODES_NUM_SCORNERS][P4EST_DIM];
  p4est_tnodes_eindex_t *snodes, nedge;
  p4est_tnodes_eindex_t dist;

  /* access element node coordinates */
  P4EST_ASSERT (p4est_tnodes_simplex_is_valid (sim));
  snodes = sim->nodes;
  for (i = 0; i < P4EST_TNODES_NUM_SCORNERS; ++i) {
    p4est_tnodes_ein_to_eco (snodes[i], enode[i]);
    for (j = 0; j < P4EST_DIM; ++j) {
      P4EST_ASSERT (P4EST_TNODES_IS_ECO (enode[i][j]));
    }
  }

  /* maximum of three edge lengths */
  msqr = 0;
  mind = -1;
  for (i = 0; i < P4EST_TNODES_NUM_SEDGES; ++i) {
    esum = 0;
    /* loop over edge endpoints j */
    for (j = 0; j < P4EST_DIM; ++j) {
      /* compute squared Euklidian distance */
      dist =
        enode[p4est_tnodes_sedge[i][0]][j] -
        enode[p4est_tnodes_sedge[i][1]][j];
      P4EST_ASSERT (-P4EST_TNODES_ESHIFT <= dist &&
                    dist <= P4EST_TNODES_ESHIFT);
      esum += dist * (int32_t) dist;
    }
    P4EST_ASSERT (esum > 0);
    if (esum > msqr) {
      msqr = esum;
      mind = i;
    }
  }
  P4EST_ASSERT (mind >= 0);

  /* assign simplex vertices of longest edge sorted by cube */
  loedge = sim->ledge;
  for (j = 0; j < 2; ++j) {
    loedge[j] = (int8_t) p4est_tnodes_sedge[mind][j];
  }
  if (snodes[loedge[0]] > snodes[loedge[1]]) {
    /* swap edge corners since they are backwards */
    lo = loedge[0];
    loedge[0] = loedge[1];
    loedge[1] = lo;
  }

  /* compute element index of longest edge midpoint */
  nedge = snodes[loedge[0]] + snodes[loedge[1]];
  P4EST_ASSERT (!(nedge & 1));
  nedge >>= 1;
  P4EST_ASSERT (P4EST_TNODES_IS_EIN (nedge));
  sim->lemnode = nedge;

#ifdef P4EST_ENABLE_DEBUG
  for (i = 0; i < P4EST_TNODES_NUM_SCORNERS; ++i) {
    P4EST_ASSERT (sim->lemnode != sim->nodes[i]);
  }
#endif
}

/* make room in the simplex list for one more */
static p4est_tnodes_simplex_t *
p4est_tnodes_simplex_init (sc_array_t *ttree, int *tind,
                           p4est_tnodes_simplex_t *parent)
{
  p4est_tnodes_simplex_t *sim;

  P4EST_ASSERT (ttree != NULL);
  P4EST_ASSERT (ttree->elem_size == sizeof (p4est_tnodes_simplex_t));
  P4EST_ASSERT (tind != NULL);
  P4EST_ASSERT (*tind >= 0 && *tind < P4EST_TNODES_NUM_SIMPLICES);
  P4EST_ASSERT (parent == NULL || p4est_tnodes_simplex_is_valid (parent));

  sim = (p4est_tnodes_simplex_t *) sc_array_index_int (ttree, *tind);
  sim->parent = parent;
  sim->index = (*tind)++;
  sim->level = (parent == NULL ? 0 : parent->level + 1);
  sim->key = P4EST_TNODES_SIMPLEX_ENDKEY;

  return sim;
}

/* create root simplex by precomputed table */
static void
p4est_tnodes_simplex_root (p4est_tnodes_simplex_t *sim, int d)
{
  int                 i;

  /* basic checks */
  P4EST_ASSERT (sim != NULL);
  P4EST_ASSERT (sim->parent == NULL);
  P4EST_ASSERT (sim->level == 0);
  P4EST_ASSERT (0 <= d && d < P4EST_TNODES_NUM_SROOTS);

  /* lookup parent, vertices, and longest edge */
  sim->sibid = d;
  for (i = 0; i < P4EST_TNODES_NUM_SCORNERS; ++i) {
    sim->nodes[i] = P4EST_TNODES_CTOEIN (p4est_tnodes_rsim[d][i]);
  }
  p4est_tnodes_longest_edge (sim);

  /* check and return */
  P4EST_ASSERT (p4est_tnodes_simplex_is_valid (sim));
}

/* create child simplex by making the longest edge midpoint a corner */
static p4est_tnodes_eindex_t *
p4est_tnodes_simplex_child (p4est_tnodes_simplex_t *sim, int d)
{
  int                 i;
  p4est_tnodes_simplex_t *parent;

  /* basic checks */
  P4EST_ASSERT (sim != NULL);
  P4EST_ASSERT (sim->parent != NULL);
  P4EST_ASSERT (sim->level == sim->parent->level + 1);
  P4EST_ASSERT (sim->level <= P4EST_DIM);
  P4EST_ASSERT (P4EST_TNODES_IS_EIN (sim->parent->lemnode));
  P4EST_ASSERT (0 <= d && d < 2);

  /* replace longest edge corners, ascending in d, by edge midpoint */
  sim->sibid = d;
  parent = sim->parent;
  for (i = 0; i < P4EST_TNODES_NUM_SCORNERS; ++i) {
    sim->nodes[i] =
      (parent->ledge[1 - d] == i) ? parent->lemnode : parent->nodes[i];
  }
  if (sim->level < P4EST_DIM) {
    p4est_tnodes_longest_edge (sim);
  }

  /* check and return */
  P4EST_ASSERT (p4est_tnodes_simplex_is_valid (sim));
  return sim->nodes;
}

static void
p4est_tnodes_simplex_verify (p4est_tnodes_simplex_t *sim,
                             int di[P4EST_DIM + 1])
{
#ifdef P4EST_ENABLE_DEBUG
  int                 j;

  P4EST_ASSERT (p4est_tnodes_simplex_is_valid (sim));
  for (j = 0; j <= sim->level; ++j) {
    P4EST_ASSERT (0 <= di[j] && di[j] <
                  (j == 0 ? P4EST_TNODES_NUM_SROOTS : 2));
  }
  for (; j <= P4EST_DIM; ++j) {
    P4EST_ASSERT (di[j] == -1);
  }

  switch (sim->level) {
  case 0:
    P4EST_ASSERT (sim->ledge[0] == 0);
#ifndef P4_TO_P8
    P4EST_ASSERT (sim->ledge[1] == (di[0] == 0 ? 2 : 1));
#else
    P4EST_ASSERT (sim->ledge[1] == 3);
#endif
    for (j = 0; j < 2; ++j) {
      P4EST_ASSERT (sim->nodes[sim->ledge[j]] == P4EST_TNODES_CTOEIN
                    (p4est_tnodes_cdiag[j]));
    }
    break;
  case 1:
    for (j = 0; j < 2; ++j) {
#ifndef P4_TO_P8
      P4EST_ASSERT (sim->nodes[sim->ledge[j]] == P4EST_TNODES_CTOEIN
                    (p4est_face_corners[p4est_tnodes_cface[di[0]][di[1]]]
                     [j]));
#else
      P4EST_ASSERT (sim->nodes[sim->ledge[j]] == P4EST_TNODES_CTOEIN
                    (p4est_face_corners[p4est_tnodes_cface[di[0]][di[1]]]
                     [j * 3]));
#endif
    }
    break;
#ifdef P4_TO_P8
  case 2:
    for (j = 0; j < 2; ++j) {
      P4EST_ASSERT (sim->nodes[sim->ledge[j]] == P4EST_TNODES_CTOEIN
                    (p8est_edge_corners
                     [p4est_tnodes_cedge[di[0]][di[1]][di[2]]][j]));
    }
    break;
#endif
  case P4EST_DIM:
    break;
  default:
    SC_ABORT_NOT_REACHED ();
  }
#endif

#ifndef P4_TO_P8
  P4EST_LDEBUGF ("Tri %d %d %d\n", sim->nodes[0], sim->nodes[1],
                 sim->nodes[2]);
#else
  P4EST_LDEBUGF ("Tet %d %d %d %d\n",
                 sim->nodes[0], sim->nodes[1], sim->nodes[2], sim->nodes[3]);
#endif
}

static void
p4est_tnodes_key_propagate (p4est_tnodes_simplex_t *sim)
{
  P4EST_ASSERT (sim != NULL);

  /* we have established our sort key */
  P4EST_ASSERT (0 <= sim->key && sim->key < P4EST_TNODES_SIMPLEX_ENDKEY);
  P4EST_LDEBUGF ("Depth %d simplex key 0x%04x %d\n",
                 sim->level, sim->key, sim->key);

  /* the minimum sort key over siblings occurs in the first child */
  if (sim->parent != NULL) {
    if (sim->parent->key == P4EST_TNODES_SIMPLEX_ENDKEY) {
      sim->parent->key = sim->key;
    }
    else {
      /* This is a bit of a fluke.  We actually construct in order. */
      P4EST_ASSERT (sim->key > sim->parent->key);
    }
  }
}

static sc_array_t  *
p4est_tnodes_eforest_refine (const p4est_tnodes_eind_code_t *eic)
{
  int                 i;
  int                 d0, d1, d2;
  int                 di[P4EST_DIM + 1];
#ifdef P4_TO_P8
  int                 de;
  p4est_tnodes_simplex_t *simd;
#endif
  int                 tind;
  int                 codims[P4EST_TNODES_NUM_SCORNERS], cd;
  sc_array_t         *ttree;
  p4est_tnodes_eindex_t *snodes;
  p4est_tnodes_simplex_t *sim, *sim1, *sim2;
  p4est_tnodes_eind_code_t eind_code;

  /* basic checks */
  P4EST_ASSERT (eic != NULL);
  P4EST_ASSERT (eic[0] == P4EST_DIM << 4);
  P4EST_ASSERT (eic[P4EST_TNODES_ERANGE - 1] ==
                (P4EST_DIM << 4) + (P4EST_CHILDREN - 1));

  /* create a simplex refinement forest for the reference cubic element */
  ttree = sc_array_new_count (sizeof (p4est_tnodes_simplex_t),
                              P4EST_TNODES_NUM_SIMPLICES);
  sc_array_memset (ttree, -1);
  tind = 0;
  memset (di, -1, (P4EST_DIM + 1) * sizeof (int));

  /* create the root simplices in the cube */
  for (d0 = 0; d0 < P4EST_TNODES_NUM_SROOTS; ++d0) {
    di[0] = d0;
    P4EST_LDEBUGF ("Tree %d simplex %d\n", d0, tind);

    /* construct the root simplex */
    sim = p4est_tnodes_simplex_init (ttree, &tind, NULL);
    p4est_tnodes_simplex_root (sim, d0);
    p4est_tnodes_simplex_verify (sim, di);

    /* split simplex at cube diagonal midpoint */
    P4EST_ASSERT (sim->key == P4EST_TNODES_SIMPLEX_ENDKEY);
    sim1 = sim;
    for (d1 = 0; d1 < 2; ++d1) {
      di[1] = d1;
      P4EST_LDEBUGF ("Tree %d branch %d simplex %d\n", d0, d1, tind);

      /* construct the child simplex */
      sim = p4est_tnodes_simplex_init (ttree, &tind, sim1);
      snodes = p4est_tnodes_simplex_child (sim, d1);
      p4est_tnodes_simplex_verify (sim, di);

#ifdef P4_TO_P8
      /* split simplex at face diagonal midpoint */
      P4EST_ASSERT (sim->key == P4EST_TNODES_SIMPLEX_ENDKEY);
      simd = sim;
      for (de = 0; de < 2; ++de) {
        di[2] = de;
        P4EST_LDEBUGF ("Tree %d branch %d %d simplex %d\n", d0, d1, de, tind);
#if 0
      }
#endif

      /* construct the child simplex */
      sim = p4est_tnodes_simplex_init (ttree, &tind, simd);
      snodes = p4est_tnodes_simplex_child (sim, de);
      p4est_tnodes_simplex_verify (sim, di);
#endif

      /* split simplex at face/edge diagonal midpoint */
      P4EST_ASSERT (sim->key == P4EST_TNODES_SIMPLEX_ENDKEY);
      sim2 = sim;
      for (d2 = 0; d2 < 2; ++d2) {
        di[P4EST_DIM] = d2;
#ifndef P4_TO_P8
        P4EST_LDEBUGF ("Tree %d branch %d %d simplex %d\n", d0, d1, d2, tind);
#else
        P4EST_LDEBUGF ("Tree %d branch %d %d %d simplex %d\n", d0, d1, de, d2,
                       tind);
#endif

        /* construct the child simplex */
        sim = p4est_tnodes_simplex_init (ttree, &tind, sim2);
        snodes = p4est_tnodes_simplex_child (sim, d2);
        p4est_tnodes_simplex_verify (sim, di);

        /* compute codimension of each corner to build sort key */
        P4EST_ASSERT (sim->key == P4EST_TNODES_SIMPLEX_ENDKEY);
        sim->key = 0;
        for (i = 0; i < P4EST_TNODES_NUM_SCORNERS; ++i) {
          codims[i] = -1;
        }
        for (i = 0; i < P4EST_TNODES_NUM_SCORNERS; ++i) {
          cd = (eind_code = eic[snodes[i]]) >> 4;
          P4EST_ASSERT (0 <= cd && cd <= P4EST_DIM);
          P4EST_ASSERT (cd > 0 || eind_code == 0);
          P4EST_ASSERT (codims[cd] == -1);
          codims[cd] = i;
          if (cd > 0) {
            sim->key |= (eind_code & 0x0F) << p4est_tnodes_codim_bits[cd - 1];
          }
        }

        /* propagate sort key up to parent */
        p4est_tnodes_key_propagate (sim);
      }
      di[P4EST_DIM] = -1;
      p4est_tnodes_key_propagate (sim2);
#ifdef P4_TO_P8
#if 0
      {
#endif
      }
      di[2] = -1;
      p4est_tnodes_key_propagate (simd);
#endif
    }
    di[1] = -1;
    p4est_tnodes_key_propagate (sim1);
  }
  di[0] = -1;

  /* we are done populating the elementary forest */
  P4EST_ASSERT (tind == P4EST_TNODES_NUM_SIMPLICES);
  return ttree;
}

#ifdef P4_TO_P8

static const int p4est_tnodes_third_dim[3][3] = {
  { -1,  2,  1 },
  {  2, -1,  0 },
  {  1,  0, -1 }
};

#endif /* !P4_TO_P8 */

static              void
p4est_tnodes_push_simplex (p4est_tnodes_t *tnodes,
                           int eindex[P4EST_TNODES_NUM_SCORNERS])
{
  p4est_locidx_t     *snodes;

  snodes = (p4est_locidx_t *) sc_array_push (tnodes->simplex_lnodes);
  snodes[0] = -1;
  snodes[1] = -1;
#ifdef P4_TO_P8
  snodes[2] = -1;
#endif
  snodes[P4EST_DIM] = -1;
}

static              void
p4est_tnodes_simplex_counts (p4est_tnodes_t *tnodes)
{
  int                 i;
  int                 mpiret;
  int                 mpisize, mpirank;
  p4est_locidx_t      local_tcount;
  p4est_gloidx_t      global_tcount;

  P4EST_ASSERT (tnodes != NULL);
  P4EST_ASSERT (tnodes->simplex_lnodes != NULL);
  P4EST_ASSERT (tnodes->local_tcount == NULL);

  /* allocate space for local simplex counts */
  mpiret = sc_MPI_Comm_rank (tnodes->lnodes->mpicomm, &mpirank);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_size (tnodes->lnodes->mpicomm, &mpisize);
  SC_CHECK_MPI (mpiret);
  tnodes->local_tcount = P4EST_ALLOC (p4est_locidx_t, mpisize);

  /* collect statistics on simplex counts */
  local_tcount = (p4est_locidx_t) tnodes->simplex_lnodes->elem_count;
  P4EST_ASSERT (local_tcount ==
                tnodes->local_toffset[tnodes->lnodes->num_local_elements]);
  mpiret = sc_MPI_Allgather (&local_tcount, 1, P4EST_MPI_LOCIDX,
                             tnodes->local_tcount, 1, P4EST_MPI_LOCIDX,
                             tnodes->lnodes->mpicomm);
  SC_CHECK_MPI (mpiret);
  P4EST_ASSERT (local_tcount == tnodes->local_tcount[mpirank]);

  /* count simplices globally */
  global_tcount = 0;
  for (i = 0; i < mpisize; ++i) {
    if (i == mpirank) {
      tnodes->global_toffset = global_tcount;
    }
    global_tcount += (p4est_gloidx_t) tnodes->local_tcount[i];
  }
  tnodes->global_tcount = global_tcount;
}

p4est_tnodes_t     *
p4est_tnodes_new_Q2 (p4est_lnodes_t * lnodes, int lnodes_take_ownership,
                     int construction_flags)
{
  int                 c;
  int                 f;
  int                 hi, i, k;
  int                 c_face_hanging;
#ifdef P4_TO_P8
  int                 e;
  int                 hj, j;
  int                 c_edge_hanging;
#endif
  int                 eindex[P4EST_TNODES_NUM_SCORNERS];
  p4est_locidx_t      el, ne;
  p4est_locidx_t     *enodes;
  p4est_locidx_t      enode[P4EST_TNODES_NUM_SCORNERS];
  p4est_lnodes_code_t fc;
  p4est_tnodes_t     *tnodes;

  P4EST_GLOBAL_PRODUCTIONF ("Into " P4EST_STRING "_tnodes_new_Q2 flags %x\n",
                            construction_flags);

  P4EST_ASSERT (lnodes != NULL);
  P4EST_ASSERT (lnodes->degree == 2 && lnodes->vnodes == P4EST_INSUL);
  P4EST_ASSERT (construction_flags == 0);

  /* remember lnodes in tnodes */
  tnodes = P4EST_ALLOC_ZERO (p4est_tnodes_t, 1);
  tnodes->lnodes = lnodes;
  tnodes->lnodes_owned = lnodes_take_ownership;

  /* prepare simplex arraw to grow on demand */
  tnodes->simplex_lnodes = sc_array_new
    (P4EST_TNODES_NUM_SCORNERS * sizeof (p4est_locidx_t));

  /* loop through local p4est elements */
  eindex[0] = p4est_tnodes_volume_index;
  enodes = lnodes->element_nodes;
  ne = lnodes->num_local_elements;
  tnodes->local_toffset = P4EST_ALLOC (p4est_locidx_t, ne + 1);
  tnodes->local_toffset[0] = 0;
  for (el = 0; el < ne; enodes += P4EST_INSUL, ++el) {
    fc = lnodes->face_code[el];

#if defined P4EST_ENABLE_DEBUG && defined P4_TO_P8
    /* verify that the hanging face surrounding edges are also hanging */
    if (fc) {
      for (i = 0; i < P4EST_DIM; ++i) {
        if (fc & (1 << (P4EST_DIM + i))) {
          for (j = 0; j < P4EST_DIM; ++j) {
            P4EST_ASSERT (i == j || (fc & (1 << (2 * P4EST_DIM + j))));
          }
        }
      }
    }
#endif

    /* loop through corners of element */
    for (c = 0; c < P4EST_CHILDREN; ++c) {

      /* prepare node indices */
      eindex[P4EST_DIM] = p4est_tnodes_corner_index[c];

      /* determine whether the element is hanging */
      hi = P4EST_DIM;
      c_face_hanging = 0;
#ifdef P4_TO_P8
      hj = P4EST_DIM;
      c_edge_hanging = 0;
#endif
      if (fc) {
        int                 cid, cxor;

        /* determine child id and child-relative corner id */
        cxor = (cid = fc & (P4EST_CHILDREN - 1)) ^ c;
        fc >>= P4EST_DIM;

        /* determine whether this corner is hanging */
        if (cxor != 0 && cxor != (P4EST_CHILDREN - 1)) {

          /* determine whether this corner is face hanging */
          for (hi = 0; hi < P4EST_DIM; ++hi) {
            if (cxor == ((P4EST_CHILDREN - 1) ^ (1 << hi)) && fc & (1 << hi)) {
              c_face_hanging = 1;

              /* replace corner with face node index */
              f = p4est_corner_faces[cid][hi];
              P4EST_ASSERT (f == p4est_corner_faces[c][hi]);
              eindex[P4EST_DIM] = p4est_tnodes_face_index[f];
              break;
            }
          }

#ifdef P4_TO_P8
          if (hi == P4EST_DIM) {
            /* determine whether this corner is edge hanging */
            for (hj = 0; hj < P4EST_DIM; ++hj) {
              if (cxor == (1 << hj) && (fc >> P4EST_DIM) & (1 << hj)) {
                c_edge_hanging = 1;

                /* replace corner with edge node index */
                e = p8est_corner_edges[cid][hj];
                P4EST_ASSERT (e == p8est_corner_edges[c][hj]);
                eindex[P4EST_DIM] = p4est_tnodes_edge_index[e];
                break;
              }
            }
          }
          P4EST_ASSERT (!(c_face_hanging && c_edge_hanging));
          P4EST_ASSERT ((hj != P4EST_DIM) == c_edge_hanging);
#endif
          P4EST_ASSERT ((hi != P4EST_DIM) == c_face_hanging);
        }
      }
      /* now add simplices touching this corner by edge and face */

#ifdef P4_TO_P8
      /* loop through the edges touching this corner */
      for (j = 0; j < P4EST_DIM; ++j) {
        if (c_face_hanging && j != hi) {
          /* face hanging corner: ignore all edges in the face plane */
          continue;
        }
        if (c_edge_hanging && j == hj) {
          /* edge hanging corner: ignore that same edge */
          continue;
        }
        eindex[2] = p8est_corner_edges[c][j];
#if 0
      }
#endif
#endif

      /* loop through the faces touching this corner/edge */
      for (k = 0; k < 2; ++k) {
#ifndef P4_TO_P8
        i = k;
        if (c_face_hanging && i == hi) {
          /* face hanging corner: ignore that same face */
          continue;
        }
#else
        /* compute face normal direction i */
        i = p8est_edge_faces[j << 2][k] >> 1;
        P4EST_ASSERT (i != j);
        if (c_edge_hanging) {
          int                 l = p4est_tnodes_third_dim[j][hj];

          /* the corner is a hanging edge midpoint */
          P4EST_ASSERT (j != hj);
          P4EST_ASSERT (l != -1 && l != j && l != hj);

          if (fc & (1 << l)) {
            /* the edge-hanging corner borders a hanging face j x hj */

            if (i != hj) {
              /* ignore all edges in the hanging face plane */
              P4EST_ASSERT (i == l);
              continue;
            }
            else {
              /* set edge corner to hanging face midpoint */
              P4EST_ASSERT (i != l);
              eindex[2] = p4est_corner_faces[c][l];
            }
          }
        }
#endif
        eindex[1] = p4est_corner_faces[c][i];

        /* push simplex to local list */
        p4est_tnodes_push_simplex (tnodes, eindex);

      }                         /* end face loop */
#ifdef P4_TO_P8
#if 0
      {
#endif
      }                         /* end edge loop */
#endif
    }                           /* end corner loop */

    /* update element simplex offset list */
    tnodes->local_toffset[el + 1] =
      (p4est_locidx_t) tnodes->simplex_lnodes->elem_count;

  }                             /* end element loop */
  P4EST_INFOF ("Created %ld local simplices\n",
               (long) tnodes->local_toffset[ne]);

  p4est_tnodes_simplex_counts (tnodes);
  P4EST_GLOBAL_PRODUCTIONF
    ("Done " P4EST_STRING "_tnodes_new_Q2 with %lld global simplices\n",
     (long long) tnodes->global_tcount);

  return tnodes;
}

p4est_tnodes_context_t *
p4est_tnodes_context_new (void)
{
  p4est_tnodes_context_t *econ;

  econ = P4EST_ALLOC_ZERO (p4est_tnodes_context_t, 1);
  econ->eind_code = p4est_tnodes_eind_code_new ();
  econ->eforest = p4est_tnodes_eforest_refine (econ->eind_code);

  return econ;
}

void
p4est_tnodes_context_destroy (p4est_tnodes_context_t *econ)
{
  P4EST_ASSERT (econ != NULL);

  sc_array_destroy (econ->eforest);
  P4EST_FREE (econ->eind_code);
  P4EST_FREE (econ);
}

typedef struct p4est_tnodes_private
{
  p4est_t            *p4est;        /**< For verification not use. */
}
p4est_tnodes_private_t;

/** A single contributor process to a node under construction. */
typedef struct tnodes_contr
{
  int                 nodene;       /**< Relative to element. */
  int                 rank;         /**< The referring process. */
  p4est_locidx_t      le;           /**< Element/ghost number. */
}
tnodes_contr_t;

/** A node under construction may have several contributors. */
typedef struct tnodes_cnode
{
  p4est_locidx_t      runid;        /**< Running count of node. */
  p4est_connect_type_t bcon;        /**< Codimension of node. */
  tnodes_contr_t     *owner;        /**< Owning contributor. */
  sc_array_t          contr;        /**< Contributing processes. */
}
tnodes_cnode_t;

#ifdef P4EST_ENABLE_MPI

/** Record one communication partner and/or node sharer. */
typedef struct tnodes_peer
{
  int                 rank;         /**< Rank of the peer process */
  int                 done;
  int                 sharind;      /**< Index of corresponding sharer */
  int                 passive;      /**< Number of passively shared nodes */
  p4est_locidx_t      lastadd;      /**< Most recently added node number */
  p4est_locidx_t      bufcount;     /**< Number items in message buffer */
  p4est_locidx_t      shacumul;     /**< Number owned nodes before peer */
  sc_array_t          sharedno;     /**< Remember local node with query */
  sc_array_t          querypos;     /**< Send/receive buffer for messages */
  sc_array_t          remosort;     /**< Pointer array to sort peer nodes */
}
tnodes_peer_t;

#endif /* P4EST_ENABLE_MPI */

/** Global control structure for the tnodes algorithm. */
typedef struct tnodes_meta
{
  int                 full_style;
  int                 with_faces;
#ifdef P4_TO_P8
  int                 with_edges;
#endif
  int                 mpisize, mpirank;
  int                *ghost_rank;
  int                 emptypeers;
  int                 locsharer;    /**< Index of local sharer in sharers */
  uint8_t            *chilev;
  sc_MPI_Comm         mpicomm;
  sc_array_t          construct;    /**< Collect nodes during traversal */
  sc_array_t          ownsort;      /**< Sorted owned nodes of local process */
  p4est_locidx_t      lenum;
  p4est_locidx_t      num_owned;
  p4est_locidx_t      num_owned_shared;         /**< Nodes we both own and share */
  p4est_locidx_t      num_shared;               /**< Nodes we don't own but share */
  p4est_locidx_t      num_all_shared;           /**< Nodes we share, owned or not */
  p4est_locidx_t      num_triangles;            /**< Number of local triangles */
  p4est_locidx_t      szero[P4EST_TNODES_MAXNE];
  p4est_locidx_t      smone[P4EST_TNODES_MAXNE];
  p4est_gloidx_t     *goffset;      /**< Global offsets for owned nodes */
  p4est_gloidx_t      num_global_triangles;     /**< Global number of triangles */
  p4est_t            *p4est;
  p4est_ghost_t      *ghost;
  p4est_tnodes_t     *tm;
#ifdef P4EST_ENABLE_MPI
  int                *proc_peer;
  sc_array_t          sortp;        /**< Sorted array pointing to peers */
  sc_array_t          peers;        /**< Unsorted peer storage */
  sc_array_t          pereq;        /**< Requests for storage */
#endif
}
tnodes_meta_t;

#ifndef P4_TO_P8

static int
config_cind (p4est_tnodes_config_t config)
{
  int                 cind;

  if (config <= 16) {
    cind = config;
  }
  else {
    P4EST_ASSERT (config == 32);
    cind = 17;
  }
  P4EST_ASSERT (0 <= cind && cind < 18);
  return cind;
}

#endif /* !P4_TO_P8 */

#ifdef P4EST_ENABLE_MPI

static p4est_lnodes_rank_t *
peer_sharer (tnodes_meta_t * me, int q)
{
  int                 pi;
  tnodes_peer_t      *peer;

  P4EST_ASSERT (me != NULL);
  P4EST_ASSERT (me->ghost != NULL);
  P4EST_ASSERT (me->ghost_rank != NULL);
  P4EST_ASSERT (me->proc_peer != NULL);
  P4EST_ASSERT (0 <= q && q < me->mpisize);

  /* currently we do not store a peer for the local process */
  P4EST_ASSERT (q != me->mpirank);

  pi = me->proc_peer[q];
  P4EST_ASSERT (0 < pi && pi <= me->mpisize);
  peer = (tnodes_peer_t *) sc_array_index_int (&me->peers, pi - 1);
  P4EST_ASSERT (peer->rank == q);
  return (p4est_lnodes_rank_t *) sc_array_index_int (me->tm->lnodes->sharers,
                                                     peer->sharind);
}

static tnodes_peer_t *
peer_access (tnodes_meta_t * me, int q)
{
  int                 pi;
  tnodes_peer_t      *peer;

  P4EST_ASSERT (me != NULL);
  P4EST_ASSERT (me->ghost != NULL);
  P4EST_ASSERT (me->ghost_rank != NULL);
  P4EST_ASSERT (me->proc_peer != NULL);
  P4EST_ASSERT (0 <= q && q < me->mpisize);

  /* currently we do not store a peer for the local process */
  P4EST_ASSERT (q != me->mpirank);

  if ((pi = me->proc_peer[q]) == 0) {
    peer = (tnodes_peer_t *) sc_array_push (&me->peers);
    me->proc_peer[q] = (int) me->peers.elem_count;
    peer->rank = q;
    peer->done = 0;
    peer->sharind = -1;
    peer->lastadd = -1;
    peer->bufcount = 0;
    peer->passive = 0;
    sc_array_init (&peer->sharedno, sizeof (p4est_locidx_t));
    sc_array_init (&peer->querypos, sizeof (p4est_locidx_t));
    sc_array_init (&peer->remosort, sizeof (tnodes_cnode_t *));
  }
  else {
    P4EST_ASSERT (0 < pi && pi <= me->mpisize);
    peer = (tnodes_peer_t *) sc_array_index_int (&me->peers, pi - 1);
    P4EST_ASSERT (peer->rank == q);
  }
  return peer;
}

/** The local owner process will receive a query for a node number. */
static void
peer_add_reply (tnodes_meta_t * me, tnodes_peer_t * peer, p4est_locidx_t lni)
{
  P4EST_ASSERT (me != NULL);
  P4EST_ASSERT (peer != NULL);
  P4EST_ASSERT (peer->rank > me->mpirank);
  P4EST_ASSERT (peer->lastadd < lni);
  P4EST_ASSERT (0 <= lni && lni < (p4est_locidx_t) me->construct.elem_count);

  ++peer->bufcount;
  peer->lastadd = lni;
}

/** The local process queries a remote owner for its node number. */
static void
peer_add_query (tnodes_meta_t * me, tnodes_peer_t * peer,
                p4est_locidx_t lni, p4est_locidx_t epos)
{
#ifdef P4EST_ENABLE_DEBUG
  p4est_gloidx_t      gdiff;
  p4est_lnodes_t     *ln;

  P4EST_ASSERT (me != NULL);
  P4EST_ASSERT (me->tm != NULL);
  ln = me->tm->lnodes;
  P4EST_ASSERT (ln != NULL);
  P4EST_ASSERT (peer != NULL);
  P4EST_ASSERT (peer->rank < me->mpirank);
  P4EST_ASSERT (peer->lastadd < lni);
  P4EST_ASSERT (0 <= lni && lni < (p4est_locidx_t) me->construct.elem_count);
  gdiff = (me->p4est->global_first_quadrant[peer->rank + 1] -
           me->p4est->global_first_quadrant[peer->rank]);
  P4EST_ASSERT (0 <= epos && epos < ln->vnodes * (p4est_locidx_t) gdiff);
#endif

  ++peer->bufcount;
  *(p4est_locidx_t *) sc_array_push (&peer->querypos) = epos;
  *(p4est_locidx_t *) sc_array_push (&peer->sharedno) = peer->lastadd = lni;
}

#endif /* P4EST_ENABLE_MPI */

static              p4est_locidx_t
tree_quad_to_le (p4est_t * p4est, p4est_topidx_t which_tree,
                 p4est_locidx_t quadid)
{
  p4est_tree_t       *tree;

  P4EST_ASSERT (p4est != NULL);
  P4EST_ASSERT (p4est->trees != NULL);
  tree = p4est_tree_array_index (p4est->trees, which_tree);

  P4EST_ASSERT (0 <= quadid &&
                quadid < (p4est_locidx_t) tree->quadrants.elem_count);
  return tree->quadrants_offset + quadid;
}

static void
check_node (tnodes_meta_t * me, p4est_locidx_t lni)
{
#ifdef P4EST_ENABLE_DEBUG
  tnodes_cnode_t     *cnode;
  tnodes_contr_t     *contr;
  int                 owner_rank;
  size_t              zz, siz;

  cnode = (tnodes_cnode_t *) sc_array_index (&me->construct, (size_t) lni);
  P4EST_ASSERT (cnode->runid == lni);
  owner_rank = cnode->owner->rank;
  siz = cnode->contr.elem_count;
  for (zz = 0; zz < siz; ++zz) {
    contr = (tnodes_contr_t *) sc_array_index (&cnode->contr, zz);
    P4EST_ASSERT (owner_rank <= contr->rank);
    if (owner_rank == contr->rank) {
      P4EST_ASSERT (contr == cnode->owner);
    }
  }
#endif /* P4EST_ENABLE_DEBUG */
}

/** Register a node position relative to an element.
 * The element is either process local or a ghost.
 * Multiple positions may reference the same local node.
 * We store only the smallest referrer for each process.
 */
static void
node_register (tnodes_meta_t * me, p4est_locidx_t * lni,
               int rank, p4est_locidx_t le, int nodene,
               p4est_connect_type_t bcon)
{
  tnodes_cnode_t     *cnode;
  tnodes_contr_t     *contr;
  p4est_lnodes_t     *ln;
  p4est_locidx_t      lnis;
  int                 owner_rank;
  size_t              zz, siz;

  /* basic checks */
  P4EST_ASSERT (me != NULL);
  P4EST_ASSERT (me->tm != NULL);
  ln = me->tm->lnodes;

  /* a new node is to be created or an existing one is passed in */
  if (lni == NULL) {
    lnis = -1;
    lni = &lnis;
  }
  P4EST_ASSERT (-1 <= *lni &&
                *lni < (p4est_locidx_t) me->construct.elem_count);

  /* abbreviate local rank */
  P4EST_ASSERT (rank >= -1 && rank != me->mpirank);
  if (rank == -1) {
    rank = me->mpirank;
  }
  P4EST_ASSERT (0 <= rank && rank < me->mpisize);

  /* check remaining arguments */
  P4EST_ASSERT (0 <= le && le < (p4est_locidx_t)
                (me->p4est->global_first_quadrant[rank + 1] -
                 me->p4est->global_first_quadrant[rank]));
  P4EST_ASSERT (0 <= nodene && nodene < ln->vnodes);
#ifndef P4_TO_P8
  P4EST_ASSERT (bcon == P4EST_CONNECT_FACE || bcon == P4EST_CONNECT_CORNER);
  P4EST_ASSERT (bcon != P4EST_CONNECT_FACE || (nodene >= 4));
  P4EST_ASSERT (bcon != P4EST_CONNECT_CORNER || (0 <= nodene && nodene < 9));
#else
  /* extend this as further progress is made */
  P4EST_ASSERT (bcon == P4EST_CONNECT_CORNER);
  P4EST_ASSERT (0 <= nodene && nodene < 8);
#endif

  if (*lni == -1) {
    /* create a new node with one instance */
    *lni = (p4est_locidx_t) me->construct.elem_count;
    cnode = (tnodes_cnode_t *) sc_array_push (&me->construct);
    cnode->runid = *lni;
    cnode->bcon = bcon;
    cnode->owner = NULL;
    sc_array_init (&cnode->contr, sizeof (tnodes_contr_t));
  }
  else {
    /* create a new instance of an existing node */
    cnode = (tnodes_cnode_t *) sc_array_index (&me->construct, (size_t) *lni);
    P4EST_ASSERT (cnode->runid == *lni);
    P4EST_ASSERT (cnode->bcon == bcon);
    P4EST_ASSERT (cnode->contr.elem_size == sizeof (tnodes_contr_t));
    P4EST_ASSERT (cnode->contr.elem_count > 0);
    P4EST_ASSERT (cnode->owner != NULL);
    check_node (me, *lni);
  }

  /* assign node to the local element position */
  if (rank == me->mpirank) {
    P4EST_ASSERT (ln->element_nodes[le * ln->vnodes + nodene] == -1);
    ln->element_nodes[le * ln->vnodes + nodene] = *lni;
  }

  /* iterate through instances to find matching process */
  siz = cnode->contr.elem_count;
  P4EST_ASSERT (siz == 0 || cnode->owner != NULL);
  for (zz = 0; zz < siz; ++zz) {
    contr = (tnodes_contr_t *) sc_array_index (&cnode->contr, zz);
    P4EST_ASSERT (cnode->owner->rank <= contr->rank);
    if (contr->rank == rank) {
      /* rank is found and we remember the smallest node position */
      if (le < contr->le || (le == contr->le && nodene < contr->nodene)) {
        contr->nodene = nodene;
        contr->le = le;
      }
      check_node (me, *lni);
      return;
    }
  }

  /* add new node process to the list for this node */
  owner_rank = cnode->owner == NULL ? -1 : cnode->owner->rank;
  contr = (tnodes_contr_t *) sc_array_push (&cnode->contr);
  contr->nodene = nodene;
  contr->rank = rank;
  contr->le = le;
  if (cnode->owner == NULL || rank <= owner_rank) {
    /* the array was empty before or we know to set a new owner */
    P4EST_ASSERT (rank != owner_rank);
    cnode->owner = contr;
  }
  else {
    /* pushing to the array has invalidated the previous owner pointer */
    cnode->owner = (tnodes_contr_t *) sc_array_index (&cnode->contr, 0);
    siz = cnode->contr.elem_count;
    for (zz = 1; zz < siz - 1; ++zz) {
      contr = (tnodes_contr_t *) sc_array_index (&cnode->contr, zz);
      if (contr->rank < cnode->owner->rank) {
        cnode->owner = contr;
      }
    }
  }
  check_node (me, *lni);
}

static void
node_lregister (tnodes_meta_t * me, p4est_locidx_t * lni,
                p4est_locidx_t le, int nodene, p4est_connect_type_t bcon)
{
  node_register (me, lni, -1, le, nodene, bcon);
}

#ifndef P4_TO_P8

static void
node_lfacetocorner (tnodes_meta_t * me, p4est_locidx_t le, int nodene)
{
  tnodes_cnode_t     *cnode;
  p4est_lnodes_t     *ln;
  p4est_locidx_t      lni;

  P4EST_ASSERT (me != NULL);
  P4EST_ASSERT (me->tm != NULL);
  ln = me->tm->lnodes;
  P4EST_ASSERT (0 <= le && le < ln->num_local_elements);

  P4EST_ASSERT (nodene == n_center);

  /* access node in local element */
  lni = ln->element_nodes[le * ln->vnodes + nodene];
  P4EST_ASSERT (lni >= 0);
  cnode = (tnodes_cnode_t *) sc_array_index (&me->construct, lni);
  P4EST_ASSERT (cnode->owner != NULL);
  P4EST_ASSERT (cnode->runid == lni);
  P4EST_ASSERT (cnode->bcon == P4EST_CONNECT_FACE);
  P4EST_ASSERT (cnode->contr.elem_size == sizeof (tnodes_contr_t));
  P4EST_ASSERT (cnode->contr.elem_count == 1);

  /* change instance of an existing node */
  cnode->bcon = P4EST_CONNECT_CORNER;
}

#endif /* !P4_TO_P8 */

static void
node_gregister (tnodes_meta_t * me, p4est_locidx_t * lni,
                p4est_locidx_t ghostid, int nodene, p4est_connect_type_t bcon)
{
  p4est_quadrant_t   *gquad;

  P4EST_ASSERT (me != NULL);
  P4EST_ASSERT (me->tm != NULL);
  P4EST_ASSERT (0 <= nodene && nodene < me->tm->lnodes->vnodes);
#ifndef P4_TO_P8
  P4EST_ASSERT (!alwaysowned[nodene]);
#else
  /* extend this as further progress is made */
  P4EST_ASSERT (0 <= nodene && nodene < 8);
#endif

  if (me->ghost != NULL) {
    P4EST_ASSERT (me->ghost_rank != NULL);
    P4EST_ASSERT (0 <= ghostid &&
                  ghostid < (p4est_locidx_t) me->ghost->ghosts.elem_count);

    /* extract remote element number from ghost quadrant */
    gquad = (p4est_quadrant_t *) sc_array_index (&me->ghost->ghosts, ghostid);
    node_register (me, lni, me->ghost_rank[ghostid],
                   gquad->p.piggy3.local_num, nodene, bcon);
  }
}

static void
iter_volume1 (p4est_iter_volume_info_t * vi, void *user_data)
{
  tnodes_meta_t      *me = (tnodes_meta_t *) user_data;
  p4est_tnodes_t     *tm = me->tm;
  p4est_locidx_t      le;
#ifndef P4_TO_P8
  int                 j;
#endif
  int                 childid;
  int8_t              level;
#ifdef P4EST_ENABLE_DEBUG
  p4est_lnodes_t     *ln = tm->lnodes;
  p4est_tree_t       *tree;

  /* initial checks  */
  P4EST_ASSERT (vi->p4est == me->p4est);
  tree = p4est_tree_array_index (vi->p4est->trees, vi->treeid);
  P4EST_ASSERT (tree->quadrants_offset + vi->quadid == me->lenum);
#endif

  /* store quadrant level and child id */
  le = me->lenum++;
  level = vi->quad->level;
  childid = p4est_quadrant_child_id (vi->quad);
  me->chilev[le] = (((uint8_t) level) << 3) | ((uint8_t) childid);
  P4EST_ASSERT (tm->configuration[le] == 0);
  P4EST_ASSERT (ln->face_code[le] == 0);
  P4EST_ASSERT (!memcmp (ln->element_nodes + le * ln->vnodes,
                         me->smone, sizeof (p4est_locidx_t) * ln->vnodes));

#ifndef P4_TO_P8
  /* add nodes as required */
  if (me->full_style || level == 0) {
    tm->configuration[le] = (((p4est_tnodes_config_t) 1) << 5);
    node_lregister (me, NULL, le, n_center, P4EST_CONNECT_CORNER);
    if (me->with_faces) {
      for (j = 0; j < 4; ++j) {
        node_lregister (me, NULL, le, n_cface[j], P4EST_CONNECT_FACE);
      }
    }
  }
  else {
    if (childid == 1 || childid == 2) {
      tm->configuration[le] = (((p4est_tnodes_config_t) 1) << 4);
    }
    if (me->with_faces) {
      node_lregister (me, NULL, le, n_center, P4EST_CONNECT_FACE);
    }
  }
#endif /* !P4_TO_P8 */
}

static void
iter_face1 (p4est_iter_face_info_t * fi, void *user_data)
{
#ifndef P4_TO_P8
  tnodes_meta_t      *me = (tnodes_meta_t *) user_data;
  int                 i, j;
  int                 face;
  int                 nodene;
  int                 childid;
  int                 swapi;
  p4est_locidx_t      le;               /**< local element number */
  p4est_locidx_t      lni;              /**< local node number */
  p4est_locidx_t      lnh[2];           /**< hanging face mids */
  p4est_iter_face_side_t *fs;
  p4est_iter_face_side_t *fss[2];
  p4est_iter_face_side_full_t *fu;
  p4est_iter_face_side_hanging_t *fh;
  p4est_lnodes_t     *ln = me->tm->lnodes;

  /* initial checks  */
  P4EST_ASSERT (fi->p4est == me->p4est);

  /* a boundary face is the easiest case */
  fs = fss[0] = fss[1] = NULL;
  if (fi->sides.elem_count == 1) {
    P4EST_ASSERT (fi->orientation == 0);
    P4EST_ASSERT (fi->tree_boundary == P4EST_CONNECT_FACE);
    fs = (p4est_iter_face_side_t *) sc_array_index_int (&fi->sides, 0);
    P4EST_ASSERT (!fs->is_hanging);
    fu = &fs->is.full;
    P4EST_ASSERT (!fu->is_ghost);
    /* a boundary face does not contribute to the configuration */
    if (me->with_faces) {
      le = tree_quad_to_le (fi->p4est, fs->treeid, fu->quadid);
      node_lregister (me, NULL, le, n_mface[fs->face], P4EST_CONNECT_FACE);
    }
    return;
  }

  /* we have two sides to the face connection */
  P4EST_ASSERT (fi->sides.elem_count == 2);
  fss[0] = (p4est_iter_face_side_t *) sc_array_index_int (&fi->sides, 0);
  fss[1] = (p4est_iter_face_side_t *) sc_array_index_int (&fi->sides, 1);
  P4EST_ASSERT (!fss[0]->is_hanging || !fss[1]->is_hanging);
  if (!fss[0]->is_hanging && !fss[1]->is_hanging) {
    /* same size face connection does not contribute to configuration */
    if (me->with_faces) {
      lni = -1;
      for (i = 0; i < 2; ++i) {
        fu = &fss[i]->is.full;
        nodene = n_mface[fss[i]->face];
        if (!fu->is_ghost) {
          le = tree_quad_to_le (fi->p4est, fss[i]->treeid, fu->quadid);
          node_lregister (me, &lni, le, nodene, P4EST_CONNECT_FACE);
        }
        else {
          node_gregister (me, &lni, fu->quadid, nodene, P4EST_CONNECT_FACE);
        }
      }
    }
    return;
  }

  /* one of the two sides is hanging */
  lni = lnh[0] = lnh[1] = -1;
  for (i = 0; i < 2; ++i) {
    swapi = (i == 0 || fi->orientation == 0) ? 0 : 1;
    if (!fss[i]->is_hanging) {
      fu = &fss[i]->is.full;
      face = fss[i]->face;
      nodene = n_mface[face];

      if (!fu->is_ghost) {
        /* this is a large local quadrant which must insert the face midpoint */
        le = tree_quad_to_le (fi->p4est, fss[i]->treeid, fu->quadid);
        if ((me->tm->configuration[le] & ~(1 << 4)) == 0) {
          /* this is a half refinement, which must be promoted to full */
          P4EST_ASSERT (!me->full_style && fu->quad->level > 0);
          if (!me->with_faces) {
            node_lregister (me, NULL, le, n_center, P4EST_CONNECT_CORNER);
          }
          else {
            node_lfacetocorner (me, le, n_center);
            for (j = 0; j < 4; ++j) {
              node_lregister (me, NULL, le, n_cface[j], P4EST_CONNECT_FACE);
            }
          }
        }
        me->tm->configuration[le] &= ~((1 << 4) | (1 << 5));
        me->tm->configuration[le] |= (1 << face);
        node_lregister (me, &lni, le, nodene, P4EST_CONNECT_CORNER);
        if (me->with_faces) {
          node_lregister (me, NULL, le, n_split[face], P4EST_CONNECT_FACE);
          for (j = 0; j < 2; ++j) {
            node_lregister (me, &lnh[j ^ swapi], le, n_hface[face][j],
                            P4EST_CONNECT_FACE);
          }
        }
      }
      else {
        node_gregister (me, &lni, fu->quadid, nodene, P4EST_CONNECT_CORNER);
        if (me->with_faces) {
          for (j = 0; j < 2; ++j) {
            node_gregister (me, &lnh[j ^ swapi], fu->quadid, n_hface[face][j],
                            P4EST_CONNECT_FACE);
          }
        }
      }
    }
    else {
      P4EST_ASSERT (fss[i]->is_hanging);

      /* for each small local quadrant contribute to its face code */
      fh = &fss[i]->is.hanging;
      face = fss[i]->face;
      for (j = 0; j < P4EST_HALF; ++j) {
        nodene = n_cornr[p4est_face_corners[face][j ^ (P4EST_HALF - 1)]];
        if (!fh->is_ghost[j]) {
          le = tree_quad_to_le (fi->p4est, fss[i]->treeid, fh->quadid[j]);
          node_lregister (me, &lni, le, nodene, P4EST_CONNECT_CORNER);
          if (me->with_faces) {
            node_lregister (me, &lnh[j ^ swapi], le, n_mface[face],
                            P4EST_CONNECT_FACE);
          }

          /* update the face code */
          childid = p4est_quadrant_child_id (fh->quad[j]);
          P4EST_ASSERT (p4est_face_corners[face][j] == childid);
          P4EST_ASSERT (face == p4est_corner_faces[childid][face >> 1]);
          P4EST_ASSERT ((ln->face_code[le] &
                         (1 << (P4EST_DIM + (face >> 1)))) == 0);
          ln->face_code[le] |= (1 << (P4EST_DIM + (face >> 1))) | childid;
        }
        else {
          node_gregister (me, &lni, fh->quadid[j], nodene,
                          P4EST_CONNECT_CORNER);
          if (me->with_faces) {
            node_gregister (me, &lnh[j ^ swapi], fh->quadid[j], n_mface[face],
                            P4EST_CONNECT_FACE);
          }
        }
      }
    }
  }
#endif /* !P4_TO_P8 */
}

#ifdef P4_TO_P8

static void
iter_edge1 (p8est_iter_edge_info_t * ei, void *user_data)
{
}

#endif /* P4_TO_P8 */

static void
iter_corner1 (p4est_iter_corner_info_t * ci, void *user_data)
{
  tnodes_meta_t      *me = (tnodes_meta_t *) user_data;
  p4est_iter_corner_side_t *cs;

  p4est_locidx_t      le;               /**< local element number */
  p4est_locidx_t      lni;              /**< local node number */
  int                 nodene;
  size_t              zz;

  /* initial checks  */
  P4EST_ASSERT (ci->p4est == me->p4est);

  lni = -1;
  for (zz = 0; zz < ci->sides.elem_count; ++zz) {
    cs = (p4est_iter_corner_side_t *) sc_array_index (&ci->sides, zz);
    P4EST_ASSERT (0 <= cs->corner && cs->corner < P4EST_CHILDREN);
    nodene = n_cornr[cs->corner];
    if (!cs->is_ghost) {
      le = tree_quad_to_le (ci->p4est, cs->treeid, cs->quadid);
      node_lregister (me, &lni, le, nodene, P4EST_CONNECT_CORNER);
    }
    else {
      node_gregister (me, &lni, cs->quadid, nodene, P4EST_CONNECT_CORNER);
    }
  }
}

static int
cnode_compare (const void *v1, const void *v2)
{
  const tnodes_cnode_t **cc1 = (const tnodes_cnode_t **) v1;
  const tnodes_cnode_t **cc2 = (const tnodes_cnode_t **) v2;
  const tnodes_contr_t *o1, *o2;
  p4est_locidx_t      ldiff;

  P4EST_ASSERT (cc1 != NULL && *cc1 != NULL);
  P4EST_ASSERT (cc2 != NULL && *cc2 != NULL);

  P4EST_ASSERT ((*cc1)->contr.elem_size == sizeof (tnodes_contr_t));
  P4EST_ASSERT ((*cc2)->contr.elem_size == sizeof (tnodes_contr_t));

  /* we sort within the same owner process */
  o1 = (*cc1)->owner;
  o2 = (*cc2)->owner;
  P4EST_ASSERT (o1->rank == o2->rank);

  /* if the elements are the same, compare node position */
  if (o1->le == o2->le) {
    return o1->nodene - o2->nodene;
  }

  /* nodes are sorted according to their element number */
  ldiff = o1->le - o2->le;
  return ldiff == 0 ? 0 : ldiff < 0 ? -1 : 1;
}

static void
owned_query_reply (tnodes_meta_t * me)
{
  tnodes_cnode_t     *cnode, **ccn;
  tnodes_contr_t     *owner;
#ifdef P4EST_ENABLE_MPI
  tnodes_contr_t     *contr;
  tnodes_peer_t      *peer;
  p4est_lnodes_t     *ln = me->tm->lnodes;
  int                 withloc;
  size_t              zc, sic;
#endif
  size_t              zz, siz;

  /* lookup nodes separately per process */
  P4EST_ASSERT (me->num_owned == 0);
  P4EST_ASSERT (me->num_owned_shared == 0);
  P4EST_ASSERT (me->num_shared == 0);
  P4EST_ASSERT (me->num_all_shared == 0);
  siz = me->construct.elem_count;
  for (zz = 0; zz < siz; ++zz) {
    cnode = (tnodes_cnode_t *) sc_array_index (&me->construct, zz);
    P4EST_ASSERT (cnode->runid == (p4est_locidx_t) zz);
    check_node (me, cnode->runid);

    owner = cnode->owner;
    if (owner->rank == me->mpirank) {
      ccn = (tnodes_cnode_t **) sc_array_push (&me->ownsort);
      *ccn = cnode;
      ++me->num_owned;

#ifdef P4EST_ENABLE_MPI
      /* post replies for all queries to self */
      sic = cnode->contr.elem_count;
      for (zc = 0; zc < sic; ++zc) {
        contr = (tnodes_contr_t *) sc_array_index (&cnode->contr, zc);
        if (contr->rank != me->mpirank) {
          P4EST_ASSERT (contr->rank > me->mpirank);
          peer = peer_access (me, contr->rank);
          peer_add_reply (me, peer, cnode->runid);
        }
        else {
          P4EST_ASSERT (owner == contr);
        }
      }
      if (sic > 1) {
        ++me->num_owned_shared;
      }
#endif
    }
    else {
#ifdef P4EST_ENABLE_MPI
      /* weed out remote-only nodes */
      withloc = 0;
      sic = cnode->contr.elem_count;
      for (zc = 0; zc < sic; ++zc) {
        contr = (tnodes_contr_t *) sc_array_index (&cnode->contr, zc);
        if (contr->rank == me->mpirank) {
          withloc = 1;
          break;
        }
      }
      if (!withloc) {
        cnode->runid = -1;
        continue;
      }
      P4EST_ASSERT (owner->rank < me->mpirank);

      /* check for passively shared nodes */
      for (zc = 0; zc < sic; ++zc) {
        contr = (tnodes_contr_t *) sc_array_index (&cnode->contr, zc);
        if (contr->rank != me->mpirank && contr->rank != owner->rank) {
          /* passively share a remotely owned node */
          P4EST_ASSERT (contr->rank > owner->rank);
          peer = peer_access (me, contr->rank);
          ++peer->passive;
        }
      }

      /* post query to remote owner */
      peer = peer_access (me, owner->rank);
      peer_add_query (me, peer, cnode->runid,
                      owner->le * ln->vnodes + owner->nodene);
      ccn = (tnodes_cnode_t **) sc_array_push (&peer->remosort);
      *ccn = cnode;
      ++me->num_shared;
#else
      SC_ABORT_NOT_REACHED ();
#endif
    }

    /* the running id will be replaced by the owner's node number */
    cnode->runid = -1;
  }
  me->num_all_shared = me->num_owned_shared + me->num_shared;
}

static void
sort_allgather (tnodes_meta_t * me)
{
  tnodes_cnode_t    **ccn;
  p4est_tnodes_t     *tm = me->tm;
  p4est_lnodes_t     *ln = tm->lnodes;
  p4est_locidx_t      lel, le, lc;
  p4est_locidx_t     *localboth, lb[2];
  p4est_gloidx_t      gc;
  const int           s = me->mpisize;
  int                 mpiret;
#ifndef P4_TO_P8
  int                 cind, lookup;
#endif
  int                 q;
  size_t              zz;

  /* sort local node list */
  sc_array_sort (&me->ownsort, cnode_compare);
  for (zz = 0; zz < me->ownsort.elem_count; ++zz) {
    ccn = (tnodes_cnode_t **) sc_array_index (&me->ownsort, zz);
    (*ccn)->runid = zz;
  }

  /* share owned count */
  ln->num_local_nodes = (ln->owned_count = me->num_owned) + me->num_shared;
  ln->nonlocal_nodes = P4EST_ALLOC (p4est_gloidx_t, me->num_shared);
  ln->global_owned_count = P4EST_ALLOC (p4est_locidx_t, s);
  lb[0] = ln->owned_count;

  /* establish local triangle count */
  lel = ln->num_local_elements;
  lc = tm->local_toffset[0] = 0;
  for (le = 0; le < lel; ++le) {
#ifndef P4_TO_P8
    cind = config_cind (me->tm->configuration[le]);
    lookup = p4est_tnodes_config_lookup[cind];
    P4EST_ASSERT (0 <= lookup && lookup < 6);
    lc = tm->local_toffset[le + 1] =
      lc + p4est_tnodes_lookup_counts[lookup][2];
#else
    /* extend this as further progress is made */
    lc = tm->local_toffset[le + 1] = lc + 0;
#endif
  }
  lb[1] = me->num_triangles = lc;

  /* parallel sharing of owned node and element counts */
  localboth = P4EST_ALLOC (p4est_locidx_t, 2 * me->mpisize);
  mpiret = sc_MPI_Allgather (lb, 2, P4EST_MPI_LOCIDX,
                             localboth, 2, P4EST_MPI_LOCIDX,
                             me->p4est->mpicomm);
  SC_CHECK_MPI (mpiret);
  me->goffset = P4EST_ALLOC (p4est_gloidx_t, s + 1);
  gc = me->goffset[0] = 0;
  P4EST_ASSERT (me->num_global_triangles == 0);
  for (q = 0; q < s; ++q) {
    gc = me->goffset[q + 1] =
      gc + (ln->global_owned_count[q] = localboth[2 * q + 0]);
    if (q == me->mpirank) {
      tm->global_toffset = me->num_global_triangles;
    }
    me->num_global_triangles += (tm->local_tcount[q] = localboth[2 * q + 1]);
  }
  ln->global_offset = me->goffset[me->mpirank];
  P4EST_FREE (localboth);
}

#ifdef P4EST_ENABLE_MPI

static int
peer_compare (const void *v1, const void *v2)
{
  const tnodes_peer_t **p1 = (const tnodes_peer_t **) v1;
  const tnodes_peer_t **p2 = (const tnodes_peer_t **) v2;
  return (*p1)->rank - (*p2)->rank;
}

static int
rnode_compare (const void *v1, const void *v2)
{
  const tnodes_cnode_t **cc1 = (const tnodes_cnode_t **) v1;
  const tnodes_cnode_t **cc2 = (const tnodes_cnode_t **) v2;
#ifdef P4EST_ENABLE_DEBUG
  const tnodes_contr_t *o1, *o2;
#endif
  p4est_locidx_t      ldiff;

  P4EST_ASSERT (cc1 != NULL && *cc1 != NULL);
  P4EST_ASSERT (cc2 != NULL && *cc2 != NULL);

  P4EST_ASSERT ((*cc1)->contr.elem_size == sizeof (tnodes_contr_t));
  P4EST_ASSERT ((*cc2)->contr.elem_size == sizeof (tnodes_contr_t));

#ifdef P4EST_ENABLE_DEBUG
  /* we sort within the same owner process */
  o1 = (*cc1)->owner;
  o2 = (*cc2)->owner;
  P4EST_ASSERT (o1->rank == o2->rank);
#endif

  /* nodes are sorted according to their runid member */
  P4EST_ASSERT ((*cc1)->runid >= 0);
  P4EST_ASSERT ((*cc2)->runid >= 0);
  ldiff = (*cc1)->runid - (*cc2)->runid;
  return ldiff == 0 ? 0 : ldiff < 0 ? -1 : 1;
}

static void
push_sharer (tnodes_meta_t * me, int *sindex, int rank)
{
  p4est_lnodes_rank_t *sharer;
  p4est_lnodes_t     *ln = me->tm->lnodes;

  P4EST_ASSERT (me != NULL);
  P4EST_ASSERT (sindex != NULL);
  P4EST_ASSERT (0 <= rank && rank < me->mpisize);

  /* push empty sharer structure */
  *sindex = (int) ln->sharers->elem_count;
  sharer = (p4est_lnodes_rank_t *) sc_array_push (ln->sharers);
  memset (sharer, -1, sizeof (p4est_lnodes_rank_t));
  sharer->rank = rank;
  sc_array_init (&sharer->shared_nodes, sizeof (p4est_locidx_t));
}

#endif /* P4EST_ENABLE_MPI */

static void
sort_peers (tnodes_meta_t * me)
{
#ifndef P4EST_ENABLE_MPI
  P4EST_ASSERT (me != NULL);
  P4EST_ASSERT (me->num_all_shared == 0);
#else
  int                 i;
  int                 num_peers;
  p4est_locidx_t      nonlofs;
  tnodes_peer_t      *tp;

  /* explicitly do nothing without a ghost layer */
  num_peers = (int) me->peers.elem_count;
  if (me->ghost == NULL || num_peers == 0) {
    P4EST_ASSERT (me->num_all_shared == 0);
    return;
  }
  P4EST_ASSERT (me->num_all_shared > 0);

  /* make it possible to iterate through peers in rank order */
  sc_array_resize (&me->sortp, num_peers);
  for (i = 0; i < num_peers; ++i) {
    *(tnodes_peer_t **) sc_array_index_int (&me->sortp, i) =
      (tnodes_peer_t *) sc_array_index_int (&me->peers, i);
  }
  sc_array_sort (&me->sortp, peer_compare);
  nonlofs = 0;
  for (i = 0; i < num_peers; ++i) {
    tp = *(tnodes_peer_t **) sc_array_index_int (&me->sortp, i);
    tp->shacumul = nonlofs;
    if (tp->rank < me->mpirank) {
      nonlofs += tp->bufcount;
    }
  }
  P4EST_ASSERT (nonlofs == me->num_shared);

  /* initialize sharers array */
  for (i = 0; i < num_peers; ++i) {
    tp = *(tnodes_peer_t **) sc_array_index_int (&me->sortp, i);
    P4EST_ASSERT (tp->rank != me->mpirank);
    if (tp->rank > me->mpirank) {
      break;
    }
    push_sharer (me, &tp->sharind, tp->rank);
  }
  push_sharer (me, &me->locsharer, me->mpirank);
  for (; i < num_peers; ++i) {
    tp = *(tnodes_peer_t **) sc_array_index_int (&me->sortp, i);
    P4EST_ASSERT (tp->rank > me->mpirank);
    push_sharer (me, &tp->sharind, tp->rank);
  }
  P4EST_ASSERT (num_peers + 1 == (int) me->tm->lnodes->sharers->elem_count);
  P4EST_ASSERT (me->locsharer >= 0);
#endif /* P4EST_ENABLE_MPI */
}

static void
post_query_reply (tnodes_meta_t * me)
{
#ifndef P4EST_ENABLE_MPI
  P4EST_ASSERT (me != NULL);
  P4EST_ASSERT (me->num_all_shared == 0);
#else
  int                 mpiret;
  size_t              zp, iz;
  sc_MPI_Request     *preq;
  tnodes_peer_t      *peer;

  /* explicitly do nothing without a ghost layer */
  zp = me->peers.elem_count;
  if (me->ghost == NULL || zp == 0) {
    P4EST_ASSERT (me->num_all_shared == 0);
    return;
  }
  P4EST_ASSERT (me->num_all_shared >= 0);

  /* go through peers (unsorted) and post messages */
  P4EST_ASSERT (me->emptypeers == 0);
  sc_array_resize (&me->pereq, zp);
  for (iz = 0; iz < zp; ++iz) {
    peer = (tnodes_peer_t *) sc_array_index (&me->peers, iz);
    preq = (sc_MPI_Request *) sc_array_index (&me->pereq, iz);
    if (peer->bufcount == 0) {
      /* purely passive peers do not send messages */
      P4EST_ASSERT (peer->passive > 0);
      *preq = sc_MPI_REQUEST_NULL;
      ++me->emptypeers;
      continue;
    }
    if (peer->rank > me->mpirank) {
      /* expect query from higher rank */
      P4EST_ASSERT (peer->querypos.elem_count == 0);
      sc_array_resize (&peer->querypos, peer->bufcount);
      mpiret = sc_MPI_Irecv (sc_array_index (&peer->querypos, 0),
                             peer->bufcount, P4EST_MPI_LOCIDX, peer->rank,
                             P4EST_COMM_TNODES_QUERY, me->mpicomm, preq);
      SC_CHECK_MPI (mpiret);
      peer->done = 1;
    }
    else {
      /* address query to lower rank */
      P4EST_ASSERT (peer->rank < me->mpirank);
      P4EST_ASSERT (peer->bufcount ==
                    (p4est_locidx_t) peer->querypos.elem_count);
      mpiret =
        sc_MPI_Isend (sc_array_index (&peer->querypos, 0), peer->bufcount,
                      P4EST_MPI_LOCIDX, peer->rank, P4EST_COMM_TNODES_QUERY,
                      me->mpicomm, preq);
      SC_CHECK_MPI (mpiret);
      peer->done = 3;
    }
  }
#endif /* P4EST_ENABLE_MPI */
}

static void
wait_query_reply (tnodes_meta_t * me)
{
#ifndef P4EST_ENABLE_MPI
  P4EST_ASSERT (me != NULL);
  P4EST_ASSERT (me->num_all_shared == 0);
#else
  int                 i, j;
  int                 mpiret;
  int                 nwalloc;
  int                 nwtotal;
  int                 nwaited;
  int                *waitind;
  sc_MPI_Request     *preq;
  p4est_locidx_t      lbc, lcl, lni, nonloc;
  p4est_locidx_t      epos, oind;
  p4est_gloidx_t      gof, gni;
  p4est_lnodes_t     *ln = me->tm->lnodes;
  tnodes_cnode_t     *cnode;
  tnodes_peer_t      *peer;

  /* explicitly do nothing without a ghost layer */
  nwalloc = (int) me->peers.elem_count;
  if (me->ghost == NULL || nwalloc == 0) {
    P4EST_ASSERT (me->num_all_shared == 0);
    return;
  }
  P4EST_ASSERT (me->num_all_shared >= 0);

  /* currently the local process does not count as peer */
  nwtotal = nwalloc - me->emptypeers;
  P4EST_ASSERT (nwtotal > 0);
  waitind = P4EST_ALLOC (int, nwalloc);
  while (nwtotal > 0) {
    mpiret = sc_MPI_Waitsome
      (nwalloc, (sc_MPI_Request *) sc_array_index (&me->pereq, 0),
       &nwaited, waitind, sc_MPI_STATUSES_IGNORE);
    SC_CHECK_MPI (mpiret);
    SC_CHECK_ABORT (nwaited > 0, "Invalid count after MPI_Waitsome");
    for (i = 0; i < nwaited; ++i) {
      j = waitind[i];
      peer = (tnodes_peer_t *) sc_array_index (&me->peers, j);
      P4EST_ASSERT (peer->rank != me->mpirank);
      preq = (sc_MPI_Request *) sc_array_index (&me->pereq, j);
      P4EST_ASSERT (*preq == sc_MPI_REQUEST_NULL);
      if (peer->rank > me->mpirank) {
        P4EST_ASSERT (peer->shacumul == me->num_shared);
        P4EST_ASSERT (peer->sharedno.elem_count == 0);
        if (peer->done == 1) {
          /* we have received a request and shall send a reply */
          lbc = peer->bufcount;
          for (lcl = 0; lcl < lbc; ++lcl) {
            epos = *(p4est_locidx_t *) sc_array_index (&peer->querypos, lcl);
#if 0
            P4EST_LDEBUGF ("Got %d gquad %d pos %d\n from %d\n", lcl,
                           epos / ln->vnodes, epos % ln->vnodes, peer->rank);
#endif
            P4EST_ASSERT (0 <= epos && epos < ln->vnodes * ln->owned_count);
#ifndef P4_TO_P8
            P4EST_ASSERT (!alwaysowned[epos % ln->vnodes]);
#else
            /* extend this as further progress is made */
            P4EST_ASSERT (epos % ln->vnodes < 8);
#endif
            lni = ln->element_nodes[epos];
            P4EST_ASSERT (0 <= lni && lni <
                          (p4est_locidx_t) me->construct.elem_count);
            cnode = (tnodes_cnode_t *) sc_array_index (&me->construct, lni);
            oind = cnode->runid;
            P4EST_ASSERT (0 <= oind && oind < ln->owned_count);

            /* send back the number of node owned by the local process */
            *(p4est_locidx_t *) sc_array_index (&peer->querypos, lcl) = oind;
          }
          mpiret = sc_MPI_Isend (sc_array_index (&peer->querypos, 0),
                                 peer->bufcount, P4EST_MPI_LOCIDX, peer->rank,
                                 P4EST_COMM_TNODES_REPLY, me->mpicomm, preq);
          SC_CHECK_MPI (mpiret);
          peer->done = 2;
        }
        else {
          /* our reply has been received */
          P4EST_ASSERT (peer->done == 2);
          peer->done = 0;
          --nwtotal;
        }
      }
      else {
        P4EST_ASSERT (peer->rank < me->mpirank);
        if (peer->done == 3) {
          /* our request has been sent and we await the reply */
          mpiret = sc_MPI_Irecv (sc_array_index (&peer->querypos, 0),
                                 peer->bufcount, P4EST_MPI_LOCIDX, peer->rank,
                                 P4EST_COMM_TNODES_REPLY, me->mpicomm, preq);
          SC_CHECK_MPI (mpiret);
          peer->done = 4;
        }
        else {
          /* process owner's node numbers in reply received */
          P4EST_ASSERT (peer->done == 4);
          lbc = peer->bufcount;
          for (lcl = 0; lcl < lbc; ++lcl) {
            oind = *(p4est_locidx_t *) sc_array_index (&peer->querypos, lcl);
            lni = *(p4est_locidx_t *) sc_array_index (&peer->sharedno, lcl);
            cnode = (tnodes_cnode_t *) sc_array_index (&me->construct, lni);
            P4EST_ASSERT (cnode->owner->rank == peer->rank);
            P4EST_ASSERT (0 <= oind &&
                          oind < ln->global_owned_count[peer->rank]);
            cnode->runid = oind;
          }
          sc_array_sort (&peer->remosort, rnode_compare);

          /* store shared node's global index */
          gof = me->goffset[peer->rank];
          for (lcl = 0; lcl < lbc; ++lcl) {
            cnode =
              *(tnodes_cnode_t **) sc_array_index (&peer->remosort, lcl);
            P4EST_ASSERT (cnode->owner->rank == peer->rank);
            P4EST_ASSERT (lcl <= cnode->runid);
            nonloc = peer->shacumul + lcl;
            P4EST_ASSERT (nonloc < me->num_shared);
            gni = gof + cnode->runid;
            P4EST_ASSERT (me->goffset[peer->rank] <= gni &&
                          gni < me->goffset[peer->rank + 1]);
            ln->nonlocal_nodes[nonloc] = gni;

            /* now the runid of each node is the local node number */
            cnode->runid = me->num_owned + nonloc;
          }
          peer->done = 0;
          --nwtotal;
        }
      }
    }
  }
  P4EST_FREE (waitind);
#ifdef P4EST_ENABLE_DEBUG
  gof = -1;
  for (lcl = 0; lcl < me->num_shared; ++lcl) {
    gni = ln->nonlocal_nodes[lcl];
    P4EST_ASSERT (0 <= gni && gni < me->goffset[me->mpisize]);
    P4EST_ASSERT (gni < me->goffset[me->mpirank] ||
                  gni >= me->goffset[me->mpirank + 1]);
    P4EST_ASSERT (gni > gof);
    gof = gni;
  }
#endif /* P4EST_ENABLE_DEBUG */
#endif /* P4EST_ENABLE_MPI */
}

#ifndef P4_TO_P8

static void
set_element_node (tnodes_meta_t * me, p4est_locidx_t le, int nodene)
{
  p4est_lnodes_t     *ln = me->tm->lnodes;
  p4est_locidx_t      lni, runid;
  tnodes_cnode_t     *cnode;

  P4EST_ASSERT (0 <= le && le < ln->num_local_elements);
  P4EST_ASSERT (0 <= nodene && nodene < ln->vnodes);
  lni = ln->element_nodes[le * ln->vnodes + nodene];
#if 0
  P4EST_LDEBUGF ("lni for %ld, %d: %ld\n", (long) le, nodene, (long) lni);
#endif
  P4EST_ASSERT (0 <= lni && lni < (p4est_locidx_t) me->construct.elem_count);

  cnode = (tnodes_cnode_t *) sc_array_index (&me->construct, lni);
  runid = cnode->runid;
#if 0
  P4EST_LDEBUGF ("Run id %ld, %d: %ld\n", nodene, (long) lni, (long) runid);
#endif
  P4EST_ASSERT (0 <= runid && runid < ln->num_local_nodes);
  P4EST_ASSERT ((runid < me->num_owned && cnode->owner->rank == me->mpirank)
                || (runid >= me->num_owned
                    && cnode->owner->rank < me->mpirank));
#ifdef P4EST_ENABLE_DEBUG
  if (runid >= me->num_owned) {
    lni = (p4est_locidx_t) (ln->nonlocal_nodes[runid - me->num_owned] -
                            me->goffset[cnode->owner->rank]);
    P4EST_ASSERT (0 <= lni &&
                  lni < ln->global_owned_count[cnode->owner->rank]);
  }
#endif
  ln->element_nodes[le * ln->vnodes + nodene] = runid;
}

static void
assign_element_nodes (tnodes_meta_t * me)
{
  int                 nodene, lookup;
  int                 ncorner, nface, cind;
  int                 ci, fi;
#ifdef P4EST_ENABLE_DEBUG
  int                 poswhich[P4EST_TNODES_MAXNE];
#endif
  p4est_lnodes_t     *ln = me->tm->lnodes;
  p4est_locidx_t      le, lel;

  /* assign final numbers of element nodes */
  lel = ln->num_local_elements;
  for (le = 0; le < lel; ++le) {
    cind = config_cind (me->tm->configuration[le]);
#ifdef P4EST_ENABLE_DEBUG
    memset (poswhich, -1, P4EST_TNODES_MAXNE * sizeof (int));
#endif
    lookup = p4est_tnodes_config_lookup[cind];
    P4EST_ASSERT (0 <= lookup && lookup < 6);
    ncorner = p4est_tnodes_lookup_counts[lookup][0];
    P4EST_ASSERT (4 <= ncorner && ncorner <= 9);
    for (ci = 0; ci < ncorner; ++ci) {
      nodene = p4est_tnodes_config_corners[cind][ci];
      P4EST_ASSERT (0 <= nodene && nodene <= 8);
      P4EST_ASSERT (poswhich[nodene] == -1);
#ifdef P4EST_ENABLE_DEBUG
      poswhich[nodene] = P4EST_DIM;
#endif
      set_element_node (me, le, nodene);
    }
#ifdef P4EST_ENABLE_DEBUG
    for (; ci < 9; ++ci) {
      P4EST_ASSERT (p4est_tnodes_config_corners[cind][ci] == -1);
    }
#endif
    if (me->with_faces) {
      nface = p4est_tnodes_lookup_counts[lookup][1];
      P4EST_ASSERT (5 <= nface && nface <= 16);
      for (fi = 0; fi < nface; ++fi) {
        nodene = p4est_tnodes_config_faces[cind][fi];
        P4EST_ASSERT (4 <= nodene && nodene <= 24);
        P4EST_ASSERT (poswhich[nodene] == -1);
#ifdef P4EST_ENABLE_DEBUG
        poswhich[nodene] = 1;
#endif
        set_element_node (me, le, nodene);
      }
#ifdef P4EST_ENABLE_DEBUG
      for (; fi < 16; ++fi) {
        P4EST_ASSERT (p4est_tnodes_config_faces[cind][fi] == -1);
      }
#endif
    }
#ifdef P4EST_ENABLE_DEBUG
    for (nodene = 0; nodene < ln->vnodes; ++nodene) {
      if (poswhich[nodene] == -1) {
        P4EST_ASSERT (ln->element_nodes[le * ln->vnodes + nodene] == -1);
      }
    }
#endif
  }
}

#endif /* !P4_TO_P8 */

static void
populate_sharers (tnodes_meta_t * me)
{
#ifdef P4EST_ENABLE_MPI
  int                 i;
  int                 num_peers;
  size_t              zz;
  p4est_locidx_t      lni;
  p4est_locidx_t      lbc, lcl;
  p4est_lnodes_t     *ln = me->tm->lnodes;
  p4est_lnodes_rank_t *sharer, *locshare;
  tnodes_peer_t      *tp;
  tnodes_cnode_t     *cnode;
  tnodes_contr_t     *contr;

  /* populate sharers array */
  num_peers = (int) me->peers.elem_count;
  if (me->ghost == NULL || num_peers == 0) {
    P4EST_ASSERT (me->num_all_shared == 0);
    return;
  }
  P4EST_ASSERT (me->num_all_shared >= 0);
  P4EST_ASSERT (num_peers + 1 == (p4est_locidx_t) ln->sharers->elem_count);
  locshare =
    (p4est_lnodes_rank_t *) sc_array_index_int (ln->sharers, me->locsharer);

  /* first iterate through owned nodes in order */
  lbc = (p4est_locidx_t) me->ownsort.elem_count;
  for (lcl = 0; lcl < lbc; ++lcl) {
    cnode = *(tnodes_cnode_t **) sc_array_index (&me->ownsort, lcl);
    P4EST_ASSERT (cnode->owner->rank == me->mpirank);
    P4EST_ASSERT (lcl == cnode->runid);
    if (cnode->contr.elem_count == 1) {
      /* this node is purely local */
      continue;
    }

    /* this node has sharers: iterate through all of them */
    for (zz = 0; zz < cnode->contr.elem_count; ++zz) {
      contr = (tnodes_contr_t *) sc_array_index (&cnode->contr, zz);
      if (contr->rank == me->mpirank) {
        /* local process is owner */
        P4EST_ASSERT (cnode->owner == contr);
        sharer = locshare;
      }
      else {
        /* remote process is sharer */
        sharer = peer_sharer (me, contr->rank);
      }
      P4EST_ASSERT (sharer->rank == contr->rank);
      *(p4est_locidx_t *) sc_array_push (&sharer->shared_nodes) = lcl;
    }
  }
  P4EST_ASSERT (me->num_owned == (p4est_locidx_t) me->ownsort.elem_count);
  P4EST_ASSERT (me->num_owned_shared ==
                (p4est_locidx_t) locshare->shared_nodes.elem_count);

  /* determine the sharer offset and count variables */
  locshare->shared_mine_offset = locshare->owned_offset = 0;
  locshare->shared_mine_count = me->num_owned_shared;
  locshare->owned_count = me->num_owned;
  for (i = 0; i < num_peers; ++i) {
    tp = *(tnodes_peer_t **) sc_array_index_int (&me->sortp, i);
    sharer =
      (p4est_lnodes_rank_t *) sc_array_index_int (ln->sharers, tp->sharind);
    P4EST_ASSERT (tp->rank == sharer->rank);
    sharer->shared_mine_offset = 0;
    sharer->shared_mine_count =
      (p4est_locidx_t) sharer->shared_nodes.elem_count;
    sharer->owned_offset = me->num_owned + tp->shacumul;
    if (tp->rank < me->mpirank) {
      P4EST_ASSERT (tp->bufcount > 0 || tp->passive > 0);
      sharer->owned_count = tp->bufcount;
    }
    else {
      P4EST_ASSERT (tp->rank > me->mpirank);
      sharer->owned_count = 0;
    }
  }

  /* iterate through the remote local nodes in order */
  lni = me->num_owned;
  for (i = 0; i < num_peers; ++i) {
    tp = *(tnodes_peer_t **) sc_array_index_int (&me->sortp, i);
    if (tp->rank < me->mpirank) {
      lbc = tp->bufcount;
      P4EST_ASSERT (lbc == (p4est_locidx_t) tp->remosort.elem_count);
      for (lcl = 0; lcl < lbc; ++lcl, ++lni) {
        cnode = *(tnodes_cnode_t **) sc_array_index (&tp->remosort, lcl);
        P4EST_ASSERT (cnode->owner->rank == tp->rank);
        P4EST_ASSERT (cnode->runid == lni);

        /* this node has sharers: iterate through all of them */
        P4EST_ASSERT (cnode->contr.elem_count > 1);
        for (zz = 0; zz < cnode->contr.elem_count; ++zz) {
          contr = (tnodes_contr_t *) sc_array_index (&cnode->contr, zz);
          if (contr->rank == me->mpirank) {
            /* local process is sharer */
            P4EST_ASSERT (cnode->owner != contr);
            sharer = locshare;
          }
          else {
            /* remote process is owner or not */
            sharer = peer_sharer (me, contr->rank);
          }
        }
        P4EST_ASSERT (sharer->rank == contr->rank);
        *(p4est_locidx_t *) sc_array_push (&sharer->shared_nodes) = lni;
      }
    }
  }
  P4EST_ASSERT (lni == ln->num_local_nodes);
#endif /* P4EST_ENABLE_MPI */
}

static void
clean_construct (tnodes_meta_t * me)
{
  tnodes_cnode_t     *cnode;
  size_t              zz;

  for (zz = 0; zz < me->construct.elem_count; ++zz) {
    cnode = (tnodes_cnode_t *) sc_array_index (&me->construct, zz);
    sc_array_reset (&cnode->contr);
  }
}

p4est_tnodes_t     *
p4est_tnodes_new (p4est_t * p4est, p4est_ghost_t * ghost, int full_style,
                  int with_faces
#ifdef P4_TO_P8
                  , int with_edges
#endif
  )
{
  int                 q, s;
  int                 vn;
  p4est_locidx_t      lel, lg, ng;
  p4est_tnodes_t     *tm;
  p4est_lnodes_t     *ln;
  tnodes_meta_t       tmeta, *me = &tmeta;
#ifdef P4EST_ENABLE_DEBUG
  p4est_tnodes_config_t configure;
#if 0
  p4est_lnodes_t     *testnodes;
#endif
  p4est_locidx_t      le;
#endif
#ifdef P4EST_ENABLE_MPI
  size_t              nz, zi;
  tnodes_peer_t      *peer;
#endif

  P4EST_ASSERT (p4est_is_balanced (p4est, P4EST_CONNECT_ALMOST));

  /* basic assignment of members */
  memset (me, 0, sizeof (tnodes_meta_t));
  memset (me->smone, -1, P4EST_TNODES_MAXNE * sizeof (p4est_locidx_t));
  me->p4est = p4est;
  me->mpicomm = p4est->mpicomm;
  s = me->mpisize = p4est->mpisize;
  me->mpirank = p4est->mpirank;
  tm = me->tm = P4EST_ALLOC_ZERO (p4est_tnodes_t, 1);
  tm->full_style = me->full_style = full_style;
  tm->with_faces = me->with_faces = with_faces;
#ifdef P4_TO_P8
  tm->with_edges = me->with_edges = with_edges;
#endif
  ln = tm->lnodes = P4EST_ALLOC_ZERO (p4est_lnodes_t, 1);
  tm->lnodes_owned = 1;
  me->locsharer = -1;
  tm->pri = P4EST_ALLOC_ZERO (p4est_tnodes_private_t, 1);
  tm->pri->p4est = p4est;

  /* lookup structure for ghost owner rank */
  if ((me->ghost = ghost) != NULL) {
    P4EST_ASSERT (ghost->proc_offsets[0] == 0);
    P4EST_ASSERT (ghost->proc_offsets[s] ==
                  (p4est_locidx_t) ghost->ghosts.elem_count);
    me->ghost_rank = P4EST_ALLOC (int, ghost->ghosts.elem_count);
    lg = 0;
    for (q = 0; q < s; ++q) {
      ng = ghost->proc_offsets[q + 1];
      for (; lg < ng; ++lg) {
        me->ghost_rank[lg] = q;
      }
    }
    P4EST_ASSERT (lg == (p4est_locidx_t) ghost->ghosts.elem_count);
#ifdef P4EST_ENABLE_MPI
    me->proc_peer = P4EST_ALLOC_ZERO (int, s);
    sc_array_init (&me->sortp, sizeof (tnodes_peer_t *));
    sc_array_init (&me->peers, sizeof (tnodes_peer_t));
    sc_array_init (&me->pereq, sizeof (sc_MPI_Request));
#endif
  }
  sc_array_init (&me->construct, sizeof (tnodes_cnode_t));
  sc_array_init (&me->ownsort, sizeof (tnodes_cnode_t *));

  /* prepare node information */
  ln->mpicomm = p4est->mpicomm;
  ln->sharers = sc_array_new (sizeof (p4est_lnodes_rank_t));
  ln->degree = 0;
#ifndef P4_TO_P8
  vn = ln->vnodes = 9 + (with_faces ? 16 : 0);
#else
  /* extend this as further progress is made */
  vn = ln->vnodes = 8;
#endif
  lel = ln->num_local_elements = p4est->local_num_quadrants;
  P4EST_ASSERT ((size_t) lel * (size_t) vn <= (size_t) P4EST_LOCIDX_MAX);
  me->chilev = P4EST_ALLOC_ZERO (uint8_t, lel);
  ln->face_code = P4EST_ALLOC_ZERO (p4est_lnodes_code_t, lel);
  ln->element_nodes = P4EST_ALLOC (p4est_locidx_t, lel * vn);
  memset (ln->element_nodes, -1, lel * vn * sizeof (p4est_locidx_t));

  /* allocate arrays for node encoding */
  tm->configuration = P4EST_ALLOC_ZERO (p4est_tnodes_config_t, lel);
  tm->local_toffset = P4EST_ALLOC (p4est_locidx_t, lel + 1);
  tm->local_tcount = P4EST_ALLOC (p4est_locidx_t, s);

  /* determine triangle configuration of each element */
  me->lenum = 0;
  p4est_iterate (p4est, ghost, me, iter_volume1, iter_face1,
#ifdef P4_TO_P8
                 iter_edge1,
#endif
                 iter_corner1);
  P4EST_ASSERT (me->lenum == lel);

  owned_query_reply (me);
  P4EST_INFOF ("p4est_tnodes_new: nodes owned %ld shared %ld\n",
               (long) me->num_owned, (long) me->num_shared);

  /* post first round of messages */
  post_query_reply (me);

  /* sort local nodes and allgather owned counts */
  sort_allgather (me);
  P4EST_INFOF ("p4est_tnodes_new: triangles owned %ld\n",
               (long) me->num_triangles);
  P4EST_GLOBAL_PRODUCTIONF
    ("p4est_tnodes_new: global triangles %lld nodes %lld\n",
     (long long) me->num_global_triangles, (long long) me->goffset[s]);

  /* sort peers by process */
  sort_peers (me);

  /* receive query messages and send replies */
  wait_query_reply (me);

  /* finalize element node assignment */
#ifndef P4_TO_P8
  assign_element_nodes (me);
#endif /* P4_TO_P8 */

  /* finalize sharer information */
  populate_sharers (me);

  /* free memory */
  P4EST_FREE (me->goffset);
  if (me->ghost != NULL) {
#ifdef P4EST_ENABLE_MPI
    nz = me->peers.elem_count;
    for (zi = 0; zi < nz; ++zi) {
      peer = (tnodes_peer_t *) sc_array_index (&me->peers, zi);
      P4EST_ASSERT (!peer->done);
      sc_array_reset (&peer->sharedno);
      sc_array_reset (&peer->querypos);
      sc_array_reset (&peer->remosort);
    }
    P4EST_FREE (me->proc_peer);
    sc_array_reset (&me->sortp);
    sc_array_reset (&me->peers);
    sc_array_reset (&me->pereq);
#endif
    P4EST_FREE (me->ghost_rank);
  }
  clean_construct (me);
  sc_array_reset (&me->construct);
  sc_array_reset (&me->ownsort);
  P4EST_FREE (me->chilev);

#ifdef P4EST_ENABLE_DEBUG
  for (le = 0; le < lel; ++le) {
    configure = tm->configuration[le];
#ifndef P4_TO_P8
    P4EST_ASSERT (configure <= 16 || configure == 32);
#else
    P4EST_ASSERT (configure == 0);
#endif
  }
  if (me->ghost != NULL) {
#if 0
    testnodes = p4est_lnodes_new (p4est, me->ghost, 2);
    P4EST_ASSERT (testnodes->num_local_elements == lel);
    for (le = 0; le < lel; ++le) {
      P4EST_ASSERT (testnodes->face_code[le] == ln->face_code[le]);
    }
    p4est_lnodes_destroy (testnodes);
#endif
  }
#endif

  return tm;
}

void
p4est_tnodes_destroy (p4est_tnodes_t * tm)
{
  P4EST_ASSERT (tm != NULL);
  P4EST_ASSERT (tm->lnodes != NULL);

  if (tm->lnodes_owned) {
    p4est_lnodes_destroy (tm->lnodes);
  }
  if (tm->simplex_lnodes != NULL) {
    sc_array_destroy (tm->simplex_lnodes);
  }
  P4EST_FREE (tm->configuration);
  P4EST_FREE (tm->local_toffset);
  P4EST_FREE (tm->local_tcount);
  P4EST_FREE (tm->pri);
  P4EST_FREE (tm);
}

#ifndef P4_TO_P8

static void
iter_triangle_properties (p4est_tnodes_iter_t * it)
{
  int                 i;
  int                 tindex;
  int                 conode, fanode;
  const int          *tnodin;
  p4est_tnodes_t     *tnodes;
  p4est_tnodes_iter_private_t *pri;
#ifdef P4EST_ENABLE_DEBUG
  int                 lookup;
  p4est_locidx_t      lni;
  p4est_lnodes_t     *ln;

  P4EST_ASSERT (it != NULL);
#endif
  tnodes = it->tnodes;
  P4EST_ASSERT (tnodes != NULL);
  pri = it->pri;
  P4EST_ASSERT (pri != NULL);

#ifdef P4EST_ENABLE_DEBUG
  /* verify iterator state */
  P4EST_ASSERT (p4est_quadrant_is_valid (it->quadrant));
  P4EST_ASSERT (it->which_quad ==
                pri->tree->quadrants_offset + pri->treequad);
  P4EST_ASSERT (0 <= pri->quadtri && pri->quadtri < pri->numqtri);
  P4EST_ASSERT (0 <= it->triangle && it->triangle < pri->numtris);
  P4EST_ASSERT (0 <= pri->cind && pri->cind < 18);
  lookup = p4est_tnodes_config_lookup[pri->cind];
  P4EST_ASSERT (p4est_tnodes_lookup_counts[lookup][2] == pri->numqtri);
  ln = tnodes->lnodes;
  P4EST_ASSERT (ln != NULL);
#endif

  /* lookup triangle */
  tindex = p4est_tnodes_config_triangles[pri->cind][pri->quadtri];
  P4EST_ASSERT (0 <= tindex && tindex < 16);
  tnodin = p4est_tnodes_triangle_nodes[tindex];
#ifdef P4EST_ENABLE_DEBUG
  for (i = pri->numqtri; i < 8; ++i) {
    P4EST_ASSERT (p4est_tnodes_config_triangles[pri->cind][i] == -1);
  }
#endif

  /* loop through triangle nodes */
  for (i = 0; i < 3; ++i) {
    /* corner node */
    conode = tnodin[i];
#ifdef P4EST_ENABLE_DEBUG
    P4EST_ASSERT (0 <= conode && conode < 9);
    lni = ln->element_nodes[it->which_quad * ln->vnodes + conode];
    P4EST_ASSERT (0 <= lni && lni < ln->num_local_nodes);
#endif
    it->corner_nodes[i] = conode;

    /* face node */
    if (tnodes->with_faces) {
      fanode = tnodin[i + 3];
#ifdef P4EST_ENABLE_DEBUG
      P4EST_ASSERT (4 <= fanode && fanode < 25);
      lni = ln->element_nodes[it->which_quad * ln->vnodes + fanode];
      P4EST_ASSERT (0 <= lni && lni < ln->num_local_nodes);
#endif
      it->face_nodes[i] = fanode;
    }
  }
}

p4est_tnodes_iter_t *
p4est_tnodes_iter_new (p4est_t * p4est, p4est_tnodes_t * tnodes)
{
  p4est_lnodes_t     *ln;
  p4est_tnodes_iter_t *it;
  p4est_tnodes_iter_private_t *pri;
  int                 lookup;

  P4EST_ASSERT (p4est != NULL);
  P4EST_ASSERT (tnodes != NULL);
  P4EST_ASSERT (tnodes->pri != NULL);
  P4EST_ASSERT (tnodes->pri->p4est == p4est);

  ln = tnodes->lnodes;
  P4EST_ASSERT (ln != NULL);
  P4EST_ASSERT (ln->degree == 0);
  P4EST_ASSERT (ln->vnodes == (tnodes->with_faces ? 25 : 9));
  P4EST_ASSERT (ln->num_local_elements == p4est->local_num_quadrants);

  /* return for an empty process */
  if (ln->num_local_elements == 0) {
    return NULL;
  }
  P4EST_ASSERT (p4est->first_local_tree <= p4est->last_local_tree);

  /* create iterator context */
  it = P4EST_ALLOC (p4est_tnodes_iter_t, 1);
  memset (it, -1, sizeof (p4est_tnodes_iter_t));
  it->p4est = p4est;
  it->tnodes = tnodes;
  pri = it->pri = P4EST_ALLOC (p4est_tnodes_iter_private_t, 1);

  /* populate iterator state */
  pri->numtris = it->tnodes->local_tcount[p4est->mpirank];
  pri->tree = p4est_tree_array_index (p4est->trees, it->which_tree =
                                      p4est->first_local_tree);
  pri->numtreeq = (p4est_locidx_t) pri->tree->quadrants.elem_count;
  it->quadrant = p4est_quadrant_array_index (&pri->tree->quadrants,
                                             pri->treequad = 0);
  pri->cind = config_cind (tnodes->configuration[it->which_quad = 0]);
  lookup = p4est_tnodes_config_lookup[pri->cind];
  pri->numqtri = p4est_tnodes_lookup_counts[lookup][2];
  pri->quadtri = 0;
  it->triangle = 0;

  /* access current triangle properties */
  iter_triangle_properties (it);

  /* iterator now points to first local triangle */
  return it;
}

void
p4est_tnodes_iter_next (p4est_tnodes_iter_t ** pit)
{
  p4est_tnodes_iter_t *it;
  p4est_tnodes_iter_private_t *pri;
  p4est_t            *p4est;
  int                 lookup;

  /* access iterator state */
  P4EST_ASSERT (pit != NULL);
  it = *pit;
  P4EST_ASSERT (it != NULL);
  p4est = it->p4est;
  P4EST_ASSERT (p4est != NULL);
  P4EST_ASSERT (it->tnodes != NULL);
  pri = it->pri;
  P4EST_ASSERT (pri != NULL);

  /* access and store next triangle */
  if ((++it->triangle, ++pri->quadtri) == pri->numqtri) {
    /* we leave the current quadrant */
    if ((++it->which_quad, ++pri->treequad) == pri->numtreeq) {
      /* we leave the current tree */
      if (it->which_tree++ == p4est->last_local_tree) {
        /* we are done iterating und free the iterator state */
        P4EST_ASSERT (it->which_quad == p4est->local_num_quadrants);
        P4EST_ASSERT (pri->quadtri == pri->numqtri);
        P4EST_ASSERT (it->triangle == pri->numtris);
        P4EST_FREE (it->pri);
        P4EST_FREE (it);
        *pit = NULL;
        return;
      }
      pri->tree = p4est_tree_array_index (p4est->trees, it->which_tree);
      pri->numtreeq = (p4est_locidx_t) pri->tree->quadrants.elem_count;
      pri->treequad = 0;
    }
    it->quadrant = p4est_quadrant_array_index (&pri->tree->quadrants,
                                               pri->treequad);
    P4EST_ASSERT (0 <= it->which_quad &&
                  it->which_quad < p4est->local_num_quadrants);
    pri->cind = config_cind (it->tnodes->configuration[it->which_quad]);
    lookup = p4est_tnodes_config_lookup[pri->cind];
    pri->numqtri = p4est_tnodes_lookup_counts[lookup][2];
    pri->quadtri = 0;
  }

  /* access current triangle properties */
  iter_triangle_properties (it);
}

#endif /* !P4_TO_P8 */
