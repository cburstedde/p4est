/*
 * Usage: p4est_delaunay <configuration> <level> <vtk basename>
 *        possible configurations:
 *        o unit      Refinement on the unit square.
 *        o three     Refinement on a forest with three trees.
 *        o moebius   Refinement on a 5-tree Moebius band.
 *        o star      Refinement on a 6-tree star shaped domain.
 *        o periodic  Refinement on the unit square with all-periodic b.c.
 *        o rotwrap   Refinement on the unit square with weird periodic b.c.
 *        o disk      Refinement on a 5-tree flat disk or square.
 *        o pdisk     Refinement on 5-tree flat disk or square, periodic b.c.
 */

#ifndef P4_TO_P8
#include <p4est_lnodes.h>
#include <p4est_extended.h>
#include <p4est_connectivity.h>
#include <p4est_bits.h>
#include <p4est_geometry.h>
#include <p4est_vtk.h>
#else
#include <p8est_lnodes.h>
#include <p8est_extended.h>
#include <p8est_connectivity.h>
#include <p8est_bits.h>
#include <p8est_geometry.h>
#include <p8est_vtk.h>
#endif /* P4_TO_P8 */
#include <sc_options.h>

typedef struct _opts
{
  int minlevel;
  int maxlevel;
  const char *vtk;
  const char *conn;
}
opts_t;

static int refine_fn (p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t *quadrant)
{
  opts_t *opts = (opts_t *) p4est->user_pointer;
  int refine_level = opts->maxlevel;
  if ((int) quadrant->level >= (refine_level - (int) (which_tree % 3))) {
    return 0;
  }
  if (quadrant->level == 1 && p4est_quadrant_child_id (quadrant) == 3) {
    return 1;
  }
  if (quadrant->x == P4EST_LAST_OFFSET (2) &&
      quadrant->y == P4EST_LAST_OFFSET (2)) {
    return 1;
  }
#ifndef P4_TO_P8
  if (quadrant->x >= P4EST_QUADRANT_LEN (2)) {
    return 0;
  }
#else
  if (quadrant->z >= P8EST_QUADRANT_LEN (2)) {
    return 0;
  }
#endif

  return 1;
}

static p4est_t *create_p4est_from_opts(opts_t *opts)
{
  p4est_connectivity_t *conn = p4est_connectivity_new_byname (opts->conn);
  p4est_t * p4est = p4est_new_ext (sc_MPI_COMM_WORLD, conn, 0, opts->minlevel, 1, 0, NULL, NULL);
  p4est->user_pointer = (void *) opts;
  p4est_refine_ext (p4est, 1, opts->maxlevel, refine_fn, NULL, NULL);
  p4est_partition (p4est, 0, NULL);
  p4est_balance (p4est, P4EST_CONNECT_FULL, NULL);
  return p4est;
}

#ifndef P4_TO_P8
#define P4EST_SIMPLICES 2
#else
#define P4EST_SIMPLICES 6
#endif

typedef struct _p4est_simplex_nodes
{
  sc_array_t *vertices;
  sc_array_t *simplices;
  p4est_lnodes_t *lnodes;
} p4est_simplex_nodes_t;

static void
p4est_quadrant_get_ordered_simplices (p4est_locidx_t corner_nodes[], p4est_gloidx_t corner_nodes_global[], p4est_locidx_t (*simplices)[P4EST_DIM + 1])
{
  p4est_gloidx_t min_global = corner_nodes_global[0];
  int min_loc = 0;
  int o; /* positive oriented symmetry? */
  for (int i = 1; i < P4EST_CHILDREN; i++) {
    if (corner_nodes_global[i] < min_global) {
      min_global = corner_nodes_global[i];
      min_loc = i;
    }
  }

  o = (min_loc == 0 || min_loc == 3 || min_loc == 5 || min_loc == 6);
#ifndef P4_TO_P8
  simplices[0][0] = corner_nodes[min_loc]; simplices[0][1] = corner_nodes[min_loc ^ (o ? 1 : 3)]; simplices[0][2] = corner_nodes[min_loc ^ (o ? 3 : 1)];
  simplices[1][0] = corner_nodes[min_loc]; simplices[1][1] = corner_nodes[min_loc ^ (o ? 3 : 2)]; simplices[1][2] = corner_nodes[min_loc ^ (o ? 2 : 3)];
#else
  for (int d = 0; d < P4EST_DIM; d++) {
    int positive = o ^ (d == 1);
    int d2 = SC_MIN (((d + 1) % 3), ((d + 2) % 3));
    int d3 = SC_MAX (((d + 1) % 3), ((d + 2) % 3));
    int a = min_loc ^ (1 << d);
    int b = a ^ (1 << d2);
    int c = a ^ (1 << d3);
    int e = min_loc ^ 7;
    p4est_gloidx_t min_ae = SC_MIN (corner_nodes_global[a], corner_nodes_global[e]);
    p4est_gloidx_t min_bc = SC_MIN (corner_nodes_global[b], corner_nodes_global[c]);

    if (min_ae < min_bc) {
      simplices[d * 2 + 0][0] = corner_nodes[min_loc]; simplices[d * 2 + 0][1] = corner_nodes[a]; simplices[d * 2 + 0][2] = corner_nodes[positive ? b : e]; simplices[d * 2 + 0][3] = corner_nodes[positive ? e : b];
      simplices[d * 2 + 1][0] = corner_nodes[min_loc]; simplices[d * 2 + 1][1] = corner_nodes[a]; simplices[d * 2 + 1][2] = corner_nodes[positive ? e : c]; simplices[d * 2 + 1][3] = corner_nodes[positive ? c : e];
    } else {
      simplices[d * 2 + 0][0] = corner_nodes[min_loc]; simplices[d * 2 + 0][1] = corner_nodes[a]; simplices[d * 2 + 0][2] = corner_nodes[positive ? b : c]; simplices[d * 2 + 0][3] = corner_nodes[positive ? c : b];
      simplices[d * 2 + 1][0] = corner_nodes[min_loc]; simplices[d * 2 + 1][1] = corner_nodes[e]; simplices[d * 2 + 1][2] = corner_nodes[positive ? c : b]; simplices[d * 2 + 1][3] = corner_nodes[positive ? b : c];
    }
  }
#endif
}

/* create simplicial nodes from fully 2:1 balanced forest */
p4est_simplex_nodes_t *
p4est_simplex_nodes_create (p4est_t *p4est, p4est_ghost_t *ghost_layer, p4est_geometry_t *geometry)
{
  p4est_ghost_t *ghost = NULL;
  p4est_geometry_t *geom = NULL;
  p4est_lnodes_t   *lnodes;
  p4est_simplex_nodes_t *snodes = P4EST_ALLOC_ZERO(p4est_simplex_nodes_t, 1);
  sc_array_t *simplices;
  sc_array_t *vertices;

  ghost = ghost_layer ? ghost_layer : p4est_ghost_new (p4est, P4EST_CONNECT_FULL);
  geom = geometry ? geometry : p4est_geometry_new_connectivity (p4est->connectivity);
  snodes->lnodes = lnodes = p4est_lnodes_new (p4est, ghost, 1);
  snodes->simplices = simplices = sc_array_new ((P4EST_DIM + 1) * sizeof(p4est_locidx_t));
  snodes->vertices = vertices = sc_array_new_size (3 * sizeof(double), lnodes->num_local_nodes);

  {
    p4est_locidx_t elem = 0;
    for (p4est_topidx_t t = p4est->first_local_tree; t <= p4est->last_local_tree; t++) {
      p4est_tree_t *tree = p4est_tree_array_index(p4est->trees, t);
      sc_array_t *quads = &tree->quadrants;
      for (p4est_locidx_t q = 0; q < quads->elem_count; q++, elem++) {
        p4est_quadrant_t *quad = p4est_quadrant_array_index (quads, q);
        p4est_locidx_t E[P4EST_CHILDREN];

        for (int i = 0; i < P4EST_CHILDREN; i++) {
          E[i] = lnodes->element_nodes[elem * P4EST_CHILDREN + i];
        }

        for (int i = 0; i < P4EST_CHILDREN; i++) {
          p4est_quadrant_t node;
          double XYZ[3] = {0};

          p4est_quadrant_corner_node (quad, i, &node);
          XYZ[0] = ((double) node.x) / ((double) P4EST_ROOT_LEN);
          XYZ[1] = ((double) node.y) / ((double) P4EST_ROOT_LEN);
#ifdef P4_TO_P8
          XYZ[2] = ((double) node.z) / ((double) P4EST_ROOT_LEN);
#endif
          geom->X (geom, t, XYZ, sc_array_index (vertices, E[i]));
        }

        if (quad->level == 0) {
          /* unstructured reference-space delaunay of connectivity mesh */
          p4est_gloidx_t E_global[P4EST_CHILDREN];

          for (int i = 0; i < P4EST_CHILDREN; i++) {
            E_global[i] = p4est_lnodes_global_index (lnodes, E[i]);
          }
          p4est_quadrant_get_ordered_simplices (E, E_global, sc_array_push_count (simplices, P4EST_SIMPLICES));
        } else {
          int sims[P4EST_SIMPLICES][P4EST_DIM + 1];
          int8_t c = p4est_quadrant_child_id (quad);
          p4est_lnodes_code_t fc = lnodes->face_code[elem];
          int corner_is_hanging[P4EST_CHILDREN] = {0};
          p4est_gloidx_t cidx;
          p4est_quadrant_t parent;
          int8_t pc;
          int o = (c == 1 || c == 2 || c == 4 || c == 7);

          p4est_quadrant_parent (quad, &parent);
          pc = p4est_quadrant_child_id (&parent);
          cidx = p4est_lnodes_global_index (lnodes, E[c]);

#ifndef P4_TO_P8
          sims[0][0] = c ^ 0; sims[0][1] = c ^ 1; sims[0][2] = c ^ 3;
          sims[1][0] = c ^ 0; sims[1][1] = c ^ 2; sims[1][2] = c ^ 3;
#else
          sims[0][0] = c ^ 0; sims[0][1] = c ^ 1; sims[0][2] = c ^ 3; sims[0][3] = c ^ 7;
          sims[1][0] = c ^ 0; sims[1][1] = c ^ 1; sims[1][2] = c ^ 5; sims[1][3] = c ^ 7;
          sims[2][0] = c ^ 0; sims[2][1] = c ^ 2; sims[2][2] = c ^ 6; sims[2][3] = c ^ 7;
          sims[3][0] = c ^ 0; sims[3][1] = c ^ 2; sims[3][2] = c ^ 3; sims[3][3] = c ^ 7;
          sims[4][0] = c ^ 0; sims[4][1] = c ^ 4; sims[4][2] = c ^ 5; sims[4][3] = c ^ 7;
          sims[5][0] = c ^ 0; sims[5][1] = c ^ 4; sims[5][2] = c ^ 6; sims[5][3] = c ^ 7;
#endif

          if (fc) {
            p4est_lnodes_code_t work = fc >> P4EST_DIM;

            for (int d = 0; d < P4EST_DIM; d++, work >>= 1) {
              if (work & 1) {
                int f = p4est_corner_faces[c][d];
                int fcorner = p4est_corner_face_corners[c][f];
                int opp_fc = fcorner ^ (P4EST_HALF - 1);
                int opp = p4est_face_corners[f][opp_fc];
                double XYZ[3] = {0};
                p4est_quadrant_t node;

                corner_is_hanging[opp] = 1;

                p4est_quadrant_corner_node (&parent, opp, &node);
                XYZ[0] = ((double) node.x) / ((double) P4EST_ROOT_LEN);
                XYZ[1] = ((double) node.y) / ((double) P4EST_ROOT_LEN);
#ifdef P4_TO_P8
                XYZ[2] = ((double) node.z) / ((double) P4EST_ROOT_LEN);
#endif
                geom->X (geom, t, XYZ, sc_array_index (vertices, E[opp]));
              }
            }
#ifdef P4_TO_P8
            for (int d = 0; d < P4EST_DIM; d++, work >>= 1) {
              if (work & 1) {
                int e = p8est_corner_edges[c][d];
                int ec = p8est_corner_edge_corners[c][e];
                int opp_ec = ec ^ 1;
                int opp = p8est_edge_corners[e][opp_ec];
                double XYZ[3] = {0};
                p4est_quadrant_t node;

                corner_is_hanging[opp] = 1;

                p4est_quadrant_corner_node (&parent, opp, &node);
                XYZ[0] = ((double) node.x) / ((double) P4EST_ROOT_LEN);
                XYZ[1] = ((double) node.y) / ((double) P4EST_ROOT_LEN);
                XYZ[2] = ((double) node.z) / ((double) P4EST_ROOT_LEN);
                geom->X (geom, t, XYZ, sc_array_index (vertices, E[opp]));
              }
            }
#endif
          }
          for (int s = 0; s < P4EST_SIMPLICES; s++) {
            p4est_locidx_t *new_simplex;
            P4EST_ASSERT (!corner_is_hanging[sims[s][0]]);
            P4EST_ASSERT (!corner_is_hanging[sims[s][P4EST_DIM]]);
            if (corner_is_hanging[sims[s][1]]) {
              if (corner_is_hanging[sims[s][P4EST_DIM-1]]) {
                /* simplex on a hanging facet */
                if (quad->level >= 2) {
                  if ((sims[s][0] != pc) && (sims[s][P4EST_DIM-1] != (pc ^ (P4EST_CHILDREN - 1)))) {
                    /* only one child will have a simplex that does not
                       satisfy this condiiton */
                    continue;
                  }
                } else {
#ifdef P4_TO_P8
                  int p = s ^ 1;
                  P4EST_ASSERT (sims[p][1] == sims[s][1]);
                  if ((cidx > p4est_lnodes_global_index (lnodes, E[sims[s][1]]))
                      || (cidx > p4est_lnodes_global_index (lnodes, E[sims[s][2]]))
                      || (cidx > p4est_lnodes_global_index (lnodes, E[sims[p][2]]))
                      ) {
                    continue;
                  }
#else
                  if (cidx > p4est_lnodes_global_index (lnodes, E[sims[s][1]])) {
                    continue;
                  }
#endif
                }
              }
              else {
                /* simplex on a hanging edge */
                if (quad->level >= 2) {
                  if ((sims[s][0] != pc) && (sims[s][1] != (pc ^ (P4EST_CHILDREN - 1)))
                      && ((sims[s][0] ^ sims[s][1]) & (sims[s][0] ^ pc))) {
                    /* only one child will have a simplex that does not
                       satisfy this condiiton */
                    continue;
                  }
                } else {
                  if (cidx > p4est_lnodes_global_index (lnodes, E[sims[s][1]])) {
                    continue;
                  }
                }
              }
            }
            /* if we did not continue above, push the simplex */
            new_simplex = (p4est_locidx_t *) sc_array_push(simplices);
            for (int i = 0; i < P4EST_DIM + 1; i++) {
              new_simplex[i] = E[sims[s][i]];
            }
            if (o ^ (s & 1)) {
              /* simplex is inverted, swap the last two for correct order */
              p4est_locidx_t tmp = new_simplex[P4EST_DIM];

              new_simplex[P4EST_DIM] = new_simplex[P4EST_DIM - 1];
              new_simplex[P4EST_DIM - 1] = tmp;
            }
          }
        }
      }
    }
  }
  if (!ghost_layer) {
    p4est_ghost_destroy (ghost);
  }
  if (!geometry) {
    p4est_geometry_destroy (geom);
  }
  return snodes;
}

void
p4est_simplex_nodes_destroy (p4est_simplex_nodes_t *snodes)
{
  sc_array_destroy (snodes->simplices);
  sc_array_destroy (snodes->vertices);
  p4est_lnodes_destroy (snodes->lnodes);
  P4EST_FREE (snodes);
}

void
p4est_simplex_nodes_test_orientation (p4est_simplex_nodes_t *snodes)
{
  for (size_t elem = 0; elem < snodes->simplices->elem_count; elem++) {
    const p4est_locidx_t *E = sc_array_index(snodes->simplices, elem);
    double v[P4EST_DIM+1][P4EST_DIM];
    double A[P4EST_DIM][P4EST_DIM];
    double det;

    for (size_t i = 0; i < P4EST_DIM + 1; i++) {
      const double *w = (const double *) sc_array_index (snodes->vertices, E[i]);
      for (size_t d = 0; d < P4EST_DIM; d++) {
        v[i][d] = w[d];
      }
    }

    for (size_t i = 0; i < P4EST_DIM; i++) {
      for (size_t j = 0; j < P4EST_DIM; j++) {
        A[i][j] = v[i+1][j] - v[0][j];
      }
    }
#ifndef P4_TO_P8
    det = A[0][0] * A[1][1] - A[1][0] * A[0][1];
    SC_CHECK_ABORTF (det > 0, "Simplex %zu (vertices %zu, %zu, %zu, coords [%g, %g], [%g, %g], [%g, %g], det %g\n", elem, (size_t) E[0], (size_t) E[1], (size_t) E[2], v[0][0], v[0][1], v[1][0], v[1][1], v[2][0], v[2][1], det);
#else
    det = A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1])
        - A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0])
        + A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);
    SC_CHECK_ABORTF (det > 0, "Simplex %zu (vertices %zu, %zu, %zu, %zu, coords [%g, %g, %g], [%g, %g, %g], [%g, %g, %g], [%g, %g, %g], det %g\n", elem, (size_t) E[0], (size_t) E[1], (size_t) E[2], (size_t) E[3],
        v[0][0], v[0][1], v[0][2],
        v[1][0], v[1][1], v[1][2],
        v[2][0], v[2][1], v[2][2],
        v[3][0], v[3][1], v[3][2],
        det);
#endif
  }
}

int
main (int argc, char **argv)
{
  sc_options_t *opt;
  opts_t       opts = {0};
  p4est_t      *p4est;
  p4est_simplex_nodes_t *snodes;
  p4est_ghost_t *ghost_layer;
  p4est_geometry_t *geom = NULL;

  SC_CHECK_MPI (sc_MPI_Init (&argc, &argv));
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  opts.minlevel = 2;
  opts.maxlevel = 5;
  opts.conn = "unit";

  opt = sc_options_new (argv[0]);
  sc_options_add_int (opt, 'l', "minlevel", &opts.minlevel, opts.minlevel, "Lowest level");
  sc_options_add_int (opt, 'L', "maxlevel", &opts.maxlevel, opts.maxlevel, "Highest level");
  sc_options_add_string (opt, 'v', "vtk", &opts.vtk, opts.vtk, "VTK basename");
  sc_options_add_string (opt, 'c', "conn", &opts.conn, opts.conn, "Name of the connectivity");

  if (sc_options_parse(p4est_package_id, SC_LP_DEFAULT, opt, argc, argv) != argc) {
    P4EST_GLOBAL_LERROR ("Error parsing options");
    sc_options_print_usage (p4est_package_id, SC_LP_ERROR, opt, NULL);
    sc_options_destroy (opt);
    sc_finalize ();
    SC_CHECK_MPI (sc_MPI_Finalize ());
    return 1;
  }

  p4est = create_p4est_from_opts (&opts);
  ghost_layer = p4est_ghost_new (p4est, P4EST_CONNECT_FULL);
#ifndef P4_TO_P8
  if (opts.conn && !strcmp(opts.conn, "disk2d")) {
    geom = p4est_geometry_new_disk2d (p4est->connectivity, 1.0, 2.0);
  } else if (opts.conn && !strcmp(opts.conn, "icosahedron")) {
    geom = p4est_geometry_new_icosahedron (p4est->connectivity, 1.0);
  } else if (opts.conn && !strcmp(opts.conn, "shell2d")) {
    geom = p4est_geometry_new_shell2d (p4est->connectivity, 2.0, 1.0);
  }
#else
  if (opts.conn && !strcmp(opts.conn, "shell")) {
    geom = p8est_geometry_new_shell (p4est->connectivity, 2.0, 1.0);
  } else if (opts.conn && !strcmp(opts.conn, "sphere")) {
    geom = p8est_geometry_new_sphere (p4est->connectivity, 3.0, 2.0, 1.0);
  } else if (opts.conn && !strcmp(opts.conn, "torus")) {
    geom = p8est_geometry_new_torus (p4est->connectivity, 1.0, 2.0, 4.0);
  }
#endif

  snodes = p4est_simplex_nodes_create (p4est, ghost_layer, geom);

#ifndef P4_TO_P8
  if (strcmp(opts.conn, "corner")
      && strcmp(opts.conn, "cubed")
      && strcmp(opts.conn, "icosahedron")
      && strcmp(opts.conn, "moebius")
      && strcmp(opts.conn, "periodic")
      && strcmp(opts.conn, "pillow")
      && strcmp(opts.conn, "rotwrap")
      && strcmp(opts.conn, "shell2d")
      && strcmp(opts.conn, "star")) {
    p4est_simplex_nodes_test_orientation (snodes);
  }
#else
  if (strcmp(opts.conn, "periodic")
      && strcmp(opts.conn, "rotwrap")
      && strcmp(opts.conn, "twowrap")) {
    p4est_simplex_nodes_test_orientation (snodes);
  }
#endif

  if (opts.vtk != NULL) {
    p4est_vtk_context_t *vtk = p4est_vtk_context_new (p4est, opts.vtk);

    vtk = p4est_vtk_write_header_simplices (vtk, snodes->simplices, snodes->vertices);
    p4est_vtk_write_footer (vtk);
  }

  p4est_simplex_nodes_destroy (snodes);
  if (geom) {
    p4est_geometry_destroy (geom);
  }

  {
    p4est_connectivity_t *conn = p4est->connectivity;

    p4est_ghost_destroy (ghost_layer);
    p4est_destroy (p4est);
    p4est_connectivity_destroy (conn);
  }

  sc_options_destroy (opt);
  sc_finalize ();
  SC_CHECK_MPI (sc_MPI_Finalize ());
  return 0;
}
