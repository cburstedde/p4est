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

/*
 * Usage: p4est_tnodes <connectivity> <level> [<options>]
 *        possible connectivities:
 *        o unit      Refinement on the unit square.
 *        o three     Refinement on a forest with three trees.
 *        o moebius   Refinement on a 5-tree Moebius band.
 *        o star      Refinement on a 6-tree star shaped domain.
 *        o periodic  Refinement on the unit square with all-periodic b.c.
 *        o rotwrap   Refinement on the unit square with weird periodic b.c.
 *        o cubed     Refinement on the 2D cubed sphere.
 *        o disk      Refinement on a 5-tree flat disk or square.
 *        o pdisk     Refinement on 5-tree flat disk or square, periodic b.c.
 *        o icosahedron   Refine on an icosahedron embedded in 3D space.
 *        options can be a string containing "N" for omitting VTK output.
 *        options can be a string containing "U" for uniform refinement.
 */

#ifndef P4_TO_P8
#include <p4est_bits.h>
#include <p4est_extended.h>
#include <p4est_tnodes.h>
#include <p4est_vtk.h>
#else
#include <p8est_bits.h>
#include <p8est_extended.h>
#include <p8est_tnodes.h>
#include <p8est_vtk.h>
#endif
#include <sc_statistics.h>

typedef enum
{
  P4EST_CONFIG_NULL,
#ifndef P4_TO_P8
  P4EST_CONFIG_UNIT,
  P4EST_CONFIG_THREE,
  P4EST_CONFIG_MOEBIUS,
  P4EST_CONFIG_STAR,
  P4EST_CONFIG_PERIODIC,
  P4EST_CONFIG_ROTWRAP,
  P4EST_CONFIG_CUBED,
  P4EST_CONFIG_DISK,
  P4EST_CONFIG_PDISK,
  P4EST_CONFIG_ICOSAHEDRON
#else
  P8EST_CONFIG_UNIT,
  P8EST_CONFIG_PERIODIC,
  P8EST_CONFIG_ROTWRAP,
  P8EST_CONFIG_TWOCUBES,
  P8EST_CONFIG_TWOWRAP,
  P8EST_CONFIG_ROTCUBES,
  P8EST_CONFIG_SHELL,
  P8EST_CONFIG_SPHERE,
  P8EST_CONFIG_TORUS
#endif
}
simple_config_t;

enum tnodes_stats_abbr
{
  MESH_TNODES_ABBR_MEM = 0,
  MESH_TNODES_ABBR_NELEM = 1,
  MESH_TNODES_ABBR_LNODES = 2,
  MESH_TNODES_ABBR_EXPORT = 3,
  MESH_TNODES_ABBR_EXPORT2 = 4
};

enum tnodes_stats_names
{
  MESH_TNODES_1Q1_MEM,
  MESH_TNODES_1Q1_NELEM,
  MESH_TNODES_1Q1_LNODES,
  MESH_TNODES_1Q1_EXPORT,
  MESH_TNODES_Q2I_MEM,
  MESH_TNODES_Q2I_NELEM,
  MESH_TNODES_Q2I_LNODES,
  MESH_TNODES_Q2I_EXPORT,
  MESH_TNODES_Q2D_EXPORT,
  MESH_TNODES_2Q1_MEM,
  MESH_TNODES_2Q1_NELEM,
  MESH_TNODES_2Q1_LNODES,
  MESH_TNODES_2Q1_EXPORT,
  MESH_TNODES_STATS_COUNT
};

static sc_statinfo_t stats[MESH_TNODES_STATS_COUNT];

static int          novtk = 0;
static int          uniform = 0;

typedef struct
{
  int                 dummy;
}
user_data_t;

typedef struct
{
  sc_MPI_Comm         mpicomm;
  int                 mpisize;
  int                 mpirank;
}
mpi_context_t;

/* Global variable is not recommended, yes.
 * To clean this up in your own application, please put your context
 * into the p4est->user_pointer field and access it from the callbacks.
 */
static int          refine_level = 0;
static const char  *configuration = NULL;

/* copy variable string */
void
tnodes_stats_set1 (sc_statinfo_t *stats, double value, const char *variable)
{
  sc_stats_set1_ext (stats, value, variable, 1, -2, -3);
}

static void
init_fn (p4est_t *p4est, p4est_topidx_t which_tree,
         p4est_quadrant_t *quadrant)
{
  user_data_t        *data = (user_data_t *) quadrant->p.user_data;
  data->dummy = -1;
}

static int
refine_uniform (p4est_t *p4est, p4est_topidx_t which_tree,
                p4est_quadrant_t *quadrant)
{
  return (int) quadrant->level < refine_level;
}

static int
refine_once (p4est_t *p4est, p4est_topidx_t which_tree,
             p4est_quadrant_t *quadrant)
{
  return 1;
}

static int
refine_normal (p4est_t *p4est, p4est_topidx_t which_tree,
               p4est_quadrant_t *quadrant)
{
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

static void
tnodes_run_Q1 (p4est_t *p4est, p4est_geometry_t *geom, p4est_ghost_t *ghost,
               const char *name, sc_statinfo_t *lstats)
{
  int                 retval;
  double              lntime, Q1time;
  size_t              mem_lnodes, mem_tnodes;
  char                concat[BUFSIZ];
  p4est_lnodes_t     *ln;
  p4est_tnodes_t     *tm;
  p4est_vtk_context_t *cont;

  P4EST_GLOBAL_PRODUCTIONF ("Example tnodes run %s\n", name);

  P4EST_ASSERT (p4est != NULL);
  P4EST_ASSERT (ghost != NULL);

  /* remember local element count */
  snprintf (concat, BUFSIZ, "%s_%s", name, "NELEM");
  tnodes_stats_set1 (lstats + MESH_TNODES_ABBR_NELEM,
                     (double) p4est->local_num_quadrants, concat);

  /* generate lnodes of degree 1 */
  lntime = sc_MPI_Wtime ();
  ln = p4est_lnodes_new (p4est, ghost, 1);
  lntime = sc_MPI_Wtime () - lntime;
  mem_lnodes = p4est_lnodes_memory_used (ln);
  P4EST_INFOF ("Memory used by %s lnodes: %lld bytes\n",
               name, (long long) mem_lnodes);
  snprintf (concat, BUFSIZ, "%s_%s", name, "LNODES");
  tnodes_stats_set1 (lstats + MESH_TNODES_ABBR_LNODES, lntime, concat);

  /* export Q1-based simplex mesh */
  Q1time = sc_MPI_Wtime ();
  tm = p4est_tnodes_new_Q1_P1 (p4est, ln);
  Q1time = sc_MPI_Wtime () - Q1time;
  mem_tnodes = p4est_tnodes_memory_used (tm);
  P4EST_INFOF ("Memory used by %s tnodes: %lld bytes\n",
               name, (long long) mem_tnodes);
  snprintf (concat, BUFSIZ, "%s_%s", name, "EXPORT");
  tnodes_stats_set1 (lstats + MESH_TNODES_ABBR_EXPORT, Q1time, concat);

  /* remember memory usage */
  snprintf (concat, BUFSIZ, "%s_%s", name, "MEM");
  tnodes_stats_set1 (lstats + MESH_TNODES_ABBR_MEM,
                     mem_lnodes + mem_tnodes, concat);

  if (!novtk) {
    /* write VTK output */

    snprintf (concat, BUFSIZ, "%s_%s_%s_%s_%02d%s", P4EST_STRING, "tnodes",
              configuration, name, refine_level, uniform ? "U" : "R");
    cont = p4est_vtk_context_new (p4est, concat);
    SC_CHECK_ABORT (cont != NULL, "Open VTK context");
    p4est_vtk_context_set_lnodes (cont, ln);
    p4est_vtk_context_set_geom (cont, geom);
    p4est_vtk_context_set_continuous (cont, 1);

    /* beware: values < 1. cause a lot more mesh nodes */
    p4est_vtk_context_set_scale (cont, 1.);

    cont = p4est_vtk_write_header_tnodes (cont, tm);
    SC_CHECK_ABORT (cont != NULL, "Write tnodes VTK header");
    cont = p4est_vtk_write_cell_dataf (cont, 1, 1, 1, 0, 0, 0, cont);
    SC_CHECK_ABORT (cont != NULL, "Write tnodes VTK cells");
    retval = p4est_vtk_write_footer (cont);
    SC_CHECK_ABORT (!retval, "Close VTK context");
  }

  /* free triangle mesh */
  p4est_tnodes_destroy (tm);
  p4est_lnodes_destroy (ln);
}

static void
compare_both_Q2_constructions (p4est_t *p4est,
                               p4est_tnodes_t *tm, p4est_tnodes_t *tl)
{
  int                 k;
  int8_t              lvm, lvl;
  int8_t             *tms, *tls;
  p4est_locidx_t      ns, s;

  P4EST_ASSERT (p4est != NULL);
  P4EST_ASSERT (tm != NULL);
  P4EST_ASSERT (tl != NULL);

  P4EST_ASSERT (tm->global_tcount == tl->global_tcount);

  ns = tm->local_tcount[p4est->mpirank];
  P4EST_ASSERT (ns == tl->local_tcount[p4est->mpirank]);
  P4EST_ASSERT (ns == tm->local_element_offset[p4est->local_num_quadrants]);
  P4EST_ASSERT (ns == tl->local_element_offset[p4est->local_num_quadrants]);

  for (s = 0; s < ns; ++s) {
    tms = (int8_t *) sc_array_index (tm->simplices, s);
    tls = (int8_t *) sc_array_index (tl->simplices, s);
    for (k = 0; k <= P4EST_DIM; ++k) {
      SC_CHECK_ABORTF (tms[k] == tls[k],
                       "Simplex mismatch %ld at %d: %d, %d\n",
                       (long) s, k, tms[k], tls[k]);
    }
    lvm = *(int8_t *) sc_array_index (tm->simplex_level, s);
    lvl = *(int8_t *) sc_array_index (tl->simplex_level, s);
    SC_CHECK_ABORTF (lvm == lvl,
                     "Level mismatch %ld: %d, %d\n", (long) s, lvm, lvl);

  }
  SC_CHECK_ABORT (sc_array_is_equal (tm->element_bits,
                                     tl->element_bits), "Bits mismatch");
}

static void
tnodes_run_Q2_both (p4est_t *p4est, p4est_geometry_t *geom,
                    p4est_ghost_t *ghost)
{
  int                 retval;
  double              lntime, Q2Itime, Q2Dtime;
  size_t              mem_lnodes, mem_tdnodes, mem_tinodes;
  char                concat[BUFSIZ];
  const char         *name;
  p4est_lnodes_t     *ln;
  p4est_tnodes_t     *tnd;
  p4est_tnodes_t     *tni;
  p4est_vtk_context_t *cont;
  sc_statinfo_t      *lstats;

  P4EST_ASSERT (p4est != NULL);
  P4EST_ASSERT (ghost != NULL);

  /* set reporting context for indirect construction */
  lstats = stats + MESH_TNODES_Q2I_MEM;
  name = "Q2I";
  P4EST_GLOBAL_PRODUCTIONF ("Example tnodes run %s\n", name);

  /* remember local element count */
  snprintf (concat, BUFSIZ, "%s_%s", name, "NELEM");
  tnodes_stats_set1 (lstats + MESH_TNODES_ABBR_NELEM,
                     (double) p4est->local_num_quadrants, concat);

  /* generate lnodes of degree 2 */
  lntime = sc_MPI_Wtime ();
  ln = p4est_lnodes_new (p4est, ghost, 2);
  lntime = sc_MPI_Wtime () - lntime;
  mem_lnodes = p4est_lnodes_memory_used (ln);
  P4EST_INFOF ("Memory used by %s lnodes: %lld bytes\n",
               name, (long long) mem_lnodes);
  snprintf (concat, BUFSIZ, "%s_%s", name, "LNODES");
  tnodes_stats_set1 (lstats + MESH_TNODES_ABBR_LNODES, lntime, concat);

  /* export Q2 simplex mesh by indirect method */
  Q2Itime = sc_MPI_Wtime ();
  tni = p4est_tnodes_new_Q2_P1_ind (p4est, ln);
  Q2Itime = sc_MPI_Wtime () - Q2Itime;
  mem_tinodes = p4est_tnodes_memory_used (tni);
  P4EST_INFOF ("Memory used by %s tnodes: %lld bytes\n",
               name, (long long) mem_tinodes);
  snprintf (concat, BUFSIZ, "%s_%s", name, "EXPORT");
  tnodes_stats_set1 (lstats + MESH_TNODES_ABBR_EXPORT, Q2Itime, concat);

  /* remember memory usage */
  snprintf (concat, BUFSIZ, "%s_%s", name, "MEM");
  tnodes_stats_set1 (lstats + MESH_TNODES_ABBR_MEM,
                     mem_lnodes + mem_tinodes, concat);

  /* set reporting context for direct construction */
  name = "Q2D";
  P4EST_GLOBAL_PRODUCTIONF ("Example tnodes run %s\n", name);

  /* export Q2 simplex mesh by direct method */
  Q2Dtime = sc_MPI_Wtime ();
  tnd = p4est_tnodes_new_Q2_P1 (p4est, ln);
  Q2Dtime = sc_MPI_Wtime () - Q2Dtime;
  mem_tdnodes = p4est_tnodes_memory_used (tnd);
  SC_CHECK_ABORTF (mem_tinodes == mem_tdnodes,
                   "Memory Q2I %lld vs. Q2D %lld",
                   (long long) mem_tinodes, (long long) mem_tdnodes);
  snprintf (concat, BUFSIZ, "%s_%s", name, "EXPORT");
  tnodes_stats_set1 (lstats + MESH_TNODES_ABBR_EXPORT2, Q2Dtime, concat);

  /* verify these two Q2 constructions are indeed identical */
  compare_both_Q2_constructions (p4est, tnd, tni);

  if (!novtk) {
    /* write VTK output */

    snprintf (concat, BUFSIZ, "%s_%s_%s_%s_%02d%s", P4EST_STRING, "tnodes",
              configuration, name, refine_level, uniform ? "U" : "R");
    cont = p4est_vtk_context_new (p4est, concat);
    SC_CHECK_ABORT (cont != NULL, "Open VTK context");
    p4est_vtk_context_set_lnodes (cont, ln);
    p4est_vtk_context_set_geom (cont, geom);
    p4est_vtk_context_set_continuous (cont, 1);

    /* beware: values < 1. cause a lot more mesh nodes */
    p4est_vtk_context_set_scale (cont, 1.);

    cont = p4est_vtk_write_header_tnodes (cont, tnd);
    SC_CHECK_ABORT (cont != NULL, "Write tnodes VTK header");
    cont = p4est_vtk_write_cell_dataf (cont, 1, 1, 1, 0, 0, 0, cont);
    SC_CHECK_ABORT (cont != NULL, "Write tnodes VTK cells");
    retval = p4est_vtk_write_footer (cont);
    SC_CHECK_ABORT (!retval, "Close VTK context");
  }

  /* free triangle mesh */
  p4est_tnodes_destroy (tnd);
  p4est_tnodes_destroy (tni);
  p4est_lnodes_destroy (ln);
}

static void
forest_run (mpi_context_t *mpi,
            p4est_connectivity_t *connectivity, p4est_geometry_t *geom)
{
  int                 l;
  unsigned            crc;
  char                msg[BUFSIZ];
  p4est_t            *p4est;
  p4est_ghost_t      *ghost;

  P4EST_GLOBAL_PRODUCTIONF ("Example forest run uniform %d\n", uniform);

  /* create new coarse p4est from specified connectivity */
  p4est = p4est_new_ext (mpi->mpicomm, connectivity, 0, 0, 1,
                         sizeof (user_data_t), init_fn, NULL);
  snprintf (msg, BUFSIZ, P4EST_STRING "_tnodes_partitioned_%02d", 0);
  if (!novtk) {
    p4est_vtk_write_file (p4est, geom, msg);
  }

  /* non-recursive refinement loop */
  for (l = 1; l <= refine_level; ++l) {
    /* refine */
    p4est_refine (p4est, 0, uniform ? refine_uniform : refine_normal,
                  init_fn);
    snprintf (msg, BUFSIZ, P4EST_STRING "_tnodes_refined_%02d", l);
    if (!novtk) {
      p4est_vtk_write_file (p4est, geom, msg);
    }

    if (!uniform) {
      /* balance */
      p4est_balance (p4est, P4EST_CONNECT_FULL, init_fn);
      snprintf (msg, BUFSIZ, P4EST_STRING "_tnodes_balanced_%02d", l);
      if (!novtk) {
        p4est_vtk_write_file (p4est, geom, msg);
      }
    }

    /* partition */
    p4est_partition (p4est, 0, NULL);
    snprintf (msg, BUFSIZ, P4EST_STRING "_tnodes_partitioned_%02d", l);
    if (!novtk) {
      p4est_vtk_write_file (p4est, geom, msg);
    }
  }
  crc = p4est_checksum (p4est);

  /* print and verify forest checksum */
  P4EST_GLOBAL_STATISTICSF ("Example forest %s checksum 0x%08x\n",
                            uniform ? "uniform" : "adapted", crc);

  /* create ghost layer and triangle mesh from Q2 nodes */
  ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FULL);
  tnodes_run_Q1 (p4est, geom, ghost, "1Q1", stats + MESH_TNODES_1Q1_MEM);
  tnodes_run_Q2_both (p4est, geom, ghost);
  p4est_ghost_destroy (ghost);

  /* refine forest uniformly by one level */
  p4est_refine (p4est, 0, refine_once, init_fn);
  P4EST_GLOBAL_STATISTICSF ("Example forest %s checksum 0x%08x\n",
                            "again", crc);

  /* create ghost layer and triangle mesh from Q1 nodes */
  ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FULL);
  tnodes_run_Q1 (p4est, geom, ghost, "2Q1", stats + MESH_TNODES_2Q1_MEM);
  p4est_ghost_destroy (ghost);

  /* destroy the p4est structure */
  p4est_destroy (p4est);
}

static void
verify_aux (void)
{
#ifdef P4EST_ENABLE_DEBUG
  int                 p, c;
  int                 kp, kc;
  int                 nfc;
  int                 hf;
#ifdef P4_TO_P8
  int                 i, j, k;
#endif
  p4est_lnodes_code_t fc, pfc, cfc[P4EST_CHILDREN];

  /* verify computation of parent simplex index */
  for (p = 0; p < P4EST_CHILDREN; ++p) {
    for (c = 0; c < P4EST_CHILDREN; ++c) {
      for (kc = 0; kc < (P4EST_DIM - 1) * P4EST_DIM; ++kc) {
        kp = p4est_tnodes_simplex_parent (p, c, kc);
        P4EST_ASSERT (p4est_tnodes_simplex_parent_is_valid (p, kp, c, kc));
        P4EST_ASSERT (p4est_tnodes_simplex_parent_is_valid
                      (0, kp, c ^ p, kc));
      }
    }
  }

  /* verify consistency of Q1 and Q2 simplex counts */
  nfc = 0;
  for (hf = 0; hf < P4EST_CHILDREN; ++hf) {

    /* iterate over all possible face codes */
#ifndef P4_TO_P8
    fc = hf << P4EST_DIM;
#else
    for (k = 0; k < 2; ++k) {
      if (!k && (hf & 3)) {
        continue;
      }
      for (j = 0; j < 2; ++j) {
        if (!j && (hf & 5)) {
          continue;
        }
        for (i = 0; i < 2; ++i) {
          if (!i && (hf & 6)) {
            continue;
          }
#if 0
        }
      }
    }
#endif
    fc = ((((k << 2) | (j << 1) | i) << P4EST_DIM) | hf) << P4EST_DIM;
#endif
    ++nfc;

    /* loop through all possible parent positions */
    for (p = 0; p < P4EST_CHILDREN; ++p) {
      pfc = fc ? (fc | p) : 0;
      p4est_lnodes_derive_child_codes (p, pfc, cfc);

      /* loop through all children of p */
      kc = 0;
      for (c = 0; c < P4EST_CHILDREN; ++c) {
        kc += p4est_tnodes_quadrant_Q1_simplices (cfc[c], p);
      }
      kp = p4est_tnodes_quadrant_Q2_simplices (pfc);
      P4EST_ASSERT (kc == kp);
    }

#ifdef P4_TO_P8
#if 0
    {
      {
        {
#endif
        }
      }
    }
#endif
  }

  /* verify total number of possible face codes */
  P4EST_ASSERT (nfc == (P4EST_DIM == 2 ? 4 : 18));
#endif
}

int
main (int argc, char **argv)
{
  int                 i;
  int                 mpiret;
  int                 wrongusage;
  const char         *usage;
  mpi_context_t       mpi_context, *mpi = &mpi_context;
  p4est_connectivity_t *connectivity;
  p4est_geometry_t   *geometry;
  simple_config_t     config;

  /* initialize MPI and p4est internals */
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  mpi->mpicomm = sc_MPI_COMM_WORLD;
  mpiret = sc_MPI_Comm_size (mpi->mpicomm, &mpi->mpisize);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (mpi->mpicomm, &mpi->mpirank);
  SC_CHECK_MPI (mpiret);

  sc_init (mpi->mpicomm, 1, 1, NULL, SC_LP_APPLICATION);
  p4est_init (NULL, SC_LP_APPLICATION);
  memset (stats, 0, sizeof (stats));

  /* process command line arguments */
  usage =
    "Arguments: <connectivity> <level> [<options>]\n"
    "   The connectivity can be any of\n"
#ifndef P4_TO_P8
    "      unit|three|moebius|star|periodic|rotwrap|\n"
    "         cubed|disk|pdisk|icosahedron\n"
#else
    "      unit|periodic|rotwrap|twocubes|twowrap|rotcubes|\n"
    "         shell|sphere|torus\n"
#endif
    "   Level controls the maximum depth of refinement\n"
    "   Options may be empty or contain N for no VTK output\n"
    "   Options may be empty or contain U for uniform refinement\n";
  wrongusage = 0;
  config = P4EST_CONFIG_NULL;
  if (!wrongusage && (argc < 3 || argc > 4)) {
    wrongusage = 1;
  }
  if (!wrongusage) {
    configuration = argv[1];
    if (!strcmp (argv[1], "unit")) {
#ifndef P4_TO_P8
      config = P4EST_CONFIG_UNIT;
#else
      config = P8EST_CONFIG_UNIT;
#endif
    }
#ifndef P4_TO_P8
    else if (!strcmp (argv[1], "three")) {
      config = P4EST_CONFIG_THREE;
    }
    else if (!strcmp (argv[1], "moebius")) {
      config = P4EST_CONFIG_MOEBIUS;
    }
    else if (!strcmp (argv[1], "star")) {
      config = P4EST_CONFIG_STAR;
    }
    else if (!strcmp (argv[1], "periodic")) {
      config = P4EST_CONFIG_PERIODIC;
    }
    else if (!strcmp (argv[1], "rotwrap")) {
      config = P4EST_CONFIG_ROTWRAP;
    }
    else if (!strcmp (argv[1], "cubed")) {
      config = P4EST_CONFIG_CUBED;
    }
    else if (!strcmp (argv[1], "disk")) {
      config = P4EST_CONFIG_DISK;
    }
    else if (!strcmp (argv[1], "pdisk")) {
      config = P4EST_CONFIG_PDISK;
    }
    else if (!strcmp (argv[1], "icosahedron")) {
      config = P4EST_CONFIG_ICOSAHEDRON;
    }
#else
    else if (!strcmp (argv[1], "periodic")) {
      config = P8EST_CONFIG_PERIODIC;
    }
    else if (!strcmp (argv[1], "rotwrap")) {
      config = P8EST_CONFIG_ROTWRAP;
    }
    else if (!strcmp (argv[1], "twocubes")) {
      config = P8EST_CONFIG_TWOCUBES;
    }
    else if (!strcmp (argv[1], "twowrap")) {
      config = P8EST_CONFIG_TWOWRAP;
    }
    else if (!strcmp (argv[1], "rotcubes")) {
      config = P8EST_CONFIG_ROTCUBES;
    }
    else if (!strcmp (argv[1], "shell")) {
      config = P8EST_CONFIG_SHELL;
    }
    else if (!strcmp (argv[1], "sphere")) {
      config = P8EST_CONFIG_SPHERE;
    }
    else if (!strcmp (argv[1], "torus")) {
      config = P8EST_CONFIG_TORUS;
    }
#endif
    else {
      wrongusage = 1;
      P4EST_GLOBAL_LERROR ("Unknown connectivity\n");
    }
  }
  if (!wrongusage) {
    refine_level = atoi (argv[2]);
    if (refine_level < 0 || refine_level > P4EST_QMAXLEVEL) {
      wrongusage = 1;
      P4EST_GLOBAL_LERROR ("Refinement level out of range\n");
    }
  }
  if (!wrongusage && argc >= 4) {
    if (strchr (argv[3], 'N') != NULL) {
      novtk = 1;
      P4EST_GLOBAL_PRODUCTION ("Option: deactivating VTK output\n");
    }
    if (strchr (argv[3], 'U') != NULL) {
      uniform = 1;
      P4EST_GLOBAL_PRODUCTION ("Option: selecting uniform refinement\n");
    }
  }
  if (wrongusage) {
    P4EST_GLOBAL_LERROR (usage);
    sc_abort_collective ("Usage error");
  }

  /* create connectivity and forest structures */
  connectivity = NULL;
  geometry = NULL;
  if (0) {
  }
#ifndef P4_TO_P8
  else if (config == P4EST_CONFIG_THREE) {
    connectivity = p4est_connectivity_new_corner ();
  }
  else if (config == P4EST_CONFIG_MOEBIUS) {
    connectivity = p4est_connectivity_new_moebius ();
  }
  else if (config == P4EST_CONFIG_STAR) {
    connectivity = p4est_connectivity_new_star ();
  }
  else if (config == P4EST_CONFIG_PERIODIC) {
    connectivity = p4est_connectivity_new_periodic ();
  }
  else if (config == P4EST_CONFIG_ROTWRAP) {
    connectivity = p4est_connectivity_new_rotwrap ();
  }
  else if (config == P4EST_CONFIG_CUBED) {
    connectivity = p4est_connectivity_new_cubed ();
  }
  else if (config == P4EST_CONFIG_DISK) {
    connectivity = p4est_connectivity_new_disk (0, 0);
  }
  else if (config == P4EST_CONFIG_PDISK) {
    connectivity = p4est_connectivity_new_disk (1, 1);
  }
  else if (config == P4EST_CONFIG_ICOSAHEDRON) {
    connectivity = p4est_connectivity_new_icosahedron ();
    geometry = p4est_geometry_new_icosahedron (connectivity, 1.);
  }
#else
  else if (config == P8EST_CONFIG_PERIODIC) {
    connectivity = p8est_connectivity_new_periodic ();
  }
  else if (config == P8EST_CONFIG_ROTWRAP) {
    connectivity = p8est_connectivity_new_rotwrap ();
  }
  else if (config == P8EST_CONFIG_TWOCUBES) {
    connectivity = p8est_connectivity_new_twocubes ();
  }
  else if (config == P8EST_CONFIG_TWOWRAP) {
    connectivity = p8est_connectivity_new_twowrap ();
  }
  else if (config == P8EST_CONFIG_ROTCUBES) {
    connectivity = p8est_connectivity_new_rotcubes ();
  }
  else if (config == P8EST_CONFIG_SHELL) {
    connectivity = p8est_connectivity_new_shell ();
    geometry = p8est_geometry_new_shell (connectivity, 1., 1.5);
  }
  else if (config == P8EST_CONFIG_SPHERE) {
    connectivity = p8est_connectivity_new_sphere ();
    geometry = p8est_geometry_new_sphere (connectivity, 1., .6, .3);
  }
  else if (config == P8EST_CONFIG_TORUS) {
    connectivity = p8est_connectivity_new_torus (8);
    geometry = p8est_geometry_new_torus (connectivity, 0.44, 1.0, 3.0);
  }
#endif
  else {
#ifndef P4_TO_P8
    connectivity = p4est_connectivity_new_unitsquare ();
#else
    connectivity = p8est_connectivity_new_unitcube ();
#endif
  }
  if (geometry == NULL) {
    geometry = p4est_geometry_new_connectivity (connectivity);
  }

  /* verify auxiliary functions */
  verify_aux ();

  /* run mesh tests */
  forest_run (mpi,              /* mpi context */
              connectivity,     /* p4est connectivity */
              geometry);        /* used for VTK output */

  /* compute and report statistics */
  sc_stats_compute (mpi->mpicomm, MESH_TNODES_STATS_COUNT, stats);
  sc_stats_print (p4est_get_package_id (), SC_LP_PRODUCTION,
                  MESH_TNODES_STATS_COUNT, stats, 1, 1);

  /* clean up and exit */
  if (geometry != NULL) {
    p4est_geometry_destroy (geometry);
  }
  p4est_connectivity_destroy (connectivity);
  for (i = 0; i < MESH_TNODES_STATS_COUNT; ++i) {
    sc_stats_reset (stats + i, 1);
  }
  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
