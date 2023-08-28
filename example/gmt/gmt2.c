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

#include <p4est_extended.h>
#include <p4est_vtk.h>
#include <sc_options.h>
#include "gmt_models.h"

typedef struct global
{
  int                 minlevel;
  int                 maxlevel;
  int                 resolution;
  int                 synthetic;
  int                 latlongno;
  sc_MPI_Comm         mpicomm;
  p4est_t            *p4est;
  p4est_gmt_model_t  *model;
}
global_t;

static int
setup_model (global_t * g)
{
  /* this function populates the model on successful initialization */
  P4EST_ASSERT (g->model == NULL);

  /* initialize model depending on command line options */
  if (g->synthetic >= 0) {
    switch (g->synthetic) {
    case 0:
      g->model = p4est_gmt_model_synth_new (g->synthetic);
      break;
    default:
      P4EST_GLOBAL_LERROR ("Synthetic model number exceeded\n");
    }
  }
  else if (g->latlongno >= 0) {
    model_latlong_params_t ap;

    switch (g->latlongno) {
    case 0:
      ap.latitude[0] = -50.;
      ap.latitude[1] = 0.;
      ap.longitude[0] = 0.;
      ap.longitude[1] = 60.;
      ap.resolution = g->resolution;
      ap.load_filename = "africa.gmt.data";
      ap.output_prefix = "africa";

      /* load data (possibly GMT, or file, or synthetic) */
      g->model = p4est_gmt_model_latlong_new (&ap);
      break;
    default:
      P4EST_GLOBAL_LERROR ("Latitute-longitude model number exceeded\n");
    }
  }

  /* on successful initalization the global model is set */
  return g->model == NULL ? -1 : 0;
}

void
run_program (global_t * g)
{
  const size_t        quad_data_size = 0;

  /* create mesh */
  g->p4est = p4est_new_ext (g->mpicomm, g->model->conn, 0, g->minlevel, 1,
                            quad_data_size, NULL, g);

  /* run mesh refinement based on data */

  /* output refined mesh */

  p4est_vtk_write_file (g->p4est, g->model->model_geom,
                        g->model->output_prefix);

  /* cleanup */
  p4est_destroy (g->p4est);
}

static int
usagerrf (sc_options_t * opt, const char *fmt, ...)
{
  va_list             ap;
  char                msg[BUFSIZ];

  va_start (ap, fmt);
  vsnprintf (msg, BUFSIZ, fmt, ap);
  va_end (ap);

  P4EST_GLOBAL_LERROR ("ERROR/\n");
  P4EST_GLOBAL_LERRORF ("ERROR: %s\n", msg);
  P4EST_GLOBAL_LERROR ("ERROR\\\n");
  return 1;
}

static int
usagerr (sc_options_t * opt, const char *msg)
{
  return usagerrf (opt, "%s", msg);
}

int
main (int argc, char **argv)
{
  int                 mpiret;
  int                 ue, fa;
  sc_options_t       *opt;
  global_t            sg, *g = &sg;

  /* initialize MPI */
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  /* initialize global context */
  memset (g, 0, sizeof (*g));
  g->mpicomm = sc_MPI_COMM_WORLD;

  /* set global logging options for p4est */
  sc_init (g->mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  /* initialize global application state */
  opt = sc_options_new (argv[0]);
  sc_options_add_int (opt, 'l', "minlevel", &g->minlevel, 0,
                      "Minimum refinement level");
  sc_options_add_int (opt, 'L', "maxlevel", &g->maxlevel, P4EST_QMAXLEVEL,
                      "Maximum refinement level");
  sc_options_add_int (opt, 'r', "resolution", &g->resolution, 0,
                      "Level of resolution (model specific)");
  sc_options_add_int (opt, 'S', "synthetic", &g->synthetic, -1,
                      "Choose specific synthetic model");
  sc_options_add_int (opt, 'M', "latlongno", &g->latlongno, -1,
                      "Choose specific latitude-longitude model");

  /* proceed in run-once loop for cleaner error checking */
  ue = 0;
  do {
    /* parse command line and assign configuration variables */
    fa = sc_options_parse (p4est_package_id, SC_LP_DEFAULT, opt, argc, argv);
    if (fa < 0 || fa != argc) {
      ue = usagerr (opt, "invalid option format or non-option argument");
      break;
    }
    P4EST_GLOBAL_PRODUCTIONF ("Manifold dimension is %d\n", P4EST_DIM);
    sc_options_print_summary (p4est_package_id, SC_LP_PRODUCTION, opt);

    /* check consistency of parameters */
    if (g->minlevel < 0 || g->minlevel > P4EST_QMAXLEVEL) {
      ue = usagerrf (opt, "minlevel not between 0 and %d", P4EST_QMAXLEVEL);
    }
    if (g->maxlevel < g->minlevel || g->maxlevel > P4EST_QMAXLEVEL) {
      ue = usagerrf (opt, "maxlevel not between minlevel and %d",
                     P4EST_QMAXLEVEL);
    }
    if (g->synthetic >= 0 && g->latlongno >= 0) {
      ue = usagerrf (opt, "set only one of the synthetic and latlong models");
    }
  }
  while (0);
  if (ue) {
    sc_options_print_usage (p4est_package_id, SC_LP_ERROR, opt, NULL);
  }

  /* setup appplication model */
  if (!ue && setup_model (g)) {
    P4EST_ASSERT (g->model == NULL);
    ue = usagerr (opt, "model-specific initialization error");
  }

  /* execute application model */
  if (!ue) {
    P4EST_ASSERT (g->model != NULL);
    run_program (g);
  }

  /* cleanup application model */
  if (g->model != NULL) {
    p4est_connectivity_destroy (g->model->conn);
    p4est_gmt_model_destroy (g->model);
  }

  /* deinit main program */
  sc_options_destroy (opt);
  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return ue ? EXIT_FAILURE : EXIT_SUCCESS;
}
