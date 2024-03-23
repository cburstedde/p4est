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

/** \file gmt2.c
 *
 * Search based refinement. There are 3 models: synthetic, sphere and latlong.
 *
 * Usage of the sphere model follows the following pipeline.
 *  -# Prepare a csv file of geodesics, following the convention described
 *     in \ref sphere_preprocessing.c . Note that world-map datasets frequently
 *     come in .shp files. The raw line geodesic data can be extracted from
 *     these easily using the python geopandas library. However, there is
 *     currently only limited library support for .shp files in C.
 *  -# Run the preprocessing script as described in \ref sphere_preprocessing.c
 *  -# Run the sphere model:
 *     p4est_gmt --sphere -r <max refinement> -F <output of preprocessing>
*/

#include <p4est_extended.h>
#include <p4est_search.h>
#include <p4est_vtk.h>
#include <p4est_communication.h>
#include <sc_options.h>
#include "gmt_models.h"
#include "gmt_global.h"

static const double irootlen = 1. / (double) P4EST_ROOT_LEN;

static int
setup_model (global_t * g)
{
  /* this function populates the model on successful initialization */
  P4EST_ASSERT (g->model == NULL);

  /* initialize model depending on command line options */
  if (g->synthetic >= 0) {
    switch (g->synthetic) {
    case 0:
      g->model = p4est_gmt_model_synth_new (g->synthetic, g->resolution);
      break;
    default:
      P4EST_GLOBAL_LERROR ("Synthetic model number exceeded\n");
    }
  }
  else if (g->latlongno >= 0) {
    p4est_gmt_model_latlong_params_t ap;

    switch (g->latlongno) {
    case 0:
      /* NOTE: We can make even this a command line input... */
      ap.latitude[0] = -50.;
      ap.latitude[1] = 0.;
      ap.longitude[0] = 0.;
      ap.longitude[1] = 60.;
      ap.resolution = g->resolution;
      ap.load_filename = g->input_filename;
      ap.output_prefix = g->output_prefix;

      /* load data (possibly GMT, or file, or synthetic) */
      g->model = p4est_gmt_model_latlong_new (&ap);
      break;
    default:
      P4EST_GLOBAL_LERROR ("Latitute-longitude model number exceeded\n");
    }
  }
  else if (g->sphere) {
    g->model =
      p4est_gmt_model_sphere_new (g->resolution, g->input_filename,
                                  g->output_prefix, g->distributed,
                                  g->mpicomm);
  }

  /* on successful initalization the global model is set */
  return g->model == NULL ? -1 : 0;
}

static void
quad_init (p4est_t * p4est,
           p4est_topidx_t which_tree, p4est_quadrant_t * quadrant)
{
  quadrant->p.user_int = 0;
}

static int
quad_refine (p4est_t * p4est,
             p4est_topidx_t which_tree, p4est_quadrant_t * quadrant)
{
  return quadrant->p.user_int;
}

static int
quad_point (p4est_t * p4est,
            p4est_topidx_t which_tree, p4est_quadrant_t * quadrant,
            p4est_locidx_t local_num, void *point_index)
{
  int                 result;
  p4est_gmt_model_t  *model;
  global_t           *g = (global_t *) p4est->user_pointer;
  p4est_locidx_t      pi; 

  /* sanity checks */
  P4EST_ASSERT (g != NULL);
  P4EST_ASSERT (g->p4est == p4est);
  model = g->model;
  P4EST_ASSERT (model != NULL);
  P4EST_ASSERT (model->intersect != NULL);

  /* retrieve point index */
  P4EST_ASSERT (point_index != NULL);
  pi = *(p4est_locidx_t*)point_index;
  P4EST_ASSERT (pi < model->M);

  /* execute intersection test */
  if ((result = g->model->intersect(which_tree, quadrant, 
        sc_array_index(model->c->points, pi) , g)) &&
        local_num >= 0 && quadrant->level < g->maxlevel) 
  {
    /* set refinement indicator for a leaf quadrant */
    quadrant->p.user_int = 1;
  }
  return result;
}

int
run_program (global_t * g)
{
  int                 refiter;
  p4est_locidx_t      il;
  char                filename[BUFSIZ];
  sc_array_t         *points = NULL;
  p4est_gloidx_t      gnq_before;
  const size_t        quad_data_size = 0;
  int                 err = 0;

  /* create mesh */
  P4EST_GLOBAL_PRODUCTION ("Create initial mesh\n");
  g->p4est = p4est_new_ext (g->mpicomm, g->model->conn, 0, g->minlevel, 1,
                            quad_data_size, quad_init, g);
  /* in non-distributed mode set up (permanent) search objects */
  if (!g->distributed) {
    P4EST_GLOBAL_PRODUCTIONF ("Setting up %lld search objects\n",
                              (long long) g->model->M);
    points = sc_array_new_count (sizeof (p4est_locidx_t), g->model->M);
    for (il = 0; il < g->model->M; ++il) {
      *(p4est_locidx_t *) sc_array_index (points, il) = il;
    }
  }
  for (refiter = 0;; ++refiter) {
    P4EST_GLOBAL_PRODUCTIONF ("Into refinement iteration %d\n", refiter);
    snprintf (filename, BUFSIZ, "p4est_gmt_%s_%02d",
              g->model->output_prefix, refiter);
    P4EST_ASSERT (g->model != NULL);
    P4EST_ASSERT (g->model->model_geom != NULL);
    p4est_vtk_write_file (g->p4est, g->model->model_geom, filename);
    gnq_before = g->p4est->global_num_quadrants;

    if (g->distributed) {
      /* communicate points */
      err = p4est_transfer_search(g->p4est, g->model->c, 
                    g->model->intersect, 0);
      g->model->M = g->model->c->points->elem_count;

      /* break on communication error */
      if (err) {
        break;
      }

      /* set up search objects for this iteration */
      P4EST_PRODUCTIONF ("Setting up %lld search objects\n",
                            (long long) g->model->M);
      points = sc_array_new_count (sizeof (p4est_locidx_t), g->model->M);
      for (il = 0; il < g->model->M; ++il) {
        *(p4est_locidx_t *) sc_array_index (points, il) = il;
      }
    }

    P4EST_GLOBAL_PRODUCTION ("Run object search\n");
    p4est_search_reorder (g->p4est, 1, NULL, NULL, NULL, quad_point, points);

    /* destroy search objects */
    if (g->distributed) {
      sc_array_destroy_null (&points);
    }

    P4EST_GLOBAL_PRODUCTION ("Run mesh refinement\n");
    p4est_refine (g->p4est, 0, quad_refine, quad_init);

    if (g->balance) {
      P4EST_GLOBAL_PRODUCTION ("Run 2:1 mesh balance\n");
      p4est_balance (g->p4est, P4EST_CONNECT_FULL, quad_init);
    }

    if (gnq_before < g->p4est->global_num_quadrants) {
      P4EST_GLOBAL_PRODUCTION ("Run mesh repartition\n");
      p4est_partition (g->p4est, 0, NULL);
    }
    else {
      P4EST_GLOBAL_PRODUCTION ("Done refinement iterations\n");
      break;
    }
  }
  /* cleanup */
  if (points != NULL) {
    sc_array_destroy_null (&points);
  }
  p4est_destroy (g->p4est);

  return err;
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
  sc_options_add_bool (opt, 'b', "balance", &g->balance, 0,
                       "Execute 2:1 balance algorithm");
  sc_options_add_int (opt, 'r', "resolution", &g->resolution, 0,
                      "Level of resolution (model specific)");
  sc_options_add_int (opt, 'S', "synthetic", &g->synthetic, -1,
                      "Choose specific synthetic model");
  sc_options_add_int (opt, 'M', "latlongno", &g->latlongno, -1,
                      "Choose specific latitude-longitude model");
  sc_options_add_bool (opt, 'W', "sphere", &g->sphere, 0, "Use sphere model");
  sc_options_add_string (opt, 'F', "in-filename", &g->input_filename, NULL,
                         "Choose model-specific input file name");
  sc_options_add_string (opt, 'O', "out-prefix", &g->output_prefix, NULL,
                         "Choose prefix for output file(s)");
  sc_options_add_bool (opt, 'd', "distributed", &g->distributed, 0,
                       "Distributed read mode");

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
    if (g->synthetic >= 0 ? (g->latlongno >= 0 || g->sphere)
        : (g->latlongno >= 0 && g->sphere)
      ) {
      ue =
        usagerrf (opt,
                  "set only one of the synthetic, sphere and latlong models");
    }
    if (g->synthetic < 0 && g->latlongno < 0 && !g->sphere) {
      ue =
        usagerrf (opt,
                  "set one of the synthetic, sphere, and latlong models");
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
    ue = run_program (g);
  }

  /* cleanup application model */
  if (g->model != NULL) {
    p4est_gmt_model_destroy (g->model);
  }

  /* deinit main program */
  sc_options_destroy (opt);
  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return ue ? EXIT_FAILURE : EXIT_SUCCESS;
}
