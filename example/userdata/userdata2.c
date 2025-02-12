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
 * This example program demonstrates how to manage application data.
 *
 *   p4est_userdata <OPTIONS> [<configuration> [<level>]]
 *
 * The following options are recognized:
 *   --help          Display a usage and help message and exit successfully.
 *   --level         The level may alternatively be specified as an option.
 *                   The second command line argument takes precedence.
 *
 * Invalid options or arguments result in an error message and exit status.
 */
#ifndef P4_TO_P8
static const char  *p4est_userdata_usage =
  "<configuration> is the first optional argument.\n"
  "  The following values are legal (default is \"unit\"):\n"
  "  o unit          Refinement on the unit square.\n"
  "  o periodic      Unit square with all-periodic boundary conditions.\n"
  "  o brick         Refinement on a 2x3 rectangle of quadtrees.\n"
  "  o disk          Refinement on a spherical disk of five trees.\n"
  "  o corner        Refinement on a non-planar hexagon made of three trees.\n"
  "  o moebius       Refinement on a 5-tree Moebius band embedded in 3D.\n"
  "  o icosahedron   Refinement on the icosahedron sphere with geometry.\n"
  "<level> is the second optional argument (default is 4).\n"
  "  It is clamped into the range of [0, P4EST_QMAXLEVEL].\n"
  "  This argument takes precedence over the option of the same name.\n"
  "No more than two non-option arguments may be specified.\n";
#endif

/* This internal header contains application data for both 2D and 3D. */
#include "userdata_global.h"

/* process the command line */
static int
p4est_userdata_process (p4est_userdata_global_t *g)
{
  P4EST_ASSERT (g != NULL && g->options != NULL);
  P4EST_ASSERT (g->configuration != NULL);
  P4EST_ASSERT (g->conn == NULL);
  P4EST_ASSERT (g->geom == NULL);

  /* choose one of the available mesh configurations */
  if (!strcmp (g->configuration, "unit")) {
#ifndef P4_TO_P8
    g->conn = p4est_connectivity_new_unitsquare ();
#else
    g->conn = p8est_connectivity_new_unitcube ();
#endif
  }
#ifndef P4_TO_P8
  else if (!strcmp (g->configuration, "periodic")) {
    g->conn = p4est_connectivity_new_periodic ();
  }
  else if (!strcmp (g->configuration, "brick")) {
    g->conn = p4est_connectivity_new_brick (2, 3, 0, 0);
  }
  else if (!strcmp (g->configuration, "disk")) {
    g->conn = p4est_connectivity_new_disk_nonperiodic ();
    g->geom = p4est_geometry_new_disk2d (g->conn, .4, 1.);
  }
  else if (!strcmp (g->configuration, "corner")) {
    g->conn = p4est_connectivity_new_corner ();
  }
  else if (!strcmp (g->configuration, "moebius")) {
    g->conn = p4est_connectivity_new_moebius ();
  }
  else if (!strcmp (g->configuration, "icosahedron")) {
    g->conn = p4est_connectivity_new_icosahedron ();
    g->geom = p4est_geometry_new_icosahedron (g->conn, 1.);
  }
#else
  else if (!strcmp (g->configuration, "periodic")) {
    g->conn = p8est_connectivity_new_periodic ();
  }
  else if (!strcmp (g->configuration, "brick")) {
    g->conn = p8est_connectivity_new_brick (2, 3, 5, 0, 0, 0);
  }
  else if (!strcmp (g->configuration, "rotcubes")) {
    g->conn = p8est_connectivity_new_rotcubes ();
  }
  else if (!strcmp (g->configuration, "sphere")) {
    g->conn = p8est_connectivity_new_sphere ();
    g->geom = p8est_geometry_new_sphere (g->conn, .3, .6, 1.);
  }
  else if (!strcmp (g->configuration, "shell")) {
    g->conn = p8est_connectivity_new_shell ();
    g->geom = p8est_geometry_new_shell (g->conn, .6, 1.);
  }
  else if (!strcmp (g->configuration, "torus")) {
    g->conn = p8est_connectivity_new_torus (6);
    g->geom = p8est_geometry_new_torus (g->conn, .1, .5, 1.);
  }
#endif

  /* clamp level into legal range */
  if (g->maxlevel < 0) {
    g->maxlevel = 0;
  }
  else if (g->maxlevel > P4EST_QMAXLEVEL) {
    g->maxlevel = P4EST_QMAXLEVEL;
  }

  /* this is the only error condition of this function */
  if (g->conn == NULL) {
    P4EST_GLOBAL_LERROR ("ERROR: Invalid configuration argument\n");
    return -1;
  }

  /* if no geometry is specified, default to vertex information */
  if (g->geom == NULL) {
    g->geom = p4est_geometry_new_connectivity (g->conn);
  }

  /* successful return! */
  return 0;
}

/* process command line options */
static int
p4est_userdata_options (int argc, char **argv, p4est_userdata_global_t *g)
{
  int                 erres;
  int                 firstarg;
  sc_options_t       *o;

  P4EST_ASSERT (argc >= 1);
  P4EST_ASSERT (argv != NULL);
  P4EST_ASSERT (g != NULL);
  P4EST_ASSERT (g->options == NULL);

  /* allocate new options processor */
  o = g->options = sc_options_new (argv[0]);
  sc_options_add_switch (o, 'h', "help", &g->help,
                         "Print help message and exit cleanly");
  sc_options_add_int (o, 'L', "maxlevel", &g->maxlevel, 4,
                      "Maximum refinement level");
  g->configuration = "unit";

  /* error condition of this function */
  erres = 0;

  /* parse command line options and log eventual errors with priority */
  if (!erres && (erres = ((firstarg = sc_options_parse
                           (p4est_get_package_id (), SC_LP_ERROR,
                            o, argc, argv)) < 0))) {
    P4EST_GLOBAL_LERROR ("ERROR: processing options\n");
  }
  P4EST_ASSERT (erres || (1 <= firstarg && firstarg <= argc));

  /* process first non-option command line argument */
  if (!erres && firstarg < argc) {
    g->configuration = argv[firstarg++];
  }

  /* process second non-option command line argument */
  if (!erres && firstarg < argc) {
    g->maxlevel = atoi (argv[firstarg++]);
  }

  /* we do not permit more than two non-option arguments */
  if (!erres && (erres = (firstarg < argc))) {
    P4EST_GLOBAL_LERROR ("ERROR: no more than two arguments are allowed\n");
  }

  /* initialize variables based on command line */
  if (!erres && (erres = p4est_userdata_process (g))) {
    P4EST_GLOBAL_LERROR ("ERROR: processing the command line\n");
  }

  if (g->help || erres) {
    /* print a usage message to explain the command line arguments */
    sc_options_print_usage (p4est_get_package_id (), SC_LP_PRODUCTION, o,
                            p4est_userdata_usage);
  }
  else {
    /* on normal operation print options for posteriority */
    sc_options_print_summary (p4est_get_package_id (), SC_LP_PRODUCTION, o);
  }

  /* this function has processed the command line */
  return erres;
}

/* free allocated application memory */
static void
p4est_userdata_cleanup (p4est_userdata_global_t *g)
{
  P4EST_ASSERT (g != NULL);
  P4EST_ASSERT (g->options != NULL);

  /* this data may or may not have been initialized */
  if (g->geom != NULL) {
    P4EST_ASSERT (g->conn != NULL);
    p4est_geometry_destroy (g->geom);
  }
  if (g->conn != NULL) {
    p4est_connectivity_destroy (g->conn);
  }

  /* options are always defined at this point */
  sc_options_destroy (g->options);
}

/* the main function of the program */
int
main (int argc, char **argv)
{
  int                 erres;
  int                 mpiret;
  p4est_userdata_global_t sglobal, *global = &sglobal;

  /* initialize MPI subsystem */
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  /* initialize global application data context */
  memset (global, 0, sizeof (*global));
  global->mpicomm = sc_MPI_COMM_WORLD;

  /*
   * The options 1, 1 catch signals and set an abort handler.
   * This is NOT what you want if (a) you are using p4est purely as a
   * library or (b) your main programs refer to a different framework for
   * such lowlevel functionality.  In these cases please use 0, 0.
   *
   * The log level SC_LP_APPLICATION is rather conservative.  If you'd
   * prefer total silence under normal operating conditions you may use
   * SC_LP_ERROR.  p4est does not trigger the error priority by itself, but
   * it can be used by the application developer, for example to issue log
   * messages on any usage or file I/O errors detected.
   */
  sc_init (global->mpicomm, 1, 1, NULL, SC_LP_APPLICATION);

  /*
   * The setting SC_LP_APPLICATION will log levels from SC_LP_PRODUCTION
   * upwards.  Thus, if your program should print performance metrics, for
   * example, that can be accomplished (for information available on rank 0)
   * using the logging macros P4EST_PRODUCTION and P4EST_PRODUCTIONF.
   */
  p4est_init (NULL, SC_LP_APPLICATION);

  /* we identify an error status of the program with a nonzero value */
  erres = 0;

  /* process command line options */
  if (!erres && (erres = p4est_userdata_options (argc, argv, global)))
  {
    P4EST_GLOBAL_LERROR ("ERROR: Usage/options\n");
  }

  /*
   * Run actual demo (except when a help message has been requested).
   * We have moved the code for this function into a separate file.
   * The reason is that the present file is an excellent template
   * for your own p4est application.  Just copy it and hack away.
   */
  if (!erres && !global->help && (erres = p4est_userdata_run (global)))
  {
    P4EST_GLOBAL_LERROR ("ERROR: running the program\n");
  }

  /* free all data regardless of the error condition */
  p4est_userdata_cleanup (global);

  /* check memory balance and clean up internal registrations */
  sc_finalize ();

  /* release the MPI subsytem */
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  /* return failure or success to the calling shell */
  return erres ? EXIT_FAILURE : EXIT_SUCCESS;
}
