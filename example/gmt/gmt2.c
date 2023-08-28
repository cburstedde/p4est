
#include <p4est_extended.h>
#include <p4est_vtk.h>
#include "gmt_models.h"

typedef struct global
{
  int                 minlev;
  int                 resolution;
  p4est_t            *p4est;
  p4est_gmt_model_t  *model;
}
global_t;

int
main (int argc, char **argv)
{
  sc_MPI_Comm         mpicomm;
  int                 modelno = 0;
  size_t              data_size = 0;
  global_t            sg, *g = &sg;

  mpicomm = sc_MPI_COMM_WORLD;
  memset (g, 0, sizeof (*g));

  /* these options come from the command line */
  g->minlev = 1;
  g->resolution = 1;

  /* parse command line for the choice of model */
  /* here just the unit square */

  if (modelno == 0) {
    model_latlong_params_t ap;
    ap.latitude[0] = -50.;
    ap.latitude[1] = 0.;
    ap.longitude[0] = 0.;
    ap.longitude[1] = 60.;
    ap.resolution = g->resolution;
    ap.load_filename = "africa.gmt.data";
    ap.output_prefix = "africa";

    /* load data (possibly GMT, or file, or synthetic) */
    g->model = p4est_gmt_model_latlong_new (&ap);
  }
  else if (modelno == 1) {
    /* "norway" instead of "africa" */
    /* etc. */
  }

  /* create mesh */
  g->p4est = p4est_new_ext (mpicomm, g->model->conn, 0, g->minlev, 1,
                            data_size, NULL, g);

  /* run mesh refinement based on data */

  /* output refined mesh */

  p4est_vtk_write_file (g->p4est, g->model->model_geom,
                        g->model->output_prefix);

  /* cleanup data */
  p4est_destroy (g->p4est);
  p4est_connectivity_destroy (g->model->conn);

  p4est_gmt_model_destroy (g->model);

  return EXIT_SUCCESS;
}
