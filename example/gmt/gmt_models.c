
#include "gmt_models.h"

static int
model_latlong_intersect (int blockno, const double *coord, int m,
                         void *context)
{
  p4est_gmt_model_t  *model = (p4est_gmt_model_t *) context;

  return 0;
}

p4est_gmt_model_t  *
p4est_gmt_model_latlong_new (model_latlong_params_t * params)
{
  p4est_gmt_model_t  *model = P4EST_ALLOC_ZERO (p4est_gmt_model_t, 1);

  model->conn = p4est_connectivity_new_unitsquare ();

  /* load model properties */

  model->model_data = NULL;     /* <- load something from params->load_filename */

  model->intersect = model_latlong_intersect;
  model->destroy_data = NULL;

  model->output_prefix = params->output_prefix;

  model->model_geom = NULL;     /* <- please populate with a function */

  return model;
}

void
p4est_gmt_model_destroy (p4est_gmt_model_t * model)
{
  if (model->destroy_data != NULL) {
    model->destroy_data (model);
  }
  P4EST_FREE (model);
}
