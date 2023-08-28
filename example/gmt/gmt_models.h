#ifndef P4EST_GMT_MODELS_H
#define P4EST_GMT_MODELS_H

#include <p4est_geometry.h>

typedef struct model_latlong_params
{
  int                 latitude[2];
  int                 longitude[2];
  int                 resolution;
  const char         *load_filename;
  const char         *output_prefix;
}
model_latlong_params_t;

typedef struct p4est_gmt_model p4est_gmt_model_t;

struct p4est_gmt_model
{
  int                 M;
  void               *model_data;

  /** Used only to free whatever is stored in model->model_data */
  void                (*destroy_data) (p4est_gmt_model_t * model);
  int                 (*intersect) (int blockno, const double *coord, int m,
                                    void *vmodel);

  const char         *output_prefix;

  p4est_connectivity_t *conn;
  p4est_geometry_t   *model_geom, sgeom;
};

p4est_gmt_model_t  *p4est_gmt_model_synth_new (int synthno);

p4est_gmt_model_t  *p4est_gmt_model_latlong_new (model_latlong_params_t *
                                                 params);

void                p4est_gmt_model_destroy (p4est_gmt_model_t * model);

#endif /* P4EST_GMT_MODELS_H */
