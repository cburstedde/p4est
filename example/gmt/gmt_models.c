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

#include "gmt_models.h"

static int
model_latlong_intersect (int blockno, const double *coord, int m,
                         void *vmodel)
{
  p4est_gmt_model_t  *model = (p4est_gmt_model_t *) vmodel;

  return 0;
}

static void
model_latlong_geom_X (p4est_geometry_t * geom, p4est_topidx_t which_tree,
                      const double abc[3], double xyz[3])
{
  p4est_gmt_model_t  *model = (p4est_gmt_model_t *) geom->user;

#if LATLONG_DATA_HAS_BEEN_PROGRAMMED
  /* put the parameters latitude, longitude into the model data */
  longitude =
    ((typecast into gmt lanlong model data *) model->model_data)->longitude;
  latitude = ...;

  xyz[0] = longitude[0] + (longitude->[1] - longitude[0]) * abc[0];
  xyz[1] = latitude[0] + (latitude->[1] - latitude[0]) * abc[1];
#else
  xyz[0] = abc[0];
  xyz[1] = abc[1];
#endif
  xyz[2] = 0.;
}

p4est_gmt_model_t  *
p4est_gmt_model_latlong_new (model_latlong_params_t * params)
{
  p4est_gmt_model_t  *model = P4EST_ALLOC_ZERO (p4est_gmt_model_t, 1);

  model->conn = p4est_connectivity_new_unitsquare ();

  /* load model properties */

  model->model_data = NULL;     /* <- Load something from params->load_filename,
                                   also deep copy the parameters into it. */

  model->intersect = model_latlong_intersect;
  model->destroy_data = NULL;

  model->output_prefix = params->output_prefix;

  model->sgeom.name = params->output_prefix;
  model->sgeom.user = model;
  model->sgeom.X = model_latlong_geom_X;
  model->sgeom.destroy = NULL;
  model->model_geom = &model->sgeom;

  return model;
}

p4est_gmt_model_t  *
p4est_gmt_model_synth_new (int synthno)
{
  return NULL;
}

void
p4est_gmt_model_destroy (p4est_gmt_model_t * model)
{
  if (model->destroy_data != NULL) {
    model->destroy_data (model);
  }
  P4EST_FREE (model);
}
