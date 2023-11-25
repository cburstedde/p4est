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

#ifndef P4EST_GMT_MODELS_H
#define P4EST_GMT_MODELS_H

#include <p4est_geometry.h>
#include <p4est.h>

/** Used to free private model data. */
typedef void        (*p4est_gmt_destroy_data_t) (void *vmodel_data);

/** Check intersection of a quadrant with an object. */
typedef int         (*p4est_gmt_intersect_t) (p4est_topidx_t which_tree,
                                              const double coord[4],
                                              size_t m, void *vmodel);

/** Find the processes that might own quadrants containing the point. */
typedef sc_array_t* (*p4est_gmt_owners_fn_t) (void *point, p4est_t * p4est);

/** Select one owner process to be responsible for propagating a point.
 * 
 * \param[in] point  The point to be propagated
 * \param[in] owners The owners of the point
 */
typedef int (*p4est_gmt_resp_fn_t) (void *point, sc_array_t* owners); 

/** General, application specific model data */
typedef struct p4est_gmt_model
{
  size_t              M;
  const char         *output_prefix;
  p4est_connectivity_t *conn;
  p4est_geometry_t   *model_geom;
  void               *model_data;
  size_t              point_size;

  /** When not NULL, free whatever is stored in model->model_data. */
  p4est_gmt_destroy_data_t destroy_data;

  /** Intersect a given rectangle with a model object. */
  p4est_gmt_intersect_t intersect;

  /** Private geometry data. */
  p4est_geometry_t    sgeom;

  /** Determine owners of point. Can be NULL if not running distributed  */
  p4est_gmt_owners_fn_t owners_fn;
  /** Determine who is responsible for propagating point */
  p4est_gmt_resp_fn_t resp_fn;
}
p4est_gmt_model_t;

/** Create a specific synthetic model */
p4est_gmt_model_t  *p4est_gmt_model_synth_new (int synthno, int resolution);

/** Parameter type for latitude-longitude model */
typedef struct p4est_gmt_model_latlong_params
{
  int                 latitude[2];
  int                 longitude[2];
  int                 resolution;
  const char         *load_filename;
  const char         *output_prefix;
}
p4est_gmt_model_latlong_params_t;

/** Create a specific latlong model */
p4est_gmt_model_t  *p4est_gmt_model_latlong_new
  (p4est_gmt_model_latlong_params_t * params);

/** Represents a segment of a geodesic in the sphere model. Includes some
 * additional data for quicker computations with this segment.
 * 
 * Segments are restricted to lying on a single face of the cube-sphere.
 * A segment is represented by its endpoints, given in tree-local
 * reference coordinates.
 * 
 * bb1 and bb2 are the lower-left coordinates of the first and last atom
 * of the smallest quadrant containing the segment. These are used to
 * determine processes whose domain possibly intersects the segment.
 */
typedef struct p4est_gmt_sphere_geodesic_seg
{
  int which_tree;
  double p1x, p1y, p2x, p2y; /* Geodesic endpoints */
  double bb1x, bb1y, bb2x, bb2y; /* Bounding quadrant start and end atoms */
} p4est_gmt_sphere_geodesic_seg_t;

/** Data used by the sphere model */
typedef struct p4est_gmt_model_sphere
{
  int resolution;
  size_t num_geodesics;
  p4est_gmt_sphere_geodesic_seg_t *points;
} p4est_gmt_model_sphere_t;

/** Create a specific sphere model.
 * 
 * The sphere model refines a spherical mesh based on geodesics. More specifically,
 * squares in the mesh are recursively refined as long as they intersect a geodesic and
 * have refinement level less than the desired resolution. An example application is
 * refining a map of the globe based on coastlines.
 * 
 * \warning Before running this function the preprocessing script
 * \ref sphere_preprocessing.c must be called.
 *
 * \param[in] resolution maximum refinement level
 */
p4est_gmt_model_t  * p4est_gmt_model_sphere_new (int resolution);

/** Destroy model */
void                p4est_gmt_model_destroy (p4est_gmt_model_t * model);

#endif /* P4EST_GMT_MODELS_H */
