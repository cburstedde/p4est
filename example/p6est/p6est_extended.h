/*
  This file is part of p4est.
  p4est is a C library to manage a collection (a forest) of multiple
  connected adaptive quadtrees or octrees in parallel.

  Copyright (C) 2013 The University of Texas System
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

#ifndef P6EST_EXTENDED_H
#define P6EST_EXTENDED_H

p6est_t            *p6est_new_ext (MPI_Comm mpicomm,
                                   p6est_connectivity_t * connectivity,
                                   p4est_locidx_t min_quadrants,
                                   int min_level, int min_zlevel,
                                   int fill_uniform, size_t data_size,
                                   p6est_init_t init_fn, void *user_pointer);
void                p6est_save_ext (const char *filename, p6est_t * p6est,
                                    int save_data, int save_partition);
p6est_t            *p6est_load_ext (const char *filename, MPI_Comm mpicomm,
                                    size_t data_size, int load_data,
                                    int autopartition, int broadcasthead,
                                    void *user_pointer,
                                    p6est_connectivity_t ** connectivity);
void                p6est_refine_columns_ext (p6est_t * p6est,
                                              int refine_recursive,
                                              int allowed_level,
                                              p6est_refine_column_t refine_fn,
                                              p6est_init_t init_fn,
                                              p6est_replace_t replace_fn);
void                p6est_refine_layers_ext (p6est_t * p6est,
                                             int refine_recursive,
                                             int allowed_level,
                                             p6est_refine_layer_t refine_fn,
                                             p6est_init_t init_fn,
                                             p6est_replace_t replace_fn);
void                p6est_coarsen_columns_ext (p6est_t * p6est,
                                               int coarsen_recursive,
                                               int callback_orphans,
                                               p6est_coarsen_column_t
                                               coarsen_fn,
                                               p6est_init_t init_fn,
                                               p6est_replace_t replace_fn);
void                p6est_coarsen_layers_ext (p6est_t * p6est,
                                              int coarsen_recursive,
                                              int callback_orphans,
                                              p6est_coarsen_layer_t
                                              coarsen_fn,
                                              p6est_init_t init_fn,
                                              p6est_replace_t replace_fn);
p4est_gloidx_t      p6est_partition_ext (p6est_t * p6est,
                                         int partition_for_coarsening,
                                         p6est_weight_t weight_fn);
void                p6est_balance_ext (p6est_t * p6est,
                                       p8est_connect_type_t btype,
                                       int max_diff, int min_diff,
                                       p6est_init_t init_fn,
                                       p6est_replace_t replace_fn);

#endif
