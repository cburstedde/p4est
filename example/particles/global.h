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
 * This is an internal header file.
 * It is included from a c file that may be preceded by p4est_to_p8est.h.
 */

#ifndef PART_GLOBAL_H
#define PART_GLOBAL_H

#ifndef P4_TO_P8
#include <p4est.h>
#else
#include <p8est.h>
#endif /* P4_TO_P8 */

SC_EXTERN_C_BEGIN;

typedef double      (*part_init_density_t) (double x, double y, double z,
                                            void *data);

typedef struct part_global
{
  sc_MPI_Comm         mpicomm;
  int                 mpisize, mpirank;

  int                 minlevel;
  int                 maxlevel;
  int                 bricklev;
  int                 order;
  int                 vtk;
  int                 check;
  double              num_particles;
  double              elem_particles;
  double              deltat;
  double              finaltime;
  const char         *prefix;

  int                 bricklength;
  int                 stage;
  long long           gpnum, gplost;
  double              global_density;
  double              t;
  sc_array_t         *padata;
  sc_array_t         *pfound;
  sc_array_t         *recevs;
  sc_array_t         *sendes;
  sc_array_t         *send_req;
  sc_array_t         *prebuf;
  sc_hash_t          *psend;
  sc_hash_t          *precv;
  sc_mempool_t       *psmem;
  part_init_density_t pidense;
  void               *piddata;

  p4est_connectivity_t *conn;
  p4est_t            *p4est;
}
part_global_t;

SC_EXTERN_C_END;

#endif /* !PART_GLOBAL_H */
