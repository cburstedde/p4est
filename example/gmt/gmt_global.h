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

typedef struct global
{
  int                 minlevel;
  int                 maxlevel;
  int                 balance;
  int                 resolution;
  int                 synthetic;
  int                 latlongno;
  int                 sphere;   /* globe sphere model */
  int                 distributed; /* distributed file read */
  const char         *input_filename;
  const char         *output_prefix;
  sc_MPI_Comm         mpicomm;
  p4est_t            *p4est;
  p4est_gmt_model_t  *model;
}
global_t;