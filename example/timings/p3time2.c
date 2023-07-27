/*
  This file is part of p4est, version 3.
  p4est is a C library to manage a collection (a forest) of multiple
  connected adaptive quadtrees or octrees in parallel.

  Copyright (C) 2019 individual authors
  Originally written by Carsten Burstedde, Lucas C. Wilcox, and Tobin Isaac

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

#ifndef P4_TO_P8
#include <p4est3_quadrant_mort2d.h>
#include <p4est3_quadrant_yx.h>
#include <p4est3_p4est.h>
#else
#include <p4est3_quadrant_mort3d.h>
#include <p4est3_quadrant_zyx.h>
#include <p4est3_p8est.h>
#endif

#include <sc_statistics.h>
#include <sc_flops.h>

#ifdef P4_TO_P8
static const int neighbor_order3d[8] = {2, 1, 3, 5, 0, 4, 0, 1};
#else
static const int neighbor_order2d[4] = {0, 1, 3, 2};
#endif

typedef enum time_function
{
  CHILD,
  MORTON_F,
  PARENT,
  SIBLING,
  FNEIGBOR,
  TBOUND,
  LAST_FUNC
}
time_function_t;

typedef enum time_qtype
{
  STANDARD,
  AVX,
  MORTON,
  LAST_QTYPE
}
time_qtype_t;

typedef struct p4est3_time
{
  int                 mpirank;
  int                 max_level;
  int                 mpisize;
  int                 ninit_quads;
  int                 max_cpu_ref;
  char               *heading;
  time_function_t     func;
  time_qtype_t        qtype;
  sc3_allocator_t    *alloc;
  sc3_array_t        *qarr;
  const p4est3_quadrant_vtable_t *qvt;
}
p4est3_time_t;

static sc3_error_t *
time_array_new (sc3_allocator_t * alloc, size_t esize,
                p4est3_locidx n, sc3_array_t ** arr)
{
  SC3E_RETVAL (arr, NULL);
  SC3A_IS (sc3_allocator_is_setup, alloc);
  SC3A_CHECK (n >= 0);

  SC3E (sc3_array_new (alloc, arr));
  SC3E (sc3_array_set_elem_size (*arr, esize));
  SC3E (sc3_array_set_elem_count (*arr, n));
  SC3E (sc3_array_setup (*arr));

  return NULL;
}

char               *
set_heading (int argc, char **argv)
{
  char               *heading;
  if (argc > 2) {
    heading = (char *) malloc (strlen (argv[1]) + 1 + strlen (argv[2]) + 1);
    strcpy (heading, argv[1]);
    strcat (heading, " ");
    strcat (heading, argv[2]);
  }
  else {
    heading = (char *) malloc (2);
    strcpy (heading, "-");
  }
  return heading;
}

static sc3_error_t *
interpret_command_line (const int argc, char **argv, p4est3_time_t * t)
{
  int                 scale;
  SC3E (sc3_MPI_Comm_size (SC3_MPI_COMM_WORLD, &t->mpisize));
  t->func = CHILD;
  t->qtype = STANDARD;
  t->max_cpu_ref = t->mpisize;
  t->max_level = 1;
  t->ninit_quads = t->max_cpu_ref / t->mpisize;
  if (argc < 5) {
    printf ("Execution without some parameters. "
            "Default parameters are applied.\n"
            "Parameter's format: "
            "<FUNCTION> <QUADRANT TYPE> <MAX CPU REF> <SCALE> <MAX LEVEL>\n");
  }
  if (argc >= 2) {
    t->func = LAST_FUNC;
    if (strcmp (argv[1], "CHILD") == 0) {
      t->func = CHILD;
    }
    else if (strcmp (argv[1], "PARENT") == 0) {
      t->func = PARENT;
    }
    else if (strcmp (argv[1], "SIBLING") == 0) {
      t->func = SIBLING;
    }
    else if (strcmp (argv[1], "FNEIGBOR") == 0) {
      t->func = FNEIGBOR;
    }
    else if (strcmp (argv[1], "TBOUND") == 0) {
      t->func = TBOUND;
    }
    else if (strcmp (argv[1], "MORTON_F") == 0) {
      t->func = MORTON_F;
    }
    SC3E_DEMAND (t->func != LAST_FUNC, "Wrong name of the test function");
  }
  if (argc >= 3) {
    t->qtype = LAST_QTYPE;
    if (strcmp (argv[2], "STANDARD") == 0) {
      t->qtype = STANDARD;
    }
    else if (strcmp (argv[2], "AVX") == 0) {
      t->qtype = AVX;
    }
    else if (strcmp (argv[2], "MORTON") == 0) {
      t->qtype = MORTON;
    }
    SC3E_DEMAND (t->qtype != LAST_QTYPE, "Wrong name of the quadrant type");
  }
  if (argc >= 4) {
    t->max_cpu_ref = atoi (argv[3]);
    SC3E_DEMAND (t->max_cpu_ref != 0,
                 "MAX CPU REF is interpreted wrong or == 0");
    t->ninit_quads = t->max_cpu_ref / t->mpisize;
  }
  if (argc >= 5) {
    scale = atoi (argv[4]);
    SC3E_DEMAND (scale != 0, "SCALE is interpreted wrong or == 0");
    t->max_cpu_ref *= scale;
    t->ninit_quads = t->max_cpu_ref / t->mpisize;
  }
  if (argc >= 6) {
    t->max_level = atoi (argv[5]);
    SC3E_DEMAND (t->max_level != 0, "MAX LEVEL is interpreted wrong or == 0");
  }

  t->heading = set_heading (argc, argv);
  return NULL;
}

static sc3_error_t *
time_fill_test_array (sc3_array_t * a, p4est3_quadrant_vtable_t * qvt)
{
  int                 i;
  p4est3_locidx       quad, child_place, n_quads;
  void               *cp, *q;

  SC3E (sc3_array_get_elem_count (a, &n_quads));
  for (quad = 0, child_place = 1; child_place < n_quads;
       ++quad, child_place += P4EST_CHILDREN) {
    for (i = 0; i < P4EST_CHILDREN; ++i) {
      SC3E (sc3_array_index (a, quad, &q));
      SC3E (sc3_array_index (a, child_place + i, &cp));
      SC3E (p4est3_quadrant_child (qvt, q, i, cp));
    }
  }
  return NULL;
}

static sc3_error_t *
time_prepare (p4est3_time_t * t, int *retval)
{
  void               *p;
  p4est3_locidx       n_quads = 0;

  SC3E_RETVAL (retval, -1);
  SC3A_CHECK (t != NULL);

  switch (t->qtype) {
  case STANDARD:
    /* the standard p4est2 virtual table always exists */
    SC3E (p4est3_quadrant_vtable_p4est (&t->qvt));
    break;

  case AVX:
    /* the AVX virtual table can only be set with hardware support */
    SC3E (p4est3_quadrant_yx_vtable (&t->qvt));
    break;

  case MORTON:
    /* the morton virtual table always exists */
    SC3E (p4est3_quadrant_mort2d_vtable (&t->qvt));
    break;
  default:
    SC3E_UNREACH ("Wrong quadrant type is chosen");
    break;
  }

  SC3E_DEMAND (t->qvt != NULL, "AVX is not supported by hardware or"
               "p4est is not build neither in 2D nor 3D");

  /* create a toplevel allocator */
  SC3E (sc3_allocator_new (sc3_allocator_nocount (), &t->alloc));
  SC3E (sc3_allocator_setup (t->alloc));

  /* allocate array of test quadrants */
  for (int i = 0; i <= t->max_level; ++i) {
    n_quads += (1 << (i * P4EST_DIM));
  }
  SC3E (p4est3_quadrant_array_new (t->alloc, t->qvt, n_quads, &t->qarr));

  /* initialize first element */
  SC3E (sc3_array_index (t->qarr, 0, &p));
  SC3E (p4est3_quadrant_root (t->qvt, p));

  /* fill in test arrays */
  SC3E (time_fill_test_array (t->qarr, t->qvt));

  /* clean and successful return */
  *retval = 0;
  return NULL;
}

static sc3_error_t *
test_morton (const int ninit_quads, sc3_array_t * a,
             p4est3_quadrant_vtable_t * qvt, int max_level)
{
  void *tmp;
  int i;
  p4est3_locidx nquad2level, q;

  SC3E (sc3_array_index (a, 0, &tmp));
  for (i = 0; i < ninit_quads; ++i) {
    for (int i = 0; i <= max_level; ++i) {
      nquad2level = (1 << (i * P4EST_DIM));
      for (q = 0; q < nquad2level; ++q) {
        SC3E (p4est3_quadrant_morton (qvt, i, q, tmp));
      }
    }
  }
  return NULL;
}

static sc3_error_t *
test_child (const int ninit_quads, sc3_array_t * a,
            p4est3_quadrant_vtable_t * qvt)
{
  p4est3_locidx       quad, n_quads;
  void               *q, *tmp;
  int i, c;

  SC3E (sc3_array_get_elem_count (a, &n_quads));
  SC3E (sc3_array_index (a, 0, &tmp));
  for (i = 0; i < ninit_quads; ++i) {
    for (quad = 1; quad < n_quads; ++quad) {
      SC3E (sc3_array_index (a, quad, &q));
      for (c = 0; c < P4EST_CHILDREN; ++c) {
        SC3E (p4est3_quadrant_child (qvt, q, c, tmp));
      }
    }
  }
  return NULL;
}

static sc3_error_t *
test_parent (const int ninit_quads, sc3_array_t * a,
             p4est3_quadrant_vtable_t * qvt)
{
  p4est3_locidx       quad, n_quads;
  void               *q, *tmp;
  int i;

  SC3E (sc3_array_get_elem_count (a, &n_quads));
  SC3E (sc3_array_index (a, 0, &tmp));
  for (i = 0; i < ninit_quads * P4EST_CHILDREN; ++i) {
    for (quad = 1; quad < n_quads; ++quad) {
      SC3E (sc3_array_index (a, quad, &q));
      SC3E (p4est3_quadrant_parent (qvt, q, tmp));
    }
  }
  return NULL;
}

static sc3_error_t *
test_sibling (const int ninit_quads, sc3_array_t * a,
              p4est3_quadrant_vtable_t * qvt)
{
  p4est3_locidx       quad, n_quads;
  void               *q, *tmp;
  int i, c;

  SC3E (sc3_array_get_elem_count (a, &n_quads));
  SC3E (sc3_array_index (a, 0, &tmp));
  for (i = 0; i < ninit_quads; ++i) {
    for (quad = 1; quad < n_quads; ++quad) {
      SC3E (sc3_array_index (a, quad, &q));
      for (c = 0; c < P4EST_CHILDREN; ++c) {
        SC3E (p4est3_quadrant_sibling (qvt, q, c, tmp));
      }
    }
  }
  return NULL;
}

static sc3_error_t *
test_face_neighbor (const int ninit_nquads, sc3_array_t * a,
                    p4est3_quadrant_vtable_t * qvt)
{
  p4est3_locidx quad, n_quads;
  void *q, *tmp;
  int i, c;
  const int *order
#ifdef P4_TO_P8
  = neighbor_order3d
#else
  = neighbor_order2d
#endif
  ;

  SC3E (sc3_array_get_elem_count (a, &n_quads));
  SC3E (sc3_array_index (a, 0, &tmp));
  for (i = 0; i < ninit_nquads; ++i) {
    for (quad = 1; quad < n_quads; ++quad) {
      SC3E (sc3_array_index (a, quad, &q));
      for (c = 0; c < P4EST_CHILDREN; ++c) {
        SC3E (p4est3_quadrant_face_neighbor
              (qvt, q, order[c], tmp));
      }
    }
  }
  return NULL;
}

static sc3_error_t *
test_tree_boundaries (const int ninit_nquads, sc3_array_t * a,
                      p4est3_quadrant_vtable_t * qvt, sc3_array_t * tmp)
{
  p4est3_locidx quad, n_quads;
  void *q;
  int i;

  SC3E (sc3_array_get_elem_count (a, &n_quads));
  for (i = 0; i < ninit_nquads * P4EST_CHILDREN; ++i) {
    for (quad = 1; quad < n_quads; ++quad) {
      SC3E (sc3_array_index (a, quad, &q));
      SC3E (p4est3_quadrant_tree_boundaries (qvt, q, tmp));
    }
  }
  return NULL;
}

static sc3_error_t *
timeavx2_measure (p4est3_time_t * t)
{
  sc_flopinfo_t       fi, snapshot;
  sc_statinfo_t       stats;
  sc3_array_t *tmp;
  switch (t->func) {
  case CHILD:
    sc_flops_snap (&fi, &snapshot);
    SC3E (test_child (t->ninit_quads, t->qarr, t->qvt));
    sc_flops_shot (&fi, &snapshot);
    break;

  case MORTON_F:
    sc_flops_snap (&fi, &snapshot);
    SC3E (test_morton (t->ninit_quads, t->qarr, t->qvt, t->max_level));
    sc_flops_shot (&fi, &snapshot);
    break;

  case PARENT:
    sc_flops_snap (&fi, &snapshot);
    SC3E (test_parent (t->ninit_quads, t->qarr, t->qvt));
    sc_flops_shot (&fi, &snapshot);
    break;

  case SIBLING:
    sc_flops_snap (&fi, &snapshot);
    SC3E (test_sibling (t->ninit_quads, t->qarr, t->qvt));
    sc_flops_shot (&fi, &snapshot);
    break;

  case FNEIGBOR:
    sc_flops_snap (&fi, &snapshot);
    SC3E (test_face_neighbor (t->ninit_quads, t->qarr, t->qvt));
    sc_flops_shot (&fi, &snapshot);
    break;

  case TBOUND:
    SC3E (time_array_new (t->alloc, sizeof (int), P4EST_DIM, &tmp));
    sc_flops_snap (&fi, &snapshot);
    SC3E (test_tree_boundaries (t->ninit_quads, t->qarr, t->qvt, tmp));
    sc_flops_shot (&fi, &snapshot);
    SC3E (sc3_array_unref (&tmp));
    break;

  default:
    SC3E_UNREACH ("Wrong function to test");
    break;
  }
  sc_stats_set1 (&stats, snapshot.iwtime, t->heading);
  sc_stats_compute (SC3_MPI_COMM_WORLD, 1, &stats);
  if (t->mpirank == 0) {
    sc_stats_print (p4est_package_id, SC_LP_ESSENTIAL, 1, &stats, 1, 1);
  }
  return NULL;
}

static sc3_error_t *
timeavx2_cleanup (p4est3_time_t * t)
{
  SC3A_CHECK (t != NULL);

  free (t->heading);
  SC3E (sc3_array_destroy (&t->qarr));
  SC3E (sc3_allocator_destroy (&t->alloc));
  return NULL;
}

int
main (int argc, char **argv)
{
  int                 retval;
  p4est3_time_t       st, *t = &st;

  /* MPI_Init comes first in a program.  We abort should this go wrong. */
  SC3X (sc3_MPI_Init (&argc, &argv));
  SC3X (sc3_MPI_Comm_rank (SC3_MPI_COMM_WORLD, &t->mpirank));

  SC3X (interpret_command_line (argc, argv, t));

  /* choose virtual tables and initialize resources */
  SC3X (time_prepare (t, &retval));

  if (!retval) {
    SC3X (timeavx2_measure (t));
    SC3X (timeavx2_cleanup (t));
  }
  /* MPI_Finalize comes last in a program.  We abort should this go wrong. */
  SC3X (sc3_MPI_Finalize ());
  return 0;
}
