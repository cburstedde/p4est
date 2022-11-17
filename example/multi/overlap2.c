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
 * Usage: p4est_overlap
 */

#ifndef P4_TO_P8
#include <p4est_extended.h>
#else
#include <p8est_extended.h>
#endif

typedef struct overlap_prodata
{
  int                 pdummy;
}
overlap_prodata_t;

typedef struct overlap_producer
{
  sc_MPI_Comm         procomm;
  p4est_connectivity_t *proconn;
  p4est_t            *pro4est;
  int                 pminl;
}
overlap_producer_t;

typedef struct overlap_condata
{
  int                 cdummy;
}
overlap_condata_t;

typedef struct overlap_consumer
{
  sc_MPI_Comm         concomm;
  p4est_connectivity_t *conconn;
  p4est_t            *con4est;
  int                 cminl;
}
overlap_consumer_t;

typedef struct overlap_global
{
  sc_MPI_Comm         glocomm;
  overlap_producer_t  pro, *p;
  overlap_consumer_t  con, *c;
}
overlap_global_t;

static void
overlap_apps_init (overlap_global_t * g, sc_MPI_Comm mpicomm)
{
  overlap_producer_t *p = g->p = &g->pro;
  overlap_consumer_t *c = g->c = &g->con;
  int                 mpiret;

  /* initialization of global data */
  g->glocomm = mpicomm;

  /* still hardwired configuration */
  p->pminl = 1;
  c->cminl = 2;

  /* setup producer app */
  mpiret = sc_MPI_Comm_dup (g->glocomm, &p->procomm);
  SC_CHECK_MPI (mpiret);
  p->proconn = p4est_connectivity_new_brick (2, 2
#ifdef P4_TO_P8
                                             , 2
#endif
                                             , 0, 0
#ifdef P4_TO_P8
                                             , 0
#endif
    );
  p->pro4est = p4est_new_ext (p->procomm, p->proconn, 0, p->pminl, 1,
                              sizeof (overlap_prodata_t), NULL, p);

  /* setup consumer app */
  mpiret = sc_MPI_Comm_dup (g->glocomm, &c->concomm);
  SC_CHECK_MPI (mpiret);
  c->conconn = p4est_connectivity_new_brick (3, 2
#ifdef P4_TO_P8
                                             , 1
#endif
                                             , 0, 0
#ifdef P4_TO_P8
                                             , 0
#endif
    );
  c->con4est = p4est_new_ext (c->concomm, c->conconn, 0, c->cminl, 1,
                              sizeof (overlap_condata_t), NULL, c);
}

static void
overlap_apps_reset (overlap_global_t * g)
{
  overlap_producer_t *p = g->p;
  overlap_consumer_t *c = g->c;
  int                mpiret;

  /* destroy producer */
  p4est_destroy (p->pro4est);
  p4est_connectivity_destroy (p->proconn);
  mpiret = sc_MPI_Comm_free (&p->procomm);
  SC_CHECK_MPI (mpiret);

  /* destroy consumer */
  p4est_destroy (c->con4est);
  p4est_connectivity_destroy (c->conconn);
  mpiret = sc_MPI_Comm_free (&c->concomm);
  SC_CHECK_MPI (mpiret);
}

int
main (int argc, char **argv)
{
  int                 mpiret;
  sc_MPI_Comm         mpicomm;
  overlap_global_t global, *g = &global;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  mpicomm = sc_MPI_COMM_WORLD;
  sc_init (mpicomm, 1, 1, NULL, SC_LP_DEFAULT);

  overlap_apps_init (g, mpicomm);

  overlap_apps_reset (g);

  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
