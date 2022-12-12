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
 *
 * Create two forest workflow apps in the same main program.
 * One app requires data from the other, which we search in parallel.
 * In this example, both apps use a duplicate of MPI_COMM_WORLD.
 * Thus, they execute their respective program alternatingly.
 *
 * The two apps are labeled producer and consumer.
 *
 * In 2D, the producer map is extruded in the third dimension by 1.
 * In 3D, producer and consumer maps are fully 3D.
 * All coordinate transforms are affine-linear.
 */

#include <sc_notify.h>
#ifndef P4_TO_P8
#include <p4est_extended.h>
#include <p4est_search.h>
#include <p4est_vtk.h>
#else
#include <p8est_extended.h>
#include <p8est_search.h>
#include <p8est_vtk.h>
#endif

#define P4EST_CON_TOLERANCE SC_1000_EPS
#define P4EST_PRO_TOLERANCE (2 * SC_1000_EPS)

typedef struct overlap_prodata
{
  double              myvalue;
  int                 isset;
}
overlap_prodata_t;

typedef struct overlap_point
{
  int                 rank;
  p4est_locidx_t      lnum;
  double              xyz[3];
  overlap_prodata_t   prodata;
}
overlap_point_t;

typedef struct overlap_send_buf
{
  int                 rank;
  sc_array_t          ops;
}
overlap_send_buf_t;

typedef overlap_send_buf_t overlap_recv_buf_t;

typedef enum comm_tag
{
  COMM_TAG_CONSDATA = P4EST_COMM_TAG_LAST,
  COMM_TAG_PRODATA
}
comm_tag_t;

typedef struct overlap_producer
{
  /* mesh constituents */
  sc_MPI_Comm         procomm;
  p4est_connectivity_t *proconn;
  p4est_t            *pro4est;
  p4est_geometry_t   *progeom, producer_geometry;

  /* parameters */
  int                 pminl;

  /* communication */
  int                 prorank;
  sc_array_t         *recv_buffer;
}
overlap_producer_t;

typedef struct overlap_condata
{
  int                 cdummy;
}
overlap_condata_t;

typedef struct overlap_consumer
{
  /* mesh constituents */
  sc_MPI_Comm         concomm;
  p4est_connectivity_t *conconn;
  p4est_t            *con4est;
  p4est_geometry_t   *congeom, consumer_geometry;

  /* parameters */
  int                 cminl;

  /* work variables */
  p4est_locidx_t      lquad_idx;
  sc_array_t         *query_xyz;

  /* minimal knowledge of the producer's mesh */
  p4est_connectivity_t *producer_conn;
  const p4est_gloidx_t *producer_gfq;
  const p4est_quadrant_t *producer_gfp;
  int                 pronum_procs;
  p4est_topidx_t      pronum_trees;

  /* communication */
  int                 conrank;
  sc_array_t         *send_buffer;
  sc_array_t         *recv_buffer;
}
overlap_consumer_t;

typedef struct overlap_global
{
  sc_MPI_Comm         glocomm;
  int                 rounds;
  overlap_producer_t  pro, *p;
  overlap_consumer_t  con, *c;
}
overlap_global_t;

#define OVERLAP_IROOTLEN (1. / P4EST_ROOT_LEN)

static double
overlap_producer_evaluate (overlap_producer_t *p, double pxyz[3])
{
  double              r[3];

  P4EST_ASSERT (p != NULL);
  P4EST_ASSERT (pxyz != NULL);

  r[0] = (pxyz[0] - .4) / 1.6;
  r[1] = (pxyz[1] + .3) / 1.1;
  r[2] = (pxyz[2] - .2) / 0.5;
  return
    .1 + .9 * exp (-.5 * (SC_SQR (r[0]) + SC_SQR (r[1]) + SC_SQR (r[2])));
}

static void         overlap_producer_invmap
  (p4est_connectivity_t *proconn, p4est_topidx_t which_tree,
   const double xyz[3], double abc[3]);

static void
overlap_producer_map (p4est_geometry_t *geom, p4est_topidx_t which_tree,
                      const double abc[3], double xyz[3])
{
  double              a, co, si, x;
#ifdef P4EST_ENABLE_DEBUG
  double              def[3];
#endif
  double             *vert;
  p4est_topidx_t      vind;
  overlap_producer_t *p;

  P4EST_ASSERT (geom != NULL);
  P4EST_ASSERT (geom->X == overlap_producer_map);

  /* preimage domain is [0, 2]^3 */

  /* access origin vertex of given tree */
  P4EST_ASSERT (geom != NULL);
  P4EST_ASSERT (geom->X == overlap_producer_map);
  p = (overlap_producer_t *) geom->user;
  P4EST_ASSERT (p->progeom == geom);
  P4EST_ASSERT (p->proconn != NULL && p->proconn->vertices != NULL);
  vind = p->proconn->tree_to_vertex[P4EST_CHILDREN * which_tree + 0];
  vert = &p->proconn->vertices[3 * vind + 0];

  /* move brick back towards origin and scale slightly */
  xyz[0] = (vert[0] + abc[0] - .5) * 1.1;
  xyz[1] = (vert[1] + abc[1] - .5) * 1.2;
  xyz[2] = (vert[2] + abc[2] - .5) * 1.3;

  /* rotate 20 degrees around the y axis */
  a = 20. * M_PI / 180.;
  co = cos (a);
  si = sin (a);
  x = xyz[0];
  xyz[0] = co * x - si * xyz[2];
  xyz[2] = si * x + co * xyz[2];

#ifdef P4EST_ENABLE_DEBUG
  /* verify inverse mapping */
  overlap_producer_invmap (p->proconn, which_tree, xyz, def);
  P4EST_ASSERT (fabs (abc[0] - def[0]) < SC_1000_EPS &&
                fabs (abc[1] - def[1]) < SC_1000_EPS &&
                fabs (abc[2] - def[2]) < SC_1000_EPS);
#endif
}

static void
overlap_producer_invmap (p4est_connectivity_t *proconn,
                         p4est_topidx_t which_tree,
                         const double xyz[3], double abc[3])
{
  double              a, co, si;
  double             *vert;
  p4est_topidx_t      vind;

  P4EST_ASSERT (proconn != NULL && proconn->vertices != NULL);
  P4EST_ASSERT (0 <= which_tree && which_tree < proconn->num_trees);
  vind = proconn->tree_to_vertex[P4EST_CHILDREN * which_tree + 0];
  vert = &proconn->vertices[3 * vind + 0];

  /* rotate 20 degrees backwards around the y axis */
  a = -20. * M_PI / 180.;
  co = cos (a);
  si = sin (a);
  abc[0] = co * xyz[0] - si * xyz[2];
  abc[2] = si * xyz[0] + co * xyz[2];
  abc[1] = xyz[1];

  /* scale slightly and move away from origin into brick */
  abc[0] = abc[0] / 1.1 + .5 - vert[0];
  abc[1] = abc[1] / 1.2 + .5 - vert[1];
  abc[2] = abc[2] / 1.3 + .5 - vert[2];
}

static void
overlap_producer_compute (p4est_iter_volume_info_t *info, void *user_data)
{
  p4est_qcoord_t      h2;
  p4est_quadrant_t   *q;
  overlap_prodata_t  *d;
  overlap_producer_t *p;
  double              qxyz[3], phys[3];

  P4EST_ASSERT (info != NULL && info->p4est != NULL);
  P4EST_ASSERT (info->p4est->user_pointer == user_data);

  p = (overlap_producer_t *) info->p4est->user_pointer;
  P4EST_ASSERT (p->pro4est == info->p4est);
  P4EST_ASSERT (info->quad != NULL);
  q = info->quad;
  d = (overlap_prodata_t *) q->p.user_data;
  P4EST_ASSERT (d != NULL);

  /* transform producer quadrant center to physical using map */
  h2 = P4EST_QUADRANT_LEN (q->level) >> 1;
  qxyz[0] = OVERLAP_IROOTLEN * (q->x + h2);
  qxyz[1] = OVERLAP_IROOTLEN * (q->y + h2);
#ifndef P4_TO_P8
  qxyz[2] = 0.;
#else
  qxyz[2] = OVERLAP_IROOTLEN * (q->z + h2);
#endif
  overlap_producer_map (p->progeom, info->treeid, qxyz, phys);

  /* interpolate prescribed field at that point */
  d->myvalue = overlap_producer_evaluate (p, phys);
  d->isset = 1;

  P4EST_LDEBUGF ("Producer input tree %d level %d quad %g %g %g\n",
                 (int) info->treeid, q->level, qxyz[0], qxyz[1], qxyz[2]);
  P4EST_LDEBUGF ("Producer compute %g %g %g\n", phys[0], phys[1], phys[2]);
}

static void
overlap_consumer_map (p4est_geometry_t *geom, p4est_topidx_t which_tree,
                      const double abc[3], double xyz[3])
{
  double              a, co, si, x;
  double             *vert;
  p4est_topidx_t      vind;
  overlap_consumer_t *c;

  /* preimage domain is 3x2x1 with origin in the lower left front */

  /* access origin vertex of given tree */
  P4EST_ASSERT (geom != NULL);
  P4EST_ASSERT (geom->X == overlap_consumer_map);
  c = (overlap_consumer_t *) geom->user;
  P4EST_ASSERT (c->congeom == geom);
  P4EST_ASSERT (c->conconn != NULL && c->conconn->vertices != NULL);
  vind = c->conconn->tree_to_vertex[P4EST_CHILDREN * which_tree + 0];
  vert = &c->conconn->vertices[3 * vind + 0];

  /* center brick around origin and scale */
  xyz[0] = (vert[0] + abc[0] - 1.5) * 1.4;
  xyz[1] = (vert[1] + abc[1] - 1.0) * 1.1;
  xyz[2] = (vert[2] + abc[2] - 0.5) * 0.7;

  /* rotate 30 degrees around the z axis */
  a = 30. * M_PI / 180.;
  co = cos (a);
  si = sin (a);
  x = xyz[0];
  xyz[0] = co * x - si * xyz[1];
  xyz[1] = si * x + co * xyz[1];
}

static void
overlap_consumer_compute (p4est_iter_volume_info_t *info, void *user_data)
{
  p4est_quadrant_t   *q;
  overlap_consumer_t *c;
  overlap_point_t    *op;
  p4est_qcoord_t      h2;
  double              qxyz[3], *phys;

  P4EST_ASSERT (info != NULL && info->p4est != NULL);
  P4EST_ASSERT (info->p4est->user_pointer == user_data);

  /* access quadrant */
  c = (overlap_consumer_t *) info->p4est->user_pointer;
  P4EST_ASSERT (c->con4est == info->p4est);
  P4EST_ASSERT (c->lquad_idx >= 0 &&
                c->lquad_idx < c->con4est->local_num_quadrants);
  P4EST_ASSERT (info->quad != NULL);
  q = info->quad;

  /* transform consumer quadrant center to physical using map */
  h2 = P4EST_QUADRANT_LEN (q->level) >> 1;
  qxyz[0] = OVERLAP_IROOTLEN * (q->x + h2);
  qxyz[1] = OVERLAP_IROOTLEN * (q->y + h2);
#ifndef P4_TO_P8
  qxyz[2] = 0.;
#else
  qxyz[2] = OVERLAP_IROOTLEN * (q->z + h2);
#endif
  phys = (op = (overlap_point_t *)
          sc_array_index (c->query_xyz, (size_t) c->lquad_idx))->xyz;
  overlap_consumer_map (c->congeom, info->treeid, qxyz, phys);
  op->lnum = c->lquad_idx++;
  op->rank = -1;
  op->prodata.myvalue = 0.;
  op->prodata.isset = 0;

  P4EST_LDEBUGF ("Consumer input tree %d level %d quad %g %g %g\n",
                 (int) info->treeid, q->level, qxyz[0], qxyz[1], qxyz[2]);
  P4EST_LDEBUGF ("Consumer point %ld compute %g %g %g\n",
                 (long) op->lnum, phys[0], phys[1], phys[2]);

  /* optimize: ignore this point if not intersecting producer domain */
}

static void
overlap_apps_init (overlap_global_t *g, sc_MPI_Comm mpicomm)
{
  overlap_producer_t *p = g->p = &g->pro;
  overlap_consumer_t *c = g->c = &g->con;
  int                 mpiret;

  /* initialization of global data */
  g->glocomm = mpicomm;
  g->rounds = 0;

  /* still hardwired configuration */
  p->pminl = 0;
  c->cminl = 0;

  /***************************** PRODUCER ****************************/

  P4EST_GLOBAL_PRODUCTION ("OVERLAP: init producer\n");

  /* setup producer geometry */
  p->progeom = &p->producer_geometry;
  p->progeom->name = "producer";
  p->progeom->user = p;
  p->progeom->X = overlap_producer_map;
  p->progeom->destroy = (p4est_geometry_destroy_t) 0;

  /* setup producer app with communicator and mesh */
  mpiret = sc_MPI_Comm_dup (g->glocomm, &p->procomm);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (p->procomm, &p->prorank);
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
  p4est_vtk_write_file (p->pro4est, p->progeom, P4EST_STRING "_producer_new");

  /* do some refinement */

  /* generate a local set of cell values by interpolating a function */
  p4est_iterate (p->pro4est, NULL, p, overlap_producer_compute, NULL
#ifdef P4_TO_P8
                 , NULL
#endif
                 , NULL);

  /* here we would need to make the global partition encoding available to
     the consumer, if the communicators were not congruent */

  /***************************** CONSUMER ****************************/

  P4EST_GLOBAL_PRODUCTION ("OVERLAP: init consumer\n");

  /* setup consumer geometry */
  c->congeom = &c->consumer_geometry;
  c->congeom->name = "consumer";
  c->congeom->user = c;
  c->congeom->X = overlap_consumer_map;
  c->congeom->destroy = (p4est_geometry_destroy_t) 0;

  /* setup consumer app with communicator and mesh */
  mpiret = sc_MPI_Comm_dup (g->glocomm, &c->concomm);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (c->concomm, &c->conrank);
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
  p4est_vtk_write_file (c->con4est, c->congeom, P4EST_STRING "_consumer_new");

  /* do some refinement */

  /* generate a local set of query points */
  c->lquad_idx = 0;
  c->query_xyz = sc_array_new_count (sizeof (overlap_point_t),
                                     c->con4est->local_num_quadrants);
  p4est_iterate (c->con4est, NULL, c, overlap_consumer_compute, NULL
#ifdef P4_TO_P8
                 , NULL
#endif
                 , NULL);
  P4EST_ASSERT (c->lquad_idx == c->con4est->local_num_quadrants);

  /* overlap works by intra-communication on the global communicator */

  P4EST_GLOBAL_PRODUCTION ("OVERLAP: init done\n");
}

static int
consumer_quadrant (p4est_t *p4est, p4est_topidx_t which_tree,
                   p4est_quadrant_t *quadrant, int pfirst, int plast,
                   void *point)
{
#ifdef P4EST_ENABLE_DEBUG
  overlap_consumer_t *c;

  /* Tree, quadrant, pfirst and plast refer to the producer. */

  P4EST_ASSERT (p4est != NULL);
  P4EST_ASSERT (quadrant != NULL);
  P4EST_ASSERT (point == NULL);

  c = (overlap_consumer_t *) p4est->user_pointer;
  P4EST_ASSERT (c != NULL);
#endif

  /* don't mess with the point search */
  return 1;
}

static void
overlap_consumer_add (overlap_consumer_t *c, overlap_point_t *op, int rank)
{
  size_t              bcount;
  overlap_send_buf_t *sb;

  P4EST_ASSERT (c != NULL && c->send_buffer != NULL);
  P4EST_ASSERT (op != NULL && op->rank == -1);
  P4EST_ASSERT (0 <= rank && rank < c->pronum_procs);
  op->rank = rank;

  /* if we have a new rank, push new send buffer */
  bcount = c->send_buffer->elem_count;
  sb = NULL;
  if (bcount > 0) {
    sb = (overlap_send_buf_t *)
      sc_array_index (c->send_buffer, bcount - 1);
    P4EST_ASSERT (sb->rank <= rank);
    P4EST_ASSERT (sb->ops.elem_count > 0);
  }
  if (bcount == 0 || sb->rank < rank) {
    sb = (overlap_send_buf_t *) sc_array_push (c->send_buffer);
    sb->rank = rank;
    sc_array_init (&sb->ops, sizeof (overlap_point_t));
  }
  memcpy (sc_array_push (&sb->ops), op, sizeof (overlap_point_t));
}

static int
producer_intersect (p4est_connectivity_t *pro_conn, p4est_topidx_t which_tree,
                    p4est_quadrant_t *quadrant, overlap_point_t *op, double tol)
{
  const double       *phys;
  double              abc[3], dh, dhz;
  double              qxyz[3];

  phys = op->xyz;

  /* transform point back to producer reference */
  overlap_producer_invmap (pro_conn, which_tree, phys, abc);

  P4EST_LDEBUGF ("Point %ld is %g %g %g\n",
                 (long) op->lnum, phys[0], phys[1], phys[2]);
  P4EST_LDEBUGF ("Tree %d level %d invert to %g %g %g\n",
                 (int) which_tree, quadrant->level, abc[0], abc[1], abc[2]);

  /* check for tree intersection */
  if ((abc[0] <= -tol || abc[0] >= 1. + tol) ||
      (abc[1] <= -tol || abc[1] >= 1. + tol) ||
      (abc[2] <= -tol || abc[2] >= 1. + tol)) {
    return 0;
  }

  P4EST_LDEBUGF ("Point %ld survive tree\n", (long) op->lnum);

  /* check for quadrant intersection */
  dh = OVERLAP_IROOTLEN * P4EST_QUADRANT_LEN (quadrant->level);
  qxyz[0] = OVERLAP_IROOTLEN * quadrant->x;
  qxyz[1] = OVERLAP_IROOTLEN * quadrant->y;
#ifndef P4_TO_P8
  qxyz[2] = 0.;
  dhz = 1.;
#else
  qxyz[2] = OVERLAP_IROOTLEN * quadrant->z;
  dhz = dh;
#endif
  if ((abc[0] < qxyz[0] - tol || abc[0] > qxyz[0] + dh + tol) ||
      (abc[1] < qxyz[1] - tol || abc[1] > qxyz[1] + dh + tol) ||
      (abc[2] < qxyz[2] - tol || abc[2] > qxyz[2] + dhz + tol)) {
    return 0;
  }

  P4EST_LDEBUGF ("Point %ld survive quadrant\n", (long) op->lnum);
  return 1;
}

static int
consumer_point (p4est_t *p4est, p4est_topidx_t which_tree,
                p4est_quadrant_t *quadrant, int pfirst, int plast,
                void *point)
{
  overlap_point_t    *op = (overlap_point_t *) point;
  overlap_consumer_t *c;
  int                 intersects;

  /* The point is owned by the consumer.
     Tree, quadrant, pfirst and plast refer to the producer. */

  P4EST_ASSERT (p4est != NULL);
  P4EST_ASSERT (quadrant != NULL);
  P4EST_ASSERT (op != NULL);
  if (op->rank >= 0) {
    /* skip a point of multiple intersections */
    P4EST_INFOF ("Skip point %ld rank %d on multiple match\n",
                 (long) op->lnum, op->rank);
    return 0;
  }

  c = (overlap_consumer_t *) p4est->user_pointer;
  P4EST_ASSERT (c != NULL);
  P4EST_ASSERT (0 <= which_tree && which_tree < c->pronum_trees);

  /* check if the point intersects the quadrant */
  P4EST_LDEBUGF ("Consumer point %ld intersection test\n", (long) op->lnum);
  intersects =
    producer_intersect (c->producer_conn, which_tree, quadrant, op,
                        P4EST_CON_TOLERANCE);
  if (!intersects) {
    return 0;
  }

  /* we have located the point in the intersection quadrant */
  if (pfirst == plast) {
    /* we have intersected with a leaf quadrant */
    P4EST_INFOF ("Point %ld leaf intersect tree %d rank %d\n",
                 (long) op->lnum, (int) which_tree, pfirst);
    overlap_consumer_add (c, op, pfirst);
  }
  return 1;
}

static int
producer_point (p4est_t *p4est, p4est_topidx_t which_tree,
                p4est_quadrant_t *quadrant, p4est_locidx_t local_num,
                void *point)
{
  int                 isleaf, intersects;
  overlap_point_t    *op = (overlap_point_t *) point;
  overlap_producer_t *p = (overlap_producer_t *) p4est->user_pointer;

  P4EST_ASSERT (p4est != NULL);
  P4EST_ASSERT (quadrant != NULL);
  P4EST_ASSERT (op != NULL);
  P4EST_ASSERT (p != NULL);
  if (op->prodata.isset) {
    /* skip a point of multiple intersections */
    P4EST_INFOF ("Skip point %ld on multiple match\n", (long) op->lnum);
    return 0;
  }

  /* check if the point intersects the quadrant */
  P4EST_LDEBUGF ("Producer point %ld intersection test\n", (long) op->lnum);
  intersects =
    producer_intersect (p->proconn, which_tree, quadrant, op,
                        P4EST_PRO_TOLERANCE);
  if (!intersects) {
    return 0;
  }

  isleaf = local_num >= 0;
  if (isleaf) {
    overlap_prodata_t  *d = (overlap_prodata_t *) quadrant->p.user_data;
    P4EST_ASSERT (d != NULL);
    op->prodata.myvalue = d->myvalue;
    op->prodata.isset = 1;
    P4EST_LDEBUGF ("Producer point %ld prodata set to %f.\n", (long) op->lnum,
                   op->prodata.myvalue);
  }

  return 1;
}

static void
consumer_update_query (sc_array_t *query_xyz, sc_array_t *buffer, int bi)
{
  overlap_recv_buf_t *rb;
  overlap_point_t    *op, *qp;
  int                 j;

  /* obtain the array of points we want to update query_xyz with */
  P4EST_ASSERT (0 <= bi && bi < (int) buffer->elem_size);
  rb = (overlap_recv_buf_t *) sc_array_index_int (buffer, bi);

  /* copy prodata into the query-point array */
  for (j = 0; j < (int) rb->ops.elem_count; ++j) {
    op = (overlap_point_t *) sc_array_index_int (&(rb->ops), j);
    qp = (overlap_point_t *) sc_array_index_int (query_xyz,
                                                 (size_t) op->lnum);
    qp->prodata.isset = op->prodata.isset;
    qp->prodata = op->prodata;
  }
}

static void
overlap_exchange (overlap_global_t * g)
{
  overlap_producer_t *p = g->p;
  overlap_consumer_t *c = g->c;
#ifdef P4EST_ENABLE_DEBUG
  overlap_point_t    *qp;
#endif
  int                 iconrank;
#if defined(P4EST_ENABLE_DEBUG) || defined(P4EST_ENABLE_MPI)
  int                 i;
#endif
  overlap_send_buf_t *sb;
#ifdef P4EST_ENABLE_MPI
  overlap_recv_buf_t *rb;
  size_t              bcount, bz;
#ifdef P4EST_ENABLE_DEBUG
  int                 prev_rank;
#endif
  int                 mpiret;
  int                 remaining, received;
  int                *prod_indices, *cons_indices;
  int                 num_receivers, num_senders, num_ops;
  int                 iprorank;
  int                 same_rank;
  sc_array_t         *receivers, *senders;
  sc_array_t         *payload_in, *payload_out;
  sc_array_t         *prod_recv_reqs, *cons_recv_reqs;
  sc_array_t         *prod_send_reqs, *cons_send_reqs;
#endif

  P4EST_GLOBAL_PRODUCTION ("OVERLAP: exchange partition\n");

  /* consumer receives global partition encoding from producer */
  /* since their communicators are congruent, this is a copy */
  c->producer_gfq = p->pro4est->global_first_quadrant;
  c->producer_gfp = p->pro4est->global_first_position;
  c->pronum_procs = p->pro4est->mpisize;
  c->pronum_trees = (c->producer_conn = p->proconn)->num_trees;

  P4EST_GLOBAL_PRODUCTION ("OVERLAP: customer partition search\n");

  c->send_buffer = sc_array_new (sizeof (overlap_send_buf_t));
  p4est_search_partition_gfx (c->producer_gfq, c->producer_gfp,
                              c->pronum_procs, c->pronum_trees, 0,
                              c, consumer_quadrant, consumer_point,
                              c->query_xyz);

#ifdef P4EST_ENABLE_MPI
  /* assemble and execute receiver and payload query */
  num_receivers = (int) (bcount = c->send_buffer->elem_count);
  receivers = sc_array_new_count (sizeof (int), bcount);
  senders = sc_array_new (sizeof (p4est_locidx_t));
  payload_in = sc_array_new_count (sizeof (int), bcount);
  payload_out = sc_array_new (sizeof (p4est_locidx_t));
  for (bz = 0; bz < bcount; ++bz) {
    sb = (overlap_send_buf_t *) sc_array_index (c->send_buffer, bz);
    *(int *) sc_array_index (receivers, bz) = sb->rank;
    *(p4est_locidx_t *) sc_array_index (payload_in, bz) =
      (p4est_locidx_t) sb->ops.elem_count;
  }
  sc_notify_ext (receivers, senders, payload_in, payload_out, g->glocomm);
  num_senders = (int) senders->elem_count;
  P4EST_INFOF ("Overlap exchange receivers %d senders %d\n",
               num_receivers, num_senders);

  /* post nonblocking receives for the point data of the consumer side */
  p->recv_buffer = sc_array_new_count (sizeof (overlap_recv_buf_t),
                                       num_senders);
  prod_recv_reqs = sc_array_new_count (sizeof (sc_MPI_Request), num_senders);
  iprorank = iconrank = -1;
  for (i = 0; i < num_senders; ++i) {
    /* initalize and allocate the buffer according to the payload */
    rb = (overlap_recv_buf_t *) sc_array_index_int (p->recv_buffer, i);
    rb->rank = *(int *) sc_array_index_int (senders, i);
    same_rank = (rb->rank == p->prorank);
    num_ops = same_rank ? 0 : *(int *) sc_array_index_int (payload_out, i);
    sc_array_init_size (&(rb->ops), sizeof (overlap_point_t),
                        (size_t) num_ops);

    if (same_rank) {
      iprorank = i;             /* save the index in the producer buffer */
      *(sc_MPI_Request *) sc_array_index_int (prod_recv_reqs, i) =
        sc_MPI_REQUEST_NULL;
      continue;
    }

    /* receive the array of overlap_point_t data and store it in the buffer */
    mpiret =
      sc_MPI_Irecv (rb->ops.array, num_ops * sizeof (overlap_point_t),
                    sc_MPI_BYTE, rb->rank, COMM_TAG_CONSDATA, g->glocomm,
                    (sc_MPI_Request *) sc_array_index_int (prod_recv_reqs,
                                                           i));
    SC_CHECK_MPI (mpiret);
  }

  /* send the point data to the producer side in a nonblocking way */
  cons_send_reqs =
    sc_array_new_count (sizeof (sc_MPI_Request), num_receivers);
  for (i = 0; i < num_receivers; ++i) {
    sb = (overlap_send_buf_t *) sc_array_index_int (c->send_buffer, i);

    if (sb->rank == c->conrank) {
      iconrank = i;             /* save the index in the consumer buffer */
      *(sc_MPI_Request *) sc_array_index_int (cons_send_reqs, i) =
        sc_MPI_REQUEST_NULL;
      continue;
    }

    mpiret =
      sc_MPI_Isend (sb->ops.array,
                    sb->ops.elem_count * sizeof (overlap_point_t),
                    sc_MPI_BYTE, sb->rank, COMM_TAG_CONSDATA, g->glocomm,
                    (sc_MPI_Request *) sc_array_index_int (cons_send_reqs,
                                                           i));
    SC_CHECK_MPI (mpiret);
  }

  /* recv the updated point data from the producer side in a nonblocking way */
  cons_recv_reqs =
    sc_array_new_count (sizeof (sc_MPI_Request), num_receivers);
  c->recv_buffer = sc_array_new_size (sizeof (overlap_recv_buf_t),
                                      c->send_buffer->elem_count);
  for (i = 0; i < num_receivers; ++i) {
    rb = (overlap_recv_buf_t *) sc_array_index_int (c->recv_buffer, i);
    sb = (overlap_send_buf_t *) sc_array_index_int (c->send_buffer, i);
    rb->rank = sb->rank;
    same_rank = (rb->rank == c->conrank);
    num_ops = same_rank ? 0 : (int) sb->ops.elem_count;
    sc_array_init_size (&(rb->ops), sizeof (overlap_point_t),
                        (size_t) num_ops);

    if (same_rank) {
      *(sc_MPI_Request *) sc_array_index_int (cons_recv_reqs, iconrank) =
        sc_MPI_REQUEST_NULL;
      continue;
    }

    /* receive the array of overlap_point_t data and store it in the buffer */
    mpiret =
      sc_MPI_Irecv (rb->ops.array,
                    rb->ops.elem_count * sizeof (overlap_point_t),
                    sc_MPI_BYTE, rb->rank, COMM_TAG_PRODATA, g->glocomm,
                    (sc_MPI_Request *) sc_array_index_int (cons_recv_reqs,
                                                           i));
    SC_CHECK_MPI (mpiret);
  }

  /* compute producer data for all incoming messages as soon as they come in */
  P4EST_GLOBAL_PRODUCTION ("OVERLAP: producer local search\n");
  prod_indices = P4EST_ALLOC (int, prod_recv_reqs->elem_count);
  prod_send_reqs = sc_array_new_count (sizeof (sc_MPI_Request), num_senders);
  remaining = num_senders;
  if (iconrank >= 0) {
    *(sc_MPI_Request *) sc_array_index_int (prod_send_reqs, iprorank) =
      sc_MPI_REQUEST_NULL;
    remaining--;      /* since we set the iprorank-th request to null earlier */
  }
  while (remaining > 0) {
    mpiret =
      sc_MPI_Waitsome (num_senders, (sc_MPI_Request *) prod_recv_reqs->array,
                       &received, prod_indices, sc_MPI_STATUSES_IGNORE);
    SC_CHECK_MPI (mpiret);
    P4EST_ASSERT (received != sc_MPI_UNDEFINED);
    P4EST_ASSERT (received > 0);

    for (i = 0; i < received; ++i) {
      /* compute the prodata for the points sent by process prod_indices[i] */
      P4EST_ASSERT (0 <= prod_indices[i] && prod_indices[i] < num_senders);
      rb = (overlap_recv_buf_t *) sc_array_index_int (p->recv_buffer,
                                                      prod_indices[i]);
      p4est_search_local (p->pro4est, 0, NULL, producer_point, &(rb->ops));

      /* send the requested producer data back in a nonblocking way */
      mpiret =
        sc_MPI_Isend (rb->ops.array,
                      rb->ops.elem_count * sizeof (overlap_point_t),
                      sc_MPI_BYTE, rb->rank, COMM_TAG_PRODATA, g->glocomm,
                      (sc_MPI_Request *) sc_array_index_int (prod_send_reqs,
                                                             prod_indices
                                                             [i]));
      SC_CHECK_MPI (mpiret);
    }

    remaining -= received;
  }

  /* compute producer data for all incoming messages as soon as they come in */
  P4EST_GLOBAL_PRODUCTION ("OVERLAP: consumer query point update\n");
  cons_indices = P4EST_ALLOC (int, cons_recv_reqs->elem_count);
  remaining = (iconrank >= 0) ? num_receivers - 1 : num_receivers;
  while (remaining > 0) {
    mpiret =
      sc_MPI_Waitsome (num_receivers,
                       (sc_MPI_Request *) cons_recv_reqs->array, &received,
                       cons_indices, sc_MPI_STATUSES_IGNORE);
    SC_CHECK_MPI (mpiret);
    P4EST_ASSERT (received != sc_MPI_UNDEFINED);
    P4EST_ASSERT (received > 0);

    for (i = 0; i < received; ++i) {
      consumer_update_query (c->query_xyz, c->recv_buffer, cons_indices[i]);
    }

    remaining -= received;
  }

  /* wait for the nonblocking sends to complete */
  mpiret =
    sc_MPI_Waitall (num_receivers, (sc_MPI_Request *) cons_send_reqs->array,
                    sc_MPI_STATUSES_IGNORE);
  SC_CHECK_MPI (mpiret);
  mpiret =
    sc_MPI_Waitall (num_senders, (sc_MPI_Request *) prod_send_reqs->array,
                    sc_MPI_STATUSES_IGNORE);
  SC_CHECK_MPI (mpiret);
#else /* !P4EST_ENABLE_MPI */
  iconrank = 0;      /* indicate that the send buffer can be updated directly */
#endif

  if (iconrank >= 0) {
    /* Interpolate point-data of local points. Instead of copying to the
     * producer buffer, we update the points in-place. */
    sb = (overlap_send_buf_t *) sc_array_index_int (c->send_buffer, iconrank);
    p4est_search_local (p->pro4est, 0, NULL, producer_point, &(sb->ops));
    consumer_update_query (c->query_xyz, c->send_buffer, iconrank);
  }

#ifdef P4EST_ENABLE_DEBUG
  for (i = 0; i < (int) c->con4est->local_num_quadrants; i++) {
    qp = (overlap_point_t *) sc_array_index_int (c->query_xyz, i);
    if (qp->prodata.isset) {
      P4EST_LDEBUGF ("Consumer query point %ld assigned prodata %f\n",
                     (long) qp->lnum, qp->prodata.myvalue);
    }
    else {
      P4EST_LDEBUGF ("Consumer query point %ld not in producer domain.\n",
                     (long) qp->lnum);
    }
  }
#endif

#ifdef P4EST_ENABLE_MPI
  /* free remaining communication data */
  P4EST_FREE (prod_indices);
  P4EST_FREE (cons_indices);
  sc_array_destroy (prod_recv_reqs);
  sc_array_destroy (cons_recv_reqs);
  sc_array_destroy (prod_send_reqs);
  sc_array_destroy (cons_send_reqs);
  sc_array_destroy (receivers);
  sc_array_destroy (senders);
  sc_array_destroy (payload_in);
  sc_array_destroy (payload_out);
#ifdef P4EST_ENABLE_DEBUG
  prev_rank = -1;
#endif
  for (bz = 0; bz < bcount; ++bz) {
    sb = (overlap_send_buf_t *) sc_array_index (c->send_buffer, bz);
    rb = (overlap_recv_buf_t *) sc_array_index (c->recv_buffer, bz);
    SC_ASSERT (sb->rank == rb->rank);
    SC_ASSERT (prev_rank < sb->rank);
#ifdef P4EST_ENABLE_DEBUG
    prev_rank = sb->rank;
    if (bz == (size_t) iconrank) {
      P4EST_ASSERT (rb->ops.elem_count == 0);
    }
    else {
      P4EST_ASSERT (rb->ops.elem_count == sb->ops.elem_count);
    }
#endif
    P4EST_ASSERT (sb->ops.elem_count > 0);
    sc_array_reset (&sb->ops);
    sc_array_reset (&rb->ops);
  }
  sc_array_destroy_null (&c->recv_buffer);
  for (i = 0; i < num_senders; ++i) {
    rb = (overlap_recv_buf_t *) sc_array_index_int (p->recv_buffer, i);
#ifdef P4EST_ENABLE_DEBUG
    if (i == iprorank) {
      P4EST_ASSERT (rb->ops.elem_count == 0);
    }
    else {
      P4EST_ASSERT (rb->ops.elem_count > 0);
    }
#endif
    sc_array_reset (&rb->ops);
  }
  sc_array_destroy_null (&p->recv_buffer);
#else /* !P4EST_ENABLE_MPI */
  sb = (overlap_send_buf_t *) sc_array_index_int (c->send_buffer, 0);
  sc_array_reset (&sb->ops);
#endif
  sc_array_destroy_null (&c->send_buffer);
}

static void
overlap_update (overlap_global_t *g)
{
}

static void
overlap_apps_reset (overlap_global_t *g)
{
  overlap_producer_t *p = g->p;
  overlap_consumer_t *c = g->c;
  int                 mpiret;

  /* destroy producer */
  p4est_destroy (p->pro4est);
  p4est_connectivity_destroy (p->proconn);
  mpiret = sc_MPI_Comm_free (&p->procomm);
  SC_CHECK_MPI (mpiret);

  /* destroy consumer */
  sc_array_destroy (c->query_xyz);
  p4est_destroy (c->con4est);
  p4est_connectivity_destroy (c->conconn);
  mpiret = sc_MPI_Comm_free (&c->concomm);
  SC_CHECK_MPI (mpiret);
}

int
main (int argc, char **argv)
{
  int                 i;
  int                 mpiret;
  sc_MPI_Comm         mpicomm;
  overlap_global_t global, *g = &global;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  mpicomm = sc_MPI_COMM_WORLD;
  sc_init (mpicomm, 1, 1, NULL, SC_LP_DEFAULT);

  overlap_apps_init (g, mpicomm);

  overlap_exchange (g);
  for (i = 0; i < g->rounds; ++i) {
    P4EST_GLOBAL_PRODUCTIONF ("Into round %d/%d\n", i, g->rounds);
    overlap_update (g);
    overlap_exchange (g);
  }

  overlap_apps_reset (g);

  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
