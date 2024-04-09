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
#include <sc_options.h>
#include <sc_flops.h>
#ifndef P4_TO_P8
#include <p4est_extended.h>
#include <p4est_search.h>
#include <p4est_vtk.h>
#include <p4est_bits.h>
#else
#include <p8est_extended.h>
#include <p8est_search.h>
#include <p8est_vtk.h>
#include <p8est_bits.h>
#endif

/* ---------------------------------------------------------------------- */
///               TIMING- AND STAT-CONTEXT FOR OVERLAP EXCHANGE
/* ---------------------------------------------------------------------- */
/* Context and auxiliary functions for in-depth timing of exchange. */

/* If set to 1, the time spent in the partition and the local search callbacks
 * is measured. By default this is disabled, since it strongly affects the total
 * run time. */
#define MEASURE_CALLBACKS 0

enum
{
  OVERLAP_EXCHANGE,
  OVERLAP_SEARCH_PARTITION,
#ifdef P4EST_ENABLE_MPI
  OVERLAP_NOTIFY,
  OVERLAP_PARTITION_NOTIFY,
  OVERLAP_NUM_PROCS_SENT,
  OVERLAP_NUM_PROCS_RECEIVED,
  OVERLAP_POST_MESSAGES,
  OVERLAP_INTERPOLATE,
  OVERLAP_UPDATE_QUERY_POINTS,
  OVERLAP_WAITALL,
  OVERLAP_UPDATE_NONLOCAL,
#endif
  OVERLAP_UPDATE_LOCAL,
  OVERLAP_UPDATE_TOTAL,
  OVERLAP_FREE_COMMUNICATION_DATA,
  OVERLAP_NUM_QP,
  OVERLAP_NUM_QP_SENT,
  OVERLAP_NUM_QP_RECEIVED,
  OVERLAP_NUM_QP_SENTRECVD,
#if MEASURE_CALLBACKS
  OVERLAP_NUM_PROD_SEARCH_OPS,
  OVERLAP_NUM_CONS_SEARCH_OPS,
  OVERLAP_NUM_SEARCH_OPS,
  OVERLAP_CONS_SEARCH_CALLBACK,
  OVERLAP_PROD_SEARCH_CALLBACK,
  OVERLAP_PROD_INTERPOLATION_CALLBACK,
#endif
  OVERLAP_SEARCH_LOCAL,
  OVERLAP_NUM_STATS
};

/* data types of the different stats: 0 - double, 1 - integer */
static int          overlap_stats_type[OVERLAP_NUM_STATS] = { 0, 0,
#ifdef P4EST_ENABLE_MPI
  0, 0, 1, 1, 0, 0, 0, 0, 0,
#endif
  0, 0, 0, 1, 1, 1, 1,
#if MEASURE_CALLBACKS
  1, 1, 1, 0, 0, 0,
#endif
  0
};

/* basic timing context to pass around between the overlap-functions */
typedef struct overlap_tstats
{
  sc_flopinfo_t       fi;
  sc_statinfo_t       stats[OVERLAP_NUM_STATS];
}
overlap_tstats_t;

/* print stats of different data types (integer or double) */
static void
sc_stats_print_x (int package_id, int log_priority, int nvars,
                  sc_statinfo_t *stats, int *stats_type, int full,
                  int summary)
{
  int                 i, ti, count;
  sc_statinfo_t      *si;
  char                buffer[BUFSIZ];

  if (full) {
    for (i = 0; i < nvars; ++i) {
      si = &stats[i];
      ti = stats_type[i];
      /* begin printing */
      if (ti) {                 /* the stat is integer */
        if (si->variable != NULL) {
          SC_GEN_LOGF (package_id, SC_LC_GLOBAL, log_priority,
                       "Statistics for   %s\n", si->variable);
        }
        else {
          SC_GEN_LOGF (package_id, SC_LC_GLOBAL, log_priority,
                       "Statistics for %d\n", i);
        }
        SC_GEN_LOGF (package_id, SC_LC_GLOBAL, log_priority,
                     "   Global number of values: %7ld\n", si->count);
        if (!si->count) {
          continue;
        }
        if (si->average != 0.) {        /* ignore the comparison warning */
          SC_GEN_LOGF (package_id, SC_LC_GLOBAL, log_priority,
                       "   Mean value (std. dev.):           %.0f (%.3g = %.3g%%)\n",
                       si->average, si->standev,
                       100. * si->standev / fabs (si->average));
        }
        else {
          SC_GEN_LOGF (package_id, SC_LC_GLOBAL, log_priority,
                       "   Mean value (std. dev.):           %.0f (%.3g)\n",
                       si->average, si->standev);
        }
        SC_GEN_LOGF (package_id, SC_LC_GLOBAL, log_priority,
                     "   Minimum attained at rank %7d: %.0f\n",
                     si->min_at_rank, si->min);
        SC_GEN_LOGF (package_id, SC_LC_GLOBAL, log_priority,
                     "   Maximum attained at rank %7d: %.0f\n",
                     si->max_at_rank, si->max);
      }
      else {
        if (si->variable != NULL) {
          SC_GEN_LOGF (package_id, SC_LC_GLOBAL, log_priority,
                       "Statistics for   %s\n", si->variable);
        }
        else {
          SC_GEN_LOGF (package_id, SC_LC_GLOBAL, log_priority,
                       "Statistics for %d\n", i);
        }
        SC_GEN_LOGF (package_id, SC_LC_GLOBAL, log_priority,
                     "   Global number of values: %7ld\n", si->count);
        if (!si->count) {
          continue;
        }
        if (si->average != 0.) {        /* ignore the comparison warning */
          SC_GEN_LOGF (package_id, SC_LC_GLOBAL, log_priority,
                       "   Mean value (std. dev.):           %g (%.3g = %.3g%%)\n",
                       si->average, si->standev,
                       100. * si->standev / fabs (si->average));
        }
        else {
          SC_GEN_LOGF (package_id, SC_LC_GLOBAL, log_priority,
                       "   Mean value (std. dev.):           %g (%.3g)\n",
                       si->average, si->standev);
        }
        SC_GEN_LOGF (package_id, SC_LC_GLOBAL, log_priority,
                     "   Minimum attained at rank %7d: %g\n",
                     si->min_at_rank, si->min);
        SC_GEN_LOGF (package_id, SC_LC_GLOBAL, log_priority,
                     "   Maximum attained at rank %7d: %g\n",
                     si->max_at_rank, si->max);
      }
    }
  }
  else {
    for (i = 0; i < nvars; ++i) {
      si = &stats[i];
      ti = stats_type[i];
      if (ti) {
        /* print just the average */
        if (si->variable != NULL) {
          snprintf (buffer, BUFSIZ, "for %s:", si->variable);
        }
        else {
          snprintf (buffer, BUFSIZ, "for %3d:", i);
        }
        if (si->average != 0.) {        /* ignore the comparison warning */
          SC_GEN_LOGF (package_id, SC_LC_GLOBAL, log_priority,
                       "Mean (sigma) %-23s %.0f (%.3g = %.3g%%)\n",
                       buffer, si->average, si->standev,
                       100. * si->standev / fabs (si->average));
        }
        else {
          SC_GEN_LOGF (package_id, SC_LC_GLOBAL, log_priority,
                       "Mean (sigma) %-23s %.0f (%.3g)\n", buffer,
                       si->average, si->standev);
        }
      }
      else {
        /* print just the average */
        if (si->variable != NULL) {
          snprintf (buffer, BUFSIZ, "for %s:", si->variable);
        }
        else {
          snprintf (buffer, BUFSIZ, "for %3d:", i);
        }
        if (si->average != 0.) {        /* ignore the comparison warning */
          SC_GEN_LOGF (package_id, SC_LC_GLOBAL, log_priority,
                       "Mean (sigma) %-23s %g (%.3g = %.3g%%)\n",
                       buffer, si->average, si->standev,
                       100. * si->standev / fabs (si->average));
        }
        else {
          SC_GEN_LOGF (package_id, SC_LC_GLOBAL, log_priority,
                       "Mean (sigma) %-23s %g (%.3g)\n", buffer,
                       si->average, si->standev);
        }
      }
    }
  }

  /* the summary always contains all variables */
  if (summary) {
    count = snprintf (buffer, BUFSIZ, "Summary = ");
    for (i = 0; i < nvars && count >= 0 && (size_t) count < BUFSIZ; ++i) {
      si = &stats[i];
      ti = stats_type[i];
      if (ti) {
        count += snprintf (buffer + count, BUFSIZ - count, "%s%.0f",
                           i == 0 ? "[ " : " ", si->average);
      }
      else {
        count += snprintf (buffer + count, BUFSIZ - count, "%s%g",
                           i == 0 ? "[ " : " ", si->average);
      }
    }
    if (count >= 0 && (size_t) count < BUFSIZ) {
      snprintf (buffer + count, BUFSIZ - count, "%s", " ];\n");
      SC_GEN_LOG (package_id, SC_LC_GLOBAL, log_priority, buffer);
    }
    else {
      SC_GEN_LOG (package_id, SC_LC_GLOBAL, log_priority,
                  "Summary overflow\n");
    }
    count = snprintf (buffer, BUFSIZ, "Maximum = ");
    for (i = 0; i < nvars && count >= 0 && (size_t) count < BUFSIZ; ++i) {
      si = &stats[i];
      ti = stats_type[i];
      if (ti) {
        count += snprintf (buffer + count, BUFSIZ - count, "%s%.0f",
                           i == 0 ? "[ " : " ", si->max);
      }
      else {
        count += snprintf (buffer + count, BUFSIZ - count, "%s%g",
                           i == 0 ? "[ " : " ", si->max);
      }
    }
    if (count >= 0 && (size_t) count < BUFSIZ) {
      snprintf (buffer + count, BUFSIZ - count, "%s", " ];\n");
      SC_GEN_LOG (package_id, SC_LC_GLOBAL, log_priority, buffer);
    }
    else {
      SC_GEN_LOG (package_id, SC_LC_GLOBAL, log_priority,
                  "Maximum overflow\n");
    }
  }
}

/* ---------------------------------------------------------------------- */
///                            OVERLAP EXCHANGE
/* ---------------------------------------------------------------------- */
/* Abstract exchange routine for unknown point type (passed around as void
 * pointer and handeĺed only by user-defined callbacks). */

/** Callback to be passed to \ref overlap_exchange.
 *
 * This callback is supposed to determine if a given user-defined \a point is
 * contained in a given \a quadrant and return the result.
 * It will be called both in a local search of the actual p4est as well as in a
 * partition search of the artificially reconstructed p4est.
 * If the p4est is artifical can be determined using \ref overlap_p4est_is_meta.
 * For an artifical p4est, the intersection test should only rely on the
 * geometrical information provided by \a quadrant, since any additional
 * information might not be available.
 */
typedef int         (*overlap_intersect_point_t) (p4est_t *p4est,
                                                  p4est_topidx_t which_tree,
                                                  p4est_quadrant_t *quadrant,
                                                  p4est_locidx_t lnum,
                                                  void *point, void *user);

/** Callback to be passed to \ref overlap_exchange.
 *
 * This callback is supposed to evaluate a given user-defined \a point that is
 * contained in a given leaf \a quadrant.
 * Otherwise, similar to \ref overlap_intersect_point_t, but only called for
 * a real \a p4est.
 */
typedef overlap_intersect_point_t overlap_interpolate_point_t;

/** Determine if a given p4est is real or artifical.
 *
 * This is an auxiliary function for the user intended to be used inside a
 * \ref overlap_intersect_point_t to determine if the p4est passed to the
 * callback is a real, producer-side p4est or an artifical, consumer-side p4est.
 * \param [in] p4est             A potentially artifical p4est.
 * \return True, iff \a p4est is artifical.
 */
static int
overlap_p4est_is_meta (p4est_t *p4est)
{
  P4EST_ASSERT (p4est != NULL);
  return (p4est->local_num_quadrants == -1);
}

typedef struct overlap_producer
{
  /* mesh constituents */
  p4est_t            *pro4est;

  /* work variables */
  p4est_locidx_t      lquad_idx;
  size_t              point_size;

  /* communication */
  sc_MPI_Comm         glocomm;
  sc_MPI_Comm         procomm;
  int                 prorank;
  int                 iprorank;
  sc_array_t         *recv_buffer;
  sc_array_t         *recv_reqs;
  sc_array_t         *send_reqs;

  /* timings */
  overlap_tstats_t   *tstats;

  /* callbacks */
  overlap_intersect_point_t intersect;
  overlap_interpolate_point_t interpolate;
  void               *user_pointer;     /* provided to callbacks */
}
overlap_producer_t;

typedef struct overlap_consumer
{
  /* work variables */
  p4est_locidx_t      lquad_idx;
  sc_array_t         *query_xyz;
  size_t              point_size;

  /* minimal knowledge of the producer's mesh */
  p4est_t            *pro4est;
  p4est_connectivity_t *producer_conn;
  const p4est_quadrant_t *producer_gfp;
  int                 pronum_procs;
  p4est_topidx_t      pronum_trees;

  /* communication */
  sc_MPI_Comm         glocomm;
  sc_MPI_Comm         concomm;
  int                 conrank;
  int                 iconrank;
  sc_array_t         *send_buffer;
  sc_array_t         *send_reqs;
  sc_array_t         *recv_buffer;
  sc_array_t         *recv_reqs;

  /* timings */
  overlap_tstats_t   *tstats;

  /* callbacks */
  overlap_intersect_point_t intersect;
  void               *user_pointer;     /* provided to callbacks */
}
overlap_consumer_t;

typedef struct overlap_global
{
  sc_MPI_Comm         glocomm;
  overlap_tstats_t   *tstats;
  overlap_producer_t  pro, *p;
  overlap_consumer_t  con, *c;
}
overlap_global_t;

typedef struct overlap_consumer_point
{
  void               *point;
  int                 rank;
  p4est_locidx_t      lnum;
}
overlap_consumer_point_t;

typedef struct overlap_producer_point
{
  void               *point;
  int                 isset;
}
overlap_producer_point_t;

typedef struct overlap_buf
{
  int                 rank;
  sc_array_t          ops;
  sc_array_t          lnums;
}
overlap_buf_t;

typedef enum overlap_comm_tag
{
  COMM_TAG_CONSDATA = P4EST_COMM_TAG_LAST,
  COMM_TAG_PRODATA
}
overlap_comm_tag_t;

static int
overlap_consumer_quadrant_fn (p4est_t *p4est, p4est_topidx_t which_tree,
                              p4est_quadrant_t *quadrant, int pfirst,
                              int plast, void *point)
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
overlap_consumer_add (overlap_consumer_t *c, overlap_consumer_point_t *op,
                      int rank)
{
  size_t              bcount;
  overlap_buf_t      *sb;

  P4EST_ASSERT (c != NULL && c->send_buffer != NULL);
  P4EST_ASSERT (op != NULL && op->rank == -1);
  P4EST_ASSERT (0 <= rank && rank < c->pronum_procs);
  op->rank = rank;

  /* if we have a new rank, push new send buffer */
  bcount = c->send_buffer->elem_count;
  sb = NULL;
  if (bcount > 0) {
    sb = (overlap_buf_t *)
      sc_array_index (c->send_buffer, bcount - 1);
    P4EST_ASSERT (sb->rank <= rank);
    P4EST_ASSERT (sb->ops.elem_count > 0);
  }
  if (bcount == 0 || sb->rank < rank) {
    sb = (overlap_buf_t *) sc_array_push (c->send_buffer);
    sb->rank = rank;
    sc_array_init (&sb->ops, c->point_size);
    sc_array_init (&sb->lnums, sizeof (p4est_locidx_t));
  }
  memcpy (sc_array_push (&sb->ops), op->point, c->point_size);
  memcpy (sc_array_push (&sb->lnums), &op->lnum, sizeof (p4est_locidx_t));
}

static int
overlap_consumer_point_fn (p4est_t *p4est, p4est_topidx_t which_tree,
                           p4est_quadrant_t *quadrant, int pfirst, int plast,
                           void *point)
{
  overlap_consumer_t *c;
  overlap_consumer_point_t *op;
  int                 intersects;

  /* The point is owned by the consumer.
     Tree, quadrant, pfirst and plast refer to the producer. */

  P4EST_ASSERT (p4est != NULL);
  P4EST_ASSERT (p4est->user_pointer != NULL);
  c = (overlap_consumer_t *) p4est->user_pointer;
  P4EST_ASSERT (quadrant != NULL);
  P4EST_ASSERT (point != NULL);
  op = (overlap_consumer_point_t *) point;

#if MEASURE_CALLBACKS
  sc_flopinfo_t       snapshot;
  c->tstats->stats[OVERLAP_NUM_CONS_SEARCH_OPS].sum_values++;
  sc_flops_snap (&c->tstats->fi, &snapshot);
#endif

  if (op->rank >= 0) {
    /* skip a point of multiple intersections */
#if MEASURE_CALLBACKS
    sc_flops_shot (&c->tstats->fi, &snapshot);
    c->tstats->stats[OVERLAP_CONS_SEARCH_CALLBACK].sum_values +=
      snapshot.iwtime;
#endif
    return 0;
  }

  P4EST_ASSERT (c != NULL);
  P4EST_ASSERT (0 <= which_tree && which_tree < c->pronum_trees);

  /* check if the point intersects the quadrant */
  P4EST_LDEBUGF ("Consumer point %ld intersection test\n", (long) op->lnum);
  c->pro4est->mpirank = (pfirst == plast) ? pfirst : -1;
  intersects =
    c->intersect (c->pro4est, which_tree, quadrant, -1, op->point,
                  c->user_pointer);
  if (!intersects) {
#if MEASURE_CALLBACKS
    sc_flops_shot (&c->tstats->fi, &snapshot);
    c->tstats->stats[OVERLAP_CONS_SEARCH_CALLBACK].sum_values +=
      snapshot.iwtime;
#endif
    return 0;
  }

  /* we have located the point in the intersection quadrant */
  if (pfirst == plast) {
    /* we have intersected with a leaf quadrant */
    overlap_consumer_add (c, op, pfirst);
  }
#if MEASURE_CALLBACKS
  sc_flops_shot (&c->tstats->fi, &snapshot);
  c->tstats->stats[OVERLAP_CONS_SEARCH_CALLBACK].sum_values +=
    snapshot.iwtime;
#endif
  return 1;
}

static void
overlap_consumer_search_partition (overlap_consumer_t *c)
{
  size_t              iz, nipz;
  sc_array_t         *query_points;
  overlap_consumer_point_t *op;

  nipz = c->query_xyz->elem_count;
  query_points = sc_array_new_count (sizeof (overlap_consumer_point_t), nipz);
  for (iz = 0; iz < nipz; ++iz) {
    /* wrap anonymous input points in struct with rank and local numbering */
    op = (overlap_consumer_point_t *) sc_array_index (query_points, iz);
    op->point = sc_array_index (c->query_xyz, iz);
    op->rank = -1;
    op->lnum = iz;
  }

  c->send_buffer = sc_array_new (sizeof (overlap_buf_t));
  p4est_search_partition_gfp (c->producer_gfp, c->pronum_procs, c->pronum_trees,
                              0, c, overlap_consumer_quadrant_fn,
                              overlap_consumer_point_fn, query_points);

  sc_array_destroy (query_points);
}

#ifdef P4EST_ENABLE_MPI
static void
overlap_consumer_producer_notify (overlap_global_t *g)
{
  overlap_producer_t *p = g->p;
  overlap_consumer_t *c = g->c;
  size_t              bz, bcount;
  sc_array_t         *receivers, *senders;
  sc_array_t         *payload_in, *payload_out;
  int                 num_senders;
  overlap_buf_t      *sb;
  overlap_buf_t      *rb;
  int                 same_rank, num_ops, i;
  int                 mpiret;
  sc_notify_t        *notifyc;
  sc_flopinfo_t       snapshot;

  /* assemble and execute receiver and payload query */
  bcount = c->send_buffer->elem_count;
  receivers = sc_array_new_count (sizeof (int), bcount);
  senders = sc_array_new (sizeof (p4est_locidx_t));
  payload_in = sc_array_new_count (sizeof (int), bcount);
  payload_out = sc_array_new (sizeof (p4est_locidx_t));
  for (bz = 0; bz < bcount; ++bz) {
    sb = (overlap_buf_t *) sc_array_index (c->send_buffer, bz);
    *(int *) sc_array_index (receivers, bz) = sb->rank;
    *(p4est_locidx_t *) sc_array_index (payload_in, bz) =
      (p4est_locidx_t) sb->ops.elem_count;
    g->tstats->stats[OVERLAP_NUM_QP_SENT].sum_values += sb->ops.elem_count;
  }

  sc_flops_snap (&g->tstats->fi, &snapshot);
  notifyc = sc_notify_new (g->glocomm);
  sc_notify_set_type (notifyc, SC_NOTIFY_NBX);
  sc_notify_payload (receivers, senders, payload_in, payload_out, 1, notifyc);
  sc_notify_destroy (notifyc);
  sc_flops_shot (&g->tstats->fi, &snapshot);
  sc_stats_set1 (&g->tstats->stats[OVERLAP_NOTIFY], snapshot.iwtime,
                 "Consumer producer notify");

  num_senders = (int) senders->elem_count;

  /* post nonblocking receives for the point data of the consumer side */
  p->recv_buffer = sc_array_new_count (sizeof (overlap_buf_t), num_senders);
  p->recv_reqs = sc_array_new_count (sizeof (sc_MPI_Request), num_senders);
  p->iprorank = c->iconrank = -1;
  for (i = 0; i < num_senders; ++i) {
    /* initalize and allocate the buffer according to the payload */
    rb = (overlap_buf_t *) sc_array_index_int (p->recv_buffer, i);
    rb->rank = *(int *) sc_array_index_int (senders, i);
    same_rank = (rb->rank == p->prorank);
    g->tstats->stats[OVERLAP_NUM_QP_RECEIVED].sum_values +=
      *(int *) sc_array_index_int (payload_out, i);
    num_ops = same_rank ? 0 : *(int *) sc_array_index_int (payload_out, i);
    sc_array_init_size (&(rb->ops), p->point_size, (size_t) num_ops);

    if (same_rank) {
      p->iprorank = i;          /* save the index in the producer buffer */
      *(sc_MPI_Request *) sc_array_index_int (p->recv_reqs, i) =
        sc_MPI_REQUEST_NULL;
      continue;
    }

    /* receive the array of simple_point_t data and store it in the buffer */
    mpiret =
      sc_MPI_Irecv (rb->ops.array, num_ops * p->point_size,
                    sc_MPI_BYTE, rb->rank, COMM_TAG_CONSDATA, p->glocomm,
                    (sc_MPI_Request *) sc_array_index_int (p->recv_reqs, i));
    SC_CHECK_MPI (mpiret);
  }

  sc_array_destroy (receivers);
  sc_array_destroy (senders);
  sc_array_destroy (payload_in);
  sc_array_destroy (payload_out);

}

static void
overlap_consumer_post_messages (overlap_consumer_t *c)
{
  overlap_buf_t      *sb;
  overlap_buf_t      *rb;
  int                 num_receivers, same_rank, num_ops, i;
  int                 mpiret;

  /* send the point data to the producer side in a nonblocking way */
  num_receivers = (int) c->send_buffer->elem_count;
  c->send_reqs = sc_array_new_count (sizeof (sc_MPI_Request), num_receivers);
  for (i = 0; i < num_receivers; ++i) {
    sb = (overlap_buf_t *) sc_array_index_int (c->send_buffer, i);

    if (sb->rank == c->conrank) {
      c->iconrank = i;          /* save the index in the consumer buffer */
      *(sc_MPI_Request *) sc_array_index_int (c->send_reqs, i) =
        sc_MPI_REQUEST_NULL;
      continue;
    }

    mpiret =
      sc_MPI_Isend (sb->ops.array,
                    sb->ops.elem_count * c->point_size,
                    sc_MPI_BYTE, sb->rank, COMM_TAG_CONSDATA, c->glocomm,
                    (sc_MPI_Request *) sc_array_index_int (c->send_reqs, i));
    SC_CHECK_MPI (mpiret);
  }

  /* recv the updated point data from the producer side in a nonblocking way */
  c->recv_reqs = sc_array_new_count (sizeof (sc_MPI_Request), num_receivers);
  c->recv_buffer = sc_array_new_size (sizeof (overlap_buf_t), num_receivers);
  for (i = 0; i < num_receivers; ++i) {
    rb = (overlap_buf_t *) sc_array_index_int (c->recv_buffer, i);
    sb = (overlap_buf_t *) sc_array_index_int (c->send_buffer, i);
    rb->rank = sb->rank;
    same_rank = (rb->rank == c->conrank);
    num_ops = same_rank ? 0 : (int) sb->ops.elem_count;
    sc_array_init_size (&(rb->ops), c->point_size, (size_t) num_ops);
    rb->lnums = sb->lnums;

    if (same_rank) {
      *(sc_MPI_Request *) sc_array_index_int (c->recv_reqs, c->iconrank) =
        sc_MPI_REQUEST_NULL;
      continue;
    }

    /* receive the array of simple_point_t data and store it in the buffer */
    mpiret =
      sc_MPI_Irecv (rb->ops.array,
                    rb->ops.elem_count * c->point_size,
                    sc_MPI_BYTE, rb->rank, COMM_TAG_PRODATA, c->glocomm,
                    (sc_MPI_Request *) sc_array_index_int (c->recv_reqs, i));
    SC_CHECK_MPI (mpiret);
  }
}
#endif /* P4EST_ENABLE_MPI */

static int
overlap_producer_point_fn (p4est_t *p4est, p4est_topidx_t which_tree,
                           p4est_quadrant_t *quadrant,
                           p4est_locidx_t local_num, void *point)
{
  int                 isleaf, intersects;
  overlap_producer_point_t *op;
  overlap_producer_t *p;

  P4EST_ASSERT (p4est != NULL);
  P4EST_ASSERT (quadrant != NULL);
  P4EST_ASSERT (p4est->user_pointer != NULL);
  p = (overlap_producer_t *) p4est->user_pointer;
  P4EST_ASSERT (point != NULL);
  op = (overlap_producer_point_t *) point;

#if MEASURE_CALLBACKS
  sc_flopinfo_t       snapshot, snapshot2;
  p->tstats->stats[OVERLAP_NUM_PROD_SEARCH_OPS].sum_values++;
  sc_flops_snap (&p->tstats->fi, &snapshot);
#endif

  if (op->isset) {
    /* skip a point of multiple intersections */
#if MEASURE_CALLBACKS
    sc_flops_shot (&p->tstats->fi, &snapshot);
    p->tstats->stats[OVERLAP_PROD_SEARCH_CALLBACK].sum_values +=
      snapshot.iwtime;
#endif
    return 0;
  }

  /* check if the point intersects the quadrant */
  P4EST_ASSERT (p->pro4est->connectivity != NULL);
  intersects =
    p->intersect (p->pro4est, which_tree, quadrant, local_num,
                  op->point, p->user_pointer);
  if (!intersects) {
#if MEASURE_CALLBACKS
    sc_flops_shot (&p->tstats->fi, &snapshot);
    p->tstats->stats[OVERLAP_PROD_SEARCH_CALLBACK].sum_values +=
      snapshot.iwtime;
#endif
    return 0;
  }

  isleaf = local_num >= 0;
  if (isleaf) {
#if MEASURE_CALLBACKS
    sc_flops_snap (&p->tstats->fi, &snapshot2);
#endif
    op->isset = 1;
    p->interpolate (p4est, which_tree, quadrant, local_num, op->point,
                    p->user_pointer);
#if MEASURE_CALLBACKS
    sc_flops_shot (&p->tstats->fi, &snapshot2);
    p->tstats->stats[OVERLAP_PROD_INTERPOLATION_CALLBACK].sum_values +=
      snapshot.iwtime;
#endif
  }

#if MEASURE_CALLBACKS
  sc_flops_shot (&p->tstats->fi, &snapshot);
  p->tstats->stats[OVERLAP_PROD_SEARCH_CALLBACK].sum_values +=
    snapshot.iwtime;
#endif
  return 1;
}

static void
overlap_producer_search_local (overlap_producer_t *p, sc_array_t *points)
{
  size_t              iz, nipz;
  sc_array_t         *query_points;
  overlap_producer_point_t *op;

  nipz = points->elem_count;
  query_points = sc_array_new_count (sizeof (overlap_producer_point_t), nipz);
  for (iz = 0; iz < nipz; ++iz) {
    /* wrap anonymous input points in struct with isset-marker */
    op = (overlap_producer_point_t *) sc_array_index (query_points, iz);
    op->point = sc_array_index (points, iz);
    op->isset = 0;
  }

  p4est_search_local (p->pro4est, 0, NULL, overlap_producer_point_fn,
                      query_points);

  sc_array_destroy (query_points);
}

#ifdef P4EST_ENABLE_MPI
static void
overlap_producer_interpolate (overlap_producer_t *p)
{
  overlap_buf_t      *rb;
  int                 num_senders, i;
  int                 remaining, received;
  int                *prod_indices;
  int                 mpiret;
  sc_flopinfo_t       snapshot;

  /* compute producer data for all incoming messages as soon as they come in */
  num_senders = (int) p->recv_reqs->elem_count;
  prod_indices = P4EST_ALLOC (int, num_senders);
  p->send_reqs = sc_array_new_count (sizeof (sc_MPI_Request), num_senders);
  remaining = num_senders;
  if (p->iprorank >= 0) {
    *(sc_MPI_Request *) sc_array_index_int (p->send_reqs, p->iprorank) =
      sc_MPI_REQUEST_NULL;
    remaining--;                /* since we set the iprorank-th request to null earlier */
  }
  while (remaining > 0) {
    mpiret =
      sc_MPI_Waitsome (num_senders, (sc_MPI_Request *) p->recv_reqs->array,
                       &received, prod_indices, sc_MPI_STATUSES_IGNORE);
    SC_CHECK_MPI (mpiret);
    P4EST_ASSERT (received != sc_MPI_UNDEFINED);
    P4EST_ASSERT (received > 0);

    for (i = 0; i < received; ++i) {
      /* compute the prodata for the points sent by process prod_indices[i] */
      P4EST_ASSERT (0 <= prod_indices[i] && prod_indices[i] < num_senders);
      rb = (overlap_buf_t *) sc_array_index_int (p->recv_buffer,
                                                 prod_indices[i]);

      sc_flops_snap (&p->tstats->fi, &snapshot);
      overlap_producer_search_local (p, &(rb->ops));
      sc_flops_shot (&p->tstats->fi, &snapshot);
      p->tstats->stats[OVERLAP_SEARCH_LOCAL].sum_values += snapshot.iwtime;
      /* send the requested producer data back in a nonblocking way */
      mpiret =
        sc_MPI_Isend (rb->ops.array,
                      rb->ops.elem_count * p->point_size,
                      sc_MPI_BYTE, rb->rank, COMM_TAG_PRODATA, p->glocomm,
                      (sc_MPI_Request *) sc_array_index_int (p->send_reqs,
                                                             prod_indices
                                                             [i]));
      SC_CHECK_MPI (mpiret);
    }

    remaining -= received;
  }

  P4EST_FREE (prod_indices);
}

#endif /* P4EST_ENABLE_MPI */

static void
overlap_consumer_update_from_buffer
  (overlap_consumer_t *c, sc_array_t *buffer, int bi)
{
  overlap_buf_t      *rb;
  void               *updated_point, *point;
  int                 j;
  p4est_locidx_t      lnum;

  /* obtain the array of points we want to update query_xyz with */
  P4EST_ASSERT (0 <= bi && bi < (int) buffer->elem_size);
  rb = (overlap_buf_t *) sc_array_index_int (buffer, bi);

  /* copy prodata into the query-point array */
  for (j = 0; j < (int) rb->ops.elem_count; ++j) {
    lnum = *(p4est_locidx_t *) sc_array_index_int (&(rb->lnums), j);
    point = sc_array_index_int (c->query_xyz, (size_t) lnum);
    updated_point = sc_array_index_int (&(rb->ops), j);
    memcpy (point, updated_point, c->point_size);
  }
}

#ifdef P4EST_ENABLE_MPI

static void
overlap_consumer_update_query_points (overlap_consumer_t *c)
{
  int                 num_receivers, i;
  int                 remaining, received;
  int                *cons_indices;
  int                 mpiret;

  /* compute producer data for all incoming messages as soon as they come in */
  num_receivers = (int) c->recv_reqs->elem_count;
  cons_indices = P4EST_ALLOC (int, num_receivers);
  remaining = (c->iconrank >= 0) ? num_receivers - 1 : num_receivers;
  while (remaining > 0) {
    mpiret =
      sc_MPI_Waitsome (num_receivers,
                       (sc_MPI_Request *) c->recv_reqs->array, &received,
                       cons_indices, sc_MPI_STATUSES_IGNORE);
    SC_CHECK_MPI (mpiret);
    P4EST_ASSERT (received != sc_MPI_UNDEFINED);
    P4EST_ASSERT (received > 0);

    for (i = 0; i < received; ++i) {
      overlap_consumer_update_from_buffer (c, c->recv_buffer,
                                           cons_indices[i]);
    }

    remaining -= received;
  }

  P4EST_FREE (cons_indices);
}

static void
overlap_consumer_waitall (overlap_consumer_t *c)
{
  int                 mpiret;
  int                 num_receivers;

  /* wait for the nonblocking sends to complete */
  num_receivers = (int) c->send_reqs->elem_count;
  mpiret =
    sc_MPI_Waitall (num_receivers, (sc_MPI_Request *) c->send_reqs->array,
                    sc_MPI_STATUSES_IGNORE);
  SC_CHECK_MPI (mpiret);
}

static void
overlap_producer_waitall (overlap_producer_t *p)
{
  int                 mpiret;
  int                 num_senders;

  /* wait for the nonblocking sends to complete */
  num_senders = (int) p->send_reqs->elem_count;
  mpiret =
    sc_MPI_Waitall (num_senders, (sc_MPI_Request *) p->send_reqs->array,
                    sc_MPI_STATUSES_IGNORE);
  SC_CHECK_MPI (mpiret);
}

#endif /* P4EST_ENABLE_MPI */

static void
overlap_consumer_producer_update_local (overlap_global_t *g)
{
  overlap_consumer_t *c = g->c;
  overlap_producer_t *p = g->p;
  overlap_buf_t      *sb;
  sc_flopinfo_t       snapshot;

  if (c->iconrank >= 0 && c->send_buffer->elem_count) {
    /* Interpolate point-data of local points. Instead of copying to the
     * producer buffer, we update the points in-place. */
    sb = (overlap_buf_t *) sc_array_index_int (c->send_buffer, c->iconrank);
    sc_flops_snap (&g->tstats->fi, &snapshot);
    overlap_producer_search_local (p, &(sb->ops));
    sc_flops_shot (&g->tstats->fi, &snapshot);
    g->tstats->stats[OVERLAP_SEARCH_LOCAL].sum_values += snapshot.iwtime;
    overlap_consumer_update_from_buffer (c, c->send_buffer, c->iconrank);
  }
}

static void
overlap_consumer_free_communication_data (overlap_consumer_t *c)
{
  overlap_buf_t      *sb;
#ifdef P4EST_ENABLE_MPI
  overlap_buf_t      *rb;
#ifdef P4EST_ENABLE_DEBUG
  int                 prev_rank;
#endif
  size_t              bz, bcount;

  sc_array_destroy (c->recv_reqs);
  sc_array_destroy (c->send_reqs);
#ifdef P4EST_ENABLE_DEBUG
  prev_rank = -1;
#endif
  bcount = c->send_buffer->elem_count;
  for (bz = 0; bz < bcount; ++bz) {
    sb = (overlap_buf_t *) sc_array_index (c->send_buffer, bz);
    rb = (overlap_buf_t *) sc_array_index (c->recv_buffer, bz);
    SC_ASSERT (sb->rank == rb->rank);
    SC_ASSERT (prev_rank < sb->rank);
#ifdef P4EST_ENABLE_DEBUG
    prev_rank = sb->rank;
    if (bz == (size_t) c->iconrank) {
      P4EST_ASSERT (rb->ops.elem_count == 0);
    }
    else {
      P4EST_ASSERT (rb->ops.elem_count == sb->ops.elem_count);
    }
#endif
    P4EST_ASSERT (sb->ops.elem_count > 0);
    sc_array_reset (&sb->ops);
    sc_array_reset (&sb->lnums);        /* rb->lnums == sb->lnums */
    sc_array_reset (&rb->ops);
  }
  sc_array_destroy_null (&c->recv_buffer);
#else /* !P4EST_ENABLE_MPI */
  if (c->send_buffer->elem_count) {
    sb = (overlap_buf_t *) sc_array_index_int (c->send_buffer, 0);
    sc_array_reset (&sb->ops);
    sc_array_reset (&sb->lnums);
  }
#endif
  sc_array_destroy_null (&c->send_buffer);
  P4EST_FREE (c->pro4est);
}

static void
overlap_producer_free_communication_data (overlap_producer_t *p)
{
#ifdef P4EST_ENABLE_MPI
  overlap_buf_t      *rb;
  int                 num_senders, i;

  sc_array_destroy (p->recv_reqs);
  sc_array_destroy (p->send_reqs);
  num_senders = (int) p->recv_buffer->elem_count;
  for (i = 0; i < num_senders; ++i) {
    rb = (overlap_buf_t *) sc_array_index_int (p->recv_buffer, i);
#ifdef P4EST_ENABLE_DEBUG
    if (i == p->iprorank) {
      P4EST_ASSERT (rb->ops.elem_count == 0);
    }
    else {
      P4EST_ASSERT (rb->ops.elem_count > 0);
    }
#endif
    sc_array_reset (&rb->ops);
  }
  sc_array_destroy_null (&p->recv_buffer);
#endif
}

/** Exchange data between a p4est and another mesh discretized by query points.
 *
 * The \a p4est is considered a producer, which provides data to another mesh.
 * The other mesh is considered a consumer and is represented by a parallel
 * distributed set of user-defined query points \a points.
 * The query points are searched in the producer p4est using the \a intersect
 * callback and sent to the respective producer process containing them.
 * On the producer process they are searche in the local part of \a p4est and
 * evaluated using the \a interpolate callback when found in a leaf quadrant.
 * Finally, the potentially updated query points are returned to their original
 * consumer process.
 * \param [in,out] pro4est  The producer p4est.
 * \param [in,out] points   A parallel distributed set of query points created
 *                          based on the consumer mesh (e.g. stemming from
 *                          quadrature). The points are not changed by this
 *                          function, they just get copied, sent and passed to
 *                          \a intersect and \a interpolate.
 * \param [in] concomm      The consumer side communicator.
 * \param [in] glocomm      Global communicator for communication between
 *                          producer and consumer.
 * \param [in] intersect    Callback function that decides for a given quadrant
 *                          and a query point from \a points, if the point lies
 *                          inside the quadrant.
 * \param [in] interpolate  Callback function that evaluates a query point for
 *                          a given producer side quadrant it is contained in.
 * \param [in,out] user     User pointer provided to all calls of \a intersect
 *                          and \a interpolate.
 */
static void
overlap_exchange (p4est_t *pro4est, sc_array_t *points,
                  sc_MPI_Comm concomm, sc_MPI_Comm glocomm,
                  overlap_intersect_point_t intersect,
                  overlap_interpolate_point_t interpolate, void *user)
{
  int                 istat;
  sc_flopinfo_t       snapshot, snapshot2, snapshot3, snapshot_total, *fi;
  int                 mpiret;
  void               *pro4est_user_pointer;
  sc_statinfo_t      *stats;
  overlap_tstats_t    tstats;
  overlap_global_t global, *g = &global;
  overlap_producer_t *p = g->p = &g->pro;
  overlap_consumer_t *c = g->c = &g->con;

  /* start overall timing */
  mpiret = sc_MPI_Barrier (glocomm);
  SC_CHECK_MPI (mpiret);
  sc_flops_start (&tstats.fi);
  p->tstats = c->tstats = g->tstats = &tstats;
  stats = g->tstats->stats;

  P4EST_ASSERT (pro4est != NULL);       /* currently we assume congruent comms */
  P4EST_ASSERT (points != NULL);

  /* initialize global_context */
  g->glocomm = glocomm;

  /* initialize producer context */
  p->glocomm = glocomm;
  p->procomm = pro4est->mpicomm;
  p->prorank = pro4est->mpirank;
  p->pro4est = pro4est;
  pro4est_user_pointer = pro4est->user_pointer; /* save previous user_pointer */
  p->pro4est->user_pointer = p;
  p->lquad_idx = 0;
  p->point_size = points->elem_size;
  p->intersect = intersect;
  p->interpolate = interpolate;
  p->user_pointer = user;

  /*initialize consumer context */
  c->glocomm = glocomm;
  c->concomm = concomm;
  sc_MPI_Comm_rank (concomm, &c->conrank);
  c->query_xyz = points;
  c->point_size = points->elem_size;
  c->intersect = intersect;
  c->user_pointer = user;

  /* initialize counters to zero */
  c->tstats->stats[OVERLAP_NUM_QP_SENT].sum_values = 0;
  c->tstats->stats[OVERLAP_NUM_QP_RECEIVED].sum_values = 0;
#if MEASURE_CALLBACKS
  c->tstats->stats[OVERLAP_NUM_CONS_SEARCH_OPS].sum_values = 0;
  c->tstats->stats[OVERLAP_NUM_PROD_SEARCH_OPS].sum_values = 0;
  c->tstats->stats[OVERLAP_CONS_SEARCH_CALLBACK].sum_values = 0;
  c->tstats->stats[OVERLAP_PROD_SEARCH_CALLBACK].sum_values = 0;
  c->tstats->stats[OVERLAP_PROD_INTERPOLATION_CALLBACK].sum_values = 0;
#endif
  c->tstats->stats[OVERLAP_SEARCH_LOCAL].sum_values = 0;

  /* make sure all processes entered the function */
  mpiret = sc_MPI_Barrier (g->glocomm);
  SC_CHECK_MPI (mpiret);

  /* total time of the exchange function */
  fi = &g->tstats->fi;
  sc_flops_snap (fi, &snapshot_total);
  P4EST_GLOBAL_PRODUCTION ("OVERLAP: exchange partition\n");

  /* consumer receives global partition encoding from producer */
  /* since their communicators are congruent, this is a copy */
  c->producer_gfp = p->pro4est->global_first_position;
  c->pronum_procs = p->pro4est->mpisize;
  c->pronum_trees = (c->producer_conn = p->pro4est->connectivity)->num_trees;
  c->pro4est = P4EST_ALLOC_ZERO (p4est_t, 1);
  c->pro4est->connectivity = c->producer_conn;
  c->pro4est->local_num_quadrants = -1; /* marks p4est as meta */

  P4EST_GLOBAL_PRODUCTION ("OVERLAP: consumer partition search\n");

  /* search for the query points in the producer-partition and create a buffer
   * to send them to the respective producer ranks */
  sc_flops_snap (fi, &snapshot2);
  sc_flops_snap (fi, &snapshot);
  overlap_consumer_search_partition (c);
  sc_flops_shot (fi, &snapshot);
  sc_stats_set1 (&stats[OVERLAP_SEARCH_PARTITION], snapshot.iwtime,
                 "Search partition");

#ifdef P4EST_ENABLE_MPI
  /* global, communication-based part of the interpolation */
  /* during this process we will mark messages that allow for a local, in-place
   * solution by setting c->iconrank */

  /* notify the producer about the point-array-messages it will receive,
   * allocate an receive buffer according to the transmitted payloads and
   * post Irecvs for the point-arrays */
  overlap_consumer_producer_notify (g);
  sc_flops_shot (fi, &snapshot2);
  sc_stats_set1 (&stats[OVERLAP_PARTITION_NOTIFY], snapshot2.iwtime,
                 "Search partition and notify");
  sc_stats_set1 (&stats[OVERLAP_NUM_PROCS_SENT],
                 (double) g->c->send_buffer->elem_count,
                 "Consumer number processes sent to");
  sc_stats_set1 (&stats[OVERLAP_NUM_PROCS_RECEIVED],
                 (double) g->p->recv_buffer->elem_count,
                 "Producer number processes received from");

  /* post Isends for the point-arrays as well as Irecvs for the updated
   * point-arrays containing the interpolation prodata */
  sc_flops_snap (fi, &snapshot3);
  sc_flops_snap (fi, &snapshot2);
  sc_flops_snap (fi, &snapshot);
  overlap_consumer_post_messages (c);
  sc_flops_shot (fi, &snapshot);
  sc_stats_set1 (&stats[OVERLAP_POST_MESSAGES], snapshot.iwtime,
                 "Consumer post messages");

  P4EST_GLOBAL_PRODUCTION ("OVERLAP: producer local search\n");

  /* interpolate the point-arrays as soon as they arrive and send them back to
   * the consumer side in a non-blocking way */
  sc_flops_snap (fi, &snapshot);
  overlap_producer_interpolate (p);
  sc_flops_shot (fi, &snapshot);
  sc_stats_set1 (&stats[OVERLAP_INTERPOLATE], snapshot.iwtime,
                 "Producer interpolate");

#else /* !P4EST_ENABLE_MPI */
  sc_flops_snap (fi, &snapshot3);
  c->iconrank = 0;              /* indicate that the send buffer can be updated directly */
#endif

  /* local, in-place part of the interpolation */
  sc_flops_snap (fi, &snapshot);
  overlap_consumer_producer_update_local (g);
  sc_flops_shot (fi, &snapshot);
  sc_stats_set1 (&stats[OVERLAP_UPDATE_LOCAL], snapshot.iwtime,
                 "Consumer producer update local");

#ifdef P4EST_ENABLE_MPI
  P4EST_GLOBAL_PRODUCTION ("OVERLAP: consumer query point update\n");
  /* compute the interpolation data of the query points based on the
   * updated point-arrays */
  sc_flops_snap (fi, &snapshot);
  overlap_consumer_update_query_points (c);
  sc_flops_shot (fi, &snapshot);
  sc_stats_set1 (&stats[OVERLAP_UPDATE_QUERY_POINTS],
                 snapshot.iwtime, "Consumer update query points");

  /* wait for the communication to complete */
  sc_flops_snap (fi, &snapshot);
  overlap_consumer_waitall (c);
  overlap_producer_waitall (p);
  sc_flops_shot (fi, &snapshot);
  sc_flops_shot (fi, &snapshot2);
  sc_stats_set1 (&stats[OVERLAP_UPDATE_NONLOCAL], snapshot2.iwtime,
                 "Consumer producer update nonlocal");
  sc_stats_set1 (&stats[OVERLAP_WAITALL], snapshot.iwtime,
                 "Consumer producer waitall");
#endif /* P4EST_ENABLE_MPI */

  P4EST_GLOBAL_PRODUCTION ("OVERLAP: local interpolation\n");

  /* local, in-place part of the interpolation */
  sc_flops_shot (fi, &snapshot3);
  sc_stats_set1 (&stats[OVERLAP_UPDATE_TOTAL], snapshot3.iwtime,
                 "Consumer producer update total");

  /* free remaining communication data */
  sc_flops_snap (fi, &snapshot);
  overlap_consumer_free_communication_data (c);
  overlap_producer_free_communication_data (p);
  sc_flops_shot (fi, &snapshot);
  sc_stats_set1 (&stats[OVERLAP_FREE_COMMUNICATION_DATA],
                 snapshot.iwtime,
                 "Consumer producer free communication data");

  /* finish timings and stats */
  sc_flops_shot (&g->tstats->fi, &snapshot_total);
  sc_stats_set1 (&stats[OVERLAP_EXCHANGE], snapshot_total.iwtime, "Exchange");
  sc_stats_set1 (&stats[OVERLAP_NUM_QP],
                 g->c->query_xyz->elem_count,
                 "Number local consumer quadrants");
  sc_stats_set1 (&stats[OVERLAP_NUM_QP_SENTRECVD],
                 stats[OVERLAP_NUM_QP_SENT].sum_values
                 + stats[OVERLAP_NUM_QP_RECEIVED].sum_values,
                 "Number query point sent and received");
  sc_stats_set1 (&stats[OVERLAP_NUM_QP_SENT],
                 stats[OVERLAP_NUM_QP_SENT].sum_values,
                 "Number query points sent to producer side");
  sc_stats_set1 (&stats[OVERLAP_NUM_QP_RECEIVED],
                 stats[OVERLAP_NUM_QP_RECEIVED].sum_values,
                 "Number query points received on producer side");
#if MEASURE_CALLBACKS
  sc_stats_set1 (&stats[OVERLAP_NUM_SEARCH_OPS],
                 stats[OVERLAP_NUM_CONS_SEARCH_OPS].sum_values
                 + stats[OVERLAP_NUM_PROD_SEARCH_OPS].sum_values,
                 "Number callback calls in all searches");
  sc_stats_set1 (&stats[OVERLAP_NUM_CONS_SEARCH_OPS],
                 stats[OVERLAP_NUM_CONS_SEARCH_OPS].sum_values,
                 "Number callback calls in partition search");
  sc_stats_set1 (&stats[OVERLAP_NUM_PROD_SEARCH_OPS],
                 stats[OVERLAP_NUM_PROD_SEARCH_OPS].sum_values,
                 "Number callback calls in local search");
  sc_stats_set1 (&stats[OVERLAP_CONS_SEARCH_CALLBACK],
                 stats[OVERLAP_CONS_SEARCH_CALLBACK].sum_values,
                 "Time spent in partition search callback");
  sc_stats_set1 (&stats[OVERLAP_PROD_SEARCH_CALLBACK],
                 stats[OVERLAP_PROD_SEARCH_CALLBACK].sum_values,
                 "Time spent in local search callback");
  sc_stats_set1 (&stats[OVERLAP_PROD_INTERPOLATION_CALLBACK],
                 stats[OVERLAP_PROD_INTERPOLATION_CALLBACK].sum_values,
                 "Time spent in interpolation callback");
#endif
  sc_stats_set1 (&stats[OVERLAP_SEARCH_LOCAL],
                 stats[OVERLAP_SEARCH_LOCAL].sum_values, "Search local");
  sc_stats_compute (g->glocomm, OVERLAP_NUM_STATS, stats);
  /* sc_stats_print_x works the same as sc_stats_print, but takes an array
   * that indicates, if the stat is a double or an integer, to decide between
   * %g and %f. We use a hardcoded overlap_stats_type for all OVERLAP_NUM_STATS
   * stats */
  sc_stats_print_x (p4est_package_id, SC_LP_ESSENTIAL, OVERLAP_NUM_STATS,
                    stats, overlap_stats_type, 1, 1);

  for (istat = 0; istat < OVERLAP_NUM_STATS; istat++) {
    sc_stats_reset (&stats[istat], 0);
  }

  /* reset user pointer of producer p4est */
  p->pro4est->user_pointer = pro4est_user_pointer;
}

/* ---------------------------------------------------------------------- */
///                     Application Producer and Consumer 
/* ---------------------------------------------------------------------- */
/* Structs to handle the producer and the consumer mesh of a rather general
 * overlap setting. Used both for the simple and the adaptive example. */

typedef struct producer
{
  /* mesh constituents */
  sc_MPI_Comm         procomm;
  p4est_connectivity_t *proconn;
  p4est_t            *pro4est;
  p4est_geometry_t   *progeom, producer_geometry;

  /* parameters */
  int                 pminl;

  /* work variables */
  p4est_locidx_t      lquad_idx;

  /* communication */
  sc_MPI_Comm         glocomm;
  int                 prorank;

  /* vtk cell data */
  sc_array_t         *interpolation_data;
}
producer_t;

typedef struct consumer
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

  /* communication */
  sc_MPI_Comm         glocomm;
  int                 conrank;

  /* vtk cell data */
  sc_array_t         *interpolation_data;
  sc_array_t         *isset_data;
  sc_array_t         *xyz_data;
}
consumer_t;

typedef struct global
{
  /* mesh overset context */
  sc_MPI_Comm         glocomm;
  producer_t          pro, *p;
  consumer_t          con, *c;
  void               *usr_ctx;

  /* application settings */
  int                 refinement_method;
  int                 refinement_maxlevel;
  int                 example;
  int                 output_vtk;
  int                 output_text;
  int                 proprocs;
  int                 conprocs;
}
global_t;

/* ---------------------------------------------------------------------- */
///                     Mappings and Intersections
/* ---------------------------------------------------------------------- */
/* A general point struct and its intersect callback. Various mappings as well
 * as corresponding inverse mappings are supplied as well. Used as a basis for
 * both the simple and the adaptive example. */

typedef struct coordinate_point
{
  /* coordinates and tree index */
  double              xyz[3];
  int                 which_tree;
  double              inv[3];

  /* local numbering for debugging */
  p4est_locidx_t      lnum;
}
coordinate_point_t;

typedef void        (*coordinate_invmap_t) (p4est_connectivity_t *conn,
                                            p4est_topidx_t which_tree,
                                            coordinate_point_t * cp);

#define COORDINATE_IROOTLEN (1. / P4EST_ROOT_LEN)

#ifdef OVERLAP_WITH_CUBE_MAP    /* cube map not currently used */
#ifdef P4EST_ENABLE_DEBUG

static void
coordinate_cube_invmap (p4est_connectivity_t *conn,
                        p4est_topidx_t which_tree, const double xyz[3],
                        double abc[3]);

#endif /* P4EST_ENABLE_DEBUG */

static void
coordinate_cube_map (p4est_geometry_t *geom, p4est_topidx_t which_tree,
                     const double abc[3], double xyz[3])
{
  double              a, co, si, x;
#ifdef P4EST_ENABLE_DEBUG
  coordinate_point_t  coordinate_point, *cp = &coordinate_point;
#endif
  double             *vert;
  p4est_topidx_t      vind;
  p4est_connectivity_t *conn;

  P4EST_ASSERT (geom != NULL);
  P4EST_ASSERT (geom->X == coordinate_cube_map);

  /* preimage domain is [0, 2]^3 */

  /* access origin vertex of given tree */
  P4EST_ASSERT (geom != NULL);
  P4EST_ASSERT (geom->X == coordinate_cube_map);
  conn = (p4est_connectivity_t *) geom->user;
  P4EST_ASSERT (conn != NULL && conn->vertices != NULL);
  vind = conn->tree_to_vertex[P4EST_CHILDREN * which_tree + 0];
  vert = &conn->vertices[3 * vind + 0];

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
  memset (cp, -1, sizeof (coordinate_point_t));
  memcpy (cp->xyz, xyz, 3 * sizeof (double));
  coordinate_cube_invmap (conn, which_tree, cp);
  P4EST_ASSERT (fabs (abc[0] - cp->inv[0]) < SC_1000_EPS &&
                fabs (abc[1] - cp->inv[1]) < SC_1000_EPS &&
                fabs (abc[2] - cp->inv[2]) < SC_1000_EPS);
#endif
}

static void
coordinate_cube_invmap (p4est_connectivity_t *conn,
                        p4est_topidx_t which_tree, coordinate_point_t *cp)
{
  double              a, co, si;
  double              abc[3];
  double             *vert;
  p4est_topidx_t      vind;

  P4EST_ASSERT (conn != NULL && conn->vertices != NULL);
  P4EST_ASSERT (0 <= which_tree && which_tree < conn->num_trees);
  vind = conn->tree_to_vertex[P4EST_CHILDREN * which_tree + 0];
  vert = &conn->vertices[3 * vind + 0];

  /* rotate 20 degrees backwards around the y axis */
  a = -20. * M_PI / 180.;
  co = cos (a);
  si = sin (a);
  abc[0] = co * cp->xyz[0] - si * cp->xyz[2];
  abc[2] = si * cp->xyz[0] + co * cp->xyz[2];
  abc[1] = cp->xyz[1];

  /* scale slightly and move away from origin into brick */
  cp->inv[0] = abc[0] / 1.1 + .5 - vert[0];
  cp->inv[1] = abc[1] / 1.2 + .5 - vert[1];
  cp->inv[2] = abc[2] / 1.3 + .5 - vert[2];
}

#endif /* 0 */

static void
coordinate_curved_invmap
  (p4est_connectivity_t *conn, p4est_topidx_t which_tree,
   coordinate_point_t *cp);

static int          nbricks[6] = { 1, 1, 1, 1, 1, 1 };

static void
coordinate_curved_map (p4est_geometry_t *geom, p4est_topidx_t which_tree,
                       const double abc[3], double xyz[3])
{
  double             *vert;
  p4est_topidx_t      vind;
#ifdef P4EST_ENABLE_DEBUG
  coordinate_point_t  coordinate_point, *cp = &coordinate_point;
#endif
  p4est_connectivity_t *conn;

  /* map to quadrant in the brick */
  P4EST_ASSERT (geom != NULL);
  P4EST_ASSERT (geom->X == coordinate_curved_map);
  conn = (p4est_connectivity_t *) geom->user;
  P4EST_ASSERT (conn != NULL && conn->vertices != NULL);
  vind = conn->tree_to_vertex[P4EST_CHILDREN * which_tree + 0];
  vert = &conn->vertices[3 * vind + 0];

  xyz[0] = vert[0] + abc[0];
  xyz[1] = vert[1] + abc[1];
  xyz[2] = vert[2] + abc[2];

  /* scale down x and y to avoid octants getting to elongated */
  xyz[0] *= 2. / (double) nbricks[0];
  xyz[1] *= 2. / (double) nbricks[1];
  xyz[2] *= 2. / (double) nbricks[2];

  /* shift y and z and according to a curve of x */
  xyz[1] += 0.2;
  xyz[1] *= 1. / 4.;
  xyz[1] += (0.75 - xyz[0]) * (0.75 - xyz[0]);

  xyz[2] -= 0.2;
  xyz[2] *= 1. / 4.;
  xyz[2] += (0.75 - xyz[0]) * (0.75 - xyz[0]);

#ifdef P4EST_ENABLE_DEBUG
  memset (cp, -1, sizeof (coordinate_point_t));
  memcpy (cp->xyz, xyz, 3 * sizeof (double));
  coordinate_curved_invmap (conn, which_tree, cp);
  P4EST_ASSERT (fabs (abc[0] - cp->inv[0]) < SC_1000_EPS &&
                fabs (abc[1] - cp->inv[1]) < SC_1000_EPS &&
                fabs (abc[2] - cp->inv[2]) < SC_1000_EPS);
#endif
}

static void
coordinate_curved_invmap (p4est_connectivity_t *conn,
                          p4est_topidx_t which_tree, coordinate_point_t *cp)
{
  double             *vert;
  double              abc[3];
  p4est_topidx_t      vind;

  abc[0] = cp->xyz[0];
  abc[1] = cp->xyz[1];
  abc[2] = cp->xyz[2];

  /* invert shifting of y and z according to a curve of x */
  abc[1] -= (0.75 - abc[0]) * (0.75 - abc[0]);
  abc[1] *= 4.;
  abc[1] -= 0.2;

  abc[2] -= (0.75 - abc[0]) * (0.75 - abc[0]);
  abc[2] *= 4.;
  abc[2] += 0.2;

  /* invert scaling of x and y */
  abc[0] *= (double) nbricks[0] / 2.;
  abc[1] *= (double) nbricks[1] / 2.;
  abc[2] *= (double) nbricks[2] / 2.;

  /* invert mapping to brick quadrant */
  P4EST_ASSERT (conn != NULL && conn->vertices != NULL);
  P4EST_ASSERT (0 <= which_tree && which_tree < conn->num_trees);
  vind = conn->tree_to_vertex[P4EST_CHILDREN * which_tree + 0];
  vert = &conn->vertices[3 * vind + 0];

  cp->inv[0] = abc[0] - vert[0];
  cp->inv[1] = abc[1] - vert[1];
  cp->inv[2] = abc[2] - vert[2];
}

static void
coordinate_brick_invmap
  (p4est_connectivity_t *conn, p4est_topidx_t which_tree,
   coordinate_point_t *cp);

static void
coordinate_brick_map (p4est_geometry_t *geom, p4est_topidx_t which_tree,
                      const double abc[3], double xyz[3])
{
  double              a, co, si, x;
#ifdef P4EST_ENABLE_DEBUG
  coordinate_point_t  coordinate_point, *cp = &coordinate_point;
#endif
  double             *vert;
  p4est_topidx_t      vind;
  p4est_connectivity_t *conn;

  /* preimage domain is 3x2x1 with origin in the lower left front */

  /* access origin vertex of given tree */
  P4EST_ASSERT (geom != NULL);
  P4EST_ASSERT (geom->X == coordinate_brick_map);
  conn = (p4est_connectivity_t *) geom->user;
  P4EST_ASSERT (conn != NULL && conn->vertices != NULL);
  vind = conn->tree_to_vertex[P4EST_CHILDREN * which_tree + 0];
  vert = &conn->vertices[3 * vind + 0];

  /* center brick around origin and scale */
  xyz[0] = (vert[0] + abc[0] - 1.5) * 0.7;
  xyz[1] = (vert[1] + abc[1] - 1.0) * 0.6;
  xyz[2] = (vert[2] + abc[2] - 0.5) * 0.5;

  /* rotate 30 degrees around the z axis and shift */
  a = 30. * M_PI / 180.;
  co = cos (a);
  si = sin (a);
  x = xyz[0];
  xyz[0] = co * x - si * xyz[1] + 0.8;
  xyz[1] = si * x + co * xyz[1] + 0.8;
  xyz[2] += 0.5;

#ifdef P4EST_ENABLE_DEBUG
  memset (cp, -1, sizeof (coordinate_point_t));
  memcpy (cp->xyz, xyz, 3 * sizeof (double));
  coordinate_brick_invmap (conn, which_tree, cp);
  P4EST_ASSERT (fabs (abc[0] - cp->inv[0]) < SC_1000_EPS &&
                fabs (abc[1] - cp->inv[1]) < SC_1000_EPS &&
                fabs (abc[2] - cp->inv[2]) < SC_1000_EPS);
#endif
}

static void
coordinate_brick_invmap (p4est_connectivity_t *conn,
                         p4est_topidx_t which_tree, coordinate_point_t *cp)
{
  double              a, co, si, x;
  double              abc[3];
  double             *vert;
  p4est_topidx_t      vind;

  P4EST_ASSERT (conn != NULL && conn->vertices != NULL);
  P4EST_ASSERT (0 <= which_tree && which_tree < conn->num_trees);
  vind = conn->tree_to_vertex[P4EST_CHILDREN * which_tree + 0];
  vert = &conn->vertices[3 * vind + 0];

  /* shift back */
  abc[0] = cp->xyz[0] - 0.8;
  abc[1] = cp->xyz[1] - 0.8;
  abc[2] = cp->xyz[2] - 0.5;

  /* complete 360 degree rotation */
  a = 330. * M_PI / 180.;
  co = cos (a);
  si = sin (a);
  x = abc[0];
  abc[0] = co * x - si * abc[1];
  abc[1] = si * x + co * abc[1];

  /* invert scaling and centering */
  cp->inv[0] = abc[0] / 0.7 + 1.5 - vert[0];
  cp->inv[1] = abc[1] / 0.6 + 1.0 - vert[1];
  cp->inv[2] = abc[2] / 0.5 + 0.5 - vert[2];
}

static void
coordinate_producer_unit_invmap
  (p4est_connectivity_t *conn, p4est_topidx_t which_tree,
   coordinate_point_t *cp);

static void
coordinate_producer_unit_map (p4est_geometry_t *geom,
                              p4est_topidx_t which_tree, const double abc[3],
                              double xyz[3])
{
#ifdef P4EST_ENABLE_DEBUG
  coordinate_point_t  coordinate_point, *cp = &coordinate_point;
#endif
  double             *vert;
  p4est_topidx_t      vind;
  p4est_connectivity_t *conn;

  /* access origin vertex of given tree */
  P4EST_ASSERT (geom != NULL);
  P4EST_ASSERT (geom->X == coordinate_producer_unit_map);
  conn = (p4est_connectivity_t *) geom->user;
  P4EST_ASSERT (conn != NULL && conn->vertices != NULL);
  vind = conn->tree_to_vertex[P4EST_CHILDREN * which_tree + 0];
  vert = &conn->vertices[3 * vind + 0];

  /* scale to unit square */
  xyz[0] = (vert[0] + abc[0]) / (double) nbricks[0];
  xyz[1] = (vert[1] + abc[1]) / (double) nbricks[1];
  xyz[2] = (vert[2] + abc[2]) / (double) nbricks[2];

#ifdef P4EST_ENABLE_DEBUG
  memset (cp, -1, sizeof (coordinate_point_t));
  memcpy (cp->xyz, xyz, 3 * sizeof (double));
  coordinate_producer_unit_invmap (conn, which_tree, cp);
  P4EST_ASSERT (fabs (abc[0] - cp->inv[0]) < SC_1000_EPS &&
                fabs (abc[1] - cp->inv[1]) < SC_1000_EPS &&
                fabs (abc[2] - cp->inv[2]) < SC_1000_EPS);
#endif
}

static void
coordinate_producer_unit_invmap (p4est_connectivity_t *conn,
                                 p4est_topidx_t which_tree,
                                 coordinate_point_t *cp)
{
  double             *vert;
  p4est_topidx_t      vind;

  P4EST_ASSERT (conn != NULL && conn->vertices != NULL);
  P4EST_ASSERT (0 <= which_tree && which_tree < conn->num_trees);
  vind = conn->tree_to_vertex[P4EST_CHILDREN * which_tree + 0];
  vert = &conn->vertices[3 * vind + 0];

  /* invert scaling and centering */
  cp->inv[0] = cp->xyz[0] * (double) nbricks[0] - vert[0];
  cp->inv[1] = cp->xyz[1] * (double) nbricks[1] - vert[1];
  cp->inv[2] = cp->xyz[2] * (double) nbricks[2] - vert[2];
}

static void
coordinate_consumer_unit_map (p4est_geometry_t *geom,
                              p4est_topidx_t which_tree, const double abc[3],
                              double xyz[3])
{
  double             *vert;
  double              def[3];
  p4est_topidx_t      vind;
  p4est_connectivity_t *conn;

  /* access origin vertex of given tree */
  P4EST_ASSERT (geom != NULL);
  P4EST_ASSERT (geom->X == coordinate_consumer_unit_map);
  conn = (p4est_connectivity_t *) geom->user;
  P4EST_ASSERT (conn != NULL && conn->vertices != NULL);
  vind = conn->tree_to_vertex[P4EST_CHILDREN * which_tree + 0];
  vert = &conn->vertices[3 * vind + 0];

  /* scale to unit square */
  def[0] = (vert[0] + abc[0]) / (double) nbricks[3];
  def[1] = (vert[1] + abc[1]) / (double) nbricks[4];
  def[2] = (vert[2] + abc[2]) / (double) nbricks[5];

  /* rotate */
  xyz[0] = 1. - def[1];
  xyz[1] = 1. - def[0];
  xyz[2] = def[2];
}

static int
coordinate_intersect_fn (p4est_t *p4est, p4est_topidx_t which_tree,
                         p4est_quadrant_t *quadrant, p4est_locidx_t lnum,
                         coordinate_point_t *cp, void *user)
{
  const double       *phys;
  double              dh, dhz;
  double              qxyz[3];
  double              tol;
  coordinate_invmap_t invmap;

  P4EST_ASSERT (p4est != NULL);
  P4EST_ASSERT (quadrant != NULL);
  P4EST_ASSERT (cp != NULL);
  P4EST_ASSERT (user != NULL);
  invmap = (coordinate_invmap_t) user;

  /* transform point back to producer reference */
  if (cp->which_tree != which_tree) {
    /* we enter a new tree in the search and have a new inverse mapping */
    invmap (p4est->connectivity, which_tree, cp);
    cp->which_tree = which_tree;
  }

  phys = cp->xyz;
  P4EST_LDEBUGF ("Point %ld is %g %g %g\n",
                 (long) cp->lnum, phys[0], phys[1], phys[2]);
  P4EST_LDEBUGF ("Tree %d level %d invert to %g %g %g\n",
                 (int) which_tree, quadrant->level, cp->inv[0], cp->inv[1],
                 cp->inv[2]);

  /* we choose a stricter tolerance on the consumer side to losing points
   * in the local search */
  tol = overlap_p4est_is_meta (p4est) ? SC_1000_EPS : (2 * SC_1000_EPS);

  /* check for quadrant intersection */
  dh = COORDINATE_IROOTLEN * P4EST_QUADRANT_LEN (quadrant->level);
  qxyz[0] = COORDINATE_IROOTLEN * quadrant->x;
  qxyz[1] = COORDINATE_IROOTLEN * quadrant->y;
#ifndef P4_TO_P8
  qxyz[2] = 0.;
  dhz = 1.;
#else
  qxyz[2] = COORDINATE_IROOTLEN * quadrant->z;
  dhz = dh;
#endif
  if ((cp->inv[0] < qxyz[0] - tol || cp->inv[0] > qxyz[0] + dh + tol) ||
      (cp->inv[1] < qxyz[1] - tol || cp->inv[1] > qxyz[1] + dh + tol) ||
      (cp->inv[2] < qxyz[2] - tol || cp->inv[2] > qxyz[2] + dhz + tol)) {
    return 0;
  }

  P4EST_LDEBUGF ("Point %ld survive quadrant\n", (long) cp->lnum);
  return 1;
}

static void
coordinate_get_quadrant_center (p4est_quadrant_t *q, double qxyz[3])
{
  p4est_qcoord_t      h2;

  /* get quadrant center reference coordinates and store them in qxyz */
  h2 = P4EST_QUADRANT_LEN (q->level) >> 1;
  qxyz[0] = COORDINATE_IROOTLEN * (q->x + h2);
  qxyz[1] = COORDINATE_IROOTLEN * (q->y + h2);
#ifndef P4_TO_P8
  qxyz[2] = 0.;
#else
  qxyz[2] = COORDINATE_IROOTLEN * (q->z + h2);
#endif
}

static void
coordinate_get_quadrant_corner (p4est_quadrant_t *q, int cid, double qxyz[3])
{
  p4est_qcoord_t      h, qcoords[3];

  /* get quadrant corner reference coordinates and store them in qxyz */
  h = P4EST_QUADRANT_LEN (q->level);
  qcoords[0] = q->x + ((cid % 2) ? h : 0);
  qcoords[1] = q->y + (((cid % 4) / 2) ? h : 0);
#ifndef P4_TO_P8
  qcoords[2] = 0;
#else
  qcoords[2] = q->z + ((cid / 4) ? h : 0);
#endif
  qxyz[0] = COORDINATE_IROOTLEN * qcoords[0];
  qxyz[1] = COORDINATE_IROOTLEN * qcoords[1];
  qxyz[2] = COORDINATE_IROOTLEN * qcoords[2];
}

/* ---------------------------------------------------------------------- */
///                             Simple Example 
/* ---------------------------------------------------------------------- */
/* Simple example of a overlap_exchange application. It computes artifical
 * solution data on the producer side and queries it for the center of every
 * cell of a consumer p4est by use of constant interpolation. */

typedef struct simple_data
{
  double              myvalue;
  int                 isset;
}
simple_data_t;

typedef struct simple_point
{
  /* data for intersection tests */
  coordinate_point_t  cp;

  /* interpolation related data */
  simple_data_t       data;
}
simple_point_t;

static void
simple_consumer_query_centers_fn (p4est_iter_volume_info_t *info,
                                  void *user_data)
{
  p4est_quadrant_t   *q;
  consumer_t         *c;
  simple_point_t     *sp;
  double              qxyz[3], *phys;

  P4EST_ASSERT (info != NULL && info->p4est != NULL);
  P4EST_ASSERT (info->p4est->user_pointer == user_data);

  /* access quadrant */
  c = (consumer_t *) info->p4est->user_pointer;
  P4EST_ASSERT (c->con4est == info->p4est);
  P4EST_ASSERT (c->lquad_idx >= 0 &&
                c->lquad_idx < c->con4est->local_num_quadrants);
  P4EST_ASSERT (info->quad != NULL);
  q = info->quad;

  /* transform consumer quadrant center to physical using map */
  coordinate_get_quadrant_center (q, qxyz);
  phys = (sp = (simple_point_t *)
          sc_array_index (c->query_xyz, (size_t) c->lquad_idx))->cp.xyz;
  memset (sp, 0, sizeof (simple_point_t));
  c->congeom->X (c->congeom, info->treeid, qxyz, phys);
  sp->cp.lnum = c->lquad_idx++;
  sp->cp.which_tree = -1;
  sp->data.myvalue = 0.;
  sp->data.isset = 0;

  P4EST_LDEBUGF ("Consumer input tree %d level %d quad %g %g %g\n",
                 (int) info->treeid, q->level, qxyz[0], qxyz[1], qxyz[2]);
  P4EST_LDEBUGF ("Consumer point %ld compute %g %g %g\n",
                 (long) sp->cp.lnum, phys[0], phys[1], phys[2]);

  /* optimize: ignore this point if not intersecting producer domain */
}

static double
simple_producer_evaluate (producer_t *p, double pxyz[3])
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

static void
simple_producer_init_quadrants_fn (p4est_iter_volume_info_t *info,
                                   void *user_data)
{
  p4est_quadrant_t   *q;
  simple_data_t      *d;
  producer_t         *p;
  double              qxyz[3], phys[3];

  P4EST_ASSERT (info != NULL && info->p4est != NULL);
  P4EST_ASSERT (info->p4est->user_pointer == user_data);

  p = (producer_t *) info->p4est->user_pointer;
  P4EST_ASSERT (p->pro4est == info->p4est);
  P4EST_ASSERT (info->quad != NULL);
  q = info->quad;
  d = (simple_data_t *) q->p.user_data;
  P4EST_ASSERT (d != NULL);

  /* transform consumer quadrant center to physical using map */
  coordinate_get_quadrant_center (q, qxyz);
  p->progeom->X (p->progeom, info->treeid, qxyz, phys);

  /* interpolate prescribed field at that point */
  d->myvalue = simple_producer_evaluate (p, phys);
  d->isset = 1;

  P4EST_LDEBUGF ("Producer input tree %d level %d quad %g %g %g\n",
                 (int) info->treeid, q->level, qxyz[0], qxyz[1], qxyz[2]);
  P4EST_LDEBUGF ("Producer compute %g %g %g\n", phys[0], phys[1], phys[2]);
}

static void
simple_producer_init_quadrants (producer_t *p)
{
  /* generate a local set of cell values by interpolating a function */
  p->lquad_idx = 0;
  p4est_iterate (p->pro4est, NULL, p, simple_producer_init_quadrants_fn, NULL
#ifdef P4_TO_P8
                 , NULL
#endif
                 , NULL);
}

static void
simple_consumer_query_centers (global_t *g)
{
  consumer_t         *c = g->c;

  /* generate a query point for every local quadrant center */
  c->lquad_idx = 0;
  c->query_xyz =
    sc_array_new_count (sizeof (simple_point_t),
                        c->con4est->local_num_quadrants);

  p4est_iterate (c->con4est, NULL, c, simple_consumer_query_centers_fn, NULL
#ifdef P4_TO_P8
                 , NULL
#endif
                 , NULL);
  P4EST_ASSERT (c->lquad_idx == c->con4est->local_num_quadrants);
}

static int
simple_intersect_fn (p4est_t *p4est, p4est_topidx_t which_tree,
                     p4est_quadrant_t *quadrant, p4est_locidx_t lnum,
                     void *point, void *user)
{
  simple_point_t     *sp;

  P4EST_ASSERT (point != NULL);
  sp = (simple_point_t *) point;
  return coordinate_intersect_fn (p4est, which_tree, quadrant, lnum, &sp->cp,
                                  user);
}

static int
simple_interpolate_fn (p4est_t *p4est, p4est_topidx_t which_tree,
                       p4est_quadrant_t *quadrant, p4est_locidx_t local_num,
                       void *point, void *user)
{
  simple_point_t     *sp;
  simple_data_t      *d;

  P4EST_ASSERT (point != NULL);
  sp = (simple_point_t *) point;
  P4EST_ASSERT (quadrant != NULL);
  P4EST_ASSERT (quadrant->p.user_data != NULL);
  d = (simple_data_t *) quadrant->p.user_data;

  /* apply producer interpolation data to consumer point */
  sp->data.myvalue = d->myvalue;
  sp->data.isset = 1;
  P4EST_LDEBUGF ("Producer point %ld prodata set to %f.\n",
                 (long) sp->cp.lnum, sp->data.myvalue);
  return 1;
}

static void
simple_exchange (global_t *g)
{
  P4EST_ASSERT (g != NULL);
  P4EST_ASSERT (g->p != NULL);
  P4EST_ASSERT (g->p->pro4est != NULL);
  P4EST_ASSERT (g->c != NULL);
  P4EST_ASSERT (g->c->query_xyz != NULL);

  overlap_exchange (g->p->pro4est, g->c->query_xyz, g->c->concomm, g->glocomm,
                    simple_intersect_fn, simple_interpolate_fn, g->usr_ctx);
}

static void
simple_verify (global_t *g)
{
  double              err, err_rel, sol, sol_norm;
  size_t              qi;
  size_t              set_qpoints;
  simple_point_t     *sp;
  consumer_t         *c = g->c;
  producer_t         *p = g->p;

  P4EST_GLOBAL_PRODUCTION ("OVERLAP: result verification\n");

  /* compute absolute error and norm of solution */
  sol_norm = 0.;
  err = 0.;
  set_qpoints = 0;
  for (qi = 0; qi < c->query_xyz->elem_count; qi++) {
    sp = (simple_point_t *) sc_array_index (c->query_xyz, qi);
    if (sp->data.isset) {
      set_qpoints++;
      sol = simple_producer_evaluate (p, sp->cp.xyz);
      sol_norm += sol * sol;
      sol -= sp->data.myvalue;
      err += sol * sol;
    }
  }

  /* compute and check relative error */
  err_rel = sqrt (err / sol_norm);
  if (err_rel > SC_1000_EPS) {
    printf
      ("We have a relative interpolation error of %f on %ld qpoints in the intersection area.\n",
       err_rel, set_qpoints);
  }
}

/* ---------------------------------------------------------------------- */
///                         VTK-Output of Simple Example 
/* ---------------------------------------------------------------------- */
/* Visualiziation of the results of the overlap_exchange. */

static void
simple_consumer_print_interpolation_data (consumer_t *c)
{
  simple_point_t     *sp;
  p4est_locidx_t      cind;

  for (cind = 0; cind < c->con4est->local_num_quadrants; cind++) {
    sp = (simple_point_t *) sc_array_index (c->query_xyz, cind);
    printf
      ("%d: Query point [%f,%f,%f] in tree %d obtained interpolated data %f.\n",
       c->conrank, sp->cp.xyz[0], sp->cp.xyz[1], sp->cp.xyz[2],
       sp->cp.which_tree, sp->data.myvalue);
  }
}

static void
simple_producer_extract_vtk_fn (p4est_iter_volume_info_t *info,
                                void *user_data)
{
  p4est_quadrant_t   *q;
  simple_data_t      *d;
  producer_t         *p;

  P4EST_ASSERT (info != NULL && info->p4est != NULL);
  P4EST_ASSERT (info->p4est->user_pointer == user_data);
  p = (producer_t *) info->p4est->user_pointer;
  P4EST_ASSERT (p->pro4est == info->p4est);
  P4EST_ASSERT (info->quad != NULL);
  q = info->quad;
  d = (simple_data_t *) q->p.user_data;
  P4EST_ASSERT (d != NULL);
  P4EST_ASSERT (d->isset = 1);

  *(double *) sc_array_index (p->interpolation_data, p->lquad_idx++) =
    d->myvalue;
}

/* write consumer p4est with interpolation data into vtk */
static void
simple_consumer_write_vtk (consumer_t *c)
{
  int                 retval;
  p4est_vtk_context_t *con_context;

  /* allocate context and set parameters */
  con_context =
    p4est_vtk_context_new (c->con4est, P4EST_STRING "_consumer_new");
  p4est_vtk_context_set_geom (con_context, c->congeom);
  p4est_vtk_context_set_continuous (con_context, 1);

  /* write header */
  con_context = p4est_vtk_write_header (con_context);
  SC_CHECK_ABORT (con_context != NULL,
                  P4EST_STRING "_vtk: Error writing header");

  /* write cell_data */
  con_context =
    p4est_vtk_write_cell_dataf (con_context, 1, 1, 1, 0, 2, 1,
                                "interpolation", c->interpolation_data,
                                "is_set", c->isset_data,
                                "xyz", c->xyz_data, con_context);
  SC_CHECK_ABORT (con_context != NULL,
                  P4EST_STRING "_vtk: Error writing celldata");

  /* properly write rest of the files' contents */
  retval = p4est_vtk_write_footer (con_context);
  SC_CHECK_ABORT (!retval, P4EST_STRING "_vtk: Error writing footer");
}

/* write producer p4est with interpolation data into vtk */
static void
simple_producer_write_vtk (producer_t *p)
{
  int                 retval;
  p4est_vtk_context_t *pro_context;

  /* allocate context and set parameters */
  pro_context =
    p4est_vtk_context_new (p->pro4est, P4EST_STRING "_producer_new");
  p4est_vtk_context_set_geom (pro_context, p->progeom);
  p4est_vtk_context_set_continuous (pro_context, 1);

  /* write header */
  pro_context = p4est_vtk_write_header (pro_context);
  SC_CHECK_ABORT (pro_context != NULL,
                  P4EST_STRING "_vtk: Error writing header");

  /* write cell_data */
  pro_context =
    p4est_vtk_write_cell_dataf (pro_context, 1, 1, 1, 0, 1, 0,
                                "interpolation", p->interpolation_data,
                                pro_context);
  SC_CHECK_ABORT (pro_context != NULL,
                  P4EST_STRING "_vtk: Error writing celldata");

  /* properly write rest of the files' contents */
  retval = p4est_vtk_write_footer (pro_context);
  SC_CHECK_ABORT (!retval, P4EST_STRING "_vtk: Error writing footer");
}

/* Output the results of the last exchange.
 * If text is true, iterates over all query points in g->c->query_xyz and prints
 * their data in a convenient format.
 * If vtk is true, vtk output of the consumer mesh (including interpolated data)
 * and the producer mesh (including quadrant user_data) is created.
 * Currently only works for one query point per quadrant. */
static void
simple_output_results (global_t *g, int text, int vtk)
{
  producer_t         *p = g->p;
  consumer_t         *c = g->c;
  simple_point_t     *sp;
  size_t              cind;
  size_t              plnq, clnq;

  /* output the interpolation data of all query points */
  if (text) {
    simple_consumer_print_interpolation_data (c);
  }

  if (vtk) {
    /* consumer side vtk output */
    /* allocate arrays needed for consumer side output */
    clnq = c->con4est->local_num_quadrants;
    c->interpolation_data = sc_array_new_count (sizeof (double), clnq);
    c->xyz_data = sc_array_new_count (sizeof (double), 3 * clnq);
    c->isset_data = sc_array_new_count (sizeof (double), clnq);

    /* extract interpolated data from query point array */
    P4EST_ASSERT (c->query_xyz != NULL);
    P4EST_ASSERT (c->query_xyz->elem_count == clnq);

    for (cind = 0; cind < (size_t) c->con4est->local_num_quadrants; cind++) {
      sp = (simple_point_t *) sc_array_index (c->query_xyz, cind);
      *(double *) sc_array_index (c->interpolation_data, cind) =
        sp->data.myvalue;
      *(double *) sc_array_index (c->isset_data, cind) =
        (double) sp->data.isset;
      *(double *) sc_array_index (c->xyz_data, 3 * cind) = sp->cp.xyz[0];
      *(double *) sc_array_index (c->xyz_data, 3 * cind + 1) = sp->cp.xyz[1];
      *(double *) sc_array_index (c->xyz_data, 3 * cind + 2) = sp->cp.xyz[2];
    }

    /* write vtk output files */
    simple_consumer_write_vtk (c);

    /* destroy vtk specific arrays */
    sc_array_destroy (c->interpolation_data);
    sc_array_destroy (c->isset_data);
    sc_array_destroy (c->xyz_data);

    /* producer side vtk output */
    /* allocate arrays needed for producer side output */
    plnq = p->pro4est->local_num_quadrants;
    p->interpolation_data = sc_array_new_count (sizeof (double), plnq);

    /* extract producer data from p4est */
    p->lquad_idx = 0;
    p4est_iterate (p->pro4est, NULL, p, simple_producer_extract_vtk_fn,
#ifdef P4_TO_P8
                   NULL,
#endif
                   NULL, NULL);

    /* write vtk output files */
    simple_producer_write_vtk (p);

    /* destroy vtk specific arrays */
    sc_array_destroy (p->interpolation_data);
  }
}

/* ---------------------------------------------------------------------- */
///                          Refinement Methods
/* ---------------------------------------------------------------------- */
/* Several refinement methods for the general example setup. */

typedef struct refine_context
{
  p4est_geometry_t   *geom;
  int                 maxlevel;
  double              geom_radius;
  sc_array_t         *polygon;
}
refine_ctx_t;

static int
refine_geometrical_fn (p4est_t *p4est, p4est_topidx_t which_tree,
                       p4est_quadrant_t *quadrant)
{
  refine_ctx_t       *r;
  double              qxyz[3];
  double              phys[3] = { 0, 0, 0 };
  double              dist;
  P4EST_ASSERT (p4est != NULL);
  P4EST_ASSERT (p4est->user_pointer != NULL);
  r = (refine_ctx_t *) p4est->user_pointer;
  P4EST_ASSERT (r->geom != NULL);
  P4EST_ASSERT (r->geom->X != NULL);

  /* transform producer quadrant center to physical using map */
  coordinate_get_quadrant_center (quadrant, qxyz);
  r->geom->X (r->geom, which_tree, qxyz, phys);

  /* compute distance from point of interest */
  phys[0] -= 0.3;
  phys[1] -= 0.5;
  phys[2] -= 0.4;
  dist = sqrt (phys[0] * phys[0] + phys[1] * phys[1] + phys[2] * phys[2]);

  /* refine quadrants that are close enough to point of interest */
  if (quadrant->level < r->maxlevel - floor (dist / r->geom_radius)) {
    return 1;
  }

  return 0;
}

static int
refine_childid_fn (p4est_t *p4est, p4est_topidx_t which_tree,
                   p4est_quadrant_t *quadrant)
{
  refine_ctx_t       *r;

  P4EST_ASSERT (p4est != NULL);
  P4EST_ASSERT (p4est->user_pointer != NULL);
  r = (refine_ctx_t *) p4est->user_pointer;

  if ((int) quadrant->level >= (r->maxlevel - (int) (which_tree % 3))) {
    return 0;
  }
  if (quadrant->level == 1 && p4est_quadrant_child_id (quadrant) == 3) {
    return 1;
  }
  if (quadrant->x == P4EST_LAST_OFFSET (2) &&
      quadrant->y == P4EST_LAST_OFFSET (2)) {
    return 1;
  }
  if (quadrant->x >= P4EST_QUADRANT_LEN (2)) {
    return 0;
  }

  return 1;
}

static int
refine_rank_fn (p4est_t *p4est, p4est_topidx_t which_tree,
                p4est_quadrant_t *quadrant)
{
  refine_ctx_t       *r;

  P4EST_ASSERT (p4est != NULL);
  P4EST_ASSERT (p4est->user_pointer != NULL);
  r = (refine_ctx_t *) p4est->user_pointer;

  if ((int) quadrant->level >= r->maxlevel) {
    return 0;
  }
  if (p4est->mpirank >= 5 && p4est->mpirank <= 15) {
    return 1;
  }

  return 0;
}

typedef struct overlap_polygon_normal
{
  double              normal[2];
  double              prod;
}
overlap_polygon_normal_t;

static sc_array_t  *
refine_get_polygon_context (global_t *g, double *xcoords,
                            double *ycoords, int ncoords)
{
  sc_array_t         *polygon;
  int                 i;
  double              edge[2], nprod;
  overlap_polygon_normal_t *opn;

  polygon = sc_array_new_count (sizeof (overlap_polygon_normal_t), ncoords);

  for (i = 0; i < ncoords; i++) {
    /* compute edge of the polygon */
    edge[0] = xcoords[(i + 1) % ncoords] - xcoords[i];
    edge[1] = ycoords[(i + 1) % ncoords] - ycoords[i];

    /* compute normal and scalar product of normal with edge */
    opn = (overlap_polygon_normal_t *) sc_array_index_int (polygon, i);
    opn->normal[0] = edge[1];
    opn->normal[1] = -edge[0];
    opn->prod = opn->normal[0] * xcoords[i] + opn->normal[1] * ycoords[i];
    P4EST_ASSERT (fabs
                  (opn->normal[0] * xcoords[(i + 1) % ncoords] +
                   opn->normal[1] * ycoords[(i + 1) % ncoords] - opn->prod) <
                  SC_1000_EPS);

    /* make sure that the polygon lies in the upper half space */
    nprod =
      opn->normal[0] * xcoords[(i + 2) % ncoords] +
      opn->normal[1] * ycoords[(i + 2) % ncoords];
    if (nprod < opn->prod) {
      opn->normal[0] = -opn->normal[0];
      opn->normal[1] = -opn->normal[1];
      opn->prod = -opn->prod;
    }
    P4EST_ASSERT (opn->normal[0] * xcoords[(i + 2) % ncoords] +
                  opn->normal[1] * ycoords[(i + 2) % ncoords] >= opn->prod);
    P4EST_ASSERT (fabs (edge[0] * opn->normal[0] + edge[1] * opn->normal[1]) <
                  SC_1000_EPS);
  }
  return polygon;
}

static int
refine_polygon_fn (p4est_t *p4est, p4est_topidx_t which_tree,
                   p4est_quadrant_t *quadrant)
{
  refine_ctx_t       *r;
  sc_array_t         *polygon;
  overlap_polygon_normal_t *opn;
  double              qxyz[3];
  double              phys[3] = { 0, 0, 0 };
  int                 cid, corners_inside;
  size_t              iz;

  P4EST_ASSERT (p4est != NULL);
  P4EST_ASSERT (p4est->user_pointer != NULL);
  r = (refine_ctx_t *) p4est->user_pointer;
  P4EST_ASSERT (r->geom != NULL);
  P4EST_ASSERT (r->geom->X != NULL);
  P4EST_ASSERT (r->polygon != NULL);
  polygon = (sc_array_t *) r->polygon;
  P4EST_ASSERT (polygon->elem_size == sizeof (overlap_polygon_normal_t));

  if (quadrant->level >= r->maxlevel) {
    /* we already refined up to the maximum allowed */
    return 0;
  }

  /* check if any corner lies inside the polygon */
  corners_inside = 0;
  for (cid = 0; cid < P4EST_CHILDREN; cid++) {
    /* transform producer quadrant corner to physical using map */
    coordinate_get_quadrant_corner (quadrant, cid, qxyz);
    r->geom->X (r->geom, which_tree, qxyz, phys);

    for (iz = 0; iz < polygon->elem_count; iz++) {
      opn = (overlap_polygon_normal_t *) sc_array_index (polygon, iz);
      if (opn->normal[0] * phys[0] + opn->normal[1] * phys[1] < opn->prod) {
        /* the point lies in the wrong half space */
        break;
      }
      if (iz == polygon->elem_count - 1) {
        /* the corner passed all tests, so it lies inside the polygon */
        corners_inside++;
      }
    }
  }
  if (corners_inside > 0 && corners_inside < P4EST_CHILDREN) {
    /* the quadrant intersects the polygon boundary, refine up to refine_level */
    return 1;
  }

  return 0;
}

/* ---------------------------------------------------------------------- */
///                          Adaptive Refinement
/* ---------------------------------------------------------------------- */
/* Adaptive refinement method based on overlap_exchange. Every consumer query
 * point is discretized by a 3x3(x3) tensor of query points. These are used to
 * mark and refine the boundary of the intersection area of the meshes on both
 * the consumer as well as the producer side. */

#ifndef P4_TO_P8
#define ADAPTIVE_NUM_TENSOR_POINTS 9
#else
#define ADAPTIVE_NUM_TENSOR_POINTS 27
#endif

typedef struct adaptive_point
{
  /* data for intersection tests */
  coordinate_point_t  cp;

  /* adaptive overlap refinement related data */
  int                 isboundary;
  int                 isset;
}
adaptive_point_t;

typedef struct adaptive_data
{
  int                 refine;
  int                 isboundary;
}
adaptive_data_t;

/* set isboundary flag of a given producer quadrant */
static void
adaptive_producer_init_quadrants_fn (p4est_iter_volume_info_t *info,
                                     void *user_data)
{
  p4est_quadrant_t   *q;
  adaptive_data_t    *d;
  producer_t         *p;
  p4est_qcoord_t      h;
  p4est_topidx_t     *ttt, tid;

  P4EST_ASSERT (info != NULL && info->p4est != NULL
                && info->p4est->user_pointer != NULL);
  P4EST_ASSERT (info->p4est->user_pointer == user_data);
  p = (producer_t *) info->p4est->user_pointer;
  P4EST_ASSERT (info->quad != NULL);
  q = info->quad;
  P4EST_ASSERT (q->p.user_data != NULL);
  d = (adaptive_data_t *) q->p.user_data;

  /* check if quadrant lies on the domain boundary */
  memset (d, 0, sizeof (adaptive_data_t));
  h = P4EST_QUADRANT_LEN (q->level);
  ttt = p->proconn->tree_to_tree;
  tid = info->treeid;
  if ((q->x == 0 && ttt[tid * P4EST_FACES] == tid)
      || (q->x + h == P4EST_ROOT_LEN && ttt[tid * P4EST_FACES + 1] == tid)
      || (q->y == 0 && ttt[tid * P4EST_FACES + 2] == tid)
      || (q->y + h == P4EST_ROOT_LEN && ttt[tid * P4EST_FACES + 3] == tid)
#ifdef P4_TO_P8
      || (q->z == 0 && ttt[tid * P4EST_FACES + 4] == tid)
      || (q->z + h == P4EST_ROOT_LEN && ttt[tid * P4EST_FACES + 5] == tid)
#endif
    ) {
    d->isboundary = 1;
  }
}

/* set isboundary flag of the producer quadrants */
static void
adaptive_producer_init_quadrants (producer_t *p)
{
  p4est_iterate (p->pro4est, NULL, p, adaptive_producer_init_quadrants_fn,
                 NULL
#ifdef P4_TO_P8
                 , NULL
#endif
                 , NULL);

}

/* create a 3x3(x3) tensor of query points for a given quadrant and set its
 * isboundary flag */
static void
adaptive_consumer_query_tensors_fn (p4est_iter_volume_info_t *info,
                                    void *user_data)
{
  p4est_quadrant_t   *q;
  adaptive_data_t    *d;
  consumer_t         *c;
  adaptive_point_t   *ap;
  p4est_qcoord_t      h, hhalf;
  p4est_topidx_t     *ttt, tid;
  double              qxyz[3], *phys;
  int                 i, j;
#ifdef P4_TO_P8
  int                 k;
#endif

  P4EST_ASSERT (info != NULL && info->p4est != NULL);
  P4EST_ASSERT (info->p4est->user_pointer == user_data);

  /* access quadrant */
  c = (consumer_t *) info->p4est->user_pointer;
  P4EST_ASSERT (c->con4est == info->p4est);
  P4EST_ASSERT (c->lquad_idx >= 0 &&
                c->lquad_idx <
                ADAPTIVE_NUM_TENSOR_POINTS * c->con4est->local_num_quadrants);
  P4EST_ASSERT (info->quad != NULL);
  q = info->quad;
  P4EST_ASSERT (q->p.user_data != NULL);
  d = (adaptive_data_t *) q->p.user_data;

  memset (d, 0, sizeof (adaptive_data_t));
  h = P4EST_QUADRANT_LEN (q->level);
  ttt = c->conconn->tree_to_tree;
  tid = info->treeid;
  if ((q->x == 0 && ttt[tid * P4EST_FACES] == tid)
      || (q->x + h == P4EST_ROOT_LEN && ttt[tid * P4EST_FACES + 1] == tid)
      || (q->y == 0 && ttt[tid * P4EST_FACES + 2] == tid)
      || (q->y + h == P4EST_ROOT_LEN && ttt[tid * P4EST_FACES + 3] == tid)
#ifdef P4_TO_P8
      || (q->z == 0 && ttt[tid * P4EST_FACES + 4] == tid)
      || (q->z + h == P4EST_ROOT_LEN && ttt[tid * P4EST_FACES + 5] == tid)
#endif
    ) {
    d->isboundary = 1;
  }

  /* create query points according to a 3x3x3 tensor product */
  hhalf = P4EST_QUADRANT_LEN (q->level + 1);    /* half of quadrant length */
  for (i = -1; i <= 1; i++) {
    for (j = -1; j <= 1; j++) {
#ifdef P4_TO_P8
      for (k = -1; k <= 1; k++) {
#endif
        /* When overlapping two identical meshes, almost all query points of a
         * 3x3 tensor would be found on the boundary of a producer quadrant.
         * Since the exchange allows for only one match, the refinement markers
         * depend on the boundary-status of the first quadrant finding the point.
         * This results in an assymetric refinement scheme for symmetric meshes.
         * To avoid this, we move all surface query points inside their quadrant
         * by a step of size 0.25 * COORDINATE_IROOTLEN, so a quarter of the
         * minimal quadrant size. */
        qxyz[0] = q->x + (1 + i) * (double) hhalf - (double) i * 0.25;
        qxyz[1] = q->y + (1 + j) * (double) hhalf - (double) j * 0.25;
#ifndef P4_TO_P8
        qxyz[2] = 0.;
#else
        qxyz[2] = q->z + (1 + k) * (double) hhalf - (double) k * 0.25;
#endif
        qxyz[0] *= COORDINATE_IROOTLEN;
        qxyz[1] *= COORDINATE_IROOTLEN;
        qxyz[2] *= COORDINATE_IROOTLEN;

        phys = (ap = (adaptive_point_t *)
                sc_array_index (c->query_xyz, (size_t) c->lquad_idx))->cp.xyz;
        memset (ap, 0, sizeof (adaptive_point_t));
        c->congeom->X (c->congeom, info->treeid, qxyz, phys);

        ap->cp.lnum = c->lquad_idx++;
        ap->cp.which_tree = -1;
        ap->isboundary = d->isboundary;
        ap->isset = 0;
#ifdef P4_TO_P8
      }
#endif
    }
  }
}

/* create a 3x3(x3) tensor of query points for all quadrants and set their
 * isboundary flag */
static void
adaptive_consumer_query_tensors (global_t *g)
{
  consumer_t         *c = g->c;

  /* generate a query point for every local quadrant center */
  c->lquad_idx = 0;
  c->query_xyz = sc_array_new_count (sizeof (adaptive_point_t),
                                     ADAPTIVE_NUM_TENSOR_POINTS *
                                     c->con4est->local_num_quadrants);
  p4est_iterate (c->con4est, NULL, c, adaptive_consumer_query_tensors_fn, NULL
#ifdef P4_TO_P8
                 , NULL
#endif
                 , NULL);
  P4EST_ASSERT (c->lquad_idx ==
                ADAPTIVE_NUM_TENSOR_POINTS * c->con4est->local_num_quadrants);
}

static int
adaptive_intersect_fn (p4est_t *p4est, p4est_topidx_t which_tree,
                       p4est_quadrant_t *quadrant, p4est_locidx_t lnum,
                       void *point, void *user)
{
  adaptive_point_t   *ap;

  P4EST_ASSERT (point != NULL);
  ap = (adaptive_point_t *) point;
  return coordinate_intersect_fn (p4est, which_tree, quadrant, lnum, &ap->cp,
                                  user);
}

static int
adaptive_interpolate_fn (p4est_t *p4est, p4est_topidx_t which_tree,
                         p4est_quadrant_t *quadrant,
                         p4est_locidx_t local_num, void *point, void *user)
{
  adaptive_point_t   *ap;
  adaptive_data_t    *d;

  P4EST_ASSERT (point != NULL);
  ap = (adaptive_point_t *) point;
  P4EST_ASSERT (quadrant != NULL);
  P4EST_ASSERT (quadrant->p.user_data != NULL);
  d = (adaptive_data_t *) quadrant->p.user_data;

  ap->isset = 1;
  if ((ap->isboundary == 1 || d->isboundary == 1)) {
    /* since a leaf intersects, we are in the intersection area
     *   ap->isboundary == 1 => we are on the consumer mesh boundary
     *   d->isboundary == 1 => we are on the producer mesh boundary */
    d->refine = 1;
    if (d->isboundary == 1) {
      ap->isset = 2;
    }
  }

  return 1;
}

/* set a given quadrants refine flag based on its updated query points */
static void
adaptive_consumer_evaluate_tensors_fn (p4est_iter_volume_info_t *info,
                                       void *user_data)
{
  p4est_quadrant_t   *q;
  consumer_t         *c;
  adaptive_point_t   *ap;
  int                 i, npin, npout;
  adaptive_data_t    *d;

  P4EST_ASSERT (info != NULL && info->p4est != NULL);
  P4EST_ASSERT (info->p4est->user_pointer == user_data);

  /* access quadrant */
  c = (consumer_t *) info->p4est->user_pointer;
  P4EST_ASSERT (c->con4est == info->p4est);
  P4EST_ASSERT (c->lquad_idx >= 0 &&
                c->lquad_idx <
                ADAPTIVE_NUM_TENSOR_POINTS * c->con4est->local_num_quadrants);
  P4EST_ASSERT (info->quad != NULL);
  q = info->quad;

  /* iterate over all children */
  npin = npout = 0;
  for (i = 0; i < ADAPTIVE_NUM_TENSOR_POINTS; i++) {
    ap = (adaptive_point_t *) sc_array_index (c->query_xyz, c->lquad_idx++);
    P4EST_ASSERT (ap->isset == 0 || ap->isset == 1 || ap->isset == 2);
    if (ap->isset) {
      npin++;
    }
    if (ap->isset == 2) {
      npout++;
    }
  }

  d = (adaptive_data_t *) q->p.user_data;
  if ((npin && d->isboundary == 1) || npout) {
    /* npin && d->isboundary == 1 => on consumer boundary + in producer mesh
     * npout                      => in consumer mesh + on producer boundary */
    d->refine = 1;
  }
  else {
    d->refine = 0;
  }
}

/* set the quadrants refine flags based on the updated query points */
static void
adaptive_consumer_evaluate_tensors (global_t *g)
{
  consumer_t         *c = g->c;

  /* generate a query point for every local quadrant center */
  c->lquad_idx = 0;
  p4est_iterate (c->con4est, NULL, c, adaptive_consumer_evaluate_tensors_fn,
                 NULL
#ifdef P4_TO_P8
                 , NULL
#endif
                 , NULL);
  P4EST_ASSERT (c->lquad_idx ==
                ADAPTIVE_NUM_TENSOR_POINTS * c->con4est->local_num_quadrants);
}

static int
adaptive_refine_fn (p4est_t *p4est, p4est_topidx_t which_tree,
                    p4est_quadrant_t *quadrant)
{
  refine_ctx_t       *r;
  adaptive_data_t    *d = (adaptive_data_t *) quadrant->p.user_data;

  P4EST_ASSERT (p4est != NULL);
  P4EST_ASSERT (p4est->user_pointer != NULL);
  r = (refine_ctx_t *) p4est->user_pointer;

  if (d->refine == 1 && quadrant->level < r->maxlevel) {
    return 1;                   /* the producer quadrant contains a boundary consumer point */
  }

  return 0;
}

/* ---------------------------------------------------------------------- */
///                          Example Applications
/* ---------------------------------------------------------------------- */
/* General example workflow creating a consumer and a producer mesh, refining
 * them and performing a simple overlap before evaluating the results. */

static void
apps_init (global_t *g, sc_MPI_Comm mpicomm)
{
  producer_t         *p = g->p = &g->pro;
  consumer_t         *c = g->c = &g->con;
  int                 mpiret;
  int                 glorank, glosize;
  p4est_connectivity_t *conns[2];
  refine_ctx_t        consumer_refine_context, *cref_ctx =
    &consumer_refine_context;
  refine_ctx_t        producer_refine_context, *pref_ctx =
    &producer_refine_context;

  /* initialization of communicators */
  g->glocomm = mpicomm;
  mpiret = sc_MPI_Comm_rank (g->glocomm, &glorank);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_size (g->glocomm, &glosize);
  SC_CHECK_MPI (mpiret);

  sc_MPI_Comm         concomm, procomm;
  mpiret =
    sc_MPI_Comm_split (g->glocomm,
                       (glorank < g->conprocs) ? 0 : sc_MPI_UNDEFINED,
                       glorank, &concomm);
  mpiret =
    sc_MPI_Comm_split (g->glocomm,
                       (glorank >= glosize - g->proprocs) ?
                       0 : sc_MPI_UNDEFINED, glorank, &procomm);

  if (glorank < g->conprocs) {
    sc_MPI_Comm_free (&concomm);
  }
  if (glorank >= glosize - g->proprocs) {
    sc_MPI_Comm_free (&procomm);
  }

  /* Create two connectivities. They will be assigned to the producer and the
   * consumer based on the value of g->example. */
  if (g->example <= 1) {
    nbricks[0] = 10;
    nbricks[1] = 2;
#ifdef P4_TO_P8
    nbricks[2] = 2;
#else
    nbricks[2] = 1;
#endif
    nbricks[3] = 3;
    nbricks[4] = 2;
    nbricks[5] = 1;
  }
  else {
    nbricks[0] = 1;
    nbricks[1] = 1;
    nbricks[2] = 1;
    nbricks[3] = 1;
    nbricks[4] = 1;
    nbricks[5] = 1;
  }
  conns[0] = p4est_connectivity_new_brick (nbricks[0], nbricks[1]
#ifdef P4_TO_P8
                                           , nbricks[2]
#endif
                                           , 0, 0
#ifdef P4_TO_P8
                                           , 0
#endif
    );
  conns[1] = p4est_connectivity_new_brick (nbricks[3], nbricks[4]
#ifdef P4_TO_P8
                                           , nbricks[5]
#endif
                                           , 0, 0
#ifdef P4_TO_P8
                                           , 0
#endif
    );

  /***************************** PRODUCER ****************************/

  P4EST_GLOBAL_PRODUCTION ("OVERLAP: init producer\n");

  /* setup producer geometry */
  p->progeom = &p->producer_geometry;
  p->progeom->name = "producer";
  p->progeom->destroy = (p4est_geometry_destroy_t) 0;

  /* setup producer communicator */
  p->glocomm = g->glocomm;
  mpiret = sc_MPI_Comm_dup (g->glocomm, &p->procomm);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (p->procomm, &p->prorank);
  SC_CHECK_MPI (mpiret);

  /* assign the geometry depending on the value of g->example */
  if (g->example == 0) {
    p->progeom->user = p->proconn = conns[1];
    p->progeom->X = coordinate_brick_map;
    g->usr_ctx = (void *) coordinate_brick_invmap;
  }
  else if (g->example == 1) {
    p->progeom->user = p->proconn = conns[0];
    p->progeom->X = coordinate_curved_map;
    g->usr_ctx = (void *) coordinate_curved_invmap;
  }
  else {
    P4EST_ASSERT ((g->example == 2) || (g->example == 3));
    p->progeom->user = p->proconn = conns[0];
    p->progeom->X = coordinate_producer_unit_map;
    g->usr_ctx = (void *) coordinate_producer_unit_invmap;
  }

  /* setup producer mesh */
  p->pro4est = p4est_new_ext (p->procomm, p->proconn, 0, p->pminl, 1,
                              sizeof (simple_data_t), NULL, p);

  /***************************** CONSUMER ****************************/

  P4EST_GLOBAL_PRODUCTION ("OVERLAP: init consumer\n");

  /* setup consumer geometry */
  c->congeom = &c->consumer_geometry;
  c->congeom->name = "consumer";
  c->congeom->destroy = (p4est_geometry_destroy_t) 0;

  /* setup consumer communicator */
  c->glocomm = g->glocomm;
  mpiret = sc_MPI_Comm_dup (g->glocomm, &c->concomm);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (c->concomm, &c->conrank);
  SC_CHECK_MPI (mpiret);

  /* assign the geometry that was not assigned to the producer side */
  if (g->example == 0) {
    c->congeom->user = c->conconn = conns[0];
    c->congeom->X = coordinate_curved_map;
  }
  else if (g->example == 1) {
    c->congeom->user = c->conconn = conns[1];
    c->congeom->X = coordinate_brick_map;
  }
  else if (g->example == 2) {
    c->congeom->user = c->conconn = conns[1];
    c->congeom->X = coordinate_consumer_unit_map;
  }
  else {
    /* we want to create a duplicate of the producer connectivity */
    P4EST_ASSERT ((nbricks[0] == nbricks[3]) && (nbricks[1] == nbricks[4])
                  && (nbricks[2] == nbricks[5]));
    c->congeom->user = c->conconn = conns[1];
    c->congeom->X = coordinate_producer_unit_map;
  }

  /* setup consumer mesh */
  c->con4est = p4est_new_ext (c->concomm, c->conconn, 0, c->cminl, 1,
                              0, NULL, c);

  /**************************** REFINEMENT ***************************/
  /* initialize consumer and producer refinement context */
  memset (cref_ctx, 0, sizeof (refine_ctx_t));
  memset (pref_ctx, 0, sizeof (refine_ctx_t));
  cref_ctx->maxlevel = pref_ctx->maxlevel = g->refinement_maxlevel;
  cref_ctx->geom = c->congeom;
  pref_ctx->geom = p->progeom;
  c->con4est->user_pointer = cref_ctx;
  p->pro4est->user_pointer = pref_ctx;

  if (g->refinement_method == 0) {
    /* Adaptively refine the boundary of the mesh intersection area.
     * We query all corners of all consumer quadrants and refine all quadrants,
     * that contains at least one point that was found in the exchange and at
     * least one that was not found.
     * We mark all producer quadrants that contain a boundary query point for
     * refinement. */
    p4est_locidx_t      old_num_proquads, old_num_consquads;

    /* we refine the brick to a higher level than the arch to make up
     * for the difference in the total tree count (6 vs. 40 in 3D example 1) */
    if (g->example == 0) {
      pref_ctx->maxlevel++;
    }
    else if (g->example == 1) {
      cref_ctx->maxlevel++;
    }

    /* prepare p4est to store quadrant data */
    p4est_reset_data (p->pro4est, sizeof (adaptive_data_t), NULL, p);
    p4est_reset_data (c->con4est, sizeof (adaptive_data_t), NULL, c);

    /* compute the maximum numbers of refinements to stay below refine_level */
    old_num_proquads = -1;
    old_num_consquads = -1;
    while (old_num_proquads != p->pro4est->global_num_quadrants ||
           old_num_consquads != c->con4est->global_num_quadrants) {
      old_num_proquads = p->pro4est->global_num_quadrants;
      old_num_consquads = c->con4est->global_num_quadrants;

      /* overlap_exchange may not touch all quadrants */
      adaptive_producer_init_quadrants (p);
      adaptive_consumer_query_tensors (g);

      /* query consumer corners and set p->refine_quadrant during the process */
      overlap_exchange (p->pro4est, c->query_xyz, c->concomm, g->glocomm,
                        adaptive_intersect_fn, adaptive_interpolate_fn,
                        g->usr_ctx);

      /* evaluate which consumer quadrants have to be refined */
      adaptive_consumer_evaluate_tensors (g);

      /* actual refinement based on the exchange results */
      c->con4est->user_pointer = cref_ctx;
      p->pro4est->user_pointer = pref_ctx;
      p4est_refine (p->pro4est, 0, adaptive_refine_fn, NULL);
      p4est_refine (c->con4est, 0, adaptive_refine_fn, NULL);
      c->con4est->user_pointer = c;     /* reset user-pointers for next exchange */
      p->pro4est->user_pointer = p;

      /* cleanup */
      sc_array_destroy (g->c->query_xyz);
    }

    p4est_reset_data (p->pro4est, sizeof (simple_data_t), NULL, p);
    p4est_reset_data (c->con4est, 0, NULL, c);
  }
  else if (g->refinement_method == 1) {
    pref_ctx->geom_radius = cref_ctx->geom_radius = 0.2;
    /* refine arch less than brick, to obtain similar quadrant counts */
    if (g->example == 0) {
      cref_ctx->geom_radius = 0.1;
    }
    else if (g->example == 1) {
      pref_ctx->geom_radius = 0.1;
    }
    p4est_refine (p->pro4est, 1, refine_geometrical_fn, NULL);
    p4est_refine (c->con4est, 1, refine_geometrical_fn, NULL);
  }
  else if (g->refinement_method == 2) {
    p4est_refine (p->pro4est, 1, refine_childid_fn, NULL);
    p4est_refine (c->con4est, 1, refine_rank_fn, NULL);
  }
  else {
    /* refine producer and consumer mesh inside a convex polygon */
    p4est_locidx_t      old_num_quads;
    double              xcoords[5] = { 0.25, 0.6, 0.8, 0.6, 0.3 };
    double              ycoords[5] = { 0.25, 0.1, 0.55, 0.8, 0.6 };
    cref_ctx->polygon = pref_ctx->polygon =
      refine_get_polygon_context (g, xcoords, ycoords, 5);

    /* refinement inside the polygon */
    old_num_quads = -1;
    while (old_num_quads != p->pro4est->global_num_quadrants) {
      old_num_quads = p->pro4est->global_num_quadrants;
      p4est_refine (p->pro4est, 0, refine_polygon_fn, NULL);
      p4est_balance (p->pro4est, P4EST_CONNECT_FACE, NULL);
    }

    old_num_quads = -1;
    while (old_num_quads != c->con4est->global_num_quadrants) {
      old_num_quads = c->con4est->global_num_quadrants;
      p4est_refine (c->con4est, 0, refine_polygon_fn, NULL);
      p4est_balance (c->con4est, P4EST_CONNECT_FACE, NULL);
    }

    /* delete polygon context */
    sc_array_destroy (pref_ctx->polygon);
  }
  p->pro4est->user_pointer = p;
  c->con4est->user_pointer = c;

  p4est_partition (p->pro4est, 0, NULL);
  p4est_partition (c->con4est, 0, NULL);

  P4EST_GLOBAL_PRODUCTION ("OVERLAP: init done\n");
}

static void
apps_run (global_t *g)
{
  /* prepare consumer and producer for exchange */
  simple_consumer_query_centers (g);
  simple_producer_init_quadrants (g->p);

  /*run the actual exchange */
  simple_exchange (g);

  /* evaluate the results of the exchange and cleanup */
  simple_output_results (g, g->output_text, g->output_vtk);
  simple_verify (g);
  sc_array_destroy (g->c->query_xyz);
}

static void
apps_reset (global_t *g)
{
  producer_t         *p = g->p;
  consumer_t         *c = g->c;
  int                 mpiret;

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

static int
usagerr (sc_options_t *opt, const char *msg)
{
  SC_GLOBAL_LERRORF ("Usage required: %s\n", msg);
  return 1;
}

int
main (int argc, char **argv)
{
  int                 mpiret;
  int                 first_argc;
  int                 ue;
  sc_MPI_Comm         mpicomm;
  int                 mpisize;
  sc_options_t       *opt;
  global_t global    , *g = &global;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  mpicomm = sc_MPI_COMM_WORLD;
  sc_init (mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
#ifndef P4EST_ENABLE_DEBUG
  sc_set_log_defaults (NULL, NULL, SC_LP_STATISTICS);
#endif
  p4est_init (NULL, SC_LP_DEFAULT);

  mpiret = sc_MPI_Comm_size (mpicomm, &mpisize);
  SC_CHECK_MPI (mpiret);

  /* process command line arguments */
  opt = sc_options_new (argv[0]);
  /* consumer communicator will consist of processes [0,a)
   * producer communicator will consist of processes [mpisize-b, mpisize) */
  sc_options_add_int (opt, 'a', "cons_processes", &g->conprocs, mpisize,
                      "Number consumer processes");
  sc_options_add_int (opt, 'b', "prod_processes", &g->proprocs, mpisize,
                      "Number producer processes");
  sc_options_add_int (opt, 'c', "cons_minlevel", &g->con.cminl, 0,
                      "Lowest consumer level");
  sc_options_add_int (opt, 'p', "prod_minlevel", &g->pro.pminl, 0,
                      "Lowest producer level");
  /* examples:
   *   0,1 - arch and brick (best suited for 3D)
   *   2,3 - identical unit squares/cubes */
  sc_options_add_int (opt, 'e', "example", &g->example, 0,
                      "Example mapping index");
  /* refinement methods:
   *   0 - refinement around mesh intersection area
   *   1 - increasing refinement around specific point in space
   *   2 - arbitrary refinement based on child ids and ranks
   *   3 - refinement on boundary of pentagon (best suited for 2D) */
  sc_options_add_int (opt, 'r', "refine_option", &g->refinement_method, 0,
                      "Refinement pattern");
  sc_options_add_int (opt, 'm', "max_level", &g->refinement_maxlevel, 3,
                      "Maximum refinement level");
  sc_options_add_bool (opt, 'v', "output_vtk", &g->output_vtk, 0,
                       "VTK output");
  sc_options_add_bool (opt, 't', "output_text", &g->output_text, 0,
                       "Text output");

  /* proceed in run-once loop for clean abort */
  ue = 0;
  do {
    first_argc = sc_options_parse (p4est_package_id, SC_LP_DEFAULT,
                                   opt, argc, argv);
    if (first_argc < 0) {
      ue = usagerr (opt, "Invalid option format");
      break;
    }
    sc_options_print_summary (p4est_package_id, SC_LP_ESSENTIAL, opt);

    /* check options for consistency */
    if (g->conprocs < 1 || g->conprocs > mpisize) {
      ue = usagerr (opt, "Consumer process count between 1 and mpisize");
    }
    if (g->proprocs < 1 || g->proprocs > mpisize) {
      ue = usagerr (opt, "Producer process count between 1 and mpisize");
    }
    if (g->con.cminl < 0 || g->con.cminl > P4EST_OLD_QMAXLEVEL) {
      ue = usagerr (opt, "Consumer minlevel between 0 and P4EST_QMAXLEVEL");
    }
    if (g->pro.pminl < 0 || g->pro.pminl > P4EST_OLD_QMAXLEVEL) {
      ue =
        usagerr (opt,
                 "Producer minlevel between minlevel and P4EST_QMAXLEVEL");
    }
    if (g->example < 0 || g->example > 3) {
      ue = usagerr (opt, "Example between 0 and 3.");
    }
    if (g->refinement_method < 0 || g->refinement_method > 3) {
      ue = usagerr (opt, "Refinement method between 0 and 3.");
    }
    if (g->refinement_maxlevel < 0
        || g->refinement_maxlevel > P4EST_OLD_QMAXLEVEL) {
      ue = usagerr (opt, "Maxlevel between minlevel and P4EST_QMAXLEVEL");
    }
    if (ue) {
      break;
    }

    /* create, refine and partition a consumer and a producer mesh */
    apps_init (g, mpicomm);

    /* create query points, run an exchange and evaluate for simple case */
    apps_run (g);

    /* destroy consumer and producer mesh */
    apps_reset (g);
  }
  while (0);
  if (ue) {
    sc_options_print_usage (p4est_package_id, SC_LP_ERROR, opt, NULL);
  }

  sc_options_destroy (opt);
  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return EXIT_SUCCESS;
}
