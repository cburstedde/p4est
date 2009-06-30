/*
  This file is part of p4est.
  p4est is a C library to manage a parallel collection of quadtrees and/or
  octrees.

  Copyright (C) 2009 Carsten Burstedde, Lucas Wilcox,
                     Toby Isaac.

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef P4_TO_P8
#include <p4est_bits.h>
#include <p4est_algorithms.h>
#include <p4est_ghost.h>
#include <p4est_iterator.h>
#else
#include <p8est_bits.h>
#include <p8est_algorithms.h>
#include <p8est_ghost.h>
#include <p8est_iterator.h>
#endif

typedef struct
{
  p4est_topidx_t      a;
  int64_t             sum;
}
user_data_t;

#ifndef P4_TO_P8
static int          refine_level = 5;
#else
static int          refine_level = 3;
#endif

static void
init_fn (p4est_t * p4est, p4est_topidx_t which_tree,
         p4est_quadrant_t * quadrant)
{
  user_data_t        *data = quadrant->p.user_data;

  data->a = which_tree;
  data->sum = quadrant->x + quadrant->y + quadrant->level;
}

static int
refine_fn (p4est_t * p4est, p4est_topidx_t which_tree,
           p4est_quadrant_t * quadrant)
{
  int                 cid;

  if (which_tree == 2 || which_tree == 3) {
    return 0;
  }

  cid = p4est_quadrant_child_id (quadrant);

  if (cid == P4EST_CHILDREN - 1 ||
      (quadrant->x >= P4EST_LAST_OFFSET (P4EST_MAXLEVEL - 2) &&
       quadrant->y >= P4EST_LAST_OFFSET (P4EST_MAXLEVEL - 2)
#ifdef P4_TO_P8
       && quadrant->z >= P4EST_LAST_OFFSET (P4EST_MAXLEVEL - 2)
#endif
      )) {
    return 1;
  }
  if ((int) quadrant->level >= (refine_level - (int) (which_tree % 3))) {
    return 0;
  }
  if (quadrant->level == 1 && cid == 2) {
    return 1;
  }
  if (quadrant->x == P4EST_QUADRANT_LEN (2) &&
      quadrant->y == P4EST_LAST_OFFSET (2)) {
    return 1;
  }
  if (quadrant->y >= P4EST_QUADRANT_LEN (2)) {
    return 0;
  }

  return 1;
}

static int
weight_one (p4est_t * p4est, p4est_topidx_t which_tree,
            p4est_quadrant_t * quadrant)
{
  return 1;
}

static void
face_do_nothing (p4est_fcb_info_t * info, void *data)
{
};

#ifdef P4_TO_P8
static void
edge_do_nothing (p4est_ecb_info_t * info, void *data)
{
};

static void
edge_test_adjacency (p4est_ecb_info_t * info, void *data)
{
  p8est_quadrant_t    temp, temp2;
  p8est_quadrant_t   *q, *p;
  int                 n = info->quads->elem_count;
  int                 i, j, k, *e, *c, f;
  bool                failure;
  int8_t              level;
  int8_t              min_level = P4EST_MAXLEVEL;
  int8_t              max_level = 0;
  int                 id;
  p4est_topidx_t     *t, *t2, nt;
  p8est_edge_info_t   ei;
  p8est_connectivity_t *conn = info->p8est->connectivity;
  int                 transform[9];
  for (i = 0; i < n; i++) {
    t = sc_array_index (info->tree, (size_t) i);
    q = sc_array_index (info->quads, (size_t) i);
    e = sc_array_index (info->edge_in_zorder, (size_t) i);
    level = q->level;
    min_level = (level < min_level) ? level : min_level;
    max_level = (level > max_level) ? level : max_level;
    for (j = 0; j < 2; j++) {
      f = p8est_edge_faces[*e][j];
      p8est_quadrant_face_neighbor (q, f, &temp);
      if (p8est_quadrant_is_outside_face (&temp)) {
        P4EST_ASSERT (info->intra_tree == false);
        p8est_find_face_transform (conn, *t, f, transform);
        p8est_quadrant_transform_face (&temp, &temp2, transform);
        temp = temp2;
        nt = conn->tree_to_tree[(*t) * 6 + f];
      }
      else {
        nt = *t;
      }
      failure = true;
      for (k = 0; k < n; k++) {
        t2 = sc_array_index (info->tree, (size_t) k);
        if (*t2 == nt) {
          c = sc_array_index (info->common_corner, (size_t) k);
          p = sc_array_index (info->quads, k);
          if (p8est_quadrant_is_equal (&temp, p)) {
            failure = false;
          }
          else if (p8est_quadrant_is_parent (p, &temp)) {
            id = p8est_quadrant_child_id (&temp);
            if (id == *c) {
              failure = false;
            }
          }
          else if (p8est_quadrant_is_parent (&temp, p)) {
            id = p8est_quadrant_child_id (p);
            if (id == *c) {
              failure = false;
            }
          }
        }
      }
      if (failure == true) {
        printf ("\ni %d j %d f %d *t %d nt %d\n", i, j, f, *t, nt);
        printf ("temp x %x y %x z %x level %d\n", temp.x, temp.y, temp.z,
                temp.level);
        for (k = 0; k < n; k++) {
          p = sc_array_index (info->quads, k);
          P4EST_ASSERT (p4est_quadrant_is_inside_root (p));
          printf ("x %x y %x z %x level %d\n", p->x, p->y, p->z, p->level);
          e = sc_array_index (info->edge_in_zorder, k);
          printf ("edge %d\n", *e);
          c = sc_array_index (info->common_corner, k);
          printf ("corner %d\n", *c);
        }
      }
      P4EST_ASSERT (failure == false);
    }
  }
  P4EST_ASSERT (max_level - min_level <= 1);
};
#endif

static void
face_test_adjacency (p4est_fcb_info_t * info, void *data)
{
  p4est_quadrant_t    temp, temp2, temp3, tempv[P4EST_CHILDREN];
  p4est_quadrant_t   *left = info->left_quad;
  p4est_quadrant_t   *right = info->right_quad;
  int8_t              left_outgoing_face = info->left_outgoing_face;
  int8_t              right_outgoing_face = info->right_outgoing_face;
  p4est_locidx_t      left_tree = info->left_tree;
  p4est_locidx_t      right_tree = info->right_tree;
#ifndef P4_TO_P8
  p4est_connectivity_t *conn = info->p4est->connectivity;
  int                 transform_left;
  int                 transform_right;
  left_outgoing_face = p4est_zface_to_rface[left_outgoing_face];
  right_outgoing_face = p4est_zface_to_rface[right_outgoing_face];
  transform_left = p4est_find_face_transform (conn, left_tree,
                                              left_outgoing_face);
  transform_right = p4est_find_face_transform (conn, right_tree,
                                               right_outgoing_face);
#else
  p4est_connectivity_t *conn = info->p8est->connectivity;
  int                 transform_left[9];
  int                 transform_right[9];
  p8est_find_face_transform (conn, left_tree, left_outgoing_face,
                             transform_left);
  p8est_find_face_transform (conn, right_tree, right_outgoing_face,
                             transform_right);
#endif
  if (!info->hanging_flag) {
    if (info->intra_tree) {
      p4est_quadrant_face_neighbor (left, left_outgoing_face, &temp);
      P4EST_ASSERT (p4est_quadrant_is_equal (&temp, right));
      p4est_quadrant_face_neighbor (right, right_outgoing_face, &temp);
      P4EST_ASSERT (p4est_quadrant_is_equal (&temp, left));
    }
    else {
      p4est_quadrant_face_neighbor (left, left_outgoing_face, &temp);
#ifndef P4_TO_P8
      P4EST_ASSERT (conn->tree_to_tree[left_tree * 4 + left_outgoing_face] ==
                    right_tree);
      if (!
          (conn->tree_to_tree[right_tree * 4 + right_outgoing_face] ==
           left_tree)) {
        printf ("%d %d %d %d\n", left_tree, right_tree, left_outgoing_face,
                right_outgoing_face);
      }
      P4EST_ASSERT (conn->
                    tree_to_tree[right_tree * 4 + right_outgoing_face] ==
                    left_tree);
      p4est_quadrant_translate_face (&temp, left_outgoing_face);
      p4est_quadrant_transform_face (&temp, &temp2, transform_left);
#else
      p8est_quadrant_transform_face (&temp, &temp2, transform_left);
#endif
      P4EST_ASSERT (p4est_quadrant_is_equal (&temp2, right));
      p4est_quadrant_face_neighbor (right, right_outgoing_face, &temp);
#ifndef P4_TO_P8
      p4est_quadrant_translate_face (&temp, right_outgoing_face);
      p4est_quadrant_transform_face (&temp, &temp2, transform_right);
#else
      p8est_quadrant_transform_face (&temp, &temp2, transform_right);
#endif
      P4EST_ASSERT (p4est_quadrant_is_equal (&temp2, left));
    }
  }
  else {
    if (info->intra_tree) {
      p4est_quadrant_childrenv (left, tempv);
      temp = tempv[info->left_corner];
      temp3 = temp;
      p4est_quadrant_face_neighbor (&temp, left_outgoing_face, &temp2);
      P4EST_ASSERT (p4est_quadrant_is_equal (&temp2, right));
      p4est_quadrant_face_neighbor (right, right_outgoing_face, &temp);
      P4EST_ASSERT (p4est_quadrant_is_equal (&temp, &temp3));
    }
    else {
      p4est_quadrant_childrenv (left, tempv);
      temp = tempv[info->left_corner];
      temp3 = temp;
      p4est_quadrant_face_neighbor (&temp, left_outgoing_face, &temp2);
#ifndef P4_TO_P8
      p4est_quadrant_translate_face (&temp2, left_outgoing_face);
      p4est_quadrant_transform_face (&temp2, &temp, transform_left);
#else
      p8est_quadrant_transform_face (&temp2, &temp, transform_left);
#endif
      P4EST_ASSERT (p4est_quadrant_is_equal (&temp, right));
      p4est_quadrant_face_neighbor (right, right_outgoing_face, &temp);
#ifndef P4_TO_P8
      p4est_quadrant_translate_face (&temp, right_outgoing_face);
      p4est_quadrant_transform_face (&temp, &temp2, transform_right);
#else
      p8est_quadrant_transform_face (&temp, &temp2, transform_right);
#endif
      P4EST_ASSERT (p4est_quadrant_is_equal (&temp2, &temp3));
    }
  }
}

int
main (int argc, char **argv)
{
  MPI_Comm            mpicomm;
  int                 mpiret;
  int                 mpisize, mpirank;
  p4est_t            *p4est;
  p4est_connectivity_t *connectivity;
  p4est_locidx_t      num_quads, array_size, num_ghosts, ghost_size;
  sc_array_t         *QtoQ, *QtoFO, *QtoH;
  sc_array_t         *gQtoQ, *gQtoFO, *gQtoH;
  sc_array_t          ghost_layer;
  int                *ghost_owner;
  bool                success;
  int                 ntests;
  int                 i;

#ifndef P4_TO_P8
  ntests = 5;
# else
  ntests = 6;
#endif

  /* initialize MPI */
  mpiret = MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  mpicomm = MPI_COMM_WORLD;
  mpiret = MPI_Comm_size (mpicomm, &mpisize);
  SC_CHECK_MPI (mpiret);
  mpiret = MPI_Comm_rank (mpicomm, &mpirank);
  SC_CHECK_MPI (mpiret);

  return 0;

  sc_init (mpicomm, true, true, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  for (i = 0; i < ntests; i++) {
    P4EST_GLOBAL_PRODUCTIONF ("Begin adjacency test %d\n", i);
    /* create connectivity and forest structures */
    switch (i) {
#ifndef P4_TO_P8
    case 0:
      connectivity = p4est_connectivity_new_unitsquare ();
      break;
    case 1:
      connectivity = p4est_connectivity_new_corner ();
      break;
    case 2:
      connectivity = p4est_connectivity_new_moebius ();
      break;
    case 3:
      connectivity = p4est_connectivity_new_star ();
      break;
    default:
      connectivity = p4est_connectivity_new_periodic ();
      break;
#else
    case 0:
      connectivity = p8est_connectivity_new_unitcube ();
      break;
    case 1:
      connectivity = p8est_connectivity_new_periodic ();
      break;
    case 2:
      connectivity = p8est_connectivity_new_rotwrap ();
      break;
    case 3:
      connectivity = p8est_connectivity_new_twocubes ();
      break;
    case 4:
      connectivity = p8est_connectivity_new_rotcubes ();
      break;
    default:
      connectivity = p8est_connectivity_new_shell ();
      break;
#endif
    }
    p4est = p4est_new (mpicomm, connectivity, 15,
                       sizeof (user_data_t), init_fn, NULL);

    /* refine to make the number of elements interesting */
    p4est_refine (p4est, true, refine_fn, init_fn);

    /* balance the forest */
    p4est_balance (p4est, P4EST_BALANCE_DEFAULT, init_fn);

    /* do a uniform partition, include the weight function for testing */
    p4est_partition (p4est, weight_one);

    /* create ghost layer */
    sc_array_init (&ghost_layer, sizeof (p4est_quadrant_t));
    success = p4est_build_ghost_layer (p4est, P4EST_BALANCE_DEFAULT,
                                       &ghost_layer, &ghost_owner);
    P4EST_ASSERT (success);

    /* create and fill adjacency arrays */
    num_quads = p4est->local_num_quadrants;
    array_size = num_quads * P4EST_DIM * 2;
    QtoQ = sc_array_new (sizeof (p4est_gloidx_t));
    sc_array_resize (QtoQ, (size_t) array_size);
    QtoFO = sc_array_new (sizeof (uint8_t));
    sc_array_resize (QtoFO, (size_t) array_size);
    QtoH = sc_array_new (sizeof (int8_t));
    sc_array_resize (QtoH, (size_t) array_size);
    num_ghosts = ghost_layer.elem_count;
    ghost_size = num_ghosts * P4EST_DIM * 2;
    gQtoQ = sc_array_new (sizeof (p4est_gloidx_t));
    sc_array_resize (gQtoQ, (size_t) ghost_size);
    gQtoFO = sc_array_new (sizeof (uint8_t));
    sc_array_resize (gQtoFO, (size_t) ghost_size);
    gQtoH = sc_array_new (sizeof (int8_t));
    sc_array_resize (gQtoH, (size_t) ghost_size);

#ifndef P4_TO_P8
    p4est_iterator (p4est, &ghost_layer, NULL, NULL, face_test_adjacency,
                    NULL);
#else
    printf ("test %d\n", i);
    p8est_iterator (p4est, &ghost_layer, NULL,
                    NULL, face_test_adjacency, edge_test_adjacency, NULL);
#endif

    /* clean up */
    sc_array_destroy (QtoQ);
    sc_array_destroy (QtoFO);
    sc_array_destroy (QtoH);
    sc_array_destroy (gQtoQ);
    sc_array_destroy (gQtoFO);
    sc_array_destroy (gQtoH);
    sc_array_reset (&ghost_layer);
    P4EST_FREE (ghost_owner);

    p4est_destroy (p4est);
    p4est_connectivity_destroy (connectivity);
    P4EST_GLOBAL_PRODUCTIONF ("End adjacency test %d\n", i);
  }

  /* exit */
  sc_finalize ();

  mpiret = MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}

/* EOF test_iterator2.c */
