/*
  This file is part of p4est.
  p4est is a C library to manage a collection (a forest) of multiple
  connected adaptive quadtrees or octrees in parallel.

  Copyright (C) 2010 The University of Texas System
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

/** \file p4est_step4.c
 *
 * This 2D example program solves the Poisson equation using finite elements.
 * Currently, it works on the unit square.  For more general domains, a
 * coordinate transformation to physical space would need to be implemented.
 * This will usually entail using a quadrature instead of exact integration.
 * The check for boundary nodes would need to be adapted to the new geometry.
 */

/* p4est has two separate interfaces for 2D and 3D, p4est*.h and p8est*.h.
 * Most API functions are available for both dimensions.  The header file
 * p4est_to_p8est.h #define's the 2D names to the 3D names such that most code
 * only needs to be written once.  In this example, we rely on this. */
#ifndef P4_TO_P8
#include <p4est_bits.h>
#include <p4est_ghost.h>
#include <p4est_lnodes.h>
#include <p4est_vtk.h>
#else
#include <p8est_bits.h>
#include <p8est_ghost.h>
#include <p8est_lnodes.h>
#include <p8est_vtk.h>
#endif

/** List number of possible independent nodes for each hanging node. */
static const int    corner_num_hanging[P4EST_CHILDREN] =
#ifndef P4_TO_P8
{ 1, 2, 2, 1 }
#else
{ 1, 2, 2, 4, 2, 4, 4, 1 }
#endif
;
static const int    zero = 0;           /**< Constant zero. */
static const int    ones = P4EST_CHILDREN - 1;  /**< One bit per dimension. */

/** For each node i of the reference quadrant, corner_num_hanging[i] many. */
static const int   *corner_to_hanging[P4EST_CHILDREN];

/** Callback function to decide on refinement.
 *
 * This function is called for every processor-local quadrant in order; its
 * return value is understood as a boolean refinement flag.  We refine around a
 * h = 1/8 block with left front lower corner (5/8, 2/8, 6/8).
 */
static int
refine_fn (p4est_t * p4est, p4est_topidx_t which_tree,
           p4est_quadrant_t * quadrant)
{
  /* Compute the integer coordinate extent of a quadrant of length 2^(-3). */
  const p4est_qcoord_t eighth = P4EST_QUADRANT_LEN (3);

  /* Compute the length of the current quadrant in integer coordinates. */
  const p4est_qcoord_t length = P4EST_QUADRANT_LEN (quadrant->level);

  /* Refine if the quadrant intersects the block in question. */
  return ((quadrant->x + length > 5 * eighth && quadrant->x < 6 * eighth) &&
          (quadrant->y + length > 2 * eighth && quadrant->y < 3 * eighth) &&
#ifdef P4_TO_P8
          (quadrant->z + length > 6 * eighth && quadrant->z < 7 * eighth) &&
#endif
          1);
}

/** Right hand side function for the 2D Poisson problem.
 *
 * This is the negative Laplace operator acting on the function \a uexact.
 * \param [in] vxyz    x, y, and z coordinates in physical space.
 *                     Even for the 2D case z is well defined, but ignored here.
 * \return             Scalar function value at vxyz.
 */
static double
func_rhs (const double vxyz[3])
{
  const double        x = vxyz[0];
  const double        y = vxyz[1];
#ifdef P4_TO_P8
  const double        z = vxyz[2];

  return 128. * (x * (1. - x) * y * (1. - y) +
                 y * (1. - y) * z * (1. - z) + z * (1. - z) * x * (1. - x));
#else
  return 32. * (x * (1. - x) + y * (1. - y));
#endif
}

/** Exact solution for the 2D Poisson problem.
 *
 * We pick a function with zero Dirichlet boundary conditions on the unit square.
 * \param [in] vxyz    x, y, and z coordinates in physical space.
 *                     Even for the 2D case z is well defined, but ignored here.
 * \return             Scalar function value at vxyz.
 */
static double
func_uexact (const double vxyz[3])
{
  const double        x = vxyz[0];
  const double        y = vxyz[1];
#ifdef P4_TO_P8
  const double        z = vxyz[2];
#endif

  return 4. * 4.
#ifdef P4_TO_P8
    * 4 * z * (1. - z)
#endif
    * x * (1. - x) * y * (1. - y);
}

/** Determine the boundary status on the unit square/cube.
 * \param [in] p4est    Can be used to access the connectivity.
 * \param [in] tt       The tree number (always zero for the unit square).
 * \param [in] node     The corner node of an element to be examined.
 * \return              True for Dirichlet boundary, false otherwise.
 */
static int
is_boundary_unitsquare (p4est_t * p4est, p4est_topidx_t tt,
                        p4est_quadrant_t * node)
{
  /* For this simple connectivity it is sufficient to check x, y (and z). */
  return (node->x == 0 || node->x == P4EST_ROOT_LEN ||
          node->y == 0 || node->y == P4EST_ROOT_LEN ||
#ifdef P4_TO_P8
          node->z == 0 || node->z == P4EST_ROOT_LEN ||
#endif
          0);
}

/** Decode the information from p{4,8}est_lnodes_t for a given element.
 *
 * \see p4est_lnodes.h for an in-depth discussion of the encoding.
 * \param [in] face_code         Bit code as defined in p{4,8}est_lnodes.h.
 * \param [out] hanging_corner   Undefined if no node is hanging.
 *                               If any node is hanging, this contains
 *                               one integer per corner, which is -1
 *                               for corners that are not hanging,
 *                               and the number of the non-hanging
 *                               corner on the hanging face/edge otherwise.
 *                               For faces in 3D, it is diagonally opposite.
 * \return true if any node is hanging, false otherwise.
 */
static int
lnodes_decode2 (p4est_lnodes_code_t face_code,
                int hanging_corner[P4EST_CHILDREN])
{
  if (face_code) {
    const int           c = (int) (face_code & ones);
    int                 i, h;
    int                 work = (int) (face_code >> P4EST_DIM);

    /* These two corners are never hanging by construction. */
    hanging_corner[c] = hanging_corner[c ^ ones] = -1;
    for (i = 0; i < P4EST_DIM; ++i) {
      /* Process face hanging corners. */
      h = c ^ (1 << i);
      hanging_corner[h ^ ones] = (work & 1) ? c : -1;
#ifdef P4_TO_P8
      /* Process edge hanging corners. */
      hanging_corner[h] = (work & P4EST_CHILDREN) ? c : -1;
#endif
      work >>= 1;
    }
    return 1;
  }
  return 0;
}

/** Compute values at hanging nodes by interpolation.
 * A face hanging node in 3D depends on the four corner nodes of the face,
 * edge hanging nodes or face hanging nodes in 2D depend on two nodes.
 * This function works in place, we have to be careful about the ordering.
 * Face hanging node values are not reused, so they are overwritten first.
 * \param [in] face_code    This number encodes the child id of the quadrant
 *                          and the hanging status of faces and edges.
 * \param [in,out] inplace  On input, the values at the independent nodes.
 *                          On output, interpolated to hanging node locations.
 */
static void
interpolate_hanging_nodes (p4est_lnodes_code_t face_code,
                           double inplace[P4EST_CHILDREN])
{
  const int           c = (int) (face_code & ones);
  int                 i, j;
  int                 ef;
  int                 work = (int) (face_code >> P4EST_DIM);
  double              sum;
  const double        factor = 1. / P4EST_HALF;

  /* Compute face hanging nodes first (this is all there is in 2D). */
  for (i = 0; i < P4EST_DIM; ++i) {
    if (work & 1) {
      ef = p4est_corner_faces[c][i];
      sum = 0.;
      for (j = 0; j < P4EST_HALF; ++j) {
        sum += inplace[p4est_face_corners[ef][j]];
      }
      inplace[c ^ ones ^ (1 << i)] = factor * sum;
    }
    work >>= 1;
  }

#ifdef P4_TO_P8
  /* Compute edge hanging nodes afterwards */
  for (i = 0; i < P4EST_DIM; ++i) {
    if (work & 1) {
      ef = p8est_corner_edges[c][i];
      inplace[c ^ (1 << i)] = .5 * (inplace[p8est_edge_corners[ef][0]] +
                                    inplace[p8est_edge_corners[ef][1]]);
    }
    work >>= 1;
  }
#endif
}

/** Parallel sum of values in node vector across all sharers.
 *
 * This function is necessary in the matrix-vector product since elements
 * from multiple processors can contribute to any given node value.
 *
 * \param [in] p4est      The mesh is not changed.
 * \param [in] lnodes     The node numbering is not changed.
 * \param [in,out] v      On input, vector with local contributions.
 *                        On output, the node-wise sum across all sharers.
 */
static void
share_sum (p4est_t * p4est, p4est_lnodes_t * lnodes, double *v)
{
  const int           nloc = lnodes->num_local_nodes;
  const int           npeers = (int) lnodes->sharers->elem_count;
  int                 iq, jn;
  int                 gl;
  sc_array_t          node_data;
  p4est_lnodes_rank_t *lrank;
  p4est_lnodes_buffer_t *buffer;

  sc_array_init_data (&node_data, v, sizeof (double), nloc);
  buffer = p4est_lnodes_share_all (&node_data, lnodes);

  for (iq = 0; iq < npeers; ++iq) {
    lrank = (p4est_lnodes_rank_t *) sc_array_index_int (lnodes->sharers, iq);
    sc_array_t         *recv_data =
      (sc_array_t *) sc_array_index_int (buffer->recv_buffers, iq);
    P4EST_ASSERT (recv_data->elem_size == node_data.elem_size);

    if (lrank->rank != p4est->mpirank) {
      const int           nshared = (int) lrank->shared_nodes.elem_count;
      const double       *w = (const double *) recv_data->array;

      P4EST_ASSERT ((int) recv_data->elem_count == nshared);

      for (jn = 0; jn < nshared; ++jn) {
        gl = (int)
          *(p4est_locidx_t *) sc_array_index_int (&lrank->shared_nodes, jn);
        P4EST_ASSERT (0 <= gl && gl < nloc);
        v[gl] += w[jn];
      }
    }
#ifdef P4EST_ENABLE_DEBUG
    else {
      P4EST_ASSERT (recv_data->elem_count == 0);
      P4EST_ASSERT (lrank->owned_offset == 0);
      P4EST_ASSERT (lrank->owned_count == lnodes->owned_count);
    }
#endif
  }

  p4est_lnodes_buffer_destroy (buffer);
  sc_array_reset (&node_data);
}

/** Compute the inner product of two node vectors in parallel.
 *
 * \param [in] p4est          The forest is not changed.
 * \param [in] lnodes         The node numbering is not changed.
 * \param [in] v1             First node vector.
 * \param [in] v2             Second node vector.
 * \return                    Parallel l_2 inner product over the domain.
 */
static double
vector_dot (p4est_t * p4est, p4est_lnodes_t * lnodes,
            const double *v1, const double *v2)
{
  const int           nown = lnodes->owned_count;
  int                 mpiret;
  int                 lnid;
  double              lsum, gsum;

  /* We only sum the locally owned values to avoid double counting. */
  lsum = 0.;
  for (lnid = 0; lnid < nown; ++lnid) {
    lsum += v1[lnid] * v2[lnid];
  }

  /* The result is made available on all processors. */
  mpiret = sc_MPI_Allreduce (&lsum, &gsum, 1, sc_MPI_DOUBLE, sc_MPI_SUM,
                             p4est->mpicomm);
  SC_CHECK_MPI (mpiret);

  return gsum;
}

/** Compute y := y + a * x.
 *
 * \param [in] p4est          The forest is not changed.
 * \param [in] lnodes         The node numbering is not changed.
 * \param [in] a              The scalar.
 * \param [in] x              First node vector.
 * \param [in,out] y          Second node vector.
 */
static void
vector_axpy (p4est_t * p4est, p4est_lnodes_t * lnodes, double a,
             const double *x, double *y)
{
  const int           nloc = lnodes->num_local_nodes;
  int                 lnid;

  for (lnid = 0; lnid < nloc; ++lnid) {
    y[lnid] += a * x[lnid];
  }
}

/** Compute y := x + b * y.
 *
 * \param [in] p4est          The forest is not changed.
 * \param [in] lnodes         The node numbering is not changed.
 * \param [in] x              First node vector.
 * \param [in] b              The scalar.
 * \param [in,out] y          Second node vector.
 */
static void
vector_xpby (p4est_t * p4est, p4est_lnodes_t * lnodes, const double *x,
             double b, double *y)
{
  const int           nloc = lnodes->num_local_nodes;
  int                 lnid;
  double              yy;

  for (lnid = 0; lnid < nloc; ++lnid) {
    yy = y[lnid];
    y[lnid] = x[lnid] + b * yy;
  }
}

/** Zero all entries of a vector.
 *
 * \param [in] p4est          The forest is not changed.
 * \param [in] lnodes         The node numbering is not changed.
 * \param [out] x             The vector.
 */
static void
vector_zero (p4est_t * p4est, p4est_lnodes_t * lnodes, double *x)
{
  const int           nloc = lnodes->num_local_nodes;

  memset (x, 0., nloc * sizeof (double));
}

/** copy a vector.
 *
 * \param [in] p4est          The forest is not changed.
 * \param [in] lnodes         The node numbering is not changed.
 * \param [in] x              Input node vector.
 * \param [out] y             output node vector.
 */
static void
vector_copy (p4est_t * p4est, p4est_lnodes_t * lnodes, const double *x,
             double *y)
{
  const int           nloc = lnodes->num_local_nodes;

  memcpy (y, x, nloc * sizeof (double));
}

/** Allocate storage for processor-relevant nodal degrees of freedom.
 *
 * \param [in] lnodes   This structure is queried for the node count.
 * \return              Allocated double array; must be freed with P4EST_FREE.
 */
static double      *
allocate_vector (p4est_lnodes_t * lnodes)
{
  return P4EST_ALLOC (double, lnodes->num_local_nodes);
}

/** Interpolate right hand side and exact solution onto mesh nodes.
 *
 * \param [in] p4est          The forest is not changed.
 * \param [in] lnodes         The node numbering is not changed.
 * \param [out] rhs_eval      Is allocated and filled with function values.
 * \param [out] uexact_eval   Is allocated and filled with function values.
 * \param [out] pbc           Boolean flags for Dirichlet boundary nodes.
 */
static void
interpolate_functions (p4est_t * p4est, p4est_lnodes_t * lnodes,
                       double **rhs_eval, double **uexact_eval, int8_t ** pbc)
{
  const p4est_locidx_t nloc = lnodes->num_local_nodes;
  int                 anyhang, hanging_corner[P4EST_CHILDREN];
  int                 i;        /* We use plain int for small loops. */
  double             *rhs, *uexact;
  double              vxyz[3];  /* We embed the 2D vertices into 3D space. */
  int8_t             *bc;
  p4est_topidx_t      tt;       /* Connectivity variables have this type. */
  p4est_locidx_t      k, q, Q;  /* Process-local counters have this type. */
  p4est_locidx_t      lni;      /* Node index relative to this processor. */
  p4est_tree_t       *tree;     /* Pointer to one octree */
  p4est_quadrant_t   *quad, *parent, sp, node;
  sc_array_t         *tquadrants;       /* Quadrant array for one tree */

  rhs = *rhs_eval = allocate_vector (lnodes);
  uexact = *uexact_eval = allocate_vector (lnodes);
  bc = *pbc = P4EST_ALLOC (int8_t, nloc);
  memset (bc, -1, sizeof (int8_t) * nloc);      /* Indicator for visiting. */

  /* We need to compute the xyz locations of non-hanging nodes to evaluate the
   * given functions.  For hanging nodes, we have to look at the corresponding
   * independent nodes.  Usually we would cache this information, here we only
   * need it once and throw it away again.
   * We also compute boundary status of independent nodes. */
  for (tt = p4est->first_local_tree, k = 0;
       tt <= p4est->last_local_tree; ++tt) {
    tree = p4est_tree_array_index (p4est->trees, tt);   /* Current tree */
    tquadrants = &tree->quadrants;
    Q = (p4est_locidx_t) tquadrants->elem_count;
    for (q = 0; q < Q; ++q, ++k) {
      /* This is now a loop over all local elements.
       * Users might aggregate the above code into a more compact iterator. */
      quad = p4est_quadrant_array_index (tquadrants, q);

      /* We need to determine whether any node on this element is hanging. */
      anyhang = lnodes_decode2 (lnodes->face_code[q], hanging_corner);
      if (!anyhang) {
        parent = NULL;          /* Defensive programming. */
      }
      else {
        /* At least one node is hanging.  We need the parent quadrant to
         * find the location of the corresponding non-hanging node. */
        parent = &sp;
        p4est_quadrant_parent (quad, parent);
      }
      for (i = 0; i < P4EST_CHILDREN; ++i) {
        lni = lnodes->element_nodes[P4EST_CHILDREN * k + i];
        P4EST_ASSERT (lni >= 0 && lni < nloc);
        if (bc[lni] < 0) {
          if (anyhang && hanging_corner[i] >= 0) {
            /* This node is hanging; access the referenced node instead. */
            p4est_quadrant_corner_node (parent, i, &node);
          }
          else {
            p4est_quadrant_corner_node (quad, i, &node);
          }

          /* Determine boundary status of independent node. */
          bc[lni] = is_boundary_unitsquare (p4est, tt, &node);

          /* Transform per-tree reference coordinates into physical space. */
          p4est_qcoord_to_vertex (p4est->connectivity, tt, node.x, node.y,
#ifdef P4_TO_P8
                                  node.z,
#endif
                                  vxyz);

          /* Use physical space coordinates to evaluate functions */
          rhs[lni] = func_rhs (vxyz);
          uexact[lni] = func_uexact (vxyz);
        }
      }
    }
  }
}

/** Apply a finite element matrix to a node vector, y = Mx.
 * \param [in] p4est          The forest is not changed.
 * \param [in] lnodes         The node numbering is not changed.
 * \param [in] bc             Boolean flags for Dirichlet boundary nodes.
 *                            If NULL, no special action is taken.
 * \param [in] stiffness      If false use scaling for the mass matrix,
 *                            if true use the scaling for stiffness matrix.
 * \param [in] matrix         A 4x4 matrix computed on the reference element.
 * \param [in] in             Input vector x.
 * \param [out] out           Output vector y = Mx.
 */
static void
multiply_matrix (p4est_t * p4est, p4est_lnodes_t * lnodes, const int8_t * bc,
                 int stiffness,
                 double (*matrix)[P4EST_CHILDREN][P4EST_CHILDREN],
                 const double *in, double *out)
{
  const int           nloc = lnodes->num_local_nodes;
  const int           nown = lnodes->owned_count;
  int                 c, h;
  int                 i, j, k;
  int                 q, Q;
  int                 anyhang, hanging_corner[P4EST_CHILDREN];
  int                 isboundary[P4EST_CHILDREN];
  int                 ncontrib;
  const int          *contrib_corner;
  double              factor, sum, inloc[P4EST_CHILDREN];
  sc_array_t         *tquadrants;
  p4est_topidx_t      tt;
  p4est_locidx_t      lni, all_lni[P4EST_CHILDREN];
  p4est_tree_t       *tree;
  p4est_quadrant_t   *quad;

  /* Initialize the output vector. */
  if (bc == NULL) {
    /* No boundary values, output vector has all zero values. */
    for (lni = 0; lni < nloc; ++lni) {
      out[lni] = 0.;
    }
  }
  else {
    /* We have boundary conditions, switch on the locally owned nodes. */
    for (lni = 0; lni < nown; ++lni) {
      out[lni] = bc[lni] ? in[lni] : 0.;
    }
    /* The node values owned by other processors will be added there. */
    for (lni = nown; lni < nloc; ++lni) {
      out[lni] = 0.;
    }
  }

  /* Loop over local quadrants to apply the element matrices. */
  for (tt = p4est->first_local_tree, k = 0;
       tt <= p4est->last_local_tree; ++tt) {
    tree = p4est_tree_array_index (p4est->trees, tt);
    tquadrants = &tree->quadrants;
    Q = (p4est_locidx_t) tquadrants->elem_count;
    for (q = 0; q < Q; ++q, ++k) {
      quad = p4est_quadrant_array_index (tquadrants, q);

      /* If we are on the 2D unit square, the Jacobian determinant is h^2;
       * for the stiffness matrix this cancels with the inner derivatives.
       * On the 3D unit cube we have an additional power of h. */
#ifndef P4_TO_P8
      factor = !stiffness ? pow (.5, 2. * quad->level) : 1.;
#else
      factor = pow (.5, (!stiffness ? 3. : 1.) * quad->level);
#endif
      for (i = 0; i < P4EST_CHILDREN; ++i) {
        /* Cache some information on corner nodes. */
        lni = lnodes->element_nodes[P4EST_CHILDREN * k + i];
        isboundary[i] = (bc == NULL ? 0 : bc[lni]);
        inloc[i] = !isboundary[i] ? in[lni] : 0.;
        all_lni[i] = lni;
      }

      /* Figure out the hanging corners on this element, if any. */
      anyhang = lnodes_decode2 (lnodes->face_code[k], hanging_corner);

      if (!anyhang) {
        /* No hanging nodes on this element; just apply the matrix. */
        for (i = 0; i < P4EST_CHILDREN; ++i) {
          if (!isboundary[i]) {
            sum = 0.;
            for (j = 0; j < P4EST_CHILDREN; ++j) {
              sum += (*matrix)[i][j] * inloc[j];
            }
            out[all_lni[i]] += factor * sum;
          }
        }
      }
      else {
        /* Compute input values at hanging nodes by interpolation. */
        interpolate_hanging_nodes (lnodes->face_code[k], inloc);

        /* Apply element matrix and then the transpose interpolation. */
        for (i = 0; i < P4EST_CHILDREN; ++i) {
          sum = 0.;
          for (j = 0; j < P4EST_CHILDREN; ++j) {
            sum += (*matrix)[i][j] * inloc[j];
          }
          if (hanging_corner[i] == -1) {
            /* This node is not hanging, there is no need to transform. */
            c = 0;
            ncontrib = 1;
            contrib_corner = &i;
            sum *= factor;
          }
          else {
            /* This node is hanging.  Use hanging node relations from the
             * reference quadrant by transforming corner numbers with ^ c. */
            c = hanging_corner[i];      /* Child id of quadrant. */
            ncontrib = corner_num_hanging[i ^ c];
            contrib_corner = corner_to_hanging[i ^ c];
            sum *= factor / (double) ncontrib;
          }
          /* The result is added into the output vector. */
          for (j = 0; j < ncontrib; ++j) {
            h = contrib_corner[j] ^ c;  /* Inverse transform of node number. */
            if (!isboundary[h]) {
              out[all_lni[h]] += sum;
            }
          }
        }
      }
    }
  }

  /* Parallel sum of result. */
  share_sum (p4est, lnodes, out);
}

/** Set Dirichlet boundary values of a node vector to zero.
 *
 * \param [in] lnodes         The node numbering is not changed.
 * \param [in] bc             Boolean flags for Dirichlet boundary nodes.
 *                            If NULL, this function does nothing.
 * \param [in,out] v          Dirichlet nodes are overwritten with zero.
 */
static void
set_dirichlet (p4est_lnodes_t * lnodes, const int8_t * bc, double *v)
{
  const int           nloc = lnodes->num_local_nodes;
  p4est_locidx_t      lni;

  if (bc == NULL) {
    return;
  }
  for (lni = 0; lni < nloc; ++lni) {
    if (bc[lni]) {
      v[lni] = 0.;
    }
  }
}

/** Multiply the mass matrix with a vector of ones to compute the volume.
 * \param [in] p4est          The forest is not changed.
 * \param [in] lnodes         The node numbering is not changed.
 * \param [in] matrix         The mass matrix should be passed in here.
 * \param [in] tmp            Must be allocated, entries are undefined.
 * \param [in,out] lump       Must be allocated, receives matrix * ones.
 */
static void
test_area (p4est_t * p4est, p4est_lnodes_t * lnodes,
           double (*matrix)[P4EST_CHILDREN][P4EST_CHILDREN],
           double *tmp, double *lump)
{
  const int           nloc = lnodes->num_local_nodes;
  double              dot;
  p4est_locidx_t      lni;

  for (lni = 0; lni < nloc; ++lni) {
    tmp[lni] = 1.;
  }
  multiply_matrix (p4est, lnodes, NULL, 0, matrix, tmp, lump);
  dot = vector_dot (p4est, lnodes, tmp, lump);
  dot = sqrt (dot);

  P4EST_GLOBAL_PRODUCTIONF ("Area of domain: %g\n", dot);
}

/** Execute the conjugate gradient method.
 * \param [in] p4est     The forest is not changed.
 * \param [in] lnodes    The node numbering is not changed.
 * \param [in] bc        Boolean flags for Dirichlet boundary nodes.
 * \param [in] stiffness If false use scaling for the mass matrix,
 *                       if true use the scaling for stiffness matrix.
 * \param [in] matrix    The mass matrix should be passed in here.
 * \param [in] b         The right hand side vector.
 * \param [out] x        The result; we use an initial value of zero.
 */
static void
solve_by_cg (p4est_t * p4est, p4est_lnodes_t * lnodes, const int8_t * bc,
             int stiffness,
             double (*matrix)[P4EST_CHILDREN][P4EST_CHILDREN],
             const double *b, double *x)
{
  const int           imax = 100;
  const double        tol = 1.e-6;

  int                 i;
  double              alpha, beta, pAp;
  double              rr, rrnew, rrorig;
  double             *aux[4];
  double             *r, *p, *Ap;

  for (i = 0; i < 3; ++i) {
    aux[i] = allocate_vector (lnodes);
  }
  r = aux[0];
  p = aux[1];
  Ap = aux[2];

  /* Initialize the solution vector to zero. */
  vector_zero (p4est, lnodes, x);
  vector_copy (p4est, lnodes, b, r);
  vector_copy (p4est, lnodes, b, p);
  rr = rrorig = vector_dot (p4est, lnodes, r, r);

  for (i = 0; i < imax && rr > rrorig * tol * tol; i++) {
    multiply_matrix (p4est, lnodes, bc, stiffness, matrix, p, Ap);
    pAp = vector_dot (p4est, lnodes, p, Ap);
    alpha = rr / pAp;
    vector_axpy (p4est, lnodes, alpha, p, x);
    vector_axpy (p4est, lnodes, -alpha, Ap, r);
    rrnew = vector_dot (p4est, lnodes, r, r);
    beta = rrnew / rr;
    vector_xpby (p4est, lnodes, r, beta, p);
    P4EST_GLOBAL_VERBOSEF ("%03d: r'r %g alpha %g beta %g\n",
                           i, rr, alpha, beta);
    rr = rrnew;
  }
  if (i < imax) {
    P4EST_GLOBAL_PRODUCTIONF ("cg converged to %g in %d iterations\n",
                              sqrt (rr), i);
  }
  else {
    P4EST_GLOBAL_PRODUCTIONF ("cg did not converge (%g) in %d iterations\n",
                              sqrt (rr), imax);
  }

  /* Free temporary storage. */
  for (i = 0; i < 3; ++i) {
    P4EST_FREE (aux[i]);
  }
}

/** Execute the numerical part of the example: Solve Poisson's equation.
 * \param [in] p4est    Solve the PDE with the given mesh refinement.
 */
static void
solve_poisson (p4est_t * p4est)
{
  /* 1D mass matrix on the reference element [0, 1]. */
  static const double m_1d[2][2] = {
    {1 / 3., 1 / 6.},
    {1 / 6., 1 / 3.},
  };
  /* 1D stiffness matrix on the reference element [0, 1]. */
  static const double s_1d[2][2] = {
    {1., -1.},
    {-1., 1.},
  };
  int                 i, j, k, l;
#ifdef P4_TO_P8
  int                 m, n;
#endif
  double              mass_dd[P4EST_CHILDREN][P4EST_CHILDREN];
  double              stiffness_dd[P4EST_CHILDREN][P4EST_CHILDREN];
  double             *rhs_eval, *uexact_eval, *lump;
  double             *rhs_fe, *u_fe, *u_diff, *diff_mass;
  double              err2, err;
  int8_t             *bc;
  p4est_ghost_t      *ghost;
  p4est_lnodes_t     *lnodes;

  /* Create the ghost layer to learn about parallel neighbors. */
  ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FULL);

  /* Create a node numbering for continuous linear finite elements. */
  lnodes = p4est_lnodes_new (p4est, ghost, 1);

  /* Destroy the ghost structure -- no longer needed after node creation. */
  p4est_ghost_destroy (ghost);
  ghost = NULL;

  /* Compute entries of reference mass and stiffness matrices in 2D.
   * In this example we can proceed without numerical integration. */
  for (l = 0; l < 2; ++l) {
    for (k = 0; k < 2; ++k) {
      for (j = 0; j < 2; ++j) {
        for (i = 0; i < 2; ++i) {
#ifndef P4_TO_P8
          mass_dd[2 * j + i][2 * l + k] = m_1d[i][k] * m_1d[j][l];
          stiffness_dd[2 * j + i][2 * l + k] =
            s_1d[i][k] * m_1d[j][l] + m_1d[i][k] * s_1d[j][l];
#else
          for (n = 0; n < 2; ++n) {
            for (m = 0; m < 2; ++m) {
              mass_dd[4 * i + 2 * n + m][4 * l + 2 * k + j] =
                m_1d[m][j] * m_1d[n][k] * m_1d[i][l];
              stiffness_dd[4 * i + 2 * n + m][4 * l + 2 * k + j] =
                s_1d[m][j] * m_1d[n][k] * m_1d[i][l] +
                m_1d[m][j] * s_1d[n][k] * m_1d[i][l] +
                m_1d[m][j] * m_1d[n][k] * s_1d[i][l];
            }
          }
#endif
        }
      }
    }
  }

  /* Assign independent nodes for hanging nodes. */
  corner_to_hanging[0] = &zero;
#ifdef P4_TO_P8
  corner_to_hanging[1] = p8est_edge_corners[0];
  corner_to_hanging[2] = p8est_edge_corners[4];
  corner_to_hanging[3] = p8est_face_corners[4];
  corner_to_hanging[4] = p8est_edge_corners[8];
#endif
  corner_to_hanging[ones - 2] = p4est_face_corners[2];
  corner_to_hanging[ones - 1] = p4est_face_corners[0];
  corner_to_hanging[ones] = &ones;

  /* Test mass matrix multiplication by computing the area of the domain. */
  rhs_fe = allocate_vector (lnodes);
  lump = allocate_vector (lnodes);
  test_area (p4est, lnodes, &mass_dd, rhs_fe, lump);

  /* Interpolate right hand side and exact solution onto mesh nodes. */
  interpolate_functions (p4est, lnodes, &rhs_eval, &uexact_eval, &bc);

  /* Apply mass matrix to create right hand side FE vector. */
  multiply_matrix (p4est, lnodes, bc, 0, &mass_dd, rhs_eval, rhs_fe);
  set_dirichlet (lnodes, bc, rhs_fe);

  /* Run conjugate gradient method with initial value zero. */
  u_fe = allocate_vector (lnodes);
  solve_by_cg (p4est, lnodes, bc, 1, &stiffness_dd, rhs_fe, u_fe);

  /* Compute the pointwise difference with the exact vector. */
  u_diff = allocate_vector (lnodes);
  vector_copy (p4est, lnodes, u_fe, u_diff);
  vector_axpy (p4est, lnodes, -1., uexact_eval, u_diff);

  /* Compute the L2 difference with the exact vector.
   * We know that this is over-optimistic: Quadrature will be sharper.
   * We could also reuse another vector instead of allocating a new one. */
  diff_mass = allocate_vector (lnodes);
  multiply_matrix (p4est, lnodes, bc, 0, &mass_dd, u_diff, diff_mass);
  err2 = vector_dot (p4est, lnodes, diff_mass, u_diff);
  err = sqrt (err2);
  P4EST_GLOBAL_PRODUCTIONF ("||u_fe - u_exact||_L2 = %g\n", err);

  /* Free finite element vectors */
  P4EST_FREE (diff_mass);
  P4EST_FREE (u_diff);
  P4EST_FREE (u_fe);
  P4EST_FREE (lump);
  P4EST_FREE (rhs_fe);
  P4EST_FREE (rhs_eval);
  P4EST_FREE (uexact_eval);
  P4EST_FREE (bc);

  /* We are done with the FE node numbering. */
  p4est_lnodes_destroy (lnodes);
}

/** The main function of the step4 example program.
 *
 * It creates a connectivity and forest, refines it, and solves the Poisson
 * equation with piecewise linear finite elements.
 */
int
main (int argc, char **argv)
{
  int                 mpiret;
  int                 startlevel, endlevel, level;
  sc_MPI_Comm         mpicomm;
  p4est_t            *p4est;
  p4est_connectivity_t *conn;

  /* Initialize MPI; see sc_mpi.h.
   * If configure --enable-mpi is given these are true MPI calls.
   * Else these are dummy functions that simulate a single-processor run. */
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  mpicomm = sc_MPI_COMM_WORLD;

  /* These functions are optional.  If called they store the MPI rank as a
   * static variable so subsequent global p4est log messages are only issued
   * from processor zero.  Here we turn off most of the logging; see sc.h. */
  sc_init (mpicomm, 1, 1, NULL, SC_LP_ESSENTIAL);
  p4est_init (NULL, SC_LP_PRODUCTION);  /* SC_LP_ERROR for silence. */
  P4EST_GLOBAL_PRODUCTIONF
    ("This is the p4est %dD demo example/steps/%s_step4\n",
     P4EST_DIM, P4EST_STRING);

  /* Create a forest that consists of multiple quadtrees/octrees.
   * This file is compiled for both 2D and 3D: the macro P4_TO_P8 can be
   * checked to execute dimension-dependent code. */
#ifndef P4_TO_P8
  conn = p4est_connectivity_new_unitsquare ();
  /* More complex domains would require a couple changes. */
  /* conn = p4est_connectivity_new_moebius (); */
#else
  conn = p8est_connectivity_new_unitcube ();
  /* For 3D computation, the matrix-vector product would need to be extended
   * by interpolating on both hanging edges and faces. */
  /* More complex domains would require a couple changes. */
  /* conn = p8est_connectivity_new_rotcubes (); */
#endif

  /* Create a forest that is not refined; it consists of the root octant.
   * The p4est_new_ext function can take a startlevel for a load-balanced
   * initial uniform refinement.  Here we refine adaptively instead. */
  startlevel = 0;
  p4est = p4est_new (mpicomm, conn, 0, NULL, NULL);

  /* Refine the forest iteratively, load balancing at each iteration.
   * This is important when starting with an unrefined forest */
  endlevel = 5;
  for (level = startlevel; level < endlevel; ++level) {
    p4est_refine (p4est, 0, refine_fn, NULL);
    /* Refinement has lead to up to 8x more elements; redistribute them. */
    p4est_partition (p4est, 0, NULL);
  }
  if (startlevel < endlevel) {
    /* For finite elements this corner balance is not strictly required.
     * We call balance only once since it's more expensive than both
     * coarsen/refine and partition. */
    p4est_balance (p4est, P4EST_CONNECT_FULL, NULL);
    /* We repartition with coarsening in mind to allow for
     * partition-independent a-posteriori adaptation (not covered in step4). */
    p4est_partition (p4est, 1, NULL);
  }

  /* Write the forest to disk for visualization, one file per processor. */
  p4est_vtk_write_file (p4est, NULL, P4EST_STRING "_step4");

  /* Execute the numerical mathematics part of the example. */
  solve_poisson (p4est);

  /* Destroy the p4est and the connectivity structure. */
  p4est_destroy (p4est);
  p4est_connectivity_destroy (conn);

  /* Verify that allocations internal to p4est and sc do not leak memory.
   * This should be called if sc_init () has been called earlier. */
  sc_finalize ();

  /* This is standard MPI programs.  Without --enable-mpi, this is a dummy. */
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;
}
