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
#if 0                           /* The 3D example is not yet implemented */
#include <p8est_bits.h>
#include <p8est_ghost.h>
#include <p8est_lnodes.h>
#endif
#include <p8est_vtk.h>
#endif

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

/* This example only does something numerical for 2D. */
#ifndef P4_TO_P8

/** Right hand side function for the 2D Poisson problem.
 *
 * This is the negative Laplace operator acting on the function \a uexact.
 * \param [in] vxyz    x, y, and z coordinates in physical space.
 * \return             Scalar function value at vxyz.
 */
static double
func_rhs (const double vxyz[3])
{
  const double        x = vxyz[0];
  const double        y = vxyz[1];

  return 32. * (x * (1. - x) + y * (1. - y));
}

/** Exact solution for the 2D Poisson problem.
 *
 * We pick a function with zero Dirichlet boundary conditions on the unit square.
 * \param [in] vxyz    x, y, and z coordinates in physical space.
 * \return             Scalar function value at vxyz.
 */
static double
func_uexact (const double vxyz[3])
{
  const double        x = vxyz[0];
  const double        y = vxyz[1];

  return 16. * x * (1. - x) * y * (1. - y);
}

/** Determine the boundary status on the unit square.
 * \param [in] p4est    Can be used to access the connectivity.
 * \param [in] tt       The tree number (always zero for the unit square).
 * \param [in] node     The corner node of an element to be examined.
 * \return              True for Dirichlet boundary, false otherwise.
 */
static int
is_boundary_unitsquare (p4est_t * p4est, p4est_topidx_t tt,
                        p4est_quadrant_t * node)
{
  /* For this simple connectivity it is sufficient to check x and y. */
  return (node->x == 0 || node->x == P4EST_ROOT_LEN ||
          node->y == 0 || node->y == P4EST_ROOT_LEN);
}

/** Decode the information from p4est_lnodes_t for a given element.
 *
 * \see p4est_lnodes.h for an in-depth discussion of the encoding.
 * \param [in] face_code         Bit code as defined in p4est_lnodes.h.
 * \param [out] hanging_corner   Undefined if no node is hanging.
 *                               If any node is hanging, this contains
 *                               one integer per corner, which is -1
 *                               for corners that are not hanging,
 *                               and the number of the non-hanging
 *                               corner on the hanging face otherwise.
 * \return true if any node is hanging, false otherwise.
 */
static int
lnodes_decode2 (p4est_lnodes_code_t face_code,
                int hanging_corner[P4EST_CHILDREN])
{
  if (face_code) {
    const int           ones = P4EST_CHILDREN - 1;
    const int           c = (int) (face_code & ones);
    int                 i, h;
    int                 work = (int) (face_code >> P4EST_DIM);

    hanging_corner[c] = hanging_corner[c ^ ones] = -1;
    for (i = 0; i < P4EST_DIM; ++i) {
      h = c ^ ones ^ (1 << i);
      hanging_corner[h] = (work & 1) ? c : -1;
      work >>= 1;
    }
    return 1;
  }
  return 0;
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

/** Compute y := y + a * x
 *
 * \param [in] p4est          The forest is not changed.
 * \param [in] lnodes         The node numbering is not changed.
 * \param [in] a              The scalar.
 * \param [in] x              First node vector.
 * \param [in] y              Second node vector.
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

/** Compute y := x + b * y
 *
 * \param [in] p4est          The forest is not changed.
 * \param [in] lnodes         The node numbering is not changed.
 * \param [in] x              First node vector.
 * \param [in] b              The scalar.
 * \param [in] y              Second node vector.
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

/** zero a vector.
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
          p4est_qcoord_to_vertex (p4est->connectivity, tt,
                                  node.x, node.y, vxyz);

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
  int                 i, j, k;
  int                 q, Q;
  int                 anyhang, hanging_corner[P4EST_CHILDREN];
  int                 isboundary[P4EST_CHILDREN];
  int                 ncontrib, contrib_corner[P4EST_HALF];
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

      /* We are on the 2D unit square.  The Jacobian determinant is h^2.
       * For the stiffness matrix this cancels with the inner derivatives. */
      factor = !stiffness ? pow (.5, 2. * quad->level) : 1.;
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
        for (j = 0; j < P4EST_CHILDREN; ++j) {
          if (hanging_corner[j] >= 0) {
            P4EST_ASSERT (hanging_corner[j] != j &&
                          hanging_corner[j] < P4EST_CHILDREN);
            P4EST_ASSERT (hanging_corner[hanging_corner[j]] == -1);
            /* (3D: process face hanging nodes before edge hanging nodes.) */
            inloc[j] = .5 * (inloc[j] + inloc[hanging_corner[j]]);
          }
        }

        /* Apply element matrix and then the transpose interpolation. */
        for (i = 0; i < P4EST_CHILDREN; ++i) {
          sum = 0.;
          for (j = 0; j < P4EST_CHILDREN; ++j) {
            sum += (*matrix)[i][j] * inloc[j];
          }
          contrib_corner[0] = i;
          if (hanging_corner[i] == -1) {
            ncontrib = 1;
            sum *= factor;
          }
          else {
            ncontrib = 2;
            contrib_corner[1] = hanging_corner[i];
            sum *= .5 * factor;
          }
          for (j = 0; j < ncontrib; ++j) {
            if (!isboundary[contrib_corner[j]]) {
              out[all_lni[contrib_corner[j]]] += sum;
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
set_dirichlet (p4est_lnodes_t * lnodes, int8_t * bc, double *v)
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
 * \param [in] tmp1           Must be allocated, entries are undefined.
 * \param [in] tmp2           Must be allocated, entries are undefined.
 */
static void
test_area (p4est_t * p4est, p4est_lnodes_t * lnodes,
           double (*matrix)[P4EST_CHILDREN][P4EST_CHILDREN],
           double *tmp1, double *tmp2)
{
  const int           nloc = lnodes->num_local_nodes;
  double              dot;
  p4est_locidx_t      lni;

  for (lni = 0; lni < nloc; ++lni) {
    tmp1[lni] = 1.;
  }
  multiply_matrix (p4est, lnodes, NULL, 0, matrix, tmp1, tmp2);
  dot = vector_dot (p4est, lnodes, tmp1, tmp2);
  dot = sqrt (dot);

  P4EST_GLOBAL_PRODUCTIONF ("Area of domain: %g\n", dot);
}

static void
solve_by_cg (p4est_t * p4est, p4est_lnodes_t * lnodes, const int8_t * bc,
             int stiffness,
             double (*matrix)[P4EST_CHILDREN][P4EST_CHILDREN],
             const double *b, double *x)
{
  int                 i;
  double              alpha, beta, rr, pAp, rrnew, rrorig;
  double             *aux[4];
  double             *r, *p, *Ap, *z;
  double              tol = 1.e-6;
  int                 imax = 100;

  for (i = 0; i < 4; ++i) {
    aux[i] = allocate_vector (lnodes);
  }
  r = aux[0];
  p = aux[1];
  Ap = aux[2];
  z = aux[3];

  /* initialize the temporary vector */
  vector_zero (p4est, lnodes, x);
  vector_copy (p4est, lnodes, b, r);
  vector_copy (p4est, lnodes, b, p);
  rr = rrorig = vector_dot (p4est, lnodes, r, r);

  for (i = 0; i < imax; i++) {
    multiply_matrix (p4est, lnodes, bc, stiffness, matrix, p, Ap);
    pAp = vector_dot (p4est, lnodes, p, Ap);
    alpha = rr / pAp;
    vector_axpy (p4est, lnodes, alpha, p, x);
    vector_axpy (p4est, lnodes, -alpha, Ap, r);
    rrnew = vector_dot (p4est, lnodes, r, r);
    if (rrnew <= rrorig * tol) {
      break;
    }
    beta = rrnew / rr;
    vector_xpby (p4est, lnodes, r, beta, p);
    P4EST_GLOBAL_VERBOSEF ("%03d: r'r %e alpha %e beta %e\n", i, rr, alpha,
                           beta);
    rr = rrnew;
  }
  if (i < imax) {
    P4EST_GLOBAL_PRODUCTIONF ("cg converged in %d iterations\n", i);
  }
  else {
    P4EST_GLOBAL_PRODUCTION ("cg did not converge\n");
  }

  for (i = 0; i < 4; ++i) {
    P4EST_FREE (aux[i]);
  }
}

/** Execute the numerical part of the example: Solve Poisson's equation.
 * \param [in] p4est    Solve the PDE with the given mesh refinement.
 */
static void
solve_poisson (p4est_t * p4est)
{
  /** 1D mass matrix on the reference element [0, 1]. */
  static const double m_1d[2][2] = {
    {1 / 3., 1 / 6.},
    {1 / 6., 1 / 3.},
  };
  /** 1D stiffness matrix on the reference element [0, 1]. */
  static const double s_1d[2][2] = {
    {1., -1.},
    {-1., 1.},
  };
  int                 i, j, k, l;
  double              mass_2d[4][4], stiffness_2d[4][4];
  double             *rhs_eval, *uexact_eval;
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
          mass_2d[2 * j + i][2 * l + k] = m_1d[i][k] * m_1d[j][l];
          stiffness_2d[2 * j + i][2 * l + k] =
            m_1d[i][k] * s_1d[j][l] + s_1d[i][k] * m_1d[j][l];
        }
      }
    }
  }

  /* Test mass matrix multiplication by computing the area of the domain. */
  rhs_fe = allocate_vector (lnodes);
  u_fe = allocate_vector (lnodes);
  test_area (p4est, lnodes, &mass_2d, rhs_fe, u_fe);

  /* Interpolate right hand side and exact solution onto mesh nodes. */
  interpolate_functions (p4est, lnodes, &rhs_eval, &uexact_eval, &bc);

  /* Apply mass matrix to create right hand side FE vector. */
  multiply_matrix (p4est, lnodes, bc, 0, &mass_2d, rhs_eval, rhs_fe);
  set_dirichlet (lnodes, bc, rhs_fe);

  /* Run conjugate gradient method with initial value zero. */
  solve_by_cg (p4est, lnodes, bc, 1, &stiffness_2d, rhs_fe, u_fe);

  /* Compute the pointwise difference with the exact vector */
  u_diff = allocate_vector (lnodes);
  vector_copy (p4est, lnodes, u_fe, u_diff);
  vector_axpy (p4est, lnodes, -1., uexact_eval, u_diff);

  /* Compute the L2 difference with the exact vector */
  diff_mass = allocate_vector (lnodes);
  multiply_matrix (p4est, lnodes, bc, 0, &mass_2d, u_diff, diff_mass);
  err2 = vector_dot (p4est, lnodes, diff_mass, u_diff);
  err = sqrt (err2);
  P4EST_GLOBAL_PRODUCTIONF ("|u_fe - u_exact| = %e\n", err);

  /* Free finite element vectors */
  P4EST_FREE (diff_mass);
  P4EST_FREE (u_diff);
  P4EST_FREE (u_fe);
  P4EST_FREE (rhs_fe);
  P4EST_FREE (rhs_eval);
  P4EST_FREE (uexact_eval);
  P4EST_FREE (bc);

  /* We are done with the FE node numbering. */
  p4est_lnodes_destroy (lnodes);
}

#endif /* !P4_TO_P8 */

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
  /* conn = p4est_connectivity_new_moebius (); */
#else
  conn = p8est_connectivity_new_rotcubes ();
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

#ifndef P4_TO_P8
  /* Execute the numerical mathematics part of the example. */
  solve_poisson (p4est);
#else
  P4EST_GLOBAL_PRODUCTION
    ("This example does not do anything in 3D.  Just writes a VTK file.\n");
#endif

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
