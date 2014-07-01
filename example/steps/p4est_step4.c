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
#if 0 /* The 3D example is not yet implemented */
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
  int8_t             *bc, *ndone;
  p4est_topidx_t      tt;       /* Connectivity variables have this type. */
  p4est_locidx_t      k, q, Q;  /* Process-local counters have this type. */
  p4est_locidx_t      lni;      /* Node index relative to this processor. */
  p4est_tree_t       *tree;     /* Pointer to one octree */
  p4est_quadrant_t   *quad, *parent, sp, node;
  sc_array_t         *tquadrants;       /* Quadrant array for one tree */

  rhs = *rhs_eval = P4EST_ALLOC (double, nloc);
  uexact = *uexact_eval = P4EST_ALLOC (double, nloc);
  bc = *pbc = P4EST_ALLOC_ZERO (int8_t, nloc);
  ndone = P4EST_ALLOC_ZERO (int8_t, nloc);

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
        parent = NULL;
      }
      else {
        /* At least one node is hanging.  We need the parent quadrant to
         * find the location of the corresponding non-hanging node. */
        parent = &sp;
        p4est_quadrant_parent (quad, parent);
      }
      for (i = 0; i < P4EST_CHILDREN; ++i) {
        lni = lnodes->element_nodes[P4EST_CHILDREN * k + i];
        P4EST_ASSERT (lni >= 0);
        if (!ndone[lni]) {
          if (anyhang && hanging_corner[i] >= 0) {
            /* This node is hanging; access the referenced node instead. */
            p4est_quadrant_corner_node (parent, i, &node);
          }
          else {
            p4est_quadrant_corner_node (quad, i, &node);
          }
          /* Transform per-tree reference coordinates into physical space. */
          p4est_qcoord_to_vertex (p4est->connectivity, tt,
                                  quad->x, quad->y, vxyz);

          /* TODO: we can now evaluate function values at vxyz. */
          /* TODO: determine boundary status of independent node. */

          /* We are done computing for this node; we only need this once. */
          ndone[lni] = 1;
        }
      }
    }
  }
  P4EST_FREE (ndone);
}

/** Execute the numerical part of the example: Solve Poisson's equation.
 * \param [in] p4est    Solve the PDE with the given mesh refinement.
 */
static void
solve_poisson (p4est_t * p4est)
{
  int8_t             *bc;
  double             *rhs_eval, *uexact_eval;
  p4est_ghost_t      *ghost;
  p4est_lnodes_t     *lnodes;

  /* Create the ghost layer to learn about parallel neighbors. */
  ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FULL);

  /* Create a node numbering for continuous linear finite elements. */
  lnodes = p4est_lnodes_new (p4est, ghost, 1);

  /* Destroy the ghost structure -- no longer needed after node creation. */
  p4est_ghost_destroy (ghost);
  ghost = NULL;

  /* Interpolate right hand side and exact solution onto mesh nodes. */
  interpolate_functions (p4est, lnodes, &rhs_eval, &uexact_eval, &bc);

  /* Free finite element vectors */
  P4EST_FREE (rhs_eval);
  P4EST_FREE (uexact_eval);
  P4EST_FREE (bc);

  /* We are done with the FE node numbering. */
  p4est_lnodes_destroy (lnodes);
  lnodes = NULL;
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
  p4est_init (NULL, SC_LP_PRODUCTION);
  P4EST_GLOBAL_PRODUCTIONF
    ("This is the p4est %dD demo example/steps/%s_step4\n",
     P4EST_DIM, P4EST_STRING);

  /* Create a forest that consists of multiple quadtrees/octrees.
   * This file is compiled for both 2D and 3D: the macro P4_TO_P8 can be
   * checked to execute dimension-dependent code. */
#ifndef P4_TO_P8
  conn = p4est_connectivity_new_moebius ();
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
