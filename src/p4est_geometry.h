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

/** \file p4est_geometry.h Transform from tree-local "reference" coordinate system
 * to global "physical space" coordinates. These are used in \ref p4est_vtk.h to write
 * global coordinate meshes to disk.
 *
 * We provide several example geometries for use. You may also implement your own
 * geometry as you see fit.
 *
 * \note For geometry purposes, each tree has the local coordinate system
 * \f$[0,1]^d\f$. For legacy/\ref p8est compatibility reasons the local
 * coordinates are always represented as a triple abc[3].
 * For a 2D quadtree mesh the local coordinates are abc[0] and abc[1] and
 * the third coordinate abc[2] should be ignored.
 *
 * \ingroup p4est
 */

#ifndef P4EST_GEOMETRY_H
#define P4EST_GEOMETRY_H

#include <p4est_lnodes.h>

SC_EXTERN_C_BEGIN;

/** This object encapsulates a custom geometry transformation. */
typedef struct p4est_geometry p4est_geometry_t;

/** Forward transformation from the tree-local coordinates to physical space.
 *
 * \note The two-dimensional connectivities built into p4est have 3D vertex coordinates
 * that can be used in the transformation if so desired. However, connectivities are
 * not in general required to have vertex coordinate information.
 *
 * \param[in]  geom       associated geometry
 * \param[in]  which_tree tree id inside forest
 * \param[in]  abc        tree-local coordinates: \f$[0,1]^d\f$.
 *                        For 2D meshes abc[2] should never be accessed.
 * \param[out] xyz        cartesian coordinates in physical space after geometry
 */
typedef void        (*p4est_geometry_X_t) (p4est_geometry_t * geom,
                                           p4est_topidx_t which_tree,
                                           const double abc[3],
                                           double xyz[3]);

/** Destructor prototype for a user-allocated \ref p4est_geometry_t.
 * It is invoked by \ref p4est_geometry_destroy.  If the user chooses to
 * reserve the structure statically, there is no need to provide it.
 */
typedef void        (*p4est_geometry_destroy_t) (p4est_geometry_t * geom);

/** Transform a quadrant reference coordinate into the geometry.
 * \param [in] geom     Properly initialized geometry object.
 * \param [in] which_tree   Valid tree number relative to the
 *                          connectivity that is underlying the geometry.
 * \param [in] coords_in    Valid quadrant reference coordinates.
 *                          They must be in [0, P4EST_ROOT_LEN]^2.
 * \param [out] coords_out  Coordinates in the physical geometry.
 */
void                p4est_geometry_transform_coordinates
  (p4est_geometry_t *geom, p4est_topidx_t which_tree,
   p4est_qcoord_t coords_in[2], double coords_out[3]);

/** Encapsulates a custom transformation from tree-local coordinates to
 * user defined physical space.
 *
 * Used in \ref p4est_vtk.h to write global-coordinate meshes.
 *
 * Some internal p4est functions assume that *user points to a
 * \ref p4est_connectivity. However, in general it can be used as the user wishes.
 *
 * This structure can be filled or allocated by the user.
 * p4est will never change its contents.
 */
struct p4est_geometry
{
  const char         *name;     /**< User's choice is arbitrary. */
  void               *user;     /**< User's choice is arbitrary. */
  p4est_geometry_X_t  X;        /**< Coordinate transformation. */
  p4est_geometry_destroy_t destroy;     /**< Destructor called by
                                             \ref p4est_geometry_destroy.  If
                                             NULL, P4EST_FREE is called. */
};

/** Can be used to conveniently destroy a geometry structure.
 * The user is free not to call this function at all if they handle the
 * memory of the \ref p4est_geometry_t in their own way.
 */
void                p4est_geometry_destroy (p4est_geometry_t * geom);

/** Create a geometry structure based on the vertices in a connectivity.
 * The transformation is constructed using bilinear interpolation.
 * \param [in] conn A connectivity with vertex coordinate information.
 *                  We do \a not take ownership and expect this structure to stay alive.
 * \return          Geometry structure; use with \ref p4est_geometry_destroy.
 */
p4est_geometry_t   *p4est_geometry_new_connectivity (p4est_connectivity_t *
                                                     conn);

/** Geometric coordinate transformation for geometry created with
 * \ref p4est_geometry_new_connectivity. This is defined by
 * tri/binlinear interpolation from vertex coordinates.
 *
 * May also be used as a building block in custom geometric coordinate transforms.
 * See for example \ref p4est_geometry_new_sphere2d or \ref p4est_geometry_new_disk2d.
 *
 * \param[in]  geom       associated geometry
 * \param[in]  which_tree tree id inside forest
 * \param[in]  abc        tree-local reference coordinates : [0,1]^3.
 *                        Note: abc[2] is only accessed by the P4_TO_P8 version
 * \param[out] xyz        Cartesian coordinates in physical space after geometry
 *
 * \warning The associated geometry is assumed to have a connectivity
 * as its *user field, and this connectivity is assumed to have vertex
 * information in its *tree_to_vertex field.
 */
void                p4est_geometry_connectivity_X (p4est_geometry_t * geom,
                                                   p4est_topidx_t which_tree,
                                                   const double abc[3],
                                                   double xyz[3]);

/** Create a geometry for mapping the sphere using 2d connectivity icosahedron.
 *
 * \param[in] conn      The result of \ref p4est_connectivity_new_icosahedron.
 * \param[in] R         The radius of the sphere.
 */
p4est_geometry_t   *p4est_geometry_new_icosahedron (p4est_connectivity_t *
                                                    conn, double R);

/** Create a geometry for mapping the annulus.
 *  This a direct adaptation of geometric shell in 3d.
 *
 * \param[in] conn      The result of \ref p4est_connectivity_new_shell2d.
 * \param[in] R1        radius of the inner circle (internal border).
 * \param[in] R2        radius of the outer circle (external border).
 */
p4est_geometry_t   *p4est_geometry_new_shell2d (p4est_connectivity_t * conn,
                                                double R2, double R1);

/**
 * Create disk2d geometry associated to disk2d connectivity.
 *
 * \param[in] conn      The result of \ref p4est_connectivity_new_disk2d.
 * \param[in] R0 radius of the inner circle.
 * \param[in] R1 radius of the outer circle (external border).
 *
 * This geometry is meant to be used with the disk2d connectivity,
 * which is a 5-tree connectivity to map the spherical disk.
 */
p4est_geometry_t   *p4est_geometry_new_disk2d (p4est_connectivity_t * conn,
                                               double R0, double R1);

/**
 * Create sphere geometry associated to cubed connectivity.
 *
 * \param[in] conn The result of \ref p4est_connectivity_new_cubed.
 * \param[in] R radius of the sphere
 *
 * This geometry is meant to be used with the cubed connectivity
 * \ref p4est_connectivity_new_cubed, which is a 6-tree connectivity,
 * to map the sphere.
 */
p4est_geometry_t   *p4est_geometry_new_sphere2d (p4est_connectivity_t * conn,
                                                 double R);

/** Create a geometry for mapping the sphere using 2d connectivity pillow.
 *
 * \param[in] conn      The result of \ref p4est_connectivity_new_pillow.
 * \param[in] R         The radius of the sphere.
 */
p4est_geometry_t   *p4est_geometry_new_pillow (p4est_connectivity_t * conn,
                                               double R);

/** Characterize different mappings of the disk using a 1-tree connectivity.
 *
 * The different mappings correspond to the ones used to produce figure 3.2 in the
 * following publication:
 *
 * "Logically rectangular grids and finite volume methods for PDEs in circular
 * and spherical domains", Calhoun et al, SIAM Review, volume 50, Issue 4, January 2008.
 * https://doi.org/10.1137/060664094
 */
typedef enum
{
  FIG32A = 0,
  FIG32B = 1,
  FIG32C = 2,
  FIG32D = 3,
}
pillow_disk_config_t;

/** Create a geometry for mapping the disk using 2d connectivity unit.
 *
 * See companion routine \ref p8est_geometry_new_pillow_sphere which maps the 3d solid sphere
 * using 1-tree unit connectivity.
 *
 * \param[in] conn      The result of \ref p4est_connectivity_new_unitsquare.
 * \param[in] R         The radius of the disk.
 * \param[in] config    The configuration to identify a mapping variant.
 */
p4est_geometry_t   *p4est_geometry_new_pillow_disk (p4est_connectivity_t *
                                                    conn, double R,
                                                    pillow_disk_config_t
                                                    config);

/** Compute node coordinates for a \ref p4est_lnodes structure.
 * Presently we allow for an lnodes degree of 1 or 2.  Cubic
 * or higher degrees may be transparently enabled in the future.
 *
 * The simple mode assigns one tree reference coordinate to each lnode.
 * This may not be suitable for visualizing periodic connectivities.
 *
 * In a more advanced mode indicated by NULL \c element_coordinates input,
 * the coordinates are made unique by reference location:  If a tree is
 * periodic, for example, its corners reference the same lnode but will
 * generate separate coordinate entries for proper visualization.
 * There will be more coordinates generated than there are lnodes.
 *
 * The coordinate numbers generated by the present version of the function
 * are partition-dependent.  This may be seen as a flaw.  Looking into it.
 *
 * \param [in] p4est    A valid forest structure.
 * \param [in] lnodes   A valid \ref p4est_lnodes structure of degree
 *                      1 or 2.  Higher degrees not presently allowed.
 *                      Must be derived from the \c p4est.
 * \param [in] refloc   Eventually used for cubic and upwards degrees.
 *                      We will expect degree + 1 many values for the
 *                      one-dimensional reference node spacing.  Out of
 *                      these, the indices from 1 to (degree - 1) / 2
 *                      inclusive will be accessed by this function.
 *                      The others default by symmetry considerations.
 * \param [in] geom     May be NULL for generating the tree-reference
 *                      coordinates, or a valid geometry object for
 *                      transforming the reference into mapped space.
 * \param [in,out] coordinates  On input, an array with entries of
 *                      3 double variables each.  Resized in this
 *                      function and populated with coordinate tuples.
 *                      With a NULL geometry, these are in [0, 1]**2.
 *                      Otherwise, they are mapped by the geometry.
 * \param [in,out] element_coordinates  This may be NULL, in which case
 *                      we generate one coordinate tuple for each lnode.
 *                      Otherwise, this must be an array with entries of
 *                      type \ref p4est_locidx_t.  Is resized to the same
 *                      number of entries as \c lnodes->element_nodes.
 *                      Its entries point into the coordinates array.
 *                      The tree index of any given entry is implicit
 *                      in that this array is derived from a p4est,
 *                      where sets of (degree + 1)**2 entries each
 *                      correspond to the forest elements in order.
 */
void                p4est_geometry_coordinates_lnodes
  (p4est_t * p4est, p4est_lnodes_t * lnodes,
   const double *refloc, p4est_geometry_t * geom,
   sc_array_t * coordinates, sc_array_t * element_coordinates);

SC_EXTERN_C_END;

#endif /* !P4EST_GEOMETRY_H */
