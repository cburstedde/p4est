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

/** \file p8est_geometry.h Transform from tree-local "reference" coordinate system
 * to global "physical space" coordinates. These are used in \ref p8est_vtk.h to write
 * global coordinate meshes to disk.
 *
 * We provide several example geometries for use. You may also implement your own
 * geometry as you see fit.
 *
 * \note For geometry purposes, each tree has the local coordinate system
 * \f$[0,1]^3\f$.
 *
 * \ingroup p8est
 */

#ifndef P8EST_GEOMETRY_H
#define P8EST_GEOMETRY_H

#include <p8est_lnodes.h>

SC_EXTERN_C_BEGIN;

/** This object encapsulates a custom geometry transformation. */
typedef struct p8est_geometry p8est_geometry_t;

/** Forward transformation from the reference unit square to physical space.
 * The physical space "xyz" is user-defined, currently used for VTK output.
 */
typedef void        (*p8est_geometry_X_t) (p8est_geometry_t * geom,
                                           p4est_topidx_t which_tree,
                                           const double abc[3],
                                           double xyz[3]);

/** Destructor prototype for a user-allocated \a p8est_geometry_t.
 * It is invoked by p8est_geometry_destroy.  If the user chooses to
 * reserve the structure statically, there is no need to provide it.
 */
typedef void        (*p8est_geometry_destroy_t) (p8est_geometry_t * geom);

/** Transform a quadrant reference coordinate into the geometry.
 * \param [in] geom     Properly initialized geometry object.
 * \param [in] which_tree   Valid tree number relative to the
 *                          connectivity that is underlying the geometry.
 * \param [in] coords_in    Valid quadrant reference coordinates.
 *                          They must be in [0, P4EST_ROOT_LEN]^3.
 * \param [out] coords_out  Coordinates in the physical geometry.
 */
void                p8est_geometry_transform_coordinates
  (p8est_geometry_t *geom, p4est_topidx_t which_tree,
   p4est_qcoord_t coords_in[3], double coords_out[3]);

/** This structure can be created by the user,
 * p4est will never change its contents.
 */
struct p8est_geometry
{
  const char         *name;     /**< User's choice is arbitrary. */
  void               *user;     /**< User's choice is arbitrary. */
  p8est_geometry_X_t  X;        /**< Coordinate transformation. */
  p8est_geometry_destroy_t destroy;     /**< Destructor called by
                                             p8est_geometry_destroy.  If
                                             NULL, P4EST_FREE is called. */
};

/** Can be used to conveniently destroy a geometry structure.
 * The user is free not to call this function at all if they handle the
 * memory of the \ref p8est_geometry_t in their own way.
 */
void                p8est_geometry_destroy (p8est_geometry_t * geom);

/** Create a geometry structure based on the vertices in a connectivity.
 * The transformation is constructed using trilinear interpolation.
 * \param [in] conn A p8est_connectivity_t with valid vertices.  We do NOT
 *                  take ownership and expect this structure to stay alive.
 * \return          Geometry structure; use with \ref p4est_geometry_destroy.
 */
p8est_geometry_t   *p8est_geometry_new_connectivity (p8est_connectivity_t *
                                                     conn);

/** Geometric coordinate transformation for geometry created with
 * \ref p8est_geometry_new_connectivity. This is defined by
 * tri/binlinear interpolation from vertex coordinates.
 *
 * May also be used as a building block in custom geometric coordinate transforms.
 * See for example \ref p8est_geometry_new_shell or \ref p8est_geometry_new_sphere.
 *
 * \param[in]  geom       associated geometry
 * \param[in]  which_tree tree id inside forest
 * \param[in]  abc        tree-local reference coordinates : [0,1]^3.
 * \param[out] xyz        Cartesian coordinates in physical space after geometry
 *
 * \warning The associated geometry is assumed to have a connectivity
 * as its *user field, and this connectivity is assumed to have vertex
 * information in its *tree_to_vertex field.
 */
void                p8est_geometry_connectivity_X (p8est_geometry_t * geom,
                                                   p4est_topidx_t which_tree,
                                                   const double abc[3],
                                                   double xyz[3]);

/** Create a geometry structure for the spherical shell of 2 trees.
 * \param [in] conn Result of p8est_connectivity_new_pillow.
 *                  We do NOT take ownership and expect it to stay alive.
 * \param [in] R2   The outer radius of the shell.
 * \param [in] R1   The inner radius of the shell.
 * \return          Geometry structure; use with \ref p4est_geometry_destroy.
 *
 * \note this coordinate transformation is describe in "Logically rectangular
 * grids and finite volume methods for PDEs in circular and spherical domains",
 * Calhoun et al., https://doi.org/10.1137/060664094
 */
p8est_geometry_t   *p8est_geometry_new_pillow (p8est_connectivity_t * conn,
                                               double R2, double R1);

/** Characterize different mappings of the solid sphere using a 1-tree connectivity.
 *
 * The different mappings correspond to the ones used to produce figure 5.2 in the
 * following publication:
 *
 * "Logically rectangular grids and finite volume methods for PDEs in circular
 * and spherical domains", Calhoun et al, SIAM Review, volume 50, Issue 4, January 2008.
 * https://doi.org/10.1137/060664094
 */
typedef enum
{
  FIG52B = 0,
  FIG52C = 1
}
pillow_sphere_config_t;

/** Create a geometry for mapping the solid sphere using the 1-tree unit connectivity.
 *
 * See companion routine \ref p4est_geometry_new_pillow_disk which maps the 2d disk
 * using 1-tree unit connectivity.
 *
 * \param[in] conn      The result of \ref p8est_connectivity_new_unitcube.
 * \param[in] R         The radius of the solid sphere.
 * \param[in] config    The configuration to identify a mapping variant.
 */
p8est_geometry_t   *p8est_geometry_new_pillow_sphere (p8est_connectivity_t *
                                                      conn, double R,
                                                      pillow_sphere_config_t
                                                      config);

/** Create a geometry structure for the spherical shell of 24 trees.
 * \param [in] conn Result of p8est_connectivity_new_shell or equivalent.
 *                  We do NOT take ownership and expect it to stay alive.
 * \param [in] R2   The outer radius of the shell.
 * \param [in] R1   The inner radius of the shell.
 * \return          Geometry structure; use with \ref p4est_geometry_destroy.
 */
p8est_geometry_t   *p8est_geometry_new_shell (p8est_connectivity_t * conn,
                                              double R2, double R1);

/** Create a geometry structure for the solid sphere of 13 trees.
 * \param [in] conn Result of p8est_connectivity_new_sphere or equivalent.
 *                  We do NOT take ownership and expect it to stay alive.
 * \param [in] R2   The outer radius of the sphere.
 * \param [in] R1   The outer radius of the inner shell.
 * \param [in] R0   The inner radius of the inner shell.
 * \return          Geometry structure; use with \ref p4est_geometry_destroy.
 */
p8est_geometry_t   *p8est_geometry_new_sphere (p8est_connectivity_t * conn,
                                               double R2, double R1,
                                               double R0);

/** Create a geometry structure for the torus.
 *
 * This geometry maps a revolution torus, obtained using
 * \ref p8est_connectivity_new_torus.
 *
 * The torus is divided into into segments around the revolution axis,
 * each segment is made of 5 trees; so here we provided the geometric
 * transformation in a piecewise manner for each tree of the connectivity.
 *
 * \param [in] conn Result of p8est_connectivity_new_torus or equivalent.
 *                  We do NOT take ownership and expect it to stay alive.
 *
 * \param [in] R0   The inner radius of the 2d disk slice as cross section.
 * \param [in] R1   The outer radius of the 2d disk slice as cross section.
 * \param [in] R2   The radius of the center circle of the torus.
 *                  The outer radius of the torus is thus \a R1 + \a R2.
 * \return          Geometry structure; use with \ref p4est_geometry_destroy.
 */
p8est_geometry_t   *p8est_geometry_new_torus (p8est_connectivity_t * conn,
                                              double R0, double R1,
                                              double R2);

/** Compute node coordinates for a \ref p8est_lnodes structure.
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
 * \param [in] p8est    A valid forest structure.
 * \param [in] lnodes   A valid \ref p4est_lnodes structure of degree
 *                      1 or 2.  Higher degrees not presently allowed.
 *                      Must be derived from the \c p8est.
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
 *                      With a NULL geometry, these are in [0, 1]**3.
 *                      Otherwise, they are mapped by the geometry.
 * \param [in,out] element_coordinates  This may be NULL, in which case
 *                      we generate one coordinate tuple for each lnode.
 *                      Otherwise, this must be an array with entries of
 *                      type \ref p4est_locidx_t.  Is resized to the same
 *                      number of entries as \c lnodes->element_nodes.
 *                      Its entries point into the coordinates array.
 *                      The tree index of any given entry is implicit
 *                      in that this array is derived from a p4est,
 *                      where sets of (degree + 1)**3 entries each
 *                      correspond to the forest elements in order.
 */
void                p8est_geometry_coordinates_lnodes
  (p8est_t * p8est, p8est_lnodes_t * lnodes,
   const double *refloc, p8est_geometry_t * geom,
   sc_array_t * coordinates, sc_array_t * element_coordinates);

SC_EXTERN_C_END;

#endif /* !P8EST_GEOMETRY_H */
