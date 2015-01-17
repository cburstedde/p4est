/*
  This file is part of p4est.
  p4est is a C library to manage a collection (a forest) of multiple
  connected adaptive quadtrees or octrees in parallel.

  Copyright (C) 2013 The University of Texas System
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

/** \file p6est.h
 *
 * A hybrid 2D+1D AMR extension.
 *
 * \ingroup p6est
 */

/** \defgroup p6est p6est
 *
 * A hybrid 2D+1D AMR extension.
 *
 * To include this component of the p4est library, configure p4est with the
 * --enable-p6est option given.  This module provides a specific kind of
 * anisotropic adaptive mesh refinement for 3D meshes: it organizes the
 * hexahedral cells into columns.  Each column has the footprint of a 2D p4est
 * quadrant.  The cells within a column can be individually refined
 * vertically, and a whole column of cells can be refined horizontally.  When
 * the forest is partitioned, each column is assigned to a single MPI process,
 * i.e., a column cannot be split between processes.  Most of the main 2D / 3D
 * interface is available for p6est: refinement and coarsening, balance, ghost
 * layers, and i/o.  Finite element nodes can also be created using
 * p6est_lnodes_new(): this creates nodes in the same data structure as
 * p8est_lnodes_new().
 */

#ifndef P6EST_H
#define P6EST_H

/* 2+1D refinement is based on the 2D p4est datatypes */
#include <p4est.h>
/* We need p8est_connect_type_t typedef from p8est_connectivity */
#include <p8est_connectivity.h>

SC_EXTERN_C_BEGIN;

/** This structure holds the 2D+1D inter-tree connectivity information.
 * It is essentially a wrapper of the 2D p4est_connecitivity_t datatype, with
 * some additional information about how the third dimension is embedded.
 */
typedef struct p6est_connectivity
{
  p4est_connectivity_t *conn4;  /**< the 2D connecitvity; owned; vertices
                                  interpreted as the vertices of the bottom of
                                  the sheet */
  double             *top_vertices;     /**< if NULL, uniform vertical profile,
                                           otherwise the vertices of the top of
                                           the sheet: should be the same size
                                           as \a conn4->tree_to_vertex; owned. */
  double              height[3];        /**< if \a top_vertices == NULL, this gives the
                                           offset from the bottom of the sheet to
                                           the top */
}
p6est_connectivity_t;

/** Create a p6est_connectivity_t from a p4est_connectivity_t.  All fields
 * are copied, so all inputs can be safey destroyed.
 *
 * \param[in] conn4         the 2D connectivity
 * \param[in] top_vertices  if NULL, then the sheet has a uniform vertical
 *                          profile; otherwise, \a top_vertices gives teh
 *                          vertices of the top of the sheet; should be the
 *                          same size as \a conn4->tree_to_vertex
 * \param[in] height        if \a top_vertices == NULL, then this gives the
 *                          offset fro the bottom of the sheet to the top.
 *
 * \return the 2D+1D connectivity information.
 */
p6est_connectivity_t *p6est_connectivity_new (p4est_connectivity_t * conn4,
                                              double *top_vertices,
                                              double height[3]);

/** Destroy a p6est_connectivity structure */
void                p6est_connectivity_destroy (p6est_connectivity_t * conn);

/** Get the vertices of the corners of a tree.
 *
 * \param[in]  conn         the 2D+1D connectivity structure
 * \param[in]  which_tree   a tree in the forest
 * \param[out] vertices     the coordinates of the corners of the tree
 */
void                p6est_tree_get_vertices (p6est_connectivity_t * conn,
                                             p4est_topidx_t which_tree,
                                             double vertices[24]);

/** Transform a quadrant coordinate into the space spanned by tree vertices.
 *
 * \param [in] connectivity     Connectivity must provide the vertices.
 * \param [in] treeid           Identify the tree that contains x, y.
 * \param [in] x, y             Quadrant coordinates relative to treeid.
 * \param [out] vxy             Transformed coordinates in vertex space.
 */
void                p6est_qcoord_to_vertex (p6est_connectivity_t *
                                            connectivity,
                                            p4est_topidx_t treeid,
                                            p4est_qcoord_t x,
                                            p4est_qcoord_t y,
                                            p4est_qcoord_t z, double vxyz[3]);

/** A 1D quadrant datatype: this is used to encode a "layer" of a column in
 * the 2D+1D AMR scheme.
 */
typedef struct p2est_quadrant
{
  p4est_qcoord_t      z;                /**< vertical coordinate */
  int8_t              level,            /**< level of refinement */
                      pad8;             /**< padding */
  int16_t             pad16;            /**< padding */
  union p6est_quadrant_data
  {
    void               *user_data;      /**< never changed by p4est */
    long                user_long;      /**< never changed by p4est */
    int                 user_int;       /**< never changed by p4est */
    p4est_topidx_t      which_tree;     /**< the tree containing the quadrant */
    struct
    {
      p4est_topidx_t      which_tree;
      int                 owner_rank;
    }
    piggy1; /**< of ghost layer, store the tree and owner rank */
    struct
    {
      p4est_topidx_t      which_tree;
      p4est_topidx_t      from_tree;
    }
    piggy2; /**< of transformed layers, store the original tree and the
                 target tree */
    struct
    {
      p4est_topidx_t      which_tree;
      p4est_locidx_t      local_num;
    }
    piggy3; /**< of ghost layers, store the tree and index in the owner's
                 numbering */
  }
  p; /**< a union of additional data attached to a layer */
}
p2est_quadrant_t;

/** The p6est forest datatype */
typedef struct p6est
{
  sc_MPI_Comm         mpicomm;          /**< MPI communicator */
  int                 mpisize,          /**< number of MPI processes */
                      mpirank;          /**< this process's MPI rank */
  size_t              data_size;        /**< size of per-quadrant p.user_data
                     (see p2est_quadrant_t::p2est_quadrant_data::user_data) */
  void               *user_pointer;     /**< convenience pointer for users,
                                             never touched by p4est */
  p6est_connectivity_t *connectivity;   /**< topology of sheet, not owned. */
  p4est_t            *columns;  /**< 2D description of column layout
                                     built from \a connectivity */
  sc_array_t         *layers;   /**< single array that stores
                                     p2est_quadrant_t layers within columns */
  sc_mempool_t       *user_data_pool;   /**< memory allocator for user data */
                                        /* WARNING: This is NULL if data size
                                         *          equals zero.  */
  sc_mempool_t       *layer_pool;       /**< memory allocator
                                             for temporary layers */
  p4est_gloidx_t     *global_first_layer; /**< first global quadrant index for
                                               each process and 1 beyond */
  p4est_qcoord_t      root_len; /**< height of the domain */
}
p6est_t;

/** Callback function prototype to initialize the layers's user data.
 *
 * \param[in] p6est        the forest
 * \param[in] which_tree   the tree in the forest
 * \param[in] column       the column in the tree in the forest
 * \param[in] layer        the layer in the column in the tree in the
 *                         forest, whose \a user_data is to be initialized
 */
typedef void        (*p6est_init_t) (p6est_t * p6est,
                                     p4est_topidx_t which_tree,
                                     p4est_quadrant_t * column,
                                     p2est_quadrant_t * layer);

/** Callback function prototype to transfer information from outgoing layers
 * to incoming layers.
 *
 * This is used by extended routines when the layers of an existing, valid
 * p6est are changed.  The callback allows the user to make changes to newly
 * initialized layers before the layers that they replace are destroyed.
 *
 * \param [in] num_outcolumns  The number of columns that contain the outgoing
 *                             layers: will be either 1 or 4.
 * \param [in] num_outlayers   The number of outgoing layers: will be either 1
 *                             (a single layer is being refined), 2 (two
 *                             layers are being vertically coarsened), or 4
 *                             (four layers are being horizontally coarsened).
 * \param [in] outcolumns      The columns of the outgoing layers
 * \param [in] outlayers       The outgoing layers: after the callback, the
 *                             user_data, if \a p6est->data_size is nonzero,
 *                             will be destroyed.
 * \param [in] num_incolumns   The number of columns that contain the outgoing
 *                             layers: will be either 1 or 4.
 * \param [in] num_inlayers    The number of incoming layers: will be either 1
 *                             (coarsening), 2 (vertical refinement), or 4
 *                             (horizontal refinement)
 * \param [in] incolumns       The columns of the incoming layers
 * \param [in,out] inlayers    The incoming layers: prior to the callback,
 *                             the user_data, if \a p6est->data_size is nonzero,
 *                             is allocated, and the p6est_init_t callback,
 *                             if it has been provided, will be called.
 */
typedef void        (*p6est_replace_t) (p6est_t * p6est,
                                        p4est_topidx_t which_tree,
                                        int num_outcolumns,
                                        int num_outlayers,
                                        p4est_quadrant_t * outcolumns[],
                                        p2est_quadrant_t * outlayers[],
                                        int num_incolumns,
                                        int num_inlayers,
                                        p4est_quadrant_t * incolumns[],
                                        p2est_quadrant_t * inlayers[]);

/** Callback function prototype to decide whether to horizontally refine a
 * column, i.e., horizontally refine all of the layers in the column.
 * \return nonzero if the layer shall be refined.
 */
typedef int         (*p6est_refine_column_t) (p6est_t * p6est,
                                              p4est_topidx_t which_tree,
                                              p4est_quadrant_t * column);

/** Callback function prototype to decide whether to vertically refine a
 * layer.
 * \return nonzero if the layer shall be refined.
 */
typedef int         (*p6est_refine_layer_t) (p6est_t * p6est,
                                             p4est_topidx_t which_tree,
                                             p4est_quadrant_t * column,
                                             p2est_quadrant_t * layer);

/** Callback function prototype to decide for horizontal coarsening.
 * \param [in] columns      Pointers to 4 sibling columns.
 * \return nonzero if the columns shall be replaced with their parent.
 */
typedef int         (*p6est_coarsen_column_t) (p6est_t * p6est,
                                               p4est_topidx_t which_tree,
                                               p4est_quadrant_t * columns[]);

/** Callback function prototype to decide for vertical coarsening.
 * \param [in] layers      Pointers to 2 vertical sibling layers.
 * \return nonzero if the layers shall be replaced with their parent.
 */
typedef int         (*p6est_coarsen_layer_t) (p6est_t * p6est,
                                              p4est_topidx_t which_tree,
                                              p4est_quadrant_t * column,
                                              p2est_quadrant_t * layers[]);

/** Callback function prototype to calculate weights for partitioning.
 * \return a 32bit integer >= 0 as the quadrant weight.
 * \note    Global sum of weights must fit into a 64bit integer.
 */
typedef int         (*p6est_weight_t) (p6est_t * p6est,
                                       p4est_topidx_t which_tree,
                                       p4est_quadrant_t * column,
                                       p2est_quadrant_t * layer);

extern void        *P2EST_DATA_UNINITIALIZED;

/** set statically allocated quadrant to defined values */
#define P2EST_QUADRANT_INIT(q) \
  ((void) memset ((q), -1, sizeof (p2est_quadrant_t)))

/** Create a new forest.
 *
 * The new forest consists of equi-partitioned root quadrants.
 * When there are more processors than trees, some processors are empty.
 *
 * \param [in] mpicomm       A valid MPI communicator.
 * \param [in] connectivity  This is the connectivity information that
 *                           the forest is built with.  Note the p6est
 *                           does not take ownership of the memory.
 * \param [in] data_size     This is the size of data for each quadrant which
 *                           can be zero.  Then user_data_pool is set to NULL.
 * \param [in] init_fn       Callback function to initialize the user_data
 *                           which is already allocated automatically.
 * \param [in] user_pointer  Assign to the user_pointer member of the p6est
 *                           before init_fn is called the first time.
 *
 * \return This returns a valid forest.
 *
 * \note The connectivity structure must not be destroyed
 *       during the lifetime of this forest.
 */
p6est_t            *p6est_new (sc_MPI_Comm mpicomm,
                               p6est_connectivity_t * connectivity,
                               size_t data_size,
                               p6est_init_t init_fn, void *user_pointer);

/** Create a new forest from an already created p4est that represents
 * columns.
 *
 * \param [in] p4est         A valid p4est.  A deep copy will be created, so
 *                           this can be destroyed without affectin the new
 *                           p6est object.
 * \param [in] top_vertices  the same as in p6est_conectivity_new()
 * \param [in] height        the same as in p6est_conectivity_new()
 * \param [in] min_zlevel    the same as in p6est_new()
 * \param [in] data_size     the same as in p6est_new()
 * \param [in] init_fn       the same as in p6est_new()
 * \param [in] user_pointer  the same as in p6est_new()
 *
 * \return This returns a valid forest.  The user must destroy the
 * connectivity for the new p6est independently.
 */
p6est_t            *p6est_new_from_p4est (p4est_t * p4est,
                                          double *top_vertices,
                                          double height[3], int min_zlevel,
                                          size_t data_size,
                                          p6est_init_t init_fn,
                                          void *user_pointer);

/** Destroy a p6est.
 *
 * \note The connectivity structure is not destroyed with the p6est.
 */
void                p6est_destroy (p6est_t * p6est);

/** Make a deep copy of a p6est.
 * The connectivity is not duplicated.
 * Copying of quadrant user data is optional.
 * If old and new data sizes are 0, the user_data field is copied regardless.
 *
 * \param [in]  copy_data  If true, data are copied.
 *                         If false, data_size is set to 0.
 * \return  Returns a valid p6est that does not depend on the input.
 */
p6est_t            *p6est_copy (p6est_t * input, int copy_data);

/** Reset user pointer and element data.
 * When the data size is changed the quadrant data is freed and allocated.
 * The initialization callback is invoked on each quadrant.
 * Old user_data content is disregarded.
 *
 * \param [in] data_size     This is the size of data for each quadrant which
 *                           can be zero.  Then user_data_pool is set to NULL.
 * \param [in] init_fn       Callback function to initialize the user_data
 *                           which is already allocated automatically.
 * \param [in] user_pointer  Assign to the user_pointer member of the p6est
 *                           before init_fn is called the first time.
 */
void                p6est_reset_data (p6est_t * p6est, size_t data_size,
                                      p6est_init_t init_fn,
                                      void *user_pointer);

/** Refine the columns of a sheet.
 *
 * \param [in,out] p6est The forest is changed in place.
 * \param [in] refine_recursive Boolean to decide on recursive refinement.
 * \param [in] refine_fn Callback function that must return true if a column
 *                       shall be refined into smaller columns.  If
 *                       refine_recursive is true, refine_fn is called for
 *                       every existing and newly created column.
 *                       Otherwise, it is called for every existing column.
 *                       It is possible that a refinement request made by the
 *                       callback is ignored.  To catch this case, you can
 *                       examine whether init_fn gets called, or use
 *                       p6est_refine_columns_ext in p6est_extended.h and examine
 *                       whether replace_fn gets called.
 * \param [in] init_fn   Callback function to initialize the user_data of newly
 *                       created layers within columns, which are already
 *                       allocated.  This function pointer may be NULL.
 */
void                p6est_refine_columns (p6est_t * p6est,
                                          int refine_recursive,
                                          p6est_refine_column_t refine_fn,
                                          p6est_init_t init_fn);

/** Refine the layers within the columns of a sheet.
 *
 * \param [in,out] p6est The forest is changed in place.
 * \param [in] refine_recursive Boolean to decide on recursive refinement.
 * \param [in] refine_fn Callback function that must return true if a layer
 *                       shall be refined into smaller layers.  If
 *                       refine_recursive is true, refine_fn is called for
 *                       every existing and newly created layer.
 *                       Otherwise, it is called for every existing layer.
 *                       It is possible that a refinement request made by the
 *                       callback is ignored.  To catch this case, you can
 *                       examine whether init_fn gets called, or use
 *                       p6est_refine_layers_ext in p6est_extended.h and examine
 *                       whether replace_fn gets called.
 * \param [in] init_fn   Callback function to initialize the user_data of newly
 *                       created layers, which are already allocated.  This
 *                       function pointer may be NULL.
 */
void                p6est_refine_layers (p6est_t * p6est,
                                         int refine_recursive,
                                         p6est_refine_layer_t refine_fn,
                                         p6est_init_t init_fn);

/** Coarsen the columns of a sheet.
 *
 * \param [in,out] p6est  The forest is changed in place.
 * \param [in] coarsen_recursive Boolean to decide on recursive coarsening.
 * \param [in] coarsen_fn Callback function that returns true if a
 *                        family of columns shall be coarsened
 * \param [in] init_fn    Callback function to initialize the user_data
 *                        which is already allocated automatically.
 */
void                p6est_coarsen_columns (p6est_t * p6est,
                                           int coarsen_recursive,
                                           p6est_coarsen_column_t coarsen_fn,
                                           p6est_init_t init_fn);

/** Coarsen the layers of a sheet.
 *
 * \param [in,out] p6est  The forest is changed in place.
 * \param [in] coarsen_recursive Boolean to decide on recursive coarsening.
 * \param [in] coarsen_fn Callback function that returns true if a
 *                        family of layers shall be coarsened
 * \param [in] init_fn    Callback function to initialize the user_data
 *                        which is already allocated automatically.
 */
void                p6est_coarsen_layers (p6est_t * p6est,
                                          int coarsen_recursive,
                                          p6est_coarsen_layer_t coarsen_fn,
                                          p6est_init_t init_fn);

/** Balance a forest.
 *
 * \param [in] p6est     The p6est to be worked on.
 * \param [in] btype     Balance type (face, corner or default, full).
 * \param [in] init_fn   Callback function to initialize the user_data
 *                       which is already allocated automatically.
 */
void                p6est_balance (p6est_t * p6est,
                                   p8est_connect_type_t btype,
                                   p6est_init_t init_fn);

typedef enum
{
  P6EST_COMM_PARTITION = 1,
  P6EST_COMM_GHOST,
  P6EST_COMM_BALANCE
}
p6est_comm_tag_t;

/** Equally partition the forest.
 *
 * The forest will be partitioned between processors where they each
 * have an approximately equal number of quadrants.
 *
 * Note that \a p6est->layers and \a p6est->global_first_layers may change
 * during this call.  Address pointers referencing these objects from before
 * \a p6est_partition is called become invalid.
 *
 * \param [in,out] p6est      The forest that will be partitioned.
 * \param [in]     weight_fn  A weighting function or NULL
 *                            for uniform partitioning.
 */
p4est_gloidx_t      p6est_partition (p6est_t * p6est,
                                     p6est_weight_t weight_fn);

/** Compute the checksum for a forest.
 * Based on quadrant arrays only. It is independent of partition and mpisize.
 * \return  Returns the checksum on processor 0 only. 0 on other processors.
 */
unsigned            p6est_checksum (p6est_t * p6est);

/** Save the complete connectivity/p6est data to disk.  This is a collective
 *
 * operation that all MPI processes need to call.  All processes write
 * into the same file, so the filename given needs to be identical over
 * all parallel invocations.
 * \param [in] filename    Name of the file to write.
 * \param [in] p6est       Valid forest structure.
 * \param [in] save_data   If true, the element data is saved.
 *                         Otherwise, a data size of 0 is saved.
 * \note            Aborts on file errors.
 */
void                p6est_save (const char *filename, p6est_t * p6est,
                                int save_data);

/** Load the complete connectivity/p6est structure from disk.
 *
 * \param [in] filename         Name of the file to read.
 * \param [in] mpicomm          A valid MPI communicator.
 * \param [in] data_size        Size of data for each quadrant which can be
 *                              zero.  Then user_data_pool is set to NULL.
 *                              If data_size is zero, load_data is ignored.
 * \param [in] load_data        If true, the element data is loaded.  This is
 *                              only permitted if the saved data size matches.
 *                              If false, the stored data size is ignored.
 * \param [in] user_pointer     Assign to the user_pointer member of the p6est
 *                              before init_fn is called the first time.
 * \param [out] connectivity    Connectivity must be destroyed separately.
 * \return          Returns a valid forest structure. A pointer to a valid
 *                  connectivity structure is returned through the last
 *                  argument.
 * \note            Aborts on file errors or invalid file contents.
 */
p6est_t            *p6est_load (const char *filename, sc_MPI_Comm mpicomm,
                                size_t data_size, int load_data,
                                void *user_pointer,
                                p6est_connectivity_t ** connectivity);

/** Return a pointer to a quadrant array element indexed by a size_t. */
/*@unused@*/
static inline p2est_quadrant_t *
p2est_quadrant_array_index (sc_array_t * array, size_t it)
{
  P4EST_ASSERT (array->elem_size == sizeof (p2est_quadrant_t));
  P4EST_ASSERT (it < array->elem_count);

  return (p2est_quadrant_t *) (array->array + sizeof (p2est_quadrant_t) * it);
}

/** Call sc_array_push for a quadrant array. */
/*@unused@*/
static inline p2est_quadrant_t *
p2est_quadrant_array_push (sc_array_t * array)
{
  P4EST_ASSERT (array->elem_size == sizeof (p2est_quadrant_t));

  return (p2est_quadrant_t *) sc_array_push (array);
}

/** Call sc_mempool_alloc for a mempool creating quadrants. */
/*@unused@*/
static inline p2est_quadrant_t *
p2est_quadrant_mempool_alloc (sc_mempool_t * mempool)
{
  P4EST_ASSERT (mempool->elem_size == sizeof (p2est_quadrant_t));

  return (p2est_quadrant_t *) sc_mempool_alloc (mempool);
}

/** Call sc_list pop for a quadrant array. */
/*@unused@*/
static inline p2est_quadrant_t *
p2est_quadrant_list_pop (sc_list_t * list)
{
  return (p2est_quadrant_t *) sc_list_pop (list);
}

#define P6EST_COLUMN_GET_RANGE(q,f,l)                \
  do {                                               \
    *(f) = (size_t) (q)->p.piggy3.local_num;         \
    *(l) = *(f) + (size_t) (q)->p.piggy3.which_tree; \
  } while (0);

#define P6EST_COLUMN_SET_RANGE(q,f,l)                        \
  do {                                                       \
    (q)->p.piggy3.local_num = (p4est_locidx_t) (f);          \
    (q)->p.piggy3.which_tree = (p4est_topidx_t) ((l) - (f)); \
  } while (0);

/*@unused@*/
static inline void
p6est_layer_init_data (p6est_t * p6est, p4est_topidx_t which_tree,
                       p4est_quadrant_t * column,
                       p2est_quadrant_t * layer, p6est_init_t init_fn)
{
  if (p6est->data_size > 0) {
    layer->p.user_data = sc_mempool_alloc (p6est->user_data_pool);
  }
  else {
    layer->p.user_data = NULL;
  }
  if (init_fn != NULL) {
    init_fn (p6est, which_tree, column, layer);
  }
}

/*@unused@*/
static inline void
p6est_layer_free_data (p6est_t * p6est, p2est_quadrant_t * layer)
{
  if (p6est->data_size > 0) {
    sc_mempool_free (p6est->user_data_pool, layer->p.user_data);
  }
  layer->p.user_data = NULL;
}

void                p6est_compress_columns (p6est_t * p6est);
void                p6est_update_offsets (p6est_t * p6est);

SC_EXTERN_C_END;

#endif /* !P6EST_H */
