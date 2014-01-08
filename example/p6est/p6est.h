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

#ifndef P6EST_H
#define P6EST_H

/* 2+1D refinement is based on the 2D p4est datatypes */
#include <p4est.h>
/* We need p8est_connect_type_t typedef from p8est_connectivity */
#include <p8est_connectivity.h>

SC_EXTERN_C_BEGIN;

typedef struct p6est_connectivity_t
{
  p4est_connectivity_t *conn4; /* owned: vertices interpreted as the vertices
                                  of the bottom of the sheet */
  double               *top_to_vertex; /* if NULL, uniform vertical profile,
                                          otherwise the vertices of the top of
                                          the sheet: should be the same size
                                          as \a conn4->tree_to_vertex; owned. */
  double                height[3]; /* if top_to_vertex == NULL, this gives the
                                      offset from the bottom of the sheet to
                                      the top */
}
p6est_connectivity_t;

/** Create a p6est_connectivity_t from a  p4est_connectivity_t.  All fields
 * are copied.
 */
p6est_connectivity_t *p6est_connectivity_new (p4est_connectivity_t *conn4,
                                              double *top_to_vertex,
                                              double height[3]);

void                p6est_connectivity_destroy (p6est_connectivity_t *conn);

void                p6est_tree_get_vertices (p6est_connectivity_t *conn,
                                             p4est_topidx_t which_tree,
                                             double vertices[24]);

/** Transform a quadrant coordinate into the space spanned by tree vertices.
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

typedef struct p2est_quadrant
{
  p4est_qcoord_t      z;
  int8_t              level, pad8;
  int16_t             pad16;
  union p6est_quadrant_data
  {
    void               *user_data;      /* never changed by p6est */
    long                user_long;      /* never changed by p6est */
    int                 user_int;       /* never changed by p6est */
    p4est_topidx_t      which_tree;
    struct
    {
      p4est_topidx_t      which_tree;
      int                 owner_rank;
    }
    piggy1;
    struct
    {
      p4est_topidx_t      which_tree;
      p4est_topidx_t      from_tree;
    }
    piggy2;
    struct
    {
      p4est_topidx_t      which_tree;
      p4est_locidx_t      local_num;
    }
    piggy3;
  }
  p;
}
p2est_quadrant_t;

typedef struct p6est
{
  MPI_Comm            mpicomm;
  int                 mpisize, mpirank;

  size_t              data_size;        /* size of per-quadrant user_data */
  void               *user_pointer;     /* convenience pointer for users,
                                           will never be touched by p4est */
  p6est_connectivity_t *connectivity;   /* topology of sheet, not owned. */
  double              height[3];        /* offset from bottom to top,
                                           if uniform height profile */
  double             *top_to_vertices; /* if this is not NULL, it should be an
                                          array of the same size as \a
                                          connectivity->tree_to_vertex: points
                                          in the interior of the sheet
                                          trilinearly interpolate between them. */
  p4est_t            *columns;          /* 2D description of column layout
                                           built from \a connectivity */
  sc_array_t         *layers;           /* single array that stores
                                           p2est_quadrant_t layers within columns */
  sc_mempool_t       *user_data_pool;   /* memory allocator for user data
                                         * WARNING: This is NULL if data size
                                         *          equals zero.  */
  sc_mempool_t       *quadrant_pool;    /* memory allocator
                                           for temporary quadrants */
  p4est_gloidx_t     *global_first_quadrant;
}
p6est_t;

/** Callback function prototype to initialize the quadrant's user data.
 */
typedef void        (*p6est_init_t) (p6est_t * p6est,
                                     p4est_topidx_t which_tree,
                                     p4est_quadrant_t * column,
                                     p2est_quadrant_t * layer);

/** Callback function prototype to decide whether to horizontally refine a
 * layer.
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
 * \param [in] layers      Pointers to 2 vertical siblings.
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
p6est_t            *p6est_new (MPI_Comm mpicomm,
                               p6est_connectivity_t * connectivity,
                               size_t data_size,
                               p6est_init_t init_fn, void *user_pointer);

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
 * \param [in] p6est     The p6est to be worked on.
 * \param [in] btype     Balance type (face, corner or default, full).
 * \param [in] init_fn   Callback function to initialize the user_data
 *                       which is already allocated automatically.
 */
void                p6est_balance (p6est_t * p6est,
                                   p8est_connect_type_t btype,
                                   p6est_init_t init_fn);

/** Equally partition the forest.
 *
 * The forest will be partitioned between processors where they each
 * have an approximately equal number of quadrants.
 *
 * \param [in,out] p6est      The forest that will be partitioned.
 * \param [in]     weight_fn  A weighting function or NULL
 *                            for uniform partitioning.
 */
void                p6est_partition (p6est_t * p6est,
                                     p6est_weight_t weight_fn);

/** Compute the checksum for a forest.
 * Based on quadrant arrays only. It is independent of partition and mpisize.
 * \return  Returns the checksum on processor 0 only. 0 on other processors.
 */
unsigned            p6est_checksum (p6est_t * p6est);

/** Save the complete connectivity/p6est data to disk.  This is a collective
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
p6est_t            *p6est_load (const char *filename, MPI_Comm mpicomm,
                                size_t data_size, int load_data,
                                void *user_pointer,
                                p4est_connectivity_t ** connectivity);

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


SC_EXTERN_C_END;

#endif /* !P6EST_H */
