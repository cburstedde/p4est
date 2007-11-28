
#include <p4est.h>
#include <p4est_algorithms.h>
#include <p4est_base.h>

static const int64_t initial_quadrants_per_processor = 15;
static const int number_toread_quadrants = 32;

p4est_t            *
p4est_new (MPI_Comm mpicomm, FILE * nout, p4est_connectivity_t * connectivity,
           int data_size, p4est_init_t init_fn)
{
#ifdef HAVE_MPI
  int                 mpiret;
#endif
  int                 i, must_remove_last_quadrant;
  int8_t              level;
  int32_t             j, num_trees;
  int64_t             tree_num_quadrants, global_num_quadrants;
  int64_t             first_tree, first_quadrant, first_tree_quadrant;
  int64_t             last_tree, last_quadrant, last_tree_quadrant;
  int64_t             quadrant_index;
  p4est_t            *p4est;
  p4est_tree_t       *tree;
  p4est_quadrant_t   *quad;

  p4est = P4EST_ALLOC_ZERO (p4est_t, 1);
  P4EST_CHECK_ALLOC (p4est);

  /* assign some data members */
  p4est->nout = nout;
  p4est->data_size = data_size;
  p4est->connectivity = connectivity;
  num_trees = connectivity->num_trees;

  p4est->mpicomm = mpicomm;
  p4est->mpisize = 1;
  p4est->mpirank = 0;
#ifdef HAVE_MPI
  if (p4est->mpicomm != MPI_COMM_NULL) {
    mpiret = MPI_Comm_size (p4est->mpicomm, &p4est->mpisize);
    P4EST_CHECK_MPI (mpiret);

    mpiret = MPI_Comm_rank (p4est->mpicomm, &p4est->mpirank);
    P4EST_CHECK_MPI (mpiret);
  }
#endif

  /* allocate memory pools */
  if (p4est->data_size > 0) {
    p4est->user_data_pool = p4est_mempool_new (p4est->data_size);
  }
  else {
    p4est->user_data_pool = NULL;
  }
  p4est->quadrant_pool = p4est_mempool_new (sizeof (p4est_quadrant_t));

  /* determine uniform level of initial tree */
  tree_num_quadrants = 1;
  for (level = 0; level < 16; ++level) {
    if (tree_num_quadrants >=
        (p4est->mpisize * initial_quadrants_per_processor) / num_trees) {
      break;
    }
    tree_num_quadrants *= 4;
  }
  P4EST_ASSERT (level < 16 && tree_num_quadrants <= INT32_MAX);

  /* compute index of first tree for this processor */
  global_num_quadrants = tree_num_quadrants * num_trees;
  first_quadrant = (global_num_quadrants * p4est->mpirank) / p4est->mpisize;
  first_tree = first_quadrant / tree_num_quadrants;
  first_tree_quadrant = first_quadrant - first_tree * tree_num_quadrants;
  last_quadrant =
    (global_num_quadrants * (p4est->mpirank + 1)) / p4est->mpisize - 1;
  last_tree = last_quadrant / tree_num_quadrants;
  last_tree_quadrant = last_quadrant - last_tree * tree_num_quadrants;

  /* check ranges of various integers to be 32bit compatible */
  P4EST_ASSERT (first_tree <= last_tree && last_tree < num_trees);
  P4EST_ASSERT (0 <= first_tree_quadrant && 0 <= last_tree_quadrant);
  P4EST_ASSERT (first_tree_quadrant < tree_num_quadrants);
  P4EST_ASSERT (last_tree_quadrant < tree_num_quadrants);
  if (first_tree == last_tree) {
    P4EST_ASSERT (first_tree_quadrant < last_tree_quadrant);
  }

  /* print some diagnostics */
  if (p4est->nout != NULL) {
    if (p4est->mpirank == 0) {
      fprintf (p4est->nout, "New forest: %d trees %d processors\n",
               num_trees, p4est->mpisize);
      fprintf (p4est->nout,
               "   initial level %d global quadrants %lld per tree %lld\n",
               level, (long long) global_num_quadrants,
               (long long) tree_num_quadrants);
    }
    fprintf (p4est->nout, "[%d] first tree %lld first quadrant %lld"
             " global quadrant %lld\n", p4est->mpirank,
             (long long) first_tree, (long long) first_tree_quadrant,
             (long long) first_quadrant);
    fprintf (p4est->nout, "[%d] last tree %lld last quadrant %lld"
             " global quadrant %lld\n", p4est->mpirank,
             (long long) last_tree, (long long) last_tree_quadrant,
             (long long) last_quadrant);
  }

  /* allocate trees and quadrants */
  p4est->trees = p4est_array_new (sizeof (p4est_tree_t));
  p4est_array_resize (p4est->trees, num_trees);
  for (j = 0; j < num_trees; ++j) {
    tree = p4est_array_index (p4est->trees, j);
    tree->quadrants = p4est_array_new (sizeof (p4est_quadrant_t));
    for (i = 0; i <= P4EST_MAXLEVEL; ++i) {
      tree->quadrants_per_level[i] = 0;
    }
    tree->maxlevel = 0;
  }

  /* for every locally non-empty tree fill first and last quadrant */
  for (j = first_tree; j <= last_tree; ++j) {
    tree = p4est_array_index (p4est->trees, j);
    must_remove_last_quadrant = 0;

    /* set morton id of first quadrant and initialize user data */
    p4est_array_resize (tree->quadrants, 1);
    quad = p4est_array_index (tree->quadrants, 0);
    if (j == first_tree) {
      p4est_quadrant_set_morton (quad, level, first_tree_quadrant);
    }
    else {
      p4est_quadrant_set_morton (quad, level, 0);
    }
    if (p4est->nout != NULL) {
      fprintf (p4est->nout, "[%d] tree %d first morton 0x%x 0x%x\n",
               p4est->mpirank, j, quad->x, quad->y);
    }
    p4est_quadrant_init_data (p4est, j, quad, init_fn);

    /* set morton id of last quadrant */
    if (tree_num_quadrants == 1 ||
        (j == first_tree && first_tree_quadrant == tree_num_quadrants - 1)) {
      /* nothing to do for this tree, we will have only one quadrant */
    }
    else {
      p4est_array_resize (tree->quadrants, 2);
      quad = p4est_array_index (tree->quadrants, 1);
      if (j == last_tree) {
        if (last_tree_quadrant == tree_num_quadrants - 1) {
          quadrant_index = last_tree_quadrant;
        }
        else {
          quadrant_index = last_tree_quadrant + 1;
          must_remove_last_quadrant = 1;
        }
        p4est_quadrant_set_morton (quad, level, quadrant_index);
      }
      else {
        p4est_quadrant_set_morton (quad, level, tree_num_quadrants - 1);
      }
      if (p4est->nout != NULL) {
        fprintf (p4est->nout, "[%d] tree %d last morton 0x%x 0x%x\n",
                 p4est->mpirank, j, quad->x, quad->y);
      }
      p4est_quadrant_init_data (p4est, j, quad, init_fn);

      /* now run algorithm CompleteRegion (&tree->quadrants) here */
    }
  }

  p4est->first_local_tree = first_tree;
  p4est->last_local_tree = last_tree;
  p4est->local_num_trees = last_tree - first_tree + 1;

  return p4est;
}

void
p4est_destroy (p4est_t * p4est)
{
  int32_t             j;
  p4est_tree_t       *tree;

  for (j = 0; j < p4est->connectivity->num_trees; ++j) {
    tree = p4est_array_index (p4est->trees, j);
    p4est_array_destroy (tree->quadrants);
  }
  p4est_array_destroy (p4est->trees);

  if (p4est->user_data_pool != NULL) {
    p4est_mempool_destroy (p4est->user_data_pool);
  }
  p4est_mempool_destroy (p4est->quadrant_pool);

  P4EST_FREE (p4est);
}

void
p4est_refine (p4est_t * p4est,
              p4est_refine_t refine_fn, p4est_init_t init_fn)
{
  int                 quadrant_pool_size;
  int                 dorefine;
  int32_t             j, movecount;
  int32_t             current, restpos, incount;
  p4est_list_t       *list;
  p4est_tree_t       *tree;
  p4est_quadrant_t   *q, *qalloc, *qpop;
  p4est_quadrant_t   *c0, *c1, *c2, *c3;
  int                *key;

  /*
    q points to a quadrant that is an array member
    qalloc is a quadrant that has been allocated through quadrant_pool
    qpop is q quadrant that has been allocated through quadrant_pool
    never mix these two types of quadrant pointers
  */
  key = &dorefine;      /* use this to create a unique user_data pointer */
  list = p4est_list_new (NULL);

  /* loop over all local trees */
  for (j = p4est->first_local_tree; j <= p4est->last_local_tree; ++j) {
    tree = p4est_array_index (p4est->trees, j);
    quadrant_pool_size = p4est->quadrant_pool->elem_count;

    if (p4est->nout != NULL) {
      fprintf (p4est->nout, "[%d] Into refine tree %d with %d\n",
               p4est->mpirank, j, tree->quadrants->elem_count);
    }

    /* run through the array to find first quadrant to be refined */
    q = NULL;
    dorefine = 0;
    incount = tree->quadrants->elem_count;
    for (current = 0; current < incount; ++current) {
      q = p4est_array_index (tree->quadrants, current);
      dorefine = refine_fn (p4est, j, q);
      if (dorefine) {
        break;
      }
    }
    if (!dorefine) {
      continue;
    }

    /* now we have a quadrant to refine, prepend it to the list */
    qalloc = p4est_mempool_alloc (p4est->quadrant_pool);
    *qalloc = *q;                       /* never prepend array members */
    p4est_list_prepend (list, qalloc);  /* only newly allocated quadrants */

    /*
      current points to the next array member to write
      restpos points to the next array member to read
    */
    restpos = current + 1;

    /* run through the list and refine recursively */
    while (list->elem_count > 0) {
      qpop = p4est_list_pop (list);
      if (dorefine || refine_fn (p4est, j, qpop)) {
        dorefine = 0;   /* a marker so that refine_fn is never called twice */
        p4est_array_resize (tree->quadrants, tree->quadrants->elem_count + 3);

        /* compute children and prepend them to the list */
        if (qpop->user_data != key) {
          p4est_quadrant_free_data (p4est, qpop);
        }
        c0 = qpop;
        c1 = p4est_mempool_alloc (p4est->quadrant_pool);
        c2 = p4est_mempool_alloc (p4est->quadrant_pool);
        c3 = p4est_mempool_alloc (p4est->quadrant_pool);
        p4est_quadrant_children (qpop, c0, c1, c2, c3);
        c0->user_data = key;
        c1->user_data = key;
        c2->user_data = key;
        c3->user_data = key;
        p4est_list_prepend (list, c3);
        p4est_list_prepend (list, c2);
        p4est_list_prepend (list, c1);
        p4est_list_prepend (list, c0);
      }
      else {
        /* need to make room in the array to store this new quadrant */
        if (restpos < incount && current == restpos) {
          movecount = P4EST_MIN (incount - restpos, number_toread_quadrants);
          while (movecount > 0) {
            q = p4est_array_index (tree->quadrants, restpos);
            qalloc = p4est_mempool_alloc (p4est->quadrant_pool);
            *qalloc = *q;       /* never append array members */
            p4est_list_append (list, qalloc);   /* only newly allocated quadrants */
            --movecount;
            ++restpos;
          }
        }

        /* store new quadrant and update counters */
        if (qpop->user_data == key) {
          p4est_quadrant_init_data (p4est, j, qpop, init_fn);
        }
        q = p4est_array_index (tree->quadrants, current);
        *q = *qpop;
        ++current;
        p4est_mempool_free (p4est->quadrant_pool, qpop);
      }
    }

    P4EST_ASSERT (restpos == incount);
    P4EST_ASSERT (current == tree->quadrants->elem_count);
    P4EST_ASSERT (list->first == NULL && list->last == NULL);
    P4EST_ASSERT (quadrant_pool_size == p4est->quadrant_pool->elem_count);
    P4EST_ASSERT (p4est_tree_is_sorted (tree));

    if (p4est->nout != NULL) {
      fprintf (p4est->nout, "[%d] Done refine tree %d now %d\n",
               p4est->mpirank, j, tree->quadrants->elem_count);
    }
  }

  p4est_list_destroy (list);
}

/* EOF p4est.c */
