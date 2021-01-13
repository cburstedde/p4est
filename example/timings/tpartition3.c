#include <p4est_to_p8est.h>
#ifndef P4_TO_P8
#error "This program is currently intended for 3D only"
#else
#include <p8est_extended.h>
#include <p8est_vtk.h>
#endif
#include <sc_flops.h>
#include <sc_options.h>
#include <sc_statistics.h>

static int          refine_level = 0;
static double       refinement_fraction = 0.25;

enum
{
  TBSEARCH_PARTITION,
  TBSEARCH_NUM_STATS
};

static int
refine_fraction (p4est_t * p4est, p4est_topidx_t which_tree,
                p4est_quadrant_t * q)
{
  /* The formula in the line below implies
   * quadrant_fraction = 7 * refinement_fraction + 1.
   * In particular we need for doubling the quadrants for a level increment
   * refinment_factor = 1 / 7.
   */
  /* p4est->global_num_quadrants is the old number of quadrants */
  const p4est_locidx_t quad_count_refinement_threshold =
    (p4est_locidx_t) (refinement_fraction * p4est->global_num_quadrants);
  p4est_locidx_t      local_refined_quads;
  p4est_gloidx_t     *refinement_counter =
    (p4est_gloidx_t *) p4est->user_pointer;

#if 0
  if ((int) q->level >= refine_level) {
    ++(*refinement_counter);
    return 0;
  }
#endif

  if (p4est->global_first_quadrant[p4est->mpirank + 1] - 1 <
      quad_count_refinement_threshold) {
    /* This process owns quadrants that we want to refine. */
    /* On this process is no refinement conting needed. */
    return 1;
  }

  if (p4est->global_first_quadrant[p4est->mpirank] <=
      quad_count_refinement_threshold
      && p4est->global_first_quadrant[p4est->mpirank + 1] - 1 >=
      quad_count_refinement_threshold) {
    /* This process owns the threshold quadrant. */
    local_refined_quads =
      quad_count_refinement_threshold -
      p4est->global_first_quadrant[p4est->mpirank];
    if (*refinement_counter < local_refined_quads) {
      /* refine only a defined fraction of the quadrants */
      ++(*refinement_counter);
      return 1;
    }
  }

  /* If this return is reached we already know that we reached the threshold.
   * That is why we do not have to increment refinement_counter.
   */
  return 0;
}

int
main (int argc, char **argv)
{
  int                 i;
  int                 mpiret;
  sc_MPI_Comm         mpicomm;
  int                 rank, size;
  int                 start_level;
  int                 write_vtk;
  p4est_gloidx_t      refinement_counter;
  p8est_connectivity_t *connectivity;
  p8est_t            *p8est;
  sc_options_t       *opt;
  sc_flopinfo_t       fi, snapshot;
  sc_statinfo_t       stats[TBSEARCH_NUM_STATS];

  /* initialize MPI and p4est internals */
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  mpicomm = sc_MPI_COMM_WORLD;
  mpiret = sc_MPI_Comm_size (mpicomm, &size);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (mpicomm, &rank);
  SC_CHECK_MPI (mpiret);

  sc_init (mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
#ifndef P4EST_ENABLE_DEBUG
  sc_set_log_defaults (NULL, NULL, SC_LP_STATISTICS);
#endif
  p4est_init (NULL, SC_LP_DEFAULT);

  opt = sc_options_new (argv[0]);
  sc_options_add_int (opt, 'l', "start level", &start_level, 4,
                      "Level of the uniform starting grid");
  sc_options_add_int (opt, 'L', "refine level", &refine_level, 10,
                      "Highest level");
  sc_options_add_double (opt, 'f', "refinement fraction",
                         &refinement_fraction, 0.25,
                         "fraction of quadrants that will be refined \
                        in each refinement iteration");
  sc_options_add_switch (opt, 'V', "write-vtk", &write_vtk,
                         "write vtk output");
  sc_options_parse (p4est_package_id, SC_LP_DEFAULT, opt, argc, argv);

  connectivity = p8est_connectivity_new_unitcube ();

  p8est =
    p4est_new_ext (mpicomm, connectivity, 0, start_level, 1, 0, NULL,
                   &refinement_counter);

  /* start overall timing */
  mpiret = sc_MPI_Barrier (p8est->mpicomm);
  SC_CHECK_MPI (mpiret);
  sc_flops_start (&fi);

#if 1
  for (i = start_level; i < refine_level; ++i) {
    refinement_counter = 0;
    p4est_refine (p8est, 0, refine_fraction, NULL);
    if (i == refine_level - 1) {
      /* We only partition for the last iteration to be more demanding 
       * for p4est_partition.
       */
      sc_flops_shot (&fi, &snapshot);
      p4est_partition (p8est, 0, NULL);
      sc_flops_shot (&fi, &snapshot);
      sc_stats_set1 (&stats[TBSEARCH_PARTITION], snapshot.iwtime,
                     "Partition");
    }
  }
#endif
#if 0
  /* Recursive refinment version that corresponds to the preprocessor if
   * in refine_fraction.
   */
  refinement_counter = 0;
  p4est_refine (p8est, 1, refine_fraction, NULL);
#endif

  sc_stats_compute (p8est->mpicomm, TBSEARCH_NUM_STATS, stats);
  sc_stats_print (p4est_package_id, SC_LP_ESSENTIAL,
                  TBSEARCH_NUM_STATS, stats, 1, 1);

  if (write_vtk) {
    p4est_vtk_write_file (p8est, NULL, "tbsearch");
  }

  p4est_destroy (p8est);
  p4est_connectivity_destroy (connectivity);

  sc_options_destroy (opt);
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
