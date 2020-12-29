#ifndef P4_TO_P8
#include <p4est_extended.h>
#else
#include <p8est_extended.h>
#endif
#include <sc_flops.h>
#include <sc_statistics.h>
#include <sc_options.h>

enum
{
  TSUCCESSOR_NEW,
  TSUCCESSOR_NUM_STATS
};

int
main (int argc, char **argv)
{
  int                 mpiret;
  sc_MPI_Comm         mpicomm;
  int                 rank, size;
  int                 refine_level;
  p4est_connectivity_t *connectivity;
  p4est_t            *p4est;
  sc_options_t       *opt;
  sc_flopinfo_t       fi, snapshot;
  sc_statinfo_t       stats[TSUCCESSOR_NUM_STATS];

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

  sc_options_add_int (opt, 'l', "refine level", &refine_level, 5,
                      "Level of the uniform grid");
  sc_options_parse (p4est_package_id, SC_LP_DEFAULT, opt, argc, argv);

#ifdef P4_TO_P8
  connectivity = p8est_connectivity_new_unitcube ();
#else
  connectivity = p4est_connectivity_new_unitsquare ();
#endif

  /* start overall timing */
  mpiret = sc_MPI_Barrier (mpicomm);
  SC_CHECK_MPI (mpiret);
  sc_flops_start (&fi);

  sc_flops_shot (&fi, &snapshot);
  p4est =
    p4est_new_ext (mpicomm, connectivity, 0, refine_level, 1, 0, NULL, NULL);
  sc_flops_shot (&fi, &snapshot);
  sc_stats_set1 (&stats[TSUCCESSOR_NEW], snapshot.iwtime,
                 "New");

  sc_stats_compute (p4est->mpicomm, TSUCCESSOR_NUM_STATS, stats);
  sc_stats_print (p4est_package_id, SC_LP_ESSENTIAL,
                  TSUCCESSOR_NUM_STATS, stats, 1, 1);

  p4est_destroy (p4est);
  p4est_connectivity_destroy (connectivity);

  sc_options_destroy (opt);
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
