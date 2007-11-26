
#include <p4est.h>
#include <p4est_file.h>
#include <p4est_base.h>

int
main (void)
{
  int                 retval;

  p4est_connectivity_t *connectivity;

  retval = p4est_connectivity_read ("mesh.p4t", &connectivity);
  P4EST_CHECK_ABORT (!retval, "Unable to read the mesh file.");

  return 0;
}

/* EOF read_forest.c */
