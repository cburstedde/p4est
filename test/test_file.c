/*
  This file is part of p4est.
  p4est is a C library to manage a parallel collection of quadtrees and/or
  octrees.

  Copyright (C) 2007 Carsten Burstedde, Lucas Wilcox.

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <p4est.h>
#include <p4est_base.h>
#include <p4est_file.h>

int
main (int argc, char **argv)
{
  int                 retval;
  int                 rank = 0;

#ifdef HAVE_MPI
  int                 use_mpi = 1;
  int                 mpiret;
#endif
  MPI_Comm            mpicomm;
  p4est_connectivity_t *connectivity;

  mpicomm = MPI_COMM_NULL;
#ifdef HAVE_MPI
  if (use_mpi) {
    mpiret = MPI_Init (&argc, &argv);
    P4EST_CHECK_MPI (mpiret);
    mpicomm = MPI_COMM_WORLD;
  }
#endif

  char                template[] = "p4est_meshXXXXXX";
  char                mesh[] = "		[Forest Info] # ]] [[ ]]\n"
    "ver = 0.0.1  # Version of the forest file\n"
    "Nk  = 3      # Number of elements\n"
    "Nv  = 7      # Number of mesh vertices\n"
    "Net = 0      # Number of element tags\n"
    "Nft = 0      # Number of face tags\n"
    "Ncf = 0      # Number of curved faces\n"
    "Nct = 0      # Number of curved types\n"
    "\n"
    "                          [Coordinates of Element Vertices]\n"
    "1 -1.00000000000e+00 -1.00000000000e+00  0.00000000000e+00\n"
    "2  0.00000000000e+00 -1.00000000000e+00  0.00000000000e+00\n"
    "3  0.00000000000e+00  0.00000000000e+00  0.00000000000e+00\n"
    "4  1.00000000000e+00  0.00000000000e+00  0.00000000000e+00\n"
    "5  1.00000000000e+00  1.00000000000e+00  0.00000000000e+00\n"
    "6  0.00000000000e+00  1.00000000000e+00  0.00000000000e+00\n"
    "7 -1.00000000000e+00  0.00000000000e+00  0.00000000000e+00\n"
    "   [Element to Vertex]\n"
    "1     1   2   4   3\n"
    "2     1   3   6   7\n"
    "3     3   4   5   6\n"
    "  [Element to Element]\n"
    "1     1   1   3   2\n"
    "2     1   3   2   2\n"
    "3     1   3   3   2\n"
    "\n"
    "[Element to Face]\n"
    "1     1   2   1   1\n"
    "2     4   4   3   4\n"
    "3     3   2   3   2\n"
    "\n"
    "[Element Tags]\n" "[Face Tags]\n" "[Curved Faces]\n" "[Curved Types]\n";

  int                 fd;
  FILE               *outfile;
#ifdef HAVE_MPI
  size_t              templatelength;
#endif
  const double        EPS = 2.22045e-16;

  int32_t             i;
  const int32_t       num_trees = 3;
  const int32_t       num_vertices = 7;
  const int32_t       tree_to_vertex[] = {
    0, 1, 3, 2, 0, 2, 5, 6, 2, 3, 4, 5
  };
  const int32_t       tree_to_tree[] = {
    0, 0, 2, 1, 0, 2, 1, 1, 0, 2, 2, 1
  };
  const int32_t       tree_to_face[] = {
    0, 1, 0, 0, 3, 3, 2, 3, 2, 1, 2, 1
  };
  const double        vertices[] = {
    -1.00000000000e+00, -1.00000000000e+00, 0.00000000000e+00,
    0.00000000000e+00, -1.00000000000e+00, 0.00000000000e+00,
    0.00000000000e+00, 0.00000000000e+00, 0.00000000000e+00,
    1.00000000000e+00, 0.00000000000e+00, 0.00000000000e+00,
    1.00000000000e+00, 1.00000000000e+00, 0.00000000000e+00,
    0.00000000000e+00, 1.00000000000e+00, 0.00000000000e+00,
    -1.00000000000e+00, 0.00000000000e+00, 0.00000000000e+00
  };

#ifdef HAVE_MPI
  mpiret = MPI_Comm_rank (mpicomm, &rank);
  P4EST_CHECK_MPI (mpiret);
#endif

  if (rank == 0) {
    /* Make a temporary file to hold the mesh */
    fd = mkstemp (template);
    P4EST_CHECK_ABORT (fd != -1, "Unable to create temp mesh file.");

    /* Promote the file descriptor to a FILE stream */
    outfile = fdopen (fd, "wb");
    P4EST_CHECK_ABORT (outfile != NULL, "Unable to fdopen temp mesh file.");

    /* Write out to the mesh to the temporary file */
    retval = fputs (mesh, outfile);
    P4EST_CHECK_ABORT (retval != EOF, "Unable to fputs temp mesh file.");

    /* Close the temporary file */
    retval = fclose (outfile);
    P4EST_CHECK_ABORT (!retval, "Unable to fclose the temp mesh file.");
  }

#ifdef HAVE_MPI
  templatelength = strlen (template);
  mpiret = MPI_Bcast (template, templatelength + 1, MPI_CHAR, 0, mpicomm);
  P4EST_CHECK_MPI (mpiret);
#endif

  /* Read in the mesh into connectivity information */
  retval = p4est_connectivity_read (template, &connectivity);
  P4EST_CHECK_ABORT (!retval, "Unable to read the mesh file.");

  /* Check what was read in */
  P4EST_CHECK_ABORT (connectivity->num_trees == num_trees, "num_trees");
  P4EST_CHECK_ABORT (connectivity->num_vertices == num_vertices,
                     "num_vertices");
  for (i = 0; i < num_trees * 4; ++i)
    P4EST_CHECK_ABORT (connectivity->tree_to_vertex[i] == tree_to_vertex[i],
                       "tree_to_vertices");
  for (i = 0; i < num_trees * 4; ++i)
    P4EST_CHECK_ABORT (connectivity->tree_to_tree[i] == tree_to_tree[i],
                       "tree_to_tree");
  for (i = 0; i < num_trees * 4; ++i)
    P4EST_CHECK_ABORT (connectivity->tree_to_face[i] == tree_to_face[i],
                       "tree_to_face");
  for (i = 0; i < num_vertices * 3; ++i)
    P4EST_CHECK_ABORT (connectivity->vertices[i] - vertices[i] < EPS,
                       "vertices");

  /* destroy the p4est and its connectivity structure */
  p4est_connectivity_destroy (connectivity);

  if (rank == 0) {
    /* remove the temporary file */
    retval = remove (template);
    P4EST_CHECK_ABORT (!retval, "Unable to remove the temp mesh file.");
  }

  /* clean up and exit */
  p4est_memory_check ();

#ifdef HAVE_MPI
  if (use_mpi) {
    mpiret = MPI_Finalize ();
    P4EST_CHECK_MPI (mpiret);
  }
#endif

  return 0;
}

/* EOF test_file.c */
