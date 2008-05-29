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

#ifndef P4EST_FILE_H
#define P4EST_FILE_H

#include <p4est.h>

/** Read in forest connectivity info from a file
 *
 * This function fills an allocated connectivity structure
 * with the data from a p4t file.
 *
 * The p4t file has the following structure.
 *
 * \verbatim
 * [Forest Info]
 * ver = 0.0.1  # Version of the forest file
 * Nk  = 3      # Number of elements
 * Nv  = 7      # Number of mesh vertices
 * Nve = 12     # Number of trees in the vertex to element list
 * Net = 0      # Number of element tags
 * Nft = 0      # Number of face tags
 * Ncf = 0      # Number of curved faces
 * Nct = 0      # Number of curved types
 * [Coordinates of Element Vertices]
 * 1 -1.00000000000e+00 -1.00000000000e+00
 * 2  0.00000000000e+00 -1.00000000000e+00
 * 3  0.00000000000e+00  0.00000000000e+00
 * 4  1.00000000000e+00  0.00000000000e+00
 * 5  1.00000000000e+00  1.00000000000e+00
 * 6  0.00000000000e+00  1.00000000000e+00
 * 7 -1.00000000000e+00  0.00000000000e+00
 * [Element to Vertex]
 * 1     1   2   4   3
 * 2     1   3   6   7
 * 3     3   4   5   6
 * [Element to Element]
 * 1     1   1   3   2
 * 2     1   3   2   2
 * 3     1   3   3   2
 * [Element to Face]
 * 1     1   2   1   1
 * 2     4   4   3   4
 * 3     3   2   3   2
 * [Vertex to Element]
 * 1     2   1   2
 * 2     1   1
 * 3     3   1   3   2
 * 4     2   1   3
 * 5     1   3
 * 6     2   2   3
 * 7     1   2
 * [Vertex to Vertex]
 * 1     2   1   1
 * 2     1   2
 * 3     3   3   3   3
 * 4     2   4   4
 * 5     1   5
 * 6     2   6   6
 * 7     1   7
 * [Element Tags]
 * [Face Tags]
 * [Curved Faces]
 * [Curved Types]
 * \endverbatim
 *
 * Note that the mesh file is 1-based which will get converted into
 * the 0-based connectivity structure.
 *
 * Note that for elements that do not have neighbors
 * the neighboring element is itself.  This is also true
 * for the neighboring face.
 *
 * \param [in]  filename     name of the file to read from
 * \param [out] connectivity allocates a connectivity structure
 *                           that is returned filled with the info
 *                           from the file \a filename
 *
 * \return Upon sucess returns 0.
 *
 * \note The user of this function is responsible for calling
 *       \c p4est_connectivity_destroy when they are finished
 *       with \a connectivity.
 *
 */
int                 p4est_connectivity_read (const char *filename,
                                             p4est_connectivity_t **
                                             connectivity);
/** Prints a forest connectivity file
 *
 * This function prints a connectivity structure to a file stream.
 *
 * \param [in]  connectivity allocates a connectivity structure
 *                           that is returned filled with the info
 *                           from the file \a filename
 * \param [out] nout         stream that the connectivity file is printed to.
 *
 */
void                p4est_connectivity_print (p4est_connectivity_t *
                                              connectivity, FILE * nout);

#endif /* !P4EST_FILE_H */
