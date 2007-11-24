
#ifndef __P4EST_FILE_H__
#define __P4EST_FILE_H__

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

#endif /* !__P4EST_FILE_H__ */
