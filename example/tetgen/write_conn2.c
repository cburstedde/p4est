/*
  This file is part of p4est.
  p4est is a C library to manage a collection (a forest) of multiple
  connected adaptive quadtrees or octrees in parallel.

  Copyright (C) 2012 Carsten Burstedde

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

/* This program is NOT collective (i.e., usually NOT called with mpirun). */

#ifndef P4_TO_P8
#include <p4est_connectivity.h>
#else
#include <p8est_connectivity.h>
#endif

int
main (int argc, char **argv)
{
  p4est_connectivity_t *conn;
  int                 retval;
  char                buf[BUFSIZ];

  if (argc != 2) {
    char               *cp, *bn;

    cp = strdup (argv[0]);
    bn = basename (cp);
    fprintf (stderr, "Usage: %s <connectivity name>\n", bn);
    free (cp);
    exit (1);
  }

  conn = p4est_connectivity_new_byname (argv[1]);
  if (conn == NULL) {
    P4EST_PRODUCTIONF ("Failed to identify connectivity \"%s\"\n", argv[1]);
  }
  else {
    snprintf (buf, BUFSIZ, "%s_%s.p%dc",
              P4EST_STRING, argv[1], P4EST_CHILDREN);
    P4EST_INFOF ("Write connectivity file \"%s\"\n", buf);
    retval = p4est_connectivity_save (buf, conn);
    if (retval) {
      P4EST_PRODUCTION ("Error writing connectivity file\n");
    }
    p4est_connectivity_destroy (conn);
  }

  return 0;
}
