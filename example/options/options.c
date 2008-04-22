/*
  This file is part of p4est.
  p4est is a C library to manage a parallel collection of quadtrees and/or
  octrees.

  Copyright (C) 2008 Carsten Burstedde, Lucas Wilcox.

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

#include <p4est_base.h>
#include <p4est_options.h>

int
main (int argc, char **argv)
{
  int                 first_arg;
  int                 w;
  int                 i1, i2;
  double              d;
  const char         *s1, *s2;
  p4est_options_t    *opt;

  opt = p4est_options_new (argv[0]);
  p4est_options_add_switch (opt, 'w', "switch", &w, "Switch");
  p4est_options_add_int (opt, 'i', "integer1", &i1, 0, "Integer 1");
  p4est_options_add_double (opt, 'd', "double", &d, 0., "Double");
  p4est_options_add_string (opt, 's', "string", &s1, NULL, "String 1");
  p4est_options_add_string (opt, 't', NULL, &s2, NULL, "String 2");
  p4est_options_add_int (opt, '\0', "integer2", &i2, 7, "Integer 2");

  first_arg = p4est_options_parse (opt, argc, argv, stdout);
  if (first_arg < 0) {
    p4est_options_print_help (opt, 1, stdout);
  }
  else {
    p4est_options_print_summary (opt, stdout);
    p4est_options_print_arguments (opt, first_arg, argc, argv, stdout);
  }

  p4est_options_destroy (opt);

  p4est_memory_check ();

  return 0;
}

/* EOF options.c */
