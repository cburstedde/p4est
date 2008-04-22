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

#ifndef __P4EST_OPTIONS_H__
#define __P4EST_OPTIONS_H__

#include <p4est_memory.h>

typedef enum
{
  P4EST_OPTION_SWITCH,
  P4EST_OPTION_INT,
  P4EST_OPTION_DOUBLE,
  P4EST_OPTION_STRING,
}
p4est_option_type_t;

typedef struct
{
  p4est_option_type_t opt_type;
  int                 opt_char;
  const char         *opt_name;
  void               *opt_var;
  const char         *help_string;
}
p4est_option_item_t;

typedef struct
{
  char                program_path[BUFSIZ];
  const char         *program_name;
  p4est_array_t      *option_items;
}
p4est_options_t;

/**
 * Create an empty options structure.
 * \param [in] program_path   Name or path name of the program.
 */
p4est_options_t    *p4est_options_new (const char *program_path);

/**
 * Destroy the options structure.
 */
void                p4est_options_destroy (p4est_options_t * opt);

/**
 * Add a switch option. This option is used without option arguments.
 * Every use increments the variable by one. Its initial value is 0.
 * Either opt_char or opt_long must be valid, that is, not '\0'/NULL.
 * \param [in] opt_char      Short option character, may be '\0'.
 * \param [in] opt_long      Option name without initial dashes, may be NULL.
 * \param [in] variable      Address of the variable to store the option value.
 * \param [in] help_string   Help string for usage message, may be NULL.
 */
void                p4est_options_add_switch (p4est_options_t * opt,
                                              int opt_char,
                                              const char *opt_long,
                                              int *variable,
                                              const char *help_string);

/**
 * Add an option that takes an integer argument.
 * \param [in] init_value   The initial value of the variable.
 */
void                p4est_options_add_int (p4est_options_t * opt,
                                           int opt_char,
                                           const char *opt_long,
                                           int *variable, int init_value,
                                           const char *help_string);

/**
 * Add an option that takes a double argument.
 */
void                p4est_options_add_double (p4est_options_t * opt,
                                              int opt_char,
                                              const char *opt_long,
                                              double *variable,
                                              double init_value,
                                              const char *help_string);

/**
 * Add a string option.
 * \param [in] init_value  The default value of the option may be NULL.
 */
void                p4est_options_add_string (p4est_options_t * opt,
                                              int opt_char,
                                              const char *opt_long,
                                              const char **variable,
                                              const char *init_value,
                                              const char *help_string);

/**
 * Print a usage message.
 * \param [in] include_args   Include the string <ARGS> in the message.
 * \param [in] nout           The stream to print to, may be NULL.
 */
void                p4est_options_print_help (p4est_options_t * opt,
                                              int include_args, FILE * nout);

/**
 * Print a summary of all option values.
 * Produces the title "Summary:", the subtitle "OPTIONS:"
 *          and a line for every option.
 */
void                p4est_options_print_summary (p4est_options_t * opt,
                                                 FILE * nout);

/**
 * Print a summary of all non-option arguments.
 * Produces the subtitle "ARGUMENTS:" and a line for every argument.
 * \param [in] first_arg   Position of first non-option argument.
 */
void                p4est_options_print_arguments (p4est_options_t * opt,
                                                   int first_arg,
                                                   int argc, char **argv,
                                                   FILE * nout);

/**
 * Parse command line options.
 * \param [in,out] argv   The options may be moved before the non-options.
 * \param [in] nerr       Stream for error output, may be NULL.
 * \return                Returns -1 on an invalid option, otherwise
 *                                the position of the first non-option argument.
 */
int                 p4est_options_parse (p4est_options_t * opt,
                                         int argc, char **argv, FILE * nerr);

#endif /* !__P4EST_OPTIONS_H__ */
