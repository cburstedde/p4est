
#include <p4est_base.h>

#ifdef P4EST_BACKTRACE
#ifdef HAVE_EXECINFO_H
#include <execinfo.h>
#endif
#endif

void
p4est_abort (void)
{
#ifdef P4EST_BACKTRACE
  int                 i;
  size_t              bt_size;
  void               *bt_buffer[64];
  char              **bt_strings;
  const char         *str;

  bt_size = backtrace (bt_buffer, 64);
  bt_strings = backtrace_symbols (bt_buffer, bt_size);

  fprintf (stderr, "Abort: Obtained %ld stack frames\n", (long int) bt_size);

#ifdef P4EST_ADDRTOLINE
  /* implement pipe connection to addr2line */
#endif

  for (i = 0; i < bt_size; i++) {
    str = strrchr (bt_strings[i], '/');
    if (str != NULL) {
      ++str;
    }
    else {
      str = bt_strings[i];
    }
    /* fprintf (stderr, "   %p %s\n", bt_buffer[i], str); */
    fprintf (stderr, "   %s\n", str);
  }
  free (bt_strings);
#endif

  fflush (stdout);
  fflush (stderr);

  abort ();
}

/* EOF p4est_base.c */
