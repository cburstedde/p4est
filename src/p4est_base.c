
#include <p4est_base.h>

#ifdef P4EST_BACKTRACE
#ifdef HAVE_EXECINFO_H
#include <execinfo.h>
#endif
#endif

static int          malloc_count = 0;
static int          free_count = 0;

void               *
p4est_malloc (size_t size)
{
  void             *ret;

  ret = malloc (size);

  if (size > 0) {
    ++malloc_count;
  }
  else {
    malloc_count += ((ret == NULL) ? 0 : 1);
  }
  
  return ret;
}

void               *
p4est_calloc (size_t nmemb, size_t size)
{
  void             *ret;

  ret = calloc (nmemb, size);

  if (nmemb * size > 0) {
    ++malloc_count;
  }
  else {
    malloc_count += ((ret == NULL) ? 0 : 1);
  }
  
  return ret;  
}

void               *
p4est_realloc (void * ptr, size_t size)
{
  void             *ret;

  ret = realloc (ptr, size);

  if (ptr == NULL) {
    if (size > 0) {
      ++malloc_count;
    }
    else {
      malloc_count += ((ret == NULL) ? 0 : 1);
    }
  }
  else if (size == 0) {
    free_count += ((ret == NULL) ? 1 : 0);
  }

  return ret;
}

void
p4est_free (void * ptr)
{
  if (ptr != NULL) {
    ++free_count;
    free (ptr);
  }
}

void
p4est_memory_check (void) 
{
  P4EST_CHECK_ABORT (malloc_count == free_count, "Memory balance");
}

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
