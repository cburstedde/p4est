
#include <p4est_file.h>
#include <p4est_base.h>

enum Section
{
  INFO,                         /* [Forest Info] */
  COORD,                        /* [Coordinates of Element Vertices] */
  ETOV,                         /* [Element to Vertex] */
  ETOE,                         /* [Element to Element] */
  ETOF,                         /* [Element to Face] */
  ET,                           /* [Element Tags] */
  FT,                           /* [Face Tags] */
  CF,                           /* [Curved Faces] */
  CT                            /* [Curved Types] */
};

int
p4est_connectivity_read (const char *filename,
                         p4est_connectivity_t ** connectivity)
{
  FILE               *file;
  int                 retval;
  char                line[BUFSIZ];
  char               *begin, *end;
  int                 i;

  file = fopen (filename, "rb");
  if (!file) {
    fprintf (stderr, "Failed to open p4est mesh file %s\n", filename);
    return 1;
  }

  i = 0;
  /* loop through the lines of the file */
  while (fgets (line, BUFSIZ, file)) {
    fprintf (stdout, "[%04d] %s", ++i, line);

    /* Look for a section */
  }

  retval = fclose (file);
  if (retval) {
    fprintf (stderr, "Failed to close p4est mesh file %s (%d:%d)\n", filename,
             retval, EOF);
    return 1;
  }

  return 0;
}

/* EOF p4est_file.h> */
