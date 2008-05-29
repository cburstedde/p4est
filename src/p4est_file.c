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

#include <p4est_file.h>

enum Section
{
  NONE,
  INFO,                         /* [Forest Info] */
  COORD,                        /* [Coordinates of Element Vertices] */
  ETOV,                         /* [Element to Vertex] */
  ETOE,                         /* [Element to Element] */
  ETOF,                         /* [Element to Face] */
  VTOE,                         /* [Vertex to Element] */
  ET,                           /* [Element Tags] */
  FT,                           /* [Face Tags] */
  CF,                           /* [Curved Faces] */
  CT                            /* [Curved Types] */
};

static void
p4est_trim_comments (char *line)
{
  /* Truncate comments by setting them to '\0' */

  for (; *line != '\0'; ++line) {
    if (*line == '#') {
      *line = '\0';
      break;
    }
  }
}

static void
p4est_trim_ending_whitespace (char *line)
{
  int                 n;

  /* Trim whitespace from the end of the line */
  for (n = strlen (line) - 1; n >= 0; n--)
    if (!isspace (line[n]))
      break;
  line[n + 1] = '\0';
}

static char        *
p4est_trim_beginning_whitespace (char *line)
{
  /* Trim whitespace from the beginning of the line */
  for (; isspace (*line); ++line) {
  }

  return line;
}

int
p4est_connectivity_read (const char *filename,
                         p4est_connectivity_t ** connectivity)
{
  FILE               *file;
  size_t              length;
  int                 retval;
  char                buf[BUFSIZ];
  char               *line;
  enum Section        section = NONE;
  int                 section_lines_read = 0;
  int                 set_num_trees = 0;
  int                 set_num_vertices = 0;
  int                 set_num_vtt = 0;
  char               *key = NULL;
  char               *value = NULL;
  int32_t             vertex, corner_trees, ctree;
  int32_t             num_trees = 0;
  int32_t             num_vertices = 0;
  int32_t             num_vtt = 0;
  int32_t             k, k0, k1, k2, k3, f0, f1, f2, f3, v0, v1, v2, v3;
  int32_t             i, Nnn;
  int32_t            *tree_to_vertex = NULL;
  int32_t            *tree_to_tree = NULL;
  int32_t            *vtt_offset = NULL;
  int32_t            *vertex_to_tree = NULL;
  int32_t            *vertex_to_vertex = NULL;
  int8_t             *tree_to_face = NULL;
  double             *vertices = NULL;
  double              vx, vy, vz;

  *connectivity = NULL;

  file = fopen (filename, "rb");
  if (!file) {
    fprintf (stderr, "Failed to open p4est mesh file %s\n", filename);
    return 1;
  }

  /* loop through the lines of the file */
  while (fgets (buf, BUFSIZ, file)) {
    line = buf;

    p4est_trim_comments (line);

    p4est_trim_ending_whitespace (line);

    line = p4est_trim_beginning_whitespace (line);

    if (*line == '\0') {
      /* skip empty lines */
    }
    else if (*line == '[') {
      /* call any checks before leaving a section */
      switch (section) {
      case ETOV:
        SC_CHECK_ABORT (section_lines_read == num_trees,
                        "Not enough entries in [Element to Vertex]");
        break;
      case ETOE:
        SC_CHECK_ABORT (section_lines_read == num_trees,
                        "Not enough entries in [Element to Element]");
        break;
      case ETOF:
        SC_CHECK_ABORT (section_lines_read == num_trees,
                        "Not enough entries in [Element to Face]");
        break;
      case VTOE:
        SC_CHECK_ABORT (section_lines_read == num_vertices,
                        "Not enough entries in [Vertex to Element]");
        break;
      default:
        ;
      }

      /* set section */
      length = strlen (line);

      /* Sections must end with ']' */
      SC_CHECK_ABORT (line[length - 1] == ']', "Sections must end with ']'");

      line[length - 1] = '\0';
      ++line;

      if (strcmp (line, "Forest Info") == 0) {
        section = INFO;
      }
      else if (strcmp (line, "Coordinates of Element Vertices") == 0) {
        section = COORD;
      }
      else if (strcmp (line, "Element to Vertex") == 0) {
        section = ETOV;
      }
      else if (strcmp (line, "Element to Element") == 0) {
        section = ETOE;
      }
      else if (strcmp (line, "Element to Face") == 0) {
        section = ETOF;
      }
      else if (strcmp (line, "Vertex to Element") == 0) {
        section = VTOE;
      }
      else if (strcmp (line, "Element Tags") == 0) {
        section = ET;
      }
      else if (strcmp (line, "Face Tags") == 0) {
        section = FT;
      }
      else if (strcmp (line, "Curved Faces") == 0) {
        section = CF;
      }
      else if (strcmp (line, "Curved Types") == 0) {
        section = CT;
      }
      else {
        SC_CHECK_ABORT (0, "Unknown section in mesh file");
      }

      SC_CHECK_ABORT (section == INFO || *connectivity != NULL,
                      "The [Forest Info] section must come first"
                      " and set Nk, Nv, and Nve.");

      section_lines_read = 0;
    }
    else {
      switch (section) {
      case INFO:
        key = strtok (line, "=");
        value = strtok (NULL, "=");

        SC_CHECK_ABORT (key != NULL && value != NULL,
                        "entries in the [Forest Info] setion must be\n"
                        "key value pairs, i.e. key=value");

        p4est_trim_ending_whitespace (key);

        if (strcmp (key, "Nk") == 0) {
          num_trees = atoi (value);
          set_num_trees = 1;
        }

        else if (strcmp (key, "Nv") == 0) {
          num_vertices = atoi (value);
          set_num_vertices = 1;
        }

        else if (strcmp (key, "Nve") == 0) {
          num_vtt = atoi (value);
          set_num_vtt = 1;
        }

        if (set_num_vertices && set_num_trees && set_num_vtt
            && *connectivity == NULL) {
          *connectivity = p4est_connectivity_new (num_trees, num_vertices,
                                                  num_vtt);
          tree_to_vertex = (*connectivity)->tree_to_vertex;
          tree_to_tree = (*connectivity)->tree_to_tree;
          tree_to_face = (*connectivity)->tree_to_face;
          vertices = (*connectivity)->vertices;
          vtt_offset = (*connectivity)->vtt_offset;
          vtt_offset[0] = 0;
          vertex_to_tree = (*connectivity)->vertex_to_tree;
          vertex_to_vertex = (*connectivity)->vertex_to_vertex;
        }

        break;
      case COORD:
        sscanf (line, "%d %lf %lf %lf", &k, &vx, &vy, &vz);
        --k;

        SC_CHECK_ABORT (k >= 0 && k < num_vertices, "Bad [] entry");

        vertices[k * 3 + 0] = vx;
        vertices[k * 3 + 1] = vy;
        vertices[k * 3 + 2] = vz;

        break;
      case ETOV:
        sscanf (line, "%d %d %d %d %d", &k, &v0, &v1, &v2, &v3);
        --k;
        --v0;
        --v1;
        --v2;
        --v3;

        SC_CHECK_ABORT (k >= 0 && k < num_trees &&
                        v0 >= 0 && v0 < num_vertices &&
                        v1 >= 0 && v1 < num_vertices &&
                        v2 >= 0 && v2 < num_vertices &&
                        v3 >= 0 && v3 < num_vertices,
                        "Bad [Element to Vertex] entry");

        tree_to_vertex[k * 4 + 0] = v0;
        tree_to_vertex[k * 4 + 1] = v1;
        tree_to_vertex[k * 4 + 2] = v2;
        tree_to_vertex[k * 4 + 3] = v3;

        break;
      case ETOE:
        sscanf (line, "%d %d %d %d %d", &k, &k0, &k1, &k2, &k3);
        --k;
        --k0;
        --k1;
        --k2;
        --k3;

        SC_CHECK_ABORT (k >= 0 && k < num_trees &&
                        k0 >= 0 && k0 < num_trees &&
                        k1 >= 0 && k1 < num_trees &&
                        k2 >= 0 && k2 < num_trees &&
                        k3 >= 0 && k3 < num_trees,
                        "Bad [Element to Element] entry");

        tree_to_tree[k * 4 + 0] = k0;
        tree_to_tree[k * 4 + 1] = k1;
        tree_to_tree[k * 4 + 2] = k2;
        tree_to_tree[k * 4 + 3] = k3;

        break;
      case ETOF:
        sscanf (line, "%d %d %d %d %d", &k, &f0, &f1, &f2, &f3);
        --k;
        --f0;
        --f1;
        --f2;
        --f3;

        SC_CHECK_ABORT (k >= 0 && k < num_trees &&
                        f0 >= 0 && f0 < 4 &&
                        f1 >= 0 && f1 < 4 &&
                        f2 >= 0 && f2 < 4 &&
                        f3 >= 0 && f3 < 4, "Bad [Element to Face] entry");

        tree_to_face[k * 4 + 0] = (int8_t) f0;
        tree_to_face[k * 4 + 1] = (int8_t) f1;
        tree_to_face[k * 4 + 2] = (int8_t) f2;
        tree_to_face[k * 4 + 3] = (int8_t) f3;

        break;
      case VTOE:
        value = strtok (line, " \t");
        v0 = atoi (value);
        --v0;
        value = strtok (NULL, " \t");
        Nnn = atoi (value);
        vtt_offset[v0 + 1] = vtt_offset[v0] + Nnn;
        for (i = 0; i < Nnn; ++i) {
          value = strtok (NULL, " \t");
          vertex_to_tree[vtt_offset[v0] + i] = atoi (value) - 1;
        }

        break;
      case ET:
        break;
      case FT:
        break;
      case CF:
        break;
      case CT:
        break;
      case NONE:
        SC_CHECK_ABORT (0, "Mesh file must start with a section");
      default:
        SC_CHECK_ABORT (0, "Unknown section in mesh file");
      }

      ++section_lines_read;
    }
  }

  /* populate vertex_to_vertex array since it is not read in currently */
  for (vertex = 0; vertex < num_vertices; ++vertex) {
    corner_trees = vtt_offset[vertex + 1] - vtt_offset[vertex];
    for (ctree = 0; ctree < corner_trees; ++ctree) {
      vertex_to_vertex[vtt_offset[vertex] + ctree] = vertex;
    }
  }

  retval = fclose (file);
  if (retval) {
    fprintf (stderr, "Failed to close p4est mesh file %s (%d:%d)\n", filename,
             retval, EOF);
    p4est_connectivity_destroy (*connectivity);
    *connectivity = NULL;
    return 1;
  }

  if (!p4est_connectivity_is_valid (*connectivity)) {
    fprintf (stderr, "Mesh file %s connectivity strucure is invalid\n",
             filename);
    p4est_connectivity_destroy (*connectivity);
    *connectivity = NULL;
    return 1;
  }

  return 0;
}

void
p4est_connectivity_print (p4est_connectivity_t * connectivity, FILE * nout)
{
  int                 i, j, k, Nnn, num_trees, num_vertices, num_vtt;
  int32_t            *tree_to_vertex, *tree_to_tree, *vtt_offset,
    *vertex_to_tree;
  int8_t             *tree_to_face;

  P4EST_ASSERT (p4est_connectivity_is_valid (connectivity));

  num_trees = connectivity->num_trees;
  num_vertices = connectivity->num_vertices;
  num_vtt = connectivity->vtt_offset[num_vertices];

  tree_to_vertex = connectivity->tree_to_vertex;
  tree_to_tree = connectivity->tree_to_tree;
  tree_to_face = connectivity->tree_to_face;

  vtt_offset = connectivity->vtt_offset;
  vertex_to_tree = connectivity->vertex_to_tree;

  fprintf (nout, "[Forest Info]\n");
  fprintf (nout, "ver = 0.0.1  # Version of the forest file\n");
  fprintf (nout, "Nk  = %d     # Number of elements\n", num_trees);
  fprintf (nout, "Nv  = %d     # Number of mesh vertices\n", num_vertices);
  fprintf (nout,
           "Nve = %d     # Number of trees in the vertex to element list\n",
           num_vtt);
  fprintf (nout, "Net = 0      # Number of element tags\n");
  fprintf (nout, "Nft = 0      # Number of face tags\n");
  fprintf (nout, "Ncf = 0      # Number of curved faces\n");
  fprintf (nout, "Nct = 0      # Number of curved types\n");
  fprintf (nout, "\n");
  fprintf (nout, "[Coordinates of Element Vertices]\n");
  fprintf (nout, "[Element to Vertex]\n");
  for (k = 0; k < num_trees; ++k)
    printf ("    %d    %d    %d    %d    %d\n",
            k + 1, tree_to_vertex[4 * k + 0] + 1,
            tree_to_vertex[4 * k + 1] + 1, tree_to_vertex[4 * k + 2] + 1,
            tree_to_vertex[4 * k + 3] + 1);
  fprintf (nout, "[Element to Element]\n");
  for (k = 0; k < num_trees; ++k)
    printf ("    %d    %d    %d    %d    %d\n",
            k + 1, tree_to_tree[4 * k + 0] + 1, tree_to_tree[4 * k + 1] + 1,
            tree_to_tree[4 * k + 2] + 1, tree_to_tree[4 * k + 3] + 1);
  fprintf (nout, "[Element to Face]\n");
  for (k = 0; k < num_trees; ++k)
    printf ("    %d    %d    %d    %d    %d\n",
            k + 1, tree_to_face[4 * k + 0] + 1, tree_to_face[4 * k + 1] + 1,
            tree_to_face[4 * k + 2] + 1, tree_to_face[4 * k + 3] + 1);
  for (i = 0; i < num_vertices; ++i) {
    Nnn = vtt_offset[i + 1] - vtt_offset[i];
    printf ("    %d   %d", i + 1, Nnn);
    for (j = 0; j < Nnn; ++j) {
      printf ("    %d", vertex_to_tree[vtt_offset[i] + j] + 1);
    }
    printf ("\n");
  }
  fprintf (nout, "[Element Tags]\n");
  fprintf (nout, "[Face Tags]\n");
  fprintf (nout, "[Curved Faces]\n");
  fprintf (nout, "[Curved Types]\n");
}

/* EOF p4est_file.h> */
