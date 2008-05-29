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

#include <p4est_base.h>
#include <p4est.h>

/* #define THEBIGTEST */

static int
compar (const void *p1, const void *p2)
{
  int                 i1 = *(int *) p1;
  int                 i2 = *(int *) p2;

  return i1 - i2;
}

int
main (int argc, char **argv)
{
  int                 rank;
  int                 i, i1, i2, i3, i3last, i4, i4last, temp, count;
  int                 s, swaps1, swaps2, swaps3, total1, total2, total3;
  int                 index;
  int                *pi;
  sc_array_t         *a1, *a2, *a3, *a4;
  int                 mpiret;
  double              start, elapsed_pqueue, elapsed_qsort;

  mpiret = MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  mpiret = MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  SC_CHECK_MPI (mpiret);

  sc_init (rank, NULL, NULL, NULL, SC_LP_DEFAULT);

  a1 = sc_array_new (sizeof (int));
  a2 = sc_array_new (sizeof (int));
  a3 = sc_array_new (sizeof (int));
  a4 = sc_array_new (sizeof (int));

#ifdef THEBIGTEST
  count = 325323;
#else
  count = 3251;
#endif
  printf ("Test pqueue with count %d\n", count);

  start = -MPI_Wtime ();

  swaps1 = swaps2 = swaps3 = 0;
  total1 = total2 = total3 = 0;
  for (i = 0; i < count; ++i) {
    sc_array_resize (a1, i + 1);
    pi = sc_array_index (a1, i);
    *pi = i;
    s = sc_array_pqueue_add (a1, &temp, compar);
    swaps1 += ((s > 0) ? 1 : 0);
    total1 += s;

    sc_array_resize (a2, i + 1);
    pi = sc_array_index (a2, i);
    *pi = count - i - 1;
    s = sc_array_pqueue_add (a2, &temp, compar);
    swaps2 += ((s > 0) ? 1 : 0);
    total2 += s;

    sc_array_resize (a3, i + 1);
    pi = sc_array_index (a3, i);
    *pi = (15 * i) % 172;
    s = sc_array_pqueue_add (a3, &temp, compar);
    swaps3 += ((s > 0) ? 1 : 0);
    total3 += s;
  }
  SC_CHECK_ABORT (swaps1 == 0 && total1 == 0, "pqueue_add");

  printf ("   Swaps %d %d %d Total %d %d %d\n",
          swaps1, swaps2, swaps3, total1, total2, total3);

  temp = 52;
  index = sc_array_bsearch (a1, &temp, compar);
  SC_CHECK_ABORT (index != -1, "array_bsearch_index");
  pi = sc_array_index (a1, index);
  SC_CHECK_ABORT (*pi == temp, "array_bsearch");

  i3last = -1;
  swaps1 = swaps2 = swaps3 = 0;
  total1 = total2 = total3 = 0;
  for (i = 0; i < count; ++i) {
    s = sc_array_pqueue_pop (a1, &i1, compar);
    swaps1 += ((s > 0) ? 1 : 0);
    total1 += s;

    s = sc_array_pqueue_pop (a2, &i2, compar);
    swaps2 += ((s > 0) ? 1 : 0);
    total2 += s;

    s = sc_array_pqueue_pop (a3, &i3, compar);
    swaps3 += ((s > 0) ? 1 : 0);
    total3 += s;

    SC_CHECK_ABORT (i == i1 && i == i2, "pqueue_pop");
    SC_CHECK_ABORT (i3 >= i3last, "pqueue_pop");
    i3last = i3;
  }
  printf ("   Swaps %d %d %d Total %d %d %d\n",
          swaps1, swaps2, swaps3, total1, total2, total3);

  elapsed_pqueue = start + MPI_Wtime ();

  sc_array_destroy (a1);
  sc_array_destroy (a2);
  sc_array_destroy (a3);

  printf ("Test array sort with count %d\n", count);

  start = -MPI_Wtime ();

  /* the resize is done to be comparable with the above procedure */
  for (i = 0; i < count; ++i) {
    sc_array_resize (a4, i + 1);
    pi = sc_array_index (a4, i);
    *pi = (15 * i) % 172;
  }
  sc_array_sort (a4, compar);

  i4last = -1;
  for (i = 0; i < count; ++i) {
    pi = sc_array_index (a4, i);
    i4 = *pi;

    SC_CHECK_ABORT (i4 >= i4last, "array_sort");
    i4last = i4;
  }
  sc_array_resize (a4, 0);

  elapsed_qsort = start + MPI_Wtime ();
  printf ("Test timings pqueue %g qsort %g\n",
          elapsed_pqueue, 3. * elapsed_qsort);

  sc_array_destroy (a4);
  sc_finalize ();

  mpiret = MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}

/* EOF test_pqueue.c */
