/*

double_vec.c

is part of

Open3DQSAR
----------

An open-source software aimed at high-throughput
chemometric analysis of molecular interaction fields

Copyright (C) 2009-2018 Paolo Tosco, Thomas Balle

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.

For further information, please contact:

Paolo Tosco, PhD
Dipartimento di Scienza e Tecnologia del Farmaco
Universita' degli Studi di Torino
Via Pietro Giuria, 9
10125 Torino (Italy)
Phone:  +39 011 670 7680
Mobile: +39 348 553 7206
Fax:    +39 011 670 7687
E-mail: paolo.tosco@unito.it

*/

#include <include/o3header.h>


DoubleVec *double_vec_resize(DoubleVec *double_vec, int size)
{
  if (!double_vec) {
    double_vec = (DoubleVec *)malloc(sizeof(DoubleVec));
    if (!double_vec) {
      return NULL;
    }
    memset(double_vec, 0, sizeof(DoubleVec));
  }
  if (size > double_vec->max_size) {
    double_vec->ve = (double *)realloc
      (double_vec->ve, size * sizeof(double));
    if (!(double_vec->ve)) {
      return NULL;
    }
    memset(&(double_vec->ve[double_vec->size]), 0,
      (size - double_vec->size) * sizeof(double));
    double_vec->max_size = size;
  }
  double_vec->size = size;
  
  return double_vec;
}


void double_vec_free(DoubleVec *double_vec)
{
  if (double_vec) {
    if (double_vec->ve) {
      free(double_vec->ve);
    }
    free(double_vec);
  }
}


DoubleVec *double_vec_sort(DoubleVec *x, IntPerm *order)
{
  int i = 0;
  int j = 0;
  int l;
  int r;
  int tmp_i = 0;
  int sp;
  int stack[MAX_STACK];
  double tmp = 0.0;
  double v;


  if ((!x) || (!order)) {
    return NULL;
  }
  if (order->size != x->size) {
    return NULL;
  } 
  for (i = 0; i < order->size; ++i) {
    order->pe[i] = i;
  }
  if (x->size <= 1) {
    return x;
  }
  /*
  adapted from
  quicksort algorithm developed by Robert Sedgewick,
       "Algorithms in C", 1990, Chapter 9, pp. 118-122.
       */
  sp = 0;
  l = 0;
  r = x->size - 1;
  v = x->ve[0];
  for (;;) {
    while (r > l) {
      v = x->ve[r];
      i = l - 1;
      j = r;
      for (;;) {
        while (x->ve[++i] < v);
        --j;
        while ((x->ve[j] > v) && (j != 0)) {
          --j;
        }
        if (i >= j) {
          break;
        }

        tmp = x->ve[i];
        x->ve[i] = x->ve[j];
        x->ve[j] = tmp;
        if (order) {
          tmp_i = order->pe[i];
          order->pe[i] = order->pe[j];
          order->pe[j] = tmp_i;
        }
      }
      tmp = x->ve[i];
      x->ve[i] = x->ve[r];
      x->ve[r] = tmp;
      if (order) {
        tmp_i = order->pe[i];
        order->pe[i] = order->pe[r];
        order->pe[r] = tmp_i;
      }

      if ((i - l) > (r - i)) {
        stack[sp++] = l;
        stack[sp++] = i - 1;
        l = i + 1;
      }
      else {
        stack[sp++] = i + 1;
        stack[sp++] = r;
        r = i - 1;
      }
    }

    if (sp == 0) {
      break;
    }
    r = stack[--sp];
    l = stack[--sp];
  }
  /*
  to make the sort stable,
  reorder the permutation vector so that
  in the presence of identical elements
  the one with the lower index comes first
  */
  i = 0;
  while (i < x->size) {
    if (i) {
      /*
      if the difference between this vector element
      and the previous one is less than ALMOST_ZERO,
      then look for similar ones, then sort
      the corresponding segment in the permutation
      vector
      */
      j = 0;
      while (fabs(x->ve[i] - tmp) < ALMOST_ZERO) {
        if (!j) {
          tmp_i = i - 1;
        }
        ++j;
        ++i;
        if (i == x->size) {
          break;
        }
      }
      if (j) {
        qsort(&(order->pe[tmp_i]), j + 1,
          sizeof(int), compare_integers);
      }
    }
    if ((!j) || (!i)) {
      tmp = x->ve[i];
      ++i;
    }
  }

  return x;
}
