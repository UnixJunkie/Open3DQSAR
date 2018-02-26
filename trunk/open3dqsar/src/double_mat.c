/*

double_mat.c

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


DoubleMat *double_mat_resize(DoubleMat *double_mat, int m, int n)
{
  if (!double_mat) {
    double_mat = (DoubleMat *)malloc(sizeof(DoubleMat));
    if (!double_mat) {
      return NULL;
    }
    memset(double_mat, 0, sizeof(DoubleMat));
  }
  if ((m > double_mat->max_m) || (n > double_mat->max_n)) {
    double_mat->base = (double *)realloc(double_mat->base,
      m * n * sizeof(double));
    if (!(double_mat->base)) {
      return NULL;
    }
    memset(double_mat->base, 0, m * n * sizeof(double));
    double_mat->max_m = m;
    double_mat->max_n = n;
  }
  double_mat->m = m;
  double_mat->n = n;
  
  return double_mat;
}


void double_mat_free(DoubleMat *double_mat)
{
  if (double_mat) {
    if (double_mat->base) {
      free(double_mat->base);
    }
    free(double_mat);
  }
}
