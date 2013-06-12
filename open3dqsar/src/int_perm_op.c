/*

int_perm_op.c

is part of

Open3DQSAR
----------

An open-source software aimed at high-throughput
chemometric analysis of molecular interaction fields

Copyright (C) 2009-2013 Paolo Tosco, Thomas Balle

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


DoubleMat *int_perm_rows(IntPerm *perm,
  DoubleMat *double_mat1, DoubleMat *double_mat2)
{
  int x;
  int y;
  
  
  if ((!perm) || (!double_mat1) || (!double_mat2)) {
      return NULL;
  }
  if ((double_mat1->m != double_mat2->m)
    || (double_mat1->n != double_mat2->n)
    || (perm->size != double_mat1->m)) {
    return NULL;
  }
  for (y = 0; y < perm->size; ++y) {
    for (x = 0; x < double_mat2->n; ++x) {
      M_POKE(double_mat2, y, x,
        M_PEEK(double_mat1, perm->pe[y], x));
    }
  }
  
  return double_mat2;
}


DoubleVec *int_perm_vec(IntPerm *perm,
  DoubleVec *double_vec1, DoubleVec *double_vec2)
{
  int i;
  
  
  if ((!perm) || (!double_vec1) || (!double_vec2)) {
    return NULL;
  }
  if ((double_vec1->size < perm->size)
    || (double_vec2->size < perm->size)) {
    return NULL;
  }
  for (i = 0; i < perm->size; ++i) {
    double_vec2->ve[i] = double_vec1->ve[perm->pe[i]];
  }
  
  return double_vec2;
}
