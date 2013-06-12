/*

int_perm.c

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


void int_perm_free(IntPerm *int_perm)
{
  if (int_perm) {
    if (int_perm->pe) {
      free(int_perm->pe);
    }
    free(int_perm);
  }
}


IntPerm *int_perm_resize(IntPerm *int_perm, int size)
{
  if (!int_perm) {
    int_perm = (IntPerm *)malloc(sizeof(IntPerm));
    if (!int_perm) {
      return NULL;
    }
    memset(int_perm, 0, sizeof(IntPerm));
  }
  if (size > int_perm->max_size) {
    int_perm->pe = (int *)realloc
      (int_perm->pe, size * sizeof(int));
    if (!(int_perm->pe)) {
      return NULL;
    }
    memset(&(int_perm->pe[int_perm->size]), 0,
      (size - int_perm->size) * sizeof(int));
    int_perm->max_size = size;
  }
  int_perm->size = size;
  
  return int_perm;
}
