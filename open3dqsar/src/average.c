/*

average.c

is part of

Open3DQSAR
----------

An open-source software aimed at high-throughput
chemometric analysis of molecular interaction fields

Copyright (C) 2009-2014 Paolo Tosco, Thomas Balle

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


int average_x_var(O3Data *od, int field_num)
{
  int i;
  int j;
  int result;
  double value;
  double ave;
  double n;
  
  
  for (j = 0, n = 0.0, ave = 0.0; j < od->object_num; ++j) {
    if (get_object_attr(od, j, ACTIVE_BIT)) {
      n += od->mel.object_weight[j];
    }
  }
  for (i = 0; i < od->x_vars; ++i) {
    for (j = 0, ave = 0.0; j < od->object_num; ++j) {
      if (get_object_attr(od, j, ACTIVE_BIT)) {
        result = get_x_value(od, field_num, j, i,
          &value, 0);
        if (result) {
          return result;
        }
        if (!MISSING(value)) {
          result = get_x_value(od, field_num, j, i,
            &value, CUTOFF_BIT);
          if (result) {
            return result;
          }
          ave += (value * od->mel.object_weight[j]);
        }
      }
    }
    if (n > 0.0) {
      ave /= n;
    }
    set_x_var_buf(od, field_num, i, AVE_BUF, ave);
  }
  
  return 0;
}


void average_y_var(O3Data *od)
{
  int i;
  int j;
  double ave;
  double n;
  
  
  for (j = 0, n = 0.0, ave = 0.0; j < od->object_num; ++j) {
    if (get_object_attr(od, j, ACTIVE_BIT)) {
      n += od->mel.object_weight[j];
    }
  }
  for (i = 0; i < od->y_vars; ++i) {
    for (j = 0, ave = 0.0; j < od->object_num; ++j) {
      if (get_object_attr(od, j, ACTIVE_BIT)) {
        ave += (get_y_value(od, j, i, WEIGHT_BIT)
          * od->mel.object_weight[j]);
      }
    }
    if (n > 0.0) {
      ave /= n;
    }
    set_y_var_buf(od, i, AVE_BUF, ave);
  }
}
