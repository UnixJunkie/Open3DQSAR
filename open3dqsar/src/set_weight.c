/*

set_weight.c

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


void set_field_weight(O3Data *od, double weight)
{
  int i;
  

  for (i = 0; i < od->field_num; ++i) {
    if (get_field_attr(od, i, OPERATE_BIT)) {
      od->mel.x_data[i].x_weight_coefficient = weight;
    }
  }
}


int set_object_weight(O3Data *od, double weight, int list_type, int options)
{
  char buffer[BUF_LEN];
  int object_num;
  int struct_num;
  int conf_num;
  int j;
  int found;
  

  memset(buffer, 0, BUF_LEN);
  if (!(list_type & (1 << FROM_FILE))) {
    for (object_num = 0; object_num < od->object_num; ++object_num) {
      struct_num = od->al.mol_info[object_num]->struct_num;
      conf_num = 1;
      if (get_object_attr(od, object_num, OPERATE_BIT)) {
        switch (options) {
          case RANDOM_WEIGHTS:
          od->mel.object_weight[object_num] = genrand_real(od);
          break;
          
          case EVEN_WEIGHTS:
          od->mel.object_weight[object_num] = 1.0 / (double)conf_num;
          break;
          
          default:
          od->mel.object_weight[object_num] = weight;
          break;
        }
      }
    }
  }
  else {
    while (fgets(buffer, BUF_LEN, od->file[ASCII_IN]->handle)) {
      buffer[BUF_LEN - 1] = '\0';
      sscanf(buffer, "%d %lf", &object_num, &weight);
      if (object_num < 1) {
        return INVALID_LIST_RANGE;
      }
      if ((list_type & (1 << OBJECT_LIST)) && (object_num > od->object_num)) {
         return INVALID_LIST_RANGE;
      }
      if (list_type & (1 << ID_LIST)) {
        for (j = 0, found = 0; (!found) && (j < od->grid.object_num); ++j) {
          if (object_num == od->al.mol_info[j]->object_id) {
            object_num = j + 1;
            found = 1;
          }
        }
        if (!found) {
          return INVALID_LIST_RANGE;
        }
      }
      if (weight < 0) {
        return WRONG_DATA_FORMAT;
      }
      od->mel.object_weight[object_num - 1] = weight;
    }
  }
  
  return 0;
}


void set_y_var_weight(O3Data *od, double weight)
{
  int i;
  

  for (i = 0; i < od->y_vars; ++i) {
    od->mel.y_data[i].y_weight_coefficient = weight;
  }
}
