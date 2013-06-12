/*

get_value.c

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


int get_x_value(O3Data *od, int field_num,
  int object_num, int x_var, double *value, int flag)
{
  int cutoff_done = 0;
  double double_value;
  double weight = 1.0;
  
  
  if (flag & WEIGHT_BIT) {
    weight = od->mel.x_data[field_num].x_weight_coefficient;
  }
  if (get_x_var_attr(od, field_num, x_var, ACTIVE_BIT)
    || (!(flag & CHECK_IF_ACTIVE_BIT))) {
    if (od->save_ram) {
      if (check_mmap(od, field_num)) {
        return OUT_OF_MEMORY;
      }
    }
    double_value = (double)(od->mel.x_var_array
      [field_num][object_num][x_var]);
    if (flag & CUTOFF_BIT) {
      if (double_value > od->mel.x_data[field_num].max_cutoff) {
        *value = od->mel.x_data[field_num].max_cutoff * weight;
        cutoff_done = 1;
      }
      else if (double_value < od->mel.x_data[field_num].min_cutoff) {
        *value = od->mel.x_data[field_num].min_cutoff * weight;
        cutoff_done = 1;
      }
    }
    if (!cutoff_done) {
      *value = double_value * weight;
    }
    /*
    if (fabs(*value) < SMALL_ENERGY_VALUE) {
      *value = 0.0;
    }
    */
  }
  else {
    *value = 0.0;
  }

  return 0;
}


double get_y_value(O3Data *od, int object_num, int y_var, int flag)
{
  double weight = 1.0;

  
  if (flag & WEIGHT_BIT) {
    weight = od->mel.y_data[y_var].y_weight_coefficient;
  }
  return ((double)(od->mel.y_var_array
    [od->y_vars * object_num + y_var]) * weight);
}
