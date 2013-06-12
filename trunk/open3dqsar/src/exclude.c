/*

exclude.c

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


int exclude(O3Data *od, int type, int ref_field)
{
  int i;
  int j;
  int k;
  int n;
  int start;
  int end;
  int excluded_x_value_count;
  int result;
  double value;
  

  tee_printf(od, "TYPE      AFFECTED\n");
  for (i = 0; i < od->field_num; ++i) {
    excluded_x_value_count = 0;
    if ((i != ref_field) && get_field_attr(od, i, OPERATE_BIT)) {
      for (k = 0; k < od->x_vars; ++k) {
        for (j = 0; j < od->object_num; ++j) {
          if (type == ANY_OBJECT) {
            start = 0;
            end = od->object_num;
          }
          else {
            start = j;
            end = j + 1;
          }
          for (n = start; n < end; ++n) {
            if (get_object_attr(od, n, ACTIVE_BIT)) {
              result = get_x_value(od, ref_field, n, k, &value, 0);
              if (result) {
                return result;
              }
              if ((value < od->mel.x_data[ref_field].min_cutoff)
                || (value > od->mel.x_data[ref_field].max_cutoff)) {
                result = set_x_value(od, i, j, k, MISSING_VALUE);
                if (result) {
                  return result;
                }
                ++excluded_x_value_count;
                break;
              }
            }
          }
        }
      }
      tee_printf(od, "x%-9d%-9d (%3d%%)\n", i + 1, excluded_x_value_count,
        (int)safe_rint((double)excluded_x_value_count /
        ((double)(od->mel.x_data[i].active_x_vars) *
        (double)(od->active_object_num)) * 100.0));
     }
  }
  result = calc_active_vars(od, FULL_MODEL);
  
  return result;
}
