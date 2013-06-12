/*

zero.c

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


int zero(O3Data *od, int type, double level)
{
  int i;
  int j;
  int k;
  int zeroed_x_value_count;
  int result;
  double threshold;
  double value;
  

  threshold = fabs(level);
  result = calc_active_vars(od, CV_MODEL);
  if (result) {
    return result;
  }
  tee_printf(od, "TYPE      AFFECTED\n");
  for (i = 0; i < od->field_num; ++i) {
    zeroed_x_value_count = 0;
    if (get_field_attr(od, i, OPERATE_BIT)) {
      for (j = 0; j < od->object_num; ++j) {
        if (get_object_attr(od, j, ACTIVE_BIT)) {
          for (k = 0; k < od->x_vars; ++k) {
            if (get_x_var_attr(od, i, k, ACTIVE_BIT)) {
              result = get_x_value(od, i, j, k, &value, 0);
              if (result) {
                return result;
              }
              if (MISSING(value)) {
                continue;
              }
              result = get_x_value(od, i, j, k,
                &value, CUTOFF_BIT);
              if (result) {
                return result;
              }
              switch (type) {
                case ZERO_NEG:
                if ((value < 0.0) && (value > -threshold)) {
                  result = set_x_value(od, i, j, k, 0.0);
                  if (result) {
                    return result;
                  }
                  ++zeroed_x_value_count;
                }
                break;
                
                case ZERO_POS:
                if ((value > 0.0) && (value < threshold)) {
                  result = set_x_value(od, i, j, k, 0.0);
                  if (result) {
                    return result;
                  }
                  ++zeroed_x_value_count;
                }
                break;
                
                default:
                if ((fabs(value) > ALMOST_ZERO) && (fabs(value) < threshold)) {
                  result = set_x_value(od, i, j, k, 0.0);
                  if (result) {
                    return result;
                  }
                  ++zeroed_x_value_count;
                }
                break;
              }            
            }
          }
        }
      }
      value = (double)(od->mel.x_data[i].active_x_vars *
        od->active_object_num);
      tee_printf(od, "x%-9d%-9d (%3d%%)\n", i + 1, zeroed_x_value_count,
        (value ? (int)safe_rint((double)zeroed_x_value_count
        / value * 100.0) : 0));
    }
  }

  result = calc_active_vars(od, FULL_MODEL);
  
  return result;
}
