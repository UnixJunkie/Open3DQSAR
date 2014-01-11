/*

autoscale.c

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


int autoscale_field(O3Data *od)
{
  int i;
  int x;
  int y;
  int result;
  double value;
  

  for (i = 0; i < od->field_num; ++i) {
    if (get_field_attr(od, i, OPERATE_BIT)) {
      result = stddev_x_var(od, i);
      if (result) {
        return result;
      }
      for (y = 0; y < od->object_num; ++y) {
        for (x = 0; x < od->x_vars; ++x) {
          result = get_x_value(od, i, y, x, &value, 0);
          if (result) {
            return result;
          }
          if (MISSING(value)) {
            value = 0.0;
          }
          else {
            result = get_x_value(od, i, y, x, &value, CUTOFF_BIT);
            if (result) {
              return result;
            }
            value -= get_x_var_buf(od, i, x, AVE_BUF);
            value /= get_x_var_buf(od, i, x, STDDEV_BUF);
          }
          result = set_x_value(od, i, y, x, value);
          if (result) {
            return result;
          }
        }
      }
    }
  }
  
  return 0;
}


int autoscale_y_var(O3Data *od)
{
  int x;
  int y;
  double value;
  

  stddev_y_var(od);
  for (y = 0; y < od->object_num; ++y) {
    for (x = 0; x < od->y_vars; ++x) {
      if (get_y_var_attr(od, x, OPERATE_BIT)) {
        value = get_y_value(od, y, x, 0);
        value -= get_y_var_buf(od, x, AVE_BUF);
        value /= get_y_var_buf(od, x, STDDEV_BUF);
        set_y_value(od, y, x, value);
      }
    }
  }
  
  return 0;
}
