/*

transform.c

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


double perform_operation(int operation, double value, double factor)
{
  switch (operation) {
    case MULTIPLY:
    value *= factor;
    break;
    
    case DIVIDE:
    value /= factor;
    break;
    
    case SUM:
    value += factor;
    break;
    
    case SUBTRACT:
    value -= factor;
    break;
    
    case OPPOSITE:
    value = -value;
    break;
    
    case ABSOLUTE_VALUE:
    value = fabs(value);
    break;
    
    case POWER:
    value = pow(value, factor);
    break;
    
    case LOGARITHM_BASE_10:
    value = log10(value);
    break;
    
    case LOGARITHM_BASE_E:
    value = log(value);
    break;
    
    case LOGARITHM_BASE_N:
    value = log(value);
    value /= log(factor);
    break;
  }
  
  return value;
}


int transform(O3Data *od, int type, int operation, double factor)
{
  int i;
  int j;
  int k;
  int result;
  double value;


  if (type == 'X') {
    for (i = 0; i < od->field_num; ++i) {
      if (get_field_attr(od, i, OPERATE_BIT)) {
        for (j = 0; j < od->object_num; ++j) {
          if (get_object_attr(od, j, ACTIVE_BIT)) {
            for (k = 0; k < od->x_vars; ++k) {
              result = get_x_value(od, i, j, k, &value, 0);
              if (result) {
                return result;
              }
              if (!MISSING(value)) {
                result = get_x_value(od, i, j, k,
                  &value, CUTOFF_BIT);
                if (result) {
                  return result;
                }
                value = perform_operation(operation, value, factor);
                result = set_x_value(od, i, j, k, value);
                if (result) {
                  return result;
                }
              }
            }
          }
        }
      }
    }
  }
  else {
    for (j = 0; j < od->object_num; ++j) {
      if (get_object_attr(od, j, ACTIVE_BIT)) {
        for (k = 0; k < od->y_vars; ++k) {
          if (get_y_var_attr(od, k, OPERATE_BIT)) {
            value = get_y_value(od, j, k, 0);
            value = perform_operation(operation, value, factor);
            set_y_value(od, j, k, value);
          }
        }
      }
    }
  }
  result = calc_active_vars(od, FULL_MODEL);
  
  return result;
}
