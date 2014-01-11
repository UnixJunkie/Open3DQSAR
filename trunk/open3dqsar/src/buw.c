/*

buw.c

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


int x_var_buw(O3Data *od)
{
  int i;
  int j;
  int k;
  int result;
  int ref_field = -1;
  double value;
  double sum_x_var_buw;
  

  od->vel.ss = double_vec_resize
    (od->vel.ss, od->field_num);
  if (!(od->vel.ss)) {
    return OUT_OF_MEMORY;
  }
  memset(od->vel.ss->ve, 0,
    od->vel.ss->size * sizeof(double));
  
  sum_x_var_buw = 0.0;
  for (i = 0; i < od->field_num; ++i) {
    if (get_field_attr(od, i, ACTIVE_BIT)) {
      if (ref_field == -1) {
        ref_field = i;
      }
      result = average_x_var(od, i);
      if (result) {
        return result;
      }
      for (j = 0; j < od->object_num; ++j) {
        if (get_object_attr(od, j, ACTIVE_BIT)) {
          for (k = 0; k < od->x_vars; ++k) {
            result = get_x_value(od, i, j, k,
              &value, 0);
            if (result) {
              return result;
            }
            if (MISSING(value)) {
              value = 0.0;
            }
            else {
              result = get_x_value(od, i, j, k,
                &value, CUTOFF_BIT);
              if (result) {
                return result;
              }
              value -= get_x_var_buf(od, i, k, AVE_BUF);
            }
            od->vel.ss->ve[i] += (square(value)
              * od->mel.object_weight[j]);
          }
        }
      }
      od->mel.x_data[i].temp_x_weight_coefficient = ((od->vel.ss->ve[i] > 0.0)
        ? sqrt(od->vel.ss->ve[ref_field] / od->vel.ss->ve[i]) : 0.0);
      sum_x_var_buw +=
        od->mel.x_data[i].temp_x_weight_coefficient;
    }
  }
  tee_printf(od, "\n%5s%25s\n", "Field", "BUW coefficient");
  tee_printf(od, "------------------------------\n");
  od->mel.x_data[ref_field].x_weight_coefficient =
    (double)(od->active_field_num) / sum_x_var_buw;
  for (i = 0; i < od->field_num; ++i) {
    if (get_field_attr(od, i, ACTIVE_BIT)) {
      od->mel.x_data[i].x_weight_coefficient =
        od->mel.x_data[ref_field].x_weight_coefficient
        * od->mel.x_data[i].temp_x_weight_coefficient;
      tee_printf(od, "%5d%25.4lf\n", i + 1,
        od->mel.x_data[i].x_weight_coefficient);
    }
  }
  tee_printf(od, "------------------------------\n");

  return 0;
}


int y_var_buw(O3Data *od)
{
  int i;
  int j;
  double value;
  double sum_y_var_buw;
  

  od->vel.ss = double_vec_resize
    (od->vel.ss, od->y_vars);
  if (!(od->vel.ss)) {
    return OUT_OF_MEMORY;
  }
  memset(od->vel.ss->ve, 0,
    od->vel.ss->size * sizeof(double));
  
  sum_y_var_buw = 0.0;
  average_y_var(od);
  for (i = 0; i < od->y_vars; ++i) {
    for (j = 0; j < od->object_num; ++j) {
      if (get_object_attr(od, j, ACTIVE_BIT)) {
        value = get_y_value(od, j, i, 0);
        value -= get_y_var_buf(od, i, AVE_BUF);
        od->vel.ss->ve[i] += (square(value)
          * od->mel.object_weight[j]);
      }
    }
    od->mel.y_data[i].temp_y_weight_coefficient =
      sqrt(od->vel.ss->ve[0] / od->vel.ss->ve[i]);
    sum_y_var_buw +=
      od->mel.y_data[i].temp_y_weight_coefficient;
  }
  tee_printf(od, "\n%5s%25s\n", "Y var", "BUW coefficient");
  tee_printf(od, "------------------------------\n");
  od->mel.y_data[0].y_weight_coefficient =
    (double)(od->y_vars) / sum_y_var_buw;
  for (i = 0; i < od->y_vars; ++i) {
    od->mel.y_data[i].y_weight_coefficient =
      od->mel.y_data[0].y_weight_coefficient
      * od->mel.y_data[i].temp_y_weight_coefficient;
    tee_printf(od, "%5d%25.4lf\n", i + 1,
      od->mel.y_data[i].y_weight_coefficient);
  }
  tee_printf(od, "------------------------------\n");

  return 0;
}
