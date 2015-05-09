/*

calc_active_vars.c

is part of

Open3DQSAR
----------

An open-source software aimed at high-throughput
chemometric analysis of molecular interaction fields

Copyright (C) 2009-2015 Paolo Tosco, Thomas Balle

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


int calc_active_vars(O3Data *od, int model_type)
{
  char format[BUF_LEN];
  int i;
  int j;
  int k;
  int first_max = 0;
  int first_min = 0;
  int max_scientific;
  int min_scientific;
  int active_count;
  int result;
  double value;
  double stddev;


  od->overall_active_x_vars = 0;
  od->overall_zero_x_values = 0;
  od->overall_zero_y_values = 0;
  memset(format, 0, BUF_LEN);
  for (i = 0; i < od->field_num; ++i) {
    if (get_field_attr(od, i, ACTIVE_BIT)) {
      active_count = 0;
      result = stddev_x_var(od, i);
      if (result) {
        return result;
      }
      for (j = 0; j < od->x_vars; ++j) {
        if (get_x_var_attr(od, i, j, DELETE_BIT)) {
          set_x_var_attr(od, i, j, ACTIVE_BIT, 0);
          continue;
        }
        stddev = get_x_var_buf(od, i, j, STDDEV_BUF);
        if (stddev > od->mel.x_data[i].sdcut_x_var) {
          set_x_var_attr(od, i, j, ACTIVE_BIT, 1);
          ++active_count;
        }
        else {
          set_x_var_attr(od, i, j, ACTIVE_BIT, 0);
        }
      }
      od->mel.x_data[i].active_x_vars = active_count;
      od->overall_active_x_vars += active_count;
    }
  }
  if (od->y_vars) {
    stddev_y_var(od);
    for (i = 0; i < od->y_vars; ++i) {
      stddev = get_y_var_buf(od, i, STDDEV_BUF);
      /*
      if (stddev < MIN_Y_VAR_SD) {
        return Y_VAR_LOW_SD;
      }
      */
    }
  }
  if (model_type & FULL_MODEL) {
    for (i = 0; i < od->field_num; ++i) {
      if (get_field_attr(od, i, ACTIVE_BIT)) {
        od->mel.x_data[i].zero_x_values = 0;
        od->mel.x_data[i].max_x_value = 0.0;
        od->mel.x_data[i].min_x_value = 0.0;
        first_max = 0;
        first_min = 0;
        for (k = 0; k < od->object_num; ++k) {
          if (get_object_attr(od, k, ACTIVE_BIT)) {
            for (j = 0; j < od->x_vars; ++j) {
              result = get_x_value(od, i, k, j, &value, 0);
              if (result) {
                return result;
              }
              if (MISSING(value)) {
                value = 0.0;
              }
              else {
                result = get_x_value(od, i, k, j, &value, CUTOFF_BIT);
                if (result) {
                  return result;
                }
              }
              if ((!first_max) || (value > od->mel.x_data[i].max_x_value)) {
                first_max = 1;
                od->mel.x_data[i].max_x_value = value;
              }
              if ((!first_min) || (value < od->mel.x_data[i].min_x_value)) {
                first_min = 1;
                od->mel.x_data[i].min_x_value = value;
              }
              if (get_x_var_attr(od, i, j, ACTIVE_BIT)
                && (fabs(value) < SMALL_ENERGY_VALUE)) {
                ++(od->mel.x_data[i].zero_x_values);
              }
            }
          }
        }
        od->overall_zero_x_values +=
          od->mel.x_data[i].zero_x_values;
      }
    }
    for (i = 0; i < od->y_vars; ++i) {
      od->mel.y_data[i].zero_y_values = 0;
      od->mel.y_data[i].max_y_value = 0.0;
      od->mel.y_data[i].min_y_value = 0.0;
      first_max = 0;
      first_min = 0;
      for (j = 0; j < od->grid.object_num; ++j) {
        if (get_object_attr(od, j, ACTIVE_BIT)) {
          value = get_y_value(od, j, i, 0);
          if ((!first_max) || (value > od->mel.y_data[i].max_y_value)) {
            first_max = 1;
            od->mel.y_data[i].max_y_value = value;
          }
          if ((!first_min) || (value < od->mel.y_data[i].min_y_value)) {
            first_min = 1;
            od->mel.y_data[i].min_y_value = value;
          }
          if (fabs(value) < ALMOST_ZERO) {
            ++(od->mel.y_data[i].zero_y_values);
          }
        }
      }
      od->overall_zero_y_values +=
        od->mel.y_data[i].zero_y_values;
    }
    tee_printf(od, "\n");
    tee_printf(od, "x variables:           %d\n", od->x_vars * od->field_num);
    tee_printf(od, "y variables:           %d\n", od->y_vars);
    tee_printf(od, "\n");
    tee_printf(od, "TYPE    STATUS      TOTAL         ACTIVE               MIN           MAX           WEIGHT\n");
    for (i = 0; i < od->field_num; ++i) {
      sprintf(format, "%lf", od->mel.x_data[i].max_x_value);
      max_scientific = (strlen(format) > 12) ? 'e' : 'f';
      sprintf(format, "%lf", od->mel.x_data[i].min_x_value);
      min_scientific = (strlen(format) > 12) ? 'e' : 'f';
      sprintf(format, "x%%-7d%%-12s%%-14d%%-8d (%%3d%%%%)      %%-14.4l%c%%-14.4l%c%%-14.4lf\n",
        min_scientific, max_scientific);
      tee_printf(od, format, i + 1,
        (get_field_attr(od, i, ACTIVE_BIT) ? "INCLUDED" : "EXCLUDED"),
        od->x_vars, od->mel.x_data[i].active_x_vars, od->x_vars
        ? (int)safe_rint((double)(od->mel.x_data[i].active_x_vars) /
        (double)(od->x_vars) * (double)100) : 0,
        od->mel.x_data[i].min_x_value,
        od->mel.x_data[i].max_x_value,
        od->mel.x_data[i].x_weight_coefficient);
    }
    for (i = 0; i < od->y_vars; ++i) {
      max_scientific = (od->mel.y_data[i].max_y_value >= 1000.0) ? 'e' : 'f';
      min_scientific = (od->mel.y_data[i].min_y_value <= -1000.0) ? 'e' : 'f';
      sprintf(format, "y%%-7d%%-12s%%-14d%%-8d (%%3d%%%%)      %%-14.4l%c%%-14.4l%c%%-14.4lf\n",
        min_scientific, max_scientific);
      tee_printf(od, format, i + 1, "INCLUDED",
        1, 1, 100,
        od->mel.y_data[i].min_y_value,
        od->mel.y_data[i].max_y_value, (double)1);
    }
    tee_printf(od, "\n");
    tee_printf(od, "Active x variables:    %8d (%3d%%)\n", od->overall_active_x_vars,
      (od->field_num && od->x_vars) ? (int)safe_rint
      ((double)((double)(od->overall_active_x_vars) /
      ((double)(od->field_num) * (double)(od->x_vars)) * (double)100)) : 0);
    if (od->field_num) {
      tee_printf(od, "Zero x values:         %8d (%3d%%)\n", od->overall_zero_x_values,
        (od->active_object_num && od->x_vars)
        ? (int)safe_rint((double)((double)(od->overall_zero_x_values) /
        ((double)(od->field_num) * (double)(od->active_object_num) *
        (double)(od->x_vars)) * (double)100)) : 0);
    }
    tee_printf(od, "Active y variables:    %8d\n", od->y_vars);
    if (od->y_vars) {
      tee_printf(od, "Zero y values:         %8d (%3d%%)\n",
        od->overall_zero_y_values,
        od->active_object_num
        ? (int)safe_rint((double)((double)(od->overall_zero_y_values) /
        ((double)(od->active_object_num) *
        (double)(od->y_vars)) * (double)100)) : 0);
    }
    tee_printf(od, "\n");
  }

  return 0;
}
