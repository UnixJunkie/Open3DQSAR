/*

cutoff.c

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


int cutoff(O3Data *od, int type, double cutoff)
{
  int i;
  int j;
  int k;
  int cut_x_value_count;
  int result;
  int first_cond;
  int second_cond;
  double value;
  double j_cutoff = 0.0;
  double k_cutoff = 0.0;
  double sign;
  

  result = calc_active_vars(od, CV_MODEL);
  if (result) {
    return result;
  }
  j_cutoff = cutoff;
  if (type & CUTOFF_QUADRATIC) {
    j_cutoff = (((type & CUTOFF_MIN) && (cutoff < 0))
      || ((type & CUTOFF_MAX) && (cutoff > 0)))
      ? 1.2 * cutoff : 0.8 * cutoff;
    k_cutoff = (((type & CUTOFF_MIN) && (cutoff < 0))
      || ((type & CUTOFF_MAX) && (cutoff > 0)))
      ? 0.8 * cutoff : 1.2 * cutoff;
  }
  tee_printf(od, "TYPE      AFFECTED\n");
  for (i = 0; i < od->field_num; ++i) {
    cut_x_value_count = 0;
    if (get_field_attr(od, i, OPERATE_BIT)) {
      for (j = 0; j < od->object_num; ++j) {
        if (get_object_attr(od, j, ACTIVE_BIT)) {
          for (k = 0; k < od->x_vars; ++k) {
            result = get_x_value(od, i, j, k, &value, 0);
            if (result) {
              return result;
            }
            if (MISSING(value)) {
              continue;
            }
            if (type & CUTOFF_MIN) {
              od->mel.x_data[i].min_cutoff = cutoff;
              sign = -1.0;
              first_cond = (value < j_cutoff);
              second_cond = ((value > j_cutoff) && (value < k_cutoff));
            }
            else {
              od->mel.x_data[i].max_cutoff = cutoff;
              sign = 1.0;
              first_cond = (value > j_cutoff);
              second_cond = ((value < j_cutoff) && (value > k_cutoff));
            }
            if (first_cond) {
              result = set_x_value(od, i, j, k,
                cutoff * (1.0 + sign * cutoff / fabs(cutoff) * 0.1));
              if (result) {
                return result;
              }
              ++cut_x_value_count;
            }
            else if ((type & CUTOFF_QUADRATIC) && second_cond) {
              result = set_x_value(od, i, j, k,
                0.5 / (k_cutoff - j_cutoff)
                * (square(value) - 2.0 * j_cutoff * value
                + square(k_cutoff)));
              if (result) {
                return result;
              }
              ++cut_x_value_count;
            }
          }
        }
      }
      value = (double)(od->mel.x_data[i].active_x_vars *
        od->active_object_num);
      tee_printf(od, "x%-9d%-9d (%3d%%)\n", i + 1, cut_x_value_count,
        (value ? (int)safe_rint((double)cut_x_value_count /
        value * 100.0) : 0));
    }
  }
  result = calc_active_vars(od, FULL_MODEL);

  return result;
}
