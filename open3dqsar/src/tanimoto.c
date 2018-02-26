/*

tanimoto.c

is part of

Open3DQSAR
----------

An open-source software aimed at high-throughput
chemometric analysis of molecular interaction fields

Copyright (C) 2009-2018 Paolo Tosco, Thomas Balle

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


int tanimoto(O3Data *od, int ref_struct)
{
  int object_num = 0;
  int ref_object_num = 0;
  int struct_num;
  int conf_num;
  int ref_conf_num;
  int n_conf;
  int operate_field_num;
  int i;
  int j;
  int k;
  int n;
  int result;
  double value;
  double ref_value;
  double sumweight;
  double sum_xy = 0.0;
  double sum_x2 = 0.0;
  double sum_y2 = 0.0;
  double tanimoto = 0.0;

  
  for (i = 0, operate_field_num = 0; i < od->field_num; ++i) {
    if (!get_field_attr(od, i, OPERATE_BIT)) {
      continue;
    }
    ++operate_field_num;
  }
  tee_printf(od, "Tanimoto similarity values against reference structure %d for field%s ",
    ref_struct + 1, (operate_field_num > 1) ? "s" : "");
  for (i = 0, j = 0; i < od->field_num; ++i) {
    if (!get_field_attr(od, i, OPERATE_BIT)) {
      continue;
    }
    tee_printf(od, "%d%s", i + 1, (j == (operate_field_num - 1)) ? "\n" : ", ");
    ++j;
  }
  tee_printf(od, "\n-----");
  for (i = 0; i < operate_field_num; ++i) {
    tee_printf(od, "------------");
  }
  tee_printf(od, "\n%5s", "Str");
  for (i = 0; i < od->field_num; ++i) {
    if (!get_field_attr(od, i, OPERATE_BIT)) {
      continue;
    }
    tee_printf(od, "%12d", i + 1);
  }
  tee_printf(od, "\n-----");
  for (i = 0; i < operate_field_num; ++i) {
    tee_printf(od, "------------");
  }
  ref_object_num = 0;
  while ((ref_object_num < od->object_num)
    && (od->al.mol_info[ref_object_num]->struct_num < ref_struct)) {
    ++ref_object_num;
  }
  ref_conf_num = 1;
  struct_num = od->al.mol_info[object_num]->struct_num;
  while (object_num < od->object_num) {
    struct_num = od->al.mol_info[object_num]->struct_num;
    conf_num = 1;
    for (n_conf = 0, sumweight = 0.0, n = 0; n_conf < conf_num; ++n_conf) {
      if (get_object_attr(od, object_num + n_conf, OPERATE_BIT)) {
        sumweight += od->mel.object_weight[object_num + n_conf];
        ++n;
      }
    }
    if (!n) {
      object_num += n_conf;
      continue;
    }
    tee_printf(od, "\n%5d", struct_num + 1);
    for (i = 0; i < od->field_num; ++i) {
      if (!get_field_attr(od, i, OPERATE_BIT)) {
        continue;
      }
      sum_xy = 0.0;
      sum_x2 = 0.0;
      sum_y2 = 0.0;
      for (k = 0; k < od->x_vars; ++k) {
        for (n_conf = 0, sumweight = 0.0; n_conf < conf_num; ++n_conf) {
          if (!get_object_attr(od, object_num + n_conf, OPERATE_BIT)) {
            continue;
          }
          sumweight += od->mel.object_weight[object_num + n_conf];
          result = get_x_value(od, i, object_num + n_conf, k, &value, 0);
          if (result) {
            return result;
          }
          if (MISSING(value)) {
            value = 0.0;
          }
          else {
            result = get_x_value(od, i, object_num + n_conf, k, &value, CUTOFF_BIT);
            if (result) {
              return result;
            }
          }
          value *= od->mel.object_weight[object_num + n_conf];
        }
        value /= sumweight;
        for (n_conf = 0, sumweight = 0.0; n_conf < ref_conf_num; ++n_conf) {
          if (!get_object_attr(od, ref_object_num + n_conf, OPERATE_BIT)) {
            continue;
          }
          sumweight += od->mel.object_weight[ref_object_num + n_conf];
          result = get_x_value(od, i, ref_object_num + n_conf, k, &ref_value, 0);
          if (result) {
            return result;
          }
          if (MISSING(ref_value)) {
            ref_value = 0.0;
          }
          else {;
            result = get_x_value(od, i, ref_object_num + n_conf, k, &ref_value, CUTOFF_BIT);
            if (result) {
              return result;
            }
          }
          ref_value *= od->mel.object_weight[ref_object_num + n_conf];
        }
        ref_value /= sumweight;
        sum_xy += (value * ref_value);
        sum_x2 += square(value);
        sum_y2 += square(ref_value);
      }
      tanimoto = sum_xy / (sum_x2 + sum_y2 - sum_xy);
      tee_printf(od, "%12.4lf", tanimoto);
    }
    object_num += n_conf;
  }
  tee_printf(od, "\n\n");
  
  return 0;
}
