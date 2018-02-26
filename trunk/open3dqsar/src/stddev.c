/*

stddev.c

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


int stddev_x_var(O3Data *od, int field_num)
{
  int i;
  int j;
  int flag;
  int result;
  int struct_num;
  int conf_num;
  int n_conf;
  double value;
  double weighted_value;
  double ave;
  double n;
  double sum1 = 0.0;
  double temp;
  double sumweight;
  double sumweight1 = 0.0;
  double stddev;
  
  
  for (i = 0; i < od->x_vars; ++i) {
    ave = 0.0;
    for (j = 0, sumweight = 0.0; j < od->object_num; ++j) {
      if (!get_object_attr(od, j, ACTIVE_BIT)) {
        continue;
      }
      result = get_x_value(od, field_num, j, i, &value, 0);
      if (result) {
        return result;
      }
      if (!MISSING(value)) {
        result = get_x_value(od, field_num, j, i,
          &value, CUTOFF_BIT);
        if (result) {
          return result;
        }
        ave += (value * od->mel.object_weight[j]);
        sumweight += od->mel.object_weight[j];
      }
    }
    if (sumweight > 0.0) {
      ave /= sumweight;
    }
    set_x_var_buf(od, field_num, i, AVE_BUF, ave);
    flag = 0;
    n = 0.0;
    j = 0;
    while (j < od->object_num) {
      struct_num = od->al.mol_info[j]->struct_num;
      conf_num = 1;
      for (n_conf = 0, sumweight = 0.0, weighted_value = 0.0;
        n_conf < conf_num; ++n_conf, ++j) {
        if (!get_object_attr(od, j, ACTIVE_BIT)) {
          continue;
        }
        result = get_x_value(od, field_num, j, i, &value, 0);
        if (result) {
          return result;
        }
        if (!MISSING(value)) {
          result = get_x_value(od, field_num, j, i,
            &value, CUTOFF_BIT);
          if (result) {
            return result;
          }
          sumweight += od->mel.object_weight[j];
          weighted_value += od->mel.object_weight[j] * value;
        }
      }
      if (sumweight > 0.0) {
        weighted_value /= sumweight;
        n += 1.0;
        if (!flag) {
          flag = 1;
          ave = weighted_value;
          sum1 = 0.0;
          sumweight1 = sumweight;
        }
        else {  
          temp = sumweight + sumweight1;
          sum1 += (sumweight1 * square(weighted_value - ave) / temp);
          ave += (weighted_value - ave) / temp;
          sumweight1 = temp;
        }
      }
    }
    stddev = 0.0;
    if (n > 1.0) {
      stddev = sqrt(sum1 * n / ((n - 1.0) * sumweight1));
    }
    set_x_var_buf(od, field_num, i, STDDEV_BUF, stddev);
  }
  
  return 0;
}


void stddev_y_var(O3Data *od)
{
  int i;
  int j;
  int flag;
  int struct_num;
  int conf_num;
  int n_conf;
  double weighted_value;
  double ave;
  double n;
  double sum1 = 0.0;
  double temp;
  double sumweight;
  double sumweight1 = 0.0;
  double stddev;
  
  
  for (i = 0; i < od->y_vars; ++i) {
    ave = 0.0;
    for (j = 0, sumweight = 0.0; j < od->grid.object_num; ++j) {
      if (!get_object_attr(od, j, ACTIVE_BIT)) {
        continue;
      }
      ave += (get_y_value(od, j, i, 0)
        * od->mel.object_weight[j]);
      sumweight += od->mel.object_weight[j];
    }
    if (sumweight > 0.0) {
      ave /= sumweight;
    }
    set_y_var_buf(od, i, AVE_BUF, ave);
    flag = 0;
    n = 0.0;
    j = 0;
    while (j < od->grid.object_num) {
      struct_num = od->al.mol_info[j]->struct_num;
      conf_num = 1;
      for (n_conf = 0, sumweight = 0.0, weighted_value = 0.0;
        n_conf < conf_num; ++n_conf, ++j) {
        if (!get_object_attr(od, j, ACTIVE_BIT)) {
          continue;
        }
        sumweight += od->mel.object_weight[j];
        weighted_value += get_y_value(od, j, i, 0)
          * od->mel.object_weight[j];
      }
      if (sumweight > 0.0) {
        weighted_value /= sumweight;
        n += 1.0;
        if (!flag) {
          flag = 1;
          ave = weighted_value;
          sum1 = 0.0;
          sumweight1 = sumweight;
        }
        else {  
          temp = sumweight + sumweight1;
          sum1 += (sumweight1 * square(weighted_value - ave) / temp);
          ave += (weighted_value - ave) / temp;
          sumweight1 = temp;
        }
      }
    }
    stddev = 0.0;
    if (n > 1.0) {
      stddev = sqrt(sum1 * n / ((n - 1.0) * sumweight1));
    }
    set_y_var_buf(od, i, STDDEV_BUF, stddev);
  }
}
