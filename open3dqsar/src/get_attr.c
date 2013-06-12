/*

get_attr.c

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


unsigned short get_field_attr(O3Data *od,
  int field_num, unsigned short attr)
{
  return (unsigned short)((od->mel.field_attr
    [field_num]) & attr);
}


unsigned short get_object_attr(O3Data *od,
  int object_num, unsigned short attr)
{
  return (unsigned short)((od->mel.object_attr
    [object_num]) & attr);
}


unsigned short get_x_var_attr(O3Data *od,
  int field_num, int x_var, unsigned short attr)
{
  return (unsigned short)((od->mel.x_var_attr
    [field_num][x_var]) & attr);
}


unsigned short get_y_var_attr(O3Data *od,
  int y_var, unsigned short attr)
{
  return (unsigned short)
    ((od->mel.y_var_attr[y_var]) & attr);
}


void get_attr_struct_ave(O3Data *od, int y_var, unsigned short attr, int *attr_struct_num, double *attr_value_ave)
{
  int object_num = 0;
  int struct_num;
  int conf_num;
  int i = 0;
  int y;
  int n_conf = 0;
  double sumweight = 0.0;
  
  
  *attr_struct_num = 0;
  if (attr_value_ave) {
    *attr_value_ave = 0.0;
  }
  while (object_num < od->object_num) {
    struct_num = od->al.mol_info[object_num]->struct_num;
    conf_num = ((od->valid & COSMOTHERM_BIT)
      ? od->al.cosmo_list[struct_num]->n_conf[BOUND] : 1);
    for (n_conf = 0, sumweight = 0.0, y = 0;
      n_conf < conf_num; ++n_conf, ++object_num) {
      if (get_object_attr(od, object_num, attr)) {
        sumweight += od->mel.object_weight[object_num];
      }
    }
    if (sumweight > 0.0) {
      if (attr_value_ave) {
        *attr_value_ave += get_y_value(od, object_num - conf_num, y_var, WEIGHT_BIT);
      }
      ++(*attr_struct_num);
    }
    i += y;
  }
  if (attr_value_ave) {
    *attr_value_ave /= (double)(*attr_struct_num);
  }
}
