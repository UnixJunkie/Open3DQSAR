/*

update_field_object_attr.c

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


void update_field_object_attr(O3Data *od, int verbose)
{
  char attribute[MAX_NAME_LEN];
  char format[MAX_NAME_LEN];
  char buffer[BUF_LEN];
  int i;
  int j;
  int active_field_num;
  int active_object_num;
  int ext_pred_object_num;
  double value;
  

  active_field_num = 0;
  memset(attribute, 0, MAX_NAME_LEN);
  memset(format, 0, MAX_NAME_LEN);
  memset(buffer, 0, BUF_LEN);
  if (verbose) {
    tee_printf(od, "\n"
      "Number of fields:      %d\n", od->field_num);
    if (od->field_num) {
      tee_printf(od,
        "\n"
        "-------------------------\n"
        "%5s    %-16s\n"
        "-------------------------\n",
        "Field", "Attribute");
    }
  }
  for (i = 0; i < od->field_num; ++i) {
    if (get_field_attr(od, i, ACTIVE_BIT)) {
      strcpy(attribute, "INCLUDED");
      ++active_field_num;
    }
    else {
      strcpy(attribute, "EXCLUDED");
    }
    if (verbose) {
      tee_printf(od, "%5d    %-16s\n", i + 1, attribute);
    }
  }

  active_object_num = 0;
  ext_pred_object_num = 0;
  if (verbose) {
    tee_printf(od,
      "\n"
      "Number of objects:   %d\n", od->grid.object_num);
    tee_printf(od,
      "\n"
      "Number of Y variables:   %d\n", od->y_vars);
    if (od->y_vars) {
      tee_printf(od,
        "\n"
        "Y variable full names:\n");
      for (i = 0; i < od->y_vars; ++i) {
        sprintf(buffer, "y%d", i + 1);
        tee_printf(od, "%4s: %s\n", buffer, od->cimal.y_var_name->me[i]);
      }
    }
    tee_printf(od,
      "\n"
      "-----------------------------------------------------------------------------------");
    for (i = 0; i < od->y_vars; ++i) {
      tee_printf(od, "------------");
    }
    tee_printf(od,
      "\n"
      "%5s%5s%5s    %-36s%-16s%12s", "", "", "", "", "", "");
    for (i = 0; i < od->y_vars; ++i) {
      sprintf(buffer, "y%d", i + 1);
      tee_printf(od, "%12s", buffer);
    }
    tee_printf(od,
      "\n"
      "%5s%5s%5s    %-36s%-16s%12s", "N", "ID", "Str", "Object name", "Attribute", "Weight");
    for (i = 0; i < od->y_vars; ++i) {
      strcpy(buffer, od->cimal.y_var_name->me[i]);
      if (strlen(buffer) > 11) {
        strcpy(&buffer[10], "*");
      }
      tee_printf(od, "%12s", buffer);
    }
    tee_printf(od,
      "\n"
      "-----------------------------------------------------------------------------------");
    for (i = 0; i < od->y_vars; ++i) {
      tee_printf(od, "------------");
    }
    tee_printf(od, "\n");
  }
  for (i = 0; i < od->grid.object_num; ++i) {
    if (get_object_attr(od, i, ACTIVE_BIT)) {
      strcpy(attribute, "TRAINING SET");
      ++active_object_num;
    }
    else if (get_object_attr(od, i, PREDICT_BIT)) {
      strcpy(attribute, "TEST SET");
      ++ext_pred_object_num;
    }
    else {
      strcpy(attribute, "EXCLUDED");
    }
    if (verbose) {
      tee_printf(od, "%5d%5d%5d    %-36s%-16s", i + 1,
        od->al.mol_info[i]->object_id,
        od->al.mol_info[i]->struct_num + 1,
        od->al.mol_info[i]->object_name, attribute);
      if (get_object_attr(od, i, ACTIVE_BIT | PREDICT_BIT)) {
        value = od->mel.object_weight[i];
        sprintf(buffer, "%lf", value);
        strcpy(format, ((strlen(buffer) > 11) ? "%12.4le" : "%12.4lf"));
        tee_printf(od, format, value);
      }
      else {
        tee_printf(od, "%12s", "-");
      }
      for (j = 0; j < od->y_vars; ++j) {
        value = get_y_value(od, i, j, 0);
        sprintf(buffer, "%lf", value);
        strcpy(format, ((strlen(buffer) > 11) ? "%12.4le" : "%12.4lf"));
        tee_printf(od, format, value);
      }
      tee_printf(od, "\n");
    }
  }
  tee_flush(od);
  od->active_field_num = active_field_num;
  od->active_object_num = active_object_num;
  od->ext_pred_object_num = ext_pred_object_num;
}
