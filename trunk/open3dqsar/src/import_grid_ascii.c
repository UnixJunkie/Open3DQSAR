/*

import_grid_ascii.c

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


double parse_grid_ascii_line(char *line, char *parsed_line, VarCoord *varcoord)
{
  int i = 0;
  int j = 0;
  double value;
  
  
  while (line[i]) {
    if (strchr("0123456789-+eE. \t", (int)line[i])) {
      parsed_line[j] = line[i];
      ++j;
    }
    ++i;
  }
  parsed_line[j] = '\0';
  sscanf(parsed_line, "%lf %lf %lf %lf", &(varcoord->cart[0]),
    &(varcoord->cart[1]), &(varcoord->cart[2]), &value);
    
  return value;
}


int import_grid_ascii(O3Data *od, char *regex_name)
{
  char line[BUF_LEN];
  char parsed_line[BUF_LEN];
  int i;
  int n_grid_x_vars;
  int object_num;
  int result;
  double value;
  double node;
  VarCoord varcoord;
  FileDescriptor inp_fd;
    

  memset(&inp_fd, 0, sizeof(FileDescriptor));
  /*
  allocate one more field
  */
  result = alloc_x_var_array(od, 1);
  if (result) {
    return result;
  }
  /*
  first, check out that all files contain od->x_vars lines
  */
  for (object_num = 0; object_num < od->object_num; ++object_num) {
    sprintf(inp_fd.name, regex_name, object_num + 1);
    inp_fd.handle = fopen(inp_fd.name, "rb");
    if (!(inp_fd.handle)) {
      return NOT_ENOUGH_OBJECTS;
    }
    memcpy(&(od->newgrid), &(od->grid), sizeof(GridInfo));
    n_grid_x_vars = 0;
    while (fgets(line, BUF_LEN, inp_fd.handle)) {
      ++n_grid_x_vars;
    }
    if (n_grid_x_vars != od->x_vars) {
      od->newgrid.x_vars = n_grid_x_vars;
      od->newgrid.object_num = object_num;
      return GRID_NOT_MATCHING;
    }
    fclose(inp_fd.handle);
  }
  for (object_num = 0; object_num < od->object_num; ++object_num) {
    sprintf(inp_fd.name, regex_name, object_num + 1);
    inp_fd.handle = fopen(inp_fd.name, "rb");
    if (!(inp_fd.handle)) {
      return NOT_ENOUGH_OBJECTS;
    }
    for (n_grid_x_vars = 0; n_grid_x_vars < od->x_vars; ++n_grid_x_vars) {
      /*
      read data points one by one from each object file
      */
      if (fgets(line, BUF_LEN, inp_fd.handle)) {
        line[BUF_LEN - 1] = '\0';
        value = parse_grid_ascii_line(line, parsed_line, &varcoord);
        result = 0;
        for (i = 0; i < 3; ++i) {
          if (((safe_rint(varcoord.cart[i] * 1.0e03) / 1.0e03)
            < (safe_rint((double)(od->newgrid.start_coord[i]) * 1.0e03) / 1.0e03))
            || ((safe_rint(varcoord.cart[i] * 1.0e03) / 1.0e03)
            > (safe_rint((double)(od->newgrid.end_coord[i]) * 1.0e03) / 1.0e03))) {
            od->newgrid.object_num = object_num;
            result = GRID_NOT_MATCHING_OUT_OF_BOUNDS;
          }
          node = (varcoord.cart[i] - (double)(od->newgrid.start_coord[i]))
            / (double)(od->newgrid.step[i]);
          if ((node - safe_rint(node)) > 1.0e-03) {
            od->newgrid.object_num = object_num;
            result = GRID_NOT_MATCHING_OFF_CENTER;
          }
          if (result) {
            od->newgrid.start_coord[i] = (float)(varcoord.cart[i]);
            od->newgrid.x_vars = n_grid_x_vars + 1;
          }
          varcoord.node[i] = (int)safe_rint(node);
        }
        if (result) {
          fclose(inp_fd.handle);
          return result;
        }
        /*
        if the current value is a MISSING value,
        then replace it with MISSING_VALUE
        */
        result = set_x_value_xyz(od, od->field_num - 1, object_num, &varcoord,
          (strstr(line, COMFA_MISSING_VALUE) ? MISSING_VALUE : value));
        if (result) {
          fclose(inp_fd.handle);
          return result;
        }
      }
      else {
        od->newgrid.x_vars = n_grid_x_vars;
        od->newgrid.object_num = object_num;
        return CANNOT_READ_GRID_DATA;
      }
    }
    fclose(inp_fd.handle);
  }
  update_field_object_attr(od, VERBOSE_BIT);

  result = calc_active_vars(od, FULL_MODEL);

  return 0;
}
