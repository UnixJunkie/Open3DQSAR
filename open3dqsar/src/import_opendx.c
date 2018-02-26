/*

import_opendx.c

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


int read_dx_header(O3Data *od, FileDescriptor *inp_fd, int object_num)
{
  char line[BUF_LEN];
  char *data_follows;
  char *gridpositions;
  char *origin;
  char *delta;
  float d[3];
  int data_follows_found;
  int delta_found;
  int origin_found;
  int gridpositions_found;
  int n_spaces = 0;
  int x_vars = 0;
  int i;
  
  memset(line, 0, BUF_LEN);
  memset(&(od->newgrid), 0, sizeof(GridInfo));
  data_follows_found = 0;
  delta_found = 0;
  origin_found = 0;
  gridpositions_found = 0;
  while ((!(data_follows_found && delta_found && origin_found
    && gridpositions_found)) && fgets(line, BUF_LEN, inp_fd->handle)) {
    /*
    skip comments
    */
    if (line[0] == '#') {
      continue;
    }
    data_follows = strstr(line, "data follows");
    if (data_follows) {
      data_follows_found = 1;
      n_spaces = 0;
      while ((data_follows != line) && (n_spaces < 2)) {
        --data_follows;
        if (isspace(*data_follows)) {
          ++n_spaces;
          while (((data_follows - 1) != line) && isspace(*(data_follows - 1))) {
            --data_follows;
          }
        }
      }
      sscanf(data_follows, "%d", &(od->newgrid.x_vars));
    }
    gridpositions = strstr(line, "gridpositions counts");
    if (gridpositions) {
      gridpositions_found = 1;
      sscanf(gridpositions, "%*s %*s %d %d %d", &(od->newgrid.nodes[0]),
        &(od->newgrid.nodes[1]), &(od->newgrid.nodes[2]));
    }
    origin = strstr(line, "origin");
    if (origin) {
      origin_found = 1;
      sscanf(origin, "%*s %f %f %f", &(od->newgrid.start_coord[0]),
        &(od->newgrid.start_coord[1]), &(od->newgrid.start_coord[2]));
    }
    delta = strstr(line, "delta");
    if (delta) {
      delta_found = 1;
      sscanf(delta, "%*s %f %f %f", &d[0], &d[1], &d[2]);
      for (i = 0; i < 3; ++i) {
        od->newgrid.step[i] += d[i];
      }
    }
  }
  if (!(data_follows_found && delta_found
    && origin_found && gridpositions_found)) {
    return BAD_DX_HEADER;
  }
  if (od->newgrid.x_vars != (od->newgrid.nodes[0]
    * od->newgrid.nodes[1] * od->newgrid.nodes[2])) {
    return BAD_DX_HEADER;
  }
  for (i = 0; i < 3; ++i) {
    od->newgrid.end_coord[i] = od->newgrid.start_coord[i]
      + (od->newgrid.nodes[i] - 1) * od->newgrid.step[i];
  }
  if (!(od->grid.nodes[0])) {
    od->newgrid.object_num = od->grid.object_num;
    memcpy(&(od->grid), &(od->newgrid), sizeof(GridInfo));
    od->x_vars = od->grid.x_vars;
    od->object_num = od->grid.object_num;
  }
  if (!match_grids(od)) {
    return GRID_NOT_MATCHING;
  }
  
  return 0;
}


int import_opendx(O3Data *od, char *regex_name)
{
  char line[BUF_LEN];
  char parsed_line[BUF_LEN];
  char alloc_field = 0;
  int i;
  int n_grid_x_vars;
  int n_read;
  int object_num;
  int result;
  double value[3];
  double node;
  VarCoord varcoord;
  FileDescriptor inp_fd;
    

  memset(&inp_fd, 0, sizeof(FileDescriptor));
  /*
  first, check that all files have the correct header
  and number of data points
  */
  for (object_num = 0; object_num < od->grid.object_num; ++object_num) {
    sprintf(inp_fd.name, regex_name, object_num + 1);
    inp_fd.handle = fopen(inp_fd.name, "rb");
    if (!(inp_fd.handle)) {
      return NOT_ENOUGH_OBJECTS;
    }
    result = read_dx_header(od, &inp_fd, object_num);
    if (result) {
      fclose(inp_fd.handle);
      od->newgrid.object_num = object_num;
      return result;
    }
    if (!alloc_field) {
      /*
      allocate one more field
      */
      result = alloc_x_var_array(od, 1);
      if (result) {
        fclose(inp_fd.handle);
        return result;
      }
      alloc_field = 1;
    }
    n_grid_x_vars = 0;
    n_read = 3;
    memset(&varcoord, 0, sizeof(VarCoord));
    while ((n_grid_x_vars < od->x_vars) && fgets(line, BUF_LEN, inp_fd.handle)) {
      n_grid_x_vars += 3;
      if (n_grid_x_vars > od->x_vars) {
        n_read -= (n_grid_x_vars - od->x_vars);
      }
      sscanf(line, "%lf %lf %lf", &value[0], &value[1], &value[2]);
      for (i = 0; i < n_read; ++i) {
        result = set_x_value_xyz(od, od->field_num - 1,
          object_num, &varcoord, value[i]);
        if (result) {
          fclose(inp_fd.handle);
          return result;
        }
        ++varcoord.node[2];
        if (varcoord.node[2] == od->grid.nodes[2]) {
          varcoord.node[2] = 0;
          ++varcoord.node[1];
          if (varcoord.node[1] == od->grid.nodes[1]) {
            varcoord.node[1] = 0;
            ++varcoord.node[0];
          }
        }
      }
    }
    fclose(inp_fd.handle);
  }
  update_field_object_attr(od, VERBOSE_BIT);

  result = calc_active_vars(od, FULL_MODEL);

  return 0;
}
