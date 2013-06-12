/*

import_gridkont.c

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


int get_gridkont_data_points(O3Data *od, int new_model,
  int replace_object_name, int endianness_switch, int dry_run)
{
  char *header_start = NULL;
  char header[BUF_LEN];
  int data_len_start;
  int data_len_end;
  int actual_len;
  int i;
  int n;
  int x = 0;
  int y = 0;
  int field_num = 0;
  int object_num = 0;
  int eof;
  float float_buf;
  FILE *kont_in;


  kont_in = od->file[BINARY_IN]->handle;
  eof = 0;
  data_len_start = 68;
  while ((!eof) && (data_len_start == 68)) {
    /*
    read the header (64 bytes)
    */
    memset(header, 0, BUF_LEN);
    actual_len = fread(header, 1, 64, kont_in);
    if (actual_len != 64) {
      return PREMATURE_DAT_EOF;
    }
    /*
    read the object number
    */
    actual_len = fread(&i, sizeof(int), 1, kont_in);
    if (actual_len != 1) {
      return PREMATURE_DAT_EOF;
    }
    fix_endianness(&i, sizeof(int), 1, endianness_switch);
    object_num = i - 1;
    if (!object_num) {
      /*
      if this is the first object, then
      we have one more field; if this is not the first field
      to be found, then the previous object number
      is the total object number
      */
      if (field_num) {
        if (!(od->newgrid.object_num)) {
          od->newgrid.object_num = object_num;
          od->newgrid.struct_num = object_num;
        }
        else {
          /*
          all fields should have the same number of objects
          */
          if (object_num != od->newgrid.object_num) {
            return PREMATURE_DAT_EOF;
          }
        }
      }
      ++field_num;
      if (dry_run) {
        od->newgrid.field_num = field_num;
      }
    }
    if (!dry_run) {
      if (od->save_ram) {
        if (fseek(od->file
          [TEMP_FIELD_DATA + od->field_num - od->newgrid.field_num + field_num - 1]->handle,
          od->object_pagesize * object_num, SEEK_SET)) {
          return CANNOT_WRITE_TEMP_FILE;
        }
      }
      if (!(header_start = strstr(header, DAT_HEADER))) {
        return PREMATURE_DAT_EOF;
      }
      *header_start = '\0';
      i = 0;
      while (header[i] && (!isspace((int)header[i]))) {
        if (strchr(" \t,:;()[]{}\"!$%&/\\|*+=?'^<>", (int)header[i])) {
          header[i] = '_';
        }
        ++i;
      }
      header[i] = '\0';
      if (new_model) {
        od->al.mol_info[object_num]->object_id = object_num + 1;
        tee_printf(od, "%5d%5d    %-36s%10d\n", object_num + 1, object_num + 1,
          od->al.mol_info[object_num]->object_name,
          od->field_num - od->newgrid.field_num + field_num);
      }
      else {
        tee_printf(od, "%5d%5d    %-36s%-36s%10d\n", object_num + 1, object_num + 1,
          header, od->al.mol_info[object_num]->object_name,
          od->field_num - od->newgrid.field_num + field_num);
      }
      if (replace_object_name) {
        strcpy(od->al.mol_info[object_num]->object_name, header);
      }
    }
    actual_len = fread(&data_len_end, sizeof(int), 1, kont_in);
    if (actual_len != 1) {
      return PREMATURE_DAT_EOF;
    }
    fix_endianness(&data_len_end, sizeof(int), 1, endianness_switch);
    if (data_len_start != data_len_end) {
      return PREMATURE_DAT_EOF;
    }
    actual_len = fread(&data_len_start, sizeof(int), 1, kont_in);
    if (actual_len != 1) {
      return PREMATURE_DAT_EOF;
    }
    fix_endianness(&data_len_start, sizeof(int), 1, endianness_switch);
    /*
    now there should be x_vars data points
    */
    n = 0;
    while (data_len_start == sizeof(float)) {
      if (n == od->newgrid.x_vars) {
        return PREMATURE_DAT_EOF;
      }
      /*
      read the energy value
      */
      actual_len = fread(&float_buf, sizeof(float), 1, kont_in);
      if (actual_len != 1) {
        return PREMATURE_DAT_EOF;
      }
      if (!dry_run) {
        /*
        store the energy value
        */
        fix_endianness(&float_buf, sizeof(float), 1, endianness_switch);
        od->mel.float_xy_mat[y * od->grid.nodes[0] + x] = float_buf;
        ++y;
        if (y == od->grid.nodes[1]) {
          y = 0;
          ++x;
          if (x == od->grid.nodes[0]) {
            x = 0;
            /*
            store the temporary matrix
            */
            if (od->save_ram) {
              actual_len = fwrite(od->mel.float_xy_mat, sizeof(float),
                od->grid.nodes[0] * od->grid.nodes[1], od->file
                [TEMP_FIELD_DATA + od->field_num - od->newgrid.field_num + field_num - 1]->handle);
              if (actual_len != (od->grid.nodes[0] * od->grid.nodes[1])) {
                return CANNOT_WRITE_TEMP_FILE;
              }
            }
            else {
              memcpy(&(od->mel.x_var_array
                [od->field_num - od->newgrid.field_num + field_num - 1]
                [object_num][n]), od->mel.float_xy_mat,
                sizeof(float) * od->grid.nodes[0] * od->grid.nodes[1]);
            }
            n += (od->grid.nodes[0] * od->grid.nodes[1]);
          }
        }
      }
      else {
        ++n;
      }
      actual_len = fread(&data_len_end, sizeof(int), 1, kont_in);
      if (actual_len != 1) {
        return PREMATURE_DAT_EOF;
      }
      fix_endianness(&data_len_end, sizeof(int), 1, endianness_switch);
      if (data_len_start != data_len_end) {
        return PREMATURE_DAT_EOF;
      }
      actual_len = fread(&data_len_start, sizeof(int), 1, kont_in);
      if (actual_len != 1) {
        /*
        have we reached EOF?
        */
        if (field_num > 1) {
          if ((object_num + 1) == od->newgrid.object_num) {
            eof = 1;
          }
          else {
            return PREMATURE_DAT_EOF;
          }
        }
        else {
          if (!(od->newgrid.object_num)) {
            od->newgrid.object_num = object_num + 1;
            od->newgrid.struct_num = object_num + 1;
          }
          eof = 1;
        }
      }
      fix_endianness(&data_len_start, sizeof(int), 1, endianness_switch);
    }
    if (n != od->newgrid.x_vars) {
      return PREMATURE_DAT_EOF;
    }
  }
  if (!eof) {
    return PREMATURE_DAT_EOF;
  }

  return 0;
}


int import_gridkont(O3Data *od, int replace_object_name)
{
  char dummy;
  int endianness_switch;
  int new_model;
  int n;
  int i;
  int actual_len;
  int result;
  int data_len_start;
  int data_len_end;
  int pos1;
  float coord[3];
  FILE *kont_in;
  

  kont_in = od->file[BINARY_IN]->handle;
  memset(&(od->newgrid), 0, sizeof(GridInfo));
  od->newgrid.x_vars = 1;
  /*
  desume endianness
  from the first int at the very beginning
  of the header
  also take into consideration
  the endianness of the machine we are running on
  */
  actual_len = fread(&data_len_start, sizeof(int), 1, kont_in);
  if (actual_len != 1) {
    return PREMATURE_DAT_EOF;
  }
  fix_endianness(&data_len_start, sizeof(int), 1, machine_type());
  
  if (data_len_start <= 0xFFFF) {
    endianness_switch = machine_type();
  }
  else {
    endianness_switch = 1 - machine_type();
  }
  fix_endianness(&data_len_start, sizeof(int), 1, endianness_switch);
  /*
  start reading data from .kont file
  */
  while (data_len_start == 17) {
    /*
    skip a byte whose meaning is unclear
    */
    actual_len = fread(&dummy, 1, 1, kont_in);
    if (actual_len != 1) {
      return PREMATURE_DAT_EOF;
    }
    /*
    read the grid point number
    */
    actual_len = fread(&n, sizeof(int), 1, kont_in);
    if (actual_len != 1) {
      return PREMATURE_DAT_EOF;
    }
    fix_endianness(&n, sizeof(int), 1, endianness_switch);
    /*
    read the grid point coordinates
    */
    actual_len = fread(coord, sizeof(float), 3, kont_in);
    if (actual_len != 3) {
      return PREMATURE_DAT_EOF;
    }
    fix_endianness(coord, sizeof(float), 3, endianness_switch);
    /*
    if this is the first grid point,
    then these are the starting grid coordinates
    */
    if (n == 1) {
      memcpy(od->newgrid.start_coord, coord, 3 * sizeof(float));
      memcpy(od->newgrid.end_coord, coord, 3 * sizeof(float));
    }
    else {
      for (i = 0; i < 3; ++i) {
        if (coord[i] > od->newgrid.end_coord[i]) {
          od->newgrid.end_coord[i] = coord[i];
          ++(od->newgrid.nodes[i]);
        }
        if (coord[i] < od->newgrid.start_coord[i]) {
          od->newgrid.start_coord[i] = coord[i];
          ++(od->newgrid.nodes[i]);
        }
      }
    }
    actual_len = fread(&data_len_end, sizeof(int), 1, kont_in);
    if (actual_len != 1) {
      return PREMATURE_DAT_EOF;
    }
    fix_endianness(&data_len_end, sizeof(int), 1, endianness_switch);
    if (data_len_start != data_len_end) {
      return PREMATURE_DAT_EOF;
    }
    actual_len = fread(&data_len_start, sizeof(int), 1, kont_in);
    if (actual_len != 1) {
      return PREMATURE_DAT_EOF;
    }
    fix_endianness(&data_len_start, sizeof(int), 1, endianness_switch);
  }
  for (i = 0; i < 3; ++i) {
    ++(od->newgrid.nodes[i]);
    od->newgrid.step[i] =
      (od->newgrid.end_coord[i] - od->newgrid.start_coord[i])
      / (float)(od->newgrid.nodes[i] - 1);
    od->newgrid.x_vars *= od->newgrid.nodes[i];
  }
  if ((n != od->newgrid.x_vars) || (data_len_start != 68)) {
    return PREMATURE_DAT_EOF;
  }
  /*
  record where we are because we will need to rewind
  to this point
  */
  pos1 = ftell(kont_in);
  result = get_gridkont_data_points(od, 0, 0, endianness_switch, DRY_RUN);
  if (result) {
    return result;
  }
  tee_printf(od, "Number of fields:    %d\n", od->newgrid.field_num);
  tee_printf(od, "Number of objects:   %d\n\n", od->newgrid.object_num);
  if (od->grid.object_num != od->newgrid.object_num) {
    return OBJECTS_NOT_MATCHING;
  }
  new_model = !(od->grid.nodes[0]);
  if (new_model) {
    memcpy(&(od->grid), &(od->newgrid), sizeof(GridInfo));
    od->x_vars = od->grid.x_vars;
    od->object_num = od->newgrid.object_num;
    od->active_object_num = od->object_num;
    od->object_pagesize = (od->x_vars * sizeof(float) / od->mmap_pagesize +
      (od->x_vars * sizeof(float) % od->mmap_pagesize ? 1 : 0))
      * od->mmap_pagesize;
  }
  else {
    if (!match_grids(od)) {
      return GRID_NOT_MATCHING;
    }
  }
  print_grid_coordinates(od, &(od->newgrid));
  if (alloc_x_var_array(od, od->newgrid.field_num)) {
    return OUT_OF_MEMORY;
  }
  if (new_model) {
    tee_printf(od, "------------------------------------------------------------\n"
      "%5s%5s    %-36s%10s\n"
      "------------------------------------------------------------\n",
      "N", "ID", "Name", "Field");
  }
  else {
    tee_printf(od, "-------------------------------------------------------------"
      "-----------------------------------\n");
    if (replace_object_name) {
      tee_printf(od, "%5s%5s    %-36s%-36s%10s\n", "N", "ID",
        "New name (replaces the old one)", "Old name", "Field");
    }
    else {
      tee_printf(od, "%5s%5s    %-36s%-36s%10s\n", "N", "ID",
        "New name", "Old name (retained)", "Field");
    }
    tee_printf(od, "--------------------------------------------------------------"
      "----------------------------------\n");
  }
  od->mel.float_xy_mat = (float *)realloc(od->mel.float_xy_mat,
    sizeof(float) * od->grid.nodes[0] * od->grid.nodes[1]);
  if (!(od->mel.float_xy_mat)) {
    return OUT_OF_MEMORY;
  }
  /*
  go back to where field data start
  */
  if (fseek(kont_in, pos1, SEEK_SET)) {
    return PREMATURE_DAT_EOF;
  }
  result = get_gridkont_data_points(od, new_model,
    replace_object_name, endianness_switch, 0);
  if (result) {
    return result;
  }
  free(od->mel.float_xy_mat);
  od->mel.float_xy_mat = NULL;
  fflush(stdout);
  if (od->save_ram) {
    sync_field_mmap(od);
  }
  
  update_field_object_attr(od, VERBOSE_BIT);

  result = calc_active_vars(od, FULL_MODEL);

  return result;
}
