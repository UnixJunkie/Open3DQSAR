/*

import_grid_formatted_cube.c

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


int import_grid_formatted_cube(O3Data *od, int mo)
{
  char buffer[BUF_LEN];
  char *context = NULL;
  char *value;
  char *eof = NULL;
  int i;
  int j;
  int mo_cube = 0;
  int found = 0;
  int object_num;
  int field_num;
  int initial_field_num;
  int result;
  int i_data[4];
  double d_data[4];
  VarCoord varcoord;
    

  memset(buffer, 0, BUF_LEN);
  initial_field_num = od->field_num;
  rewind(od->file[TEMP_SORTED_MATCH]->handle);
  for (object_num = 0; object_num < od->object_num; ++object_num) {
    memset(&varcoord, 0, sizeof(VarCoord));
    if (!fgets(od->file[ASCII_IN]->name, BUF_LEN, od->file[TEMP_SORTED_MATCH]->handle)) {
      return CANNOT_READ_TEMP_FILE;
    }
    remove_newline(od->file[ASCII_IN]->name);
    od->file[ASCII_IN]->handle = fopen(od->file[ASCII_IN]->name, "rb");
    if (!(od->file[ASCII_IN]->handle)) {
      return NOT_ENOUGH_OBJECTS;
    }
    if (!(found = fgrep(od->file[ASCII_IN]->handle, buffer, "$CUBE"))) {
      found = fgrep(od->file[ASCII_IN]->handle, buffer, "START OF CUBE FORMAT");
    }
    if (!found) {
      rewind(od->file[ASCII_IN]->handle);
    }
    /*
    skip two lines
    */
    found = 1;
    for (i = 0; (i < 2) && found; ++i) {
      found = (fgets(buffer, BUF_LEN, od->file[ASCII_IN]->handle) ? 1 : 0);
    }
    if (!found) {
      return PREMATURE_DAT_EOF;
    }
    /*
    now read 4 blocks constituted by:
    - 1 int
    - 3 doubles (origin x, y, z)
    */
    for (i = 0; i < 4; ++i) {
      if (!fgets(buffer, BUF_LEN, od->file[ASCII_IN]->handle)) {
        return PREMATURE_DAT_EOF;
      }
      buffer[BUF_LEN - 1] = '\0';
      sscanf(buffer, "%d %lf %lf %lf", &i_data[i], &d_data[0], &d_data[1], &d_data[2]);
      if (!i) {
        /*
        if NAtoms is negative, this is a MO cube file
        */
        mo_cube = (i_data[0] < 0);
        i_data[0] = absval(i_data[0]);
        for (j = 0; j < 3; ++j) {
          od->newgrid.start_coord[j] = (float)
            (safe_rint(d_data[j] * BOHR_RADIUS * 1.0e04) / 1.0e04);
        }
      }
      else {
        od->newgrid.nodes[i - 1] = i_data[i];
        od->newgrid.step[i - 1] = (float)
          (safe_rint(d_data[i - 1] * BOHR_RADIUS * 1.0e04) / 1.0e04);
        od->newgrid.end_coord[i - 1] = (float)(safe_rint
          (((double)(od->newgrid.start_coord[i - 1]) + (double)(od->newgrid.step[i - 1])
          * (double)(i_data[i] - 1)) * 1.0e04) / 1.0e04);
      }
    }
    if (!match_grids(od)) {
      od->newgrid.x_vars = od->newgrid.nodes[0]
        *  od->newgrid.nodes[1] *  od->newgrid.nodes[2];
      od->newgrid.object_num = object_num;
      return GRID_NOT_MATCHING;
    }
    /*
    now skip natoms lines constituted by:
    - 1 int (atom number)
    - 4 doubles (charge, x_coord, y_coord, z_coord)
    */
    for (i = 0; i < i_data[0]; ++i) {
      if (!fgets(buffer, BUF_LEN - 1, od->file[ASCII_IN]->handle)) {
        return PREMATURE_DAT_EOF;
      }
    }
    /*
    now if this is a MO cube file there is the list of MOs
    */
    found = -1;
    i_data[1] = 1;
    if (mo_cube) {
      i = 0;
      j = 0;
      found = 0;
      while ((j < i_data[1]) && (eof = fgets(buffer, BUF_LEN - 1, od->file[ASCII_IN]->handle))) {
        buffer[BUF_LEN - 1] = '\0';
        value = strtok_r(buffer, " \t\n\r\0", &context);
        if (!i) {
          /*
          read how many MOs we have
          */
          sscanf(buffer, "%d", &i_data[1]);
          if ((!value) || (!i_data[1])) {
            eof = NULL;
            break;
          }
          value = strtok_r(NULL, " \t\n\r\0", &context);
        }
        /*
        now go through the MO list
        */
        while (value && (j < i_data[1])) {
          ++j;
          /*
          if the user wants all MOs, mo = 0, then
          found will remain set to 0 as it was originally;
          however we need to go through the whole list
          all the same, to skip all the lines
          which constitute the list itself
          */
          if ((!found) && (mo == j)) {
            /*
            if we found the MO we want, we save it
            */
            found = j;
          }
          value = strtok_r(NULL, " \t\n\r\0", &context);
        }
        ++i;
      }
      if ((!found) && (!mo)) {
        /*
        we want all MOs, so found = -1
        */
        found = -1;
      }
      if ((!eof) || (!found) || (found > i_data[1])) {
        od->newgrid.object_num = object_num;
        return CANNOT_FIND_MO;
      }
    }
    /*
    now there should be x_nodes * y_nodes * z_nodes (* n_mo)
    data points
    z: fastest varying coordinate
    x: slowest varying coordinate
    */
    i = 0;
    while ((i < od->grid.x_vars)
      && fgets(buffer, BUF_LEN, od->file[ASCII_IN]->handle)) {
      buffer[BUF_LEN - 1] = '\0';
      if (strstr(buffer, "$END")) {
        break;
      }
      j = 0;
      field_num = initial_field_num;
      value = strtok_r(buffer, " \t\n\r\0", &context);
      while (value && (i < od->grid.x_vars)) {
        ++j;
        /*
        either we want all MOs or this is is the MO we want
        */
        if ((found == -1) || (j == found)) {
          if ((!object_num) && (!i)) {
            /*
            allocate one more field
            */
            result = alloc_x_var_array(od, 1);
            if (result) {
              return result;
            }
          }
          sscanf(value, "%lf", &d_data[0]);
          result = (mo ? set_x_value_xyz(od, field_num,
            object_num, &varcoord, d_data[0])
            : set_x_value_xyz_unbuffered(od, field_num,
            object_num, &varcoord, d_data[0]));
          if (result) {
            return result;
          }
          ++field_num;
        }
        if (j == i_data[1]) {
          /*
          if this is the last MO, go back to the first field
          and increment x_var count
          */
          j = 0;
          field_num = initial_field_num;
          ++i;
          ++(varcoord.node[2]);
          if (varcoord.node[2] == od->grid.nodes[2]) {
            varcoord.node[2] = 0;
            ++(varcoord.node[1]);
          }
          if (varcoord.node[1] == od->grid.nodes[1]) {
            varcoord.node[1] = 0;
            ++(varcoord.node[0]);
          }
        }
        value = strtok_r(NULL, " \t\n\r\0", &context);
      }
    }
    if (od->file[ASCII_IN]->handle) {
      fclose(od->file[ASCII_IN]->handle);
      od->file[ASCII_IN]->handle = NULL;
    }
  }
  if (od->save_ram) {
    sync_field_mmap(od);
  }
  update_field_object_attr(od, VERBOSE_BIT);

  result = calc_active_vars(od, FULL_MODEL);

  return 0;
}
