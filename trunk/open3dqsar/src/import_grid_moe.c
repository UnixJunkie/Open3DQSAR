/*

import_grid_moe.c

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


int import_grid_moe(O3Data *od, char *regex_name)
{
  char buffer[MAX_NAME_LEN];
  char header[MAX_NAME_LEN];
  int i;
  int n;
  int title_len;
  int actual_len;
  int object_num;
  int field_num;
  int initial_field_num;
  int dimensions;
  int n_grids;
  int swap_endianness;
  int result;
  float coord;
  double value;
  double node;
  VarCoord varcoord;
  FILE *moe_grid_in = NULL;


  initial_field_num = od->field_num;
  swap_endianness = (int)fabs(machine_type() - O3Q_LITTLE_ENDIAN);
  memset(od->newgrid.start_coord, 0, 3 * sizeof(float));
  for (object_num = 0; object_num < od->object_num; ++object_num) {
    sprintf(od->file[BINARY_IN]->name, regex_name, object_num + 1);
    moe_grid_in = fopen(od->file[BINARY_IN]->name, "rb");
    if (!moe_grid_in) {
      return NOT_ENOUGH_OBJECTS;
    }
    od->file[BINARY_IN]->handle = moe_grid_in;
    sprintf(header, "%s%.2f", MOE_ID_TOKEN, MOE_FORMAT_VERSION);
    actual_len = fread(buffer, sizeof(char), 12, moe_grid_in);
    if ((actual_len != 12) || strncmp(buffer, header, 12)) {
      return PREMATURE_DAT_EOF;
    }
    /*
    skip 2 int values
    */
    if (fseek(moe_grid_in, 2 * sizeof(int), SEEK_CUR)) {
      return PREMATURE_DAT_EOF;
    }
    /*
    read title length in order to skip it
    */
    actual_len = fread(&title_len, sizeof(int), 1, moe_grid_in);
    if (actual_len != 1) {
      return PREMATURE_DAT_EOF;
    }
    fix_endianness(&title_len, sizeof(int), 1, swap_endianness);
    if (title_len) {
      if (fseek(moe_grid_in, title_len, SEEK_CUR)) {
        return PREMATURE_DAT_EOF;
      }
    }
    /*
    there must be 3 dimensions
    */
    actual_len = fread(&dimensions, sizeof(int), 1, moe_grid_in);
    fix_endianness(&dimensions, sizeof(int), 1, swap_endianness);
    if ((actual_len != 1) || (dimensions != 3)) {
      return PREMATURE_DAT_EOF;
    }
    /*
    now there should be 3 blocks constituted by:
    - 1 int or long int (n = number of ticks)
    - n floats (tick coordinates, we skip them)
    */
    for (i = 0; i < 3; ++i) {
      actual_len = fread(&(od->newgrid.nodes[i]),
        sizeof(int), 1, moe_grid_in);
      fix_endianness(&(od->newgrid.nodes[i]),
        sizeof(int), 1, swap_endianness);
      if ((actual_len != 1)
        || (od->newgrid.nodes[i] != od->grid.nodes[i])) {
        return PREMATURE_DAT_EOF;
      }
      for (n = 0; n < od->newgrid.nodes[i]; ++n) {
        actual_len = fread(&coord,
          sizeof(float), 1, moe_grid_in);
        if (actual_len != 1) {
          return PREMATURE_DAT_EOF;
        }
        fix_endianness(&(od->newgrid.nodes[i]),
          sizeof(float), 1, swap_endianness);
        value = (double)coord;
        node = (value - (double)(od->grid.start_coord[i]))
          / (double)(od->grid.step[i]);
        if ((node - safe_rint(node)) > 1.0e-03) {
          od->newgrid.start_coord[i] = value;
          od->newgrid.object_num = object_num;
          return GRID_NOT_MATCHING_OFF_CENTER;
        }
      }
    }
    actual_len = fread(&n_grids, sizeof(int), 1, moe_grid_in);
    if (actual_len != 1) {
      return PREMATURE_DAT_EOF;
    }
    fix_endianness(&n_grids, sizeof(int), 1, swap_endianness);
    for (field_num = 0; field_num < n_grids; ++field_num) {
      /*
      allocate one more field
      */
      if (!object_num) {
        result = alloc_x_var_array(od, 1);
        if (result) {
          return result;
        }
      }
      /*
      read label length in order to skip it
      */
      actual_len = fread(&n, sizeof(int), 1, moe_grid_in);
      if (actual_len != 1) {
        return PREMATURE_DAT_EOF;
      }
      fix_endianness(&n, sizeof(int), 1, swap_endianness);
      if (fseek(moe_grid_in, n, SEEK_CUR)) {
        return PREMATURE_DAT_EOF;
      }
      /*
      now there should be x_nodes * y_nodes arrays,
      each constituted by z_nodes doubles
      z: fastest varying coordinate
      x: slowest varying coordinate
      */
      for (varcoord.node[0] = 0; varcoord.node[0] < od->newgrid.nodes[0]; ++varcoord.node[0]) {
        for (varcoord.node[1] = 0; varcoord.node[1] < od->newgrid.nodes[1]; ++varcoord.node[1]) {
          for (varcoord.node[2] = 0; varcoord.node[2] < od->newgrid.nodes[2]; ++varcoord.node[2]) {
            actual_len = fread(&value, sizeof(double), 1, moe_grid_in);
            if (actual_len != 1) {
              return PREMATURE_DAT_EOF;
            }
            fix_endianness(&value, sizeof(double), 1, swap_endianness);
            result = set_x_value_xyz(od, field_num + initial_field_num,
              object_num, &varcoord, value);
            if (result) {
              return result;
            }
          }
        }
      }
    }
    if (moe_grid_in) {
      fclose(moe_grid_in);
      od->file[BINARY_IN]->handle = NULL;
    }
  }
  if (od->save_ram) {
    sync_field_mmap(od);
  }
  update_field_object_attr(od, VERBOSE_BIT);

  result = calc_active_vars(od, FULL_MODEL);

  return result;
}
