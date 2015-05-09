/*

import_grid_unformatted_cube.c

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


int import_grid_unformatted_cube(O3Data *od, int mo)
{
  int i;
  int j;
  int mo_cube = 0;
  int actual_len;
  int object_num;
  int field_num;
  int initial_field_num;
  int data_len_start;
  int data_len_end;
  int endianness_switch;
  int result;
  int n_atoms;
  int n_mo;
  int found = 0;
  int cube_word_size;
  long int long_n_atoms;
  long int long_node;
  long int long_n_mo;
  long int long_i;
  double value;
  double origin[3];
  double step[3][3];
  VarCoord varcoord;
  FILE *gaussian_cube_in = NULL;
  


  initial_field_num = od->field_num;
  rewind(od->file[TEMP_SORTED_MATCH]->handle);
  for (object_num = 0; object_num < od->object_num; ++object_num) {
    if (!fgets(od->file[BINARY_IN]->name, BUF_LEN, od->file[TEMP_SORTED_MATCH]->handle)) {
      return CANNOT_READ_TEMP_FILE;
    }
    remove_newline(od->file[BINARY_IN]->name);
    gaussian_cube_in = fopen(od->file[BINARY_IN]->name, "rb");
    if (!gaussian_cube_in) {
      return NOT_ENOUGH_OBJECTS;
    }
    od->file[BINARY_IN]->handle = gaussian_cube_in;
    /*
    desume endianness
    from the first int at the very beginning
    of the GAUSSIAN cube file
    also take into consideration
    the endianness of the machine we are running on
    */

    actual_len = fread(&data_len_start, sizeof(int), 1, gaussian_cube_in);
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
    at the beginning there are two 80-byte headers
    that we can skip
    */
    if (fseek(gaussian_cube_in, data_len_start, SEEK_CUR)) {
      return PREMATURE_DAT_EOF;
    }
    actual_len = fread(&data_len_end, sizeof(int), 1, gaussian_cube_in);
    fix_endianness(&data_len_end, sizeof(int), 1, endianness_switch);
    if ((actual_len != 1) || (data_len_start != data_len_end)) {
      return PREMATURE_DAT_EOF;
    }
    actual_len = fread(&data_len_start, sizeof(int), 1, gaussian_cube_in);
    if (actual_len != 1) {
      return PREMATURE_DAT_EOF;
    }
    fix_endianness(&data_len_start, sizeof(int), 1, endianness_switch);
    if (fseek(gaussian_cube_in, data_len_start, SEEK_CUR)) {
      return PREMATURE_DAT_EOF;
    }
    actual_len = fread(&data_len_end, sizeof(int), 1, gaussian_cube_in);
    fix_endianness(&data_len_end, sizeof(int), 1, endianness_switch);
    if ((actual_len != 1) || (data_len_start != data_len_end)) {
      return PREMATURE_DAT_EOF;
    }
    /*
    now there should be a block constituted by:
    - 1 int or long int (number of atoms)
    - 3 doubles (origin x, y, z)
    */
    actual_len = fread(&data_len_start, sizeof(int), 1, gaussian_cube_in);
    if (actual_len != 1) {
      return PREMATURE_DAT_EOF;
    }
    fix_endianness(&data_len_start, sizeof(int), 1, endianness_switch);
    data_len_start -= 3 * sizeof(double);
    if ((data_len_start != sizeof(int))
      && (data_len_start != sizeof(long int))) {
      return PREMATURE_DAT_EOF;
    }
    cube_word_size = data_len_start;
    actual_len = fread(&long_n_atoms, cube_word_size, 1, gaussian_cube_in);
    if (actual_len != 1) {
      return PREMATURE_DAT_EOF;
    }
    fix_endianness(&long_n_atoms, cube_word_size, 1, endianness_switch);
    n_atoms = (int)long_n_atoms;
    mo_cube = (n_atoms < 0);
    n_atoms = absval(n_atoms);
    for (i = 0; i < 3; ++i) {
      actual_len = fread(&origin[i], sizeof(double), 1, gaussian_cube_in);
      if (actual_len != 1) {
        return PREMATURE_DAT_EOF;
      }
      fix_endianness(&origin[i], sizeof(double), 1, endianness_switch);
      od->newgrid.start_coord[i] = (float)(safe_rint(origin[i] * BOHR_RADIUS * 1.0e04) / 1.0e04);
    }
    actual_len = fread(&data_len_end, sizeof(int), 1, gaussian_cube_in);
    if (actual_len != 1) {
      return PREMATURE_DAT_EOF;
    }
    fix_endianness(&data_len_end, sizeof(int), 1, endianness_switch);
    data_len_end -= 3 * sizeof(double);
    if (data_len_start != data_len_end) {
      return PREMATURE_DAT_EOF;
    }
    /*
    now there should be 3 blocks constituted by:
    - 1 int or long int (number of [X,Y,Z]_nodes)
    - 3 doubles (step size Xn,Yn,Zn)
    */
    for (i = 0; i < 3; ++i) {
      actual_len = fread(&data_len_start, sizeof(int), 1, gaussian_cube_in);
      if (actual_len != 1) {
        return PREMATURE_DAT_EOF;
      }
      fix_endianness(&data_len_start, sizeof(int), 1, endianness_switch);
      data_len_start -= 3 * sizeof(double);
      if (data_len_start != cube_word_size) {
        return PREMATURE_DAT_EOF;
      }
      actual_len = fread(&long_node, cube_word_size, 1, gaussian_cube_in);
      if (actual_len != 1) {
        return PREMATURE_DAT_EOF;
      }
      fix_endianness(&long_node, cube_word_size, 1, endianness_switch);
      od->newgrid.nodes[i] = (int)long_node;
      for (j = 0; j < 3; ++j) {
        actual_len = fread(&step[i][j], sizeof(double), 1, gaussian_cube_in);
        if (actual_len != 1) {
          return PREMATURE_DAT_EOF;
        }
        fix_endianness(&step[i][j], sizeof(double), 1, endianness_switch);
        step[i][j] = safe_rint(step[i][j] * BOHR_RADIUS * 1.0e04) / 1.0e04;
      }
      actual_len = fread(&data_len_end, sizeof(int), 1, gaussian_cube_in);
      if (actual_len != 1) {
        return PREMATURE_DAT_EOF;
      }
      fix_endianness(&data_len_end, sizeof(int), 1, endianness_switch);
      data_len_end -= 3 * sizeof(double);
      if (data_len_start != data_len_end) {
        return PREMATURE_DAT_EOF;
      }
    }
    for (i = 0; i < 3; ++i) {
      od->newgrid.step[i] = (float)step[i][i];
      od->newgrid.end_coord[i] = (float)
        (safe_rint(((double)(od->newgrid.start_coord[i])
        + (double)(od->newgrid.nodes[i] - 1)
        * step[i][i]) * 1.0e04) / 1.0e04);
    }
    /*
    now there should be n_atoms blocks constituted by:
    - 1 int or long int (atom number)
    - 4 doubles (charge, x_coord, y_coord, z_coord)
    */
    actual_len = fread(&data_len_start, sizeof(int), 1, gaussian_cube_in);
    fix_endianness(&data_len_start, sizeof(int), 1, endianness_switch);
    if ((actual_len != 1) || ((data_len_start != (n_atoms * (sizeof(int) + 4 * sizeof(double))))
      && (data_len_start != (n_atoms * (sizeof(long int) + 4 * sizeof(double)))))) {
      return PREMATURE_DAT_EOF;
    }
    if (fseek(gaussian_cube_in, data_len_start, SEEK_CUR)) {
      return PREMATURE_DAT_EOF;
    }
    actual_len = fread(&data_len_end, sizeof(int), 1, gaussian_cube_in);
    fix_endianness(&data_len_end, sizeof(int), 1, endianness_switch);
    if ((actual_len != 1) || (data_len_start != data_len_end)) {
      return PREMATURE_DAT_EOF;
    }
    if (!match_grids(od)) {
      return GRID_NOT_MATCHING;
    }
    /*
    now if this is a MO cube file there is the list of MOs
    */
    found = -1;
    n_mo = 1;
    if (mo_cube) {
      actual_len = fread(&data_len_start, sizeof(int), 1, gaussian_cube_in);
      if (actual_len != 1) {
        return PREMATURE_DAT_EOF;
      }
      fix_endianness(&data_len_start, sizeof(int), 1, endianness_switch);
      actual_len = fread(&long_n_mo, cube_word_size, 1, gaussian_cube_in);
      if (actual_len != 1) {
        return PREMATURE_DAT_EOF;
      }
      fix_endianness(&long_n_mo, cube_word_size, 1, endianness_switch);
      n_mo = (int)long_n_mo;
      if (data_len_start != ((n_mo + 1) * cube_word_size)) {
        return PREMATURE_DAT_EOF;
      }
      j = 0;
      found = 0;
      while (j < n_mo) {
        actual_len = fread(&long_i, cube_word_size, 1, gaussian_cube_in);
        if (actual_len != 1) {
          return PREMATURE_DAT_EOF;
        }
        ++j;
        if ((!found) && (mo == j)) {
          found = j;
        }
      }
      if ((!found) && (!mo)) {
        found = -1;
      }
      if ((!found) || (found > n_mo)) {
        od->newgrid.object_num = object_num;
        return CANNOT_FIND_MO;
      }
      actual_len = fread(&data_len_end, sizeof(int), 1, gaussian_cube_in);
      fix_endianness(&data_len_end, sizeof(int), 1, endianness_switch);
      if ((actual_len != 1) || (data_len_end != data_len_start)) {
        return PREMATURE_DAT_EOF;
      }
    }
    /*
    now there should be x_nodes * y_nodes arrays,
    each constituted by z_nodes doubles
    z: fastest varying coordinate
    x: slowest varying coordinate
    */
    for (varcoord.node[0] = 0; varcoord.node[0] < od->newgrid.nodes[0]; ++varcoord.node[0]) {
      for (varcoord.node[1] = 0; varcoord.node[1] < od->newgrid.nodes[1]; ++varcoord.node[1]) {
        actual_len = fread(&data_len_start, sizeof(int), 1, gaussian_cube_in);
        fix_endianness(&data_len_start, sizeof(int), 1, endianness_switch);
        if ((actual_len != 1) || (data_len_start != (od->newgrid.nodes[2] * n_mo * sizeof(double)))) {
          return PREMATURE_DAT_EOF;
        }
        for (varcoord.node[2] = 0; varcoord.node[2] < od->newgrid.nodes[2]; ++varcoord.node[2]) {
          for (j = 1, field_num = initial_field_num; j <= n_mo; ++j) {
            actual_len = fread(&value, sizeof(double), 1, gaussian_cube_in);
            if (actual_len != 1) {
              return PREMATURE_DAT_EOF;
            }
            fix_endianness(&value, sizeof(double), 1, endianness_switch);
            if ((found == -1) || (found == j)) {
              if ((!object_num) && (!varcoord.node[0])
                && (!varcoord.node[1]) && (!varcoord.node[2])) {
                /*
                allocate one more field
                */
                result = alloc_x_var_array(od, 1);
                if (result) {
                  return result;
                }
              }
              result = (mo ? set_x_value_xyz(od, field_num,
                object_num, &varcoord, value)
                : set_x_value_xyz_unbuffered(od, field_num,
                object_num, &varcoord, value));
              if (result) {
                return result;
              }
              ++field_num;
            }
          }
        }  
        actual_len = fread(&data_len_end, sizeof(int), 1, gaussian_cube_in);
        fix_endianness(&data_len_end, sizeof(int), 1, endianness_switch);
        if ((actual_len != 1) || (data_len_start != data_len_end)) {
          return PREMATURE_DAT_EOF;
        }
      }
    }
    if (gaussian_cube_in) {
      fclose(gaussian_cube_in);
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
