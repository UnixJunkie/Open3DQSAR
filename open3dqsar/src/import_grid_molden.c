/*

import_grid_molden.c

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
#ifdef WIN32
#include <windows.h>
#endif


int import_grid_molden(O3Data *od)
{
  int i;
  int x;
  int y;
  int actual_len;
  int object_num;
  int n_atoms;
  int data_len_start;
  int data_len_end;
  int endianness_switch;
  int iplat;
  int result;
  int n_xy_temp;
  float *float_xy_mat = NULL;
  double value;
  double adjus;
  double p[3];
  double c[3];
  double r[3];
  FILE *molden_in;
  


  /*
  allocate one more field
  */
  if (alloc_x_var_array(od, 1)) {
    return OUT_OF_MEMORY;
  }
  rewind(od->file[TEMP_SORTED_MATCH]->handle);
  for (object_num = 0; object_num < od->object_num; ++object_num) {
    if (od->save_ram) {
      if (fseek(od->file
        [TEMP_FIELD_DATA + od->field_num - 1]->handle,
        od->object_pagesize * object_num, SEEK_SET)) {
        return CANNOT_WRITE_TEMP_FILE;
      }
    }
    if (!fgets(od->file[BINARY_IN]->name, BUF_LEN,
      od->file[TEMP_SORTED_MATCH]->handle)) {
      return CANNOT_READ_TEMP_FILE;
    }
    remove_newline(od->file[BINARY_IN]->name);
    molden_in = fopen(od->file[BINARY_IN]->name, "rb");
    if (!molden_in) {
      return NOT_ENOUGH_OBJECTS;
    }
    od->file[BINARY_IN]->handle = molden_in;
    /*
    desume endianness
    from the first int at the very beginning
    of the header
    also take into consideration
    the endianness of the machine we are running on
    */

    actual_len = fread(&data_len_start, sizeof(int), 1, molden_in);
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
    read number of atoms
    */
    actual_len = fread(&n_atoms, sizeof(int), 1, molden_in);
    if (actual_len != 1) {
      return PREMATURE_DAT_EOF;
    }
    fix_endianness(&n_atoms, sizeof(int), 1, endianness_switch);
    actual_len = fread(&data_len_end, sizeof(int), 1, molden_in);
    if (actual_len != 1) {
      return PREMATURE_DAT_EOF;
    }
    fix_endianness(&data_len_end, sizeof(int), 1, endianness_switch);
    if (data_len_start != data_len_end) {
      return PREMATURE_DAT_EOF;
    }
    actual_len = fread(&data_len_start, sizeof(int), 1, molden_in);
    if (actual_len != 1) {
      return PREMATURE_DAT_EOF;
    }
    fix_endianness(&data_len_start, sizeof(int), 1, endianness_switch);
    /*
    skip the atom numbers (not interesting)
    */
    if (fseek(molden_in, data_len_start, SEEK_CUR)) {
      return PREMATURE_DAT_EOF;
    }
    actual_len = fread(&data_len_end, sizeof(int), 1, molden_in);
    fix_endianness(&data_len_end, sizeof(int), 1, endianness_switch);
    if (actual_len != 1) {
      return PREMATURE_DAT_EOF;
    }
    fix_endianness(&data_len_end, sizeof(int), 1, endianness_switch);
    if (data_len_start != data_len_end) {
      return PREMATURE_DAT_EOF;
    }
    actual_len = fread(&data_len_start, sizeof(int), 1, molden_in);
    if (actual_len != 1) {
      return PREMATURE_DAT_EOF;
    }
    fix_endianness(&data_len_start, sizeof(int), 1, endianness_switch);
    /*
    we expect a double; it should be the bohr-to-angstrom conversion factor
    */
    if (data_len_start != sizeof(double)) {
      return PREMATURE_DAT_EOF;
    }
    actual_len = fread(&adjus, sizeof(double), 1, molden_in);
    if (actual_len != 1) {
      return PREMATURE_DAT_EOF;
    }
    fix_endianness(&adjus, sizeof(double), 1, endianness_switch);
    actual_len = fread(&data_len_end, sizeof(int), 1, molden_in);
    if (actual_len != 1) {
      return PREMATURE_DAT_EOF;
    }
    fix_endianness(&data_len_end, sizeof(int), 1, endianness_switch);
    if (data_len_start != data_len_end) {
      return PREMATURE_DAT_EOF;
    }
    /*
    we expect n_atoms * 3 doubles, corresponding to atomic coordinates
    we are not interested, so let's skip them
    */
    for (i = 0; i < n_atoms; ++i) {
      actual_len = fread(&data_len_start, sizeof(int), 1, molden_in);
      if (actual_len != 1) {
        return PREMATURE_DAT_EOF;
      }
      fix_endianness(&data_len_start, sizeof(int), 1, endianness_switch);
      if (data_len_start != (3 * sizeof(double))) {
        return PREMATURE_DAT_EOF;
      }
      /*
      skip the 3 doubles
      */
      if (fseek(molden_in, data_len_start, SEEK_CUR)) {
        return PREMATURE_DAT_EOF;
      }
      actual_len = fread(&data_len_end, sizeof(int), 1, molden_in);
      if (actual_len != 1) {
        return PREMATURE_DAT_EOF;
      }
      fix_endianness(&data_len_end, sizeof(int), 1, endianness_switch);
      if (data_len_start != data_len_end) {
        return PREMATURE_DAT_EOF;
      }
    }
    /*
    we expect 9 doubles (px, py, pz, cx, cy, cz, rx, ry, rz)
    followed by 4 ints (nptsx, nptsy, nptsz, iplat)
    */
    actual_len = fread(&data_len_start, sizeof(int), 1, molden_in);
    if (actual_len != 1) {
      return PREMATURE_DAT_EOF;
    }
    fix_endianness(&data_len_start, sizeof(int), 1, endianness_switch);
    if (data_len_start != (9 * sizeof(double) + 4 * sizeof(int))) {
      return PREMATURE_DAT_EOF;
    }
    for (i = 0; i < 3; ++i) {
      actual_len = fread(&p[i], sizeof(double), 1, molden_in);
      if (actual_len != 1) {
        return PREMATURE_DAT_EOF;
      }
      fix_endianness(&p[i], sizeof(double), 1, endianness_switch);
    }
    for (i = 0; i < 3; ++i) {
      actual_len = fread(&c[i], sizeof(double), 1, molden_in);
      if (actual_len != 1) {
        return PREMATURE_DAT_EOF;
      }
      fix_endianness(&c[i], sizeof(double), 1, endianness_switch);
    }
    for (i = 0; i < 3; ++i) {
      actual_len = fread(&r[i], sizeof(double), 1, molden_in);
      if (actual_len != 1) {
        return PREMATURE_DAT_EOF;
      }
      fix_endianness(&r[i], sizeof(double), 1, endianness_switch);
    }
    for (i = 0; i < 3; ++i) {
      actual_len = fread(&(od->newgrid.nodes[i]), sizeof(int), 1, molden_in);
      if (actual_len != 1) {
        return PREMATURE_DAT_EOF;
      }
      fix_endianness(&(od->newgrid.nodes[i]), sizeof(int), 1, endianness_switch);
    }
    actual_len = fread(&iplat, sizeof(int), 1, molden_in);
    if (actual_len != 1) {
      return PREMATURE_DAT_EOF;
    }
    fix_endianness(&iplat, sizeof(int), 1, endianness_switch);

    actual_len = fread(&data_len_end, sizeof(int), 1, molden_in);
    if (actual_len != 1) {
      return PREMATURE_DAT_EOF;
    }
    fix_endianness(&data_len_end, sizeof(int), 1, endianness_switch);
    if (data_len_start != data_len_end) {
      return PREMATURE_DAT_EOF;
    }
    od->newgrid.x_vars = 1;
    for (i = 0; i < 3; ++i) {
                  od->newgrid.start_coord[i] = (float)((p[i] - r[i] / (double)2) * adjus);
                  od->newgrid.end_coord[i] = (float)((p[i] + r[i] / (double)2) * adjus);
      od->newgrid.step[i] =
        (od->newgrid.end_coord[i] - od->newgrid.start_coord[i])
        / (float)(od->newgrid.nodes[i] - 1);
      od->newgrid.x_vars *= od->newgrid.nodes[i];
    }
    if (!match_grids(od)) {
      return GRID_NOT_MATCHING;
                }
    /*
    read binary data from the .kont file until EOF
    and store it in a temporary x*y matrix; after
    reading a whole xy plane, store it into memory
    */
    float_xy_mat = (float *)realloc(od->mel.float_xy_mat,
      od->newgrid.nodes[0] *
      od->newgrid.nodes[1] * sizeof(float));
    if (!float_xy_mat) {
      return OUT_OF_MEMORY;
    }
    od->mel.float_xy_mat = float_xy_mat;
    for (i = 0; i < od->newgrid.nodes[2]; ++i) {
      /*
      now we expect od->newgrid.z_nodes blocks each containing
      od->newgrid.x_nodes * od->newgrid.y_nodes doubles
      */
      actual_len = fread(&data_len_start, sizeof(int), 1, molden_in);
      if (actual_len != 1) {
        return PREMATURE_DAT_EOF;
      }
      fix_endianness(&data_len_start, sizeof(int), 1, endianness_switch);
      if (data_len_start != (od->newgrid.nodes[0]
        * od->newgrid.nodes[1] * sizeof(double))) {
        return PREMATURE_DAT_EOF;
      }
      
      x = 0;
      y = 0;
      n_xy_temp = 0;
      while (n_xy_temp < (od->grid.nodes[0] * od->grid.nodes[1])) {
        actual_len = fread(&value, sizeof(double), 1, molden_in);
        if (actual_len != 1) {
          return PREMATURE_DAT_EOF;
        }
        fix_endianness(&value, sizeof(double), 1, endianness_switch);
        /*
        store the xy plane
        */
        float_xy_mat[y * od->grid.nodes[0] + x] =
          (float)value;
        ++n_xy_temp;
        ++y;
        if (y == od->grid.nodes[1]) {
          y = 0;
          ++x;
        }
      }
      /*
      store the temporary matrix
      */
      if (od->save_ram) {
        actual_len = fwrite(float_xy_mat, sizeof(float),
        od->grid.nodes[0] * od->grid.nodes[1],
        od->file[TEMP_FIELD_DATA + od->field_num - 1]->handle);
        if (actual_len != (od->grid.nodes[0] * od->grid.nodes[1])) {
          return CANNOT_WRITE_TEMP_FILE;
        }
      }
      else {
        memcpy(&(od->mel.x_var_array[od->field_num - 1][object_num]
          [i * od->grid.nodes[0] * od->grid.nodes[1]]), float_xy_mat,
          sizeof(float) * od->grid.nodes[0] * od->grid.nodes[1]);
      }

      actual_len = fread(&data_len_end, sizeof(int), 1, molden_in);
      if (actual_len != 1) {
        return PREMATURE_DAT_EOF;
      }
      fix_endianness(&data_len_end, sizeof(int), 1, endianness_switch);
      if (data_len_start != data_len_end) {
        return PREMATURE_DAT_EOF;
      }
    }
    if (molden_in) {
      fclose(molden_in);
      od->file[BINARY_IN]->handle = NULL;
    }
  }
  if (float_xy_mat) {
    free(float_xy_mat);
    od->mel.float_xy_mat = NULL;
  }
  if (od->save_ram) {
    sync_field_mmap(od);
  }
  update_field_object_attr(od, VERBOSE_BIT);

  result = calc_active_vars(od, FULL_MODEL);

  return result;
}
