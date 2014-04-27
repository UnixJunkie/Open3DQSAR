/*

load_dat.c

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
#ifdef WIN32
#include <windows.h>
#endif

#define FIELD_NUM    0
#define OBJECT_NUM    1
#define OBJECT_ID    2
#define X_VARS_NUM    2
#define Y_VARS_NUM    3


int load_dat(O3Data *od, int file_id, int options)
{
  char buffer[LARGE_BUF_LEN];
  char header[BUF_LEN];
  char o3_header[TITLE_LEN + 1];
  char use_dummy_y_vars = 0;
  int actual_len;
  int object_num = 0;
  int old_object_num = 0;
  int temp_out = 0;
  int mol_len;
  int to_be_read;
  int fields_objects_x_y_vars[4] = { 0, 0, 0, 0 };
  int field_object_num_id[3] = { -1, -1, -1 };
  int endianness_switch;
  int i;
  int j;
  int result;
  int valid_size = 0;
  int all_predict = 1;
  uint64_t valid[1];
  float o3_version;
  float value;
  fzPtr *dat_in;
  int dummy;
  

  dat_in = (fzPtr *)(od->file[file_id]->handle);
  /*
  the first 4-byte word of the .dat file
  is an endianness indicator
  */
  actual_len = fzread(&i, sizeof(int), 1, dat_in);
  if (actual_len != 1) {
    O3_ERROR_LOCATE(&(od->task));
    return PREMATURE_DAT_EOF;
  }
  fix_endianness(&i, sizeof(int), 1, machine_type());
  
  if (i <= 0xFFFF) {
    endianness_switch = machine_type();
  }
  else {
    endianness_switch = 1 - machine_type();
  }
  /*
  then, read a short header about software version
  */
  memset(o3_header, ' ', TITLE_LEN);
  o3_header[TITLE_LEN] = '\0';
  actual_len = fzread(o3_header, 1, TITLE_LEN, dat_in);
  if (actual_len != TITLE_LEN) {
    O3_ERROR_LOCATE(&(od->task));
    return PREMATURE_DAT_EOF;
  }
  i = 13;
  while (o3_header[i] != ' ') {
    ++i;
  }
  o3_header[i] = '\0';
  o3_version = 1.5;
  if (!strncmp(o3_header, "Open3DQSAR", 10)) {
    sscanf(&o3_header[13], "%f", &o3_version);
  }
  temp_out = ((o3_version > 1.25) ? TEMP_MOLFILE : TEMP_OUT);
  /*
  then, read 4 int:
  - number of fields;
  - number of objects;
  - number of grid x_vars;
  - number of y variables.
  */
  actual_len = fzread(fields_objects_x_y_vars, sizeof(int), 4, dat_in);
  if (actual_len != 4) {
    O3_ERROR_LOCATE(&(od->task));
    return PREMATURE_DAT_EOF;
  }
  fix_endianness(fields_objects_x_y_vars, sizeof(int), 4, endianness_switch);
  if (options & VERBOSE_BIT) {
    tee_printf(od, "Header:                %s\n", o3_header);
    tee_printf(od, "Number of fields:      %d\n", fields_objects_x_y_vars[FIELD_NUM]);
    tee_printf(od, "Number of objects:     %d\n", fields_objects_x_y_vars[OBJECT_NUM]);
    tee_printf(od, "Number of variables:   %d\n", fields_objects_x_y_vars[X_VARS_NUM]
      * fields_objects_x_y_vars[FIELD_NUM] + fields_objects_x_y_vars[Y_VARS_NUM]);
    tee_printf(od, "Number of x variables: %d\n",
      fields_objects_x_y_vars[X_VARS_NUM] * fields_objects_x_y_vars[FIELD_NUM]);
    tee_printf(od, "Number of y variables: %d\n\n", fields_objects_x_y_vars[Y_VARS_NUM]);
  }
  if (options & APPEND_BIT) {
    if (od->field_num != fields_objects_x_y_vars[FIELD_NUM]) {
      O3_ERROR_LOCATE(&(od->task));
      return WRONG_NUMBER_OF_FIELDS;
    }
    if (od->x_vars != fields_objects_x_y_vars[X_VARS_NUM]) {
      O3_ERROR_LOCATE(&(od->task));
      return WRONG_NUMBER_OF_X_VARS;
    }
    if (od->y_vars != fields_objects_x_y_vars[Y_VARS_NUM]) {
      if (!fields_objects_x_y_vars[Y_VARS_NUM]) {
        /*
        if there are no Y vars in the loaded .dat file,
        we assume that the user does not know their value
        so we place a fake activity of 1.00
        (suggested by Claus Erhardt)
        */
        tee_printf(od,
          "Since y variables are not present in the appended dataset, while they were present\n"
          "in the pre-existing dataset, a dummy value of %.4lf will be assumed\n\n", DUMMY_Y_VAR_VALUE);
        use_dummy_y_vars = 1;
      }
      else {
        O3_ERROR_LOCATE(&(od->task));
        return WRONG_NUMBER_OF_Y_VARS;
      }
    }
  }
  /* then, read 6 float:
  - x, y, z start coordinates of 3D-grid;
  - x, y, z end coordinates of 3D-grid;
  */
  memset(&(od->newgrid), 0, sizeof(GridInfo));
  od->newgrid.x_vars = fields_objects_x_y_vars[X_VARS_NUM];
  actual_len = fzread(od->newgrid.start_coord, sizeof(float), 3, dat_in);
  if (actual_len != 3) {
    O3_ERROR_LOCATE(&(od->task));
    return PREMATURE_DAT_EOF;
  }
  fix_endianness(od->newgrid.start_coord, sizeof(float), 3, endianness_switch);
  actual_len = fzread(od->newgrid.end_coord, sizeof(float), 3, dat_in);
  if (actual_len != 3) {
    O3_ERROR_LOCATE(&(od->task));
    return PREMATURE_DAT_EOF;
  }
  fix_endianness(od->newgrid.end_coord, sizeof(float), 3, endianness_switch);
  
  /* then, read 3 int:
  - number of nodes in x, y, z directions;
  */
  actual_len = fzread(od->newgrid.nodes, sizeof(float), 3, dat_in);
  if (actual_len != 3) {
    O3_ERROR_LOCATE(&(od->task));
    return PREMATURE_DAT_EOF;
  }
  fix_endianness(od->newgrid.nodes, sizeof(float), 3, endianness_switch);
  /* then, read 3 float:
  - step size in x, y, z directions;
  */
  actual_len = fzread(od->newgrid.step, sizeof(float), 3, dat_in);
  if (actual_len != 3) {
    O3_ERROR_LOCATE(&(od->task));
    return PREMATURE_DAT_EOF;
  }
  fix_endianness(od->newgrid.step, sizeof(float), 3, endianness_switch);
  if (options & APPEND_BIT) {
    if (!match_grids(od)) {
      O3_ERROR_LOCATE(&(od->task));
      return GRID_NOT_MATCHING;
    }
    old_object_num = od->grid.object_num;
  }
  else {
    memcpy(&(od->grid), &(od->newgrid), sizeof(GridInfo));
  }
  if (!fields_objects_x_y_vars[OBJECT_NUM]) {
    O3_ERROR_LOCATE(&(od->task));
    return PREMATURE_DAT_EOF;
  }
  od->grid.object_num += fields_objects_x_y_vars[OBJECT_NUM];
  od->grid.struct_num += fields_objects_x_y_vars[OBJECT_NUM];
  /*
  alloc MolInfo structure array
  */
  if (options & APPEND_BIT) {
    if (!(od->al.mol_info = (MolInfo **)realloc
      (od->al.mol_info, (od->grid.object_num + 1) * sizeof(MolInfo *)))) {
      O3_ERROR_LOCATE(&(od->task));
      return OUT_OF_MEMORY;
    }
    memset(&(od->al.mol_info[old_object_num]), 0,
      (fields_objects_x_y_vars[OBJECT_NUM] + 1) * sizeof(MolInfo *));
    for (i = old_object_num; i < od->grid.object_num; ++i) {
      if (!(od->al.mol_info[i] = (MolInfo *)malloc(sizeof(MolInfo)))) {
        O3_ERROR_LOCATE(&(od->task));
        return OUT_OF_MEMORY;
      }
      memset(od->al.mol_info[i], 0, sizeof(MolInfo));
    }
  }
  else {
    free_array(od->al.mol_info);
    if (!(od->al.mol_info = (MolInfo **)
      alloc_array(od->grid.object_num, sizeof(MolInfo)))) {
      O3_ERROR_LOCATE(&(od->task));
      return OUT_OF_MEMORY;
    }
  }
  if (!use_dummy_y_vars) {
    od->y_vars = fields_objects_x_y_vars[Y_VARS_NUM];
  }
  if (alloc_object_attr(od, old_object_num)) {
    O3_ERROR_LOCATE(&(od->task));
    return OUT_OF_MEMORY;
  }
  if (!IS_O3A(od)) {
    od->x_vars = fields_objects_x_y_vars[X_VARS_NUM];
    od->object_num = od->grid.object_num;
    od->object_pagesize = (od->x_vars * sizeof(float) / od->mmap_pagesize +
      (od->x_vars * sizeof(float) % od->mmap_pagesize ? 1 : 0))
      * od->mmap_pagesize;
    if (fields_objects_x_y_vars[FIELD_NUM]) {
      if (options & APPEND_BIT) {
        if (realloc_x_var_array(od, old_object_num)) {
          O3_ERROR_LOCATE(&(od->task));
          return OUT_OF_MEMORY;
        }
      }
      else {
        if (alloc_x_var_array(od, fields_objects_x_y_vars[FIELD_NUM])) {
          O3_ERROR_LOCATE(&(od->task));
          return OUT_OF_MEMORY;
        }
      }
    }
  }
  
  /*
  then, read an uint16_t/long:
  - valid analyses
  */
  valid_size = ((o3_version < 2.105)
    ? sizeof(uint16_t) : sizeof(uint64_t));
  actual_len = fzread(valid, valid_size, 1, dat_in);
  if (actual_len != 1) {
    O3_ERROR_LOCATE(&(od->task));
    return PREMATURE_DAT_EOF;
  }
  if (!(options & APPEND_BIT)) {
    fix_endianness(valid, valid_size, 1, endianness_switch);
    od->valid = ((o3_version < 2.105)
      ? (unsigned long)(*((uint16_t *)valid)) : *valid);
  }

  /*
  then, read a number of XData structures equal to active_field_num
  */
  if (fields_objects_x_y_vars[FIELD_NUM]) {
    if (IS_O3A(od) || (options & APPEND_BIT)) {
      if (fzseek(dat_in, sizeof(XData)
        * fields_objects_x_y_vars[FIELD_NUM], SEEK_CUR)) {
        O3_ERROR_LOCATE(&(od->task));
        return PREMATURE_DAT_EOF;
      }
    }
    else {
      actual_len = fzread(od->mel.x_data, sizeof(XData),
        fields_objects_x_y_vars[FIELD_NUM], dat_in);
      if (actual_len != fields_objects_x_y_vars[FIELD_NUM]) {
        O3_ERROR_LOCATE(&(od->task));
        return PREMATURE_DAT_EOF;
      }
    }
  }
  /*
  start reading data from .dat file
  */
  while (!((field_object_num_id[FIELD_NUM] ==
    (fields_objects_x_y_vars[FIELD_NUM] - 1))
    && (field_object_num_id[OBJECT_NUM] == (od->grid.object_num - 1)))) {
    /*
    in the first block of MAX_NAME_LEN bytes
    there must be the DAT_HEADER header
    */
    actual_len = fzread(header, 1, MAX_NAME_LEN, dat_in);
    if (actual_len != MAX_NAME_LEN) {
      O3_ERROR_LOCATE(&(od->task));
      return PREMATURE_DAT_EOF;
    }
    if (strncmp(header, DAT_HEADER, 6)) {
      O3_ERROR_LOCATE(&(od->task));
      return PREMATURE_DAT_EOF;
    }
    /*
    in the second block of MAX_NAME_LEN bytes
    there must be the actual object name, from
    which we strip trailing spaces
    */
    actual_len = fzread(header, 1, MAX_NAME_LEN, dat_in);
    if (actual_len != MAX_NAME_LEN) {
      O3_ERROR_LOCATE(&(od->task));
      return PREMATURE_DAT_EOF;
    }
    i = 0;
    while ((header[i] != ' ') && (i < MAX_NAME_LEN)) {
      ++i;
    }
    header[i] = '\0';
    /*
    read 3 int:
    - field number
    - object number
    - object ID
    */
    actual_len = fzread(field_object_num_id, sizeof(int), 3, dat_in);
    if (actual_len != 3) {
      O3_ERROR_LOCATE(&(od->task));
      return PREMATURE_DAT_EOF;
    }
    fix_endianness(field_object_num_id, sizeof(int), 3, endianness_switch);
    if (options & APPEND_BIT) {
      field_object_num_id[OBJECT_NUM] += old_object_num;
      for (i = 0; i < od->grid.object_num; ++i) {
      }
      field_object_num_id[OBJECT_ID] += od->al.mol_info[old_object_num - 1]->object_id;
    }
    if ((field_object_num_id[OBJECT_NUM] >= od->grid.object_num)
      || (field_object_num_id[FIELD_NUM] >= fields_objects_x_y_vars[FIELD_NUM])) {
      O3_ERROR_LOCATE(&(od->task));
      return PREMATURE_DAT_EOF;
    }
    od->al.mol_info[field_object_num_id[OBJECT_NUM]]->object_id =
      field_object_num_id[OBJECT_ID];
    od->al.mol_info[field_object_num_id[OBJECT_NUM]]->struct_num =
      field_object_num_id[OBJECT_NUM];
    /*
    read MOL information
    */
    actual_len = fzread(&mol_len, sizeof(int), 1, dat_in);
    fix_endianness(&mol_len, sizeof(int), 1, endianness_switch);
    if (actual_len != 1) {
      O3_ERROR_LOCATE(&(od->task));
      return PREMATURE_DAT_EOF;
    }
    if (mol_len) {
      if (!(od->file[temp_out]->handle)) {
        if (!(od->file[temp_out]->handle = fopen
          (od->file[temp_out]->name, (options & APPEND_BIT) ? "ab" : "wb"))) {
          O3_ERROR_LOCATE(&(od->task));
          return CANNOT_WRITE_TEMP_FILE;
        }
      }
      while (mol_len > 0) {
        to_be_read = (((mol_len - LARGE_BUF_LEN) >= 0)
          ? LARGE_BUF_LEN : mol_len);
        actual_len = fzread(buffer, 1, to_be_read, dat_in);
        if (actual_len != to_be_read) {
          O3_ERROR_LOCATE(&(od->task));
          return PREMATURE_DAT_EOF;
        }
        actual_len = fwrite(buffer, 1, to_be_read,
          od->file[temp_out]->handle);
        if (actual_len != to_be_read) {
          O3_ERROR_LOCATE(&(od->task));
          return CANNOT_WRITE_TEMP_FILE;
        }
        mol_len -= to_be_read;
      }
      /*
      if molecules are stored as MOL files,
      add a $$$$ delimiter at the end of file
      */
      if (o3_version > 1.25) {
        fputs(SDF_DELIMITER"\n", od->file[temp_out]->handle);
      }
    }
    if (field_object_num_id[FIELD_NUM] >= 0) {
      if (IS_O3A(od)) {
        if (fzseek(dat_in, sizeof(float)
          * fields_objects_x_y_vars[X_VARS_NUM], SEEK_CUR)) {
          O3_ERROR_LOCATE(&(od->task));
          return PREMATURE_DAT_EOF;
        }
      }
      else {
        /*
        read x_vars
        */
        if (od->save_ram) {
          if (od->mmap_field_num != -1) {
            for (j = 0; j < od->object_num; ++j) {
              #ifndef WIN32
              munmap(od->mel.x_var_array[od->mmap_field_num][j], od->object_pagesize);
              #else
              UnmapViewOfFile(od->mel.x_var_array[od->mmap_field_num][j]);
              CloseHandle(od->al.mol_info[j]->hMapHandle);
              #endif
            }
            od->mmap_field_num = -1;
          }
          if (fseek(od->file[TEMP_FIELD_DATA + field_object_num_id[FIELD_NUM]]->handle,
            od->object_pagesize * field_object_num_id[OBJECT_NUM], SEEK_SET)) {
            O3_ERROR_LOCATE(&(od->task));
            return PREMATURE_DAT_EOF;
          }
          to_be_read = LARGE_BUF_LEN;
          for (j = (sizeof(float) * od->x_vars / LARGE_BUF_LEN); j >= 0; --j) {
            if (!j) {
              to_be_read = sizeof(float) * od->x_vars % LARGE_BUF_LEN;
            }
            actual_len = fzread(buffer, 1, to_be_read, dat_in);
            if (actual_len != to_be_read) {
              O3_ERROR_LOCATE(&(od->task));
              return PREMATURE_DAT_EOF;
            }
            actual_len = fwrite(buffer, 1, to_be_read,
              od->file[TEMP_FIELD_DATA + field_object_num_id[FIELD_NUM]]->handle);
            if (actual_len != to_be_read) {
              O3_ERROR_LOCATE(&(od->task));
              return PREMATURE_DAT_EOF;
            }
          }
        }
        else {
          actual_len = fzread(od->mel.x_var_array
            [field_object_num_id[FIELD_NUM]][field_object_num_id[OBJECT_NUM]],
            sizeof(float), od->x_vars, dat_in);
          if (actual_len != od->x_vars) {
            O3_ERROR_LOCATE(&(od->task));
            return PREMATURE_DAT_EOF;
          }
        }
      }
    }
  }
  if (od->file[temp_out]->handle) {
    fclose(od->file[temp_out]->handle);
    od->file[temp_out]->handle = NULL;
    if (o3_version < 1.25) {
      /*
      if the file was saved by Open3DQSAR version
      < 1.3, structure was stored as a MOL2 file
      and needs to be converted to MOL format
      */
      result = convert_mol(od, od->file[TEMP_OUT]->name,
        od->file[TEMP_MOLFILE]->name, "mol2", "sdf", "-xl");
      if (result) {
        O3_ERROR_LOCATE(&(od->task));
        return result;
      }
    }
    strcpy(od->file[MOLFILE_IN]->name, od->file[TEMP_MOLFILE]->name);
    if (!(od->file[MOLFILE_IN]->handle = fopen
      (od->file[MOLFILE_IN]->name, "rb"))) {
      O3_ERROR_LOCATE(&(od->task));
      return CANNOT_READ_TEMP_FILE;
    }
    result = parse_sdf(od, 0, NULL);
    fclose(od->file[MOLFILE_IN]->handle);
    od->file[MOLFILE_IN]->handle = NULL;
    if (result) {
      O3_ERROR_LOCATE(&(od->task));
      return result;
    }
  }
  /*
  Now load attributes:
  - field_attr
  - object_attr
  - object_weight (if o3_version > 2.03)
  - x_var_attr
  */
  if (fields_objects_x_y_vars[FIELD_NUM]) {
    /*
    if the APPEND_BIT is set, then skip
    field_attr and keep the current ones
    */
    if (IS_O3A(od) || (options & APPEND_BIT)) {
      if (fzseek(dat_in, sizeof(uint16_t)
        * fields_objects_x_y_vars[FIELD_NUM], SEEK_CUR)) {
        O3_ERROR_LOCATE(&(od->task));
        return PREMATURE_DAT_EOF;
      }
    }
    else {
      actual_len = fzread(od->mel.field_attr, sizeof(uint16_t),
        fields_objects_x_y_vars[FIELD_NUM], dat_in);
      if (actual_len != fields_objects_x_y_vars[FIELD_NUM]) {
        O3_ERROR_LOCATE(&(od->task));
        return PREMATURE_DAT_EOF;
      }
      fix_endianness(od->mel.field_attr, sizeof(uint16_t),
        fields_objects_x_y_vars[FIELD_NUM], endianness_switch);
    }
  }
  actual_len = fzread(&(od->mel.object_attr[old_object_num]),
    sizeof(uint16_t), fields_objects_x_y_vars[OBJECT_NUM], dat_in);
  if (actual_len != fields_objects_x_y_vars[OBJECT_NUM]) {
    O3_ERROR_LOCATE(&(od->task));
    return PREMATURE_DAT_EOF;
  }
  fix_endianness(&(od->mel.object_attr[old_object_num]),
    sizeof(uint16_t), fields_objects_x_y_vars[OBJECT_NUM],
    endianness_switch);
  if (o3_version > 2.03) {
    actual_len = fzread(&(od->mel.object_weight[old_object_num]),
      sizeof(double), fields_objects_x_y_vars[OBJECT_NUM], dat_in);
    if (actual_len != fields_objects_x_y_vars[OBJECT_NUM]) {
      O3_ERROR_LOCATE(&(od->task));
      return PREMATURE_DAT_EOF;
    }
    fix_endianness(&(od->mel.object_weight[old_object_num]),
      sizeof(double), fields_objects_x_y_vars[OBJECT_NUM],
      endianness_switch);
  }
  if (options & APPEND_BIT) {
    for (i = old_object_num; all_predict && (i < od->grid.object_num); ++i) {
      all_predict = get_object_attr(od, i, PREDICT_BIT);
    }
    if (!all_predict) {
      od->valid = SDF_BIT;
    }
  }
  for (i = 0; i < fields_objects_x_y_vars[FIELD_NUM]; ++i) {
    /*
    if the APPEND_BIT is set, then skip
    x_var_attr and keep the current ones
    */
    if (IS_O3A(od) || (options & APPEND_BIT)) {
      if (fzseek(dat_in, sizeof(uint16_t)
        * fields_objects_x_y_vars[X_VARS_NUM], SEEK_CUR)) {
        O3_ERROR_LOCATE(&(od->task));
        return PREMATURE_DAT_EOF;
      }
      /*
      if one or more imported objects are part of the
      training set, then set all x_vars as ACTIVE
      */
      if ((!IS_O3A(od)) && (!all_predict)) {
        for (j = 0; j < od->x_vars; ++j) {
          od->mel.x_var_attr[i][j] = ACTIVE_BIT;
        }
      }
    }
    else {
      actual_len = fzread(od->mel.x_var_attr[i], sizeof(uint16_t),
        fields_objects_x_y_vars[X_VARS_NUM], dat_in);
      if (actual_len != fields_objects_x_y_vars[X_VARS_NUM]) {
        O3_ERROR_LOCATE(&(od->task));
        return PREMATURE_DAT_EOF;
      }
      fix_endianness(od->mel.x_var_attr[i], sizeof(uint16_t),
        fields_objects_x_y_vars[X_VARS_NUM], endianness_switch);
    }
  }

  /*
  if a SRD analysis is present, then load
  relevant data
  */
  if ((od->valid & SEED_BIT) && ((!(options & APPEND_BIT))
    || ((options & APPEND_BIT) && all_predict))) {
    if (IS_O3Q(od)) {
      free_array(od->al.voronoi_composition);
      od->al.voronoi_composition = NULL;
      if ((options & VERBOSE_BIT)) {
        tee_printf(od, "SRD analysis present:\n\n"
          "%5s%20s\n", "Field", "Groups");
        tee_printf(od, "-------------------------\n");
      }
    }
    actual_len = fzread(&(od->voronoi_num), sizeof(int), 1, dat_in);
    if (actual_len != 1) {
      O3_ERROR_LOCATE(&(od->task));
      return PREMATURE_DAT_EOF;
    }
    if (IS_O3Q(od)) {
      if (alloc_voronoi(od, od->voronoi_num)) {
        O3_ERROR_LOCATE(&(od->task));
        return OUT_OF_MEMORY;
      }
      actual_len = fzread(od->mel.voronoi_active, sizeof(int),
        od->voronoi_num, dat_in);
      if (actual_len != od->voronoi_num) {
        O3_ERROR_LOCATE(&(od->task));
        return PREMATURE_DAT_EOF;
      }
    }
    else {
      if (fzseek(dat_in, od->voronoi_num * sizeof(int), SEEK_CUR)) {
        O3_ERROR_LOCATE(&(od->task));
        return PREMATURE_DAT_EOF;
      }
    }
    for (i = 0; i < od->voronoi_num; ++i) {
      if (IS_O3Q(od)) {
        actual_len = fzread(&(od->mel.voronoi_fill[i]),
          sizeof(int), 1, dat_in);
      }
      else {
        actual_len = fzread(&dummy, sizeof(int), 1, dat_in);
      }
      if (actual_len != 1) {
        O3_ERROR_LOCATE(&(od->task));
        return PREMATURE_DAT_EOF;
      }
      if (IS_O3Q(od)) {
        od->al.voronoi_composition[i] =
          alloc_int_array(NULL, od->mel.voronoi_fill[i]);
        if (!(od->al.voronoi_composition[i])) {
          O3_ERROR_LOCATE(&(od->task));
          return OUT_OF_MEMORY;
        }
        actual_len = fzread(od->al.voronoi_composition[i],
          sizeof(int), od->mel.voronoi_fill[i],
          dat_in);
        if (actual_len != od->mel.voronoi_fill[i]) {
          O3_ERROR_LOCATE(&(od->task));
          return PREMATURE_DAT_EOF;
        }
      }
      else {
        if (fzseek(dat_in, dummy * sizeof(int), SEEK_CUR)) {
          O3_ERROR_LOCATE(&(od->task));
          return PREMATURE_DAT_EOF;
        }
      }
    }
    for (i = 0; i < od->field_num; ++i) {
      if (IS_O3Q(od)) {
        actual_len = fzread(&(od->mel.seed_count[i]),
          sizeof(int), 1, dat_in);
        if (actual_len != 1) {
          O3_ERROR_LOCATE(&(od->task));
          return PREMATURE_DAT_EOF;
        }
        actual_len = fzread(&(od->mel.seed_count_before_collapse[i]),
          sizeof(int), 1, dat_in);
        if (actual_len != 1) {
          O3_ERROR_LOCATE(&(od->task));
          return PREMATURE_DAT_EOF;
        }
        actual_len = fzread(&(od->mel.group_zero[i]),
          sizeof(int), 1, dat_in);
        if (actual_len != 1) {
          O3_ERROR_LOCATE(&(od->task));
          return PREMATURE_DAT_EOF;
        }
        actual_len = fzread(od->cimal.voronoi_buf->me[i],
          sizeof(int), od->x_vars, dat_in);
        if (actual_len != od->x_vars) {
          O3_ERROR_LOCATE(&(od->task));
          return PREMATURE_DAT_EOF;
        }
        if ((options & VERBOSE_BIT)) {
          tee_printf(od, "%5d%20d\n", i + 1, od->mel.seed_count[i]);
        }
      }
      else {
        if (fzseek(dat_in, (od->x_vars + 3) * sizeof(int), SEEK_CUR)) {
          O3_ERROR_LOCATE(&(od->task));
          return PREMATURE_DAT_EOF;
        }
      }
    }
    if (IS_O3Q(od) && (options & VERBOSE_BIT)) {
      tee_printf(od, "\n");
    }
  }
  if (od->y_vars) {
    if (options & APPEND_BIT) {
      if (realloc_y_var_array(od, old_object_num)) {
        O3_ERROR_LOCATE(&(od->task));
        return OUT_OF_MEMORY;
      }
    }
    else {
      if (alloc_y_var_array(od)) {
        O3_ERROR_LOCATE(&(od->task));
        return OUT_OF_MEMORY;
      }
    }
    if (!(options & APPEND_BIT)) {
      free_char_matrix(od->cimal.y_var_name);
      od->cimal.y_var_name = NULL;
      od->cimal.y_var_name = alloc_char_matrix
        (od->cimal.y_var_name, od->y_vars, MAX_NAME_LEN);
      if (!(od->cimal.y_var_name)) {
        O3_ERROR_LOCATE(&(od->task));
        return OUT_OF_MEMORY;
      }
    }
    if (!use_dummy_y_vars) {
      /*
      read the DAT_HEADER string
      */
      actual_len = fzread(header, 1, MAX_NAME_LEN, dat_in);
      if (actual_len != MAX_NAME_LEN) {
        O3_ERROR_LOCATE(&(od->task));
        return PREMATURE_DAT_EOF;
      }
      if (strncmp(header, DAT_HEADER, 6)) {
        O3_ERROR_LOCATE(&(od->task));
        return PREMATURE_DAT_EOF;
      }
      for (i = 0; i < od->y_vars; ++i) {
        /*
        read dependent variable names
        */
        actual_len = fzread(header, 1, MAX_NAME_LEN, dat_in);
        if (actual_len != MAX_NAME_LEN) {
          O3_ERROR_LOCATE(&(od->task));
          return PREMATURE_DAT_EOF;
        }
        if (!(options & APPEND_BIT)) {
          j = 0;
          while ((header[j] != ' ') && (j < MAX_NAME_LEN)) {
            ++j;
          }
          header[j] = '\0';
          strcpy(od->cimal.y_var_name->me[i], header);
        }
      }
      /*
      then, read a number of YData structures equal to y_vars
      */
      if (options & APPEND_BIT) {
        if (fzseek(dat_in, sizeof(YData) * od->y_vars, SEEK_CUR)) {
          O3_ERROR_LOCATE(&(od->task));
          return PREMATURE_DAT_EOF;
        }
      }
      else {
        actual_len = fzread(od->mel.y_data,
          sizeof(YData), od->y_vars, dat_in);
        if (actual_len != od->y_vars) {
          O3_ERROR_LOCATE(&(od->task));
          return PREMATURE_DAT_EOF;
        }
      }
      object_num = -1;
      while (object_num != (od->grid.object_num - 1)) {
        /*
        read 1 int:
        - object number
        */
        actual_len = fzread(&object_num, sizeof(int), 1, dat_in);
        fix_endianness(&object_num, sizeof(int), 1, endianness_switch);
        if (actual_len != 1) {
          O3_ERROR_LOCATE(&(od->task));
          return PREMATURE_DAT_EOF;
        }
        object_num += old_object_num;
        if ((object_num < old_object_num) || (object_num >= od->grid.object_num)) {
          O3_ERROR_LOCATE(&(od->task));
          return PREMATURE_DAT_EOF;
        }
        /*
        read y_vars
        */
        for (j = 0; j < od->y_vars; ++j) {
          actual_len = fzread(&value, sizeof(float), 1, dat_in);
          if (actual_len != 1) {
            O3_ERROR_LOCATE(&(od->task));
            return PREMATURE_DAT_EOF;
          }
          fix_endianness(&value, sizeof(float), 1, endianness_switch);
          set_y_value(od, object_num, j, (double)value);
        }
      }
      /*
      Now load y_var_attr attributes
      */
      if (options & APPEND_BIT) {
        if (fzseek(dat_in, sizeof(uint16_t) * od->y_vars, SEEK_CUR)) {
          O3_ERROR_LOCATE(&(od->task));
          return PREMATURE_DAT_EOF;
        }
      }
      else {
        actual_len = fzread(od->mel.y_var_attr, sizeof(uint16_t),
          od->y_vars, dat_in);
        fix_endianness(od->mel.y_var_attr, sizeof(float),
          od->y_vars, endianness_switch);
        if (actual_len != od->y_vars) {
          O3_ERROR_LOCATE(&(od->task));
          return PREMATURE_DAT_EOF;
        }
      }
    }
    else {
      for (object_num = old_object_num;
        object_num < od->grid.object_num; ++object_num) {
        for (j = 0; j < od->y_vars; ++j) {
          set_y_value(od, object_num, j, DUMMY_Y_VAR_VALUE);
        }
      }
    }
  }
  if (dat_in) {
    fzclose(dat_in);
  }
  od->file[file_id]->handle = NULL;
  if (od->save_ram) {
    sync_field_mmap(od);
  }
  
  update_field_object_attr(od, VERBOSE_BIT);

  result = calc_active_vars(od, FULL_MODEL);

  return result;
}
