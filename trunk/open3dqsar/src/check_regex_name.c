/*

check_regex_name.c

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
#include <shlwapi.h>
#endif


int check_regex_name(char *regex_name, int n_regex)
{
  char number[BUF_LEN];
  int i;
  int n;
  int n_percent;
  int n_zero;
  
  
  i = 0;
  n_percent = 0;
  while (regex_name[i]) {
    if (regex_name[i] == '%') {
      ++n_percent;
      if (n_percent > n_regex) {
        return PARSE_INPUT_ERROR;
      }
      ++i;
      if (i >= BUF_LEN) {
        return PARSE_INPUT_ERROR;
      }
      if ((regex_name[i] != '0') && (regex_name[i] != 'd')) {
        return PARSE_INPUT_ERROR;
      }
      ++i;
      if (i >= BUF_LEN) {
        return PARSE_INPUT_ERROR;
      }
      if (regex_name[i - 1] == '0') {
        n = 0;
        while (regex_name[i] != 'd') {
          number[n] = regex_name[i];
          ++n;
          ++i;
          if (i >= BUF_LEN) {
            return PARSE_INPUT_ERROR;
          }
        }
        number[n] = '\0';
        ++i;
        if (i >= BUF_LEN) {
          return PARSE_INPUT_ERROR;
        }
        sscanf(number, "%d", &n_zero);
        if ((n_zero <= 0) || (n_zero > 9)) {
          return PARSE_INPUT_ERROR;
        }
      }
    }
    else {
      ++i;
      if (i >= BUF_LEN) {
        return PARSE_INPUT_ERROR;
      }
    }
  }

  return 0;
}


int get_datafile_coord(O3Data *od, FileDescriptor *data_fd, int n_atom,
  int total_n_atoms, int *cube_word_size, double *data_coord, int datafile_type)
{
  char atombuf[8];
  char buffer[BUF_LEN];
  int n;
  int n_skip = 0;
  int data_pos;
  int found;
  int n_data_atoms;
  int actual_len;
  int data_len_start;
  int data_len_end;
  int endianness_switch = 0;
  long int long_n_atoms;


  if (!(data_fd->handle)) {
    if (!fgets(data_fd->name, BUF_LEN, od->file[TEMP_OBJECT_MATCH]->handle)) {
      O3_ERROR_LOCATE(&(od->task));
      O3_ERROR_STRING(&(od->task), od->file[TEMP_OBJECT_MATCH]->name);
      return CANNOT_READ_TEMP_FILE;
    }
    remove_newline(data_fd->name);
    if (!(data_fd->handle = fopen(data_fd->name, "rb"))) {
      O3_ERROR_LOCATE(&(od->task));
      O3_ERROR_STRING(&(od->task), data_fd->name);
      return PREMATURE_DAT_EOF;
    }
    switch (datafile_type) {
      case COSMO_INPUT_FILE:
      found = 0;
      while ((!found) && fgets(buffer, BUF_LEN, data_fd->handle)) {
        found = (!strncasecmp(buffer, "$coord_car", 10));
        if (!found) {
          strncpy(atombuf, buffer, 4);
          atombuf[4] = '\0';
        }
      }
      if (!found) {
        O3_ERROR_LOCATE(&(od->task));
        O3_ERROR_STRING(&(od->task), data_fd->name);
        return PREMATURE_DAT_EOF;
      }
      sscanf(atombuf, "%d", &n_data_atoms);
      if (n_data_atoms != total_n_atoms) {
        return OBJECTS_NOT_MATCHING;
      }
      data_pos = ftell(data_fd->handle);
      found = 0;
      n = 0;
      while ((!found) && fgets(buffer, BUF_LEN, data_fd->handle)) {
        found = (!strncasecmp(buffer, "end", 3));
        if (!found) {
          ++n;
        }
      }
      if (!found) {
        O3_ERROR_LOCATE(&(od->task));
        O3_ERROR_STRING(&(od->task), data_fd->name);
        return PREMATURE_DAT_EOF;
      }
      n_skip = n - n_data_atoms;
      if (n_skip < 0) {
        O3_ERROR_LOCATE(&(od->task));
        O3_ERROR_STRING(&(od->task), data_fd->name);
        return PREMATURE_DAT_EOF;
      }
      if (fseek(data_fd->handle, data_pos, SEEK_SET)) {
        O3_ERROR_LOCATE(&(od->task));
        O3_ERROR_STRING(&(od->task), data_fd->name);
        return PREMATURE_DAT_EOF;
      }
      break;
      
      case FORMATTED_CUBE_INPUT_FILE:
      if (!(found = fgrep(data_fd->handle, buffer, "$CUBE"))) {
        found = fgrep(data_fd->handle, buffer, "START OF CUBE FORMAT");
      }
      if (!found) {
        rewind(data_fd->handle);
      }
      /*
      skip two lines
      */
      found = 1;
      for (n = 0; (n < 2) && found; ++n) {
        found = (fgets(buffer, BUF_LEN, data_fd->handle) ? 1 : 0);
      }
      if (!found) {
        return PREMATURE_DAT_EOF;
      }
      /*
      read number of atoms
      */
      if (!fgets(buffer, BUF_LEN, data_fd->handle)) {
        return PREMATURE_DAT_EOF;
      }
      buffer[BUF_LEN - 1] = '\0';
      sscanf(buffer, "%d", &n_data_atoms);
      n_data_atoms = absval(n_data_atoms);
      if (n_data_atoms != total_n_atoms) {
        return OBJECTS_NOT_MATCHING;
      }
      /*
      now skip 3 lines
      */
      n_skip = 3;
      break;
      
      case UNFORMATTED_CUBE_INPUT_FILE:
      /*
      desume endianness from the first int at the
      very beginning of the cube file
      also take into consideration
      the endianness of the machine we are running on
      */
      actual_len = fread(&data_len_start, sizeof(int), 1, data_fd->handle);
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
      if (fseek(data_fd->handle, data_len_start, SEEK_CUR)) {
        return PREMATURE_DAT_EOF;
      }
      actual_len = fread(&data_len_end, sizeof(int), 1, data_fd->handle);
      fix_endianness(&data_len_end, sizeof(int), 1, endianness_switch);
      if ((actual_len != 1) || (data_len_start != data_len_end)) {
        return PREMATURE_DAT_EOF;
      }
      actual_len = fread(&data_len_start, sizeof(int), 1, data_fd->handle);
      if (actual_len != 1) {
        return PREMATURE_DAT_EOF;
      }
      fix_endianness(&data_len_start, sizeof(int), 1, endianness_switch);
      if (fseek(data_fd->handle, data_len_start, SEEK_CUR)) {
        return PREMATURE_DAT_EOF;
      }
      actual_len = fread(&data_len_end, sizeof(int), 1, data_fd->handle);
      fix_endianness(&data_len_end, sizeof(int), 1, endianness_switch);
      if ((actual_len != 1) || (data_len_start != data_len_end)) {
        return PREMATURE_DAT_EOF;
      }
      /*
      now there should be a block constituted by:
      - 1 int or long int (number of atoms)
      - 3 doubles (origin x, y, z; we may skip those)
      */
      actual_len = fread(&data_len_start, sizeof(int), 1, data_fd->handle);
      if (actual_len != 1) {
        return PREMATURE_DAT_EOF;
      }
      fix_endianness(&data_len_start, sizeof(int), 1, endianness_switch);
      data_len_start -= 3 * sizeof(double);
      if ((data_len_start != sizeof(int))
        && (data_len_start != sizeof(long int))) {
        return PREMATURE_DAT_EOF;
      }
      *cube_word_size = data_len_start;
      actual_len = fread(&long_n_atoms, *cube_word_size, 1, data_fd->handle);
      if (actual_len != 1) {
        return PREMATURE_DAT_EOF;
      }
      fix_endianness(&long_n_atoms, *cube_word_size, 1, endianness_switch);
      n_data_atoms = (int)long_n_atoms;
      n_data_atoms = absval(n_data_atoms);
      if (n_data_atoms != total_n_atoms) {
        return OBJECTS_NOT_MATCHING;
      }
      if (fseek(data_fd->handle, 3 * sizeof(double), SEEK_CUR)) {
        return PREMATURE_DAT_EOF;
      }
      actual_len = fread(&data_len_end, sizeof(int), 1, data_fd->handle);
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
      we may skip all of them
      */
      for (n = 0; n < 3; ++n) {
        actual_len = fread(&data_len_start, sizeof(int), 1, data_fd->handle);
        if (actual_len != 1) {
          return PREMATURE_DAT_EOF;
        }
        fix_endianness(&data_len_start, sizeof(int), 1, endianness_switch);
        data_len_start -= 3 * sizeof(double);
        if (data_len_start != *cube_word_size) {
          return PREMATURE_DAT_EOF;
        }
        if (fseek(data_fd->handle, *cube_word_size + 3 * sizeof(double), SEEK_CUR)) {
          return PREMATURE_DAT_EOF;
        }
        actual_len = fread(&data_len_end, sizeof(int), 1, data_fd->handle);
        if (actual_len != 1) {
          return PREMATURE_DAT_EOF;
        }
        fix_endianness(&data_len_end, sizeof(int), 1, endianness_switch);
        data_len_end -= 3 * sizeof(double);
        if (data_len_start != data_len_end) {
          return PREMATURE_DAT_EOF;
        }
      }
      n_skip = 0;
      break;
      
      case MOLDEN_INPUT_FILE:
      /*
      desume endianness from the first int at the
      very beginning of the cube file
      also take into consideration
      the endianness of the machine we are running on
      */
      actual_len = fread(&data_len_start, sizeof(int), 1, data_fd->handle);
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
      actual_len = fread(&n_data_atoms, sizeof(int), 1, data_fd->handle);
      if (actual_len != 1) {
        return PREMATURE_DAT_EOF;
      }
      fix_endianness(&n_data_atoms, sizeof(int), 1, endianness_switch);
      actual_len = fread(&data_len_end, sizeof(int), 1, data_fd->handle);
      if (actual_len != 1) {
        return PREMATURE_DAT_EOF;
      }
      if (n_data_atoms != total_n_atoms) {
        return OBJECTS_NOT_MATCHING;
      }
      fix_endianness(&data_len_end, sizeof(int), 1, endianness_switch);
      if (data_len_start != data_len_end) {
        return PREMATURE_DAT_EOF;
      }
      actual_len = fread(&data_len_start, sizeof(int), 1, data_fd->handle);
      if (actual_len != 1) {
        return PREMATURE_DAT_EOF;
      }
      fix_endianness(&data_len_start, sizeof(int), 1, endianness_switch);
      /*
      skip the atom numbers (not interesting)
      */
      if (fseek(data_fd->handle, data_len_start, SEEK_CUR)) {
        return PREMATURE_DAT_EOF;
      }
      actual_len = fread(&data_len_end, sizeof(int), 1, data_fd->handle);
      fix_endianness(&data_len_end, sizeof(int), 1, endianness_switch);
      if (actual_len != 1) {
        return PREMATURE_DAT_EOF;
      }
      fix_endianness(&data_len_end, sizeof(int), 1, endianness_switch);
      if (data_len_start != data_len_end) {
        return PREMATURE_DAT_EOF;
      }
      actual_len = fread(&data_len_start, sizeof(int), 1, data_fd->handle);
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
      if (fseek(data_fd->handle, sizeof(double), SEEK_CUR)) {
        return PREMATURE_DAT_EOF;
      }
      actual_len = fread(&data_len_end, sizeof(int), 1, data_fd->handle);
      if (actual_len != 1) {
        return PREMATURE_DAT_EOF;
      }
      fix_endianness(&data_len_end, sizeof(int), 1, endianness_switch);
      if (data_len_start != data_len_end) {
        return PREMATURE_DAT_EOF;
      }
      /*
      now we expect n_atoms * 3 doubles, corresponding to atomic coordinates
      */
      n_skip = 0;
      break;
    }
    n = 0;
    while ((n < n_skip) && fgets(buffer, BUF_LEN, data_fd->handle)) {
      ++n;
    }
    if (n != n_skip) {
      O3_ERROR_LOCATE(&(od->task));
      O3_ERROR_STRING(&(od->task), data_fd->name);
      return PREMATURE_DAT_EOF;
    }
  }
  switch (datafile_type) {
    case COSMO_INPUT_FILE:
    case FORMATTED_CUBE_INPUT_FILE:
    if (!fgets(buffer, BUF_LEN, data_fd->handle)) {
      O3_ERROR_LOCATE(&(od->task));
      O3_ERROR_STRING(&(od->task), data_fd->name);
      return PREMATURE_DAT_EOF;
    }
    buffer[BUF_LEN - 1] = '\0';
    if (datafile_type == COSMO_INPUT_FILE) {
      sscanf(buffer, "%*s %lf %lf %lf",
        &data_coord[0], &data_coord[1], &data_coord[2]);
    }
    else {
      sscanf(buffer, "%*s %*s %lf %lf %lf",
        &data_coord[0], &data_coord[1], &data_coord[2]);
      for (n = 0; n < 3; ++n) {
        data_coord[n] *= BOHR_RADIUS;
      }
    }
    break;
    
    case UNFORMATTED_CUBE_INPUT_FILE:
    /*
    - 1 int or long int (atom number)
    - 4 doubles (charge, x_coord, y_coord, z_coord)
    */
    if (!n_atom) {
      actual_len = fread(&data_len_start, sizeof(int), 1, data_fd->handle);
      fix_endianness(&data_len_start, sizeof(int), 1, endianness_switch);
      if ((actual_len != 1) || (data_len_start !=
        (total_n_atoms * (*cube_word_size + 4 * sizeof(double))))) {
        O3_ERROR_LOCATE(&(od->task));
        O3_ERROR_STRING(&(od->task), data_fd->name);
        return PREMATURE_DAT_EOF;
      }
    }
    if (fseek(data_fd->handle, *cube_word_size + sizeof(double), SEEK_CUR)) {
      O3_ERROR_LOCATE(&(od->task));
      O3_ERROR_STRING(&(od->task), data_fd->name);
      return PREMATURE_DAT_EOF;
    }
    actual_len = fread(data_coord, sizeof(double), 3, data_fd->handle);
    if (actual_len != 3) {
      O3_ERROR_LOCATE(&(od->task));
      O3_ERROR_STRING(&(od->task), data_fd->name);
      return PREMATURE_DAT_EOF;
    }
    fix_endianness(data_coord, sizeof(double), 3, endianness_switch);
    for (n = 0; n < 3; ++n) {
      data_coord[n] *= BOHR_RADIUS;
    }
    if (n_atom == total_n_atoms) {
      actual_len = fread(&data_len_end, sizeof(int), 1, data_fd->handle);
      fix_endianness(&data_len_end, sizeof(int), 1, endianness_switch);
      if ((actual_len != 1) || (data_len_start != data_len_end)) {
        O3_ERROR_LOCATE(&(od->task));
        O3_ERROR_STRING(&(od->task), data_fd->name);
        return PREMATURE_DAT_EOF;
      }
    }
    break;
    
    case MOLDEN_INPUT_FILE:
    actual_len = fread(&data_len_start, sizeof(int), 1, data_fd->handle);
    if (actual_len != 1) {
      return PREMATURE_DAT_EOF;
    }
    fix_endianness(&data_len_start, sizeof(int), 1, endianness_switch);
    if (data_len_start != (3 * sizeof(double))) {
      return PREMATURE_DAT_EOF;
    }
    /*
    read the coordinates
    */
    actual_len = fread(data_coord, sizeof(double), 3, data_fd->handle);
    if (actual_len != 3) {
      O3_ERROR_LOCATE(&(od->task));
      O3_ERROR_STRING(&(od->task), data_fd->name);
      return PREMATURE_DAT_EOF;
    }
    fix_endianness(data_coord, sizeof(double), 3, endianness_switch);
    for (n = 0; n < 3; ++n) {
      data_coord[n] *= BOHR_RADIUS;
    }
    actual_len = fread(&data_len_end, sizeof(int), 1, data_fd->handle);
    if (actual_len != 1) {
      return PREMATURE_DAT_EOF;
    }
    fix_endianness(&data_len_end, sizeof(int), 1, endianness_switch);
    if (data_len_start != data_len_end) {
      return PREMATURE_DAT_EOF;
    }
    break;
  }
  
  return 0;
}


int match_objects_with_datafile(O3Data *od, char *file_pattern, int datafile_type)
{
  char buffer[BUF_LEN];
  char *datafile_base;
  int object_num;
  int n_sdf_atoms;
  int sdf_pos;
  int cube_word_size;
  int attempts;
  int found;
  int result;
  double sdf_coord[3];
  double data_coord[3];
  FileDescriptor mol_fd;
  FileDescriptor data_fd;
  
  
  memset(&mol_fd, 0, sizeof(FileDescriptor));
  memset(&data_fd, 0, sizeof(FileDescriptor));
  memset(buffer, 0, BUF_LEN);
  result = check_file_pattern(od->file[TEMP_OBJECT_MATCH]->handle, file_pattern, od->object_num);
  if (result) {
    return result;
  }
  for (object_num = 0; object_num < od->object_num; ++object_num) {
    rewind(od->file[TEMP_OBJECT_MATCH]->handle);
    sprintf(mol_fd.name, "%s%c%04d.mol", od->field.mol_dir,
      SEPARATOR, od->al.mol_info[object_num]->object_id);
    if (!(mol_fd.handle = fopen(mol_fd.name, "rb"))) {
      O3_ERROR_LOCATE(&(od->task));
      O3_ERROR_STRING(&(od->task), mol_fd.name);
      return CANNOT_READ_ORIGINAL_SDF;
    }
    if (find_conformation_in_sdf(mol_fd.handle, NULL, 0)) {
      O3_ERROR_LOCATE(&(od->task));
      O3_ERROR_STRING(&(od->task), mol_fd.name);
      fclose(mol_fd.handle);
      return CANNOT_READ_ORIGINAL_SDF;
    }
    sdf_pos = ftell(mol_fd.handle);
    n_sdf_atoms = 0;
    attempts = 0;
    while ((attempts < od->object_num)
      && (n_sdf_atoms < od->al.mol_info[object_num]->n_atoms)
      && fgets(buffer, BUF_LEN, mol_fd.handle)) {
      buffer[BUF_LEN - 1] = '\0';
      parse_sdf_coord_line(od->al.mol_info[object_num]->sdf_version,
        buffer, NULL, sdf_coord, NULL);
      result = get_datafile_coord(od, &data_fd, n_sdf_atoms,
        od->al.mol_info[object_num]->n_atoms,
        &cube_word_size, data_coord, datafile_type);
      if (result && (result != OBJECTS_NOT_MATCHING)) {
        fclose(data_fd.handle);
        fclose(mol_fd.handle);
        return result;
      }
      if ((result == OBJECTS_NOT_MATCHING)
        || (fabs(data_coord[0] - sdf_coord[0]) > 0.0005)
        || (fabs(data_coord[1] - sdf_coord[1]) > 0.0005)
        || (fabs(data_coord[2] - sdf_coord[2]) > 0.0005)) {
        fclose(data_fd.handle);
        data_fd.handle = NULL;
        n_sdf_atoms = 0;
        ++attempts;
        if (fseek(mol_fd.handle, sdf_pos, SEEK_SET)) {
          O3_ERROR_LOCATE(&(od->task));
          O3_ERROR_STRING(&(od->task), data_fd.name);
          fclose(mol_fd.handle);
          return PREMATURE_DAT_EOF;
        }
        continue;
      }
      ++n_sdf_atoms;
    }
    if (data_fd.handle) {
      fclose(data_fd.handle);
      data_fd.handle = NULL;
    }
    fclose(mol_fd.handle);
    if (!n_sdf_atoms) {
      found = 0;
      rewind(od->file[TEMP_OBJECT_MATCH]->handle);
      while ((!found) && fgets(data_fd.name,
        BUF_LEN, od->file[TEMP_OBJECT_MATCH]->handle)) {
        remove_newline(data_fd.name);
        strcpy(buffer, data_fd.name);
        remove_extension(buffer);
        datafile_base = get_basename(buffer);
        found = (!strcasecmp(datafile_base, od->al.mol_info[object_num]->object_name));
      }
      if (!found) {
        od->task.data[DATA_OBJECT_NUM] = object_num;
        O3_ERROR_LOCATE(&(od->task));
        return OBJECTS_NOT_MATCHING;
      }
      n_sdf_atoms = od->al.mol_info[object_num]->n_atoms;
    }
    if (n_sdf_atoms < od->al.mol_info[object_num]->n_atoms) {
      O3_ERROR_LOCATE(&(od->task));
      O3_ERROR_STRING(&(od->task), mol_fd.name);
      return CANNOT_READ_ORIGINAL_SDF;
    }
    fprintf(od->file[TEMP_SORTED_MATCH]->handle, "%s\n", data_fd.name);
  }
  
  return 0;
}


int check_file_pattern(FILE *handle, char *file_pattern, int object_num)
{
  char dir_name[BUF_LEN];
  char *base_pattern;
  char *nextfile;
  int n = 0;
  #ifndef WIN32
  int match;
  DIR *dir = NULL;
  struct dirent dir_entry;
  struct dirent *result;
  #else
  char buffer[BUF_LEN];
  BOOL match;
  HANDLE dir;
  WIN32_FIND_DATA filedata;
  #endif
  
  
  memset(dir_name, 0, BUF_LEN);
  strncpy(dir_name, file_pattern, BUF_LEN - 1);
  dir_name[BUF_LEN - 1] = '\0';
  get_dirname(dir_name);
  base_pattern = get_basename(file_pattern);
  if (!dexist(dir_name)) {
    return CANNOT_OPEN_DIRECTORY;
  }
  
  #ifndef WIN32
  dir = opendir(dir_name);
  if (dir) {
    while (!readdir_r(dir, &dir_entry, &result)) {
      if (!result) {
        break;
      }
      nextfile = dir_entry.d_name;
  #else
  sprintf(buffer, "%s\\*", dir_name);
  dir = FindFirstFileA(buffer, &filedata);
  if (dir != INVALID_HANDLE_VALUE) {
    while (FindNextFileA(dir, &filedata)) {
      nextfile = filedata.cFileName;
  #endif
      if ((!strcmp(nextfile, "."))
        || (!strcmp(nextfile, ".."))) {
        continue;
      }
      #ifndef WIN32
      match = (fnmatch((const char *)base_pattern, (const char *)nextfile, 0) ? 0 : 1);
      #else
      match = PathMatchSpecA(nextfile, base_pattern);
      #endif
      if (match) {
        ++n;
        fprintf(handle, "%s%c%s\n", dir_name, SEPARATOR, nextfile);
      }
    }
  }
  #ifndef WIN32
  closedir(dir);
  #else
  CloseHandle(dir);
  #endif
  
  return ((n == object_num) ? 0 : ((n < object_num)
    ? NOT_ENOUGH_OBJECTS : TOO_MANY_OBJECTS));
}
