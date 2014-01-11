/*

parse_sdf.c

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
#include <include/basis_set.h>
#include <include/proc_env.h>
#include <include/ff_parm.h>
#include <include/prog_exe_info.h>
#include <include/extern.h>


extern EnvList babel_env[];
extern EnvList turbomole_env[];

char *get_y_var_name(char *buffer, char *y_name)
{
  char *y_name_start = NULL;
  char *y_name_end = NULL;
  int i;
  int found;


  if ((found = (buffer[0] == '>'))) {
    memset(y_name, 0, MAX_NAME_LEN);
    if ((y_name_start = strchr(buffer, '<'))
      && (y_name_end = strchr(y_name_start, '>'))) {
      ++y_name_start;
      --y_name_end;
      i = 0;
      while ((&y_name_start[i] <= y_name_end)
        && (i < (MAX_NAME_LEN - 1))) {
        y_name[i] = y_name_start[i];
        ++i;
      }
    }
  }
  
  return (((!y_name[0]) || (!found)) ? NULL : y_name);
}


void parse_sdf_coord_line(int sdf_version, char *buffer, char *element, double *coord, int *sdf_charge)
{
  char b;
  char *ptr = NULL;
  int i;
  
  
  if (sdf_version == V2000) {
    for (i = 0; i < 3; ++i) {
      b = buffer[10 * (i + 1)];
      buffer[10 * (i + 1)] = '\0';
      sscanf(&buffer[10 * i], "%lf", &coord[i]);
      buffer[10 * (i + 1)] = b;
    }
    if (element) {
      sscanf(&buffer[31], "%s", element);
    }
  }
  else {
    /*
    get formal charges (V3000 only)
    */
    if (element) {
      sscanf(buffer, "%*s %*s %*s %s %lf %lf %lf", element,
        &coord[0], &coord[1], &coord[2]);
    }
    else {
      sscanf(buffer, "%*s %*s %*s %*s %lf %lf %lf",
        &coord[0], &coord[1], &coord[2]);
    }
    if (sdf_charge && (ptr = strstr(buffer, "CHG="))) {
      sscanf(&ptr[4], "%d", sdf_charge);
    }
  }
}


int get_n_atoms_bonds(MolInfo *mol_info, FILE *handle, char *buffer)
{
  char buffer2[BUF_LEN];
  char b;
  char *counts = NULL;
  int pos;
  
  
  if ((strlen(buffer) <= 34) || (!strncasecmp(&buffer[34], "V2000", 5))) {
    mol_info->sdf_version = V2000;
    /*
    get number of atoms
    */
    b = buffer[3];
    buffer[3] = '\0';
    sscanf(buffer, "%d", &(mol_info->n_atoms));
    buffer[3] = b;
    /*
    get number of bonds
    */
    b = buffer[6];
    buffer[6] = '\0';
    sscanf(&buffer[3], "%d", &(mol_info->n_bonds));
    buffer[6] = b;
  }
  else if ((strlen(buffer) > 34) && (!strncasecmp(&buffer[34], "V3000", 5))) {
    mol_info->sdf_version = V3000;
    /*
    get number of atoms and bonds
    */
    counts = NULL;
    pos = ftell(handle);
    while ((!counts) && fgets(buffer2, BUF_LEN, handle)) {
      buffer2[BUF_LEN - 1] = '\0';
      counts = strstr(buffer2, "COUNTS");
    }
    if (!counts) {
      return PREMATURE_EOF;
    }
    if (fseek(handle, pos, SEEK_SET)) {
      return PREMATURE_EOF;
    }  
    sscanf(counts, "%*s %d %d",
      &(mol_info->n_atoms), &(mol_info->n_bonds));
  }
  else {
    return PREMATURE_EOF;
  }
  
  return 0;
}


int break_sdf_to_mol(O3Data *od, TaskInfo *task, FileDescriptor *from_fd, char *to_dir)
{
  char buffer[BUF_LEN];
  int object_num;
  int found;
  FileDescriptor mol_fd;


  memset(&mol_fd, 0, sizeof(FileDescriptor));
  memset(buffer, 0, BUF_LEN);
  if (!(from_fd->handle = fopen(from_fd->name, "rb"))) {
    O3_ERROR_LOCATE(task);
    O3_ERROR_STRING(task, from_fd->name);
    return FL_CANNOT_READ_SDF_FILE;
  }
  for (object_num = 0; object_num < od->grid.object_num; ++object_num) {
    sprintf(mol_fd.name, "%s%c%04d.mol", to_dir, SEPARATOR,
      od->al.mol_info[object_num]->object_id);
    if (!(mol_fd.handle = fopen(mol_fd.name, "wb"))) {
      O3_ERROR_LOCATE(task);
      O3_ERROR_STRING(task, mol_fd.name);
      fclose(from_fd->handle);
      return FL_CANNOT_WRITE_TEMP_FILE;
    }
    found = 0;
    while ((!found) && fgets(buffer, BUF_LEN, from_fd->handle)) {
      buffer[BUF_LEN - 1] = '\0';
      remove_newline(buffer);
      fprintf(mol_fd.handle, "%s\n", buffer);
      found = (!strncmp(buffer, MOL_DELIMITER, strlen(MOL_DELIMITER)));
    }
    if (!found) {
      O3_ERROR_LOCATE(task);
      O3_ERROR_STRING(task, from_fd->name);
      fclose(from_fd->handle);
      fclose(mol_fd.handle);
      return FL_CANNOT_READ_SDF_FILE;
    }
    found = 0;
    while ((!found) && fgets(buffer, BUF_LEN, from_fd->handle)) {
      buffer[BUF_LEN - 1] = '\0';
      found = (!strncmp(buffer, SDF_DELIMITER, 4));
    }
    fclose(mol_fd.handle);
    if (!found) {
      O3_ERROR_LOCATE(task);
      O3_ERROR_STRING(task, from_fd->name);
      fclose(from_fd->handle);
      return FL_CANNOT_READ_SDF_FILE;
    }
  }
  fclose(from_fd->handle);
  
  return 0;
}


int break_sdf_to_sdf(O3Data *od, TaskInfo *task, FileDescriptor *from_fd, char *to_dir)
{
  char buffer[BUF_LEN];
  int object_num;
  int found;
  FileDescriptor sdf_fd;


  memset(&sdf_fd, 0, sizeof(FileDescriptor));
  memset(buffer, 0, BUF_LEN);
  if (!(from_fd->handle = fopen(from_fd->name, "rb"))) {
    O3_ERROR_LOCATE(task);
    O3_ERROR_STRING(task, from_fd->name);
    return FL_CANNOT_READ_SDF_FILE;
  }
  for (object_num = 0; object_num < od->grid.object_num; ++object_num) {
    if (od->al.mol_info[object_num]->done & OBJECT_ALREADY_DONE) {
      continue;
    }
    sprintf(sdf_fd.name, "%s%c%04d.sdf", to_dir, SEPARATOR,
      od->al.mol_info[object_num]->object_id);
    if (!(sdf_fd.handle = fopen(sdf_fd.name, "wb"))) {
      O3_ERROR_LOCATE(task);
      O3_ERROR_STRING(task, sdf_fd.name);
      fclose(from_fd->handle);
      return FL_CANNOT_WRITE_TEMP_FILE;
    }
    found = 0;
    while ((!found) && fgets(buffer, BUF_LEN, from_fd->handle)) {
      buffer[BUF_LEN - 1] = '\0';
      remove_newline(buffer);
      fprintf(sdf_fd.handle, "%s\n", buffer);
      found = (!strncmp(buffer, SDF_DELIMITER, strlen(SDF_DELIMITER)));
    }
    fclose(sdf_fd.handle);
    if (!found) {
      O3_ERROR_LOCATE(task);
      O3_ERROR_STRING(task, from_fd->name);
      fclose(from_fd->handle);
      return FL_CANNOT_READ_SDF_FILE;
    }
  }
  fclose(from_fd->handle);
  
  return 0;
}


int join_mol_to_sdf(O3Data *od, TaskInfo *task, FileDescriptor *to_fd, char *from_dir)
{
  char buffer[BUF_LEN];
  int object_num;
  FileDescriptor mol_fd;


  memset(&mol_fd, 0, sizeof(FileDescriptor));
  memset(buffer, 0, BUF_LEN);
  if (!(to_fd->handle = fopen(to_fd->name, "wb"))) {
    O3_ERROR_LOCATE(task);
    O3_ERROR_STRING(task, to_fd->name);
    return FL_CANNOT_WRITE_TEMP_FILE;
  }
  for (object_num = 0; object_num < od->grid.object_num; ++object_num) {
    if (od->al.mol_info[object_num]->done & OBJECT_ALREADY_DONE) {
      continue;
    }
    sprintf(mol_fd.name, "%s%c%04d.mol",
      from_dir, SEPARATOR, od->al.mol_info[object_num]->object_id);
    if (!(mol_fd.handle = fopen(mol_fd.name, "rb"))) {
      O3_ERROR_LOCATE(task);
      O3_ERROR_STRING(task, mol_fd.name);
      fclose(to_fd->handle);
      return FL_CANNOT_READ_SDF_FILE;
    }
    while (fgets(buffer, BUF_LEN, mol_fd.handle)) {
      buffer[BUF_LEN - 1] = '\0';
      remove_newline(buffer);
      fprintf(to_fd->handle, "%s\n", buffer);
    }
    fprintf(to_fd->handle, SDF_DELIMITER"\n");
    fclose(mol_fd.handle);
  }
  fclose(to_fd->handle);
  
  return 0;
}


int convert_mol(O3Data *od, char *from_filename, char *to_filename, char *from_ext, char *to_ext, char *flags)
{
  char buffer[BUF_LEN];
  int result;
  int pid;
  ProgExeInfo prog_exe_info;
  
  
  memset(buffer, 0, BUF_LEN);
  memset(&prog_exe_info, 0, sizeof(ProgExeInfo));
  prog_exe_info.proc_env = fill_env
    (od, babel_env, od->field.babel_exe_path, 0);
  if (!(prog_exe_info.proc_env)) {
    return OUT_OF_MEMORY;
  }
  prog_exe_info.need_stdin = NEED_STDIN_NORMAL;
  prog_exe_info.stdout_fd = 0;
  prog_exe_info.stderr_fd = od->file[TEMP_LOG];
  prog_exe_info.exedir = od->field.babel_exe_path;
  prog_exe_info.sep_proc_grp = 1;
  sprintf(prog_exe_info.command_line,
    "%s%c"BABEL_EXE" %s -i%s %s -o%s %s",
    od->field.babel_exe_path, SEPARATOR, flags,
    from_ext, from_filename, to_ext, to_filename);
  pid = ext_program_exe(&prog_exe_info, &result);
  if (result) {
    return result;
  }
  ext_program_wait(&prog_exe_info, pid);
  free_proc_env(prog_exe_info.proc_env);
  if (!(od->file[TEMP_LOG]->handle = fopen
    (od->file[TEMP_LOG]->name, "rb"))) {
    return OPENBABEL_ERROR;
  }
  result = fgrep(od->file[TEMP_LOG]->handle, buffer, "error");
  fclose(od->file[TEMP_LOG]->handle);
  od->file[TEMP_LOG]->handle = NULL;
  if (result) {  
    return OPENBABEL_ERROR;
  }

  return 0;
}


int update_mol(O3Data *od)
{
  char buffer[BUF_LEN];
  int i;
  int n;
  int object_num;
  double atom_coord[3];
  FileDescriptor mol_fd;


  memset(buffer, 0, BUF_LEN);
  memset(&mol_fd, 0, sizeof(FileDescriptor));
  od->field.max_n_atoms = 0;
  od->field.max_n_bonds = 0;
  if ((od->pymol.use_pymol) || (od->jmol.use_jmol)) {
    if (!convert_mol(od, od->file[MOLFILE_IN]->name,
      od->file[TEMP_OUT]->name, "sdf", "mol2", "-xl")) {
      object_num = 0;
      if ((od->file[TEMP_OUT]->handle =
        fopen(od->file[TEMP_OUT]->name, "rb"))) {
        while (fgets(buffer, BUF_LEN, od->file[TEMP_OUT]->handle)) {
          buffer[BUF_LEN - 1] = '\0';
          remove_newline(buffer);
          if (strstr(buffer, "@<TRIPOS>MOLECULE")) {
            if (od->file[BINARY_IN]->handle) {
              fclose(od->file[BINARY_IN]->handle);
              od->file[BINARY_IN]->handle = NULL;
            }
            sprintf(od->file[BINARY_IN]->name, "%s%c%04d.mol2", od->field.mol_dir,
              SEPARATOR, od->al.mol_info[object_num]->object_id);
            od->file[BINARY_IN]->handle = fopen
              (od->file[BINARY_IN]->name, "wb");
            ++object_num;
          }
          if (od->file[BINARY_IN]->handle) {
            fprintf(od->file[BINARY_IN]->handle, "%s\n", buffer);
          }
        }
        if (od->file[BINARY_IN]->handle) {
          fclose(od->file[BINARY_IN]->handle);
          od->file[BINARY_IN]->handle = NULL;
        }
        fclose(od->file[TEMP_OUT]->handle);
        od->file[TEMP_OUT]->handle = NULL;
      }
    }
  }
  for (object_num = 0; object_num < od->grid.object_num; ++object_num) {
    sprintf(mol_fd.name, "%s%c%04d.mol", od->field.mol_dir,
      SEPARATOR, od->al.mol_info[object_num]->object_id);
    if (!(mol_fd.handle = fopen(mol_fd.name, "rb"))) {
      O3_ERROR_LOCATE(&(od->task));
      O3_ERROR_STRING(&(od->task), mol_fd.name);
      return CANNOT_READ_TEMP_FILE;
    }
    /*
    get number of atoms
    */
    n = 0;
    while ((n < 4) && fgets(buffer, BUF_LEN, mol_fd.handle)) {
      ++n;
    }
    if (n != 4) {
      fclose(mol_fd.handle);
      O3_ERROR_LOCATE(&(od->task));
      O3_ERROR_STRING(&(od->task), mol_fd.name);
      return CANNOT_READ_TEMP_FILE;
    }
    remove_newline(buffer);
    if (get_n_atoms_bonds(od->al.mol_info[object_num], mol_fd.handle, buffer)) {
      fclose(mol_fd.handle);
      O3_ERROR_LOCATE(&(od->task));
      O3_ERROR_STRING(&(od->task), mol_fd.name);
      return CANNOT_READ_TEMP_FILE;
    }
    if (od->al.mol_info[object_num]->n_atoms > od->field.max_n_atoms) {
      od->field.max_n_atoms = od->al.mol_info[object_num]->n_atoms;
    }
    if (od->al.mol_info[object_num]->n_bonds > od->field.max_n_bonds) {
      od->field.max_n_bonds = od->al.mol_info[object_num]->n_bonds;
    }
    if (od->al.mol_info[object_num]->sdf_version == V3000) {
      if (!fgrep(mol_fd.handle, buffer, "BEGIN ATOM")) {
        fclose(mol_fd.handle);
        O3_ERROR_LOCATE(&(od->task));
        O3_ERROR_STRING(&(od->task), mol_fd.name);
        return CANNOT_READ_TEMP_FILE;
      }
    }    
    n = 0;
    while ((n < od->al.mol_info[object_num]->n_atoms)
      && fgets(buffer, BUF_LEN, mol_fd.handle)) {
      buffer[BUF_LEN - 1] = '\0';
      parse_sdf_coord_line(od->al.mol_info[object_num]->sdf_version,
        buffer, NULL, atom_coord, NULL);
      for (i = 0; i < 3; ++i) {
        if (((!n) && (!object_num)) || (atom_coord[i]) < od->min_coord[i]) {
          od->min_coord[i] = atom_coord[i];
        }
        if (((!n) && (!object_num)) || (atom_coord[i]) > od->max_coord[i]) {
          od->max_coord[i] = atom_coord[i];
        }
      }
      ++n;
    }
    fclose(mol_fd.handle);
    if (n != od->al.mol_info[object_num]->n_atoms) {
      O3_ERROR_LOCATE(&(od->task));
      O3_ERROR_STRING(&(od->task), mol_fd.name);
      return CANNOT_READ_TEMP_FILE;
    }
    if (od->al.mol_info[object_num]->name_count) {
      continue;
    }
    for (i = object_num + 1, n = 1; i < od->grid.object_num; ++i) {
      if (!strcasecmp(od->al.mol_info[object_num]->object_name,
        od->al.mol_info[i]->object_name)) {
        od->al.mol_info[i]->name_count = n;
        ++n;
      }
    }
  }
  
  return 0;
}


int parse_sdf(O3Data *od, int options, char *name_list)
{
  char *ptr;
  char *context = NULL;
  char buffer[BUF_LEN];
  char y_name[MAX_NAME_LEN];
  int i;
  int n = 0;
  int read_y_value = 0;
  int found = 0;
  int y_vars = 0;
  int look_for_sdf_delimiter = 0;
  int line = 0;
  int molecule_num;
  int object_num = 0;
  double value;
  FileDescriptor mol_fd;


  memset(&mol_fd, 0, sizeof(FileDescriptor));
  memset(buffer, 0, BUF_LEN);
  molecule_num = 0;
  /*
  count how many molecules are there in the sdf file
  */
  while (fgets(buffer, BUF_LEN, od->file[MOLFILE_IN]->handle)) {
    buffer[BUF_LEN - 1] = '\0';
    if (strstr(buffer, SDF_DELIMITER)) {
      ++molecule_num;
    }
  }
  rewind(od->file[MOLFILE_IN]->handle);
  if (!molecule_num) {
    return PREMATURE_EOF;
  }
  /*
  fill grid info
  */
  od->grid.object_num = molecule_num;
  od->grid.struct_num = molecule_num;
  if (options & ALLOC_MOL_INFO_BIT) {
    /*
    alloc MolInfo structure array
    */
    free_array(od->al.mol_info);
    if (!(od->al.mol_info = (MolInfo **)
      alloc_array(molecule_num, sizeof(MolInfo)))) {
      return OUT_OF_MEMORY;
    }
    for (i = 0; i < molecule_num; ++i) {
      od->al.mol_info[i]->object_id = i + 1;
      od->al.mol_info[i]->struct_num = i;
      od->al.mol_info[i]->conf_num = 1;
    }
  }
  if (options & IMPORT_Y_VARS_BIT) {
    /*
    tokenize string using
    comma as separator
    */
    free_char_matrix(od->cimal.y_var_name);
    od->cimal.y_var_name = NULL;
    y_vars = 0;
    if (strcasecmp(name_list, "all")) {
      strcpy(buffer, name_list);
      ptr = strtok_r(buffer, ",\0", &context);
      while (ptr) {
        ++y_vars;
        ptr = strtok_r(NULL, ",\0", &context);
      }
      if (!(od->cimal.y_var_name = alloc_char_matrix
        (od->cimal.y_var_name, y_vars, MAX_NAME_LEN))) {
        return OUT_OF_MEMORY;
      }
      strcpy(buffer, name_list);
      ptr = strtok_r(buffer, ",\0", &context);
      y_vars = 0;
      while (ptr) {
        for (i = 0, found = 0; ((!found) && (i < y_vars)); ++i) {
          found = (!strcmp(od->cimal.y_var_name->me[i], ptr));
        }
        if (!found) {
          strcpy(od->cimal.y_var_name->me[y_vars], ptr);
          ++y_vars;
        }
        ptr = strtok_r(NULL, ",\0", &context);
      }
      if (!(od->pel.numberlist[Y_VAR_LIST] = int_perm_resize
        (od->pel.numberlist[Y_VAR_LIST], y_vars))) {
        return OUT_OF_MEMORY;
      }
      memset(od->pel.numberlist[Y_VAR_LIST]->pe, 0, y_vars * sizeof(int));
      while (fgets(buffer, BUF_LEN, od->file[MOLFILE_IN]->handle)) {
        buffer[BUF_LEN - 1] = '\0';
        if (get_y_var_name(buffer, y_name)) {
          for (i = 0; i < y_vars; ++i) {
            if ((!(od->pel.numberlist[Y_VAR_LIST]->pe[i]))
              && (!strcmp(od->cimal.y_var_name->me[i], y_name))) {
               od->pel.numberlist[Y_VAR_LIST]->pe[i] = 1;
            }
          }
        }
      }
      for (i = 0, found = 1; (found && (i < y_vars)); ++i) {
        found = od->pel.numberlist[Y_VAR_LIST]->pe[i];
      }
      if (!found) {
        return CANNOT_FIND_Y_VAR_NAME;
      }
    }
    else {
      while (fgets(buffer, BUF_LEN, od->file[MOLFILE_IN]->handle)) {
        buffer[BUF_LEN - 1] = '\0';
        if (get_y_var_name(buffer, y_name)) {
          for (i = 0, found = 0; ((!found) && (i < y_vars)); ++i) {
            found = (!strcmp(od->cimal.y_var_name->me[i], y_name));
          }
          if (!found) {
            ++y_vars;
            if (!(od->cimal.y_var_name = alloc_char_matrix
              (od->cimal.y_var_name, y_vars, MAX_NAME_LEN))) {
              return OUT_OF_MEMORY;
            }
            strcpy(od->cimal.y_var_name->me[y_vars - 1], y_name);
          }
        }
      }
      if (!(od->pel.numberlist[Y_VAR_LIST] = int_perm_resize
        (od->pel.numberlist[Y_VAR_LIST], y_vars))) {
        return OUT_OF_MEMORY;
      }
    }
    memset(od->pel.numberlist[Y_VAR_LIST]->pe, 0, y_vars * sizeof(int));
    rewind(od->file[MOLFILE_IN]->handle);
    od->y_vars = y_vars;
    if (alloc_y_var_array(od)) {
      return OUT_OF_MEMORY;
    }
  }
  while (fgets(buffer, BUF_LEN, od->file[MOLFILE_IN]->handle)
    && (object_num < molecule_num)) {
    buffer[BUF_LEN - 1] = '\0';
    remove_newline(buffer);
    if (look_for_sdf_delimiter) {
      if (mol_fd.handle) {
        fclose(mol_fd.handle);
        memset(&mol_fd, 0, sizeof(FileDescriptor));
      }
      if (read_y_value) {
        /*
        get y var value
        */
        read_y_value = 0;
        value = 0.0;
        sscanf(buffer, "%lf", &value);
        set_y_value(od, object_num, y_vars, value);
      }
      if ((options & IMPORT_Y_VARS_BIT) && (get_y_var_name(buffer, y_name))) {
        for (i = 0, read_y_value = 0; ((!read_y_value) && (i < od->y_vars)); ++i) {
          read_y_value = (!strcmp(od->cimal.y_var_name->me[i], y_name));
          if (read_y_value) {
            od->pel.numberlist[Y_VAR_LIST]->pe[i] = 1;
            y_vars = i;
            ++n;
          }
        }
      }
      else if (!strncmp(buffer, SDF_DELIMITER, 4)) {
        look_for_sdf_delimiter = 0;
        line = 0;
        if (options & IMPORT_Y_VARS_BIT) {
          if (od->y_vars - n) {
            tee_printf(od, "Object %4d: undefined variable%s ",
              object_num + 1, (((od->y_vars - n) > 1) ? "s" : ""));
            for (i = 0; i < od->y_vars; ++i) {
              if (!(od->pel.numberlist[Y_VAR_LIST]->pe[i])) {
                tee_printf(od, "%s%s", od->cimal.y_var_name->me[i],
                  ((i == (od->y_vars - 1)) ? "" : ", "));
                set_y_value(od, object_num, i, 0.0);
              }
            }
            tee_printf(od, "; setting to 0.0\n");
          }
          /*
          initialize found variable vector and missing
          variable count for the current object
          */
          y_vars = 0;
          n = 0;
          memset(od->pel.numberlist[Y_VAR_LIST]->pe, 0, od->y_vars * sizeof(int));
        }
        ++object_num;
      }
      continue;
    }
    ++line;
    if (line == 1) {
      sprintf(mol_fd.name, "%s%c%04d.mol",
        od->field.mol_dir, SEPARATOR, od->al.mol_info[object_num]->object_id);
      if (!(mol_fd.handle = fopen(mol_fd.name, "wb+"))) {
        O3_ERROR_LOCATE(&(od->task));
        O3_ERROR_STRING(&(od->task), mol_fd.name);
        return CANNOT_WRITE_TEMP_FILE;
      }
      /*
      get molecule name
      */
      i = 0;
      while ((i < MAX_NAME_LEN) && buffer[i]) {
        /*
        replace eventual special chars with underscores
        */
        if (strchr(" \t,:;()[]{}\"!$%&/\\|*+=?'^<>", (int)buffer[i])) {
          buffer[i] = '_';
        }
        /*
        when linefeed is met, end-of-string
        */
        else if (strchr("\r\n", (int)buffer[i])) {
          buffer[i] = '\0';
          break;
        }
        ++i;
      }
      buffer[MAX_NAME_LEN - 1] = '\0';
      if (!buffer[0]) {
        strcpy(buffer, "NONAME");
      }
      strcpy(od->al.mol_info[object_num]->object_name, buffer);
    }
    if (line == 4) {
      if (get_n_atoms_bonds(od->al.mol_info[object_num],
        od->file[MOLFILE_IN]->handle, buffer)) {
        return PREMATURE_EOF;
      }
    }
    fprintf(mol_fd.handle, "%s\n", buffer);
    look_for_sdf_delimiter = (!strncmp(buffer, MOL_DELIMITER, strlen(MOL_DELIMITER)));
  }
  if (options & IMPORT_Y_VARS_BIT) {
    tee_printf(od, "\n");
  }
  if (object_num != molecule_num) {
    return PREMATURE_EOF;
  }
  
  return 0;
}


int call_obenergy(O3Data *od, int force_field)
{
  char buffer[BUF_LEN];
  int error;
  int pid;
  int object_num;
  ProgExeInfo prog_exe_info;


  memset(buffer, 0, BUF_LEN);
  for (object_num = 0; object_num < od->grid.object_num; ++object_num) {
    sprintf(buffer, "%s%c%04d.%s", od->field.mol_dir, SEPARATOR,
      od->al.mol_info[object_num]->object_id, "mmff");
    if (!fexist(buffer)) {
      break;
    }
  }
  /*
  obenergy has already been called previously
  */
  if (object_num == od->grid.object_num) {
    return 0;
  }
  /*
  Call obenergy (OpenBabel) to obtain force-field-specific
  atom types and charges, as well as formal charges
  (the latter only for MMFF94)
  */
  memset(&prog_exe_info, 0, sizeof(ProgExeInfo));
  prog_exe_info.proc_env = fill_env
    (od, babel_env, od->field.babel_exe_path, 0);
  if (!(prog_exe_info.proc_env)) {
    return OUT_OF_MEMORY;
  }
  prog_exe_info.need_stdin = 0;
  prog_exe_info.stdout_fd = od->file[TEMP_OUT];
  prog_exe_info.stderr_fd = od->file[TEMP_LOG];
  prog_exe_info.exedir = od->field.babel_exe_path;
  prog_exe_info.sep_proc_grp = 1;
  sprintf(prog_exe_info.command_line,
    "%s%c"OBENERGY_EXE" -ff %s %s",
    od->field.babel_exe_path, SEPARATOR, "MMFF94",
    od->file[MOLFILE_IN]->name);
  pid = ext_program_exe(&prog_exe_info, &error);
  if (error) {
    return error;
  }
  ext_program_wait(&prog_exe_info, pid);
  free_proc_env(prog_exe_info.proc_env);
  if (!(od->file[TEMP_LOG]->handle = fopen
    (od->file[TEMP_LOG]->name, "rb"))) {
    return OPENBABEL_ERROR;
  }
  fclose(od->file[TEMP_LOG]->handle);
  od->file[TEMP_LOG]->handle = NULL;
  /*
  split obenergy output into multiple files
  */
  if (!(od->file[TEMP_OUT]->handle = fopen
    (od->file[TEMP_OUT]->name, "rb"))) {
    return OPENBABEL_ERROR;
  }
  object_num = 0;
  while (fgets(buffer, BUF_LEN, od->file[TEMP_OUT]->handle)) {
    buffer[BUF_LEN - 1] = '\0';
    remove_newline(buffer);
    if (strstr(buffer, "A T O M   T Y P E S")) {
      if (od->file[BINARY_IN]->handle) {
        fclose(od->file[BINARY_IN]->handle);
        od->file[BINARY_IN]->handle = NULL;
      }
      sprintf(od->file[BINARY_IN]->name, "%s%c%04d.%s", od->field.mol_dir,
        SEPARATOR, od->al.mol_info[object_num]->object_id, "mmff");
      if (!(od->file[BINARY_IN]->handle = fopen
        (od->file[BINARY_IN]->name, "wb"))) {
        return OPENBABEL_ERROR;
      }
      ++object_num;
    }
    if (od->file[BINARY_IN]->handle) {
      fprintf(od->file[BINARY_IN]->handle, "%s\n", buffer);
    }
  }
  if (od->file[BINARY_IN]->handle) {
    fclose(od->file[BINARY_IN]->handle);
    od->file[BINARY_IN]->handle = NULL;
  }
  fclose(od->file[TEMP_OUT]->handle);
  od->file[TEMP_OUT]->handle = NULL;
  if (object_num != od->grid.object_num) {
    return OPENBABEL_ERROR;
  }
  
  return 0;
}


int fill_atom_info(O3Data *od, TaskInfo *task, AtomInfo **atom, BondList **bond_list, int object_num, char force_field)
{
  char b = 0;
  char buffer[BUF_LEN];
  char ring_type[MAX_FF_TYPE_LEN];
  char *ptr;
  char *context = NULL;
  int found = 0;
  int unknown;
  int i;
  int j;
  int n;
  int n_atoms;
  int n_bonds;
  int n_charges;
  int sdf_charge;
  int a[2];
  int order;
  int result = 0;
  FileDescriptor mol_fd;
  FileDescriptor out_fd;


  memset(&mol_fd, 0, sizeof(FileDescriptor));
  memset(&out_fd, 0, sizeof(FileDescriptor));
  sprintf(mol_fd.name, "%s%c%04d.mol", od->field.mol_dir,
    SEPARATOR, od->al.mol_info[object_num]->object_id);
  sprintf(out_fd.name, "%s%c%04d.%s", od->field.mol_dir,
    SEPARATOR, od->al.mol_info[object_num]->object_id, "mmff");
  /*
  open MOL file and get number of atoms, bonds
  */
  if (!(mol_fd.handle = fopen(mol_fd.name, "rb"))) {
    O3_ERROR_LOCATE(task);
    O3_ERROR_STRING(task, mol_fd.name);
    return FL_CANNOT_READ_MOL_FILE;
  }
  n = 0;
  while ((n < 4) && fgets(buffer, BUF_LEN, mol_fd.handle)) {
    ++n;
  }
  if (n != 4) {
    fclose(mol_fd.handle);
    O3_ERROR_LOCATE(task);
    O3_ERROR_STRING(task, mol_fd.name);
    return FL_CANNOT_READ_MOL_FILE;
  }
  if (od->al.mol_info[object_num]->sdf_version == V3000) {
    if (!fgrep(mol_fd.handle, buffer, "BEGIN ATOM")) {
      fclose(mol_fd.handle);
      O3_ERROR_LOCATE(task);
      O3_ERROR_STRING(task, mol_fd.name);
      return FL_CANNOT_READ_MOL_FILE;
    }
  }    
  n_atoms = od->al.mol_info[object_num]->n_atoms;
  n_bonds = od->al.mol_info[object_num]->n_bonds;
  for (i = 0; i < n_atoms; ++i) {
    memset(atom[i], 0, sizeof(AtomInfo));
  }
  for (i = 0; bond_list && (i < n_bonds); ++i) {
    memset(bond_list[i], 0, sizeof(BondList));
  }
  if (!(out_fd.handle = fopen(out_fd.name, "rb"))) {
    fclose(mol_fd.handle);
    O3_ERROR_LOCATE(task);
    O3_ERROR_STRING(task, out_fd.name);
    return FL_CANNOT_READ_OB_OUTPUT;
  }
  if (!fgrep(out_fd.handle, buffer,
    "F O R M A L   C H A R G E S")) {
    fclose(out_fd.handle);
    fclose(mol_fd.handle);
    O3_ERROR_LOCATE(task);
    O3_ERROR_STRING(task, out_fd.name);
    return FL_CANNOT_READ_OB_OUTPUT;
  }
  found = 0;
  while ((!found) && fgets(buffer, BUF_LEN, out_fd.handle)) {
    buffer[BUF_LEN - 1] = '\0';
    found = (strstr(buffer, "CHARGE") ? 1 : 0);
  }
  if (!found) {
    fclose(out_fd.handle);
    fclose(mol_fd.handle);
    O3_ERROR_LOCATE(task);
    O3_ERROR_STRING(task, out_fd.name);
    return FL_CANNOT_READ_OB_OUTPUT;
  }
  found = 0;
  n = 0;
  while ((n < n_atoms) && fgets(buffer, BUF_LEN, out_fd.handle)) {
    buffer[BUF_LEN - 1] = '\0';
    if (isspace((int)buffer[0])) {
      break;
    }
    sscanf(buffer, "%*s %lf", &(atom[n]->formal_charge));
    ++n;
  }
  if (n != n_atoms) {
    fclose(out_fd.handle);
    fclose(mol_fd.handle);
    O3_ERROR_LOCATE(task);
    O3_ERROR_STRING(task, out_fd.name);
    return FL_CANNOT_READ_OB_OUTPUT;
  }
  /*
  find atom types for this molecule
  */
  if (!fgrep(out_fd.handle, buffer, "A T O M   T Y P E S")) {
    fclose(out_fd.handle);
    fclose(mol_fd.handle);
    O3_ERROR_LOCATE(task);
    O3_ERROR_STRING(task, out_fd.name);
    return FL_CANNOT_READ_OB_OUTPUT;
  }
  found = 0;
  while ((!found) && fgets(buffer, BUF_LEN, out_fd.handle)) {
    buffer[BUF_LEN - 1] = '\0';
    found = (strstr(buffer, "TYPE") ? 1 : 0);
  }
  if (!found) {
    fclose(out_fd.handle);
    fclose(mol_fd.handle);
    O3_ERROR_LOCATE(task);
    O3_ERROR_STRING(task, out_fd.name);
    return FL_CANNOT_READ_OB_OUTPUT;
  }
  found = 0;
  n = 0;
  while ((n < n_atoms) && fgets(buffer, BUF_LEN, out_fd.handle)) {
    buffer[BUF_LEN - 1] = '\0';
    if (isspace((int)buffer[0])) {
      break;
    }
    unknown = 1;
    memset(ring_type, 0, MAX_FF_TYPE_LEN);
    sscanf(buffer, "%*s %d %s", &j, ring_type);
    i = 0;
    while (ff_parm[O3_MMFF94][i].type_num
      && (unknown = (j != ff_parm[O3_MMFF94][i].type_num))) {
      ++i;
    }
    if (unknown) {
      --i;
      result = FL_UNKNOWN_ATOM_TYPE;
    }
    else {
      atom[n]->atom_type = ff_parm[(int)force_field][i].type_num;
      atom[n]->atom_num = n + 1;
      switch (ring_type[1]) {
        case 'R':
        atom[n]->ring = RING_BIT | AROMATIC;
        break;
        
        case 'L':
        atom[n]->ring = RING_BIT;
        break;
        
        default:
        atom[n]->ring = 0;
        break;
      }
      strcpy(atom[n]->atom_name, ff_parm[(int)force_field][i].type_chr);
    }
    ++n;
  }
  /*
  this is to indicate that the atom list is over
  */
  atom[n]->atom_type = -1;
  if (n != n_atoms) {
    fclose(out_fd.handle);
    fclose(mol_fd.handle);
    O3_ERROR_LOCATE(task);
    O3_ERROR_STRING(task, out_fd.name);
    return FL_CANNOT_READ_OB_OUTPUT;
  }
  /*
  find partial charges for this molecule
  */
  if (!fgrep(out_fd.handle, buffer, "P A R T I A L   C H A R G E S")) {
    fclose(out_fd.handle);
    fclose(mol_fd.handle);
    O3_ERROR_LOCATE(task);
    O3_ERROR_STRING(task, out_fd.name);
    return FL_CANNOT_READ_OB_OUTPUT;
  }
  found = 0;
  while ((!found) && fgets(buffer, BUF_LEN, out_fd.handle)) {
    buffer[BUF_LEN - 1] = '\0';
    found = (strstr(buffer, "CHARGE") ? 1 : 0);
  }
  if (!found) {
    fclose(out_fd.handle);
    fclose(mol_fd.handle);
    O3_ERROR_LOCATE(task);
    O3_ERROR_STRING(task, out_fd.name);
    return FL_CANNOT_READ_OB_OUTPUT;
  }
  found = 0;
  n = 0;
  while ((n < n_atoms) && fgets(buffer, BUF_LEN, out_fd.handle)) {
    buffer[BUF_LEN - 1] = '\0';
    if (isspace((int)buffer[0])) {
      break;
    }
    sscanf(buffer, "%*s %lf", &(atom[n]->charge));
    ++n;
  }
  if (n != n_atoms) {
    fclose(out_fd.handle);
    fclose(mol_fd.handle);
    O3_ERROR_LOCATE(task);
    O3_ERROR_STRING(task, out_fd.name);
    return FL_CANNOT_READ_OB_OUTPUT;
  }
  fclose(out_fd.handle);
  out_fd.handle = NULL;
  /*
  read atom coordinates from the mol file
  */
  n = 0;
  od->al.mol_info[object_num]->n_heavy_atoms = 0;
  while ((n < n_atoms) && fgets(buffer, BUF_LEN, mol_fd.handle)) {
    buffer[BUF_LEN - 1] = '\0';
    parse_sdf_coord_line(od->al.mol_info[object_num]->sdf_version,
      buffer, atom[n]->element, atom[n]->coord, &(atom[n]->sdf_charge));
    if (strcmp(atom[n]->element, "H")) {
      ++(od->al.mol_info[object_num]->n_heavy_atoms);
    }
    ++n;
  }
  if (n != n_atoms) {
    fclose(mol_fd.handle);
    O3_ERROR_LOCATE(task);
    O3_ERROR_STRING(task, mol_fd.name);
    return FL_CANNOT_READ_MOL_FILE;
  }
  /*
  build connectivity table
  */
  n = 0;
  if (od->al.mol_info[object_num]->sdf_version == V3000) {
    if (!fgrep(mol_fd.handle, buffer, "BEGIN BOND")) {
      fclose(mol_fd.handle);
      O3_ERROR_LOCATE(task);
      O3_ERROR_STRING(task, mol_fd.name);
      return FL_CANNOT_READ_MOL_FILE;
    }
  }    
  while ((n < n_bonds) && fgets(buffer, BUF_LEN, mol_fd.handle)) {
    buffer[BUF_LEN - 1] = '\0';
    if (od->al.mol_info[object_num]->sdf_version == V2000) {
      if (!strncmp(buffer, MOL_DELIMITER, strlen(MOL_DELIMITER))) {
        break;
      }
      for (i = 0; i <= 1; ++i) {
        b = buffer[3 * (i + 1)];
        buffer[3 * (i + 1)] = '\0';
        sscanf(&buffer[3 * i], "%d", &a[i]);
        buffer[3 * (i + 1)] = b;
      }
      sscanf(&buffer[8], "%d", &order);
    }
    else {
      sscanf(buffer, "%*s %*s %*s %d %d %d", &order, &a[0], &a[1]);
    }
    for (i = 0; i <= 1; ++i) {
      if ((a[i] < 1) || (a[i] > n_atoms)) {
        fclose(mol_fd.handle);
        O3_ERROR_LOCATE(task);
        O3_ERROR_STRING(task, mol_fd.name);
        return FL_CANNOT_READ_MOL_FILE;
      }
      atom[a[i] - 1]->bonded[atom[a[i] - 1]->n_bonded].num = a[1 - i] - 1;
      atom[a[i] - 1]->bonded[atom[a[i] - 1]->n_bonded].order = order;
      if (bond_list) {
        bond_list[n]->a[i] = a[i] - 1;
        if (i) {
          bond_list[n]->order = (((atom[a[0] - 1]->ring & AROMATIC)
            && (atom[a[1] - 1]->ring & AROMATIC)) ? AROMATIC : order);
        }
      }
      ++(atom[a[i] - 1]->n_bonded);
    }
    ++n;
  }
  if (n != n_bonds) {
    fclose(mol_fd.handle);
    O3_ERROR_LOCATE(task);
    O3_ERROR_STRING(task, mol_fd.name);
    return FL_CANNOT_READ_MOL_FILE;
  }
  /*
  get formal charges (V2000 only)
  */
  if (od->al.mol_info[object_num]->sdf_version == V2000) {
    while (fgets(buffer, BUF_LEN, mol_fd.handle)) {
      buffer[BUF_LEN - 1] = '\0';
      if (!strncmp(buffer, MOL_DELIMITER, strlen(MOL_DELIMITER))) {
        break;
      }
      if (strstr(buffer, "M  CHG")) {
        ptr = strtok_r(buffer, " \t\n\r\0", &context);
        j = 2;
        while (ptr && j) {
          --j;
          ptr = strtok_r(NULL, " \t\n\r\0", &context);
        }
        /*
        read how many charges are indicated
        on this line
        */
        if (!ptr) {
          continue;
        }
        sscanf(ptr, "%d", &n_charges);
        j = 1;
        while (ptr && j) {
          --j;
          ptr = strtok_r(NULL, " \t\n\r\0", &context);
        }
        for (i = 0; ptr && (i < n_charges); ++i) {
          sscanf(ptr, "%d", &n);
          j = 1;
          while (ptr && j) {
            --j;
            ptr = strtok_r(NULL, " \t\n\r\0", &context);
          }
          sscanf(ptr, "%d", &sdf_charge);
          if ((n > 0) && (n <= n_atoms)) {
            atom[n - 1]->sdf_charge = sdf_charge;
          }
          j = 1;
          while (ptr && j) {
            --j;
            ptr = strtok_r(NULL, " \t\n\r\0", &context);
          }
        }
      }
    }
  }
  for (n = 0; n < n_atoms; ++n) {
    if (atom[n]->n_bonded) {
      qsort(&(atom[n]->bonded), atom[n]->n_bonded,
        sizeof(BondInfo), compare_bond_info);
    }
  }
  fclose(mol_fd.handle);
  
  return result;
}


int fill_md_grid_types(AtomInfo **atom)
{
  char *univalent_cations[] = {
    "B", "Cs", "Ga", "K", "Li", "Na", "Rb", "Ti", NULL };
  char *divalent_cations[] = {
    "Ba", "Be", "Ca", "Co", "Cr", "Cu", "Fe", "Mg",
    "Mn", "Ni", "Pb", "Sn", "Sr", "Ti", "V", "Zn", NULL };
  char *trivalent_cations[] = {
    "Co", "Cr", "Fe", "Ga", "In", "Mn", "Sc", "Ti",
    "V", "Y", NULL };
  char *tetravalent_cations[] = {
    "Cr", "Ge", "Mn", "Mo", "Nb", "Pb", "Sn", "Ti",
    "V", "W", "Zr", NULL };
  int i;
  int j;
  int k;
  int n_alkyl;
  int n_aryl_vinyl;
  int n_carb_sulf_phos;
  int n_atoms;
  int n_hydro;
  int n_nitro;
  int n_oxy;
  int n_sulf;
  int n_fluorine;
  int n_free;
  int oxy_bound_to_hydro;
  int error = 0;
  int sdf_charge;
  double formal_charge;
  double dist;
  AtomInfo *carbon;
  AtomInfo *neighbor;
  AtomInfo *oxy;
  CationList valence[4];    


  valence[0].cations = univalent_cations;
  valence[1].cations = divalent_cations;
  valence[2].cations = trivalent_cations;
  valence[3].cations = tetravalent_cations;
  n_atoms = 0;
  while (atom[n_atoms]->atom_type != -1) {
    ++n_atoms;
  }
  for (i = 0; i < n_atoms; ++i) {
    strcpy(atom[i]->atom_name, "UNK");
  }
  for (i = 0; i < n_atoms; ++i) {
    for (j = 0, n_hydro = 0, n_nitro = 0, n_oxy = 0, n_sulf = 0;
      j < atom[i]->n_bonded; ++j) {
      if (!strcmp(atom[atom[i]->bonded[j].num]->element, "H")) {
        ++n_hydro;
      }
      else if (!strcmp(atom[atom[i]->bonded[j].num]->element, "N")) {
        ++n_nitro;
      }
      else if (!strcmp(atom[atom[i]->bonded[j].num]->element, "O")) {
        ++n_oxy;
      }
      else if (!strcmp(atom[atom[i]->bonded[j].num]->element, "S")) {
        ++n_sulf;
      }
    }
    if ((!strcmp(atom[i]->element, "C")) && (!strcmp(atom[i]->atom_name, "UNK"))) {
      /*
      CARBON
      */
      switch (n_hydro) {
        case 3:
        if (atom[i]->atom_type == 1) {
          strcpy(atom[i]->atom_name, " C3");
        }
        break;

        case 2:
        if ((atom[i]->atom_type == 1) || (atom[i]->atom_type == 20)
          || (atom[i]->atom_type == 22)) {
          strcpy(atom[i]->atom_name, " C2");
        }
        else if ((atom[i]->atom_type == 2) || (atom[i]->atom_type == 3)) {
          strcpy(atom[i]->atom_name, " C2=");
        }
        break;

        case 1:
        if ((atom[i]->atom_type == 1) || (atom[i]->atom_type == 20)
          || (atom[i]->atom_type == 22)) {
          strcpy(atom[i]->atom_name, " C1");
        }
        else if ((atom[i]->atom_type == 37) || (atom[i]->atom_type == 2)
          || (atom[i]->atom_type == 63) || (atom[i]->atom_type == 64)
          || (atom[i]->atom_type == 78) || (atom[i]->atom_type == 30)
          || (atom[i]->atom_type == 80)) {
          strcpy(atom[i]->atom_name, " C1=");
        }
        else if (atom[i]->atom_type == 3) {
          if (atom[i]->n_bonded == 3) {
            for (j = 0; j < 3; ++j) {
              if (atom[atom[i]->bonded[j].num]->atom_type == 7) {
                /*
                aldehyde
                */
                strcpy(atom[i]->atom_name, " CH");
                break;
              }
              else if ((atom[atom[i]->bonded[j].num]->atom_type == 9)
                || (atom[atom[i]->bonded[j].num]->atom_type == 16)) {
                /*
                oxime or thioaldehyde
                */
                strcpy(atom[i]->atom_name, " C1=");
                break;
              }
            }
          }
        }
        else if (atom[i]->atom_type == 41) {
          strcpy(atom[i]->atom_name, " C-H");
        }
        else if (atom[i]->atom_type == 4) {
          strcpy(atom[i]->atom_name, " C1#");
        }
        break;

        case 0:
        if ((atom[i]->atom_type == 1) || (atom[i]->atom_type == 20)
          || (atom[i]->atom_type == 22)) {
          strcpy(atom[i]->atom_name, " C0");
        }
        else if ((atom[i]->atom_type == 37) || (atom[i]->atom_type == 2)
          || (atom[i]->atom_type == 63) || (atom[i]->atom_type == 64)
          || (atom[i]->atom_type == 78) || (atom[i]->atom_type == 30)
          || (atom[i]->atom_type == 80)) {
          strcpy(atom[i]->atom_name, " C=");
        }
        else if (atom[i]->atom_type == 57) {
          /*
          guanidinium carbon
          */
          strcpy(atom[i]->atom_name, " C+1");
        }
        else if ((atom[i]->atom_type == 3) || (atom[i]->atom_type == 41)) {
          if (!n_nitro) {
            if (n_oxy == 2) {
              for (j = 0, oxy_bound_to_hydro = 0, formal_charge = 0.0,
                sdf_charge = 0; j < atom[i]->n_bonded; ++j) {
                if (!strcmp(atom[atom[i]->bonded[j].num]->element, "O")) {
                  oxy = atom[atom[i]->bonded[j].num];
                  if ((oxy->sdf_charge == 0) || ((oxy->formal_charge < 0.05) && (oxy->formal_charge > -0.05))){
                    for (k = 0; (k < oxy->n_bonded) && (!oxy_bound_to_hydro); ++k) {
                      /*
                      if charge is zero and hydrogen is bound to
                      oxygen, it is an unionised carboxylic acid
                      */
                      oxy_bound_to_hydro = (!strcmp(atom[oxy->bonded[k].num]->element, "H"));
                    }
                  }
                  sdf_charge += oxy->sdf_charge;
                  formal_charge += oxy->formal_charge;
                }
              }
              if (((formal_charge < 0.05) && (formal_charge > -0.05)) || (sdf_charge == 0)) {
                /*
                unionised carboxylic acid, ester or anhydride
                */
                strcpy(atom[i]->atom_name, (oxy_bound_to_hydro ? " C" : " C="));
              }
              else if (((formal_charge < -0.95) && (formal_charge > -1.05)) || (sdf_charge == -1)) {
                /*
                ionised carboxylic acid
                */
                strcpy(atom[i]->atom_name, " C-1");
              }
              else if ((formal_charge < -1.95) || (sdf_charge == -2)) {
                /*
                ionised carbonate
                */
                strcpy(atom[i]->atom_name, " C-2");
              }
            }
            else if ((n_oxy == 1) && (!n_sulf)) {
              /*
              unionised carbonyl or acyl halide
              */
              strcpy(atom[i]->atom_name, " C");
            }
            else if (n_sulf >= 1) {
              /*
              thioester or thioketone
              */
              strcpy(atom[i]->atom_name, " C=");
            }
          }
          else {
            if (!n_sulf) {
              /*
              amide or carbamate
              */
              strcpy(atom[i]->atom_name, " C");
            }
            else {
              for (j = 0; j < atom[i]->n_bonded; ++j) {
                if ((atom[atom[i]->bonded[j].num]->atom_type == 9)
                  || (atom[atom[i]->bonded[j].num]->atom_type == 16)) {
                  /*
                  thioamide, thio or dithiocarbamate (>N-C(H,C,O,S)=S)
                  */
                  strcpy(atom[i]->atom_name, " C=");
                  break;
                }
                else if (atom[atom[i]->bonded[j].num]->atom_type == 7) {
                  /*
                  thiocarbamate (>N-C(S)=O)
                  */
                  strcpy(atom[i]->atom_name, " C ");
                  break;
                }
              }
            }
          }
        }
        else if (atom[i]->atom_type == 4) {
          /*
          acetylene
          */
          strcpy(atom[i]->atom_name, " C#");
        }
        else if (atom[i]->atom_type == 60) {
          /*
          acetylene anion
          */
          strcpy(atom[i]->atom_name, " C:#");
        }
        break;
      }
    }
    else if (!strcmp(atom[i]->element, "N")) {
      /*
      NITROGEN
      */
      switch (n_hydro) {
        case 3:
        if (atom[i]->atom_type == 34) {
          /*
          protonated primary amine
          */
          strcpy(atom[i]->atom_name, " N3+");
        }
        break;

        case 2:
        if (atom[i]->atom_type == 34) {
          /*
          protonated secondary amine
          */
          strcpy(atom[i]->atom_name, " N2+");
        }
        else if (atom[i]->atom_type == 8) {
          /*
          primary amine
          */
          strcpy(atom[i]->atom_name, " N2:");
        }
        else if (atom[i]->atom_type == 54) {
          /*
          protonated primary imine
          */
          strcpy(atom[i]->atom_name, " N2=");
        }
        else if ((atom[i]->atom_type == 10) || (atom[i]->atom_type == 40)
          || (atom[i]->atom_type == 43) || (atom[i]->atom_type == 55)
          || (atom[i]->atom_type == 56)) {
          /*
          amide, amidinium, sulfonamide, aniline or guanidinium
          */
          strcpy(atom[i]->atom_name, " N2");
        }
        break;

        case 1:
        if (atom[i]->atom_type == 34) {
          /*
          protonated tertiary amine
          */
          strcpy(atom[i]->atom_name, " N1+");
        }
        else if ((atom[i]->atom_type == 54) || (atom[i]->atom_type == 58)
          || (atom[i]->atom_type == 81)) {
          /*
          protonated secondary imine or imidazolium
          */
          strcpy(atom[i]->atom_name, " N1=");
        }
        else if (atom[i]->atom_type == 8) {
          /*
          secondary amine, we are not making any difference
          between "N1:" and "N1<", but also GREATER
          does not seem to distinguish between the two
          */
          strcpy(atom[i]->atom_name, " N1:");
        }
        else if ((atom[i]->atom_type == 10) || (atom[i]->atom_type == 40)
          || (atom[i]->atom_type == 43) || (atom[i]->atom_type == 56)
          || (atom[i]->atom_type == 39)) {
          /*
          N-alkylamide, N-alkylsulfonamide, N-alkylaniline, guanidium, pyrrole
          */
          strcpy(atom[i]->atom_name, " N1");
        }
        else if ((atom[i]->atom_type == 9) || (atom[i]->atom_type == 62)) {
          /*
          primary imine or sulfonamide anion
          */
          strcpy(atom[i]->atom_name, " NH=");
        }
        else if (atom[i]->atom_type == 61) {
          /*
          isonitrile nitrogen bound to hydrogen
          */
          strcpy(atom[i]->atom_name, " N1#");
        }
        break;

        case 0:
        if (atom[i]->atom_type == 34) {
          /*
          quaternary ammonium
          */
          strcpy(atom[i]->atom_name, " N+");
        }
        else if (atom[i]->atom_type == 8) {
          /*
          tertiary amine
          */
          strcpy(atom[i]->atom_name, " N:");
        }
        else if ((atom[i]->atom_type == 55) && (atom[i]->atom_type == 58)) {
          /*
          nitrogen in alkylpyridinium, alkylimidazolium, etc.
          */
          strcpy(atom[i]->atom_name, " N=");
        }
        else if (((atom[i]->atom_type == 10) || (atom[i]->atom_type == 40)
          || (atom[i]->atom_type == 43) || (atom[i]->atom_type == 56)
          || (atom[i]->atom_type == 39) || (atom[i]->atom_type == 45)
          || (atom[i]->atom_type == 46) || (atom[i]->atom_type == 67)
          || (atom[i]->atom_type == 68) || (atom[i]->atom_type == 69)
          || (atom[i]->atom_type == 82))
          && (n_oxy < 3)) {
          /*
          N,N-dialkylamide, N,N-dialkylsulfonamide, N,N-dialkylaniline,
          nitroso, nitro, N-alkylpyrrole, alkylguanidinium,
          pyridine N-oxide, trialkylamine N-oxide,
          N-alkyl alkylimine N-oxide
          */
          strcpy(atom[i]->atom_name, " N0");
        }
        else if (atom[i]->atom_type == 61) {
          /*
          isonitrile nitrogen
          */
          strcpy(atom[i]->atom_name, " N#+");
        }
        else if (atom[i]->atom_type == 62) {
          /*
          ionized N-acyl/alkylsulfonamide
          */
          strcpy(atom[i]->atom_name, " N0:");
        }
        else if ((atom[i]->atom_type == 9) || (atom[i]->atom_type == 38)
          || (atom[i]->atom_type == 65) || (atom[i]->atom_type == 66)
          || (atom[i]->atom_type == 79)) {
          /*
          secondary imine, =N- in aromatic rings
          */
          strcpy(atom[i]->atom_name, " N:=");
        }
        else if (atom[i]->atom_type == 47) {
          /*
          azide -1 nitrogens
          */
          strcpy(atom[i]->atom_name, " N::");
        }
        else if (atom[i]->atom_type == 53) {
          /*
          azide +1 nitrogens
          */
          strcpy(atom[i]->atom_name, " N#");
        }
        else if (atom[i]->atom_type == 42) {
          /*
          triple-bonded nitrogen
          */
          strcpy(atom[i]->atom_name, " N:#");
        }
        else if ((atom[i]->atom_type == 45) && (n_oxy == 3)) {
          /*
          nitrate anion
          */
          strcpy(atom[i]->atom_name, " N-1");
        }
        else if (atom[i]->atom_type == 76) {
          /*
          anionic tetrazole nitrogen
          */
          for (j = 0, k = 0; j < atom[i]->n_bonded; ++j) {
            if (!strcmp(atom[atom[i]->bonded[j].num]->element, "N")) {
              ++k;
            }
          }
          if (k == 2) {
            strcpy(atom[i]->atom_name, " N:=");
          }
          else if (k == 1) {
            for (j = 0; j < atom[i]->n_bonded; ++j) {
              if (!strcmp(atom[atom[i]->bonded[j].num]->element, "C")) {
                strcpy(atom[atom[i]->bonded[j].num]->atom_name, " C+1");
              }
              else if (!strcmp(atom[atom[i]->bonded[j].num]->element, "N")) {
                strcpy(atom[i]->atom_name, " N-:");
              }
            }
          }
        }
      }
    }
    else if (!strcmp(atom[i]->element, "O")) {
      /*
      OXYGEN
      */
      switch (n_hydro) {
        case 1:
        if (atom[i]->atom_type == 6) {
          if (atom[i]->n_bonded == 2) {
            for (j = 0; j < 2; ++j) {
              /*
              hydroxyl oxygen not belonging to phenol, enol,
              carboxylic acid
              */
              if (strcmp(atom[atom[i]->bonded[j].num]->element, "H")
                && (atom[atom[i]->bonded[j].num]->atom_type != 37)
                && (atom[atom[i]->bonded[j].num]->atom_type != 2)
                && (atom[atom[i]->bonded[j].num]->atom_type != 63)
                && (atom[atom[i]->bonded[j].num]->atom_type != 64)
                && (atom[atom[i]->bonded[j].num]->atom_type != 78)
                && (atom[atom[i]->bonded[j].num]->atom_type != 30)
                && (atom[atom[i]->bonded[j].num]->atom_type != 3)) {
                strcpy(atom[i]->atom_name, " O1");
                break;
              }
            }
          }
          if (strcmp(atom[i]->atom_name, " O1")) {
            /*
            hydroxyl oxygen belonging to phenol, enol,
            carboxylic acid
            */
            strcpy(atom[i]->atom_name, " OH");
          }
        }
        break;

        case 0:
        if (atom[i]->atom_type == 7) {
          if (atom[i]->n_bonded == 1) {
            if (atom[atom[i]->bonded[0].num]->atom_type == 46) {
              /*
              nitroso group oxygen
              */
              strcpy(atom[i]->atom_name, " ON");
            }
            else if (atom[atom[i]->bonded[0].num]->atom_type == 17) {
              /*
              sulfoxide group oxygen
              */
              strcpy(atom[i]->atom_name, " OS");
            }
            else {
              /*
              carbonyl oxygen
              */
              strcpy(atom[i]->atom_name, " O");
            }
          }
        }
        else if (atom[i]->atom_type == 32) {
          if (atom[i]->n_bonded == 1) {
            if ((atom[atom[i]->bonded[0].num]->atom_type == 41)
              || (atom[atom[i]->bonded[0].num]->atom_type == 67)
              || (atom[atom[i]->bonded[0].num]->atom_type == 68)
              || (atom[atom[i]->bonded[0].num]->atom_type == 69)
              || (atom[atom[i]->bonded[0].num]->atom_type == 82)) {
              /*
              oxygen in pyridine-N-oxide, trialkylamine N-oxide,
              N-alkyl alkylimine N-oxide
              in carboxylate anion
              */
              strcpy(atom[i]->atom_name, " O::");
            }
            else if (atom[atom[i]->bonded[0].num]->atom_type == 25) {
              /*
              oxygen in phosphate anion
              */
              strcpy(atom[i]->atom_name, " O=");
            }
            else if (atom[atom[i]->bonded[0].num]->atom_type == 45) {
              /*
              nitro group oxygen
              */
              strcpy(atom[i]->atom_name, " ON");
            }
            else if ((atom[atom[i]->bonded[0].num]->atom_type == 18)
              || (atom[atom[i]->bonded[0].num]->atom_type == 73)) {
              /*
              sulfonyl group oxygen or ionized sulfinate
              */
              strcpy(atom[i]->atom_name, " O=S");
            }
            else if (atom[atom[i]->bonded[0].num]->atom_type == 25) {
              /*
              phosphate oxygen
              */
              strcpy(atom[i]->atom_name, " O=");
            }
          }
        }
        else if (atom[i]->atom_type == 35) {
          if (atom[i]->n_bonded == 1) {
            if (atom[atom[i]->bonded[0].num]->atom_type == 4) {
              /*
              oxygen in cyanate
              */
              strcpy(atom[i]->atom_name, " O#");
            }
            else if (atom[atom[i]->bonded[0].num]->atom_type == 1) {
              /*
              oxygen in alkoxide
              */
              strcpy(atom[i]->atom_name, " O-1");
            }
            else if ((atom[atom[i]->bonded[j].num]->atom_type == 37)
              || (atom[atom[i]->bonded[j].num]->atom_type == 2)
              || (atom[atom[i]->bonded[j].num]->atom_type == 63)
              || (atom[atom[i]->bonded[j].num]->atom_type == 64)
              || (atom[atom[i]->bonded[j].num]->atom_type == 78)
              || (atom[atom[i]->bonded[j].num]->atom_type == 30)
              || (atom[atom[i]->bonded[j].num]->atom_type == 80)) {
              /*
              oxygen in phenoxide/enolate
              */
              strcpy(atom[i]->atom_name, " O-");
            }
          }
        }
        else if (atom[i]->atom_type == 6) {
          /*
          count how many alkyl/aryl-vinyl/carbonyl substituents are there
          */
          if (atom[i]->n_bonded == 2) {
            for (j = 0, n_aryl_vinyl = 0, n_alkyl = 0, n_carb_sulf_phos = 0, n_nitro = 0; j < 2; ++j) {
              if (atom[atom[i]->bonded[j].num]->atom_type == 1) {
                ++n_alkyl;
              }
              if (!strcmp(atom[atom[i]->bonded[j].num]->element, "N")) {
                ++n_nitro;
              }
              else if ((atom[atom[i]->bonded[j].num]->atom_type == 37)
                || (atom[atom[i]->bonded[j].num]->atom_type == 2)
                || (atom[atom[i]->bonded[j].num]->atom_type == 63)
                || (atom[atom[i]->bonded[j].num]->atom_type == 64)
                || (atom[atom[i]->bonded[j].num]->atom_type == 78)
                || (atom[atom[i]->bonded[j].num]->atom_type == 30)
                || (atom[atom[i]->bonded[j].num]->atom_type == 80)) {
                ++n_aryl_vinyl;
              }
              else if ((atom[atom[i]->bonded[j].num]->atom_type == 3)
                || (atom[atom[i]->bonded[j].num]->atom_type == 18)
                || (atom[atom[i]->bonded[j].num]->atom_type == 17)
                || (atom[atom[i]->bonded[j].num]->atom_type == 25)) {
                ++n_carb_sulf_phos;
              }
            }
            if ((n_alkyl == 2) && (n_aryl_vinyl == 0)) {
              /*
              dialkylether
              */
              strcpy(atom[i]->atom_name, " OC2");
            }
            else if (((n_alkyl == 1) || (n_nitro == 1)) && (n_aryl_vinyl == 1)) {
              /*
              alkylaryl(vinyl)ether
              */
              strcpy(atom[i]->atom_name, " OC1");
            }
            else if ((n_alkyl == 0) && (n_aryl_vinyl == 2)) {
              /*
              diaryl(vinyl)ether
              */
              strcpy(atom[i]->atom_name, " OC=");
            }
            else if (n_carb_sulf_phos) {
              /*
              carboxy/sulfinate/sulfonate/phosphate ester or anhydride
              */
              strcpy(atom[i]->atom_name, " OES");
            }
          }

        }
        else if (atom[i]->atom_type == 59) {
          strcpy(atom[i]->atom_name, " OFU");
        }
        break;
      }
    }
    else if (!strcmp(atom[i]->element, "S")) {
      /*
      SULFUR
      */
      switch (n_hydro) {
        case 1:
        if (atom[i]->atom_type == 15) {
          /*
          thiol
          */
          strcpy(atom[i]->atom_name, " S1");
        }
        break;

        case 0:
        if (atom[i]->atom_type == 44) {
          /*
          sulfide/thiophene/thiazole
          */
          strcpy(atom[i]->atom_name, " STH");
        }
        if (atom[i]->atom_type == 15) {
          if (atom[i]->n_bonded == 2) {
            for (j = 0; j < 2; ++j) {
              if (atom[atom[i]->bonded[j].num]->atom_type == 15) {
                /*
                disulfide
                */
                strcpy(atom[i]->atom_name, " S");
                break;
              }
            }
            if (strcmp(atom[i]->atom_name, " S")) {
              /*
              sulfide
              */
              strcpy(atom[i]->atom_name, " STH");
            }
          }
        }
        else if (atom[i]->atom_type == 16) {
          /*
          thioketone sulfur
          */
          strcpy(atom[i]->atom_name, " S");
        }
        else if (atom[i]->atom_type == 18) {
          for (j = 0, formal_charge = 0.0, sdf_charge = 0; j < atom[i]->n_bonded; ++j) {
            sdf_charge += atom[atom[i]->bonded[j].num]->sdf_charge;
            formal_charge += atom[atom[i]->bonded[j].num]->formal_charge;
          }
          if (((formal_charge < 0.05) && (formal_charge > -0.05)) || (sdf_charge == 0)) {
            /*
            unionized sulfone/sulfonamide
            */
            strcpy(atom[i]->atom_name, " S");
          }
          else if (((formal_charge < -0.95) && (formal_charge > -1.05)) || (sdf_charge == -1)) {
            /*
            ionized sulfonate, sulfonamide
            or sulfamate
            */
            strcpy(atom[i]->atom_name, " S-1");
          }
          else {
            /*
            sulfate
            */
            strcpy(atom[i]->atom_name, " S-2");
          }
        }
        else if (atom[i]->atom_type == 17) {
          if (atom[i]->n_bonded == 3) {
            for (j = 0; j < 3; ++j) {
              if (atom[atom[i]->bonded[j].num]->atom_type == 7) {
                /*
                sulfoxide sulfur
                */
                strcpy(atom[i]->atom_name, " SO");
                break;
              }
            }
            if (strcmp(atom[i]->atom_name, " SO")) {
              /*
              sulfonium salt
              */
              strcpy(atom[i]->atom_name, " S+1");
            }
          }
        }
        else if (atom[i]->atom_type == 72) {
          /*
          ionized thiol
          */
          strcpy(atom[i]->atom_name, " S-");
        }
        else if (atom[i]->atom_type == 73) {
          /*
          ionized sulfinate
          */
          strcpy(atom[i]->atom_name, " S-1");
        }
      }
    }
    else if ((!strcmp(atom[i]->element, "P")) && (!n_hydro)) {
      /*
      PHOSPHORUS
      */
      if (atom[i]->atom_type == 25) {
        for (j = 0, n_free = 0; j < atom[i]->n_bonded; ++j) {
          /*
          check if bound oxygens/sulfurs are bonded
          to another atom or not
          */
          if (((!strcmp(atom[atom[i]->bonded[j].num]->element, "O"))
            || (!strcmp(atom[atom[i]->bonded[j].num]->element, "S")))
            && (atom[atom[i]->bonded[j].num]->n_bonded == 1)) {
            ++n_free;
          }
        }
        if (n_free == 1) {
          strcpy(atom[i]->atom_name, " P");
        }
        else if ((n_free > 1) && (n_free < 5)) {
          sprintf(atom[i]->atom_name, " P-%d", n_free - 1);
        }
      }
      else if (atom[i]->atom_type == 26) {
        if (atom[i]->n_bonded > 3) {
          strcpy(atom[i]->atom_name, " P+1");
        }
        else {
          strcpy(atom[i]->atom_name, " P");
        }
      }
    }
    else if ((!strcmp(atom[i]->element, "Si")) && (!n_hydro)) {
      /*
      SILICON
      */
      for (j = 0, n_free = 0; j < atom[i]->n_bonded; ++j) {
        /*
        check if bound oxygens/sulfurs are bonded
        to another atom or not
        */
        if (((!strcmp(atom[atom[i]->bonded[j].num]->element, "O"))
          || (!strcmp(atom[atom[i]->bonded[j].num]->element, "S")))
          && (atom[atom[i]->bonded[j].num]->n_bonded == 1)) {
          ++n_free;
        }
      }
      if (n_free == 0) {
        strcpy(atom[i]->atom_name, "SI");
      }
      else if ((n_free > 0) && (n_free < 5)) {
        sprintf(atom[i]->atom_name, "SI-%d", n_free);
      }
    }
    else if ((!strcmp(atom[i]->element, "As")) && (!n_hydro)) {
      /*
      ARSENIC
      */
      strcpy(atom[i]->atom_name, "AS");
    }
    else if ((!strcmp(atom[i]->element, "Sb")) && (!n_hydro)) {
      /*
      ANTIMONY
      */
      strcpy(atom[i]->atom_name, "SB");
    }
    else if (!strcmp(atom[i]->element, "Se")) {
      /*
      SELENIUM
      */
      switch (n_hydro) {
        case 1:
        strcpy(atom[i]->atom_name, "SE1");
        break;

        case 0:
        for (j = 0, formal_charge = 0.0, sdf_charge = 0; j < atom[i]->n_bonded; ++j) {
          sdf_charge += atom[i]->sdf_charge;
          formal_charge += atom[i]->formal_charge;
        }
        if (((formal_charge > -0.05) && (formal_charge < 0.05)) || (sdf_charge == 0)) {
          strcpy(atom[i]->atom_name, "SE");
        }
        else if (((formal_charge > -1.05) && (formal_charge < -0.95)) || (sdf_charge == -1)) {
          strcpy(atom[i]->atom_name, "SE-1");
        }
        else if (((formal_charge > -2.05) && (formal_charge < -1.95)) || (sdf_charge == -2)) {
          strcpy(atom[i]->atom_name, "SE-2");
        }
        break;
      }
    }
    else if (!strcmp(atom[i]->element, "F")) {
      /*
      FLUORINE
      */
      if (atom[i]->n_bonded == 1) {
        if (!strcmp(atom[atom[i]->bonded[0].num]->element, "C")) {
          strcpy(atom[i]->atom_name, " F");
          carbon = atom[atom[i]->bonded[0].num];
          for (j = 0, n_fluorine = 0; j < carbon->n_bonded; ++j) {
            if (!strcmp(atom[carbon->bonded[j].num]->element, "F")) {
              ++n_fluorine;
            }
          }
          if (n_fluorine >= 2) {
            /*
            CF2 or CF3
            */
            strcpy(atom[i]->atom_name, " F3");
          }
          else {
            /*
            look at the nearest carbon(s)
            if they are aromatic/olefinic carbons,
            check if they bear fluorine substituents
            and, if so, if they are in cis position
            */
            for (j = 0; j < carbon->n_bonded; ++j) {
              if ((atom[carbon->bonded[j].num]->atom_type == 37)
                || (atom[carbon->bonded[j].num]->atom_type == 2)
                || (atom[carbon->bonded[j].num]->atom_type == 63)
                || (atom[carbon->bonded[j].num]->atom_type == 64)
                || (atom[carbon->bonded[j].num]->atom_type == 78)
                || (atom[carbon->bonded[j].num]->atom_type == 30)
                || (atom[carbon->bonded[j].num]->atom_type == 80)) {
                neighbor = atom[carbon->bonded[j].num];
                for (k = 0; k < neighbor->n_bonded; ++k) {
                  if (!strcmp(atom[neighbor->bonded[k].num]->element, "F")) {
                    dist = squared_euclidean_distance
                      (atom[neighbor->bonded[k].num]->coord, atom[i]->coord);
                    if ((dist > 6.76) && (dist < 9.61)) {
                      /*
                      cis configuration
                      */
                      strcpy(atom[i]->atom_name, " F3");
                      break;
                    }
                  }
                }
              }
            }
          }
        }
        else if ((atom[atom[i]->bonded[0].num]->atom_type == 17)
          || (atom[atom[i]->bonded[0].num]->atom_type == 18)) {
          strcpy(atom[i]->atom_name, " FS");
        }
      }
      else if ((atom[i]->n_bonded == 0) && ((atom[i]->sdf_charge == -1)
        || ((atom[i]->formal_charge < -0.05) && (atom[i]->formal_charge > -1.05)))) {
        strcpy(atom[i]->atom_name, " F-");
      }
    }
    else if (!strcmp(atom[i]->element, "Cl")) {
      /*
      CHLORINE
      */
      if (atom[i]->n_bonded == 1) {
        strcpy(atom[i]->atom_name, "CL");
      }
      else if ((atom[i]->n_bonded == 0) && ((atom[i]->sdf_charge == -1)
        || ((atom[i]->formal_charge < -0.05) && (atom[i]->formal_charge > -1.05)))) {
        strcpy(atom[i]->atom_name, "CL-");
      }
    }
    else if (!strcmp(atom[i]->element, "Br")) {
      /*
      BROMINE
      */
      if (atom[i]->n_bonded == 1) {
        strcpy(atom[i]->atom_name, "BR");
      }
      else if ((atom[i]->n_bonded == 0) && ((atom[i]->sdf_charge == -1)
        || ((atom[i]->formal_charge < -0.05) && (atom[i]->formal_charge > -1.05)))) {
        strcpy(atom[i]->atom_name, "BR-");
      }
    }
    else if (!strcmp(atom[i]->element, "I")) {
      /*
      IODINE
      */
      if (atom[i]->n_bonded == 1) {
        strcpy(atom[i]->atom_name, " I");
      }
      else if ((atom[i]->n_bonded == 0) && ((atom[i]->sdf_charge == -1)
        || ((atom[i]->formal_charge < -0.05) && (atom[i]->formal_charge > -1.05)))) {
        strcpy(atom[i]->atom_name, " I-");
      }
    }
    else if ((atom[i]->formal_charge > 0.95) || (atom[i]->sdf_charge > 0)) {
      /*
      UNI, DI, TRI, TETRA -VALENT CATIONS
      */
      for (k = 0; k < 4; ++k) {
        if (((atom[i]->formal_charge > (0.95 + (double)k))
          && (atom[i]->formal_charge < (1.05 + (double)k)))
          || (atom[i]->sdf_charge == (k + 1))) {
          j = 0;
          while (valence[k].cations[j]
            && strcmp(atom[i]->element, valence[k].cations[j])) {
            ++j;
          }
          if (valence[k].cations[j]) {
            sprintf(atom[i]->atom_name, "%2s", valence[k].cations[j]);
            atom[i]->atom_name[1] = (char)toupper((int)(atom[i]->atom_name[1]));
            strcat(atom[i]->atom_name, "+");
            if (k) {
              sprintf(&(atom[i]->atom_name[strlen(atom[i]->atom_name)]), "%d", k + 1);
            }
          }
          break;
        }
      }
    }
    else if (!strcmp(atom[i]->element, "H")) {
      /*
      HYDROGEN
      */
      strcpy(atom[i]->atom_name, " H");
    }
    if (!strcmp(atom[i]->atom_name, "UNK")) {
      error = FL_UNKNOWN_ATOM_TYPE;
    }
  }
  
  return error;
}


int prep_cosmo_input(O3Data *od, TaskInfo *task, AtomInfo **atom, int object_num)
{
  char buffer[BUF_LEN];
  char input_dir[BUF_LEN];
  char basis_set_name[MAX_NAME_LEN];
  char lowercase_elem[MAX_NAME_LEN];
  unsigned char used_elem_array[BUF_LEN];
  unsigned char elem_i = 0;
  int i;
  int n;
  int n_atoms;
  int found = 0;
  int result;
  int pos;
  int error = 0;
  int pid = 0;
  double total_formal_charge = 0.0;
  ProgExeInfo prog_exe_info;
  FileDescriptor inp_fd;
  FileDescriptor log_fd;
  FileDescriptor coord_fd;
  #ifndef WIN32
  FILE *pipe_handle = NULL;
  #else
  DWORD n_chr = 0;
  HANDLE pipe_handle;
  #endif


  memset(&inp_fd, 0, sizeof(FileDescriptor));
  memset(&log_fd, 0, sizeof(FileDescriptor));
  memset(&coord_fd, 0, sizeof(FileDescriptor));
  memset(buffer, 0, BUF_LEN);
  memset(input_dir, 0, BUF_LEN);
  memset(basis_set_name, 0, MAX_NAME_LEN);
  memset(lowercase_elem, 0, MAX_NAME_LEN);
  result = fill_atom_info(od, task, atom, NULL, object_num, O3_MMFF94);
  if (result) {
    return result;
  }
  sprintf(input_dir, "%s%c%s_%04d", od->field.qm_dir, SEPARATOR,
    cosmo_label, od->al.mol_info[object_num]->object_id);
  if (dexist(input_dir)) {
    remove_recursive(input_dir);
  }
  #ifndef WIN32
  result = mkdir(input_dir, S_IRWXU | S_IRGRP | S_IROTH);
  #else
  result = mkdir(input_dir);
  #endif
  if (result == -1) {
    return FL_CANNOT_CREATE_SCRDIR;
  }
  sprintf(inp_fd.name, "%s%cdefine_%04d.inp", input_dir,
    SEPARATOR, od->al.mol_info[object_num]->object_id);
  /*
  find total formal charge for this molecule
  */
  n_atoms = 0;
  total_formal_charge = 0.0;
  while (atom[n_atoms]->atom_type != -1) {
    total_formal_charge += atom[n_atoms]->formal_charge;
    ++n_atoms;
  }
  if (!(inp_fd.handle = fopen(inp_fd.name, "wb+"))) {
    O3_ERROR_LOCATE(task);
    O3_ERROR_STRING(task, inp_fd.name);
    return FL_CANNOT_WRITE_INP_FILE;
  }
  /*
  prepare job file in TURBOMOLE format
  */
  switch (od->field.basis_set) {
    case O3_SVP:
    strcpy(basis_set_name, "def-SVP");
    break;

    case O3_TZVP:
    strcpy(basis_set_name, "def-TZVP");
    break;
  }
  fprintf(inp_fd.handle,
    "\n"
    "%04d_cosmo\n"
    "a coord\n"
    "*\n"
    "no\n"
    "b all %s\n",
    od->al.mol_info[object_num]->object_id,
    basis_set_name);
  sprintf(coord_fd.name, "%s%ccoord", input_dir, SEPARATOR);
  if (!(coord_fd.handle = fopen(coord_fd.name, "wb+"))) {
    O3_ERROR_LOCATE(task);
    O3_ERROR_STRING(task, coord_fd.name);
    return FL_CANNOT_WRITE_INP_FILE;
  }
  fprintf(coord_fd.handle, "$coord\n");
  memset(used_elem_array, 0xFF, BUF_LEN);
  for (n = 0; n < n_atoms; ++n) {
    elem_i = 0;
    while (element_data[(int)elem_i].atomic_number) {
      if (!strcasecmp(atom[n]->element, element_data[(int)elem_i].atomic_symbol)) {
        break;
      }
      ++elem_i;
    }
    if (!(element_data[(int)elem_i].atomic_number)) {
      O3_ERROR_LOCATE(task);
      return FL_UNKNOWN_ATOM_TYPE;
    }
    for (i = 0, found = 0; ((used_elem_array[i] != 0xFF) && (!found)); ++i) {
      found = (used_elem_array[i] == elem_i);
    }
    if (!found) {
      used_elem_array[i] = elem_i;
    }
    strcpy(lowercase_elem, atom[n]->element);
    string_to_lowercase(lowercase_elem);
    fprintf(coord_fd.handle,
      "%20.14lf  %20.14lf  %20.14lf      %s\n",
      atom[n]->coord[0] / BOHR_RADIUS,
      atom[n]->coord[1] / BOHR_RADIUS,
      atom[n]->coord[2] / BOHR_RADIUS,
      lowercase_elem);
  }
  fprintf(coord_fd.handle,
    "$user-defined bonds\n"
    "$end\n");
  fclose(coord_fd.handle);
  fprintf(inp_fd.handle,
    "*\n"
    "eht\n"
    "y\n"
    "%d\n"
    "%s"
    "dft\n"
    "on\n"
    "func b-p\n"
    "*\n"
    "ri\n"
    "on\n"
    "m %d\n"
    "*\n"
    "*\n",
    (int)safe_rint(total_formal_charge),
    ((od->field.spin == O3_UNRESTRICTED)
    ? "n\n"
    "uf 0\n"
    "*\n"
    "n\n"
    : "y\n"),
    MAX_TURBOMOLE_RI_MEMORY);
  rewind(inp_fd.handle);
  memset(&prog_exe_info, 0, sizeof(ProgExeInfo));
  prog_exe_info.proc_env = fill_env
    (od, turbomole_env, od->field.qm_exe_path, object_num);
  if (!(prog_exe_info.proc_env)) {
    O3_ERROR_LOCATE(task);
    return FL_OUT_OF_MEMORY;
  }
  prog_exe_info.need_stdin = NEED_STDIN_NORMAL;
  prog_exe_info.stdout_fd = &log_fd;
  prog_exe_info.stderr_fd = &log_fd;
  prog_exe_info.exedir = input_dir;
  prog_exe_info.sep_proc_grp = 1;
  sprintf(log_fd.name, "%s%cdefine_%04d.log", input_dir,
    SEPARATOR, od->al.mol_info[object_num]->object_id);
  sprintf(prog_exe_info.command_line, "%s%c%s",
    od->field.qm_exe_path, SEPARATOR, TURBOMOLE_DEFINE_EXE);
  pid = ext_program_exe(&prog_exe_info, &error);
  if (error) {
    O3_ERROR_LOCATE(task);
    return error;
  }
  #ifndef WIN32
  if (!(pipe_handle = fdopen(prog_exe_info.pipe_des[1], "w"))) {
    O3_ERROR_LOCATE(task);
    return FL_CANNOT_CREATE_CHANNELS;
  }
  #else
  pipe_handle = prog_exe_info.stdin_wr;
  #endif
  /*
  pipe job information into define
  */
  while (fgets(buffer, BUF_LEN, inp_fd.handle)) {
    buffer[BUF_LEN - 1] = '\0';
    FWRITE_WRAP(pipe_handle, buffer, &n_chr);
  }
  fclose(inp_fd.handle);
  FFLUSH_WRAP(pipe_handle);
  FCLOSE_WRAP(pipe_handle);
  ext_program_wait(&prog_exe_info, pid);
  if (!(log_fd.handle = fopen(log_fd.name, "rb"))) {
    O3_ERROR_LOCATE(task);
    O3_ERROR_STRING(task, log_fd.name);
    return FL_CANNOT_READ_TEMP_FILE;
  }
  if (!fgrep(log_fd.handle, buffer, TURBOMOLE_NORMAL_TERMINATION)) {
    O3_ERROR_LOCATE(task);
    error = FL_ABNORMAL_TERMINATION;
  }
  fclose(log_fd.handle);
  if (error) {
    return error;
  }
  sprintf(inp_fd.name, "%s%ccontrol", input_dir, SEPARATOR);
  if (!(inp_fd.handle = fopen(inp_fd.name, "rb+"))) {
    O3_ERROR_LOCATE(task);
    O3_ERROR_STRING(task, inp_fd.name);
    return FL_CANNOT_READ_TEMP_FILE;
  }
  found = 0;
  while (!found) {
    pos = ftell(inp_fd.handle);
    if (!fgets(buffer, BUF_LEN, inp_fd.handle)) {
      break;
    }
    buffer[BUF_LEN - 1] = '\0';
    remove_newline(buffer);
    found = (!strncmp(buffer, "$scfiterlimit", 13));
  }
  if (!found) {
    fclose(inp_fd.handle);
    O3_ERROR_LOCATE(task);
    O3_ERROR_STRING(task, inp_fd.name);
    return FL_CANNOT_READ_TEMP_FILE;
  }
  if (fseek(inp_fd.handle, pos, SEEK_SET)) {
    fclose(inp_fd.handle);
    O3_ERROR_LOCATE(task);
    O3_ERROR_STRING(task, inp_fd.name);
    return FL_CANNOT_READ_TEMP_FILE;
  }
  n = strlen(buffer) - 13;
  sprintf(buffer, "$scfiterlimit%%%dd\n", n);
  fprintf(inp_fd.handle, buffer, MAX_TURBOMOLE_SCFITER);
  found = 0;
  while (!found) {
    pos = ftell(inp_fd.handle);
    if (!fgets(buffer, BUF_LEN, inp_fd.handle)) {
      break;
    }
    buffer[BUF_LEN - 1] = '\0';
    found = (!strncmp(buffer, "$last step", 10));
  }
  if (!found) {
    fclose(inp_fd.handle);
    O3_ERROR_LOCATE(task);
    O3_ERROR_STRING(task, inp_fd.name);
    return FL_CANNOT_READ_TEMP_FILE;
  }
  if (fseek(inp_fd.handle, pos, SEEK_SET)) {
    fclose(inp_fd.handle);
    O3_ERROR_LOCATE(task);
    O3_ERROR_STRING(task, inp_fd.name);
    return FL_CANNOT_READ_TEMP_FILE;
  }
  fprintf(inp_fd.handle,
    "$marij\n"
    "$last step     define\n"
    "$end\n");
  rewind(inp_fd.handle);
  if (!fgrep(inp_fd.handle, buffer, "$atoms")) {
    fclose(inp_fd.handle);
    O3_ERROR_LOCATE(task);
    O3_ERROR_STRING(task, inp_fd.name);
    return FL_CANNOT_READ_TEMP_FILE;
  }
  pos = ftell(inp_fd.handle);
  /*
  prepare input for cosmo
  */
  sprintf(coord_fd.name, "%s%ccosmoprep", input_dir, SEPARATOR);
  if (!(coord_fd.handle = fopen(coord_fd.name, "wb"))) {
    O3_ERROR_LOCATE(task);
    O3_ERROR_STRING(task, inp_fd.name);
    return FL_CANNOT_WRITE_INP_FILE;
  }
  fprintf(coord_fd.handle,
    "$cosmo\n"
    "  rsolv = %.4lf\n"
    "$cosmo_atoms\n"
    "# radii in Angstrom units\n", COSMO_DEFAULT_RSOLV);
  i = 0;
  while ((i < BUF_LEN) && (used_elem_array[i] != 0xFF)) {
    if (fseek(inp_fd.handle, pos, SEEK_SET)) {
      fclose(inp_fd.handle);
      O3_ERROR_LOCATE(task);
      O3_ERROR_STRING(task, inp_fd.name);
      return FL_CANNOT_READ_TEMP_FILE;
    }
    strcpy(lowercase_elem, element_data[(int)(used_elem_array[i])].atomic_symbol);
    string_to_lowercase(lowercase_elem);
    strcat(lowercase_elem, " ");
    found = 0;
    while ((!found) && fgets(buffer, BUF_LEN - 1, inp_fd.handle)) {
      buffer[BUF_LEN - 1] = '\0';
      remove_newline(buffer);
      found = (!strncmp(buffer, lowercase_elem, strlen(lowercase_elem)));
    }
    if (!found) {
      fclose(inp_fd.handle);
      O3_ERROR_LOCATE(task);
      O3_ERROR_STRING(task, inp_fd.name);
      return FL_CANNOT_READ_TEMP_FILE;
    }
    fprintf(coord_fd.handle,
      "%s\n"
      "   radius=%8.4lf\n", buffer,
      element_data[(int)(used_elem_array[i])].cosmo_radius);
    ++i;
  }
  fclose(inp_fd.handle);
  fprintf(coord_fd.handle, "$cosmo_out file=%s_%04d"TURBOMOLE_COSMO_EXT"\n",
    cosmo_label, od->al.mol_info[object_num]->object_id);
  fclose(coord_fd.handle);
  if (!fcopy(inp_fd.name, coord_fd.name, "ab")) {
    O3_ERROR_LOCATE(task);
    O3_ERROR_STRING(task, inp_fd.name);
    return FL_CANNOT_READ_TEMP_FILE;
  }
  remove(inp_fd.name);
  if (rename(coord_fd.name, inp_fd.name) == -1) {
    O3_ERROR_LOCATE(task);
    O3_ERROR_STRING(task, inp_fd.name);
    return FL_CANNOT_READ_TEMP_FILE;
  }

  return 0;
}


int prep_cs3d_input(O3Data *od)
{
  char buffer[BUF_LEN];
  char fdir[BUF_LEN];
  char *cosmo_base;
  int object_num = 0;
  FileDescriptor inp_fd;
  
  
  memset(buffer, 0, BUF_LEN);
  memset(fdir, 0, BUF_LEN);
  memset(&inp_fd, 0, sizeof(FileDescriptor));
  sprintf(inp_fd.name, "%s%ccs3d.inp", od->field.qm_dir, SEPARATOR);
  if (!(inp_fd.handle = fopen(inp_fd.name, "wb"))) {
    O3_ERROR_LOCATE(&(od->task));
    O3_ERROR_STRING(&(od->task), inp_fd.name);
    return CANNOT_WRITE_QM_INP_FILE;
  }
  fprintf(inp_fd.handle, "1 NSCSTART=1 MIN1D=0 IDELSIG=%d NRAND=-1 TDIR=\"%s\"\n"
    "SETBOX=%.4f,%.4f,%.4f,%d,%d,%d GRIDS=%.4f NOLSANAM WRTLSPA COMPRESS=%s\n",
    od->field.idelsig,
    od->field.qm_dir,
    od->grid.start_coord[0],
    od->grid.start_coord[1],
    od->grid.start_coord[2],
    od->grid.nodes[0],
    od->grid.nodes[1],
    od->grid.nodes[2],
    od->grid.step[0],
    ((od->field.compress == O3_COMPRESS_GZIP)
    ? "GZ" : ((od->field.compress == O3_COMPRESS_ZIP)
    ? "ZIP" : "NONE")));
  rewind(od->file[TEMP_SORTED_MATCH]->handle);
  while ((object_num < od->object_num)
    && fgets(buffer, BUF_LEN, od->file[TEMP_SORTED_MATCH]->handle)) {
    buffer[BUF_LEN - 1] = '\0';
    remove_newline(buffer);
    if (!object_num) {
      strcpy(fdir, buffer);
      get_dirname(fdir);
      absolute_path(fdir);
      fprintf(inp_fd.handle, "FDIR=\"%s\"\n", fdir);
    }
    cosmo_base = get_basename(buffer);
    fprintf(inp_fd.handle, "%s\n", cosmo_base);
    ++object_num;
  }
  fclose(inp_fd.handle);
  if (object_num != od->object_num) {
    O3_ERROR_LOCATE(&(od->task));
    O3_ERROR_STRING(&(od->task), od->file[TEMP_SORTED_MATCH]->name);
    return CANNOT_READ_TEMP_FILE;
  }
  
  return 0;
}


int prep_qm_input(O3Data *od, TaskInfo *task, AtomInfo **atom, int object_num)
{
  char buffer[BUF_LEN];
  char input_dir[BUF_LEN];
  char basis_set_name[MAX_NAME_LEN];
  char lowercase_elem[MAX_NAME_LEN];
  unsigned char used_elem_array[BUF_LEN];
  unsigned char elem_i = 0;
  int i;
  int n;
  int basis_set_i = 0;
  int n_atoms;
  int found = 0;
  int result;
  int pos;
  int error = 0;
  int pid = 0;
  double total_formal_charge = 0.0;
  ProgExeInfo prog_exe_info;
  FileDescriptor inp_fd;
  FileDescriptor log_fd;
  FileDescriptor basis_fd;
  FileDescriptor coord_fd;
  #ifndef WIN32
  FILE *pipe_handle = NULL;
  #else
  DWORD n_chr = 0;
  HANDLE pipe_handle;
  #endif


  memset(&inp_fd, 0, sizeof(FileDescriptor));
  memset(&log_fd, 0, sizeof(FileDescriptor));
  memset(&basis_fd, 0, sizeof(FileDescriptor));
  memset(&coord_fd, 0, sizeof(FileDescriptor));
  memset(buffer, 0, BUF_LEN);
  memset(input_dir, 0, BUF_LEN);
  memset(basis_set_name, 0, MAX_NAME_LEN);
  memset(lowercase_elem, 0, MAX_NAME_LEN);
  result = fill_atom_info(od, task, atom, NULL, object_num, O3_MMFF94);
  if (result) {
    return result;
  }
  if (od->field.type & PREP_TURBOMOLE_INPUT) {
    sprintf(input_dir, "%s%cturbomole_%04d", od->field.qm_dir,
      SEPARATOR, od->al.mol_info[object_num]->object_id);
    if (dexist(input_dir)) {
      remove_recursive(input_dir);
    }
    #ifndef WIN32
    result = mkdir(input_dir, S_IRWXU | S_IRGRP | S_IROTH);
    #else
    result = mkdir(input_dir);
    #endif
    if (result == -1) {
      return FL_CANNOT_CREATE_SCRDIR;
    }
    sprintf(inp_fd.name, "%s%cdefine_%04d.inp", input_dir,
      SEPARATOR, od->al.mol_info[object_num]->object_id);
  }
  else {
    strcpy(input_dir, od->field.qm_dir);
    sprintf(inp_fd.name, "%s%c%s_%04d.inp", input_dir, SEPARATOR,
      od->field.qm_software, od->al.mol_info[object_num]->object_id);
  }
  /*
  find total formal charge for this molecule
  */
  n_atoms = 0;
  total_formal_charge = 0.0;
  while (atom[n_atoms]->atom_type != -1) {
    total_formal_charge += atom[n_atoms]->formal_charge;
    ++n_atoms;
  }
  if (!(inp_fd.handle = fopen(inp_fd.name, "wb+"))) {
    O3_ERROR_LOCATE(task);
    O3_ERROR_STRING(task, inp_fd.name);
    return FL_CANNOT_WRITE_INP_FILE;
  }
  /*
  an appropriate QM input file has to be generated
  */
  if (od->field.f_func) {
    if (!(od->field.d_func)) {
      od->field.d_func = 1;
    }
    if (!(od->field.p_func)) {
      od->field.p_func = 1;
    }
  }
  if (od->field.p_func) {
    if (od->field.d_func < od->field.p_func) {
      od->field.d_func = od->field.p_func;
    }
  }
  if (od->field.basis_set & EMSL_BASIS_SET) {
    switch (od->field.basis_set & (~EMSL_BASIS_SET)) {
      case O3_3_21G:
      basis_set_i = 0;
      strcpy(basis_set_name, "EMSL_3-21G");
      break;

      case O3_6_311G:
      basis_set_i = 1;
      strcpy(basis_set_name, "EMSL_6-311G");
      break;

      case O3_6_311GXX:
      basis_set_i = 2;
      strcpy(basis_set_name, "EMSL_6-311Gxx");
      break;
    }
  }
  if (od->field.type & PREP_GAUSSIAN_INPUT) {
    /*
    prepare job file in GAUSSIAN format
    */
    fprintf(inp_fd.handle,
      "%%Chk=%s%c%s_%04d.chk\n"
      "# %c%s", input_dir, SEPARATOR, od->field.qm_software,
      od->al.mol_info[object_num]->object_id,
      od->field.spin, od->field.theory);
    if (!(od->field.basis_set & EMSL_BASIS_SET)) {
      switch (od->field.basis_set) {
        case O3_STO_3G:
        fprintf(inp_fd.handle, "/STO-3G");
        break;

        case O3_3_21G:
        fprintf(inp_fd.handle, "/3-21");
        break;

        case O3_6_31G:
        fprintf(inp_fd.handle, "/6-31");
        break;

        case O3_6_311G:
        fprintf(inp_fd.handle, "/6-311");
        break;
      }
      if (od->field.basis_set != O3_STO_3G) {
        if (od->field.diff_sp) {
          fprintf(inp_fd.handle, "+");
        }
        fprintf(inp_fd.handle, "G");
        if (od->field.d_func) {
          if (od->field.d_func > 1) {
            fprintf(inp_fd.handle, "(%dd", od->field.d_func);
          }
          else {
            fprintf(inp_fd.handle, "(d");
          }
          if (od->field.f_func) {
            fprintf(inp_fd.handle, "f");
          }
          if (od->field.p_func) {
            if (od->field.p_func > 1) {
              fprintf(inp_fd.handle, ",%dp", od->field.p_func);
            }
            else {
              fprintf(inp_fd.handle, ",p");
            }
          }
          fprintf(inp_fd.handle, ")");
        }
      }

    }
    else {
      fprintf(inp_fd.handle, "/Gen");
    }
    fprintf(inp_fd.handle,
      " Density=SCF NoSymm Test\n\n"
      "%04d_energy\n\n"
      "%d  1\n",
      od->al.mol_info[object_num]->object_id,
      (int)safe_rint(total_formal_charge));
  }
  else if ((od->field.type & PREP_FIREFLY_INPUT)
    || (od->field.type & PREP_GAMESS_INPUT)) {
    /*
    prepare job file in GAMESS/FIREFLY format
    */
    fprintf(inp_fd.handle,
      " $CONTRL UNITS=ANGS COORD=UNIQUE RUNTYP=ENERGY SCFTYP=%cHF\n"
      " ICHARG=%d MULT=1 MAXIT=199 ", od->field.spin,
      (int)safe_rint(total_formal_charge));
    if (!strcmp(od->field.theory, O3_DFT_B3LYP)) {
      fprintf(inp_fd.handle, "DFTTYP=B3LYP ");
    }
    fprintf(inp_fd.handle,
      "$END\n"
      " $SCF DIRSCF=.TRUE. ");
    if (od->field.diff_sp) {
      fprintf(inp_fd.handle, "FDIFF=.FALSE. ");
    }
    fprintf(inp_fd.handle,
      "$END\n"
      " $SYSTEM TIMLIM=999999 MWORDS=24 $END\n"
      " $GUESS GUESS=HUCKEL $END\n");
    if ((od->field.type & PREP_FIREFLY_INPUT)
      || (od->field.type & PREP_GAMESS_INPUT)) {
      if (od->field.type & QM_ELE_FIELD) {
        fprintf(inp_fd.handle,
          " $ELPOT IEPOT=1 %s $END\n",
          (od->field.type & PREP_GAMESS_INPUT)
          ? "WHERE=GRID OUTPUT=PUNCH" : "");
      }
      else {
        fprintf(inp_fd.handle,
          " $ELDENS IEDEN=1 %s $END\n",
          (od->field.type & PREP_GAMESS_INPUT)
          ? "WHERE=GRID OUTPUT=PUNCH" : "");
      }
    }
    if (!(od->field.basis_set & EMSL_BASIS_SET)) {
      fprintf(inp_fd.handle, " $BASIS GBASIS=");
      switch (od->field.basis_set) {
        case O3_STO_3G:
        fprintf(inp_fd.handle, "STO NGAUSS=3 ");
        break;

        case O3_3_21G:
        fprintf(inp_fd.handle, "N21 NGAUSS=3 ");
        break;

        case O3_6_31G:
        fprintf(inp_fd.handle, "N31 NGAUSS=6 ");
        break;

        case O3_6_311G:
        fprintf(inp_fd.handle, "N311 NGAUSS=6 ");
        break;
      }
      if (od->field.basis_set != O3_STO_3G) {
        if (od->field.d_func) {
          fprintf(inp_fd.handle, "NDFUNC=%d ", od->field.d_func);
        }
        if (od->field.p_func) {
          fprintf(inp_fd.handle, "NPFUNC=%d ", od->field.p_func);
        }
        if (od->field.f_func) {
          fprintf(inp_fd.handle, "NFFUNC=%d ", od->field.f_func);
        }
        if (od->field.diff_sp) {
          fprintf(inp_fd.handle, "DIFFSP=.TRUE. ");
        }
      }
      fprintf(inp_fd.handle, "$END\n");
    }
    fprintf(inp_fd.handle,
      " $DATA\n"
      "%04d_energy\n"
      "C1\n", od->al.mol_info[object_num]->object_id);
  }
  else if (od->field.type & PREP_TURBOMOLE_INPUT) {
    /*
    prepare job file in TURBOMOLE format
    */
    if (!(od->field.basis_set & EMSL_BASIS_SET)) {
      switch (od->field.basis_set) {
        case O3_STO_3G:
        strcpy(basis_set_name, "sto-3g hondo");
        break;

        case O3_SV:
        strcpy(basis_set_name, "SV");
        break;

        case O3_SVP:
        strcpy(basis_set_name, "SVP");
        break;

        case O3_TZVP:
        strcpy(basis_set_name, "TZVP");
        break;
      }
    }
    fprintf(inp_fd.handle,
      "\n"
      "%04d_energy\n"
      "a coord\n"
      "*\n"
      "no\n"
      "b all %s\n",
      od->al.mol_info[object_num]->object_id,
      basis_set_name);
    sprintf(coord_fd.name, "%s%ccoord", input_dir, SEPARATOR);
    if (!(coord_fd.handle = fopen(coord_fd.name, "wb+"))) {
      O3_ERROR_LOCATE(task);
      O3_ERROR_STRING(task, coord_fd.name);
      return FL_CANNOT_WRITE_INP_FILE;
    }
    fprintf(coord_fd.handle, "$coord\n");
    if (od->field.basis_set & EMSL_BASIS_SET) {
      sprintf(basis_fd.name, "%s%cemsl_%04d.txt", input_dir,
        SEPARATOR, od->al.mol_info[object_num]->object_id);
      if (!(basis_fd.handle = fopen(basis_fd.name, "wb+"))) {
        O3_ERROR_LOCATE(task);
        O3_ERROR_STRING(task, basis_fd.name);
        return FL_CANNOT_WRITE_INP_FILE;
      }
      fprintf(basis_fd.handle,
        "$basis\n"
        "*\n");
    }
  }
  memset(used_elem_array, 0xFF, BUF_LEN);
  for (n = 0; n < n_atoms; ++n) {
    elem_i = 0;
    while (element_data[(int)elem_i].atomic_number) {
      if (!strcasecmp(atom[n]->element, element_data[(int)elem_i].atomic_symbol)) {
        break;
      }
      ++elem_i;
    }
    if (!(element_data[(int)elem_i].atomic_number)) {
      O3_ERROR_LOCATE(task);
      return FL_UNKNOWN_ATOM_TYPE;
    }
    for (i = 0, found = 0; ((used_elem_array[i] != 0xFF) && (!found)); ++i) {
      found = (used_elem_array[i] == elem_i);
    }
    if (!found) {
      used_elem_array[i] = elem_i;
    }
    if (od->field.type & PREP_GAUSSIAN_INPUT) {
      fprintf(inp_fd.handle,
        "%-3s%16.5lf%16.5lf%16.5lf\n", atom[n]->element,
        atom[n]->coord[0], atom[n]->coord[1], atom[n]->coord[2]);
    }
    else if ((od->field.type & PREP_FIREFLY_INPUT)
      || (od->field.type & PREP_GAMESS_INPUT)) {
      fprintf(inp_fd.handle,
        "%-3s%7.1lf  %16.10lf%16.10lf%16.10lf\n", atom[n]->element,
        (double)(element_data[(int)elem_i].atomic_number),
        atom[n]->coord[0], atom[n]->coord[1], atom[n]->coord[2]);
      if (od->field.basis_set & EMSL_BASIS_SET) {
        fprintf(inp_fd.handle,
          "%s\n", element_data[(int)elem_i].gamess_basis_set[basis_set_i]);
      }
    }
    else if (od->field.type & PREP_TURBOMOLE_INPUT) {
      strcpy(lowercase_elem, atom[n]->element);
      string_to_lowercase(lowercase_elem);
      fprintf(coord_fd.handle,
        "%20.14lf  %20.14lf  %20.14lf      %s\n",
        atom[n]->coord[0] / BOHR_RADIUS,
        atom[n]->coord[1] / BOHR_RADIUS,
        atom[n]->coord[2] / BOHR_RADIUS,
        lowercase_elem);
    }
  }
  if (od->field.type & PREP_GAUSSIAN_INPUT) {
    fprintf(inp_fd.handle, "\n");
    if (od->field.basis_set & EMSL_BASIS_SET) {
      for (i = 0; used_elem_array[i] != 0xFF; ++i) {
        fprintf(inp_fd.handle,
          "%s 0\n"
          "%s"
          "****\n",
          element_data[(int)(used_elem_array[i])].atomic_symbol,
          element_data[(int)(used_elem_array[i])].gaussian_basis_set[basis_set_i]);
      }
      fprintf(inp_fd.handle, "\n");
    }
    fprintf(inp_fd.handle,
      "!!! PLEASE DO NOT REMOVE/MODIFY THE FOLLOWING LINES !!!\n"
      "! "FORMCHK_EXE" %s_%04d.chk %s_%04d"FORMCHK_EXT"\n"
      "! "CUBEGEN_EXE" 0 %s=SCF %s_%04d"FORMCHK_EXT
      " %s_%04d"GAUSSIAN_CUBE_EXT" -1 h\n"
      "! %5d%12.6f%12.6f%12.6f\n"
      "! %5d%12.6f%12.6f%12.6f\n"
      "! %5d%12.6f%12.6f%12.6f\n"
      "! %5d%12.6f%12.6f%12.6f\n",
      od->field.qm_software,
      od->al.mol_info[object_num]->object_id,
      od->field.qm_software,
      od->al.mol_info[object_num]->object_id,
      ((od->field.type & QM_ELE_FIELD) ? "potential" : "density"),
      od->field.qm_software,
      od->al.mol_info[object_num]->object_id,
      od->field.qm_software,
      od->al.mol_info[object_num]->object_id,
      GAUSSIAN_UNIT_NUMBER,
      od->grid.start_coord[0],
      od->grid.start_coord[1],
      od->grid.start_coord[2],
      od->grid.nodes[0],
      od->grid.step[0], 0.0, 0.0,
      od->grid.nodes[1],
      0.0, od->grid.step[1], 0.0,
      od->grid.nodes[2],
      0.0, 0.0, od->grid.step[2]);
  }
  else if ((od->field.type & PREP_FIREFLY_INPUT)
    || (od->field.type & PREP_GAMESS_INPUT)) {
    fprintf(inp_fd.handle, " $END\n");
    if (od->field.type & PREP_FIREFLY_INPUT) {
      fprintf(inp_fd.handle,
        " $CUBE CUBE=.TRUE. NXCUBE=%d NYCUBE=%d NZCUBE=%d\n"
        " X0CUBE=%.4lf Y0CUBE=%.4lf Z0CUBE=%.4lf\n"
        " XINCUB=%.4lf YINCUB=%.4lf ZINCUB=%.4lf $END\n\n",
        od->grid.nodes[0], od->grid.nodes[1], od->grid.nodes[2],
        (double)(od->grid.start_coord[0]) / BOHR_RADIUS,
        (double)(od->grid.start_coord[1]) / BOHR_RADIUS,
        (double)(od->grid.start_coord[2]) / BOHR_RADIUS,
        (double)(od->grid.step[0]) / BOHR_RADIUS,
        (double)(od->grid.step[1]) / BOHR_RADIUS,
        (double)(od->grid.step[2]) / BOHR_RADIUS);
    }
    else {
      fprintf(inp_fd.handle,
        " $GRID MODGRD=1 SIZE=%.4lf UNITS=BOHR\n"
        " ORIGIN(1)=%.4lf ORIGIN(2)=%.4lf ORIGIN(3)=%.4lf\n"
        " XVEC(1)=%.4lf XVEC(2)=%.4lf XVEC(3)=%.4lf\n"
        " YVEC(1)=%.4lf YVEC(2)=%.4lf YVEC(3)=%.4lf\n"
        " ZVEC(1)=%.4lf ZVEC(2)=%.4lf ZVEC(3)=%.4lf $END\n\n",
        (double)(od->grid.step[0]) / BOHR_RADIUS,
        (double)(od->grid.start_coord[0]) / BOHR_RADIUS,
        (double)(od->grid.start_coord[1]) / BOHR_RADIUS,
        (double)(od->grid.start_coord[2]) / BOHR_RADIUS,
        (double)(od->grid.end_coord[0] + 0.3 * od->grid.step[0]) / BOHR_RADIUS,
        (double)(od->grid.start_coord[1]) / BOHR_RADIUS,
        (double)(od->grid.start_coord[2]) / BOHR_RADIUS,
        (double)(od->grid.start_coord[0]) / BOHR_RADIUS,
        (double)(od->grid.end_coord[1] + 0.3 * od->grid.step[0]) / BOHR_RADIUS,
        (double)(od->grid.start_coord[2]) / BOHR_RADIUS,
        (double)(od->grid.start_coord[0]) / BOHR_RADIUS,
        (double)(od->grid.start_coord[1]) / BOHR_RADIUS,
        (double)(od->grid.end_coord[2] + 0.3 * od->grid.step[0]) / BOHR_RADIUS);
    }
  }
  else if (od->field.type & PREP_TURBOMOLE_INPUT) {
    fprintf(coord_fd.handle,
      "$user-defined bonds\n"
      "$end\n");
    fclose(coord_fd.handle);
    if (od->field.basis_set & EMSL_BASIS_SET) {
      for (i = 0; used_elem_array[i] != 0xFF; ++i) {
        fprintf(inp_fd.handle,
          "file\n"
          "%s\n",
          get_basename(basis_fd.name));
        strcpy(lowercase_elem, element_data[(int)(used_elem_array[i])].atomic_symbol);
        string_to_lowercase(lowercase_elem);
        fprintf(basis_fd.handle,
          "%s %s\n"
          "# %s     %s\n"
          "*\n"
          "%s"
          "*\n",
          lowercase_elem, basis_set_name,
          lowercase_elem, basis_set_name,
          element_data[(int)(used_elem_array[i])].turbomole_basis_set[basis_set_i]);
      }
      fprintf(basis_fd.handle, "$end\n");
      fclose(basis_fd.handle);
    }
    fprintf(inp_fd.handle,
      "*\n"
      "eht\n"
      "y\n"
      "%d\n"
      "%s",
      (int)safe_rint(total_formal_charge),
      ((od->field.spin == O3_UNRESTRICTED)
      ? "n\n"
      "uf 0\n"
      "*\n"
      "n\n"
      : "y\n"));
    if (!strcmp(od->field.theory, O3_DFT_B3LYP)) {
      fprintf(inp_fd.handle,
      "dft\n"
      "on\n"
      "func b3-lyp_Gaussian\n"
      "*\n");
    }
    fprintf(inp_fd.handle,
      "*\n");
    rewind(inp_fd.handle);
    memset(&prog_exe_info, 0, sizeof(ProgExeInfo));
    prog_exe_info.proc_env = fill_env
      (od, turbomole_env, od->field.qm_exe_path, object_num);
    if (!(prog_exe_info.proc_env)) {
      O3_ERROR_LOCATE(task);
      return FL_OUT_OF_MEMORY;
    }
    prog_exe_info.need_stdin = NEED_STDIN_NORMAL;
    prog_exe_info.stdout_fd = &log_fd;
    prog_exe_info.stderr_fd = &log_fd;
    prog_exe_info.exedir = input_dir;
    prog_exe_info.sep_proc_grp = 1;
    sprintf(log_fd.name, "%s%cdefine_%04d.log", input_dir,
      SEPARATOR, od->al.mol_info[object_num]->object_id);
    sprintf(prog_exe_info.command_line, "%s%c%s",
      od->field.qm_exe_path, SEPARATOR, TURBOMOLE_DEFINE_EXE);
    pid = ext_program_exe(&prog_exe_info, &error);
    if (error) {
      O3_ERROR_LOCATE(task);
      return error;
    }
    #ifndef WIN32
    if (!(pipe_handle = fdopen(prog_exe_info.pipe_des[1], "w"))) {
      O3_ERROR_LOCATE(task);
      return FL_CANNOT_CREATE_CHANNELS;
    }
    #else
    pipe_handle = prog_exe_info.stdin_wr;
    #endif
    /*
    pipe job information into define
    */
    while (fgets(buffer, BUF_LEN, inp_fd.handle)) {
      buffer[BUF_LEN - 1] = '\0';
      FWRITE_WRAP(pipe_handle, buffer, &n_chr);
    }
    fclose(inp_fd.handle);
    FFLUSH_WRAP(pipe_handle);
    FCLOSE_WRAP(pipe_handle);
    ext_program_wait(&prog_exe_info, pid);
    if (!(log_fd.handle = fopen(log_fd.name, "rb"))) {
      O3_ERROR_LOCATE(task);
      O3_ERROR_STRING(task, log_fd.name);
      return FL_CANNOT_READ_TEMP_FILE;
    }
    if (!fgrep(log_fd.handle, buffer, TURBOMOLE_NORMAL_TERMINATION)) {
      O3_ERROR_LOCATE(task);
      error = FL_ABNORMAL_TERMINATION;
    }
    fclose(log_fd.handle);
    if (error) {
      return error;
    }
    sprintf(inp_fd.name, "%s%ccontrol", input_dir, SEPARATOR);
    if (!(inp_fd.handle = fopen(inp_fd.name, "rb+"))) {
      O3_ERROR_LOCATE(task);
      O3_ERROR_STRING(task, inp_fd.name);
      return FL_CANNOT_READ_TEMP_FILE;
    }
    found = 0;
    while (!found) {
      pos = ftell(inp_fd.handle);
      if (!fgets(buffer, BUF_LEN, inp_fd.handle)) {
        break;
      }
      buffer[BUF_LEN - 1] = '\0';
      remove_newline(buffer);
      found = (!strncmp(buffer, "$scfiterlimit", 13));
    }
    if (!found) {
      fclose(inp_fd.handle);
      O3_ERROR_LOCATE(task);
      O3_ERROR_STRING(task, inp_fd.name);
      return FL_CANNOT_READ_TEMP_FILE;
    }
    if (fseek(inp_fd.handle, pos, SEEK_SET)) {
      O3_ERROR_LOCATE(task);
      O3_ERROR_STRING(task, inp_fd.name);
      return FL_CANNOT_READ_TEMP_FILE;
    }
    n = strlen(buffer) - 13;
    sprintf(buffer, "$scfiterlimit%%%dd\n", n);
    fprintf(inp_fd.handle, buffer, MAX_TURBOMOLE_SCFITER);
    found = 0;
    while (!found) {
      pos = ftell(inp_fd.handle);
      if (!fgets(buffer, BUF_LEN, inp_fd.handle)) {
        break;
      }
      buffer[BUF_LEN - 1] = '\0';
      found = (!strncmp(buffer, "$last step", 10));
    }
    if (!found) {
      fclose(inp_fd.handle);
      O3_ERROR_LOCATE(task);
      O3_ERROR_STRING(task, inp_fd.name);
      return FL_CANNOT_READ_TEMP_FILE;
    }
    if (fseek(inp_fd.handle, pos, SEEK_SET)) {
      O3_ERROR_LOCATE(task);
      O3_ERROR_STRING(task, inp_fd.name);
      return FL_CANNOT_READ_TEMP_FILE;
    }
    fprintf(inp_fd.handle,
      "$pointval %s fmt=cub\n"
      "   grid1 vector 1 0 0 range %.6lf,%.6lf points %d\n"
      "   grid2 vector 0 1 0 range %.6lf,%.6lf points %d\n"
      "   grid3 vector 0 0 1 range %.6lf,%.6lf points %d\n"
      "   origin 0.0 0.0 0.0\n"
      "$end\n",
      ((od->field.type & QM_ELE_FIELD) ? "pot" : ""),
      (double)(od->grid.start_coord[0]) / BOHR_RADIUS,
      (double)(od->grid.end_coord[0]) / BOHR_RADIUS,
      od->grid.nodes[0],
      (double)(od->grid.start_coord[1]) / BOHR_RADIUS,
      (double)(od->grid.end_coord[1]) / BOHR_RADIUS,
      od->grid.nodes[1],
      (double)(od->grid.start_coord[2]) / BOHR_RADIUS,
      (double)(od->grid.end_coord[2]) / BOHR_RADIUS,
      od->grid.nodes[2]);
  }
  if (inp_fd.handle) {
    fclose(inp_fd.handle);
    inp_fd.handle = NULL;
  }

  return 0;
}


int prep_molden_input(O3Data *od, int object_num)
{
  char buffer[BUF_LEN];
  int i;
  int result;
  double c[3];
  double ed[3];
  FileDescriptor temp_fd;
  
  
  memset(buffer, 0, BUF_LEN);
  memset(&temp_fd, 0, sizeof(FileDescriptor));
  sprintf(buffer, "%s%cmolden_%04d", od->field.qm_dir,
    SEPARATOR, od->al.mol_info[object_num]->object_id);
  #ifndef WIN32
  result = mkdir(buffer, S_IRWXU | S_IRGRP | S_IROTH);
  #else
  result = mkdir(buffer);
  #endif
  if (result == -1) {
    return FL_CANNOT_CREATE_MDNDIR;
  }
  sprintf(temp_fd.name,
    "%s%cmolden_%04d"MOLDEN_INP_EXT, buffer,
      SEPARATOR, od->al.mol_info[object_num]->object_id);
  if (!(temp_fd.handle =
    fopen(temp_fd.name, "wb+"))) {
    O3_ERROR_LOCATE(&(od->task));
    O3_ERROR_STRING(&(od->task), temp_fd.name);
    return FL_CANNOT_WRITE_INP_FILE;
  }
  for (i = 0; i < 3; ++i) {
    c[i] = (double)(od->grid.start_coord[i]
      + od->grid.end_coord[i])
      / (2.0 * BOHR_RADIUS);
    ed[i] = (double)(od->grid.end_coord[i]
      - od->grid.start_coord[i])
      / BOHR_RADIUS;
  }
  fprintf(temp_fd.handle,
    "molden_%04d %s\n"
    "FILE=..%cmolden_%04d.out SPACE=0.01 %sWRBAS\n"
    "CENTER=(%.6lf,%.6lf,%.6lf) LINE=(0.0,0.0,-1.0)\n"
    "EDX=%.6lf EDY=%.6lf EDZ=%.6lf "
    "NPTSX=%d NPTSY=%d NPTSZ=%d\n",
    od->al.mol_info[object_num]->object_id,
    ((od->field.type & QM_ELE_FIELD) ? "esp" : "density"),
    SEPARATOR, od->al.mol_info[object_num]->object_id,
    ((od->field.type & QM_ELE_FIELD) ? "ELPO " : ""),
    c[0], c[1], c[2],
    ed[0], ed[1], ed[2],
    od->grid.nodes[0],
    od->grid.nodes[1],
    od->grid.nodes[2]);
  fclose(temp_fd.handle);
  
  return 0;
}


void prep_sybyl_input(O3Data *od)
{
  char current_time[BUF_LEN];


  if (get_current_time(current_time)) {
    strcpy(current_time, "Date unavailable");
  }
  fprintf(od->file[PREPINP_OUT]->handle,
    "Filename:  %s\n"
    "Created.:  %s"
    "Points:  %d\n"
    "Boxes:  1\n"
    "\n"
    "  Box 1            _____ X _____  _____ Y _____  _____ Z _____\n"
    "  Lower Corner:    %13.6f  %13.6f  %13.6f\n"
    "  High Corner.:    %13.6f  %13.6f  %13.6f\n"
    "  Step Size...:    %13.6f  %13.6f  %13.6f\n"
    "  Number steps:%17d%17d%17d\n"
    "  Probe Atom..:  C.3\n"
    "  Charge......:  1.000000\n",
    od->file[PREPINP_OUT]->name,
    current_time,
    od->grid.nodes[0] * od->grid.nodes[1]
    * od->grid.nodes[2],
    od->grid.start_coord[0],
    od->grid.start_coord[1],
    od->grid.start_coord[2],
    od->grid.end_coord[0],
    od->grid.end_coord[1],
    od->grid.end_coord[2],
    od->grid.step[0],
    od->grid.step[1],
    od->grid.step[2],
    od->grid.nodes[0],
    od->grid.nodes[1],
    od->grid.nodes[2]);
}


void prep_moe_grid_input(O3Data *od)
{
  int i;
  int n;
  
  
  fprintf(od->file[PREPINP_OUT]->handle,
    "function o3q_define_grid_shape []\n"
    "\n"
    "  return\n"
    "  [\n");
  for (i = 0; i < 3; ++i) {
    fprintf(od->file[PREPINP_OUT]->handle, "    [ ");
    for (n = 0; n < od->grid.nodes[i]; ++n) {
      fprintf(od->file[PREPINP_OUT]->handle, "%.4lf%s",
        safe_rint((double)n
        * (double)(od->grid.step[i]) * 1.0e04) / 1.0e04
        + (double)(od->grid.start_coord[i]),
        ((n == (od->grid.nodes[i] - 1)) ? " ]" : ", "));
    }
    fprintf(od->file[PREPINP_OUT]->handle,
      ((i == 2) ? "\n" : ",\n"));
  }
  fprintf(od->file[PREPINP_OUT]->handle,
    "  ];\n"
    "\n"
    "endfunction\n");
}


int replace_coord(int sdf_version, char *buffer, double *coord)
{
  char buffer2[BUF_LEN];
  char *ptr = NULL;
  char *context = NULL;
  int i;
  int x;

  
  if (sdf_version == V2000) {
    sprintf(buffer, "%10.4f%10.4f%10.4f", coord[0], coord[1], coord[2]);
    buffer[30] = ' ';
  }
  else {
    memset(buffer2, 0, BUF_LEN);
    ptr = strtok_r(buffer, " \t\0", &context);
    i = 4;
    while (ptr && i) {
      strcat(buffer2, ptr);
      strcat(buffer2, ((i == 4) ? "  " : " "));
      --i;
      ptr = strtok_r(NULL, " \t\0", &context);
    }
    if (!ptr) {
      return FL_CANNOT_READ_MOL_FILE;
    }
    for (x = 0; x < 3; ++x) {
      i = strlen(buffer2);
      sprintf(&buffer2[i], "%.4f ", coord[x]);
    }
    i = 3;
    while (ptr && i) {
      --i;
      ptr = strtok_r(NULL, " \t\0", &context);
    }
    if (!ptr) {
      return FL_CANNOT_READ_MOL_FILE;
    }
    while (ptr) {
      strcat(buffer2, ptr);
      strcat(buffer2, " ");
      ptr = strtok_r(NULL, " \t\0", &context);
    }
    strcpy(buffer, buffer2);
    i = strlen(buffer);
    buffer[i - 1] = '\0';
  }
  
  return 0;
}


int find_conformation_in_sdf(FILE *handle_in, FILE *handle_out, int conf_num)
{
  char buffer[BUF_LEN];
  char *ptr = NULL;
  int i = 0;
  int n = 0;
  int skip_mol_name = 0;
  
  
  memset(buffer, 0, BUF_LEN);
  if (conf_num == SKIP_MOL_NAME) {
    skip_mol_name = 1;
    conf_num = 0;
  }
  while (n <= conf_num) {
    i = 0;
    while ((i < 4) && fgets(buffer, BUF_LEN, handle_in)) {
      buffer[BUF_LEN - 1] = '\0';
      remove_newline(buffer);
      if (handle_out && (n == conf_num)
        && ((!skip_mol_name) || (i == 3))) {
        fprintf(handle_out, "%s\n", buffer);
      }
      ++i;
    }
    if (i < 4) {
      return FL_CANNOT_READ_SDF_FILE;
    }
    if (n == conf_num) {
      break;
    }
    i = 0;
    while ((!i) && fgets(buffer, BUF_LEN, handle_in)) {
      buffer[BUF_LEN - 1] = '\0';
      i = (!strncmp(buffer, SDF_DELIMITER, 4));
    }
    ++n;
  }
  if (n != conf_num) {
    return FL_CANNOT_FIND_CONF;
  }
  /*
  move to the beginning of the atom coordinate section
  */
  if (!strncasecmp(&buffer[34], "V3000", 5)) {
    ptr = NULL;
    while ((!ptr) && fgets(buffer, BUF_LEN, handle_in)) {
      buffer[BUF_LEN - 1] = '\0';
      remove_newline(buffer);
      if (handle_out) {
        fprintf(handle_out, "%s\n", buffer);
      }
      ptr = strstr(buffer, "BEGIN ATOM");
    }
    if (!ptr) {
      return FL_CANNOT_READ_SDF_FILE;
    }
  }
  
  return 0;
}


void prepare_rototrans_matrix(double *rt_mat, double *t_mat1, double *t_mat2, double *rad)
{
  double rt_mat1[RT_MAT_SIZE];
  double rt_mat2[RT_MAT_SIZE];


  /*
  objects are translated to the origin, rotated along x, y, z axes
  then re-translated to their original position; the transformation
  is not applied object-wise but dataset-wise
  */
  memset(rt_mat, 0, RT_MAT_SIZE * sizeof(double));
  rt_mat[0]                   = cos(rad[2]);
  rt_mat[1]                   = sin(rad[2]);
  rt_mat[RT_VEC_SIZE]         = -sin(rad[2]);
  rt_mat[RT_VEC_SIZE + 1]     = cos(rad[2]);
  rt_mat[RT_VEC_SIZE * 2 + 2] = 1.0;
  rt_mat[RT_VEC_SIZE * 3 + 3] = 1.0;
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
    RT_VEC_SIZE, RT_VEC_SIZE, RT_VEC_SIZE,
    1.0, t_mat2, RT_VEC_SIZE, rt_mat, RT_VEC_SIZE,
    0.0, rt_mat1, RT_VEC_SIZE);
  memset(rt_mat2, 0, RT_MAT_SIZE * sizeof(double));
  rt_mat2[0]                   = cos(rad[1]);
  rt_mat2[2]                   = -sin(rad[1]);
  rt_mat2[RT_VEC_SIZE + 1]     = 1.0;
  rt_mat2[RT_VEC_SIZE * 2]     = sin(rad[1]);
  rt_mat2[RT_VEC_SIZE * 2 + 2] = cos(rad[1]);
  rt_mat2[RT_VEC_SIZE * 3 + 3] = 1.0;
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
    RT_VEC_SIZE, RT_VEC_SIZE, RT_VEC_SIZE,
    1.0, rt_mat1, RT_VEC_SIZE, rt_mat2, RT_VEC_SIZE,
    0.0, rt_mat, RT_VEC_SIZE);
  memset(rt_mat2, 0, RT_MAT_SIZE * sizeof(double));
  rt_mat2[0]                   = 1.0;
  rt_mat2[RT_VEC_SIZE + 1]     = cos(rad[0]);
  rt_mat2[RT_VEC_SIZE + 2]     = sin(rad[0]);
  rt_mat2[RT_VEC_SIZE * 2 + 1] = -sin(rad[0]);
  rt_mat2[RT_VEC_SIZE * 2 + 2] = cos(rad[0]);
  rt_mat2[RT_VEC_SIZE * 3 + 3] = 1.0;
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
    RT_VEC_SIZE, RT_VEC_SIZE, RT_VEC_SIZE,
    1.0, rt_mat, RT_VEC_SIZE, rt_mat2, RT_VEC_SIZE,
    0.0, rt_mat1, RT_VEC_SIZE);
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
    RT_VEC_SIZE, RT_VEC_SIZE, RT_VEC_SIZE,
    1.0, rt_mat1, RT_VEC_SIZE, t_mat1, RT_VEC_SIZE,
    0.0, rt_mat, RT_VEC_SIZE);
}


int rototrans(O3Data *od, char *out_sdf_name, double *trans, double *rot)
{
  char buffer[BUF_LEN];
  int i;
  int x;
  int object_num;
  int n_atoms = 0;
  int result;
  double centroid[3];
  double t_mat1[RT_MAT_SIZE];
  double t_mat2[RT_MAT_SIZE];
  double rt_mat[RT_MAT_SIZE];
  double t_vec1[RT_VEC_SIZE];
  double t_vec2[RT_VEC_SIZE];
  FileDescriptor mol_fd;
  FileDescriptor out_sdf_fd;
  
  
  memset(buffer, 0, BUF_LEN);
  memset(&mol_fd, 0, sizeof(FileDescriptor));
  memset(&out_sdf_fd, 0, sizeof(FileDescriptor));
  strcpy(out_sdf_fd.name, out_sdf_name);
  if (!(out_sdf_fd.handle = fopen(out_sdf_fd.name, "wb"))) {
    O3_ERROR_LOCATE(&(od->task));
    O3_ERROR_STRING(&(od->task), out_sdf_fd.name);
    return CANNOT_WRITE_ROTOTRANSED_SDF;
  }
  memset(centroid, 0, 3 * sizeof(double));
  for (object_num = 0; object_num < od->grid.object_num; ++object_num) {
    if (!(od->al.mol_info[object_num]->atom = (AtomInfo **)
      alloc_array(od->al.mol_info[object_num]->n_atoms + 1, sizeof(AtomInfo)))) {
      O3_ERROR_LOCATE(&(od->task));
      free_atom_array(od);
      return OUT_OF_MEMORY;
    }
    result = fill_atom_info(od, &(od->task),
      od->al.mol_info[object_num]->atom, NULL, object_num, O3_MMFF94);
    if (result) {
      free_atom_array(od);
      return result;
    }
    /*
    add to centroid
    */
    if (get_object_attr(od, object_num, OPERATE_BIT)) {
      for (x = 0; x < 3; ++x) {
        for (i = 0; i < od->al.mol_info[object_num]->n_atoms; ++i) {
          centroid[x] += od->al.mol_info[object_num]->atom[i]->coord[x];
        }
      }
      n_atoms += od->al.mol_info[object_num]->n_atoms;
    }
  }
  if (n_atoms) {
    for (x = 0; x < 3; ++x) {
      centroid[x] /= (double)n_atoms;
    }
  }
  for (object_num = 0; object_num < od->grid.object_num; ++object_num) {
    sprintf(mol_fd.name, "%s%c%04d.mol", od->field.mol_dir,
      SEPARATOR, od->al.mol_info[object_num]->object_id);
    if (!(mol_fd.handle = fopen(mol_fd.name, "rb"))) {
      O3_ERROR_LOCATE(&(od->task));
      O3_ERROR_STRING(&(od->task), mol_fd.name);
      fclose(out_sdf_fd.handle);
      free_atom_array(od);
      return CANNOT_READ_ORIGINAL_SDF;
    }
    if (find_conformation_in_sdf(mol_fd.handle, out_sdf_fd.handle, 0)) {
      O3_ERROR_LOCATE(&(od->task));
      O3_ERROR_STRING(&(od->task), mol_fd.name);
      fclose(mol_fd.handle);
      fclose(out_sdf_fd.handle);
      free_atom_array(od);
      return CANNOT_READ_ORIGINAL_SDF;
    }
    if (get_object_attr(od, object_num, OPERATE_BIT)) {
      /*
      initialize t_mat1 and t_mat2
      */
      memset(t_mat1, 0, RT_MAT_SIZE * sizeof(double));
      for (i = 0; i < RT_MAT_SIZE; i += 5) {
        t_mat1[i] = 1.0;
      }
      memcpy(t_mat2, t_mat1, RT_MAT_SIZE * sizeof(double));
      /*
      translate molecule to origin
      */
      cblas_daxpy(3, -1.0, centroid, 1, &t_mat1[3 * RT_VEC_SIZE], 1);
      /*
      translate molecule by the desired amount
      */
      for (x = 0; x < 3; ++x) {
        t_mat2[3 * RT_VEC_SIZE + x] = centroid[x] + trans[x];
      }
      memset(t_vec1, 0, RT_VEC_SIZE * sizeof(double));
      t_vec1[3] = 1.0;
      /*
      rototranslate molecule by the desired amount
      */
      prepare_rototrans_matrix(rt_mat, t_mat1, t_mat2, rot);
    }
    for (i = 0; i < od->al.mol_info[object_num]->n_atoms; ++i) {
      if (get_object_attr(od, object_num, OPERATE_BIT)) {
        cblas_dcopy(3, od->al.mol_info[object_num]->atom[i]->coord, 1, t_vec1, 1);
        cblas_dgemv(CblasColMajor, CblasNoTrans,
          RT_VEC_SIZE, RT_VEC_SIZE, 1.0, rt_mat, RT_VEC_SIZE,
          t_vec1, 1, 0.0, t_vec2, 1);
      }
      if (!fgets(buffer, BUF_LEN, mol_fd.handle)) {
        O3_ERROR_LOCATE(&(od->task));
        O3_ERROR_STRING(&(od->task), mol_fd.name);
        fclose(mol_fd.handle);
        fclose(out_sdf_fd.handle);
        free_atom_array(od);
        return CANNOT_READ_ORIGINAL_SDF;
      }
      buffer[BUF_LEN - 1] = '\0';
      remove_newline(buffer);
      if (get_object_attr(od, object_num, OPERATE_BIT)) {
        if (replace_coord(od->al.mol_info[object_num]->sdf_version, buffer, t_vec2)) {
          O3_ERROR_LOCATE(&(od->task));
          O3_ERROR_STRING(&(od->task), mol_fd.name);
          fclose(mol_fd.handle);
          fclose(out_sdf_fd.handle);
          free_atom_array(od);
          return CANNOT_READ_ORIGINAL_SDF;
        }
      }
      fprintf(out_sdf_fd.handle, "%s\n", buffer);
    }
    while (fgets(buffer, BUF_LEN, mol_fd.handle)) {
      buffer[BUF_LEN - 1] = '\0';
      remove_newline(buffer);
      fprintf(out_sdf_fd.handle, "%s\n", buffer);
    }
    fprintf(out_sdf_fd.handle, SDF_DELIMITER"\n");
    fclose(mol_fd.handle);
  }
  free_atom_array(od);
  fclose(out_sdf_fd.handle);

  return 0;
}
