/*

import_dependent.c

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


int parse_o3_line(char *buffer)
{
  char *ptr = NULL;
  int i = 0;
  
  
  /*
  remove comments
  */
  if ((ptr = strchr(buffer, (int)'#'))) {
    *ptr = '\0';
  }
  /*
  skip blank lines
  */
  while (buffer[i] && isspace(buffer[i])) {
    ++i;
  }
  return (buffer[i] ? 1 : 0);
}


int import_dependent(O3Data *od, char *name_list)
{
  char buffer[BUF_LEN];
  char buffer_copy[BUF_LEN];
  char *ptr;
  char *context = NULL;
  int y_vars = 0;
  int y_vars_file;
  int lines;
  int i;
  int found = 0;
  int pos;
  int object_num;
  int struct_num;
  int result;
  double value;


  memset(buffer, 0, BUF_LEN);
  memset(buffer_copy, 0, BUF_LEN);
  /*
  there must be at least struct_num lines in addition
  to the first line, which is reserved for variable names
  */
  free_char_matrix(od->cimal.y_var_name);
  od->cimal.y_var_name = NULL;
  if (strcasecmp(name_list, "all")) {
    /*
    copy into buffer the list of variables
    requested by user
    */
    strcpy(buffer, name_list);
  }
  else {
    /*
    user requested ALL available variables,
    so just get the first line with labels
    */
    while (fgets(buffer_copy, BUF_LEN, od->file[DEP_IN]->handle)) {
      buffer_copy[BUF_LEN - 1] = '\0';
      if (parse_o3_line(buffer_copy)) {
        strcpy(buffer, buffer_copy);
        break;
      }
    }
    rewind(od->file[DEP_IN]->handle);
  }
  strcpy(buffer_copy, buffer);
  y_vars = 0;
  ptr = strtok_r(buffer, ", \t\r\n\0", &context);
  while (ptr) {
    for (i = 0, found = 0; ((!found) && (i < y_vars)); ++i) {
      found = (!strcmp(od->cimal.y_var_name->me[i], ptr));
    }
    if (!found) {
      ++y_vars;
      if (!(od->cimal.y_var_name = alloc_char_matrix
        (od->cimal.y_var_name, y_vars, MAX_NAME_LEN))) {
        return OUT_OF_MEMORY;
      }
      strcpy(od->cimal.y_var_name->me[y_vars - 1], ptr);
    }
    ptr = strtok_r(NULL, ", \t\r\n\0", &context);
  }
  strcpy(buffer, buffer_copy);
  if (!(od->pel.numberlist[Y_VAR_LIST] = int_perm_resize
    (od->pel.numberlist[Y_VAR_LIST], y_vars))) {
    return OUT_OF_MEMORY;
  }
  memset(od->pel.numberlist[Y_VAR_LIST]->pe, 0, y_vars * sizeof(int));
  od->y_vars = y_vars;
  /*
  get the first line with labels
  */
  while (fgets(buffer_copy, BUF_LEN, od->file[DEP_IN]->handle)) {
    buffer_copy[BUF_LEN - 1] = '\0';
    if (parse_o3_line(buffer_copy)) {
      strcpy(buffer, buffer_copy);
      break;
    }
  }
  y_vars_file = 0;
  y_vars = 0;
  ptr = strtok_r(buffer, " \t\n\r\0", &context);
  while (ptr) {
    ++y_vars_file;
    for (i = 0; i < od->pel.numberlist[Y_VAR_LIST]->size; ++i) {
      if (!strcmp(od->cimal.y_var_name->me[i], ptr)) {
        od->pel.numberlist[Y_VAR_LIST]->pe[y_vars] = y_vars_file;
        ++y_vars;
        break;
      }
    }
    ptr = strtok_r(NULL, " \t\n\r\0", &context);
  }
  for (i = 0; i < od->y_vars; ++i) {
    if (!(od->pel.numberlist[Y_VAR_LIST]->pe[i])) {
      return CANNOT_FIND_Y_VAR_NAME;
    }
  }
  lines = 0;
  pos = ftell(od->file[DEP_IN]->handle);
  while (fgets(buffer_copy, BUF_LEN, od->file[DEP_IN]->handle)) {
    buffer[BUF_LEN - 1] = '\0';
    if (parse_o3_line(buffer)) {
      ++lines;
    }
  }
  if (lines != od->grid.struct_num) {
    return NOT_ENOUGH_OBJECTS;
  }
  if (fseek(od->file[DEP_IN]->handle, pos, SEEK_SET)) {
    return PREMATURE_DEP_IN_EOF;
  }
  if (alloc_y_var_array(od)) {
    return OUT_OF_MEMORY;
  }
  for (struct_num = 0; struct_num < od->grid.struct_num; ++struct_num) {
    found = 0;
    while ((!found) && fgets(buffer, BUF_LEN, od->file[DEP_IN]->handle)) {
      buffer[BUF_LEN - 1] = '\0';
      found = parse_o3_line(buffer);
    }
    if (!found) {
      return PREMATURE_DEP_IN_EOF;
    }
    strcpy(buffer_copy, buffer);
    /*
    read variables; on each line the number of variables
    must be equal to the number of variable names which
    have been read on the first line
    */
    ptr = strtok_r(buffer, " \t\n\r\0", &context);
    y_vars = 0;
    while (ptr) {
      ++y_vars;
      ptr = strtok_r(NULL, " \t\n\r\0", &context);
    }
    if (y_vars != y_vars_file) {
      return WRONG_NUMBER_OF_Y_VARS;
    }
    ptr = strtok_r(buffer_copy, " \t\n\r\0", &context);
    y_vars = 0;
    i = 0;
    while (ptr && (y_vars <= y_vars_file)) {
      ++y_vars;
      if (od->pel.numberlist[Y_VAR_LIST]->pe[i] == y_vars) {
        /*
        record variables in RAM
        */
        sscanf(ptr, "%lf", &value);
        object_num = 0;
        while (object_num < od->grid.object_num) {
          if (od->al.mol_info[object_num]->struct_num == struct_num) {
            set_y_value(od, object_num, i, value);
          }
          ++object_num;
        }
        ++i;
      }
      ptr = strtok_r(NULL, " \t\n\r\0", &context);
    }
  }
  update_field_object_attr(od, VERBOSE_BIT);
  result = calc_active_vars(od, FULL_MODEL);

  return result;
}
