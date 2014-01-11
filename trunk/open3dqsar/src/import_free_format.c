/*

import_free_format.c

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


int check_duplicate_parameter(int **max_vary, int *parameter)
{
  int i = 0;
  int found = 0;
  
  
  while ((!found) && max_vary[i]) {
    found = (parameter == max_vary[i]);
    ++i;
  }
  
  return found;
}


int find_vary_speed(O3Data *od, char *name_list, int **max_vary, int **vary,
  int *field_num, int *object_num, VarCoord *varcoord)
{
  char buffer[BUF_LEN];
  char *ptr;
  char *context = NULL;
  int i;
  int result;
  
  
  memset(buffer, 0, BUF_LEN);
  memset(max_vary, 0, (MAX_FREE_FORMAT_PARAMETERS + 1) * sizeof(int *));
  memset(vary, 0, (MAX_FREE_FORMAT_PARAMETERS + 1) * sizeof(int *));
  strcpy(buffer, name_list);
  i = 0;
  ptr = strtok_r(buffer, ",\0", &context);
  while ((i < MAX_FREE_FORMAT_PARAMETERS) && ptr) {
    if (!strncasecmp(ptr, "N_FIELD", 7)) {
      if (!(result = check_duplicate_parameter(max_vary, &(od->field_num)))) {
        max_vary[i] = &(od->field_num);
        vary[i] = field_num;
      }
    }
    else if (!strncasecmp(ptr, "N_OBJECT", 8)) {
      if (!(result = check_duplicate_parameter(max_vary, &(od->object_num)))) {
        max_vary[i] = &(od->object_num);
        vary[i] = object_num;
      }
    }
    else if (!strncasecmp(ptr, "X_COORD", 7)) {
      if (!(result = check_duplicate_parameter(max_vary, &(od->grid.nodes[0])))) {
        max_vary[i] = &(od->grid.nodes[0]);
        vary[i] = &(varcoord->node[0]);
      }
    }
    else if (!strncasecmp(ptr, "Y_COORD", 7)) {
      if (!(result = check_duplicate_parameter(max_vary, &(od->grid.nodes[1])))) {
        max_vary[i] = &(od->grid.nodes[1]);
        vary[i] = &(varcoord->node[1]);
      }
    }
    else if (!strncasecmp(ptr, "Z_COORD", 7)) {
      if (!(result = check_duplicate_parameter(max_vary, &(od->grid.nodes[2])))) {
        max_vary[i] = &(od->grid.nodes[2]);
        vary[i] = &(varcoord->node[2]);
      }
    }
    else {
      return WRONG_PARAMETER_NAME;
    }
    if (result) {
      return DUPLICATE_PARAMETER_NAME;
    }
    ++i;
    ptr = strtok_r(NULL, ",\0", &context);
  }
  
  return ((i < MAX_FREE_FORMAT_PARAMETERS) ? NOT_ENOUGH_PARAMETERS : 0);
}


int get_double_from_ascii_file(FILE *handle, char *read_buffer, int read_buffer_size,
  char **context, char **ptr, int *back, double *value)
{
  int actual_len = 0;
  int loop = 1;
  
  
  while (loop) {
    if (!(*context)) {
      if (*back) {
        memmove(read_buffer, &read_buffer[read_buffer_size - 1 - (*back)], *back);
      }
      if (!(actual_len = fzread(&read_buffer[*back], 1,
        read_buffer_size - 1 - (*back), (fzPtr *)handle))) {
        return EOF_FOUND;
      }
      if (actual_len == (read_buffer_size - 1 - (*back))) {
        *back = 0;
        actual_len = read_buffer_size - 2;
        while (actual_len && read_buffer[actual_len]
          && (!strchr("\t ,;\r\n", (int)read_buffer[actual_len]))) {
          --actual_len;
          ++(*back);
        }
        read_buffer[actual_len] = '\0';
      }
      else {
        read_buffer[actual_len + (*back)] = '\0';
      }
      *ptr = strtok_r(read_buffer, "\t ,;\r\n", context);
    }
    else {
      *ptr = strtok_r(NULL, "\t ,;\r\n", context);
    }
    if (*ptr) {
      sscanf(*ptr, "%lf", value);
      loop = 0;
    }
    else {
      *context = NULL;
    }
  }
  
  return 0;
}


int import_free_format(O3Data *od, char *name_list, int skip_header, int *n_values)
{
  char buffer[LARGE_BUF_LEN];
  char *context = NULL;
  char *ptr = NULL;
  int i;
  int num_fields;
  int field_num = 0;
  int object_num = 0;
  int result = 0;
  int back = 0;
  int *vary[MAX_FREE_FORMAT_PARAMETERS + 1];
  int *max_vary[MAX_FREE_FORMAT_PARAMETERS + 1];
  double value;
  VarCoord varcoord;
  
  
  memset(buffer, 0, LARGE_BUF_LEN);
  memset(&varcoord, 0, sizeof(VarCoord));
  memset(vary, 0, (MAX_FREE_FORMAT_PARAMETERS + 1) * sizeof(int));
  memset(max_vary, 0, (MAX_FREE_FORMAT_PARAMETERS + 1) * sizeof(int));
  *n_values = 0;
  i = skip_header;
  while (!result) {
    if (!(result = get_double_from_ascii_file(od->file[ASCII_IN]->handle,
      buffer, LARGE_BUF_LEN, &context, &ptr, &back, &value))) {
      if (i) {
        --i;
      }
      else {
        ++(*n_values);
      }
    }
  }
  if (result != EOF_FOUND) {
    O3_ERROR_LOCATE(&(od->task));
    O3_ERROR_STRING(&(od->task), od->file[ASCII_IN]->name);
    return result;
  }
  if (*n_values < (od->object_num * od->x_vars)) {
    O3_ERROR_LOCATE(&(od->task));
    O3_ERROR_STRING(&(od->task), od->file[ASCII_IN]->name);
    return NOT_ENOUGH_VALUES;
  }
  if (*n_values % (od->object_num * od->x_vars)) {
    O3_ERROR_LOCATE(&(od->task));
    O3_ERROR_STRING(&(od->task), od->file[ASCII_IN]->name);
    return INCORRECT_NUMBER_OF_VALUES;
  }
  num_fields = *n_values / (od->object_num * od->x_vars);
  fzrewind((fzPtr *)(od->file[ASCII_IN]->handle));
  ptr = NULL;
  context = NULL;
  result = alloc_x_var_array(od, num_fields);
  if (result) {
    O3_ERROR_LOCATE(&(od->task));
    return result;
  }
  find_vary_speed(od, name_list, max_vary, vary, &field_num, &object_num, &varcoord);
  result = 0;
  i = skip_header;
  while (i) {
    result = get_double_from_ascii_file(od->file[ASCII_IN]->handle,
      buffer, LARGE_BUF_LEN, &context, &ptr, &back, &value);
    if (result) {
      return result;
    }
    --i;
  }
  for (*vary[4] = 0; (!result) && (*vary[4] < *max_vary[4]); ++(*vary[4])) {
    for (*vary[3] = 0; (!result) && (*vary[3] < *max_vary[3]); ++(*vary[3])) {
      for (*vary[2] = 0; (!result) && (*vary[2] < *max_vary[2]); ++(*vary[2])) {
        for (*vary[1] = 0; (!result) && (*vary[1] < *max_vary[1]); ++(*vary[1])) {
          for (*vary[0] = 0; (!result) && (*vary[0] < *max_vary[0]); ++(*vary[0])) {
            result = get_double_from_ascii_file(od->file[ASCII_IN]->handle,
              buffer, LARGE_BUF_LEN, &context, &ptr, &back, &value);
            if (!result) {
              result = set_x_value_xyz_unbuffered
                (od, field_num, object_num, &varcoord, value);
            }
          }
        }
      }
    }
  }
  if (od->save_ram) {
    sync_field_mmap(od);
  }
  if (!result) {
    update_field_object_attr(od, VERBOSE_BIT);
    result = calc_active_vars(od, FULL_MODEL);
  }
  else {
    O3_ERROR_LOCATE(&(od->task));
  }
  
  return result;
}
