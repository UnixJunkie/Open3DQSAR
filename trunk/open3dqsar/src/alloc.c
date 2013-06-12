/*

alloc.c

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


char **alloc_array(int n, int size)
{
  char **array;
  int i;
  
  
  if (!(array = (char **)malloc((n + 1) * sizeof(char *)))) {
    return NULL;
  }
  memset(array, 0, (n + 1) * sizeof(char *));
  for (i = 0; i < n; ++i) {
    if (!(array[i] = (char *)malloc(size))) {
      free_array(array);
      return NULL;
    }
    memset(array[i], 0, size);
  }

  return array;
}


CharMat *alloc_char_matrix(CharMat *old_char_mat, int m, int n)
{
  int i;
  int j;
  CharMat *char_mat;


  char_mat = old_char_mat;
  if (!char_mat) {
    /*
    if the CharMat is being created from scratch
    */
    char_mat = (CharMat *)malloc(sizeof(CharMat));
    if (!char_mat) {
      return NULL;
    }
    memset(char_mat, 0, sizeof(CharMat));
  }
  if (m < char_mat->m) {
    /*
    if the new array of pointers is shorter
    than the previous one, free the old pointers
    */
    for (i = (char_mat->m - 1); i >= m; --i) {
      if (char_mat->me[i]) {
        free(char_mat->me[i]);
      }
    }
  }
  /*
  allocate the array of pointers
  */
  char_mat->me = (char **)realloc(char_mat->me, m * sizeof(char *));
  if (!(char_mat->me)) {
    return NULL;
  }
  if (m > char_mat->m) {
    /*
    if the new array of pointers is longer
    than the previous one, initialize to 0 all new pointers
    */
    memset(&(char_mat->me[char_mat->m]), 0,
      (m - char_mat->m) * sizeof(char *));
  }
  /*
  now allocate the arrays of chars
  */
  for (i = 0; i < m; ++i) {
    if ((i >= char_mat->m) || (n > char_mat->n)) {
      char_mat->me[i] = (char *)realloc(char_mat->me[i], n);
      if (!(char_mat->me[i])) {
        /*
        if allocation fails, free the arrays of chars
        already allocated, as well as the
        array of pointers
        */
        for (j = 0; j < i; ++j) {
          if (char_mat->me[j]) {
            free(char_mat->me[j]);
            char_mat->me[j] = NULL;
          }
          if (char_mat->me) {
            free(char_mat->me);
            char_mat->me = NULL;
          }
          return NULL;
        }
      }
      /*
      if needed, initialize to zero the arrays of chars
      */
      if (i >= char_mat->m) {
        memset(char_mat->me[i], 0, n);
      }
      else if (n > char_mat->n) {
        memset(&char_mat->me[i][char_mat->n], 0, n - char_mat->n);
      }
    }
  }
  /*
  update CharMat information
  */
  char_mat->m = m;
  char_mat->n = n;
  
  return char_mat;
}


int alloc_file_descriptor(O3Data *od, int file_num)
{
  int i;
  
  
  od->mel.file_descriptor = (FileDescriptor **)
    realloc(od->mel.file_descriptor, 
    sizeof(FileDescriptor *) * (od->file_num + file_num));
  if (!(od->mel.file_descriptor)) {
    return OUT_OF_MEMORY;
  }
  memset(&(od->mel.file_descriptor[od->file_num]), 0,
    sizeof(FileDescriptor *) * file_num);
  od->file = od->mel.file_descriptor;
  for (i = 0; i < file_num; ++i) {
    od->mel.file_descriptor[od->file_num + i] =
      (FileDescriptor *)malloc(sizeof(FileDescriptor));
    if (!(od->mel.file_descriptor[od->file_num + i])) {
      return OUT_OF_MEMORY;
    }
    memset(od->mel.file_descriptor[od->file_num + i],
      0, sizeof(FileDescriptor));
  }
  od->file_num += file_num;
  
  return 0;
}


int *alloc_int_array(int *old_ptr, int places)
{
  int *membuf;

  
  membuf = (int *)realloc(old_ptr, sizeof(int) * places);
  if (!membuf) {
    return NULL;
  }
  memset(membuf, 0, sizeof(int) * places);

  return membuf;
}


IntMat *alloc_int_matrix(IntMat *old_int_mat, int m, int n)
{
  int i;
  int j;
  IntMat *int_mat;


  int_mat = old_int_mat;
  if (!int_mat) {
    /*
    if the IntMat is being created from scratch
    */
    int_mat = (IntMat *)malloc(sizeof(IntMat));
    if (!int_mat) {
      return NULL;
    }
    memset(int_mat, 0, sizeof(IntMat));
  }
  if (m < int_mat->m) {
    /*
    if the new array of pointers is shorter
    than the previous one, free the old pointers
    */
    for (i = (int_mat->m - 1); i >= m; --i) {
      if (int_mat->me[i]) {
        free(int_mat->me[i]);
      }
    }
  }
  /*
  allocate the array of pointers
  */
  int_mat->me = (int **)realloc(int_mat->me, m * sizeof(int *));
  if (!(int_mat->me)) {
    return NULL;
  }
  if (m > int_mat->m) {
    /*
    if the new array of pointers is longer
    than the previous one, initialize to 0 all new pointers
    */
    memset(&(int_mat->me[int_mat->m]), 0,
      (m - int_mat->m) * sizeof(int *));
  }
  /*
  now allocate the arrays of ints
  */
  for (i = 0; i < m; ++i) {
    if ((i >= int_mat->m) || (n > int_mat->n)) {
      int_mat->me[i] = (int *)realloc(int_mat->me[i], n * sizeof(int));
      if (!(int_mat->me[i])) {
        /*
        if allocation fails, free the arrays of ints
        already allocated, as well as the
        array of pointers
        */
        for (j = 0; j < i; ++j) {
          if (int_mat->me[j]) {
            free(int_mat->me[j]);
            int_mat->me[j] = NULL;
          }
        }
        if (int_mat->me) {
          free(int_mat->me);
          int_mat->me = NULL;
        }
        return NULL;
      }
      /*
      if needed, initialize to zero the arrays of ints
      */
      if (i >= int_mat->m) {
        memset(int_mat->me[i], 0, n * sizeof(int));
      }
      else if (n > int_mat->n) {
        memset(&int_mat->me[i][int_mat->n], 0, (n - int_mat->n) * sizeof(int));
      }
    }
  }
  /*
  update IntMat information
  */
  int_mat->m = m;
  int_mat->n = n;
  
  return int_mat;
}


int alloc_threads(O3Data *od)
{
  int i;


  memset(od->mel.thread_info, 0, MAX_THREADS * sizeof(ThreadInfo *));
  for (i = 0; i < od->n_proc; ++i) {
    od->mel.thread_info[i] = (ThreadInfo *)malloc(sizeof(ThreadInfo));
    if (!(od->mel.thread_info[i])) {
      return OUT_OF_MEMORY;
    }
    memset(od->mel.thread_info[i], 0, sizeof(ThreadInfo));
  }
  #ifndef WIN32
  od->mel.mutex = (pthread_mutex_t *)realloc
    (od->mel.mutex, sizeof(pthread_mutex_t));
  #else
  od->mel.mutex = (HANDLE *)realloc
    (od->mel.mutex, sizeof(HANDLE));
  #endif
  if (!(od->mel.mutex)) {
    return OUT_OF_MEMORY;
  }

  return 0;
}


int alloc_x_var_array(O3Data *od, int num_fields)
{
  char buffer[LARGE_BUF_LEN];
  char field_dir_name[BUF_LEN];
  int i;
  int j;
  int n;
  int len;
  int actual_len;
  
  
  memset(buffer, 0, LARGE_BUF_LEN);
  memset(field_dir_name, 0, BUF_LEN);
  od->valid &= (SDF_BIT | COSMOTHERM_BIT);
  od->mel.x_var_array = (float ***)
    realloc(od->mel.x_var_array,
    sizeof(float *) * (od->field_num + num_fields));
  if (!(od->mel.x_var_array)) {
    return OUT_OF_MEMORY;
  }
  memset(&(od->mel.x_var_array[od->field_num]), 0,
    sizeof(float *) * num_fields);
  for (i = 0; i < num_fields; ++i) {
    od->mel.x_var_array[od->field_num + i] =
      (float **)malloc(sizeof(float *) * od->object_num);
    if (!(od->mel.x_var_array[od->field_num + i])) {
      return OUT_OF_MEMORY;
    }
    memset(od->mel.x_var_array[od->field_num + i], 0,
      sizeof(float *) * od->object_num);
  }
  od->mel.x_var_attr = (unsigned short **)realloc
    (od->mel.x_var_attr, sizeof(unsigned short *)
    * (od->field_num + num_fields));
  if (!(od->mel.x_var_attr)) {
    return OUT_OF_MEMORY;
  }
  memset(&(od->mel.x_var_attr[od->field_num]), 0,
    sizeof(unsigned short *) * num_fields);
  for (i = 0; i < MAX_VAR_BUF; ++i) {
    od->mel.x_var_buf[i] = (double **)realloc
    (od->mel.x_var_buf[i], sizeof(double *)
      * (od->field_num + num_fields));
    if (!(od->mel.x_var_buf[i])) {
      return OUT_OF_MEMORY;
    }
    memset(&(od->mel.x_var_buf[i][od->field_num]), 0,
      sizeof(double *) * num_fields);
  }
  for (i = 0; i < num_fields; ++i) {
    if (od->save_ram) {
      if (!i) {
        if (alloc_file_descriptor(od, num_fields)) {
          return OUT_OF_MEMORY;
        }
        if (open_temp_dir(od, NULL, "field_dir", field_dir_name)) {
          return CANNOT_WRITE_TEMP_FILE;
        }
      }
      if (!(od->file[TEMP_FIELD_DATA + od->field_num + i]->handle)) {
        sprintf(od->file[TEMP_FIELD_DATA + od->field_num + i]->name,
          "%s%cfield-%02d_data", field_dir_name, SEPARATOR,
          od->field_num + i + 1);
        if (!(od->file[TEMP_FIELD_DATA + od->field_num + i]->handle = fopen
          (od->file[TEMP_FIELD_DATA + od->field_num + i]->name, "wb+"))) {
          return CANNOT_WRITE_TEMP_FILE;
        }
      }
      rewind(od->file[TEMP_FIELD_DATA + od->field_num + i]->handle);
      n = od->object_pagesize * od->object_num / LARGE_BUF_LEN + 1;
      for (j = 0; j < n; ++j) {
        len = ((j == (n - 1))
          ? od->object_pagesize * od->object_num % LARGE_BUF_LEN : LARGE_BUF_LEN);
        if (len) {
          actual_len = fwrite(buffer, 1, len,
            od->file[TEMP_FIELD_DATA + od->field_num + i]->handle);
          if (actual_len != len) {
            return CANNOT_WRITE_TEMP_FILE;
          }
        }
      }
      rewind(od->file[TEMP_FIELD_DATA + od->field_num + i]->handle);
    }
    else {
      /*
      reserve RAM for the x variable array
      */
      for (j = 0; j < od->object_num; ++j) {
        od->mel.x_var_array[od->field_num + i][j] =
          (float *)malloc(sizeof(float) * od->x_vars);
        if (!(od->mel.x_var_array[od->field_num + i][j])) {
          return OUT_OF_MEMORY;
        }
      }
    }
    od->mel.x_var_attr[od->field_num + i] =
      (unsigned short *)malloc(od->x_vars * sizeof(unsigned short));
    if (!(od->mel.x_var_attr[od->field_num + i])) {
      return OUT_OF_MEMORY;
    }
    for (j = 0; j < od->x_vars; ++j) {
      od->mel.x_var_attr[od->field_num + i][j] = ACTIVE_BIT;
    }
    for (j = 0; j < MAX_VAR_BUF; ++j) {
      od->mel.x_var_buf[j][od->field_num + i] =
        (double *)malloc(sizeof(double) * od->x_vars);
      if (!(od->mel.x_var_buf[j][od->field_num + i])) {
        return OUT_OF_MEMORY;
      }
      memset(od->mel.x_var_buf[j][od->field_num + i], 0,
        sizeof(double) * od->x_vars);
    }
  }
  
  od->mel.field_attr = (unsigned short *)realloc
    (od->mel.field_attr,
    (od->field_num + num_fields)
    * sizeof(unsigned short));
  if (!(od->mel.field_attr)) {
    return OUT_OF_MEMORY;
  }
  for (i = 0; i < num_fields; ++i) {
    od->mel.field_attr[od->field_num + i] = ACTIVE_BIT;
  }
  od->mel.x_data = (XData *)realloc(od->mel.x_data,
    (od->field_num + num_fields) * sizeof(XData));
  if (!(od->mel.x_data)) {
    return OUT_OF_MEMORY;
  }
  memset(&(od->mel.x_data[od->field_num]), 0,
    num_fields * sizeof(XData));
  for (i = 0; i < num_fields; ++i) {
    od->mel.x_data
      [od->field_num + i].x_weight_coefficient = 1.0;
    od->mel.x_data
      [od->field_num + i].min_cutoff = -MAX_CUTOFF;
    od->mel.x_data
      [od->field_num + i].max_cutoff = MAX_CUTOFF;
  }
  
  od->field_num += num_fields;
  od->active_field_num += num_fields;

  return 0;
}


int alloc_object_attr(O3Data *od)
{
  int i;
  
  
  od->mel.object_attr =
    (unsigned short *)realloc(od->mel.object_attr,
    od->grid.object_num * sizeof(unsigned short));
  if (!(od->mel.object_attr)) {
    return OUT_OF_MEMORY;
  }
  od->mel.object_weight =
    (double *)realloc(od->mel.object_weight,
    od->grid.object_num * sizeof(double));
  if (!(od->mel.object_attr)) {
    return OUT_OF_MEMORY;
  }
  for (i = 0; i < od->grid.object_num; ++i) {
    od->mel.object_attr[i] = ACTIVE_BIT;
    od->mel.object_weight[i] = 1.0;
  }
  
  return 0;
}


int alloc_y_var_array(O3Data *od)
{
  int i;
  
  
  /*
  reserve RAM for the y variable array
  */
  od->mel.y_var_array =
    (float *)realloc(od->mel.y_var_array,
    sizeof(float) * od->grid.object_num
    * od->y_vars);
  if (!(od->mel.y_var_array)) {
    return OUT_OF_MEMORY;
  }
  memset(od->mel.y_var_array, 0,
    sizeof(float) * od->grid.object_num * od->y_vars);
  od->mel.y_var_attr =
    realloc(od->mel.y_var_attr,
    od->y_vars * sizeof(unsigned short));
  if (!(od->mel.y_var_attr)) {
    return OUT_OF_MEMORY;
  }
  memset(od->mel.y_var_attr, ACTIVE_BIT,
    od->y_vars * sizeof(unsigned short));
  for (i = 0; i < MAX_VAR_BUF; ++i) {
    od->mel.y_var_buf[i] =
      (double *)realloc(od->mel.y_var_buf[i],
      sizeof(double) * od->y_vars);
    if (!(od->mel.y_var_buf[i])) {
      return OUT_OF_MEMORY;
    }
  }
  od->mel.y_data = (YData *)realloc
    (od->mel.y_data, od->y_vars * sizeof(YData));
  if (!(od->mel.y_data)) {
    return OUT_OF_MEMORY;
  }
  memset(od->mel.y_data, 0, od->y_vars * sizeof(YData));
  for (i = 0; i < od->y_vars; ++i) {
    od->mel.y_data[i].y_weight_coefficient = 1.0;
  }

  return 0;
}


int fill_thread_info(O3Data *od, int n_tasks)
{
  int i;
  int n_threads;
  int n_calc_per_thread;
  int exceeding;
  ThreadInfo **ti;


  ti = od->mel.thread_info;
  n_calc_per_thread = n_tasks / od->n_proc;
  exceeding = n_tasks % od->n_proc;
  n_threads = ((n_tasks < od->n_proc) ? n_tasks : od->n_proc);
  for (i = 0; i < n_threads; ++i) {
    memcpy(&(ti[i]->od), od, sizeof(O3Data));
    ti[i]->thread_num = i;
    ti[i]->n_calc = n_calc_per_thread;
    if (exceeding) {
      ++(ti[i]->n_calc);
      --exceeding;
    }
    if (i) {
      ti[i]->start = ti[i - 1]->end + 1;
    }
    else {
      ti[i]->start = 0;
    }
    ti[i]->end = ti[i]->start + ti[i]->n_calc - 1;
  }
  
  return n_threads;
}


int alloc_average_mat(O3Data *od, int model_type,
  int cv_type, int groups, int runs)
{
  int i;
  int run_count;
  int group_num;
  int struct_count;
  int object_num = 0;
  int object_num2 = 0;
  int struct_num = 0;
  int struct_num2 = 0;
  int active_struct_num;
  int conf_num;
  int conf_num2;
  int n_conf = 0;
  double sumweight = 0.0;


  /*
  allocate matrices to store average values to be used
  by trim_mean_center_[xy]_matrix_hp
  */
  od->mal.large_e_mat_ave = double_mat_resize(od->mal.large_e_mat_ave,
    od->cv.overall_cv_runs, od->mal.large_e_mat->n);
  if (!(od->mal.large_e_mat_ave)) {
    return OUT_OF_MEMORY;
  }
  memset(od->mal.large_e_mat_ave->base, 0,
    od->mal.large_e_mat_ave->m * od->mal.large_e_mat_ave->n * sizeof(double));
  if (model_type != SCRAMBLE_CV_MODEL) {
    od->mal.large_f_mat_ave = double_mat_resize(od->mal.large_f_mat_ave,
      od->cv.overall_cv_runs, od->mal.large_f_mat->n);
    if (!(od->mal.large_f_mat_ave)) {
      return OUT_OF_MEMORY;
    }
    memset(od->mal.large_f_mat_ave->base, 0,
      od->mal.large_f_mat_ave->m * od->mal.large_f_mat_ave->n * sizeof(double));
  }
  run_count = 0;
  switch (cv_type) {
    case LEAVE_ONE_OUT:
    int_perm_resize(od->pel.out_structs, 1);
    while (object_num < od->object_num) {
      struct_num = od->al.mol_info[object_num]->struct_num;
      conf_num = ((od->valid & COSMOTHERM_BIT)
        ? od->al.cosmo_list[struct_num]->n_conf[BOUND] : 1);
      for (n_conf = 0, sumweight = 0.0;
        n_conf < conf_num; ++n_conf, ++object_num) {
        if (get_object_attr(od, object_num, ACTIVE_BIT)) {
          sumweight += od->mel.object_weight[object_num];
        }
      }
      if (sumweight > 0.0) {
        /*
        one struct at a time is going to be left out
        */
        od->pel.out_structs->pe[0] = struct_num;
        /*
        precalculate averages for large_e_mat
        and large_f_mat, one for each LOO CV model
        */
        calc_large_mat_ave(od, od->mal.large_e_mat,
          od->mal.large_e_mat_ave, run_count);
        if (model_type != SCRAMBLE_CV_MODEL) {
          calc_large_mat_ave(od, od->mal.large_f_mat,
            od->mal.large_f_mat_ave, run_count);
        }
        ++run_count;
      }
    }
    break;

    case LEAVE_TWO_OUT:
    int_perm_resize(od->pel.out_structs, 2);
    get_attr_struct_ave(od, 0, ACTIVE_BIT, &active_struct_num, NULL);
    while (object_num < od->object_num) {
      struct_num = od->al.mol_info[object_num]->struct_num;
      conf_num = ((od->valid & COSMOTHERM_BIT)
        ? od->al.cosmo_list[struct_num]->n_conf[BOUND] : 1);
      for (n_conf = 0, sumweight = 0.0;
        n_conf < conf_num; ++n_conf, ++object_num) {
        if (get_object_attr(od, object_num, ACTIVE_BIT)) {
          sumweight += od->mel.object_weight[object_num];
        }
      }
      if (sumweight > 0.0) {
        object_num2 = object_num;
        while (object_num2 < od->object_num) {
          struct_num2 = od->al.mol_info[object_num2]->struct_num;
          conf_num2 = ((od->valid & COSMOTHERM_BIT)
            ? od->al.cosmo_list[struct_num2]->n_conf[BOUND] : 1);
          for (n_conf = 0, sumweight = 0.0;
            n_conf < conf_num2; ++n_conf, ++object_num2) {
            if (get_object_attr(od, object_num2, ACTIVE_BIT)) {
              sumweight += od->mel.object_weight[object_num2];
            }
          }
          if (sumweight > 0.0) {
            /*
            two objects at a time are going to be left out
            */
            od->pel.out_structs->pe[0] = struct_num;
            od->pel.out_structs->pe[1] = struct_num2;
            /*
            precalculate averages for large_e_mat
            and large_f_mat, one for each LTO CV model
            */
            calc_large_mat_ave(od, od->mal.large_e_mat,
              od->mal.large_e_mat_ave, run_count);
            if (model_type != SCRAMBLE_CV_MODEL) {
              calc_large_mat_ave(od, od->mal.large_f_mat,
                od->mal.large_f_mat_ave, run_count);
            }
            ++run_count;
          }
        }
      }
    }
    break;

    case LEAVE_MANY_OUT:
    for (i = 0; i < runs; ++i) {
      for (group_num = 0; group_num < groups; ++group_num) {
        /*
        initialize the left-out objects vector before PLS
        */
        int_perm_resize(od->pel.out_structs, od->mel.struct_per_group[group_num]);
        for (struct_count = 0; struct_count
          < od->mel.struct_per_group[group_num]; ++struct_count) {
          od->pel.out_structs->pe[struct_count] =
            od->cimal.group_composition_list[i]->me
            [group_num][struct_count];
        }
        qsort(od->pel.out_structs->pe,
          od->pel.out_structs->size, sizeof(int),
          compare_integers);
        /*
        precalculate averages for large_e_mat
        and large_f_mat, one for each LMO CV model
        */
        calc_large_mat_ave(od, od->mal.large_e_mat,
          od->mal.large_e_mat_ave, run_count);
        if (model_type != SCRAMBLE_CV_MODEL) {
          calc_large_mat_ave(od, od->mal.large_f_mat,
            od->mal.large_f_mat_ave, run_count);
        }
        ++run_count;
      }
    }
    break;
  }
  
  return 0;
}


int alloc_cv_sdep(O3Data *od, int pc_num, int runs)
{
  od->mel.struct_list = alloc_int_array
    (od->mel.struct_list, od->grid.struct_num);
  if (!(od->mel.struct_list)) {
    return OUT_OF_MEMORY;
  }
  od->mal.sdep_mat = double_mat_resize
    (od->mal.sdep_mat, runs, pc_num + 1);
  if (!(od->mal.sdep_mat)) {
    return OUT_OF_MEMORY;
  }
  od->mal.ave_press = double_mat_resize
    (od->mal.ave_press, pc_num + 1, od->y_vars);
  if (!(od->mal.ave_press)) {
    return OUT_OF_MEMORY;
  }
  memset(od->mal.ave_press->base, 0,
    od->mal.ave_press->m * od->mal.ave_press->n * sizeof(double));
  od->vel.sumweight = double_vec_resize(od->vel.sumweight, pc_num + 1);
  if (!(od->vel.sumweight)) {
    return OUT_OF_MEMORY;
  }
  memset(od->vel.sumweight->ve, 0, od->vel.sumweight->size * sizeof(double));
  od->vel.sum1 = double_vec_resize(od->vel.sum1, pc_num + 1);
  if (!(od->vel.sum1)) {
    return OUT_OF_MEMORY;
  }
  memset(od->vel.sum1->ve, 0, od->vel.sum1->size * sizeof(double));
  od->vel.q2 = double_vec_resize(od->vel.q2, pc_num + 1);
  if (!(od->vel.q2)) {
    return OUT_OF_MEMORY;
  }
  od->vel.secv = double_vec_resize(od->vel.secv, pc_num + 1);
  if (!(od->vel.secv)) {
    return OUT_OF_MEMORY;
  }
  od->pel.sdep_rank = int_perm_resize(od->pel.sdep_rank, pc_num + 1);
  if (!(od->pel.sdep_rank)) {
    return OUT_OF_MEMORY;
  }

  return 0;
}


int alloc_pls(O3Data *od, int x_vars, int pc_num, int model_type)
{
  int coeff_size;
  
  
  coeff_size = 1;
  if (model_type & (UVEPLS_FULL_MODEL | UVEPLS_CV_MODEL)
    && (!(od->uvepls.ive))) {
    coeff_size = 2;
  }
  od->mal.e_mat = double_mat_resize(od->mal.e_mat,
    od->object_num, x_vars * coeff_size);
  if (!(od->mal.e_mat)) {
    return OUT_OF_MEMORY;
  }
  od->vel.e_mat_ave = double_vec_resize
    (od->vel.e_mat_ave, x_vars * coeff_size);
  if (!(od->vel.e_mat_ave)) {
    return OUT_OF_MEMORY;
  }
  if (od->y_vars) {
    od->mal.f_mat = double_mat_resize
      (od->mal.f_mat, od->object_num, od->y_vars);
    if (!(od->mal.f_mat)) {
      return OUT_OF_MEMORY;
    }
    od->vel.f_mat_ave = double_vec_resize
      (od->vel.f_mat_ave, od->y_vars);
    if (!(od->vel.f_mat_ave)) {
      return OUT_OF_MEMORY;
    }
    od->mal.pred_f_mat = double_mat_resize
      (od->mal.pred_f_mat, od->object_num, od->y_vars);
    if (!(od->mal.pred_f_mat)) {
      return OUT_OF_MEMORY;
    }
    od->mal.y_scores = double_mat_resize
      (od->mal.y_scores, od->object_num, pc_num + 1);
    if (!(od->mal.y_scores)) {
      return OUT_OF_MEMORY;
    }
    od->mal.y_loadings = double_mat_resize
      (od->mal.y_loadings, od->y_vars, pc_num + 1);
    if (!(od->mal.y_loadings)) {
      return OUT_OF_MEMORY;
    }
    od->mal.b_coefficients = double_mat_resize
      (od->mal.b_coefficients, x_vars * coeff_size, od->y_vars);
    if (!(od->mal.b_coefficients)) {
      return OUT_OF_MEMORY;
    }
  }
  od->mal.x_scores = double_mat_resize
    (od->mal.x_scores, od->object_num, pc_num + 1);
  if (!(od->mal.x_scores)) {
    return OUT_OF_MEMORY;
  }
  od->mal.temp = double_mat_resize
    (od->mal.temp, pc_num + 1, pc_num + 1);
  if (!(od->mal.temp)) {
    return OUT_OF_MEMORY;
  }
  od->mel.ipiv = alloc_int_array(od->mel.ipiv, pc_num + 1);
  if (!(od->mel.ipiv)) {
    return OUT_OF_MEMORY;
  }
  #if (!defined HAVE_LIBLAPACK_ATLAS) && (!defined HAVE_LIBSUNPERF)
  od->mel.work = realloc(od->mel.work,
    (pc_num + 1) * LWORK_BLOCK_SIZE * sizeof(double));
  if (!(od->mel.work)) {
    return OUT_OF_MEMORY;
  }
  #endif
  od->mal.x_weights = double_mat_resize
    (od->mal.x_weights, x_vars * coeff_size, pc_num + 1);
  if (!(od->mal.x_weights)) {
    return OUT_OF_MEMORY;
  }
  od->mal.x_weights_star = double_mat_resize
    (od->mal.x_weights_star, x_vars * coeff_size, pc_num + 1);
  if (!(od->mal.x_weights_star)) {
    return OUT_OF_MEMORY;
  }
  od->mal.x_loadings = double_mat_resize
    (od->mal.x_loadings, x_vars * coeff_size, pc_num + 1);
  if (!(od->mal.x_loadings)) {
    return OUT_OF_MEMORY;
  }
  od->vel.v = double_vec_resize(od->vel.v, od->object_num);
  if (!(od->vel.v)) {
    return OUT_OF_MEMORY;
  }
  od->vel.v_new = double_vec_resize(od->vel.v_new, od->object_num);
  if (!(od->vel.v_new)) {
    return OUT_OF_MEMORY;
  }
  od->vel.ro = double_vec_resize(od->vel.ro, pc_num + 1);
  if (!(od->vel.ro)) {
    return OUT_OF_MEMORY;
  }
  od->vel.explained_s2_x = double_vec_resize(od->vel.explained_s2_x, pc_num + 1);
  if (!(od->vel.explained_s2_x)) {
    return OUT_OF_MEMORY;
  }
  od->vel.explained_s2_y = double_vec_resize(od->vel.explained_s2_y, pc_num + 1);
  if (!(od->vel.explained_s2_y)) {
    return OUT_OF_MEMORY;
  }
  od->vel.ave_sdep = double_vec_resize(od->vel.ave_sdep, pc_num + 1);
  if (!(od->vel.ave_sdep)) {
    return OUT_OF_MEMORY;
  }
  od->vel.ave_sdec = double_vec_resize(od->vel.ave_sdec, pc_num + 1);
  if (!(od->vel.ave_sdec)) {
    return OUT_OF_MEMORY;
  }
  od->vel.r2 = double_vec_resize(od->vel.r2, pc_num + 1);
  if (!(od->vel.r2)) {
    return OUT_OF_MEMORY;
  }
  od->vel.r2_pred = double_vec_resize(od->vel.r2_pred, pc_num + 1);
  if (!(od->vel.r2_pred)) {
    return OUT_OF_MEMORY;
  }
  memset(od->vel.ave_sdec->ve, 0, od->vel.ave_sdec->size * sizeof(double));
  od->pel.out_structs = int_perm_resize(od->pel.out_structs, od->grid.struct_num);
  if (!(od->pel.out_structs)) {
    return OUT_OF_MEMORY;
  }

  return 0;
}


int alloc_voronoi(O3Data *od, int seed_num)
{
  int x;
  int y;
  
  
  od->mel.seed_count =
    alloc_int_array(od->mel.seed_count, od->field_num);
  if (!(od->mel.seed_count)) {
    return OUT_OF_MEMORY;
  }
  od->mel.seed_count_before_collapse = alloc_int_array
    (od->mel.seed_count_before_collapse, od->field_num);
  if (!(od->mel.seed_count_before_collapse)) {
    return OUT_OF_MEMORY;
  }
  od->mel.group_zero = alloc_int_array
    (od->mel.group_zero, od->field_num);
  if (!(od->mel.group_zero)) {
    return OUT_OF_MEMORY;
  }
  /*
  allocate an array of pointers to arrays of integers
  we need as many pointers as seeds for the current field
  */
  od->al.voronoi_composition = (int **)realloc
    (od->al.voronoi_composition, (seed_num + 1) * sizeof(int *));
  if (!(od->al.voronoi_composition)) {
    return OUT_OF_MEMORY;
  }
  memset(od->al.voronoi_composition, 0, (seed_num + 1) * sizeof(int *));
  od->mel.voronoi_active = alloc_int_array
    (od->mel.voronoi_active, seed_num);
  if (!(od->mel.voronoi_active)) {
    return OUT_OF_MEMORY;
  }
  /*
  allocate an array which keeps track of the number of x_vars
  assigned to the polyhedron identified by each seed
  */
  od->mel.voronoi_fill = alloc_int_array
    (od->mel.voronoi_fill, seed_num);
  if (!(od->mel.voronoi_fill)) {
    return OUT_OF_MEMORY;
  }
  od->cimal.voronoi_buf = alloc_int_matrix
    (od->cimal.voronoi_buf, od->field_num, od->x_vars);
  if (!(od->cimal.voronoi_buf)) {
    return OUT_OF_MEMORY;
  }
  for (y = 0; y < od->field_num; ++y) {
    for (x = 0; x < od->x_vars; ++x) {
      set_voronoi_buf(od, y, x, -1);
    }
  }
  
  return 0;
}
