/*

cv_thread.c

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


#ifndef WIN32
void *lmo_cv_thread(void *pointer)
#else
DWORD lmo_cv_thread(void *pointer)
#endif
{
  int i;
  int j;
  int x;
  int num_predictions;
  int group_num;
  int object_count;
  int struct_count;
  int result;
  double cum_press;
  ThreadInfo *ti;


  ti = (ThreadInfo *)pointer;
  for (j = ti->start; j <= ti->end; ++j) {
    memset(ti->od.mal.press->base, 0,
      ti->od.mal.press->m
      * ti->od.mal.press->n
      * sizeof(double));
    num_predictions = 0;
    for (group_num = 0; group_num < ti->groups; ++group_num) {
      /*
      initialize the left-out objects vector before PLS
      */
      int_perm_resize(ti->od.pel.out_structs,
        ti->od.mel.struct_per_group[group_num]);
      for (struct_count = 0, object_count = 0;
        struct_count < ti->od.mel.struct_per_group[group_num];
        ++struct_count) {
        ti->od.pel.out_structs->pe[struct_count] =
          ti->od.cimal.group_composition_list[j]
          ->me[group_num][struct_count];
        ++object_count;
      }
      qsort(ti->od.pel.out_structs->pe,
        ti->od.pel.out_structs->size, sizeof(int), compare_integers);
      if (ti->model_type & UVEPLS_CV_MODEL) {
        trim_mean_center_x_matrix_hp(&(ti->od), ti->model_type,
          ti->od.active_object_num - object_count,
          j * ti->groups + group_num);
        trim_mean_center_y_matrix_hp(&(ti->od),
          ti->od.active_object_num - object_count,
          j * ti->groups + group_num);
      }
      else if (ti->model_type & SCRAMBLE_CV_MODEL) {
        trim_mean_center_x_matrix_hp(&(ti->od), ti->model_type,
          ti->od.active_object_num - object_count,
          j * ti->groups + group_num);
        trim_mean_center_matrix(&(ti->od), ti->od.mal.large_f_mat,
          &(ti->od.mal.f_mat), &(ti->od.vel.f_mat_ave),
          ti->model_type, ti->od.active_object_num - object_count);
      }
      else {
        trim_mean_center_matrix(&(ti->od), ti->od.mal.large_e_mat,
          &(ti->od.mal.e_mat), &(ti->od.vel.e_mat_ave),
          ti->model_type, ti->od.active_object_num - object_count);
        trim_mean_center_matrix(&(ti->od), ti->od.mal.large_f_mat,
          &(ti->od.mal.f_mat), &(ti->od.vel.f_mat_ave),
          ti->model_type, ti->od.active_object_num - object_count);
      }
      pls(&(ti->od), ti->pc_num, ti->model_type);
      result = pred_y_values(&(ti->od), ti,
        ti->od.pc_num, ti->model_type,
        j * ti->groups + group_num);
      if (result) {
        ti->cannot_write_temp_file = 1;
      }
      num_predictions += (ti->od.y_vars
        * ti->od.mel.struct_per_group[group_num]);
    }
    for (i = 0; i <= ti->od.pc_num; ++i) {
      cum_press = (double)0;
      #ifndef WIN32
      pthread_mutex_lock(ti->od.mel.mutex);
      #else
      WaitForSingleObject(ti->od.mel.mutex, INFINITE);
      #endif
      for (x = 0; x < ti->od.y_vars; ++x) {
        M_POKE(ti->od.mal.ave_press, i, x,
          M_PEEK(ti->od.mal.ave_press, i, x)
          + M_PEEK(ti->od.mal.press, i, x));
        cum_press += M_PEEK(ti->od.mal.press, i, x);

      }
      M_POKE(ti->od.mal.sdep_mat, j, i,
        sqrt(cum_press / (double)num_predictions));
      #ifndef WIN32
      pthread_mutex_unlock(ti->od.mel.mutex);
      #else
      ReleaseMutex(ti->od.mel.mutex);
      #endif
    }
  }
  #ifndef WIN32
  pthread_exit(pointer);
  #else
  return 0;
  #endif
}


#ifndef WIN32
void *loo_cv_thread(void *pointer)
#else
DWORD loo_cv_thread(void *pointer)
#endif
{
  int i;
  int object_count;
  int result;
  ThreadInfo *ti;


  ti = (ThreadInfo *)pointer;
  for (i = ti->start; i <= ti->end; ++i) {
    /*
    initialize the left-out objects vector before PLS
    */
    int_perm_resize(ti->od.pel.out_structs, 1);
    ti->od.pel.out_structs->pe[0] =
      ti->od.mel.struct_list[i];
    object_count = 1;
    if (ti->model_type & UVEPLS_CV_MODEL) {
      trim_mean_center_x_matrix_hp(&(ti->od), ti->model_type,
        ti->od.active_object_num - object_count, i);
      trim_mean_center_y_matrix_hp(&(ti->od),
        ti->od.active_object_num - object_count, i);
    }
    else if (ti->model_type & SCRAMBLE_CV_MODEL) {
      trim_mean_center_x_matrix_hp(&(ti->od), ti->model_type,
        ti->od.active_object_num - object_count, i);
      trim_mean_center_matrix(&(ti->od), ti->od.mal.large_f_mat,
        &(ti->od.mal.f_mat), &(ti->od.vel.f_mat_ave),
        ti->model_type, ti->od.active_object_num - object_count);
    }
    else {
      trim_mean_center_matrix(&(ti->od), ti->od.mal.large_e_mat,
        &(ti->od.mal.e_mat), &(ti->od.vel.e_mat_ave),
        ti->model_type, ti->od.active_object_num - object_count);
      trim_mean_center_matrix(&(ti->od), ti->od.mal.large_f_mat,
        &(ti->od.mal.f_mat), &(ti->od.vel.f_mat_ave),
        ti->model_type, ti->od.active_object_num - object_count);
    }
    pls(&(ti->od), ti->pc_num, ti->model_type);
    result = pred_y_values(&(ti->od), ti,
      ti->od.pc_num, ti->model_type, i);
    if (result) {
      ti->cannot_write_temp_file = 1;
    }
  }
  #ifndef WIN32
  pthread_exit(pointer);
  #else
  return 0;
  #endif
}


#ifndef WIN32
void *lto_cv_thread(void *pointer)
#else
DWORD lto_cv_thread(void *pointer)
#endif
{
  int i;
  int object_count;
  int result;
  ThreadInfo *ti;


  ti = (ThreadInfo *)pointer;
  for (i = (ti->start * 2); i <= (ti->end * 2); i += 2) {
    /*
    initialize the left-out objects vector before PLS
    */
    int_perm_resize(ti->od.pel.out_structs, 2);
    ti->od.pel.out_structs->pe[0] =
      ti->od.mel.struct_list[i];
    ti->od.pel.out_structs->pe[1] =
      ti->od.mel.struct_list[i + 1];
    object_count = 2;
    if (ti->model_type & UVEPLS_CV_MODEL) {
      trim_mean_center_x_matrix_hp(&(ti->od), ti->model_type,
        ti->od.active_object_num - object_count, i / 2);
      trim_mean_center_y_matrix_hp(&(ti->od),
        ti->od.active_object_num - object_count, i / 2);
    }
    if (ti->model_type & SCRAMBLE_CV_MODEL) {
      trim_mean_center_x_matrix_hp(&(ti->od), ti->model_type,
        ti->od.active_object_num - object_count, i / 2);
      trim_mean_center_matrix(&(ti->od), ti->od.mal.large_f_mat,
        &(ti->od.mal.f_mat), &(ti->od.vel.f_mat_ave),
        ti->model_type, ti->od.active_object_num - object_count);
    }
    else {
      trim_mean_center_matrix(&(ti->od), ti->od.mal.large_e_mat,
        &(ti->od.mal.e_mat), &(ti->od.vel.e_mat_ave),
        ti->model_type, ti->od.active_object_num - object_count);
      trim_mean_center_matrix(&(ti->od), ti->od.mal.large_f_mat,
        &(ti->od.mal.f_mat), &(ti->od.vel.f_mat_ave),
        ti->model_type, ti->od.active_object_num - object_count);
    }
    pls(&(ti->od), ti->pc_num, ti->model_type);
    result = pred_y_values(&(ti->od), ti,
      ti->od.pc_num, ti->model_type, i / 2);
    if (result) {
      ti->cannot_write_temp_file = 1;
    }
  }
  #ifndef WIN32
  pthread_exit(pointer);
  #else
  return 0;
  #endif
}
