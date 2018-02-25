/*

pred_y_values.c

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


int pred_y_values(O3Data *od, ThreadInfo *ti, int pc_num, int model_type, int cv_run)
{
  int i;
  int j;
  int x;
  int y;
  int actual_len;
  int object_num;
  int struct_num;
  int conf_num;
  int n_conf;
  int info;
  #if (!defined HAVE_LIBLAPACK_ATLAS) && (!defined HAVE_LIBSUNPERF)
  int lwork;
  #endif
  double value;
  double ave_res;
  double sqrt_weight;
  double sumweight;
  FileDescriptor *temp_pred;
  FileDescriptor *temp_cv_coeff;
  
  
  temp_pred = (ti ? &(ti->temp_pred) : od->file[TEMP_PRED]);
  temp_cv_coeff = (ti ? &(ti->temp_cv_coeff) : od->file[TEMP_CV_COEFF]);
  /*
  fill left_out_values with values for the
  left out struct
  */
  object_num = 0;
  i = 0;
  y = 0;
  while ((object_num < od->object_num)
    && (i < od->pel.out_structs->size)) {
    struct_num = od->al.mol_info[object_num]->struct_num;
    conf_num = 1;
    if (struct_num != od->pel.out_structs->pe[i]) {
      object_num += conf_num;
      continue;
    }
    for (n_conf = 0; n_conf < conf_num; ++n_conf, ++object_num) {
      if (!get_object_attr(od, object_num, ACTIVE_BIT)) {
        continue;
      }
      fill_y_vector(od, object_num, y, model_type, cv_run);
      fill_x_vector(od, object_num, y, model_type, cv_run);
      ++y;
    }
    ++i;
  }
  /*
  no need to check the return value, since this is just a logical
  resizing: all these matrices have been allocated
  with maximal size with this purpose
  */
  double_mat_resize(od->mal.f_mat, y, od->y_vars);
  double_mat_resize(od->mal.e_mat, y, od->mal.e_mat->n);
  double_mat_resize(od->mal.pred_f_mat, y, od->y_vars);
  if (model_type & (CV_MODEL | SILENT_PLS)) {
    /*
    write number of PCs
    */
    actual_len = fwrite(&pc_num, sizeof(int), 1, temp_pred->handle);
    if (actual_len != 1) {
      return CANNOT_WRITE_TEMP_FILE;
    }
    /*
    write number of y_vars
    */
    actual_len = fwrite(&(od->y_vars), sizeof(int), 1, temp_pred->handle);
    if (actual_len != 1) {
      return CANNOT_WRITE_TEMP_FILE;
    }
    /*
    write number of objects in this block
    */
    actual_len = fwrite(&y, sizeof(int), 1, temp_pred->handle);
    if (actual_len != 1) {
      return CANNOT_WRITE_TEMP_FILE;
    }
    /*
    write object numbers in this block
    */
    object_num = 0;
    i = 0;
    while ((object_num < od->object_num)
      && (i < od->pel.out_structs->size)) {
      struct_num = od->al.mol_info[object_num]->struct_num;
      conf_num = 1;
      if (struct_num != od->pel.out_structs->pe[i]) {
        object_num += conf_num;
        continue;
      }
      for (n_conf = 0; n_conf < conf_num; ++n_conf, ++object_num) {
        if (!get_object_attr(od, object_num, ACTIVE_BIT)) {
          continue;
        }
        actual_len = fwrite(&object_num, sizeof(int), 1, temp_pred->handle);
        if (actual_len != 1) {
          return CANNOT_WRITE_TEMP_FILE;
        }
      }
      ++i;
    }
  }
  /*
  predicted y values are stored in this way:
  pc_num
  number of variables (equal to x_max_y, that is y_vars)
  value_num(here, y_max), var_num(here, 0), var-0_value-0_pc-pc_num, var-0_value-1_pc-pc_num, var-0_value-2_pc-pc_num, ..., var-0_value-y_max_pc-pc_num
  value_num(here, y_max), var_num(here, 1), var-1_value-0_pc-pc_num, var-1_value-1_pc-pc_num, var-1_value-2_pc-pc_num, ..., var-1_value-y_max_pc-pc_num
  value_num(here, y_max), var_num(here, 2), var-2_value-0_pc-pc_num, var-2_value-1_pc-pc_num, var-2_value-2_pc-pc_num, ..., var-2_value-y_max_pc-pc_num
  ...
  value_num(here, y_max), var_num(here, x_max_y), var-x_max_y_value-0_pc-pc_num, var-x_max_y_value-1_pc-pc_num, var-x_max_y_value-2_pc-pc_num, ..., var-x_max_y_value-y_max_pc-pc_num

  value_num(here, y_max), var_num(here, 0), var-0_value-0_pc-(pc_num-1), var-0_value-1_pc-(pc_num-1), var-0_value-2_pc-(pc_num-1), ..., var-0_value-y_max_pc-(pc_num-1)
  value_num(here, y_max), var_num(here, 1), var-1_value-0_pc-(pc_num-1), var-1_value-1_pc-(pc_num-1), var-1_value-2_pc-(pc_num-1), ..., var-1_value-y_max_pc-(pc_num-1)
  value_num(here, y_max), var_num(here, 2), var-2_value-0_pc-(pc_num-1), var-2_value-1_pc-(pc_num-1), var-2_value-2_pc-(pc_num-1), ..., var-2_value-y_max_pc-(pc_num-1)
  ...
  value_num(here, y_max), var_num(here, x_max_y), var-x_max_y_value-0_pc-(pc_num-1), var-x_max_y_value-1_pc-(pc_num-1), var-x_max_y_value-2_pc-(pc_num-1), ..., var-x_max_y_value-y_max_pc-(pc_num-1)

  value_num(here, y_max), var_num(here, 0), var-0_value-0_pc-(pc_num-2), var-0_value-1_pc-(pc_num-2), var-0_value-2_pc-(pc_num-2), ..., var-0_value-y_max_pc-(pc_num-2)
  value_num(here, y_max), var_num(here, 1), var-1_value-0_pc-(pc_num-2), var-1_value-1_pc-(pc_num-2), var-1_value-2_pc-(pc_num-2), ..., var-1_value-y_max_pc-(pc_num-2)
  value_num(here, y_max), var_num(here, 2), var-2_value-0_pc-(pc_num-2), var-2_value-1_pc-(pc_num-2), var-2_value-2_pc-(pc_num-2), ..., var-2_value-y_max_pc-(pc_num-2)
  ...
  value_num(here, y_max), var_num(here, x_max_y), var-x_max_y_value-0_pc-(pc_num-2), var-x_max_y_value-1_pc-(pc_num-2), var-x_max_y_value-2_pc-(pc_num-2), ..., var-x_max_y_value-y_max_pc-(pc_num-2)

  ...
  ...
  ...

  value_num(here, y_max), var_num(here, 0), var-0_value-0_pc-0, var-0_value-1_pc-0, var-0_value-2_pc-0, ..., var-0_value-y_max_pc-0
  value_num(here, y_max), var_num(here, 1), var-1_value-0_pc-0, var-1_value-1_pc-0, var-1_value-2_pc-0, ..., var-1_value-y_max_pc-0
  value_num(here, y_max), var_num(here, 2), var-2_value-0_pc-0, var-2_value-1_pc-0, var-2_value-2_pc-0, ..., var-2_value-y_max_pc-0
  ...
  value_num(here, y_max), var_num(here, x_max_y), var-x_max_y_value-0_pc-0, var-x_max_y_value-1_pc-0, var-x_max_y_value-2_pc-0, ..., var-x_max_y_value-y_max_pc-0
  */    
  /*
  no need to check the return value, since this is just a logical
  resizing: all these matrices have been allocated
  with maximal size with this purpose
  */
  double_mat_resize(od->mal.x_weights, od->mal.e_mat->n, pc_num);
  double_mat_resize(od->mal.b_coefficients, od->mal.e_mat->n, od->y_vars);
  double_mat_resize(od->mal.x_loadings, od->mal.e_mat->n, pc_num);
  double_mat_resize(od->mal.temp, pc_num, pc_num);
  
  cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
    od->mal.x_loadings->n,
    od->mal.x_weights->n,
    od->mal.x_loadings->m, 1.0,
    od->mal.x_loadings->base,
    od->mal.x_loadings->max_m,
    od->mal.x_weights->base,
    od->mal.x_weights->max_m, 0.0,
    od->mal.temp->base,
    od->mal.temp->max_m);
  
  /*
  if the matrix is singular, replace weights_star with weights
  */
  #ifdef HAVE_LIBMKL
  dgetrf(&(od->mal.temp->m),
    &(od->mal.temp->m),
    od->mal.temp->base,
    &(od->mal.temp->max_m),
    od->mel.ipiv, &info);
  #elif HAVE_LIBSUNPERF
  dgetrf(od->mal.temp->m,
    od->mal.temp->m,
    od->mal.temp->base,
    od->mal.temp->max_m,
    od->mel.ipiv, &info);
  #elif HAVE_LIBLAPACK_ATLAS
  info = clapack_dgetrf(CblasColMajor,
    od->mal.temp->m,
    od->mal.temp->m,
    od->mal.temp->base,
    od->mal.temp->max_m,
    od->mel.ipiv);
  #elif HAVE_LIBLAPACKE
  info = LAPACKE_dgetrf(LAPACK_COL_MAJOR,
    od->mal.temp->m,
    od->mal.temp->m,
    od->mal.temp->base,
    od->mal.temp->max_m,
    od->mel.ipiv);
  #else
  dgetrf_(&(od->mal.temp->m),
    &(od->mal.temp->m),
    od->mal.temp->base,
    &(od->mal.temp->max_m),
    od->mel.ipiv, &info);
  #endif
  if (!info) {
    #ifdef HAVE_LIBMKL
    lwork = (pc_num + 1) * LWORK_BLOCK_SIZE * sizeof(double);
    dgetri(&(od->mal.temp->m),
      od->mal.temp->base,
      &(od->mal.temp->max_m),
      od->mel.ipiv,
      od->mel.work, &lwork, &info);
    #elif HAVE_LIBSUNPERF
    dgetri(od->mal.temp->m,
      od->mal.temp->base,
      od->mal.temp->max_m,
      od->mel.ipiv, &info);
    #elif HAVE_LIBLAPACK_ATLAS
    info = clapack_dgetri(CblasColMajor,
      od->mal.temp->m,
      od->mal.temp->base,
      od->mal.temp->max_m,
      od->mel.ipiv);
    #elif HAVE_LIBLAPACKE
    info = LAPACKE_dgetri(LAPACK_COL_MAJOR,
      od->mal.temp->m,
      od->mal.temp->base,
      od->mal.temp->max_m,
      od->mel.ipiv);
    #else
    lwork = (pc_num + 1) * LWORK_BLOCK_SIZE * sizeof(double);
    dgetri_(&(od->mal.temp->m),
      od->mal.temp->base,
      &(od->mal.temp->max_m),
      od->mel.ipiv,
      od->mel.work, &lwork, &info);
    #endif
  }
  if (!info) {
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
      od->mal.x_weights->m,
      od->mal.temp->n,
      pc_num, 1.0,
      od->mal.x_weights->base,
      od->mal.x_weights->max_m,
      od->mal.temp->base,
      od->mal.temp->max_m, 0.0,
      od->mal.x_weights_star->base,
      od->mal.x_weights_star->max_m);
  }
  else {
    memcpy(od->mal.x_weights_star->base, od->mal.x_weights->base,
      od->mal.x_weights->m * od->mal.x_weights->n * sizeof(double));
  }
  for (i = pc_num; i >= 0; --i) {
    /*
    no need to check the return value, since this is just a logical
    resizing: all these matrices have been allocated
    with maximal size with this purpose
    */
    double_mat_resize(od->mal.x_weights_star, od->mal.e_mat->n, i);
    double_mat_resize(od->mal.y_loadings, od->y_vars, i);
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,
      od->mal.x_weights_star->m,
      od->mal.y_loadings->m,
      i, 1.0,
      od->mal.x_weights_star->base,
      od->mal.x_weights_star->max_m,
      od->mal.y_loadings->base,
      od->mal.y_loadings->max_m, 0.0,
      od->mal.b_coefficients->base,
      od->mal.b_coefficients->max_m);
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
      od->mal.e_mat->m,
      od->mal.b_coefficients->n,
      od->mal.e_mat->n, 1.0,
      od->mal.e_mat->base,
      od->mal.e_mat->max_m,
      od->mal.b_coefficients->base,
      od->mal.b_coefficients->max_m, 0.0,
      od->mal.pred_f_mat->base,
      od->mal.pred_f_mat->max_m);
    for (x = 0; x < od->y_vars; ++x) {
      object_num = 0;
      j = 0;
      y = 0;
      while ((object_num < od->object_num)
        && (j < od->pel.out_structs->size)) {
        struct_num = od->al.mol_info[object_num]->struct_num;
        conf_num = 1;
        if (struct_num != od->pel.out_structs->pe[j]) {
          object_num += conf_num;
          continue;
        }
        for (n_conf = 0, sumweight = 0.0, ave_res = 0.0;
          n_conf < conf_num; ++n_conf, ++object_num) {
          if (!get_object_attr(od, object_num, ACTIVE_BIT)) {
            continue;
          }
          sumweight += od->mel.object_weight[object_num];
          sqrt_weight = sqrt(od->mel.object_weight[object_num]);
          ave_res += (M_PEEK(od->mal.f_mat, y, x)
            - M_PEEK(od->mal.pred_f_mat, y, x))
            * sqrt_weight;
          if (model_type & (CV_MODEL | SILENT_PLS)) {
            /*
            write predicted value to file
            */
            value = ((sqrt_weight > 0.0)
              ? M_PEEK(od->mal.pred_f_mat, y, x) / sqrt_weight
              + od->vel.f_mat_ave->ve[x] : INFINITY);
            actual_len = fwrite(&value, sizeof(double), 1,
              temp_pred->handle);
            if (actual_len != 1) {
              return CANNOT_WRITE_TEMP_FILE;
            }
          }
          ++y;
        }
        if (sumweight > 0.0) {
          ave_res /= sumweight;
          if (model_type & PARALLEL_CV) {
            #ifndef WIN32
            pthread_mutex_lock(od->mel.mutex);
            #else
            WaitForSingleObject(*(od->mel.mutex), INFINITE);
            #endif
          }
          M_POKE(od->mal.press, i, x,
            M_PEEK(od->mal.press, i, x) + square(ave_res));
          if (model_type & PARALLEL_CV) {
            #ifndef WIN32
            pthread_mutex_unlock(od->mel.mutex);
            #else
            ReleaseMutex(*(od->mel.mutex));
            #endif
          }
        }
        ++j;
      }
    }
    if ((model_type & UVEPLS_CV_MODEL) && (i == pc_num)) {
      /*
      store coefficients matrix
      */
      for (x = 0; x < od->mal.b_coefficients->n; ++x) {
        if (od->uvepls.save_ram) {
          if (model_type & PARALLEL_CV) {
            #ifndef WIN32
            pthread_mutex_lock(od->mel.mutex);
            #else
            WaitForSingleObject(*(od->mel.mutex), INFINITE);
            #endif
          }
          actual_len = fwrite
            (&M_PEEK(od->mal.b_coefficients, 0, x),
            sizeof(double), od->mal.b_coefficients->m,
            temp_cv_coeff->handle);
          if (actual_len != od->mal.b_coefficients->m) {
            return CANNOT_WRITE_TEMP_FILE;
          }
          if (model_type & PARALLEL_CV) {
            #ifndef WIN32
            pthread_mutex_unlock(od->mel.mutex);
            #else
            ReleaseMutex(*(od->mel.mutex));
            #endif
          }
        }
        else {
          cblas_dcopy(od->mal.b_coefficients->m,
            &M_PEEK(od->mal.b_coefficients, 0, x), 1,
            &M_PEEK(od->mal.b_coefficients_store,
            cv_run * od->mal.b_coefficients->n + x, 0),
            od->mal.b_coefficients_store->max_m);
        }
      }
    }
  }

  return 0;
}


int pred_ext_y_values(O3Data *od, int pc_num, int model_type)
{
  int i;
  int object_num;
  int struct_num;
  int conf_num;
  int n_conf;
  int actual;
  int x;
  int y;
  int result;
  double value;
  double ave_res;
  double sqrt_weight;
  double sumweight;
  

  /*
  no need to check the return value, since this is just a logical
  resizing: all these matrices have been allocated
  with maximal size with this purpose
  */
  double_mat_resize(od->mal.f_mat,
    od->ext_pred_object_num, od->y_vars);
  double_mat_resize(od->mal.e_mat,
    od->ext_pred_object_num, od->mal.e_mat->n);
  double_mat_resize(od->mal.pred_f_mat,
    od->ext_pred_object_num, od->y_vars);
  for (object_num = 0, y = 0; object_num < od->object_num; ++object_num) {
    if (get_object_attr(od, object_num, PREDICT_BIT)) {
      fill_y_vector(od, object_num, y, model_type, 0);
      fill_x_vector(od, object_num, y, model_type, 0);
      ++y;
    }
  }
  if (model_type & FULL_MODEL) {
    /*
    write number of PCs
    */
    actual = fwrite(&pc_num, sizeof(int), 1,
      od->file[TEMP_EXT_PRED]->handle);
    if (actual != 1) {
      return CANNOT_WRITE_TEMP_FILE;
    }
    /*
    write number of y_vars
    */
    actual = fwrite(&(od->y_vars), sizeof(int), 1,
      od->file[TEMP_EXT_PRED]->handle);
    if (actual != 1) {
      return CANNOT_WRITE_TEMP_FILE;
    }
    /*
    write number of objects in this block
    */
    actual = fwrite(&(od->ext_pred_object_num), sizeof(int), 1,
      od->file[TEMP_EXT_PRED]->handle);
    if (actual != 1) {
      return CANNOT_WRITE_TEMP_FILE;
    }
    /*
    write object numbers in this block
    */
    for (object_num = 0; object_num < od->object_num; ++object_num) {
      if (get_object_attr(od, object_num, PREDICT_BIT)) {
        actual = fwrite(&object_num, sizeof(int), 1,
          od->file[TEMP_EXT_PRED]->handle);
        if (actual != 1) {
          return CANNOT_WRITE_TEMP_FILE;
        }
      }
    }
  }
  /*
  predicted y values are stored in this way:
  pc_num
  number of variables (equal to x_max_y, that is y_vars)
  value_num(here, y_max), var_num(here, 0), var-0_value-0_pc-pc_num, var-0_value-1_pc-pc_num, var-0_value-2_pc-pc_num, ..., var-0_value-y_max_pc-pc_num
  value_num(here, y_max), var_num(here, 1), var-1_value-0_pc-pc_num, var-1_value-1_pc-pc_num, var-1_value-2_pc-pc_num, ..., var-1_value-y_max_pc-pc_num
  value_num(here, y_max), var_num(here, 2), var-2_value-0_pc-pc_num, var-2_value-1_pc-pc_num, var-2_value-2_pc-pc_num, ..., var-2_value-y_max_pc-pc_num
  ...
  value_num(here, y_max), var_num(here, x_max_y), var-x_max_y_value-0_pc-pc_num, var-x_max_y_value-1_pc-pc_num, var-x_max_y_value-2_pc-pc_num, ..., var-x_max_y_value-y_max_pc-pc_num

  value_num(here, y_max), var_num(here, 0), var-0_value-0_pc-(pc_num-1), var-0_value-1_pc-(pc_num-1), var-0_value-2_pc-(pc_num-1), ..., var-0_value-y_max_pc-(pc_num-1)
  value_num(here, y_max), var_num(here, 1), var-1_value-0_pc-(pc_num-1), var-1_value-1_pc-(pc_num-1), var-1_value-2_pc-(pc_num-1), ..., var-1_value-y_max_pc-(pc_num-1)
  value_num(here, y_max), var_num(here, 2), var-2_value-0_pc-(pc_num-1), var-2_value-1_pc-(pc_num-1), var-2_value-2_pc-(pc_num-1), ..., var-2_value-y_max_pc-(pc_num-1)
  ...
  value_num(here, y_max), var_num(here, x_max_y), var-x_max_y_value-0_pc-(pc_num-1), var-x_max_y_value-1_pc-(pc_num-1), var-x_max_y_value-2_pc-(pc_num-1), ..., var-x_max_y_value-y_max_pc-(pc_num-1)

  value_num(here, y_max), var_num(here, 0), var-0_value-0_pc-(pc_num-2), var-0_value-1_pc-(pc_num-2), var-0_value-2_pc-(pc_num-2), ..., var-0_value-y_max_pc-(pc_num-2)
  value_num(here, y_max), var_num(here, 1), var-1_value-0_pc-(pc_num-2), var-1_value-1_pc-(pc_num-2), var-1_value-2_pc-(pc_num-2), ..., var-1_value-y_max_pc-(pc_num-2)
  value_num(here, y_max), var_num(here, 2), var-2_value-0_pc-(pc_num-2), var-2_value-1_pc-(pc_num-2), var-2_value-2_pc-(pc_num-2), ..., var-2_value-y_max_pc-(pc_num-2)
  ...
  value_num(here, y_max), var_num(here, x_max_y), var-x_max_y_value-0_pc-(pc_num-2), var-x_max_y_value-1_pc-(pc_num-2), var-x_max_y_value-2_pc-(pc_num-2), ..., var-x_max_y_value-y_max_pc-(pc_num-2)

  ...
  ...
  ...

  value_num(here, y_max), var_num(here, 0), var-0_value-0_pc-0, var-0_value-1_pc-0, var-0_value-2_pc-0, ..., var-0_value-y_max_pc-0
  value_num(here, y_max), var_num(here, 1), var-1_value-0_pc-0, var-1_value-1_pc-0, var-1_value-2_pc-0, ..., var-1_value-y_max_pc-0
  value_num(here, y_max), var_num(here, 2), var-2_value-0_pc-0, var-2_value-1_pc-0, var-2_value-2_pc-0, ..., var-2_value-y_max_pc-0
  ...
  value_num(here, y_max), var_num(here, x_max_y), var-x_max_y_value-0_pc-0, var-x_max_y_value-1_pc-0, var-x_max_y_value-2_pc-0, ..., var-x_max_y_value-y_max_pc-0
  */
  memset(od->mal.press->base, 0,
    od->mal.press->m * od->mal.press->n * sizeof(double));
  for (i = pc_num; i >= 0; --i) {
    if (i) {
      result = reload_coefficients(od, i);
      if (result) {
        return CANNOT_READ_TEMP_FILE;
      }
    }
    else {
      memset(od->mal.b_coefficients->base, 0,
        od->mal.b_coefficients->m
        * od->mal.b_coefficients->n * sizeof(double));
    }
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
      od->mal.e_mat->m,
      od->mal.b_coefficients->n,
      od->mal.e_mat->n, 1.0,
      od->mal.e_mat->base,
      od->mal.e_mat->max_m,
      od->mal.b_coefficients->base,
      od->mal.b_coefficients->max_m, 0.0,
      od->mal.pred_f_mat->base,
      od->mal.pred_f_mat->max_m);
    for (x = 0; x < od->y_vars; ++x) {
      object_num = 0;
      y = 0;
      while (object_num < od->object_num) {
        struct_num = od->al.mol_info[object_num]->struct_num;
        conf_num = 1;
        for (n_conf = 0, sumweight = 0.0, ave_res = 0.0;
          n_conf < conf_num; ++n_conf, ++object_num) {
          if (!get_object_attr(od, object_num, PREDICT_BIT)) {
            continue;
          }
          sumweight += od->mel.object_weight[object_num];
          sqrt_weight = sqrt(od->mel.object_weight[object_num]);
          ave_res += (M_PEEK(od->mal.f_mat, y, x)
            - M_PEEK(od->mal.pred_f_mat, y, x))
            * sqrt_weight;
          value = ((sqrt_weight > 0.0)
            ? M_PEEK(od->mal.pred_f_mat, y, x) / sqrt_weight
            + od->vel.f_mat_full_ave->ve[x] : INFINITY);
          if (model_type & FULL_MODEL) {
            /*
            write predicted value to file
            */
            actual = fwrite(&value, sizeof(double), 1,
              od->file[TEMP_EXT_PRED]->handle);
            if (actual != 1) {
              return CANNOT_WRITE_TEMP_FILE;
            }
          }
          ++y;
        }
        if (sumweight > 0.0) {
          ave_res /= sumweight;
          M_POKE(od->mal.press, i, x,
            M_PEEK(od->mal.press, i, x) + square(ave_res));
        }
      }
    }
  }

  return 0;
}
