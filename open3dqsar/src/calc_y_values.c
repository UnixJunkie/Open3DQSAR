/*

calc_y_values.c

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


int calc_y_values(O3Data *od)
{
  int i;
  int x;
  int y;
  int num_blocks;
  int x_max_x;
  int x_max_y;
  int y_max;
  int object_num;
  int struct_num;
  int conf_num;
  int n_conf;
  int pc_num;
  int actual_len;
  int result;
  int info;
  #if (!defined HAVE_LIBLAPACK_ATLAS) && (!defined HAVE_LIBSUNPERF)
  int lwork;
  #endif
  double sqrt_weight;
  
  
  /*
  calculated y values are stored in this way:
  pc_num
  total number of variables (equal to x_max_y, that is y_vars)

  *** BLOCK 1 ***
  number of objects in this block (cmp_blk)
  list of object numbers in this block (length of the list: cmp_blk)
  var-0_value-0_pc-pc_num, var-0_value-1_pc-pc_num, var-0_value-2_pc-pc_num, ..., var-0_value-cmp_blk-pc_num
  var-1_value-0_pc-pc_num, var-1_value-1_pc-pc_num, var-1_value-2_pc-pc_num, ..., var-1_value-cmp_blk-pc_num
  var-2_value-0_pc-pc_num, var-2_value-1_pc-pc_num, var-2_value-2_pc-pc_num, ..., var-2_value-cmp_blk-pc_num
  ...
  var-x_max_y_value-0_pc-pc_num, var-x_max_y_value-1_pc-pc_num, var-x_max_y_value-2_pc-pc_num, ..., var-x_max_y_value-cmp_blk-pc_num

  var-0_value-0_pc-(pc_num-1), var-0_value-1_pc-(pc_num-1), var-0_value-2_pc-(pc_num-1), ..., var-0_value-cmp_blk-(pc_num-1)
  var-1_value-0_pc-(pc_num-1), var-1_value-1_pc-(pc_num-1), var-1_value-2_pc-(pc_num-1), ..., var-1_value-cmp_blk-(pc_num-1)
  var-2_value-0_pc-(pc_num-1), var-2_value-1_pc-(pc_num-1), var-2_value-2_pc-(pc_num-1), ..., var-2_value-cmp_blk-(pc_num-1)
  ...
  var-x_max_y_value-0_pc-(pc_num-1), var-x_max_y_value-1_pc-(pc_num-1), var-x_max_y_value-2_pc-(pc_num-1), ..., var-x_max_y_value-cmp_blk-(pc_num-1)

  var-0_value-0_pc-(pc_num-2), var-0_value-1_pc-(pc_num-2), var-0_value-2_pc-(pc_num-2), ..., var-0_value-cmp_blk-(pc_num-2)
  var-1_value-0_pc-(pc_num-2), var-1_value-1_pc-(pc_num-2), var-1_value-2_pc-(pc_num-2), ..., var-1_value-cmp_blk-(pc_num-2)
  var-2_value-0_pc-(pc_num-2), var-2_value-1_pc-(pc_num-2), var-2_value-2_pc-(pc_num-2), ..., var-2_value-cmp_blk-(pc_num-2)
  ...
  var-x_max_y_value-0_pc-(pc_num-2), var-x_max_y_value-1_pc-(pc_num-2), var-x_max_y_value-2_pc-(pc_num-2), ..., var-x_max_y_value-cmp_blk-(pc_num-2)

  ...
  ...
  ...

  var-0_value-0_pc-(pc-0), var-0_value-1_pc-(pc-0), var-0_value-2_pc-(pc-0), ..., var-0_value-cmp_blk-(pc-0)
  var-1_value-0_pc-(pc-0), var-1_value-1_pc-(pc-0), var-1_value-2_pc-(pc-0), ..., var-1_value-cmp_blk-(pc-0)
  var-2_value-0_pc-(pc-0), var-2_value-1_pc-(pc-0), var-2_value-2_pc-(pc-0), ..., var-2_value-cmp_blk-(pc-0)
  ...
  var-x_max_y_value-0_pc-(pc-0), var-x_max_y_value-1_pc-(pc-0), var-x_max_y_value-2_pc-(pc-0), ..., var-x_max_y_value-cmp_blk-(pc-0)

  *** BLOCK 2 ***

  ...
  ...
  ...

  */    
  x_max_x = od->overall_active_x_vars;
  x_max_y = od->y_vars;
  y_max = od->active_object_num;
  pc_num = od->pc_num;
  /*
  write number of PCs
  */
  actual_len = fwrite(&(od->pc_num), sizeof(int), 1,
    od->file[TEMP_CALC]->handle);
  if (actual_len != 1) {
    return CANNOT_WRITE_TEMP_FILE;
  }
  /*
  write number of y_vars
  */
  actual_len = fwrite(&x_max_y, sizeof(int), 1,
    od->file[TEMP_CALC]->handle);
  if (actual_len != 1) {
    return CANNOT_WRITE_TEMP_FILE;
  }
  /*
  write number of objects in this block
  */
  actual_len = fwrite(&y_max, sizeof(int), 1,
    od->file[TEMP_CALC]->handle);
  if (actual_len != 1) {
    return CANNOT_WRITE_TEMP_FILE;
  }
  /*
  write object numbers in this block
  */
  for (i = 0; i < od->object_num; ++i) {
    if (get_object_attr(od, i, ACTIVE_BIT)) {
      actual_len = fwrite(&i, sizeof(int), 1,
        od->file[TEMP_CALC]->handle);
      if (actual_len != 1) {
        return CANNOT_WRITE_TEMP_FILE;
      }
    }
  }
  
  /*
  save coefficient matrix on a temporary file
  the first integer is the number of blocks
  */
  num_blocks = pc_num + 1;
  actual_len = fwrite(&num_blocks, sizeof(int), 1,
    od->file[TEMP_PLS_COEFF]->handle);
  if (actual_len != 1) {
    return CANNOT_WRITE_TEMP_FILE;
  }
  double_mat_resize(od->mal.x_weights, x_max_x, pc_num);
  double_mat_resize(od->mal.x_weights_star, x_max_x, pc_num);
  double_mat_resize(od->mal.x_loadings, x_max_x, pc_num);
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
    memcpy(od->mal.x_weights_star->base,
      od->mal.x_weights->base,
      od->mal.x_weights->m
      * od->mal.x_weights->n
      * sizeof(double));
  }
  for (i = pc_num; i >= 0; --i) {
    double_mat_resize(od->mal.x_weights_star, x_max_x, i);
    double_mat_resize(od->mal.y_loadings, x_max_y, i);
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
    
    double_mat_resize(od->mal.x_scores, y_max, i);
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,
      od->mal.x_scores->m,
      od->mal.y_loadings->m,
      i, 1.0,
      od->mal.x_scores->base,
      od->mal.x_scores->max_m,
      od->mal.y_loadings->base,
      od->mal.y_loadings->max_m, 0.0,
      od->mal.pred_f_mat->base,
      od->mal.pred_f_mat->max_m);
    
    /*
    save coefficient matrix on a temporary file
    from the highest to the lowest PC (num_blocks overall)
    */
    for (x = 0; x < x_max_y; ++x) {
      for (y = 0; y < x_max_x; ++y) {
        actual_len = fwrite(&M_PEEK(od->mal.b_coefficients, y, x),
          sizeof(double), 1, od->file[TEMP_PLS_COEFF]->handle);
        if (actual_len != 1) {
          return CANNOT_WRITE_TEMP_FILE;
        }
      }
    }
    for (x = 0; x < x_max_y; ++x) {
      object_num = 0;
      y = 0;
      while (object_num < od->object_num) {
        struct_num = od->al.mol_info[object_num]->struct_num;
        conf_num = 1;
        for (n_conf = 0; n_conf < conf_num; ++n_conf, ++object_num) {
          if (get_object_attr(od, object_num, ACTIVE_BIT)) {
            sqrt_weight = sqrt(od->mel.object_weight[object_num]);
            M_POKE(od->mal.pred_f_mat, y, x,
              ((sqrt_weight > 0.0)
              ? M_PEEK(od->mal.pred_f_mat, y, x) / sqrt_weight
              + od->vel.f_mat_full_ave->ve[x] : INFINITY));
            /*
            write calculated pred_f_mat to file
            */
            actual_len = fwrite(&M_PEEK(od->mal.pred_f_mat, y, x),
              sizeof(double), 1, od->file[TEMP_CALC]->handle);
            if (actual_len != 1) {
              return CANNOT_WRITE_TEMP_FILE;
            }
            ++y;
          }
        }
      }
    }
  }
  rewind(od->file[TEMP_CALC]->handle);
  if (od->file[TEMP_PLS_COEFF]->handle) {
    fclose(od->file[TEMP_PLS_COEFF]->handle);
    od->file[TEMP_PLS_COEFF]->handle = NULL;
  }
  result = reload_weights_loadings(od);

  return result;
}
