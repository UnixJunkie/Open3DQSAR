/*

trim_mean_center_matrix.c

is part of

Open3DQSAR
----------

An open-source software aimed at high-throughput
chemometric analysis of molecular interaction fields

Copyright (C) 2009-2018 Paolo Tosco, Thomas Balle

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


void trim_mean_center_matrix(O3Data *od, DoubleMat *large_mat, DoubleMat **mat,
  DoubleVec **mat_ave, int model_type, int active_object_num)
{
  int y;
  int x;
  int real_y;
  int n;
  int object_num;
  int struct_num;
  int conf_num;
  int n_conf;
  double value;
  double sumweight = 0.0;
  
  
  /*
  allocate a active_object_num * large_mat->n matrix
  */
  *mat = double_mat_resize(*mat, active_object_num, (*mat)->n);
  /*
  allocate one vector for averages
  */
  *mat_ave = double_vec_resize(*mat_ave, (*mat)->n);
  memset((*mat_ave)->ve, 0, (*mat_ave)->size * sizeof(double));
  /*
  copy values from large-E matrix,
  mean-center them and store them into the E matrix
  */
  for (x = 0; x < large_mat->n; ++x) {
    n = 0;
    object_num = 0;
    y = 0;
    sumweight = 0.0;
    while (object_num < od->object_num) {
      struct_num = od->al.mol_info[object_num]->struct_num;
      conf_num = 1;
      /*
      if it is a CV model some compounds shall be left out
      */
      if (model_type & (CV_MODEL | SCRAMBLE_CV_MODEL)) {
        if (n < od->pel.out_structs->size) {
          if (struct_num == od->pel.out_structs->pe[n]) {
            ++n;
            object_num += conf_num;
            y += conf_num;
            continue;
          }
        }
      }
      for (n_conf = 0; n_conf < conf_num; ++n_conf, ++object_num) {
        if (get_object_attr(od, object_num, ACTIVE_BIT)) {
          value = M_PEEK(large_mat, y, x);
          if (!MISSING(value)) {
            (*mat_ave)->ve[x] += (value * od->mel.object_weight[object_num]);
            sumweight += od->mel.object_weight[object_num];
          }
          ++y;
        }
      }
    }
    if (sumweight > 0.0) {
      (*mat_ave)->ve[x] /= sumweight;
    }
    real_y = 0;
    n = 0;
    object_num = 0;
    y = 0;
    while (object_num < od->object_num) {
      struct_num = od->al.mol_info[object_num]->struct_num;
      conf_num = 1;
      /*
      if it is a CV model some compounds shall be left out
      */
      if (model_type & (CV_MODEL | SCRAMBLE_CV_MODEL)) {
        if (n < od->pel.out_structs->size) {
          if (struct_num == od->pel.out_structs->pe[n]) {
            ++n;
            object_num += conf_num;
            y += conf_num;
            continue;
          }
        }
      }
      for (n_conf = 0; n_conf < conf_num; ++n_conf, ++object_num) {
        if (get_object_attr(od, object_num, ACTIVE_BIT)) {
          value = M_PEEK(large_mat, y, x);
          M_POKE(*mat, real_y, x, (MISSING(value)
            ? 0.0 : (value - (*mat_ave)->ve[x]))
            * sqrt(od->mel.object_weight[object_num]));
          ++real_y;
          ++y;
        }
      }
    }
  }
}


void trim_mean_center_x_matrix_hp(O3Data *od, int model_type, int active_object_num, int run)
{
  char sel_one_zero;
  int y;
  int x;
  int real_y;
  int real_x;
  int size_coeff;
  int n;
  int object_num;
  int struct_num;
  int conf_num;
  int n_conf;
  double value;
  
  
  size_coeff = 1;
  if ((model_type & UVEPLS_CV_MODEL)
    && (!(od->uvepls.ive))) {
    size_coeff = 2;
  }
  if (model_type & FFDSEL_CV_MODEL) {
    od->mal.e_mat->n = od->mal.large_e_mat->n;
  }
  /*
  allocate a active_object_num * od->mal.large_e_mat->n E matrix
  */
  od->mal.e_mat = double_mat_resize(od->mal.e_mat,
    active_object_num, od->mal.e_mat->n);
  /*
  copy values from large-E matrix,
  mean-center them and store them into the E matrix
  */
  real_x = 0;
  real_y = 0;
  sel_one_zero = '1';
  for (x = 0; x < (od->mal.large_e_mat->n / size_coeff); ++x) {
    if (model_type & (FFDSEL_FULL_MODEL | FFDSEL_CV_MODEL)) {
      sel_one_zero = od->mel.ffdsel_included[x];
    }
    if (model_type & (UVEPLS_FULL_MODEL | UVEPLS_CV_MODEL)) {
      sel_one_zero = od->mel.uvepls_included[x];
    }
    if (sel_one_zero == '1') {
      real_y = 0;
      n = 0;
      object_num = 0;
      y = 0;
      while (object_num < od->object_num) {
        struct_num = od->al.mol_info[object_num]->struct_num;
        conf_num = 1;
        if (n < od->pel.out_structs->size) {
          if (struct_num == od->pel.out_structs->pe[n]) {
            ++n;
            object_num += conf_num;
            y += conf_num;
            continue;
          }
        }
        for (n_conf = 0; n_conf < conf_num; ++n_conf, ++object_num) {
          if (get_object_attr(od, object_num, ACTIVE_BIT)) {
            value = M_PEEK(od->mal.large_e_mat, y, x);
            M_POKE(od->mal.e_mat, real_y, real_x, (MISSING(value)
              ? 0.0 : (value - M_PEEK(od->mal.large_e_mat_ave, run, x))
              * sqrt(od->mel.object_weight[object_num])));
            if ((model_type & UVEPLS_CV_MODEL) && (!(od->uvepls.ive))) {
              value = M_PEEK(od->mal.large_e_mat, y,
                x + od->mal.large_e_mat->n / 2);
              M_POKE(od->mal.e_mat, real_y,
                real_x + od->mal.e_mat->n / 2, (MISSING(value)
                ? 0.0 : (value - M_PEEK(od->mal.large_e_mat_ave, run,
                x + od->mal.large_e_mat->n / 2))
                * sqrt(od->mel.object_weight[object_num])));
            }
            ++real_y;
            ++y;
          }
        }
      }
      ++real_x;
    }
  }
  /*
  the following loop removes the empty gap between
  the block of real variables and the block of dummy
  variables
  */
  if ((model_type & UVEPLS_CV_MODEL) && (!(od->uvepls.ive))) {
    for (y = 0; y < real_y; ++y) {
      for (x = 0; x < real_x; ++x) {
        M_POKE(od->mal.e_mat, y, x + real_x,
          M_PEEK(od->mal.e_mat, y,
          x + od->mal.e_mat->n / 2));
      }
    }
  }
  /*
  trim the E matrix size
  */
  od->mal.e_mat = double_mat_resize(od->mal.e_mat,
    active_object_num, real_x * size_coeff);
}


void trim_mean_center_y_matrix_hp(O3Data *od, int active_object_num, int run)
{
  int y;
  int x;
  int real_y;
  int n;
  int object_num;
  int struct_num;
  int conf_num;
  int n_conf;
  double value;
  
  
  /*
  allocate a active_object_num * od->mal.large_f_mat->n F matrix
  */
  od->mal.f_mat = double_mat_resize(od->mal.f_mat,
    active_object_num, od->mal.large_f_mat->n);
  /*
  copy values from large-F matrix,
  mean-center them and store them into the F matrix
  */
  for (x = 0; x < od->mal.large_f_mat->n; ++x) {
    real_y = 0;
    n = 0;
    object_num = 0;
    y = 0;
    while (object_num < od->object_num) {
      struct_num = od->al.mol_info[object_num]->struct_num;
      conf_num = 1;
      if (n < od->pel.out_structs->size) {
        if (struct_num == od->pel.out_structs->pe[n]) {
          ++n;
          object_num += conf_num;
          y += conf_num;
          continue;
        }
      }
      for (n_conf = 0; n_conf < conf_num; ++n_conf, ++object_num) {
        if (get_object_attr(od, object_num, ACTIVE_BIT)) {
          value = M_PEEK(od->mal.large_f_mat, y, x);
          M_POKE(od->mal.f_mat, real_y, x,
            (value - M_PEEK(od->mal.large_f_mat_ave, run, x))
            * sqrt(od->mel.object_weight[object_num]));
          ++real_y;
          ++y;
        }
      }
    }
  }
}


void trim_mean_center_x_matrix_pca(O3Data *od)
{
  int y;
  int x;
  int object_num;
  double value;
  
  
  /*
  allocate a od->active_object_num * od->mal.large_e_mat->n E matrix
  */
  od->mal.e_mat = double_mat_resize(od->mal.e_mat,
    od->active_object_num, od->mal.e_mat->n);
  /*
  allocate one vector for averages
  */
  od->vel.e_mat_ave = double_vec_resize
    (od->vel.e_mat_ave, od->mal.e_mat->n);
  memset(od->vel.e_mat_ave->ve, 0,
    od->vel.e_mat_ave->size * sizeof(double));
  /*
  copy values from large-E matrix,
  mean-center them and store them into the E matrix
  */
  for (x = 0; x < od->mal.large_e_mat->n; ++x) {
    for (object_num = 0, y = 0; object_num < od->object_num; ++object_num) {
      if (get_object_attr(od, object_num, ACTIVE_BIT)) {
        value = M_PEEK(od->mal.large_e_mat, y, x);
        if (!MISSING(value)) {
          od->vel.e_mat_ave->ve[x] += value;
          ++y;
        }
      }
    }
    if (y) {
      od->vel.e_mat_ave->ve[x] /= (double)y;
    }
    for (object_num = 0, y = 0; object_num < od->object_num; ++object_num) {
      if (get_object_attr(od, object_num, ACTIVE_BIT)) {
        value = M_PEEK(od->mal.large_e_mat, y, x);
        M_POKE(od->mal.e_mat, y, x, (MISSING(value)
          ? 0.0 : (value - od->vel.e_mat_ave->ve[x])));
        ++y;
      }
    }
  }
}
