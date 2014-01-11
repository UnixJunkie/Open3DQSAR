/*

fill_vector.c

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


void fill_x_vector(O3Data *od, int object_num, int row, int model_type, int cv_run)
{
  char sel_one_zero;
  int n_active;
  int n_predict;
  int x;
  int i;
  int real_x;
  int size_coeff;
  double value;
  
  
  /*
  copy x_vars values for a certain object,
  mean-center and store them into the E matrix
  */
  real_x = 0;
  sel_one_zero = '1';
  size_coeff = 1;
  if ((model_type & UVEPLS_CV_MODEL)
    && (!(od->uvepls.ive))) {
    size_coeff = 2;
  }
  for (i = 0, n_active = 0, n_predict = od->active_object_num; i < object_num; ++i) {
    if (get_object_attr(od, i, ACTIVE_BIT)) {
      ++n_active;
    }
    if (get_object_attr(od, i, PREDICT_BIT)) {
      ++n_predict;
    }
  }
  for (x = 0; x < (od->mal.large_e_mat->n / size_coeff); ++x) {
    if (model_type & FFDSEL_CV_MODEL) {
      sel_one_zero = od->mel.ffdsel_included[x];
    }
    if (model_type & UVEPLS_CV_MODEL) {
      sel_one_zero = od->mel.uvepls_included[x];
    }
    if (sel_one_zero == '1') {
      value = M_PEEK(od->mal.large_e_mat, (get_object_attr
        (od, object_num, ACTIVE_BIT) ? n_active : n_predict), x);
      if (model_type & (FFDSEL_CV_MODEL | UVEPLS_CV_MODEL | SCRAMBLE_CV_MODEL)) {
        M_POKE(od->mal.e_mat, row, real_x,
          (MISSING(value) ? 0.0 : (value
          - M_PEEK(od->mal.large_e_mat_ave, cv_run, x))
          * sqrt(od->mel.object_weight[object_num])));
      }
      else {
        M_POKE(od->mal.e_mat, row, real_x,
          (MISSING(value) ? 0.0 : (value
          - od->vel.e_mat_ave->ve[real_x])
          * sqrt(od->mel.object_weight[object_num])));
      }
      if ((model_type & UVEPLS_CV_MODEL) && (!(od->uvepls.ive))) {
        value = M_PEEK(od->mal.large_e_mat, (get_object_attr
        (od, object_num, ACTIVE_BIT) ? n_active : n_predict),
        x + od->mal.large_e_mat->n / 2);
        M_POKE(od->mal.e_mat, row,
          real_x + od->mal.e_mat->n / 2,
          (MISSING(value) ? 0.0 : (value
          - M_PEEK(od->mal.large_e_mat_ave,
          cv_run, x + od->mal.large_e_mat->n / 2))
          * sqrt(od->mel.object_weight[object_num])));
      }
      ++real_x;
    }
  }
}


void fill_y_vector(O3Data *od, int object_num, int row, int model_type, int cv_run)
{
  int n_active;
  int n_predict;
  int x;
  int i;
  double value;
  
  
  /*
  copy y_vars values for a certain object,
  mean-center and store them into the F matrix
  */
  for (i = 0, n_active = 0, n_predict = od->active_object_num; i < object_num; ++i) {
    if (get_object_attr(od, i, ACTIVE_BIT)) {
      ++n_active;
    }
    if (get_object_attr(od, i, PREDICT_BIT)) {
      ++n_predict;
    }
  }
  for (x = 0; x < od->mal.f_mat->n; ++x) {
    value = M_PEEK(od->mal.large_f_mat, (get_object_attr
      (od, object_num, ACTIVE_BIT) ? n_active : n_predict), x);
    if (model_type & (FFDSEL_CV_MODEL | UVEPLS_CV_MODEL)) {
      M_POKE(od->mal.f_mat, row, x,
        (value - M_PEEK(od->mal.large_f_mat_ave, cv_run, x))
        * sqrt(od->mel.object_weight[object_num]));
    }
    else {
      M_POKE(od->mal.f_mat, row, x,
        (value - od->vel.f_mat_ave->ve[x])
        * sqrt(od->mel.object_weight[object_num]));
    }
  }
}
