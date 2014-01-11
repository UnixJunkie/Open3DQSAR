/*

pls.c

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


void pls(O3Data *od, int suggested_pc_num, int model_type)
{
  int i;
  int x;
  int y;
  int pc_num;
  int conv;
  double norm_c;
  double norm_d;
  double ss_u;
  double ss_v_diff;
  double s2_x_tot;
  double s2_x_tot_start = 0.0;
  double s2_y_tot;
  double s2_y_tot_start = 0.0;
  double sum_object_weight = 0.0;
  double cum_explained_s2_x = 0.0;
  double cum_explained_s2_y = 0.0;
    

  pc_num = suggested_pc_num;
  /*
  this check is useful for FFD variable selection:
  if the user wishes to extract n PCs but in this particular
  combination there are not enough active variables, the number
  of PCs must be reduced
  */
  if (suggested_pc_num > od->mal.e_mat->n) {
    pc_num = od->mal.e_mat->n;
  }
  /*
  no need to check the return value, since this is just a logical
  resizing: all these matrices and vectors have been allocated
  with maximal size with this purpose
  */
  double_mat_resize(od->mal.x_scores,
    od->mal.e_mat->m, pc_num + 1);
  double_mat_resize(od->mal.x_weights,
    od->mal.e_mat->n, pc_num + 1);
  double_mat_resize(od->mal.x_loadings,
    od->mal.e_mat->n, pc_num + 1);
  double_mat_resize(od->mal.pred_f_mat,
    od->mal.e_mat->m, od->mal.f_mat->n);
  double_vec_resize(od->vel.ave_sdep, pc_num + 1);
  double_vec_resize(od->vel.v, od->mal.e_mat->m);
  double_vec_resize(od->vel.v_new, od->mal.e_mat->m);
  double_mat_resize(od->mal.y_loadings,
    od->mal.f_mat->n, pc_num + 1);
  double_mat_resize(od->mal.y_scores,
    od->mal.f_mat->m, pc_num + 1);
  double_vec_resize(od->vel.ro, pc_num + 1);
  od->pc_num = pc_num;

  /*
  Copy first dependent variable into vector v
  */
  cblas_dcopy(od->mal.f_mat->m,
    od->mal.f_mat->base, 1,
    od->vel.v->ve, 1);
  /*
  calculate total variances of the E, F matrices
  */
  if (model_type & FULL_MODEL) {
    for (i = 0; i < od->object_num; ++i) {
      if (get_object_attr(od, i, ACTIVE_BIT)) {
        sum_object_weight += od->mel.object_weight[i];
      }
    }
    for (x = 0, s2_x_tot_start = 0.0; x < od->mal.e_mat->n; ++x) {
      s2_x_tot_start += cblas_ddot(od->mal.e_mat->m,
        &M_PEEK(od->mal.e_mat, 0, x), 1, &M_PEEK(od->mal.e_mat, 0, x), 1);
    }
    for (x = 0, s2_y_tot_start = 0.0; x < od->mal.f_mat->n; ++x) {
      s2_y_tot_start += cblas_ddot(od->mal.f_mat->m,
        &M_PEEK(od->mal.f_mat, 0, x), 1, &M_PEEK(od->mal.f_mat, 0, x), 1);
    }
    tee_printf(od, "  %12s%12s%12s%12s\n", "Exp.", "Cum. exp.",
      "Exp.", "Cum. exp.");
    tee_printf(od, "PC%12s%12s%12s%12s%12s%12s\n", "var. X %", "var. X %",
      "var. Y %", "var. Y %", "SDEC", "r2");
    tee_printf(od, "--------------------------------------------------------------------------\n");
  }
        for (i = 0; i <= pc_num; ++i) {
    if (model_type & FULL_MODEL) {
      for (x = 0, s2_x_tot = 0.0; x < od->mal.e_mat->n; ++x) {
        s2_x_tot += cblas_ddot(od->mal.e_mat->m,
          &M_PEEK(od->mal.e_mat, 0, x), 1, &M_PEEK(od->mal.e_mat, 0, x), 1);
      }
      for (x = 0, s2_y_tot = 0.0; x < od->mal.f_mat->n; ++x) {
        s2_y_tot += cblas_ddot(od->mal.f_mat->m,
          &M_PEEK(od->mal.f_mat, 0, x), 1, &M_PEEK(od->mal.f_mat, 0, x), 1);
      }
      od->vel.explained_s2_x->ve[i] =
        1.0 - cum_explained_s2_x - s2_x_tot / s2_x_tot_start;
      cum_explained_s2_x += od->vel.explained_s2_x->ve[i];
      od->vel.explained_s2_y->ve[i] =
        1.0 - cum_explained_s2_y - s2_y_tot / s2_y_tot_start;
      cum_explained_s2_y += od->vel.explained_s2_y->ve[i];
      od->vel.ave_sdec->ve[i] =
        sqrt(s2_y_tot / sum_object_weight);
      od->vel.r2->ve[i] = cum_explained_s2_y;
      tee_printf(od, "%2d%12.4lf%12.4lf%12.4lf%12.4lf%12.4lf%12.4lf\n",
        i, od->vel.explained_s2_x->ve[i] * 100,
        cum_explained_s2_x * 100,
        od->vel.explained_s2_y->ve[i] * 100,
        cum_explained_s2_y * 100,
        od->vel.ave_sdec->ve[i],
        cum_explained_s2_y);
    }
    conv = 0;
    while (!conv) {
      /*
      c = E[i]'v
      */
      cblas_dgemv(CblasColMajor, CblasTrans,
        od->mal.e_mat->m,
        od->mal.e_mat->n, 1.0,
        od->mal.e_mat->base,
        od->mal.e_mat->max_m,
        od->vel.v->ve, 1, 0.0,
        &M_PEEK(od->mal.x_weights, 0, i), 1);
      
      /*
      c = c * (c'c)^(-0.5)
      */
      norm_c = cblas_dnrm2(od->mal.x_weights->m,
        &M_PEEK(od->mal.x_weights, 0, i), 1);
      cblas_dscal(od->mal.x_weights->m, 1.0 / norm_c,
        &M_PEEK(od->mal.x_weights, 0, i), 1);
      /*
      u = E[i]c
      */
      cblas_dgemv(CblasColMajor, CblasNoTrans,
        od->mal.e_mat->m,
        od->mal.e_mat->n, 1.0,
        od->mal.e_mat->base,
        od->mal.e_mat->max_m,
        &M_PEEK(od->mal.x_weights, 0, i), 1, 0.0,
        &M_PEEK(od->mal.x_scores, 0, i), 1);
      /*
      d = F[i]'u
      */
      cblas_dgemv(CblasColMajor, CblasTrans,
        od->mal.f_mat->m,
        od->mal.f_mat->n, 1.0,
        od->mal.f_mat->base,
        od->mal.f_mat->max_m,
        &M_PEEK(od->mal.x_scores, 0, i), 1, 0.0,
        &M_PEEK(od->mal.y_loadings, 0, i), 1);
      /*
      d = d * (d'd)^(-0.5)
      */
      norm_d = cblas_dnrm2(od->mal.y_loadings->m,
        &M_PEEK(od->mal.y_loadings, 0, i), 1);
      cblas_dscal(od->mal.y_loadings->m, 1.0 / norm_d,
        &M_PEEK(od->mal.y_loadings, 0, i), 1);
      /*
      v_new = F[i]d
      */
      cblas_dgemv(CblasColMajor, CblasNoTrans,
        od->mal.f_mat->m,
        od->mal.f_mat->n, 1.0,
        od->mal.f_mat->base,
        od->mal.f_mat->max_m,
        &M_PEEK(od->mal.y_loadings, 0, i),
        1, 0.0, od->vel.v_new->ve, 1);
      conv = 1;
      if (od->mal.f_mat->n > 1) {
        /*
        if there are multiple y vars, check convergence
        */
        cblas_daxpy(od->vel.v_new->size, -1.0,
          od->vel.v_new->ve, 1,
          od->vel.v->ve, 1);
        ss_v_diff = cblas_ddot(od->mal.y_scores->m,
          od->vel.v->ve, 1,
          od->vel.v->ve, 1);
        if (fabs(ss_v_diff) >= PLS_CONV_THRESHOLD) {
          conv = 0;
        }
      }
      cblas_dcopy(od->vel.v_new->size,
        od->vel.v_new->ve, 1,
        od->vel.v->ve, 1);
    }
    /*
    ro = u'v / u'u
    */
    ss_u = cblas_ddot(od->mal.x_scores->m,
      &M_PEEK(od->mal.x_scores, 0, i), 1,
      &M_PEEK(od->mal.x_scores, 0, i), 1);
    od->vel.ro->ve[i] =
      cblas_ddot(od->mal.x_scores->m,
      &M_PEEK(od->mal.x_scores, 0, i), 1,
      od->vel.v->ve, 1) / ss_u;
    
    if (model_type & FULL_MODEL) {
      cblas_dcopy(od->mal.f_mat->m,
        od->vel.v->ve, 1,
        &M_PEEK(od->mal.y_scores, 0, i), 1);
    }
    /*
    d = ro d
    */
    cblas_dscal(od->mal.y_loadings->m,
      od->vel.ro->ve[i],
      &M_PEEK(od->mal.y_loadings, 0, i), 1);

    /*
    b = E[i]'u / u'u
    */
    cblas_dgemv(CblasColMajor, CblasTrans,
      od->mal.e_mat->m,
      od->mal.e_mat->n, 1.0,
      od->mal.e_mat->base,
      od->mal.e_mat->max_m,
      &M_PEEK(od->mal.x_scores, 0, i), 1, 0.0,
      &M_PEEK(od->mal.x_loadings, 0, i), 1);
    cblas_dscal(od->mal.x_loadings->m, 1.0 / ss_u,
      &M_PEEK(od->mal.x_loadings, 0, i), 1);

    /*
    E[i + 1] = E[i] - ub'
    F[i + 1] = F[i] - ud'
    */
    for (y = 0; y < od->mal.e_mat->m; ++y) {
      cblas_daxpy(od->mal.e_mat->n,
        - M_PEEK(od->mal.x_scores, y, i),
        &M_PEEK(od->mal.x_loadings, 0, i), 1,
        &M_PEEK(od->mal.e_mat, y, 0),
        od->mal.e_mat->max_m);
      cblas_daxpy(od->mal.f_mat->n,
        - M_PEEK(od->mal.x_scores, y, i),
        &M_PEEK(od->mal.y_loadings, 0, i), 1,
        &M_PEEK(od->mal.f_mat, y, 0),
        od->mal.f_mat->max_m);
    }
  }
}
