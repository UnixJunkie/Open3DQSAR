/*

pca.c

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


int pca(O3Data *od, int pc_num)
{
  int i;
  int j;
  int x;
  int y;
  int x_max;
  int y_max;
  int conv;
  int max_ss_x;
  int cmp;
  int actual_len;
  double norm_b;
  double max_ss;
  double ss_u;
  double ss_b;
  double ss_u_old;
  double s2_x_tot;
  double s2_x_tot_start;
  double cum_explained_s2_x;
  double score;


  ss_u_old = 0;
  x_max = od->overall_active_x_vars;
  y_max = od->active_object_num;
  od->mal.e_mat =
    double_mat_resize(od->mal.e_mat, y_max, x_max);
  if (!(od->mal.e_mat)) {
    return OUT_OF_MEMORY;
  }
  od->mal.x_weights =
    double_mat_resize(od->mal.x_weights, x_max, pc_num + 1);
  if (!(od->mal.x_weights)) {
    return OUT_OF_MEMORY;
  }
  od->vel.explained_s2_x = double_vec_resize
    (od->vel.explained_s2_x, pc_num);
  if (!(od->vel.explained_s2_x)) {
    return OUT_OF_MEMORY;
  }
  od->vel.u =
    double_vec_resize(od->vel.u, y_max);
  if (!(od->vel.u)) {
    return OUT_OF_MEMORY;
  }
  od->vel.b =
    double_vec_resize(od->vel.b, x_max);
  if (!(od->vel.b)) {
    return OUT_OF_MEMORY;
  }

  /*
  calculate total variance of this matrix
  */
  for (x = 0, s2_x_tot_start = 0.0; x < od->mal.e_mat->n; ++x) {
    s2_x_tot_start += cblas_ddot(od->mal.e_mat->m,
      &M_PEEK(od->mal.e_mat, 0, x), 1, &M_PEEK(od->mal.e_mat, 0, x), 1);
  }
  cum_explained_s2_x = 0.0;
  tee_printf(od, "  %20s%20s\n", "Explained", "Cum. explained");
  tee_printf(od, "PC%20s%20s\n", "variance X %", "variance X %");
  tee_printf(od, "------------------------------------------\n");
  fwrite(&pc_num, sizeof(int), 1,
    od->file[TEMP_PCA_LOADINGS]->handle);
  actual_len = fwrite(&(od->mal.e_mat->n), sizeof(int), 1,
    od->file[TEMP_PCA_LOADINGS]->handle);
  if (actual_len != 1) {
    return CANNOT_WRITE_TEMP_FILE;
  }
  fwrite(&pc_num, sizeof(int), 1,
    od->file[TEMP_PCA_SCORES]->handle);
  actual_len = fwrite(&(od->mal.e_mat->m), sizeof(int), 1,
    od->file[TEMP_PCA_SCORES]->handle);
  if (actual_len != 1) {
    return CANNOT_WRITE_TEMP_FILE;
  }
    
        for (i = 0; i < pc_num; ++i) {
    /*
    Find the column vector with the largest sum of squares
    */
    max_ss = 0.0;
    max_ss_x = 0;
    for (x = 0; x < x_max; ++x) {
      /*
      extract x-th column vector from X
      */
      cblas_dcopy(od->mal.e_mat->m,
        &M_PEEK(od->mal.e_mat, 0, x), 1,
        od->vel.u->ve, 1);
      ss_u = cblas_ddot(od->vel.u->size,
        od->vel.u->ve, 1,
        od->vel.u->ve, 1);
      if (ss_u > max_ss) {
        max_ss = ss_u;
        max_ss_x = x;
      }
    }
    cblas_dcopy(od->mal.e_mat->m,
      &M_PEEK(od->mal.e_mat, 0, max_ss_x), 1,
      od->vel.u->ve, 1);
    conv = 0;
    cmp = 0;
    while (!conv) {
      /*
      b = (X[i]'t) / (t't)
      */
      cblas_dgemv(CblasColMajor, CblasTrans,
        od->mal.e_mat->m,
        od->mal.e_mat->n, 1.0,
        od->mal.e_mat->base,
        od->mal.e_mat->max_m,
        od->vel.u->ve,
        1, 0.0, od->vel.b->ve, 1);
      ss_u = cblas_ddot(od->vel.u->size,
        od->vel.u->ve, 1,
        od->vel.u->ve, 1);
      cblas_dscal(od->vel.b->size,
        1.0 / ss_u, od->vel.b->ve, 1);
      
      /*
      b = b * (b'b)^(-0.5)
      */
      norm_b = cblas_dnrm2(od->vel.b->size,
        od->vel.b->ve, 1);
      cblas_dscal(od->vel.b->size,
        1.0 / norm_b, od->vel.b->ve, 1);
    
      /*
      u = (E[i]b) / (b'b)
      */
      cblas_dgemv(CblasColMajor, CblasNoTrans,
        od->mal.e_mat->m,
        od->mal.e_mat->n, 1.0,
        od->mal.e_mat->base,
        od->mal.e_mat->max_m,
        od->vel.b->ve,
        1, 0.0, od->vel.u->ve, 1);
      ss_b = cblas_ddot(od->vel.b->size,
        od->vel.b->ve, 1,
        od->vel.b->ve, 1);
      cblas_dscal(od->vel.u->size,
        1.0 / ss_b, od->vel.u->ve, 1);
      
      /*
      check convergence
      */
      ss_u = cblas_ddot(od->vel.u->size,
        od->vel.u->ve, 1,
        od->vel.u->ve, 1);

      if (!cmp) {
        cmp = 1;
      }
      else {
        if (fabs(ss_u - ss_u_old)
          < (PCA_CONV_THRESHOLD * ss_u)) {
          conv = 1;
        }
      }
      ss_u_old = ss_u;
    }
    /*
    E[i + 1] = E[i] - ub'
    */
    for (y = 0; y < y_max; ++y) {
      cblas_daxpy(od->mal.e_mat->n,
        - od->vel.u->ve[y],
        od->vel.b->ve, 1,
        &M_PEEK(od->mal.e_mat, y, 0),
        od->mal.e_mat->max_m);
    }
    actual_len = fwrite(od->vel.u->ve,
      sizeof(double), od->mal.e_mat->m,
      od->file[TEMP_PCA_SCORES]->handle);
    if (actual_len != od->mal.e_mat->m) {
      return CANNOT_WRITE_TEMP_FILE;
    }
    actual_len = fwrite(od->vel.b->ve,
      sizeof(double), od->mal.e_mat->n,
      od->file[TEMP_PCA_LOADINGS]->handle);
    if (actual_len != od->mal.e_mat->n) {
      return CANNOT_WRITE_TEMP_FILE;
    }
    for (x = 0, s2_x_tot = 0.0; x < od->mal.e_mat->n; ++x) {
      s2_x_tot += cblas_ddot(od->mal.e_mat->m,
        &M_PEEK(od->mal.e_mat, 0, x), 1, &M_PEEK(od->mal.e_mat, 0, x), 1);
    }
    od->vel.explained_s2_x->ve[i] =
      1.0 - cum_explained_s2_x
      - s2_x_tot / s2_x_tot_start;
    cum_explained_s2_x +=
      od->vel.explained_s2_x->ve[i];
    tee_printf(od, "%2d%20.4lf%20.4lf\n", i + 1,
      od->vel.explained_s2_x->ve[i] * 100,
      cum_explained_s2_x * 100);
  }
  tee_printf(od, "\nPCA X scores");
  tee_printf(od, "\n--------------------------------------------");
  for (j = 0; j < pc_num; ++j) {
    tee_printf(od, "------------");
  }
  tee_printf(od, "\n");
  tee_printf(od, "%4s    %-36s", "N", "Name");
  for (j = 0; j < pc_num; ++j) {
    tee_printf(od, "%12d", j + 1);
  }
  tee_printf(od, "\n--------------------------------------------");
  for (j = 0; j < pc_num; ++j) {
    tee_printf(od, "------------");
  }
  tee_printf(od, "\n");
  y = 0;
  for (i = 0; i < od->object_num; ++i) {
    if (get_object_attr(od, i, ACTIVE_BIT)) {
      tee_printf(od, "%4d    %-36s", i + 1, od->al.mol_info[i]->object_name);
      for (j = 0; j < pc_num; ++j) {
        if (fseek(od->file[TEMP_PCA_SCORES]->handle,
          (j * od->mal.e_mat->m + y) * sizeof(double)
          + 2 * sizeof(int), SEEK_SET)) {
          return CANNOT_READ_TEMP_FILE;
        }
        actual_len = fread(&score, sizeof(double), 1,
          od->file[TEMP_PCA_SCORES]->handle);
        if (actual_len != 1) {
          return CANNOT_READ_TEMP_FILE;
        }
        tee_printf(od, "%12.4lf", score);
      }
      tee_printf(od, "\n");
      ++y;
    }
  }
  tee_printf(od, "\n");
  if (od->mal.e_mat) {
    double_mat_free(od->mal.e_mat);
    od->mal.e_mat = NULL;
  }
  if (od->mal.x_weights) {
    double_mat_free(od->mal.x_weights);
    od->mal.x_weights = NULL;
  }
  if (od->vel.u) {
    double_vec_free(od->vel.u);
    od->vel.u = NULL;
  }
  if (od->vel.b) {
    double_vec_free(od->vel.b);
    od->vel.b = NULL;
  }
  
  return 0;
}
