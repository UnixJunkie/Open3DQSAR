/*

d_optimal.c

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


double calc_delta_ij(O3Data *od,
  DoubleMat *dispersion_mat, int i, int j)
{
  double dxi;
  double dxj;
  double dxixj;

  
  dxi = calc_dxi(od, dispersion_mat,
    od->mal.design_mat, i);
  dxj = calc_dxi(od, dispersion_mat,
    od->mal.support_mat, j);
  dxixj = calc_dxixj(od, dispersion_mat, i, j);

  return (dxj - (dxi * dxj - square(dxixj)) - dxi);
}


DoubleMat *calc_dispersion_matrix(O3Data *od,
  DoubleMat *candidates_mat)
{
  int info;
  #if (!defined HAVE_LIBLAPACK_ATLAS) && (!defined HAVE_LIBSUNPERF)
  int lwork;
  #endif
  
  
  cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
    candidates_mat->n,
    candidates_mat->n,
    candidates_mat->m, 1.0,
    candidates_mat->base,
    candidates_mat->max_m,
    candidates_mat->base,
    candidates_mat->max_m, 0.0,
    od->mal.temp->base,
    od->mal.temp->max_m);
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
    lwork = (candidates_mat->n + 1) * LWORK_BLOCK_SIZE * sizeof(double);
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
    lwork = (candidates_mat->n + 1) * LWORK_BLOCK_SIZE * sizeof(double);
    dgetri_(&(od->mal.temp->m),
      od->mal.temp->base,
      &(od->mal.temp->max_m),
      od->mel.ipiv,
      od->mel.work, &lwork, &info);
    #endif
  }
  
  return od->mal.temp;
}


double calc_dxi(O3Data *od, DoubleMat *dispersion_mat,
  DoubleMat *candidates_mat, int i)
{
  int n;
  double dxi;

  
  for (n = 0; n < candidates_mat->n; ++n) {
    od->vel.xn_vec->ve[n] = M_PEEK(candidates_mat, i, n);
  }
  cblas_dgemv(CblasColMajor, CblasNoTrans,
    dispersion_mat->m,
    dispersion_mat->n, 1.0,
    dispersion_mat->base,
    dispersion_mat->max_m,
    od->vel.xn_vec->ve,
    1, 0.0, od->vel.temp_xn_vec->ve, 1);
  dxi = cblas_ddot(od->vel.xn_vec->size,
    od->vel.xn_vec->ve, 1,
    od->vel.temp_xn_vec->ve, 1);
  
  return dxi;
}


double calc_dxixj(O3Data *od,
  DoubleMat *dispersion_mat, int i, int j)
{
  int n;
  double dxixj;
  
  
  for (n = 0; n < od->mal.design_mat->n; ++n) {
    od->vel.xn_vec->ve[n] =
      M_PEEK(od->mal.design_mat, i, n);
  }
  cblas_dgemv(CblasColMajor, CblasNoTrans,
    dispersion_mat->m,
    dispersion_mat->n, 1.0,
    dispersion_mat->base,
    dispersion_mat->max_m,
    od->vel.xn_vec->ve,
    1, 0.0, od->vel.temp_xn_vec->ve, 1);
  for (n = 0; n < od->mal.support_mat->n; ++n) {
    od->vel.xn_vec->ve[n] =
      M_PEEK(od->mal.support_mat, j, n);
  }
  dxixj = cblas_ddot(od->vel.xn_vec->size,
    od->vel.xn_vec->ve, 1,
    od->vel.temp_xn_vec->ve, 1);
  
  return dxixj;
}


int create_design_support_matrices(O3Data *od,
  DoubleMat *candidates_mat, int design_points)
{
  int i;
  int j;
  int n;
  int candidates;
  int random_candidate;
  int factors;
  int support_points;
  
  
  candidates = candidates_mat->m;
  factors = candidates_mat->n;
  support_points = candidates - design_points;
  od->mal.sorted_candidates_mat = double_mat_resize
    (od->mal.sorted_candidates_mat,
    candidates, factors);
  if (!(od->mal.sorted_candidates_mat)) {
    return OUT_OF_MEMORY;
  }
  od->mal.design_mat = double_mat_resize
    (od->mal.design_mat, design_points, factors);
  if (!(od->mal.design_mat)) {
    return OUT_OF_MEMORY;
  }
  od->mal.support_mat = double_mat_resize
    (od->mal.support_mat, support_points, factors);
  if (!(od->mal.support_mat)) {
    return OUT_OF_MEMORY;
  }
  od->pel.ordered_var_list = int_perm_resize
    (od->pel.ordered_var_list, candidates);
  if (!(od->pel.ordered_var_list)) {
    return OUT_OF_MEMORY;
  }
    
  for (i = 0; i < candidates; ++i) {
    od->pel.ordered_var_list->pe[i] = i;
  }
  n = candidates;
  set_random_seed(od, od->random_seed);
  for (i = 0; i < candidates; ++i) {
    random_candidate = (int)safe_rint
      ((double)(n - 1) * genrand_real(od));
    od->pel.dxn_perm->pe[i] =
      od->pel.ordered_var_list->pe[random_candidate];
    for (j = random_candidate; j < (candidates - 1); ++j) {
      od->pel.ordered_var_list->pe[j] =
        od->pel.ordered_var_list->pe[j + 1];
    }
    --n;
  }

  (void)int_perm_rows(od->pel.dxn_perm,
    candidates_mat, od->mal.sorted_candidates_mat);
  
  for (i = 0; i < factors; ++i) {
    memcpy(&M_PEEK(od->mal.support_mat, 0, i),
      &M_PEEK(od->mal.sorted_candidates_mat, 0, i),
      support_points * sizeof(double));
    memcpy(&M_PEEK(od->mal.design_mat, 0, i),
      &M_PEEK(od->mal.sorted_candidates_mat,
      support_points, i),
      design_points * sizeof(double));
  }
  double_mat_free(od->mal.sorted_candidates_mat);
  od->mal.sorted_candidates_mat = NULL;

  return 0;
}


int k_exchange(O3Data *od, DoubleMat *dispersion_mat)
{
  int i;
  int j;
  int k;
  int design_points;
  int support_points;
  int factors;
  int iterations;
  int delta_over_threshold;
  int old_index;
  int worst_i;
  int best_j;
  
  
  design_points = od->mal.design_mat->m;
  support_points = od->mal.support_mat->m;
  factors = dispersion_mat->m;
  iterations = 0;
  delta_over_threshold = 0;
  k = design_points / 16;
  if (!k) {
    k = 1;
  }
  set_random_seed(od, od->random_seed);
  while (iterations < MAX_K_EXCHANGE_ITERATIONS) {
    /*
    calculate dxi for all rows in design_mat
    and store the values in dxi_vec
    */
    for (i = 0; i < design_points; ++i) {
      od->vel.dxi_vec->ve[i] =
        calc_dxi(od, dispersion_mat,
        od->mal.design_mat, i);
    }
    (void)double_vec_sort(od->vel.dxi_vec, od->pel.dxi_perm);
    for (i = 0; i < k; ++i) {
      for (j = 0; j < support_points; ++j) {
        od->vel.delta_vec->ve[j] =
          calc_delta_ij(od, dispersion_mat,
          od->pel.dxi_perm->pe[i], j);
      }
      (void)double_vec_sort(od->vel.delta_vec,
        od->pel.delta_perm);
      if (od->vel.delta_vec->ve[support_points - 1]
        > MAX_DELTA_THRESHOLD) {
        delta_over_threshold = 1;
        /*
        update dxn_perm
        */
        worst_i = od->pel.dxi_perm->pe[i];
        best_j = od->pel.delta_perm->pe[support_points - 1];
        old_index = od->pel.dxn_perm->pe[worst_i + support_points];
        od->pel.dxn_perm->pe[worst_i + support_points] =
          od->pel.dxn_perm->pe[best_j];
        od->pel.dxn_perm->pe[best_j] = old_index;
        /*
        exchange vectors between design and support matrices
        */
        cblas_dcopy(factors, &M_PEEK(od->mal.design_mat,
          worst_i, 0), od->mal.design_mat->max_m,
          od->vel.xn_vec->ve, 1);
        cblas_dcopy(factors, &M_PEEK(od->mal.support_mat,
          best_j, 0), od->mal.support_mat->max_m,
          od->vel.temp_xn_vec->ve, 1);
        cblas_dcopy(factors, od->vel.xn_vec->ve, 1,
          &M_PEEK(od->mal.support_mat,
          best_j, 0),
          od->mal.support_mat->max_m);
        cblas_dcopy(factors, od->vel.temp_xn_vec->ve, 1,
          &M_PEEK(od->mal.design_mat,
          worst_i, 0),
          od->mal.design_mat->max_m);
        dispersion_mat = calc_dispersion_matrix(od,
          od->mal.design_mat);
        if (!dispersion_mat) {
          return SINGULAR_MATRIX;
        }
      }
    }
    if (!delta_over_threshold) {
      break;
    }
    else {
      delta_over_threshold = 0;
      ++iterations;
    }
  }
  
  return 0;
}


int d_optimal(O3Data *od, int factors, int design_points, int type)
{
  int candidates;
  int support_points;
  int result;
  DoubleMat *candidates_mat;
  DoubleMat *dispersion_mat;
  

  candidates = od->overall_active_x_vars;
  result = reload_weights_loadings(od);
  if (result) {
    return result;
  }
  if (type == LOADINGS) {
    od->mal.x_loadings = double_mat_resize
      (od->mal.x_loadings, candidates, factors);
    candidates_mat = od->mal.x_loadings;
  }
  else {
    od->mal.x_weights = double_mat_resize
      (od->mal.x_weights, candidates, factors);
    candidates_mat = od->mal.x_weights;
  }
  od->mal.temp = double_mat_resize
    (od->mal.temp, factors, factors);
  support_points = candidates - design_points;
  od->vel.xn_vec = double_vec_resize
    (od->vel.xn_vec, factors);
  if (!(od->vel.xn_vec)) {
    return OUT_OF_MEMORY;
  }
  od->vel.temp_xn_vec = double_vec_resize
    (od->vel.temp_xn_vec, factors);
  if (!(od->vel.temp_xn_vec)) {
    return OUT_OF_MEMORY;
  }
  od->vel.dxn_vec = double_vec_resize
    (od->vel.dxn_vec, candidates);
  if (!(od->vel.dxn_vec)) {
    return OUT_OF_MEMORY;
  }
  od->vel.dxi_vec = double_vec_resize
    (od->vel.dxi_vec, design_points);
  if (!(od->vel.dxi_vec)) {
    return OUT_OF_MEMORY;
  }
  od->vel.delta_vec = double_vec_resize
    (od->vel.delta_vec, support_points);
  if (!(od->vel.delta_vec)) {
    return OUT_OF_MEMORY;
  }
  od->pel.dxn_perm = int_perm_resize
    (od->pel.dxn_perm, candidates);
  if (!(od->pel.dxn_perm)) {
    return OUT_OF_MEMORY;
  }
  od->pel.dxi_perm = int_perm_resize
    (od->pel.dxi_perm, design_points);
  if (!(od->pel.dxi_perm)) {
    return OUT_OF_MEMORY;
  }
  od->pel.delta_perm = int_perm_resize
    (od->pel.delta_perm, support_points);
  if (!(od->pel.delta_perm)) {
    return OUT_OF_MEMORY;
  }
  result = create_design_support_matrices
    (od, candidates_mat, design_points);
  if (result) {
    return OUT_OF_MEMORY;
  }
  
  dispersion_mat = calc_dispersion_matrix
    (od, od->mal.design_mat);
  if (!dispersion_mat) {
    return SINGULAR_MATRIX;
  }
  result = k_exchange(od, dispersion_mat);

  return result;
}


int d_optimal_select(O3Data *od, int pc_num,
  int design_points, int type)
{
  int i;
  int n;
  int field_num;
  int x_var;
  int actual_len;
  int candidates;
  int support_points;
  int factors;
  int result;
  

  candidates = od->overall_active_x_vars;
  factors = pc_num;
  support_points = candidates - design_points;
  result = d_optimal(od, factors, design_points, type);
  if (result) {
    return result;
  }
  /*
  set D_OPTIMAL bit on variables in the support matrix
  */
  od->file[TEMP_X_MATRIX]->handle =
    fopen(od->file[TEMP_X_MATRIX]->name, "rb");
  if (!(od->file[TEMP_X_MATRIX]->handle)) {
    return CANNOT_READ_TEMP_FILE;
  }
  n = 0;
  while (1) {
    fread(&field_num, sizeof(int), 1,
      od->file[TEMP_X_MATRIX]->handle);
    actual_len = fread(&x_var, sizeof(int), 1,
      od->file[TEMP_X_MATRIX]->handle);
    if (!actual_len) {
      break;
    }
    if (!get_field_attr(od, field_num, ACTIVE_BIT)) {
      continue;
    }
    for (i = 0; i < support_points; ++i) {
      if (n == od->pel.dxn_perm->pe[i]) {
        set_x_var_attr(od,
          field_num, x_var, D_OPTIMAL_BIT, 1);
        break;
      }
    }
    ++n;
  }
  
  return 0;
}
