/*

scramble.c

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


int prepare_scrambling(O3Data *od)
{
  int i;
  int j;
  int k;
  int n;
  int x;
  int y;
  int object_num;
  int struct_num;
  int active_struct_num;
  int conf_num;
  int n_conf;
  int active_struct_count;
  int actual_len;
  int found;
  int top_outgap_len;
  int bottom_outgap_len;
  int bins;
  int structs_per_bin;
  int excess_structs;
  int random_index;
  int positions_to_be_assigned;
  int start;
  double sumweight;
  
  
  get_attr_struct_ave(od, 0, ACTIVE_BIT, &active_struct_num, NULL);
  od->vel.y_values_ave = double_vec_resize
    (od->vel.y_values_ave, active_struct_num);
  if (!(od->vel.y_values_ave)) {
    return OUT_OF_MEMORY;
  }
  memset(od->vel.y_values_ave->ve, 0,
    active_struct_num * sizeof(double));
  od->pel.scrambling_temp = int_perm_resize
    (od->pel.scrambling_temp, active_struct_num);
  if (!(od->pel.scrambling_temp)) {
    return OUT_OF_MEMORY;
  }
  od->pel.scrambling_order = int_perm_resize
    (od->pel.scrambling_order, active_struct_num);
  if (!(od->pel.scrambling_order)) {
    return OUT_OF_MEMORY;
  }
  /*
  allocate a od->active_object_num
  * od->y_vars ordered F matrix
  */
  od->mal.ordered_f_mat = double_mat_resize
    (od->mal.ordered_f_mat, od->active_object_num, od->y_vars);
  if (!(od->mal.ordered_f_mat)) {
    return OUT_OF_MEMORY;
  }
  /*
  if multiple y variables are present,
  then they will be ordered according to their
  average value for each object
  */
  object_num = 0;
  i = 0;
  while (object_num < od->object_num) {
    struct_num = od->al.mol_info[object_num]->struct_num;
    conf_num = 1;
    for (n_conf = 0, y = 0, sumweight = 0.0; n_conf < conf_num; ++n_conf, ++object_num) {
      if (!get_object_attr(od, object_num, ACTIVE_BIT)) {
        continue;
      }
      for (j = 0; j < od->y_vars; ++j) {
        od->vel.y_values_ave->ve[i] += (od->mel.object_weight[object_num]
          * get_y_value(od, object_num, j, WEIGHT_BIT));
        sumweight += od->mel.object_weight[object_num];
      }
      ++y;
    }
    if (y) {
      od->vel.y_values_ave->ve[i] /= sumweight;
      ++i;
    }
  }
  /*
  sort y average values in increasing order
  */
  double_vec_sort(od->vel.y_values_ave, od->pel.scrambling_temp);
  /*
  flip the permutation vector order upside down
  so to have y values in decreasing order
  as in J. Comput.-Aided Mol. Des. 2004, 18, 563-576
  */
  for (i = 0; i < od->pel.scrambling_temp->size; ++i) {
    od->pel.scrambling_order->pe[i] =
      od->pel.scrambling_temp->pe
      [od->pel.scrambling_temp->size - i - 1];
  }
  /*
  copy values from active objects, y_vars
  and store them in decreasing order in ordered_f_mat
  */
  for (n = 0, j = 0; n < od->pel.scrambling_order->size; ++n) {
    active_struct_count = 0;
    object_num = 0;
    found = 0;
    while ((!found) && (object_num < od->object_num)) {
      struct_num = od->al.mol_info[object_num]->struct_num;
      conf_num = 1;
      for (n_conf = 0, found = 0, y = 0; n_conf < conf_num; ++n_conf, ++object_num) {
        if (get_object_attr(od, object_num, ACTIVE_BIT)) {
          if (active_struct_count == od->pel.scrambling_order->pe[n]) {
            found = 1;
            break;
          }
          ++y;
        }
      }
      if (y && (!found)) {
        ++active_struct_count;
      }
    }
    if (found) {
      while (n_conf < conf_num) {
        if (get_object_attr(od, object_num, ACTIVE_BIT)) {
          for (k = 0; k < od->y_vars; ++k) {
            M_POKE(od->mal.ordered_f_mat, j, k,
              get_y_value(od, object_num, k, WEIGHT_BIT));
          }
          ++j;
        }
        ++n_conf;
        ++object_num;
      }
    }
  }
  /*
  starting number of bins
  */
  od->mel.bin_populations = alloc_int_array
    (od->mel.bin_populations, od->scramble.max_bins);
  if (!(od->mel.bin_populations)) {
    return OUT_OF_MEMORY;
  }
  od->mel.candidate_pos = alloc_int_array
    (od->mel.candidate_pos, active_struct_num);
  if (!(od->mel.candidate_pos)) {
    return OUT_OF_MEMORY;
  }
  /*
  create a temporary file to store scrambled orders
  */
  /*
  write the initial order of objects
  */
  actual_len = fwrite(od->pel.scrambling_order->pe,
    sizeof(int), active_struct_num,
    od->file[TEMP_SCRAMBLE]->handle);
  if (actual_len != active_struct_num) {
    return CANNOT_WRITE_TEMP_FILE;
  }
  set_random_seed(od, od->random_seed);
  for (bins = od->scramble.max_bins;
    bins >= od->scramble.min_bins; --bins) {
    structs_per_bin = active_struct_num / bins;
    excess_structs = active_struct_num % bins;
    top_outgap_len = (bins - excess_structs) / 2;
    bottom_outgap_len = top_outgap_len;
    if ((bins - excess_structs) % 2) {
      ++top_outgap_len;
    }
    for (j = 0; j < bins; ++j) {
      od->mel.bin_populations[j] = structs_per_bin;
    }
    for (j = top_outgap_len; j < (bins - bottom_outgap_len); ++j) {
      ++(od->mel.bin_populations[j]);
    }
    for (j = 0; j < od->scramble.scramblings; ++j) {
      x = 0;
      for (k = 0; k < bins; ++k) {
        for (n = 0; n < od->mel.bin_populations[k]; ++n) {
          od->mel.candidate_pos[n] = n;
        }
        positions_to_be_assigned = od->mel.bin_populations[k];
        for (n = 0, start = 0; n < k; ++n) {
          start += od->mel.bin_populations[n];
        }
        while (positions_to_be_assigned) {
          random_index = (int)safe_rint
            ((double)(positions_to_be_assigned - 1) * genrand_real(od));
          --positions_to_be_assigned;
          od->pel.scrambling_temp->pe[x] = od->pel.scrambling_order->pe
            [start + od->mel.candidate_pos[random_index]];
          ++x;
          for (n = random_index; n < positions_to_be_assigned; ++n) {
            od->mel.candidate_pos[n] = od->mel.candidate_pos[n + 1];
          }
        }
      }
      actual_len = fwrite(od->pel.scrambling_temp->pe,
        sizeof(int), active_struct_num,
        od->file[TEMP_SCRAMBLE]->handle);
      if (actual_len != active_struct_num) {
        return CANNOT_WRITE_TEMP_FILE;
      }
    }
  }
  rewind(od->file[TEMP_SCRAMBLE]->handle);
  
  return 0;
}


int scramble(O3Data *od, int pc_num)
{
  int bins;
  int i;
  int j;
  int k;
  int x;
  int result;
  int scramble_run;
  int info;
  #if (!defined HAVE_LIBLAPACK_ATLAS) && (!defined HAVE_LIBSUNPERF)
  int lwork;
  #endif
  double r2yy;
  double r2_num;
  double r2_den1;
  double r2_den2;
  double critical_fit[2];
  
  
  od->mal.press = double_mat_resize
    (od->mal.press, pc_num + 1, od->y_vars);
  if (!(od->mal.press)) {
    return OUT_OF_MEMORY;
  }
  od->scramble.overall_scramblings =
    (od->scramble.max_bins - od->scramble.min_bins + 1)
    * od->scramble.scramblings;
  od->mal.r2_fit_mat = double_mat_resize(od->mal.r2_fit_mat,
    od->scramble.overall_scramblings, od->scramble.fit_order + 1);
  if (!(od->mal.r2_fit_mat)) {
    return OUT_OF_MEMORY;
  }
  od->mal.fit_temp = double_mat_resize(od->mal.fit_temp,
    od->scramble.fit_order + 1, od->scramble.fit_order + 1);
  if (!(od->mal.fit_temp)) {
    return OUT_OF_MEMORY;
  }
  memset(od->mal.r2_fit_mat->base, 0,
    od->mal.r2_fit_mat->m * od->mal.r2_fit_mat->n * sizeof(double));
  for (i = 0; i < 2; ++i) {
    od->vel.fit_vec[i] = double_vec_resize
      (od->vel.fit_vec[i], od->scramble.overall_scramblings);
    if (!(od->vel.fit_vec[i])) {
      return OUT_OF_MEMORY;
    }
    od->vel.fit_coeff[i] = double_vec_resize
      (od->vel.fit_coeff[i], od->scramble.fit_order + 1);
    if (!(od->vel.fit_coeff[i])) {
      return OUT_OF_MEMORY;
    }
    od->vel.fit_prod[i] = double_vec_resize
      (od->vel.fit_prod[i], od->scramble.fit_order + 1);
    if (!(od->vel.fit_prod[i])) {
      return OUT_OF_MEMORY;
    }
  }
  result = alloc_pls(od, od->overall_active_x_vars, pc_num, CV_MODEL);
  if (result) {
    return OUT_OF_MEMORY;
  }
  result = alloc_cv_sdep(od, pc_num, od->scramble.runs);
  if (result) {
    return OUT_OF_MEMORY;
  }
  set_random_seed(od, od->random_seed);
  result = prepare_cv(od, pc_num, od->scramble.cv_type,
    od->scramble.groups, od->scramble.runs);
  if (result) {
    return result;
  }
  result = prepare_scrambling(od);
  if (result) {
    return result;
  }
  result = fill_x_matrix_scrambled(od);
  if (result) {
    return result;
  }
  result = alloc_average_mat(od, SCRAMBLE_CV_MODEL,
    od->scramble.cv_type, od->scramble.groups,
    od->scramble.runs);
  if (result) {
    return OUT_OF_MEMORY;
  }
  scramble_run = 0;
  if (od->scramble.print_runs) {
    tee_printf(od, "------------------------------------------------\n"
      "%12s%12s%12s%12s\n"
      "------------------------------------------------\n",
      "Scramble run", "r2(yy')", "q2", "SE(cv)");
  }
  for (bins = od->scramble.max_bins;
    bins >= od->scramble.min_bins; --bins) {
    for (j = 0; j < od->scramble.scramblings; ++j) {
      result = fill_y_matrix_scrambled(od);
      if (result) {
        return result;
      }
      memset(od->mal.ave_press->base, 0, od->mal.ave_press->m
        * od->mal.ave_press->n * sizeof(double));
      if (od->n_proc > 1) {
        result = parallel_cv(od, od->overall_active_x_vars,
          pc_num, PARALLEL_CV | SCRAMBLE_CV_MODEL,
          od->scramble.cv_type, od->scramble.groups,
          od->scramble.runs);
        if (result) {
          return result;
        }
      }
      else {
        result = cv(od, pc_num, SCRAMBLE_CV_MODEL,
          od->scramble.cv_type, od->scramble.groups,
          od->scramble.runs);
        if (result) {
          return result;
        }
      }
      if (od->pc_num != pc_num) {
        od->vel.fit_vec[Q2_FIT]->ve[scramble_run] = -10.0;
        od->vel.fit_vec[SECV_FIT]->ve[scramble_run] = 10.0;
      }
      else {
        od->vel.fit_vec[Q2_FIT]->ve[scramble_run] =
          od->vel.q2->ve[pc_num];
        od->vel.fit_vec[SECV_FIT]->ve[scramble_run] =
          od->vel.secv->ve[pc_num];
      }
      r2_num = 0.0;
      r2_den1 = 0.0;
      r2_den2 = 0.0;
      r2yy = 0.0;
      M_POKE(od->mal.r2_fit_mat, scramble_run, 0, 1.0);
      for (x = 0; x < od->y_vars; ++x) {
        for (k = 0; k < od->active_object_num; ++k) {
          r2_num += ((M_PEEK(od->mal.ordered_f_mat, k, x)
            - od->vel.f_mat_full_ave->ve[x])
            * (M_PEEK(od->mal.large_f_mat, k, x)
            - od->vel.f_mat_full_ave->ve[x]));
          r2_den1 += square
            (M_PEEK(od->mal.ordered_f_mat, k, x)
            - od->vel.f_mat_full_ave->ve[x]);
          r2_den2 += square
            (M_PEEK(od->mal.large_f_mat, k, x)
            - od->vel.f_mat_full_ave->ve[x]);
        }
        r2yy += square(r2_num) / (r2_den1 * r2_den2);
      }
      r2yy /= od->y_vars;
      M_POKE(od->mal.r2_fit_mat, scramble_run, 1, r2yy);
      M_POKE(od->mal.r2_fit_mat, scramble_run, 2, square(r2yy));
      if (od->scramble.fit_order == 3) {
        M_POKE(od->mal.r2_fit_mat, scramble_run, 3, pow(r2yy, 3));
      }
      if (od->scramble.print_runs) {
        tee_printf(od, "%12d%12.4lf%12.4lf%12.4lf\n", scramble_run + 1, r2yy,
          od->vel.fit_vec[Q2_FIT]->ve[scramble_run],
          od->vel.fit_vec[SECV_FIT]->ve[scramble_run]);
      }
      ++scramble_run;
    }
  }
  if (od->scramble.print_runs) {
    tee_printf(od, "\n");
  }
  cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
    od->mal.r2_fit_mat->n,
    od->mal.r2_fit_mat->n,
    od->mal.r2_fit_mat->m, 1.0,
    od->mal.r2_fit_mat->base,
    od->mal.r2_fit_mat->max_m,
    od->mal.r2_fit_mat->base,
    od->mal.r2_fit_mat->max_m, 0.0,
    od->mal.fit_temp->base,
    od->mal.fit_temp->max_m);
  #ifdef HAVE_LIBMKL
  dgetrf(&(od->mal.fit_temp->m),
    &(od->mal.fit_temp->m),
    od->mal.fit_temp->base,
    &(od->mal.fit_temp->max_m),
    od->mel.ipiv, &info);
  #elif HAVE_LIBSUNPERF
  dgetrf(od->mal.fit_temp->m,
    od->mal.fit_temp->m,
    od->mal.fit_temp->base,
    od->mal.fit_temp->max_m,
    od->mel.ipiv, &info);
  #elif HAVE_LIBLAPACK_ATLAS
  info = clapack_dgetrf(CblasColMajor,
    od->mal.fit_temp->m,
    od->mal.fit_temp->m,
    od->mal.fit_temp->base,
    od->mal.fit_temp->max_m,
    od->mel.ipiv);
  #else
  dgetrf_(&(od->mal.fit_temp->m),
    &(od->mal.fit_temp->m),
    od->mal.fit_temp->base,
    &(od->mal.fit_temp->max_m),
    od->mel.ipiv, &info);
  #endif
  if (!info) {
    #ifdef HAVE_LIBMKL
    lwork = (od->mal.fit_temp->n + 1) * LWORK_BLOCK_SIZE * sizeof(double);
    dgetri(&(od->mal.fit_temp->m),
      od->mal.fit_temp->base,
      &(od->mal.fit_temp->max_m),
      od->mel.ipiv,
      od->mel.work, &lwork, &info);
    #elif HAVE_LIBSUNPERF
    dgetri(od->mal.fit_temp->m,
      od->mal.fit_temp->base,
      od->mal.fit_temp->max_m,
      od->mel.ipiv, &info);
    #elif HAVE_LIBLAPACK_ATLAS
    info = clapack_dgetri(CblasColMajor,
      od->mal.fit_temp->m,
      od->mal.fit_temp->base,
      od->mal.fit_temp->max_m,
      od->mel.ipiv);
    #else
    lwork = (od->mal.fit_temp->n + 1) * LWORK_BLOCK_SIZE * sizeof(double);
    dgetri_(&(od->mal.fit_temp->m),
      od->mal.fit_temp->base,
      &(od->mal.fit_temp->max_m),
      od->mel.ipiv,
      od->mel.work, &lwork, &info);
    #endif
  }
  for (i = 0; i < 2; ++i) {
    cblas_dgemv(CblasColMajor, CblasTrans,
      od->mal.r2_fit_mat->m,
      od->mal.r2_fit_mat->n, 1.0,
      od->mal.r2_fit_mat->base,
      od->mal.r2_fit_mat->max_m,
      od->vel.fit_vec[i]->ve, 1, 0.0,
      od->vel.fit_prod[i]->ve, 1);
    cblas_dgemv(CblasColMajor, CblasNoTrans,
      od->mal.fit_temp->m,
      od->mal.fit_temp->n, 1.0,
      od->mal.fit_temp->base,
      od->mal.fit_temp->max_m,
      od->vel.fit_prod[i]->ve, 1, 0.0,
      od->vel.fit_coeff[i]->ve, 1);
  }
  for (i = 0; i < 2; ++i) {
    critical_fit[i] = 0.0;
    for (j = 0; j <= od->scramble.fit_order; ++j) {
      critical_fit[i] += od->vel.fit_coeff[i]->ve[j]
        * pow(od->scramble.critical_point, j);
    }
  }
  tee_printf(od, "--------------------------------------------------------------------\n"
    "%16s%16s%16s%20s\n"
    "--------------------------------------------------------------------\n"
    "%16.4lf%16d%16.4lf%20.4lf\n\n",
    "Critical r2(yy')", "Fit order", "Fitted q2", "Fitted SE(cv)",
    od->scramble.critical_point,
    od->scramble.fit_order,
    critical_fit[Q2_FIT], critical_fit[SECV_FIT]);
  if (od->vel.y_values_ave) {
    double_vec_free(od->vel.y_values_ave);
    od->vel.y_values_ave = NULL;
  }
  for (i = 0; i < 2; ++i) {
    if (od->vel.fit_prod[i]) {
      double_vec_free(od->vel.fit_prod[i]);
      od->vel.fit_prod[i] = NULL;
    }
  }
  if (od->pel.scrambling_temp) {
    int_perm_free(od->pel.scrambling_temp);
    od->pel.scrambling_temp = NULL;
  }
  if (od->pel.scrambling_order) {
    int_perm_free(od->pel.scrambling_order);
    od->pel.scrambling_order = NULL;
  }
  if (od->mal.ordered_f_mat) {
    double_mat_free(od->mal.ordered_f_mat);
    od->mal.ordered_f_mat = NULL;
  }
  if (od->mal.large_e_mat_ave) {
    double_mat_free(od->mal.large_e_mat_ave);
    od->mal.large_e_mat_ave = NULL;
  }
  if (od->mel.bin_populations) {
    free(od->mel.bin_populations);
    od->mel.bin_populations = NULL;
  }
  if (od->mel.candidate_pos) {
    free(od->mel.candidate_pos);
    od->mel.candidate_pos = NULL;
  }
  if (od->n_proc > 1) {
    free_parallel_cv(od, od->mel.thread_info,
      SCRAMBLE_CV_MODEL, od->scramble.cv_type,
      od->scramble.runs);
  }
  if (od->scramble.cv_type == LEAVE_MANY_OUT) {
    free_cv_groups(od, od->scramble.runs);
  }
  od->valid |= SCRAMBLE_BIT;

  return 0;
}
