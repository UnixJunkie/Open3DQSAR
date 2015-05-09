/*

uvepls.c

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


int get_cv_coeff(O3Data *od, int cv_run,
  int y, int x, double *cv_coeff, int save_ram)
{
  int i;
  int actual_len;
  ThreadInfo **thread_info;
  FILE *file;
  
  
  if (save_ram) {
    if (od->n_proc > 1) {
      thread_info = od->mel.thread_info;
      i = 0;
      while((cv_run - thread_info[i]->n_calc
          * od->uvepls.groups) >= 0) {
        cv_run -= thread_info[i]->n_calc
          * od->uvepls.groups;
        ++i;
      }
      file = thread_info[i]->temp_cv_coeff.handle;
    }
    else {
      file = od->file[TEMP_CV_COEFF]->handle;
    }
    fseek(file, (od->mal.b_coefficients->m
      * (y + cv_run * od->y_vars) + x)
      * sizeof(double), SEEK_SET);
    actual_len = fread(cv_coeff, sizeof(double), 1, file);
    if (actual_len != 1) {
      return CANNOT_READ_TEMP_FILE;
    }
  }
  else {
    *cv_coeff = M_PEEK(od->mal.b_coefficients_store,
      y + cv_run * od->y_vars, x);
  }
  
  return 0;
}


int uvepls(O3Data *od, int pc_num)
{
  int i;
  int j;
  int k = 0;
  int n;
  int n_excluded;
  int overall_n_excluded;
  int ive_iterations;
  int coeff_size;
  int actual_len;
  int result;
  int first_half_end;
  int second_half_start;
  int nth_quantile;
  int model_num;
  int included;
  int eliminated;
  int c_count;
  int eff_c_count = 0;
  int lowest_c_pos;
  int overall_seed_count;
  int voronoi_num;
  int cv_run;
  int limit;
  double threshold_c;
  double lower_quartile;
  double higher_quartile;
  DoubleMat *b_coefficients;
  DoubleMat *b_coefficients_ave;
  DoubleMat *b_coefficients_sd;
  DoubleMat *b_coefficients_store = NULL;
  DoubleVec *median_vec = NULL;
  DoubleVec *c;
  DoubleVec *eff_c;
  DoubleVec *eff_dummy_c = NULL;
  DoubleVec *predictivity_list = NULL;
  IntPerm *median_perm = NULL;
  IntPerm *eff_c_rank = NULL;
  IntPerm *eff_dummy_c_rank = NULL;
  IntPerm *predictivity_list_rank = NULL;
  ThreadInfo **thread_info;
  FILE *file;
  
  
  thread_info = od->mel.thread_info;
  coeff_size = 2;
  /*
  if an iterative procedure will be used,
  dummy variables are not necessary
  */
  if (od->uvepls.ive) {
    coeff_size = 1;
  }
  od->uvepls.uvepls_included_vars = set_sel_included_bit(od,
    od->uvepls.use_srd_groups);
  /*
  allocate a vector of characters indicating whether variables
  are included ('1') or not ('0')
  */
  od->mel.uvepls_included = realloc
    (od->mel.uvepls_included,
    od->uvepls.uvepls_included_vars + 4);
  if (!(od->mel.uvepls_included)) {
    return OUT_OF_MEMORY;
  }
  /*
  initially all SEL_INCLUDED variables are included in the model
  */
  for (i = 0; i < od->uvepls.uvepls_included_vars; ++i) {
    od->mel.uvepls_included[i] = '1';
  }
  strncpy(&(od->mel.uvepls_included
    [od->uvepls.uvepls_included_vars]), "\n\0", 2);
  /*
  allocate a large variable matrix where all included
  variables, eventual dummy variables
  and all active + ext_pred objects will be put
  */
  od->mal.large_e_mat = double_mat_resize(od->mal.large_e_mat,
    od->active_object_num + od->ext_pred_object_num,
    od->uvepls.uvepls_included_vars * coeff_size);
  if (!(od->mal.large_e_mat)) {
    return OUT_OF_MEMORY;
  }
  if (!(od->uvepls.ive)) {
    /*
    if UVE-PLS was chosen, initialize random_seed
    */
    set_random_seed(od, od->random_seed);
  }
  else {
    /*
    allocate a vector to store q2 values for the models obtained
    by stepwise variable elimination
    */
    od->vel.predictivity_list = double_vec_resize(od->vel.predictivity_list,
      od->uvepls.uvepls_included_vars - pc_num + 1);
    if (!(od->vel.predictivity_list)) {
      return OUT_OF_MEMORY;
    }
    predictivity_list = od->vel.predictivity_list;
    od->pel.predictivity_list_rank = int_perm_resize(od->pel.predictivity_list_rank,
      od->uvepls.uvepls_included_vars - pc_num + 1);
    if (!(od->pel.predictivity_list_rank)) {
      return OUT_OF_MEMORY;
    }
    predictivity_list_rank = od->pel.predictivity_list_rank;
  }

  /*
  fill up both x and y matrices
  */
  result = fill_x_matrix(od, UVEPLS_FULL_MODEL,
    od->uvepls.use_srd_groups);
  if (result) {
    return result;
  }
  result = fill_y_matrix(od);
  if (result) {
    return result;
  }
      
  od->mal.press = double_mat_resize(od->mal.press,
    pc_num + 1, od->y_vars);
  if (!(od->mal.press)) {
    return OUT_OF_MEMORY;
  }
  memset(od->mal.press->base, 0, od->mal.press->m
    * od->mal.press->n * sizeof(double));
  result = alloc_pls(od, od->uvepls.uvepls_included_vars,
    pc_num, UVEPLS_CV_MODEL);
  if (result) {
    return OUT_OF_MEMORY;
  }
  result = alloc_cv_sdep(od, pc_num, od->uvepls.runs);
  if (result) {
    return OUT_OF_MEMORY;
  }
  result = prepare_cv(od, pc_num,  od->uvepls.cv_type,
    od->uvepls.groups, od->uvepls.runs);
  if (result) {
    return result;
  }
  if (od->uvepls.uve_m) {
    /*
    if UVE-M was chosen, allocate a vector and
    a permutation array to compute the median
    */
    od->vel.median_vec = double_vec_resize(od->vel.median_vec,
      od->cv.overall_cv_runs);
    if (!(od->vel.median_vec)) {
      return OUT_OF_MEMORY;
    }
    median_vec = od->vel.median_vec;
    od->pel.median_perm = int_perm_resize(od->pel.median_perm,
      od->cv.overall_cv_runs);
    if (!(od->pel.median_perm)) {
      return OUT_OF_MEMORY;
    }
    median_perm = od->pel.median_perm;
  }
  result = alloc_average_mat(od, UVEPLS_CV_MODEL,
    od->uvepls.cv_type, od->uvepls.groups,
    od->uvepls.runs);
  if (result) {
    return OUT_OF_MEMORY;
  }
  model_num = 0;
  /*
  if IVE-PLS was chosen, the minimum number of retained variables
  is equal to pc_num, so there will be (uvepls_included_vars - pc_num)
  variable eliminations
  */
  i = od->uvepls.uvepls_included_vars;
  ive_iterations = 0;
  while (i >= pc_num) {
    ++ive_iterations;
    if (!(od->uvepls.save_ram)) {
      od->mal.b_coefficients_store =
        double_mat_resize(od->mal.b_coefficients_store,
        od->cv.overall_cv_runs * od->y_vars,
        i * coeff_size);
      if (!(od->mal.b_coefficients_store)) {
        return OUT_OF_MEMORY;
      }
      b_coefficients_store = od->mal.b_coefficients_store;
      memset(b_coefficients_store->base, 0, b_coefficients_store->m
        * b_coefficients_store->n * sizeof(double));
    }
    /*
    allocate a matrix where average (or median, if UVE-M was chosen)
    PLS coefficients will be stored
    */
    od->mal.b_coefficients_ave = double_mat_resize
      (od->mal.b_coefficients_ave, i * coeff_size, od->y_vars);
    if (!(od->mal.b_coefficients_ave)) {
      return OUT_OF_MEMORY;
    }
    b_coefficients_ave = od->mal.b_coefficients_ave;
    memset(b_coefficients_ave->base, 0, b_coefficients_ave->m
      * b_coefficients_ave->n * sizeof(double));
    /*
    allocate a matrix where the standard deviation
    (or the interquartile range, if UVE-M was chosen)
    of PLS coefficients will be stored
    */
    od->mal.b_coefficients_sd = double_mat_resize
      (od->mal.b_coefficients_sd, i * coeff_size, od->y_vars);
    if (!(od->mal.b_coefficients_sd)) {
      return OUT_OF_MEMORY;
    }
    b_coefficients_sd = od->mal.b_coefficients_sd;
    memset(b_coefficients_sd->base, 0, b_coefficients_sd->m
      * b_coefficients_sd->n * sizeof(double));
    od->vel.eff_c = double_vec_resize(od->vel.eff_c, i);
    if (!(od->vel.eff_c)) {
      return OUT_OF_MEMORY;
    }
    eff_c = od->vel.eff_c;
    memset(od->vel.eff_c->ve, 0, od->vel.eff_c->size * sizeof(double));
    od->pel.eff_c_rank = int_perm_resize(od->pel.eff_c_rank, i);
    if (!(od->pel.eff_c_rank)) {
      return OUT_OF_MEMORY;
    }
    eff_c_rank = od->pel.eff_c_rank;
    if (od->uvepls.ive) {
      /*
      if IVE-PLS was chosen, then write the uvepls_included
      array in the TEMP_UVEPLS temporary file,
      including the final '\n' character
      */
      actual_len = fwrite(od->mel.uvepls_included, 1,
        od->uvepls.uvepls_included_vars + 1,
        od->file[TEMP_UVEPLS]->handle);
      if (actual_len !=
        (od->uvepls.uvepls_included_vars + 1)) {
        return CANNOT_WRITE_TEMP_FILE;
      }
    }
    else {
      od->vel.eff_dummy_c = double_vec_resize(od->vel.eff_dummy_c, i);
      if (!(od->vel.eff_dummy_c)) {
        return OUT_OF_MEMORY;
      }
      eff_dummy_c = od->vel.eff_dummy_c;
      memset(od->vel.eff_dummy_c->ve, 0, od->vel.eff_dummy_c->size * sizeof(double));
      od->pel.eff_dummy_c_rank = int_perm_resize(od->pel.eff_dummy_c_rank, i);
      if (!(od->pel.eff_dummy_c_rank)) {
        return OUT_OF_MEMORY;
      }
      eff_dummy_c_rank = od->pel.eff_dummy_c_rank;
    }      
    od->mal.e_mat = double_mat_resize(od->mal.e_mat,
      od->object_num, i * coeff_size);
    od->vel.e_mat_ave = double_vec_resize(od->vel.e_mat_ave, i * coeff_size);
    memset(od->mal.ave_press->base, 0,
      od->mal.ave_press->m * od->mal.ave_press->n * sizeof(double));
    if (od->n_proc > 1) {
      result = parallel_cv(od, i, pc_num,
        PARALLEL_CV | UVEPLS_CV_MODEL,
        od->uvepls.cv_type, od->uvepls.groups,
        od->uvepls.runs);
      if (result) {
        return result;
      }
    }
    else {
      result = cv(od, pc_num, UVEPLS_CV_MODEL,
        od->uvepls.cv_type, od->uvepls.groups,
        od->uvepls.runs);
      if (result) {
        return result;
      }
    }
    b_coefficients = od->mal.b_coefficients;
    /*
    the c vector allocated for PLS will be used
    to store the reliability criterion for included variables
    */
    od->vel.c = double_vec_resize(od->vel.c, i * coeff_size);
    if (!(od->vel.c)) {
      return OUT_OF_MEMORY;
    }
    c = od->vel.c;
    memset(c->ve, 0, i * sizeof(double));
    if (od->uvepls.uve_m) {
      /*
      if UVE-M was chosen, calculate the median of coefficients
      and store it into a matrix
      */
      for (n = 0; n < (i * coeff_size); ++n) {
        for (k = 0; k < od->y_vars; ++k) {
          if (od->uvepls.save_ram) {
            for (j = 0; j < od->cv.overall_cv_runs; ++j) {
              result = get_cv_coeff(od, j, k, n,
                &(median_vec->ve[j]), od->uvepls.save_ram);
              if (result) {
                return CANNOT_READ_TEMP_FILE;
              }
            }
          }
          else {
            cblas_dcopy(od->cv.overall_cv_runs,
              &M_PEEK(b_coefficients_store, k, n),
              od->y_vars,
              median_vec->ve, 1);
          }
          (void)double_vec_sort(median_vec, median_perm);
          /*
          if there is an even number of CV runs,
          then the median will be the average
          between the two central ones
          */
          if (od->cv.overall_cv_runs % 2) {
            /*
            if there is an odd number of CV runs,
            then the median will be the central one
            */
            M_POKE(b_coefficients_ave, n, k,
              median_vec->ve[od->cv.overall_cv_runs / 2]);
          }
          else {
            /*
            if there is an even number of CV runs,
            then the median will be the average
            between the two central ones
            */
            M_POKE(b_coefficients_ave, n, k,
              (median_vec->ve[od->cv.overall_cv_runs / 2 - 1]
              + median_vec->ve[od->cv.overall_cv_runs / 2])
              / 2.0);
          }
          first_half_end = od->cv.overall_cv_runs / 2;
          second_half_start = first_half_end
            + od->cv.overall_cv_runs % 2 + 1;
          if (first_half_end % 2) {
            lower_quartile =
              median_vec->ve[first_half_end / 2];
            higher_quartile =
              median_vec->ve[first_half_end / 2 + second_half_start - 1];
          }
          else {
            lower_quartile =
              (median_vec->ve[first_half_end / 2 - 1]
              + median_vec->ve[first_half_end / 2])
              / 2.0;
            higher_quartile =
              (median_vec->ve[first_half_end / 2 + second_half_start - 2]
              + median_vec->ve[first_half_end / 2 + second_half_start - 1])
              / 2.0;
          }
          M_POKE(b_coefficients_sd, n, k, higher_quartile - lower_quartile);
        }
      }
    }
    else {
      /*
      if standard UVE was chosen, then compute average
      PLS coefficient and their SD
      */
      for (j = 0; j < od->cv.overall_cv_runs; ++j) {
        if (od->uvepls.save_ram) {
          for (k = 0; k < od->y_vars; ++k) {
            for (n = 0; n < (i * coeff_size); ++n) {
              result = get_cv_coeff(od, j, k, n,
                &M_PEEK(b_coefficients, n, k),
                od->uvepls.save_ram);
              if (result) {
                return CANNOT_READ_TEMP_FILE;
              }
            }
          }
        }
        else {
          for (k = 0; k < od->y_vars; ++k) {
            cblas_dcopy(i * coeff_size,
              &M_PEEK(b_coefficients_store,
              j * od->y_vars + k, 0),
              b_coefficients_store->max_m,
              &M_PEEK(b_coefficients, 0, k), 1);
          }
        }
        cblas_daxpy(b_coefficients->m * b_coefficients->n,
          1.0, b_coefficients->base, 1,
          b_coefficients_ave->base, 1);
      }
      for (k = 0; k < od->y_vars; ++k) {
        for (j = 0; j < od->cv.overall_cv_runs; ++j) {
          if (od->uvepls.save_ram) {
            if (od->n_proc > 1) {
              n = 0;
              cv_run = j;
              while((cv_run - thread_info[n]->n_calc) >= 0) {
                cv_run -= thread_info[n]->n_calc;
                ++n;
              }
              file = thread_info[n]->temp_cv_coeff.handle;
            }
            else {
              file = od->file[TEMP_CV_COEFF]->handle;
            }
            fseek(file, (b_coefficients->m * (k + j * od->y_vars))
              * sizeof(double), SEEK_SET);
            actual_len = fread(&M_PEEK(b_coefficients, 0, k),
              sizeof(double), b_coefficients->m, file);
            if (actual_len != b_coefficients->m) {
              return CANNOT_READ_TEMP_FILE;
            }
          }
          else {
            cblas_dcopy(b_coefficients->m,
              &M_PEEK(b_coefficients_store,
              j * od->y_vars + k, 0),
              b_coefficients_store->max_m,
              &M_PEEK(b_coefficients, 0, k), 1);
          }
          if (!j) {
            cblas_dscal(b_coefficients_ave->m * b_coefficients_ave->n,
              1.0 / od->cv.overall_cv_runs,
              b_coefficients_ave->base, 1);
          }
          cblas_daxpy(b_coefficients->m * b_coefficients->n,
            -1.0, b_coefficients_ave->base, 1,
            b_coefficients->base, 1);
          for (n = 0; n < b_coefficients->m; ++n) {
            M_POKE(b_coefficients, n, k, 
              square(M_PEEK(b_coefficients, n, k)));
          }
          cblas_daxpy(b_coefficients->m * b_coefficients->n,
            1.0, b_coefficients->base, 1,
            b_coefficients_sd->base, 1);
        }
        for (n = 0; n < b_coefficients->m; ++n) {
          M_POKE(b_coefficients_sd, n, k,
            sqrt(M_PEEK(b_coefficients_sd, n, k)
            / (od->cv.overall_cv_runs - 1)));
        }
      }
    }
    /*
    compute the reliability coefficient c for real variables
    and eventually, if UVE-PLS was chosen, for dummy ones;
    if there are multiple y vars, then the overall reliability
    coefficient will be the sum of the coefficients computed
    for each y variable
    */
    for (k = 0; k < od->y_vars; ++k) {
      for (n = 0; n < b_coefficients->m; ++n) {
        /*
        c = sum(ave_b(n) / sd[b(n)]) over all y vars
        */
        c->ve[n] += (safe_rint(fabs(M_PEEK(b_coefficients_ave, n, k)
          / M_PEEK(b_coefficients_sd, n, k)) * 1.0e08) / 1.0e08);
      }
    }
    if (od->uvepls.use_srd_groups) {
      k = 0;
      c_count = 0;
      eff_c_count = 0;
      for (j = 0; j < od->voronoi_num; ++j) {
        included = 0;
        for (n = 0; n < od->mel.voronoi_active[j]; ++n) {
          if (od->mel.uvepls_included[k] == '1') {
            included = 1;
            eff_c->ve[eff_c_count] += c->ve[c_count];
            if (!(od->uvepls.ive)) {
              eff_dummy_c->ve[eff_c_count] +=
                c->ve[c_count + i];
            }
            ++c_count;
          }
          ++k;
        }
        if (included) {
          eff_c->ve[eff_c_count] /= od->mel.voronoi_active[j];
          if (!(od->uvepls.ive)) {
            eff_dummy_c->ve[eff_c_count] /=
              od->mel.voronoi_active[j];
          }
          ++eff_c_count;
        }
      }
      (void)double_vec_resize(eff_c, eff_c_count);
      if (!(od->uvepls.ive)) {
        (void)double_vec_resize(eff_dummy_c, eff_c_count);
      }
    }
    else {
      eff_c_count = i;
      for (n = 0; n < eff_c_count; ++n) {
        eff_c->ve[n] = c->ve[n];
        if (!(od->uvepls.ive)) {
          eff_dummy_c->ve[n] = c->ve[n + eff_c_count];
        }
      }
    }
    if (!(od->uvepls.ive)) {
      (void)double_vec_sort(eff_dummy_c, eff_dummy_c_rank);
      if ((int)(od->uvepls.uve_alpha)) {
        /*
        take the appropriate quantile
        as threshold value
        */
        nth_quantile = (int)safe_rint(od->uvepls.uve_alpha
          / (double)100 * (double)eff_c_count);
        threshold_c = eff_dummy_c->ve[nth_quantile];
      }
      else {
        threshold_c = eff_dummy_c->ve[eff_c_count - 1];
      }
      if (od->uvepls.use_srd_groups) {
        k = 0;
        eff_c_count = 0;
        for (j = 0; j < od->voronoi_num; ++j) {
          included = 0;
          for (n = 0; n < od->mel.voronoi_active[j]; ++n) {
            if (od->mel.uvepls_included[k] == '1') {
              included = 1;
              if (eff_c->ve[eff_c_count] < threshold_c) {
                od->mel.uvepls_included[k] = '0';
              }
            }
            ++k;
          }
          if (included) {
            ++eff_c_count;
          }
        }
      }
      else {
        for (n = 0; n < i; ++n) {
          if (eff_c->ve[n] < threshold_c) {
            od->mel.uvepls_included[n] = '0';
          }
        }
      }
      break;
    }
    else {
      /*
      store the q2 for the requested pc_num in the predictivity_list vector
      */
      if (od->pc_num == pc_num) {
        predictivity_list->ve[model_num] = od->vel.q2->ve[pc_num];
      }
      else {
        predictivity_list->ve[model_num] = -1.0e99;
      }
      /*
      if ive_external_sdep was requested, then store
      also the SDEP on the external test set, which will
      be used as a criterion to select the best model
      instead of q2
      */
      if (od->uvepls.ive_external_sdep) {
        trim_mean_center_matrix(od, od->mal.large_e_mat,
          &(od->mal.e_mat), &(od->vel.e_mat_ave),
          UVEPLS_FULL_MODEL, od->active_object_num);
        trim_mean_center_matrix(od, od->mal.large_f_mat,
          &(od->mal.f_mat), &(od->vel.f_mat_ave),
          UVEPLS_FULL_MODEL, od->active_object_num);
        pls(od, pc_num, UVEPLS_FULL_MODEL);
        result = pred_ext_y_values(od, od->pc_num,
          UVEPLS_FULL_MODEL);
        if (result) {
          return result;
        }
        if (od->pc_num == pc_num) {
          predictivity_list->ve[model_num] =
            od->vel.ave_sdep->ve[pc_num];
        }
        else {
          predictivity_list->ve[model_num] = 1.0e99;
        }
      }
      (void)double_vec_sort(eff_c, eff_c_rank);
      lowest_c_pos = eff_c_rank->pe[0];
      if (od->uvepls.use_srd_groups) {
        k = 0;
        eff_c_count = 0;
        for (j = 0; j < od->voronoi_num; ++j) {
          included = 0;
          eliminated = 0;
          for (n = 0; n < od->mel.voronoi_active[j]; ++n) {
            if (od->mel.uvepls_included[k] == '1') {
              included = 1;
              if (eff_c_count == lowest_c_pos) {
                eliminated = 1;
                od->mel.uvepls_included[k] = '0';
                --i;
              }
            }
            ++k;
          }
          if (eliminated) {
            break;
          }
          if (included) {
            ++eff_c_count;
          }
        }
      }
      else {
        k = 0;
        eff_c_count = 0;
        while (k < od->uvepls.uvepls_included_vars) {
          if (od->mel.uvepls_included[k] == '1') {
            if (eff_c_count == lowest_c_pos) {
              break;
            }
            ++eff_c_count;
          }
          ++k;
        }
        od->mel.uvepls_included[k] = '0';
        --i;
      }
      ++model_num;
    }
  }
  if (od->uvepls.ive) {
    (void)double_vec_resize(predictivity_list, ive_iterations);
    (void)int_perm_resize(predictivity_list_rank, ive_iterations);
    (void)double_vec_sort(predictivity_list, predictivity_list_rank);
    limit = (int)safe_rint
      ((double)(od->uvepls.ive_percent_limit)  / 100.0
      * (double)(od->uvepls.uvepls_included_vars));
    if (od->uvepls.ive_external_sdep) {
      n = 0;
      while (predictivity_list_rank->pe[n] > limit) {
        ++n;
      }
    }
    else {
      n = predictivity_list_rank->size - 1;
      while (predictivity_list_rank->pe[n] > limit) {
        --n;
      }
    }
    fseek(od->file[TEMP_UVEPLS]->handle,
      (od->uvepls.uvepls_included_vars + 1)
      * predictivity_list_rank->pe[n], SEEK_SET);
    actual_len = fread(od->mel.uvepls_included, 1,
      od->uvepls.uvepls_included_vars + 1,
      od->file[TEMP_UVEPLS]->handle);
    if (actual_len != (od->uvepls.uvepls_included_vars + 1)) {
      return CANNOT_READ_TEMP_FILE;
    }
  }
  /*
  print out the results of the variable selection
  */
  tee_printf(od,
    "%-5s%24s\n"
    "-----------------------------\n",
    "Field", "Excluded");
  n = 0;
  overall_n_excluded = 0;
  if (!(od->uvepls.use_srd_groups)) {
    for (i = 0; i < od->field_num; ++i) {
      if (get_field_attr(od, i, ACTIVE_BIT)) {
        n_excluded = 0;
        for (k = 0; k < od->x_vars; ++k) {
          if (get_x_var_attr(od, i, k, SEL_INCLUDED_BIT)) {
            if (od->mel.uvepls_included[n] == '0') {
              set_x_var_attr(od, i, k, UVEPLS_BIT, 1);
              ++n_excluded;
            }
            ++n;
          }
        }
        tee_printf(od, "%5d%24d\n", i + 1, n_excluded);
        overall_n_excluded += n_excluded;
      }
    }
  }
  else {
    overall_seed_count = 0;
    voronoi_num = 0;
    i = 0;
    while (i < od->field_num) {
      if (get_field_attr(od, i, ACTIVE_BIT)) {
        n_excluded = 0;
        while ((voronoi_num - overall_seed_count) < od->mel.seed_count[i]) {
          for (k = 0; k < od->x_vars; ++k) {
            if (get_x_var_attr(od, i, k, SEL_INCLUDED_BIT)
              && (get_voronoi_buf(od, i, k) == voronoi_num)) {
              if (od->mel.uvepls_included[n] == '0') {
                set_x_var_attr(od, i, k, UVEPLS_BIT, 1);
                ++n_excluded;
              }
              ++n;
            }
          }
          ++voronoi_num;
        }
        tee_printf(od, "%5d%24d\n", i + 1, n_excluded);
        overall_n_excluded += n_excluded;
      }
      overall_seed_count += od->mel.seed_count[i];
      ++i;
    }
  }
  tee_printf(od,
    "%-5s%24d\n"
    "-----------------------------\n\n",
    "Total", overall_n_excluded);
  if (od->n_proc > 1) {
    free_parallel_cv(od, thread_info, UVEPLS_CV_MODEL,
      od->uvepls.cv_type, od->uvepls.runs);
  }
  else if (od->uvepls.cv_type == LEAVE_MANY_OUT) {
    free_cv_groups(od, od->uvepls.runs);
  }
  if (od->mal.b_coefficients_store) {
    double_mat_free(od->mal.b_coefficients_store);
    od->mal.b_coefficients_store = NULL;
  }
  if (od->mal.b_coefficients_ave) {
    double_mat_free(od->mal.b_coefficients_ave);
    od->mal.b_coefficients_ave = NULL;
  }
  if (od->mal.b_coefficients_sd) {
    double_mat_free(od->mal.b_coefficients_sd);
    od->mal.b_coefficients_sd = NULL;
  }
  if (od->mal.large_e_mat_ave) {
    double_mat_free(od->mal.large_e_mat_ave);
    od->mal.large_e_mat_ave = NULL;
  }
  if (od->mal.large_f_mat_ave) {
    double_mat_free(od->mal.large_f_mat_ave);
    od->mal.large_f_mat_ave = NULL;
  }
  if (od->vel.c) {
    double_vec_free(od->vel.c);
    od->vel.c = NULL;
  }
  if (od->vel.eff_c) {
    double_vec_free(od->vel.eff_c);
    od->vel.eff_c = NULL;
  }
  if (od->pel.eff_c_rank) {
    int_perm_free(od->pel.eff_c_rank);
    od->pel.eff_c_rank = NULL;
  }
  
  return 0;
}
