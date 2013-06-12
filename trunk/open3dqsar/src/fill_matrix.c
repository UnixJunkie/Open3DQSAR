/*

fill_matrix.c

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


int fill_x_matrix(O3Data *od, int model_type, int use_srd_groups)
{
  int x_max_x = 0;
  int y_max;
  int i;
  int j;
  int k;
  int m;
  int n;
  int n_bit;
  int overall_seed_count;
  int x;
  int y;
  int actual_len;
  int result;
  unsigned short bit[2] = {ACTIVE_BIT, PREDICT_BIT};
  double value;
  double sumweight;
  
  
  if (model_type & (FFDSEL_FULL_MODEL | FFDSEL_CV_MODEL
    | UVEPLS_FULL_MODEL | UVEPLS_CV_MODEL)) {
    x_max_x = od->mal.large_e_mat->n;
  }
  else if (model_type & CV_MODEL) {
    x_max_x = od->overall_active_x_vars;
  }
  else {
    /*
    FULL_MODEL
    */
    x_max_x = od->overall_active_x_vars;
    if (open_temp_file(od, od->file[TEMP_X_MATRIX], "x_matrix")) {
      return CANNOT_WRITE_TEMP_FILE;
    }
  }
  y_max = od->active_object_num + od->ext_pred_object_num;
  /*
  allocate a y_max * x_max_x large-E matrix
  */
  od->mal.large_e_mat = double_mat_resize
    (od->mal.large_e_mat, y_max, x_max_x);
  if (!(od->mal.large_e_mat)) {
    return OUT_OF_MEMORY;
  }
  /*
  allocate a y_max * x_max_x E matrix
  */
  od->mal.e_mat = double_mat_resize(od->mal.e_mat, y_max, x_max_x);
  if (!(od->mal.e_mat)) {
    return OUT_OF_MEMORY;
  }
  /*
  allocate a vector for full averages
  */
  od->vel.e_mat_full_ave = double_vec_resize
    (od->vel.e_mat_full_ave, x_max_x);
  if (!(od->vel.e_mat_full_ave)) {
    return OUT_OF_MEMORY;
  }
  memset(od->vel.e_mat_full_ave->ve, 0,
    od->vel.e_mat_full_ave->size * sizeof(double));
  /*
  copy values from active fields, objects, x_vars
  mean-center them and store them into a matrix
  */
  
  x = 0;
  for (j = 0, sumweight = 0.0; j < od->object_num; ++j) {
    if (get_object_attr(od, j, ACTIVE_BIT)) {
      sumweight += od->mel.object_weight[j];
    }
  }
  if ((model_type & (FFDSEL_FULL_MODEL | UVEPLS_FULL_MODEL)) && use_srd_groups) {
    overall_seed_count = 0;
    for (i = 0; i < od->field_num; ++i) {
      if (get_field_attr(od, i, ACTIVE_BIT)) {
        for (m = 0; m < od->mel.seed_count[i]; ++m) {
          for (n = 0; n < od->mel.voronoi_fill[overall_seed_count + m]; ++n) {
            k = od->al.voronoi_composition[overall_seed_count + m][n];
            if (get_x_var_attr(od, i, k, ACTIVE_BIT)) {
              y = 0;
              for (n_bit = 0; n_bit < 2; ++n_bit) {
                for (j = 0; j < od->object_num; ++j) {
                  if (get_object_attr(od, j, bit[n_bit])) {
                    result = get_x_value(od, i, j, k, &value, 0);
                    if (result) {
                      return result;
                    }
                    if (!MISSING(value)) {
                      result = get_x_value(od, i, j, k,
                        &M_PEEK(od->mal.large_e_mat, y, x),
                        CUTOFF_BIT | WEIGHT_BIT);
                      if (result) {
                        return result;
                      }
                      if ((model_type & UVEPLS_FULL_MODEL)
                        && (!(od->uvepls.ive))) {
                        M_POKE(od->mal.large_e_mat,
                          y, x + x_max_x / 2,
                          (od->uvepls.dummy_range_coefficient
                          * (od->mel.x_data[i].max_x_value
                          - od->mel.x_data[i].min_x_value) * (genrand_real(od)
                          - 0.5) + 0.5 * (od->mel.x_data[i].max_x_value
                          + od->mel.x_data[i].min_x_value))
                          * od->mel.x_data[i].x_weight_coefficient
                          * od->uvepls.dummy_value_coefficient);
                      }
                      if (get_object_attr(od, j, ACTIVE_BIT)) {
                        od->vel.e_mat_full_ave->ve[x] +=
                          (M_PEEK(od->mal.large_e_mat, y, x)
                          * od->mel.object_weight[j]);
                        if ((model_type & UVEPLS_FULL_MODEL)
                          && (!(od->uvepls.ive))) {
                          od->vel.e_mat_full_ave->ve[x + x_max_x / 2] += 
                            (M_PEEK(od->mal.large_e_mat,
                            y, x + x_max_x / 2)
                            * od->mel.object_weight[j]);
                        }
                      }
                    }
                    else {
                      M_POKE(od->mal.large_e_mat, y, x, MISSING_VALUE);
                      if ((model_type & UVEPLS_FULL_MODEL)
                        && (!(od->uvepls.ive))) {
                        M_POKE(od->mal.large_e_mat,
                          y, x + x_max_x / 2, MISSING_VALUE);
                      }
                    }
                    ++y;
                  }
                }
              }
              if (sumweight > 0.0) {
                od->vel.e_mat_full_ave->ve[x] /= sumweight;
                if ((model_type & UVEPLS_FULL_MODEL)
                  && (!(od->uvepls.ive))) {
                  od->vel.e_mat_full_ave->ve[x + x_max_x / 2] /= sumweight;
                }
              }
              ++x;
            }
          }
        }
      }
      overall_seed_count += od->mel.seed_count[i];
    }
  }
  else {
    for (i = 0; i < od->field_num; ++i) {
      if (get_field_attr(od, i, ACTIVE_BIT)) {
        for (k = 0; k < od->x_vars; ++k) {
          if (get_x_var_attr(od, i, k, ACTIVE_BIT)) {
            y = 0;
            for (n_bit = 0; n_bit < 2; ++n_bit) {
              for (j = 0; j < od->object_num; ++j) {
                if (get_object_attr(od, j, bit[n_bit])) {
                  result = get_x_value(od, i, j, k, &value, 0);
                  if (result) {
                    return result;
                  }
                  if (!MISSING(value)) {
                    result = get_x_value(od, i, j, k,
                      &M_PEEK(od->mal.large_e_mat, y, x),
                      CUTOFF_BIT | WEIGHT_BIT);
                    if (result) {
                      return result;
                    }
                    if ((model_type & UVEPLS_FULL_MODEL)
                      && (!(od->uvepls.ive))) {
                      M_POKE(od->mal.large_e_mat, y, x + x_max_x / 2,
                        (od->uvepls.dummy_range_coefficient
                        * (od->mel.x_data[i].max_x_value
                        - od->mel.x_data[i].min_x_value) * (genrand_real(od)
                        - 0.5) + 0.5 * (od->mel.x_data[i].max_x_value
                        + od->mel.x_data[i].min_x_value))
                        * od->mel.x_data[i].x_weight_coefficient
                        * od->uvepls.dummy_value_coefficient);
                    }
                    if (get_object_attr(od, j, ACTIVE_BIT)) {
                      od->vel.e_mat_full_ave->ve[x] +=
                        (M_PEEK(od->mal.large_e_mat, y, x)
                        * od->mel.object_weight[j]);
                      if ((model_type & UVEPLS_FULL_MODEL)
                        && (!(od->uvepls.ive))) {
                        od->vel.e_mat_full_ave->ve[x + x_max_x / 2] += 
                          (M_PEEK(od->mal.large_e_mat,
                          y, x + x_max_x / 2)
                          * od->mel.object_weight[j]);
                      }
                    }
                  }
                  else {
                    M_POKE(od->mal.large_e_mat, y, x, MISSING_VALUE);
                    if ((model_type & UVEPLS_FULL_MODEL)
                      && (!(od->uvepls.ive))) {
                      M_POKE(od->mal.large_e_mat,
                        y, x + x_max_x / 2, MISSING_VALUE);
                    }
                  }
                  ++y;
                }
              }
            }
            if (sumweight > 0.0) {
              od->vel.e_mat_full_ave->ve[x] /= sumweight;
              if ((model_type & UVEPLS_FULL_MODEL)
                && (!(od->uvepls.ive))) {
                od->vel.e_mat_full_ave->ve[x + x_max_x / 2] /= sumweight;
              }
            }
            if (model_type & FULL_MODEL) {
              fwrite(&i, sizeof(int), 1,
                od->file[TEMP_X_MATRIX]->handle);
              actual_len = fwrite(&k, sizeof(int), 1,
                od->file[TEMP_X_MATRIX]->handle);
              if (actual_len != 1) {
                return CANNOT_WRITE_TEMP_FILE;
              }
            }
            ++x;
          }
        }
      }
    }
  }
  if (od->file[TEMP_X_MATRIX]->handle) {
    fclose(od->file[TEMP_X_MATRIX]->handle);
    od->file[TEMP_X_MATRIX]->handle = NULL;
  }

  return 0;
}


int fill_y_matrix(O3Data *od)
{
  int j;
  int k;
  int y;
  int y_max;
  int n_bit;
  unsigned short bit[2] = {ACTIVE_BIT, PREDICT_BIT};
  double sumweight;
  
  
  y_max = od->active_object_num + od->ext_pred_object_num;
  /*
  allocate a y_max * od->y_vars large-F matrix
  */
  od->mal.large_f_mat = double_mat_resize
    (od->mal.large_f_mat, y_max, od->y_vars);
  if (!(od->mal.large_f_mat)) {
    return OUT_OF_MEMORY;
  }
  /*
  allocate a y_max * od->y_vars F matrix
  */
  od->mal.f_mat = double_mat_resize
    (od->mal.f_mat, y_max, od->y_vars);
  if (!(od->mal.f_mat)) {
    return OUT_OF_MEMORY;
  }
  /*
  allocate a vector for full averages
  */
  od->vel.f_mat_full_ave = double_vec_resize
    (od->vel.f_mat_full_ave, od->y_vars);
  if (!(od->vel.f_mat_full_ave)) {
    return OUT_OF_MEMORY;
  }
  memset(od->vel.f_mat_full_ave->ve, 0,
    od->vel.f_mat_full_ave->size * sizeof(double));
  od->vel.active_value_ave = double_vec_resize
    (od->vel.active_value_ave, od->y_vars);
  if (!(od->vel.active_value_ave)) {
    return OUT_OF_MEMORY;
  }
  memset(od->vel.active_value_ave->ve, 0,
    od->vel.active_value_ave->size * sizeof(double));
  /*
  copy values from active objects, y_vars
  mean-center them and store them into a matrix
  */
  for (j = 0, sumweight = 0.0; j < od->object_num; ++j) {
    if (get_object_attr(od, j, ACTIVE_BIT)) {
      sumweight += od->mel.object_weight[j];
    }
  }
  for (k = 0; k < od->y_vars; ++k) {
    y = 0;
    for (n_bit = 0; n_bit < 2; ++n_bit) {
      for (j = 0; j < od->object_num; ++j) {
        if (get_object_attr(od, j, bit[n_bit])) {
          M_POKE(od->mal.large_f_mat, y, k,
            get_y_value(od, j, k, WEIGHT_BIT));
          if (get_object_attr(od, j, ACTIVE_BIT)) {
            od->vel.f_mat_full_ave->ve[k] +=
              M_PEEK(od->mal.large_f_mat, y, k)
              * od->mel.object_weight[j];
          }
          ++y;
        }
      }
    }
    if (sumweight > 0.0) {
      od->vel.f_mat_full_ave->ve[k] /= sumweight;
    }
  }
  
  return 0;
}


int fill_x_matrix_scrambled(O3Data *od)
{
  int x_max_x = 0;
  int y_max;
  int i;
  int j;
  int k;
  int n;
  int offset;
  int x = 0;
  int y;
  int object_num;
  int struct_num;
  int conf_num;
  int n_conf;
  int active_struct_num;
  int active_struct_count;
  int actual_len;
  int found;
  int result;
  double value;
  double sumweight;
  
  
  get_attr_struct_ave(od, 0, ACTIVE_BIT, &active_struct_num, NULL);
  x_max_x = od->overall_active_x_vars;
  y_max = od->active_object_num + od->ext_pred_object_num;
  /*
  allocate a y_max * x_max_x large-E matrix
  */
  od->mal.large_e_mat = double_mat_resize
    (od->mal.large_e_mat, y_max, x_max_x);
  if (!(od->mal.large_e_mat)) {
    return OUT_OF_MEMORY;
  }
  memset(od->mal.large_e_mat->base, 0,
    od->mal.large_e_mat->m * od->mal.large_e_mat->n
    * sizeof(double));
  /*
  allocate a y_max * x_max_x E matrix
  */
  od->mal.e_mat = double_mat_resize(od->mal.e_mat, y_max, x_max_x);
  if (!(od->mal.e_mat)) {
    return OUT_OF_MEMORY;
  }
  /*
  allocate a vector for full averages
  */
  od->vel.e_mat_full_ave =
    double_vec_resize(od->vel.e_mat_full_ave, x_max_x);
  if (!(od->vel.e_mat_full_ave)) {
    return OUT_OF_MEMORY;
  }
  memset(od->vel.e_mat_full_ave->ve, 0,
    od->vel.e_mat_full_ave->size * sizeof(double));
  /*
  read the permutation order of objects
  when they are sorted according to decreasing
  y value (average y value if there are multiple y's)
  */
  actual_len = fread(od->pel.scrambling_order->pe,
    sizeof(int), active_struct_num,
    od->file[TEMP_SCRAMBLE]->handle);
  if (actual_len != active_struct_num) {
    return CANNOT_READ_TEMP_FILE;
  }
  /*
  copy values from active fields, objects, x_vars
  mean-center them and store them into a matrix
  */
  for (j = 0, sumweight = 0.0; j < od->object_num; ++j) {
    if (get_object_attr(od, j, ACTIVE_BIT)) {
      sumweight += od->mel.object_weight[j];
    }
  }
  for (i = 0, offset = 0; i < od->field_num; ++i) {
    if (get_field_attr(od, i, ACTIVE_BIT)) {
      for (n = 0, j = 0; n < od->pel.scrambling_order->size; ++n) {
        active_struct_count = 0;
        object_num = 0;
        found = 0;
        while ((!found) && (object_num < od->object_num)) {
          struct_num = od->al.mol_info[object_num]->struct_num;
          conf_num = ((od->valid & COSMOTHERM_BIT)
            ? od->al.cosmo_list[struct_num]->n_conf[BOUND] : 1);
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
              for (k = 0, x = 0; k < od->x_vars; ++k) {
                if (get_x_var_attr(od, i, k, ACTIVE_BIT)) {
                  result = get_x_value(od, i, object_num, k, &value, 0);
                  if (result) {
                    return result;
                  }
                  if (!MISSING(value)) {
                    result = get_x_value(od, i, object_num, k,
                      &M_PEEK(od->mal.large_e_mat, j, x + offset),
                      CUTOFF_BIT | WEIGHT_BIT);
                    if (result) {
                      return result;
                    }
                    od->vel.e_mat_full_ave->ve[x + offset] +=
                      (M_PEEK(od->mal.large_e_mat, j, x + offset)
                      * od->mel.object_weight[object_num]);
}
                  else {
                    M_POKE(od->mal.large_e_mat, j, x + offset, MISSING_VALUE);
                  }
                  ++x;
                }
              }
              ++j;
            }
            ++n_conf;
            ++object_num;
          }
        }
      }
      offset += od->mel.x_data[i].active_x_vars;
    }
  }
  if (sumweight > 0.0) {
    for (x = 0; x < od->vel.e_mat_full_ave->size; ++x) {
      od->vel.e_mat_full_ave->ve[x] /= sumweight;
    }
  }
  /*
  tee_printf(od, "LARGE_E_MAT\n");
  for (i = 0; i < 4; ++i) {
    for (j = 0; j < od->mal.large_e_mat->n; ++j) {
      tee_printf(od, "%12.4lf", M_PEEK(od->mal.large_e_mat, i, j));
    }
    tee_printf(od, "\n");
  }
  */

  return 0;
}


int fill_y_matrix_scrambled(O3Data *od)
{
  int j;
  int k;
  int n;
  int y;
  int object_num;
  int struct_num;
  int conf_num;
  int n_conf;
  int active_struct_count;
  int y_max;
  int actual_len;
  int found;
  double sumweight;
  
  
  y_max = od->active_object_num + od->ext_pred_object_num;
  /*
  allocate a y_max * od->y_vars large-F matrix
  */
  od->mal.large_f_mat = double_mat_resize
    (od->mal.large_f_mat, y_max, od->y_vars);
  if (!(od->mal.large_f_mat)) {
    return OUT_OF_MEMORY;
  }
  memset(od->mal.large_f_mat->base, 0,
    od->mal.large_f_mat->m * od->mal.large_f_mat->n
    * sizeof(double));
  /*
  allocate a y_max * od->y_vars F matrix
  */
  od->mal.f_mat = double_mat_resize
    (od->mal.f_mat, y_max, od->y_vars);
  if (!(od->mal.f_mat)) {
    return OUT_OF_MEMORY;
  }
  /*
  allocate a vector for full averages
  */
  od->vel.f_mat_full_ave = double_vec_resize
    (od->vel.f_mat_full_ave, od->y_vars);
  if (!(od->vel.f_mat_full_ave)) {
    return OUT_OF_MEMORY;
  }
  memset(od->vel.f_mat_full_ave->ve, 0,
    od->vel.f_mat_full_ave->size * sizeof(double));
  actual_len = fread(od->pel.scrambling_temp->pe, sizeof(int),
    od->active_object_num, od->file[TEMP_SCRAMBLE]->handle);
  if (actual_len != od->active_object_num) {
    return CANNOT_READ_TEMP_FILE;
  }
  /*
  copy values from active objects, y_vars
  and store them in the scrambled order in
  large_f_mat
  */
  for (j = 0, sumweight = 0.0; j < od->object_num; ++j) {
    if (get_object_attr(od, j, ACTIVE_BIT)) {
      sumweight += od->mel.object_weight[j];
    }
  }
  for (k = 0; k < od->y_vars; ++k) {
    for (n = 0, j = 0; n < od->pel.scrambling_temp->size; ++n) {
      active_struct_count = 0;
      object_num = 0;
      found = 0;
      while ((!found) && (object_num < od->object_num)) {
        struct_num = od->al.mol_info[object_num]->struct_num;
        conf_num = ((od->valid & COSMOTHERM_BIT)
          ? od->al.cosmo_list[struct_num]->n_conf[BOUND] : 1);
        for (n_conf = 0, found = 0, y = 0; n_conf < conf_num; ++n_conf, ++object_num) {
          if (get_object_attr(od, object_num, ACTIVE_BIT)) {
            if (active_struct_count == od->pel.scrambling_temp->pe[n]) {
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
            M_POKE(od->mal.large_f_mat, j, k,
              get_y_value(od, object_num, k, WEIGHT_BIT));
            od->vel.f_mat_full_ave->ve[k] +=
              M_PEEK(od->mal.large_f_mat, j, k)
              * od->mel.object_weight[object_num];
            ++j;
          }
          ++n_conf;
          ++object_num;
        }
      }
    }
    od->vel.f_mat_full_ave->ve[k] /= sumweight;
  }
  /*
  tee_printf(od, "LARGE_F_MAT\n");
  for (k = 0; k < 4; ++k) {
    for (j = 0; j < od->mal.large_f_mat->n; ++j) {
      tee_printf(od, "%12.4lf", M_PEEK(od->mal.large_f_mat, k, j));
    }
    tee_printf(od, "\n");
  }
  */
  
  return 0;
}


int fill_x_matrix_pca(O3Data *od)
{
  int x_max_x = 0;
  int y_max;
  int i;
  int j;
  int k;
  int n_active;
  int x;
  int y;
  int actual_len;
  int result;
  double value;
  
  
  /*
  FULL_MODEL
  */
  x_max_x = od->overall_active_x_vars;
  if (open_temp_file(od, od->file[TEMP_X_MATRIX], "x_matrix")) {
    return CANNOT_WRITE_TEMP_FILE;
  }
  y_max = od->active_object_num
    + od->ext_pred_object_num;
  /*
  allocate a y_max * x_max_x large-E matrix
  */
  od->mal.large_e_mat =
    double_mat_resize(od->mal.large_e_mat, y_max, x_max_x);
  if (!(od->mal.large_e_mat)) {
    return OUT_OF_MEMORY;
  }
  /*
  allocate a y_max * x_max_x E matrix
  */
  od->mal.e_mat =
    double_mat_resize(od->mal.e_mat, y_max, x_max_x);
  if (!(od->mal.e_mat)) {
    return OUT_OF_MEMORY;
  }
  /*
  allocate a vector for full averages
  */
  od->vel.e_mat_full_ave =
    double_vec_resize(od->vel.e_mat_full_ave,
    x_max_x);
  if (!(od->vel.e_mat_full_ave)) {
    return OUT_OF_MEMORY;
  }
  memset(od->vel.e_mat_full_ave->ve, 0,
    od->vel.e_mat_full_ave->size
    * sizeof(double));
  /*
  copy values from active fields, objects, x_vars
  mean-center them and store them into a matrix
  */
  
  x = 0;
  for (i = 0; i < od->field_num; ++i) {
    if (get_field_attr(od, i, ACTIVE_BIT)) {
      for (k = 0; k < od->x_vars; ++k) {
        if (get_x_var_attr(od, i, k, ACTIVE_BIT)) {
          y = 0;
          for (j = 0, n_active = 0; j < od->object_num; ++j) {
            if (get_object_attr(od, j, ACTIVE_BIT | PREDICT_BIT)) {
              result = get_x_value(od, i, j, k, &value, 0);
              if (result) {
                return result;
              }
              if (!MISSING(value)) {
                result = get_x_value(od, i, j, k,
                  &M_PEEK(od->mal.large_e_mat, y, x),
                  CUTOFF_BIT | WEIGHT_BIT);
                if (result) {
                  return result;
                }
                if (get_object_attr(od, j, ACTIVE_BIT)) {
                  od->vel.e_mat_full_ave->ve[x] +=
                    M_PEEK(od->mal.large_e_mat, y, x);
                  ++n_active;
                }
              }
              else {
                M_POKE(od->mal.large_e_mat, y, x, MISSING_VALUE);
              }
              ++y;
            }
          }
          if (n_active) {
            od->vel.e_mat_full_ave->ve[x] /= (double)n_active;
          }
          fwrite(&i, sizeof(int), 1,
            od->file[TEMP_X_MATRIX]->handle);
          actual_len = fwrite(&k, sizeof(int), 1,
            od->file[TEMP_X_MATRIX]->handle);
          if (actual_len != 1) {
            return CANNOT_WRITE_TEMP_FILE;
          }
          ++x;
        }
      }
    }
  }
  if (od->file[TEMP_X_MATRIX]->handle) {
    fclose(od->file[TEMP_X_MATRIX]->handle);
    od->file[TEMP_X_MATRIX]->handle = NULL;
  }

  return 0;
}
