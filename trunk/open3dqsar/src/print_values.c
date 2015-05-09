/*

print_values.c

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


int mol_to_sdf(O3Data *od, int object_num, double actual_value)
{
  char buffer[BUF_LEN];
  int found;
  FileDescriptor mol_fd;
  
  
  memset(buffer, 0, BUF_LEN);
  if (od->file[ASCII_IN]->handle) {
    memset(&mol_fd, 0, sizeof(FileDescriptor));
    sprintf(mol_fd.name, "%s%c%04d.mol", od->field.mol_dir,
      SEPARATOR, od->al.mol_info[object_num]->object_id);
    if (!(mol_fd.handle = fopen(mol_fd.name, "rb"))) {
      return CANNOT_READ_TEMP_FILE;
    }
    found = 0;
    while (fgets(buffer, BUF_LEN, mol_fd.handle) && (!found)) {
      buffer[BUF_LEN - 1] = '\0';
      fprintf(od->file[ASCII_IN]->handle, "%s", buffer);
      found = (!strncmp(buffer, MOL_DELIMITER, strlen(MOL_DELIMITER)));
    }
    if (!found) {
      return CANNOT_READ_TEMP_FILE;
    }
    fclose(mol_fd.handle);
    fprintf(od->file[ASCII_IN]->handle,
      ">  <N>\n%d\n\n"
      ">  <ID>\n%d\n\n"
      ">  <Name>\n%s\n\n"
      ">  <Actual>\n%.4lf\n\n",
      object_num + 1, od->al.mol_info[object_num]->object_id,
      od->al.mol_info[object_num]->object_name, actual_value);
  }
  
  return 0;
}


int print_calc_values(O3Data *od, int options)
{
  char format[MAX_NAME_LEN];
  char buffer[BUF_LEN];
  int i;
  int j;
  int k;
  int pos1;
  int pos2;
  int x;
  int y;
  int n;
  int object_num;
  int struct_num;
  int conf_num;
  int n_conf;
  int n_hat;
  int x_max_x;
  int x_max_y;
  int y_max;
  int result;
  int pc_num;
  int opt_pc_n;
  int actual_len;
  int active_struct_num = 0;
  int active_found;
  int predict_found;
  double leverage;
  double sumweight;
  double best_delta;
  double value;
  double actual_value;
  double active_value_ave;
  double tss;
  
  
  /*
  calculated y values are stored in this way:
  pc_num
  total number of variables (equal to x_max_y, that is y_vars)
  total number of blocks in this file

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
  memset(format, 0, MAX_NAME_LEN);
  memset(buffer, 0, BUF_LEN);
  opt_pc_n = 0;
  /*
  Read number of PCs
  */
  actual_len = fread(&pc_num, sizeof(int), 1,
    od->file[TEMP_CALC]->handle);
  if (actual_len != 1) {
    return CANNOT_READ_TEMP_FILE;
  }
  if (!(od->mel.weighted_value = (double *)realloc
    (od->mel.weighted_value, (pc_num + 1) * sizeof(double)))) {
    return OUT_OF_MEMORY;
  }
  /*
  Read total number of variables, that is y_vars
  */
  actual_len = fread(&x_max_y, sizeof(int), 1,
    od->file[TEMP_CALC]->handle);
  if (actual_len != 1) {
    return CANNOT_READ_TEMP_FILE;
  }
  /*
  Read number of objects in the first (and only) block
  */
  actual_len = fread(&y_max, sizeof(int), 1,
    od->file[TEMP_CALC]->handle);
  if (actual_len != 1) {
    return CANNOT_READ_TEMP_FILE;
  }
  /*
  Skip object numbers: no need to read them, we already
  know which they are
  */
  if (fseek(od->file[TEMP_CALC]->handle, y_max * sizeof(int), SEEK_CUR)) {
    return CANNOT_READ_TEMP_FILE;
  }
  
  pos1 = ftell(od->file[TEMP_CALC]->handle);
  tss = 0.0;
  memset(od->vel.r2->ve, 0, (pc_num + 1) * sizeof(double));
  memset(od->vel.ave_sdec->ve, 0, (pc_num + 1)  * sizeof(double));
  for (x = 0; x < x_max_y; ++x) {
    /*
    move to the x-th variable
    */
    if (fseek(od->file[TEMP_CALC]->handle,
      x * y_max * sizeof(double), SEEK_CUR)) {
      return CANNOT_READ_TEMP_FILE;
    }  
    pos2 = ftell(od->file[TEMP_CALC]->handle);
    tee_printf(od, "Calculated values for dependent variable %2d (%s)\n",
      x + 1, od->cimal.y_var_name->me[x]);
    tee_printf(od, "-------------------------------------------------------------------");
    for (i = 0; i < (pc_num + 1 + ((options & CALC_LEVERAGE_BIT) ? 1 : 0)); ++i) {
      tee_printf(od, "------------");
    }
    tee_printf(od, "\n%5s%5s%5s    %-36s%12s", "N", "ID", "Str", "Name", "Actual");
    for (i = 0; i < pc_num; ++i) {
      tee_printf(od, "%12d", i + 1);
    }
    tee_printf(od, "%12s", "Opt PC n");
    if (options & CALC_LEVERAGE_BIT) {
      tee_printf(od, "%12s", "Leverage");
    }
    tee_printf(od, "\n-------------------------------------------------------------------");
    for (i = 0; i < (pc_num + 1 + ((options & CALC_LEVERAGE_BIT) ? 1 : 0)); ++i) {
      tee_printf(od, "------------");
    }
    object_num = 0;
    struct_num = 0;
    i = 0;
    y = 0;
    n_hat = 0;
    active_struct_num = 0;
    active_value_ave = 0.0;
    actual_value = 0.0;
    while (object_num < od->object_num) {
      struct_num = od->al.mol_info[object_num]->struct_num;
      conf_num = 1;
      for (n_conf = 0, active_found = 0, predict_found = 0,
        leverage = 0.0; n_conf < conf_num; ++n_conf) {
        if (get_object_attr(od, object_num + n_conf, ACTIVE_BIT)) {
          if (!active_found) {
            actual_value = get_y_value(od, object_num + n_conf, x, WEIGHT_BIT);
            active_value_ave += actual_value;
            ++active_struct_num;
          }
          if (options & CALC_LEVERAGE_BIT) {
            leverage += M_PEEK(od->mal.hat_mat, n_hat, n_hat);
          }
          ++n_hat;
          ++active_found;
        }
      }
      leverage /= (double)conf_num;
      if (!active_found) {
        object_num += n_conf;
        continue;
      }
      sprintf(buffer, "%lf", actual_value);
      strcpy(format, ((strlen(buffer) > 11) ? "%12.4le" : "%12.4lf"));
      tee_printf(od, "\n%5d%5d%5d    %-36s", object_num + 1,
        od->al.mol_info[object_num]->object_id, struct_num + 1,
        od->al.mol_info[object_num]->object_name);
      tee_printf(od, format, actual_value);
      if (mol_to_sdf(od, object_num, actual_value)) {
        return CANNOT_READ_TEMP_FILE;
      }
      best_delta = 1.0e99;
      memset(od->mel.weighted_value, 0, (pc_num + 1) * sizeof(double));
      for (n_conf = 0, sumweight = 0.0; n_conf < conf_num; ++n_conf, ++object_num) {
        if (!get_object_attr(od, object_num, ACTIVE_BIT | PREDICT_BIT)) {
          continue;
        }
        sumweight += od->mel.object_weight[object_num];
        for (j = 0; j <= pc_num; ++j) {
          /*
          move to the j-th component
          */
          if (fseek(od->file[TEMP_CALC]->handle,
            (y_max * x_max_y * (pc_num - j) + y)
            * sizeof(double), SEEK_CUR)) {
            return CANNOT_READ_TEMP_FILE;
          }
          /*
          Read predicted value for y-th variable, j-th component
          */
          actual_len = fread(&value, sizeof(double), 1,
            od->file[TEMP_CALC]->handle);
          if (actual_len != 1) {
            return CANNOT_READ_TEMP_FILE;
          }
          od->mel.weighted_value[j] += (value
            * od->mel.object_weight[object_num]);
          /*
          go back to starting position
          */
          if (fseek(od->file[TEMP_CALC]->handle, pos2, SEEK_SET)) {
            return CANNOT_READ_TEMP_FILE;
          }
        }
        ++y;
        ++i;
      }
      for (j = 0; j <= pc_num; ++j) {
        od->mel.weighted_value[j] /= sumweight;
        if (j) {
          if (fabs(actual_value - od->mel.weighted_value[j]) < best_delta) {
            best_delta = fabs(actual_value - od->mel.weighted_value[j]);
            opt_pc_n = j;
          }
          if (sumweight > 0.0) {
            sprintf(buffer, "%lf", od->mel.weighted_value[j]);
            strcpy(format, ((strlen(buffer) > 11) ? "%12.4le" : "%12.4lf"));
            tee_printf(od, format, od->mel.weighted_value[j]);
            if (od->file[ASCII_IN]->handle) {
              fprintf(od->file[ASCII_IN]->handle, ">  <PC %d>\n%.4lf\n\n",
                j, od->mel.weighted_value[j]);
            }
          }
          else {
            tee_printf(od, "%12s", "-");
            if (od->file[ASCII_IN]->handle) {
              fprintf(od->file[ASCII_IN]->handle,
                ">  <PC %d>\n%s\n\n", j, "-");
            }
          }
        }
        od->vel.ave_sdec->ve[j] += square
          (actual_value - od->mel.weighted_value[j]);
      }
      tee_printf(od, "%12d", opt_pc_n);
      if (options & CALC_LEVERAGE_BIT) {
        tee_printf(od, "%12.4lf", leverage);
      }
      if (od->file[ASCII_IN]->handle) {
        fprintf(od->file[ASCII_IN]->handle,
          ">  <Opt PC n>\n%d\n\n"
          "$$$$\n", opt_pc_n);
      }
    }
    active_value_ave /= (double)active_struct_num;
    object_num = 0;
    while (object_num < od->object_num) {
      struct_num = od->al.mol_info[object_num]->struct_num;
      conf_num = 1;
      for (n_conf = 0, active_found = 0; n_conf < conf_num; ++n_conf, ++object_num) {
        if (get_object_attr(od, object_num, ACTIVE_BIT)) {
          if (!active_found) {
            actual_value = get_y_value(od, object_num, x, WEIGHT_BIT);
            tss += square(actual_value - active_value_ave);
          }
          ++active_found;
        }
      }
    }
    if (fseek(od->file[TEMP_CALC]->handle, pos1, SEEK_SET)) {
      return CANNOT_READ_TEMP_FILE;
    }
  }
  tee_printf(od, "\n\n\nPC%12s%12s%12s\n", "SDEC", "r2", "F-test");
  tee_printf(od, "--------------------------------------\n");
  for (i = 0; i <= pc_num; ++i) {
    od->vel.r2->ve[i] = 1.0 - od->vel.ave_sdec->ve[i] / tss;
    od->vel.ave_sdec->ve[i] = sqrt
      (od->vel.ave_sdec->ve[i] / (double)active_struct_num);
    if (i) {
      sprintf(buffer, "%12.4lf", (double)(od->active_object_num - i - 1)
        * od->vel.r2->ve[i] / ((double)i * (1.0 - od->vel.r2->ve[i])));
    }
    else {
      sprintf(buffer, "%12s", "-");
    }
    tee_printf(od, "%2d%12.4lf%12.4lf%12s\n", i,
      od->vel.ave_sdec->ve[i], od->vel.r2->ve[i], buffer);
  }
  tee_printf(od, "\n");
  if (options & CALC_FIELD_CONTRIB_BIT) {
    tee_printf(od, "\nRelative field contributions\n--");
    for (i = 0; i < od->field_num; ++i) {
      if (get_field_attr(od, i, ACTIVE_BIT)) {
        tee_printf(od, "------------");
      }
    }
    tee_printf(od, "\nPC");
    for (i = 0; i < od->field_num; ++i) {
      if (get_field_attr(od, i, ACTIVE_BIT)) {
        sprintf(buffer, "x%d", i + 1);
        tee_printf(od, "%12s", buffer);
      }
    }
    tee_printf(od, "\n--");
    for (i = 0; i < od->field_num; ++i) {
      if (get_field_attr(od, i, ACTIVE_BIT)) {
        tee_printf(od, "------------");
      }
    }
    for (i = 0; i <= pc_num; ++i) {
      tee_printf(od, "\n%2d", i);
      result = reload_coefficients(od, i);
      if (result) {
        return result;
      }
      memset(od->mel.field_contrib, 0, sizeof(double) * od->field_num);
      value = 0.0;
      n = 0;
      for (j = 0; j < od->field_num; ++j) {
        if (get_field_attr(od, j, ACTIVE_BIT)) {
          for (x = 0; x < od->x_vars; ++x) {
            if (get_x_var_attr(od, j, x, ACTIVE_BIT)) {
              for (k = 0; k < od->y_vars; ++k) {
                od->mel.field_contrib[j] += (fabs(M_PEEK(od->mal.b_coefficients, n, k))
                  * get_x_var_buf(od, j, x, STDDEV_BUF));
              }
              ++n;
            }
          }
        }
        value += od->mel.field_contrib[j];
      }
      for (j = 0; j < od->field_num; ++j) {
        if (get_field_attr(od, j, ACTIVE_BIT)) {
          tee_printf(od, "%12.4lf", ((value > 0.0)
            ? od->mel.field_contrib[j] / value : 0.0));
        }
      }
    }
    tee_printf(od, "\n\n\n");
  }
  if (od->mal.hat_mat) {
    double_mat_free(od->mal.hat_mat);
    od->mal.hat_mat = NULL;
  }
  
  return 0;
}


int print_ext_pred_values(O3Data *od)
{
  char format[MAX_NAME_LEN];
  char buffer[BUF_LEN];
  int i;
  int j;
  int pos1;
  int pos2;
  int x;
  int y;
  int object_num;
  int struct_num;
  int conf_num;
  int n_conf;
  int x_max_y;
  int y_max;
  int pc_num;
  int opt_pc_n;
  int actual_len;
  int active_struct_num = 0;
  int predict_struct_num = 0;
  int active_found = 0;
  int predict_found = 0;
  double sumweight;
  double best_delta;
  double value;
  double actual_value;
  double active_value_ave;
  double cum_press;
  double tss;
  
  
  /*
  calculated y values are stored in this way:
  pc_num
  total number of variables (equal to x_max_y, that is y_vars)
  total number of blocks in this file

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
  memset(format, 0, MAX_NAME_LEN);
  memset(buffer, 0, BUF_LEN);
  opt_pc_n = 0;
  /*
  Read number of PCs
  */
  actual_len = fread(&pc_num, sizeof(int), 1,
    od->file[TEMP_EXT_PRED]->handle);
  if (actual_len != 1) {
    return CANNOT_READ_TEMP_FILE;
  }
  if (!(od->mel.weighted_value = (double *)realloc
    (od->mel.weighted_value, (pc_num + 1) * sizeof(double)))) {
    return OUT_OF_MEMORY;
  }
  /*
  Read total number of variables, that is y_vars
  */
  actual_len = fread(&x_max_y, sizeof(int), 1,
    od->file[TEMP_EXT_PRED]->handle);
  if (actual_len != 1) {
    return CANNOT_READ_TEMP_FILE;
  }
  /*
  Read number of objects in the first (and only) block
  */
  actual_len = fread(&y_max, sizeof(int), 1,
    od->file[TEMP_EXT_PRED]->handle);
  if (actual_len != 1) {
    return CANNOT_READ_TEMP_FILE;
  }
  /*
  Skip object numbers: no need to read them, we already
  know which they are
  */
  if (fseek(od->file[TEMP_EXT_PRED]->handle, y_max * sizeof(int), SEEK_CUR)) {
    return CANNOT_READ_TEMP_FILE;
  }
  
  pos1 = ftell(od->file[TEMP_EXT_PRED]->handle);
  tss = 0.0;
  for (x = 0; x < x_max_y; ++x) {
    /*
    move to the x-th variable
    */
    if (fseek(od->file[TEMP_EXT_PRED]->handle,
      x * y_max * sizeof(double), SEEK_CUR)) {
      return CANNOT_READ_TEMP_FILE;
    }  
    pos2 = ftell(od->file[TEMP_EXT_PRED]->handle);
    tee_printf(od, "External predictions for dependent variable %2d (%s)\n",
      x + 1, od->cimal.y_var_name->me[x]);
    tee_printf(od, "-------------------------------------------------------------------");
    for (i = 0; i < (pc_num + 1); ++i) {
      tee_printf(od, "------------");
    }
    tee_printf(od, "\n%5s%5s%5s    %-36s%12s", "N", "ID", "Str", "Name", "Actual");
    for (i = 0; i < pc_num; ++i) {
      tee_printf(od, "%12d", i + 1);
    }
    tee_printf(od, "%12s", "Opt PC n");
    tee_printf(od, "\n-------------------------------------------------------------------");
    for (i = 0; i < (pc_num + 1); ++i) {
      tee_printf(od, "------------");
    }
    object_num = 0;
    struct_num = 0;
    y = 0;
    active_struct_num = 0;
    predict_struct_num = 0;
    active_value_ave = 0.0;
    actual_value = 0.0;
    while (object_num < od->object_num) {
      struct_num = od->al.mol_info[object_num]->struct_num;
      conf_num = 1;
      for (n_conf = 0, active_found = 0, predict_found = 0; n_conf < conf_num; ++n_conf) {
        if (get_object_attr(od, object_num + n_conf, ACTIVE_BIT)) {
          if (!active_found) {
            actual_value = get_y_value(od, object_num + n_conf, x, WEIGHT_BIT);
            active_value_ave += actual_value;
            ++active_struct_num;
          }
          ++active_found;
        }
        if (get_object_attr(od, object_num + n_conf, PREDICT_BIT)) {
          if (!predict_found) {
            actual_value = get_y_value(od, object_num + n_conf, x, WEIGHT_BIT);
            ++predict_struct_num;
          }
          ++predict_found;
        }
      }
      if (!predict_found) {
        object_num += conf_num;
        i += active_found;
        continue;
      }
      sprintf(buffer, "%lf", actual_value);
      strcpy(format, ((strlen(buffer) > 11) ? "%12.4le" : "%12.4lf"));
      tee_printf(od, "\n%5d%5d%5d    %-36s", object_num + 1,
        od->al.mol_info[object_num]->object_id, struct_num + 1,
        od->al.mol_info[object_num]->object_name);
      tee_printf(od, format, actual_value);
      if (mol_to_sdf(od, object_num, actual_value)) {
        return CANNOT_READ_TEMP_FILE;
      }
      best_delta = 1.0e99;
      memset(od->mel.weighted_value, 0, (pc_num + 1) * sizeof(double));
      for (n_conf = 0, sumweight = 0.0; n_conf < conf_num; ++n_conf, ++object_num) {
        if (!get_object_attr(od, object_num, ACTIVE_BIT | PREDICT_BIT)) {
          continue;
        }
        sumweight += od->mel.object_weight[object_num];
        for (j = 0; j <= pc_num; ++j) {
          /*
          move to the j-th component
          */
          if (fseek(od->file[TEMP_EXT_PRED]->handle,
            (y_max * x_max_y * (pc_num - j) + y)
            * sizeof(double), SEEK_CUR)) {
            return CANNOT_READ_TEMP_FILE;
          }
          /*
          Read predicted value for y-th variable, j-th component
          */
          actual_len = fread(&value, sizeof(double), 1,
            od->file[TEMP_EXT_PRED]->handle);
          if (actual_len != 1) {
            return CANNOT_READ_TEMP_FILE;
          }
          od->mel.weighted_value[j] += (value
            * od->mel.object_weight[object_num]);
          /*
          go back to starting position
          */
          if (fseek(od->file[TEMP_EXT_PRED]->handle, pos2, SEEK_SET)) {
            return CANNOT_READ_TEMP_FILE;
          }
        }
        ++y;
      }
      for (j = 0; j <= pc_num; ++j) {
        od->mel.weighted_value[j] /= sumweight;
        if (j) {
          if (fabs(actual_value - od->mel.weighted_value[j]) < best_delta) {
            best_delta = fabs(actual_value - od->mel.weighted_value[j]);
            opt_pc_n = j;
          }
          if (sumweight > 0.0) {
            sprintf(buffer, "%lf", od->mel.weighted_value[j]);
            strcpy(format, ((strlen(buffer) > 11) ? "%12.4le" : "%12.4lf"));
            tee_printf(od, format, od->mel.weighted_value[j]);
            if (od->file[ASCII_IN]->handle) {
              fprintf(od->file[ASCII_IN]->handle, ">  <PC %d>\n%.4lf\n\n",
                j, od->mel.weighted_value[j]);
            }
          }
          else {
            tee_printf(od, "%12s", "-");
            if (od->file[ASCII_IN]->handle) {
              fprintf(od->file[ASCII_IN]->handle,
                ">  <PC %d>\n%s\n\n", j, "-");
            }
          }
        }
      }
      tee_printf(od, "%12d", opt_pc_n);
      if (od->file[ASCII_IN]->handle) {
        fprintf(od->file[ASCII_IN]->handle,
          ">  <Opt PC n>\n%d\n\n"
          "$$$$\n", opt_pc_n);
      }
    }
    active_value_ave /= (double)active_struct_num;
    object_num = 0;
    while (object_num < od->object_num) {
      struct_num = od->al.mol_info[object_num]->struct_num;
      conf_num = 1;
      for (n_conf = 0, predict_found = 0; n_conf < conf_num; ++n_conf, ++object_num) {
        if (!get_object_attr(od, object_num, ACTIVE_BIT | PREDICT_BIT)) {
          continue;
        }
        if (get_object_attr(od, object_num, PREDICT_BIT)) {
          if (!predict_found) {
            actual_value = get_y_value(od, object_num, x, WEIGHT_BIT);
            tss += square(actual_value - active_value_ave);
          }
          ++predict_found;
        }
      }
    }
    if (fseek(od->file[TEMP_EXT_PRED]->handle, pos1, SEEK_SET)) {
      return CANNOT_READ_TEMP_FILE;
    }
  }
  tee_printf(od, "\n\n\nPC%12s%12s\n", "SDEP", "r2(pred)");
  tee_printf(od, "--------------------------\n");
  for (i = 0; i <= pc_num; ++i) {
    for (x = 0, cum_press = 0.0; x < od->y_vars; ++x) {
      cum_press += M_PEEK(od->mal.press, i, x);
    }
    od->vel.r2_pred->ve[i] = 1.0 - cum_press / tss;
    od->vel.ave_sdep->ve[i] = sqrt(cum_press / (double)predict_struct_num);
    tee_printf(od, "%2d%12.4lf%12.4lf\n", i,
      od->vel.ave_sdep->ve[i], od->vel.r2_pred->ve[i]);
  }
  tee_printf(od, "\n");
  
  return 0;
}


int print_pred_values(O3Data *od)
{
  char format[MAX_NAME_LEN];
  char buffer[BUF_LEN];
  int i;
  int j;
  int n;
  int object_num;
  int struct_num;
  int conf_num;
  int n_conf;
  int pos1;
  int x;
  int y_max;
  int y_vars;
  int pc_num;
  int opt_pc_n;
  int actual_len;
  int num_values;
  int object_is_in_list;
  int active_found = 0;
  int predict_found = 0;
  double sumweight;
  double value;
  double best_delta;
  double actual_value;
  
  
  /*
  predicted y values are stored in this way:
  pc_num
  total number of variables (equal to od->y_vars, that is y_vars)
  total number of blocks in this file

  *** BLOCK 1 ***
  number of objects in this block (cmp_blk)
  list of object numbers in this block (length of the list: cmp_blk)
  var-0_value-0_pc-pc_num, var-0_value-1_pc-pc_num, var-0_value-2_pc-pc_num, ..., var-0_value-cmp_blk-pc_num
  var-1_value-0_pc-pc_num, var-1_value-1_pc-pc_num, var-1_value-2_pc-pc_num, ..., var-1_value-cmp_blk-pc_num
  var-2_value-0_pc-pc_num, var-2_value-1_pc-pc_num, var-2_value-2_pc-pc_num, ..., var-2_value-cmp_blk-pc_num
  ...
  var-od->y_vars_value-0_pc-pc_num, var-od->y_vars_value-1_pc-pc_num, var-od->y_vars_value-2_pc-pc_num, ..., var-od->y_vars_value-cmp_blk-pc_num

  var-0_value-0_pc-(pc_num-1), var-0_value-1_pc-(pc_num-1), var-0_value-2_pc-(pc_num-1), ..., var-0_value-cmp_blk-(pc_num-1)
  var-1_value-0_pc-(pc_num-1), var-1_value-1_pc-(pc_num-1), var-1_value-2_pc-(pc_num-1), ..., var-1_value-cmp_blk-(pc_num-1)
  var-2_value-0_pc-(pc_num-1), var-2_value-1_pc-(pc_num-1), var-2_value-2_pc-(pc_num-1), ..., var-2_value-cmp_blk-(pc_num-1)
  ...
  var-od->y_vars_value-0_pc-(pc_num-1), var-od->y_vars_value-1_pc-(pc_num-1), var-od->y_vars_value-2_pc-(pc_num-1), ..., var-od->y_vars_value-cmp_blk-(pc_num-1)

  var-0_value-0_pc-(pc_num-2), var-0_value-1_pc-(pc_num-2), var-0_value-2_pc-(pc_num-2), ..., var-0_value-cmp_blk-(pc_num-2)
  var-1_value-0_pc-(pc_num-2), var-1_value-1_pc-(pc_num-2), var-1_value-2_pc-(pc_num-2), ..., var-1_value-cmp_blk-(pc_num-2)
  var-2_value-0_pc-(pc_num-2), var-2_value-1_pc-(pc_num-2), var-2_value-2_pc-(pc_num-2), ..., var-2_value-cmp_blk-(pc_num-2)
  ...
  var-od->y_vars_value-0_pc-(pc_num-2), var-od->y_vars_value-1_pc-(pc_num-2), var-od->y_vars_value-2_pc-(pc_num-2), ..., var-od->y_vars_value-cmp_blk-(pc_num-2)

  ...
  ...
  ...

  var-0_value-0_pc-(pc-0), var-0_value-1_pc-(pc-0), var-0_value-2_pc-(pc-0), ..., var-0_value-cmp_blk-(pc-0)
  var-1_value-0_pc-(pc-0), var-1_value-1_pc-(pc-0), var-1_value-2_pc-(pc-0), ..., var-1_value-cmp_blk-(pc-0)
  var-2_value-0_pc-(pc-0), var-2_value-1_pc-(pc-0), var-2_value-2_pc-(pc-0), ..., var-2_value-cmp_blk-(pc-0)
  ...
  var-od->y_vars_value-0_pc-(pc-0), var-od->y_vars_value-1_pc-(pc-0), var-od->y_vars_value-2_pc-(pc-0), ..., var-od->y_vars_value-cmp_blk-(pc-0)

  *** BLOCK 2 ***

  ...
  ...
  ...

  */
  memset(format, 0, MAX_NAME_LEN);
  memset(buffer, 0, BUF_LEN);
  opt_pc_n = 0;
  od->mal.pred_f_mat = double_mat_resize
    (od->mal.pred_f_mat, od->active_object_num, od->y_vars);
  od->mel.predicted_object_list = (int *)realloc
    (od->mel.predicted_object_list, sizeof(int) * od->object_num);
  if (!(od->mel.predicted_object_list)) {
    return OUT_OF_MEMORY;
  }
  od->mel.sum = (double *)realloc(od->mel.sum, (od->pc_num + 1) * sizeof(double));
  if (!(od->mel.sum)) {
    return OUT_OF_MEMORY;
  }
  if (!(od->mel.weighted_value = (double *)realloc
    (od->mel.weighted_value, (od->pc_num + 1) * sizeof(double)))) {
    return OUT_OF_MEMORY;
  }
  for (x = 0; x < od->y_vars; ++x) {
    tee_printf(od, "Predicted values for dependent variable %2d (%s)\n",
      x + 1, od->cimal.y_var_name->me[x]);
    tee_printf(od, "-------------------------------------------------------------------");
    for (i = 0; i < (od->pc_num + 1); ++i) {
      tee_printf(od, "------------");
    }
    tee_printf(od, "\n%5s%5s%5s    %-36s%12s", "N", "ID", "Str", "Name", "Actual");
    for (i = 0; i < od->pc_num; ++i) {
      tee_printf(od, "%12d", i + 1);
    }
    tee_printf(od, "%12s", "Opt PC n");
    tee_printf(od, "\n-------------------------------------------------------------------");
    for (i = 0; i < (od->pc_num + 1); ++i) {
      tee_printf(od, "------------");
    }
    object_num = 0;
    struct_num = 0;
    actual_value = 0.0;
    while (object_num < od->object_num) {
      struct_num = od->al.mol_info[object_num]->struct_num;
      conf_num = 1;
      for (n_conf = 0, active_found = 0, predict_found = 0; n_conf < conf_num; ++n_conf) {
        if (get_object_attr(od, object_num + n_conf, ACTIVE_BIT)) {
          actual_value = get_y_value(od, object_num + n_conf, x, WEIGHT_BIT);
          ++active_found;
        }
        if (get_object_attr(od, object_num + n_conf, PREDICT_BIT)) {
          ++predict_found;
        }
      }
      if (!active_found) {
        object_num += n_conf;
        continue;
      }
      sprintf(buffer, "%lf", actual_value);
      strcpy(format, ((strlen(buffer) > 11) ? "%12.4le" : "%12.4lf"));
      tee_printf(od, "\n%5d%5d%5d    %-36s", object_num + 1,
        od->al.mol_info[object_num]->object_id, struct_num + 1,
        od->al.mol_info[object_num]->object_name);
      tee_printf(od, format, actual_value);
      if (mol_to_sdf(od, object_num, actual_value)) {
        return CANNOT_READ_TEMP_FILE;
      }
      memset(od->mel.weighted_value, 0, (od->pc_num + 1) * sizeof(double));
      for (n_conf = 0, sumweight = 0.0; n_conf < conf_num; ++n_conf, ++object_num) {
        if (!get_object_attr(od, object_num, ACTIVE_BIT | PREDICT_BIT)) {
          continue;
        }
        sumweight += od->mel.object_weight[object_num];
        num_values = 0;
        memset(od->mel.sum, 0, (od->pc_num + 1) * sizeof(double));
        rewind(od->file[TEMP_PRED]->handle);
        while (!feof(od->file[TEMP_PRED]->handle)) {
          /*
          Read number of PCs
          */
          actual_len = fread(&pc_num, sizeof(int), 1,
            od->file[TEMP_PRED]->handle);
          if (actual_len != 1) {
            if (!feof(od->file[TEMP_PRED]->handle)) {
              return CANNOT_READ_TEMP_FILE;
            }
            else {
              break;
            }
          }
          /*
          Read total number of variables, that is y_vars
          */
          actual_len = fread(&y_vars, sizeof(int), 1,
            od->file[TEMP_PRED]->handle);
          if ((actual_len != 1) || (y_vars != od->y_vars)) {
            return CANNOT_READ_TEMP_FILE;
          }
          /*
          Read number of objects in the next block
          */
          actual_len = fread(&y_max, sizeof(int), 1,
            od->file[TEMP_PRED]->handle);
          if (actual_len != 1) {
            return CANNOT_READ_TEMP_FILE;
          }
          /*
          Read object numbers for this block, then check
          if the object for which we want to write
          predicted values is present in this block
          */
          object_is_in_list = -1;
          for (n = 0; n < y_max; ++n) {
            actual_len = fread(&(od->mel.predicted_object_list[n]),
              sizeof(int), 1, od->file[TEMP_PRED]->handle);
            if (actual_len != 1) {
              return CANNOT_READ_TEMP_FILE;
            }
            /*
            if the object is present, remember
            its position in the block
            */
            if (od->mel.predicted_object_list[n] == object_num) {
              object_is_in_list = n;
            }
          }
          /*
          if so, read the value for the current var and add it
          */
          if (object_is_in_list != -1) {
            /*
            remember the starting position
            for the current block
            */
            ++num_values;
            pos1 = ftell(od->file[TEMP_PRED]->handle);
            for (j = 0; j <= pc_num; ++j) {
              /*
              move to the j-th component, x-th var,
              object_is_in_list-th value
              */
              if (fseek(od->file[TEMP_PRED]->handle,
                (y_max * od->y_vars * (pc_num - j) + y_max * x +
                object_is_in_list) * sizeof(double), SEEK_CUR)) {
                return CANNOT_READ_TEMP_FILE;
              }
              /*
              Read predicted value for x-th variable,
              j-th component, object_is_in_list-th value
              */
              actual_len = fread(&value, sizeof(double), 1,
                od->file[TEMP_PRED]->handle);
              if (actual_len != 1) {
                return CANNOT_READ_TEMP_FILE;
              }
              od->mel.sum[j] += value;
              /*
              go back to starting position for this block
              */
              if (fseek(od->file[TEMP_PRED]->handle, pos1, SEEK_SET)) {
                return CANNOT_READ_TEMP_FILE;
              }
            }
          }
          /*
          skip to the next block (if there is one)
          */
          if (fseek(od->file[TEMP_PRED]->handle,
            y_max * od->y_vars * (pc_num + 1) * sizeof(double), SEEK_CUR)) {
            return CANNOT_READ_TEMP_FILE;
          }
        }
        for (j = 0; j <= pc_num; ++j) {
          od->mel.weighted_value[j] +=
            (od->mel.sum[j] / (double)num_values
            * od->mel.object_weight[object_num]);
        }
      }
      best_delta = 1.0e99;
      for (j = 0; j <= pc_num; ++j) {
        od->mel.weighted_value[j] /= sumweight;
        if (j) {
          if (fabs(actual_value - od->mel.weighted_value[j]) < best_delta) {
            best_delta = fabs(actual_value - od->mel.weighted_value[j]);
            opt_pc_n = j;
          }
          if (sumweight > 0.0) {
            sprintf(buffer, "%lf", od->mel.weighted_value[j]);
            strcpy(format, ((strlen(buffer) > 11) ? "%12.4le" : "%12.4lf"));
            tee_printf(od, format, od->mel.weighted_value[j]);
            if (od->file[ASCII_IN]->handle) {
              fprintf(od->file[ASCII_IN]->handle,
                ">  <PC %d>\n%.4lf\n\n", j, od->mel.weighted_value[j]);
            }
          }
          else {
            tee_printf(od, "%12s", "-");
            if (od->file[ASCII_IN]->handle) {
              fprintf(od->file[ASCII_IN]->handle,
                ">  <PC %d>\n%s\n\n", j, "-");
            }
          }
        }
      }
      tee_printf(od, "%12d", opt_pc_n);
      if (od->file[ASCII_IN]->handle) {
        fprintf(od->file[ASCII_IN]->handle,
          ">  <Opt PC n>\n%d\n\n"
          "$$$$\n", opt_pc_n);
      }
    }
    tee_printf(od, "\n\n");
  }
  
  return 0;
}


void print_pls_scores(O3Data *od, int options)
{
  char format[MAX_NAME_LEN];
  char buffer[BUF_LEN];
  int bit[2] = { X_SCORES, Y_SCORES };
  int i;
  int j;
  int k;
  int y;
  int object_num;
  int struct_num;
  int info;
  #if (!defined HAVE_LIBLAPACK_ATLAS) && (!defined HAVE_LIBSUNPERF)
  int lwork;
  #endif
  double actual_value;
  
  
  memset(format, 0, MAX_NAME_LEN);
  memset(buffer, 0, BUF_LEN);
  if (options & PREDICT_BIT) {
    double_mat_resize(od->mal.x_weights, od->mal.e_mat->n, od->pc_num);
    double_mat_resize(od->mal.x_loadings, od->mal.e_mat->n, od->pc_num);
    double_mat_resize(od->mal.temp, od->pc_num, od->pc_num);
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
      lwork = (od->pc_num + 1) * LWORK_BLOCK_SIZE * sizeof(double);
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
      lwork = (od->pc_num + 1) * LWORK_BLOCK_SIZE * sizeof(double);
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
        od->pc_num, 1.0,
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
    double_mat_resize(od->mal.x_scores, od->mal.e_mat->m, od->pc_num);
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
      od->mal.e_mat->m,
      od->mal.x_weights_star->n,
      od->mal.e_mat->n, 1.0,
      od->mal.e_mat->base,
      od->mal.e_mat->max_m,
      od->mal.x_weights_star->base,
      od->mal.x_weights_star->max_m, 0.0,
      od->mal.x_scores->base,
      od->mal.x_scores->max_m);
  }
  for (k = 0; k <= 1; ++k) {
    if (!(options & bit[k])) {
      continue;
    }
    tee_printf(od, "\n"
      "%s %c scores\n"
      "-------------------------------------------------------",
      (options & PREDICT_BIT) ? "Predicted" : "PLS", k ? 'Y' : 'X');
    for (i = 0; i < od->pc_num; ++i) {
      tee_printf(od, "------------");
    }
    tee_printf(od, "\n%5s%5s%5s    %-36s", "N", "ID", "Str", "Name");
    for (i = 0; i < od->pc_num; ++i) {
      tee_printf(od, "%12d", i + 1);
    }
    tee_printf(od, "\n"
      "-------------------------------------------------------");
    for (i = 0; i < od->pc_num ; ++i) {
      tee_printf(od, "------------");
    }
    object_num = 0;
    struct_num = 0;
    for (object_num = 0, y = 0; object_num < od->object_num; ++object_num) {
      if (!get_object_attr(od, object_num,
        (options & PREDICT_BIT) ? PREDICT_BIT : ACTIVE_BIT)) {
        continue;
      }
      struct_num = od->al.mol_info[object_num]->struct_num;
      tee_printf(od, "\n%5d%5d%5d    %-36s", object_num + 1,
        od->al.mol_info[object_num]->object_id, struct_num + 1,
        od->al.mol_info[object_num]->object_name);
      for (j = 0; j < od->pc_num; ++j) {
        actual_value = M_PEEK(k ? od->mal.y_scores : od->mal.x_scores, y, j);
        sprintf(buffer, "%lf", actual_value);
        strcpy(format, ((strlen(buffer) > 11) ? "%12.4le" : "%12.4lf"));
        tee_printf(od, format, actual_value);
      }
      ++y;
    }
    tee_printf(od, "\n"
      "-------------------------------------------------------");
    for (i = 0; i < od->pc_num ; ++i) {
      tee_printf(od, "------------");
    }
    tee_printf(od, "\n\n");
  }
}
