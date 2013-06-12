/*

import_cosmotherm.c

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
#ifdef WIN32
#include <windows.h>
#endif


int import_cosmotherm(O3Data *od, char regex_name[MAX_COSMOTHERM_STATES][BUF_LEN], int *state)
{
  char buffer[BUF_LEN];
  char dirname[BUF_LEN];
  char *base_regex;
  char *nextfile;
  char *eof = NULL;
  int j = 0;
  int n;
  int n_mol = 0;
  int n_mol_copy = 0;
  int n_wat = 0;
  int n_solute;
  int n_solute_copy;
  int object_num = 0;
  int struct_num = 0;
  int conf_num = 0;
  int found = 0;
  double mw;
  #ifndef WIN32
  DIR *dir = NULL;
  struct dirent dir_entry;
  struct dirent *result;
  #else
  HANDLE dir;
  WIN32_FIND_DATA filedata;
  #endif


  memset(buffer, 0, BUF_LEN);
  memset(dirname, 0, BUF_LEN);
  for (*state = 0; *state < MAX_COSMOTHERM_STATES; ++(*state)) {
    free_array(od->al.regex_list[*state]);
    od->al.regex_list[*state] = NULL;
    od->cosmotherm_object_num[*state] = 0;
    if (!regex_name[*state][0]) {
      continue;
    }
    strcpy(dirname, regex_name[*state]);
    get_dirname(dirname);
    base_regex = get_basename(regex_name[*state]);
    for (j = 0; j < 2; ++j) {
      /*
      in the first run (j = 0) objects are counted to know how much memory
      should be allocated, in the second run (j = 1) data are stored
      */
      #ifndef WIN32
      dir = opendir(dirname);
      if (!dir) {
      #else
      sprintf(buffer, "%s\\*", dirname);
      dir = FindFirstFileA(buffer, &filedata);
      if (dir == INVALID_HANDLE_VALUE) {
      #endif
        return CANNOT_OPEN_DIRECTORY;
      }
      od->cosmotherm_object_num[*state] = 0;
      #ifndef WIN32
      while (!readdir_r(dir, &dir_entry, &result)) {
        if (!result) {
          break;
        }
        nextfile = dir_entry.d_name;
      #else
      while (FindNextFileA(dir, &filedata)) {
        nextfile = filedata.cFileName;
      #endif
        n = sscanf(nextfile, base_regex, &struct_num, &conf_num);
        /*
        if the current dir_entry matches base_regex, add it to the list
        */
        if (n == 2) {
          sprintf(buffer, base_regex, struct_num, conf_num);
          #ifndef WIN32
          if (!strcmp(buffer, nextfile)) {
          #else
          if (!strcasecmp(buffer, nextfile)) {
          #endif
            /*
            for *state = 0, there must be exactly object_num files matching regex_name
            */
            if ((!(*state)) && (od->cosmotherm_object_num[*state] == od->grid.object_num)) {
              return TOO_MANY_OBJECTS;
            }
            if (j) {
              od->al.regex_list[*state][od->cosmotherm_object_num[*state]]->struct_num = struct_num - 1;
              od->al.regex_list[*state][od->cosmotherm_object_num[*state]]->conf_num = conf_num - 1;
            }
            ++(od->cosmotherm_object_num[*state]);
          }
        }
      }
      #ifndef WIN32
      closedir(dir);
      #else
      CloseHandle(dir);
      #endif
      if ((!(*state)) && (od->cosmotherm_object_num[*state] < od->grid.object_num)) {
        return NOT_ENOUGH_OBJECTS;
      }
      if ((!j) && od->cosmotherm_object_num[*state]) {
        if (!(od->al.regex_list[*state] = (RegexData **)alloc_array
          (od->cosmotherm_object_num[*state], sizeof(RegexData)))) {
          return OUT_OF_MEMORY;
        }
      }
    }
    qsort(od->al.regex_list[*state], od->cosmotherm_object_num[*state], sizeof(RegexData *), compare_regex_data);
    for (object_num = 0; object_num < od->cosmotherm_object_num[*state]; ++object_num) {
      if (!(*state)) {
        /*
        this is only needed for conformations in the bound state
        */
        od->al.mol_info[object_num]->struct_num = od->al.regex_list[*state][object_num]->struct_num;
        od->al.mol_info[object_num]->conf_num = od->al.regex_list[*state][object_num]->conf_num;
      }
      sprintf(od->file[ASCII_IN]->name, regex_name[*state],
        od->al.regex_list[*state][object_num]->struct_num + 1,
        od->al.regex_list[*state][object_num]->conf_num + 1);
      if (!(od->file[ASCII_IN]->handle = fopen(od->file[ASCII_IN]->name, "rb"))) {
        return CANNOT_OPEN_COSMOTHERM_FILE;
      }
      /*
      in the mixture we should not have more than 2 compounds
      (i.e., solute and water)
      */
      if (!fgrep(od->file[ASCII_IN]->handle, buffer,
        "Total number of processed compounds")) {
        return PREMATURE_EOF;
      }
      sscanf(buffer, "%*s %*s %*s %*s %*s %d", &n_mol);
      if ((n_mol != 2)) {
        return PREMATURE_EOF;
      }
      rewind(od->file[ASCII_IN]->handle);
      /*
      find out which molecule number we are interested in
      */
      n_wat = 0;
      n = 0;
      eof = 0;
      while ((!n_wat) && (n < n_mol) && (eof = fgets
        (buffer, BUF_LEN, od->file[ASCII_IN]->handle))) {
        buffer[BUF_LEN - 1] = '\0';
        if (strstr(buffer, "Compound Information for compound")) {
          ++n;
          if (strstr(buffer, "h2o")) {
            sscanf(buffer, "%*s %*s %*s %*s %d", &n_wat);
            while ((!found) && (eof = fgets
              (buffer, BUF_LEN, od->file[ASCII_IN]->handle))) {
              if (strchr("\n\r", (int)buffer[0])) {
                break;
              }
              if ((found = (strstr(buffer, "Molecular Weight") ? 1 : 0))) {
                sscanf(buffer, "%*s %*s %*s %lf", &mw);
                if ((int)(safe_rint(mw * 1.0e03)) == 18015) {
                  n_wat = n;
                }
              }
            }
            if ((!eof) || (!found)) {
              return PREMATURE_EOF;
            }
          }
        }
      }
      if ((!eof) || (!n_wat)) {
        return PREMATURE_EOF;
      }
      n_solute = 3 - n_wat;
      if (!fgrep(od->file[ASCII_IN]->handle, buffer, "Results for mixture")) {
        return PREMATURE_EOF;
      }
      n_mol_copy = n_mol;
      n_solute_copy = n_solute;
      while (n_mol_copy && n_solute_copy && (eof = fgets
        (buffer, BUF_LEN, od->file[ASCII_IN]->handle))) {
        if (strstr(buffer, "Free energy of molecule in mix")) {
          --n_mol_copy;
          --n_solute_copy;
          if (!n_solute_copy) {
            found = 1;
            sscanf(buffer, "%*s %*s %*s %*s %*s %*s %*s %*s %lf",
              &(od->al.regex_list[*state][object_num]->g));
            /*
            divide COSMOtherm g by RT (units: kcal, K, mol)
            */
            od->al.regex_list[*state][object_num]->g /= (R_KCAL_K_MOL * 298.15);
          }
        }
      }
      if (!eof) {
        return PREMATURE_EOF;
      }
      /*
      for bound conformations also Kow is needed
      */
      if (!(*state)) {
        rewind(od->file[ASCII_IN]->handle);
        n_mol_copy = n_mol;
        n_solute_copy = n_solute;
        while (n_mol_copy && n_solute_copy && (eof = fgets
          (buffer, BUF_LEN, od->file[ASCII_IN]->handle))) {
          if (strstr(buffer, "QSPR property log(Pow)")) {
            --n_mol_copy;
            --n_solute_copy;
            if (!n_solute_copy) {
              found = 1;
              sscanf(buffer, "%*s %*s %*s %*s %lf",
                &(od->al.mol_info[object_num]->ln_k));
              /*
              convert log(Kow) to ln(Kow)
              */
              od->al.mol_info[object_num]->ln_k /= M_LOG10E;
              /*
              printf("%04d_%06d\t%lf\t%lf\n",
                od->al.mol_info[object_num]->struct_num + 1,
                od->al.mol_info[object_num]->conf_num + 1,
                od->al.mol_info[object_num]->g, od->al.mol_info[object_num]->k);
              */
              //od->al.mol_info[object_num]->ln_k = 0.0;
            }
          }
        }
        if (!eof) {
          return PREMATURE_EOF;
        }
      }
      fclose(od->file[ASCII_IN]->handle);
      od->file[ASCII_IN]->handle = NULL;
    }
    if ((!(*state)) && object_num) {
      od->grid.struct_num = od->al.mol_info[object_num - 1]->struct_num + 1;
    }
    if (!(od->cosmotherm_object_num[*state])) {
      memset(regex_name[*state], 0, BUF_LEN);
    }
  }

  return 0;
}


void calc_conf_energies(O3Data *od)
{
  char min_g_assigned;
  int i;
  int j;
  int state;
  int flag;
  int object_num[MAX_COSMOTHERM_STATES];
  int same_conf_object_num[MAX_COSMOTHERM_STATES];
  int old_struct_num[MAX_COSMOTHERM_STATES];
  double min_g;
  
  
  memset(object_num, 0, MAX_COSMOTHERM_STATES * sizeof(int));
  /*
  compute conformer g relative to the global minimum
  for each structure
  */
  flag = 1;
  while (flag) {
    for (state = 0, min_g_assigned = 0, min_g = 0.0; state < MAX_COSMOTHERM_STATES; ++state) {
      if (!(od->al.regex_list[state])) {
        continue;
      }
      old_struct_num[state] = od->al.regex_list[state][object_num[state]]->struct_num;
      same_conf_object_num[state] = object_num[state];
      od->al.cosmo_list[old_struct_num[state]]->n_conf[state] = 0;
      while ((same_conf_object_num[state] < od->cosmotherm_object_num[state])
        && (od->al.regex_list[state][same_conf_object_num[state]]->struct_num == old_struct_num[state])) {
        ++same_conf_object_num[state];
        ++(od->al.cosmo_list[old_struct_num[state]]->n_conf[state]);
      }
      for (i = object_num[state]; i < same_conf_object_num[state]; ++i) {
        if ((!state) && (!get_object_attr(od, i, ACTIVE_BIT | PREDICT_BIT))) {
          continue;
        }
        if ((!min_g_assigned) || (od->al.regex_list[state][i]->g < min_g)) {
          min_g = od->al.regex_list[state][i]->g;
          min_g_assigned = 1;
        }
      }
    }
    for (state = 0; state < MAX_COSMOTHERM_STATES; ++state) {
      if (!(od->al.regex_list[state])) {
        continue;
      }
      for (i = object_num[state]; i < same_conf_object_num[state]; ++i) {
        od->al.regex_list[state][i]->rel_g = 0.0;
        if ((!state) && (!get_object_attr(od, i, ACTIVE_BIT | PREDICT_BIT))) {
          continue;
        }
        od->al.regex_list[state][i]->rel_g = od->al.regex_list[state][i]->g - min_g;
      }
      object_num[state] = same_conf_object_num[state];
      if (flag) {
        flag = (object_num[state] < od->cosmotherm_object_num[state]);
      }
    }
  }
  /*
  compute z_water
  */
  memset(object_num, 0, MAX_COSMOTHERM_STATES * sizeof(int));
  for (i = 0; i < od->grid.struct_num; ++i) {
    od->al.cosmo_list[i]->z_water = 0.0;
    for (state = 0; state < MAX_COSMOTHERM_STATES; ++state) {
      if (!(od->al.regex_list[state])) {
        continue;
      }
      for (j = 0; j < od->al.cosmo_list[i]->n_conf[state]; ++j, ++object_num[state]) {
        od->al.cosmo_list[i]->z_water +=
          exp(-(od->al.regex_list[state][object_num[state]]->rel_g));
      }
    }
  }
}


void calc_conf_weights_training_set(O3Data *od)
{
  int object_num;
  int struct_num;
  int conf_num;
  int n_conf;


  /*
  compute z_bound
  */
  object_num = 0;
  while (object_num < od->grid.object_num) {
    struct_num = od->al.mol_info[object_num]->struct_num;
    conf_num = od->al.cosmo_list[struct_num]->n_conf[BOUND];
    od->al.cosmo_list[struct_num]->z_bound = 0.0;
    for (n_conf = 0; n_conf < conf_num; ++n_conf, ++object_num) {
      if (!get_object_attr(od, object_num, ACTIVE_BIT)) {
        continue;
      }
      od->al.mol_info[object_num]->exp_g_minus_ln_k =
        exp(-(od->al.regex_list[BOUND][object_num]->rel_g
        - od->al.mol_info[object_num]->ln_k));
      od->al.cosmo_list[struct_num]->z_bound +=
        od->al.mol_info[object_num]->exp_g_minus_ln_k;
    }
    if (od->al.cosmo_list[struct_num]->z_bound > 0.0) {
      od->al.cosmo_list[struct_num]->thermo_ln_k =
        log(od->al.cosmo_list[struct_num]->z_bound)
        - log(od->al.cosmo_list[struct_num]->z_water);
      od->al.cosmo_list[struct_num]->delta =
        od->al.cosmo_list[struct_num]->thermo_ln_k
        - od->al.cosmo_list[struct_num]->orig_y / M_LOG10E;
    }
    /*
    tee_printf(od, "%4d) ln_K_exp = %.4lf, thermo_ln_K = %.4lf, delta = %.4lf\n",
      struct_num + 1,
      od->al.cosmo_list[struct_num]->orig_y / M_LOG10E,
      od->al.cosmo_list[struct_num]->thermo_ln_k,
      od->al.cosmo_list[struct_num]->delta);
    */
  }
  object_num = 0;
  while (object_num < od->grid.object_num) {
    struct_num = od->al.mol_info[object_num]->struct_num;
    conf_num = od->al.cosmo_list[struct_num]->n_conf[BOUND];
    for (n_conf = 0; n_conf < conf_num; ++n_conf, ++object_num) {
      if (!get_object_attr(od, object_num, ACTIVE_BIT)) {
        continue;
      }
      od->mel.object_weight[object_num] =
        od->al.mol_info[object_num]->exp_g_minus_ln_k
        / od->al.cosmo_list[struct_num]->z_bound;
      set_y_value(od, object_num, 0,
        od->al.mol_info[object_num]->ln_k
        - od->al.cosmo_list[struct_num]->delta);
    }
  }
}


int update_conf_ln_k(O3Data *od, int model_type, int pc_num, double *ln_k_rmsd, int conv_method)
{
  int n;
  int y = 0;
  int y_max;
  int y_vars;
  int pos1;
  int object_num = 0;
  int actual_len;
  int num_values;
  int object_is_in_list;
  double sum = 0.0;
  double value = 0.0;
  double new_ln_k;
  double sqrt_weight = 0.0;
  
  
  if (model_type == FULL_MODEL) {
    double_mat_resize(od->mal.y_loadings, 1, pc_num);
    double_mat_resize(od->mal.x_scores, od->active_object_num, pc_num);
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,
      od->mal.x_scores->m,
      od->mal.y_loadings->m,
      pc_num, 1.0,
      od->mal.x_scores->base,
      od->mal.x_scores->max_m,
      od->mal.y_loadings->base,
      od->mal.y_loadings->max_m, 0.0,
      od->mal.pred_f_mat->base,
      od->mal.pred_f_mat->max_m);
  }
  else {
    od->mal.pred_f_mat = double_mat_resize
      (od->mal.pred_f_mat, od->active_object_num, od->y_vars);
    od->mel.predicted_object_list = (int *)realloc
      (od->mel.predicted_object_list, sizeof(int) * od->object_num);
    if (!(od->mel.predicted_object_list)) {
      return OUT_OF_MEMORY;
    }
    for (object_num = 0, y = 0; object_num < od->object_num; ++object_num) {
      if (!get_object_attr(od, object_num, ACTIVE_BIT)) {
        continue;
      }
      num_values = 0;
      sum = 0.0;
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
          /*
          move to the pc_num-th component, 0-th var,
          object_is_in_list-th value
          */
          if (fseek(od->file[TEMP_PRED]->handle,
            object_is_in_list * sizeof(double), SEEK_CUR)) {
            return CANNOT_READ_TEMP_FILE;
          }
          /*
          Read predicted value for 0-th variable,
          j-th component, object_is_in_list-th value
          */
          actual_len = fread(&value, sizeof(double), 1,
            od->file[TEMP_PRED]->handle);
          if (actual_len != 1) {
            return CANNOT_READ_TEMP_FILE;
          }
          sum += value;
          /*
          go back to starting position for this block
          */
          if (fseek(od->file[TEMP_PRED]->handle,
            pos1, SEEK_SET)) {
            return CANNOT_READ_TEMP_FILE;
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
      M_POKE(od->mal.pred_f_mat, y, 0, sum / (double)num_values);
      ++y;
    }
  }
  *ln_k_rmsd = 0.0;
  for (object_num = 0, y = 0; object_num < od->object_num; ++object_num) {
    if (!get_object_attr(od, object_num, ACTIVE_BIT)) {
      continue;
    }
    sqrt_weight = sqrt(od->mel.object_weight[object_num]);
    new_ln_k = ((model_type == FULL_MODEL) ? ((sqrt_weight > 0.0)
      ? M_PEEK(od->mal.pred_f_mat, y, 0) / sqrt_weight
      + od->vel.f_mat_full_ave->ve[0] : INFINITY)
      : M_PEEK(od->mal.pred_f_mat, y, 0));
    /*
    tee_printf(od, "%8d, old_ln_k = %16.4lf, new_ln_k = %16.4lf\n", object_num + 1, od->al.mol_info[object_num]->ln_k, new_ln_k);
    */
    /*
    tee_printf(od, "%8d, ln_k = %16.4lf, rel_g = %16.4lf\n", object_num + 1, od->al.mol_info[object_num]->ln_k, od->al.regex_list[BOUND][object_num]->rel_g);
    */
    *ln_k_rmsd += square(od->al.mol_info[object_num]->ln_k - new_ln_k);
    od->al.mol_info[object_num]->ln_k = (conv_method
      ? (od->al.mol_info[object_num]->ln_k + new_ln_k) / 2.0 : new_ln_k);
    ++y;
  }
  *ln_k_rmsd = (y ? sqrt(*ln_k_rmsd / (double)y) : MAX_CUTOFF);
  
  return 0;
}


void replace_orig_y(O3Data *od)
{
  int i;
  int j;
  int object_num;


  for (i = 0, object_num = 0; i < od->grid.struct_num; ++i) {
    od->al.cosmo_list[i]->orig_y = get_y_value(od, object_num, 0, 0);
    for (j = 0; j < od->al.cosmo_list[i]->n_conf[BOUND]; ++j, ++object_num) {
      set_y_value(od, object_num, 0, od->al.cosmo_list[i]->orig_y / M_LOG10E);
    }
  }
}


void restore_orig_y(O3Data *od)
{
  int i;
  int j;
  int object_num;


  for (i = 0, object_num = 0; i < od->grid.struct_num; ++i) {
    for (j = 0; j < od->al.cosmo_list[i]->n_conf[BOUND]; ++j, ++object_num) {
      set_y_value(od, object_num, 0, od->al.cosmo_list[i]->orig_y);
    }
  }
}


void calc_conf_weights_test_set(O3Data *od, int pc_num)
{
  int y;
  int x_max_x;
  int x_max_y;
  int y_max;
  int object_num;
  int struct_num;
  int conf_num;
  int n_conf;
  int info;
  #if (!defined HAVE_LIBLAPACK_ATLAS) && (!defined HAVE_LIBSUNPERF)
  int lwork;
  #endif
  

  x_max_x = od->overall_active_x_vars;
  x_max_y = od->y_vars;
  y_max = od->active_object_num;
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
  double_mat_resize(od->mal.x_weights_star, x_max_x, pc_num);
  double_mat_resize(od->mal.y_loadings, x_max_y, pc_num);
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,
    od->mal.x_weights_star->m,
    od->mal.y_loadings->m,
    pc_num, 1.0,
    od->mal.x_weights_star->base,
    od->mal.x_weights_star->max_m,
    od->mal.y_loadings->base,
    od->mal.y_loadings->max_m, 0.0,
    od->mal.b_coefficients->base,
    od->mal.b_coefficients->max_m);
  double_mat_resize(od->mal.f_mat, od->ext_pred_object_num, 1);
  double_mat_resize(od->mal.e_mat, od->ext_pred_object_num, od->mal.e_mat->n);
  double_mat_resize(od->mal.pred_f_mat, od->ext_pred_object_num, 1);
  for (object_num = 0, y = 0; object_num < od->grid.object_num; ++object_num) {
    if (get_object_attr(od, object_num, PREDICT_BIT)) {
      fill_x_vector(od, object_num, y, FULL_MODEL, 0);
      ++y;
    }
  }
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
    od->mal.e_mat->m,
    od->mal.b_coefficients->n,
    od->mal.e_mat->n, 1.0,
    od->mal.e_mat->base,
    od->mal.e_mat->max_m,
    od->mal.b_coefficients->base,
    od->mal.b_coefficients->max_m, 0.0,
    od->mal.pred_f_mat->base,
    od->mal.pred_f_mat->max_m);
  /*
  compute z_bound
  */
  object_num = 0;
  y = 0;
  while (object_num < od->grid.object_num) {
    struct_num = od->al.mol_info[object_num]->struct_num;
    conf_num = od->al.cosmo_list[struct_num]->n_conf[BOUND];
    od->al.cosmo_list[struct_num]->z_bound = 0.0;
    for (n_conf = 0; n_conf < conf_num; ++n_conf, ++object_num) {
      if (!get_object_attr(od, object_num, PREDICT_BIT)) {
        continue;
      }
      od->al.mol_info[object_num]->ln_k =
        M_PEEK(od->mal.pred_f_mat, y, 0)
        + od->vel.f_mat_full_ave->ve[0];
      od->al.mol_info[object_num]->exp_g_minus_ln_k =
        exp(-(od->al.regex_list[BOUND][object_num]->rel_g
        - od->al.mol_info[object_num]->ln_k));
      od->al.cosmo_list[struct_num]->z_bound +=
        od->al.mol_info[object_num]->exp_g_minus_ln_k;
      ++y;
    }
  }
  object_num = 0;
  while (object_num < od->grid.object_num) {
    struct_num = od->al.mol_info[object_num]->struct_num;
    conf_num = od->al.cosmo_list[struct_num]->n_conf[BOUND];
    for (n_conf = 0; n_conf < conf_num; ++n_conf, ++object_num) {
      if (!get_object_attr(od, object_num, PREDICT_BIT)) {
        continue;
      }
      od->mel.object_weight[object_num] =
        od->al.mol_info[object_num]->exp_g_minus_ln_k
        / od->al.cosmo_list[struct_num]->z_bound;
    }
  }
}
