/*

free.c

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


void free_atom_array(O3Data *od)
{
  int i = 0;
  
  
  for (i = 0; i < od->grid.object_num; ++i) {
    if (od->al.mol_info[i] && od->al.mol_info[i]->atom) {
      free_array(od->al.mol_info[i]->atom);
      od->al.mol_info[i]->atom = NULL;
    }
  }
}


void free_array(void *array)
{
  int i = 0;
  
  
  if (array) {
    while (((char **)array)[i]) {
      free(((char **)array)[i]);
      ((char **)array)[i] = NULL;
      ++i;
    }
    free(array);
  }
}


void free_char_matrix(CharMat *char_mat)
{
  int i;
  
  
  if (char_mat) {
    if (char_mat->me) {
      for (i = 0; i < char_mat->m; ++i) {
        if (char_mat->me[i]) {
          free(char_mat->me[i]);
          char_mat->me[i] = NULL;
        }
      }
      free(char_mat->me);
      char_mat->me = NULL;
    }
    free(char_mat);
  }
}


void free_mem(O3Data *od)
{
  char **array;
  char **mem;
  char **vec;
  char **mat;
  char **cimat;
  char **perm;
  int i;
  int n;
  O3Data *thread_od;
  
  
  free_x_var_array(od);
  thread_od = od;
  for (n = od->n_proc - 1; n >= 0; --n) {
    if (n) {
      if (!(od->mel.thread_info[n])) {
        continue;
      }
      thread_od = &(od->mel.thread_info[n]->od);
    }
    array = (char **)&(thread_od->al);
    i = 0;
    while (i < (sizeof(ArrayList) / sizeof(char *))) {
      if (array[i]) {
        free_array(array[i]);
        array[i] = NULL;
      }
      ++i;
    }
    mem = (char **)&(thread_od->mel);
    i = 0;
    while (i < (sizeof(MemList) / sizeof(char *))) {
      if (mem[i]) {
        free(mem[i]);
        mem[i] = NULL;
      }
      ++i;
    }
    vec = (char **)&(thread_od->vel);
    i = 0;
    while (i < (sizeof(VecList) / sizeof(DoubleVec *))) {
      if (vec[i]) {
        double_vec_free((DoubleVec *)vec[i]);
        vec[i] = NULL;
      }
      ++i;
    }
    mat = (char **)&(thread_od->mal);
    i = 0;
    while (i < (sizeof(MatList) / sizeof(DoubleMat *))) {
      if (mat[i]) {
        double_mat_free((DoubleMat *)mat[i]);
        mat[i] = NULL;
      }
      ++i;
    }
    cimat = (char **)&(thread_od->cimal);
    i = 0;
    while (i < (sizeof(CIMatList) / sizeof(CharMat *))) {
      if (cimat[i]) {
        free_char_matrix((CharMat *)cimat[i]);
        cimat[i] = NULL;
      }
      ++i;
    }
    perm = (char **)&(thread_od->pel);
    i = 0;
    while (i < (sizeof(PermList) / sizeof(IntPerm *))) {
      if (perm[i]) {
        int_perm_free((IntPerm *)perm[i]);
        perm[i] = NULL;
      }
      ++i;
    }
  }
}


void free_threads(O3Data *od)
{
  int i;

  
  for (i = 0; i < od->n_proc; ++i) {
    if (od->mel.thread_info[i]) {
      free(od->mel.thread_info[i]);
      od->mel.thread_info[i] = NULL;
    }
  }
  if (od->mel.mutex) {
    free(od->mel.mutex);
    od->mel.mutex = NULL;
  }
}


void free_x_var_array(O3Data *od)
{
  int i;
  int j;
  
  
  for (i = 0; i < od->field_num; ++i) {
    if (od->mel.x_var_array) {
      if (od->mel.x_var_array[i]) {
        for (j = 0; j < od->object_num; ++j) {
          if ((!(od->save_ram)) && od->mel.x_var_array[i][j]) {
            free(od->mel.x_var_array[i][j]);
            od->mel.x_var_array[i][j] = NULL;
          }
        }
        free(od->mel.x_var_array[i]);
        od->mel.x_var_array[i] = NULL;
      }
    }
    if (od->mel.x_var_attr[i]) {
      free(od->mel.x_var_attr[i]);
      od->mel.x_var_attr[i] = NULL;
    }
    for (j = 0; j < MAX_VAR_BUF; ++j) {
      if (od->mel.x_var_buf[j][i]) {
        free(od->mel.x_var_buf[j][i]);
        od->mel.x_var_buf[j][i] = NULL;
      }
    }
  }
  od->field_num = 0;
}


void free_y_var_array(O3Data *od)
{
  int i;
  
  
  if (od->mel.y_var_array) {
    free(od->mel.y_var_array);
  }
  od->mel.y_var_array = NULL;
  if (od->mel.y_var_attr) {
    free(od->mel.y_var_attr);
  }
  od->mel.y_var_attr = NULL;
  for (i = 0; i < MAX_VAR_BUF; ++i) {
    if (od->mel.y_var_buf[i]) {
      free(od->mel.y_var_buf[i]);
      od->mel.y_var_buf[i] = NULL;
    }
  }
  if (od->mel.y_data) {
    free(od->mel.y_data);
    od->mel.y_data = NULL;
  }
  od->y_vars = 0;
}


void free_cv_groups(O3Data *od, int runs)
{
  int i;
  
  
  if (od->mel.struct_per_group) {
    free(od->mel.struct_per_group);
    od->mel.struct_per_group = NULL;
  }
  if (od->cimal.group_composition_list) {
    for (i = 0; i < runs; ++i) {
      if (od->cimal.group_composition_list[i]) {
        free_char_matrix((CharMat *)
          (od->cimal.group_composition_list[i]));
        od->cimal.group_composition_list[i] = NULL;
      }
    }
    free(od->cimal.group_composition_list);
    od->cimal.group_composition_list = NULL;
  }
  if (od->mal.sdep_mat) {
    double_mat_free(od->mal.sdep_mat);
    od->mal.sdep_mat = NULL;
  }
}


void free_cv_sdep(O3Data *od)
{
  if (od->mel.struct_list) {
    free(od->mel.struct_list);
    od->mel.struct_list = NULL;
  }
  if (od->mal.sdep_mat) {
    double_mat_free(od->mal.sdep_mat);
    od->mal.sdep_mat = NULL;
  }
  if (od->mal.ave_press) {
    double_mat_free(od->mal.ave_press);
    od->mal.ave_press = NULL;
  }
  if (od->vel.sumweight) {
    double_vec_free(od->vel.sumweight);
    od->vel.sumweight = NULL;
  }
  if (od->vel.sum1) {
    double_vec_free(od->vel.sum1);
    od->vel.sum1 = NULL;
  }
  if (od->vel.q2) {
    double_vec_free(od->vel.q2);
    od->vel.q2 = NULL;
  }
  if (od->vel.secv) {
    double_vec_free(od->vel.secv);
    od->vel.secv = NULL;
  }
  if (od->pel.sdep_rank) {
    int_perm_free(od->pel.sdep_rank);
    od->pel.sdep_rank = NULL;
  }
}


void free_parallel_cv(O3Data *od, ThreadInfo **thread_info,
  int model_type, int cv_type, int runs)
{
  int i;
  
  
  for (i = 1; i < od->cv.n_threads; ++i) {
    if (thread_info[i]->temp_pred.handle) {
      fclose(thread_info[i]->temp_pred.handle);
      thread_info[i]->temp_pred.handle = NULL;
      remove(thread_info[i]->temp_pred.name);
      thread_info[i]->temp_pred.name[0] = '\0';
    }
    if (thread_info[i]->temp_cv_coeff.handle) {
      fclose(thread_info[i]->temp_cv_coeff.handle);
      thread_info[i]->temp_cv_coeff.handle = NULL;
      remove(thread_info[i]->temp_cv_coeff.name);
      thread_info[i]->temp_cv_coeff.name[0] = '\0';
    }
    /*
    free all duplicate structures
    */
    free_pls(&(thread_info[i]->od));
    if (cv_type == LEAVE_MANY_OUT) {
      if (thread_info[i]->od.mal.press) {
        double_mat_free(thread_info[i]->od.mal.press);
        thread_info[i]->od.mal.press = NULL;
      }
    }
  }
  if (cv_type == LEAVE_MANY_OUT) {
    free_cv_groups(od, runs);
  }
  if ((!(model_type & SILENT_PLS)) && (od->file[TEMP_PRED]->handle)) {
    fclose(od->file[TEMP_PRED]->handle);
    thread_info[0]->temp_pred.handle = NULL;
    od->file[TEMP_PRED]->handle = NULL;
  }
  if (od->file[TEMP_CV_COEFF]->handle) {
    fclose(od->file[TEMP_CV_COEFF]->handle);
    thread_info[0]->temp_cv_coeff.handle = NULL;
    od->file[TEMP_CV_COEFF]->handle = NULL;
  }
  free_threads(od);
}


void free_pls(O3Data *od)
{
  if (od->mal.b_coefficients) {
    double_mat_free(od->mal.b_coefficients);
    od->mal.b_coefficients = NULL;
  }
  if (od->mal.e_mat) {
    double_mat_free(od->mal.e_mat);
    od->mal.e_mat = NULL;
  }
  if (od->mal.f_mat) {
    double_mat_free(od->mal.f_mat);
    od->mal.f_mat = NULL;
  }
  if (od->mal.pred_f_mat) {
    double_mat_free(od->mal.pred_f_mat);
    od->mal.pred_f_mat = NULL;
  }
  if (od->vel.e_mat_ave) {
    double_vec_free(od->vel.e_mat_ave);
    od->vel.e_mat_ave = NULL;
  }
  if (od->vel.f_mat_ave) {
    double_vec_free(od->vel.f_mat_ave);
    od->vel.f_mat_ave = NULL;
  }
  if (od->vel.v) {
    double_vec_free(od->vel.v);
    od->vel.v = NULL;
  }
  if (od->vel.v_new) {
    double_vec_free(od->vel.v_new);
    od->vel.v_new = NULL;
  }
  if (od->vel.ro) {
    double_vec_free(od->vel.ro);
    od->vel.ro = NULL;
  }
  if (od->vel.explained_s2_x) {
    double_vec_free(od->vel.explained_s2_x);
    od->vel.explained_s2_x = NULL;
  }
  if (od->vel.explained_s2_y) {
    double_vec_free(od->vel.explained_s2_y);
    od->vel.explained_s2_y = NULL;
  }
  if (od->vel.ave_sdep) {
    double_vec_free(od->vel.ave_sdep);
    od->vel.ave_sdep = NULL;
  }
  if (od->vel.ave_sdec) {
    double_vec_free(od->vel.ave_sdec);
    od->vel.ave_sdec = NULL;
  }
  if (od->vel.r2) {
    double_vec_free(od->vel.r2);
    od->vel.r2 = NULL;
  }
  if (od->vel.r2_pred) {
    double_vec_free(od->vel.r2_pred);
    od->vel.r2_pred = NULL;
  }
  if (od->mal.x_scores) {
    double_mat_free(od->mal.x_scores);
    od->mal.x_scores = NULL;
  }
  if (od->mal.x_weights) {
    double_mat_free(od->mal.x_weights);
    od->mal.x_weights = NULL;
  }
  if (od->mal.x_weights_star) {
    double_mat_free(od->mal.x_weights_star);
    od->mal.x_weights_star = NULL;
  }
  if (od->mal.x_loadings) {
    double_mat_free(od->mal.x_loadings);
    od->mal.x_loadings = NULL;
  }
  if (od->mal.y_loadings) {
    double_mat_free(od->mal.y_loadings);
    od->mal.y_loadings = NULL;
  }
  if (od->mal.y_scores) {
    double_mat_free(od->mal.y_scores);
    od->mal.y_scores = NULL;
  }
  if (od->mal.temp) {
    double_mat_free(od->mal.temp);
    od->mal.temp = NULL;
  }
  if (od->mel.ipiv) {
    free(od->mel.ipiv);
    od->mel.ipiv = NULL;
  }
  #ifdef HAVE_LIBMKL
  if (od->mel.work) {
    free(od->mel.work);
    od->mel.work = NULL;
  }
  #endif
  if (od->pel.out_structs) {
    free(od->pel.out_structs);
    od->pel.out_structs = NULL;
  }
}
