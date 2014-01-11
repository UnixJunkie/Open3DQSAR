/*

cv.c

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


int cv(O3Data *od, int suggested_pc_num,
  int model_type, int cv_type, int groups, int runs)
{
  int i;
  int j;
  int x;
  int y;
  int y2;
  int group_num;
  int object_num = 0;
  int object_num2 = 0;
  int object_count;
  int struct_num = 0;
  int struct_num2 = 0;
  int struct_count;
  int active_struct_num = 0;
  int conf_num = 0;
  int conf_num2 = 0;
  int n_conf;
  int result;
  IntMat *group_composition;
  double sumweight = 0.0;
  double cum_press;
  double sd_sdep;
  double temp;
  

  result = 0;
  od->cv.overall_cv_runs = 0;
  switch (cv_type) {
    case LEAVE_ONE_OUT:
    memset(od->mal.press->base, 0,
      od->mal.press->m * od->mal.press->n * sizeof(double));
    int_perm_resize(od->pel.out_structs, 1);
    while (object_num < od->object_num) {
      struct_num = od->al.mol_info[object_num]->struct_num;
      conf_num = 1;
      for (n_conf = 0, sumweight = 0.0, y = 0;
        n_conf < conf_num; ++n_conf, ++object_num) {
        if (get_object_attr(od, object_num, ACTIVE_BIT)) {
          sumweight += od->mel.object_weight[object_num];
          ++y;
        }
      }
      if (sumweight > 0.0) {
        /*
        one struct at a time is going to be left out
        */
        od->pel.out_structs->pe[0] = struct_num;
        if (model_type & (FFDSEL_CV_MODEL | UVEPLS_CV_MODEL)) {
          trim_mean_center_x_matrix_hp(od,
            model_type, od->active_object_num - y,
            od->cv.overall_cv_runs);
          trim_mean_center_y_matrix_hp(od,
            od->active_object_num - y,
            od->cv.overall_cv_runs);
        }
        else if (model_type & SCRAMBLE_CV_MODEL) {
          trim_mean_center_x_matrix_hp(od,
            model_type, od->active_object_num - y,
            od->cv.overall_cv_runs);
          trim_mean_center_matrix(od, od->mal.large_f_mat,
            &(od->mal.f_mat), &(od->vel.f_mat_ave),
            model_type, od->active_object_num - y);
        }
        else {
          trim_mean_center_matrix(od, od->mal.large_e_mat,
            &(od->mal.e_mat), &(od->vel.e_mat_ave),
            model_type, od->active_object_num - y);
          trim_mean_center_matrix(od, od->mal.large_f_mat,
            &(od->mal.f_mat), &(od->vel.f_mat_ave),
            model_type, od->active_object_num - y);
        }
        pls(od, suggested_pc_num, model_type);
        result = pred_y_values(od, NULL, od->pc_num, model_type,
          od->cv.overall_cv_runs);
        ++(od->cv.overall_cv_runs);
        if (result) {
          return CANNOT_WRITE_TEMP_FILE;
        }
      }
    }
    if (model_type & CV_MODEL) {
      result = print_pred_values(od);
      if (result) {
        return result;
      }
      tee_printf(od, "\nPC%12s%12s\n",
        "SDEP", "q2");
      tee_printf(od, "--------------------------\n");
    }
    for (i = 0; i <= od->pc_num; ++i) {
      cum_press = 0.0;
      for (x = 0; x < od->y_vars; ++x) {
        cum_press += M_PEEK(od->mal.press, i, x);
      }
      od->vel.ave_sdep->ve[i] =
        sqrt(cum_press / (double)(od->cv.num_predictions));
      od->vel.q2->ve[i] = 1.0 - cum_press / od->cv.tss;
      if (model_type & SCRAMBLE_CV_MODEL) {
        od->vel.secv->ve[i] =
          sqrt(cum_press / (od->y_vars
          * ((double)(od->cv.num_predictions) - i - 1)));
      }
      if (model_type & CV_MODEL) {
        tee_printf(od, "%2d%12.4lf%12.4lf\n", i,
          od->vel.ave_sdep->ve[i],
          od->vel.q2->ve[i]);
      }
    }
    break;

    case LEAVE_TWO_OUT:
    memset(od->mal.press->base, 0,
      od->mal.press->m * od->mal.press->n * sizeof(double));
    int_perm_resize(od->pel.out_structs, 2);
    get_attr_struct_ave(od, 0, ACTIVE_BIT, &active_struct_num, NULL);
    while (object_num < od->object_num) {
      struct_num = od->al.mol_info[object_num]->struct_num;
      conf_num = 1;
      for (n_conf = 0, sumweight = 0.0, y = 0;
        n_conf < conf_num; ++n_conf, ++object_num) {
        if (get_object_attr(od, object_num, ACTIVE_BIT)) {
          sumweight += od->mel.object_weight[object_num];
          ++y;
        }
      }
      if (sumweight > 0.0) {
        object_num2 = object_num;
        while (object_num2 < od->object_num) {
          struct_num2 = od->al.mol_info[object_num2]->struct_num;
          conf_num2 = 1;
          for (n_conf = 0, sumweight = 0.0, y2 = 0;
            n_conf < conf_num2; ++n_conf, ++object_num2) {
            if (get_object_attr(od, object_num2, ACTIVE_BIT)) {
              sumweight += od->mel.object_weight[object_num2];
              ++y2;
            }
          }
          if (sumweight > 0.0) {
            /*
            two objects at a time are going to be left out
            */
            od->pel.out_structs->pe[0] = struct_num;
            od->pel.out_structs->pe[1] = struct_num2;
            if (model_type & (FFDSEL_CV_MODEL | UVEPLS_CV_MODEL)) {
              trim_mean_center_x_matrix_hp(od,
                model_type, od->active_object_num - (y + y2),
                od->cv.overall_cv_runs);
              trim_mean_center_y_matrix_hp(od,
                od->active_object_num - (y + y2),
                od->cv.overall_cv_runs);
            }
            else if (model_type & SCRAMBLE_CV_MODEL) {
              trim_mean_center_x_matrix_hp(od,
                model_type, od->active_object_num - (y + y2),
                od->cv.overall_cv_runs);
              trim_mean_center_matrix(od, od->mal.large_f_mat,
                &(od->mal.f_mat), &(od->vel.f_mat_ave),
                model_type, od->active_object_num - (y + y2));
            }
            else {
              trim_mean_center_matrix(od, od->mal.large_e_mat,
                &(od->mal.e_mat), &(od->vel.e_mat_ave),
                model_type, od->active_object_num - (y + y2));
              trim_mean_center_matrix(od, od->mal.large_f_mat,
                &(od->mal.f_mat), &(od->vel.f_mat_ave),
                model_type, od->active_object_num - (y + y2));
            }
            pls(od, suggested_pc_num, model_type);
            result = pred_y_values(od, NULL, od->pc_num, model_type,
              od->cv.overall_cv_runs);
            ++(od->cv.overall_cv_runs);
            if (result) {
              return CANNOT_WRITE_TEMP_FILE;
            }
          }
        }
      }
    }
    if (model_type & CV_MODEL) {
      result = print_pred_values(od);
      if (result) {
        return result;
      }
      tee_printf(od, "\nPC%12s%12s\n",
        "SDEP", "q2");
      tee_printf(od, "--------------------------\n");
    }
    for (i = 0; i <= od->pc_num; ++i) {
      cum_press = 0.0;
      for (x = 0; x < od->y_vars; ++x) {
        cum_press += M_PEEK(od->mal.press, i, x);
      }
      od->vel.ave_sdep->ve[i] =
        sqrt(cum_press / (double)(od->cv.num_predictions));
      od->vel.q2->ve[i] = 1.0 - cum_press / od->cv.tss;
      if (model_type & SCRAMBLE_CV_MODEL) {
        od->vel.secv->ve[i] =
          sqrt(cum_press / (od->y_vars
          * ((double)(od->cv.num_predictions) - i - 1)));
      }
      if (model_type & CV_MODEL) {
        tee_printf(od, "%2d%12.4lf%12.4lf\n", i,
          od->vel.ave_sdep->ve[i],
          od->vel.q2->ve[i]);
      }
    }
    break;
    
    case LEAVE_MANY_OUT:
    for (i = 0; i < runs; ++i) {
      memset(od->mal.press->base, 0,
        od->mal.press->m * od->mal.press->n * sizeof(double));
      od->cv.num_predictions = 0;
      group_composition = od->cimal.group_composition_list[i];
      for (group_num = 0; group_num < groups; ++group_num) {
        /*
        initialize the left-out objects vector before PLS
        */
        int_perm_resize(od->pel.out_structs,
          od->mel.struct_per_group[group_num]);
        for (struct_count = 0, object_count = 0; struct_count
          < od->mel.struct_per_group[group_num]; ++struct_count) {
          od->pel.out_structs->pe[struct_count] =
            group_composition->me[group_num][struct_count];
          ++object_count;
        }
        qsort(od->pel.out_structs->pe,
          od->pel.out_structs->size, sizeof(int),
          compare_integers);
        if (model_type & (FFDSEL_CV_MODEL | UVEPLS_CV_MODEL)) {
          trim_mean_center_x_matrix_hp(od, model_type,
            od->active_object_num - object_count,
            od->cv.overall_cv_runs);
          trim_mean_center_y_matrix_hp(od,
            od->active_object_num - object_count,
            od->cv.overall_cv_runs);
        }
        else if (model_type & SCRAMBLE_CV_MODEL) {
          trim_mean_center_x_matrix_hp(od, model_type,
            od->active_object_num - object_count,
            od->cv.overall_cv_runs);
          trim_mean_center_matrix(od, od->mal.large_f_mat,
            &(od->mal.f_mat), &(od->vel.f_mat_ave),
            model_type, od->active_object_num - object_count);
        }
        else {
          trim_mean_center_matrix(od, od->mal.large_e_mat,
            &(od->mal.e_mat), &(od->vel.e_mat_ave),
            model_type, od->active_object_num - object_count);
          trim_mean_center_matrix(od, od->mal.large_f_mat,
            &(od->mal.f_mat), &(od->vel.f_mat_ave),
            model_type, od->active_object_num - object_count);
        }
        pls(od, suggested_pc_num, model_type);
        result = pred_y_values(od, NULL, od->pc_num, model_type,
          od->cv.overall_cv_runs);
        if (result) {
          return CANNOT_WRITE_TEMP_FILE;
        }
        ++(od->cv.overall_cv_runs);
        od->cv.num_predictions += (od->y_vars
          * od->mel.struct_per_group[group_num]);
      }
      for (j = 0; j <= od->pc_num; ++j) {
        cum_press = 0.0;
        for (x = 0; x < od->y_vars; ++x) {
          M_POKE(od->mal.ave_press, j, x,
            M_PEEK(od->mal.ave_press, j, x)
            + M_PEEK(od->mal.press, j, x));
          cum_press += M_PEEK(od->mal.press, j, x);
          
        }
        M_POKE(od->mal.sdep_mat, i, j,
          sqrt(cum_press / (double)(od->cv.num_predictions)));
      }
    }
    od->cv.num_predictions *= runs;
    if (model_type & CV_MODEL) {
      result = print_pred_values(od);
      if (result) {
        return result;
      }
      tee_printf(od, "\nPC%12s%12s%12s\n",
        "SDEP", "SD on SDEP", "Average q2");
      tee_printf(od, "--------------------------------------\n");
    }
    for (j = 0; j <= od->pc_num; ++j) {
      for (i = 0; i < runs; ++i) {
        if (!i) {
          od->vel.ave_sdep->ve[j] =
            M_PEEK(od->mal.sdep_mat, i, j);
          od->vel.sum1->ve[j] = 0.0;
          od->vel.sumweight->ve[j] = 1.0;
        }
        else {
          temp = 1.0 + od->vel.sumweight->ve[j];
          od->vel.sum1->ve[j] +=
            (od->vel.sumweight->ve[j]
            * square(M_PEEK(od->mal.sdep_mat, i, j)
            - od->vel.ave_sdep->ve[j]) / temp);
          od->vel.ave_sdep->ve[j] +=
            ((M_PEEK(od->mal.sdep_mat, i, j)
            - od->vel.ave_sdep->ve[j]) / temp);
          od->vel.sumweight->ve[j] = temp;
        }
      }
      cum_press = 0.0;
      sd_sdep = sqrt(od->vel.sum1->ve[j]
        * runs / ((runs - 1)
        * od->vel.sumweight->ve[j]));
      for (x = 0; x < od->y_vars; ++x) {
        cum_press += M_PEEK(od->mal.ave_press, j, x);
      }
      od->vel.q2->ve[j] = 1.0 - cum_press / od->cv.tss;
      if (model_type & SCRAMBLE_CV_MODEL) {
        od->vel.secv->ve[j] =
          sqrt(cum_press / (od->y_vars
          * ((double)(od->cv.num_predictions) - j - 1)));
      }
      if (model_type & CV_MODEL) {
        tee_printf(od, "%2d%12.4lf%12.4lf%12.4lf\n", j,
          od->vel.ave_sdep->ve[j],
          sd_sdep, od->vel.q2->ve[j]);
      }
    }
    break;
  }

  return 0;
}
