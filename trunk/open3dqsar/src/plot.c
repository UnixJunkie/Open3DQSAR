/*

plot.c

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


int plot(O3Data *od, char *filename, int type,
  int label, int requested_pc_num, int *pc_axis)
{
  char type_of_data[TITLE_LEN];
  char fx[BUF_LEN];
  char z_dimension[BUF_LEN];
  char plot_command[BUF_LEN];
  char with_labels[BUF_LEN];
  char y_axis_label[TITLE_LEN];
  int i;
  int j = 0;
  int n;
  int pos1;
  int pos2;
  int x;
  int y;
  int x_max_y;
  int y_max;
  int pc_num;
  int max_pc;
  int num_values;
  int actual_len;
  int object_is_in_list;
  int elem;
  int y_axis;
  unsigned short bit;
  double x_value;
  double y_value;
  double z_value;
  double abs_max_value;
  double value;
  double actual_value;
  DoubleMat *mat = NULL;
  DoubleVec *vec = NULL;
  FILE *plt_dat_out = NULL;
  FILE *plt_cmd_out = NULL;
  FILE *temp_handle = NULL;
  
  
  if (type & VS_EXP_PLOT) {
    if ((type & EXT_PRED_VS_EXP) || (type & RECALC_VS_EXP)) {
      /*
      if recalculated vs experimental
      or external predictions vs experimental
      were requested
      */
      if (type & RECALC_VS_EXP) {
        temp_handle = fopen(od->file[TEMP_CALC]->name, "rb");
        if (!temp_handle) {
          return CANNOT_READ_TEMP_FILE;
        }
        od->file[TEMP_CALC]->handle = temp_handle;
        strcpy(type_of_data, "Calculated values");
        bit = ACTIVE_BIT;
      }
      else {
        temp_handle = fopen(od->file[TEMP_EXT_PRED]->name, "rb");
        if (!temp_handle) {
          return CANNOT_READ_TEMP_FILE;
        }
        od->file[TEMP_EXT_PRED]->handle = temp_handle;
        strcpy(type_of_data, "External predictions");
        bit = PREDICT_BIT;
      }
      if (type & RESIDUALS) {
        strcat(type_of_data, " (residuals)");
      }
      /*
      Read number of PCs
      */
      actual_len = fread(&pc_num, sizeof(int), 1, temp_handle);
      if (actual_len != 1) {
        return CANNOT_READ_TEMP_FILE;
      }
      /*
      Read total number of variables, that is y_vars
      */
      actual_len = fread(&x_max_y, sizeof(int), 1, temp_handle);
      if (actual_len != 1) {
        return CANNOT_READ_TEMP_FILE;
      }
      /*
      Read number of objects in the first (and only) block
      */
      actual_len = fread(&y_max, sizeof(int), 1, temp_handle);
      if (actual_len != 1) {
        return CANNOT_READ_TEMP_FILE;
      }
      /*
      Skip object numbers: no need to read them, we already
      know which they are
      */
      if (fseek(temp_handle, y_max * sizeof(int), SEEK_CUR)) {
        return CANNOT_READ_TEMP_FILE;
      }
      pos1 = ftell(temp_handle);
      for (x = 0; x < x_max_y; ++x) {
        if (get_y_var_attr(od, x, OPERATE_BIT)) {
          sprintf(od->file[PLT_DAT_OUT]->name, "%s%s%02d%s",
            filename, Y_VAR_EXTENSION, x + 1, PLT_DAT_EXTENSION);
          if (!(plt_dat_out =
            fopen(od->file[PLT_DAT_OUT]->name, "w+"))) {
            return PREMATURE_PLT_DAT_EOF;
          }
          od->file[PLT_DAT_OUT]->handle = plt_dat_out;
          sprintf(od->file[PLT_CMD_OUT]->name, "%s%s%02d%s",
            filename, Y_VAR_EXTENSION, x + 1, PLT_CMD_EXTENSION);
          if (!(plt_cmd_out =
            fopen(od->file[PLT_CMD_OUT]->name, "w+"))) {
            return PREMATURE_PLT_CMD_EOF;
          }
          od->file[PLT_CMD_OUT]->handle = plt_cmd_out;
          /*
          move to the x-th variable
          */
          if (fseek(temp_handle, x * y_max * sizeof(double), SEEK_CUR)) {
            return CANNOT_READ_TEMP_FILE;
          }
          pos2 = ftell(temp_handle);
          if (type & RESIDUALS) {
            sprintf(fx,
              "f(x) = 0\n");
          }
          else {
            sprintf(fx,
              "f(x) = m * x + b\n"
              "fit f(x) \"%s\" using 2:3 via m, b\n",
              od->file[PLT_DAT_OUT]->name);
          }
          if (label) {
            sprintf(with_labels,
              "\"%s\" \\\n\t"
              "using 2:3:1 with labels left "
              "offset 0.5,0.2 lt -1 ,\\\n\t",
              od->file[PLT_DAT_OUT]->name);
          }
          else {
            memset(with_labels, 0, BUF_LEN);
          }
          fprintf(plt_cmd_out,
            "unset key\n"
            "unset parametric\n"
            "set title \"%s for dependent variable %2d (%s)\"\n"
            "show title\n"
            "set rmargin 16\n"
            "set xlabel \"Experimental values\" tc lt -1\n"
            "show xlabel\n"
            "set ylabel \"%s\" tc lt -1\n"
            "show ylabel\n"
            "%s"
            "plot \"%s\" \\\n\t"
            "using 2:3 with points lt rgb \"blue\" "
            "ps 0.75 pt 5, \\\n\t"
            "%s"
            "f(x) lt rgb \"red\" lw 2\n",
            type_of_data, x + 1,
            od->cimal.y_var_name->me[x],
            type_of_data, fx,
            od->file[PLT_DAT_OUT]->name,
            with_labels);
          y = 0;
          for (i = 0; i < od->object_num; ++i) {
            if (get_object_attr(od, i, bit)) {
              actual_value = M_PEEK(od->mal.large_f_mat,
                ((bit == PREDICT_BIT) ? od->active_object_num : 0) + y, x);
              /*
              move to the requested PC
              */
              if (fseek(temp_handle,
                (y_max * x_max_y * (pc_num - requested_pc_num) + y)
                * sizeof(double), SEEK_CUR)) {
                return CANNOT_READ_TEMP_FILE;
              }
              /*
              Read calculated/predicted value for y-th variable
              */
              actual_len = fread(&value, sizeof(double), 1, temp_handle);
              if (actual_len != 1) {
                return CANNOT_READ_TEMP_FILE;
              }
              if (type & RESIDUALS) {
                value -= actual_value;
              }
              if (label) {
                if (!pc_axis[2]) {
                  sprintf(with_labels,
                    "\"%s\" \\\n\t"
                    "using 2:3:1 with labels left "
                    "offset 0.5,0.2 lt -1 ,\\\n\t",
                    od->file[PLT_DAT_OUT]->name);
                }
                else {
                  sprintf(with_labels,
                    "\"%s\" \\\n\t"
                    "using 2:3:4:1 with labels left "
                    "offset 0.5,0.2 lt -1 ,\\\n\t",
                    od->file[PLT_DAT_OUT]->name);
                }
              }
              else {
                memset(with_labels, 0, BUF_LEN);
              }
              if (label == NAME_LABEL) {
                fprintf(plt_dat_out, "%-40s",
                  od->al.mol_info[i]->object_name);
              }
              else if (label == ID_LABEL) {
                fprintf(plt_dat_out, "%-10d",
                  od->al.mol_info[i]->object_id);
              }
              else {
                fprintf(plt_dat_out, "%-10d", i + 1);
              }
              fprintf(plt_dat_out, "%16.6le%16.6le\n",
                 actual_value, value);
              /*
              go back to starting position
              */
              if (fseek(temp_handle, pos2, SEEK_SET)) {
                return CANNOT_READ_TEMP_FILE;
              }
              ++y;
            }
          }
          if (fseek(temp_handle, pos1, SEEK_SET)) {
            return CANNOT_READ_TEMP_FILE;
          }
        }
        if (plt_dat_out) {
          fclose(plt_dat_out);
          od->file[PLT_DAT_OUT]->handle = NULL;
        }
        if (plt_cmd_out) {
          fclose(plt_cmd_out);
          od->file[PLT_CMD_OUT]->handle = NULL;
        }
      }
      if (temp_handle) {
        fclose(temp_handle);
        if (type & RECALC_VS_EXP) {
          od->file[TEMP_CALC]->handle = NULL;
        }
        else {
          od->file[TEMP_EXT_PRED]->handle = NULL;
        }
      }
    }
    else {
      temp_handle = fopen(od->file[TEMP_PRED]->name, "rb");
      if (!temp_handle) {
        return CANNOT_READ_TEMP_FILE;
      }
      od->file[TEMP_PRED]->handle = temp_handle;
      strcpy(type_of_data, "Predicted values");
      if (type & RESIDUALS) {
        strcat(type_of_data, " (residuals)");
      }
      /*
      if internal predictions vs experimental
      were requested
      */
      od->mal.pred_f_mat =
        double_mat_resize(od->mal.pred_f_mat,
        od->active_object_num, od->y_vars);
      od->mel.predicted_object_list = (int *)realloc
        (od->mel.predicted_object_list,
        sizeof(int) * od->object_num);
      if (!(od->mel.predicted_object_list)) {
        return OUT_OF_MEMORY;
      }
      od->mel.sum =
        (double *)realloc(od->mel.sum,
        (od->pc_num + 1) * sizeof(double));
      if (!(od->mel.sum)) {
        return OUT_OF_MEMORY;
      }
      for (x = 0; x < od->y_vars; ++x) {
        if (get_y_var_attr(od, x, OPERATE_BIT)) {
          sprintf(od->file[PLT_DAT_OUT]->name, "%s%s%02d%s",
            filename, Y_VAR_EXTENSION, x + 1, PLT_DAT_EXTENSION);
          if (!(plt_dat_out =
            fopen(od->file[PLT_DAT_OUT]->name, "w+"))) {
            return PREMATURE_PLT_DAT_EOF;
          }
          od->file[PLT_DAT_OUT]->handle = plt_dat_out;
          sprintf(od->file[PLT_CMD_OUT]->name, "%s%s%02d%s",
            filename, Y_VAR_EXTENSION, x + 1, PLT_CMD_EXTENSION);
          if (!(plt_cmd_out =
            fopen(od->file[PLT_CMD_OUT]->name, "w+"))) {
            return PREMATURE_PLT_CMD_EOF;
          }
          od->file[PLT_CMD_OUT]->handle = plt_cmd_out;
          if (type & RESIDUALS) {
            sprintf(fx,
              "f(x) = 0\n");
          }
          else {
            sprintf(fx,
              "f(x) = m * x + b\n"
              "fit f(x) \"%s\" using 2:3 via m, b\n",
              od->file[PLT_DAT_OUT]->name);
          }
          fprintf(plt_cmd_out,
            "unset key\n"
            "unset parametric\n"
            "set title \"%s for dependent variable %2d (%s)\"\n"
            "show title\n"
            "set rmargin 16\n"
            "set xlabel \"Experimental values\" tc lt -1\n"
            "show xlabel\n"
            "set ylabel \"%s\" tc lt -1\n"
            "show ylabel\n"
            "%s"
            "plot \"%s\" \\\n\t"
            "using 2:3 with points lt rgb \"blue\" "
            "ps 0.75 pt 5 ,\\\n\t"
            "%s"
            "f(x) lt rgb \"red\" lw 2\n",
            type_of_data, x + 1,
            od->cimal.y_var_name->me[x],
            type_of_data, fx,
            od->file[PLT_DAT_OUT]->name,
            with_labels);
          y = 0;
          for (i = 0; i < od->object_num; ++i) {
            if (get_object_attr(od, i, ACTIVE_BIT)) {
              num_values = 0;
              memset(od->mel.sum, 0,
                (od->pc_num + 1) * sizeof(double));
              actual_value = M_PEEK(od->mal.large_f_mat, y, x);
              rewind(temp_handle);
              while (!feof(temp_handle)) {
                /*
                Read number of PCs
                */
                actual_len = fread(&pc_num, sizeof(int), 1,
                  temp_handle);
                if (actual_len != 1) {
                  if (!feof(temp_handle)) {
                    return CANNOT_READ_TEMP_FILE;
                  }
                  else {
                    break;
                  }
                }
                /*
                Read total number of variables, that is y_vars
                */
                actual_len = fread(&od->y_vars, sizeof(int), 1,
                  temp_handle);
                if (actual_len != 1) {
                  return CANNOT_READ_TEMP_FILE;
                }
                /*
                Read number of objects in the next block
                */
                actual_len = fread(&y_max, sizeof(int), 1,
                  temp_handle);
                if (actual_len != 1) {
                  return CANNOT_READ_TEMP_FILE;
                }
                /*
                Read object numbers for this block
                */
                /*
                check if the object for which we want to write predicted values
                is present in this block
                */
                object_is_in_list = -1;
                for (n = 0; n < y_max; ++n) {
                  actual_len =
                    fread(&(od->mel.predicted_object_list[n]),
                    sizeof(int), 1, temp_handle);
                  if (actual_len != 1) {
                    return CANNOT_READ_TEMP_FILE;
                  }
                  /*
                  if the object is present, remember
                  its position in the block
                  */
                  if (od->mel.predicted_object_list[n] == i) {
                    object_is_in_list = n;
                  }
                }
                /*
                if so, read the value for the current var and add it;
                */
                if (object_is_in_list != -1) {
                  /*
                  remember the starting position
                  for the current block
                  */
                  ++num_values;
                  pos1 = ftell(temp_handle);
                  /*
                  move to the x-th var,
                  object_is_in_list-th value
                  */
                  if (fseek(temp_handle,
                    (y_max * od->y_vars * (pc_num - requested_pc_num)
                    + y_max * x + object_is_in_list) * sizeof(double), SEEK_CUR)) {
                    return CANNOT_READ_TEMP_FILE;
                  }
                  /*
                  Read predicted value for x-th variable,
                  j-th component, object_is_in_list-th value
                  */
                  actual_len = fread(&value, sizeof(double), 1, temp_handle);
                  if (actual_len != 1) {
                    return CANNOT_READ_TEMP_FILE;
                  }
                  if (type & RESIDUALS) {
                    value -= actual_value;
                  }
                  od->mel.sum[requested_pc_num] += value;
                  /*
                  go back to starting position for this block
                  */
                  if (fseek(temp_handle, pos1, SEEK_SET)) {
                    return CANNOT_READ_TEMP_FILE;
                  }
                }
                /*
                skip to the next block (if there is one)
                */
                if (fseek(temp_handle,
                  y_max * od->y_vars * (pc_num + 1)
                  * sizeof(double), SEEK_CUR)) {
                  return CANNOT_READ_TEMP_FILE;
                }
              }
              if (label) {
                if (!pc_axis[2]) {
                  sprintf(with_labels,
                    "\"%s\" \\\n\t"
                    "using 2:3:1 with labels left "
                    "offset 0.5,0.2 lt -1 ,\\\n\t",
                    od->file[PLT_DAT_OUT]->name);
                }
                else {
                  sprintf(with_labels,
                    "\"%s\" \\\n\t"
                    "using 2:3:4:1 with labels left "
                    "offset 0.5,0.2 lt -1, \\\n\t",
                    od->file[PLT_DAT_OUT]->name);
                }
              }
              else {
                memset(with_labels, 0, BUF_LEN);
              }
              if (label == NAME_LABEL) {
                fprintf(plt_dat_out, "%-40s",
                  od->al.mol_info[i]->object_name);
              }
              else if (label == ID_LABEL) {
                fprintf(plt_dat_out, "%-10d",
                  od->al.mol_info[i]->object_id);
              }
              else {
                fprintf(plt_dat_out, "%-10d", i + 1);
              }
              fprintf(plt_dat_out, "%16.6le%16.6le\n", actual_value,
                od->mel.sum[requested_pc_num] / num_values);
              ++y;
            }
          }
          if (plt_dat_out) {
            fclose(plt_dat_out);
            od->file[PLT_DAT_OUT]->handle = NULL;
          }
          if (plt_cmd_out) {
            fclose(plt_cmd_out);
            od->file[PLT_CMD_OUT]->handle = NULL;
          }
        }
      }
      if (temp_handle) {
        fclose(temp_handle);
        od->file[TEMP_PRED]->handle = NULL;
      }
    }
    return 0;
  }
  else {
    sprintf(od->file[PLT_DAT_OUT]->name,
      "%s%s", filename, PLT_DAT_EXTENSION);
    sprintf(od->file[PLT_CMD_OUT]->name,
      "%s%s", filename, PLT_CMD_EXTENSION);
    if (!(plt_dat_out =
      fopen(od->file[PLT_DAT_OUT]->name, "w+"))) {
      return PREMATURE_PLT_DAT_EOF;
    }
    od->file[PLT_DAT_OUT]->handle = plt_dat_out;
    if (!(plt_cmd_out =
      fopen(od->file[PLT_CMD_OUT]->name, "w+"))) {
      return PREMATURE_PLT_CMD_EOF;
    }
    od->file[PLT_CMD_OUT]->handle = plt_cmd_out;
  }
  if ((type & SDEC_VS_PC_PLOT)
    || (type & SDEP_VS_PC_PLOT)
    || (type & R2_VS_PC_PLOT)
    || (type & Q2_VS_PC_PLOT)) {
    max_pc = od->pc_num;
    if (type == SDEC_VS_PC_PLOT) {
      vec = od->vel.ave_sdec;
      strcpy(type_of_data, "SDEC");
    }
    else if (type == SDEP_VS_PC_PLOT) {
      vec = od->vel.ave_sdep;
      max_pc = od->cv.pc_num;
      strcpy(type_of_data, "SDEP");
    }
    else if (type == R2_VS_PC_PLOT) {
      vec = od->vel.r2;
      strcpy(type_of_data, "r2");
    }
    else if (type == Q2_VS_PC_PLOT) {
      vec = od->vel.q2;
      max_pc = od->cv.pc_num;
      strcpy(type_of_data, "q2");
    }
    for (i = 0; i <= max_pc; ++i) {
      fprintf(plt_dat_out,
        "%5d%16.6le\n", i, vec->ve[i]);
    }
    fprintf(plt_cmd_out,
      "unset key\n"
      "set parametric\n"
      "set title \"%s vs PC\"\n"
      "show title\n"
      "set rmargin 0\n"
      "set xlabel \"PCs\" tc lt -1\n"
      "show xlabel\n"
      "set ylabel \"%s\" tc lt -1\n"
      "show ylabel\n"
      "plot \"%s\" \\\n\t"
      "using 1:2 with points lt rgb \"blue\" "
      "ps 0.75 pt 5 ,\\\n\t"
      "\"%s\" \\\n\t"
      "using 1:2 lt rgb \"red\" lw 2 \\\n\t"
      "smooth bezier \\\n",
      type_of_data,
      type_of_data,
      od->file[PLT_DAT_OUT]->name,
      od->file[PLT_DAT_OUT]->name);
  }
  else if (type & PLS_X_VS_Y_SCORES_PLOT) {
    n = 0;
    for (i = 0; i < od->object_num; ++i) {
      if (get_object_attr(od, i, ACTIVE_BIT)) {
        y_value = M_PEEK(od->mal.y_scores,
          n, requested_pc_num - 1);
        x_value = M_PEEK(od->mal.x_scores,
          n, requested_pc_num - 1);
        if (label) {
          sprintf(with_labels,
            "\"%s\" \\\n\t"
            "using 2:3:1 with labels left "
            "offset 0.5,0.2 lt -1 ,\\\n\t",
            od->file[PLT_DAT_OUT]->name);
        }
        else {
          memset(with_labels, 0, BUF_LEN);
        }
        if (label == NAME_LABEL) {
          fprintf(plt_dat_out, "%-40s",
            od->al.mol_info[i]->object_name);
        }
        else if (label == ID_LABEL) {
          fprintf(plt_dat_out, "%-10d",
            od->al.mol_info[i]->object_id);
        }
        else {
          fprintf(plt_dat_out, "%-10d", i + 1);
        }
        fprintf(plt_dat_out, "%16.6le%16.6le\n",
           x_value, y_value);
        ++n;
      }
    }
    fprintf(plt_cmd_out,
      "unset key\n"
      "unset parametric\n"
      "set title \"PLS Y vs X scores\"\n"
      "show title\n"
      "set rmargin 16\n"
      "set xlabel \"X scores\" tc lt -1\n"
      "show xlabel\n"
      "set ylabel \"Y scores\" tc lt -1\n"
      "show ylabel\n"
      "f(x) = m * x + b\n"
      "fit f(x) \"%s\" using 2:3 via m, b\n"
      "plot \"%s\" \\\n\t"
      "using 2:3 with points lt rgb \"blue\" "
      "ps 0.75 pt 5 ,\\\n\t"
      "%s"
      "f(x) lt rgb \"red\" lw 2\n",
      od->file[PLT_DAT_OUT]->name,
      od->file[PLT_DAT_OUT]->name,
      with_labels);
  }
  else if ((type & SCRAMBLED_Q2_VS_R2)
    || (type & SCRAMBLED_SECV_VS_R2)) {
    y_axis = ((type & SCRAMBLED_Q2_VS_R2) ? Q2_FIT : SECV_FIT);
    if (y_axis == Q2_FIT) {
      strcpy(type_of_data, "q2");
      strcpy(y_axis_label, "q2");
    }
    else {
      strcpy(type_of_data, "SE(cv)");
      strcpy(y_axis_label, "SE(cv)");
    }
    sprintf(fx, "%.4le + (%.4le) * x + (%.4le) * x**2",
      od->vel.fit_coeff[y_axis]->ve[0],
      od->vel.fit_coeff[y_axis]->ve[1],
      od->vel.fit_coeff[y_axis]->ve[2]);
    if (od->scramble.fit_order == 3) {
      n = strlen(fx);
      sprintf(&fx[n], " + (%.4le) * x**3",
        od->vel.fit_coeff[y_axis]->ve[3]);
    }
    for (i = 0; i < od->scramble.overall_scramblings; ++i) {
      x_value = M_PEEK(od->mal.r2_fit_mat, i, 1);
      y_value = od->vel.fit_vec[y_axis]->ve[i];
      fprintf(plt_dat_out, "%16.6le%16.6le\n",
         x_value, y_value);
    }
    fprintf(plt_cmd_out,
      "unset key\n"
      "unset parametric\n"
      "set title \"Scrambled %s vs r2(yy')\"\n"
      "show title\n"
      "set rmargin 2\n"
      "set xlabel \"r2(yy')\" tc lt -1\n"
      "show xlabel\n"
      "set ylabel \"%s\" tc lt -1\n"
      "show ylabel\n"
      "f(x) = %s\n"
      "plot \"%s\" \\\n\t"
      "using 1:2 with points lt rgb \"blue\" "
      "ps 0.75 pt 5 ,\\\n\t"
      "f(x) lt rgb \"red\" lw 2\n",
      type_of_data,
      y_axis_label,
      fx,
      od->file[PLT_DAT_OUT]->name);
  }
  else if ((type & PLS_PLOT) || (type & PCA_PLOT)) {
    if (type & LOADINGS_PLOT) {
      if (type & PLS_PLOT) {
        mat = od->mal.x_loadings;
        elem = mat->m;
        strcpy(type_of_data, "PLS X loadings");
      }
      else {
        /*
        skip the PC number we are not interested in
        */
        temp_handle = od->file[TEMP_PCA_LOADINGS]->handle;
        if (fseek(temp_handle, sizeof(int), SEEK_SET)) {
          return CANNOT_READ_TEMP_FILE;
        }
        /*
        read the number of loadings present in the file
        */
        actual_len = fread(&elem, sizeof(int), 1, temp_handle);
        if (actual_len != 1) {
          return CANNOT_READ_TEMP_FILE;
        }
        strcpy(type_of_data, "PCA X loadings");
      }
    }
    else if (type & SCORES_PLOT) {
      if (type & PLS_PLOT) {
        mat = od->mal.x_scores;
        elem = mat->m;
        strcpy(type_of_data, "PLS X scores");
      }
      else {
        /*
        skip the PC number we are not interested in
        */
        temp_handle = od->file[TEMP_PCA_SCORES]->handle;
        if (fseek(temp_handle, sizeof(int), SEEK_SET)) {
          return CANNOT_READ_TEMP_FILE;
        }
        /*
        read the number of scores present in the file
        */
        actual_len = fread(&elem, sizeof(int), 1, temp_handle);
        if (actual_len != 1) {
          return CANNOT_READ_TEMP_FILE;
        }
        strcpy(type_of_data, "PCA X scores");
      }
    }
    else if (type & WEIGHTS_PLOT) {
      mat = od->mal.x_weights;
      elem = mat->m;
      strcpy(type_of_data, "PLS X weights");
    }
    abs_max_value = 0.0;
    for (i = 0; i < elem; ++i) {
      if (type & PLS_PLOT) {
        x_value = M_PEEK(mat, i, pc_axis[0] - 1);
        y_value = M_PEEK(mat, i, pc_axis[1] - 1);
      }
      else {
        if (fseek(temp_handle, ((pc_axis[0] - 1) * elem + i) * sizeof(double)
          + 2 * sizeof(int), SEEK_SET)) {
          return CANNOT_READ_TEMP_FILE;
        }  
        actual_len = fread(&x_value, sizeof(double), 1, temp_handle);
        if (actual_len != 1) {
          return CANNOT_READ_TEMP_FILE;
        }
        if (fseek(temp_handle, ((pc_axis[1] - 1) * elem
          + i) * sizeof(double)
          + 2 * sizeof(int), SEEK_SET)) {
          return CANNOT_READ_TEMP_FILE;
        }
        actual_len = fread(&y_value, sizeof(double), 1, temp_handle);
        if (actual_len != 1) {
          return CANNOT_READ_TEMP_FILE;
        }
      }
      if (fabs(x_value) > abs_max_value) {
        abs_max_value = fabs(x_value);
      }
      if (fabs(y_value) > abs_max_value) {
        abs_max_value = fabs(y_value);
      }
      if ((type & SCORES_PLOT) && label) {
        if (!pc_axis[2]) {
          sprintf(with_labels,
            "\"%s\" \\\n\t"
            "using 2:3:1 with labels left "
            "offset 0.5,0.2 lt -1, \\\n\t",
            od->file[PLT_DAT_OUT]->name);
        }
        else {
          sprintf(with_labels,
            "\"%s\" \\\n\t"
            "using 2:3:4:1 with labels left "
            "offset 0.5,0.2 lt -1, \\\n\t",
            od->file[PLT_DAT_OUT]->name);
        }
        n = 0;
        for (j = 0; j < od->object_num; ++j) {
          if (get_object_attr(od, j, ACTIVE_BIT)) {
            if (n == i) {
              break;
            }
            ++n;
          }
        }
      }
      else {
        memset(with_labels, 0, BUF_LEN);
      }
      if (label == NAME_LABEL) {
        fprintf(plt_dat_out, "%-40s",
          od->al.mol_info[j]->object_name);
      }
      else if (label == ID_LABEL) {
        fprintf(plt_dat_out, "%-10d",
          od->al.mol_info[j]->object_id);
      }
      else {
        fprintf(plt_dat_out, "%-10d", j + 1);
      }
      fprintf(plt_dat_out, "%16.6le%16.6le",
         x_value, y_value);
      if (pc_axis[2]) {
        if (type & PLS_PLOT) {
          z_value = M_PEEK(mat, i, pc_axis[2] - 1);
        }
        else {
          if (fseek(temp_handle, ((pc_axis[2] - 1) * elem
            + i) * sizeof(double)
            + 2 * sizeof(int), SEEK_SET)) {
            return CANNOT_READ_TEMP_FILE;
          }
          actual_len = fread(&z_value, sizeof(double), 1, temp_handle);
          if (actual_len != 1) {
            return CANNOT_READ_TEMP_FILE;
          }
        }
        if (fabs(z_value) > abs_max_value) {
          abs_max_value = fabs(z_value);
        }
        fprintf(plt_dat_out, "%16.6le", z_value);
        sprintf(z_dimension,  "set grid\n"
              "set urange [%.6le:%.6le]\n"
              "set ticslevel 0\n"
              "set xtics offset graph 0,graph 0, graph -0.03\n"
              "set ytics offset graph 0,graph 0, graph -0.03\n"
              "set view 60, 60, 1, 1\n"
              "set zlabel \"PC %d\" tc lt -1\n"
              "show zlabel\n"
              "set zrange [%.6le:%.6le]\n"
              "show zrange\n",
              -abs_max_value, abs_max_value,
              pc_axis[2],
              -abs_max_value, abs_max_value);
        sprintf(plot_command,  "splot \"%s\" \\\n\t"
              "using 2:3:4 with points lt rgb \"blue\" "
              "ps 0.75 pt 5, \\\n\t"
              "%s"
              "u,0,0 lt rgb \"red\" lw 2 ,\\\n\t"
              "0,u,0 lt rgb \"red\" lw 2 ,\\\n\t"
              "0,0,u lt rgb \"red\" lw 2\n",
              od->file[PLT_DAT_OUT]->name,
              with_labels);
      }
      else {
        sprintf(z_dimension, "%s", "");
        sprintf(plot_command,  "set grid\n"
              "set trange [%.6le:%.6le]\n"
              "plot \"%s\" \\\n\t"
              "using 2:3 with points lt rgb \"blue\" "
              "ps 0.75 pt 5, \\\n\t"
              "%s"
              "t,0 lt rgb \"red\" lw 2 ,\\\n\t"
              "0,t lt rgb \"red\" lw 2\n",
              -abs_max_value, abs_max_value,
              od->file[PLT_DAT_OUT]->name,
              with_labels);
      }
      fprintf(plt_dat_out, "\n");
    }
    abs_max_value *= 1.1;
    fprintf(plt_cmd_out,
      "unset key\n"
      "set parametric\n"
      "set title \"%s\"\n"
      "show title\n"
      "set rmargin 16\n"
      "set xlabel \"PC %d\" tc lt -1\n"
      "show xlabel\n"
      "set xrange [%.6le:%.6le]\n"
      "show xrange\n"
      "set ylabel \"PC %d\" tc lt -1\n"
      "show ylabel\n"
      "set yrange [%.6le:%.6le]\n"
      "show yrange\n"
      "%s"
      "%s",
      type_of_data,
      pc_axis[0],
      -abs_max_value, abs_max_value,
      pc_axis[1],
      -abs_max_value, abs_max_value,
      z_dimension,
      plot_command);
  }
  if (od->file[TEMP_PCA_LOADINGS]->handle) {
    fclose(od->file[TEMP_PCA_LOADINGS]->handle);
    od->file[TEMP_PCA_LOADINGS]->handle = NULL;
  }
  if (od->file[TEMP_PCA_SCORES]->handle) {
    fclose(od->file[TEMP_PCA_SCORES]->handle);
    od->file[TEMP_PCA_SCORES]->handle = NULL;
  }
  if (plt_dat_out) {
    fclose(plt_dat_out);
    od->file[PLT_DAT_OUT]->handle = NULL;
  }
  if (plt_cmd_out) {
    fclose(plt_cmd_out);
    od->file[PLT_CMD_OUT]->handle = NULL;
  }

  return 0;
}
