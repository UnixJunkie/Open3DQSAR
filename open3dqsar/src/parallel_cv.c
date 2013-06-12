/*

parallel_cv.c

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


int join_thread_files(O3Data *od, ThreadInfo **thread_info)
{
  char buffer[LARGE_BUF_LEN];
  int i;
  int eof;
  int actual_len_read;
  int actual_len_write;


  for (i = 0; i < od->cv.n_threads; ++i) {
    eof = 0;
    fflush(thread_info[i]->temp_pred.handle);
    if (!i) {
      continue;
    }
    rewind(thread_info[i]->temp_pred.handle);
    while (!eof) {
      actual_len_read = fread(buffer, 1, LARGE_BUF_LEN,
        thread_info[i]->temp_pred.handle);
      eof = (actual_len_read < LARGE_BUF_LEN);
      actual_len_write = fwrite(buffer, 1, actual_len_read,
        thread_info[0]->temp_pred.handle);
      if (actual_len_write != actual_len_read) {
        return CANNOT_WRITE_TEMP_FILE;
      }
    }
  }
  fflush(thread_info[0]->temp_pred.handle);
  rewind(thread_info[0]->temp_pred.handle);

  return 0;
}


int parallel_cv(O3Data *od, int x_vars, int suggested_pc_num,
  int model_type, int cv_type, int groups, int runs)
{
  char pred_y_temp_filename[TITLE_LEN];
  char cv_coeff_temp_filename[TITLE_LEN];
  int i;
  int j = 0;
  int x;
  int result;
  int first_parallel_run = 0;
  int run_count = 0;
  int object_num = 0;
  int object_num2 = 0;
  int struct_num = 0;
  int struct_num2 = 0;
  int active_struct_num = 0;
  int conf_num;
  int conf_num2;
  int n_conf = 0;
  double sumweight = 0.0;
  double cum_press;
  double sd_sdep;
  double temp;
  ThreadInfo **ti;
  #ifndef WIN32
  pthread_attr_t thread_attr;
  #endif
  

  result = 0;
  ti = od->mel.thread_info;
  if (((model_type & UVEPLS_CV_MODEL) && (!ti[0]))
    || (!(model_type & UVEPLS_CV_MODEL))) {
    alloc_threads(od);
    first_parallel_run = 1;
  }
  #ifndef WIN32
  /*
  set pthread attributes
  */
  pthread_mutex_init(od->mel.mutex, NULL);
  pthread_attr_init(&thread_attr);
  pthread_attr_setdetachstate(&thread_attr, PTHREAD_CREATE_JOINABLE);
  #else
  if (!(*(od->mel.mutex) = CreateMutex(NULL, FALSE, NULL))) {
    return CANNOT_CREATE_THREAD;
  }
  #endif

  switch (cv_type) {
    case LEAVE_ONE_OUT:
    memset(od->mal.press->base, 0,
      od->mal.press->m * od->mal.press->n * sizeof(double));
    if (((model_type & UVEPLS_CV_MODEL)
      && first_parallel_run) || (!(model_type & UVEPLS_CV_MODEL))) {
      od->cv.cv_thread = (void *)loo_cv_thread;
      run_count = 0;
      while (object_num < od->object_num) {
        struct_num = od->al.mol_info[object_num]->struct_num;
        conf_num = 1;
        for (n_conf = 0, sumweight = 0.0;
          n_conf < conf_num; ++n_conf, ++object_num) {
          if (get_object_attr(od, object_num, ACTIVE_BIT)) {
            sumweight += od->mel.object_weight[object_num];
          }
        }
        if (sumweight > 0.0) {
          /*
          one struct at a time is going to be left out
          */
          od->mel.struct_list[run_count] = struct_num;
          ++run_count;
        }
      }
    }
    break;

    case LEAVE_TWO_OUT:
    memset(od->mal.press->base, 0,
      od->mal.press->m * od->mal.press->n * sizeof(double));
    if (((model_type & UVEPLS_CV_MODEL) && first_parallel_run)
      || (!(model_type & UVEPLS_CV_MODEL))) {
      od->cv.cv_thread = (void *)lto_cv_thread;
      od->mel.struct_list =
        alloc_int_array(od->mel.struct_list,
        od->cv.overall_cv_runs * 2);
      if (!(od->mel.struct_list)) {
        return OUT_OF_MEMORY;
      }
      run_count = 0;
      get_attr_struct_ave(od, 0, ACTIVE_BIT, &active_struct_num, NULL);
      while (object_num < od->object_num) {
        struct_num = od->al.mol_info[object_num]->struct_num;
        conf_num = 1;
        for (n_conf = 0, sumweight = 0.0;
          n_conf < conf_num; ++n_conf, ++object_num) {
          if (get_object_attr(od, object_num, ACTIVE_BIT)) {
            sumweight += od->mel.object_weight[object_num];
          }
        }
        if (sumweight > 0.0) {
          struct_num2 = struct_num + 1;
          object_num2 = object_num;
          while (object_num2 < od->object_num) {
            struct_num2 = od->al.mol_info[object_num2]->struct_num;
            conf_num2 = 1;
            for (n_conf = 0, sumweight = 0.0;
              n_conf < conf_num2; ++n_conf, ++object_num2) {
              if (get_object_attr(od, object_num2, ACTIVE_BIT)) {
                sumweight += od->mel.object_weight[object_num2];
              }
            }
            if (sumweight > 0.0) {
              /*
              two objects at a time are going to be left out
              */
              od->mel.struct_list[run_count * 2] = struct_num;
              od->mel.struct_list[run_count * 2 + 1] = struct_num2;
              ++run_count;
            }
          }
        }
      }
    }
    break;

    case LEAVE_MANY_OUT:
    if (((model_type & UVEPLS_CV_MODEL)
      && first_parallel_run) || (!(model_type & UVEPLS_CV_MODEL))) {
      od->cv.cv_thread = (void *)lmo_cv_thread;
    }
    run_count = runs;
    break;
  }
  /*
  duplicate the main od structure for each ti structure
  also allocate a new array of matrices containing the group_composition_list
  peculiar of each thread
  */
  od->cv.n_threads = ((run_count < od->n_proc) ? run_count : od->n_proc);
  if (((model_type & UVEPLS_CV_MODEL) && first_parallel_run)
    || (!(model_type & UVEPLS_CV_MODEL))) {
    od->cv.n_threads = fill_thread_info(od, run_count);
  }
  for (i = 0; i < od->cv.n_threads; ++i) {
    if (((model_type & UVEPLS_CV_MODEL) && first_parallel_run)
      || (!(model_type & UVEPLS_CV_MODEL))) {
      /*
      ave_press, ave_sdep, sum1, sumweight need to be updated by all threads
      */
      ti[i]->pc_num = suggested_pc_num;
      ti[i]->groups = groups;
      ti[i]->model_type = model_type;
      /*
      it is not necessary to reallocate data structures for thread 0,
      since they are already allocated
      */
      if (i) {
        /*
        this function sets to NULL all data structures which need to be reallocated
        */
        init_pls(&(ti[i]->od));
        /*
        this function allocates new data structures
        */
        result = alloc_pls(&(ti[i]->od), x_vars, suggested_pc_num, model_type);
        if (result) {
          return OUT_OF_MEMORY;
        }
        if (cv_type == LEAVE_MANY_OUT) {
          ti[i]->od.mal.press =
            double_mat_alloc(suggested_pc_num + 1, od->y_vars);
          if (!(ti[i]->od.mal.press)) {
            return OUT_OF_MEMORY;
          }
        }
      }
      if (model_type & (CV_MODEL | SILENT_PLS)) {
        memset(pred_y_temp_filename, 0, TITLE_LEN);
        sprintf(pred_y_temp_filename, "pred_y_t%02d", i + 1);
        if (open_temp_file(&(ti[i]->od),
          &(ti[i]->temp_pred), pred_y_temp_filename)) {
          return CANNOT_WRITE_TEMP_FILE;
        }
        if (!i) {
          od->file[TEMP_PRED]->handle = ti[0]->temp_pred.handle;
          memcpy(od->file[TEMP_PRED]->name, ti[0]->temp_pred.name, BUF_LEN);
        }
      }
      else if ((model_type & UVEPLS_CV_MODEL) && od->uvepls.save_ram) {
        memset(cv_coeff_temp_filename, 0, TITLE_LEN);
        sprintf(cv_coeff_temp_filename, "cv_coeff_t%02d", i + 1);
        if (open_temp_file(&(ti[i]->od),
          &(ti[i]->temp_cv_coeff), cv_coeff_temp_filename)) {
          return CANNOT_WRITE_TEMP_FILE;
        }
        if (!i) {
          od->file[TEMP_CV_COEFF]->handle = ti[0]->temp_cv_coeff.handle;
          memcpy(od->file[TEMP_CV_COEFF]->name, ti[0]->temp_cv_coeff.name, BUF_LEN);
        }
      }
    }
    if (cv_type == LEAVE_MANY_OUT) {
      memset(ti[i]->od.mal.press->base, 0,
        ti[i]->od.mal.press->m * ti[i]->od.mal.press->n * sizeof(double));
    }
    if ((model_type & UVEPLS_CV_MODEL) && od->uvepls.save_ram) {
      rewind(ti[i]->temp_cv_coeff.handle);
    }
    #ifndef WIN32
    od->error_code = pthread_create
      (&(od->thread_id[i]), &thread_attr,
      (void *(*)(void *))(od->cv.cv_thread),
      (void *)ti[i]);
    if (od->error_code) {
      return CANNOT_CREATE_THREAD;
    }
    #else
    od->hThreadArray[i] = CreateThread(NULL, 0,
      (LPTHREAD_START_ROUTINE)(od->cv.cv_thread),
      ti[i], 0, &(od->dwThreadIdArray[i]));
    if (!(od->hThreadArray[i])) {
      return CANNOT_CREATE_THREAD;
    }
    #endif
    
  }
  #ifndef WIN32
  /*
  free the pthread attribute memory
  */
  pthread_attr_destroy(&thread_attr);
  /*
  wait for all threads to have finished
  */
  for (i = 0; i < od->cv.n_threads; ++i) {
    od->error_code = pthread_join(od->thread_id[i],
      &(od->thread_result[i]));
    if (od->error_code) {
      return CANNOT_JOIN_THREAD;
    }
  }
  pthread_mutex_destroy(od->mel.mutex);
  #else
  WaitForMultipleObjects(od->cv.n_threads, od->hThreadArray, TRUE, INFINITE);
  for (i = 0; i < od->cv.n_threads; ++i) {
    CloseHandle(od->hThreadArray[i]);
  }
  CloseHandle(*(od->mel.mutex));
  #endif
  for (i = 0; i < od->cv.n_threads; ++i) {
    if (ti[i]->cannot_write_temp_file) {
      return CANNOT_WRITE_TEMP_FILE;
    }
  }
  od->pc_num = ti[0]->od.pc_num;
  if (model_type & (CV_MODEL | SILENT_PLS)) {
    result = join_thread_files(od, ti);
    if (result) {
      return result;
    }
    if (model_type & CV_MODEL) {
      result = print_pred_values(od);
      if (result) {
        return result;
      }
    }
  }
  if (cv_type == LEAVE_MANY_OUT) {
    if (model_type & CV_MODEL) {
      tee_printf(od, "\nPC%12s%12s%12s\n",
        "SDEP", "SD on SDEP", "Average q2");
      tee_printf(od, "--------------------------------------\n");
    }
    memset(od->vel.ave_sdep->ve, 0,
      od->vel.ave_sdep->size * sizeof(double));
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
        tee_printf(od, "%2d%12.4lf%12.4lf%12.4lf\n", j, od->vel.ave_sdep->ve[j],
          sd_sdep, od->vel.q2->ve[j]);
      }
    }
  }
  else {
    if (model_type & CV_MODEL) {
      tee_printf(od, "\nPC%12s%12s\n",
        "SDEP", "q2");
      tee_printf(od, "--------------------------\n");
    }
    for (i = 0; i <= od->pc_num; ++i) {
      cum_press = 0.0;
      for (x = 0; x < od->y_vars; ++x) {
        cum_press += M_PEEK(od->mal.press, i, x);
      }
      od->vel.ave_sdep->ve[i] = sqrt(cum_press
        / (double)(od->cv.num_predictions));
      od->vel.q2->ve[i] = 1.0 - cum_press / od->cv.tss;
      if (model_type & SCRAMBLE_CV_MODEL) {
        od->vel.secv->ve[i] =
          sqrt(cum_press / (od->y_vars
          * ((double)(od->cv.num_predictions) - i - 1)));
      }
      if (model_type & CV_MODEL) {
        tee_printf(od, "%2d%12.4lf%12.4lf\n", i, od->vel.ave_sdep->ve[i],
          od->vel.q2->ve[i]);
      }
    }
  }
  
  return 0;
}
