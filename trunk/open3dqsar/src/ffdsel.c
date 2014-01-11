/*

ffdsel.c

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
#ifdef WIN32
#include <windows.h>
#endif


#ifndef WIN32
void *ffdsel_thread(void *pointer)
#else
DWORD ffdsel_thread(void *pointer)
#endif
{
  ThreadInfo *ti;
  int i;


  ti = (ThreadInfo *)pointer;
  for (i = ti->start; i <= ti->end; ++i) {
    prepare_design_model(&(ti->od), i);
    if (ti->od.ffdsel.cv_type == EXTERNAL_PREDICTION) {
      trim_mean_center_matrix(&(ti->od), ti->od.mal.large_e_mat,
        &(ti->od.mal.e_mat), &(ti->od.vel.e_mat_ave),
        FFDSEL_FULL_MODEL, ti->od.active_object_num);
      trim_mean_center_matrix(&(ti->od), ti->od.mal.large_f_mat,
        &(ti->od.mal.f_mat), &(ti->od.vel.f_mat_ave),
        FFDSEL_FULL_MODEL, ti->od.active_object_num);
      pls(&(ti->od), ti->pc_num, FFDSEL_FULL_MODEL);
      pred_ext_y_values(&(ti->od), ti->pc_num, FFDSEL_FULL_MODEL);
    }
    else {
      cv(&(ti->od), ti->pc_num, FFDSEL_CV_MODEL,
        ti->od.ffdsel.cv_type, ti->od.ffdsel.groups, ti->od.ffdsel.runs);
    }
    double_vec_sort(ti->od.vel.ave_sdep, ti->od.pel.sdep_rank);
    ti->od.vel.best_sdep->ve[i] = ti->od.vel.ave_sdep->ve[0];
  }
  #ifndef WIN32
  pthread_exit(pointer);
  #else
  return 0;
  #endif
}


void prepare_design_model(O3Data *od, int design_row)
{
  int i;
  int j;
  int x = 0;
  int n = 0;
  int voronoi_count = 0;
  
  
  /*
  loop over the width of the design matrix
  */
  for (i = 0; i < od->ffdsel.design_x; ++i) {
    /*
    if the i-th column does not contain a dummy variable...
    */
    if (od->mel.ffdsel_status[i] != FFDSEL_DUMMY) {
      /*
      if SRD groups are being used, get the number of active variables
      belonging to this Voronoi polyhedron
      */
      n = (od->ffdsel.use_srd_groups ? od->mel.voronoi_active[voronoi_count] : 1);
      /*
      if SRD groups are being used, replicate each column of the design
      matrix for n times, to account for all active variables belonging
      to this Voronoi polyhedron
      */
      for (j = 0; j < n; ++j) {
        /*
        prepare a vector with '1' for the real variables which will be included
        in this model, '0' for those which will be kept out
        */
        od->mel.ffdsel_included[x] = od->cimal.ffd_design_mat->me[design_row][i];
        ++x;
      }
      ++voronoi_count;
    }
  }
}


void write_ffd_design_matrix_col(O3Data *od,
  int first_element, int col, int decimal)
{
  CharMat *ffd_design_mat;
  int i;
  int j;
  int one_or_zero;
  int number;
  int *binary;
  
  
  number = decimal;
  binary = od->mel.binary;
  ffd_design_mat = od->cimal.ffd_design_mat;
  for (i = od->ffdsel.power_of_two - 1; i >= 0; --i) {
    one_or_zero = number / (1 << i);
    number -= (one_or_zero * (1 << i));
    binary[i] = one_or_zero * 2 - 1;
  }
  ffd_design_mat->me[0][col] = (char)(((first_element + 1) / 2) + '0');
  for (i = 0; i < od->ffdsel.power_of_two; ++i) {
    for (j = 0; j < (1 << i); ++j) {
      ffd_design_mat->me[(1 << i) + j][col] =
        (char)(((((((int)ffd_design_mat->me[j][col] -
        (int)'0') * 2 - 1) * binary[i]) + 1) / 2) + '0');
    }
  }
  if (od->ffdsel.fold_over) {
    for (i = ((1 << od->ffdsel.power_of_two) - 1); i > 0; --i) {
      ffd_design_mat->me[i * 2][col] =
        (char)('1' - ffd_design_mat->me[i][col] + '0');
      ffd_design_mat->me[i * 2 - 1][col] = ffd_design_mat->me[i][col];
    }
  }
}


int ffdsel(O3Data *od, int pc_num)
{
  CharMat *ffd_design_mat;
  int dummy_step = 0;
  int column;
  int coeff;
  int incremented;
  int new_column;
  int min_design_y;
  int *binary;
  int *exponents;
  int i;
  int j;
  int k;
  int n_fixed;
  int n_excluded;
  int n_uncertain;
  int overall_n_fixed;
  int overall_n_excluded;
  int overall_n_uncertain;
  int sign;
  int overall_seed_count;
  int voronoi_num;
  int result;
  int which;
  int status;
  int n_threads;
  double ref_value;
  double df;
  double p;
  double q;
  double t;
  double bound;
  ThreadInfo **ti;
  #ifndef WIN32
  pthread_attr_t thread_attr;
  #endif
  
  
  ti = od->mel.thread_info;
    /*
    if FFD selection is based on SRD groups, then
    as many columns as Voronoi polyhedra are needed
    in the design matrix; if it is based on normal
    variables as many columns as active x variables
    are needed
    */
  od->ffdsel.design_x = (od->ffdsel.use_srd_groups
    ? od->voronoi_num : od->overall_active_x_vars);
  od->ffdsel.n_dummies = 0;
  /*
  calculate the number of dummies according
  to the percentage chosen by the user
  */
  if (od->ffdsel.percent_dummies) {
    /*
    calculate every how many columns a dummy variable
    has to be placed
    */
    dummy_step = (int)safe_rint((double)(100 - od->ffdsel.percent_dummies)
      / (double)(od->ffdsel.percent_dummies));
    /*
    adjust the number of dummies accordingly
    */
    od->ffdsel.n_dummies = (int)((double)od->ffdsel.design_x
      / (double)dummy_step) + 1;
    ++dummy_step;
  }
  /*
  set the SEL_INCLUDED_BIT
  */
  od->ffdsel.ffdsel_included_vars = set_sel_included_bit(od, od->ffdsel.use_srd_groups);
  /*
  the minimum number of combinations is:
  (number SRD groups/variables + number of dummies) * combination/variable ratio
  then this number is adjusted to the nearest power of two p, provided that p does not
  exceed the number of SRD/groups variables, or duplicate combinations would arise
  */  
  min_design_y = (int)safe_rint(od->ffdsel.combination_variable_ratio
    * (double)(od->ffdsel.design_x + od->ffdsel.n_dummies));
  od->ffdsel.power_of_two = 0;
  od->ffdsel.design_y = 1;
  while ((od->ffdsel.design_y < min_design_y)
    && (od->ffdsel.power_of_two < od->ffdsel.design_x)) {
    ++(od->ffdsel.power_of_two);
    od->ffdsel.design_y <<= 1;
  }
  od->ffdsel.design_x += od->ffdsel.n_dummies;
  /*
  if the number of columns is exactly equal to the number of rows in the design matrix,
  then double the number of rows (combinations)
  */
  if (od->ffdsel.design_x == od->ffdsel.design_y) {
    ++(od->ffdsel.power_of_two);
    od->ffdsel.design_y <<= 1;
  }
  /*
  In this array the status of the different SRD groups/variables is stored
  possible values are FFDSEL_DUMMY, FFDSEL_UNCERTAIN, FFDSEL_FIXED, FFDSEL_EXCLUDED
  */
  od->mel.ffdsel_status = (unsigned char *)realloc
    (od->mel.ffdsel_status, od->ffdsel.design_x);
  if (!(od->mel.ffdsel_status)) {
    return OUT_OF_MEMORY;
  }
  /*
  at the beginning all SRD groups/variables are initialized to FFDSEL_UNCERTAIN
  */
  memset(od->mel.ffdsel_status,
    FFDSEL_UNCERTAIN, od->ffdsel.design_x);
  /*
  this array is used by the write_ffd_design_matrix_col function
  */
  od->mel.binary = alloc_int_array(od->mel.binary, od->ffdsel.design_y);
  if (!(od->mel.binary)) {
    return OUT_OF_MEMORY;
  }
  binary = od->mel.binary;
  /*
  if a fold-over design matrix was requested by the user, then
  double the number of combinations and subtract 1, since it would
  be a combination of no variables (all -1 columns)
  */
  od->ffdsel.design_y = od->ffdsel.design_y
    * (od->ffdsel.fold_over + 1) - od->ffdsel.fold_over;
  /*
  allocate the design matrix
  */
  od->cimal.ffd_design_mat = alloc_char_matrix(od->cimal.ffd_design_mat,
    od->ffdsel.design_y, od->ffdsel.design_x);
  if (!(od->cimal.ffd_design_mat)) {
    return OUT_OF_MEMORY;
  }
  ffd_design_mat = od->cimal.ffd_design_mat;
  /*
  allocate some exponent arrays used to generate the design matrix columns
  */
  od->mel.exponents = alloc_int_array(od->mel.exponents,
    od->ffdsel.power_of_two - 2);
  if (!(od->mel.exponents)) {
    return OUT_OF_MEMORY;
  }
  exponents = od->mel.exponents;

  if (od->ffdsel.n_dummies) {
    /*
    assign dummies positions in the design matrix
    one every dummy_step steps
    */
    for (i = 0; i < od->ffdsel.design_x; i += dummy_step) {
      od->mel.ffdsel_status[i] = FFDSEL_DUMMY;
    }
  }
  for (i = 0; i < (od->ffdsel.power_of_two - 2); ++i) {
    exponents[i] = -1;
  }
  column = 0;
  incremented = 0;
  /*
  write down the first p columns according to a full factorial design
  */
  for (i = 0; i < od->ffdsel.power_of_two; ++i) {
    coeff = (1 << od->ffdsel.power_of_two) - (1 << i) - 1;
    write_ffd_design_matrix_col(od, 1, column, coeff);
    ++column;
  }
  /*
  now write the remaining columns as in the design matrix used in
  Baroni M., Clementi S., Cruciani G., Costantino G., Riganelli D.
  J. Chemometrics 1992, 6, 347-356.
  */
  while (1) {
    coeff = 0;
    for (i = 0; i < (od->ffdsel.power_of_two - 2); ++i) {
      if (exponents[i] >= 0) {
        coeff += (1 << exponents[i]);
      }
    }
    write_ffd_design_matrix_col(od, 1, column, coeff);
    ++column;
    if (column == od->ffdsel.design_x) {
      break;
    }
    new_column = 1;
    for (i = incremented; i>= 0; --i) {
      if (exponents[i] < (od->ffdsel.power_of_two - 1 - incremented + i)) {
        ++exponents[i];
        new_column = 0;
        for (j = i + 1; j <= incremented; ++j) {
          exponents[j] = exponents[j - 1] + 1;
        }
        break;
      }
    }
    if (new_column) {
      ++incremented;
      for (i = 0; i <= incremented; ++i) {
        exponents[i] = i;
      }
    }
  }
  /*
  allocate a large variable matrix including
  all active variables and all objects
  */
  od->mal.large_e_mat =
    double_mat_resize(od->mal.large_e_mat,
    od->active_object_num + od->ext_pred_object_num,
    od->ffdsel.ffdsel_included_vars);
  if (!(od->mal.large_e_mat)) {
    return OUT_OF_MEMORY;
  }
  /*
  fill up both x and y matrices
  */
  result = fill_x_matrix(od,
    FFDSEL_FULL_MODEL, od->ffdsel.use_srd_groups);
  if (result) {
    return result;
  }
  result = fill_y_matrix(od);
  if (result) {
    return result;
  }
  /*
  print out some information about the design matrix
  */
  tee_printf(od, "Design matrix: %d combinations x %d ",
    od->ffdsel.design_y, od->ffdsel.design_x);
  if (od->ffdsel.use_srd_groups) {
    tee_printf(od, "groups of ");
  }
  tee_printf(od, "variables\n");
  if (od->ffdsel.n_dummies) {
    tee_printf(od, "(%d real, %d dummy)\n", od->ffdsel.design_x
       - od->ffdsel.n_dummies, od->ffdsel.n_dummies);
  }
  if (od->ffdsel.use_srd_groups) {
    tee_printf(od, "Active variables not in group zero: %d\n",
      od->ffdsel.ffdsel_included_vars);
  }
  tee_printf(od, "\n");
  tee_flush(od);
  /*
  allocate a vector to store best SDEP values
  for the models derived by each combination
  */
  od->vel.best_sdep = double_vec_resize
    (od->vel.best_sdep, od->ffdsel.design_y);
  if (!(od->vel.best_sdep)) {
    return OUT_OF_MEMORY;
  }
  /*
  allocate a vector to store effects of each SRD group/variable
  */
  od->vel.effect = double_vec_resize
    (od->vel.effect, od->ffdsel.design_x);
  if (!(od->vel.effect)) {
    return OUT_OF_MEMORY;
  }
  result = alloc_cv_sdep(od, pc_num, od->ffdsel.runs);
  if (result) {
    return OUT_OF_MEMORY;
  }
  result = alloc_pls(od, od->ffdsel.ffdsel_included_vars,
    pc_num, FFDSEL_CV_MODEL);
  if (result) {
    return OUT_OF_MEMORY;
  }
  set_random_seed(od, od->random_seed);
  result = prepare_cv(od, pc_num, od->ffdsel.cv_type,
    od->ffdsel.groups, od->ffdsel.runs);
  if (result) {
    return result;
  }
  /*
  allocate and pre-calculate averages for all CV models
  */
  result = alloc_average_mat(od, FFDSEL_CV_MODEL,
    od->ffdsel.cv_type, od->ffdsel.groups, od->ffdsel.runs);
  if (result) {
    return OUT_OF_MEMORY;
  }
  /*
  allocate structures which will be passed to each computational thread
  */
  alloc_threads(od);
  #ifndef WIN32
  /*
  set pthread attributes
  */
  pthread_attr_init(&thread_attr);
  pthread_attr_setdetachstate(&thread_attr, PTHREAD_CREATE_JOINABLE);
  #endif
  /*
  duplicate the main od structure for each ti structure
  also allocate a new array describing the SRD group/variable composition
  of each model for each thread
  */
  n_threads = fill_thread_info(od, od->ffdsel.design_y);
  for (i = 0; i < n_threads; ++i) {
    ti[i]->od.mel.ffdsel_included =
      malloc(od->ffdsel.ffdsel_included_vars);
    if (!(ti[i]->od.mel.ffdsel_included)) {
      return OUT_OF_MEMORY;
    }
    /*
    it is not necessary to reallocate data structures for thread 0,
    since they are already allocated
    */
    if (i) {
      init_cv_sdep(&(ti[i]->od));
      result = alloc_cv_sdep(&(ti[i]->od), pc_num, od->ffdsel.runs);
      if (result) {
        return OUT_OF_MEMORY;
      }
      /*
      this function sets to NULL all data structures
      which need to be reallocated
      */
      init_pls(&(ti[i]->od));
      /*
      this function allocates new data structures
      */
      ti[i]->od.mal.press =
        double_mat_alloc(pc_num + 1, od->y_vars);
      if (!(ti[i]->od.mal.press)) {
        return OUT_OF_MEMORY;
      }
      memset(ti[i]->od.mal.press->base, 0,
        ti[i]->od.mal.press->m
        * ti[i]->od.mal.press->n
        * sizeof(double));
      result = alloc_pls(&(ti[i]->od),
        od->ffdsel.ffdsel_included_vars,
        pc_num, FFDSEL_CV_MODEL);
      if (result) {
        return OUT_OF_MEMORY;
      }
    }
    ti[i]->pc_num = pc_num;
    /*
    create the i-th thread
    */
    #ifndef WIN32
    od->error_code = pthread_create(&(od->thread_id[i]), &thread_attr,
      ffdsel_thread, (void *)ti[i]);
    if (od->error_code) {
      return CANNOT_CREATE_THREAD;
    }
    #else
    od->hThreadArray[i] = CreateThread(NULL, 0,
      (LPTHREAD_START_ROUTINE)ffdsel_thread,
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
  for (i = 0; i < n_threads; ++i) {
    od->error_code = pthread_join(od->thread_id[i],
      &(od->thread_result[i]));
    if (od->error_code) {
      return CANNOT_JOIN_THREAD;
    }
  }
  #else
  WaitForMultipleObjects(n_threads, od->hThreadArray, TRUE, INFINITE);
  for (i = 0; i < n_threads; ++i) {
    CloseHandle(od->hThreadArray[i]);
  }
  #endif
  for (i = 0; i < n_threads; ++i) {
    /*
    free all duplicate structures
    */
    if (ti[i]->od.mel.ffdsel_included) {
      free(ti[i]->od.mel.ffdsel_included);
      ti[i]->od.mel.ffdsel_included = NULL;
    }
    if (i) {
      free_cv_sdep(&(ti[i]->od));
      if (ti[i]->od.mal.press) {
        double_mat_free(ti[i]->od.mal.press);
        ti[i]->od.mal.press = NULL;
      }
      free_pls(&(ti[i]->od));
    }
  }
  if (od->ffdsel.cv_type == LEAVE_MANY_OUT) {
    free_cv_groups(od, od->ffdsel.runs);
  }
  free_threads(od);

  /*
  if the user requested extensive printout of all SDEPs
  */
  if (od->ffdsel.print_sdep) {
    tee_printf(od,
      "-------------------------\n"
      "SDEP of each combination:\n"
      "-------------------------\n");
    for (i = 0; i < od->ffdsel.design_y; ++i) {
      tee_printf(od, "%12d%12.4lf\n", i + 1,
        od->vel.best_sdep->ve[i]);
    }
    tee_printf(od, "\n");
  }
  
  /*
  calculate effects on SDEP of each variable/dummy
  using Yates' algorithm
  */
  for (j = 0; j < od->ffdsel.design_x; ++j) {
    od->vel.effect->ve[j] =
      od->vel.best_sdep->ve[0];
    for (i = 1; i < od->ffdsel.design_y; ++i) {
      sign = (int)(ffd_design_mat->me[i][j] - '0') * 2 - 1;
      od->vel.effect->ve[j] += (double)sign
        * od->vel.best_sdep->ve[i];
    }
    od->vel.effect->ve[j]
      /= (double)(od->ffdsel.design_x);
  }

  /*
  if the user requested extensive printout of all effects
  */
  if (od->ffdsel.print_effect) {
    tee_printf(od, "-----------------------");
    if (od->ffdsel.use_srd_groups) {
      tee_printf(od, "----------");
    }
    tee_printf(od, "\n"
      "Effect of each ");
    if (od->ffdsel.use_srd_groups) {
      tee_printf(od, "group of ");
    }
    tee_printf(od, "variable");
    if (od->ffdsel.use_srd_groups) {
      tee_printf(od, "s");
    }
    tee_printf(od, "\n"
      "-----------------------");
    if (od->ffdsel.use_srd_groups) {
      tee_printf(od, "----------");
    }
    tee_printf(od, "\n");
    for (i = 0; i < od->ffdsel.design_x; ++i) {
      tee_printf(od, "%12d%12.4lf\n", i + 1,
      od->vel.effect->ve[i]);
    }
    tee_printf(od, "\n");
  }
  
  /*
  compute the mean effect of dummies according to
  Baroni M., Costantino G., Cruciani G., Riganelli D., Valigi R., Clementi S.
  Quant. Struct.-Act. Relat. 1993, 12, 9-20
  */
  df = (double)(od->ffdsel.n_dummies);
  which = 2;
  p = (double)(od->ffdsel.confidence_level) / (double)100;
  q = (double)1 - p;
  cdft(&which, &p, &q, &t, &df, &status, &bound);
  if (status) {
    return STUDENT_T_OUT_OF_BOUNDS;
  }
  ref_value = (double)0;
  for (i = 0; i < od->ffdsel.design_x; ++i) {
    if (od->mel.ffdsel_status[i] == FFDSEL_DUMMY) {
      ref_value += square
        (od->vel.effect->ve[i]);
    }

  }
  ref_value = sqrt(ref_value / od->ffdsel.n_dummies) * t;
  /*
  now compute the number of excluded, uncertain and fixed variables
  */
  for (i = 0; i < od->ffdsel.design_x; ++i) {
    if (od->mel.ffdsel_status[i] == FFDSEL_DUMMY) {
      continue;
    }
    if (((fabs(od->vel.effect->ve[i]) - ref_value) > -ALMOST_ZERO)
      && (od->vel.effect->ve[i] > -ALMOST_ZERO)) {
      od->mel.ffdsel_status[i] = FFDSEL_EXCLUDED;
    }
    else if (fabs(od->vel.effect->ve[i]) < ref_value) {
      od->mel.ffdsel_status[i] = FFDSEL_UNCERTAIN;
    }
    else {
      od->mel.ffdsel_status[i] = FFDSEL_FIXED;
    }
  }
  /*
  print out the results of the variable selection
  */
  tee_printf(od,
    "%-5s%16s%16s%16s\n"
    "-----------------------------------------------------\n",
    "Field", "Excluded", "Uncertain", "Fixed");
  j = 0;
  /*
  set the FFDSEL_BIT to 1 on all variables eligible for removal
  with the REMOVE tool
  */
  overall_n_fixed = 0;
  overall_n_excluded = 0;
  overall_n_uncertain = 0;
  if (!(od->ffdsel.use_srd_groups)) {
    for (i = 0; i < od->field_num; ++i) {
      if (get_field_attr(od, i, ACTIVE_BIT)) {
        n_fixed = 0;
        n_excluded = 0;
        n_uncertain = 0;
        for (k = 0; k < od->x_vars; ++k) {
          if (get_x_var_attr(od, i, k, SEL_INCLUDED_BIT)) {
            while (od->mel.ffdsel_status[j] == FFDSEL_DUMMY) {
              ++j;
            }
            if (od->mel.ffdsel_status[j] == FFDSEL_EXCLUDED) {
              set_x_var_attr(od, i, k, FFDSEL_BIT, 1);
              ++n_excluded;
            }
            else if (od->mel.ffdsel_status[j] == FFDSEL_UNCERTAIN) {
              if (!(od->ffdsel.retain_uncertain)) {
                set_x_var_attr(od, i, k, FFDSEL_BIT, 1);
              }
              ++n_uncertain;
            }
            else {
              ++n_fixed;
            }
            ++j;
          }
        }
        tee_printf(od,
          "%5d%16d%16d%16d\n",
          i + 1, n_excluded, n_uncertain, n_fixed);
        overall_n_fixed += n_fixed;
        overall_n_excluded += n_excluded;
        overall_n_uncertain += n_uncertain;
      }
    }
  }
  else {
    overall_seed_count = 0;
    voronoi_num = 0;
    i = 0;
    while (i < od->field_num) {
      if (get_field_attr(od, i, ACTIVE_BIT)) {
        n_fixed = 0;
        n_excluded = 0;
        n_uncertain = 0;
        while ((voronoi_num - overall_seed_count)
          < od->mel.seed_count[i]) {
          while (od->mel.ffdsel_status[j] == FFDSEL_DUMMY) {
            ++j;
          }
          for (k = 0; k < od->x_vars; ++k) {
            if (get_x_var_attr(od, i, k, SEL_INCLUDED_BIT)
              && (get_voronoi_buf(od, i, k) == voronoi_num)) {
              if (od->mel.ffdsel_status[j] == FFDSEL_EXCLUDED) {
                set_x_var_attr(od, i, k, FFDSEL_BIT, 1);
                ++n_excluded;
              }
              else if (od->mel.ffdsel_status[j] == FFDSEL_UNCERTAIN) {
                if (!(od->ffdsel.retain_uncertain)) {
                  set_x_var_attr(od, i, k, FFDSEL_BIT, 1);
                }
                ++n_uncertain;
              }
              else {
                ++n_fixed;
              }
            }
          }
          ++voronoi_num;
          ++j;
        }
        tee_printf(od,
          "%5d%16d%16d%16d\n",
          i + 1, n_excluded, n_uncertain, n_fixed);
        overall_n_fixed += n_fixed;
        overall_n_excluded += n_excluded;
        overall_n_uncertain += n_uncertain;
      }
      overall_seed_count += od->mel.seed_count[i];
      ++i;
    }
  }
  tee_printf(od,
    "%-5s%16d%16d%16d\n"
    "-----------------------------------------------------\n\n",
    "Total", overall_n_excluded, overall_n_uncertain, overall_n_fixed);
  if (od->mal.large_e_mat_ave) {
    double_mat_free(od->mal.large_e_mat_ave);
    od->mal.large_e_mat_ave = NULL;
  }
  if (od->mal.large_f_mat_ave) {
    double_mat_free(od->mal.large_f_mat_ave);
    od->mal.large_f_mat_ave = NULL;
  }
  if (od->cimal.ffd_design_mat) {
    free_char_matrix(od->cimal.ffd_design_mat);
    od->cimal.ffd_design_mat = NULL;
  }
  
  return 0;
}
