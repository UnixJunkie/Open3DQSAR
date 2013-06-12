/*

srd.c

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


int get_voronoi_buf(O3Data *od, int field_num, int x_var)
{
  return od->cimal.voronoi_buf->me[field_num][x_var];
}


void set_voronoi_buf(O3Data *od, int field_num, int x_var, int voronoi_num)
{
  od->cimal.voronoi_buf->me[field_num][x_var] = voronoi_num;
}


int calc_p_vectors(O3Data *od, int field_num, int seed_num)
{
  int i;
  int j;
  int n;
  int object_num;
  int overall_seed_count;
  int cum_value_found;
  int pos_value_found;
  int neg_value_found;
  int result;
  double value;
  
  
  object_num = 0;
  overall_seed_count = 0;
  for (i = 0; i < field_num; ++i) {
    if (get_field_attr(od, i, ACTIVE_BIT)) {
      overall_seed_count += od->mel.seed_count[i];
    }
  }
  for (j = 0; j < od->object_num; ++j) {
    if (!get_object_attr(od, j, ACTIVE_BIT)) {
      continue;
    }
    M_POKE(od->mal.cum_ave, seed_num, object_num, 0.0);
    M_POKE(od->mal.pos_ave, seed_num, object_num, 0.0);
    M_POKE(od->mal.neg_ave, seed_num, object_num, 0.0);
    cum_value_found = 0;
    pos_value_found = 0;
    neg_value_found = 0;
    for (n = 0; n < od->mel.voronoi_fill[seed_num + overall_seed_count]; ++n) {
      if (!get_x_var_attr(od, field_num, od->al.voronoi_composition
        [seed_num + overall_seed_count][n], ACTIVE_BIT)) {
        continue;
      }
      result = get_x_value(od, field_num, j,
        od->al.voronoi_composition[seed_num + overall_seed_count][n], &value, 0);
      if (result) {
        return result;
      }
      if (!MISSING(value)) {
        ++cum_value_found;
        result = get_x_value(od, field_num, j,
          od->al.voronoi_composition[seed_num + overall_seed_count][n],
          &value, CHECK_IF_ACTIVE_BIT | CUTOFF_BIT | WEIGHT_BIT);
        if (result) {
          return result;
        }
        M_POKE(od->mal.cum_ave, seed_num, object_num,
          M_PEEK(od->mal.cum_ave, seed_num, object_num) + value);
        if (value > ALMOST_ZERO) {
          ++pos_value_found;
          M_POKE(od->mal.pos_ave, seed_num, object_num,
            M_PEEK(od->mal.pos_ave, seed_num, object_num) + value);
        }
        if (value < -ALMOST_ZERO) {
          ++neg_value_found;
          M_POKE(od->mal.neg_ave, seed_num, object_num,
            M_PEEK(od->mal.neg_ave, seed_num, object_num) + value);
        }
      }
    }
    M_POKE(od->mal.cum_ave, seed_num, object_num,
      M_PEEK(od->mal.cum_ave, seed_num, object_num) / (double)cum_value_found);
    if (pos_value_found) {
      M_POKE(od->mal.pos_ave, seed_num, object_num,
        M_PEEK(od->mal.pos_ave, seed_num, object_num) / (double)pos_value_found);
    }
    else {
      /*
      a negative value where only a positive value makes sense
      actually indicates a missing value
      */
      M_POKE(od->mal.pos_ave, seed_num, object_num, -1.0);
    }
    if (neg_value_found) {
      M_POKE(od->mal.neg_ave, seed_num, object_num,
        M_PEEK(od->mal.neg_ave, seed_num, object_num) / (double)neg_value_found);
    }
    else {
      /*
      a positive value where only a negative value makes sense
      actually indicates a missing value
      */
      M_POKE(od->mal.neg_ave, seed_num, object_num, 1.0);
    }
    ++object_num;
  }
  
  return 0;
}


double pearson_r(DoubleMat *mat, int *x,
  int check_missing, double missing)
{
  int i;
  int initial;
  int k;
  int n;
  int actual_n;
  double sum_x[2];
  double sum_sq_x[2];
  double sum_x0x1;
  double mean_x[2];
  double sweep;
  double delta_x[2];
  double pop_sd_x[2];
  double cov_x0x1;
  double den;
  double correlation;
  
  n = mat->n;
  for (k = 0; k <= 1; ++k) {
    sum_x[k] = 0.0;
    sum_sq_x[k] = 0.0;
    actual_n = 0;
    initial = 0;
    if (check_missing) {
      while ((int)M_PEEK(mat, x[k], initial) == (int)missing) {
        ++initial;
        if (initial == n) {
          return 1.0;
        }
      }
    }
    mean_x[k] = M_PEEK(mat, x[k], initial);
    for (i = 0; i < n; ++i) {
      if (check_missing) {
        if ((int)M_PEEK(mat, x[k], i) == (int)missing) {
          continue;
        }
      }
      ++actual_n;
      sum_x[k] += M_PEEK(mat, x[k], i);
    }
  }
  sum_x0x1 = 0.0;
  actual_n = 1;
  for (i = initial + 1; i < n; ++i) {
    if (check_missing) {
      if ((int)M_PEEK(mat, x[0], i) == (int)missing) {
        continue;
      }
    }
    ++actual_n;
    sweep = (double)(actual_n - 1) / (double)actual_n;
    for (k = 0; k <= 1; ++k) {
      delta_x[k] = M_PEEK(mat, x[k], i) - mean_x[k];
      sum_sq_x[k] += (delta_x[k] * delta_x[k] * sweep);
      mean_x[k] += (delta_x[k] / (double)actual_n);
    }
    sum_x0x1 += (delta_x[0] * delta_x[1] * sweep);
  }
  for (k = 0; k <= 1; ++k) {
    pop_sd_x[k] = sqrt(sum_sq_x[k] / (double)actual_n);
  }
  cov_x0x1 = sum_x0x1 / (double)actual_n;
  den = pop_sd_x[0] * pop_sd_x[1];
  if (den < ALMOST_ZERO) {
    correlation = 9.9999;
  }
  else {
    correlation = cov_x0x1 / den;
  }

  return correlation;
}


void pseudo_seed_coord(O3Data *od, int field_num, int *seed)
{
  int i;
  int k;
  int n;
  int fill;
  int n_active;
  int overall_seed_count;
  double cart[3];

  
  memset(cart, 0, 3 * sizeof(double));
  for (i = 0, overall_seed_count = 0; i < field_num; ++i) {
    overall_seed_count += od->mel.seed_count[i];
  }
  for (i = 0; i < 3; ++i) {
    for (k = 0, fill = 0; k <= 1; ++k) {
      for (n = 0, n_active = 0; n < od->mel.voronoi_fill[overall_seed_count + seed[k]]; ++n) {
        if ((!get_x_var_attr(od, field_num, od->al.voronoi_composition[overall_seed_count + seed[k]][n], SEED_BIT))
          && get_x_var_attr(od, field_num, od->al.voronoi_composition[overall_seed_count + seed[k]][n], ACTIVE_BIT)) {
          ++n_active;
        }
      }
      cart[i] += (od->al.seed_coord[seed[k]]->cart[i]
        * (n_active + 1));
      fill += (n_active + 1);
    }
    cart[i] /= (double)fill;
  }
  memcpy(od->al.seed_coord[seed[1]]->cart, cart, 3 * sizeof(double));
}


int srd(O3Data *od, int pc_num, int seed_num, int type,
  int collapse, double critical_distance, double collapse_distance)
{
  int i;
  int j;
  int k;
  int n;
  int field_num;
  int result;
  int x_var;
  int seed;
  int factors;
  int candidates;
  int support_points;
  int actual_len;
  int nearest;
  int max_elem;
  int elem;
  int collapse_checked_len;
  int is_checked;
  int attempt;
  int overall_seed_count;
  int overall_group_zero;
  int overall_seed_count_before_collapse;
  int voronoi_population;
  int different_pattern;
  int candidate_seed[2];
  int **voronoi_composition;
  int *new_voronoi_composition;
  IntMat *collapse_checked = NULL;
  IntMat *seed_list;
  int *seed_count;
  int *seed_count_before_collapse;
  int *voronoi_fill;
  int *voronoi_active;
  int *group_zero;
  double distance;
  double min_distance;
  double cum_r;
  double pos_r;
  double neg_r;
  VarCoord var_coord;
  VarCoord current_seed_coord;
  IntPerm *x_vars_list;
  FILE *x_matrix_handle;
  
  /*
  start extracting seed_num seeds with d-optimal design
  */
  factors = pc_num;
  candidates = od->overall_active_x_vars;
  support_points = candidates - seed_num;
  
  if ((result = d_optimal(od, factors, seed_num, type))) {
    return result;
  }
  
  x_matrix_handle = fopen(od->file[TEMP_X_MATRIX]->name, "rb");
  if (!x_matrix_handle) {
    return CANNOT_READ_TEMP_FILE;
  }
  od->file[TEMP_X_MATRIX]->handle = x_matrix_handle;
  x_vars_list = od->pel.dxn_perm;
  free_char_matrix((CharMat *)(od->cimal.seed_list));
  od->cimal.seed_list = NULL;
  od->cimal.seed_list = alloc_int_matrix
    (od->cimal.seed_list, od->field_num, seed_num);
  if (!(od->cimal.seed_list)) {
    return OUT_OF_MEMORY;
  }
  seed_list = od->cimal.seed_list;
  /*
  free eventual arrays previously allocated,
  then allocate new ones
  */
  free_array(od->al.voronoi_composition);
  od->al.voronoi_composition = NULL;
  if (alloc_voronoi(od, seed_num)) {
    return OUT_OF_MEMORY;
  }
  od->voronoi_num = seed_num;
  seed_count = od->mel.seed_count;
  seed_count_before_collapse = od->mel.seed_count_before_collapse;
  group_zero = od->mel.group_zero;
  voronoi_fill = od->mel.voronoi_fill;
  voronoi_active = od->mel.voronoi_active;
  voronoi_composition = od->al.voronoi_composition;
  /*
  Read from the TEMP_X_MATRIX file field_num
  and x_var numbers corresponding to active x_vars
  */
  n = 0;
  while (1) {
    fread(&field_num, sizeof(int), 1, x_matrix_handle);
    actual_len = fread(&x_var, sizeof(int), 1, x_matrix_handle);
    if (!actual_len) {
      break;
    }
    /*
    find the field and the x_var corresponding
    to each selected seed and put its progressive number in
    an int matrix with field_num rows and seed_num columns
    */
    for (i = 0; i < seed_num; ++i) {
      if (n == x_vars_list->pe[i + support_points]) {
        seed_list->me[field_num][seed_count[field_num]] = x_var;
        ++seed_count[field_num];
        break;
      }
    }
    ++n;
  }
  tee_printf(od, "%5s%12s%20s\n", "Field", "Groups", "Vars in group 0");
  tee_printf(od, "-------------------------------------\n");
  /*
  For each field build Voronoi polyhedra getting the
  nearest variables around each seed;
  then store the composition of each polyhedron in a
  temporary file
  */
  overall_seed_count = 0;
  overall_group_zero = 0;
  for (i = 0; i < od->field_num; ++i) {
    if (!get_field_attr(od, i, ACTIVE_BIT)) {
      continue;
    }
    /*
    sort x_var numbers corresponding to seeds in seed_list->me[i]
    */
    qsort(seed_list->me[i], seed_count[i], sizeof(int), compare_integers);
    if (open_temp_file(od, od->file[TEMP_VORONOI], "voronoi")) {
      return CANNOT_WRITE_TEMP_FILE;
    }
    /*
    initialize SEED_BIT and GROUP_BIT in x_var_attr:
    turn them off for all vars, then turn SEED_BIT on
    for seeds, since they must not be assigned to any polyhedron
    */
    for (j = 0; j < od->x_vars; ++j) {
      set_x_var_attr(od, i, j, SEED_BIT | GROUP_BIT, 0);
    }
    for (j = 0; j < seed_count[i]; ++j) {
      set_x_var_attr(od, i, seed_list->me[i][j], SEED_BIT, 1);
    }
    for (k = 0; k < od->x_vars; ++k) {
      /*
      if the current x_var is also a seed, skip it
      since it will have 0.00000 distance
      */
      if (get_x_var_attr(od, i, k, SEED_BIT)) {
        continue;
      }
      /*
      get Cartesian coordinates of the k-th variable
      */
      var_to_xyz(od, k, &var_coord);
      /*
      set min_distance to be higher than critical_distance
      */
      min_distance = critical_distance + 1.0;
      /*
      loop over all seeds for the current field
      */
      for (j = 0; j < seed_count[i]; ++j) {
        var_to_xyz(od, seed_list->me[i][j], &current_seed_coord);
        /*
        get distance from this seed to the current x_var
        */
        distance = sqrt(squared_euclidean_distance
          (var_coord.cart, current_seed_coord.cart)) / od->grid.step[0];
        /*
        if the distance is the lowest until now record the
        distance and the overall seed number
        */
        if (distance < min_distance) {
          min_distance = distance;
          nearest = j + overall_seed_count;
        }
      }
      /*
      if min_distance is lower than critical_distance
      then assign this x_var to the Voronoi polyhedron
      identified by this seed
      */
      if ((min_distance - critical_distance) < ALMOST_ZERO) {
        /*
        format of this binary file is:
        - x_var number
        - seed number
        */
        fwrite(&k, sizeof(int), 1,
          od->file[TEMP_VORONOI]->handle);
        actual_len = fwrite(&nearest, sizeof(int), 1,
          od->file[TEMP_VORONOI]->handle);
        if (actual_len != 1) {
          return CANNOT_WRITE_TEMP_FILE;
        }
        /*
        mark this variable as already assigned to a polyhedron
        */
        set_x_var_attr(od, i, k, GROUP_BIT, 1);
      }
    }
    /*
    loop over all seeds for the current field
    */
    for (j = 0; j < seed_count[i]; ++j) {
      /*
      count how many x_vars are there in each polyhedron
      */
      rewind(od->file[TEMP_VORONOI]->handle);
      /*
      starting population is 1 because we count also the seed
      */
      voronoi_population = 1;
      while (1) {
        fread(&x_var, sizeof(int), 1,
          od->file[TEMP_VORONOI]->handle);
        actual_len = fread(&seed, sizeof(int), 1,
          od->file[TEMP_VORONOI]->handle);
        if (actual_len != 1) {
          break;
        }
        if (seed == (j + overall_seed_count)) {
          ++voronoi_population;
        }
      }
      /*
      for each seed allocate an array of integers to store
      x_var numbers belonging to the Voronoi polyhedron
      */
      voronoi_composition[j + overall_seed_count] =
        alloc_int_array(NULL, voronoi_population);
      if (!voronoi_composition[j + overall_seed_count]) {
        return OUT_OF_MEMORY;
      }
      /*
      initialize the population of each group by setting
      the x_var corresponding to the respective seed
      as the first element
      */
      voronoi_composition[j + overall_seed_count][0] = seed_list->me[i][j];
      voronoi_fill[j + overall_seed_count] = 1;
    }
    /*
    now re-read the file once and assign x_vars to each seed
    */
    rewind(od->file[TEMP_VORONOI]->handle);
    while (1) {
      fread(&x_var, sizeof(int), 1,
        od->file[TEMP_VORONOI]->handle);
      actual_len = fread(&seed, sizeof(int), 1,
        od->file[TEMP_VORONOI]->handle);
      if (actual_len != 1) {
        break;
      }
      voronoi_composition[seed][voronoi_fill[seed]] = x_var;
      ++voronoi_fill[seed];
    }
    group_zero[i] = 0;
    for (k = 0; k < od->x_vars; ++k) {
      if (!get_x_var_attr(od, i, k, SEED_BIT | GROUP_BIT)) {
        ++group_zero[i];
      }
    }
    for (j = 0; j < seed_count[i]; ++j) {
      qsort(voronoi_composition[overall_seed_count + j],
        voronoi_fill[overall_seed_count + j],
        sizeof(int), compare_integers);
    }
    overall_seed_count += seed_count[i];
    if (seed_count[i]) {
      overall_group_zero += group_zero[i];
      tee_printf(od, "%5d%12d%20d\n", i + 1, seed_count[i], group_zero[i]);
    }
    else {
      tee_printf(od, "%5d%12d%20s\n", i + 1, seed_count[i], "-");
    }
    if (od->file[TEMP_VORONOI]->handle) {
      fclose(od->file[TEMP_VORONOI]->handle);
      od->file[TEMP_VORONOI]->handle = NULL;
    }
  }
  if (overall_seed_count) {
    tee_printf(od, "%5s%12d%20d\n", "Total", overall_seed_count, overall_group_zero);
  }
  else {
    tee_printf(od, "%5s%12d%20s\n", "Total", overall_seed_count, "-");
  }
  tee_printf(od, "-------------------------------------\n\n");
  tee_flush(od);
  for (i = 0; i < od->field_num; ++i) {
    for (j = 0; j < od->x_vars; ++j) {
      set_voronoi_buf(od, i, j, -1);
    }
  }
  overall_seed_count = 0;
  for (i = 0; i < od->field_num; ++i) {
    if (get_field_attr(od, i, ACTIVE_BIT)) {
      for (j = 0; j < seed_count[i]; ++j) {
        voronoi_active[overall_seed_count + j] = 0;
        for (k = 0; k < voronoi_fill[overall_seed_count + j]; ++k) {
          set_voronoi_buf(od, i, voronoi_composition
            [overall_seed_count + j][k], overall_seed_count + j);
          if (get_x_var_attr(od, i, voronoi_composition
            [overall_seed_count + j][k], ACTIVE_BIT)) {
            ++voronoi_active[overall_seed_count + j];
          }
        }
      }
    }
    overall_seed_count += seed_count[i];
  }
  if (collapse) {
    tee_printf(od, "Collapsing Voronoi polyhedra.\n\n");
    tee_flush(od);
    overall_seed_count = 0;
    for (i = 0; i < od->field_num; ++i) {
      if (!get_field_attr(od, i, ACTIVE_BIT)) {
        continue;
      }
      seed_count_before_collapse[i] = seed_count[i];
      /*
      allocate and initialize seed coordinate vector
      */
      if (!(od->al.seed_coord = (VarCoord **)
        alloc_array(seed_count[i], sizeof(VarCoord)))) {
        return OUT_OF_MEMORY;
      }
      for (j = 0; j < seed_count[i]; ++j) {
        var_to_xyz(od, seed_list->me[i][j], od->al.seed_coord[j]);
      }
      /*
      resize P vector matrices according to the number of seeds
      in the current field
      */
      od->mal.cum_ave = double_mat_resize(od->mal.cum_ave,
        seed_count[i], od->active_object_num);
      if (!(od->mal.cum_ave)) {
        return OUT_OF_MEMORY;;
      }
      od->mal.pos_ave = double_mat_resize(od->mal.pos_ave,
        seed_count[i], od->active_object_num);
      if (!(od->mal.pos_ave)) {
        return OUT_OF_MEMORY;;
      }
      od->mal.neg_ave = double_mat_resize(od->mal.neg_ave,
        seed_count[i], od->active_object_num);
      if (!(od->mal.neg_ave)) {
        return OUT_OF_MEMORY;;
      }
      /*
      compute P vectors for all Voronoi polyhedra
      */
      for (j = 0; j < seed_count[i]; ++j) {
        result = calc_p_vectors(od, i, j);
        if (result) {
          return result;
        }
      }
      for (j = 0, max_elem = 0; j < (seed_count[i] - 1); ++j) {
        max_elem += (seed_count[i] - (j + 1));
      }
      free_char_matrix((CharMat *)(od->cimal.collapse_checked));
      od->cimal.collapse_checked = NULL;
      od->cimal.collapse_checked = alloc_int_matrix
        (od->cimal.collapse_checked, max_elem, 2);
      if (!(od->cimal.collapse_checked)) {
        return OUT_OF_MEMORY;
      }
      collapse_checked = od->cimal.collapse_checked;
      collapse_checked_len = 0;
      od->al.nearest_mat = (SeedDistMat **)alloc_array(max_elem, sizeof(SeedDistMat));
      if (!(od->al.nearest_mat)) {
        return OUT_OF_MEMORY;
      }
      collapse = 1;
      while (collapse) {
        /*
        compute pairwise distances between all seeds
        in the current field and find the two nearest ones
        */
        elem = 0;
        for (j = 0; j < (seed_count[i] - 1); ++j) {
          if (voronoi_fill[overall_seed_count + j] != -1) {
            for (k = j + 1; k < seed_count[i]; ++k) {
              if (voronoi_fill[overall_seed_count + k] != -1) {
                distance = sqrt(squared_euclidean_distance
                  (od->al.seed_coord[j]->cart, od->al.seed_coord[k]->cart))
                   / od->grid.step[0];
                if (distance <= collapse_distance) {
                  od->al.nearest_mat[elem]->dist = distance;
                  od->al.nearest_mat[elem]->seed[0] = j;
                  od->al.nearest_mat[elem]->seed[1] = k;
                  ++elem;
                }
              }
            }
          }
        }
        if (!elem) {
          break;
        }
        qsort(od->al.nearest_mat, elem, sizeof(SeedDistMat *), compare_seed_dist);
        /*
        if the nearest seeds are farther than collapse distance,
        then abort collapsing procedure
        */
        collapse = 0;
        for (attempt = 0; attempt < elem; ++attempt) {
          for (k = 0; k <= 1; ++k) {
            candidate_seed[k] = od->al.nearest_mat[attempt]->seed[k];
          }
          for (n = 0, is_checked = 0; (!is_checked) && (n < collapse_checked_len); ++n) {
            is_checked = ((candidate_seed[0] == collapse_checked->me[n][0])
              && (candidate_seed[1] == collapse_checked->me[n][1]));
          }
          if (is_checked) {
            continue;
          }
          collapse_checked->me[collapse_checked_len][0] = candidate_seed[0];
          collapse_checked->me[collapse_checked_len][1] = candidate_seed[1];
          ++collapse_checked_len;
          /*
          compare the missing value patterns between the
          Voronoi polyhedra identified by the two nearest seeds
          */
          different_pattern = 0;
          for (j = 0; (!different_pattern) && (j < od->active_object_num); ++j) {
            for (k = 0, different_pattern = 0; (!different_pattern) && (k <= 1); ++k) {
              if (((int)M_PEEK(od->mal.pos_ave, candidate_seed[k], j) == -1) &&
                ((int)M_PEEK(od->mal.pos_ave, candidate_seed[1 - k], j) != -1)) {
                different_pattern = 1;
              }
              if (((int)M_PEEK(od->mal.neg_ave, candidate_seed[k], j) == 1) &&
                ((int)M_PEEK(od->mal.neg_ave, candidate_seed[1 - k], j) != 1)) {
                different_pattern = 1;
              }
            }
          }
          if (different_pattern) {
            continue;
          }
          cum_r = pearson_r(od->mal.cum_ave, candidate_seed, 0, 0.0);
          pos_r = pearson_r(od->mal.pos_ave, candidate_seed, 1, -1.0);
          neg_r = pearson_r(od->mal.neg_ave, candidate_seed, 1, 1.0);
          if ((cum_r > 0.8) && (pos_r > 0.5)  && (neg_r > 0.5)) {
            collapse = 1;
            /*
            find coordinates of the new pseudo-seed
            */
            pseudo_seed_coord(od, i, candidate_seed);
            /*
            collapse the two corresponding Voronoi polyhedra into one
            */
            new_voronoi_composition =
              (int *)malloc((voronoi_fill[overall_seed_count + candidate_seed[0]] +
              voronoi_fill[overall_seed_count + candidate_seed[1]]) * sizeof(int));
            if (!new_voronoi_composition) {
              return OUT_OF_MEMORY;
            }
            /*
            populate the new Voronoi polyhedron
            with the populations from the two merged ones;
            also free the old arrays
            */
            n = 0;
            for (k = 0; k <= 1; ++k) {
              for (j = 0; j < voronoi_fill[overall_seed_count + candidate_seed[k]]; ++j) {
                new_voronoi_composition[n] = voronoi_composition
                  [overall_seed_count + candidate_seed[k]][j];
                ++n;
              }
              if (voronoi_composition[overall_seed_count + candidate_seed[k]]) {
                free(voronoi_composition[overall_seed_count + candidate_seed[k]]);
              }
            }
            qsort(new_voronoi_composition, n, sizeof(int), compare_integers);
            voronoi_composition[overall_seed_count + candidate_seed[1]] =
              new_voronoi_composition;
            voronoi_composition[overall_seed_count + candidate_seed[0]] = NULL;
            if (od->al.seed_coord[candidate_seed[0]]) {
              free(od->al.seed_coord[candidate_seed[0]]);
              od->al.seed_coord[candidate_seed[0]] = NULL;
            }
            /*
            update the population of the new seed
            */
            voronoi_fill[overall_seed_count + candidate_seed[1]] +=
              voronoi_fill[overall_seed_count + candidate_seed[0]];
            voronoi_active[overall_seed_count + candidate_seed[1]] +=
              voronoi_active[overall_seed_count + candidate_seed[0]];
            /*
            update the P vectors for the new seed
            */
            result = calc_p_vectors(od, i, candidate_seed[1]);
            if (result) {
              return result;
            }
            /*
            mark the lowest number Voronoi polyhedron as non-existing
            */
            voronoi_fill[overall_seed_count + candidate_seed[0]] = -1;
            voronoi_active[overall_seed_count + candidate_seed[0]] = -1;
            /*
            blacklist also the corresponding elements in the matrix
            */
            /*
            for (j = attempt; j < elem; ++j) {
              for (k = 0; (!(od->al.nearest_mat[j]->blacklisted)) && (k <= 1); ++k) {
                od->al.nearest_mat[j]->blacklisted
                  = ((od->al.nearest_mat[j]->seed[k] == candidate_seed[0])
                  || (od->al.nearest_mat[j]->seed[k] == candidate_seed[1]));
              }
            }
            */
            break;
          }
        }
      }
      /*
      update seed and Voronoi lists
      */
      j = 0;
      while (j < seed_count[i]) {
        if (voronoi_fill[overall_seed_count + j] == -1) {
          for (k = overall_seed_count + j; k < (od->voronoi_num - 1); ++k) {
            voronoi_fill[k] = voronoi_fill[k + 1];
            voronoi_active[k] = voronoi_active[k + 1];
            voronoi_composition[k] = voronoi_composition[k + 1];
          }
          voronoi_composition[k] = NULL;
          --(od->voronoi_num);
          --seed_count[i];
        }
        else {
          ++j;
        }
      }
      for (j = 0; j < seed_count[i]; ++j) {
        qsort(voronoi_composition[overall_seed_count + j],
          voronoi_fill[overall_seed_count + j],
          sizeof(int), compare_integers);
      }
      free_array(od->al.seed_coord);
      od->al.seed_coord = NULL;
      free_array(od->al.nearest_mat);
      od->al.nearest_mat = NULL;
      overall_seed_count += seed_count[i];
    }
    tee_printf(od, "%5s%20s%30s\n", "Field", "Initial groups", "Groups after collapsing");
    tee_printf(od, "-------------------------------------------------------\n");
    overall_seed_count_before_collapse = 0;
    for (i = 0; i < od->field_num; ++i) {
      if (get_field_attr(od, i, ACTIVE_BIT)) {
        overall_seed_count_before_collapse += seed_count_before_collapse[i];
        if (seed_count_before_collapse[i]) {
          tee_printf(od, "%5d%20d%30d\n", i + 1,
            seed_count_before_collapse[i], seed_count[i]);
        }
        else {
          tee_printf(od, "%5d%20s%30s\n", i + 1, " -", "-");
        }
      }
    }
    if (overall_seed_count_before_collapse) {
      tee_printf(od, "%5s%20d%30d\n", "Total",
        overall_seed_count_before_collapse, overall_seed_count);
    }
    else {
      tee_printf(od, "%5s%20s%30s\n", "Total", "-", "-");
    }
    tee_printf(od, "-------------------------------------------------------\n\n");
    if (od->mal.cum_ave) {
      double_mat_free(od->mal.cum_ave);
      od->mal.cum_ave = NULL;
    }
    if (od->mal.pos_ave) {
      double_mat_free(od->mal.pos_ave);
      od->mal.pos_ave = NULL;
    }
    if (od->mal.neg_ave) {
      double_mat_free(od->mal.neg_ave);
      od->mal.neg_ave = NULL;
    }
    if (collapse_checked) {
      free_char_matrix((CharMat *)collapse_checked);
      od->cimal.collapse_checked = NULL;
    }
  }
  for (i = 0; i < od->field_num; ++i) {
    for (j = 0; j < od->x_vars; ++j) {
      set_voronoi_buf(od, i, j, -1);
    }
  }
  overall_seed_count = 0;
  for (i = 0; i < od->field_num; ++i) {
    if (get_field_attr(od, i, ACTIVE_BIT)) {
      for (j = 0; j < seed_count[i]; ++j) {
        for (k = 0; k < voronoi_fill[overall_seed_count + j]; ++k) {
          set_voronoi_buf(od, i, voronoi_composition
            [overall_seed_count + j][k], overall_seed_count + j);
        }
      }
    }
    overall_seed_count += seed_count[i];
  }
  if (x_matrix_handle) {
    fclose(x_matrix_handle);
    od->file[TEMP_X_MATRIX]->handle = NULL;
  }
  
  return 0;  
}
