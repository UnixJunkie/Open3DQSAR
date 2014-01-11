/*

prepare_cv.c

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


int prepare_cv(O3Data *od, int pc_num, int cv_type, int groups, int runs)
{
  int i;
  int j;
  int n;
  int x;
  int object_num;
  int object_num2;
  int struct_num;
  int struct_num2;
  int active_struct_num;
  int conf_num;
  int conf_num2;
  int n_conf;
  int group_num;
  int excess;
  int structs_to_be_assigned;
  int struct_count;
  int position;
  int random_struct;
  double sumweight = 0.0;
  double active_value_ave = 0.0;
  
  
  od->cv.tss = 0.0;
  od->cv.num_predictions = 0;
  od->cv.overall_cv_runs = 0;
  switch (cv_type) {
    case LEAVE_ONE_OUT:
    for (x = 0; x < od->y_vars; ++x) {
      get_attr_struct_ave(od, x, ACTIVE_BIT, &active_struct_num, &active_value_ave);
      if (active_struct_num < 2) {
        return NOT_ENOUGH_OBJECTS;
      }
      object_num = 0;
      i = 0;
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
          od->cv.tss += square(get_y_value
            (od, object_num - conf_num, x, WEIGHT_BIT) - active_value_ave);
          if (!x) {
            ++(od->cv.overall_cv_runs);
            ++(od->cv.num_predictions);
          }
        }
      }
    }
    break;

    case LEAVE_TWO_OUT:
    for (x = 0; x < od->y_vars; ++x) {
      get_attr_struct_ave(od, x, ACTIVE_BIT, &active_struct_num, &active_value_ave);
      if (active_struct_num < 3) {
        return NOT_ENOUGH_OBJECTS;
      }
      struct_num = 0;
      object_num = 0;
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
              od->cv.tss +=
                square(get_y_value(od, object_num - 1, x, WEIGHT_BIT) - active_value_ave)
                + square(get_y_value(od, object_num2 - 1, x, WEIGHT_BIT) - active_value_ave);
              if (!x) {
                ++(od->cv.overall_cv_runs);
                od->cv.num_predictions += 2;
              }
            }
          }
        }
      }
    }
    break;
    
    case LEAVE_MANY_OUT:
    /*
    initialize structure list with all active structure numbers
    */
    od->cimal.group_composition_list =
      (IntMat **)malloc(runs * sizeof(IntMat *));
    if (!(od->cimal.group_composition_list)) {
      return OUT_OF_MEMORY;
    }
    od->mel.struct_per_group =
      alloc_int_array(od->mel.struct_per_group, groups);
    if (!(od->mel.struct_per_group)) {
      return OUT_OF_MEMORY;
    }
    for (x = 0; x < od->y_vars; ++x) {
      get_attr_struct_ave(od, 0, ACTIVE_BIT,
        &active_struct_num, &(od->vel.active_value_ave->ve[x]));
    }
    excess = active_struct_num % groups;
    for (i = 0; i < groups; ++i) {
      od->mel.struct_per_group[i] = active_struct_num / groups;
      if (excess) {
        ++(od->mel.struct_per_group[i]);
        --excess;
      }
    }
    od->mal.sdep_mat = double_mat_alloc(runs, pc_num + 1);
    if (!(od->mal.sdep_mat)) {
      return OUT_OF_MEMORY;
    }
    for (i = 0; i < runs; ++i) {
      od->cimal.group_composition_list[i] = 
        alloc_int_matrix(NULL, groups,
          active_struct_num / groups + 1);
      if (!(od->cimal.group_composition_list[i])) {
        return OUT_OF_MEMORY;
      }
      group_num = 0;
      struct_count = 0;
      structs_to_be_assigned = active_struct_num;
      object_num = 0;
      n = 0;
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
          od->mel.struct_list[n] = struct_num;
          ++n;
        }
      }
      while (structs_to_be_assigned) {
        position = (int)safe_rint((double)(structs_to_be_assigned - 1)
          * genrand_real(od));
        random_struct = od->mel.struct_list[position];
        object_num = 0;
        while (1) {
          struct_num = od->al.mol_info[object_num]->struct_num;
          if (struct_num == random_struct) {
            break;
          }
          conf_num = 1;
          for (n_conf = 0;
            n_conf < conf_num; ++n_conf, ++object_num) {
            if (!get_object_attr(od, object_num, ACTIVE_BIT | PREDICT_BIT)) {
              continue;
            }
          }
        }
        for (x = 0; x < od->y_vars; ++x) {
          od->cv.tss += square(get_y_value
            (od, object_num, x, WEIGHT_BIT) - od->vel.active_value_ave->ve[x]);
        }
        ++(od->cv.num_predictions);
        /*
        take random objects out of the list shifting
        those which follow up one place
        */
        --structs_to_be_assigned;
        for (j = position; j < structs_to_be_assigned; ++j) {
          od->mel.struct_list[j] =
            od->mel.struct_list[j + 1];
        }
        /*
        assign random object to current group; if object run_count
        exceeds size of current group, then reset the counter
        and pass on to the next group
        */
        od->cimal.group_composition_list[i]->me
          [group_num][struct_count] = random_struct;
        ++struct_count;
        if (struct_count >=
          od->mel.struct_per_group[group_num]) {
          struct_count = 0;
          ++group_num;
        }
      }
    }
    od->cv.overall_cv_runs = groups * runs;
    break;
  }
  od->cv.num_predictions *= od->y_vars;
  
  return 0;
}
