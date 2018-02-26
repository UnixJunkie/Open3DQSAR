/*

print_variables.c

is part of

Open3DQSAR
----------

An open-source software aimed at high-throughput
chemometric analysis of molecular interaction fields

Copyright (C) 2009-2018 Paolo Tosco, Thomas Balle

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


int print_variables(O3Data *od, int type)
{
  uint16_t bit = 0;
  uint16_t state = 0;
  int i;
  int k;
  VarCoord varcoord;
  int j;
  int n;
  int group_num = 0;
  int overall_seed_count = 0;
  int field_num = 0;
  IntPerm *field_list = NULL;
  IntPerm *group_list = NULL;


  switch (type) {
    case TWO_LEVEL:
    bit = TWO_LEVEL_BIT;
    state = 1;
    tee_printf(od, "2-level variables\n"
      "-----------------\n\n");
    break;

    case THREE_LEVEL:
    bit = THREE_LEVEL_BIT;
    state = 1;
    tee_printf(od, "3-level variables\n"
      "-----------------\n\n");
    break;

    case FOUR_LEVEL:
    bit = FOUR_LEVEL_BIT;
    state = 1;
    tee_printf(od, "4-level variables\n"
      "-----------------\n\n");
    break;

    case D_OPTIMAL:
    bit = D_OPTIMAL_BIT;
    state = 0;
    tee_printf(od, "D-optimal variables\n"
      "-------------------\n\n");
    break;

    case FFDSEL:
    bit = FFDSEL_BIT;
    state = 0;
    tee_printf(od, "FFD-selected variables\n"
      "----------------------\n\n");
    break;

    case SEEDS:
    bit = SEED_BIT;
    state = 0;
    tee_printf(od, "SRD-selected seeds\n"
      "----------------------\n\n");
    
    break;

    case GROUPS:
    field_list = od->pel.numberlist[FIELD_LIST];
    field_num = field_list->size;
    group_list = od->pel.numberlist[GROUP_LIST];
    group_num = group_list->size;
    if (!group_num) {
      group_num = od->mel.seed_count[field_list->pe[0] - 1];
      group_list = int_perm_resize(group_list, group_num);
      if (!group_list) {
        od->pel.numberlist[GROUP_LIST] = NULL;
        return OUT_OF_MEMORY;
      }
      od->pel.numberlist[GROUP_LIST] = group_list;
      for (i = 0; i < group_num; ++i) {
        group_list->pe[i] = i + 1;
      }
    }
    if (field_num > 1) {
      if ((group_num > 1)  ||
        ((group_num == 1) && (group_list->pe[0] != 0))) {
        return ONE_FIELD_AT_A_TIME;
      }
    }
    else {
      if ((group_list->pe[0] < 0) || (group_list->pe[group_num - 1]
        > od->mel.seed_count[field_list->pe[0] - 1])) {
        return INVALID_LIST_RANGE;
      }
    }
    bit = GROUP_BIT | SEED_BIT;
    tee_printf(od, "SRD-selected variables\n");
    tee_printf(od, "----------------------\n\n");
    break;
  }
  for (i = 0; i < od->field_num; ++i) {
    if (!get_field_attr(od, i, ACTIVE_BIT)) {
      continue;
    }
    if (get_field_attr(od, i, OPERATE_BIT)) {
      tee_printf(od, "\nField %2d\n", i + 1);
      tee_printf(od, "--------\n\n");
      if (type == GROUPS) {
        for (j = 0; j < group_num; ++j) {
          tee_printf(od, "Group %5d: %d variables (%d active)\n",
            group_list->pe[j + 1],
            od->mel.voronoi_fill
            [overall_seed_count + group_list->pe[j] - 1],
            od->mel.voronoi_active
            [overall_seed_count + group_list->pe[j] - 1]);
          tee_printf(od, "-------------------------------\n\n");
          if (group_list->pe[j] == 0) {
            for (k = 0; k < od->x_vars; ++k) {
              if (!get_x_var_attr(od, i, k, bit)) {
                var_to_xyz(od, k, &varcoord);
                tee_printf(od, "%8d%12.4lf%12.4lf%12.4lf      ", k + 1,
                  varcoord.cart[0], varcoord.cart[1], varcoord.cart[2]);
                if (get_x_var_attr(od, i, k, ACTIVE_BIT)) {
                  tee_printf(od, "%s\n", "ACTIVE");
                }
                else {
                  tee_printf(od, "%s\n", "INACTIVE");
                }
              }
            }
          }
          else {
            for (k = 0; k < od->mel.voronoi_fill
              [overall_seed_count + group_list->pe[j + 1] - 1]; ++k) {
              n = od->al.voronoi_composition
                [overall_seed_count + group_list->pe[j + 1] - 1][k];
              var_to_xyz(od, n, &varcoord);
              tee_printf(od, "%8d%12.4lf%12.4lf%12.4lf      ", n + 1,
                varcoord.cart[0], varcoord.cart[1], varcoord.cart[2]);
              if (get_x_var_attr(od, i, n, ACTIVE_BIT)) {
                tee_printf(od, "%s\n", "ACTIVE");
              }
              else {
                tee_printf(od, "%s\n", "INACTIVE");
              }
            }
          }
          tee_printf(od, "\n");
        }
      }
      else {
        for (k = 0; k < od->x_vars; ++k) {
          if (get_x_var_attr(od, i, k, ACTIVE_BIT)
            && ((get_x_var_attr(od, i, k, bit) ? 1: 0) == state)) {
            var_to_xyz(od, k, &varcoord);
            tee_printf(od, "%8d%12.4lf%12.4lf%12.4lf\n", k + 1,
              varcoord.cart[0], varcoord.cart[1], varcoord.cart[2]);
          }
        }
      }
      tee_printf(od, "\n");
    }
    overall_seed_count += od->mel.seed_count[i];
  }
  tee_printf(od, "\n");

  return 0;
}
