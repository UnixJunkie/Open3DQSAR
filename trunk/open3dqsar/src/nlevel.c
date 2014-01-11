/*

nlevel.c

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


int nlevel(O3Data *od)
{
  uint16_t level_bit;
  int i;
  int j;
  int k;
  int m;
  int n;
  int p;
  int len;
  int different_values;
  int equal;
  int result;
  int nlevel_op[NLEVEL_INIT];
  int count[NLEVEL_INIT];
  int two_level_count[NLEVEL_INIT];
  int three_level_count[NLEVEL_INIT][NLEVEL_INIT];
  int four_level_count[NLEVEL_INIT][NLEVEL_INIT];
  IntPerm *numberlist;
  double value;
  double nth_value[NLEVEL_INIT];
  

  /*
  initialize to zero all counters
  */
  for (i = 0; i < NLEVEL_INIT; ++i) {
    nlevel_op[i] = 0;
    two_level_count[i] = 0;
    for (j = 0; j < NLEVEL_INIT; ++j) {
      three_level_count[i][j] = 0;
      four_level_count[i][j] = 0;
    }
  }
  numberlist = od->pel.numberlist[NLEVEL_LIST];
  len = numberlist->size;
  if (!len) {
    /*
    if nlevel=all
    */
    for (n = 2; n <= 4; ++n) {
      nlevel_op[n] = 1;
    }
  }
  else {
    /*
    read which n-level variables the user
    wishes to eliminate
    */
    for (i = 0; i < len; ++i) {
      if ((numberlist->pe[i] < 2) ||
        (numberlist->pe[i] > 4)) {
        return PARSE_NLEVEL_TYPE_ERROR;
      }
      else {  
        nlevel_op[numberlist->pe[i]] = 1;
      }
    }
  }
  
  /*
  cycle over all active fields
  */
  for (i = 0; i < od->field_num; ++i) {
    if (get_field_attr(od, i, OPERATE_BIT)) {
      /*
      cycle over all active vars
      */
      for (j = 0; j < od->x_vars; ++j) {
        if (get_x_var_attr(od, i, j, ACTIVE_BIT)) {
          /*
          initialize counters
          */
          for (m = 0; m < NLEVEL_INIT; ++m) {
            count[m] = 0;
          }
          different_values = 0;
          /*
          cycle over all active objects
          */
          for (k = 0; k < od->object_num; ++k) {
            if (get_object_attr(od, k, ACTIVE_BIT)) {
              /*
              get the value of the j-th var on the k-th
              object on the i-th field
              */
              result = get_x_value(od, i, k, j, &value,
                CUTOFF_BIT);
              if (result) {
                return result;
              }
              
              equal = 0;
              /*
              check if this value is equal to one
              of the values (which are at most n) that
              have already been found previously
              */
              for (m = 0; m < different_values; ++m) {
                if (fabs((value - nth_value[m])
                  / value) < 0.0005) {
                  equal = 1;
                  /*
                  if it is almost equal to a value previously found,
                  then increase the counter for that value
                  */
                  ++count[m];
                  break;
                }
              }
              /*
              if it is different, then store it
              */
              if (!equal) {
                /*
                if 4 different values have already
                been found for this var, then break
                */
                if (different_values == 4) {
                  ++different_values;
                  break;
                }
                nth_value[different_values] = value;
                count[different_values] = 1;
                ++different_values;
              }
            }
          }
          /*
          check for 2-, 3-, 4-level vars
          */
          level_bit = 0;
          if ((different_values >= 2) && (different_values <= 4)) {
            qsort(count, different_values, sizeof(int), compare_integers);
            if ((different_values == 2) && nlevel_op[2]) {
              /*
              if one of the two different values found
              for this var appears 3 times or less, then
              remove this var
              */
              if (count[0] <= 3) {
                level_bit = TWO_LEVEL_BIT;
                ++two_level_count[count[0]];
              }
            }
            else if ((different_values == 3) && nlevel_op[3]) {
              if ((count[0] <= 2) && (count[1] <= 2)) {
                level_bit = THREE_LEVEL_BIT;
                ++three_level_count[count[0]][count[1]];
              }
            }
            else if ((different_values == 4) && nlevel_op[4]) {
              if ((count[0] == 1) && (count[1] <= 2) && (count[2] <= 2)) {
                level_bit = FOUR_LEVEL_BIT;
                ++four_level_count[count[1]][count[2]];
              }
            }
            if (level_bit) {
              set_x_var_attr(od, i, j, level_bit, 1);
            }
          }
        }
      }
    }
  }
  if (nlevel_op[2]) {
    tee_printf(od, "\n2-level variables found:\n");
    for (m = 1; m <= 3; ++m) {
      if (two_level_count[m]) {
        tee_printf(od, "%d|N        (%d)\n", m,
          two_level_count[m]);
      }
    }
    od->valid |= TWO_LEVEL_BIT;
  }
  if (nlevel_op[3]) {
    tee_printf(od, "\n3-level variables found:\n");
    for (m = 1; m <= 2; ++m) {
      for (p = 1; p <= 2; ++p) {
        if (three_level_count[m][p]) {
          tee_printf(od, "%d|%d|N      (%d)\n", m, p,
            three_level_count[m][p]);
        }
      }
    }
    od->valid |= THREE_LEVEL_BIT;
  }
  if (nlevel_op[4]) {
    tee_printf(od, "\n4-level variables found:\n");
    for (m = 1; m <= 2; ++m) {
      for (p = 1; p <= 2; ++p) {
        if (four_level_count[m][p]) {
          tee_printf(od, "1|%d|%d|N    (%d)\n", m, p,
            four_level_count[m][p]);
        }
      }
    }
    od->valid |= FOUR_LEVEL_BIT;
  }
  result = calc_active_vars(od, FULL_MODEL);

  return result;
}
