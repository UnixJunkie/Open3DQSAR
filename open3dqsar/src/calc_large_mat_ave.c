/*

calc_large_mat_ave.c

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


void calc_large_mat_ave(O3Data *od, DoubleMat *large_mat, DoubleMat *large_mat_ave, int run)
{
  int y;
  int x;
  int n;
  int object_num;
  int struct_num;
  int conf_num;
  int n_conf;
  double value;
  double sumweight;
  
  
  /*
  calculate averages from columns of large_(e,f)_mat,
  and store them into large_(e,f)_mat_ave
  */
  for (x = 0; x < large_mat->n; ++x) {
    n = 0;
    object_num = 0;
    y = 0;
    sumweight = 0.0;
    while (object_num < od->object_num) {
      struct_num = od->al.mol_info[object_num]->struct_num;
      conf_num = 1;
      if (n < od->pel.out_structs->size) {
        if (struct_num == od->pel.out_structs->pe[n]) {
          ++n;
          object_num += conf_num;
          y += conf_num;
          continue;
        }
      }
      for (n_conf = 0; n_conf < conf_num; ++n_conf, ++object_num) {
        if (get_object_attr(od, object_num, ACTIVE_BIT)) {
          value = M_PEEK(large_mat, y, x);
          if (!MISSING(value)) {
            M_POKE(large_mat_ave, run, x,
              M_PEEK(large_mat_ave, run, x)
              + value * od->mel.object_weight[object_num]);
            sumweight += od->mel.object_weight[object_num];
          }
          ++y;
        }
      }
    }
    if (sumweight > 0.0) {
      M_POKE(large_mat_ave, run, x,
        M_PEEK(large_mat_ave, run, x) / sumweight);
    }
  }
}
