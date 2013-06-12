/*

predict.c

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


int predict(O3Data *od, int pc_num)
{
  int result;
  
  
  trim_mean_center_matrix(od, od->mal.large_e_mat,
    &(od->mal.e_mat), &(od->vel.e_mat_ave),
    FULL_MODEL, od->active_object_num);
  trim_mean_center_matrix(od, od->mal.large_f_mat,
    &(od->mal.f_mat), &(od->vel.f_mat_ave),
    FULL_MODEL, od->active_object_num);
  result = pred_ext_y_values(od, pc_num, FULL_MODEL);
  if (result) {
    return result;
  }
  rewind(od->file[TEMP_EXT_PRED]->handle);
  result = print_ext_pred_values(od);
  if (result) {
    return CANNOT_READ_TEMP_FILE;
  }

  return 0;
}
