/*

store_weights_loadings.c

is part of

Open3DQSAR
----------

An open-source software aimed at high-throughput
chemometric analysis of molecular interaction fields

Copyright (C) 2009-2015 Paolo Tosco, Thomas Balle

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


int store_weights_loadings(O3Data *od)
{
  int result;
  int actual_len;
  
  
  result = open_temp_file(od, od->file[TEMP_WLS], "wls");
  if (result) {
    return CANNOT_WRITE_TEMP_FILE;
  }
  actual_len = fwrite(&(od->pc_num), sizeof(int), 1,
    od->file[TEMP_WLS]->handle);
  if (actual_len != 1) {
    return CANNOT_WRITE_TEMP_FILE;
  }
  actual_len = fwrite(od->mal.x_weights->base,
    sizeof(double), od->mal.e_mat->n * (od->pc_num + 1),
    od->file[TEMP_WLS]->handle);
  if (actual_len != (od->mal.e_mat->n * (od->pc_num + 1))) {
    return CANNOT_WRITE_TEMP_FILE;
  }

  actual_len = fwrite(&(od->pc_num), sizeof(int), 1,
    od->file[TEMP_WLS]->handle);
  if (actual_len != 1) {
    return CANNOT_WRITE_TEMP_FILE;
  }
  actual_len = fwrite(od->mal.x_loadings->base,
    sizeof(double), od->mal.e_mat->n * (od->pc_num + 1),
    od->file[TEMP_WLS]->handle);
  if (actual_len != (od->mal.e_mat->n * (od->pc_num + 1))) {
    return CANNOT_WRITE_TEMP_FILE;
  }

  actual_len = fwrite(&(od->pc_num), sizeof(int), 1,
    od->file[TEMP_WLS]->handle);
  if (actual_len != 1) {
    return CANNOT_WRITE_TEMP_FILE;
  }
  actual_len = fwrite(od->mal.y_loadings->base,
    sizeof(double), (od->mal.f_mat->n * (od->pc_num + 1)),
    od->file[TEMP_WLS]->handle);
  if (actual_len != (od->mal.f_mat->n * (od->pc_num + 1))) {
    return CANNOT_WRITE_TEMP_FILE;
  }
  fclose(od->file[TEMP_WLS]->handle);
  od->file[TEMP_WLS]->handle = NULL;
  
  return 0;
}
