/*

reload_coeff_weights_loadings.c

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


int reload_coefficients(O3Data *od, int pc_num)
{
  int actual_len;
  int x_max_x;
  int x_max_y;
  int x;
  int y;
  int num_blocks;
  double coeff;
  
  od->file[TEMP_PLS_COEFF]->handle = fopen(od->file[TEMP_PLS_COEFF]->name, "rb");
  if (!(od->file[TEMP_PLS_COEFF]->handle)) {
    return CANNOT_READ_TEMP_FILE;
  }
  /*
  Read the number of blocks in the coefficients file
  */
  actual_len = fread(&num_blocks, sizeof(int), 1, od->file[TEMP_PLS_COEFF]->handle);
  if (!actual_len) {
    return CANNOT_READ_TEMP_FILE;
  }
  /*
  Skip the blocks containing the PCs we are not interested in
  */
  x_max_x = od->overall_active_x_vars;
  x_max_y = od->y_vars;
  actual_len = fseek(od->file[TEMP_PLS_COEFF]->handle, x_max_y * x_max_x
    * (num_blocks - pc_num - 1) * sizeof(double), SEEK_CUR);
  if (actual_len == -1) {
    return CANNOT_READ_TEMP_FILE;
  }
  /*
  Read coefficients for the requested PC
  */
  for (x = 0; x < x_max_y; ++x) {
    for (y = 0; y < x_max_x; ++y) {
      actual_len = fread(&coeff,
        sizeof(double), 1, od->file[TEMP_PLS_COEFF]->handle);
      M_POKE(od->mal.b_coefficients, y, x, coeff);
      if (!actual_len) {
        return CANNOT_READ_TEMP_FILE;
      }
    }
  }
  fclose(od->file[TEMP_PLS_COEFF]->handle);
  od->file[TEMP_PLS_COEFF]->handle = NULL;
  
  return 0;
}


int reload_weights_loadings(O3Data *od)
{
  int x_max_x;
  int x_max_y;
  int actual_len;
  
  
  /*
  reload original weights and loadings
  */
  x_max_x = od->overall_active_x_vars;
  x_max_y = od->y_vars;
  od->file[TEMP_WLS]->handle =
    fopen(od->file[TEMP_WLS]->name, "rb");
  if (!(od->file[TEMP_WLS]->handle)) {
    return CANNOT_READ_TEMP_FILE;
  }
  
  actual_len = fread(&(od->pc_num), sizeof(int), 1,
    od->file[TEMP_WLS]->handle);
  if (actual_len != 1) {
    return CANNOT_READ_TEMP_FILE;
  }
  if (od->mal.x_loadings) {
    double_mat_free(od->mal.x_loadings);
  }
  od->mal.x_loadings = double_mat_alloc(x_max_x, od->pc_num + 1);
  if (!(od->mal.x_loadings)) {
    return OUT_OF_MEMORY;
  }
  if (od->mal.y_loadings) {
    double_mat_free(od->mal.y_loadings);
  }
  od->mal.y_loadings = double_mat_alloc(x_max_x, od->pc_num + 1);
  if (!(od->mal.y_loadings)) {
    return OUT_OF_MEMORY;
  }
  if (od->mal.x_weights) {
    double_mat_free(od->mal.x_weights);
  }
  od->mal.x_weights = double_mat_alloc(x_max_x, od->pc_num + 1);
  if (!(od->mal.x_weights)) {
    return OUT_OF_MEMORY;
  }
  if (od->mal.x_weights_star) {
    double_mat_free(od->mal.x_weights_star);
  }
  od->mal.x_weights_star = double_mat_alloc(x_max_x, od->pc_num + 1);
  if (!(od->mal.x_weights_star)) {
    return OUT_OF_MEMORY;
  }
  actual_len = fread(od->mal.x_weights->base,
    sizeof(double), x_max_x * (od->pc_num + 1),
    od->file[TEMP_WLS]->handle);
  if (actual_len != (x_max_x * (od->pc_num + 1))) {
    return CANNOT_READ_TEMP_FILE;
  }

  actual_len = fread(&(od->pc_num), sizeof(int), 1,
    od->file[TEMP_WLS]->handle);
  if (actual_len != 1) {
    return CANNOT_READ_TEMP_FILE;
  }
  actual_len = fread(od->mal.x_loadings->base,
    sizeof(double), x_max_x * (od->pc_num + 1),
    od->file[TEMP_WLS]->handle);
  if (actual_len != (x_max_x * (od->pc_num + 1))) {
    return CANNOT_READ_TEMP_FILE;
  }

  actual_len = fread(&(od->pc_num), sizeof(int), 1,
    od->file[TEMP_WLS]->handle);
  if (actual_len != 1) {
    return CANNOT_READ_TEMP_FILE;
  }
  actual_len = fread(od->mal.y_loadings->base,
    sizeof(double), x_max_y * (od->pc_num + 1),
    od->file[TEMP_WLS]->handle);
  if (actual_len != (x_max_y * (od->pc_num + 1))) {
    return CANNOT_READ_TEMP_FILE;
  }
  fclose(od->file[TEMP_WLS]->handle);
  od->file[TEMP_WLS]->handle = NULL;

  return 0;
}
