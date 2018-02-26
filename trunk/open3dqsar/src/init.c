/*

init.c

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


void init_cv_sdep(O3Data *od)
{
  od->mel.struct_list = NULL;
  od->mal.sdep_mat = NULL;
  od->mal.ave_press = NULL;
  od->vel.sumweight = NULL;
  od->vel.sum1 = NULL;
  od->vel.q2 = NULL;
  od->vel.secv = NULL;
  od->pel.sdep_rank = NULL;
}


void init_pls(O3Data *od)
{
  od->mal.e_mat = NULL;
  od->vel.e_mat_ave = NULL;
  od->mal.f_mat = NULL;
  od->vel.f_mat_ave = NULL;
  od->mal.pred_f_mat = NULL;
  od->mal.x_scores = NULL;
  od->mal.y_loadings = NULL;
  od->mal.y_scores = NULL;
  od->mal.temp = NULL;
  od->mel.ipiv = NULL;
  #ifdef HAVE_LIBMKL
  od->mel.work = NULL;
  #endif
  od->mal.x_weights = NULL;
  od->mal.x_weights_star = NULL;
  od->mal.x_loadings = NULL;
  od->mal.b_coefficients = NULL;
  od->vel.v = NULL;
  od->vel.v_new = NULL;
  od->vel.ro = NULL;
  od->vel.explained_s2_x = NULL;
  od->vel.explained_s2_y = NULL;
  od->vel.ave_sdep = NULL;
  od->vel.ave_sdec = NULL;
  od->vel.r2 = NULL;
  od->vel.r2_pred = NULL;
  od->pel.out_structs = NULL;
}
