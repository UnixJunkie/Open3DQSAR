/*

set_sel_included_bit.c

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


int set_sel_included_bit(O3Data *od, int use_srd_groups)
{
  int i;
  int j;
  int k;
  int voronoi_num;
  int sel_included_vars;
  
  
  /*
  initialize to zero the SEL_INCLUDED_BIT
  */
  for (i = 0; i < od->field_num; ++i) {
    for (j = 0; j < od->x_vars; ++j) {
      set_x_var_attr(od, i, j, SEL_INCLUDED_BIT, 0);
    }
  }
  sel_included_vars = 0;
  for (i = 0; i < od->field_num; ++i) {
    if (get_field_attr(od, i, ACTIVE_BIT)) {
      for (k = 0; k < od->x_vars; ++k) {
        if (get_x_var_attr(od, i, k, ACTIVE_BIT)) {
          if (use_srd_groups) {
            /*
            if SRD groups are being used, find out to which
            Voronoi polyhedron this variable belongs
            */
            voronoi_num = get_voronoi_buf(od, i, k);
            /*
            if it is not in group zero, then include it
            in the design matrix and set the SEL_INCLUDED_BIT
            accordingly
            */
            if (voronoi_num >= 0) {
              set_x_var_attr(od, i, k, SEL_INCLUDED_BIT, 1);
              ++sel_included_vars;
            }
          }
          else {
            /*
            if SRD groups are not being used, include all active
            x variables in the design matrix and set the
            SEL_INCLUDED_BIT accordingly;
            */
            set_x_var_attr(od, i, k, SEL_INCLUDED_BIT, 1);
            ++sel_included_vars;
          }
        }
      }
    }
  }
  
  return sel_included_vars;
}
