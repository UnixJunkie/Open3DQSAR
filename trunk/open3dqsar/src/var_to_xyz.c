/*

var_to_xyz.c

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


void var_to_xyz(O3Data *od, int x_var, VarCoord *varcoord)
{
  varcoord->node[2] = x_var
    / (od->grid.nodes[0] * od->grid.nodes[1]);
  varcoord->cart[2] = safe_rint((double)(varcoord->node[2])
    * (double)(od->grid.step[2]) * 1.0e04) / 1.0e04
    + (double)(od->grid.start_coord[2]);
  varcoord->node[1] = (x_var - varcoord->node[2]
    * od->grid.nodes[0] * od->grid.nodes[1])
    / od->grid.nodes[0];
  varcoord->cart[1] = safe_rint((double)(varcoord->node[1])
    * (double)(od->grid.step[1]) * 1.0e04) / 1.0e04
    + (double)(od->grid.start_coord[1]);
  varcoord->node[0] = x_var - varcoord->node[2]
    * od->grid.nodes[0] * od->grid.nodes[1]
    - varcoord->node[1] * od->grid.nodes[0];
  varcoord->cart[0] = safe_rint((double)(varcoord->node[0])
    * (double)(od->grid.step[0]) * 1.0e04) / 1.0e04
    + (double)(od->grid.start_coord[0]);
}


int xyz_to_var(O3Data *od, VarCoord *varcoord)
{
  return ((int)(varcoord->node[2]) * od->grid.nodes[0]
    * od->grid.nodes[1] + (int)(varcoord->node[1])
    * od->grid.nodes[0] + (int)(varcoord->node[0]));
}
