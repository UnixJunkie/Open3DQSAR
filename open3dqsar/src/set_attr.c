/*

set_attr.c

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


void set_field_attr(O3Data *od,
  int field_num, uint16_t attr, int onoff)
{
  if (onoff) {
    od->mel.field_attr[field_num] =
      (uint16_t)(od->mel.
      field_attr[field_num] | attr);
  }
  else {
    od->mel.field_attr[field_num] =
      (uint16_t)(od->mel.
      field_attr[field_num] & (~attr));
  }
  update_field_object_attr(od, SILENT);
}


void set_object_attr(O3Data *od, int object_num, uint16_t attr, int onoff)
{
  if (onoff) {
    od->mel.object_attr[object_num] =
      (uint16_t)(od->mel.
      object_attr[object_num] | attr);
  }
  else {
    od->mel.object_attr[object_num] =
      (uint16_t)(od->mel.
      object_attr[object_num] & (~attr));
  }
  update_field_object_attr(od, SILENT);
}


void set_x_var_attr(O3Data *od, int field_num,
  int x_var, uint16_t attr, int onoff)
{
  if (onoff) {
    od->mel.x_var_attr[field_num][x_var] =
      (uint16_t)(od->mel.
      x_var_attr[field_num][x_var] | attr);
  }
  else {
    od->mel.x_var_attr[field_num][x_var] =
      (uint16_t)(od->mel.
      x_var_attr[field_num][x_var] & (~attr));
  }
}


void set_y_var_attr(O3Data *od, int y_var, uint16_t attr, int onoff)
{
  if (onoff) {
    od->mel.y_var_attr[y_var] =
      (uint16_t)(od->mel.
      y_var_attr[y_var] | attr);
  }
  else {
    od->mel.y_var_attr[y_var] =
      (uint16_t)(od->mel.
      y_var_attr[y_var] & (~attr));  
  }
}
