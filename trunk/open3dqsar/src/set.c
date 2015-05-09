/*

set.c

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


int set(O3Data *od, int type,
  uint16_t attr, int state, int verbose)
{
  int i;
  int j;
  int found;
  int result;
  int object_num;
  int struct_num;
  int conf_num;
  int len;
  int overall_conf_num;
  IntPerm *numberlist;
  

  numberlist = od->pel.numberlist[type];
  len = numberlist->size;
  
  if (type == ID_LIST) {
    od->pel.numberlist[OBJECT_LIST] = int_perm_resize
      (od->pel.numberlist[OBJECT_LIST], len);
    if (!len) {
      len = od->grid.object_num;
      fill_numberlist(od, len, OBJECT_LIST);
    }
    else {
      for (i = 0, j = 0; i < len; ++i) {
        for (found = 0; (!found) && (j < od->grid.object_num); ++j) {
          if (numberlist->pe[i] == od->al.mol_info[j]->object_id) {
            od->pel.numberlist[OBJECT_LIST]->pe[i] = j + 1;
            found = 1;
          }
        }
        if (!found) {
          return INVALID_LIST_RANGE;
        }
      }
    }
    type = OBJECT_LIST;
    numberlist = od->pel.numberlist[OBJECT_LIST];
  }
  else if (type == STRUCT_LIST) {
    object_num = 0;
    overall_conf_num = 0;
    i = 0;
    while ((object_num < od->grid.object_num) && (i < len)) {
      struct_num = od->al.mol_info[object_num]->struct_num;
      conf_num = 1;
      if (numberlist->pe[i] == (struct_num + 1)) {
        overall_conf_num += conf_num;
        ++i;
      }
      object_num += conf_num;
    }
    od->pel.numberlist[OBJECT_LIST] = int_perm_resize
      (od->pel.numberlist[OBJECT_LIST], overall_conf_num);
    if (!len) {
      len = od->grid.object_num;
      fill_numberlist(od, len, OBJECT_LIST);
    }
    else {
      i = 0;
      j = 0;
      object_num = 0;
      while (i < len) {
        found = 0;
        while ((!found) && (object_num < od->grid.object_num)) {
          struct_num = od->al.mol_info[object_num]->struct_num;
          conf_num = 1;
          if (numberlist->pe[i] == (struct_num + 1)) {
            while (conf_num) {
              od->pel.numberlist[OBJECT_LIST]->pe[j] = object_num + 1;
              ++j;
              ++object_num;
              --conf_num;
            }
            found = 1;
            ++i;
          }
          else {
            object_num += conf_num;
          }
        }
        if (!found) {
          return INVALID_LIST_RANGE;
        }
      }
    }
    len = overall_conf_num;
    type = OBJECT_LIST;
    numberlist = od->pel.numberlist[OBJECT_LIST];
  }
  if (type == OBJECT_LIST) {
    if (!len) {
      len = od->grid.object_num;
      fill_numberlist(od, len, type);
      numberlist = od->pel.numberlist[type];
    }
    for (i = 0; i < len; ++i) {
      if ((numberlist->pe[i] < 1)
        || (numberlist->pe[i] > od->grid.object_num)) {
        return INVALID_LIST_RANGE;
      }
    }
    if (attr == OPERATE_BIT) {
      /*
      initialize to zero the OPERATE_BIT on all objects
      */
      for (i = 0; i < od->grid.object_num; ++i) {
        set_object_attr(od, i, OPERATE_BIT, 0);
      }
    }
    for (i = 0; i < len; ++i) {
      set_object_attr(od, numberlist->pe[i] - 1, attr, state);
      if (attr == PREDICT_BIT) {
        /*
        the PREDICT_BIT is mutually exclusive with the
        ACTIVE_BIT and DELETE_BIT, so when PREDICT_BIT
        is switched on, ACTIVE_BIT and DELETE_BIT are
        switched off and viceversa
        */
        set_object_attr(od, numberlist->pe[i] - 1,
          ACTIVE_BIT | DELETE_BIT, 0);
      }
      else if (attr == ACTIVE_BIT) {
        /*
        either the ACTIVE_BIT is switched on or off,
        the PREDICT_BIT and DELETE_BIT must be switched off
        */
        set_object_attr(od, numberlist->pe[i] - 1,
          PREDICT_BIT | DELETE_BIT, 0);
      }
      else if (attr == DELETE_BIT) {
        /*
        if the DELETE_BIT is switched on,
        the PREDICT_BIT and ACTIVE_BIT must be switched off
        */
        set_object_attr(od, numberlist->pe[i] - 1,
          PREDICT_BIT | ACTIVE_BIT, 0);
      }
    }
  }
  else if (type == FIELD_LIST) {
    if (!len) {
      len = od->field_num;
      fill_numberlist(od, len, type);
      numberlist = od->pel.numberlist[type];
    }
    for (i = 0; i < len; ++i) {
      if ((numberlist->pe[i] < 1)
        || (numberlist->pe[i] > od->field_num)) {
        return INVALID_LIST_RANGE;
      }
    }
    if (attr == OPERATE_BIT) {
      /*
      initialize to zero the OPERATE_BIT on all fields
      */
      for (i = 0; i < od->field_num; ++i) {
        set_field_attr(od, i, OPERATE_BIT, 0);
      }
    }
    for (i = 0; i < len; ++i) {
      set_field_attr(od, numberlist->pe[i] - 1, attr, state);
      if (attr == ACTIVE_BIT) {
        /*
        either the ACTIVE_BIT is switched on or off,
        the DELETE_BIT must be switched off
        */
        set_field_attr(od,
          numberlist->pe[i] - 1, DELETE_BIT, 0);
      }
      else if (attr == DELETE_BIT) {
        /*
        if the DELETE_BIT is switched on,
        the ACTIVE_BIT must be switched off
        */
        set_field_attr(od,
          numberlist->pe[i] - 1, ACTIVE_BIT, 0);
      }
    }
  }
  else if (type == Y_VAR_LIST) {
    if (!len) {
      len = od->y_vars;
      fill_numberlist(od, len, type);
      numberlist = od->pel.numberlist[type];
    }
    for (i = 0; i < len; ++i) {
      if ((numberlist->pe[i] < 1)
        || (numberlist->pe[i] > od->y_vars)) {
        return INVALID_LIST_RANGE;
      }
    }
    if (attr == OPERATE_BIT) {
      /*
      intialize to zero the OPERATE_BIT on all y vars
      */
      for (i = 0; i < od->y_vars; ++i) {
        set_y_var_attr(od, i, OPERATE_BIT, 0);
      }
    }
    for (i = 0; i < len; ++i) {
      set_y_var_attr(od,
        numberlist->pe[i] - 1, attr, state);
    }
  }
  update_field_object_attr(od, verbose);
  result = calc_active_vars(od, verbose ? FULL_MODEL : CV_MODEL);
  
  return result;
}
