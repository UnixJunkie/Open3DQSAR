/*

parse_comma_hyphen_list_to_array.c

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


int parse_comma_hyphen_list_to_array(O3Data *od, char *list, int type)
{
  char *ptr;
  char *context = NULL;
  char *list_copy;
  char *list_orig;
  char from[BUF_LEN];
  char *to;
  IntPerm *numberlist;
  int start;
  int end;
  int count;
  int i;
  int j;
  

  memset(from, 0, BUF_LEN);
  if (((strlen(list) >= 3) && (!strncasecmp(list, "all", 3)))
    || (!strlen(list))) {
    od->pel.numberlist[type] = int_perm_resize
      (od->pel.numberlist[type], 0);
    if (!(od->pel.numberlist[type])) {
      return OUT_OF_MEMORY;
    }
    numberlist = od->pel.numberlist[type];
    return 0;
  }
  od->mel.list_orig = realloc(od->mel.list_orig, strlen(list) + 1);
  if (!(od->mel.list_orig)) {
    return OUT_OF_MEMORY;
  }
  list_orig = od->mel.list_orig;
  memset(list_orig, 0, strlen(list) + 1);
  od->mel.list_copy = realloc(od->mel.list_copy, strlen(list) + 1);
  if (!(od->mel.list_copy)) {
    return OUT_OF_MEMORY;
  }
  list_copy = od->mel.list_copy;
  memset(list_copy, 0, strlen(list) + 1);
  strcpy(list_orig, list);
  strcpy(list_copy, list);
  /*
  tokenize list string using
  comma as separator
  */
  ptr = strtok_r(list_copy, ",\0", &context);
  count = 0;
  while (ptr) {
    strcpy(from, ptr);
    to = from;
    /*
    was a hyphen-separated
    start-end range specified in the n-th token?
    */
    while ((*to != '-') && *to) {
      ++to;
    };
    if (*to) {
      *to = '\0';
      ++to;
      sscanf(to, "%d", &end);
    }
    sscanf(from, "%d", &start);
    /*
    otherwise set a start-start range
    */
    if (!(*to)) {
      end = start;
    }
    /*
    if start > end in the specified range
    bomb out
    */
    if (start > end) {
      return INVALID_LIST_RANGE;
    }
    count += (end - start + 1);
    ptr = strtok_r(NULL, ",\0", &context);
  }
  od->pel.numberlist[type] = int_perm_resize
    (od->pel.numberlist[type], count);
  if (!(od->pel.numberlist[type])) {
    return OUT_OF_MEMORY;
  }
  numberlist = od->pel.numberlist[type];
  
  ptr = strtok_r(list_orig, ",\0", &context);
  i = 0;
  while (ptr) {
    strcpy(from, ptr);
    to = from;
    /*
    was a hyphen-separated
    start-end range specified in the n-th token?
    */
    while ((*to != '-') && *to) {
      ++to;
    };
    if (*to) {
      *to = '\0';
      ++to;
      sscanf(to, "%d", &end);
    }
    sscanf(from, "%d", &start);
    /*
    otherwise set a start-start range
    */
    if (!(*to)) {
      end = start;
    }
    for (j = start; j <= end; ++j) {  
      numberlist->pe[i] = j;
      ++i;
    }
    ptr = strtok_r(NULL, ",\0", &context);
  }
  qsort(numberlist->pe, count, sizeof(int), compare_integers);

  return 0;
}
