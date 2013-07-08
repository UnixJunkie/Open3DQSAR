/*

get_args.c

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


char *get_args(O3Data *od, char *parameter_name)
{
  int i;
  int j;
  int len;
  int parameter_found;
  #ifndef WIN32
  char *home;
  int prefix_len;
  #endif

  
  /*
  if parameter_name == NULL
  it is a sanity check;
  just check if syntax is ok
  */
  /*
  loop over all parameters
  */
  for (i = 1, parameter_found = 0; i < od->argn; ++i) {
    j = 0;
    /*
    until end-of-string is not reached and parameter length
    is below MAX_NAME_LEN keep searching for '='
    */
    while (od->cimal.arg->me[i][j] && (j < MAX_NAME_LEN)
      && (od->cimal.arg->me[i][j] != '=')) {
      ++j;
    }
    /*
    if we can't find '='
    */
    if (!(od->cimal.arg->me[i][j])) {
      return od->cimal.arg->me[i];
    }
    /*
    if this is not a sanity check,
    check if it is the parameter we were looking for;
    if it is, then break
    */
    if (parameter_name) {
      len = strlen(parameter_name);
      parameter_found = ((len == j) && (!strncasecmp(parameter_name,
        od->cimal.arg->me[i], len)));
    }
    if (parameter_found) {
      break;
    }
  }
  /*
  if no parameters matches the requested one, return NULL
  */
  if (!parameter_found) {
    return NULL;
  }
  od->cimal.arg->me[i][BUF_LEN - 1] = '\0';
  /*
  If the argument is between quotes, remove them
  */
  if (od->cimal.arg->me[i][j + 1] == '\"') {
    ++j;
    len = strlen(&(od->cimal.arg->me[i][j + 1]));
    if (od->cimal.arg->me[i][j + len] == '\"') {
      od->cimal.arg->me[i][j + len] = '\0';
    }
  }
  #ifndef WIN32
  /*
  if we are under *NIX and the argument
  begins with ~/, then expand it
  */
  if ((od->cimal.arg->me[i][j + 1] == '~')
    && ((od->cimal.arg->me[i][j + 2] == SEPARATOR)
    || (!(od->cimal.arg->me[i][j + 2])))) {
    home = getenv("HOME");
    if (home) {
      len = strlen(&(od->cimal.arg->me[i][j + 1]));
      prefix_len = strlen(home) + 1;
      while (len) {
        od->cimal.arg->me[i][j + len + prefix_len] =
          od->cimal.arg->me[i][j + len + 1];
        --len;
      }
      memcpy(&(od->cimal.arg->me[i][j + 1]), home, prefix_len);
      od->cimal.arg->me[i][j + prefix_len] = SEPARATOR;
    }
  }
  #endif
  
  return (&(od->cimal.arg->me[i][j + 1]));
}

