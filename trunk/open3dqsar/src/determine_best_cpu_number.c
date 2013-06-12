/*

determine_best_cpu_number.c

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


void determine_best_cpu_number(O3Data *od, char *parameter)
{
  int max_procs;
  
  
  od->n_proc = 0;
  max_procs = get_number_of_procs();
  if (max_procs) {
    if (max_procs > MAX_THREADS) {
      max_procs = MAX_THREADS;
    }
    if (!strncasecmp(parameter, "all", 3)) {
      od->n_proc = max_procs;
    }
  }
  if (!(od->n_proc)) {
    sscanf(parameter, "%d", &(od->n_proc));
    if (od->n_proc < 1) {
      od->n_proc = 1;
    }
    if (max_procs) {
      if (od->n_proc > max_procs) {
        tee_printf(od, "It would not be convenient to start "
          "more threads than available CPUs.\n");
        od->n_proc = max_procs;
      }
    }
  }
}
