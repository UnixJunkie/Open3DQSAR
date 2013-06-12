/*

set_nice_value.c

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


#ifndef WIN32
void set_nice_value(O3Data *od, int nice_value)
{
  int current_nice;
  

  errno = 0;
  current_nice = getpriority(PRIO_PROCESS, 0);
  if (errno) {
    current_nice = -99;
  }
  errno = 0;
  (void)setpriority(PRIO_PROCESS, 0, nice_value);
  if (errno) {
    tee_printf(od, "The nice value could not be changed "
      "from the current value");
    if (current_nice == -99) {
      tee_printf(od, ".\n\n");
    }
    else {
      tee_printf(od, " of %d.\n\n", current_nice);
    }
  }
  else {
    tee_printf(od, "The nice value has been set to %d.\n\n", nice_value);
  }
}
#else
#include <include/nice_windows.h>


void set_nice_value(O3Data *od, int nice_value)
{
  int current_nice_index;
  BOOL result;
  DWORD current_nice;
  HANDLE process_handle;
  
  
  process_handle = GetCurrentProcess();
  current_nice = GetPriorityClass(process_handle);
  current_nice_index = 0;
  while (current_nice_index < 6) {
    if (current_nice == nice_dword[current_nice_index]) {
      break;
    }
    ++current_nice_index;
  }
  result = 0;
  if (nice_value < 6) {
    result = SetPriorityClass(process_handle, nice_dword[nice_value]);
  }
  if (!result) {
    tee_printf(od, "The nice value could not be changed "
      "from the current value (%s).\n", nice_name[current_nice_index]);
  }
  else {
    tee_printf(od, "The nice value has been "
      "set to %s.\n\n", nice_name[nice_value]);
  }
}
#endif
