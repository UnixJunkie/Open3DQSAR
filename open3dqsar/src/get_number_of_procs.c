/*

get_number_of_procs.c

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

#ifdef WIN32

#define O3_KNOWN_ARCHITECTURE 1
#include <windows.h>

int get_number_of_procs()
{
  SYSTEM_INFO system_info;

  GetSystemInfo(&system_info);
  return (int)(system_info.dwNumberOfProcessors);
}

#endif

#ifndef __ICC
#if defined(__APPLE__) || defined(__FreeBSD__)
#define O3_KNOWN_ARCHITECTURE 1
#include <sys/param.h>
#include <sys/sysctl.h>

int get_number_of_procs()
{
  int mib[2];
  size_t len;
  int n_proc = 1;

  mib[0] = CTL_HW;
  mib[1] = HW_NCPU;
  len = sizeof(n_proc);
  (void)sysctl(mib, 2, &n_proc, &len, NULL, 0);
  return n_proc;
}

#elif sun

#define O3_KNOWN_ARCHITECTURE 1
#include <unistd.h>

int get_number_of_procs()
{
  return (int)sysconf(_SC_NPROCESSORS_CONF);
}

#endif
#endif

#ifdef linux
#define O3_KNOWN_ARCHITECTURE 1

#include <sys/sysinfo.h>

int get_number_of_procs()
{
  return get_nprocs();
}

/*
if we are running on different systems
the choice of the number of threaads to run will
be completely left to the user independently of
the number of physical CPUs/cores
*/
#endif

#ifndef O3_KNOWN_ARCHITECTURE

int get_number_of_procs()
{
  return 0;
}

#endif
