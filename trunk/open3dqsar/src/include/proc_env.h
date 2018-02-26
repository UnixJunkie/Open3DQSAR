/*

proc_env.h

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


#define EL_QM_EXE_PATH      (1<<0)
#define EL_PLAIN_SET      (1<<1)
#define EL_QM_PUNCH      (1<<2)
#define EL_QM_GAUSCR      (1<<3)
#define EL_QM_GMSSCR      (1<<4)
#ifdef WIN32
#define EL_EXE_PATH      (1<<5)
#define free_proc_env(x)    if (x) free(x)
#else
#define EL_DYN_LIBRARY_PATH    (1<<5)
#define EL_DYN_32_LIBRARY_PATH    (1<<6)
#define free_proc_env(x)    if (x) free_array(x)
#endif
#define EL_TEMP_DIR      (1<<7)
#define EL_BABEL_DATADIR    (1<<8)
#define EL_BABEL_LIBDIR      (1<<9)


#ifndef WIN32
extern char **environ;
#endif


extern EnvList babel_env[];
extern EnvList minimal_env[];
extern EnvList gaussian_env[];
extern EnvList gamess_env[];
extern EnvList turbomole_env[];
