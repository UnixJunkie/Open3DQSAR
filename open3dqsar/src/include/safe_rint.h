/*

safe_rint.h

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


#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <math.h>


#ifdef HAVE_SAFE_RINT
#ifdef __x86_64__
#include <endian.h>
#include <stdint.h>
#include <sys/types.h>

#if __FLOAT_WORD_ORDER == BIG_ENDIAN

typedef union
{
  double value;
  struct
  {
    u_int32_t msw;
    u_int32_t lsw;
  } parts;
  uint64_t word;
} ieee_double_shape_type;

#endif

#if __FLOAT_WORD_ORDER == LITTLE_ENDIAN

typedef union
{
  double value;
  struct
  {
    u_int32_t lsw;
    u_int32_t msw;
  } parts;
  uint64_t word;
} ieee_double_shape_type;

#endif

#define EXTRACT_WORDS64(i,d)          \
do {                \
  ieee_double_shape_type gh_u;          \
  gh_u.value = (d);            \
  (i) = gh_u.word;            \
} while (0)

#define INSERT_WORDS64(d,i)          \
do {                \
  ieee_double_shape_type iw_u;          \
  iw_u.word = (i);            \
  (d) = iw_u.value;            \
} while (0)


#ifdef __STDC__
static const double
#else
static double
#endif
TWO52[2]={
  4.50359962737049600000e+15, /* 0x43300000, 0x00000000 */
 -4.50359962737049600000e+15, /* 0xC3300000, 0x00000000 */
};

double safe_rint(double x);
#endif
#endif
#ifndef HAVE_SAFE_RINT
#define safe_rint(x)      rint(x)
#endif
#ifndef __x86_64__
#define safe_rint(x)      rint(x)
#endif
