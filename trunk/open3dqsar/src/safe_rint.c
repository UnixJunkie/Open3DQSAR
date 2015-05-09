/*

safe_rint.c

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


#include <include/safe_rint.h>


#ifdef HAVE_SAFE_RINT
#ifdef __x86_64__
double safe_rint(double x)
{
  int64_t i0,sx;
  int32_t j0;
  EXTRACT_WORDS64(i0,x);
  sx = (i0>>63)&1;
  j0 = ((i0>>52)&0x7ff)-0x3ff;
  if(j0<52) {
      if(j0<0) {
    if((i0 & UINT64_C(0x7fffffffffffffff))==0) return x;
    uint64_t i = i0 & UINT64_C(0xfffffffffffff);
    i0 &= UINT64_C(0xfffe000000000000);
    i0 |= (((i|-i) >> 12) & UINT64_C(0x8000000000000));
    INSERT_WORDS64(x,i0);
    double w = TWO52[sx]+x;
    double t =  w-TWO52[sx];
    EXTRACT_WORDS64(i0,t);
    INSERT_WORDS64(t,(i0&UINT64_C(0x7fffffffffffffff))|(sx<<63));
    return t;
      } else {
    uint64_t i = UINT64_C(0x000fffffffffffff)>>j0;
    if((i0&i)==0) return x; /* x is integral */
    i>>=1;
    if((i0&i)!=0)
        i0 = (i0&(~i))|(UINT64_C(0x4000000000000)>>j0);
      }
  } else {
      if(j0==0x400) return x+x;  /* inf or NaN */
      else return x;    /* x is integral */
  }
  INSERT_WORDS64(x,i0);
  double w = TWO52[sx]+x;
  return w-TWO52[sx];
}
#endif
#endif
