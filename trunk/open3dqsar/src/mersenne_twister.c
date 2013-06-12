/* 
   A C-program for MT19937, with initialization improved 2002/1/26.
   Coded by Takuji Nishimura and Makoto Matsumoto.

   Before using, initialize the state by using init_genrand(seed)  
   or init_by_array(init_key, key_length).

   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.                          

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. * Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote 
        products derived from this software without specific prior written 
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


   Any feedback is very welcome.
   http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
   email: m-mat @ math.sci.hiroshima-u.ac.jp (remove space)
*/

#include <include/o3header.h>


/* initializes mt[MERSENNE_N] with a seed */
void init_genrand(O3Data *od, unsigned long s)
{
  od->mt[0]= s & 0xffffffffUL;
  for (od->mti = 1; od->mti < MERSENNE_N; ++(od->mti)) {
    od->mt[od->mti] = 
      (1812433253UL * (od->mt[od->mti-1]
      ^ (od->mt[od->mti-1] >> 30)) + od->mti); 
    /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for coefficient. */
    /* In the previous versions, MSBs of the seed affect   */
    /* only MSBs of the array mt[].                        */
    /* 2002/01/09 modified by Makoto Matsumoto             */
    od->mt[od->mti] &= 0xffffffffUL;
    /* for >32 bit machines */
  }
}

/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32(O3Data *od)
{
  unsigned long y;
  static unsigned long mag01[2] = { 0x0UL, MATRIX_A };
  int kk;

  /* mag01[x] = x * MATRIX_A  for x=0,1 */

  if (od->mti >= MERSENNE_N) { /* generate MERSENNE_N words at one time */
    if (od->mti == (MERSENNE_N + 1))   /* if init_genrand() has not been called, */
      init_genrand(od, od->random_seed); /* a default initial seed is used */

    for (kk = 0; kk < (MERSENNE_N - MERSENNE_M); ++kk) {
      y = (od->mt[kk] & UPPER_MASK) | (od->mt[kk + 1] & LOWER_MASK);
      od->mt[kk] = od->mt[kk + MERSENNE_M] ^ (y >> 1) ^ mag01[y & 0x1UL];
    }
    for (; kk < (MERSENNE_N - 1); ++kk) {
      y = (od->mt[kk] & UPPER_MASK) | (od->mt[kk + 1] & LOWER_MASK);
      od->mt[kk] = od->mt[kk + (MERSENNE_M - MERSENNE_N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
    }
    y = (od->mt[MERSENNE_N - 1] & UPPER_MASK) | (od->mt[0] & LOWER_MASK);
    od->mt[MERSENNE_N - 1] = od->mt[MERSENNE_M - 1] ^ (y >> 1) ^ mag01[y & 0x1UL];

    od->mti = 0;
  }
  
  y = od->mt[(od->mti)++];

  /* Tempering */
  y ^= (y >> 11);
  y ^= (y << 7) & 0x9d2c5680UL;
  y ^= (y << 15) & 0xefc60000UL;
  y ^= (y >> 18);

  return y;
}

/* generates a random number on [0,1]-real-interval */
double genrand_real(O3Data *od)
{
  return genrand_int32(od) * (1.0 / 4294967295.0); 
  /* divided by 2^32-1 */ 
}
