/*

ff_parm.h

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


#define MMFF94_ALPHA    0
#define MMFF94_N    1
#define MMFF94_A    2
#define MMFF94_G    3
#define MMFF94_RIJ    0
#define MMFF94_RIJ7    1
#define MMFF94_EIJ    2
#define MMFF94_COUL    332.0716
#define MMFF94_ELEC_BUFF  0.05
#define MMFF94_POWER    0.25
#define MMFF94_B    0.2
#define MMFF94_BETA    12.0
#define MMFF94_DARAD    0.8
#define MMFF94_DAEPS    0.5


extern const char cosmo_label[];

typedef struct FFParm FFParm;

struct FFParm {
  short type_num;
  char type_chr[MAX_FF_TYPE_LEN];
  char da;
  double vdw_parm[MAX_FF_PARM];
};

extern FFParm ff_parm[MAX_FF_N][100];
FFParm *get_mmff_parm(unsigned int num);
