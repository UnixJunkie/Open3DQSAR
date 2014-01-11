/*

tee.c

is part of

Open3DQSAR
----------

An open-source software aimed at high-throughput
chemometric analysis of molecular interaction fields

Copyright (C) 2009-2014 Paolo Tosco, Thomas Balle

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
#include <stdarg.h>


void tee_error(O3Data *od, int run_type, int overall_line_num, char *fmt, ...)
{
  int i;
  FILE *out_stream[2];
  FileDescriptor *source;
  
  out_stream[0] =  ((od->prompt || (!(od->out))) ? stdout : NULL);
  out_stream[1] = od->out;
  source = od->mel.source;
  while (source && source->next) {
    source = source->next;
  }
  va_list arg;
  for (i = 0; i < 2; ++i) {
    if  (out_stream[i]) {
      if (!(run_type & INTERACTIVE_RUN)) {
        fprintf(out_stream[i], "%s - Line %d\n",
          (source ? source->name :
          (od->file[MAIN_INPUT]->name[0] ?
          od->file[MAIN_INPUT]->name : "Standard input")),
          overall_line_num);
      }
      va_start(arg, fmt);
      vfprintf(out_stream[i], fmt, arg);
      va_end(arg);
    }
  }
}


void tee_flush(O3Data *od)
{
  if  (!(od->out) || od->prompt) {
    fflush(stdout);
  }
  if (od->out) {
    fflush(od->out);
  }
}


void tee_printf(O3Data *od, char *fmt, ...)
{
  int i;
  FILE *out_stream[2];
  
  out_stream[0] =  ((od->prompt || (!(od->out))) ? stdout : NULL);
  out_stream[1] = od->out;
  va_list arg;
  for (i = 0; i < 2; ++i) {
    if  (out_stream[i]) {
      va_start(arg, fmt);
      vfprintf(out_stream[i], fmt, arg);
      va_end(arg);
    }
  }
}
