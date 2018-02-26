/*

time.c

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

#include <include/o3header.h>
#ifdef WIN32
#include <windows.h>
#endif


void elapsed_time(O3Data *od, struct timeval *start, struct timeval *end)
{
  double seconds;
  struct timeval elapsed;
  
  memset(&elapsed, 0, sizeof(struct timeval));
  seconds = (double)0;
  if (!timeval_subtract(&elapsed, end, start)) {
    seconds = (double)(elapsed.tv_sec) +
      (double)(elapsed.tv_usec) / (double)1.0e06;
    tee_printf(od, "\nElapsed time: %.4lf seconds.\n\n", seconds);
  }
  else {
    tee_printf(od, "\nElapsed time appears to be negative. "
      "This CPU must be really fast... :-)\n\n");
  }
  tee_flush(od);
}


int get_current_time(char *time_string)
{
  #ifndef WIN32
  time_t current_time;
  #else
  char time[BUF_LEN];
  char date[BUF_LEN];
  SYSTEMTIME current_time;
  #endif
  
  
  #ifndef WIN32
  if (time(&current_time) == -1) {
    return TIME_NOT_AVAILABLE;
  }
  if (!ctime_r(&current_time, time_string)) {
    return TIME_NOT_AVAILABLE;
  }
  #else
  GetLocalTime(&current_time);
  GetTimeFormat(LOCALE_USER_DEFAULT, 0,
    (const SYSTEMTIME *)&current_time,
    NULL, time, BUF_LEN);
  GetDateFormat(LOCALE_USER_DEFAULT, 0,
    (const SYSTEMTIME *)&current_time,
    NULL, date, BUF_LEN);
  snprintf(time_string, BUF_LEN - 1, "%s at %s\n", date, time);
  #endif
  
  return 0;
}


int timeval_subtract(struct timeval *result, struct timeval *x, struct timeval *y)
{
  /*
  Perform the carry for the later subtraction by updating y.
  */
  if (x->tv_usec < y->tv_usec) {
    int nsec = (y->tv_usec - x->tv_usec) / 1000000 + 1;
    y->tv_usec -= 1000000 * nsec;
    y->tv_sec += nsec;
  }
  if (x->tv_usec - y->tv_usec > 1000000) {
    int nsec = (x->tv_usec - y->tv_usec) / 1000000;
    y->tv_usec += 1000000 * nsec;
    y->tv_sec -= nsec;
  }

  /*
  Compute the time remaining to wait
  tv_usec is certainly positive
  */
  result->tv_sec = x->tv_sec - y->tv_sec;
  result->tv_usec = x->tv_usec - y->tv_usec;

  /*
  Return 1 if result is negative.
  */
  return x->tv_sec < y->tv_sec;
}
