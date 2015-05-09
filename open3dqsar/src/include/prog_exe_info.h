/*

prog_exe_info.h

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


#ifdef WIN32
#include <windows.h>
#endif


struct ProgExeInfo {
  char command_line[BUF_LEN];
  char *exedir;
  int need_stdin;
  int sep_proc_grp;
  FileDescriptor *stdout_fd;
  FileDescriptor *stderr_fd;
  #ifndef WIN32
  char **proc_env;
  int pipe_des[2];
  #else
  char *proc_env;
  HANDLE stdin_rd;
  HANDLE stdin_wr;
  HANDLE des;
  HANDLE out;
  HANDLE log;
  SECURITY_ATTRIBUTES out_sa_attr;
  SECURITY_ATTRIBUTES log_sa_attr;
  SECURITY_ATTRIBUTES pipe_sa_attr;
  PROCESS_INFORMATION proc_info;
  STARTUPINFO startup_info;
  #endif
};
