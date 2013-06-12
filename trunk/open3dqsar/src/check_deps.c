/*

check_deps.c

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
#include <include/prog_exe_info.h>
#include <include/proc_env.h>


extern EnvList babel_env[];

int check_babel(O3Data *od, char *bin)
{
  char buffer[BUF_LEN];
  char obenergy_dir[BUF_LEN];
  char cu1pw1[] =
    "CU1PW1\n"
    " OpenBabel          3D\n"
    "\n"
    "  4  2  0  0  0  0  0  0  0  0999 V2000\n"
    "   -1.9861   -0.7088    4.1281 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "   -2.5426   -1.0556    3.3951 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "   -2.5672   -0.9578    4.8815 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "   -0.0086    0.3270    4.0927 Cu  0  3  0  0  0  0  0  0  0  0  0  0\n"
    "  1  3  1  0  0  0  0\n"
    "  1  2  1  0  0  0  0\n"
    "M  CHG  1   4   1\n"
    "M  END\n"
    "$$$$\n";
  int i;
  int n;
  int error;
  int result;
  int pid = 0;
  ProgExeInfo prog_exe_info;


  memset(buffer, 0, BUF_LEN);
  memset(obenergy_dir, 0, BUF_LEN);
  for (i = 0, result = 0; i < 3; ++i) {
    memset(&prog_exe_info, 0, sizeof(ProgExeInfo));
    if (open_temp_file(od, od->file[TEMP_LOG], "babel_log")) {
      return CANNOT_WRITE_TEMP_FILE;
    }
    fclose(od->file[TEMP_LOG]->handle);
    od->file[TEMP_LOG]->handle = NULL;
    prog_exe_info.proc_env = fill_env(od, babel_env, bin, 0);
    if (!(prog_exe_info.proc_env)) {
      return OUT_OF_MEMORY;
    }
    prog_exe_info.need_stdin = 0;
    prog_exe_info.stdout_fd = od->file[TEMP_LOG];
    prog_exe_info.stderr_fd = od->file[TEMP_LOG];
    prog_exe_info.exedir = od->field.babel_exe_path;
    prog_exe_info.sep_proc_grp = 0;
    if (i < 2) {
      sprintf(prog_exe_info.command_line,
        "%s%c"BABEL_EXE" -L sdf",
        od->field.babel_exe_path, SEPARATOR);
    }
    else {
      memset(od->file[TEMP_BABEL]->name, 0, BUF_LEN);
      if (open_temp_dir(od, NULL, "obenergy_test", obenergy_dir)) {
        return CANNOT_WRITE_TEMP_FILE;
      }
      sprintf(od->file[TEMP_BABEL]->name, "%s%ccu1pw1.sdf", obenergy_dir, SEPARATOR);
      if (!(od->file[TEMP_BABEL]->handle = fopen(od->file[TEMP_BABEL]->name, "wb"))) {
        return CANNOT_WRITE_TEMP_FILE;
      }
      fputs(cu1pw1, od->file[TEMP_BABEL]->handle);
      fputs(cu1pw1, od->file[TEMP_BABEL]->handle);
      fclose(od->file[TEMP_BABEL]->handle);
      od->file[TEMP_BABEL]->handle = NULL;
      sprintf(prog_exe_info.command_line,
        "%s%c"OBENERGY_EXE" -ff MMFF94 %s",
        od->field.babel_exe_path, SEPARATOR, od->file[TEMP_BABEL]->name);
    }
    pid = ext_program_exe(&prog_exe_info, &error);
    ext_program_wait(&prog_exe_info, pid);
    free_proc_env(prog_exe_info.proc_env);
    if (!(od->file[TEMP_LOG]->handle = fopen
      (od->file[TEMP_LOG]->name, "rb"))) {
      return CANNOT_READ_TEMP_FILE;
    }
    if (i < 2) {
      result = 0;
      if (!fgrep(od->file[TEMP_LOG]->handle, buffer, "MDL MOL format")) {
        result = BABEL_PLUGINS_NOT_FOUND;
        sprintf(od->field.babel_libdir, "%s%cplugins", od->field.babel_exe_path, SEPARATOR);
        sprintf(od->field.babel_datadir, "%s%c..%cshare%copenbabel",
          od->field.babel_exe_path, SEPARATOR, SEPARATOR, SEPARATOR);
      }
    }
    else if (!result) {
      n = 0;
      while (fgets(buffer, BUF_LEN, od->file[TEMP_LOG]->handle)) {
        buffer[BUF_LEN - 1] = '\0';
        if (!strncmp(buffer, "4\t97", 4)) {
          ++n;
        }
      }
      result = ((n < 2) ? BABEL_TOO_OLD : 0);
      remove_recursive(obenergy_dir);
    }
    else {
      result = BABEL_NOT_WORKING;
    }
    fclose(od->file[TEMP_LOG]->handle);
    od->file[TEMP_LOG]->handle = NULL;
    if ((result == BABEL_NOT_WORKING) || (result == BABEL_TOO_OLD)) {
      break;
    }
    remove(od->file[TEMP_LOG]->name);
  }

  return result;
}


int check_pharao(O3Data *od, char *bin)
{
  char buffer[BUF_LEN];
  char pharao_dir[BUF_LEN];
  char methane[] =
    "METHANE\n"
    " OpenBabel          3D\n"
    "\n"
    "  5  4  0  0  0  0  0  0  0  0999 V2000\n"
    "    0.0000   -0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    1.0922   -0.0000   -0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "   -0.3640    0.7281    0.7281 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "   -0.3641   -0.9946    0.2665 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "   -0.3641    0.2665   -0.9946 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "  1  2  1  0  0  0  0\n"
    "  1  3  1  0  0  0  0\n"
    "  1  4  1  0  0  0  0\n"
    "  1  5  1  0  0  0  0\n"
    "M  END\n"
    "$$$$\n";
  int i;
  int error;
  int result;
  int pid = 0;
  ProgExeInfo prog_exe_info;


  memset(buffer, 0, BUF_LEN);
  memset(pharao_dir, 0, BUF_LEN);
  for (i = 0, result = 0; i < 2; ++i) {
    memset(&prog_exe_info, 0, sizeof(ProgExeInfo));
    if (open_temp_file(od, od->file[TEMP_LOG], "pharao_log")) {
      return CANNOT_WRITE_TEMP_FILE;
    }
    fclose(od->file[TEMP_LOG]->handle);
    od->file[TEMP_LOG]->handle = NULL;
    prog_exe_info.proc_env = fill_env(od, babel_env, bin, 0);
    if (!(prog_exe_info.proc_env)) {
      return OUT_OF_MEMORY;
    }
    prog_exe_info.need_stdin = 0;
    prog_exe_info.stdout_fd = od->file[TEMP_LOG];
    prog_exe_info.stderr_fd = od->file[TEMP_LOG];
    prog_exe_info.exedir = od->field.babel_exe_path;
    prog_exe_info.sep_proc_grp = 0;
    memset(od->file[TEMP_BABEL]->name, 0, BUF_LEN);
    if (open_temp_dir(od, NULL, "pharao_test", pharao_dir)) {
      return CANNOT_WRITE_TEMP_FILE;
    }
    sprintf(od->file[TEMP_BABEL]->name, "%s%cmethane.sdf", pharao_dir, SEPARATOR);
    if (!(od->file[TEMP_BABEL]->handle = fopen(od->file[TEMP_BABEL]->name, "wb"))) {
      return CANNOT_WRITE_TEMP_FILE;
    }
    fputs(methane, od->file[TEMP_BABEL]->handle);
    fclose(od->file[TEMP_BABEL]->handle);
    od->file[TEMP_BABEL]->handle = NULL;
    sprintf(prog_exe_info.command_line,
      "%s --noHybrid -d %s -p %s%cmethane.phar", od->align.pharao_exe,
      od->file[TEMP_BABEL]->name, pharao_dir, SEPARATOR);
    pid = ext_program_exe(&prog_exe_info, &error);
    ext_program_wait(&prog_exe_info, pid);
    free_proc_env(prog_exe_info.proc_env);
    if (!(od->file[TEMP_LOG]->handle = fopen
      (od->file[TEMP_LOG]->name, "rb"))) {
      return CANNOT_READ_TEMP_FILE;
    }
    sprintf(od->file[TEMP_BABEL]->name, "%s%cmethane.phar", pharao_dir, SEPARATOR);
    if (!(od->file[TEMP_BABEL]->handle = fopen
      (od->file[TEMP_BABEL]->name, "rb"))) {
      return CANNOT_READ_TEMP_FILE;
    }
    if (!i) {
      if (fgrep(od->file[TEMP_LOG]->handle, buffer, "Warning")
        || fgrep(od->file[TEMP_LOG]->handle, buffer, "Processed 0 molecules")) {
        result = BABEL_PLUGINS_NOT_FOUND;
        sprintf(od->field.babel_libdir, "%s%cplugins", od->field.babel_exe_path, SEPARATOR);
        sprintf(od->field.babel_datadir, "%s%c..%cshare%copenbabel",
          od->field.babel_exe_path, SEPARATOR, SEPARATOR, SEPARATOR);
      }
    }
    if ((!result) && ((!fgrep(od->file[TEMP_LOG]->handle, buffer, "Processed 1 molecules"))
      || (!fgrep(od->file[TEMP_BABEL]->handle, buffer, "LIPO")))) {
      result = BABEL_NOT_WORKING;
    }
    fclose(od->file[TEMP_LOG]->handle);
    od->file[TEMP_LOG]->handle = NULL;
    fclose(od->file[TEMP_BABEL]->handle);
    od->file[TEMP_BABEL]->handle = NULL;
    remove_recursive(pharao_dir);
    if (result != BABEL_PLUGINS_NOT_FOUND) {
      break;
    }
    remove(od->file[TEMP_LOG]->name);
  }

  return result;
}


int check_define(O3Data *od, char *bin)
{
  char buffer[BUF_LEN];
  char define_dir[BUF_LEN];
  char basis_set_name[MAX_NAME_LEN];
  char lowercase_elem[MAX_NAME_LEN];
  char methane[] =
    "$coord\n"
    "    0.00000000000000     -0.00000000000000      0.00000000000000      c\n"
    "    2.06395887371800     -0.00000000000000     -0.00000000000000      h\n"
    "   -0.68786030949767      1.37590959160783      1.37590959160783      h\n"
    "   -0.68804928211017     -1.87952160391862      0.50361201231079      h\n"
    "   -0.68804928211017      0.50361201231079     -1.87952160391862      h\n"
    "$user-defined bonds\n"
    "$end\n";
  int error = 0;
  int pid = 0;
  ProgExeInfo prog_exe_info;
  FileDescriptor log_fd;
  FileDescriptor coord_fd;
  #ifndef WIN32
  FILE *pipe_handle = NULL;
  #else
  DWORD n_chr = 0;
  HANDLE pipe_handle;
  #endif


  memset(define_dir, 0, BUF_LEN);
  memset(&prog_exe_info, 0, sizeof(ProgExeInfo));
  if (open_temp_dir(od, NULL, "define_test", define_dir)) {
    O3_ERROR_LOCATE(&(od->task));
    return CANNOT_WRITE_TEMP_FILE;
  }
  memset(&log_fd, 0, sizeof(FileDescriptor));
  memset(&coord_fd, 0, sizeof(FileDescriptor));
  memset(buffer, 0, BUF_LEN);
  memset(basis_set_name, 0, MAX_NAME_LEN);
  memset(lowercase_elem, 0, MAX_NAME_LEN);
  sprintf(coord_fd.name, "%s%ccoord", define_dir, SEPARATOR);
  if (!(coord_fd.handle = fopen(coord_fd.name, "wb"))) {
    O3_ERROR_LOCATE(&(od->task));
    O3_ERROR_STRING(&(od->task), coord_fd.name);
    return CANNOT_WRITE_TEMP_FILE;
  }
  fprintf(coord_fd.handle, "%s", methane);
  fclose(coord_fd.handle);
  memset(&prog_exe_info, 0, sizeof(ProgExeInfo));
  prog_exe_info.proc_env = fill_env
    (od, turbomole_env, od->field.qm_exe_path, 0);
  if (!(prog_exe_info.proc_env)) {
    O3_ERROR_LOCATE(&(od->task));
    return FL_OUT_OF_MEMORY;
  }
  prog_exe_info.need_stdin = 1;
  prog_exe_info.stdout_fd = &log_fd;
  prog_exe_info.stderr_fd = &log_fd;
  prog_exe_info.exedir = define_dir;
  prog_exe_info.sep_proc_grp = 1;
  sprintf(log_fd.name, "%s%cdefine.log", define_dir, SEPARATOR);
  sprintf(prog_exe_info.command_line, "%s%c%s",
    od->field.qm_exe_path, SEPARATOR, TURBOMOLE_DEFINE_EXE);
  pid = ext_program_exe(&prog_exe_info, &error);
  if (error) {
    O3_ERROR_LOCATE(&(od->task));
    return error;
  }
  #ifndef WIN32
  if (!(pipe_handle = fdopen(prog_exe_info.pipe_des[1], "w"))) {
    O3_ERROR_LOCATE(&(od->task));
    return FL_CANNOT_CREATE_CHANNELS;
  }
  #else
  pipe_handle = prog_exe_info.stdin_wr;
  #endif
  /*
  pipe job information into define
  */
  FWRITE_WRAP(pipe_handle,
    "\n"
    "define_test\n"
    "a coord\n"
    "*\n"
    "no\n"
    "b all def-SVP\n"
    "*\n"
    "*\n"
    "*\n", &n_chr);
  FFLUSH_WRAP(pipe_handle);
  FCLOSE_WRAP(pipe_handle);
  ext_program_wait(&prog_exe_info, pid);
  if (!(log_fd.handle = fopen(log_fd.name, "rb"))) {
    O3_ERROR_LOCATE(&(od->task));
    O3_ERROR_STRING(&(od->task), log_fd.name);
    return FL_CANNOT_READ_TEMP_FILE;
  }
  if (!fgrep(log_fd.handle, buffer, TURBOMOLE_NORMAL_TERMINATION)) {
    O3_ERROR_LOCATE(&(od->task));
    error = FL_ABNORMAL_TERMINATION;
  }
  fclose(log_fd.handle);
  if (error) {
    return error;
  }
  remove_recursive(define_dir);

  return 0;
}
