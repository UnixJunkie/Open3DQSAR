/*

calc_field.c

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
#include <include/proc_env.h>
#include <include/ff_parm.h>
#include <include/prog_exe_info.h>


extern EnvList minimal_env[];
extern EnvList gaussian_env[];
extern EnvList gamess_env[];
extern EnvList turbomole_env[];


int calc_field(O3Data *od, void *thread_func, int prep_or_calc)
{
  char buffer[BUF_LEN];
  char buffer2[BUF_LEN];
  int i;
  int n_threads;
  #ifndef WIN32
  pthread_attr_t thread_attr;
  #endif
  ThreadInfo **ti;


  ti = od->mel.thread_info;
  memset(buffer, 0, BUF_LEN);
  if (alloc_threads(od)) {
    return OUT_OF_MEMORY;
  }
  if (!(od->al.task_list = (TaskInfo **)alloc_array
    (od->grid.object_num, sizeof(TaskInfo)))) {
    return OUT_OF_MEMORY;
  }
  for (i = 0; i < od->grid.object_num; ++i) {
    od->al.mol_info[i]->done = 0;
  }
  if ((void *)thread_func == (void *)calc_mm_thread) {
    /*
    allocate one more field
    */
    if (alloc_x_var_array(od, 1)) {
      return OUT_OF_MEMORY;
    }
  }
  else if ((void *)thread_func == (void *)calc_md_grid_thread) {
    /*
    copy GRUB.DAT from MD GRID directory
    to the temporary od_grid_dir directory
    */
    sprintf(buffer, "%s%c"GRUB_FILENAME, od->field.md_grid_exe_path, SEPARATOR);
    sprintf(buffer2, "%s%c"GRUB_FILENAME, od->field.md_grid_dir, SEPARATOR);
    if (!fcopy(buffer, buffer2, "wb")) {
      return CANNOT_COPY_GRUB_DAT;
    }
  }
  #ifndef WIN32
  pthread_mutex_init(od->mel.mutex, NULL);
  pthread_attr_init(&thread_attr);
  pthread_attr_setdetachstate(&thread_attr, PTHREAD_CREATE_JOINABLE);
  #else
  if (!(*(od->mel.mutex) = CreateMutex(NULL, FALSE, NULL))) {
    return CANNOT_CREATE_THREAD;
  }
  #endif
  n_threads = fill_thread_info(od, od->object_num);
  for (i = 0; i < n_threads; ++i) {
    memcpy(&(ti[i]->od), od, sizeof(O3Data));
    ti[i]->model_type = prep_or_calc;
    /*
    create the i-th thread
    */
    #ifndef WIN32
    od->error_code = pthread_create(&(od->thread_id[i]),
      &thread_attr, (void *(*)(void *))thread_func, ti[i]);
    if (od->error_code) {
      return CANNOT_CREATE_THREAD;
    }
    #else
    od->hThreadArray[i] = CreateThread(NULL, 0,
      (LPTHREAD_START_ROUTINE)thread_func,
      ti[i], 0, &(od->dwThreadIdArray[i]));
    if (!(od->hThreadArray[i])) {
      return CANNOT_CREATE_THREAD;
    }
    #endif
  }
  #ifndef WIN32
  /*
  free the pthread attribute memory
  */
  pthread_attr_destroy(&thread_attr);
  /*
  wait for all threads to have finished
  */
  for (i = 0; i < n_threads; ++i) {
    od->error_code = pthread_join(od->thread_id[i],
      &(od->thread_result[i]));
    if (od->error_code) {
      return CANNOT_JOIN_THREAD;
    }
  }
  pthread_mutex_destroy(od->mel.mutex);
  #else
  WaitForMultipleObjects(n_threads, od->hThreadArray, TRUE, INFINITE);
  for (i = 0; i < n_threads; ++i) {
    CloseHandle(od->hThreadArray[i]);
  }
  CloseHandle(*(od->mel.mutex));
  #endif
  free_threads(od);
  
  return 0;
}


#ifndef WIN32
void *calc_cosmo_thread(void *pointer)
#else
DWORD calc_cosmo_thread(void *pointer)
#endif
{
  char buffer[BUF_LEN];
  char input_dir[BUF_LEN];
  int len;
  int object_num;
  int pid = 0;
  int result;
  int assigned = 1;
  ProgExeInfo prog_exe_info;
  AtomInfo **atom;
  ThreadInfo *ti;
  FileDescriptor inp_fd;
  FileDescriptor out_fd;
  FileDescriptor log_fd;


  ti = (ThreadInfo *)pointer;
  memset(buffer, 0, BUF_LEN);
  memset(input_dir, 0, BUF_LEN);
  /*
  allocate memory for AtomInfo structure array
  */
  if (!(atom = (AtomInfo **)alloc_array(ti->od.field.max_n_atoms + 1, sizeof(AtomInfo)))) {
    for (object_num = ti->start; object_num <= ti->end; ++object_num) {
      O3_ERROR_LOCATE(ti->od.al.task_list[object_num]);
      ti->od.al.task_list[object_num]->code = FL_OUT_OF_MEMORY;
    }
    #ifndef WIN32
    pthread_exit(pointer);
    #else
    return 0;
    #endif
  }
  while (assigned) {
    object_num = 0;
    assigned = 0;
    while ((!assigned) && (object_num < ti->od.grid.object_num)) {
      if (!(ti->od.al.mol_info[object_num]->done)) {
        #ifndef WIN32
        pthread_mutex_lock(ti->od.mel.mutex);
        #else
        WaitForSingleObject(ti->od.mel.mutex, INFINITE);
        #endif
        if (!(ti->od.al.mol_info[object_num]->done)) {
          ti->od.al.mol_info[object_num]->done = 1;
          assigned = 1;
        }
        #ifndef WIN32
        pthread_mutex_unlock(ti->od.mel.mutex);
        #else
        ReleaseMutex(ti->od.mel.mutex);
        #endif
      }
      else {
        ++object_num;
      }
    }
    if (!assigned) {
      break;
    }
    sprintf(buffer, "%s%c%s_%04d"TURBOMOLE_COSMO_EXT, ti->od.field.qm_dir, SEPARATOR,
      cosmo_label, ti->od.al.mol_info[object_num]->object_id);
    if ((ti->model_type == 'c') && fexist(buffer)) {
      continue;
    }
    memset(&inp_fd, 0, sizeof(FileDescriptor));
    memset(&out_fd, 0, sizeof(FileDescriptor));
    memset(&log_fd, 0, sizeof(FileDescriptor));
    sprintf(inp_fd.name, "%s%c%s_%04d.inp", ti->od.field.qm_dir, SEPARATOR,
      cosmo_label, ti->od.al.mol_info[object_num]->object_id);
    sprintf(out_fd.name, "%s%c%s_%04d.out", ti->od.field.qm_dir, SEPARATOR,
      cosmo_label, ti->od.al.mol_info[object_num]->object_id);
    sprintf(log_fd.name, "%s%c%s_%04d.log", ti->od.field.qm_dir, SEPARATOR,
      cosmo_label, ti->od.al.mol_info[object_num]->object_id);
    ti->od.al.task_list[object_num]->code = prep_cosmo_input
      (&(ti->od), ti->od.al.task_list[object_num], atom, object_num);
    if (ti->od.al.task_list[object_num]->code) {
      continue;
    }
    if (ti->model_type == 'c') {
      prog_exe_info.proc_env = fill_env
        (&(ti->od), turbomole_env, ti->od.field.qm_exe_path, object_num);
      if (!(prog_exe_info.proc_env)) {
        O3_ERROR_LOCATE(ti->od.al.task_list[object_num]);
        ti->od.al.task_list[object_num]->code = FL_OUT_OF_MEMORY;
        continue;
      }
      sprintf(buffer, "%s%c%s_%04d.0", ti->od.field.qm_scratch, SEPARATOR,
        cosmo_label, ti->od.al.mol_info[object_num]->object_id);
      #ifndef WIN32
      result = mkdir(buffer, S_IRWXU | S_IRGRP | S_IROTH);
      #else
      result = mkdir(buffer);
      #endif
      if (result) {
        O3_ERROR_LOCATE(ti->od.al.task_list[object_num]);
        O3_ERROR_STRING(ti->od.al.task_list[object_num], buffer);
        ti->od.al.task_list[object_num]->code = FL_CANNOT_CREATE_SCRDIR;
        continue;
      }
      sprintf(input_dir, "%s%c%s_%04d", ti->od.field.qm_dir, SEPARATOR,
        cosmo_label, ti->od.al.mol_info[object_num]->object_id);
      prog_exe_info.need_stdin = 0;
      prog_exe_info.stdout_fd = &out_fd;
      prog_exe_info.stderr_fd = &log_fd;
      prog_exe_info.exedir = input_dir;
      prog_exe_info.sep_proc_grp = 1;
      sprintf(prog_exe_info.command_line, "%s%c"TURBOMOLE_RIDFT_EXE,
        ti->od.field.qm_exe_path, SEPARATOR);
      pid = ext_program_exe(&prog_exe_info, &(ti->od.al.task_list[object_num]->code));
      ext_program_wait(&prog_exe_info, pid);
      /*
      check if the calculation underwent normal termination
      */
      if (ti->od.al.task_list[object_num]->code) {
        O3_ERROR_LOCATE(ti->od.al.task_list[object_num]);
        continue;
      }
      if (!(out_fd.handle = fopen(out_fd.name, "rb"))) {
        O3_ERROR_LOCATE(ti->od.al.task_list[object_num]);
        O3_ERROR_STRING(ti->od.al.task_list[object_num], out_fd.name);
        ti->od.al.task_list[object_num]->code = FL_CANNOT_READ_OUT_FILE;
        continue;
      }
      if (!fgrep(out_fd.handle, buffer, TURBOMOLE_NORMAL_TERMINATION)) {
        O3_ERROR_LOCATE(ti->od.al.task_list[object_num]);
        ti->od.al.task_list[object_num]->code = FL_ABNORMAL_TERMINATION;
      }
      fclose(out_fd.handle);
      out_fd.handle = NULL;
      sprintf(buffer, "%s%ccontrol", input_dir, SEPARATOR);
      rename(buffer, inp_fd.name);
      sprintf(buffer, "%s%c%s_%04d"TURBOMOLE_COSMO_EXT, input_dir, SEPARATOR,
        cosmo_label, ti->od.al.mol_info[object_num]->object_id);
      len = strlen(inp_fd.name);
      strcpy(&(inp_fd.name[len - 4]), TURBOMOLE_COSMO_EXT);
      rename(buffer, inp_fd.name);
      free_proc_env(prog_exe_info.proc_env);
    }
  }
  #ifndef WIN32
  pthread_exit(pointer);
  #else
  return 0;
  #endif
}


int call_cs3d_program(O3Data *od)
{
  char buffer[BUF_LEN];
  int object_num;
  int error = 0;
  int pid;
  FileDescriptor inp_fd;
  FileDescriptor log_fd;
  FileDescriptor temp_fd;
  ProgExeInfo prog_exe_info;

  
  memset(buffer, 0, BUF_LEN);
  memset(&inp_fd, 0, sizeof(FileDescriptor));
  memset(&log_fd, 0, sizeof(FileDescriptor));
  memset(&temp_fd, 0, sizeof(FileDescriptor));
  sprintf(inp_fd.name, "%s%ccs3d.inp", od->field.qm_dir, SEPARATOR);
  sprintf(log_fd.name, "%s%ccs3d.log", od->field.qm_dir, SEPARATOR);
  prog_exe_info.proc_env = fill_env
    (od, turbomole_env, od->field.qm_exe_path, 0);
  if (!(prog_exe_info.proc_env)) {
    O3_ERROR_LOCATE(&(od->task));
    return OUT_OF_MEMORY;
  }
  strcpy(buffer, od->field.cs3d_exe);
  get_dirname(buffer);
  prog_exe_info.need_stdin = 0;
  prog_exe_info.stdout_fd = &log_fd;
  prog_exe_info.stderr_fd = &log_fd;
  prog_exe_info.exedir = buffer;
  prog_exe_info.sep_proc_grp = 1;
  sprintf(prog_exe_info.command_line, "%s %s",
    od->field.cs3d_exe, inp_fd.name);
  pid = ext_program_exe(&prog_exe_info, &error);
  ext_program_wait(&prog_exe_info, pid);
  free_proc_env(prog_exe_info.proc_env);
  /*
  check if the calculation underwent normal termination
  */
  if (!(log_fd.handle = fopen(log_fd.name, "rb"))) {
    O3_ERROR_LOCATE(&(od->task));
    O3_ERROR_STRING(&(od->task), log_fd.name);
    return CS3D_ERROR;
  }
  error = (fgrep(log_fd.handle, buffer, "COSMOsim3D stops after") ? 0 : CS3D_ERROR);
  fclose(log_fd.handle);
  if (!error) {
    sscanf(buffer, "%*s %*s %*s %d", &object_num);
    error = ((object_num == od->object_num) ? 0 : CS3D_ERROR);
  }
  if (error) {
    O3_ERROR_LOCATE(&(od->task));
  }
  
  return error;
}


#ifndef WIN32
void *calc_qm_thread(void *pointer)
#else
DWORD calc_qm_thread(void *pointer)
#endif
{
  char buffer[LARGE_BUF_LEN];
  char input_dir[BUF_LEN];
  int i;
  int len;
  int object_num;
  int pid = 0;
  int result;
  int assigned = 1;
  #ifndef WIN32
  FILE *cubegen_handle = NULL;
  #else
  DWORD n_chr = 0;
  HANDLE cubegen_handle;
  #endif
  ProgExeInfo prog_exe_info;
  FILE *handle = NULL;
  AtomInfo **atom;
  ThreadInfo *ti;
  FileDescriptor inp_fd;
  FileDescriptor out_fd;
  FileDescriptor log_fd;


  ti = (ThreadInfo *)pointer;
  memset(buffer, 0, LARGE_BUF_LEN);
  memset(input_dir, 0, BUF_LEN);
  /*
  allocate memory for AtomInfo structure array
  */
  if (!(atom = (AtomInfo **)alloc_array(ti->od.field.max_n_atoms + 1, sizeof(AtomInfo)))) {
    for (object_num = ti->start; object_num <= ti->end; ++object_num) {
      O3_ERROR_LOCATE(ti->od.al.task_list[object_num]);
      ti->od.al.task_list[object_num]->code = FL_OUT_OF_MEMORY;
    }
    #ifndef WIN32
    pthread_exit(pointer);
    #else
    return 0;
    #endif
  }
  while (assigned) {
    object_num = 0;
    assigned = 0;
    while ((!assigned) && (object_num < ti->od.grid.object_num)) {
      #ifndef WIN32
      pthread_mutex_lock(ti->od.mel.mutex);
      #else
      WaitForSingleObject(ti->od.mel.mutex, INFINITE);
      #endif
      if (!(ti->od.al.mol_info[object_num]->done)) {
        ti->od.al.mol_info[object_num]->done = 1;
        assigned = 1;
      }
      #ifndef WIN32
      pthread_mutex_unlock(ti->od.mel.mutex);
      #else
      ReleaseMutex(ti->od.mel.mutex);
      #endif
      if (!assigned) {
        ++object_num;
      }
    }
    if (!assigned)  {
      break;
    }
    memset(&inp_fd, 0, sizeof(FileDescriptor));
    memset(&out_fd, 0, sizeof(FileDescriptor));
    memset(&log_fd, 0, sizeof(FileDescriptor));
    sprintf(inp_fd.name, "%s%c%s_%04d.inp", ti->od.field.qm_dir, SEPARATOR,
      ti->od.field.qm_software, ti->od.al.mol_info[object_num]->object_id);
    sprintf(out_fd.name, "%s%c%s_%04d.out", ti->od.field.qm_dir, SEPARATOR,
      ti->od.field.qm_software, ti->od.al.mol_info[object_num]->object_id);
    sprintf(log_fd.name, "%s%c%s_%04d.log", ti->od.field.qm_dir, SEPARATOR,
      ti->od.field.qm_software, ti->od.al.mol_info[object_num]->object_id);
    ti->od.al.task_list[object_num]->code = prep_qm_input
      (&(ti->od), ti->od.al.task_list[object_num], atom, object_num);
    if (ti->od.al.task_list[object_num]->code || (ti->model_type == 'p')) {
      continue;
    }
    ti->od.al.task_list[object_num]->code = 0;
    if (ti->od.field.type & PREP_GAUSSIAN_INPUT) {
      sprintf(input_dir, "%s%c%s_%04d.0", ti->od.field.qm_scratch,
        SEPARATOR, ti->od.field.qm_software, ti->od.al.mol_info[object_num]->object_id);
      #ifndef WIN32
      result = mkdir(input_dir, S_IRWXU | S_IRGRP | S_IROTH);
      #else
      result = mkdir(input_dir);
      #endif
      if (result) {
        O3_ERROR_LOCATE(ti->od.al.task_list[object_num]);
        O3_ERROR_STRING(ti->od.al.task_list[object_num], input_dir);
        ti->od.al.task_list[object_num]->code = FL_CANNOT_CREATE_SCRDIR;
        continue;
      }
      memset(&prog_exe_info, 0, sizeof(ProgExeInfo));
      prog_exe_info.proc_env = fill_env
        (&(ti->od), gaussian_env, ti->od.field.qm_exe_path, object_num);
      if (!(prog_exe_info.proc_env)) {
        O3_ERROR_LOCATE(ti->od.al.task_list[object_num]);
        ti->od.al.task_list[object_num]->code = FL_OUT_OF_MEMORY;
        continue;
      }
      prog_exe_info.need_stdin = 0;
      prog_exe_info.stdout_fd = &log_fd;
      prog_exe_info.stderr_fd = &log_fd;
      prog_exe_info.exedir = input_dir;
      prog_exe_info.sep_proc_grp = 1;
      sprintf(prog_exe_info.command_line, "%s%c%s %s %s",
        ti->od.field.qm_exe_path, SEPARATOR, ti->od.field.qm_exe,
        inp_fd.name, out_fd.name);
      pid = ext_program_exe(&prog_exe_info, &(ti->od.al.task_list[object_num]->code));
      ext_program_wait(&prog_exe_info, pid);
      if (ti->od.al.task_list[object_num]->code) {
        O3_ERROR_LOCATE(ti->od.al.task_list[object_num]);
        continue;
      }
      if (!(out_fd.handle = fopen(out_fd.name, "rb"))) {
        O3_ERROR_LOCATE(ti->od.al.task_list[object_num]);
        O3_ERROR_STRING(ti->od.al.task_list[object_num], out_fd.name);
        ti->od.al.task_list[object_num]->code = FL_CANNOT_READ_OUT_FILE;
        continue;
      }
      if (!fgrep(out_fd.handle, buffer, GAUSSIAN_NORMAL_TERMINATION)) {
        O3_ERROR_LOCATE(ti->od.al.task_list[object_num]);
        ti->od.al.task_list[object_num]->code = FL_ABNORMAL_TERMINATION;
      }
      if (out_fd.handle) {
        fclose(out_fd.handle);
        out_fd.handle = NULL;
      }
      if (ti->od.al.task_list[object_num]->code) {
        continue;
      }
      /*
      get formchk arguments from the input file
      */
      if (!(inp_fd.handle = fopen(inp_fd.name, "rb"))) {
        O3_ERROR_LOCATE(ti->od.al.task_list[object_num]);
        O3_ERROR_STRING(ti->od.al.task_list[object_num], inp_fd.name);
        ti->od.al.task_list[object_num]->code = FL_CANNOT_READ_INP_FILE;
        continue;
      }
      if (!fgrep(inp_fd.handle, buffer, "! "FORMCHK_EXE)) {
        O3_ERROR_LOCATE(ti->od.al.task_list[object_num]);
        O3_ERROR_STRING(ti->od.al.task_list[object_num], inp_fd.name);
        ti->od.al.task_list[object_num]->code = FL_CANNOT_READ_INP_FILE;
      }
      if (inp_fd.handle) {
        fclose(inp_fd.handle);
        inp_fd.handle = NULL;
      }
      if (ti->od.al.task_list[object_num]->code) {
        continue;
      }
      remove_newline(buffer);
      /*
      run formchk
      */
      prog_exe_info.need_stdin = 0;
      prog_exe_info.stdout_fd = &log_fd;
      prog_exe_info.stderr_fd = &log_fd;
      prog_exe_info.exedir = ti->od.field.qm_dir;
      prog_exe_info.sep_proc_grp = 1;
      sprintf(prog_exe_info.command_line, "%s%c%s",
        ti->od.field.qm_exe_path, SEPARATOR, &buffer[2]);
      pid = ext_program_exe(&prog_exe_info, &(ti->od.al.task_list[object_num]->code));
      ext_program_wait(&prog_exe_info, pid);
      /*
      check that the .fchk file exists and is readable
      */
      if (ti->od.al.task_list[object_num]->code) {
        O3_ERROR_LOCATE(ti->od.al.task_list[object_num]);
        continue;
      }
      sprintf(buffer, "%s%c%s_%04d"FORMCHK_EXT, ti->od.field.qm_dir, SEPARATOR,
        ti->od.field.qm_software, ti->od.al.mol_info[object_num]->object_id);
      if (!(handle = fopen(buffer, "rb"))) {
        O3_ERROR_LOCATE(ti->od.al.task_list[object_num]);
        O3_ERROR_STRING(ti->od.al.task_list[object_num], buffer);
        ti->od.al.task_list[object_num]->code = FL_CANNOT_READ_FCHK_FILE;
        continue;
      }
      if (handle) {
        fclose(handle);
        handle = NULL;
      }
      /*
      get cubegen arguments from the input file
      */
      if (!(inp_fd.handle = fopen(inp_fd.name, "rb"))) {
        O3_ERROR_LOCATE(ti->od.al.task_list[object_num]);
        O3_ERROR_STRING(ti->od.al.task_list[object_num], inp_fd.name);
        ti->od.al.task_list[object_num]->code = FL_CANNOT_READ_INP_FILE;
        continue;
      }
      if (!fgrep(inp_fd.handle, buffer, "! "CUBEGEN_EXE)) {
        O3_ERROR_LOCATE(ti->od.al.task_list[object_num]);
        O3_ERROR_STRING(ti->od.al.task_list[object_num], inp_fd.name);
        ti->od.al.task_list[object_num]->code = FL_CANNOT_READ_INP_FILE;
        fclose(inp_fd.handle);
        inp_fd.handle = NULL;
        continue;
      }
      remove_newline(buffer);
      prog_exe_info.need_stdin = 1;
      prog_exe_info.stdout_fd = &log_fd;
      prog_exe_info.stderr_fd = &log_fd;
      prog_exe_info.exedir = ti->od.field.qm_dir;
      prog_exe_info.sep_proc_grp = 1;
      sprintf(prog_exe_info.command_line, "%s%c%s",
        ti->od.field.qm_exe_path, SEPARATOR, &buffer[2]);
      pid = ext_program_exe(&prog_exe_info, &(ti->od.al.task_list[object_num]->code));
      if (ti->od.al.task_list[object_num]->code) {
        O3_ERROR_LOCATE(ti->od.al.task_list[object_num]);
        continue;
      }
      #ifndef WIN32
      if (!(cubegen_handle = fdopen(prog_exe_info.pipe_des[1], "w"))) {
        O3_ERROR_LOCATE(ti->od.al.task_list[object_num]);
        ti->od.al.task_list[object_num]->code = FL_CANNOT_CREATE_CHANNELS;
        continue;
      }
      #else
      cubegen_handle = prog_exe_info.stdin_wr;
      #endif
      /*
      pipe grid information into cubegen
      */
      i = 0;
      fgrep(inp_fd.handle, buffer, "! "CUBEGEN_EXE);
      while (fgets(buffer, BUF_LEN, inp_fd.handle) && (i < 4)) {
        buffer[BUF_LEN - 1] = '\0';
        FWRITE_WRAP(cubegen_handle, &buffer[2], &n_chr);
        ++i;
      }
      FFLUSH_WRAP(cubegen_handle);
      FCLOSE_WRAP(cubegen_handle);
      fclose(inp_fd.handle);
      inp_fd.handle = NULL;
      /*
      check that the .gcube file exists and is readable
      */
      ext_program_wait(&prog_exe_info, pid);
      sprintf(buffer, "%s%c%s_%04d"GAUSSIAN_CUBE_EXT, ti->od.field.qm_dir, SEPARATOR,
        ti->od.field.qm_software, ti->od.al.mol_info[object_num]->object_id);
      if (!(handle = fopen(buffer, "rb"))) {
        O3_ERROR_LOCATE(ti->od.al.task_list[object_num]);
        O3_ERROR_STRING(ti->od.al.task_list[object_num], buffer);
        ti->od.al.task_list[object_num]->code = FL_CANNOT_READ_GCUBE_FILE;
      }
      if (handle) {
        fclose(handle);
      }
    }
    else if (ti->od.field.type & PREP_FIREFLY_INPUT) {
      prog_exe_info.proc_env = fill_env
        (&(ti->od), minimal_env, ti->od.field.qm_exe_path, object_num);
      if (!(prog_exe_info.proc_env)) {
        O3_ERROR_LOCATE(ti->od.al.task_list[object_num]);
        ti->od.al.task_list[object_num]->code = FL_OUT_OF_MEMORY;
        continue;
      }
      prog_exe_info.need_stdin = 0;
      prog_exe_info.stdout_fd = &log_fd;
      prog_exe_info.stderr_fd = &log_fd;
      prog_exe_info.exedir = ti->od.field.qm_exe_path;
      prog_exe_info.sep_proc_grp = 1;
      #ifdef linux
      sprintf(prog_exe_info.command_line, "%s%c%s "FIREFLY_OPTIONS
        " -nompi -ncores 1 -nthreads 1 -ex %s -t %s%c%s_%04d -i %s -o %s ",
        ti->od.field.qm_exe_path, SEPARATOR, ti->od.field.qm_exe,
        ti->od.field.qm_exe_path, ti->od.field.qm_dir, SEPARATOR,
        ti->od.field.qm_software, ti->od.al.mol_info[object_num]->object_id,
        inp_fd.name, out_fd.name);
      #elif __APPLE__
      sprintf(prog_exe_info.command_line, "%s%cWINE%cbin%cwine %s%c%s -osx "FIREFLY_OPTIONS
        " -t %s%c%s_%04d -i %s -o %s -np 1",
        ti->od.field.qm_exe_path, SEPARATOR, SEPARATOR, SEPARATOR, 
        ti->od.field.qm_exe_path, SEPARATOR, ti->od.field.qm_exe,
        ti->od.field.qm_dir, SEPARATOR, ti->od.field.qm_software,
        ti->od.al.mol_info[object_num]->object_id, inp_fd.name, out_fd.name);
      #elif WIN32
      sprintf(prog_exe_info.command_line, "%s%c%s "FIREFLY_OPTIONS
        " -nompi -t %s%c%s_%04d -i %s -o %s ",
        ti->od.field.qm_exe_path, SEPARATOR, ti->od.field.qm_exe,
        ti->od.field.qm_dir, SEPARATOR, ti->od.field.qm_software,
        ti->od.al.mol_info[object_num]->object_id, inp_fd.name, out_fd.name);
      #else
      sprintf(prog_exe_info.command_line, "%s%c%s "FIREFLY_OPTIONS
        " -t %s%c%s_%04d -i %s -o %s "
        "-nompi -ncores 1 -nthreads 1",
        ti->od.field.qm_exe_path, SEPARATOR, ti->od.field.qm_exe,
        ti->od.field.qm_dir, SEPARATOR, ti->od.field.qm_software,
        ti->od.al.mol_info[object_num]->object_id, inp_fd.name, out_fd.name);
      #endif
      pid = ext_program_exe(&prog_exe_info, &(ti->od.al.task_list[object_num]->code));
      ext_program_wait(&prog_exe_info, pid);
      /*
      check if the calculation underwent normal termination
      */
      if (ti->od.al.task_list[object_num]->code) {
        O3_ERROR_LOCATE(ti->od.al.task_list[object_num]);
        continue;
      }
      if (!(out_fd.handle = fopen(out_fd.name, "rb"))) {
        O3_ERROR_LOCATE(ti->od.al.task_list[object_num]);
        O3_ERROR_STRING(ti->od.al.task_list[object_num], out_fd.name);
        ti->od.al.task_list[object_num]->code = FL_CANNOT_READ_OUT_FILE;
        continue;
      }
      if (!fgrep(out_fd.handle, buffer, FIREFLY_NORMAL_TERMINATION)) {
        O3_ERROR_LOCATE(ti->od.al.task_list[object_num]);
        ti->od.al.task_list[object_num]->code = FL_ABNORMAL_TERMINATION;
      }
      fclose(out_fd.handle);
      out_fd.handle = NULL;
    }
    else if (ti->od.field.type & PREP_GAMESS_INPUT) {
      prog_exe_info.proc_env = fill_env
        (&(ti->od), gamess_env, ti->od.field.qm_exe_path, object_num);
      if (!(prog_exe_info.proc_env)) {
        O3_ERROR_LOCATE(ti->od.al.task_list[object_num]);
        ti->od.al.task_list[object_num]->code = FL_OUT_OF_MEMORY;
        continue;
      }
      sprintf(buffer, "%s%c%s_%04d.0", ti->od.field.qm_scratch, SEPARATOR,
        ti->od.field.qm_software, ti->od.al.mol_info[object_num]->object_id);
      #ifndef WIN32
      result = mkdir(buffer, S_IRWXU | S_IRGRP | S_IROTH);
      #else
      result = mkdir(buffer);
      #endif
      if (result) {
        O3_ERROR_LOCATE(ti->od.al.task_list[object_num]);
        O3_ERROR_STRING(ti->od.al.task_list[object_num], buffer);
        ti->od.al.task_list[object_num]->code = FL_CANNOT_CREATE_SCRDIR;
        continue;
      }
      /*
      copy the .inp file in the qm_scratch
      dir with the .F05 extension
      */
      sprintf(buffer, "%s%c%s_%04d.0%c%s_%04d.F05", ti->od.field.qm_scratch,
        SEPARATOR, ti->od.field.qm_software, ti->od.al.mol_info[object_num]->object_id,
        SEPARATOR, ti->od.field.qm_software, ti->od.al.mol_info[object_num]->object_id);
      if (!fcopy(inp_fd.name, buffer, "wb")) {
        O3_ERROR_LOCATE(ti->od.al.task_list[object_num]);
        O3_ERROR_STRING(ti->od.al.task_list[object_num], buffer);
        ti->od.al.task_list[object_num]->code = FL_CANNOT_WRITE_INP_FILE;
        continue;
      }
      prog_exe_info.need_stdin = 0;
      prog_exe_info.stdout_fd = &out_fd;
      prog_exe_info.stderr_fd = &log_fd;
      prog_exe_info.exedir = ti->od.field.qm_exe_path;
      prog_exe_info.sep_proc_grp = 1;
      prog_exe_info.command_line[0] = '\0';
      #ifdef WIN32
      if (ti->od.field.mpiexec_exe[0]) {
        sprintf(inp_fd.name, "%s%c%s_%04d.0%cenvlist",
          ti->od.field.qm_scratch, SEPARATOR, ti->od.field.qm_software,
          ti->od.al.mol_info[object_num]->object_id, SEPARATOR);
        if (!(inp_fd.handle = fopen(inp_fd.name, "wb"))) {
          O3_ERROR_LOCATE(ti->od.al.task_list[object_num]);
          O3_ERROR_STRING(ti->od.al.task_list[object_num], inp_fd.name);
          ti->od.al.task_list[object_num]->code = FL_CANNOT_WRITE_TEMP_FILE;
          continue;
        }
        i = 0;
        while (prog_exe_info.proc_env[i]) {
          fprintf(inp_fd.handle, "%s\n", &(prog_exe_info.proc_env[i]));
          i += (strlen(&(prog_exe_info.proc_env[i])) + 1);
        }
        fclose(inp_fd.handle);
        inp_fd.handle = NULL;
        sprintf(prog_exe_info.command_line, "%s -env ENVFIL %s -n 2 %s%c%s",
          ti->od.field.mpiexec_exe, inp_fd.name, ti->od.field.qm_exe_path,
          SEPARATOR, ti->od.field.qm_exe);
      }
      #endif
      if (!(prog_exe_info.command_line[0])) {
        sprintf(prog_exe_info.command_line, "%s%c"GAMESS_DDIKICK_EXE" %s%c%s "
          "%s_%04d -ddi 1 1 127.0.0.1 -scr %s%c%s_%04d.0",
          ti->od.field.qm_exe_path, SEPARATOR,
          ti->od.field.qm_exe_path, SEPARATOR, ti->od.field.qm_exe,
          ti->od.field.qm_software, ti->od.al.mol_info[object_num]->object_id,
          ti->od.field.qm_scratch, SEPARATOR,
          ti->od.field.qm_software, ti->od.al.mol_info[object_num]->object_id);
      }
      pid = ext_program_exe(&prog_exe_info, &(ti->od.al.task_list[object_num]->code));
      ext_program_wait(&prog_exe_info, pid);
      /*
      check if the calculation underwent normal termination
      */
      if (ti->od.al.task_list[object_num]->code) {
        O3_ERROR_LOCATE(ti->od.al.task_list[object_num]);
        continue;
      }
      if (!(out_fd.handle = fopen(out_fd.name, "rb"))) {
        O3_ERROR_LOCATE(ti->od.al.task_list[object_num]);
        O3_ERROR_STRING(ti->od.al.task_list[object_num], out_fd.name);
        ti->od.al.task_list[object_num]->code = FL_CANNOT_READ_OUT_FILE;
        continue;
      }
      if (!fgrep(out_fd.handle, buffer, GAMESS_NORMAL_TERMINATION)) {
        O3_ERROR_LOCATE(ti->od.al.task_list[object_num]);
        ti->od.al.task_list[object_num]->code = FL_ABNORMAL_TERMINATION;
      }
      fclose(out_fd.handle);
      out_fd.handle = NULL;
    }
    else if (ti->od.field.type & PREP_TURBOMOLE_INPUT) {
      prog_exe_info.proc_env = fill_env
        (&(ti->od), turbomole_env, ti->od.field.qm_exe_path, object_num);
      if (!(prog_exe_info.proc_env)) {
        O3_ERROR_LOCATE(ti->od.al.task_list[object_num]);
        ti->od.al.task_list[object_num]->code = FL_OUT_OF_MEMORY;
        continue;
      }
      sprintf(buffer, "%s%c%s_%04d.0", ti->od.field.qm_scratch, SEPARATOR,
        ti->od.field.qm_software, ti->od.al.mol_info[object_num]->object_id);
      #ifndef WIN32
      result = mkdir(buffer, S_IRWXU | S_IRGRP | S_IROTH);
      #else
      result = mkdir(buffer);
      #endif
      if (result) {
        O3_ERROR_LOCATE(ti->od.al.task_list[object_num]);
        O3_ERROR_STRING(ti->od.al.task_list[object_num], buffer);
        ti->od.al.task_list[object_num]->code = FL_CANNOT_CREATE_SCRDIR;
        continue;
      }
      sprintf(input_dir, "%s%cturbomole_%04d", ti->od.field.qm_dir,
        SEPARATOR, ti->od.al.mol_info[object_num]->object_id);
      prog_exe_info.need_stdin = 0;
      prog_exe_info.stdout_fd = &out_fd;
      prog_exe_info.stderr_fd = &log_fd;
      prog_exe_info.exedir = input_dir;
      prog_exe_info.sep_proc_grp = 1;
      sprintf(prog_exe_info.command_line, "%s%c"TURBOMOLE_DSCF_EXE,
        ti->od.field.qm_exe_path, SEPARATOR);
      pid = ext_program_exe(&prog_exe_info, &(ti->od.al.task_list[object_num]->code));
      ext_program_wait(&prog_exe_info, pid);
      /*
      check if the calculation underwent normal termination
      */
      if (ti->od.al.task_list[object_num]->code) {
        O3_ERROR_LOCATE(ti->od.al.task_list[object_num]);
        continue;
      }
      if (!(out_fd.handle = fopen(out_fd.name, "rb"))) {
        O3_ERROR_LOCATE(ti->od.al.task_list[object_num]);
        O3_ERROR_STRING(ti->od.al.task_list[object_num], out_fd.name);
        ti->od.al.task_list[object_num]->code = FL_CANNOT_READ_OUT_FILE;
        continue;
      }
      if (!fgrep(out_fd.handle, buffer, TURBOMOLE_NORMAL_TERMINATION)) {
        O3_ERROR_LOCATE(ti->od.al.task_list[object_num]);
        ti->od.al.task_list[object_num]->code = FL_ABNORMAL_TERMINATION;
      }
      fclose(out_fd.handle);
      out_fd.handle = NULL;
      sprintf(buffer, "%s%ccontrol", input_dir, SEPARATOR);
      rename(buffer, inp_fd.name);
      sprintf(buffer, "%s%ct%c.cub", input_dir, SEPARATOR,
        (ti->od.field.type & QM_ELE_FIELD ? 'p' : 'd'));
      len = strlen(inp_fd.name);
      strcpy(&(inp_fd.name[len - 4]), GAMESS_PUNCH_EXT);
      rename(buffer, inp_fd.name);
    }
    free_proc_env(prog_exe_info.proc_env);
  }
  #ifndef WIN32
  pthread_exit(pointer);
  #else
  return 0;
  #endif
}


#ifndef WIN32
void *calc_mm_thread(void *pointer)
#else
DWORD calc_mm_thread(void *pointer)
#endif
{
  int i;
  int j;
  int n;
  int n_atoms;
  int result;
  int shift[3];
  int object_num;
  int assigned = 1;
  double R_i = 0.0;
  double R_j = 0.0;
  double R_ij2 = 0.0;
  double R_ij6 = 0.0;
  double gamma_ij = 0.0;
  double r_ij = 0.0;
  double r_ij2 = 0.0;
  double r_ij6 = 0.0;
  double r_ij7 = 0.0;
  double energy = 0.0;
  double probe_coord[3];
  VarCoord pc;
  AtomInfo **atom = NULL;
  ThreadInfo *ti;
  

  ti = (ThreadInfo *)pointer;
  R_j = ff_parm[O3_MMFF94][ti->od.field.probe.atom_type].vdw_parm[MMFF94_A]
    * pow(ff_parm[O3_MMFF94][ti->od.field.probe.atom_type].vdw_parm[MMFF94_ALPHA],
    MMFF94_POWER);
  memset(&pc, 0, sizeof(VarCoord));
  /*
  allocate memory for AtomInfo structure array
  */
  if (!(atom = (AtomInfo **)alloc_array(ti->od.field.max_n_atoms + 1, sizeof(AtomInfo)))) {
    for (object_num = ti->start; object_num <= ti->end; ++object_num) {
      O3_ERROR_LOCATE(ti->od.al.task_list[object_num]);
      ti->od.al.task_list[object_num]->code = FL_OUT_OF_MEMORY;
    }
    #ifndef WIN32
    pthread_exit(pointer);
    #else
    return 0;
    #endif
  }
  while (assigned) {
    object_num = 0;
    assigned = 0;
    while ((!assigned) && (object_num < ti->od.grid.object_num)) {
      #ifndef WIN32
      pthread_mutex_lock(ti->od.mel.mutex);
      #else
      WaitForSingleObject(ti->od.mel.mutex, INFINITE);
      #endif
      if (!(ti->od.al.mol_info[object_num]->done)) {
        ti->od.al.mol_info[object_num]->done = 1;
        assigned = 1;
      }
      #ifndef WIN32
      pthread_mutex_unlock(ti->od.mel.mutex);
      #else
      ReleaseMutex(ti->od.mel.mutex);
      #endif
      if (!assigned) {
        ++object_num;
      }
    }
    if (!assigned)  {
      break;
    }
    ti->od.al.task_list[object_num]->code =
      fill_atom_info(&(ti->od), ti->od.al.task_list[object_num],
      atom, NULL, object_num, ti->od.field.force_field);
    if (ti->od.al.task_list[object_num]->code) {
      continue;
    }
    n_atoms = ti->od.al.mol_info[object_num]->n_atoms;
    for (i = 0; i < n_atoms; ++i) {
      if (ti->od.field.type & VDW_FIELD) {
        R_i = ff_parm[O3_MMFF94][atom[i]->atom_type].vdw_parm[MMFF94_A]
          * pow(ff_parm[O3_MMFF94][atom[i]->atom_type].vdw_parm[MMFF94_ALPHA],
          MMFF94_POWER);
        gamma_ij = (R_i - R_j) / (R_i + R_j);
        atom[i]->parm[MMFF94_RIJ] =
          MMFF94_DAEPS * (R_i + R_j) * (1.0 + MMFF94_B
          * (1 - exp(-MMFF94_BETA * square(gamma_ij))));
        R_ij2 = square(atom[i]->parm[MMFF94_RIJ]);
        R_ij6 = R_ij2 * R_ij2 * R_ij2;
        atom[i]->parm[MMFF94_RIJ7] =
          R_ij6 * atom[i]->parm[MMFF94_RIJ];
        atom[i]->parm[MMFF94_EIJ] =
          181.16 * ff_parm[O3_MMFF94][atom[i]->atom_type].vdw_parm[MMFF94_G]
          * ff_parm[O3_MMFF94][ti->od.field.probe.atom_type].vdw_parm[MMFF94_G]
          * ff_parm[O3_MMFF94][atom[i]->atom_type].vdw_parm[MMFF94_ALPHA]
          * ff_parm[O3_MMFF94][ti->od.field.probe.atom_type].vdw_parm[MMFF94_ALPHA]
          / ((pow(ff_parm[O3_MMFF94][atom[i]->atom_type].vdw_parm[MMFF94_ALPHA]
          / ff_parm[O3_MMFF94][atom[i]->atom_type].vdw_parm[MMFF94_N], 0.5)
          + pow(ff_parm[O3_MMFF94][ti->od.field.probe.atom_type].vdw_parm[MMFF94_ALPHA]
          / ff_parm[O3_MMFF94][ti->od.field.probe.atom_type].vdw_parm[MMFF94_N], 0.5)) * R_ij6);
      }
    }
    result = 0;
    for (pc.node[2] = 0; (pc.node[2] < ti->od.grid.nodes[2]) && (!result); ++(pc.node[2])) {
      for (pc.node[1] = 0; (pc.node[1] < ti->od.grid.nodes[1]) && (!result); ++(pc.node[1])) {
        for (pc.node[0] = 0; (pc.node[0] < ti->od.grid.nodes[0]) && (!result); ++(pc.node[0])) {
          energy = 0.0;
          n = 0;
          for (j = ti->od.field.smooth_probe_flag; j <= 0; ++j) {
            for (shift[2] = j; shift[2] <= 1; shift[2] += 2) {
              for (shift[1] = j; shift[1] <= 1; shift[1] += 2) {
                for (shift[0] = j; shift[0] <= 1; shift[0] += 2) {
                  for (i = 0; i < 3; ++i) {
                    probe_coord[i] =
                      safe_rint(((pc.node[i] + (double)shift[i] / 3.0)
                      * (double)(ti->od.grid.step[i])) * 1.0e04) / 1.0e04
                      + (double)(ti->od.grid.start_coord[i]);
                  }
                  if (ti->od.field.type & VDW_FIELD) {
                    for (i = 0; i < n_atoms; ++i) {
                      r_ij2 = squared_euclidean_distance
                        (probe_coord, atom[i]->coord);
                      r_ij6 = r_ij2 * r_ij2 * r_ij2;
                      r_ij = sqrt(r_ij2);
                      r_ij7 = r_ij6 * r_ij;
                      energy += (atom[i]->parm[MMFF94_EIJ]
                        * pow(1.07 * atom[i]->parm[MMFF94_RIJ]
                        / (r_ij + 0.07 * atom[i]->parm[MMFF94_RIJ]), 7.0)
                        * (1.12 * atom[i]->parm[MMFF94_RIJ7]
                        / (r_ij7 + 0.12 * atom[i]->parm[MMFF94_RIJ]) - 2.0));
                    }
                  }
                  else {
                    for (i = 0; i < n_atoms; ++i) {
                      r_ij2 = squared_euclidean_distance
                        (probe_coord, atom[i]->coord);
                      if (ti->od.field.diel_dep == CONST_DIELECTRIC) {
                        r_ij2 = sqrt(r_ij2);
                      }
                      energy += (MMFF94_COUL * atom[i]->charge
                        / (ti->od.field.diel_const * (r_ij2 + MMFF94_ELEC_BUFF)));
                    }
                  }
                  ++n;
                }
              }
            }
          }
          /*
          if the probe_smooth option is used, then the energy value
          in each grid point is the average of the energy values
          on the vertexes of a cube with 2 * probe_smooth * step_size 
          side in angstrom, centred on the grid point, and in the
          centre of the cube itself
          */
          energy /= (double)n;
          result = set_x_value_xyz(&(ti->od), ti->od.field_num - 1, object_num, &pc, energy);
          if (result) {
            O3_ERROR_LOCATE(ti->od.al.task_list[object_num]);
            ti->od.al.task_list[object_num]->code = FL_OUT_OF_MEMORY;
          }
        }
      }
    }
  }
  free_array(atom);

  #ifndef WIN32
  pthread_exit(pointer);
  #else
  return 0;
  #endif
}


#ifndef WIN32
void *calc_md_grid_thread(void *pointer)
#else
DWORD calc_md_grid_thread(void *pointer)
#endif
{
  char buffer[BUF_LEN];
  int i;
  int n_atoms;
  int n;
  int object_num;
  int assigned = 1;
  int pid = 0;
  AtomInfo **atom = NULL;
  ProgExeInfo prog_exe_info;
  ThreadInfo *ti;
  FileDescriptor inp_fd;
  FileDescriptor log_fd;
  #ifndef WIN32
  FILE *pipe_handle = NULL;
  #else
  DWORD n_chr;
  HANDLE pipe_handle = NULL;
  #endif
  

  ti = (ThreadInfo *)pointer;
  memset(&inp_fd, 0, sizeof(FileDescriptor));
  memset(&log_fd, 0, sizeof(FileDescriptor));
  /*
  allocate memory for AtomInfo structure array
  */
  if (!(atom = (AtomInfo **)alloc_array(ti->od.field.max_n_atoms + 1, sizeof(AtomInfo)))) {
    for (object_num = ti->start; object_num <= ti->end; ++object_num) {
      O3_ERROR_LOCATE(ti->od.al.task_list[object_num]);
      ti->od.al.task_list[object_num]->code = FL_OUT_OF_MEMORY;
    }
    #ifndef WIN32
    pthread_exit(pointer);
    #else
    return 0;
    #endif
  }
  while (assigned) {
    object_num = 0;
    assigned = 0;
    while ((!assigned) && (object_num < ti->od.grid.object_num)) {
      #ifndef WIN32
      pthread_mutex_lock(ti->od.mel.mutex);
      #else
      WaitForSingleObject(ti->od.mel.mutex, INFINITE);
      #endif
      if (!(ti->od.al.mol_info[object_num]->done)) {
        ti->od.al.mol_info[object_num]->done = 1;
        assigned = 1;
      }
      #ifndef WIN32
      pthread_mutex_unlock(ti->od.mel.mutex);
      #else
      ReleaseMutex(ti->od.mel.mutex);
      #endif
      if (!assigned) {
        ++object_num;
      }
    }
    if (!assigned)  {
      break;
    }
    ti->od.al.task_list[object_num]->code =
      fill_atom_info(&(ti->od), ti->od.al.task_list[object_num],
      atom, NULL, object_num, O3_MMFF94);
    if (ti->od.al.task_list[object_num]->code) {
      continue;
    }
    ti->od.al.task_list[object_num]->code = fill_md_grid_types(atom);
    n_atoms = ti->od.al.mol_info[object_num]->n_atoms;
    sprintf(inp_fd.name, "%s%c%04d.pdb",
      ti->od.field.md_grid_dir, SEPARATOR, ti->od.al.mol_info[object_num]->object_id);
    if (!(inp_fd.handle = fopen(inp_fd.name, "wb+"))) {
      O3_ERROR_LOCATE(ti->od.al.task_list[object_num]);
      O3_ERROR_STRING(ti->od.al.task_list[object_num], inp_fd.name);
      ti->od.al.task_list[object_num]->code = FL_CANNOT_WRITE_TEMP_FILE;
      continue;
    }
    fprintf(inp_fd.handle, MD_GRID_PDB_FILE_HEADER"\n");
    for (i = 0; i < n_atoms; ++i) {
      fprintf(inp_fd.handle,
        "HETATM%5d %-4s <1>     1    %8.3lf%8.3lf%8.3lf  0.00  0.00\n",
        i + 1, atom[i]->atom_name, atom[i]->coord[0],
        atom[i]->coord[1], atom[i]->coord[2]);
    }
    fprintf(inp_fd.handle, "END\n");
    fclose(inp_fd.handle);
    inp_fd.handle = NULL;
    if (ti->od.al.task_list[object_num]->code) {
      continue;
    }

    sprintf(log_fd.name, "%s%c%04d.log", ti->od.field.md_grid_dir,
      SEPARATOR, ti->od.al.mol_info[object_num]->object_id);
    sprintf(inp_fd.name, "%s%c%04d.grin", ti->od.field.md_grid_dir,
      SEPARATOR, ti->od.al.mol_info[object_num]->object_id);
    if (!(inp_fd.handle = fopen(inp_fd.name, "wb+"))) {
      O3_ERROR_LOCATE(ti->od.al.task_list[object_num]);
      O3_ERROR_STRING(ti->od.al.task_list[object_num], inp_fd.name);
      ti->od.al.task_list[object_num]->code = FL_CANNOT_WRITE_TEMP_FILE;
      continue;
    }
    fprintf(inp_fd.handle,
      "LOUT %04d.lout\n"
      "KOUT %04d.kout\n"
      "INAT "GRUB_FILENAME"\n"
      "INKO %04d.pdb\n"
      "IEND\n",
      ti->od.al.mol_info[object_num]->object_id,
      ti->od.al.mol_info[object_num]->object_id,
      ti->od.al.mol_info[object_num]->object_id);
    fclose(inp_fd.handle);
    if (!(inp_fd.handle = fopen(inp_fd.name, "rb"))) {
      O3_ERROR_LOCATE(ti->od.al.task_list[object_num]);
      O3_ERROR_STRING(ti->od.al.task_list[object_num], inp_fd.name);
      ti->od.al.task_list[object_num]->code = FL_CANNOT_READ_TEMP_FILE;
      continue;
    }
    memset(&prog_exe_info, 0, sizeof(ProgExeInfo));
    prog_exe_info.proc_env = fill_env
      (&(ti->od), minimal_env, ti->od.field.md_grid_exe_path, object_num);
    if (!(prog_exe_info.proc_env)) {
      O3_ERROR_LOCATE(ti->od.al.task_list[object_num]);
      ti->od.al.task_list[object_num]->code = FL_OUT_OF_MEMORY;
      continue;
    }
    prog_exe_info.need_stdin = 1;
    prog_exe_info.stdout_fd = &log_fd;
    prog_exe_info.stderr_fd = &log_fd;
    prog_exe_info.exedir = ti->od.field.md_grid_dir;
    prog_exe_info.sep_proc_grp = 1;
    sprintf(prog_exe_info.command_line, "%s%c%s",
      ti->od.field.md_grid_exe_path, SEPARATOR, GRIN_EXE);
    pid = ext_program_exe(&prog_exe_info, &(ti->od.al.task_list[object_num]->code));
    if (!(ti->od.al.task_list[object_num]->code)) {
      #ifndef WIN32
      if (!(pipe_handle = fdopen(prog_exe_info.pipe_des[1], "w"))) {
        O3_ERROR_LOCATE(ti->od.al.task_list[object_num]);
        ti->od.al.task_list[object_num]->code = FL_CANNOT_CREATE_CHANNELS;
      }
      #else
      pipe_handle = prog_exe_info.stdin_wr;
      #endif
    }
    if (!(ti->od.al.task_list[object_num]->code)) {
      while (fgets(buffer, BUF_LEN, inp_fd.handle)) {
        buffer[BUF_LEN - 1] = '\0';
        FWRITE_WRAP(pipe_handle, buffer, &n_chr);
      }
      FFLUSH_WRAP(pipe_handle);
      FCLOSE_WRAP(pipe_handle);
    }
    if (inp_fd.handle) {
      fclose(inp_fd.handle);
      inp_fd.handle = NULL;
    }
    if (!(ti->od.al.task_list[object_num]->code)) {
      ext_program_wait(&prog_exe_info, pid);
      if (!(log_fd.handle = fopen(log_fd.name, "rb"))) {
        O3_ERROR_LOCATE(ti->od.al.task_list[object_num]);
        O3_ERROR_STRING(ti->od.al.task_list[object_num], log_fd.name);
        ti->od.al.task_list[object_num]->code = FL_CANNOT_READ_TEMP_FILE;
      }
    }
    if (!(ti->od.al.task_list[object_num]->code)) {
      n = 0;
      while ((!feof(log_fd.handle))
        && fgets(buffer, BUF_LEN, log_fd.handle)) {
        ++n;
      }
      if (n == 1) {
        n = (strncasecmp(buffer, "FORTRAN STOP", 12) ? 2 : 1);
      }
      if (n > 1) {
        O3_ERROR_LOCATE(ti->od.al.task_list[object_num]);
        ti->od.al.task_list[object_num]->code = FL_GRIN_ERROR;
      }
    }
    if (log_fd.handle) {
      fclose(log_fd.handle);
      log_fd.handle = NULL;
    }
    free_proc_env(prog_exe_info.proc_env);
  }
  free_array(atom);

  #ifndef WIN32
  pthread_exit(pointer);
  #else
  return 0;
  #endif
}  
  

int call_md_grid_program(O3Data *od)
{
  char buffer[BUF_LEN];
  char probe_name[MAX_FF_TYPE_LEN];
  int object_num;
  int n;
  int error;
  int pid = 0;
  int retry;
  FILE *handle;
  ProgExeInfo prog_exe_info;
  #ifndef WIN32
  FILE *pipe_handle = NULL;
  #else
  DWORD n_chr;
  HANDLE pipe_handle = NULL;
  #endif


  sprintf(buffer, "%s%cfile.list", od->field.md_grid_dir, SEPARATOR);
  if (!(handle = fopen(buffer, "wb+"))) {
    O3_ERROR_LOCATE(&(od->task));
    O3_ERROR_STRING(&(od->task), buffer);
    return CANNOT_WRITE_TEMP_FILE;
  }
  for (object_num = 0; object_num < od->grid.object_num; ++object_num) {
    fprintf(handle, "%04d\n", od->al.mol_info[object_num]->object_id);
  }
  fclose(handle);
  sprintf(buffer, "%s%cgrid.inp", od->field.md_grid_dir, SEPARATOR);
  if (!(handle = fopen(buffer, "wb+"))) {
    O3_ERROR_LOCATE(&(od->task));
    O3_ERROR_STRING(&(od->task), buffer);
    return CANNOT_WRITE_TEMP_FILE;
  }
  /*
  if the probe is CONHR or AR.CONHR GRID requires that information
  whether the hydrogen atom is cis or trans is provided; for this purpose,
  let's change the underscore in the probe name into a newline character
  */
  strcpy(probe_name, od->field.probe.atom_name);
  if ((od->field.probe.atom_type >= 59) && (od->field.probe.atom_type <= 62)) {
    n = 0;
    while (probe_name[n] && (probe_name[n] != '_')) {
      ++n;
    }
    if (probe_name[n]) {
      probe_name[n] = '\n';
    }
  }
  fprintf(handle,
    "LONT grid.lont\n"
    "KONT grid.kont\n"
    "INPT file.list\n"
    "BOTX %.4f\n"
    "BOTY %.4f\n"
    "BOTZ %.4f\n"
    "TOPX %.4f\n"
    "TOPY %.4f\n"
    "TOPZ %.4f\n"
    "NPLA %.4lf\n"
    "DWAT %.4lf\n"
    "EMAX %.4lf\n"
    "FARH 5.0\n"
    "FARR 8.0\n"
    "KWIK 0\n"
    "LEAU 0\n"
    "NUMB 0\n"
    "LENG 0\n"
    "LEVL 1\n"
    "MOVE 0\n"
    "LIST -2\n"
    "ALMD 0\n"
    "NETA 0\n"
    "%s\n"
    "IEND\n"
    PACKAGE_NAME" driven run\n"
    "0 1\n",
    od->grid.start_coord[0], od->grid.start_coord[1], od->grid.start_coord[2],
    od->grid.end_coord[0], od->grid.end_coord[1], od->grid.end_coord[2],
    (double)(od->grid.nodes[0] - 1) / ((double)(od->grid.end_coord[0])
    - (double)(od->grid.start_coord[0])), od->field.diel_const,
    od->field.md_grid_cutoff, probe_name);
  fclose(handle);
  retry = 0;
  while (retry < MAX_ATTEMPTS_FILE) {
    sprintf(buffer, "%s%cgrid.inp", od->field.md_grid_dir, SEPARATOR);
    if (!(handle = fopen(buffer, "rb"))) {
      O3_ERROR_LOCATE(&(od->task));
      O3_ERROR_STRING(&(od->task), buffer);
      return CANNOT_READ_TEMP_FILE;
    }
    memset(&prog_exe_info, 0, sizeof(ProgExeInfo));
    prog_exe_info.proc_env = fill_env
      (od, minimal_env, od->field.md_grid_exe_path, 0);
    if (!(prog_exe_info.proc_env)) {
      O3_ERROR_LOCATE(&(od->task));
      return OUT_OF_MEMORY;
    }
    prog_exe_info.need_stdin = 1;
    prog_exe_info.stdout_fd = od->file[TEMP_LOG];
    prog_exe_info.stderr_fd = od->file[TEMP_LOG];
    prog_exe_info.exedir = od->field.md_grid_dir;
    prog_exe_info.sep_proc_grp = 1;
    sprintf(prog_exe_info.command_line, "%s%c%s",
      od->field.md_grid_exe_path, SEPARATOR, GRID_EXE);
    pid = ext_program_exe(&prog_exe_info, &error);
    if (!error) {
      #ifndef WIN32
      if (!(pipe_handle = fdopen(prog_exe_info.pipe_des[1], "w"))) {
        O3_ERROR_LOCATE(&(od->task));
        return FL_CANNOT_CREATE_CHANNELS;
      }
      #else
      pipe_handle = prog_exe_info.stdin_wr;
      #endif
      while (fgets(buffer, BUF_LEN, handle)) {
        buffer[BUF_LEN - 1] = '\0';
        FWRITE_WRAP(pipe_handle, buffer, &n_chr);
      }
      FFLUSH_WRAP(pipe_handle);
      FCLOSE_WRAP(pipe_handle);
      fclose(handle);
      handle = NULL;
      ext_program_wait(&prog_exe_info, pid);
      if (!(od->file[TEMP_LOG]->handle = fopen
        (od->file[TEMP_LOG]->name, "rb"))) {
        O3_ERROR_LOCATE(&(od->task));
        O3_ERROR_STRING(&(od->task), od->file[TEMP_LOG]->name);
        return CANNOT_READ_TEMP_FILE;
      }
      error = 0;
      if (fgets(buffer, BUF_LEN, od->file[TEMP_LOG]->handle)) {
        buffer[BUF_LEN - 1] = '\0';
        if (strncasecmp(buffer, "FORTRAN STOP", 12)) {
          if (strstr(buffer, "date")) {
            /*
            this is a workaround for what appears to be a GRID bug:
            basically, GRID seems to write, and subsequently attempt
            to read, a file named /var/tmp/date0785. When working
            with multiple Open3DGRID processes running on the same
            multi-CPU machine, the file may just be erased/replaced
            by one of the running GRID processes, causing scattered,
            though unlikely, failures. The workaround is very simple:
            just retry until MAX_ATTEMPTS_FILE attempts have been made
            before acknowledging a failure of probably different nature
            */
            error = POSSIBLE_GRID_BUG;
          }
          else {
            error = GRID_ERROR;
          }
        }
      }
      fclose(od->file[TEMP_LOG]->handle);
      od->file[TEMP_LOG]->handle = NULL;
      switch (error) {
        case POSSIBLE_GRID_BUG:
        ++retry;
        break;

        case GRID_ERROR:
        O3_ERROR_LOCATE(&(od->task));
        return GRID_ERROR;

        default:
        retry = MAX_ATTEMPTS_FILE;
        break;
      }
    }
    if (handle) {
      fclose(handle);
    }
  }
  free_proc_env(prog_exe_info.proc_env);
  sprintf(buffer, "%s%c"GRUB_FILENAME, od->field.md_grid_dir, SEPARATOR);
  remove(buffer);

  return 0;
}
