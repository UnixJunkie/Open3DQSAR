/*

utils.c

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

#include <include/o3header.h>
#include <include/prog_exe_info.h>


int intlog2(int n)
{
  int i = 0;
  
  
  while (n >>= 1) ++i;
  
  return i;
}

  
int fcopy(char *from_filename, char *to_filename, char *mode)
{
  char buffer[LARGE_BUF_LEN];
  int eof = 0;
  int actual_len_read;
  int actual_len_write;
  int result;
  FILE *from_handle;
  FILE *to_handle;
  
  
  from_handle = fopen(from_filename, "rb");
  if (!from_handle) {
    return 0;
  }
  to_handle = fopen(to_filename, mode);
  if (!to_handle) {
    fclose(from_handle);
    return 0;
  }
  result = 1;
  while ((!eof) && result) {
    actual_len_read = fread(buffer, 1, LARGE_BUF_LEN, from_handle);
    eof = (actual_len_read < LARGE_BUF_LEN);
    actual_len_write = fwrite(buffer, 1, actual_len_read, to_handle);
    result = ((actual_len_write == actual_len_read) ? 1 : 0);
  }
  fclose(from_handle);
  fclose(to_handle);
  
  return result;
}


int fmove(char *filename1, char *filename2)
{
  if (!fcopy(filename1, filename2, "wb")) {
    return 0;
  }
  remove(filename1);
  
  return 1;
}


int fgrep(FILE *handle, char *buffer, char *grep_key)
{
  char last_occurrence[BUF_LEN];
  int found = 0;
  int line_n = 0;
  int last_occurrence_n = 0;
  

  rewind(handle);
  while (fgets(buffer, BUF_LEN, handle)) {
    ++line_n;
    buffer[BUF_LEN - 1] = '\0';
    if (strstr(buffer, grep_key)) {
      found = ftell(handle);
      memset(last_occurrence, 0, BUF_LEN);
      strcpy(last_occurrence, buffer);
      last_occurrence_n = line_n;
    }
  }
  if (found) {
    strcpy(buffer, last_occurrence);
    (void)fseek(handle, found, SEEK_SET);
  }

  return last_occurrence_n;
}


double squared_euclidean_distance(double *coord1, double *coord2)
{
  return (square(coord1[0] - coord2[0]) +
    square(coord1[1] - coord2[1]) +
    square(coord1[2] - coord2[2]));
}


void slash_to_backslash(char *string)
{
  int i = 0;
  
  
  if (string) {
    while (string[i]) {
      if (string[i] == '/') {
        string[i] = '\\';
      }
      ++i;
    }
  }
}


void remove_exe(char *string)
{
  int len;
  
  
  if (string) {
    len = strlen(string);
    if (len > 4) {
      if (!strncasecmp(&string[len - 4], ".exe", 4)) {
        string[len - 4] = '\0';
      }
    }
  }
}


void remove_extension(char *string)
{
  int len;
  int found = 0;
  
  
  if (string) {
    len = strlen(string);
    while (len && (!(found = (string[len - 1] == '.')))) {
      --len;
    }
    if (found) {
      string[len - 1] = '\0';
    }
  }
}


void remove_newline(char *string)
{
  int len;


  if (string) {
    len = strlen(string) - 1;
    while ((len >= 0) && strchr("\r\n", (int)string[len])) {
      string[len] = '\0';
      --len;
    }
  }
}


void string_to_lowercase(char *string)
{
  int i = 0;
  
  
  while (string[i]) {
    string[i] = (char)tolower((int)string[i]);
    ++i;
  }
}


void fix_endianness(void *chunk, int word_size, int chunk_len, int swap_endianness)
{

  /*
  This routine reverses word_size-byte chunks
  in a sequence of chunk_len bytes to comply
  with different binary formats
  */

  char temp_chunk[8];
  int i;
  int j;
  
  if (swap_endianness) {
    j = 0;
    while (j < (chunk_len * word_size)) {
      memcpy(temp_chunk, (char *)chunk + j, word_size);
      for (i = 0; i < word_size; ++i) {
        *((char *)chunk + j + i) =
          temp_chunk[word_size - i - 1];
      }
      j += word_size;
    }
  }
}


int machine_type()
{
  long one = 1;
  
  return (int)(!(*((char *)(&one))));
}


int ext_program_exe(ProgExeInfo *prog_exe_info, int *error)
{
  #ifndef WIN32
  char current_dir[BUF_LEN];
  char *execve_args[MAX_ARG];
  char *context = NULL;
  int i;
  int quote;
  int pid;
  int des = -1;
  int out = -1;
  int log = -1;
  int temp_des;
  struct rlimit rlp;
  #endif
  
  
  *error = 0;
  #ifndef WIN32
  memset(execve_args, 0, MAX_ARG * sizeof(char *));
  if (!getcwd(current_dir, BUF_LEN)) {
    current_dir[0] = '\0';
  }
  current_dir[BUF_LEN - 1] = '\0';
  if (prog_exe_info->need_stdin) {
    if (pipe(prog_exe_info->pipe_des)) {
      *error = FL_CANNOT_CREATE_CHANNELS;
      return 0;
    }
  }
  pid = fork();
  if (!pid) {
    /*
    we are the child process; close all inherited
    file descriptors
    */
    if (prog_exe_info->need_stdin) {
      /*
      close the write end of the pipe
      */
      close(prog_exe_info->pipe_des[1]);
    }
    if (prog_exe_info->sep_proc_grp) {
      setpgid(0, 0);
    }
    memset(&rlp, 0, sizeof(struct rlimit));
    if (!getrlimit(RLIMIT_NOFILE, &rlp)) {
      temp_des = (int)(rlp.rlim_cur);
      while (temp_des >= 0) {
        if ((!(prog_exe_info->need_stdin))
          || (temp_des != prog_exe_info->pipe_des[0])) {
          close(temp_des);
        }
        --temp_des;
      }
    }
    /*
    if needed, open a /dev/null file descriptor
    */
    if ((!(prog_exe_info->need_stdin)) || (!(prog_exe_info->stdout_fd))
      || (!(prog_exe_info->stderr_fd))) {
      des = open("/dev/null", O_RDWR);
      if (des == -1) {
        *error = FL_CANNOT_CREATE_CHANNELS;
      }
    }
    if (!(*error)) {
      if (prog_exe_info->need_stdin) {
        /*
        duplicate the read end of the pipe to
        CHILD_STDIN, then close it
        */
        dup2(prog_exe_info->pipe_des[0], STDIN_FILENO);
        close(prog_exe_info->pipe_des[0]);
      }
      else {
        dup2(des, STDIN_FILENO);
      }
      if (prog_exe_info->stdout_fd) {
        remove(prog_exe_info->stdout_fd->name);
        out = open
          (prog_exe_info->stdout_fd->name, O_CREAT | O_RDWR,
          S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
        if (out == -1) {
          *error = FL_CANNOT_CREATE_CHANNELS;
        }
        else {
          dup2(out, STDOUT_FILENO);
        }
      }
      else {
        dup2(des, STDOUT_FILENO);
      }
    }
    if (!(*error)) {
      if (prog_exe_info->stderr_fd) {
        if (prog_exe_info->stdout_fd != prog_exe_info->stderr_fd) {
          remove(prog_exe_info->stderr_fd->name);
          log = open(prog_exe_info->stderr_fd->name,
            O_CREAT | O_RDWR, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
          if (log == -1) {
            *error = FL_CANNOT_CREATE_CHANNELS;
          }
          else {
            dup2(log, STDERR_FILENO);
          }
        }
        else {
          dup2(out, STDERR_FILENO);
        }
      }
      else {
        dup2(des, STDERR_FILENO);
      }
    }
    if (!(*error)) {
      i = 0;
      quote = 0;
      while (prog_exe_info->command_line[i]) {
        if (prog_exe_info->command_line[i] == '\"') {
          prog_exe_info->command_line[i] = '\t';
          quote = 1 - quote;
        }
        else if ((prog_exe_info->command_line[i] == ' ') && (!quote)) {
          prog_exe_info->command_line[i] = '\t';
        }
        ++i;
      }
      i = 0;
      execve_args[i] = strtok_r(prog_exe_info->command_line, "\t\n\r\0", &context);
      while (execve_args[i] && (i < (MAX_ARG - 1))) {
        ++i;
        execve_args[i] = strtok_r(NULL, "\t\n\r\0", &context);
      }
      if (chdir(prog_exe_info->exedir)) {
        *error = FL_CANNOT_CHDIR;
      }
      if (!(*error)) {
        if (prog_exe_info->proc_env) {
          execve((const char *)execve_args[0],
            (char *const *)execve_args,
            (char *const *)(prog_exe_info->proc_env));
        }
        else {
          execv((const char *)execve_args[0],
            (char *const *)execve_args);
        }
        if (des >= 0) {
          close(des);
        }
        if (out >= 0) {
          close(out);
        }
        if (log >= 0) {
          close(log);
        }
      }
      if (chdir(current_dir)) {
        *error = FL_CANNOT_CHDIR;
      }
    }
    exit(0);
  }
  else if (pid > 0) {
    if (prog_exe_info->need_stdin
      && (prog_exe_info->need_stdin & NEED_STDIN_NORMAL)) {
      /*
      close the read end of the pipe
      */
      close(prog_exe_info->pipe_des[0]);
    }
    return pid;
  }
  
  return pid;
  #else
  memset(&(prog_exe_info->proc_info), 0, sizeof(PROCESS_INFORMATION));
  memset(&(prog_exe_info->startup_info), 0, sizeof(STARTUPINFO));
  memset(&(prog_exe_info->out_sa_attr), 0, sizeof(SECURITY_ATTRIBUTES));
  prog_exe_info->out_sa_attr.nLength = (DWORD)sizeof(SECURITY_ATTRIBUTES);
  prog_exe_info->out_sa_attr.bInheritHandle = TRUE;
  memcpy(&(prog_exe_info->pipe_sa_attr), &(prog_exe_info->out_sa_attr),
    sizeof(SECURITY_ATTRIBUTES));
  memcpy(&(prog_exe_info->log_sa_attr), &(prog_exe_info->out_sa_attr),
    sizeof(SECURITY_ATTRIBUTES));
  if ((!(prog_exe_info->need_stdin)) || (!(prog_exe_info->stdout_fd))
    || (!(prog_exe_info->stderr_fd))) {
    if ((prog_exe_info->des = CreateFile("nul",
      GENERIC_READ | GENERIC_WRITE,
      FILE_SHARE_READ | FILE_SHARE_WRITE, NULL,
      OPEN_EXISTING, 0, NULL)) == INVALID_HANDLE_VALUE) {
      *error = FL_CANNOT_CREATE_CHANNELS;
      return 0;
    }
  }
  if (prog_exe_info->need_stdin) {
    if (!CreatePipe(&(prog_exe_info->stdin_rd),
      &(prog_exe_info->stdin_wr),
      &(prog_exe_info->pipe_sa_attr), 0)) {
      *error = FL_CANNOT_CREATE_CHANNELS;
      return 0;
    }
    if (!(*error)) {
      if (!SetHandleInformation(prog_exe_info->stdin_wr,
        HANDLE_FLAG_INHERIT, 0)) {
        *error = FL_CANNOT_CREATE_CHANNELS;
        return 0;
      }
    }
  }
  else {
    prog_exe_info->stdin_rd = NULL;
    prog_exe_info->stdin_wr = NULL;
  }
  if (prog_exe_info->stdout_fd) {
    remove(prog_exe_info->stdout_fd->name);
    if ((prog_exe_info->out = CreateFile
      (prog_exe_info->stdout_fd->name,
      GENERIC_WRITE, FILE_SHARE_READ | FILE_SHARE_WRITE,
      &(prog_exe_info->out_sa_attr), CREATE_ALWAYS,
      FILE_ATTRIBUTE_NORMAL, NULL)) == INVALID_HANDLE_VALUE) {
      *error = FL_CANNOT_CREATE_CHANNELS;
      return 0;
    }
  }
  else {
    if (!DuplicateHandle(GetCurrentProcess(),
      prog_exe_info->des, GetCurrentProcess(),
      &(prog_exe_info->out), GENERIC_WRITE,
      TRUE, DUPLICATE_SAME_ACCESS)) {
      *error = FL_CANNOT_CREATE_CHANNELS;
      return 0;
    }
  }
  if (prog_exe_info->stdout_fd != prog_exe_info->stderr_fd) {
    if (prog_exe_info->stderr_fd) {
      remove(prog_exe_info->stderr_fd->name);
      if ((prog_exe_info->log = CreateFile
        (prog_exe_info->stderr_fd->name,
        GENERIC_WRITE, FILE_SHARE_READ | FILE_SHARE_WRITE,
        &(prog_exe_info->log_sa_attr), CREATE_ALWAYS,
        FILE_ATTRIBUTE_NORMAL, NULL)) == INVALID_HANDLE_VALUE) {
        *error = FL_CANNOT_CREATE_CHANNELS;
        return 0;
      }
    }
    else {
      if (!DuplicateHandle(GetCurrentProcess(), prog_exe_info->des,
        GetCurrentProcess(), &(prog_exe_info->log), GENERIC_WRITE,
        TRUE, DUPLICATE_SAME_ACCESS)) {
        *error = FL_CANNOT_CREATE_CHANNELS;
        return 0;
      }
    }
  }
  else {
    if (!DuplicateHandle(GetCurrentProcess(), prog_exe_info->out,
      GetCurrentProcess(), &(prog_exe_info->log), GENERIC_WRITE,
      TRUE, DUPLICATE_SAME_ACCESS)) {
      *error = FL_CANNOT_CREATE_CHANNELS;
      return 0;
    }
  }
  prog_exe_info->startup_info.cb = sizeof(STARTUPINFO);
  prog_exe_info->startup_info.hStdInput =
    (prog_exe_info->need_stdin ? prog_exe_info->stdin_rd : prog_exe_info->des);
  prog_exe_info->startup_info.hStdOutput = prog_exe_info->out;
  prog_exe_info->startup_info.hStdError = prog_exe_info->log;
  prog_exe_info->startup_info.dwFlags = STARTF_USESTDHANDLES;
  if (!CreateProcess(NULL, prog_exe_info->command_line,
    NULL, NULL, TRUE, CREATE_NO_WINDOW,
    (LPVOID)(prog_exe_info->proc_env),
    prog_exe_info->exedir, &(prog_exe_info->startup_info),
    &(prog_exe_info->proc_info))) {
    *error = FL_CANNOT_CREATE_PROCESS;
    return 0;
  }
  
  return 1;
  #endif
}


void ext_program_wait(ProgExeInfo *prog_exe_info, int pid)
{
  
  #ifndef WIN32
  int status;
  
  
  waitpid(pid, &status, 0);
  #else
  WaitForSingleObject(prog_exe_info->proc_info.hProcess, INFINITE);
  CloseHandle(prog_exe_info->proc_info.hProcess);
  CloseHandle(prog_exe_info->proc_info.hThread);
  if (prog_exe_info->des) {
    FlushFileBuffers(prog_exe_info->des);
    CloseHandle(prog_exe_info->des);
  }
  if (prog_exe_info->stdin_rd) {
    FlushFileBuffers(prog_exe_info->stdin_rd);
    CloseHandle(prog_exe_info->stdin_rd);
  }
  if (prog_exe_info->stdin_wr) {
    FlushFileBuffers(prog_exe_info->stdin_wr);
    CloseHandle(prog_exe_info->stdin_wr);
  }
  if (prog_exe_info->out) {
    FlushFileBuffers(prog_exe_info->out);
    CloseHandle(prog_exe_info->out);
  }
  if (prog_exe_info->log) {
    FlushFileBuffers(prog_exe_info->log);
    CloseHandle(prog_exe_info->log);
  }
  #endif
}


void print_debug_info(O3Data *od, TaskInfo *task)
{
  tee_printf(od, "DEBUG: error raised in %s (line %d, %s)\n",
    task->func, task->line, task->file);
}


int exe_shell_cmd(O3Data *od, char *command, char *exedir, char *shell)
{
  char buffer[BUF_LEN];
  char cwd[BUF_LEN];
  int result;
  int pid;
  ProgExeInfo prog_exe_info;
  struct stat file_stat;
  
  
  memset(buffer, 0, BUF_LEN);
  memset(cwd, 0, BUF_LEN);
  memset(&prog_exe_info, 0, sizeof(ProgExeInfo));
  prog_exe_info.stdout_fd = od->file[TEMP_OUT];
  prog_exe_info.stderr_fd = od->file[TEMP_LOG];
  if (!exedir) {
    exedir = getcwd(cwd, BUF_LEN - 1);
    if (!exedir) {
      return CANNOT_CHANGE_DIR;
    }
  }
  prog_exe_info.exedir = exedir;
  prog_exe_info.sep_proc_grp = 1;
  sprintf(prog_exe_info.command_line, "%s \"%s\"", shell
    ? shell : DEFAULT_SHELL_CMD, command);
  pid = ext_program_exe(&prog_exe_info, &result);
  if (result) {
    return result;
  }
  ext_program_wait(&prog_exe_info, pid);
  if (stat(od->file[TEMP_LOG]->name, &file_stat) == -1) {
    return CANNOT_READ_LOG_FILE;
  }
  if (file_stat.st_size) {
    return LOG_FILE_NOT_EMPTY;
  }
  if (stat(od->file[TEMP_OUT]->name, &file_stat) == -1) {
    return CANNOT_READ_OUT_FILE;
  }
  if (file_stat.st_size) {
    return OUT_FILE_NOT_EMPTY;
  }

  return 0;
  
}


#ifndef HAVE_STRTOK_R
char *strtok_r(char *s1, const char *s2, char **lasts)
{
  char *ret;

  if (s1 == NULL)
    s1 = *lasts;
  while ((*s1) && strchr(s2, *s1))
    ++s1;
  if (*s1 == '\0')
    return NULL;
  ret = s1;
  while ((*s1) && (!strchr(s2, *s1)))
    ++s1;
  if (*s1)
    *s1++ = '\0';
  *lasts = s1;

  return ret;
}
#endif


int fill_date_string(char *date_string)
{
  int ret = 0;
  #ifndef WIN32
  time_t current_time;
  struct tm broken_down_time;

  
  ret = ((time(&current_time) == -1) ? 0 : 1);
  if (ret) {
    ret = (localtime_r(&current_time, &broken_down_time) ? 1 : 0);
  }
  if (ret) {
    sprintf(date_string, "%04d%02d%02d_%02d%02d%02d",
      broken_down_time.tm_year + 1900,
      broken_down_time.tm_mon + 1,
      broken_down_time.tm_mday,
      broken_down_time.tm_hour,
      broken_down_time.tm_min,
      broken_down_time.tm_sec);
  }
  #else
  SYSTEMTIME lt;
  
  
  GetLocalTime(&lt);
  sprintf(date_string, "%04d%02d%02d_%02d%02d%02d",
    lt.wYear,
    lt.wMonth,
    lt.wDay,
    lt.wHour,
    lt.wMinute,
    lt.wSecond);
  ret = 1;
  #endif
  
  return ret;
}


#ifndef HAVE_MKDTEMP
#ifndef TMP_MAX
#define TMP_MAX      238328
#endif
#define  MAX_RND_LETTERS    62
#define  MAX_RND_PART    6
#define MKDTEMP_VALUE_INCR  7777


char *mkdtemp(char *tmpl)
{
  const char letters[] =
    "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789";
  int i;
  int j;
  int len;
  int fd;
  int old_errno;
  uint64_t random_time_bits;
  uint64_t value;
  uint64_t temp_value;
  struct timeval tv;

  
  len = strlen (tmpl);
  if ((len < MAX_RND_PART)
    || strcmp(&tmpl[len - MAX_RND_PART], "XXXXXX")) {
    errno = EINVAL;
    return NULL;
  }
  if (gettimeofday(&tv, NULL) == -1) {
    errno = ENOSYS;
    return NULL;
  }
  random_time_bits = ((uint64_t)(tv.tv_usec) << 16) ^ tv.tv_sec;
  value = random_time_bits ^ getpid();
  old_errno = errno;
  for (i = 0, fd = -1; (fd == -1) && (i < TMP_MAX);
    value += MKDTEMP_VALUE_INCR, ++i) {
    temp_value = value;
    for (j = 0; j < MAX_RND_PART; ++j) {
      tmpl[len - MAX_RND_PART + j] =
        letters[temp_value % MAX_RND_LETTERS];
      temp_value /= MAX_RND_LETTERS;
    }
    #ifndef WIN32
    fd = mkdir(tmpl, S_IRUSR | S_IWUSR | S_IXUSR);
    #else
    fd = mkdir(tmpl);
    #endif
    if ((fd == -1) && (errno != EEXIST)) {
      break;
    }
  }
  if (fd != -1) {
    errno = old_errno;
  }
  
  return ((fd == -1) ? NULL : tmpl);
}
#endif


#ifndef HAVE_MKSTEMP
#ifndef O_EXCL
#define O_EXCL      _O_EXCL
#endif
#ifndef O_RDWR
#define O_RDWR      _O_RDWR
#endif
#ifndef O_CREAT
#define O_CREAT      _O_CREAT
#endif
#ifndef O_BINARY
#ifndef _O_BINARY
#define _O_BINARY    0
#endif
#define O_BINARY    _O_BINARY
#endif
#ifndef TMP_MAX
#define TMP_MAX      238328
#endif
#define  MAX_RND_LETTERS    62
#define  MAX_RND_PART    6
#define MKSTEMP_VALUE_INCR  7777


int mkstemp(char *tmpl)
{
  const char letters[] =
    "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789";
  int i;
  int j;
  int len;
  int fd;
  int old_errno;
  uint64_t random_time_bits;
  uint64_t value;
  uint64_t temp_value;
  struct timeval tv;

  
  len = strlen (tmpl);
  if ((len < MAX_RND_PART)
    || strcmp(&tmpl[len - MAX_RND_PART], "XXXXXX")) {
    return -1;
  }
  if (gettimeofday(&tv, NULL) == -1) {
    return -1;
  }
  random_time_bits = ((uint64_t)(tv.tv_usec) << 16) ^ tv.tv_sec;
  value = random_time_bits ^ getpid();
  old_errno = errno;
  for (i = 0, fd = -1; (fd == -1) && (i < TMP_MAX);
    value += MKSTEMP_VALUE_INCR, ++i) {
    temp_value = value;
    for (j = 0; j < MAX_RND_PART; ++j) {
      tmpl[len - MAX_RND_PART + j] =
        letters[temp_value % MAX_RND_LETTERS];
      temp_value /= MAX_RND_LETTERS;
    }
    errno = 0;
    fd = open(tmpl,
      _O_EXCL | _O_RDWR | _O_BINARY | _O_CREAT, 0600);
    if ((fd == -1) && (errno != EEXIST)) {
      break;
    }
  }
  if (fd != -1) {
    errno = old_errno;
  }
  
  return fd;
}
#endif
