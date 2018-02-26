/*

close_files.c

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


void remove_recursive(char *filename)
{
  char buffer[BUF_LEN];
  char *nextfile;
  #ifndef WIN32
  DIR *dir = NULL;
  struct dirent dir_entry;
  struct dirent *result;
  #else
  HANDLE dir;
  WIN32_FIND_DATA filedata;
  #endif
  
  
  #ifndef WIN32
  dir = opendir(filename);
  if (dir) {
    while (!readdir_r(dir, &dir_entry, &result)) {
      if (!result) {
        break;
      }
      nextfile = dir_entry.d_name;
  #else
  sprintf(buffer, "%s\\*", filename);
  dir = FindFirstFileA(buffer, &filedata);
  if (dir != INVALID_HANDLE_VALUE) {
    while (FindNextFileA(dir, &filedata)) {
      nextfile = filedata.cFileName;
  #endif
      if ((!strcmp(nextfile, "."))
        || (!strcmp(nextfile, ".."))) {
        continue;
      }
      snprintf(buffer, BUF_LEN, "%s%c%s",
        filename, SEPARATOR, nextfile);
      if (unlink(buffer)) {
        remove_recursive(buffer);
      }
    }
    #ifndef WIN32
    closedir(dir);
    #else
    FindClose(dir);
    #endif
  }
  rmdir(filename);
}


void remove_temp_files(char *package_code)
{
  char *temp_dir_string;
  char *nextfile;
  char prefix[BUF_LEN];
  char buffer[BUF_LEN];
  FILE *clean_up_handle = NULL;
  #ifndef WIN32
  DIR *dir = NULL;
  struct dirent dir_entry;
  struct dirent *result;
  #else
  HANDLE dir;
  WIN32_FIND_DATA filedata;
  #endif
  
  
  memset(prefix, 0, BUF_LEN);
  temp_dir_string = getenv(TEMP_DIR_ENV);
  if (temp_dir_string) {
    #ifndef WIN32
    dir = opendir(temp_dir_string);
    if (dir) {
      while (!readdir_r(dir, &dir_entry, &result)) {
        if (!result) {
          break;
        }
        nextfile = dir_entry.d_name;
    #else
    sprintf(buffer, "%s\\*", temp_dir_string);
    dir = FindFirstFileA(buffer, &filedata);
    if (dir != INVALID_HANDLE_VALUE) {
      while (FindNextFileA(dir, &filedata)) {
        nextfile = filedata.cFileName;
    #endif
        if (strstr(nextfile, CLEAN_UP_SUFFIX)) {
          strcpy(prefix, nextfile);
          prefix[strlen(prefix) - strlen(CLEAN_UP_SUFFIX)] = '\0';
          if (!remove_with_prefix(temp_dir_string, prefix)) {
            sprintf(buffer, "%s%c%s", temp_dir_string,
              SEPARATOR, prefix);
            remove(buffer);
          }
        }
      }
      #ifndef WIN32
      closedir(dir);
      #else
      FindClose(dir);
      #endif
    }
    sprintf(prefix, "%s.%d", package_code, (int)getpid());
    if (remove_with_prefix(temp_dir_string, prefix)) {
      /*
      if some files could not be removed
      (presumably because they were in use)
      touch a clean up file which will instruct
      us to clean up the temporary folder
      next time we will be executed
      */
      sprintf(buffer, "%s%c%s.%d%s", temp_dir_string,
        SEPARATOR, package_code, (int)getpid(),
        CLEAN_UP_SUFFIX);
      if ((clean_up_handle = fopen(buffer, "wb+"))) {
        fclose(clean_up_handle);
      }
    }
  }
}


int remove_with_prefix(char *temp_dir_string, char *prefix)
{
  char *nextfile;
  char buffer[BUF_LEN];
  int n = 1;
  int attempts = 0;
  #ifndef WIN32
  DIR *dir = NULL;
  struct dirent dir_entry;
  struct dirent *result;
  #else
  HANDLE dir;
  WIN32_FIND_DATA filedata;
  #endif


  while (n && (attempts < MAX_ATTEMPTS_FILE)) {
    n = 0;
    #ifndef WIN32
    dir = opendir(temp_dir_string);
    if (dir) {
      while (!readdir_r(dir, &dir_entry, &result)) {
        if (!result) {
          break;
        }
        nextfile = dir_entry.d_name;
    #else
    sprintf(buffer, "%s\\*", temp_dir_string);
    dir = FindFirstFileA(buffer, &filedata);
    if (dir != INVALID_HANDLE_VALUE) {
      while (FindNextFileA(dir, &filedata)) {
        nextfile = filedata.cFileName;
    #endif
        if ((!strcmp(nextfile, "."))
          || (!strcmp(nextfile, ".."))) {
          continue;
        }
        if (strstr(nextfile, prefix)) {
          snprintf(buffer, BUF_LEN, "%s%c%s",
            temp_dir_string, SEPARATOR, nextfile);
          if (unlink(buffer)) {
            remove_recursive(buffer);
          }
          ++n;
        }
      }
      #ifndef WIN32
      closedir(dir);
      #else
      FindClose(dir);
      #endif
    }
    if (n) {
      #ifndef WIN32
      usleep(1000 * MS_SLEEP_BEFORE_RETRY);
      #else
      Sleep(MS_SLEEP_BEFORE_RETRY);
      #endif
      ++attempts;
    }
  }
  
  return n;
}


void close_files(O3Data *od, int from)
{
  int i;
  

  for (i = from; i < od->file_num; ++i) {
    if (od->file[i]->handle) {
      fclose(od->file[i]->handle);
      od->file[i]->handle = NULL;
    }
    if (i >= TEMP_START) {
      if (od->file[i]->name) {
        if (od->file[i]->name[0]) {
          remove(od->file[i]->name);
        }
      }
    }
    if (od->mel.file_descriptor[i]) {
      free(od->mel.file_descriptor[i]);
      od->mel.file_descriptor[i] = NULL;
    }
  }
  od->file_num = from;
}
