/*

fill_env.c

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
#ifdef WIN32
#include <windows.h>
#endif


#ifndef WIN32
char **fill_env(O3Data *od, EnvList personalized_env[], char *bin, int object_num)
#else
char *fill_env(O3Data *od, EnvList personalized_env[], char *bin, int object_num)
#endif
{
  #ifndef WIN32
  char **proc_env = NULL;
  #else
  char *curr_environ;
  char *proc_env = NULL;
  #endif
  char *equal_sign = NULL;
  char *temp_dir;
  char *babel_datadir;
  char *babel_libdir;
  char *string = NULL;
  int len;
  int max_len;
  int n_env = 0;
  int dont_append;
  int i = 0;
  int j = 0;
  int k = 0;
  int found = 0;
  #ifndef WIN32
  char *library_path;
  #else
  char *exe_path;
  #endif


  n_env = 0;
  while (personalized_env[n_env].name) {
    ++n_env;
  }
  #ifndef WIN32
  /*
  set personalized environment
  */
  i = 0;
  j = 0;
  /*
  count how many environment variables are there
  in the current environment, excluding those defined
  in our personalized environment, which we need to set
  according to our needs
  */
  max_len = 0;
  while (environ[i]) {
    len = strlen(environ[i]) + MAX_NAME_LEN;
    if (len > max_len) {
      max_len = len;
      string = realloc(string, len);
      if (!string) {
        return NULL;
      }
      memset(string, 0, len);
    }
    strcpy(string, environ[i]);
    if ((equal_sign = strchr(string, '='))) {
      *equal_sign = '\0';
    }
    k = 0;
    found = 0;
    while (personalized_env[k].name && (!found)) {
      found = (!strcmp(string, personalized_env[k].name));
      ++k;
    }
    if (!found) {
      ++j;
    }
    ++i;
  }
  proc_env = (char **)malloc((j + n_env + 1) * sizeof(char *));
  if (!proc_env) {
    free(string);
    return NULL;
  }
  memset(proc_env, 0, (j + n_env + 1) * sizeof(char *));
  i = 0;
  j = 0;
  while (environ[i]) {
    /*
    copy the current environment, excluding
    variables included in personalized_env
    */
    strcpy(string, environ[i]);
    if ((equal_sign = strchr(string, '='))) {
      *equal_sign = '\0';
    }
    k = 0;
    found = 0;
    while (personalized_env[k].name && (!found)) {
      found = (!strcmp(string, personalized_env[k].name));
      ++k;
    }
    if (!found) {
      len = strlen(environ[i]) + MAX_NAME_LEN;
      proc_env[j] = malloc(len);
      if (!(proc_env[j])) {
        free_array(proc_env);
        free(string);
        return NULL;
      }
      memset(proc_env[j], 0, len);
      strcpy(proc_env[j], environ[i]);
      ++j;
    }
    ++i;
  }
  #else
  /*
  set personalized environment
  */
  i = 0;
  j = 0;
  curr_environ = GetEnvironmentStrings();
  /*
  count how many environment variables are there
  in the current environment, excluding those defined
  in our personalized environment, which we need to set
  according to our needs
  */
  max_len = 0;
  while (curr_environ[i]) {
    len = strlen(&curr_environ[i]) + MAX_NAME_LEN;
    if (len > max_len) {
      max_len = len;
      string = realloc(string, len);
      if (!string) {
        return NULL;
      }
      memset(string, 0, len);
    }
    strcpy(string, &curr_environ[i]);
    if ((equal_sign = strchr(string, '='))) {
      *equal_sign = '\0';
    }
    k = 0;
    len = strlen(&curr_environ[i]) + MAX_NAME_LEN;
    found = 0;
    while (personalized_env[k].name && (!found)) {
      /*
      Windows environment variables are not
      case-sensitive
      */
      found = (!strcasecmp(string, personalized_env[k].name));
      ++k;
    }
    if (!found) {
      j += len;
    }
    i += len;
  }
  k = 0;
  len = 0;
  while (personalized_env[k].name) {
    len += BUF_LEN;
    ++k;
  }
  proc_env = malloc(j + len + MAX_NAME_LEN);
  if (!proc_env) {
    free(string);
    return NULL;
  }
  memset(proc_env, 0, j + len + MAX_NAME_LEN);
  /*
  copy the current environment, excluding
  variables included in personalized_env
  */
  i = 0;
  j = 0;
  while (curr_environ[i]) {
    strcpy(string, &curr_environ[i]);
    if ((equal_sign = strchr(string, '='))) {
      *equal_sign = '\0';
    }
    k = 0;
    len = strlen(&curr_environ[i]) + MAX_NAME_LEN;
    found = 0;
    while (personalized_env[k].name && (!found)) {
      /*
      Windows environment variables are not
      case-sensitive
      */
      found = (!strcasecmp(string, personalized_env[k].name));
      ++k;
    }
    if (!found) {
      strcpy(&proc_env[j], &curr_environ[i]);
      j += len;
    }
    i += len;
  }
  FreeEnvironmentStrings(curr_environ);
  #endif
  /*
  append to the current environment new values
  for selected variables
  */
  k = 0;
  while (personalized_env[k].name) {
    dont_append = 0;
    len = BUF_LEN;
    switch (personalized_env[k].flag) {
      case EL_QM_EXE_PATH:
      len = strlen(od->field.qm_exe_path) + strlen(personalized_env[k].ext) + BUF_LEN;
      if (len > max_len) {
        max_len = len;
        string = realloc(string, len);
        if (!string) {
          free(proc_env);
          return NULL;
        }
        memset(string, 0, len);
      }
      sprintf(string, "%s%c%s", od->field.qm_exe_path,
        SEPARATOR, personalized_env[k].ext);
      break;

      case EL_TEMP_DIR:
      temp_dir = getenv(TEMP_DIR_ENV);
      len = BUF_LEN;
      if (temp_dir) {
        len += strlen(temp_dir);
      }
      if (len > max_len) {
        max_len = len;
        string = realloc(string, len);
        if (!string) {
          free(proc_env);
          return NULL;
        }
        memset(string, 0, len);
      }
      if (temp_dir) {
        strcpy(string, temp_dir);
      }
      break;

      case EL_PLAIN_SET:
      len = strlen(personalized_env[k].ext) + BUF_LEN;
      if (len > max_len) {
        max_len = len;
        string = realloc(string, len);
        if (!string) {
          free(proc_env);
          return NULL;
        }
        memset(string, 0, len);
      }
      strcpy(string, personalized_env[k].ext);
      break;

      case EL_QM_PUNCH:
      len = strlen(od->field.qm_dir) + strlen(od->field.qm_software) + BUF_LEN;
      if (len > max_len) {
        max_len = len;
        string = realloc(string, len);
        if (!string) {
          free(proc_env);
          return NULL;
        }
        memset(string, 0, len);
      }
      sprintf(string, "%s%c%s_%04d"GAMESS_PUNCH_EXT,
        od->field.qm_dir, SEPARATOR, od->field.qm_software,
        od->al.mol_info[object_num]->object_id);
      break;

      case EL_QM_GAUSCR:
      len = strlen(od->field.qm_scratch) + strlen(od->field.qm_software) + BUF_LEN;
      if (len > max_len) {
        max_len = len;
        string = realloc(string, len);
        if (!string) {
          free(proc_env);
          return NULL;
        }
        memset(string, 0, len);
      }
      sprintf(string, "%s%c%s_%04d.0", od->field.qm_scratch,
        SEPARATOR, od->field.qm_software,
        od->al.mol_info[object_num]->object_id);
      break;
      
      case EL_QM_GMSSCR:
      len = strlen(od->field.qm_scratch) + 2 * strlen(od->field.qm_software) + BUF_LEN;
      if (len > max_len) {
        max_len = len;
        string = realloc(string, len);
        if (!string) {
          free(proc_env);
          return NULL;
        }
        memset(string, 0, len);
      }
      if (personalized_env[k].ext[0]) {
        sprintf(string, "%s%c%s_%04d.0%c%s_%04d%s", od->field.qm_scratch,
          SEPARATOR, od->field.qm_software, od->al.mol_info[object_num]->object_id,
          SEPARATOR, od->field.qm_software, od->al.mol_info[object_num]->object_id,
          personalized_env[k].ext);
      }
      else {
        sprintf(string, "%s%c%s_%04d.0", od->field.qm_scratch,
          SEPARATOR, od->field.qm_software, od->al.mol_info[object_num]->object_id);
      }
      break;

      #ifndef WIN32
      case EL_DYN_LIBRARY_PATH:
      library_path = getenv(DYN_LIBRARY_PATH);
      len = strlen(bin)  + BUF_LEN;
      if (library_path) {
        len += strlen(library_path);
      }
      if (len > max_len) {
        max_len = len;
        string = realloc(string, len);
        if (!string) {
          free(proc_env);
          return NULL;
        }
        memset(string, 0, len);
      }
      strcpy(string, bin);
      if (library_path && library_path[0]) {
        strcat(string, PATH_SEPARATOR);
        strcat(string, library_path);
      }
      break;
      #else
      case EL_EXE_PATH:
      exe_path = getenv(EXE_PATH);
      len = strlen(bin)  + BUF_LEN;
      if (exe_path) {
        len += strlen(exe_path);
      }
      if (len > max_len) {
        max_len = len;
        string = realloc(string, len);
        if (!string) {
          free(proc_env);
          return NULL;
        }
        memset(string, 0, len);
      }
      strcpy(string, bin);
      if (exe_path && exe_path[0]) {
        strcat(string, PATH_SEPARATOR);
        strcat(string, exe_path);
      }
      break;
      #endif
      
      #ifdef __FreeBSD__
      case EL_DYN_32_LIBRARY_PATH:
      library_path = getenv(DYN_32_LIBRARY_PATH);
      len = strlen(bin)  + BUF_LEN;
      if (library_path) {
        len += strlen(library_path);
      }
      if (len > max_len) {
        max_len = len;
        string = realloc(string, len);
        if (!string) {
          free(proc_env);
          return NULL;
        }
        memset(string, 0, len);
      }
      strcpy(string, bin);
      if (library_path && library_path[0]) {
        strcat(string, PATH_SEPARATOR);
        strcat(string, library_path);
      }
      break;
      #endif

      case EL_BABEL_DATADIR:
      if (od->field.babel_datadir[0]) {
        len = strlen(od->field.babel_datadir) + BUF_LEN;
        if (len > max_len) {
          max_len = len;
          string = realloc(string, len);
          if (!string) {
            free(proc_env);
            return NULL;
          }
          memset(string, 0, len);
        }
        strcpy(string, od->field.babel_datadir);
      }
      else {
        babel_datadir = getenv(BABEL_DATADIR_ENV);
        if (babel_datadir) {
          len = strlen(babel_datadir) + BUF_LEN;
          if (len > max_len) {
            max_len = len;
            string = realloc(string, len);
            if (!string) {
              free(proc_env);
              return NULL;
            }
            memset(string, 0, len);
          }
          strcpy(string, babel_datadir);
        }
        else {
          dont_append = 1;
        }
      }
      break;
      
      case EL_BABEL_LIBDIR:
      if (od->field.babel_libdir[0]) {
        len = strlen(od->field.babel_libdir) + BUF_LEN;
        if (len > max_len) {
          max_len = len;
          string = realloc(string, len);
          if (!string) {
            free(proc_env);
            return NULL;
          }
          memset(string, 0, len);
        }
        strcpy(string, od->field.babel_libdir);
      }
      else {
        babel_libdir = getenv(BABEL_LIBDIR_ENV);
        if (babel_libdir) {
          len = strlen(babel_libdir) + BUF_LEN;
          if (len > max_len) {
            max_len = len;
            string = realloc(string, len);
            if (!string) {
              free(proc_env);
              return NULL;
            }
            memset(string, 0, len);
          }
          strcpy(string, babel_libdir);
        }
        else {
          dont_append = 1;
        }
      }
      break;
    }
    if (!dont_append) {
      #ifndef WIN32
      proc_env[j] = malloc(len);
      if (!(proc_env[j])) {
        free_array(proc_env);
        free(string);
        return NULL;
      }
      sprintf(proc_env[j], "%s=%s", personalized_env[k].name, string);
      ++j;
      #else
      sprintf(&proc_env[j], "%s=%s", personalized_env[k].name, string);
      j += (strlen(&proc_env[j]) + MAX_NAME_LEN);
      #endif
    }
    ++k;
  }
  #ifndef WIN32
  proc_env[j] = NULL;
  #endif
  if (string) {
    free(string);
  }
  
  return proc_env;
}
