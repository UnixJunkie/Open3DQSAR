/*

main.c

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
#include <include/keywords.h>
#include <include/ff_parm.h>
#include <include/proc_env.h>
#include <include/rl_runtime.h>
#include <sys/stat.h>
#ifdef WIN32
#include <include/nice_windows.h>
#define SCREEN_BUFFER_X      1024
#define SCREEN_BUFFER_Y      5000
#ifndef ENABLE_EXTENDED_FLAGS
#define ENABLE_EXTENDED_FLAGS    0x0080
#endif
#ifndef ENABLE_INSERT_MODE
#define ENABLE_INSERT_MODE    0x0020
#endif
#ifndef ENABLE_QUICK_EDIT_MODE
#define ENABLE_QUICK_EDIT_MODE    0x0040
#endif
#endif

#ifdef O3A
#define PACKAGE_CODE      "O3A"
#elif O3G
#define PACKAGE_CODE      "O3G"
#elif O3Q
#define PACKAGE_CODE      "O3Q"
#endif
#ifdef HAVE_EDITLINE_FUNCTIONALITY
#ifdef HAVE_GNU_READLINE
#include <readline/readline.h>
#include <readline/history.h>
#else
#include <editline/readline.h>
#ifndef WIN32
#include <histedit.h>
#endif
#endif
#endif
#define EL_RC_FILE      ".editrc"


#ifndef WIN32
EditLineData el_data;
#endif
O3Data *extern_od;
void *(*_dlsym_add_history)(const char *) = NULL;
HIST_ENTRY *(*_dlsym_next_history)(void) = NULL;
HIST_ENTRY *(*_dlsym_previous_history)(void) = NULL;
int (*_dlsym_read_history)(const char *) = NULL;
int (*_dlsym_write_history)(const char *) = NULL;
void (*_dlsym_using_history)(void) = NULL;
char *(*_dlsym_readline)(const char *) = NULL;
void (*_dlsym_rl_free)(void *) = NULL;
void *(*_dlsym_rl_attempted_completion_function) = NULL;
#ifdef WIN32
void *(*_dlsym_rl_user_completion_entry_free_function) = NULL;
#endif
int *_dlsym_rl_attempted_completion_over = NULL;
int *_dlsym_rl_completion_append_character = NULL;
char **(*_dlsym_rl_completion_matches)(const char *, void *) = NULL;
char *(*_dlsym_rl_filename_completion_function)(const char *, int) = NULL;
char *(*_dlsym_rl_line_buffer) = NULL;
int *_dlsym_rl_point = NULL;
int *_dlsym_rl_catch_signals = NULL;
int *_dlsym_rl_delete_text = NULL;
char have_editline = 0;
char have_gnu_readline = 0;


char *o3_get_keyword(int *keyword_len)
{
  char *keyword = NULL;
  int i = 0;
  
  
  while (rl_line_buffer[i] && isspace(rl_line_buffer[i])) {
    ++i;
  }
  if (rl_line_buffer[i]) {
    *keyword_len = 0;
    keyword = &rl_line_buffer[i];
    while (rl_line_buffer[i] && (!isspace(rl_line_buffer[i]))) {
      ++i;
      ++(*keyword_len);
    }
  }
  
  return keyword;
}


char *o3_completion_generator(const char *text, int state)
{
  char buffer[BUF_LEN];
  char *completion = NULL;
  char *file_completion = NULL;
  char *keyword = NULL;
  char *parameter = NULL;
  char *temp_string;
  int open_quote = 0;
  int skip = 0;
  int i;
  int j;
  int k;
  int m;
  int n;
  int ff_n;
  int probe_n;
  int default_seeds[3];
  int keyword_len;
  int parameter_len;
  int text_len;
  int only_space = 0;
  int complete = 0;
  struct stat filestat;
  
  
  #ifndef WIN32
  rl_completion_append_character = '\0';
  #endif
  /*
  check where we are on the line to decide
  which completion should be accomplished
  */
  i = rl_point;
  text_len = strlen(text);
  while (i) {
    --i;
    /*
    toggle open_quote state if quotes are found
    */
    if (rl_line_buffer[i] == '\"') {
      complete ^= O3_COMPLETION_QUOTE;
    }
    else {
      if (isspace(rl_line_buffer[i])) {
        j = i;
        only_space = 1;
        while (j && only_space) {
          --j;
          only_space = isspace(rl_line_buffer[j]);
        }
        if (!(complete & O3_COMPLETION_QUOTE)) {
          if (!only_space) {
            complete |= O3_COMPLETION_PARAMETER;
          }
          else {
            complete |= O3_COMPLETION_KEYWORD;
          }
        }
        else {
          break;
        }
      }
    }
    if (rl_line_buffer[i] == '=') {
      if (((complete & O3_COMPLETION_PARAMETER)
        && (complete & O3_COMPLETION_QUOTE))
        || (complete & O3_COMPLETION_QUOTE)  //check
        || (!complete)) {
        complete |= O3_COMPLETION_VALUE;
        complete &= (~O3_COMPLETION_PARAMETER);
      }
      break;
    }
  }
  if (!complete) {
    complete = O3_COMPLETION_KEYWORD;
  }
  complete &= (~O3_COMPLETION_QUOTE);
  switch (complete) {
    case O3_COMPLETION_KEYWORD:
    j = 0;
    n = 0;
    while (keyword_data[j].keyword) {
      if (!strncasecmp(keyword_data[j].keyword,
        text, text_len)) {
        if (n == state) {
          if ((!completion)
            && (completion = malloc(strlen(keyword_data[j].keyword) + 2))) {
            strcpy(completion, keyword_data[j].keyword);
          }
        }
        ++n;
      }
      ++j;
    }
    if ((completion) && (n == 1)) {
      strcat(completion, " ");
    }
    break;
    
    case O3_COMPLETION_PARAMETER:
    /*
    get the keyword from rl_line_buffer
    */
    keyword = o3_get_keyword(&keyword_len);
    if (keyword) {
      j = 0;
      while (keyword_data[j].keyword
        && strncasecmp(keyword, keyword_data[j].keyword, keyword_len)) {
        ++j;
      }
      if (keyword_data[j].keyword) {
        k = 0;
        n = 0;
        while (keyword_data[j].parameter_data[k].parameter) {
          if (!strncasecmp
            (keyword_data[j].parameter_data[k].parameter,
            text, text_len)) {
            if (n == state) {
              if ((!completion) && (completion = malloc(strlen
                (keyword_data[j].parameter_data[k].parameter) + 2))) {
                strcpy(completion,
                  keyword_data[j].parameter_data[k].parameter);
              }
            }
            ++n;
          }
          ++k;
        }
      }
    }
    if ((completion) && (n == 1)) {
      strcat(completion, "=");
    }
    break;
    
    case O3_COMPLETION_VALUE:
    keyword = o3_get_keyword(&keyword_len);
    if (keyword) {
      j = 0;
      while (keyword_data[j].keyword
        && strncasecmp(keyword, keyword_data[j].keyword, keyword_len)) {
        ++j;
      }
      if (keyword_data[j].keyword) {
        i = rl_point;
        while (i && (rl_line_buffer[i - 1] != '=')) {
          --i;
        }
        k = i;
        while (i && (!isspace(rl_line_buffer[i - 1]))) {
          --i;
        }
        parameter = &rl_line_buffer[i];
        parameter_len = k - i - 1;
        if (!parameter_len) {
          break;
        }
        k = 0;
        while (keyword_data[j].parameter_data[k].parameter
          && strncasecmp(parameter, keyword_data[j].parameter_data
          [k].parameter, parameter_len)) {
          ++k;
        }
        if (keyword_data[j].parameter_data[k].parameter) {
          switch (keyword_data[j].parameter_data[k].type) {
            case O3_PARAM_NUMERIC:
            case O3_PARAM_STRING:
            m = 0;
            n = 0;
            while (keyword_data[j].parameter_data[k].choice[m]) {
              if (!strncasecmp
                (keyword_data[j].parameter_data[k].choice[m],
                text, text_len)) {
                if (n == state) {
                  if ((!completion) && (completion = malloc(strlen
                    (keyword_data[j].parameter_data[k].choice[m]) + 2))) {
                    strcpy(completion,
                      keyword_data[j].parameter_data[k].choice[m]);
                  }
                }
                ++n;
              }
              ++m;
            }
            if ((completion) && (n == 1)) {
              strcat(completion, " ");
            }
            break;

            case O3_PARAM_FILE:
            case O3_PARAM_DIRECTORY:
            #ifndef WIN32
            n = rl_point - text_len;
            m = 0;
            while (n && (rl_line_buffer[n - 1] != '\"')
              && (rl_line_buffer[n - 1] != '=')) {
              --n;
              ++m;
            }
            open_quote = ((rl_line_buffer[n - 1] == '\"') ? 1 : 0);
            if ((rl_point >= n) && (temp_string = malloc(rl_point - n + 1))) {
              if (rl_point > n) {
                strncpy(temp_string, &rl_line_buffer[n], rl_point - n);
              }
              temp_string[rl_point - n] = '\0';
              file_completion = rl_filename_completion_function(temp_string, state);
              completion = (file_completion ? strdup(file_completion) : NULL);
              free(temp_string);
            }
            if (completion) {
              #ifdef HAVE_EDITLINE_FUNCTIONALITY
              #ifdef HAVE_GNU_READLINE
              rl_free(file_completion);
              #else
              free(file_completion);
              #endif
              #else
              if (have_editline && _dlsym_rl_free) {
                rl_free(file_completion);
              }
              else {
                free(file_completion);
              }
              #endif
              n = strlen(completion);
              if ((completion = realloc(completion, n + 4))) {
                skip = 0;
                if (strchr(completion, ' ') && (!open_quote)) {
                  memmove(&completion[1], completion, n + 1);
                  completion[0] = '\"';
                  skip = 1;
                }
                if (stat(&completion[skip], &filestat) != -1) {
                  if (S_ISDIR(filestat.st_mode)) {
                    strcat(completion, "/");
                  }
                  else {
                    if (open_quote) {
                      strcat(completion, "\"");
                    }
                    strcat(completion, " ");
                  }
                }
              }
            }
            if (m && completion) {
              memmove(completion, &completion[m], strlen(completion) + 1);
            }
            #else
            file_completion = rl_filename_completion_function(text, state);
            completion = (file_completion ? strdup(file_completion) : NULL);
            if (completion) {
              #ifdef HAVE_EDITLINE_FUNCTIONALITY
              #ifdef HAVE_GNU_READLINE
              rl_free(file_completion);
              #else
              free(file_completion);
              #endif
              #else
              if (have_editline && _dlsym_rl_free) {
                rl_free(file_completion);
              }
              else {
                free(file_completion);
              }
              #endif
              if ((temp_string = strdup(completion))) {
                n = strlen(temp_string);
                if (n && (open_quote = (temp_string[n - 1] == '\"'))) {
                  temp_string[n - 1] = '\0';
                }
                skip = ((temp_string[0] == '\"') ? 1 : 0);
                if (stat(&temp_string[skip], &filestat) != -1) {
                  if (S_ISDIR(filestat.st_mode)) {
                    if (open_quote) {
                      completion[n - 1] ='\0';
                    }
                  }
                }
                free(temp_string);
              }
            }
            #endif
            break;
            
            case O3_PARAM_PC:
            if (extern_od->pc_num) {
              sprintf(buffer, "%d ", extern_od->pc_num);
              if ((!strncasecmp(buffer, text, text_len)) && (!state)) {
                completion = strdup(buffer);
              }
            }
            break;
            
            case O3_PARAM_DESIGN_POINTS:
            if (extern_od->mal.x_weights && extern_od->mal.x_weights->m) {
              sprintf(buffer, "%d ", extern_od->mal.x_weights->m);
              if ((!strncasecmp(buffer, text, text_len)) && (!state)) {
                completion = strdup(buffer);
              }
            }
            break;
            
            case O3_PARAM_N_CPUS:
            n = (extern_od->n_proc ? extern_od->n_proc : 1);
            sprintf(buffer, "%d ", n);
            if ((!strncasecmp(buffer, text, text_len)) && (!state)) {
              completion = strdup(buffer);
            }
            break;
            
            case O3_PARAM_PROBE:
            n = 0;
            ff_n = 0;
            probe_n = 0;
            while (ff_n < MAX_FF_N) {
              while (ff_parm[ff_n][probe_n].type_num) {
                if (!strncmp(ff_parm[ff_n][probe_n].type_chr, text, text_len)) {
                  if (n == state) {
                    if ((!completion) && (completion = malloc(strlen
                      (ff_parm[ff_n][probe_n].type_chr) + 2))) {
                      strcpy(completion,
                        ff_parm[ff_n][probe_n].type_chr);
                    }
                  }
                  ++n;
                }
                ++probe_n;
              }
              probe_n = 0;
              ++ff_n;
            }
            if (completion && (n == 1)) {
              strcat(completion, " ");
            }
            break;
            
            case O3_PARAM_SCRAMBLE_MAX_BINS:
            get_attr_struct_ave(extern_od, 0, ACTIVE_BIT, &n, NULL);
            if (n) {
              sprintf(buffer, "%d ", n / 3);
              if ((!strncasecmp(buffer, text, text_len)) && (!state)) {
                completion = strdup(buffer);
              }
            }
            break;
            
            case O3_PARAM_SRD_SEEDS:
            default_seeds[0] = extern_od->x_vars / 10;
            default_seeds[1] = extern_od->overall_active_x_vars / 2;
            default_seeds[2] = 3000;
            qsort(default_seeds, 3, sizeof(int), compare_integers);
            n = default_seeds[0];
            if (n) {
              sprintf(buffer, "%d ", n);
              if ((!strncasecmp(buffer, text, text_len)) && (!state)) {
                completion = strdup(buffer);
              }
            }
            break;
          }
        }
      }
    }
    break;
  }
  rl_attempted_completion_over = 1;
  
  return completion;
}


char **o3_completion_matches(const char *text, int start, int end)
{
  return rl_completion_matches(text, o3_completion_generator);
}


#ifdef WIN32
void o3_compentry_free(void *mem)
{
  free(mem);
}
#endif


#ifndef WIN32
void program_signal_handler(int signum)
{
  char current_time[BUF_LEN];


  if ((have_editline
    && (((signum == SIGINT) && (extern_od->prompt)
    && rl_line_buffer && (!rl_line_buffer[0]))
    || ((signum == SIGINT) && (!(extern_od->prompt)))
    || (!(rl_line_buffer)) || (signum != SIGINT)))
    || (!have_editline)) {
    SET_INK(extern_od, NORMAL_INK);
    tee_printf(extern_od, "\nJob terminated ");
    memset(current_time, 0, BUF_LEN);
    if (!get_current_time(current_time)) {
      remove_newline(current_time);
      tee_printf(extern_od, "on %s ", current_time);
    }
    tee_printf(extern_od, "by signal SIG");
    switch (signum) {
      case SIGINT:
      tee_printf(extern_od, "INT.\n");
      break;
      
      case SIGTERM:
      tee_printf(extern_od, "TERM.\n");
      break;
      
      case SIGHUP:
      tee_printf(extern_od, "HUP.\n");
      break;
      
      default:
      tee_printf(extern_od, "UNK.\n");
      break;
    }
    tee_flush(extern_od);
    reset_user_terminal(extern_od);
    remove_temp_files(PACKAGE_CODE);
    signal(signum, SIG_DFL);
    raise(signum);
  }
  else if (have_editline && rl_line_buffer) {
    pthread_cancel(el_data.thread_id);
  }
}
#else
BOOL program_signal_handler(DWORD fdwCtrlType)
{
  char current_time[BUF_LEN];


  if (have_editline && (fdwCtrlType == CTRL_C_EVENT)
    && extern_od->prompt && rl_line_buffer
    && strlen(rl_line_buffer)) {
    return TRUE;
  }
  if (extern_od->prompt) {
    SET_INK(extern_od, NORMAL_INK);
  }
  tee_printf(extern_od, "\nJob terminated ");
  memset(current_time, 0, BUF_LEN);
  if (!get_current_time(current_time)) {
    remove_newline(current_time);
    tee_printf(extern_od, "on %s ", current_time);
  }
  tee_printf(extern_od, "by signal CTRL_");
  switch (fdwCtrlType) {
    case CTRL_C_EVENT:
    tee_printf(extern_od, "C_EVENT.\n");
    break;

    case CTRL_CLOSE_EVENT:
    tee_printf(extern_od, "CLOSE_EVENT.\n");
    break;

    case CTRL_BREAK_EVENT:
    tee_printf(extern_od, "BREAK_EVENT.\n");
    break;

    case CTRL_LOGOFF_EVENT:
    tee_printf(extern_od, "LOGOFF_EVENT.\n");
    break;

    case CTRL_SHUTDOWN_EVENT:
    tee_printf(extern_od, "SHUTDOWN_EVENT.\n");
    break;
  }
  tee_flush(extern_od);
  remove_temp_files(PACKAGE_CODE);
  
  return FALSE;
}
#endif


void reset_user_terminal(O3Data *od)
{
  #ifndef WIN32
  if (have_editline && od->user_termios) {
    tcsetattr(STDOUT_FILENO, TCSADRAIN, od->user_termios);
  }
  SET_INK(od, DEFAULT_INK);
  #endif
}


int main(int argc, char **argv)
{
  char *temp_dir_string;
  char *save_ram_string;
  char *n_cpus_string;
  char *nice_string;
  char *babel_path_string;
  char *pymol_string;
  char *jmol_string;
  char current_time[BUF_LEN];
  char current_dir[BUF_LEN];
  char buffer[BUF_LEN];
  char bin[BUF_LEN];
  char temp_dir_env[BUF_LEN];
  char *appdata;
  extern char E_PROGRAM_EXIT[];
  extern char E_OUT_OF_MEMORY[];
  extern char E_FILE_CANNOT_BE_OPENED_FOR_WRITING[];
  extern char E_FILE_CANNOT_BE_OPENED_FOR_READING[];
  extern char E_OPENBABEL_NOT_WORKING[];
  extern char E_OPENBABEL_DATA_PLUGINS[];
  extern char E_OPENBABEL_MISSING_OR_TOO_OLD[];
  extern char M_NUMBER_OF_CPUS[];
  extern char O3_FAILED[];
  char *too_many_arguments =
    PACKAGE_NAME_LOWERCASE": Too many arguments\n";
  char *invalid_option =
    PACKAGE_NAME_LOWERCASE": invalid option -- %s\n";
  char *unrecognized_option =
    PACKAGE_NAME_LOWERCASE": unrecognized option `%s'\n";
  char *option_requires_argument =
    PACKAGE_NAME_LOWERCASE": option requires an argument -- %s\n";
  char *try_help_usage =
    "Try `"PACKAGE_NAME_LOWERCASE" --help' or `"
    PACKAGE_NAME_LOWERCASE" --usage' "
    "for more information.\n";
  char *package_version =
    "\n"
    PACKAGE_NAME" version " VERSION "\n"
    #ifdef O3Q
    "Copyright (C) 2009-2015 Paolo Tosco, Thomas Balle\n"
    #else
    "Copyright (C) 2010-2014 Paolo Tosco, Thomas Balle\n"
    #endif
    "Licensed under the terms of GPLv3\n"
    "Report bugs to " PACKAGE_BUGREPORT "\n"
    "\n";
  char help[] =
    "\n"
    "Usage: "PACKAGE_NAME_LOWERCASE" [OPTION...]\n"
    "\n"
    "------------\n"
    PACKAGE_NAME"\n"
    "------------\n"
    "\n"
    #ifdef O3Q
    "An open-source software aimed at high-throughput\n"
    "chemometric analysis of molecular interaction fields\n"
    #elif O3G
    "An open-source software aimed at high-throughput\n"
    "generation of molecular interaction fields (MIFs)\n"
    #elif O3A
    "An open-source software aimed at unsupervised molecular alignment\n"
    #endif
    "\n"
    "Version " VERSION "\n"
    #ifdef O3Q
    "Copyright (C) 2009-2015 Paolo Tosco, Thomas Balle\n"
    #else
    "Copyright (C) 2010-2014 Paolo Tosco, Thomas Balle\n"
    #endif
    "\n"
    "\n"
    "  -i <filein>                Input is read from <filein>\n"
    "  -o <fileout>               Output is written to <fileout>\n"
    "  -p                         Input is piped through standard input\n"
    "  -?, --help                 Give this help list\n"
    "      --usage                Give a short usage message\n"
    "  -V, --version              Print program version\n"
    "\n"
    "Report bugs to " PACKAGE_BUGREPORT "\n"
    "\n";
  char usage[] =
    "Usage: "PACKAGE_NAME_LOWERCASE" [-p?V] [-i FILE] [-o FILE] "
    "[--help] [--usage] [--version]\n";
  int i;
  int result;
  int found;
  int nice_value;
  O3Data od;
  CLIArgs cli_args;
  #ifndef WIN32
  struct sigaction setup_action;
  sigset_t block_mask;
  struct termios user_termios;
  void *dl_handle;
  #else
  char cmd_cli[BUF_LEN];
  int wincmd = 0;
  int error = 0;
  DWORD n;
  HANDLE hIcon = NULL;
  HMODULE hModule = NULL;
  COORD dwSize;
  COORD cursor_home;
  DWORD n_chars;
  WORD wVersionRequested;
  WSADATA wsaData;
  HMODULE dl_handle;
  #endif
  char el_rc[BUF_LEN];
  FILE *el_rc_handle = NULL;

  #undef HAVE_DEFAULT_EL_RC
  #ifdef WIN32
  #define HAVE_DEFAULT_EL_RC  1
  char default_el_rc[] =
    "history size 200\n";
  char el_rc_env[BUF_LEN];
  #endif
  #ifndef __ICC
  #ifdef __FreeBSD__
  #define HAVE_DEFAULT_EL_RC  1
  char default_el_rc[] =
    "history size 200\n"
    "bind ^W ed-delete-prev-word\n"
    "bind ^[[3~ ed-delete-next-char\n"
    "bind ^[[1;5D vi-prev-word\n"
    "bind ^[[1;5C vi-next-word\n";
  #elif __APPLE__
  #define HAVE_DEFAULT_EL_RC  1
  char default_el_rc[] =
    "history size 200\n"
    "bind ^W ed-delete-prev-word\n"
    "bind ^[[3~ ed-delete-next-char\n"
    "bind ^[[1;5D vi-prev-word\n"
    "bind ^[[1;5C vi-next-word\n"
    "bind ^I rl_complete\n";
  #elif sun
  #define HAVE_DEFAULT_EL_RC  1
  char default_el_rc[] =
    "history size 200\n"
    "bind ^W ed-delete-prev-word\n"
    "bind ^[[3~ ed-delete-next-char\n"
    "bind ^[[5D vi-prev-word\n"
    "bind ^[[5C vi-next-word\n";
  #endif
  #endif
  #ifdef linux
  #define HAVE_DEFAULT_EL_RC  1
  char default_el_rc[] =
    "history size 200\n"
    "bind ^W ed-delete-prev-word\n"
    "bind ^[[3~ ed-delete-next-char\n"
    "bind ^[[1;5D vi-prev-word\n"
    "bind ^[[1;5C vi-next-word\n";
  #endif
  #ifndef O3G
  char *random_seed_string;
  long random_seed;
  #endif
  #ifdef O3Q
  char *gnuplot_string;
  #endif
  #ifdef O3A
  char *pharao_string;
  char *tinker_path_string;
  extern char E_PHARAO_NOT_WORKING[];
  extern char E_PHARAO_MISSING_OR_TOO_OLD[];
  #else
  char *qm_engine_string;
  char *cs3d_string;
  char *md_grid_path_string;
  #endif


  memset(&cli_args, 0, sizeof(CLIArgs));
  memset(&od, 0, sizeof(O3Data));
  extern_od = &od;
  memset(current_dir, 0, BUF_LEN);
  memset(buffer, 0, BUF_LEN);
  memset(temp_dir_env, 0, BUF_LEN);
  memset(bin, 0, BUF_LEN);
  strcpy(od.package_code, PACKAGE_CODE);
  strncpy(bin, argv[0], BUF_LEN - 2);
  absolute_path(bin);
  get_dirname(bin);
  dl_handle = check_readline();
  if (have_editline) {
    memset(el_rc, 0, BUF_LEN);
    #ifndef HAVE_EDITLINE_FUNCTIONALITY
    rl_attempted_completion_function = (void *)o3_completion_matches;
    #else
    rl_attempted_completion_function = o3_completion_matches;
    #endif
    #ifdef WIN32
    rl_user_completion_entry_free_function = o3_compentry_free;
    #else
    #if (defined HAVE_EDITLINE_FUNCTIONALITY && defined HAVE_GNU_READLINE)
    rl_catch_signals = 0;
    #endif
    #ifndef HAVE_EDITLINE_FUNCTIONALITY
    if (have_gnu_readline) {
      rl_catch_signals = 0;
    }
    #endif
    memset(&el_data, 0, sizeof(EditLineData));
    #endif
  }
  #ifndef WIN32
  memset (&user_termios, 0, sizeof(struct termios));
  od.user_termios = &user_termios;
  #endif
  cli_args.prompt = INTERACTIVE_RUN;
  i = 1;
  while (i < argc) {
    if (!strcmp(argv[i], "-p")) {
      cli_args.prompt = 0;
      ++i;
    }
    else if (!strcmp(argv[i], "--debug")) {
      od.debug = 1;
      ++i;
    }
    else if ((!strcasecmp(argv[i], "-v"))
      || (!strcmp(argv[i], "--version"))) {
      puts(package_version);
      return 0;
    }
    else if ((!strcmp(argv[i], "-?"))
      || (!strcmp(argv[i], "--help"))) {
      puts(help);
      return 0;
    }
    else if (!strcmp(argv[i], "--usage")) {
      puts(usage);
      return 0;
    }
    else if (!strcmp(argv[i], "-i")) {
      if ((i + 1) < argc) {
        strncpy(cli_args.input_file, argv[i + 1], BUF_LEN);
        cli_args.input_file[BUF_LEN - 2] = '\0';
        cli_args.input = 1;
        cli_args.prompt = 0;
        i += 2;
        continue;
      }
      else {
        fprintf(stderr, option_requires_argument, "i");
        fputs(try_help_usage, stderr);
        return -1;
      }
    }
    else if (!strcmp(argv[i], "-o")) {
      if ((i + 1) < argc) {
        strncpy(cli_args.output_file, argv[i + 1], BUF_LEN);
        cli_args.output_file[BUF_LEN - 2] = '\0';
        cli_args.output = 1;
        i += 2;
        continue;
      }
      else {
        fprintf(stderr, option_requires_argument, "o");
        fputs(try_help_usage, stderr);
        return -1;
      }
    }
    else if (!strcmp(argv[i], "--term")) {
      od.terminal = 1;
      #ifdef WIN32
      wincmd = 1;
      #endif
      ++i;
    }
    #ifdef WIN32
    else if (!strcmp(argv[i], "--wincmd")) {
      wincmd = 1;
      ++i;
    }
    #endif
    else {
      if (!strncmp(argv[i], "--", 2)) {
        fprintf(stderr, unrecognized_option, argv[i]);
      }
      else if (argv[i][0] == '-') {
        fprintf(stderr, invalid_option, argv[i]);
      }
      else {
        fputs(too_many_arguments, stderr);
      }
      fputs(try_help_usage, stderr);
      return -1;
    }
  }
  if (!getcwd(current_dir, BUF_LEN)) {
    current_dir[0] = '\0';
  }
  current_dir[BUF_LEN - 2] = '\0';
  #ifdef WIN32
  if (cli_args.prompt && (!wincmd)) {
    /*
    if we are running under Windows
    and we have been started in interactive
    mode, restart ourselves using cmd.exe
    so we have wheel scrolling
    in the console window
    */
    strcpy(cmd_cli, " /c");
    for (i = 0; i < argc; ++i) {
      strcat(cmd_cli, " ");
      if (strchr(argv[i], ' ')) {
        strcat(cmd_cli, "\"");
      }
      strcat(cmd_cli, argv[i]);
      if (strchr(argv[i], ' ')) {
        strcat(cmd_cli, "\"");
      }
    }
    strcat(cmd_cli, " --wincmd");
    result = (INT_PTR)ShellExecute(NULL, TEXT("open"),
      TEXT(getenv("COMSPEC")),
      TEXT(cmd_cli), TEXT(current_dir), SW_SHOW);
    if (result <= 32) {
      result = (INT_PTR)ShellExecute(NULL, TEXT("open"),
        TEXT("cmd.exe"),
        TEXT(cmd_cli), TEXT(current_dir), SW_SHOW);
    }
    if (result > 32) {
      return 0;
    }
  }
  #else
  if (have_editline && tcgetattr(STDOUT_FILENO, &user_termios)) {
    od.user_termios = NULL;
  }
  #endif
  /*
  otherwise Intel MKL on IVE-PLS will be a memory hog
  */
  #ifdef HAVE_LIBMKL
  putenv("MKL_DISABLE_FAST_MM=1");
  #endif
  od.active_field_num = -1;
  od.mmap_field_num = -1;
  od.active_object_num = -1;
  /*
  trap SIGINT, SIGTERM and SIGHUP
  signals so temporary files can be
  removed if the user terminates the program
  */
  #ifndef WIN32
  sigemptyset(&block_mask);
  sigaddset(&block_mask, SIGTERM);
  sigaddset(&block_mask, SIGHUP);
  setup_action.sa_handler = program_signal_handler;
  setup_action.sa_mask = block_mask;
  setup_action.sa_flags = 0;
  sigaction(SIGINT, &setup_action, NULL);
  sigemptyset(&block_mask);
  sigaddset(&block_mask, SIGINT);
  sigaddset(&block_mask, SIGHUP);
  setup_action.sa_handler = program_signal_handler;
  setup_action.sa_mask = block_mask;
  setup_action.sa_flags = 0;
  sigaction(SIGTERM, &setup_action, NULL);
  sigemptyset(&block_mask);
  sigaddset(&block_mask, SIGTERM);
  sigaddset(&block_mask, SIGINT);
  setup_action.sa_handler = program_signal_handler;
  setup_action.sa_mask = block_mask;
  setup_action.sa_flags = 0;
  sigaction(SIGHUP, &setup_action, NULL);
  #else
  SetConsoleCtrlHandler((PHANDLER_ROUTINE)
    program_signal_handler, TRUE);
  #endif
  od.prompt = cli_args.prompt;
  od.in = stdin;
  #ifdef WIN32
  od.hInput = GetStdHandle(STD_INPUT_HANDLE);
  od.hOutput = GetStdHandle(STD_OUTPUT_HANDLE);
  #endif
  if (alloc_file_descriptor(&od, MAX_FILES)) {
    fprintf(stderr, E_OUT_OF_MEMORY, O3_FAILED);
    if (od.prompt) {
      fprintf(stderr, E_PROGRAM_EXIT);
      fflush(stderr);
      #ifndef WIN32
      fgets(buffer, BUF_LEN - 2, od.in);
      #else
      ReadConsole(od.hInput, buffer, BUF_LEN - 2, &n, NULL);
      #endif
    }
    return -1;
  }
  if (cli_args.input) {
    od.file[MAIN_INPUT]->handle =
      fopen(cli_args.input_file, "rb");
    if (!(od.file[MAIN_INPUT]->handle)) {
      fprintf(stderr, E_FILE_CANNOT_BE_OPENED_FOR_READING,
        cli_args.input_file, O3_FAILED);
      return -1;
    }
    strcpy(od.file[MAIN_INPUT]->name, cli_args.input_file);
    od.in = od.file[MAIN_INPUT]->handle;
    od.prompt = 0;
  }
  if (cli_args.output) {
    od.file[MAIN_OUTPUT]->handle =
      fopen(cli_args.output_file, "wb");
    if (!(od.file[MAIN_OUTPUT]->handle)) {
      fprintf(stderr, E_FILE_CANNOT_BE_OPENED_FOR_WRITING,
        cli_args.output_file, O3_FAILED);
      return -1;
    }
    strcpy(od.file[MAIN_OUTPUT]->name, cli_args.output_file);
    od.out = od.file[MAIN_OUTPUT]->handle;
  }
  #ifdef WIN32
  if (od.prompt) {
    memset(&cursor_home, 0, sizeof(COORD));
    memset(&dwSize, 0, sizeof(COORD));
    SetConsoleTitle(PACKAGE_NAME);
    /*
    paint the console background in white
    foreground in blue
    */
    SetConsoleMode(od.hInput,
      ENABLE_EXTENDED_FLAGS | ENABLE_PROCESSED_INPUT
      | ENABLE_LINE_INPUT | ENABLE_ECHO_INPUT
      | ENABLE_QUICK_EDIT_MODE | ENABLE_INSERT_MODE);
    SetConsoleMode(od.hOutput, ENABLE_PROCESSED_OUTPUT
      | ENABLE_WRAP_AT_EOL_OUTPUT);
    if (!(od.terminal)) {
      dwSize.X = SCREEN_BUFFER_X;
      dwSize.Y = SCREEN_BUFFER_Y;
      SetConsoleScreenBufferSize(od.hOutput, dwSize);
      SET_INK(&od, NORMAL_INK);
      FillConsoleOutputAttribute(od.hOutput, NORMAL_INK,
        dwSize.X * dwSize.Y, cursor_home, &n_chars);
    }
    /*
    set the program icon
    on the console window
    */
    hModule = GetModuleHandle(argv[0]);
    if (hModule && (!(od.terminal))) {
      hIcon = LoadImage(hModule, MAKEINTRESOURCE(3),
        IMAGE_ICON, 0, 0, LR_DEFAULTSIZE);
      if (hIcon) {
        SendMessage(GetConsoleWindow(), WM_SETICON,
          ICON_SMALL, (LPARAM)hIcon);
      }
    }  
  }
  #else
  SET_INK(&od, NORMAL_INK);
  #endif
  od.mti = MERSENNE_N + 1;
  tee_printf(&od, "%s\n", package_version);
  tee_flush(&od);
  get_system_information(&od);
  n_cpus_string = getenv("O3_N_CPUS");
  determine_best_cpu_number(&od, "all");
  if (n_cpus_string) {
    determine_best_cpu_number(&od, n_cpus_string);
  }
  tee_printf(&od, M_NUMBER_OF_CPUS, PACKAGE_NAME, od.n_proc);
  
  /*
  if the O3_NICE environment variable
  has not been set, then set default nice at PRIO_MIN
  */
  memset(od.home_dir, 0, BUF_LEN);
  appdata = getenv(APPDATA);
  if (appdata) {
    sprintf(od.home_dir, "%s%c", appdata, SEPARATOR);
    #ifndef WIN32
    if (have_editline && (!have_gnu_readline)) {
      sprintf(el_rc, "%s%c%s", appdata, SEPARATOR, EL_RC_FILE);
    }
    #endif
  }
  strcat(od.home_dir, HOMEDIR);
  #ifndef WIN32
  nice_value = PRIO_MAX;
  mkdir(od.home_dir, S_IRWXU | S_IRGRP | S_IROTH);
  #else
  nice_value = 1;
  mkdir(od.home_dir);
  if (have_editline && (!have_gnu_readline)) {
    sprintf(el_rc, "%s%c%s", od.home_dir, SEPARATOR, EL_RC_FILE);
    sprintf(el_rc_env, "EDITRC=%s", el_rc);
    putenv(el_rc_env);
  }
  #endif
  if (have_editline && (!have_gnu_readline)) {
    #ifdef HAVE_DEFAULT_EL_RC
    if (el_rc[0]) {
      if (access(el_rc, F_OK) == -1) {
        el_rc_handle = fopen(el_rc, "w+");
        if (el_rc_handle) {
          fprintf(el_rc_handle, "%s", default_el_rc);
          fclose(el_rc_handle);
        }
      }
    }
    #endif
  }
  nice_string = getenv("O3_NICE");
  if (nice_string) {
    #ifndef WIN32
    sscanf(nice_string, "%d", &nice_value);
    #else
    nice_value = 0;
    while (nice_value < 6) {
      if (!strcasecmp(nice_string, nice_name[nice_value])) {
        break;
      }
      ++nice_value;
    }
    #endif
  }
  #ifndef WIN32
  if ((nice_value < PRIO_MIN) || (nice_value > PRIO_MAX)) {
    nice_value = PRIO_MAX;
  }
  #else
  if (nice_value == 6) {
    nice_value = 1;
  }
  #endif
  set_nice_value(&od, nice_value);
  #ifndef O3G
  random_seed = DEFAULT_RANDOM_SEED;
  if ((random_seed_string = getenv("O3_RANDOM_SEED"))) {
    sscanf(random_seed_string, "%ld", &random_seed);
  }
  set_random_seed(&od, (unsigned long)absval(random_seed));
  #endif
  strcpy(od.temp_dir, current_dir);
  temp_dir_string = GET_TEMPDIR;
  if (temp_dir_string) {
    strcpy(od.temp_dir, temp_dir_string);
  }
  temp_dir_string = getenv(TEMP_DIR_ENV);
  if (temp_dir_string) {
    strcpy(od.temp_dir, temp_dir_string);
    absolute_path(od.temp_dir);
  }
  if (!(od.temp_dir[0])) {
    fprintf(stderr, "Could not find a directory "
      "for temporary storage. Please "
      "consider setting the "
      TEMP_DIR_ENV" environment "
      "variable to a suitable value.\n%s\n\n", O3_FAILED);
    if (od.prompt) {
      fprintf(stderr, E_PROGRAM_EXIT);
      fflush(stderr);
      #ifndef WIN32
      fgets(buffer, BUF_LEN - 2, od.in);
      #else
      ReadConsole(od.hInput, buffer, BUF_LEN - 2, &n, NULL);
      #endif
      reset_user_terminal(&od);
    }
    return -1;
  }
  sprintf(temp_dir_env, TEMP_DIR_ENV"=%s", od.temp_dir);
  putenv(temp_dir_env);
  tee_printf(&od,
    "The directory for temporary storage is:\n"
    "%s\n\n"
    "The current working directory is:\n"
    "%s\n\n", od.temp_dir, current_dir);
  #if (!defined HAVE_LIBMINIZIP) || (!defined HAVE_MINIZIP_ZIP_H) || (!defined HAVE_MINIZIP_UNZIP_H)
  tee_printf(&od,
      "Since "PACKAGE_NAME" was not linked against libminizip, "
      "support for ZIP files will not be available.\n\n");
  #endif
  tee_flush(&od);
  strcpy(od.field.babel_exe_path, bin);
  if ((babel_path_string = getenv(BABEL_PATH_ENV))) {
    if (babel_path_string[0]) {
      strncpy(od.field.babel_exe_path, babel_path_string, BUF_LEN - 2);
      absolute_path(od.field.babel_exe_path);
    }
  }
  found = 0;
  if (dexist(od.field.babel_exe_path)) {
    sprintf(buffer, "%s%c%s", od.field.babel_exe_path, SEPARATOR, BABEL_EXE);
    if (fexist(buffer)) {
      sprintf(buffer, "%s%c%s", od.field.babel_exe_path, SEPARATOR, OBENERGY_EXE);
      if (fexist(buffer)) {
        found = 1;
      }
    }
  }
  if (!found) {
    if (is_in_path(BABEL_EXE, buffer) && is_in_path(OBENERGY_EXE, buffer)) {
      found = 1;
      strncpy(od.field.babel_exe_path, buffer, BUF_LEN - 2);
      od.field.babel_exe_path[BUF_LEN - 1] = '\0';
      get_dirname(od.field.babel_exe_path);
    }
  }
  result = 0;
  if (found) {
    result = check_babel(&od, bin);
    switch (result) {
      case OUT_OF_MEMORY:
      fprintf(stderr, E_OUT_OF_MEMORY, O3_FAILED);
      if (od.prompt) {
        fprintf(stderr, E_PROGRAM_EXIT);
        fflush(stderr);
        #ifndef WIN32
        fgets(buffer, BUF_LEN - 2, od.in);
        #else
        ReadConsole(od.hInput, buffer, BUF_LEN - 2, &n, NULL);
        #endif
        reset_user_terminal(&od);
      }
      return -1;
    
      case CANNOT_READ_TEMP_FILE:
      fprintf(stderr, E_FILE_CANNOT_BE_OPENED_FOR_READING,
        od.file[TEMP_LOG]->name, "");
      break;
      
      case CANNOT_WRITE_TEMP_FILE:
      fprintf(stderr, E_FILE_CANNOT_BE_OPENED_FOR_WRITING,
        od.file[TEMP_LOG]->name, "");
      break;

      case BABEL_NOT_WORKING:
      fprintf(stderr, E_OPENBABEL_NOT_WORKING);
      if ((od.file[TEMP_LOG]->handle = fopen(od.file[TEMP_LOG]->name, "rb"))) {
        while (fgets(buffer, BUF_LEN, od.file[TEMP_LOG]->handle)) {
          buffer[BUF_LEN - 1] = '\0';
          fprintf(stderr, "%s", buffer);
        }
        fclose(od.file[TEMP_LOG]->handle);
        od.file[TEMP_LOG]->handle = NULL;
      }
      break;

      case BABEL_PLUGINS_NOT_FOUND:
      fprintf(stderr, E_OPENBABEL_DATA_PLUGINS, "");
      break;
    }
  }
  if (found && (!result)) {
    tee_printf(&od, "OpenBabel binaries are in:\n"
      "%s\n", od.field.babel_exe_path);
  }
  else {
    if (found && result) {
      tee_printf(&od, "Working ");
    }
    tee_printf(&od, E_OPENBABEL_MISSING_OR_TOO_OLD, "");
    memset(od.field.babel_exe_path, 0, BUF_LEN);
  }
  tee_printf(&od, "The path to OpenBabel binaries can be %s "
    "through the "BABEL_PATH_ENV" environment variable or "
    "by the \"env babel_path\" keyword.\n\n",
    (found ? "changed" : "set"));
  #ifndef O3A
  memset(od.field.qm_exe, 0, BUF_LEN);
  memset(od.field.qm_exe_path, 0, BUF_LEN);
  found = 1;
  if (!is_in_path(FIREFLY_EXE, buffer)) {
    if (!is_in_path(G09_EXE, buffer)) {
      if (!is_in_path(G03_EXE, buffer)) {
        if (!is_in_path(TURBOMOLE_DSCF_EXE, buffer)) {
          found = 0;
        }
      }
    }
  }
  if (found) {
    strncpy(od.field.qm_exe_path, buffer, BUF_LEN - 2);
  }
  if ((qm_engine_string = getenv("O3_QM_ENGINE"))) {
    if (qm_engine_string[0]) {
      strncpy(od.field.qm_exe_path, qm_engine_string, BUF_LEN - 2);
      absolute_path(od.field.qm_exe_path);
    }
  }
  if ((found = fexist(od.field.qm_exe_path))) {
    tee_printf(&od, "The QM engine is:\n"
      "%s\n", od.field.qm_exe_path);
    strcpy(od.field.qm_exe, get_basename(od.field.qm_exe_path));
    get_dirname(od.field.qm_exe_path);
  }
  else {
    tee_printf(&od, "No QM engine was found.\n");
  }
  tee_printf(&od, "The QM engine can be %s through the "
    "O3_QM_ENGINE environment variable or "
    "by the \"env qm_engine\" keyword.\n\n",
    (found ? "changed" : "chosen"));
  if ((cs3d_string = getenv("O3_CS3D"))) {
    if (cs3d_string[0]) {
      strncpy(od.field.cs3d_exe, cs3d_string, BUF_LEN - 2);
      absolute_path(od.field.cs3d_exe);
    }
  }
  else if (is_in_path(CS3D_EXE, buffer)) {
    strcpy(od.field.cs3d_exe, buffer);
  }
  if ((found = fexist(od.field.cs3d_exe))) {
    tee_printf(&od, "The CS3D executable is:\n"
      "%s\n", od.field.cs3d_exe);
  }
  else {
    tee_printf(&od, "CS3D was not found.\n");
  }
  tee_printf(&od, "The path to CS3D can be %s through the "
    "O3_CS3D environment variable or "
    "by the \"env cs3d\" keyword.\n\n",
    (found ? "changed" : "set"));
  memset(od.field.md_grid_exe_path, 0, BUF_LEN);
  if (is_in_path(GRID_EXE, buffer)) {
    strncpy(od.field.md_grid_exe_path, buffer, BUF_LEN - 2);
    get_dirname(od.field.md_grid_exe_path);
  }
  if ((md_grid_path_string = getenv("O3_MD_GRID_PATH"))) {
    if (md_grid_path_string[0]) {
      strncpy(od.field.md_grid_exe_path, md_grid_path_string, BUF_LEN - 2);
      absolute_path(od.field.md_grid_exe_path);
    }
  }
  found = 0;
  if (dexist(od.field.md_grid_exe_path)) {
    sprintf(buffer, "%s%c%s", od.field.md_grid_exe_path, SEPARATOR, GRID_EXE);
    if (fexist(buffer)) {
      sprintf(buffer, "%s%c%s", od.field.md_grid_exe_path, SEPARATOR, GRIN_EXE);
      if (fexist(buffer)) {
        found = 1;
      }
    }
  }
  if (found) {
    tee_printf(&od, "MD GRID binaries are in:\n"
      "%s\n", od.field.md_grid_exe_path);
  }
  else {
    tee_printf(&od, "MD GRID binaries could not be found.\n");
  }
  tee_printf(&od, "The path to MD GRID binaries can be %s "
    "through the O3_MD_GRID_PATH environment variable or "
    "by the \"env md_grid_path\" keyword.\n\n",
    (found ? "changed" : "set"));
  #else
  memset(od.qmd.tinker_exe_path, 0, BUF_LEN);
  if (is_in_path(TINKER_DYNAMIC_EXE, buffer)) {
    strncpy(od.qmd.tinker_exe_path, buffer, BUF_LEN - 2);
    get_dirname(od.qmd.tinker_exe_path);
  }
  else {
    strcpy(od.qmd.tinker_exe_path, bin);
  }
  if ((tinker_path_string = getenv("O3_TINKER_PATH"))) {
    if (tinker_path_string[0]) {
      strncpy(od.qmd.tinker_exe_path, tinker_path_string, BUF_LEN - 2);
      absolute_path(od.qmd.tinker_exe_path);
    }
  }
  found = 0;
  if (dexist(od.qmd.tinker_exe_path)) {
    sprintf(buffer, "%s%c%s", od.qmd.tinker_exe_path, SEPARATOR, TINKER_DYNAMIC_EXE);
    if (fexist(buffer)) {
      sprintf(buffer, "%s%c%s", od.qmd.tinker_exe_path, SEPARATOR, TINKER_OPTIMIZE_EXE);
      if (fexist(buffer)) {
        sprintf(buffer, "%s%c%s", od.qmd.tinker_exe_path, SEPARATOR, TINKER_MINIMIZE_EXE);
        if (fexist(buffer)) {
          found = 1;
        }
      }
    }
  }
  if (found) {
    tee_printf(&od, "TINKER binaries are in:\n"
      "%s\n", od.qmd.tinker_exe_path);
  }
  else {
    tee_printf(&od, "TINKER binaries could not be found.\n");
  }
  tee_printf(&od, "The path to TINKER binaries can be %s "
    "through the O3_TINKER_PATH environment variable or "
    "by the \"env tinker_path\" keyword.\n\n",
    (found ? "changed" : "set"));
  memset(od.align.pharao_exe, 0, BUF_LEN);
  found = 1;
  result = 1;
  if (od.field.babel_exe_path[0]) {
    sprintf(od.align.pharao_exe, "%s%c%s",
      od.field.babel_exe_path, SEPARATOR, PHARAO_EXE);
    found = fexist(od.align.pharao_exe);
    if (!found) {
      if ((found = is_in_path(PHARAO_EXE, buffer))) {
        strncpy(od.align.pharao_exe, buffer, BUF_LEN - 2);
      }
    }
    if ((pharao_string = getenv("O3_PHARAO"))) {
      if (pharao_string[0]) {
        strncpy(od.align.pharao_exe, pharao_string, BUF_LEN - 2);
        absolute_path(od.align.pharao_exe);
        found = fexist(od.align.pharao_exe);
      }
    }
    if (found) {
      result = check_pharao(&od, bin);
      switch (result) {
        case OUT_OF_MEMORY:
        fprintf(stderr, E_OUT_OF_MEMORY, O3_FAILED);
        if (od.prompt) {
          fprintf(stderr, E_PROGRAM_EXIT);
          fflush(stderr);
          #ifndef WIN32
          fgets(buffer, BUF_LEN - 2, od.in);
          #else
          ReadConsole(od.hInput, buffer, BUF_LEN - 2, &n, NULL);
          #endif
          reset_user_terminal(&od);
        }
        return -1;

        case CANNOT_READ_TEMP_FILE:
        fprintf(stderr, E_FILE_CANNOT_BE_OPENED_FOR_READING,
          od.file[TEMP_LOG]->name, "");
        break;

        case CANNOT_WRITE_TEMP_FILE:
        fprintf(stderr, E_FILE_CANNOT_BE_OPENED_FOR_WRITING,
          od.file[TEMP_OUT]->name, "");
        break;

        case BABEL_NOT_WORKING:
        fprintf(stderr, E_PHARAO_NOT_WORKING);
        if ((od.file[TEMP_LOG]->handle = fopen(od.file[TEMP_LOG]->name, "rb"))) {
          while (fgets(buffer, BUF_LEN, od.file[TEMP_LOG]->handle)) {
            buffer[BUF_LEN - 1] = '\0';
            fprintf(stderr, "%s", buffer);
          }
          fclose(od.file[TEMP_LOG]->handle);
          od.file[TEMP_LOG]->handle = NULL;
        }
        break;
      }
    }
  }
  if (found && (!result)) {
    tee_printf(&od, "The PHARAO binary is:\n"
      "%s\n", od.align.pharao_exe);
  }
  else {
    tee_printf(&od, "A ");
    if (found && result) {
      tee_printf(&od, "working ");
    }
    tee_printf(&od, E_PHARAO_MISSING_OR_TOO_OLD, "");
    memset(od.align.pharao_exe, 0, BUF_LEN);
  }
  tee_printf(&od, "The path to PHARAO can be %s through the "
    "O3_PHARAO environment variable or "
    "by the \"env pharao\" keyword.\n\n",
    ((found && (!result)) ? "changed" : "set"));
  #endif
  #ifdef O3Q
  od.gnuplot.use_gnuplot = od.prompt;
  memset(od.gnuplot.gnuplot_exe, 0, BUF_LEN);
  if (is_in_path(GNUPLOT_EXE, buffer)) {
    strncpy(od.gnuplot.gnuplot_exe, buffer, BUF_LEN - 2);
  }
  if ((gnuplot_string = getenv("O3_GNUPLOT"))) {
    if (gnuplot_string[0]) {
      strncpy(od.gnuplot.gnuplot_exe, gnuplot_string, BUF_LEN - 2);
      absolute_path(od.gnuplot.gnuplot_exe);
      od.gnuplot.use_gnuplot = 1;
    }
    else {
      od.gnuplot.use_gnuplot = 0;
    }
  }
  found = 0;
  if (od.gnuplot.use_gnuplot) {
    if ((found = fexist(od.gnuplot.gnuplot_exe))) {
      tee_printf(&od, "The GNUPLOT binary is:\n"
        "%s\n", od.gnuplot.gnuplot_exe);
    }
    else {
      tee_printf(&od, "The GNUPLOT binary could not be found.\n");
      od.gnuplot.use_gnuplot = 0;
    }
  }
  else {
    tee_printf(&od, "GNUPLOT will not be used.\n");
  }
  tee_printf(&od, "The path to GNUPLOT can be %s through the "
    "O3_GNUPLOT environment variable or "
    "by the \"env gnuplot\" keyword.\n\n",
    (found ? "changed" : "set"));
  #endif
  od.pymol.use_pymol = od.prompt;
  memset(od.pymol.pymol_exe, 0, BUF_LEN);
  if (is_in_path(PYMOL_EXE, buffer)) {
    strncpy(od.pymol.pymol_exe, buffer, BUF_LEN - 2);
  }
  if ((pymol_string = getenv("O3_PYMOL"))) {
    if (pymol_string[0]) {
      strncpy(od.pymol.pymol_exe, pymol_string, BUF_LEN - 2);
      absolute_path(od.pymol.pymol_exe);
      od.pymol.use_pymol = 1;
    }
    else {
      od.pymol.use_pymol = 0;
    }
  }
  found = 0;
  if (od.pymol.use_pymol) {
    if ((found = fexist(od.pymol.pymol_exe))) {
      tee_printf(&od, "The PyMOL binary is:\n"
        "%s\n", od.pymol.pymol_exe);
      strcat(od.pymol.pymol_exe, " " PYMOL_ARGS);
    }
    else {
      tee_printf(&od, "The PyMOL binary could not be found.\n");
      od.pymol.use_pymol = 0;
    }
  }
  else {
    tee_printf(&od, "PyMOL will not be used.\n");
  }
  tee_printf(&od, "The path to PyMOL can be %s through the "
    "O3_PYMOL environment variable or "
    "by the \"env pymol\" keyword.\n\n",
    (found ? "changed" : "set"));
  #ifdef WIN32
  wVersionRequested = MAKEWORD(2, 2);
  error = WSAStartup(wVersionRequested, &wsaData);
  if (!error) {
    error = ((LOBYTE(wsaData.wVersion) != 2)
      || (HIBYTE(wsaData.wVersion) != 2));
    if (error) {
      WSACleanup();
    }
  }
  if (error) {
    od.jmol.port = -1;
  }
  #endif
  if (od.prompt && (!(od.pymol.use_pymol))
    && (od.jmol.port != -1)) {
    memset(od.jmol.jmol_exe, 0, BUF_LEN);
    if (is_in_path(JMOL_EXE, buffer)) {
      strncpy(od.jmol.jmol_exe, buffer, BUF_LEN - 2);
    }
    if ((jmol_string = getenv("O3_JMOL"))) {
      if (jmol_string[0]) {
        strncpy(od.jmol.jmol_exe, jmol_string, BUF_LEN - 2);
        absolute_path(od.jmol.jmol_exe);
        od.jmol.use_jmol = 1;
        od.jmol.port = JMOL_FIRST_PORT;
      }
      else {
        od.jmol.use_jmol = 0;
        od.jmol.port = -1;
      }
    }
    found = 0;
    if (od.jmol.use_jmol) {
      if ((found = fexist(od.jmol.jmol_exe))) {
        tee_printf(&od, "The Jmol binary is:\n"
          "%s\n", od.jmol.jmol_exe);
        strcat(od.jmol.jmol_exe, " " JMOL_ARGS);
      }
      else {
        tee_printf(&od, "The Jmol binary could not be found.\n");
        od.jmol.use_jmol = 0;
      }
    }
    else {
      tee_printf(&od, "Jmol will not be used.\n");
    }
    tee_printf(&od, "The path to Jmol can be %s through the "
      "O3_JMOL environment variable or "
      "by the \"env jmol\" keyword.\n\n",
      (found ? "changed" : "set"));
  }
  save_ram_string = getenv("O3_SAVE_RAM");
  if (save_ram_string && (!strncasecmp(save_ram_string, "y", 1))) {
    od.save_ram = 1;
    tee_printf(&od, "Page files will be used to "
      "minimize physical RAM usage.\n\n"); 
  }
  tee_flush(&od);
  if (!get_current_time(current_time)) {
    tee_printf(&od, "Job started on %s\n", current_time);
  }
  tee_flush(&od);
  result = 0;
  if (!(cli_args.prompt) && (cli_args.input)) {
    /*
    if an input script is being sourced,
    make a dry run to check out its consistency
    and avoid wasting CPU time
    */
    result = parse_input(&od, od.in, DRY_RUN);
    if (result) {
      tee_printf(&od, "\n"
        "Your input did not pass the consistency check.\n"
        "No computation was performed.\n");
      od.valid = 0;
      od.object_num = 0;
      od.field_num = 0;
    }
  }
  if (!result) {
    result = parse_input(&od, od.in, cli_args.prompt);
    if (!get_current_time(current_time)) {
      tee_printf(&od, "\n\n"
        "Job finished on %s\n", current_time);
    }
    if (result) {
      tee_printf(&od, "Abnormal exit.\n");
    }  
    else {
      tee_printf(&od, "Successful completion.\n");
    }
  }
  tee_flush(&od);
  od.out = NULL;
  if (current_dir[0]) {
    chdir(current_dir);
  }
  if (od.pymol.pymol_handle) {
    strcpy(buffer, "quit\n");
    FWRITE_WRAP(od.pymol.pymol_handle, buffer, &n);
    FFLUSH_WRAP(od.pymol.pymol_handle);
    free_proc_env(od.pymol.proc_env);
  }
  if (od.jmol.use_jmol && (od.jmol.port != -1)) {
    send_jmol_command(&od, "exitJmol");
    free_proc_env(od.jmol.proc_env);
    #ifndef WIN32
    close(od.jmol.sock);
    #else
    closesocket(od.jmol.sock);
    #endif
  }
  #ifdef WIN32
  if (hIcon) {
    DestroyIcon(hIcon);
  }
  if (od.jmol.port != -1) {
    WSACleanup();
  }
  #endif
  if (od.prompt && result) {
    fprintf(stderr, E_PROGRAM_EXIT);
    fflush(stderr);
    #ifndef WIN32
    fgets(buffer, BUF_LEN - 2, od.in);
    #else
    ReadConsole(od.hInput, buffer, BUF_LEN - 2, &n, NULL);
    #endif
  }
  reset_user_terminal(&od);
  close_files(&od, 0);
  if ((!result) || (!(od.debug))) {
    remove_temp_files(PACKAGE_CODE);
  }
  if (od.mel.line) {
    #if (defined HAVE_EDITLINE_FUNCTIONALITY && defined HAVE_GNU_READLINE)
    rl_free(od.mel.line);
    od.mel.line = NULL;
    #endif
    #ifndef HAVE_EDITLINE_FUNCTIONALITY
    if (have_editline && _dlsym_rl_free) {
      rl_free(od.mel.line);
      od.mel.line = NULL;
    }
    #endif
  }
  free_mem(&od);
  #ifndef HAVE_EDITLINE_FUNCTIONALITY
  if (dl_handle) {
    #ifndef WIN32
    dlclose(dl_handle);
    #else
    FreeLibrary(dl_handle);
    #endif
  }
  #endif

  return result;
}
