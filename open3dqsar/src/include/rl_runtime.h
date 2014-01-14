/*

rl_runtime.h

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


#ifndef WIN32
#include <dlfcn.h>
#endif
#ifndef HAVE_EDITLINE_FUNCTIONALITY
#ifdef __STDC__
typedef void *histdata_t;
#else
typedef char *histdata_t;
#endif

typedef struct _hist_entry {
  char *line;
  char *timestamp;
  histdata_t data;
} HIST_ENTRY;

#define add_history (*_dlsym_add_history)
#define next_history (*_dlsym_next_history)
#define previous_history (*_dlsym_previous_history)
#define read_history (*_dlsym_read_history)
#define using_history (*_dlsym_using_history)
#define write_history (*_dlsym_write_history)
#define readline (*_dlsym_readline)
#define rl_free (*_dlsym_rl_free)
#define rl_attempted_completion_function (*_dlsym_rl_attempted_completion_function)
#define rl_user_completion_entry_free_function (*_dlsym_rl_user_completion_entry_free_function)
#define rl_attempted_completion_over (*_dlsym_rl_attempted_completion_over)
#define rl_completion_append_character (*_dlsym_rl_completion_append_character)
#define rl_completion_matches (*_dlsym_rl_completion_matches)
#define rl_filename_completion_function (*_dlsym_rl_filename_completion_function)
#define rl_line_buffer (*_dlsym_rl_line_buffer)
#define rl_point (*_dlsym_rl_point)
#define rl_catch_signals (*_dlsym_rl_catch_signals)
#define rl_delete_text (*_dlsym_rl_delete_text)
#else
#ifdef HAVE_GNU_READLINE
#include <readline/readline.h>
#include <readline/history.h>
#else
#include <editline/readline.h>
#include <histedit.h>
#endif
#endif

extern void *(*_dlsym_add_history)();
extern HIST_ENTRY *(*_dlsym_next_history)();
extern HIST_ENTRY *(*_dlsym_previous_history)();
extern int (*_dlsym_read_history)();
extern int (*_dlsym_write_history)();
extern void (*_dlsym_using_history)();
extern char *(*_dlsym_readline)();
extern void (*_dlsym_rl_free)();
extern void *(*_dlsym_rl_attempted_completion_function);
extern void *(*_dlsym_rl_user_completion_entry_free_function);
extern int *_dlsym_rl_attempted_completion_over;
extern int *_dlsym_rl_completion_append_character;
extern char **(*_dlsym_rl_completion_matches)();
extern char *(*_dlsym_rl_filename_completion_function)();
extern char *(*_dlsym_rl_line_buffer);
extern int *_dlsym_rl_point;
extern int *_dlsym_rl_catch_signals;
extern int *_dlsym_rl_delete_text;
extern char have_editline;
extern char have_gnu_readline;
