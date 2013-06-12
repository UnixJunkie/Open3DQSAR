/*

fn_complete.c

is part of:

MinGWEditLine
Copyright 2010-2012 Paolo Tosco <paolo.tosco@unito.it>
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:

    * Redistributions of source code must retain the above copyright
    notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
    notice, this list of conditions and the following disclaimer in the
    documentation and/or other materials provided with the distribution.
    * Neither the name of MinGWEditLine nor the name of its contributors
    may be used to endorse or promote products derived from this software
    without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/


#define _UNICODE
#define UNICODE

#include <editline/readline.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <wctype.h>
#include <tchar.h>


BOOL _el_replace_char(wchar_t *string, wchar_t src, wchar_t dest)
{
	int i = 0;
	BOOL replaced = FALSE;
	
	
	while (string[i]) {
		if (string[i] == src) {
			replaced = TRUE;
			string[i] = dest;
		}
		++i;
	}
	
	return replaced;
}


wchar_t *_el_get_compl_text(int *start, int *end)
{
	wchar_t *text = NULL;
	int i;
	int nearest_quote = 0;
	int open_quote = 0;


	i = rl_point - 1;
	/*
	start scanning the line backwards from
	the current logical cursor position;
	record the position of the first quotes
	found, and count all quotes
	*/
	while (i >= 0) {
		if (_el_line_buffer[i] == _T('\"')) {
			if (!open_quote) {
				nearest_quote = i;
			}
			++open_quote;
		}
		--i;
	}
	/*
	are quotes in odd number? Then the user
	has just opened quotes, so the completion
	should be considered as starting from that
	position
	*/
	if (open_quote % 2) {
		*start = nearest_quote + 1;
	}
	/*
	otherwise (no quotes or quotes in even number)
	just look for the nearest space or beginning
	of line
	*/
	else {
		i = rl_point;
		while (i && (!wcschr((_el_completer_word_break_characters
			? _el_completer_word_break_characters
			: _el_basic_word_break_characters),
			_el_line_buffer[i - 1]))) {
			--i;
		}
		*start = i;
	}
	*end = rl_point;
	i = *end - *start;
	if ((text = (wchar_t *)malloc((i + 1) * sizeof(wchar_t)))) {
		if (i) {
			memcpy(text, &_el_line_buffer[*start], i * sizeof(wchar_t));
		}
		text[i] = '\0';
	}
	
	return text;
}


char **rl_completion_matches(const char *text, char *entry_func(const char *, int))
{
	char **array = NULL;
	char *entry;
	int i = 0;
	int j = 0;
	int n_loop = 0;
	
	
	while (n_loop < 2) {
		while (1) {
			entry = (n_loop
				? (rl_completion_entry_function
				? rl_completion_entry_function(text, j)
				: rl_filename_completion_function(text, j))
				: entry_func(text, i));
			if (!(array = (char **)realloc(array, (i + 2) * sizeof(char *)))) {
				return NULL;
			}
			array[i] = NULL;
			array[i + 1] = NULL;
			if (entry) {
				if (!(array[i] = strdup(entry))) {
					return NULL;
				}
				free(entry);
			}
			else {
				break;
			}
			++i;
			if (n_loop) {
				++j;
			}
		}
		if ((!i) && (!rl_attempted_completion_over)) {
			++n_loop;
		}
		else {
			break;
		}
	}
	/*
	if (i == 1) {
		j = strlen(array[0]);
		if (!(array[0] = realloc(array[0], j + 2))) {
			return NULL;
		}
		array[0][j] = rl_completion_append_character;
		array[0][j + 1] = '\0';
	}
	*/
	
	return array;
}


char *rl_filename_completion_function(const char *text, int state)
{
	char *compl_mb = NULL;
	wchar_t *separator;
	wchar_t n;
	int i;
	int j;
	int len;
	int dir_name_len;
	int nearest_quote = 0;
	int open_quote;
	int close_quote;
	int end_of_list = 0;
	HANDLE *dHandle = NULL;
	WIN32_FIND_DATA filedata;
	
	
	/*
	if a new TAB-completion cycle has been started
	*/
	if (!state) {
		i = rl_point;
		open_quote = 0;
		/*
		start scanning the line backwards from
		the current logical cursor position;
		record the position of the first quotes
		found, and count all quotes
		*/
		while (i >= 0) {
			if (_el_line_buffer[i] == _T('\"')) {
				if (!open_quote) {
					nearest_quote = i;
				}
				++open_quote;
			}
			--i;
		}
		/*
		are quotes in odd number? Then the user
		has just opened quotes, so the filename
		should be considered as starting from that
		position
		*/
		if (open_quote % 2) {
			i = nearest_quote;
		}
		/*
		otherwise (no quotes or quotes in even number)
		just look for the nearest space or beginning
		of line
		*/
		else {
			/*
			if there are are quotes in even number
			and no spaces in between, then
			remove all of them since they are useless
			*/
			if (open_quote) {
				i = rl_point;
				j = 0;
				while ((!j) && (i >= 0) && (_el_line_buffer[i] != _T('\"'))) {
					j = iswspace(_el_line_buffer[i]);
					--i;
				}
				if (!j) {
					/*
					remove excess quotes from _el_line_buffer
					and from _el_text
					*/
					if (_el_delete_char(VK_BACK, open_quote)) {
						return NULL;
					}
					if (_el_text) {
						j = wcslen(_el_text) - open_quote;
						if (j >= 0) {
							_el_text[j] = _T('\0');
						}
					}
				}
				i = rl_point;
			}
			else {
				i = rl_point;
				while (i && (!wcschr(_el_basic_file_break_characters,
					_el_line_buffer[i - 1]))) {
					--i;
				}
			}
		}
		/*
		i is the position of the first character
		from which file completion will be attempted
		now we shall check if the
		user has already written a (partial) path
		check if the first character is a quote
		*/
		open_quote = (((_el_line_buffer[i] == _T('\"'))
			&& (i != rl_point)) ? 1 : 0);
		/*
		if the first character of the path is a
		slash or backslash, then abort completion,
		since we do not want to take into account
		UNIX-like paths beginning with slashes, but
		only native Windows paths
		*/
		/*
		if (_el_line_buffer[i + open_quote]
			&& strchr("/\\", (int)(_el_line_buffer[i + open_quote]))) {
			return NULL;
		}
		*/
		/*
		j is the position of the character just before the cursor
		*/
		j = rl_point - 1;
		/*
		allocate memory for the partially written path
		(may also be an empty string, that is just '\0')
		*/
		if (!(_el_old_arg = (wchar_t *)realloc
			(_el_old_arg, (j + 2 - i) * sizeof(wchar_t)))) {
			return NULL;
		}
		memset(_el_old_arg, 0, (j + 2 - i) * sizeof(wchar_t));
		if (j + 1 - i) {
			memcpy(_el_old_arg, &_el_line_buffer[i], (j + 1 - i) * sizeof(wchar_t));
		}
		separator = NULL;
		_el_replace_char(_el_old_arg, _T('/'), _T('\\'));
		/*
		is the partially written path pointing
		to a directory?
		*/
		memset(&filedata, 0, sizeof(WIN32_FIND_DATA));
		dHandle = FindFirstFile(_el_old_arg, &filedata);
		if ((dHandle == INVALID_HANDLE_VALUE)
			|| (filedata.dwFileAttributes != FILE_ATTRIBUTE_DIRECTORY)) {
			/*
			if not, scan backwards searching
			for the nearest slash or backslash
			*/
			while (j > i) {
				if ((separator = wcschr
					(_T("/\\"), _el_line_buffer[j]))) {
					break;
				}
				--j;
			}
			if (separator) {
				/*
				if a slash or backslash has been found,
				then j is its position
				*/
				if (!(_el_dir_name = (wchar_t *)realloc(_el_dir_name,
					(j + 2 - (i + open_quote)) * sizeof(wchar_t)))) {
					return NULL;
				}
				memset(_el_dir_name, 0, (j + 2 - (i + open_quote)) * sizeof(wchar_t));
				/*
				copy the path into dir_name
				until the last slash or backslash
				including the latter
				*/
				memcpy(_el_dir_name, &_el_line_buffer[i + open_quote],
					(j + 1 - (i + open_quote)) * sizeof(wchar_t));
				if (!(_el_file_name = (wchar_t *)realloc(_el_file_name,
					(rl_point - j) * sizeof(wchar_t)))) {
					return NULL;
				}
				memset(_el_file_name, 0, (rl_point - j) * sizeof(wchar_t));
				/*
				copy what comes after
				the last slash or backslash
				*/
				memcpy(_el_file_name, &_el_line_buffer[j + 1],
					(rl_point - j - 1) * sizeof(wchar_t));
			}
			else {
				/*
				if a slash or backslash has not been found,
				then we look in the current directory
				*/
				if (!(_el_dir_name = (wchar_t *)realloc
					(_el_dir_name, 3 * sizeof(wchar_t)))) {
					return NULL;
				}
				wcscpy(_el_dir_name, _T("./"));
				if (!(_el_file_name = (wchar_t *)realloc(_el_file_name,
					(rl_point - (i + open_quote) + 1) * sizeof(wchar_t)))) {
					return NULL;
				}
				memset(_el_file_name, 0, (rl_point - (i + open_quote) + 1) * sizeof(wchar_t));
				if (rl_point - (i + open_quote)) {
					memcpy(_el_file_name, &_el_line_buffer[i + open_quote],
						(rl_point - (i + open_quote)) * sizeof(wchar_t));
				}
			}
		}
		else {
			/*
			if the path was directly pointing to a directory
			*/
			len = wcslen(_el_old_arg);
			/*
			check if a close quote is included
			*/
			close_quote = ((_el_old_arg[len - 1] == _T('\"')) ? 1 : 0);
			/*
			check if a final slash or backslash is included
			*/
			n = (!wcschr(_T("/\\"),
				_el_old_arg[len - close_quote - 1])) ? _T('\\') : _T('\0');
			if (!(_el_dir_name = (wchar_t *)realloc
				(_el_dir_name, (len + 4) * sizeof(wchar_t)))) {
				return NULL;
			}
			memset(_el_dir_name, 0, (len + 4) * sizeof(wchar_t));
			/*
			if the slash or backslash was not included,
			then add it to the directory name
			*/
			wsprintf(_el_dir_name, _T("%s%c"), &_el_old_arg[open_quote], n);
		}
		if (dHandle != INVALID_HANDLE_VALUE) {
			CloseHandle(dHandle);
		}
		/*
		replace all slashes in the directory path
		with backslashes
		*/
		_el_replace_char(_el_dir_name, _T('/'), _T('\\'));
		/*
		if the directory does not open,
		no completion is possible
		*/
		if (!(_el_wide = (wchar_t *)realloc(_el_wide,
			(wcslen(_el_dir_name) + 3) * sizeof(wchar_t)))) {
			return NULL;
		}
		wsprintf(_el_wide, _T("%s\\*"), _el_dir_name);
		dHandle = FindFirstFile(_el_wide, &filedata);
		if (dHandle == INVALID_HANDLE_VALUE) {
			return NULL;
		}
		_el_n_compl = 0;
		/*
		read the next directory entry
		*/
		while (FindNextFile(dHandle, &filedata)) {
			/*
			skip "." and ".." if present
			(this is probably unnecessary under Windows)
			*/
			if ((!(wcscmp(filedata.cFileName, _T("."))))
				|| (!(wcscmp(filedata.cFileName, _T(".."))))) {
				continue;
			}
			/*
			if the name of this entry matches the
			partial filename written by the user,
			increase the count of matching entries
			*/
			if (_el_check_root_identity(_el_file_name, filedata.cFileName)) {
				++_el_n_compl;
			}
		}
		CloseHandle(dHandle);
		/*
		if the directory is empty,
		no completion is possible
		*/
		if (!_el_n_compl) {
			return NULL;
		}
		/*
		after counting, let's allocate a string array,
		re-read and store directory entries in the array
		*/
		dHandle = FindFirstFile(_el_wide, &filedata);
		if (dHandle == INVALID_HANDLE_VALUE) {
			return NULL;
		}
		_el_free_array(_el_compl_array);
		if (!(_el_compl_array = (wchar_t **)
			malloc((_el_n_compl + 1) * sizeof(wchar_t *)))) {
			return NULL;
		}
		memset(_el_compl_array, 0, (_el_n_compl + 1) * sizeof(wchar_t *));
		_el_n_compl = 0;
		dir_name_len = wcslen(_el_dir_name);
		while (FindNextFile(dHandle, &filedata)) {
			/*
			skip "." and ".." if present
			(this is probably unnecessary under Windows)
			*/
			if ((!wcscmp(filedata.cFileName, _T(".")))
				|| (!wcscmp(filedata.cFileName, _T("..")))) {
				continue;
			}
			if (!_el_check_root_identity(_el_file_name, filedata.cFileName)) {
				continue;
			}
			len = wcslen(filedata.cFileName);
			if (!(_el_compl_array[_el_n_compl] = (wchar_t *)
				malloc((dir_name_len + len + 3) * sizeof(wchar_t)))) {
				return NULL;
			}
			/*
			prepend the directory path to the current filename
			*/
			wsprintf(_el_compl_array[_el_n_compl], _T("%s%s"),
				_el_dir_name, filedata.cFileName);
			/*
			if the filename contains spaces, then enclose
			it between quotes
			*/
			if (wcschr(_el_compl_array[_el_n_compl], _T(' '))) {
				if (_el_compl_array[_el_n_compl][0] != _T('\"')) {
					len = wcslen(_el_compl_array[_el_n_compl]);
					memmove(&_el_compl_array[_el_n_compl][1],
						_el_compl_array[_el_n_compl], (len + 1) * sizeof(wchar_t));
					_el_compl_array[_el_n_compl][0] = _T('\"');
				}
				wcscat(_el_compl_array[_el_n_compl], _T("\""));
			}
			++_el_n_compl;
		}
		CloseHandle(dHandle);
		/*
		do an alphabetic sort on directory entries
		*/
		qsort(_el_compl_array, _el_n_compl, sizeof(wchar_t *),
			_el_fn_qsort_string_compare);
		/*
		we start cycling from the first entry
		*/
		_el_compl_index = 1;
	}
	else {
		/*
		if TAB has been pressed again after the first time,
		then we just need to cycle over directory entries
		*/
		++_el_compl_index;
		if (_el_compl_index == (_el_n_compl + 1)) {
			--_el_compl_index;
			end_of_list = 1;
		}
	}
	if (_el_file_name) {
		free(_el_file_name);
		_el_file_name = NULL;
	}
	
	return (end_of_list ? NULL : _el_w2mb
		(_el_compl_array[_el_compl_index - 1], &compl_mb));
}
