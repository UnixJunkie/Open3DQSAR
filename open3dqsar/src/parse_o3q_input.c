/*

parse_o3q_input.c

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
#include <include/error_messages.h>
#include <include/ff_parm.h>
#ifdef WIN32
#include <include/nice_windows.h>
#endif

#ifdef HAVE_EDITLINE_FUNCTIONALITY
#include <editline/readline.h>
#include <histedit.h>
#endif


int parse_input(O3Data *od, FILE *input_stream, int run_type)
{
  char *ptr;
  char *context = NULL;
  char *multi_file_name[] = {
    "",
    "FORMATTED_CUBE",
    "UNFORMATTED_CUBE",
    "MOLDEN",
    "MOE_GRID",
    "GRID_ASCII",
    "COSMO"
  };
  char prompt[TITLE_LEN];
  char grid_fill[MAX_NAME_LEN];
  char history_file[BUF_LEN];
  char line_orig[BUF_LEN];
  char last_line[BUF_LEN];
  char whole_line[BUF_LEN];
  char regex_name[MAX_STATES][BUF_LEN];
  char file_basename[BUF_LEN];
  char name_list[BUF_LEN];
  char comma_hyphen_list[BUF_LEN];
  char buffer[BUF_LEN];
  char *turbodir = NULL;
  char *parameter = NULL;
  char *point;
  char *eof;
  char *default_probe[] = { "CR", "C3" };
  unsigned short attr = 0;
  int fail = 0;
  int command = 0;
  int nesting = 0;
  int i;
  int j;
  int found;
  int space_count;
  int len;
  int actual_len;
  int argn;
  int line_num;
  int overall_line_num = 0;
  int line_complete;
  int spaces;
  int result;
  int import_y_vars = 0;
  int pc;
  int max_pc;
  int mo;
  int format;
  int type = 0;
  int ref_field;
  int operation;
  int sign;
  int state = 0;
  int synonym = 0;
  int list_type = 0;
  int interpolate;
  int requested_endianness;
  int replace_object_name;
  int pc_axis[3];
  int multi_file_type;
  int quote;
  int quote_error;
  int nice_value;
  int prep_or_calc = 0;
  int from_file = 0;
  int skip_header;
  int options;
  int n_values;
  int label = 0;
  int *max_vary[MAX_FREE_FORMAT_PARAMETERS + 1];
  int *vary[MAX_FREE_FORMAT_PARAMETERS + 1];
  char *failed = NULL;
  long random_seed;
  double factor;
  double level;
  double weight;
  double outgap = 0.0;
  double trans[3];
  double rot[3];
  CharMat *arg;
  struct timeval start;
  struct timeval end;
  GridInfo temp_grid;
  VarCoord temp_varcoord;
  FileDescriptor **source;
  FileDescriptor *current_source;
  FileDescriptor *multi_file_fd = NULL;
  FileDescriptor log_fd;
  #ifdef HAVE_EDITLINE_FUNCTIONALITY
  FILE *touch;
  HIST_ENTRY *last_entry;
  #endif
  int active_struct_num;
  int groups = 0;
  int runs = 0;
  int percent_remove;
  int start_gnuplot = 0;
  int srd_collapse;
  int default_seeds[3];
  int seeds;
  int candidates = 0;
  int design_points;
  double critical_distance;
  double collapse_distance;
  double grid_size[3];
  double grid_diagonal;
  #ifndef WIN32
  int des = 0;
  pid_t pid1;
  struct rlimit rlp;
  #endif
  

  memset(pc_axis, 0, 3 * sizeof(int));
  memset(comma_hyphen_list, 0, BUF_LEN);
  memset(history_file, 0, BUF_LEN);
  memset(last_line, 0, BUF_LEN);
  memset(&temp_grid, 0, sizeof(GridInfo));
  memset(&log_fd, 0, sizeof(FileDescriptor));
  current_source = od->mel.source;
  while (current_source) {
    ++nesting;
    current_source = current_source->next;
  }
  #ifdef HAVE_EDITLINE_FUNCTIONALITY
  if (!(run_type & INTERACTIVE_RUN)) {
    /*
    if this is not an interactive session,
    then we are going to fetch lines by fgets
    so allocate memory for buffer
    */
    od->mel.line = realloc(od->mel.line, BUF_LEN);
    if (!(od->mel.line)) {
      tee_error(od, run_type, overall_line_num,
        E_OUT_OF_MEMORY, PARSE_INPUT_FAILED);
      return PARSE_INPUT_ERROR;
    }
  }
  else {
    /*
    if this is an interactive session,
    then check if the history file exists,
    otherwise create it
    */
    sprintf(history_file, "%s%c%s", od->home_dir, SEPARATOR, HISTORY_FILE);
    if (access(history_file, F_OK) == -1) {
      touch = fopen(history_file, "w+");
      if (touch) {
        fclose(touch);
      }
    }
    /*
    read history from file and get the last entry
    */
    using_history();
    read_history(history_file);
    last_entry = previous_history();
    next_history();
    if (last_entry) {
      strncpy(last_line, last_entry->line, BUF_LEN);
      last_line[BUF_LEN - 1] = '\0';
    }
  }
  #else
  od->mel.line = realloc(od->mel.line, BUF_LEN);
  if (!(od->mel.line)) {
    tee_error(od, run_type, overall_line_num,
      E_OUT_OF_MEMORY, PARSE_INPUT_FAILED);
    return PARSE_INPUT_ERROR;
  }
  #endif
  /*
  allocate an array of pointers to string buffers
  for command line arguments and values, then store
  the pointers in the MemList substructure
  */
  od->cimal.arg = alloc_char_matrix(od->cimal.arg, MAX_ARG, BUF_LEN);
  if (!(od->cimal.arg)) {
    tee_error(od, run_type, overall_line_num,
      E_OUT_OF_MEMORY, PARSE_INPUT_FAILED);
    return PARSE_INPUT_ERROR;
  }
  arg = od->cimal.arg;
  /*
  main parse loop
  read one line at a time
  */
  while (!fail) {
    /*
    start closing file handles used by previous
    commands (if any); make sure we do not close
    MAIN_INPUT and MAIN_OUTPUT files
    */
    for (i = 0; i < MAX_FILES; ++i) {
      if ((i != MAIN_INPUT) && (i != MAIN_OUTPUT)) {
        if (od->file[i]->handle) {
          fclose(od->file[i]->handle);
          od->file[i]->handle = NULL;
        }
      }
    }
    /*
    if we are running interactively then we are using
    readline which allocates a pointer to the entered
    line; the user is supposed to free it after use
    */
    #ifdef HAVE_EDITLINE_FUNCTIONALITY
    if (run_type & INTERACTIVE_RUN) {
      if (od->mel.line) {
        free(od->mel.line);
        od->mel.line = NULL;
      }
    }
    else {
      memset(od->mel.line, 0, BUF_LEN);
    }
    #else
    memset(od->mel.line, 0, BUF_LEN);
    #endif
    memset(line_orig, 0, BUF_LEN);
    line_complete = 0;
    line_num = 0;
    memset(whole_line, 0, BUF_LEN);
    sprintf(prompt, "%s%s", PACKAGE_NAME, SHORT_PROMPT);
    /*
    read line(s) until a line feed which is not
    preceded by a backslash is met
    */
    while (!line_complete) {
      if (run_type & INTERACTIVE_RUN) {
        /*
        if it is an interactive session we use readline
        */
        SET_INK(od, HIGHLIGHT_INK);
        #ifdef HAVE_EDITLINE_FUNCTIONALITY
        if (rl_line_buffer && rl_line_buffer[0]) {
          memset(rl_line_buffer, 0, strlen(rl_line_buffer));
        }
        od->mel.line = readline(prompt);
        if (rl_line_buffer) {
          rl_line_buffer[0] = '\0';
        }
        eof = od->mel.line;
        #else
        printf("%s", prompt);
        if ((eof = fgets(od->mel.line, BUF_LEN, input_stream))) {
          od->mel.line[BUF_LEN - 1] = '\0';
          remove_newline(od->mel.line);
        }
        #endif
        SET_INK(od, NORMAL_INK);
      }
      else {
        /*
        otherwise we use fgets
        */
        if ((eof = fgets(od->mel.line, BUF_LEN, input_stream))) {
          od->mel.line[BUF_LEN - 1] = '\0';
          remove_newline(od->mel.line);
        }
      }
      if (eof) {
        /*
        if line ends with backslash, then wait for
        more input; at the end, put the line together
        */
        len = strlen(od->mel.line);
        /*
        remove spaces eventually present after the last character
        (often present when doing copy/paste from scrolled konsole)
        */
        while (len && isspace(od->mel.line[len - 1])) {
          --len;
        }
        if (len && (od->mel.line[len - 1] == '\\')) {
          ++line_num;
          strcpy(prompt, SHORT_PROMPT);
          od->mel.line[len - 1] = '\0';
          if ((strlen(whole_line) + len) < BUF_LEN) {
            strcat(whole_line, od->mel.line);
          }
          else {
            strncat(whole_line, od->mel.line,
              BUF_LEN - strlen(whole_line));
            line_complete = 1;
          }
        }
        else {
          strncat(whole_line, od->mel.line,
            BUF_LEN - strlen(whole_line));
          line_complete = 1;
        }
      }
      else {
        break;
      }
    }
    if (!eof) {
      /*
      if we are running interactively and the user
      pushes CTRL-D, print a new prompt
      this will prevent editline from quitting
      us upon window resizing
      */
      if (run_type & INTERACTIVE_RUN) {
        printf("\n");
        continue;
      }
      else {
        break;
      }
    }
    whole_line[BUF_LEN - 1] = '\0';
    i = 0;
    quote = 0;
    /*
    pretreat the line:
    - removing comments
    - converting tabs into spaces
    - converting ' ' into ';' provided it is not
      enclosed between quotes
    - converting '\ ' into ' ' provided it is not
      enclosed between quotes
    */
    len = strlen(whole_line);
    while (i < len) {
      if (isspace(whole_line[i]) && (!quote)) {
        whole_line[i] = ';';
      }
      if (!isprint(whole_line[i])) {
        for (j = i; j < len; ++j) {
          whole_line[j] = whole_line[j + 1];
        }
        --len;
        continue;
      }
      if (whole_line[i] == '#') {
        whole_line[i] = '\0';
        break;
      }
      if (whole_line[i] == '\"') {
        quote = 1 - quote;
      }
      if ((whole_line[i] == '\\') && (!quote)) {
        if ((i + 1) < len) {
          if (whole_line[i + 1] == ' ') {
            for (j = i; j < len; ++j) {
              whole_line[j] = whole_line[j + 1];
            }
            --len;
            ++i;
            continue;
          }
        }
      }
      ++i;
    }
    /*
    if quotes were opened but never closed issue an error
    */
    quote_error = 0;
    if (quote) {
      tee_error(od, run_type, overall_line_num, "Error: missing quote\n\n");
      quote_error = 1;
    }
    argn = 0;
    
    /*
    copy command line arguments in the "arg"
    array until MAX_ARG is reached
    */
    ptr = strtok_r(whole_line, ";", &context);
    while (ptr && (argn < MAX_ARG)) {
      strcpy(arg->me[argn], ptr);
      if (argn) {
        strcat(line_orig, ";");
      }
      strcat(line_orig, ptr);
      ++argn;
      ptr = strtok_r(NULL, ";", &context);
    }
    len = strlen(line_orig);
    i = 0;
    quote = 0;
    /*
    regenerate the original line:
    - converting ';' into ' ' provided it is not
      enclosed between quotes
    - converting ' ' into '\ ' provided it is not
      enclosed between quotes
    so basically we remove from the entered line
    comments and multiple spaces
    */
    while (i < len) {
      if (whole_line[i] == '\"') {
        quote = 1 - quote;
      }
      if ((line_orig[i] == ' ') && (!quote)) {
        for (j = len; j >= i; --j) {
          line_orig[j + 1] = line_orig[j];
        }
        line_orig[i] = '\\';
        ++i;
        ++len;
      }
      else if ((line_orig[i] == ';') && (!quote)) {
        line_orig[i] = ' ';
      }
      ++i;
    }
    if (line_orig[0] && (run_type & INTERACTIVE_RUN)
      && strcmp(last_line, line_orig)) {
      #ifdef HAVE_EDITLINE_FUNCTIONALITY
      add_history(line_orig);
      write_history(history_file);
      #endif
      strcpy(last_line, line_orig);
    }
    /*
    if a quote was missing it is time to bomb out,
    but at least the user finds the line in history
    so he can correct it easily
    */
    if (quote_error) {
      fail = !(run_type & INTERACTIVE_RUN);
      continue;
    }
    /*
    if there are no arguments (because the command line
    is empty or constituted by comments only), we need
    to count it or pretreatment errors will not refer to
    the correct text line
    */
    if (!argn) {
      ++overall_line_num;
      continue;
    }
    od->argn = argn;
    /*
    the first argument is the command itself
    */
    if ((parameter = get_args(od, NULL))) {
      tee_error(od, run_type, overall_line_num, "Error parsing "
        "the following command:\n\n"
        "> %s\n\n"
        "Please check your input near:\n\n", line_orig);
      point = strstr(line_orig, parameter);
      if (point) {
        spaces = 0;
        space_count = 0;
        while (((point - spaces) > line_orig)
          && (spaces <= MAX_SPACES)) {
          if (*(point - spaces) == ' ') {
            ++space_count;
          }
          if (space_count == 2) {
            break;
          }
          ++spaces;
        }
        i = 0;
        space_count = 0;
        while (point[i] && (point[i] != '\n')
          && (point[i] != '\r')
          && (i <= MAX_SPACES)) {
          if (point[i] == ' ') {
            ++space_count;
          }
          if (space_count == 2) {
            break;
          } 
          ++i;
        }
        point[i] = '\0';
        if (spaces > MAX_SPACES) {
          tee_printf(od, "[...] ");
        }
        tee_printf(od, "%s", (char *)(point - spaces));
        if (i > MAX_SPACES) {
          tee_printf(od, "[...] ");
        }
        tee_printf(od, "\n");
        for (i = 0; i < spaces; ++i) {
          tee_printf(od, " ");
        }
        tee_printf(od, "^\n");
      }
      fail = !(run_type & INTERACTIVE_RUN);
      continue;
    }
    overall_line_num += (line_num + 1);

    if (!(od->file[TEMP_OUT]->name[0])) {
      result = open_temp_file(od, od->file[TEMP_OUT], "temp_out");
      if (result) {
        tee_error(od, run_type, overall_line_num,
          E_TEMP_FILE_CANNOT_BE_OPENED_FOR_WRITING,
          od->file[TEMP_OUT]->name, IMPORT_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      if (od->file[TEMP_OUT]->handle) {
        fclose(od->file[TEMP_OUT]->handle);
        od->file[TEMP_OUT]->handle = NULL;
      }
    }
    if (!(od->file[TEMP_LOG]->name[0])) {
      result = open_temp_file(od, od->file[TEMP_LOG], "temp_log");
      if (result) {
        tee_error(od, run_type, overall_line_num,
          E_TEMP_FILE_CANNOT_BE_OPENED_FOR_WRITING,
          od->file[TEMP_LOG]->name, IMPORT_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      if (od->file[TEMP_LOG]->handle) {
        fclose(od->file[TEMP_LOG]->handle);
        od->file[TEMP_LOG]->handle = NULL;
      }
    }
    if (!strcasecmp(arg->me[0], "box")) {
      if (!(od->valid & SDF_BIT)) {
        tee_error(od, run_type, overall_line_num,
          E_IMPORT_MOLFILE_FIRST, BOX_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      if (od->field_num) {
        tee_error(od, run_type, overall_line_num,
          E_CANNOT_ALTER_GRID_BOX,
          BOX_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      memset(grid_fill, 0, MAX_NAME_LEN);
      temp_grid.step[0] = 1.0;
      outgap = 5.0;
      from_file = 0;
      if ((parameter = get_args(od, "file"))) {
        grid_fill[0] = 1;
        for (i = 0; i < 6; ++i) {
          grid_fill[i] = 1;
        }
        outgap = -1.0;
        strcpy(od->file[ASCII_IN]->name, parameter);
        if (!(od->file[ASCII_IN]->handle =
          fopen(od->file[ASCII_IN]->name, "rb"))) {
          tee_error(od, run_type, overall_line_num,
            E_FILE_CANNOT_BE_OPENED_FOR_READING,
            od->file[ASCII_IN]->name, BOX_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        from_file = 1;
      }
      else {
        for (result = 0, i = 0; i < 3; ++i) {
          sprintf(buffer, "%c_start", i + 'x');
          if ((parameter = get_args(od, buffer))) {
            grid_fill[i] = 1;
            sscanf(parameter, "%f", &(temp_grid.start_coord[i]));
          }
          synonym = 0;
          sprintf(buffer, "%c_end", i + 'x');
          if ((parameter = get_args(od, buffer))) {
            grid_fill[i + 3] = 1;
            sscanf(parameter, "%f", &(temp_grid.end_coord[i]));
            ++synonym;
          }
          sprintf(buffer, "%c_nodes", i + 'x');
          if ((parameter = get_args(od, buffer))) {
            grid_fill[i + 3] = 1;
            sscanf(parameter, "%d", &(temp_grid.nodes[i]));
            ++synonym;
          }
          if (synonym > 1) {
            tee_error(od, run_type, overall_line_num,
              "Conflicting grid definitions were detected; "
              "please specify either %c_end or %c_nodes.\n",
              i + 'x', i + 'x');
            result = 1;
          }
        }
        if (result) {
          tee_error(od, run_type, overall_line_num,
            "%s", BOX_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        if ((parameter = get_args(od, "step"))) {
          grid_fill[6] = 1;
          sscanf(parameter, "%f", &(temp_grid.step[0]));
        }
        if ((parameter = get_args(od, "outgap"))) {
          grid_fill[7] = 1;
          sscanf(parameter, "%lf", &outgap);
        }
        if (grid_fill[7]) {
          /*
          if the outgap parameter is given, no other parameters
          should be given except step (which defaults to 1.0)
          */
          for (i = 0, j = 0; ((i <= 5) && !j); ++i) {
            j += (int)grid_fill[i];
          }
          if (j) {
            tee_error(od, run_type, overall_line_num,
              "When the outgap parameter is input, "
              "only the step parameter may be "
              "present.\n%s",
              BOX_FAILED);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
          outgap = (double)safe_rint(fabs(outgap) * 1.0e04) / 1.0e04;
        }
        else {
          for (i = 0, j = 0; i < 6; ++i) {
            j += (int)grid_fill[i];
          }
          /*
          if the outgap parameter is not given, full grid
          start and end coordinates should be supplied,
          or none; in the latter case the default outgap (5.0)
          will be used
          */
          if (j && (j < 6)) {
            tee_error(od, run_type, overall_line_num,
              "Full grid box coordinates are needed.\n%s",
              BOX_FAILED);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
          if (j == 6) {
            outgap = -1.0;
          }
        }
        if (temp_grid.step[0] < ALMOST_ZERO) {
          tee_error(od, run_type, overall_line_num,
            "The step value must "
            "be greater than 0.0.\n%s",
            BOX_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
      }
      if (!(run_type & DRY_RUN)) {
        ++command;
        tee_printf(od, M_TOOL_INVOKE, nesting, command, "BOX", line_orig);
        tee_flush(od);
        result = create_box(od, &temp_grid, outgap, from_file);
        switch (result) {
          case PREMATURE_DAT_EOF:
          tee_error(od, run_type, overall_line_num,
            E_FILE_CORRUPTED_OR_IN_WRONG_FORMAT, "CoMFA region",
            od->file[ASCII_IN]->name, BOX_FAILED);
          return PARSE_INPUT_ERROR;
        }
        print_grid_coordinates(od, &(od->grid));
        if (od->file[ASCII_IN]->handle) {
          fclose(od->file[ASCII_IN]->handle);
          od->file[ASCII_IN]->handle = NULL;
        }
        update_pymol(od);
        update_jmol(od);
        tee_printf(od, M_TOOL_SUCCESS, nesting, command, "BOX");
        tee_flush(od);
      }
      else {
        /*
        fake data for DRY_RUN
        */
        od->grid.nodes[0] = 10;
        od->grid.nodes[1] = 10;
        od->grid.nodes[2] = 10;
        od->grid.object_num = 10;
        od->grid.struct_num = 10;
        od->object_num = 10;
      }
    }
    else if (!strcasecmp(arg->me[0], "rototrans")) {
      if (!(od->valid & SDF_BIT)) {
        tee_error(od, run_type, overall_line_num,
          E_IMPORT_MOLFILE_FIRST, ROTOTRANS_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      if (!(od->field.babel_exe_path[0])) {
        tee_error(od, run_type, overall_line_num,
          E_OPENBABEL_PATH, ROTOTRANS_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      if (od->field_num) {
        tee_error(od, run_type, overall_line_num,
          E_CANNOT_ALTER_OBJECT_COORDINATES,
          ROTOTRANS_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      if (!(parameter = get_args(od, "file"))) {
        tee_error(od, run_type, overall_line_num, "Please specify a file "
          "where new coordinates should be saved.\n%s",
          ROTOTRANS_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      strcpy(file_basename, parameter);
      gettimeofday(&start, NULL);
      result = parse_synonym_lists(od, "ROTOTRANS", ROTOTRANS_FAILED,
        (1 << OBJECT_LIST) | (1 << ID_LIST) | (1 << STRUCT_LIST), &list_type,
        OBJECT_LIST, run_type, overall_line_num);
      switch (result) {
        case PARSE_INPUT_RECOVERABLE_ERROR:
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
        
        case PARSE_INPUT_ERROR:
        return PARSE_INPUT_ERROR;
      }
      memset(trans, 0, 3 * sizeof(double));
      memset(rot, 0, 3 * sizeof(double));
      for (i = 0; i < 3; ++i) {
        sprintf(buffer, "%c_trans", i + 'x');
        if ((parameter = get_args(od, buffer))) {
          sscanf(parameter, "%lf", &trans[i]);
        }
        sprintf(buffer, "%c_rot", i + 'x');
        if ((parameter = get_args(od, buffer))) {
          sscanf(parameter, "%lf", &rot[i]);
        }
      }
      for (i = 0; i < 3; ++i) {
        rot[i] = angle2rad(rot[i]);
      }
      if (!(run_type & DRY_RUN)) {
        ++command;
        tee_printf(od, M_TOOL_INVOKE, nesting, command, "ROTOTRANS", line_orig);
        tee_flush(od);
        result = call_obenergy(od, O3_MMFF94);
        switch (result) {
          case FL_CANNOT_CREATE_CHANNELS:
          tee_error(od, run_type, overall_line_num,
            E_CANNOT_CREATE_PIPE, ROTOTRANS_FAILED);
          return PARSE_INPUT_ERROR;

          case FL_CANNOT_CHDIR:
          tee_error(od, run_type, overall_line_num,
            E_CANNOT_CHANGE_DIR, od->field.babel_exe_path, ROTOTRANS_FAILED);
          return PARSE_INPUT_ERROR;

          case FL_CANNOT_CREATE_PROCESS:
          tee_error(od, run_type, overall_line_num,
            E_CANNOT_CREATE_PROCESS, "OpenBabel", ROTOTRANS_FAILED);
          return PARSE_INPUT_ERROR;

          case OUT_OF_MEMORY:
          tee_error(od, run_type, overall_line_num,
            E_OUT_OF_MEMORY, ROTOTRANS_FAILED);
          return PARSE_INPUT_ERROR;

          case OPENBABEL_ERROR:
          tee_error(od, run_type, overall_line_num,
            E_PROGRAM_ERROR, "OpenBabel");
          if ((od->file[TEMP_LOG]->handle = fopen
            (od->file[TEMP_LOG]->name, "rb"))) {
            while (fgets(buffer, BUF_LEN,
              od->file[TEMP_LOG]->handle)) {
              buffer[BUF_LEN - 1] = '\0';
              tee_printf(od, "%s", buffer);
            }
            fclose(od->file[TEMP_LOG]->handle);
            od->file[TEMP_LOG]->handle = NULL;
            tee_printf(od, "\n%s", ROTOTRANS_FAILED);
          }
          else {
            tee_error(od, run_type, overall_line_num,
              E_CANNOT_READ_PROGRAM_LOG, "OpenBabel", ROTOTRANS_FAILED);
          }
          return PARSE_INPUT_ERROR;
        }
        result = rototrans(od, file_basename, trans, rot);
        switch (result) {
          case OUT_OF_MEMORY:
          tee_error(od, run_type, overall_line_num,
            E_OUT_OF_MEMORY, ROTOTRANS_FAILED);
          O3_ERROR_PRINT(&(od->task));
          return PARSE_INPUT_ERROR;

          case CANNOT_WRITE_ROTOTRANSED_SDF:
          tee_error(od, run_type, overall_line_num,
            E_ERROR_IN_WRITING_SDF_FILE,
            od->task.string, ROTOTRANS_FAILED);
          O3_ERROR_PRINT(&(od->task));
          return PARSE_INPUT_ERROR;

          case CANNOT_READ_ORIGINAL_SDF:
          tee_error(od, run_type, overall_line_num,
            E_ERROR_IN_READING_SDF_FILE,
            od->task.string, ROTOTRANS_FAILED);
          O3_ERROR_PRINT(&(od->task));
          return PARSE_INPUT_ERROR;
          
          case FL_CANNOT_READ_MOL_FILE:
          tee_error(od, run_type, overall_line_num,
            E_ERROR_IN_READING_MOL_FILE,
            od->task.string, ROTOTRANS_FAILED);
          O3_ERROR_PRINT(&(od->task));
          return PARSE_INPUT_ERROR;
          
          case FL_CANNOT_READ_OB_OUTPUT:
          tee_error(od, run_type, overall_line_num,
            E_ERROR_IN_READING_OB_OUTPUT,
            od->task.string, ROTOTRANS_FAILED);
          O3_ERROR_PRINT(&(od->task));
          return PARSE_INPUT_ERROR;
          
          case FL_UNKNOWN_ATOM_TYPE:
          tee_error(od, run_type, overall_line_num,
            E_UNKNOWN_ATOM_TYPE, "atom", ROTOTRANS_FAILED);
          O3_ERROR_PRINT(&(od->task));
          return PARSE_INPUT_ERROR;
        }
        tee_printf(od, M_TOOL_SUCCESS, nesting, command, "ROTOTRANS");
        tee_flush(od);
      }
    }
    else if (!strcasecmp(arg->me[0], "import")) {
      /*
      user wants to IMPORT something; check
      if the file type was specified
      */
      multi_file_type = 0;
      gettimeofday(&start, NULL);
      if (!(parameter = get_args(od, "type"))) {
        tee_error(od, run_type, overall_line_num,
          E_SPECIFY_FILE_TYPE, IMPORT_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      if ((!strcasecmp(parameter, "sdf"))
        || (!strcasecmp(parameter, "mol2"))) {
        multi_file_type = tolower(parameter[0]);
        if (!(od->field.babel_exe_path[0])) {
          tee_error(od, run_type, overall_line_num,
            E_OPENBABEL_PATH, IMPORT_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        if (!(parameter = get_args(od, "file"))) {
          tee_error(od, run_type, overall_line_num,
            E_SPECIFY_STRUCTURE_FILE, "import", IMPORT_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        if (!(run_type & DRY_RUN)) {
          memset(od->file[MOLFILE_IN]->name, 0, BUF_LEN);
          strcpy(od->file[MOLFILE_IN]->name, parameter);
          absolute_path(od->file[MOLFILE_IN]->name);
          if (!fexist(od->file[MOLFILE_IN]->name)) {
            tee_error(od, run_type, overall_line_num,
              E_FILE_CANNOT_BE_OPENED_FOR_READING,
              parameter, IMPORT_FAILED);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
          if (multi_file_type == 'm') {
            result = open_temp_dir(od, NULL, "sdf_dir", buffer);
            if (result) {
              tee_error(od, run_type, overall_line_num,
                E_TEMP_DIR_CANNOT_BE_CREATED, buffer,
                IMPORT_FAILED);
                fail = !(run_type & INTERACTIVE_RUN);
                continue;
            }
            sprintf(file_basename, "%s%cmolfile.sdf", buffer, SEPARATOR);
            if (convert_mol(od, od->file[MOLFILE_IN]->name,
              file_basename, "mol2", "sdf", "-xl")) {
              tee_error(od, run_type, overall_line_num,
                E_FILE_CORRUPTED_OR_IN_WRONG_FORMAT,
                "MOL2", od->file[MOLFILE_IN]->name, IMPORT_FAILED);
              fail = !(run_type & INTERACTIVE_RUN);
              continue;
            }
            strcpy(od->file[MOLFILE_IN]->name, file_basename);
          }
          if (od->valid & SDF_BIT) {
            for (i = 0; i < od->grid.object_num; ++i) {
              set_object_attr(od, i, OPERATE_BIT, 1);
            }
            remove_object(od);
            if (od->pel.pymol_old_object_id) {
              int_perm_free(od->pel.pymol_old_object_id);
              od->pel.pymol_old_object_id = NULL;
            }
            if (od->pel.jmol_old_object_id) {
              int_perm_free(od->pel.jmol_old_object_id);
              od->pel.jmol_old_object_id = NULL;
            }
          }
          import_y_vars = 0;
          memset(name_list, 0, BUF_LEN);
          if ((multi_file_type == 's')
            && (parameter = get_args(od, "y_var_name"))) {
            strcpy(name_list, parameter);
            import_y_vars |= IMPORT_Y_VARS_BIT;
          }
          ++command;
          tee_printf(od, M_TOOL_INVOKE, nesting, command,
            (multi_file_type == 's') ? "IMPORT SDF" : "IMPORT MOL2", line_orig);
          tee_flush(od);
          if (!(od->file[MOLFILE_IN]->handle =
            fopen(od->file[MOLFILE_IN]->name, "rb"))) {
            tee_error(od, run_type, overall_line_num,
              E_FILE_CANNOT_BE_OPENED_FOR_READING,
              parameter, IMPORT_FAILED);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
          result = open_temp_dir(od, NULL, "mol_dir", od->field.mol_dir);
          if (result) {
            tee_error(od, run_type, overall_line_num,
              E_TEMP_DIR_CANNOT_BE_CREATED, od->field.mol_dir,
              IMPORT_FAILED);
              fail = !(run_type & INTERACTIVE_RUN);
              continue;
          }
          memset(&(od->grid), 0, sizeof(GridInfo));
          result = parse_sdf(od, import_y_vars | ALLOC_MOL_INFO_BIT, name_list);
          if (od->file[MOLFILE_IN]->handle) {
            fclose(od->file[MOLFILE_IN]->handle);
            od->file[MOLFILE_IN]->handle = NULL;
          }
          gettimeofday(&end, NULL);
          elapsed_time(od, &start, &end);
          switch (result) {
            case PREMATURE_EOF:
            tee_error(od, run_type, overall_line_num, 
              E_FILE_CORRUPTED_OR_IN_WRONG_FORMAT, "SDF",
              od->file[MOLFILE_IN]->name, IMPORT_FAILED);
            return PARSE_INPUT_ERROR;

            case CANNOT_WRITE_TEMP_FILE:
            tee_error(od, run_type, overall_line_num,
              E_ERROR_IN_WRITING_TEMP_FILE,
              od->task.string, IMPORT_FAILED);
            O3_ERROR_PRINT(&(od->task));
            return PARSE_INPUT_ERROR;

            case OUT_OF_MEMORY:
            tee_error(od, run_type, overall_line_num,
              E_OUT_OF_MEMORY, IMPORT_FAILED);
            return PARSE_INPUT_ERROR;

            case CANNOT_FIND_Y_VAR_NAME:
            tee_error(od, run_type, overall_line_num, E_CANNOT_FIND_Y_VAR_NAME,
              od->file[MOLFILE_IN]->name, IMPORT_FAILED);
            return PARSE_INPUT_ERROR;

            case Y_VAR_LOW_SD:
            tee_error(od, run_type, overall_line_num,
              E_Y_VAR_LOW_SD, IMPORT_FAILED);
            return PARSE_INPUT_ERROR;
          }
          if (alloc_object_attr(od)) {
            tee_error(od, run_type, overall_line_num,
              E_OUT_OF_MEMORY, IMPORT_FAILED);
            return PARSE_INPUT_ERROR;
          }
          od->valid |= SDF_BIT;
          if (update_mol(od)) {
            tee_error(od, run_type, overall_line_num,
              E_ERROR_IN_READING_TEMP_FILE,
              od->task.string, IMPORT_FAILED);
            O3_ERROR_PRINT(&(od->task));
            return PARSE_INPUT_ERROR;
          }
          update_field_object_attr(od, VERBOSE_BIT);
          update_pymol(od);
          update_jmol(od);
          tee_printf(od, "\n"
            "This set of molecules fits in a grid box "
            "whose bottom left, top right cooordinates "
            "(no outgap) are at least:\n"
            "[(%.4f,%.4f,%.4f), (%.4f,%.4f,%.4f)]\n\n",
            od->min_coord[0], od->min_coord[1],
            od->min_coord[2], od->max_coord[0],
            od->max_coord[1], od->max_coord[2]);
          strcpy(od->default_folder, od->file[MOLFILE_IN]->name);
          get_dirname(od->default_folder);
          tee_printf(od, M_TOOL_SUCCESS, nesting, command,
            (multi_file_type == 's') ? "IMPORT SDF" : "IMPORT MOL2");
          tee_flush(od);
          multi_file_type = 0;
        }
        else {
          multi_file_type = 0;
          od->valid |= SDF_BIT;
          od->grid.object_num = 10;
          od->grid.struct_num = 10;
        }
      }
      else if (!strcasecmp(parameter, "gridkont")) {
        if (!(od->valid & SDF_BIT)) {
          tee_error(od, run_type, overall_line_num,
            E_IMPORT_MOLFILE_FIRST, IMPORT_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        /*
        user wants to IMPORT a GRIDKONT file
        */
        strcpy(od->file[BINARY_IN]->name, "grid.kont");
        if ((parameter = get_args(od, "file"))) {
          strcpy(od->file[BINARY_IN]->name, parameter);
        }
        replace_object_name = 0;
        if ((parameter = get_args(od, "replace_object_name"))) {
          if (!strncasecmp(parameter, "y", 1)) {
            replace_object_name = 1;
          }
        }
        /*
        we need to open the GRID.KONT file
        */
        if (!(od->file[BINARY_IN]->handle =
          fopen(od->file[BINARY_IN]->name, "rb"))) {
          tee_error(od, run_type, overall_line_num,
            E_FILE_CANNOT_BE_OPENED_FOR_READING,
            od->file[BINARY_IN]->name, IMPORT_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        if (!(run_type & DRY_RUN)) {
          /*
          read the GRIDKONT file
          */
          ++command;
          tee_printf(od, M_TOOL_INVOKE, nesting, command,
            "IMPORT GRIDKONT", line_orig);
          tee_flush(od);
          result = import_gridkont(od, replace_object_name);
          gettimeofday(&end, NULL);
          elapsed_time(od, &start, &end);
          /*
          check for errors
          */
          switch (result) {
            case PREMATURE_DAT_EOF:
            tee_error(od, run_type, overall_line_num,
              E_FILE_CORRUPTED_OR_IN_WRONG_FORMAT, ".kont",
              od->file[BINARY_IN]->name,
              IMPORT_FAILED);
            return PARSE_INPUT_ERROR;

            case CANNOT_WRITE_TEMP_FILE:
            tee_error(od, run_type, overall_line_num,
              E_ERROR_IN_WRITING_TEMP_FILE, "TEMP_FIELD", IMPORT_FAILED);
            return PARSE_INPUT_ERROR;

            case OUT_OF_MEMORY:
            tee_error(od, run_type, overall_line_num,
              E_OUT_OF_MEMORY, IMPORT_FAILED);
            return PARSE_INPUT_ERROR;

            case OBJECTS_NOT_MATCHING:
            tee_error(od, run_type, overall_line_num,
              E_GRIDKONT_NOT_MATCHING,
              od->file[BINARY_IN]->name, od->newgrid.object_num,
              od->object_num, IMPORT_FAILED);
            return PARSE_INPUT_ERROR;

            case GRID_NOT_MATCHING:
            tee_error(od, run_type, overall_line_num,
              E_GRID_NOT_MATCHING, "GRIDKONT",
              od->file[BINARY_IN]->name, "");
            print_grid_comparison(od);
            tee_error(od, run_type, overall_line_num, "%s\n",
              IMPORT_FAILED);
            return PARSE_INPUT_ERROR;

            default:
            update_pymol(od);
            update_jmol(od);
            tee_printf(od, M_TOOL_SUCCESS, nesting, command, "IMPORT GRIDKONT");
            tee_flush(od);
          }
          if (od->file[BINARY_IN]->handle) {
            fclose(od->file[BINARY_IN]->handle);
            od->file[BINARY_IN]->handle = NULL;
          }
        }
        else {
          /*
          fake data for DRY_RUN
          */
          od->grid.nodes[0] = 10;
          od->grid.nodes[1] = 10;
          od->grid.nodes[2] = 10;
          od->grid.object_num = 10;
          od->grid.struct_num = 10;
          od->object_num = 10;
          od->field_num = 2;
          od->valid &= (~PLS_BIT);
        }
      }
      else if (!strcasecmp(parameter, "free_format")) {
        strcpy(name_list, "N_FIELD,Z_COORD,Y_COORD,X_COORD,N_OBJECT");
        if ((parameter = get_args(od, "data_order"))) {
          strcpy(name_list, parameter);
        }
        skip_header = 0;
        if ((parameter = get_args(od, "skip_header"))) {
          sscanf(parameter, "%d", &skip_header);
        }
        if (skip_header < 0) {
          tee_error(od, run_type, overall_line_num, E_POSITIVE_NUMBER,
            "skip_header parameter", IMPORT_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        if (!(od->valid & SDF_BIT)) {
          tee_error(od, run_type, overall_line_num,
            E_IMPORT_MOLFILE_FIRST, IMPORT_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        if (!(od->grid.nodes[0])) {
          tee_error(od, run_type, overall_line_num,
            E_GRID_BOX_FIRST, IMPORT_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        result = find_vary_speed(od, name_list, max_vary, vary, &i, &j, &temp_varcoord);
        switch (result) {
          case WRONG_PARAMETER_NAME:
          tee_error(od, run_type, overall_line_num,
            "The only valid keywords for the data_order "
            "parameter are \"N_FIELD\", \"N_OBJECT\", "
            "\"X_COORD\", \"Y_COORD\" and \"Z_COORD\".\n%s",
            IMPORT_FAILED);
          break;

          case DUPLICATE_PARAMETER_NAME:
          tee_error(od, run_type, overall_line_num,
            "A duplicate keyword was detected in the "
            "data_order parameter.\n%s", IMPORT_FAILED);
          break;

          case NOT_ENOUGH_PARAMETERS:
          tee_error(od, run_type, overall_line_num,
            "Please enter all the five \"N_FIELD\", "
            "\"N_OBJECT\", \"X_COORD\", \"Y_COORD\" and "
            "\"Z_COORD\" keywords in the data_order "
            "parameter.\n%s", IMPORT_FAILED);
          break;
        }
        if (result) {
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        if (!(parameter = get_args(od, "file"))) {
          tee_error(od, run_type, overall_line_num,
            "Please specify a file from which you wish "
            "to import free format data.\n%s",
            IMPORT_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        strcpy(od->file[ASCII_IN]->name, parameter);
        absolute_path(od->file[ASCII_IN]->name);
        if (!fexist(od->file[ASCII_IN]->name)) {
          tee_error(od, run_type, overall_line_num,
            E_FILE_CANNOT_BE_OPENED_FOR_READING,
            od->file[ASCII_IN]->name, IMPORT_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        if (!(run_type & DRY_RUN)) {
          if (!(od->file[ASCII_IN]->handle =
            (FILE *)fzopen(od->file[ASCII_IN]->name, "rb"))) {
            tee_error(od, run_type, overall_line_num,
              E_FILE_CANNOT_BE_OPENED_FOR_READING,
              od->file[ASCII_IN]->name, IMPORT_FAILED);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
          ++command;
          tee_printf(od, M_TOOL_INVOKE, nesting, command, "IMPORT FREE_FORMAT", line_orig);
          tee_flush(od);
          result = import_free_format(od, name_list, skip_header, &n_values);
          gettimeofday(&end, NULL);
          elapsed_time(od, &start, &end);
          if (od->file[ASCII_IN]->handle) {
            fzclose((fzPtr *)(od->file[ASCII_IN]->handle));
            od->file[ASCII_IN]->handle = NULL;
          }
          if (result) {
            O3_ERROR_PRINT(&(od->task));
            switch (result) {
              case OUT_OF_MEMORY:
              tee_error(od, run_type, overall_line_num,
                E_OUT_OF_MEMORY, IMPORT_FAILED);
              return PARSE_INPUT_ERROR;

              case WRONG_DATA_FORMAT:
              tee_error(od, run_type, overall_line_num,
                "The free format ASCII file %s should "
                "consist in a series of floating-point "
                "values separated by space, tab, comma, "
                "semicolon or newline.\n%s\n",
                od->file[ASCII_IN]->name, IMPORT_FAILED);
              return PARSE_INPUT_ERROR;

              case NOT_ENOUGH_VALUES:
              tee_error(od, run_type, overall_line_num,
                "The free format ASCII file %s should "
                "consist in a series of at least %d floating-point "
                "values separated by space, tab, comma, "
                "semicolon or newline; instead, "
                "only %d values were read.\n%s\n",
                od->file[ASCII_IN]->name, od->x_vars, n_values,
                IMPORT_FAILED);
              return PARSE_INPUT_ERROR;

              case INCORRECT_NUMBER_OF_VALUES:
              tee_error(od, run_type, overall_line_num,
                "The free format ASCII file %s should "
                "consist in a series of floating-point "
                "values separated by space, tab, comma, "
                "semicolon or newline multiple of %d; "
                "instead, %d values were read, which "
                "is not an exact multiple.\n%s\n",
                od->file[ASCII_IN]->name, od->x_vars, n_values,
                IMPORT_FAILED);
              return PARSE_INPUT_ERROR;

              case EOF_FOUND:
              case PREMATURE_EOF:
              tee_error(od, run_type, overall_line_num,
                E_FILE_CORRUPTED_OR_IN_WRONG_FORMAT, "ASCII",
                od->file[ASCII_IN]->name, IMPORT_FAILED);
              return PARSE_INPUT_ERROR;

              case Y_VAR_LOW_SD:
              tee_error(od, run_type, overall_line_num,
                E_Y_VAR_LOW_SD, IMPORT_FAILED);
              return PARSE_INPUT_ERROR;
            }
          }
          tee_printf(od, M_TOOL_SUCCESS, nesting, command, "IMPORT FREE_FORMAT");
          tee_flush(od);
        }
        else {
          /*
          fake data for DRY_RUN
          */
          od->field_num = 2;
          od->valid &= (~PLS_BIT);
        }
      }
      else if (!strcasecmp(parameter, "dependent")) {
        /*
        user wants to IMPORT dependent variables
        */
        if (!(parameter = get_args(od, "file"))) {
          tee_error(od, run_type, overall_line_num,
            "Please specify a file from which you wish "
            "to import dependent variables.\n%s",
            IMPORT_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        if (!(od->valid & SDF_BIT)) {
          tee_error(od, run_type, overall_line_num,
            E_IMPORT_MOLFILE_FIRST, IMPORT_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        /*
        we need a txt file
        */
        memset(od->file[DEP_IN]->name, 0, BUF_LEN);
        strcpy(od->file[DEP_IN]->name, parameter);
        absolute_path(od->file[DEP_IN]->name);
        if (!((od->file[DEP_IN]->handle) =
          fopen(od->file[DEP_IN]->name, "rb"))) {
          tee_error(od, run_type, overall_line_num,
            E_FILE_CANNOT_BE_OPENED_FOR_READING,
            parameter, IMPORT_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        /*
        import the txt file
        */
        if (!(run_type & DRY_RUN)) {
          strcpy(name_list, "all");
          if ((parameter = get_args(od, "y_var_name"))) {
            strcpy(name_list, parameter);
          }
          ++command;
          tee_printf(od, M_TOOL_INVOKE, nesting, command,
            "IMPORT DEPENDENT", line_orig);
          tee_flush(od);
          result = import_dependent(od, name_list);
          gettimeofday(&end, NULL);
          elapsed_time(od, &start, &end);
          /*
          check for errors
          */
          switch (result) {
            case CANNOT_READ_TEMP_FILE:
            tee_error(od, run_type, overall_line_num,
              E_ERROR_IN_READING_TEMP_FILE, "TEMP_FIELD", IMPORT_FAILED);
            return PARSE_INPUT_ERROR;

            case NOT_ENOUGH_OBJECTS:
            tee_error(od, run_type, overall_line_num,
              "In the file \"%s\" there must be exactly "
              "%d lines, the first one with variable "
              "names followed by dependent variable values.\n%s",
              od->file[DEP_IN]->name, od->grid.object_num + 1,
              IMPORT_FAILED);
            return PARSE_INPUT_ERROR;

            case WRONG_NUMBER_OF_Y_VARS:
            tee_error(od, run_type, overall_line_num, E_WRONG_NUMBER_OF_Y_VARS,
              od->file[DEP_IN]->name, IMPORT_FAILED);
            return PARSE_INPUT_ERROR;

            case PREMATURE_DEP_IN_EOF:
            tee_error(od, run_type, overall_line_num,
              E_FILE_CORRUPTED_OR_IN_WRONG_FORMAT, "ASCII",
              od->file[DEP_IN]->name, IMPORT_FAILED);
            return PARSE_INPUT_ERROR;

            case CANNOT_FIND_Y_VAR_NAME:
            tee_error(od, run_type, overall_line_num, E_CANNOT_FIND_Y_VAR_NAME,
              od->file[DEP_IN]->name, IMPORT_FAILED);
            return PARSE_INPUT_ERROR;

            case Y_VAR_LOW_SD:
            tee_error(od, run_type, overall_line_num,
              E_Y_VAR_LOW_SD, IMPORT_FAILED);
            return PARSE_INPUT_ERROR;

            default:
            tee_printf(od, M_TOOL_SUCCESS, nesting, command, "IMPORT DEPENDENT");
            tee_flush(od);
          }
        }
      }
      else if ((!strcasecmp(parameter, "gamess_cube"))
        || (!strcasecmp(parameter, "formatted_cube"))) {
        multi_file_type = FORMATTED_CUBE_INPUT_FILE;
        multi_file_fd = od->file[ASCII_IN];
      }
      else if ((!strcasecmp(parameter, "gaussian_cube"))
        || (!strcasecmp(parameter, "unformatted_cube"))) {
        multi_file_type = UNFORMATTED_CUBE_INPUT_FILE;
        multi_file_fd = od->file[BINARY_IN];
      }
      else if (!strcasecmp(parameter, "molden")) {
        multi_file_type = MOLDEN_INPUT_FILE;
        multi_file_fd = od->file[BINARY_IN];
      }
      else if (!strcasecmp(parameter, "moe_grid")) {
        multi_file_type = MOE_GRID_INPUT_FILE;
        multi_file_fd = od->file[BINARY_IN];
      }
      else if (!strcasecmp(parameter, "grid_ascii")) {
        multi_file_type = GRID_ASCII_INPUT_FILE;
        multi_file_fd = NULL;
      }
      else {
        tee_error(od, run_type, overall_line_num,
          "Only \"SDF\", \"MOL2\", \"FREE_FORMAT\", \"GRIDKONT\", "
          "\"DEPENDENT\", \"GAMESS_CUBE\", \"FORMATTED_CUBE\", "
          "\"GAUSSIAN_CUBE\", \"UNFORMATTED_CUBE\", "
          "\"MOLDEN\", \"MOE_GRID\" and \"GRID_ASCII\" "
          "types are allowed for the IMPORT keyword.\n%s",
          IMPORT_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      if (multi_file_type) {
        /*
        user wants to IMPORT grid data
        */
        memset(regex_name[0], 0, BUF_LEN);
        if ((parameter = get_args(od, "file"))) {
          strncpy(regex_name[0], parameter, BUF_LEN - 1);
          regex_name[0][BUF_LEN - 1] = '\0';
          #ifdef WIN32
          slash_to_backslash(regex_name[0]);
          #endif
          strcpy(buffer, regex_name[0]);
          get_dirname(buffer);
          if (!dexist(buffer)) {
            tee_error(od, run_type, overall_line_num,
              E_DIR_NOT_EXISTING, buffer,
              IMPORT_FAILED);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
        }
        else {
          tee_error(od, run_type, overall_line_num,
            "Please specify a regex name identifying the files "
            "from which you wish to import %s data.\n%s",
            multi_file_name[multi_file_type],
            IMPORT_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        mo = 0;
        /*
        by default, if it is a MO cube file, all MOs all imported
        */
        if ((multi_file_type == UNFORMATTED_CUBE_INPUT_FILE)
          || (multi_file_type == FORMATTED_CUBE_INPUT_FILE)) {
          if ((parameter = get_args(od, "mo"))) {
            if (strncasecmp(parameter, "all", 3)) {
              sscanf(parameter, "%d", &mo);
              if (mo < 1) {
                tee_error(od, run_type, overall_line_num,
                  "The mo parameter should be an integer "
                  "greater than or equal to 1.\n%s",
                  IMPORT_FAILED);
                fail = !(run_type & INTERACTIVE_RUN);
                continue;
              }
            }
          }
        }
        if (!(od->valid & SDF_BIT)) {
          tee_error(od, run_type, overall_line_num,
            E_IMPORT_MOLFILE_FIRST, IMPORT_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        if (!(od->grid.nodes[0])) {
          tee_error(od, run_type, overall_line_num,
            E_GRID_BOX_FIRST, IMPORT_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        if ((multi_file_type == GRID_ASCII_INPUT_FILE)
          || (multi_file_type == MOE_GRID_INPUT_FILE)) {
          if (check_regex_name(regex_name[0], 1)) {
            tee_error(od, run_type, overall_line_num,
              E_ONLY_ONE_WILDCARD, IMPORT_FAILED);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
        }
        else if (!(run_type & DRY_RUN)) {
          result = open_temp_file(od, od->file[TEMP_SORTED_MATCH], "sorted_match");
          if (result) {
            tee_error(od, run_type, overall_line_num,
              E_TEMP_FILE_CANNOT_BE_OPENED_FOR_WRITING,
              od->file[TEMP_SORTED_MATCH]->name, IMPORT_FAILED);
            return PARSE_INPUT_ERROR;
          }
          result = open_temp_file(od, od->file[TEMP_OBJECT_MATCH], "object_match");
          if (result) {
            tee_error(od, run_type, overall_line_num,
              E_TEMP_FILE_CANNOT_BE_OPENED_FOR_WRITING,
              od->file[TEMP_OBJECT_MATCH]->name, IMPORT_FAILED);
            return PARSE_INPUT_ERROR;
          }
          result = match_objects_with_datafile(od, regex_name[0], multi_file_type);
          strcpy(buffer, regex_name[0]);
          get_dirname(buffer);
          switch (result) {
            case CANNOT_OPEN_DIRECTORY:
            tee_error(od, run_type, overall_line_num,
              E_DIR_NOT_EXISTING, buffer, IMPORT_FAILED);
            return PARSE_INPUT_ERROR;
            
            case NOT_ENOUGH_OBJECTS:
            case TOO_MANY_OBJECTS:
            tee_error(od, run_type, overall_line_num,
              E_NOT_ENOUGH_OBJECTS,
              od->grid.object_num, regex_name[0],
              IMPORT_FAILED);
            return PARSE_INPUT_ERROR;

            case CANNOT_READ_ORIGINAL_SDF:
            tee_error(od, run_type, overall_line_num,
              E_ERROR_IN_READING_SDF_FILE,
              od->task.string, IMPORT_FAILED);
            O3_ERROR_PRINT(&(od->task));
            return PARSE_INPUT_ERROR;
        
            case CANNOT_READ_TEMP_FILE:
            tee_error(od, run_type, overall_line_num,
              E_ERROR_IN_READING_TEMP_FILE, 
              od->task.string, IMPORT_FAILED);
            O3_ERROR_PRINT(&(od->task));
            return PARSE_INPUT_ERROR;

            case PREMATURE_DAT_EOF:
            tee_error(od, run_type, overall_line_num,
              E_FILE_CORRUPTED_OR_IN_WRONG_FORMAT,
              multi_file_name[multi_file_type],
              od->task.string, IMPORT_FAILED);
            O3_ERROR_PRINT(&(od->task));
            return PARSE_INPUT_ERROR;

            case OBJECTS_NOT_MATCHING:
            tee_error(od, run_type, overall_line_num,
              E_OBJECTS_NOT_MATCHING, regex_name[0],
              od->al.mol_info[od->task.data[DATA_OBJECT_NUM]]->object_id,
              od->al.mol_info[od->task.data[DATA_OBJECT_NUM]]->object_name,
              IMPORT_FAILED);
            O3_ERROR_PRINT(&(od->task));
            return PARSE_INPUT_ERROR;
          }
        }
        if (!(run_type & DRY_RUN)) {
          /*
          import the grid file
          */
          ++command;
          tee_printf(od, M_TOOL_INVOKE, nesting, command,
            multi_file_name[multi_file_type], line_orig);
          tee_flush(od);
          switch (multi_file_type) {
            case FORMATTED_CUBE_INPUT_FILE:
            result = import_grid_formatted_cube(od, mo);
            break;

            case UNFORMATTED_CUBE_INPUT_FILE:
            result = import_grid_unformatted_cube(od, mo);
            break;

            case MOLDEN_INPUT_FILE:
            result = import_grid_molden(od);
            break;

            case MOE_GRID_INPUT_FILE:
            result = import_grid_moe(od, regex_name[0]);
            break;

            default:
            result = import_grid_ascii(od, regex_name[0]);
            break;
          }
          gettimeofday(&end, NULL);
          elapsed_time(od, &start, &end);
          /*
          check for errors
          */
          switch (result) {
            case NOT_ENOUGH_OBJECTS:
            tee_error(od, run_type, overall_line_num,
              E_NOT_ENOUGH_OBJECTS,
              od->object_num, regex_name[0],
              IMPORT_FAILED);
            return PARSE_INPUT_ERROR;

            case PREMATURE_DAT_EOF:
            tee_error(od, run_type, overall_line_num,
              E_FILE_CORRUPTED_OR_IN_WRONG_FORMAT,
              multi_file_name[multi_file_type],
              od->file[BINARY_IN]->name, IMPORT_FAILED);
            return PARSE_INPUT_ERROR;

            case CANNOT_READ_TEMP_FILE:
            tee_error(od, run_type, overall_line_num,
              E_ERROR_IN_READING_TEMP_FILE,
              od->file[TEMP_SORTED_MATCH]->name, IMPORT_FAILED);
            return PARSE_INPUT_ERROR;

            case CANNOT_WRITE_TEMP_FILE:
            tee_error(od, run_type, overall_line_num,
              E_ERROR_IN_WRITING_TEMP_FILE, "TEMP_FIELD", IMPORT_FAILED);
            return PARSE_INPUT_ERROR;

            case OUT_OF_MEMORY:
            tee_error(od, run_type, overall_line_num,
              E_OUT_OF_MEMORY, IMPORT_FAILED);
            return PARSE_INPUT_ERROR;

            case CANNOT_READ_GRID_DATA:
            sprintf(buffer, regex_name[0], od->newgrid.object_num + 1);
            tee_error(od, run_type, overall_line_num,
              E_ERROR_IN_GRID_DATA,
              od->newgrid.x_vars, buffer,
              IMPORT_FAILED);
            return PARSE_INPUT_ERROR;

            case GRID_NOT_MATCHING:
            sprintf(buffer, regex_name[0], od->newgrid.object_num + 1);
            tee_error(od, run_type, overall_line_num,
              E_GRID_NOT_MATCHING,
              multi_file_name[multi_file_type],
              multi_file_fd ? multi_file_fd->name : buffer,
              IMPORT_FAILED);
            print_grid_comparison(od);
            return PARSE_INPUT_ERROR;

            case GRID_NOT_MATCHING_OFF_CENTER:
            sprintf(buffer, regex_name[0], od->newgrid.object_num + 1);
            tee_error(od, run_type, overall_line_num,
              E_GRID_NOT_MATCHING_DATA_POINT,
              od->newgrid.x_vars,
              od->newgrid.start_coord[0],
              od->newgrid.start_coord[1],
              od->newgrid.start_coord[2],
              buffer, "is off-center",
              IMPORT_FAILED);
            return PARSE_INPUT_ERROR;

            case GRID_NOT_MATCHING_OUT_OF_BOUNDS:
            sprintf(buffer, regex_name[0], od->newgrid.object_num + 1);
            tee_error(od, run_type, overall_line_num,
              E_GRID_NOT_MATCHING_DATA_POINT,
              od->newgrid.x_vars,
              od->newgrid.start_coord[0],
              od->newgrid.start_coord[1],
              od->newgrid.start_coord[2],
              buffer, "is out of bounds",
              IMPORT_FAILED);
            return PARSE_INPUT_ERROR;

            case CANNOT_FIND_MO:
            sprintf(buffer, regex_name[0], od->newgrid.object_num + 1);
            tee_error(od, run_type, overall_line_num,
              E_CANNOT_FIND_MO, buffer, IMPORT_FAILED);
            return PARSE_INPUT_ERROR;

            default:
            update_pymol(od);
            update_jmol(od);
            sprintf(buffer, "IMPORT %s", multi_file_name[multi_file_type]);
            tee_printf(od, M_TOOL_SUCCESS, nesting, command, buffer);
            tee_flush(od);
          }
        }
        else {
          /*
          fake data for DRY_RUN
          */
          od->grid.nodes[0] = 10;
          od->grid.nodes[1] = 10;
          od->grid.nodes[2] = 10;
          od->grid.object_num = 10;
          od->grid.struct_num = 10;
          od->object_num = 10;
          od->field_num = 2;
          od->valid &= (~PLS_BIT);
        }
      }
    }
    else if ((!strcasecmp(arg->me[0], "calc_field"))
      || (!strcasecmp(arg->me[0], "prepare"))) {
      prep_or_calc = tolower((int)(arg->me[0][0]));
      failed = ((prep_or_calc == 'c') ? CALC_FIELD_FAILED : PREPARE_FAILED);
      gettimeofday(&start, NULL);
      /*
      user wants to compute a MIF out of a SDF file;
      check if it has already been loaded
      */

      if (!(od->valid & SDF_BIT)) {
        tee_error(od, run_type, overall_line_num,
          E_IMPORT_MOLFILE_FIRST, failed);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      if (!(od->field.babel_exe_path[0])) {
        tee_error(od, run_type, overall_line_num,
          E_OPENBABEL_PATH, failed);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      if (!(od->grid.nodes[0])) {
        tee_error(od, run_type, overall_line_num,
          E_GRID_BOX_FIRST, failed);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      /*
      remove the extension from the molfile name
      */
      memset(file_basename, 0, BUF_LEN);
      strcpy(file_basename, get_basename(od->file[MOLFILE_IN]->name));
      remove_extension(file_basename);
      od->field.type = 0;
      if (!(parameter = get_args(od, "type"))) {
        tee_error(od, run_type, overall_line_num,  
          "Please specify whether a VDW, MM_ELE, CS3D, QM_ELE, "
          "QM_DEN or MD_GRID field should be computed.\n%s",
          failed);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      if (prep_or_calc == 'c') {
        if (strncasecmp(parameter, "vdw", 3)
          && strncasecmp(parameter, "mm_ele", 6)
          && strncasecmp(parameter, "cs3d", 4)
          && strncasecmp(parameter, "qm_ele", 6)
          && strncasecmp(parameter, "qm_den", 6)
          && strncasecmp(parameter, "md_grid", 7)) {
          tee_error(od, run_type, overall_line_num,
            "The allowed field types are VDW, MM_ELE, "
            "CS3D, QM_ELE, QM_DEN and MD_GRID.\n%s",
            failed);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
      }
      else {
        if (strncasecmp(parameter, "cs3d", 4)
          && strncasecmp(parameter, "qm_ele", 5)
          && strncasecmp(parameter, "qm_den", 5)
          && strncasecmp(parameter, "moe_grid", 8)
          && strncasecmp(parameter, "sybyl", 5)) {
          tee_error(od, run_type, overall_line_num,
            "The allowed field types are "
            "CS3D, MOE_GRID, SYBYL, QM_ELE and QM_DEN.\n%s",
            PREPARE_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
      }
      if (!strncasecmp(parameter, "vdw", 3)) {
        od->field.type |= VDW_FIELD;
      }
      else if (!strncasecmp(parameter, "mm_ele", 5)) {
        od->field.type |= MM_ELE_FIELD;
        od->field.diel_const = 1.0;
        if ((parameter = get_args(od, "diel_const"))) {
          sscanf(parameter, "%lf", &(od->field.diel_const));
          if (od->field.diel_const <= 0.0) {
            tee_error(od, run_type, overall_line_num,
              E_POSITIVE_NUMBER, "dielectric constant", failed);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
        }
        od->field.diel_dep = CONST_DIELECTRIC;
        if ((parameter = get_args(od, "diel_dep"))) {
          if (!strncasecmp(parameter, "const", 5)) {
            od->field.diel_dep = CONST_DIELECTRIC;
          }
          else if (!strncasecmp(parameter, "dist", 4)) {
            od->field.diel_dep = DIST_DEP_DIELECTRIC;
          }
          else {
            tee_error(od, run_type, overall_line_num,
              "The only allowed keywords for the diel_dep "
              "parameter are CONSTANT and DISTANCE.\n%s",
              failed);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
        }
      }
      else if (!strncasecmp(parameter, "qm_ele", 5)) {
        od->field.type |= QM_ELE_FIELD;
      }
      else if (!strncasecmp(parameter, "qm_den", 5)) {
        od->field.type |= QM_DEN_FIELD;
      }
      else if (!strncasecmp(parameter, "cs3d", 4)) {
        od->field.type |= CS3D_FIELD;
        if (!(od->field.cs3d_exe[0])) {
          tee_error(od, run_type, overall_line_num,
            E_CS3D_EXE, CALC_FIELD_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        memset(regex_name[0], 0, BUF_LEN);
        od->field.match = 0;
        if ((prep_or_calc == 'c') && (parameter = get_args(od, "file"))) {
          strncpy(regex_name[0], parameter, BUF_LEN - 1);
          regex_name[0][BUF_LEN - 1] = '\0';
          #ifdef WIN32
          slash_to_backslash(regex_name[0]);
          #endif
          strcpy(buffer, regex_name[0]);
          get_dirname(buffer);
          if (!dexist(buffer)) {
            tee_error(od, run_type, overall_line_num,
              E_DIR_NOT_EXISTING, buffer,
              CALC_FIELD_FAILED);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
          od->field.match = 1;
        }
        od->field.basis_set = O3_SVP;
        if ((parameter = get_args(od, "basis_set"))) {
          if (!strcasecmp(parameter, "TZVP")) {
            od->field.basis_set = O3_TZVP;
          }
          else if (strcasecmp(parameter, "SVP")) {
            tee_error(od, run_type, overall_line_num,
              "The allowed basis sets are SVP and TZVP.\n%s",
              CALC_FIELD_FAILED);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
        }
        od->field.idelsig = CS3D_DEFAULT_IDELSIG;
        if ((parameter = get_args(od, "delsig"))) {
          sscanf(parameter, "%d", &(od->field.idelsig));
          if (od->field.idelsig <= 0) {
            tee_error(od, run_type, overall_line_num,
              E_POSITIVE_NUMBER, "delsig", CALC_FIELD_FAILED);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
        }
        od->field.compress = O3_COMPRESS_GZIP;
        if ((parameter = get_args(od, "compress"))) {
          switch (tolower(parameter[0])) {
            case 'g':
            od->field.compress = O3_COMPRESS_GZIP;
            break;
            
            #ifdef HAVE_LIBMINIZIP
            case 'z':
            od->field.compress = O3_COMPRESS_ZIP;
            break;
            #endif
            
            default:
            od->field.compress = 0;
            break;
          }
        }
      }
      else if (!strncasecmp(parameter, "md_grid", 7)) {
        od->field.type |= MD_GRID_FIELD;
        od->field.force_field = O3_MD_GRID;
        if (!(od->field.md_grid_exe_path[0])) {
          tee_error(od, run_type, overall_line_num,
            E_MD_GRID_PATH, IMPORT_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        od->field.diel_const = 80.0;
        if ((parameter = get_args(od, "diel_const"))) {
          sscanf(parameter, "%lf", &(od->field.diel_const));
          if (od->field.diel_const <= 0.0) {
            tee_error(od, run_type, overall_line_num,
              E_POSITIVE_NUMBER, "dielectric constant", failed);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
        }
        od->field.md_grid_cutoff = 5.0;
        if ((parameter = get_args(od, "cutoff"))) {
          sscanf(parameter, "%lf", &(od->field.md_grid_cutoff));
          if (od->field.md_grid_cutoff <= 0.0) {
            tee_error(od, run_type, overall_line_num,
              E_POSITIVE_NUMBER, "max energy cutoff", failed);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
        }
      }
      else if (!strncasecmp(parameter, "sybyl", 5)) {
        od->field.type |= PREP_SYBYL_INPUT;
        if (!(parameter = get_args(od, "file"))) {
          tee_error(od, run_type, overall_line_num,
            "Please indicate a file where region "
            "information should be written.\n%s",
            failed);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        len = strlen(parameter);
        result = (len < 5);
        if (!result) {
          result = strncmp(&parameter[len - 4], ".rgn", 4);
        }
        if (!result) {
          for (i = 0; (i < (len - 4)) && (!result); ++i) {
            result = !isalnum((int)parameter[i]);
          }
        }
        if (result) {
          tee_error(od, run_type, overall_line_num,
            "Please use only alphanumeric characters "
            "and the .rgn extension in your filename.\n%s",
            failed);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        strcpy(od->file[PREPINP_OUT]->name, parameter);
      }
      else if (!strncasecmp(parameter, "moe_grid", 8)) {
        od->field.type |= PREP_MOE_GRID_INPUT;
        if (!(parameter = get_args(od, "file"))) {
          tee_error(od, run_type, overall_line_num,
            "Please indicate a file where MOE grid "
            "information should be written.\n%s",
            failed);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        len = strlen(parameter);
        result = (len < 5);
        if (!result) {
          result = strncmp(&parameter[len - 4], ".svl", 4);
        }
        if (!result) {
          for (i = 0; (i < (len - 4)) && (!result); ++i) {
            result = !isalnum((int)parameter[i]);
          }
        }
        if (result) {
          tee_error(od, run_type, overall_line_num,
            "Please use only alphanumeric characters "
            "and the .svl extension in your filename.\n%s",
            failed);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        strcpy(od->file[PREPINP_OUT]->name, parameter);
      }
      od->field.spin = 0;
      if (od->field.type & (VDW_FIELD | MM_ELE_FIELD | MD_GRID_FIELD)) {
        strcpy(od->field.probe.atom_name, default_probe[(int)(od->field.force_field)]);
        found = 0;
        if ((parameter = get_args(od, "probe_type"))) {
          strncpy(od->field.probe.atom_name, parameter, MAX_FF_TYPE_LEN - 1);
          od->field.probe.atom_name[MAX_FF_TYPE_LEN - 1] = '\0';
        }
        i = 0;
        while (ff_parm[(int)(od->field.force_field)][i].type_num
          && (!(found = (!strcasecmp(od->field.probe.atom_name,
          ff_parm[(int)(od->field.force_field)][i].type_chr))))) {
          ++i;
        }
        if (!found) {
          tee_error(od, run_type, overall_line_num,
            E_UNKNOWN_ATOM_TYPE,
            "probe", failed);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
      }
      if (od->field.type & (VDW_FIELD | MM_ELE_FIELD)) {
        od->field.force_field = O3_MMFF94;
        od->field.probe.atom_type = i;
        od->field.smooth_probe_flag = 0;
        if ((parameter = get_args(od, "smooth_probe"))) {
          if (!strncasecmp(parameter, "y", 1)) {
            od->field.smooth_probe_flag = -1;
          }
        }
      }
      else if (od->field.type & (QM_ELE_FIELD | QM_DEN_FIELD | CS3D_FIELD)) {
        if (!(run_type & DRY_RUN)) {
          result = open_temp_file(od, od->file[TEMP_SORTED_MATCH], "sorted_match");
          if (result) {
            tee_error(od, run_type, overall_line_num,
              E_TEMP_FILE_CANNOT_BE_OPENED_FOR_WRITING,
              od->file[TEMP_SORTED_MATCH]->name, failed);
            return PARSE_INPUT_ERROR;
          }
        }
        if (!(od->field.type & CS3D_FIELD)) {
          if ((!strncasecmp(od->field.qm_exe, "g98", 3))
            || (!strncasecmp(od->field.qm_exe, "g03", 3))
            || (!strncasecmp(od->field.qm_exe, "g09", 3))) {
            od->field.type |= PREP_GAUSSIAN_INPUT;
            strcpy(od->field.qm_software, "gaussian");
          }
          else if (!strncasecmp(od->field.qm_exe, "gamess", 6)) {
            od->field.type |= PREP_GAMESS_INPUT;
            strcpy(od->field.qm_software, "gamess");
            #ifdef WIN32
            od->field.mpiexec_exe[0] = '\0';
            sprintf(buffer, "%s%c%s", od->field.qm_exe_path,
              SEPARATOR, GAMESS_DDIKICK_EXE);
            if (!fexist(buffer)) {
              if (!is_in_path(GAMESS_MPIEXEC_EXE, od->field.mpiexec_exe)) {
                sprintf(od->field.mpiexec_exe, "%s%c%s",
                  od->field.qm_exe_path, SEPARATOR, GAMESS_MPIEXEC_EXE);
                if (!fexist(od->field.mpiexec_exe)) {
                  tee_error(od, run_type, overall_line_num,
                    E_EXECUTABLE_PATH, GAMESS_MPIEXEC_EXE,
                    od->field.qm_exe_path, failed);
                  fail = !(run_type & INTERACTIVE_RUN);
                  continue;
                }
              }
              if (!is_in_path(GAMESS_SMPD_EXE, buffer)) {
                sprintf(buffer, "%s%c%s",
                  od->field.qm_exe_path, SEPARATOR, GAMESS_SMPD_EXE);
                if (!fexist(buffer)) {
                  tee_error(od, run_type, overall_line_num,
                    E_EXECUTABLE_PATH, GAMESS_SMPD_EXE,
                    od->field.qm_exe_path, failed);
                  fail = !(run_type & INTERACTIVE_RUN);
                  continue;
                }
              }
            }
            #endif
          }
          else if (!strncasecmp(od->field.qm_exe, "firefly", 7)) {
            od->field.type |= PREP_FIREFLY_INPUT;
            strcpy(od->field.qm_software, "firefly");
          }
          else if ((!strncasecmp(od->field.qm_exe, "dscf", 4))
            || (!strncasecmp(od->field.qm_exe, "define", 6))
            || (!strncasecmp(od->field.qm_exe, "ridft", 5))) {
            od->field.type |= PREP_TURBOMOLE_INPUT;
            strcpy(od->field.qm_software, "turbomole");
          }
          else if ((prep_or_calc == 'p')
            && (!strncasecmp(od->field.qm_exe, "molden", 6))) {
            od->field.type |= PREP_MOLDEN_INPUT;
            strcpy(od->field.qm_software, "molden");
          }
          else {
            tee_error(od, run_type, overall_line_num,
              "Unknown QM engine. Please set the "
              "O3_QM_ENGINE environment variable "
              "or use the \"env qm_engine\" keyword "
              "to indicate the QM engine you wish to use; "
              "currently GAUSSIAN, FIREFLY, GAMESS-US, "
              "TURBOMOLE and MOLDEN "
              "(in the case of MOLDEN PREPARE only) "
              "are supported.\n%s", failed);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
          if (!(od->field.type & PREP_MOLDEN_INPUT)) {
            od->field.spin = O3_RESTRICTED;
            if ((parameter = get_args(od, "spin"))) {
              if (!strncasecmp(parameter, "u", 1)) {
                od->field.spin = O3_UNRESTRICTED;
              }
              else if (strncasecmp(parameter, "r", 1)) {
                tee_error(od, run_type, overall_line_num,
                  "The allowed SPIN types are "
                  "R and U.\n%s",
                  failed);
                fail = !(run_type & INTERACTIVE_RUN);
                continue;
              }
            }
            strcpy(od->field.theory, O3_HF);
            if ((parameter = get_args(od, "theory"))) {
              if (!strncasecmp(parameter, "DFT", 3)) {
                strcpy(od->field.theory, O3_DFT_B3LYP);
              }
              else if (strncasecmp(parameter, "HF", 2)) {
                tee_error(od, run_type, overall_line_num,
                  "The allowed theory levels are "
                  "HF and DFT.\n%s",
                  failed);
                fail = !(run_type & INTERACTIVE_RUN);
                continue;
              }
            }
            od->field.basis_set = O3_6_31G;
            strcpy(buffer, "qm_");
            if ((parameter = get_args(od, "basis_set"))) {
              strcat(buffer, parameter);
              if (!strcasecmp(parameter, "STO-3G")) {
                od->field.basis_set = O3_STO_3G;
              }
              else if ((!strcasecmp(parameter, "3-21G"))
                && (!(od->field.type & PREP_TURBOMOLE_INPUT))) {
                od->field.basis_set = O3_3_21G;
              }
              else if ((!strcasecmp(parameter, "6-31G"))
                && (!(od->field.type & PREP_TURBOMOLE_INPUT))) {
                od->field.basis_set = O3_6_31G;
              }
              else if ((!strcasecmp(parameter, "6-311G"))
                && (!(od->field.type & PREP_TURBOMOLE_INPUT))) {
                od->field.basis_set = O3_6_311G;
              }
              else if ((!strcasecmp(parameter, "SV"))
                && (od->field.type & PREP_TURBOMOLE_INPUT)) {
                od->field.basis_set = O3_SV;
              }
              else if ((!strcasecmp(parameter, "SVP"))
                && (od->field.type & PREP_TURBOMOLE_INPUT)) {
                od->field.basis_set = O3_SVP;
              }
              else if ((!strcasecmp(parameter, "TZVP"))
                && (od->field.type & PREP_TURBOMOLE_INPUT)) {
                od->field.basis_set = O3_TZVP;
              }
              else if (!strcasecmp(parameter, "EMSL_3-21G")) {
                od->field.basis_set = EMSL_BASIS_SET | O3_3_21G;
              }
              else if (!strcasecmp(parameter, "EMSL_6-311G")) {
                od->field.basis_set = EMSL_BASIS_SET | O3_6_311G;
              }
              else if (!strcasecmp(parameter, "EMSL_6-311Gxx")) {
                od->field.basis_set = EMSL_BASIS_SET | O3_6_311GXX;
              }
              else {
                if (!(od->field.type & PREP_TURBOMOLE_INPUT)) {
                  tee_error(od, run_type, overall_line_num,
                    "The allowed basis sets are STO-3G, "
                    "3-21G, 6-31G, 6-311G, EMSL_3-21G, "
                    "EMSL_6-311G and EMSL_6-311Gxx.\n%s",
                    failed);
                }
                else {
                  tee_error(od, run_type, overall_line_num,
                    "The allowed basis sets are STO-3G, "
                    "SV, SVP, TZVP, EMSL_3-21G, "
                    "EMSL_6-311G and EMSL_6-311Gxx.\n%s",
                    failed);
                }
                fail = !(run_type & INTERACTIVE_RUN);
                continue;
              }
            }
            else {
              strcat(buffer, ((!(od->field.type & PREP_TURBOMOLE_INPUT)) ? "6-31G" : "SVP"));
            }
            len = strlen(buffer);
            od->field.d_func = 0;
            od->field.p_func = 0;
            od->field.f_func = 0;
            od->field.diff_sp = 0;
            if ((od->field.basis_set != O3_STO_3G) && (!(od->field.type & PREP_TURBOMOLE_INPUT))) {
              if (!((od->field.type & PREP_GAUSSIAN_INPUT)
                && (od->field.basis_set == O3_3_21G))) {
                od->field.d_func = 1;
                if ((parameter = get_args(od, "d_func"))) {
                  sscanf(parameter, "%d", &(od->field.d_func));
                  if ((od->field.d_func < 0) || (od->field.d_func > 3)) {
                    tee_error(od, run_type, overall_line_num,
                      E_POLARIZATION_FUNCTIONS, "d", "0 - 3",
                      failed);
                    fail = !(run_type & INTERACTIVE_RUN);
                    continue;
                  }
                  if (od->field.d_func) {
                    sprintf(&buffer[len], "_");
                    ++len;
                    if (od->field.d_func > 1) {
                      sprintf(&buffer[len], "%d", od->field.d_func);
                      ++len;
                    }
                    sprintf(&buffer[len], "d");
                    ++len;
                  }
                }
              }
              if ((parameter = get_args(od, "p_func"))) {
                sscanf(parameter, "%d", &(od->field.p_func));
                if ((od->field.p_func < 0) || (od->field.p_func > 3)) {
                  tee_error(od, run_type, overall_line_num,
                    E_POLARIZATION_FUNCTIONS, "p", "0 - 3",
                    failed);
                  fail = !(run_type & INTERACTIVE_RUN);
                  continue;
                }
                if (od->field.p_func) {
                  sprintf(&buffer[len], "_");
                  ++len;
                  if (od->field.p_func > 1) {
                    sprintf(&buffer[len], "%d", od->field.p_func);
                    ++len;
                  }
                  sprintf(&buffer[len], "p");
                  ++len;
                }
              }
              if ((parameter = get_args(od, "f_func"))) {
                sscanf(parameter, "%d", &(od->field.f_func));
                if ((od->field.f_func < 0) || (od->field.f_func > 1)) {
                  tee_error(od, run_type, overall_line_num,
                    E_POLARIZATION_FUNCTIONS, "f", "0 - 1",
                    failed);
                  fail = !(run_type & INTERACTIVE_RUN);
                  continue;
                }
                if (od->field.f_func) {
                  sprintf(&buffer[len], "_");
                  ++len;
                  if (od->field.f_func > 1) {
                    sprintf(&buffer[len], "%d", od->field.f_func);
                    ++len;
                  }
                  sprintf(&buffer[len], "f");
                  ++len;
                }
              }
              if ((parameter = get_args(od, "diff_sp"))) {
                if (!strncasecmp(parameter, "y", 1)) {
                  od->field.diff_sp = 1;
                  strcat(buffer, "_plus");
                }
              }
            }
          }
          else {
            strcpy(buffer, "qm_molden");
          }
        }
        else if ((!(run_type & DRY_RUN)) && (od->field.match)) {
          result = open_temp_file(od, od->file[TEMP_OBJECT_MATCH], "object_match");
          if (result) {
            tee_error(od, run_type, overall_line_num,
              E_TEMP_FILE_CANNOT_BE_OPENED_FOR_WRITING,
              od->file[TEMP_OBJECT_MATCH]->name, failed);
            return PARSE_INPUT_ERROR;
          }
          result = match_objects_with_datafile(od, regex_name[0], COSMO_INPUT_FILE);
          strcpy(buffer, regex_name[0]);
          get_dirname(buffer);
          switch (result) {
            case CANNOT_OPEN_DIRECTORY:
            tee_error(od, run_type, overall_line_num,
              E_DIR_NOT_EXISTING, buffer, failed);
            return PARSE_INPUT_ERROR;
            
            case NOT_ENOUGH_OBJECTS:
            case TOO_MANY_OBJECTS:
            tee_error(od, run_type, overall_line_num,
              E_NOT_ENOUGH_OBJECTS,
              od->grid.object_num, regex_name[0],
              failed);
            return PARSE_INPUT_ERROR;

            case CANNOT_READ_ORIGINAL_SDF:
            tee_error(od, run_type, overall_line_num,
              E_ERROR_IN_READING_SDF_FILE,
              od->task.string, failed);
            O3_ERROR_PRINT(&(od->task));
            return PARSE_INPUT_ERROR;
        
            case CANNOT_READ_TEMP_FILE:
            tee_error(od, run_type, overall_line_num,
              E_ERROR_IN_READING_TEMP_FILE, 
              od->task.string, failed);
            O3_ERROR_PRINT(&(od->task));
            return PARSE_INPUT_ERROR;

            case PREMATURE_DAT_EOF:
            tee_error(od, run_type, overall_line_num,
              E_FILE_CORRUPTED_OR_IN_WRONG_FORMAT,
              "COSMO", od->task.string,
              failed);
            O3_ERROR_PRINT(&(od->task));
            return PARSE_INPUT_ERROR;

            case OBJECTS_NOT_MATCHING:
            tee_error(od, run_type, overall_line_num,
              E_OBJECTS_NOT_MATCHING, regex_name[0],
              od->al.mol_info[od->task.data[DATA_OBJECT_NUM]]->object_id,
              od->al.mol_info[od->task.data[DATA_OBJECT_NUM]]->object_name,
              failed);
            O3_ERROR_PRINT(&(od->task));
            return PARSE_INPUT_ERROR;
          }
        }
        else if ((!(od->field.match)) && (prep_or_calc == 'c')) {
          sprintf(file_basename, "%s%c"TURBOMOLE_RIDFT_EXE,
            od->field.qm_exe_path, SEPARATOR);
          found = fexist(file_basename);
          if (!found) {
            found = is_in_path(TURBOMOLE_RIDFT_EXE, od->field.qm_exe_path);
          }
          if (!found) {
            tee_error(od, run_type, overall_line_num,
              "Cannot find "TURBOMOLE_RIDFT_EXE
              ". Please set the "
              "O3_QM_ENGINE environment variable "
              "or use the \"env qm_engine\" keyword "
              "to indicate the full path to "TURBOMOLE_RIDFT_EXE
              ".\n%s", failed);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
        }
        if (od->field.type & (PREP_TURBOMOLE_INPUT | CS3D_FIELD)) {
          turbodir = getenv("TURBODIR");
          i = 1;
          if ((!turbodir) || (!turbodir[0]) || (!dexist(turbodir))) {
            strcpy(file_basename, od->field.qm_exe_path);
            i = up_n_levels(file_basename, 2);
          }
          else {
            strcpy(file_basename, turbodir);
          }
          sprintf(buffer, "%s%cbasen", file_basename, SEPARATOR);
          if ((!dexist(buffer)) || (!i)) {
            tee_error(od, run_type, overall_line_num,
              "Please set the TURBODIR environment "
              "variable to the root dir of your "
              "TURBOMOLE installation (i.e., "
              "the directory which contains the \"basen\" "
              "folder).\n%s", failed);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
          sprintf(buffer, "TURBODIR=%s", file_basename);
          putenv(buffer);
          #ifdef WIN32
          sprintf(buffer, "TURBOWINDIR=%s", file_basename);
          putenv(buffer);
          #endif
          if (!(run_type & DRY_RUN)) {
            result = check_define(od, od->field.qm_exe_path);
            if (result && strchr(od->field.qm_exe_path, ' ')) {
              tee_error(od, run_type, overall_line_num,
                E_DEFINE_SPACED_PATH, failed);
              fail = !(run_type & INTERACTIVE_RUN);
              O3_ERROR_PRINT(&(od->task));
              continue;
            }
          }
        }
        if (!(run_type & DRY_RUN)) {
          /*
          if fields are going to be computed by QM,
          create a folder for input and output files
          */
          if ((parameter = get_args(od, "qm_dir"))) {
            strcpy(od->field.qm_dir, parameter);
            absolute_path(od->field.qm_dir);
            result = 0;
            if (!dexist(od->field.qm_dir)) {
              #ifndef WIN32
              result = mkdir(od->field.qm_dir, S_IRWXU | S_IRGRP | S_IROTH);
              #else
              result = mkdir(od->field.qm_dir);
              #endif
            }
          }
          else {
            result = open_perm_dir(od, od->default_folder,
              ((od->field.type & CS3D_FIELD) ? "cs3d" : "qm_dir"),
              od->field.qm_dir);
          }
          if (result) {
            tee_error(od, run_type, overall_line_num,
              E_TEMP_DIR_CANNOT_BE_CREATED, od->field.qm_dir,
              failed);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
          tee_printf(od, M_INPUT_OUTPUT_LOG_DIR, "qm_dir", od->field.qm_dir);
          if ((od->field.type & CS3D_FIELD) && (!(od->field.match))) {
            for (i = 0; i < od->object_num; ++i) {
              fprintf(od->file[TEMP_SORTED_MATCH]->handle,
                "%s%c%s_%04d"TURBOMOLE_COSMO_EXT"\n",
                od->field.qm_dir, SEPARATOR, cosmo_label,
                od->al.mol_info[i]->object_id);
            }
          }
        }
      }
      if (!(run_type & DRY_RUN)) {
        gettimeofday(&start, NULL);
        ++command;
        if (od->field.type & (PREP_SYBYL_INPUT | PREP_MOE_GRID_INPUT)) {
          tee_printf(od, M_TOOL_INVOKE, nesting, command, "PREPARE", line_orig);
          tee_flush(od);
          if (!(od->file[PREPINP_OUT]->handle =
            fopen(od->file[PREPINP_OUT]->name, "wb"))) {
            tee_error(od, run_type, overall_line_num,
              E_FILE_CANNOT_BE_OPENED_FOR_WRITING,
              od->file[PREPINP_OUT]->name, failed);
            fail = !(run_type & INTERACTIVE_RUN);
            continue; 
          }
          if (od->field.type & PREP_SYBYL_INPUT) {
            prep_sybyl_input(od);
          }
          else {
            prep_moe_grid_input(od);
          }
          if (od->file[PREPINP_OUT]->handle) {
            fclose(od->file[PREPINP_OUT]->handle);
            od->file[PREPINP_OUT]->handle = NULL;
          }
          tee_printf(od, M_TOOL_SUCCESS, nesting, command, "PREPARE");
          tee_flush(od);
        }  
        else if (od->field.type & PREP_MOLDEN_INPUT) {
          tee_printf(od, M_TOOL_INVOKE, nesting, command, "PREPARE", line_orig);
          tee_flush(od);
          for (i = 0; i < od->object_num; ++i) {
            result = prep_molden_input(od, i);
            switch (result) {
              case FL_CANNOT_CREATE_MDNDIR:
              sprintf(buffer, "%s%cmolden_%04d", od->field.qm_dir,
                SEPARATOR, od->al.mol_info[i]->object_id);
              tee_error(od, run_type, overall_line_num,
                E_TEMP_DIR_CANNOT_BE_CREATED,
                buffer, failed);
              return PARSE_INPUT_ERROR;
              
              case FL_CANNOT_WRITE_INP_FILE:
              sprintf(buffer, "%s%cmolden_%04d%cmolden_%04d"
                MOLDEN_INP_EXT, od->field.qm_dir, SEPARATOR,
                od->al.mol_info[i]->object_id, SEPARATOR,
                od->al.mol_info[i]->object_id);
              tee_error(od, run_type, overall_line_num,
                E_FILE_CANNOT_BE_OPENED_FOR_WRITING,
                buffer, failed);
              return PARSE_INPUT_ERROR;
            }
          }
          tee_printf(od, M_TOOL_SUCCESS, nesting, command, "PREPARE");
          tee_flush(od);
        }  
        else if (od->field.type & (QM_ELE_FIELD | QM_DEN_FIELD | CS3D_FIELD)) {
          tee_printf(od, M_TOOL_INVOKE, nesting, command,
            (prep_or_calc == 'c' ? ((od->field.type & CS3D_FIELD)
            ? "CALC_CS3D_FIELD" : "CALC_QM_FIELD") : "PREPARE"), line_orig);
          if (!(od->field.match)) {
            if (prep_or_calc == 'c') {
              if ((parameter = get_args(od, "qm_scratch"))) {
                strcpy(od->field.qm_scratch, parameter);
                absolute_path(od->field.qm_scratch);
                result = 0;
                if (!dexist(od->field.qm_scratch)) {
                  #ifndef WIN32
                  result = mkdir(od->field.qm_scratch, S_IRWXU | S_IRGRP | S_IROTH);
                  #else
                  result = mkdir(od->field.qm_scratch);
                  #endif
                }
              }
              else {
                result = open_temp_dir(od, od->temp_dir, "qm_scratch", od->field.qm_scratch);
                if (result) {
                  tee_error(od, run_type, overall_line_num,
                    E_TEMP_DIR_CANNOT_BE_CREATED, od->field.qm_scratch,
                    failed);
                    fail = !(run_type & INTERACTIVE_RUN);
                    continue;
                }
              }
              tee_printf(od, "The qm_scratch directory is:\n%s\n\n", od->field.qm_scratch);
            }
            tee_flush(od);
            result = call_obenergy(od, O3_MMFF94);
            switch (result) {
              case FL_CANNOT_CREATE_CHANNELS:
              tee_error(od, run_type, overall_line_num,
                E_CANNOT_CREATE_PIPE, failed);
              return PARSE_INPUT_ERROR;
              
              case FL_CANNOT_CHDIR:
              tee_error(od, run_type, overall_line_num,
                E_CANNOT_CHANGE_DIR, od->field.babel_exe_path, failed);
              return PARSE_INPUT_ERROR;
              
              case FL_CANNOT_CREATE_PROCESS:
              tee_error(od, run_type, overall_line_num,
                E_CANNOT_CREATE_PROCESS, "OpenBabel", failed);
              return PARSE_INPUT_ERROR;

              case OUT_OF_MEMORY:
              tee_error(od, run_type, overall_line_num,
                E_OUT_OF_MEMORY, failed);
              return PARSE_INPUT_ERROR;

              case OPENBABEL_ERROR:
              tee_error(od, run_type, overall_line_num,
                E_PROGRAM_ERROR, "OpenBabel");
              if ((od->file[TEMP_LOG]->handle = fopen
                (od->file[TEMP_LOG]->name, "rb"))) {
                while (fgets(buffer, BUF_LEN,
                  od->file[TEMP_LOG]->handle)) {
                  buffer[BUF_LEN - 1] = '\0';
                  tee_printf(od, "%s", buffer);
                }
                fclose(od->file[TEMP_LOG]->handle);
                od->file[TEMP_LOG]->handle = NULL;
                tee_printf(od, "\n%s", failed);
              }
              else {
                tee_error(od, run_type, overall_line_num,
                  E_CANNOT_READ_PROGRAM_LOG, "OpenBabel", failed);
              }
              return PARSE_INPUT_ERROR;
            }
            result = calc_field(od, (void *)((od->field.type & CS3D_FIELD)
              ? calc_cosmo_thread : calc_qm_thread), prep_or_calc);
            switch (result) {
              case CANNOT_CREATE_THREAD:
              tee_error(od, run_type, overall_line_num,
                E_THREAD_ERROR, "create",
                od->error_code, failed);
              return PARSE_INPUT_ERROR;

              case CANNOT_JOIN_THREAD:
              tee_error(od, run_type, overall_line_num,
                E_THREAD_ERROR, "join",
                od->error_code, failed);
              return PARSE_INPUT_ERROR;
            }
            for (i = 0, result = 0; i < od->object_num; ++i) {
              if (od->al.task_list[i]->code) {
                result = 1;
                tee_printf(od, "Object ID %4d:\n", od->al.mol_info[i]->object_id);
                switch (od->al.task_list[i]->code) {
                  case FL_OUT_OF_MEMORY:
                  tee_printf(od, E_OUT_OF_MEMORY, "");
                  break;

                  case FL_CANNOT_CREATE_CHANNELS:
                  tee_printf(od, E_CANNOT_CREATE_PIPE, "");
                  break;

                  case FL_CANNOT_CREATE_PROCESS:
                  tee_printf(od, E_CANNOT_CREATE_PROCESS, "QM engine", "");
                  break;

                  case FL_CANNOT_READ_OUT_FILE:
                  tee_printf(od, E_CANNOT_READ_WRITE_QM_FILE, "read", ".out",
                    od->al.task_list[i]->string, "");
                  break;

                  case FL_CANNOT_READ_INP_FILE:
                  tee_printf(od, E_CANNOT_READ_WRITE_QM_FILE, "read", ".inp",
                    od->al.task_list[i]->string, "");
                  break;

                  case FL_CANNOT_WRITE_INP_FILE:
                  tee_printf(od, E_CANNOT_READ_WRITE_QM_FILE, "write", ".inp",
                    od->al.task_list[i]->string, "");
                  break;

                  case FL_CANNOT_READ_FCHK_FILE:
                  tee_printf(od, E_CANNOT_READ_WRITE_QM_FILE, "read", FORMCHK_EXT,
                    od->al.task_list[i]->string, "");
                  break;

                  case FL_CANNOT_READ_GCUBE_FILE:
                  tee_printf(od, E_CANNOT_READ_WRITE_QM_FILE, "read", GAUSSIAN_CUBE_EXT,
                    od->al.task_list[i]->string, "");
                  break;

                  case FL_CANNOT_CREATE_SCRDIR:
                  tee_printf(od, E_QM_DIR_CANNOT_BE_CREATED, "scratch",
                    od->al.task_list[i]->string, "");
                  break;

                  case FL_ABNORMAL_TERMINATION:
                  tee_printf(od, E_CALCULATION_ERROR, "QM calculations", "");
                  break;

                  case FL_CANNOT_CHDIR:
                  tee_printf(od, E_CANNOT_CHANGE_DIR, "the QM engine folder", "");
                  break;

                  case FL_CANNOT_READ_MOL_FILE:
                  tee_printf(od, E_ERROR_IN_READING_MOL_FILE, 
                    od->al.task_list[i]->string, "");
                  break;

                  case FL_CANNOT_READ_OB_OUTPUT:
                  tee_printf(od, E_ERROR_IN_READING_OB_OUTPUT,
                    od->al.task_list[i]->string, "");
                  break;

                  case FL_UNKNOWN_ATOM_TYPE:
                  tee_printf(od, E_UNKNOWN_ATOM_TYPE, "atom", "");
                  break;
                }
                O3_ERROR_PRINT(od->al.task_list[i]);
              }
            }
            if (result) {
              tee_error(od, run_type, overall_line_num,
                E_CALCULATION_ERROR, "QM calculations", failed);
              return PARSE_INPUT_ERROR;
            }
            free_array(od->al.task_list);
            od->al.task_list = NULL;
          }
          if ((prep_or_calc == 'c') && (od->field.type & PREP_GAUSSIAN_INPUT)) {
            for (i = 0; i < od->object_num; ++i) {
              fprintf(od->file[TEMP_SORTED_MATCH]->handle,
                "%s%c%s_%04d"GAUSSIAN_CUBE_EXT"\n",
                od->field.qm_dir, SEPARATOR, od->field.qm_software,
                od->al.mol_info[i]->object_id);
            }
            sprintf(regex_name[0], "%s%c%s_####"GAUSSIAN_CUBE_EXT,
              od->field.qm_dir, SEPARATOR, od->field.qm_software);
            result = import_grid_unformatted_cube(od, 1);
            gettimeofday(&end, NULL);
            elapsed_time(od, &start, &end);
            switch (result) {
              case NOT_ENOUGH_OBJECTS:
              tee_error(od, run_type, overall_line_num,
                E_NOT_ENOUGH_OBJECTS, od->object_num,
                regex_name[0], failed);
              return PARSE_INPUT_ERROR;

              case PREMATURE_DAT_EOF:
              tee_error(od, run_type, overall_line_num,
                E_FILE_CORRUPTED_OR_IN_WRONG_FORMAT,
                "GAUSSIAN CUBE", od->file[BINARY_IN]->name,
                failed);
              return PARSE_INPUT_ERROR;

              case CANNOT_WRITE_TEMP_FILE:
              tee_error(od, run_type, overall_line_num,
                E_ERROR_IN_WRITING_TEMP_FILE, "TEMP_FIELD", failed);
              return PARSE_INPUT_ERROR;

              case OUT_OF_MEMORY:
              tee_error(od, run_type, overall_line_num,
                E_OUT_OF_MEMORY, failed);
              return PARSE_INPUT_ERROR;

              case GRID_NOT_MATCHING:
              tee_error(od, run_type, overall_line_num,
                E_GRID_NOT_MATCHING, "GAUSSIAN_CUBE",
                od->file[BINARY_IN]->name, "");
              print_grid_comparison(od);
              tee_error(od, run_type, overall_line_num, "%s\n",
                failed);
              return PARSE_INPUT_ERROR;

              default:
              tee_printf(od, M_TOOL_SUCCESS, nesting, command, "CALC_QM_FIELD");
              tee_flush(od);
            }
          }
          else if ((prep_or_calc == 'c') && ((od->field.type & PREP_FIREFLY_INPUT)
            || (od->field.type & PREP_GAMESS_INPUT))) {
            for (i = 0; i < od->object_num; ++i) {
              fprintf(od->file[TEMP_SORTED_MATCH]->handle,
                "%s%c%s_%04d"GAMESS_PUNCH_EXT"\n",
                od->field.qm_dir, SEPARATOR, od->field.qm_software,
                od->al.mol_info[i]->object_id);
            }
            sprintf(regex_name[0], "%s%c%s_####"GAMESS_PUNCH_EXT,
              od->field.qm_dir, SEPARATOR, od->field.qm_software);
            result = import_grid_formatted_cube(od, 1);
            gettimeofday(&end, NULL);
            elapsed_time(od, &start, &end);
            switch (result) {
              case NOT_ENOUGH_OBJECTS:
              tee_error(od, run_type, overall_line_num,
                E_NOT_ENOUGH_OBJECTS, od->object_num,
                regex_name[0], failed);
              return PARSE_INPUT_ERROR;

              case PREMATURE_DAT_EOF:
              tee_error(od, run_type, overall_line_num,
                E_FILE_CORRUPTED_OR_IN_WRONG_FORMAT,
                "GAMESS CUBE", od->file[ASCII_IN]->name,
                failed);
              return PARSE_INPUT_ERROR;

              case CANNOT_WRITE_TEMP_FILE:
              tee_error(od, run_type, overall_line_num,
                E_ERROR_IN_WRITING_TEMP_FILE, "TEMP_FIELD", failed);
              return PARSE_INPUT_ERROR;

              case OUT_OF_MEMORY:
              tee_error(od, run_type, overall_line_num,
                E_OUT_OF_MEMORY, failed);
              return PARSE_INPUT_ERROR;

              case GRID_NOT_MATCHING:
              tee_error(od, run_type, overall_line_num,
                E_GRID_NOT_MATCHING, "GAMESS CUBE",
                od->file[ASCII_IN]->name, "");
              print_grid_comparison(od);
              tee_error(od, run_type, overall_line_num, "%s\n", failed);
              return PARSE_INPUT_ERROR;

              default:
              tee_printf(od, M_TOOL_SUCCESS, nesting, command,
                (prep_or_calc == 'c') ? "CALC_QM_FIELD" : "PREPARE");
              tee_flush(od);
            }
          }
          else if ((prep_or_calc == 'c') && (od->field.type & PREP_TURBOMOLE_INPUT)) {
            for (i = 0; i < od->object_num; ++i) {
              fprintf(od->file[TEMP_SORTED_MATCH]->handle,
                "%s%c%s_%04d"GAMESS_PUNCH_EXT"\n",
                od->field.qm_dir, SEPARATOR, od->field.qm_software,
                od->al.mol_info[i]->object_id);
            }
            sprintf(regex_name[0], "%s%c%s_####"GAMESS_PUNCH_EXT,
              od->field.qm_dir, SEPARATOR, od->field.qm_software);
            result = import_grid_formatted_cube(od, 1);
            gettimeofday(&end, NULL);
            elapsed_time(od, &start, &end);
            switch (result) {
              case NOT_ENOUGH_OBJECTS:
              tee_error(od, run_type, overall_line_num,
                E_NOT_ENOUGH_OBJECTS, od->object_num,
                regex_name[0], failed);
              return PARSE_INPUT_ERROR;

              case PREMATURE_DAT_EOF:
              tee_error(od, run_type, overall_line_num,
                E_FILE_CORRUPTED_OR_IN_WRONG_FORMAT,
                "TURBOMOLE CUBE", od->file[ASCII_IN]->name,
                failed);
              return PARSE_INPUT_ERROR;

              case CANNOT_WRITE_TEMP_FILE:
              tee_error(od, run_type, overall_line_num,
                E_ERROR_IN_WRITING_TEMP_FILE, "TEMP_FIELD", failed);
              return PARSE_INPUT_ERROR;

              case OUT_OF_MEMORY:
              tee_error(od, run_type, overall_line_num,
                E_OUT_OF_MEMORY, failed);
              return PARSE_INPUT_ERROR;

              case GRID_NOT_MATCHING:
              tee_error(od, run_type, overall_line_num,
                E_GRID_NOT_MATCHING, "TURBOMOLE CUBE",
                od->file[ASCII_IN]->name, "");
              print_grid_comparison(od);
              tee_error(od, run_type, overall_line_num, "%s\n", failed);
              return PARSE_INPUT_ERROR;

              default:
              tee_printf(od, M_TOOL_SUCCESS, nesting, command,
                (prep_or_calc == 'c') ? "CALC_QM_FIELD" : "PREPARE");
              tee_flush(od);
            }
          }
          else if (od->field.type & CS3D_FIELD) {
            result = prep_cs3d_input(od);
            switch (result) {
              case CANNOT_WRITE_QM_INP_FILE:
              sprintf(buffer, "%s%ccs3d.inp", od->field.qm_dir, SEPARATOR);
              tee_error(od, run_type, overall_line_num,
                E_FILE_CANNOT_BE_OPENED_FOR_WRITING,
                buffer, failed);
              O3_ERROR_PRINT(&(od->task));
              return PARSE_INPUT_ERROR;
              
              case CANNOT_READ_TEMP_FILE:
              tee_error(od, run_type, overall_line_num,
                E_ERROR_IN_READING_TEMP_FILE,
                od->file[TEMP_SORTED_MATCH]->name, failed);
              O3_ERROR_PRINT(&(od->task));
              return PARSE_INPUT_ERROR;
            }
            if (prep_or_calc == 'c') {
              result = call_cs3d_program(od);
              switch (result) {
                case OUT_OF_MEMORY:
                tee_error(od, run_type, overall_line_num,
                  E_OUT_OF_MEMORY, failed);
                return PARSE_INPUT_ERROR;

                case FL_CANNOT_CREATE_PROCESS:
                tee_error(od, run_type, overall_line_num,
                  E_CANNOT_CREATE_PROCESS, "cs3d", failed);
                return PARSE_INPUT_ERROR;

                case CS3D_ERROR:
                tee_error(od, run_type, overall_line_num,
                  E_PROGRAM_ERROR, "cs3d");
                sprintf(log_fd.name, "%s%ccs3d.log", od->field.qm_dir, SEPARATOR);
                if ((log_fd.handle = fopen(log_fd.name, "rb"))) {
                  while (fgets(buffer, BUF_LEN, log_fd.handle)) {
                    buffer[BUF_LEN - 1] = '\0';
                    tee_printf(od, "%s", buffer);
                  }
                  fclose(log_fd.handle);
                  log_fd.handle = NULL;
                  tee_printf(od, "\n%s", failed);
                }
                else {
                  tee_error(od, run_type, overall_line_num,
                    E_CANNOT_READ_PROGRAM_LOG, "cs3d", failed);
                }
                return PARSE_INPUT_ERROR;
              }
              strcpy(name_list, "N_FIELD,Z_COORD,Y_COORD,X_COORD,N_OBJECT");
              sprintf(od->file[ASCII_IN]->name, "%s%ccs3d.lsp%s",
                od->field.qm_dir, SEPARATOR,
                ((od->field.compress == O3_COMPRESS_GZIP)
                ? ".gz" : ((od->field.compress == O3_COMPRESS_ZIP)
                ? ".zip" : "")));
              if (!(od->file[ASCII_IN]->handle =
                (FILE *)fzopen(od->file[ASCII_IN]->name, "rb"))) {
                tee_error(od, run_type, overall_line_num,
                  E_FILE_CANNOT_BE_OPENED_FOR_READING,
                  od->file[ASCII_IN]->name, failed);
                fail = !(run_type & INTERACTIVE_RUN);
                continue;
              }
              result = import_free_format(od, name_list, CS3D_SKIP_HEADER, &n_values);
              if (od->file[ASCII_IN]->handle) {
                fzclose((fzPtr *)(od->file[ASCII_IN]->handle));
                od->file[ASCII_IN]->handle = NULL;
              }
              if (result) {
                O3_ERROR_PRINT(&(od->task));
                switch (result) {
                  case OUT_OF_MEMORY:
                  tee_error(od, run_type, overall_line_num,
                    E_OUT_OF_MEMORY, IMPORT_FAILED);
                  return PARSE_INPUT_ERROR;

                  case WRONG_DATA_FORMAT:
                  case NOT_ENOUGH_VALUES:
                  case INCORRECT_NUMBER_OF_VALUES:
                  case EOF_FOUND:
                  case PREMATURE_EOF:
                  tee_error(od, run_type, overall_line_num,
                    E_FILE_CORRUPTED_OR_IN_WRONG_FORMAT, "ASCII",
                    od->file[ASCII_IN]->name, IMPORT_FAILED);
                  return PARSE_INPUT_ERROR;

                  case Y_VAR_LOW_SD:
                  tee_error(od, run_type, overall_line_num,
                    E_Y_VAR_LOW_SD, IMPORT_FAILED);
                  return PARSE_INPUT_ERROR;
                }
              }
            }
            tee_printf(od, M_TOOL_SUCCESS, nesting, command,
              (prep_or_calc == 'c') ? "CALC_CS3D_FIELD" : "PREPARE");
            tee_flush(od);
            if (od->file[TEMP_OBJECT_MATCH]->handle) {
              fclose(od->file[TEMP_OBJECT_MATCH]->handle);
              od->file[TEMP_OBJECT_MATCH]->handle = NULL;
            }
            if (od->file[TEMP_SORTED_MATCH]->handle) {
              fclose(od->file[TEMP_SORTED_MATCH]->handle);
              od->file[TEMP_SORTED_MATCH]->handle = NULL;
            }
          }
        }
        else if ((prep_or_calc == 'c') && (od->field.type & MD_GRID_FIELD)) {
          tee_printf(od, M_TOOL_INVOKE, nesting, command, "CALC_MD_GRID", line_orig);
          /*
          if fields are going to be computed by MD GRID,
          create a folder for input and output files
          */
          memset(file_basename, 0, BUF_LEN);
          sprintf(buffer, "md_grid_%s", od->field.probe.atom_name);
          if ((parameter = get_args(od, "md_grid_dir"))) {
            strcpy(od->field.md_grid_dir, parameter);
            absolute_path(od->field.md_grid_dir);
            result = 0;
            if (!dexist(od->field.md_grid_dir)) {
              #ifndef WIN32
              result = mkdir(od->field.md_grid_dir, S_IRWXU | S_IRGRP | S_IROTH);
              #else
              result = mkdir(od->field.md_grid_dir);
              #endif
            }
          }
          else {
            result = open_perm_dir(od, od->default_folder, buffer, od->field.md_grid_dir);
          }
          if (result) {
            tee_error(od, run_type, overall_line_num,
              E_TEMP_DIR_CANNOT_BE_CREATED, od->field.md_grid_dir, failed);
              fail = !(run_type & INTERACTIVE_RUN);
              continue;
          }
          tee_printf(od, M_INPUT_OUTPUT_LOG_DIR, buffer, od->field.md_grid_dir);
          tee_flush(od);
          result = call_obenergy(od, O3_MMFF94);
          switch (result) {
            case FL_CANNOT_CREATE_CHANNELS:
            tee_error(od, run_type, overall_line_num,
              E_CANNOT_CREATE_PIPE, failed);
            return PARSE_INPUT_ERROR;
            
            case FL_CANNOT_CHDIR:
            tee_error(od, run_type, overall_line_num,
              E_CANNOT_CHANGE_DIR, od->field.babel_exe_path, failed);
            return PARSE_INPUT_ERROR;
            
            case FL_CANNOT_CREATE_PROCESS:
            tee_error(od, run_type, overall_line_num,
              E_CANNOT_CREATE_PROCESS, "OpenBabel", failed);
            return PARSE_INPUT_ERROR;

            case OUT_OF_MEMORY:
            tee_error(od, run_type, overall_line_num,
              E_OUT_OF_MEMORY, failed);
            return PARSE_INPUT_ERROR;

            case OPENBABEL_ERROR:
            tee_error(od, run_type, overall_line_num,
              E_PROGRAM_ERROR, "OpenBabel");
            if ((od->file[TEMP_LOG]->handle = fopen
              (od->file[TEMP_LOG]->name, "rb"))) {
              while (fgets(buffer, BUF_LEN,
                od->file[TEMP_LOG]->handle)) {
                buffer[BUF_LEN - 1] = '\0';
                tee_printf(od, "%s", buffer);
              }
              fclose(od->file[TEMP_LOG]->handle);
              od->file[TEMP_LOG]->handle = NULL;
              tee_printf(od, "\n%s", failed);
            }
            else {
              tee_error(od, run_type, overall_line_num,
                E_CANNOT_READ_PROGRAM_LOG, "OpenBabel", failed);
            }
            return PARSE_INPUT_ERROR;
          }
          result = calc_field(od, (void *)calc_md_grid_thread, prep_or_calc);
          switch (result) {
            case OUT_OF_MEMORY:
            tee_error(od, run_type, overall_line_num,
              E_OUT_OF_MEMORY, failed);
            return PARSE_INPUT_ERROR;

            case CANNOT_CREATE_THREAD:
            tee_error(od, run_type, overall_line_num,
              E_THREAD_ERROR, "create",
              od->error_code, failed);
            return PARSE_INPUT_ERROR;

            case CANNOT_JOIN_THREAD:
            tee_error(od, run_type, overall_line_num,
              E_THREAD_ERROR, "join",
              od->error_code, failed);
            return PARSE_INPUT_ERROR;

            case CANNOT_COPY_GRUB_DAT:
            tee_error(od, run_type, overall_line_num,
              "Cannot copy %s%c"GRUB_FILENAME" to %s%c"GRUB_FILENAME".\n%s",
              od->field.md_grid_exe_path, SEPARATOR,
              od->field.md_grid_dir, SEPARATOR, failed);
            return PARSE_INPUT_ERROR;
          }
          for (i = 0, result = 0; i < od->object_num; ++i) {
            if (od->al.task_list[i]->code) {
              result = 1;
              tee_printf(od, "Object ID %4d:\n", od->al.mol_info[i]->object_id);
              switch (od->al.task_list[i]->code) {
                case FL_OUT_OF_MEMORY:
                tee_printf(od, E_OUT_OF_MEMORY, "");
                break;

                case FL_CANNOT_READ_TEMP_FILE:
                tee_printf(od, E_ERROR_IN_READING_TEMP_FILE, 
                  od->al.task_list[i]->string, "");
                break;

                case FL_CANNOT_WRITE_TEMP_FILE:
                tee_printf(od, E_ERROR_IN_WRITING_TEMP_FILE,
                  od->al.task_list[i]->string, "");
                break;

                case FL_CANNOT_CREATE_CHANNELS:
                tee_printf(od, E_CANNOT_CREATE_PIPE, "");
                break;

                case FL_CANNOT_CREATE_PROCESS:
                tee_printf(od, E_CANNOT_CREATE_PROCESS, "GRIN", "");
                break;

                case FL_UNKNOWN_ATOM_TYPE:
                tee_printf(od, E_UNKNOWN_ATOM_TYPE, "atom", "");
                break;

                case FL_CANNOT_CHDIR:
                tee_printf(od, E_CANNOT_CHANGE_DIR, od->field.md_grid_dir, "");
                break;

                case FL_CANNOT_READ_MOL_FILE:
                tee_printf(od, E_ERROR_IN_READING_MOL_FILE,
                  od->al.task_list[i]->string, "");
                break;

                case FL_CANNOT_READ_OB_OUTPUT:
                tee_printf(od, E_ERROR_IN_READING_OB_OUTPUT,
                  od->al.task_list[i]->string, "");
                break;

                case FL_GRIN_ERROR:
                tee_printf(od, E_PROGRAM_ERROR, "GRIN");
                sprintf(log_fd.name, "%s%c%04d.log", od->field.md_grid_dir,
                  SEPARATOR, od->al.mol_info[i]->object_id);
                if ((log_fd.handle = fopen(log_fd.name, "rb"))) {
                  while (fgets(buffer, BUF_LEN, log_fd.handle)) {
                    buffer[BUF_LEN - 1] = '\0';
                    tee_printf(od, "%s", buffer);
                  }
                  fclose(log_fd.handle);
                  log_fd.handle = NULL;
                }
                else {
                  tee_printf(od, E_CANNOT_READ_PROGRAM_LOG, "GRIN", "");
                }
                break;
              }
              O3_ERROR_PRINT(od->al.task_list[i]);
            }
          }
          if (result) {
            tee_error(od, run_type, overall_line_num,
              E_CALCULATION_ERROR, "MD GRID computations", failed);
            return PARSE_INPUT_ERROR;
          }
          free_array(od->al.task_list);
          od->al.task_list = NULL;
          result = call_md_grid_program(od);
          switch (result) {
            case OUT_OF_MEMORY:
            tee_error(od, run_type, overall_line_num,
              E_OUT_OF_MEMORY, failed);
            O3_ERROR_PRINT(&(od->task));
            return PARSE_INPUT_ERROR;

            case CANNOT_READ_TEMP_FILE:
            tee_error(od, run_type, overall_line_num,
              E_ERROR_IN_READING_TEMP_FILE, 
              od->task.string, failed);
            O3_ERROR_PRINT(&(od->task));
            return PARSE_INPUT_ERROR;

            case CANNOT_WRITE_TEMP_FILE:
            tee_error(od, run_type, overall_line_num,
              E_ERROR_IN_WRITING_TEMP_FILE,
              od->task.string, failed);
            O3_ERROR_PRINT(&(od->task));
            return PARSE_INPUT_ERROR;

            case FL_CANNOT_CREATE_CHANNELS:
            tee_error(od, run_type, overall_line_num,
              E_CANNOT_CREATE_PIPE, failed);
            O3_ERROR_PRINT(&(od->task));
            return PARSE_INPUT_ERROR;

            case FL_CANNOT_CREATE_PROCESS:
            tee_error(od, run_type, overall_line_num,
              E_CANNOT_CREATE_PROCESS, "MD GRID", failed);
            O3_ERROR_PRINT(&(od->task));
            return PARSE_INPUT_ERROR;

            case GRID_ERROR:
            tee_error(od, run_type, overall_line_num,
              E_MD_GRID_ERROR, "MD GRID", failed);
            if ((od->file[TEMP_LOG]->handle = fopen(od->file[TEMP_LOG]->name, "rb"))) {
              while (fgets(buffer, BUF_LEN, od->file[TEMP_LOG]->handle)) {
                buffer[BUF_LEN - 1] = '\0';
                tee_printf(od, "%s", buffer);
              }
              fclose(od->file[TEMP_LOG]->handle);
              od->file[TEMP_LOG]->handle = NULL;
            }
            O3_ERROR_PRINT(&(od->task));
            return PARSE_INPUT_ERROR;
          }
          sprintf(od->file[BINARY_IN]->name, "%s%cgrid.kont",
            od->field.md_grid_dir, SEPARATOR);
          if (!(od->file[BINARY_IN]->handle =
            fopen(od->file[BINARY_IN]->name, "rb"))) {
            tee_error(od, run_type, overall_line_num,
              E_FILE_CANNOT_BE_OPENED_FOR_READING,
              od->file[BINARY_IN]->name, failed);
            return PARSE_INPUT_ERROR;
          }
          result = import_gridkont(od, 0);
          gettimeofday(&end, NULL);
          elapsed_time(od, &start, &end);
          /*
          check for errors
          */
          switch (result) {
            case PREMATURE_DAT_EOF:
            tee_error(od, run_type, overall_line_num,
              E_FILE_CORRUPTED_OR_IN_WRONG_FORMAT, ".kont",
              od->file[BINARY_IN]->name, failed);
            return PARSE_INPUT_ERROR;

            case CANNOT_WRITE_TEMP_FILE:
            tee_error(od, run_type, overall_line_num,
              E_ERROR_IN_WRITING_TEMP_FILE, "TEMP_FIELD", failed);
            return PARSE_INPUT_ERROR;

            case OUT_OF_MEMORY:
            tee_error(od, run_type, overall_line_num,
              E_OUT_OF_MEMORY, failed);
            return PARSE_INPUT_ERROR;

            case OBJECTS_NOT_MATCHING:
            tee_error(od, run_type, overall_line_num,
              E_GRIDKONT_NOT_MATCHING,
              od->file[BINARY_IN]->name, od->newgrid.object_num,
              od->object_num, failed);
            return PARSE_INPUT_ERROR;

            case GRID_NOT_MATCHING:
            tee_error(od, run_type, overall_line_num,
              E_GRID_NOT_MATCHING, "GRIDKONT",
              od->file[BINARY_IN]->name, "");
            print_grid_comparison(od);
            tee_error(od, run_type, overall_line_num, "%s\n", failed);
            return PARSE_INPUT_ERROR;

            default:
            tee_printf(od, M_TOOL_SUCCESS, nesting, command, "CALC_MD_GRID");
            tee_flush(od);
          }
          if (od->file[BINARY_IN]->handle) {
            fclose(od->file[BINARY_IN]->handle);
            od->file[BINARY_IN]->handle = NULL;
          }
        }
        else if (prep_or_calc == 'c') {
          tee_printf(od, M_TOOL_INVOKE, nesting, command, "CALC_MM_FIELD", line_orig);
          tee_flush(od);
          result = call_obenergy(od, (int)(od->field.force_field));
          switch (result) {
            case FL_CANNOT_CREATE_CHANNELS:
            tee_error(od, run_type, overall_line_num,
              E_CANNOT_CREATE_PIPE, failed);
            return PARSE_INPUT_ERROR;
            
            case FL_CANNOT_CHDIR:
            tee_error(od, run_type, overall_line_num,
              E_CANNOT_CHANGE_DIR, od->field.babel_exe_path, failed);
            return PARSE_INPUT_ERROR;
            
            case FL_CANNOT_CREATE_PROCESS:
            tee_error(od, run_type, overall_line_num,
              E_CANNOT_CREATE_PROCESS, "OpenBabel", failed);
            return PARSE_INPUT_ERROR;

            case OUT_OF_MEMORY:
            tee_error(od, run_type, overall_line_num,
              E_OUT_OF_MEMORY, failed);
            return PARSE_INPUT_ERROR;

            case OPENBABEL_ERROR:
            tee_error(od, run_type, overall_line_num,
              E_PROGRAM_ERROR, "OpenBabel");
            if ((od->file[TEMP_LOG]->handle = fopen
              (od->file[TEMP_LOG]->name, "rb"))) {
              while (fgets(buffer, BUF_LEN,
                od->file[TEMP_LOG]->handle)) {
                buffer[BUF_LEN - 1] = '\0';
                tee_printf(od, "%s", buffer);
              }
              fclose(od->file[TEMP_LOG]->handle);
              od->file[TEMP_LOG]->handle = NULL;
              tee_printf(od, "\n%s", failed);
            }
            else {
              tee_error(od, run_type, overall_line_num,
                E_CANNOT_READ_PROGRAM_LOG, "OpenBabel", failed);
            }
            return PARSE_INPUT_ERROR;
          }
          result = calc_field(od, (void *)calc_mm_thread, 0);
          switch (result) {
            case OUT_OF_MEMORY:
            tee_error(od, run_type, overall_line_num,
              E_OUT_OF_MEMORY, failed);
            return PARSE_INPUT_ERROR;

            case CANNOT_CREATE_THREAD:
            tee_error(od, run_type, overall_line_num,
              E_THREAD_ERROR, "create",
              od->error_code, failed);
            return PARSE_INPUT_ERROR;

            case CANNOT_JOIN_THREAD:
            tee_error(od, run_type, overall_line_num,
              E_THREAD_ERROR, "join",
              od->error_code, failed);
            return PARSE_INPUT_ERROR;
          }
          for (i = 0, result = 0; i < od->object_num; ++i) {
            if (od->al.task_list[i]->code) {
              result = 1;
              tee_printf(od, "Object ID %4d:\n", od->al.mol_info[i]->object_id);
              switch (od->al.task_list[i]->code) {
                case FL_OUT_OF_MEMORY:
                tee_printf(od, E_OUT_OF_MEMORY, "");
                break;

                case FL_CANNOT_READ_OB_OUTPUT:
                tee_printf(od, E_ERROR_IN_READING_OB_OUTPUT,
                  od->al.task_list[i]->string, "");
                break;

                case FL_CANNOT_READ_MOL_FILE:
                tee_printf(od, E_ERROR_IN_READING_MOL_FILE,
                  od->al.task_list[i]->string, "");
                break;

                case FL_CANNOT_WRITE_TEMP_FILE:
                tee_printf(od, E_ERROR_IN_WRITING_TEMP_FILE,
                  od->al.task_list[i]->string, "");
                break;

                case FL_UNKNOWN_ATOM_TYPE:
                tee_printf(od, E_UNKNOWN_ATOM_TYPE, "atom", "");
                break;
              }
              O3_ERROR_PRINT(od->al.task_list[i]);
            }
          }
          if (result) {
            tee_error(od, run_type, overall_line_num,
              E_CALCULATION_ERROR, "MM calculations", failed);
            return PARSE_INPUT_ERROR;
          }
          update_field_object_attr(od, VERBOSE_BIT);
          result = calc_active_vars(od, FULL_MODEL);
          switch (result) {
            case CANNOT_READ_TEMP_FILE:
            tee_error(od, run_type, overall_line_num,
              E_ERROR_IN_READING_TEMP_FILE, "TEMP_FIELD", failed);
            return PARSE_INPUT_ERROR;

            case Y_VAR_LOW_SD:
            tee_error(od, run_type, overall_line_num,
              E_Y_VAR_LOW_SD, failed);
            return PARSE_INPUT_ERROR;
          }
          tee_printf(od, M_TOOL_SUCCESS, nesting, command, "CALC_MM_FIELD");
          tee_flush(od);
        }
        free_array(od->al.task_list);
        od->al.task_list = NULL;
      }
      else if (prep_or_calc == 'c') {
        /*
        fake data for DRY_RUN
        */
        od->field_num = 2;
        od->valid &= (~PLS_BIT);
      }
    }
    else if (!strcasecmp(arg->me[0], "load")) {
      gettimeofday(&start, NULL);
      if (!(od->field.babel_exe_path[0])) {
        tee_error(od, run_type, overall_line_num,
          E_OPENBABEL_PATH, LOAD_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      /*
      user wants to LOAD a .dat file; check
      if the file name was specified
      */
      if (!(parameter = get_args(od, "file"))) {
        tee_error(od, run_type, overall_line_num,
          "Please specify an "PACKAGE_NAME" file "
          "from which data should be loaded.\n%s",
          LOAD_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      strcpy(od->file[DAT_IN]->name, parameter);
      /*
      the following is for prep_qm_input, so that the default dir
      for QM input files is defined also if user LOADs a DAT file
      instead of IMPORTing a SDF file
      */
      strcpy(od->default_folder, od->file[DAT_IN]->name);
      get_dirname(od->default_folder);
      if (!(run_type & DRY_RUN)) {
        ++command;
        tee_printf(od, M_TOOL_INVOKE, nesting, command, "LOAD", line_orig);
        tee_flush(od);
        close_files(od, MAX_FILES);
        remove_temp_files(od->package_code);
        free_x_var_array(od);
        if (od->pel.pymol_old_object_id) {
          int_perm_free(od->pel.pymol_old_object_id);
          od->pel.pymol_old_object_id = NULL;
        }
        if (od->pel.jmol_old_object_id) {
          int_perm_free(od->pel.jmol_old_object_id);
          od->pel.jmol_old_object_id = NULL;
        }
        if (!(od->file[DAT_IN]->handle = (FILE *)
          fzopen(od->file[DAT_IN]->name, "rb"))) {
          tee_error(od, run_type, overall_line_num,
            E_FILE_CANNOT_BE_OPENED_FOR_READING,
            parameter, LOAD_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        /*
        load .dat file
        */
        result = open_temp_dir(od, NULL, "mol_dir", od->field.mol_dir);
        if (result) {
          tee_error(od, run_type, overall_line_num,
            E_TEMP_DIR_CANNOT_BE_CREATED, od->field.mol_dir,
            LOAD_FAILED);
            fzclose((fzPtr *)(od->file[DAT_IN]->handle));
            od->file[DAT_IN]->handle = NULL;
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
        }
        sprintf(od->file[TEMP_MOLFILE]->name, "%s%cloaded_molfile.sdf",
          od->field.mol_dir, SEPARATOR);
        result = load_dat(od, DAT_IN, VERBOSE_BIT);
        gettimeofday(&end, NULL);
        elapsed_time(od, &start, &end);
        switch (result) {
          case PREMATURE_DAT_EOF:
          tee_error(od, run_type, overall_line_num,
            E_FILE_CORRUPTED_OR_IN_WRONG_FORMAT, ".dat",
            parameter, LOAD_FAILED);
          return PARSE_INPUT_ERROR;

          case STATUS_INCONSISTENCY:
          tee_error(od, run_type, overall_line_num,
            "In the .dat file \"%s\" the active/inactive "
            "status of objects is not consistent among the "
            "different fields.\n%s",
            parameter, LOAD_FAILED);
          return PARSE_INPUT_ERROR;

          case OUT_OF_MEMORY:
          tee_error(od, run_type, overall_line_num,
            E_OUT_OF_MEMORY, LOAD_FAILED);
          return PARSE_INPUT_ERROR;

          case Y_VAR_LOW_SD:
          tee_error(od, run_type, overall_line_num,
            E_Y_VAR_LOW_SD, LOAD_FAILED);
          return PARSE_INPUT_ERROR;

          case CANNOT_READ_TEMP_FILE:
          tee_error(od, run_type, overall_line_num,
            E_ERROR_IN_READING_TEMP_FILE, "TEMP_FIELD", LOAD_FAILED);
          return PARSE_INPUT_ERROR;

          case FL_CANNOT_CREATE_CHANNELS:
          tee_error(od, run_type, overall_line_num,
            E_CANNOT_CREATE_PIPE, LOAD_FAILED);
          return PARSE_INPUT_ERROR;

          case FL_CANNOT_CHDIR:
          tee_error(od, run_type, overall_line_num,
            E_CANNOT_CHANGE_DIR, od->field.babel_exe_path, LOAD_FAILED);
          return PARSE_INPUT_ERROR;

          case FL_CANNOT_CREATE_PROCESS:
          tee_error(od, run_type, overall_line_num,
            E_CANNOT_CREATE_PROCESS, "OpenBabel", LOAD_FAILED);
          return PARSE_INPUT_ERROR;

          case OPENBABEL_ERROR:
          tee_error(od, run_type, overall_line_num,
            E_PROGRAM_ERROR, "OpenBabel");
          if ((od->file[TEMP_LOG]->handle = fopen
            (od->file[TEMP_LOG]->name, "rb"))) {
            while (fgets(buffer, BUF_LEN,
              od->file[TEMP_LOG]->handle)) {
              buffer[BUF_LEN - 1] = '\0';
              tee_printf(od, "%s", buffer);
            }
            tee_printf(od, "\n%s", LOAD_FAILED);
            fclose(od->file[TEMP_LOG]->handle);
            od->file[TEMP_LOG]->handle = NULL;
          }
          else {
            tee_error(od, run_type, overall_line_num,
              E_CANNOT_READ_PROGRAM_LOG, "OpenBabel", LOAD_FAILED);
          }
          return PARSE_INPUT_ERROR;
        }
        if (update_mol(od)) {
          tee_error(od, run_type, overall_line_num,
            E_ERROR_IN_READING_TEMP_FILE,
            od->task.string, LOAD_FAILED);
          O3_ERROR_PRINT(&(od->task));
          return PARSE_INPUT_ERROR;
        }
        update_pymol(od);
        update_jmol(od);
        tee_printf(od, M_TOOL_SUCCESS, nesting, command, "LOAD");
        tee_flush(od);
      }
      else {
        /*
        fake data for DRY_RUN
        */
        od->grid.nodes[0] = 10;
        od->grid.nodes[1] = 10;
        od->grid.nodes[2] = 10;
        od->grid.object_num = 10;
        od->grid.struct_num = 10;
        od->object_num = 10;
        od->field_num = 2;
        od->valid &= (~PLS_BIT);
        od->valid |= SDF_BIT;
      }
    }
    else if (!strcasecmp(arg->me[0], "save")) {
      gettimeofday(&start, NULL);
      /*
      user wants to SAVE a .dat file; check
      if the file name was specified
      */
      if (!(parameter = get_args(od, "file"))) {
        tee_error(od, run_type, overall_line_num, "Please specify a file "
          "where data should be saved.\n%s",
          SAVE_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      strcpy(od->file[DAT_OUT]->name, parameter);
      if (!(run_type & DRY_RUN)) {
        if (!(od->grid.object_num)) {
          tee_error(od, run_type, overall_line_num,
            E_NO_OBJECTS_PRESENT, SAVE_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        /*
        save .dat file
        */
        if (!(od->file[DAT_OUT]->handle = (FILE *)
          fzopen(od->file[DAT_OUT]->name, "wb"))) {
          tee_error(od, run_type, overall_line_num,
            E_TEMP_FILE_CANNOT_BE_OPENED_FOR_WRITING,
            od->file[DAT_OUT]->name, SAVE_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        ++command;
        tee_printf(od, M_TOOL_INVOKE, nesting, command, "SAVE", line_orig);
        tee_flush(od);
        result = save_dat(od, DAT_OUT);
        gettimeofday(&end, NULL);
        elapsed_time(od, &start, &end);
        switch (result) {
          case PREMATURE_DAT_EOF:
          tee_error(od, run_type, overall_line_num,
            E_ERROR_IN_WRITING_TEMP_FILE,
            od->file[DAT_OUT]->name, SAVE_FAILED);
          return PARSE_INPUT_ERROR;

          case OUT_OF_MEMORY:
          tee_error(od, run_type, overall_line_num,
            E_OUT_OF_MEMORY, SAVE_FAILED);
          return PARSE_INPUT_ERROR;
        }
        tee_printf(od, M_TOOL_SUCCESS, nesting, command, "SAVE");
        tee_flush(od);
      }
    }
    else if (!strcasecmp(arg->me[0], "cutoff")) {
      gettimeofday(&start, NULL);
      /*
      user wants to set minimum, maximum values inside
      all blocks of x variables (all fields) or a subset of them
      */
      if (!(parameter = get_args(od, "type"))) {
        tee_error(od, run_type, overall_line_num,
          "Please specify the type of cutoff "
          "(MIN | MAX) to be operated.\n%s",
          CUTOFF_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      if (!strncasecmp(parameter, "min", 3)) {
        type = CUTOFF_MIN;
      }
      else if (!strncasecmp(parameter, "max", 3)) {
        type = CUTOFF_MAX;
      }
      else {
        tee_error(od, run_type, overall_line_num,
          "Error while parsing the type of cutoff "
          "(MIN | MAX) to be operated.\n%s",
          CUTOFF_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      if ((parameter = get_args(od, "smooth"))) {
        if (!strncasecmp(parameter, "quad", 4)) {
          type |= CUTOFF_QUADRATIC;
        }
        else if (strncasecmp(parameter, "none", 4)) {
          tee_error(od, run_type, overall_line_num,
            "Only \"NONE\" and \"QUADRATIC\" "
            "cutoff smoothing paradigms "
            "are available.\n%s",
            CUTOFF_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
      }
      if (!(parameter = get_args(od, "level"))) {
        tee_error(od, run_type, overall_line_num,
          "Please specify the level of cutoff "
          "to be operated.\n%s",
          CUTOFF_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      sscanf(parameter, "%lf", &level);
      if (!(run_type & DRY_RUN)) {
        if (!(od->field_num)) {
          tee_error(od, run_type, overall_line_num,
            E_NO_FIELDS_PRESENT, CUTOFF_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        strcpy(comma_hyphen_list, "all");
        if ((parameter = get_args(od, "field_list"))) {
          strcpy(comma_hyphen_list, parameter);
        }
        result = parse_comma_hyphen_list_to_array
          (od, comma_hyphen_list, FIELD_LIST);
        switch (result) {
          case OUT_OF_MEMORY:
          tee_error(od, run_type, overall_line_num,
            E_OUT_OF_MEMORY, CUTOFF_FAILED);
          return PARSE_INPUT_ERROR;

          case INVALID_LIST_RANGE:
          tee_error(od, run_type, overall_line_num,
            E_LIST_PARSING, "fields", "CUTOFF",
            CUTOFF_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        result = set(od, FIELD_LIST, OPERATE_BIT, 1, SILENT);
        if (result) {
          tee_error(od, run_type, overall_line_num,
            E_ALLOWED_FIELD_RANGE,
            od->field_num, CUTOFF_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        /*
        perform cutoff operation
        */
        ++command;
        tee_printf(od, M_TOOL_INVOKE, nesting, command, "CUTOFF", line_orig);
        tee_flush(od);
        result = cutoff(od, type, level);
        gettimeofday(&end, NULL);
        elapsed_time(od, &start, &end);
        switch (result) {
          case Y_VAR_LOW_SD:
          tee_error(od, run_type, overall_line_num,
            E_Y_VAR_LOW_SD, CUTOFF_FAILED);
          return PARSE_INPUT_ERROR;

          case OUT_OF_MEMORY:
          tee_error(od, run_type, overall_line_num,
            E_OUT_OF_MEMORY, CUTOFF_FAILED);
          return PARSE_INPUT_ERROR;

          default:
          tee_printf(od, M_TOOL_SUCCESS, nesting, command, "CUTOFF");
          tee_flush(od);
        }
      }
    }
    else if (!strcasecmp(arg->me[0], "transform")) {
      gettimeofday(&start, NULL);
      /*
      user wants to transform x or y variables
      */
      if (!(parameter = get_args(od, "type"))) {
        tee_error(od, run_type, overall_line_num,
          "Please specify the type of variable "
          "on which the transformation "
          "should be carried out.\n%s",
          TRANSFORM_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      if (!strncasecmp(parameter, "x", 1)) {
        type = 'X';
      }
      else if (!strncasecmp(parameter, "y", 1)) {
        type = 'Y';
      }
      else {
        tee_error(od, run_type, overall_line_num,
          "Only \"X\" and \"Y\" variable types "
          "are allowed.\n%s",
          TRANSFORM_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      factor = 0.0;
      if (!(parameter = get_args(od, "operation"))) {
        tee_error(od, run_type, overall_line_num,
          E_TRANSFORM_MISSING_OPERATION, type,
          TRANSFORM_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      if (!strncasecmp(parameter, "mul", 3)) {
        operation = MULTIPLY;
        if (!(parameter = get_args(od, "value"))) {
          tee_error(od, run_type, overall_line_num,
            E_TRANSFORM_MISSING_VALUE, type,
            "multiplied for", TRANSFORM_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        sscanf(parameter, "%lf", &factor);
      }
      else if (!strncasecmp(parameter, "div", 3)) {
        operation = DIVIDE;
        if (!(parameter = get_args(od, "value"))) {
          tee_error(od, run_type, overall_line_num,
            E_TRANSFORM_MISSING_VALUE, type,
            "divided by", TRANSFORM_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        sscanf(parameter, "%lf", &factor);
        if (factor < ALMOST_ZERO) {
          tee_error(od, run_type, overall_line_num,
            "Cannot divide by zero.\n%s.",
            TRANSFORM_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
      }
      else if (!strncasecmp(parameter, "sum", 3)) {
        operation = SUM;
        if (!(parameter = get_args(od, "value"))) {
          tee_error(od, run_type, overall_line_num,
            E_TRANSFORM_MISSING_VALUE, type,
            "added", TRANSFORM_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        sscanf(parameter, "%lf", &factor);
      }
      else if (!strncasecmp(parameter, "sub", 3)) {
        operation = SUBTRACT;
        if (!(parameter = get_args(od, "value"))) {
          tee_error(od, run_type, overall_line_num,
            E_TRANSFORM_MISSING_VALUE, type,
            "subtracted", TRANSFORM_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        sscanf(parameter, "%lf", &factor);
      }
      else if (!strncasecmp(parameter, "opp", 3)) {
        operation = OPPOSITE;
      }
      else if (!strncasecmp(parameter, "abs", 3)) {
        operation = ABSOLUTE_VALUE;
      }
      else if (!strncasecmp(parameter, "pow", 3)) {
        operation = POWER;
        if (!(parameter = get_args(od, "exponent"))) {
          tee_error(od, run_type, overall_line_num,
            E_TRANSFORM_MISSING_EXPONENT, type,
            TRANSFORM_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        sscanf(parameter, "%lf", &factor);
      }
      else if (!strcasecmp(parameter, "log10")) {
        operation = LOGARITHM_BASE_10;
      }
      else if (!strcasecmp(parameter, "ln")) {
        operation = LOGARITHM_BASE_E;
      }
      else if (!strcasecmp(parameter, "log_base_n")) {
        operation = LOGARITHM_BASE_N;
        if (!(parameter = get_args(od, "n"))) {
          tee_error(od, run_type, overall_line_num,
            E_TRANSFORM_MISSING_BASE, type,
            TRANSFORM_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        sscanf(parameter, "%lf", &factor);
        if (factor <= 0.) {
          tee_error(od, run_type, overall_line_num,
            "The logarithm base should be "
            "a positive value.\n%s.",
            TRANSFORM_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
      }
      else {
        tee_error(od, run_type, overall_line_num,
          "Only \"SUM\", \"SUBTRACT\", \"MULTIPLY\", "
          "\"DIVIDE\", \"OPPOSITE\", \"ABS\", \"POWER\", "
          "\"LOG10\", \"LN\", and \"LOG_BASE_N\" "
          "operations are allowed.\n%s",
          TRANSFORM_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }

      if (!(run_type & DRY_RUN)) {
        if (!(od->field_num)) {
          tee_error(od, run_type, overall_line_num,
            E_NO_FIELDS_PRESENT, TRANSFORM_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        strcpy(comma_hyphen_list, "all");
        if (type == 'X') {
          if ((parameter = get_args(od, "field_list"))) {
            strcpy(comma_hyphen_list, parameter);
          }
          result = parse_comma_hyphen_list_to_array
            (od, comma_hyphen_list, FIELD_LIST);
          switch (result) {
            case OUT_OF_MEMORY:
            tee_error(od, run_type, overall_line_num,
              E_OUT_OF_MEMORY, TRANSFORM_FAILED);
            return PARSE_INPUT_ERROR;

            case INVALID_LIST_RANGE:
            tee_error(od, run_type, overall_line_num,
              E_LIST_PARSING, "fields", "TRANFORM",
              TRANSFORM_FAILED);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
          result = set(od, FIELD_LIST, OPERATE_BIT, 1, SILENT);
          if (result) {
            tee_error(od, run_type, overall_line_num,
              E_ALLOWED_FIELD_RANGE,
              od->field_num, TRANSFORM_FAILED);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
        }
        else {
          if ((parameter = get_args(od, "y_var_list"))) {
            strcpy(comma_hyphen_list, parameter);
          }
          result = parse_comma_hyphen_list_to_array
            (od, comma_hyphen_list, Y_VAR_LIST);
          switch (result) {
            case OUT_OF_MEMORY:
            tee_error(od, run_type, overall_line_num,
              E_OUT_OF_MEMORY, TRANSFORM_FAILED);
            return PARSE_INPUT_ERROR;

            case INVALID_LIST_RANGE:
            tee_error(od, run_type, overall_line_num,
              E_LIST_PARSING, "Y variables",
              "TRANSFORM", TRANSFORM_FAILED);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
          result = set(od, Y_VAR_LIST, OPERATE_BIT, 1, SILENT);
          if (result) {
            tee_error(od, run_type, overall_line_num,
              E_Y_VAR_ALLOWED_RANGE,
              od->y_vars, TRANSFORM_FAILED);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
        }
        /*
        perform transformation
        */
        ++command;
        tee_printf(od, M_TOOL_INVOKE, nesting, command, "TRANSFORM", line_orig);
        tee_flush(od);
        result = transform(od, type, operation, factor);
        gettimeofday(&end, NULL);
        elapsed_time(od, &start, &end);
        switch (result) {
          case OUT_OF_MEMORY:
          tee_error(od, run_type, overall_line_num,
            E_OUT_OF_MEMORY, TRANSFORM_FAILED);
          return PARSE_INPUT_ERROR;

          default:
          tee_printf(od, M_TOOL_SUCCESS, nesting, command, "TRANSFORM");
          tee_flush(od);
        }
      }
    }
    else if (!strcasecmp(arg->me[0], "zero")) {
      gettimeofday(&start, NULL);
      /*
      user wants to zero all values below
      a certain threshold
      */
      type = ZERO_ALL;
      if ((parameter = get_args(od, "sign"))) {
        if (!strncasecmp(parameter, "pos", 3)) {
          type = ZERO_POS;
        }
        else if (!strncasecmp(parameter, "neg", 3)) {
          type = ZERO_NEG;
        }
        else if (strncasecmp(parameter, "both", 4)) {
          tee_error(od, run_type, overall_line_num,
            "Allowed values for sign are \"POS\", "
            "\"NEG\" and \"BOTH\".\n%s", ZERO_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
      }
      if (!(parameter = get_args(od, "level"))) {
        tee_error(od, run_type, overall_line_num,
          "Please specify the absolute value threshold "
          "below which ZERO should operate.\n%s",
          ZERO_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      sscanf(parameter, "%lf", &level);

      if (!(run_type & DRY_RUN)) {
        if (!(od->field_num)) {
          tee_error(od, run_type, overall_line_num,
            E_NO_FIELDS_PRESENT, ZERO_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        strcpy(comma_hyphen_list, "all");
        if ((parameter = get_args(od, "field_list"))) {
          strcpy(comma_hyphen_list, parameter);
        }
        result = parse_comma_hyphen_list_to_array(od,
          comma_hyphen_list, FIELD_LIST);
        switch (result) {
          case OUT_OF_MEMORY:
          tee_error(od, run_type, overall_line_num,
            E_OUT_OF_MEMORY, ZERO_FAILED);
          return PARSE_INPUT_ERROR;

          case INVALID_LIST_RANGE:
          tee_error(od, run_type, overall_line_num,
            E_LIST_PARSING, "fields", "ZERO",
            ZERO_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        result = set(od, FIELD_LIST, OPERATE_BIT, 1, SILENT);
        if (result) {
          tee_error(od, run_type, overall_line_num,
            E_ALLOWED_FIELD_RANGE,
            od->field_num, ZERO_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        /*
        perform zero operation
        */
        ++command;
        tee_printf(od, M_TOOL_INVOKE, nesting, command, "ZERO", line_orig);
        tee_flush(od);
        result = zero(od, type, level);
        gettimeofday(&end, NULL);
        elapsed_time(od, &start, &end);
        switch (result) {
          case Y_VAR_LOW_SD:
          tee_error(od, run_type, overall_line_num,
            E_Y_VAR_LOW_SD, ZERO_FAILED);
          return PARSE_INPUT_ERROR;

          case OUT_OF_MEMORY:
          tee_error(od, run_type, overall_line_num,
            E_OUT_OF_MEMORY, ZERO_FAILED);
          return PARSE_INPUT_ERROR;

          default:
          tee_printf(od, M_TOOL_SUCCESS, nesting, command, "ZERO");
          tee_flush(od);
        }
      }
    }
    else if (!strcasecmp(arg->me[0], "exclude")) {
      gettimeofday(&start, NULL);
      /*
      user wants to exclude from selected fields
      those values which, in a reference field, are
      outside cutoff limits; values which are outside
      the limits in any of the objects or in matching
      objects can be excluded
      */
      type = MATCHING_OBJECTS;
      if ((parameter = get_args(od, "type"))) {
        if (!strncasecmp(parameter, "any", 3)) {
          type = ANY_OBJECT;
        }
        else if (strncasecmp(parameter, "match", 5)) {
          tee_error(od, run_type, overall_line_num,
            "Allowed values for type are \"ANY\", "
            "and \"MATCH\".\n%s", EXCLUDE_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
      }
      if (!(parameter = get_args(od, "ref_field"))) {
        tee_error(od, run_type, overall_line_num,
          "Please specify the reference field "
          "whose values with respect to cutoff "
          "levels should be checked.\n%s",
          EXCLUDE_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      sscanf(parameter, "%d", &ref_field);
      if ((ref_field < 1) || (ref_field > od->field_num)) {
        tee_error(od, run_type, overall_line_num,
          E_ALLOWED_FIELD_RANGE, od->field_num, EXCLUDE_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      if (!(run_type & DRY_RUN)) {
        if (!(od->field_num)) {
          tee_error(od, run_type, overall_line_num,
            E_NO_FIELDS_PRESENT, EXCLUDE_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        strcpy(comma_hyphen_list, "all");
        if ((parameter = get_args(od, "field_list"))) {
          strcpy(comma_hyphen_list, parameter);
        }
        result = parse_comma_hyphen_list_to_array(od,
          comma_hyphen_list, FIELD_LIST);
        switch (result) {
          case OUT_OF_MEMORY:
          tee_error(od, run_type, overall_line_num,
            E_OUT_OF_MEMORY, EXCLUDE_FAILED);
          return PARSE_INPUT_ERROR;

          case INVALID_LIST_RANGE:
          tee_error(od, run_type, overall_line_num,
            E_LIST_PARSING, "fields", "EXCLUDE", EXCLUDE_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        result = set(od, FIELD_LIST, OPERATE_BIT, 1, SILENT);
        if (result) {
          tee_error(od, run_type, overall_line_num, E_ALLOWED_FIELD_RANGE,
            od->field_num, EXCLUDE_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        /*
        perform exclude operation
        */
        ++command;
        tee_printf(od, M_TOOL_INVOKE, nesting, command, "EXCLUDE", line_orig);
        tee_flush(od);
        result = exclude(od, type, ref_field - 1);
        gettimeofday(&end, NULL);
        elapsed_time(od, &start, &end);
        switch (result) {
          case Y_VAR_LOW_SD:
          tee_error(od, run_type, overall_line_num,
            E_Y_VAR_LOW_SD, EXCLUDE_FAILED);
          return PARSE_INPUT_ERROR;

          case OUT_OF_MEMORY:
          tee_error(od, run_type, overall_line_num,
            E_OUT_OF_MEMORY, EXCLUDE_FAILED);
          return PARSE_INPUT_ERROR;

          default:
          tee_printf(od, M_TOOL_SUCCESS, nesting, command, "EXCLUDE");
          tee_flush(od);
        }
      }
    }
    else if (!strcasecmp(arg->me[0], "sdcut")) {
      gettimeofday(&start, NULL);
      /*
      user wants to remove variables
      whose standard deviation lies below
      a certain threshold
      */
      /*
      perform sdcut operation
      */
      if (!(parameter = get_args(od, "level"))) {
        tee_error(od, run_type, overall_line_num,
          "Please specify the SD threshold below which "
          "SDCUT should operate.\n%s",
          SDCUT_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      sscanf(parameter, "%lf", &level);
      if (!(run_type & DRY_RUN)) {
        if (!(od->field_num)) {
          tee_error(od, run_type, overall_line_num,
            E_NO_FIELDS_PRESENT, SDCUT_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        strcpy(comma_hyphen_list, "all");
        if ((parameter = get_args(od, "field_list"))) {
          strcpy(comma_hyphen_list, parameter);
        }
        result = parse_comma_hyphen_list_to_array
          (od, comma_hyphen_list, FIELD_LIST);
        switch (result) {
          case OUT_OF_MEMORY:
          tee_error(od, run_type, overall_line_num,
            E_OUT_OF_MEMORY, SDCUT_FAILED);
          return PARSE_INPUT_ERROR;

          case INVALID_LIST_RANGE:
          tee_error(od, run_type, overall_line_num,
            E_LIST_PARSING, "fields", "SDCUT",
            SDCUT_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        result = set(od, FIELD_LIST, OPERATE_BIT, 1, SILENT);
        if (result) {
          tee_error(od, run_type, overall_line_num, E_ALLOWED_FIELD_RANGE,
            od->field_num, SDCUT_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        ++command;
        tee_printf(od, M_TOOL_INVOKE, nesting, command, "SDCUT", line_orig);
        tee_flush(od);
        result = sdcut(od, level);
        gettimeofday(&end, NULL);
        elapsed_time(od, &start, &end);
        switch (result) {
          case CANNOT_READ_TEMP_FILE:
          tee_error(od, run_type, overall_line_num,
            E_ERROR_IN_READING_TEMP_FILE, "TEMP_FIELD", SDCUT_FAILED);
          return PARSE_INPUT_ERROR;

          case Y_VAR_LOW_SD:
          tee_error(od, run_type, overall_line_num,
            E_Y_VAR_LOW_SD, SDCUT_FAILED);
          return PARSE_INPUT_ERROR;

          default:
          tee_printf(od, M_TOOL_SUCCESS, nesting, command, "SDCUT");
          tee_flush(od);
        }
      }
    }
    else if (!strcasecmp(arg->me[0], "nlevel")) {
      gettimeofday(&start, NULL);
      /*
      user wants to remove n-level variables from
      all blocks of x variables (all fields) or a subset of them
      */
      strcpy(comma_hyphen_list, "all");
      if ((parameter = get_args(od, "level"))) {
        strcpy(comma_hyphen_list, parameter);
      }
      result = parse_comma_hyphen_list_to_array
        (od, comma_hyphen_list, NLEVEL_LIST);
      switch (result) {
        case OUT_OF_MEMORY:
        tee_error(od, run_type, overall_line_num,
          E_OUT_OF_MEMORY, NLEVEL_FAILED);
        return PARSE_INPUT_ERROR;

        case INVALID_LIST_RANGE:
        tee_error(od, run_type, overall_line_num,
          "Error while parsing the range of n-level variables "
          "which NLEVEL should remove.\n%s",
          NLEVEL_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      if (!(run_type & DRY_RUN)) {
        if (!(od->field_num)) {
          tee_error(od, run_type, overall_line_num,
            E_NO_FIELDS_PRESENT, NLEVEL_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        strcpy(comma_hyphen_list, "all");
        if ((parameter = get_args(od, "field_list"))) {
          strcpy(comma_hyphen_list, parameter);
        }
        result = parse_comma_hyphen_list_to_array
          (od, comma_hyphen_list, FIELD_LIST);
        switch (result) {
          case OUT_OF_MEMORY:
          tee_error(od, run_type, overall_line_num,
            E_OUT_OF_MEMORY, NLEVEL_FAILED);
          return PARSE_INPUT_ERROR;

          case INVALID_LIST_RANGE:
          tee_error(od, run_type, overall_line_num,
            E_LIST_PARSING, "fields", "NLEVEL", NLEVEL_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        result = set(od, FIELD_LIST, OPERATE_BIT, 1, SILENT);
        if (result) {
          tee_error(od, run_type, overall_line_num, E_ALLOWED_FIELD_RANGE,
            od->field_num, NLEVEL_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        /*
        perform nlevel operation
        */
        ++command;
        tee_printf(od, M_TOOL_INVOKE, nesting, command, "NLEVEL", line_orig);
        tee_flush(od);
        result = nlevel(od);
        gettimeofday(&end, NULL);
        elapsed_time(od, &start, &end);
        switch (result) {
          case PARSE_NLEVEL_TYPE_ERROR:
          tee_error(od, run_type, overall_line_num,
            "The variable levels to be removed should lie "
            "in the 2 - 4 range.\n%s",
            NLEVEL_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;

          case OUT_OF_MEMORY:
          tee_error(od, run_type, overall_line_num,
            E_OUT_OF_MEMORY, NLEVEL_FAILED);
          return PARSE_INPUT_ERROR;

          case Y_VAR_LOW_SD:
          tee_error(od, run_type, overall_line_num,
            E_Y_VAR_LOW_SD, NLEVEL_FAILED);
          return PARSE_INPUT_ERROR;

          default:
          tee_printf(od, M_TOOL_SUCCESS, nesting, command, "NLEVEL");
          tee_flush(od);
        }
      }
      else {
        od->valid |= (TWO_LEVEL_BIT
          | THREE_LEVEL_BIT | FOUR_LEVEL_BIT);
      }
    }
    else if (!strcasecmp(arg->me[0], "pca")) {
      if (!(run_type & DRY_RUN)) {
        gettimeofday(&start, NULL);
        /*
        user wants to accomplish
        PCA analysis
        */
        if (!(od->field_num)) {
          tee_error(od, run_type, overall_line_num,
            E_NO_FIELDS_PRESENT, PCA_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        pc = 5;
        if ((parameter = get_args(od, "pc"))) {
          sscanf(parameter, "%d", &pc);
        }
        if (pc > od->overall_active_x_vars) {
          tee_error(od, run_type, overall_line_num,
            E_TOO_FEW_MANY_FOR_AVAILABLE_DATA, "many PCs",
            PCA_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        if (pc < 1) {
          tee_error(od, run_type, overall_line_num,
            E_AT_LEAST_ONE_PC, PCA_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        /*
        perform pca operation
        */
        result = calc_active_vars(od, CV_MODEL);
        switch (result) {
          case CANNOT_READ_TEMP_FILE:
          tee_error(od, run_type, overall_line_num,
            E_ERROR_IN_READING_TEMP_FILE, "TEMP_FIELD", PCA_FAILED);
          return PARSE_INPUT_ERROR;

          case Y_VAR_LOW_SD:
          tee_error(od, run_type, overall_line_num,
            E_Y_VAR_LOW_SD, PCA_FAILED);
          return PARSE_INPUT_ERROR;
        }
        result = fill_x_matrix_pca(od);
        switch (result) {
          case CANNOT_READ_TEMP_FILE:
          tee_error(od, run_type, overall_line_num,
            E_ERROR_IN_READING_TEMP_FILE, "TEMP_FIELD", PCA_FAILED);
          return PARSE_INPUT_ERROR;

          case CANNOT_WRITE_TEMP_FILE:
          tee_error(od, run_type, overall_line_num,
            E_ERROR_IN_WRITING_TEMP_FILE, "TEMP_FIELD", PCA_FAILED);
          return PARSE_INPUT_ERROR;

          case OUT_OF_MEMORY:
          tee_error(od, run_type, overall_line_num,
            E_OUT_OF_MEMORY, PCA_FAILED);
          return PARSE_INPUT_ERROR;
        }
        trim_mean_center_x_matrix_pca(od);
        result = open_temp_file(od, od->file[TEMP_PCA_SCORES], "pca_scores");
        if (result) {
          tee_error(od, run_type, overall_line_num,
            E_TEMP_FILE_CANNOT_BE_OPENED_FOR_WRITING,
            od->file[TEMP_PCA_SCORES]->name, PCA_FAILED);
          return PARSE_INPUT_ERROR;
        }
        result = open_temp_file(od, od->file[TEMP_PCA_LOADINGS], "pca_loadings");
        if (result) {
          tee_error(od, run_type, overall_line_num,
            E_TEMP_FILE_CANNOT_BE_OPENED_FOR_WRITING,
            od->file[TEMP_PCA_LOADINGS]->name, PCA_FAILED);
          return PARSE_INPUT_ERROR;
        }
        od->file[ASCII_IN]->name[0] = '\0';
        if ((parameter = get_args(od, "file"))) {
          strcpy(od->file[ASCII_IN]->name, parameter);
          if (!(od->file[ASCII_IN]->handle =
            fopen(od->file[ASCII_IN]->name, "wb+"))) {
            tee_error(od, run_type, overall_line_num,
              E_FILE_CANNOT_BE_OPENED_FOR_WRITING,
              od->file[ASCII_IN]->name, PCA_FAILED);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
        }
        ++command;
        tee_printf(od, M_TOOL_INVOKE, nesting, command, "PCA", line_orig);
        tee_flush(od);
        result = pca(od, pc);
        gettimeofday(&end, NULL);
        elapsed_time(od, &start, &end);
        switch (result) {
          case OUT_OF_MEMORY:
          tee_error(od, run_type, overall_line_num,
            E_OUT_OF_MEMORY, PCA_FAILED);
          return PARSE_INPUT_ERROR;

          case CANNOT_READ_TEMP_FILE:
          tee_error(od, run_type, overall_line_num,
            E_ERROR_IN_READING_TEMP_FILE, "TEMP_FIELD", PCA_FAILED);
          return PARSE_INPUT_ERROR;

          case CANNOT_WRITE_TEMP_FILE:
          tee_error(od, run_type, overall_line_num,
            E_ERROR_IN_WRITING_TEMP_FILE, "TEMP_FIELD", PCA_FAILED);
          return PARSE_INPUT_ERROR;

          default:
          tee_printf(od, M_TOOL_SUCCESS, nesting, command, "PCA");
          tee_flush(od);
          od->valid |= PCA_BIT;
        }
        if (od->file[TEMP_PCA_SCORES]->handle) {
          fclose(od->file[TEMP_PCA_SCORES]->handle);
          od->file[TEMP_PCA_SCORES]->handle = NULL;
        }
        if (od->file[TEMP_PCA_LOADINGS]->handle) {
          fclose(od->file[TEMP_PCA_LOADINGS]->handle);
          od->file[TEMP_PCA_LOADINGS]->handle = NULL;
        }
      }
      else {
        od->valid |= PCA_BIT;
      }
    }
    else if (!strcasecmp(arg->me[0], "set")) {
      gettimeofday(&start, NULL);
      result = parse_synonym_lists(od, "SET", SET_FAILED,
        (1 << OBJECT_LIST) | (1 << ID_LIST)
        | (1 << STRUCT_LIST) | (1 << FIELD_LIST),
        &list_type, 0, run_type, overall_line_num);
      switch (result) {
        case PARSE_INPUT_RECOVERABLE_ERROR:
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
        
        case PARSE_INPUT_ERROR:
        return PARSE_INPUT_ERROR;
      }
      if (list_type & (1 << FIELD_LIST)) {
        if ((parameter = get_args(od, "attribute"))) {
          if (!strncasecmp(parameter, "included", 8)) {
            attr = ACTIVE_BIT;
            state = 1;
          }
          else if (!strncasecmp(parameter, "excluded", 8)) {
            attr = ACTIVE_BIT;
            state = 0;
          }
          else {
            tee_error(od, run_type, overall_line_num,
              "Only \"INCLUDED\" and \"EXCLUDED\" "
              "attributes are allowed.\n%s",
              SET_FAILED);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
        }
        else {
          tee_error(od, run_type, overall_line_num,
            "Please specify an attribute which "
            "should be set.\n%s", SET_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
      }
      else {
        if ((parameter = get_args(od, "attribute"))) {
          if (!strncasecmp(parameter, "training", 8)) {
            attr = ACTIVE_BIT;
            state = 1;
          }
          else if (!strncasecmp(parameter, "excluded", 8)) {
            attr = ACTIVE_BIT;
            state = 0;
          }
          else if (!strncasecmp(parameter, "test", 4)) {
            attr = PREDICT_BIT;
            state = 1;
          }
          else {
            tee_error(od, run_type, overall_line_num,
              "Only \"TRAINING\", \"EXCLUDED\", "
              "\"TEST\" parameters are allowed.\n%s",
              SET_FAILED);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
        }
        else {
          tee_error(od, run_type, overall_line_num,
            "Please specify an attribute which should be set.\n%s",
            SET_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
      }
      if (!(run_type & DRY_RUN)) {
        ++command;
        tee_printf(od, M_TOOL_INVOKE, nesting, command, "SET", line_orig);
        tee_flush(od);
        result = set(od, intlog2(list_type), attr, state, VERBOSE_BIT);
        gettimeofday(&end, NULL);
        elapsed_time(od, &start, &end);
        switch (result) {
          case CANNOT_READ_TEMP_FILE:
          tee_error(od, run_type, overall_line_num,
            E_ERROR_IN_READING_TEMP_FILE, "TEMP_FIELD", SET_FAILED);
          return PARSE_INPUT_ERROR;

          case Y_VAR_LOW_SD:
          tee_error(od, run_type, overall_line_num,
            E_Y_VAR_LOW_SD, SET_FAILED);
          return PARSE_INPUT_ERROR;

          case INVALID_LIST_RANGE:
          tee_error(od, run_type, overall_line_num,
            "The range you selected is not valid.\n%s", SET_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;

          default:
          update_pymol(od);
          update_jmol(od);
          tee_printf(od, M_TOOL_SUCCESS, nesting, command, "SET");
          tee_flush(od);
        }
      }
    }
    else if (!strcasecmp(arg->me[0], "scale_object")) {
      weight = 1.0;
      options = 0;
      if ((parameter = get_args(od, "weight"))) {
        if (!strncasecmp(parameter, "rand", 4)) {
          options = RANDOM_WEIGHTS;
        }
        else if (!strncasecmp(parameter, "even", 4)) {
          options = EVEN_WEIGHTS;
        }
        else {
          sscanf(parameter, "%lf", &weight);
          if (weight < 0.0) {
            tee_error(od, run_type, overall_line_num, E_POSITIVE_NUMBER,
              "weight parameter", SCALE_OBJECT_FAILED);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
        }
      }
      result = parse_synonym_lists(od, "SCALE_OBJECT", SCALE_OBJECT_FAILED,
        (1 << OBJECT_LIST) | (1 << ID_LIST)
        | (1 << STRUCT_LIST) | (1 << FROM_FILE),
        &list_type, OBJECT_LIST, run_type, overall_line_num);
      switch (result) {
        case PARSE_INPUT_RECOVERABLE_ERROR:
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
        
        case PARSE_INPUT_ERROR:
        return PARSE_INPUT_ERROR;
      }
      if (!(run_type & DRY_RUN)) {
        if (!(od->field_num)) {
          tee_error(od, run_type, overall_line_num,
            E_NO_FIELDS_PRESENT, SCALE_OBJECT_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        gettimeofday(&start, NULL);
        ++command;
        tee_printf(od, M_TOOL_INVOKE, nesting, command, "SCALE_OBJECT", line_orig);
        tee_flush(od);
        result = set_object_weight(od, weight, list_type, options);
        if (od->file[ASCII_IN]->handle) {
          fclose(od->file[ASCII_IN]->handle);
          od->file[ASCII_IN]->handle = NULL;
        }
        switch (result) {
          case PREMATURE_EOF:
          tee_error(od, run_type, overall_line_num, 
            E_FILE_CORRUPTED_OR_IN_WRONG_FORMAT, "ASCII",
            od->file[ASCII_IN]->name, SCALE_OBJECT_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          break;

          case INVALID_LIST_RANGE:
          tee_error(od, run_type, overall_line_num,
            "At least one of the object/ID numbers specified "
            "in the file %s is out of range.\n%s",
            od->file[ASCII_IN]->name, SCALE_OBJECT_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          break;

          case WRONG_DATA_FORMAT:
          tee_error(od, run_type, overall_line_num,
            "Weights specified in the file %s "
            "must be positive numbers.\n%s",
            od->file[ASCII_IN]->name, SCALE_OBJECT_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          break;
        }
        if (result) {
          continue;
        }
        update_field_object_attr(od, VERBOSE_BIT);
        gettimeofday(&end, NULL);
        elapsed_time(od, &start, &end);
        tee_printf(od, M_TOOL_SUCCESS, nesting, command, "SCALE_OBJECT");
        tee_flush(od);
      }
    }
    else if (!strcasecmp(arg->me[0], "scale_x_vars")) {
      gettimeofday(&start, NULL);
      type = CUSTOM_SCALE;
      weight = 1.0;
      if ((parameter = get_args(od, "type"))) {
        if (!strncasecmp(parameter, "custom", 6)) {
          if ((parameter = get_args(od, "weight"))) {
            sscanf(parameter, "%lf", &weight);
          }
        }
        else if (!strncasecmp(parameter, "auto", 4)) {
          type = AUTO_SCALE;
        }
        else if (!strncasecmp(parameter, "buw", 3)) {
          type = BUW_SCALE;
        }
        else {
          tee_error(od, run_type, overall_line_num,
            E_ONLY_CUSTOM_AUTO_BUW_ALLOWED, SCALE_X_VARS_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
      }
      if (!(run_type & DRY_RUN)) {
        if (!(od->field_num)) {
          tee_error(od, run_type, overall_line_num,
            E_NO_FIELDS_PRESENT, SCALE_X_VARS_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        strcpy(comma_hyphen_list, "all");
        if ((parameter = get_args(od, "field_list"))) {
          strcpy(comma_hyphen_list, parameter);
        }
        result = parse_comma_hyphen_list_to_array(od,
          comma_hyphen_list, FIELD_LIST);
        switch (result) {
          case OUT_OF_MEMORY:
          tee_error(od, run_type, overall_line_num,
            E_OUT_OF_MEMORY, SCALE_X_VARS_FAILED);
          return PARSE_INPUT_ERROR;

          case INVALID_LIST_RANGE:
          tee_error(od, run_type, overall_line_num, E_LIST_PARSING, "fields",
            "SCALE_X_VARS", SCALE_X_VARS_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        gettimeofday(&start, NULL);
        result = set(od, FIELD_LIST, OPERATE_BIT, 1, SILENT);
        if (result) {
          tee_error(od, run_type, overall_line_num, E_ALLOWED_FIELD_RANGE,
            od->field_num, SCALE_X_VARS_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        ++command;
        tee_printf(od, M_TOOL_INVOKE, nesting, command, "SCALE_X_VARS", line_orig);
        tee_flush(od);
        switch (type) {
          case CUSTOM_SCALE:
          set_field_weight(od, weight);
          break;

          case BUW_SCALE:
          result = x_var_buw(od);
          switch (result) {
            case OUT_OF_MEMORY:
            tee_error(od, run_type, overall_line_num,
              E_OUT_OF_MEMORY, SCALE_X_VARS_FAILED);
            return PARSE_INPUT_ERROR;
          }
          break;

          case AUTO_SCALE:
          result = autoscale_field(od);
          switch (result) {
            case OUT_OF_MEMORY:
            tee_error(od, run_type, overall_line_num,
              E_OUT_OF_MEMORY, SCALE_X_VARS_FAILED);
            return PARSE_INPUT_ERROR;
          }
          break;
        }
        gettimeofday(&end, NULL);
        elapsed_time(od, &start, &end);
        tee_printf(od, M_TOOL_SUCCESS, nesting, command, "SCALE_X_VARS");
        tee_flush(od);
      }
    }
    else if (!strcasecmp(arg->me[0], "scale_y_vars")) {
      gettimeofday(&start, NULL);
      type = CUSTOM_SCALE;
      weight = 1.0;
      if ((parameter = get_args(od, "type"))) {
        if (!strncasecmp(parameter, "custom", 6)) {
          if ((parameter = get_args(od, "weight"))) {
            sscanf(parameter, "%lf", &weight);
          }
        }
        else if (!strncasecmp(parameter, "auto", 4)) {
          type = AUTO_SCALE;
        }
        else if (!strncasecmp(parameter, "buw", 3)) {
          type = BUW_SCALE;
        }
        else {
          tee_error(od, run_type, overall_line_num,
            E_ONLY_CUSTOM_AUTO_BUW_ALLOWED, SCALE_Y_VARS_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
      }
      if (!(run_type & DRY_RUN)) {
        if (!(od->y_vars)) {
          tee_error(od, run_type, overall_line_num,
            E_NO_Y_VARS_PRESENT, SCALE_Y_VARS_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        strcpy(comma_hyphen_list, "all");
        if ((parameter = get_args(od, "y_var_list"))) {
          strcpy(comma_hyphen_list, parameter);
        }
        result = parse_comma_hyphen_list_to_array(od,
          comma_hyphen_list, Y_VAR_LIST);
        switch (result) {
          case OUT_OF_MEMORY:
          tee_error(od, run_type, overall_line_num,
            E_OUT_OF_MEMORY, SCALE_Y_VARS_FAILED);
          return PARSE_INPUT_ERROR;

          case INVALID_LIST_RANGE:
          tee_error(od, run_type, overall_line_num,
            E_LIST_PARSING, "Y variables",
            "SCALE_Y_VARS", SCALE_Y_VARS_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        gettimeofday(&start, NULL);
        result = set(od, Y_VAR_LIST, OPERATE_BIT, 1, SILENT);
        if (result) {
          tee_error(od, run_type, overall_line_num,
            E_Y_VAR_ALLOWED_RANGE,
            od->y_vars, SCALE_Y_VARS_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        ++command;
        tee_printf(od, M_TOOL_INVOKE, nesting, command, "SCALE_Y_VARS", line_orig);
        tee_flush(od);
        switch (type) {
          case CUSTOM_SCALE:
          set_y_var_weight(od, weight);
          break;

          case BUW_SCALE:
          result = y_var_buw(od);
          if (result) {
            tee_error(od, run_type, overall_line_num,
              E_OUT_OF_MEMORY, SCALE_Y_VARS_FAILED);
            return PARSE_INPUT_ERROR;
          }
          break;

          case AUTO_SCALE:
          autoscale_y_var(od);
          break;
        }
        gettimeofday(&end, NULL);
        elapsed_time(od, &start, &end);
        tee_printf(od, M_TOOL_SUCCESS, nesting, command, "SCALE_Y_VARS");
        tee_flush(od);
      }
    }
    else if (!strcasecmp(arg->me[0], "env")) {
      if (!(run_type & DRY_RUN)) {
        ++command;
        tee_printf(od, M_TOOL_INVOKE, nesting, command, "ENV", line_orig);
        tee_flush(od);
      }
      if ((parameter = get_args(od, "n_cpus"))) {
        determine_best_cpu_number(od, "all");
        determine_best_cpu_number(od, parameter);
        tee_printf(od, M_NUMBER_OF_CPUS, PACKAGE_NAME, od->n_proc);
      }
      else if ((parameter = get_args(od, "random_seed"))) {
        if (!(run_type & DRY_RUN)) {
          sscanf(parameter, "%ld", &random_seed);
          set_random_seed(od, (unsigned long)absval(random_seed));
          tee_printf(od, "The random seed has been set to %ld.\n\n",
            random_seed);
        }
      }
      else if ((parameter = get_args(od, "temp_dir"))) {
        if (dexist(parameter)) {
          strcpy(od->temp_dir, parameter);
          absolute_path(od->temp_dir);
          if (!(run_type & DRY_RUN)) {
            tee_printf(od, "The temporary directory has been set to %s.\n\n",
              od->temp_dir);
          }
        }
        else {
          tee_error(od, run_type, overall_line_num,
            E_DIR_NOT_EXISTING, parameter, ENV_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
      }
      else if ((parameter = get_args(od, "qm_engine"))) {
        memset(od->field.qm_exe, 0, BUF_LEN);
        memset(od->field.qm_exe_path, 0, BUF_LEN);
        found = 0;
        if (parameter[0]) {
          found = fexist(parameter);
          if (!found) {
            found = is_in_path(parameter, od->field.qm_exe_path);
          }
          else {
            strcpy(od->field.qm_exe_path, parameter);
          }
        }
        if (found) {
          absolute_path(od->field.qm_exe_path);
          strcpy(od->field.qm_exe, get_basename(od->field.qm_exe_path));
          get_dirname(od->field.qm_exe_path);
          if (!(run_type & DRY_RUN)) {
            tee_printf(od, "%s was chosen as QM engine.\n\n", od->field.qm_exe);
          }
        }
        else {
          tee_error(od, run_type, overall_line_num,
            "The QM engine %s could not be found.\n%s\n",
            parameter, ENV_FAILED);
          memset(od->field.qm_exe_path, 0, BUF_LEN);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
      }
      else if ((parameter = get_args(od, "babel_path"))) {
        memset(od->field.babel_exe_path, 0, BUF_LEN);
        found = 0;
        result = 1;
        if (parameter[0]) {
          strncpy(od->field.babel_exe_path, parameter, BUF_LEN - 1);
          absolute_path(od->field.babel_exe_path);
          if (dexist(od->field.babel_exe_path)) {
            sprintf(buffer, "%s%c%s", od->field.babel_exe_path, SEPARATOR, BABEL_EXE);
            if (fexist(buffer)) {
              sprintf(buffer, "%s%c%s", od->field.babel_exe_path, SEPARATOR, OBENERGY_EXE);
              if (fexist(buffer)) {
                found = 1;
              }
            }
          }
        }
        if (!(run_type & DRY_RUN)) {
          if (found) {
            result = check_babel(od, od->field.babel_exe_path);
            switch (result) {
              case OUT_OF_MEMORY:
              tee_error(od, run_type, overall_line_num,
                E_OUT_OF_MEMORY, ENV_FAILED);
              return PARSE_INPUT_ERROR;

              case CANNOT_READ_TEMP_FILE:
              tee_printf(od, E_FILE_CANNOT_BE_OPENED_FOR_READING,
                od->file[TEMP_LOG]->name);
              break;

              case CANNOT_WRITE_TEMP_FILE:
              tee_printf(od, E_FILE_CANNOT_BE_OPENED_FOR_WRITING,
                od->file[TEMP_LOG]->name);
              break;

              case BABEL_NOT_WORKING:
              tee_printf(od, "A simple OpenBabel test failed with the following error:\n");
              if ((od->file[TEMP_LOG]->handle = fopen(od->file[TEMP_LOG]->name, "rb"))) {
                while (fgets(buffer, BUF_LEN, od->file[TEMP_LOG]->handle)) {
                  buffer[BUF_LEN - 1] = '\0';
                  tee_printf(od, "%s", buffer);
                }
                fclose(od->file[TEMP_LOG]->handle);
                od->file[TEMP_LOG]->handle = NULL;
              }
              break;

              case BABEL_PLUGINS_NOT_FOUND:
              tee_printf(od, "OpenBabel data/plugins could not be found.\n"
                "This issue can be fixed setting the "BABEL_DATADIR_ENV" and "
                BABEL_LIBDIR_ENV" environment "
                "variables to suitable values.\n");
              break;
            }
          }
          if ((!found) || (result)) {
            if (found) {
              tee_printf(od, "Working ");
            }
            tee_error(od, run_type, overall_line_num,
              "OpenBabel binaries could not be found.\n", ENV_FAILED);
            memset(od->field.babel_exe_path, 0, BUF_LEN);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
          tee_printf(od, M_EXECUTABLE_PATH, "OpenBabel",
            od->field.babel_exe_path);
        }
      }
      else if ((parameter = get_args(od, "md_grid_path"))) {
        memset(od->field.md_grid_exe_path, 0, BUF_LEN);
        found = 0;
        if (parameter[0]) {
          strncpy(od->field.md_grid_exe_path, parameter, BUF_LEN - 1);
          absolute_path(od->field.md_grid_exe_path);
          if (dexist(od->field.md_grid_exe_path)) {
            sprintf(buffer, "%s%c%s", od->field.md_grid_exe_path, SEPARATOR, GRID_EXE);
            if (fexist(buffer)) {
              sprintf(buffer, "%s%c%s", od->field.md_grid_exe_path, SEPARATOR, GRIN_EXE);
              if (fexist(buffer)) {
                found = 1;
              }
            }
          }
        }
        if (found) {
          if (!(run_type & DRY_RUN)) {
            tee_printf(od, M_EXECUTABLE_PATH, "MD GRID",
              od->field.md_grid_exe_path);
          }
        }
        else {
          tee_error(od, run_type, overall_line_num,
            "MD GRID binaries could not be found.\n%s\n", ENV_FAILED);
          memset(od->field.md_grid_exe_path, 0, BUF_LEN);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
      }
      else if ((parameter = get_args(od, "pymol"))) {
        memset(od->pymol.pymol_exe, 0, BUF_LEN);
        od->pymol.use_pymol = 0;
        found = 0;
        if (parameter[0]) {
          found = fexist(parameter);
          if (!found) {
            found = is_in_path(parameter, od->pymol.pymol_exe);
          }
          else {
            strcpy(od->pymol.pymol_exe, parameter);
          }
        }
        if (parameter[0] && (!found)) {
          tee_error(od, run_type, overall_line_num,
            "The PyMOL binary %s could not be found.\n%s\n",
            parameter, ENV_FAILED);
          memset(od->pymol.pymol_exe, 0, BUF_LEN);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        if (!found) {
          if (!(run_type & DRY_RUN)) {
            tee_printf(od, "PyMOL will not be used.\n\n");
          }
        }
        else {
          od->pymol.use_pymol = 1;
          absolute_path(od->pymol.pymol_exe);
          tee_printf(od, M_EXECUTABLE, "PyMOL",
            od->pymol.pymol_exe);
          strcat(od->pymol.pymol_exe, " " PYMOL_ARGS);
        }
      }
      else if ((parameter = get_args(od, "jmol"))) {
        memset(od->jmol.jmol_exe, 0, BUF_LEN);
        od->jmol.use_jmol = 0;
        found = 0;
        if (parameter[0]) {
          found = fexist(parameter);
          if (!found) {
            found = is_in_path(parameter, od->jmol.jmol_exe);
          }
          else {
            strcpy(od->jmol.jmol_exe, parameter);
          }
        }
        if (parameter[0] && (!found)) {
          tee_error(od, run_type, overall_line_num,
            "The Jmol start script %s could not be found.\n%s\n",
            parameter, ENV_FAILED);
          memset(od->jmol.jmol_exe, 0, BUF_LEN);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        if (!found) {
          if (!(run_type & DRY_RUN)) {
            tee_printf(od, "Jmol will not be used.\n\n");
          }
        }
        else {
          od->jmol.use_jmol = 1;
          absolute_path(od->jmol.jmol_exe);
          tee_printf(od, M_EXECUTABLE, "Jmol",
            od->jmol.jmol_exe);
          strcat(od->jmol.jmol_exe, " " JMOL_ARGS);
        }
      }
      else if ((parameter = get_args(od, "gnuplot"))) {
        memset(od->gnuplot.gnuplot_exe, 0, BUF_LEN);
        od->gnuplot.use_gnuplot = 0;
        found = 0;
        if (parameter[0]) {
          found = fexist(parameter);
          if (!found) {
            found = is_in_path(parameter, od->gnuplot.gnuplot_exe);
          }
          else {
            strcpy(od->gnuplot.gnuplot_exe, parameter);
          }
        }
        if (parameter[0] && (!found)) {
          tee_error(od, run_type, overall_line_num,
            "The gnuplot binary %s could not be found.\n%s\n",
            parameter, ENV_FAILED);
          memset(od->gnuplot.gnuplot_exe, 0, BUF_LEN);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        if (!found) {
          if (!(run_type & DRY_RUN)) {
            tee_printf(od, "gnuplot will not be used.\n\n");
          }
        }
        else {
          od->gnuplot.use_gnuplot = 1;
          absolute_path(od->gnuplot.gnuplot_exe);
          tee_printf(od, M_EXECUTABLE, "gnuplot",
            od->gnuplot.gnuplot_exe);
        }
      }
      else if ((parameter = get_args(od, "cs3d"))) {
        memset(od->field.cs3d_exe, 0, BUF_LEN);
        found = 0;
        if (parameter[0]) {
          found = fexist(parameter);
          if (!found) {
            found = is_in_path(parameter, od->field.cs3d_exe);
          }
          else {
            strcpy(od->field.cs3d_exe, parameter);
          }
        }
        if (parameter[0] && (!found)) {
          tee_error(od, run_type, overall_line_num,
            "The CS3D binary %s could not be found.\n%s\n",
            parameter, ENV_FAILED);
          memset(od->field.cs3d_exe, 0, BUF_LEN);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        if (!found) {
          if (!(run_type & DRY_RUN)) {
            tee_printf(od, "CS3D will not be used.\n\n");
          }
        }
        else {
          absolute_path(od->field.cs3d_exe);
          tee_printf(od, M_EXECUTABLE, "CS3D",
            od->field.cs3d_exe);
        }
      }
      else if ((parameter = get_args(od, "nice"))) {
        if (!(run_type & DRY_RUN)) {
          #ifndef WIN32
          sscanf(parameter, "%d", &nice_value);
          #else
          for (nice_value = 0; nice_value < 6; ++nice_value) {
            if (!strcasecmp(parameter, nice_name[nice_value])) {
              break;
            }
          }
          #endif
          set_nice_value(od, nice_value);
        }
      }
      else {
        tee_error(od, run_type, overall_line_num,
          "Allowed environmental variables which may be set are: "
          "\"random_seed\", \"temp_dir\", \"n_cpus\", \"nice\", "
          "\"babel_path\", \"md_grid_path\", "
          "\"qm_engine\", \"cs3d\", \"gnuplot\", "
          "\"jmol\" and \"pymol\".\n%s",
          ENV_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      if (!(run_type & DRY_RUN)) {
        tee_printf(od, M_TOOL_SUCCESS, nesting, command, "ENV");
        tee_flush(od);
      }
    }
    else if (!strcasecmp(arg->me[0], "remove_box")) {
      if (od->field_num) {
        tee_error(od, run_type, overall_line_num,
          E_CANNOT_ALTER_GRID_BOX,
          REMOVE_BOX_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      if (!(od->grid.nodes[0])) {
        tee_error(od, run_type, overall_line_num,
          E_NO_GRID_BOX_PRESENT,
          REMOVE_BOX_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      if (!(run_type & DRY_RUN)) {
        ++command;
        tee_printf(od, M_TOOL_INVOKE, nesting, command, "REMOVE_BOX", line_orig);
        tee_flush(od);
        remove_box(od);
        update_pymol(od);
        update_jmol(od);
        tee_printf(od, M_TOOL_SUCCESS, nesting, command, "REMOVE_BOX");
        tee_flush(od);
      }
    }
    else if (!strcasecmp(arg->me[0], "remove_x_vars")) {
      gettimeofday(&start, NULL);
      if (!(parameter = get_args(od, "type"))) {
        tee_error(od, run_type, overall_line_num,
          "Please specify which type of variables "
          "should be removed.\n%s",
          REMOVE_X_VARS_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      if (!strncasecmp(parameter, "nlevel", 6)) {
        strcpy(comma_hyphen_list, "2-4");
        if ((parameter = get_args(od, "level"))) {
          strcpy(comma_hyphen_list, parameter);
        }
        result = parse_comma_hyphen_list_to_array(od,
          comma_hyphen_list, NLEVEL_LIST);
        switch (result) {
          case OUT_OF_MEMORY:
          tee_error(od, run_type, overall_line_num,
            E_OUT_OF_MEMORY, REMOVE_X_VARS_FAILED);
          return PARSE_INPUT_ERROR;

          case INVALID_LIST_RANGE:
          tee_error(od, run_type, overall_line_num,
            E_LIST_PARSING, "NLEVEL variables",
            "REMOVE_X_VARS", REMOVE_X_VARS_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        if (!(run_type & DRY_RUN)) {
          if (!(od->field_num)) {
            tee_error(od, run_type, overall_line_num,
              E_NO_FIELDS_PRESENT, REMOVE_X_VARS_FAILED);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
          for (i = 0; i < od->pel.numberlist[NLEVEL_LIST]->size; ++i) {
            if (od->pel.numberlist[NLEVEL_LIST]->pe[i] == 2) {
              if (!(od->valid & TWO_LEVEL_BIT)) {
                tee_error(od, run_type, overall_line_num,
                  E_NO_VALID_SELECTION, "2-level",
                  REMOVE_X_VARS_FAILED);
                fail = !(run_type & INTERACTIVE_RUN);
                continue;
              }
            }
            else if (od->pel.numberlist[NLEVEL_LIST]->pe[i] == 3) {
              if (!(od->valid & THREE_LEVEL_BIT)) {
                tee_error(od, run_type, overall_line_num,
                  E_NO_VALID_SELECTION, "3-level",
                  REMOVE_X_VARS_FAILED);
                fail = !(run_type & INTERACTIVE_RUN);
                continue;
              }
            }
            else if (od->pel.numberlist[NLEVEL_LIST]->pe[i] == 4) {
              if (!(od->valid & FOUR_LEVEL_BIT)) {
                tee_error(od, run_type, overall_line_num,
                  E_NO_VALID_SELECTION, "4-level",
                  REMOVE_X_VARS_FAILED);
                fail = !(run_type & INTERACTIVE_RUN);
                continue;
              }
            }
          }
        }
        attr = TWO_LEVEL_BIT;
      }
      else if (!strncasecmp(parameter, "d_optimal", 9)) {
        if (!(run_type & DRY_RUN)) {
          if (!(od->valid & D_OPTIMAL_BIT)) {
            tee_error(od, run_type, overall_line_num,
              E_NO_VALID_SELECTION, "D-optimal",
              REMOVE_X_VARS_FAILED);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
        }
        attr = D_OPTIMAL_BIT;
      }
      else if (!strncasecmp(parameter, "ffdsel", 6)) {
        if (!(run_type & DRY_RUN)) {
          if (!(od->valid & FFDSEL_BIT)) {
            tee_error(od, run_type, overall_line_num,
              E_NO_VALID_SELECTION, "FFDSEL",
              REMOVE_X_VARS_FAILED);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
        }
        attr = FFDSEL_BIT;
      }
      else if (!strncasecmp(parameter, "uvepls", 6)) {
        if (!(run_type & DRY_RUN)) {
          if (!(od->valid & UVEPLS_BIT)) {
            tee_error(od, run_type, overall_line_num,
              E_NO_VALID_SELECTION, "UVEPLS",
              REMOVE_X_VARS_FAILED);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
        }
        attr = UVEPLS_BIT;
      }
      else if (!strncasecmp(parameter, "groups", 10)) {
        /*
        if not specified, default group to be removed is group zero
        */
        if (!(run_type & DRY_RUN)) {
          if (!(od->valid & SEED_BIT)) {
            tee_error(od, run_type, overall_line_num,
              E_NO_VALID_SELECTION, "SRD",
              REMOVE_X_VARS_FAILED);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
          strcpy(comma_hyphen_list, "0");
          if ((parameter = get_args(od, "group_list"))) {
            strcpy(comma_hyphen_list, parameter);
          }
          result = parse_comma_hyphen_list_to_array
            (od, comma_hyphen_list, GROUP_LIST);
          switch (result) {
            case OUT_OF_MEMORY:
            tee_error(od, run_type, overall_line_num,
              E_OUT_OF_MEMORY, REMOVE_X_VARS_FAILED);
            return PARSE_INPUT_ERROR;

            case INVALID_LIST_RANGE:
            tee_error(od, run_type, overall_line_num,
              E_LIST_PARSING, "SRD groups", "REMOVE_X_VARS",
              REMOVE_X_VARS_FAILED);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
        }
        attr = SEED_BIT;
      }
      else {
        tee_error(od, run_type, overall_line_num,
          "Only \"NLEVEL\", \"D_OPTIMAL\", \"GROUPS\", "
          "\"FFDSEL\" types are allowed.\n%s",
          REMOVE_X_VARS_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      if (!(run_type & DRY_RUN)) {
        /*
        if not specified, all fields are taken into consideration
        */
        strcpy(comma_hyphen_list, "all");
        if ((parameter = get_args(od, "field_list"))) {
          strcpy(comma_hyphen_list, parameter);
        }
        result = parse_comma_hyphen_list_to_array
          (od, comma_hyphen_list, FIELD_LIST);
        switch (result) {
          case OUT_OF_MEMORY:
          tee_error(od, run_type, overall_line_num,
            E_OUT_OF_MEMORY, REMOVE_X_VARS_FAILED);
          return PARSE_INPUT_ERROR;

          case INVALID_LIST_RANGE:
          tee_error(od, run_type, overall_line_num,
            E_LIST_PARSING, "fields", "REMOVE_X_VARS",
            REMOVE_X_VARS_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        result = set(od, FIELD_LIST, OPERATE_BIT, 1, SILENT);
        if (result) {
          tee_error(od, run_type, overall_line_num, E_ALLOWED_FIELD_RANGE,
            od->field_num, REMOVE_X_VARS_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        ++command;
        tee_printf(od, M_TOOL_INVOKE, nesting, command, "REMOVE_X_VARS", line_orig);
        tee_flush(od);
        result = remove_x_vars(od, attr);
        gettimeofday(&end, NULL);
        elapsed_time(od, &start, &end);
        switch (result) {

          case INVALID_LIST_RANGE:
          tee_error(od, run_type, overall_line_num,
            "The selected range is not valid.\n%s",
            REMOVE_X_VARS_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;

          case ONE_FIELD_AT_A_TIME:
          tee_error(od, run_type, overall_line_num,
            E_GROUPS_OTHER_THAN_ZERO, "removed",
            REMOVE_X_VARS_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;

          default:
          tee_printf(od, M_TOOL_SUCCESS, nesting, command, "REMOVE_X_VARS");
          tee_flush(od);
        }
      }
      else {
        od->valid &= SDF_BIT;
      }
    }  
    else if (!strcasecmp(arg->me[0], "remove_field")) {
      gettimeofday(&start, NULL);
      if (!(parameter = get_args(od, "field_list"))) {
        tee_error(od, run_type, overall_line_num,
          "Please specify which fields should be removed.\n%s",
          REMOVE_FIELD_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      if (!(run_type & DRY_RUN)) {
        if (!(od->field_num)) {
          tee_error(od, run_type, overall_line_num,
            E_NO_FIELDS_PRESENT, REMOVE_FIELD_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        result = parse_comma_hyphen_list_to_array
          (od, parameter, FIELD_LIST);
        switch (result) {
          case OUT_OF_MEMORY:
          tee_error(od, run_type, overall_line_num,
            E_OUT_OF_MEMORY, REMOVE_FIELD_FAILED);
          return PARSE_INPUT_ERROR;

          case INVALID_LIST_RANGE:
          tee_error(od, run_type, overall_line_num,
            E_LIST_PARSING, "fields", "REMOVE_FIELD",
            REMOVE_FIELD_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        result = set(od, FIELD_LIST, OPERATE_BIT, 1, SILENT);
        if (result) {
          tee_error(od, run_type, overall_line_num,
            E_ALLOWED_FIELD_RANGE,
            od->field_num, REMOVE_FIELD_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        ++command;
        tee_printf(od, M_TOOL_INVOKE, nesting, command, "REMOVE_FIELD", line_orig);
        tee_flush(od);
        result = remove_field(od);
        gettimeofday(&end, NULL);
        elapsed_time(od, &start, &end);
        switch (result) {
          case OUT_OF_MEMORY:
          tee_error(od, run_type, overall_line_num,
            E_OUT_OF_MEMORY, REMOVE_FIELD_FAILED);
          return PARSE_INPUT_ERROR;

          case CANNOT_READ_TEMP_FILE:
          case PREMATURE_DAT_EOF:
          tee_error(od, run_type, overall_line_num,
            E_ERROR_IN_READING_TEMP_FILE, "TEMP_FIELD", REMOVE_FIELD_FAILED);
          return PARSE_INPUT_ERROR;

          case CANNOT_WRITE_TEMP_FILE:
          tee_error(od, run_type, overall_line_num,
            E_ERROR_IN_WRITING_TEMP_FILE, "TEMP_FIELD", REMOVE_FIELD_FAILED);
          return PARSE_INPUT_ERROR;

          default:
          if (!(od->object_num)) {
            update_pymol(od);
            update_jmol(od);
            tee_printf(od, "Since all fields were removed, all grid "
              "information has been removed as well.\n");
          }
          tee_printf(od, M_TOOL_SUCCESS, nesting, command, "REMOVE_FIELD");
          tee_flush(od);
        }
      }
      else {
        od->valid &= SDF_BIT;
      }
    }  
    else if (!strcasecmp(arg->me[0], "remove_object")) {
      gettimeofday(&start, NULL);
      result = parse_synonym_lists(od, "REMOVE_OBJECT", REMOVE_OBJECT_FAILED,
        (1 << OBJECT_LIST) | (1 << ID_LIST) | (1 << STRUCT_LIST), &list_type,
        0, run_type, overall_line_num);
      switch (result) {
        case PARSE_INPUT_RECOVERABLE_ERROR:
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
        
        case PARSE_INPUT_ERROR:
        return PARSE_INPUT_ERROR;
      }
      if (!(run_type & DRY_RUN)) {
        ++command;
        tee_printf(od, M_TOOL_INVOKE, nesting, command, "REMOVE_OBJECT", line_orig);
        tee_flush(od);
        result = remove_object(od);
        gettimeofday(&end, NULL);
        elapsed_time(od, &start, &end);
        switch (result) {
          case OUT_OF_MEMORY:
          tee_error(od, run_type, overall_line_num,
            E_OUT_OF_MEMORY, REMOVE_OBJECT_FAILED);
          return PARSE_INPUT_ERROR;

          case CANNOT_READ_TEMP_FILE:
          case PREMATURE_DAT_EOF:
          tee_error(od, run_type, overall_line_num,
            E_ERROR_IN_READING_TEMP_FILE, "TEMP_FIELD", REMOVE_OBJECT_FAILED);
          return PARSE_INPUT_ERROR;

          case CANNOT_WRITE_TEMP_FILE:
          tee_error(od, run_type, overall_line_num,
            E_ERROR_IN_WRITING_TEMP_FILE, "TEMP_FIELD", REMOVE_OBJECT_FAILED);
          return PARSE_INPUT_ERROR;

          default:
          if (!(od->grid.object_num)) {
            tee_printf(od, "Since all objects were removed, all grid "
              "information has been removed as well.\n");
          }
          update_pymol(od);
          update_jmol(od);
          tee_printf(od, M_TOOL_SUCCESS, nesting, command, "REMOVE_OBJECT");
          tee_flush(od);
        }
      }
      else {
        od->valid &= SDF_BIT;
      }
    }  
    else if (!strcasecmp(arg->me[0], "remove_y_vars")) {
      gettimeofday(&start, NULL);
      if (od->y_vars) {
        if (!(run_type & DRY_RUN)) {
          strcpy(comma_hyphen_list, "all");
          if ((parameter = get_args(od, "y_var_list"))) {
            strcpy(comma_hyphen_list, parameter);
          }
          result = parse_comma_hyphen_list_to_array
            (od, comma_hyphen_list, Y_VAR_LIST);
          switch (result) {
            case OUT_OF_MEMORY:
            tee_error(od, run_type, overall_line_num,
              E_OUT_OF_MEMORY, REMOVE_Y_VARS_FAILED);
            return PARSE_INPUT_ERROR;

            case INVALID_LIST_RANGE:
            tee_error(od, run_type, overall_line_num,
              E_LIST_PARSING, "Y variables", "EXPORT", REMOVE_Y_VARS_FAILED);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
          result = set(od, Y_VAR_LIST, OPERATE_BIT, 1, SILENT);
          if (result) {
            tee_error(od, run_type, overall_line_num, E_Y_VAR_ALLOWED_RANGE,
              od->y_vars, REMOVE_Y_VARS_FAILED);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
          result = remove_y_vars(od);
          switch (result) {
            case CANNOT_WRITE_TEMP_FILE:
            tee_error(od, run_type, overall_line_num, E_ERROR_IN_WRITING_TEMP_FILE,
              od->file[DEP_IN]->name, REMOVE_Y_VARS_FAILED);
            return PARSE_INPUT_ERROR;

            case CANNOT_READ_TEMP_FILE:
            case NOT_ENOUGH_OBJECTS:
            case WRONG_NUMBER_OF_Y_VARS:
            case PREMATURE_DEP_IN_EOF:
            case CANNOT_FIND_Y_VAR_NAME:
            tee_error(od, run_type, overall_line_num, E_ERROR_IN_READING_TEMP_FILE,
              od->file[DEP_IN]->name, REMOVE_Y_VARS_FAILED);
            return PARSE_INPUT_ERROR;

            case Y_VAR_LOW_SD:
            tee_error(od, run_type, overall_line_num,
              E_Y_VAR_LOW_SD, REMOVE_Y_VARS_FAILED);
            return PARSE_INPUT_ERROR;
          }
        }
      }
      else {
        tee_printf(od, "No Y variables are present; nothing was done.\n\n");
      }
      tee_printf(od, M_TOOL_SUCCESS, nesting, command, "REMOVE_Y_VARS");
      tee_flush(od);
    }
    else if (!strcasecmp(arg->me[0], "print")) {
      gettimeofday(&start, NULL);
      if (!(parameter = get_args(od, "type"))) {
        tee_error(od, run_type, overall_line_num,
          "Please specify whether PRINT should print "
          "\"2_LEVEL\", \"3_LEVEL\", \"4_LEVEL\", \"D_OPTIMAL\", "
          "\"FFDSEL\", \"SEEDS\", \"UVEPLS\" or \"GROUPS\" "
          "variables.\n%s",
          PRINT_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      if ((!strncasecmp(parameter, "2_level", 7))
        || (!strncasecmp(parameter, "3_level", 7))
        || (!strncasecmp(parameter, "4_level", 7))
        || (!strncasecmp(parameter, "d_optimal", 9))
        || (!strncasecmp(parameter, "ffdsel", 6))
        || (!strncasecmp(parameter, "uvepls", 6))
        || (!strncasecmp(parameter, "seed", 4))
        || (!strncasecmp(parameter, "group", 5))
        ) {
        type = (int)tolower(parameter[0]);
        if (!(run_type & DRY_RUN)) {
          if (!(od->field_num)) {
            tee_error(od, run_type, overall_line_num,
              E_NO_FIELDS_PRESENT, PRINT_FAILED);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
          strcpy(comma_hyphen_list, "all");
          if ((parameter = get_args(od, "field_list"))) {
            strcpy(comma_hyphen_list, parameter);
          }
          result = parse_comma_hyphen_list_to_array
            (od, comma_hyphen_list, FIELD_LIST);
          switch (result) {
            case OUT_OF_MEMORY:
            tee_error(od, run_type, overall_line_num,
              E_OUT_OF_MEMORY, PRINT_FAILED);
            return PARSE_INPUT_ERROR;

            case INVALID_LIST_RANGE:
            tee_error(od, run_type, overall_line_num,
              E_LIST_PARSING, "fields",
              "PRINT", PRINT_FAILED);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
          result = set(od, FIELD_LIST,
            OPERATE_BIT, 1, SILENT);
          if (result) {
            tee_error(od, run_type, overall_line_num,
              E_ALLOWED_FIELD_RANGE,
              od->field_num, PRINT_FAILED);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
        }
        switch (type) {
          case GROUPS:
          if (!(od->valid & SEED_BIT)) {
            tee_error(od, run_type, overall_line_num,
              E_NO_VALID_SELECTION, "SRD", PRINT_FAILED);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
          strcpy(comma_hyphen_list, "all");
          if ((parameter = get_args(od, "group_list"))) {
            strcpy(comma_hyphen_list, parameter);
          }
          if (!(run_type & DRY_RUN)) {
            result = parse_comma_hyphen_list_to_array
              (od, comma_hyphen_list, GROUP_LIST);
            switch (result) {
              case OUT_OF_MEMORY:
              tee_error(od, run_type, overall_line_num,
                E_OUT_OF_MEMORY, PRINT_FAILED);
              return PARSE_INPUT_ERROR;

              case INVALID_LIST_RANGE:
              tee_error(od, run_type, overall_line_num,
                E_LIST_PARSING, "SRD groups", "PRINT",
                PRINT_FAILED);
              return PARSE_INPUT_ERROR;
            }
          }
          break;

          case SEEDS:
          if (!(od->valid & SEED_BIT)) {
            tee_error(od, run_type, overall_line_num,
              E_NO_VALID_SELECTION, "SRD seeds", PRINT_FAILED);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
          break;

          case D_OPTIMAL:
          if (!(od->valid & D_OPTIMAL_BIT)) {
            tee_error(od, run_type, overall_line_num,
              E_NO_VALID_SELECTION, "D-optimal", PRINT_FAILED);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
          break;

          case FFDSEL:
          if (!(od->valid & FFDSEL_BIT)) {
            tee_error(od, run_type, overall_line_num,
              E_NO_VALID_SELECTION, "FFDSEL", PRINT_FAILED);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
          break;

          case UVEPLS:
          if (!(od->valid & UVEPLS_BIT)) {
            tee_error(od, run_type, overall_line_num,
              E_NO_VALID_SELECTION, "UVE-PLS", PRINT_FAILED);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
          break;

          case TWO_LEVEL:
          if (!(od->valid & TWO_LEVEL_BIT)) {
            tee_error(od, run_type, overall_line_num,
              E_NO_VALID_SELECTION, "2-level", PRINT_FAILED);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
          break;

          case THREE_LEVEL:
          if (!(od->valid & THREE_LEVEL_BIT)) {
            tee_error(od, run_type, overall_line_num,
              E_NO_VALID_SELECTION, "3-level", PRINT_FAILED);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
          break;

          case FOUR_LEVEL:
          if (!(od->valid & FOUR_LEVEL_BIT)) {
            tee_error(od, run_type, overall_line_num,
              E_NO_VALID_SELECTION, "4-level", PRINT_FAILED);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
          break;
        }
        if (!(run_type & DRY_RUN)) {
          ++command;
          tee_printf(od, M_TOOL_INVOKE, nesting, command, "PRINT", line_orig);
          tee_flush(od);
          result = print_variables(od, type);
          gettimeofday(&end, NULL);
          elapsed_time(od, &start, &end);
          switch (result) {

            case OUT_OF_MEMORY:
            tee_error(od, run_type, overall_line_num,
              E_OUT_OF_MEMORY, PRINT_FAILED);
            return PARSE_INPUT_ERROR;

            case ONE_FIELD_AT_A_TIME:
            tee_error(od, run_type, overall_line_num,
              E_GROUPS_OTHER_THAN_ZERO, "listed", PRINT_FAILED);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;

            case INVALID_LIST_RANGE:
            tee_error(od, run_type, overall_line_num,
              "The range of SRD groups specified is not valid.\n%s",
              PRINT_FAILED);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;

            default:
            tee_printf(od, M_TOOL_SUCCESS, nesting, command, "PRINT");
            tee_flush(od);
          }
        }
      }
      else {
        tee_error(od, run_type, overall_line_num,
          "Only \"2_LEVEL\", \"3_LEVEL\", \"4_LEVEL\", "
          "\"D_OPTIMAL\", \"FFDSEL\", \"UVEPLS\", "
          "\"SEEDS\" and \"GROUPS\" "
          "types are allowed.\n%s",
          PRINT_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
    }
    else if (!strcasecmp(arg->me[0], "export")) {
      gettimeofday(&start, NULL);
      requested_endianness = machine_type();
      if ((parameter = get_args(od, "endianness"))) {
        if (!strncasecmp(parameter, "little_endian", 13)) {
          requested_endianness = O3Q_LITTLE_ENDIAN;
        }
        else if (!strncasecmp(parameter, "big_endian", 10)) {
          requested_endianness = O3Q_BIG_ENDIAN;
        }
        else if (strncasecmp(parameter, "native", 6)) {
          tee_error(od, run_type, overall_line_num,
            "Allowed values for endianness are "
            "\"LITTLE_ENDIAN\", \"BIG_ENDIAN\" "
            "AND \"NATIVE\".\n%s",
            EXPORT_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
      }
      format = INSIGHT_FORMAT;
      if ((parameter = get_args(od, "format"))) {
        if (!strncasecmp(parameter, "maestro", 7)) {
          format = MAESTRO_FORMAT;
        }
        else if (!strncasecmp(parameter, "moe", 3)) {
          format = MOE_FORMAT;
          requested_endianness = O3Q_LITTLE_ENDIAN;
        }
        else if (!strncasecmp(parameter, "sybyl", 5)) {
          format = SYBYL_FORMAT;
        }
        else if (!strncasecmp(parameter, "ascii", 5)) {
          format = ASCII_FORMAT;
        }
        else if (!strncasecmp(parameter, "xyz", 3)) {
          format = XYZ_FORMAT;
        }
        else if (strncasecmp(parameter, "insight", 7)) {
          tee_error(od, run_type, overall_line_num,
            "Allowed values for format are \"INSIGHT\", "
            "\"MAESTRO\", \"MOE\", \"SYBYL\", "
            "\"ASCII\" and \"XYZ\".\n%s",
            EXPORT_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
      }
      sign = 0;
      if ((parameter = get_args(od, "sign"))) {
        if (!strncasecmp(parameter, "pos", 3)) {
          sign = 1;
        }
        else if (!strncasecmp(parameter, "neg", 3)) {
          sign = -1;
        }
        else if (strncasecmp(parameter, "both", 4)) {
          tee_error(od, run_type, overall_line_num,
            "Allowed values for sign are \"POS\", "
            "\"NEG\" and \"BOTH\".\n%s",
            EXPORT_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
      }
      if ((parameter = get_args(od, "type"))) {
        if ((!strncasecmp(parameter, "2_level", 7))
          || (!strncasecmp(parameter, "3_level", 7))
          || (!strncasecmp(parameter, "4_level", 7))
          || (!strncasecmp(parameter, "pca_loading", 11))
          || (!strncasecmp(parameter, "field_sd", 8))
          || (!strncasecmp(parameter, "object_field", 12))
          || (!strncasecmp(parameter, "weight", 6))
          || (!strncasecmp(parameter, "loading", 7))
          || (!strncasecmp(parameter, "coefficient", 11))
          || (!strncasecmp(parameter, "d_optimal", 9))
          || (!strncasecmp(parameter, "ffdsel", 6))
          || (!strncasecmp(parameter, "uvepls", 6))
          || (!strncasecmp(parameter, "group", 5))
          ) {
          type = (int)tolower(parameter[0]);
          if ((type == WEIGHTS) || (type == LOADINGS)
            || (type == PCA_LOADINGS) || (type == COEFFICIENTS)) {
            if (type == PCA_LOADINGS) {
              if (!(od->valid & PCA_BIT)) {
                tee_error(od, run_type, overall_line_num,
                  E_NO_PCA_ANALYSIS, EXPORT_FAILED);
                fail = !(run_type & INTERACTIVE_RUN);
                continue;
              }
              if (!(run_type & DRY_RUN)) {
                od->file[TEMP_PCA_LOADINGS]->handle =
                  fopen(od->file[TEMP_PCA_LOADINGS]->name, "rb");
                if (!(od->file[TEMP_PCA_LOADINGS]->handle)) {
                  tee_error(od, run_type, overall_line_num,
                    E_ERROR_IN_READING_TEMP_FILE, "TEMP_FIELD", EXPORT_FAILED);
                  return PARSE_INPUT_ERROR;
                }
                actual_len = fread(&pc, sizeof(int), 1,
                  od->file[TEMP_PCA_LOADINGS]->handle);
                if (actual_len != 1) {
                  tee_error(od, run_type, overall_line_num,
                    E_ERROR_IN_READING_TEMP_FILE, "TEMP_FIELD", EXPORT_FAILED);
                  return PARSE_INPUT_ERROR;
                }
              }
            }
            else {
              pc = od->pc_num;
            }
            if (!(run_type & DRY_RUN)) {
              if (!(od->field_num)) {
                tee_error(od, run_type, overall_line_num,
                  E_NO_FIELDS_PRESENT, EXPORT_FAILED);
                fail = !(run_type & INTERACTIVE_RUN);
                continue;
              }
              max_pc = pc;
              if ((parameter = get_args(od, "pc"))) {
                sscanf(parameter, "%d", &pc);
              }
              if (pc < 1) {
                tee_error(od, run_type, overall_line_num,
                  E_AT_LEAST_ONE_PC, EXPORT_FAILED);
                fail = !(run_type & INTERACTIVE_RUN);
                continue;
              }
              if (pc > max_pc) {
                if (type == PCA_LOADINGS) {
                  tee_error(od, run_type, overall_line_num,
                    E_TOO_FEW_MANY_FOR_CURRENT_PCA_ANALYSIS,
                    "many PCs", EXPORT_FAILED);
                }
                else {
                  tee_error(od, run_type, overall_line_num,
                    E_TOO_FEW_MANY_FOR_CURRENT_PLS_MODEL,
                    "many PCs", EXPORT_FAILED);
                }
                fail = !(run_type & INTERACTIVE_RUN);
                continue;
              }
            }
          }
          if (!(parameter = get_args(od, "file"))) {
            tee_error(od, run_type, overall_line_num,
              "Please specify a file basename "
              "where grid data should be saved.\n%s",
              EXPORT_FAILED);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
          if (!(run_type & DRY_RUN)) {
            memset(file_basename, 0, BUF_LEN);
            strcpy(file_basename, parameter);
            strcpy(comma_hyphen_list, "all");
            if ((parameter = get_args(od, "field_list"))) {
              strcpy(comma_hyphen_list, parameter);
            }
            result = parse_comma_hyphen_list_to_array(od,
              comma_hyphen_list, FIELD_LIST);
            switch (result) {
              case OUT_OF_MEMORY:
              tee_error(od, run_type, overall_line_num,
                E_OUT_OF_MEMORY, EXPORT_FAILED);
              return PARSE_INPUT_ERROR;

              case INVALID_LIST_RANGE:
              tee_error(od, run_type, overall_line_num,
                E_LIST_PARSING, "fields", "EXPORT", EXPORT_FAILED);
              fail = !(run_type & INTERACTIVE_RUN);
              continue;
            }
            result = set(od, FIELD_LIST, OPERATE_BIT, 1, SILENT);
            if (result) {
              tee_error(od, run_type, overall_line_num, E_ALLOWED_FIELD_RANGE,
                od->field_num, EXPORT_FAILED);
              fail = !(run_type & INTERACTIVE_RUN);
              continue;
            }
          }
          switch (type) {
            case TWO_LEVEL:
            if (!(od->valid & TWO_LEVEL_BIT)) {
              tee_error(od, run_type, overall_line_num,
                E_NO_VALID_SELECTION, "2-level", EXPORT_FAILED);
              fail = !(run_type & INTERACTIVE_RUN);
              continue;
            }
            break;

            case THREE_LEVEL:
            if (!(od->valid & THREE_LEVEL_BIT)) {
              tee_error(od, run_type, overall_line_num,
                E_NO_VALID_SELECTION, "3-level", EXPORT_FAILED);
              fail = !(run_type & INTERACTIVE_RUN);
              continue;
            }
            break;

            case FOUR_LEVEL:
            if (!(od->valid & FOUR_LEVEL_BIT)) {
              tee_error(od, run_type, overall_line_num,
                E_NO_VALID_SELECTION, "4-level", EXPORT_FAILED);
              fail = !(run_type & INTERACTIVE_RUN);
              continue;
            }
            break;

            case SD:
            if (!(run_type & DRY_RUN)) {
              if (!(od->object_num)) {
                tee_error(od, run_type, overall_line_num,
                  E_NO_OBJECTS_PRESENT, EXPORT_FAILED);
                fail = !(run_type & INTERACTIVE_RUN);
                continue;
              }
            }
            break;

            case OBJECT_FIELD:
            if (!(run_type & DRY_RUN)) {
              if (!(od->object_num)) {
                tee_error(od, run_type, overall_line_num,
                  E_NO_OBJECTS_PRESENT, EXPORT_FAILED);
                fail = !(run_type & INTERACTIVE_RUN);
                continue;
              }
              label = 0;
              if ((parameter = get_args(od, "inactive"))) {
                if (!strncasecmp(parameter, "label", 5)) {
                  label = LABEL_INACTIVE_BIT;
                }
                else if (!strncasecmp(parameter, "zero", 4)) {
                  label = ZERO_INACTIVE_BIT;
                }
                else if (strncasecmp(parameter, "orig", 4)) {
                  tee_error(od, run_type, overall_line_num,
                    "Allowed values for inactive are "
                    "\"ORIGINAL\", \"ZERO\" and \"LABEL\".\n%s",
                    EXPORT_FAILED);
                  fail = !(run_type & INTERACTIVE_RUN);
                  continue;
                }
              }
              if ((parameter = get_args(od, "missing"))) {
                if (!strncasecmp(parameter, "label", 5)) {
                  label |= LABEL_MISSING_BIT;
                }
                else if (strncasecmp(parameter, "zero", 4)) {
                  tee_error(od, run_type, overall_line_num,
                    "Allowed values for missing are "
                    "\"ZERO\" and \"LABEL\".\n%s",
                    EXPORT_FAILED);
                  fail = !(run_type & INTERACTIVE_RUN);
                  continue;
                }
              }
              result = parse_synonym_lists(od, "EXPORT", EXPORT_FAILED,
                (1 << OBJECT_LIST) | (1 << ID_LIST) | (1 << STRUCT_LIST),
                &list_type, OBJECT_LIST, run_type, overall_line_num);
              switch (result) {
                case PARSE_INPUT_RECOVERABLE_ERROR:
                fail = !(run_type & INTERACTIVE_RUN);
                continue;

                case PARSE_INPUT_ERROR:
                return PARSE_INPUT_ERROR;
              }
            }
            break;

            case PCA_LOADINGS:
            if (!(od->valid & PCA_BIT)) {
              tee_error(od, run_type, overall_line_num,
                E_NO_PCA_ANALYSIS, EXPORT_FAILED);
              fail = !(run_type & INTERACTIVE_RUN);
              continue;
            }
            break;

            case WEIGHTS:
            case LOADINGS:
            if (!(od->valid & PLS_BIT)) {
              tee_error(od, run_type, overall_line_num,
                E_NO_PLS_MODEL, EXPORT_FAILED);
              fail = !(run_type & INTERACTIVE_RUN);
              continue;
            }
            break;

            case D_OPTIMAL:
            if (!(od->valid & D_OPTIMAL_BIT)) {
              tee_error(od, run_type, overall_line_num,
                E_NO_VALID_SELECTION, "D-optimal", EXPORT_FAILED);
              fail = !(run_type & INTERACTIVE_RUN);
              continue;
            }
            break;

            case FFDSEL:
            if (!(od->valid & FFDSEL_BIT)) {
              tee_error(od, run_type, overall_line_num,
                E_NO_VALID_SELECTION, "FFD", EXPORT_FAILED);
              fail = !(run_type & INTERACTIVE_RUN);
              continue;
            }
            break;

            case UVEPLS:
            if (!(od->valid & UVEPLS_BIT)) {
              tee_error(od, run_type, overall_line_num,
                E_NO_VALID_SELECTION, "UVE-PLS", EXPORT_FAILED);
              fail = !(run_type & INTERACTIVE_RUN);
              continue;
            }
            break;

            case COEFFICIENTS:
            if (!(od->valid & PLS_BIT)) {
              tee_error(od, run_type, overall_line_num,
                E_NO_PLS_MODEL, EXPORT_FAILED);
              fail = !(run_type & INTERACTIVE_RUN);
              continue;
            }
            if (!(run_type & DRY_RUN)) {
              strcpy(comma_hyphen_list, "all");
              if ((parameter = get_args(od, "y_var_list"))) {
                strcpy(comma_hyphen_list, parameter);
              }
              result = parse_comma_hyphen_list_to_array
                (od, comma_hyphen_list, Y_VAR_LIST);
              switch (result) {
                case OUT_OF_MEMORY:
                tee_error(od, run_type, overall_line_num,
                  E_OUT_OF_MEMORY, EXPORT_FAILED);
                return PARSE_INPUT_ERROR;

                case INVALID_LIST_RANGE:
                tee_error(od, run_type, overall_line_num,
                  E_LIST_PARSING, "Y variables", "EXPORT", EXPORT_FAILED);
                fail = !(run_type & INTERACTIVE_RUN);
                continue;
              }
              result = set(od, Y_VAR_LIST, OPERATE_BIT, 1, SILENT);
              if (result) {
                tee_error(od, run_type, overall_line_num, E_Y_VAR_ALLOWED_RANGE,
                  od->y_vars, EXPORT_FAILED);
                fail = !(run_type & INTERACTIVE_RUN);
                continue;
              }
            }
            break;

            case GROUPS:
            if (!(od->valid & SEED_BIT)) {
              tee_error(od, run_type, overall_line_num,
                E_NO_VALID_SELECTION, "SRD", EXPORT_FAILED);
              fail = !(run_type & INTERACTIVE_RUN);
              continue;
            }
            if (!(parameter = get_args(od, "group_list"))) {
              tee_error(od, run_type, overall_line_num,
                "Please specify a list of SRD groups "
                "which EXPORT should display.\n%s",
                EXPORT_FAILED);
              fail = !(run_type & INTERACTIVE_RUN);
              continue;
            }
            if (!(run_type & DRY_RUN)) {
              result = parse_comma_hyphen_list_to_array(od,
                parameter, GROUP_LIST);
            }
            break;
          }
        }
        else {
          tee_error(od, run_type, overall_line_num,
            "Only \"WEIGHTS\", \"LOADINGS\", \"PCA_LOADINGS\", "
            "\"COEFFICIENTS\", \"2_LEVEL\", \"3_LEVEL\", "
            "\"4_LEVEL\", \"D_OPTIMAL\", \"FFDSEL\", "
            "\"FIELD_SD\", \"SRD\" and \"OBJECT_FIELD\" "
            "types are allowed.\n%s",
            EXPORT_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
      }
      else {
        tee_error(od, run_type, overall_line_num,
          "Please specify which type of data "
          "should be written in .grd format.\n%s",
          EXPORT_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      interpolate = 0;
      if ((parameter = get_args(od, "interpolate"))) {
        sscanf(parameter, "%d", &interpolate);
      }
      if ((interpolate < 0) || (interpolate > 5)) {
        tee_error(od, run_type, overall_line_num,
          "The level of interpolation should "
          "lie in the range 0 - 5.\n%s",
          EXPORT_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      if (!(run_type & DRY_RUN)) {
        ++command;
        tee_printf(od, M_TOOL_INVOKE, nesting, command, "EXPORT", line_orig);
        tee_flush(od);
        result = grid_write(od, file_basename, pc, type,
          sign, format, label, interpolate,
          requested_endianness);
        gettimeofday(&end, NULL);
        elapsed_time(od, &start, &end);
        switch (result) {
          case PREMATURE_EOF:
          tee_error(od, run_type, overall_line_num,
            "Error while writing .grd file \"%s\".\n%s",
            od->file[GRD_OUT]->name, EXPORT_FAILED);
          O3_ERROR_PRINT(&(od->task));
          return PARSE_INPUT_ERROR;

          case CANNOT_READ_TEMP_FILE:
          tee_error(od, run_type, overall_line_num,
            E_ERROR_IN_READING_TEMP_FILE, od->task.string, EXPORT_FAILED);
          O3_ERROR_PRINT(&(od->task));
          return PARSE_INPUT_ERROR;

          case CANNOT_WRITE_TEMP_FILE:
          tee_error(od, run_type, overall_line_num,
            E_ERROR_IN_WRITING_TEMP_FILE, od->task.string, EXPORT_FAILED);
          O3_ERROR_PRINT(&(od->task));
          return PARSE_INPUT_ERROR;

          case OUT_OF_MEMORY:
          tee_error(od, run_type, overall_line_num,
            E_OUT_OF_MEMORY, EXPORT_FAILED);
          return PARSE_INPUT_ERROR;

          case ONE_FIELD_AT_A_TIME:
          tee_error(od, run_type, overall_line_num,
            E_GROUPS_OTHER_THAN_ZERO, "displayed", EXPORT_FAILED);
          return PARSE_INPUT_ERROR;

          case INVALID_LIST_RANGE:
          tee_error(od, run_type, overall_line_num,
            "The range of SRD groups specified is not valid.\n%s",
            EXPORT_FAILED);
          return PARSE_INPUT_ERROR;

          default:
          tee_printf(od, M_TOOL_SUCCESS, nesting, command, "EXPORT");
          tee_flush(od);
        }
      }
    }
    else if ((!strcasecmp(arg->me[0], "cd"))
      || (!strcasecmp(arg->me[0], "chdir"))) {
      if (!(run_type & DRY_RUN)) {
        memset(file_basename, 0, BUF_LEN);
        if ((parameter = get_args(od, "dir"))) {
          strncpy(file_basename, parameter, BUF_LEN - 1);
          file_basename[BUF_LEN - 1] = '\0';
        }
        ++command;
        tee_printf(od, M_TOOL_INVOKE, nesting, command, "CHDIR", line_orig);
        tee_flush(od);
        result = 0;
        if (file_basename[0]) {
          result = chdir(file_basename);
        }
        if (result) {
          tee_error(od, run_type, overall_line_num,
            E_CANNOT_CHANGE_DIR, file_basename, CHANGE_DIR_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        else {
          if (getcwd(file_basename, BUF_LEN - 1)) {
            tee_printf(od, "Current directory: %s\n\n", file_basename);
            tee_printf(od, M_TOOL_SUCCESS, nesting, command, "CHDIR");
            tee_flush(od);
          }
        }
      }
    }
    else if (!strcasecmp(arg->me[0], "source")) {
      if (!(parameter = get_args(od, "file"))) {
        tee_error(od, run_type, overall_line_num,
          "Please specify the input file you wish "
          "to source.\n%s",
          SOURCE_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      source = &(od->mel.source);
      while (*source) {
        source = &((*source)->next);
      }
      *source = (FileDescriptor *)malloc(sizeof(FileDescriptor));
      if (!(*source)) {
        tee_error(od, run_type, overall_line_num, E_OUT_OF_MEMORY,
          SOURCE_FAILED);
        return PARSE_INPUT_ERROR;
      }
      memset(*source, 0, sizeof(FileDescriptor));
      strcpy((*source)->name, parameter);
      (*source)->handle = fopen((*source)->name, "rb");
      if (!((*source)->handle)) {
        tee_error(od, run_type, overall_line_num,
          E_FILE_CANNOT_BE_OPENED_FOR_READING,
          parameter, SOURCE_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      ++command;
      tee_printf(od, M_TOOL_INVOKE, nesting, command, "SOURCE", line_orig);
      tee_flush(od);
      if (run_type & DRY_RUN) {
        result = parse_input(od, (*source)->handle, DRY_RUN);
        if (result) {
          tee_error(od, run_type, overall_line_num,
            E_WHILE_SOURCING, (*source)->name,
            SOURCE_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
      }
      result = parse_input(od, (*source)->handle, 0);
      if (result) {
        tee_error(od, run_type, overall_line_num,
          E_WHILE_SOURCING, (*source)->name,
          SOURCE_FAILED);
      }
      fclose((*source)->handle);
      (*source)->handle = NULL;
      free(*source);
      *source = NULL;
      switch (result) {
        case PARSE_INPUT_RECOVERABLE_ERROR:
        fail = !(run_type & INTERACTIVE_RUN);
        break;

        case PARSE_INPUT_ERROR:
        return PARSE_INPUT_ERROR;

        default:
        tee_printf(od, M_TOOL_SUCCESS, nesting, command, "SOURCE");
        tee_flush(od);
      }
      continue;
    }
    else if (!strcasecmp(arg->me[0], "dataset")) {
      update_field_object_attr(od, VERBOSE_BIT);
      calc_active_vars(od, FULL_MODEL);
    }
    else if ((!strcasecmp(arg->me[0], "exit"))
      || (!strcasecmp(arg->me[0], "quit"))
      || (!strcasecmp(arg->me[0], "stop"))) {
      break;
    }
    else if (!strcasecmp(arg->me[0], "pls")) {
      if (!(run_type & DRY_RUN)) {
        gettimeofday(&start, NULL);
        /*
        user wants to accomplish
        PLS analysis
        */
        if (!(od->field_num)) {
          tee_error(od, run_type, overall_line_num,
            E_NO_FIELDS_PRESENT, PLS_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        pc = 5;
        if ((parameter = get_args(od, "pc"))) {
          sscanf(parameter, "%d", &pc);
        }
        if (pc > od->overall_active_x_vars) {
          tee_error(od, run_type, overall_line_num,
            E_TOO_FEW_MANY_FOR_AVAILABLE_DATA, "many PCs",
            PLS_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        if (pc < 1) {
          tee_error(od, run_type, overall_line_num,
            E_AT_LEAST_ONE_PC, PLS_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        if (!(od->y_vars)) {
          tee_error(od, run_type, overall_line_num,
            E_NO_Y_VARS_PRESENT, PLS_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        /*
        perform pls operation
        */
        od->mal.press = double_mat_resize
          (od->mal.press, pc + 1, od->y_vars);
        if (!(od->mal.press)) {
          tee_error(od, run_type, overall_line_num,
            E_OUT_OF_MEMORY, PLS_FAILED);
          return PARSE_INPUT_ERROR;
        }
        memset(od->mal.press->base, 0,
          od->mal.press->m * od->mal.press->n * sizeof(double));
        result = alloc_pls(od, od->overall_active_x_vars,
          pc, FULL_MODEL);
        if (result) {
          tee_error(od, run_type, overall_line_num,
            E_OUT_OF_MEMORY, PLS_FAILED);
          return PARSE_INPUT_ERROR;
        }
        result = calc_active_vars(od, FULL_MODEL);
        switch (result) {
          case CANNOT_READ_TEMP_FILE:
          tee_error(od, run_type, overall_line_num,
            E_ERROR_IN_READING_TEMP_FILE, "TEMP_FIELD", PLS_FAILED);
          return PARSE_INPUT_ERROR;

          case Y_VAR_LOW_SD:
          tee_error(od, run_type, overall_line_num,
            E_Y_VAR_LOW_SD, PLS_FAILED);
          return PARSE_INPUT_ERROR;
        }
        result = fill_x_matrix(od, FULL_MODEL, 0);
        switch (result) {
          case OUT_OF_MEMORY:
          tee_error(od, run_type, overall_line_num,
            E_OUT_OF_MEMORY, PLS_FAILED);
          return PARSE_INPUT_ERROR;
        }
        result = fill_y_matrix(od);
        if (result) {
          tee_error(od, run_type, overall_line_num,
            E_OUT_OF_MEMORY, PLS_FAILED);
          return PARSE_INPUT_ERROR;
        }
        trim_mean_center_matrix(od, od->mal.large_e_mat,
          &(od->mal.e_mat), &(od->vel.e_mat_ave),
          FULL_MODEL, od->active_object_num);
        trim_mean_center_matrix(od, od->mal.large_f_mat,
          &(od->mal.f_mat), &(od->vel.f_mat_ave),
          FULL_MODEL, od->active_object_num);
        od->file[ASCII_IN]->name[0] = '\0';
        if ((parameter = get_args(od, "file"))) {
          strcpy(od->file[ASCII_IN]->name, parameter);
          if (!(od->file[ASCII_IN]->handle =
            fopen(od->file[ASCII_IN]->name, "wb+"))) {
            tee_error(od, run_type, overall_line_num,
              E_FILE_CANNOT_BE_OPENED_FOR_WRITING,
              od->file[ASCII_IN]->name, PLS_FAILED);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
        }
        ++command;
        tee_printf(od, M_TOOL_INVOKE, nesting, command, "PLS", line_orig);
        tee_flush(od);
        pls(od, pc, FULL_MODEL);
        gettimeofday(&end, NULL);
        elapsed_time(od, &start, &end);
        tee_flush(od);
        result = store_weights_loadings(od);
        if (result) {
          tee_error(od, run_type, overall_line_num,
            E_ERROR_IN_WRITING_TEMP_FILE, "TEMP_FIELD", PLS_FAILED);
          return PARSE_INPUT_ERROR;
        }
        result = open_temp_file(od, od->file[TEMP_CALC], "calc_y");
        if (result) {
          tee_error(od, run_type, overall_line_num,
            E_TEMP_FILE_CANNOT_BE_OPENED_FOR_WRITING,
            od->file[TEMP_CALC]->name, PLS_FAILED);
          return PARSE_INPUT_ERROR;
        }
        result = open_temp_file(od, od->file[TEMP_PLS_COEFF], "x_coefficients");
        if (result) {
          tee_error(od, run_type, overall_line_num,
            E_TEMP_FILE_CANNOT_BE_OPENED_FOR_WRITING,
            od->file[TEMP_PLS_COEFF]->name, PLS_FAILED);
          return PARSE_INPUT_ERROR;
        }
        result = calc_y_values(od);
        if (od->file[TEMP_PLS_COEFF]->handle) {
          fclose(od->file[TEMP_PLS_COEFF]->handle);
          od->file[TEMP_PLS_COEFF]->handle = NULL;
        }
        switch (result) {
          case CANNOT_WRITE_TEMP_FILE:
          tee_error(od, run_type, overall_line_num,
            E_ERROR_IN_WRITING_TEMP_FILE, "TEMP_FIELD", PLS_FAILED);
          return PARSE_INPUT_ERROR;

          case CANNOT_READ_TEMP_FILE:
          tee_error(od, run_type, overall_line_num,
            E_ERROR_IN_READING_TEMP_FILE, "TEMP_FIELD", PLS_FAILED);
          return PARSE_INPUT_ERROR;
        }
        result = print_calc_values(od);
        switch (result) {
          case CANNOT_READ_TEMP_FILE:
          tee_error(od, run_type, overall_line_num,
            E_ERROR_IN_READING_TEMP_FILE, "TEMP_FIELD", PLS_FAILED);
          return PARSE_INPUT_ERROR;
        }
        if (od->file[TEMP_CALC]->handle) {
          fclose(od->file[TEMP_CALC]->handle);
          od->file[TEMP_CALC]->handle = NULL;
        }
        options = 0;
        if ((parameter = get_args(od, "scores"))) {
          if (!strncmp(parameter, "x", 1)) {
            options = X_SCORES;
          }
          else if (!strncmp(parameter, "y", 1)) {
            options = Y_SCORES;
          }
          else if (!strncmp(parameter, "both", 4)) {
            options = X_SCORES | Y_SCORES;
          }
          else if (strncmp(parameter, "none", 4)) {
            tee_error(od, run_type, overall_line_num,
              "Only \"X\", \"Y\", \"BOTH\", "
              "\"NONE\" keywords are allowed "
              "as score types.\n%s", EXPORT_FAILED);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
          print_pls_scores(od, options);
        }
        tee_printf(od, M_TOOL_SUCCESS, nesting, command, "PLS");
        tee_flush(od);
        od->valid |= PLS_BIT;
      }
      else {
        od->valid |= PLS_BIT;
      }
      if (od->file[ASCII_IN]->handle) {
        fclose(od->file[ASCII_IN]->handle);
        od->file[ASCII_IN]->handle = NULL;
      }
    }
    else if (!strcasecmp(arg->me[0], "cv")) {
      gettimeofday(&start, NULL);
      if (!(od->valid & PLS_BIT)) {
        tee_error(od, run_type, overall_line_num,
          E_NO_PLS_MODEL, CV_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      type = LEAVE_ONE_OUT;
      if ((parameter = get_args(od, "type"))) {
        if (!strncasecmp(parameter, "lto", 3)) {
          type = LEAVE_TWO_OUT;
        }
        else if (!strncasecmp(parameter, "lmo", 3)) {
          type = LEAVE_MANY_OUT;
        }
        else if (strncasecmp(parameter, "loo", 3)) {
          tee_error(od, run_type, overall_line_num,
            E_ONLY_LOO_LTO_LMO_CV_ALLOWED, CV_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
      }
      if (!(run_type & DRY_RUN)) {
        pc = od->pc_num;
        if ((parameter = get_args(od, "pc"))) {
          sscanf(parameter, "%d", &pc);
        }
        if (pc < 1) {
          tee_error(od, run_type, overall_line_num,
            E_AT_LEAST_ONE_PC, CV_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        if (pc > od->pc_num) {
          tee_error(od, run_type, overall_line_num,
            E_TOO_FEW_MANY_FOR_CURRENT_PLS_MODEL,
            "many PCs", CV_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        od->cv.pc_num = pc;
        ++command;
        tee_printf(od, M_TOOL_INVOKE, nesting, command, "CV", line_orig);
        switch (type) {
          case LEAVE_ONE_OUT:
          tee_printf(od, "LOO CV was chosen.\n\n");
          break;

          case LEAVE_TWO_OUT:
          tee_printf(od, "LTO CV was chosen.\n\n");
          break;

          case LEAVE_MANY_OUT:
          tee_printf(od, "LMO CV was chosen.\n\n");
          result = check_lmo_parameters(od, CV_FAILED,
            &groups, &runs, run_type, overall_line_num);
          if (result) {
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
          break;
        }
        tee_flush(od);
        if ((pc > od->pc_num)
          || (pc > od->overall_active_x_vars)) {
          tee_error(od, run_type, overall_line_num,
            E_TOO_FEW_MANY_FOR_AVAILABLE_DATA, "many PCs",
            CV_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        od->mal.press = double_mat_resize
          (od->mal.press, pc + 1, od->y_vars);
        if (!(od->mal.press)) {
          tee_error(od, run_type, overall_line_num,
            E_OUT_OF_MEMORY, CV_FAILED);
          return PARSE_INPUT_ERROR;
        }
        memset(od->mal.press->base, 0,
          od->mal.press->m * od->mal.press->n * sizeof(double));
        result = alloc_pls(od, od->overall_active_x_vars, pc, CV_MODEL);
        if (result) {
          tee_error(od, run_type, overall_line_num,
            E_OUT_OF_MEMORY, CV_FAILED);
          return PARSE_INPUT_ERROR;
        }
        result = fill_x_matrix(od, CV_MODEL, 0);
        switch (result) {
          case OUT_OF_MEMORY:
          tee_error(od, run_type, overall_line_num,
            E_OUT_OF_MEMORY, CV_FAILED);
          return PARSE_INPUT_ERROR;
        }
        result = fill_y_matrix(od);
        if (result) {
          tee_error(od, run_type, overall_line_num,
            E_OUT_OF_MEMORY, CV_FAILED);
          return PARSE_INPUT_ERROR;
        }
        result = alloc_cv_sdep(od, pc, runs);
        if (result) {
          tee_error(od, run_type, overall_line_num,
            E_OUT_OF_MEMORY, CV_FAILED);
          return PARSE_INPUT_ERROR;
        }
        set_random_seed(od, od->random_seed);
        result = prepare_cv(od, pc, type, groups, runs);
        switch (result) {
          case OUT_OF_MEMORY:
          tee_error(od, run_type, overall_line_num,
            E_OUT_OF_MEMORY, CV_FAILED);
          return PARSE_INPUT_ERROR;
          
          case NOT_ENOUGH_OBJECTS:
          tee_error(od, run_type, overall_line_num,
            E_TOO_FEW_STRUCTURES_FOR_CV, CV_FAILED);
          return PARSE_INPUT_ERROR;
        }
        od->file[ASCII_IN]->name[0] = '\0';
        if ((parameter = get_args(od, "file"))) {
          strcpy(od->file[ASCII_IN]->name, parameter);
          if (!(od->file[ASCII_IN]->handle =
            fopen(od->file[ASCII_IN]->name, "wb+"))) {
            tee_error(od, run_type, overall_line_num,
              E_FILE_CANNOT_BE_OPENED_FOR_WRITING,
              od->file[ASCII_IN]->name, CV_FAILED);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
        }
        if (od->n_proc > 1) {
          result = parallel_cv(od, od->overall_active_x_vars,
             pc, PARALLEL_CV | CV_MODEL, type, groups, runs);
          free_parallel_cv(od, od->mel.thread_info, CV_MODEL, type, runs);
        }
        else {
          if (open_temp_file(od, od->file[TEMP_PRED], "pred_y")) {
            tee_error(od, run_type, overall_line_num,
              E_TEMP_FILE_CANNOT_BE_OPENED_FOR_WRITING,
              od->file[TEMP_PRED]->name, CV_FAILED);
            return PARSE_INPUT_ERROR;
          }
          result = cv(od, pc, CV_MODEL, type, groups, runs);
          if (od->file[TEMP_PRED]->handle) {
            fclose(od->file[TEMP_PRED]->handle);
            od->file[TEMP_PRED]->handle = NULL;
          }
          if (type == LEAVE_MANY_OUT) {
            free_cv_groups(od, runs);
          }
        }
        gettimeofday(&end, NULL);
        elapsed_time(od, &start, &end);
        switch (result) {
          case OUT_OF_MEMORY:
          tee_error(od, run_type, overall_line_num,
            E_OUT_OF_MEMORY, CV_FAILED);
          return PARSE_INPUT_ERROR;

          case CANNOT_WRITE_TEMP_FILE:
          tee_error(od, run_type, overall_line_num,
            E_ERROR_IN_WRITING_TEMP_FILE, "TEMP_FIELD", CV_FAILED);
          return PARSE_INPUT_ERROR;

          case CANNOT_READ_TEMP_FILE:
          tee_error(od, run_type, overall_line_num,
            E_ERROR_IN_READING_TEMP_FILE, "TEMP_FIELD", CV_FAILED);
          return PARSE_INPUT_ERROR;

          case CANNOT_CREATE_THREAD:
          tee_error(od, run_type, overall_line_num,
            E_THREAD_ERROR, "create",
            od->error_code, CV_FAILED);
          return PARSE_INPUT_ERROR;

          case CANNOT_JOIN_THREAD:
          tee_error(od, run_type, overall_line_num,
            E_THREAD_ERROR, "join",
            od->error_code, CV_FAILED);
          return PARSE_INPUT_ERROR;

          default:
          tee_printf(od, M_TOOL_SUCCESS, nesting, command, "CV");
          tee_flush(od);
          od->valid |= CV_BIT;
          result = reload_weights_loadings(od);
          if (result) {
            tee_error(od, run_type, overall_line_num,
              E_RELOAD_WEIGHTS_LOADINGS, CV_FAILED);
            return PARSE_INPUT_ERROR;
          }
        }
      }
      else {
        od->valid |= CV_BIT;
      }
      if (od->file[ASCII_IN]->handle) {
        fclose(od->file[ASCII_IN]->handle);
        od->file[ASCII_IN]->handle = NULL;
      }
    }
    else if (!strcasecmp(arg->me[0], "scramble")) {
      gettimeofday(&start, NULL);
      if (!(od->valid & PLS_BIT)) {
        tee_error(od, run_type, overall_line_num,
          E_NO_PLS_MODEL, SCRAMBLE_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      od->scramble.scramblings = 10;
      if ((parameter = get_args(od, "scramblings"))) {
        sscanf(parameter, "%d", &(od->scramble.scramblings));
      }
      if (od->scramble.scramblings < 1) {
        tee_error(od, run_type, overall_line_num,
          "At least one scrambling operation must be carried out "
          "at each binning level.\n%s",
          SCRAMBLE_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      od->scramble.critical_point = 0.85;
      if ((parameter = get_args(od, "critical_point"))) {
        sscanf(parameter, "%lf", &(od->scramble.critical_point));
      }
      od->scramble.print_runs = 0;
      if ((parameter = get_args(od, "print_runs"))) {
        if (!strncmp(parameter, "y", 1)) {
          od->scramble.print_runs = 1;
        }
      }
      if ((od->scramble.critical_point < 0.0)
        || (od->scramble.critical_point > 1.0))  {
        tee_error(od, run_type, overall_line_num,
          "The critical_point parameter should lie in the range "
          "0.0 - 1.0.\n%s", SCRAMBLE_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      od->scramble.fit_order = 3;
      if ((parameter = get_args(od, "fit_order"))) {
        sscanf(parameter, "%d", &(od->scramble.fit_order));
      }
      if ((od->scramble.fit_order != 2)
        && (od->scramble.fit_order != 3))  {
        tee_error(od, run_type, overall_line_num,
          "The fit_order parameter should be set "
          "either to 2 or 3.\n%s", SCRAMBLE_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      od->scramble.cv_type = LEAVE_ONE_OUT;
      if ((parameter = get_args(od, "type"))) {
        if (!strncasecmp(parameter, "lto", 3)) {
          od->scramble.cv_type = LEAVE_TWO_OUT;
        }
        else if (!strncasecmp(parameter, "lmo", 3)) {
          od->scramble.cv_type = LEAVE_MANY_OUT;
        }
        else if (strncasecmp(parameter, "loo", 3)) {
          tee_error(od, run_type, overall_line_num,
            E_ONLY_LOO_LTO_LMO_CV_ALLOWED, SCRAMBLE_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
      }
      if (!(run_type & DRY_RUN)) {
        get_attr_struct_ave(od, 0, ACTIVE_BIT, &active_struct_num, NULL);
        od->scramble.max_bins = active_struct_num / 3;
        if ((parameter = get_args(od, "max_bins"))) {
          sscanf(parameter, "%d", &(od->scramble.max_bins));
        }
        if ((od->scramble.max_bins > (active_struct_num / 3))
          || (od->scramble.max_bins < 2))  {
          tee_error(od, run_type, overall_line_num,
            "The max_bins parameter should lie in the range "
            "2 - %d.\n%s", active_struct_num / 3,
            SCRAMBLE_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        od->scramble.min_bins = 2;
        if ((parameter = get_args(od, "min_bins"))) {
          sscanf(parameter, "%d", &(od->scramble.min_bins));
        }
        if ((od->scramble.min_bins < 2)
          || (od->scramble.min_bins > od->scramble.max_bins))  {
          tee_error(od, run_type, overall_line_num,
            "The min_bins parameter should lie in the range "
            "2 - %d.\n%s", od->scramble.max_bins,
            SCRAMBLE_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        pc = od->pc_num;
        if ((parameter = get_args(od, "pc"))) {
          sscanf(parameter, "%d", &pc);
        }
        if (pc < 1) {
          tee_error(od, run_type, overall_line_num,
            E_AT_LEAST_ONE_PC, SCRAMBLE_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        if (pc > od->pc_num) {
          tee_error(od, run_type, overall_line_num,
            E_TOO_FEW_MANY_FOR_CURRENT_PLS_MODEL,
            "many PCs", SCRAMBLE_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        ++command;
        tee_printf(od, M_TOOL_INVOKE, nesting, command, "SCRAMBLE", line_orig);
        if (od->scramble.cv_type == LEAVE_MANY_OUT) {
          result = check_lmo_parameters(od, SCRAMBLE_FAILED,
            &(od->scramble.groups), &(od->scramble.runs),
            run_type, overall_line_num);
          if (result) {
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
        }
        tee_flush(od);
        if ((pc > od->pc_num)
          || (pc > od->overall_active_x_vars)) {
          tee_error(od, run_type, overall_line_num,
            E_TOO_FEW_MANY_FOR_AVAILABLE_DATA, "many PCs",
            SCRAMBLE_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        if (open_temp_file(od, od->file[TEMP_SCRAMBLE], "scramble")) {
          tee_error(od, run_type, overall_line_num,
            E_TEMP_FILE_CANNOT_BE_OPENED_FOR_WRITING,
            od->file[TEMP_SCRAMBLE]->name, SCRAMBLE_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        result = scramble(od, pc);
        gettimeofday(&end, NULL);
        elapsed_time(od, &start, &end);
        switch (result) {
          case OUT_OF_MEMORY:
          tee_error(od, run_type, overall_line_num,
            E_OUT_OF_MEMORY, SCRAMBLE_FAILED);
          return PARSE_INPUT_ERROR;

          case CANNOT_WRITE_TEMP_FILE:
          tee_error(od, run_type, overall_line_num,
            E_ERROR_IN_WRITING_TEMP_FILE, "TEMP_FIELD", SCRAMBLE_FAILED);
          return PARSE_INPUT_ERROR;

          case CANNOT_READ_TEMP_FILE:
          tee_error(od, run_type, overall_line_num,
            E_ERROR_IN_READING_TEMP_FILE, "TEMP_FIELD", SCRAMBLE_FAILED);
          return PARSE_INPUT_ERROR;

          case CANNOT_CREATE_THREAD:
          tee_error(od, run_type, overall_line_num, E_THREAD_ERROR, "create",
            od->error_code, SCRAMBLE_FAILED);
          return PARSE_INPUT_ERROR;

          case CANNOT_JOIN_THREAD:
          tee_error(od, run_type, overall_line_num, E_THREAD_ERROR, "join",
            od->error_code, SCRAMBLE_FAILED);
          return PARSE_INPUT_ERROR;

          case NOT_ENOUGH_OBJECTS:
          tee_error(od, run_type, overall_line_num,
            E_TOO_FEW_STRUCTURES_FOR_CV, SCRAMBLE_FAILED);
          return PARSE_INPUT_ERROR;

          default:
          if (od->file[TEMP_SCRAMBLE]->handle) {
            fclose(od->file[TEMP_SCRAMBLE]->handle);
            od->file[TEMP_SCRAMBLE]->handle = NULL;
          }
          result = reload_weights_loadings(od);
          if (result) {
            tee_error(od, run_type, overall_line_num,
              E_RELOAD_WEIGHTS_LOADINGS, SCRAMBLE_FAILED);
            return PARSE_INPUT_ERROR;
          }
          tee_printf(od, M_TOOL_SUCCESS, nesting, command, "SCRAMBLE");
          tee_flush(od);
        }
      }
      else {
        od->valid |= SCRAMBLE_BIT;
      }
    }
    else if (!strcasecmp(arg->me[0], "ffdsel")) {
      gettimeofday(&start, NULL);
      if (!(od->valid & PLS_BIT)) {
        tee_error(od, run_type, overall_line_num,
          E_NO_PLS_MODEL, FFDSEL_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      od->ffdsel.cv_type = LEAVE_ONE_OUT;
      if ((parameter = get_args(od, "type"))) {
        if (!strncasecmp(parameter, "ext", 3)) {
          od->ffdsel.cv_type = EXTERNAL_PREDICTION;
        }
        else if (!strncasecmp(parameter, "lto", 3)) {
          od->ffdsel.cv_type = LEAVE_TWO_OUT;
        }
        else if (!strncasecmp(parameter, "lmo", 3)) {
          od->ffdsel.cv_type = LEAVE_MANY_OUT;
        }
        else if (strncasecmp(parameter, "loo", 3)) {
          tee_error(od, run_type, overall_line_num,
            "Only \"EXTERNAL\", \"LOO\", \"LTO\", \"LMO\" "
            "keywords are allowed as validation types.\n%s",
            FFDSEL_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
      }
      od->ffdsel.percent_dummies = 20;
      if ((parameter = get_args(od, "percent_dummies"))) {
        sscanf(parameter, "%d", &(od->ffdsel.percent_dummies));
        if ((od->ffdsel.percent_dummies < 0)
          || (od->ffdsel.percent_dummies > 50)) {
          tee_error(od, run_type, overall_line_num,
            "The percentage of dummy variables "
            "should lie in the range 0 - 50%%.\n%s",
            FFDSEL_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
      }
      od->ffdsel.use_srd_groups = 0;
      if (od->voronoi_num) {
        od->ffdsel.use_srd_groups = 1;
      }
      if ((parameter = get_args(od, "use_srd_groups"))) {
        od->ffdsel.use_srd_groups = 0;
        if (!strncasecmp(parameter, "y", 1)) {
          if (!(od->valid & SEED_BIT)) {
            tee_error(od, run_type, overall_line_num,
              "No variable group information is available. "
              "Please run SRD tool first.\n%s",
              FFDSEL_FAILED);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
          od->ffdsel.use_srd_groups = 1;
        }
      }
      od->ffdsel.retain_uncertain = 1;
      if ((parameter = get_args(od, "retain_uncertain"))) {
        od->ffdsel.retain_uncertain = 0;
        if (!strncasecmp(parameter, "y", 1)) {
          od->ffdsel.retain_uncertain = 1;
        }
      }
      od->ffdsel.fold_over = 0;
      if ((parameter = get_args(od, "fold_over"))) {
        if (!strncasecmp(parameter, "y", 1)) {
          od->ffdsel.fold_over = 1;
        }
      }
      od->ffdsel.print_sdep = 0;
      if ((parameter = get_args(od, "print_sdep"))) {
        if (!strncasecmp(parameter, "y", 1)) {
          od->ffdsel.print_sdep = 1;
        }
      }
      od->ffdsel.print_effect = 0;
      if ((parameter = get_args(od, "print_effect"))) {
        if (!strncasecmp(parameter, "y", 1)) {
          od->ffdsel.print_effect = 1;
        }
      }
      od->ffdsel.combination_variable_ratio = 2.0;
      if ((parameter = get_args(od, "combination_variable_ratio"))) {
        sscanf(parameter, "%lf", &(od->ffdsel.combination_variable_ratio));
      }
      if ((od->ffdsel.combination_variable_ratio < 1.0)
        || (od->ffdsel.combination_variable_ratio > 10.0)) {
        tee_error(od, run_type, overall_line_num,
          "The combination/variable ratio "
          "should lie in the range 1.0 - 10.0.\n%s",
          FFDSEL_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      od->ffdsel.confidence_level = 99;
      if ((parameter = get_args(od, "confidence_level"))) {
        sscanf(parameter, "%d", &(od->ffdsel.confidence_level));
      }
      if ((od->ffdsel.confidence_level < 80.0)
        || (od->ffdsel.confidence_level > 99.0)) {
        tee_error(od, run_type, overall_line_num,
          "The confidence level should lie in the range 80 - 99.\n%s",
          FFDSEL_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      if (!(run_type & DRY_RUN)) {
        pc = od->pc_num;
        if ((parameter = get_args(od, "pc"))) {
          sscanf(parameter, "%d", &pc);
        }
        if (pc < 1) {
          tee_error(od, run_type, overall_line_num,
            E_AT_LEAST_ONE_PC, FFDSEL_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        if (pc > od->pc_num) {
          tee_error(od, run_type, overall_line_num,
            E_TOO_FEW_MANY_FOR_CURRENT_PLS_MODEL,
            "many PCs", FFDSEL_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        if ((od->ffdsel.cv_type == EXTERNAL_PREDICTION)
          && (od->ext_pred_object_num < 4)) {
          tee_error(od, run_type, overall_line_num,
            "A minimum of 4 test set compounds are necessary "
            "to use external validation.\n%s",
            FFDSEL_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        if (od->ffdsel.cv_type == LEAVE_MANY_OUT) {
          result = check_lmo_parameters(od, FFDSEL_FAILED,
            &(od->ffdsel.groups), &(od->ffdsel.runs),
            run_type, overall_line_num);
          if (result) {
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
        }
        ++command;
        tee_printf(od, M_TOOL_INVOKE, nesting, command, "FFDSEL", line_orig);
        tee_flush(od);
        result = ffdsel(od, pc);
        gettimeofday(&end, NULL);
        elapsed_time(od, &start, &end);
        switch (result) {

          case OUT_OF_MEMORY:
          tee_error(od, run_type, overall_line_num,
            E_OUT_OF_MEMORY, FFDSEL_FAILED);
          return PARSE_INPUT_ERROR;

          case CANNOT_CREATE_THREAD:
          tee_error(od, run_type, overall_line_num,
            E_THREAD_ERROR, "create",
            od->error_code, FFDSEL_FAILED);
          return PARSE_INPUT_ERROR;

          case CANNOT_JOIN_THREAD:
          tee_error(od, run_type, overall_line_num,
            E_THREAD_ERROR, "join",
            od->error_code, FFDSEL_FAILED);
          return PARSE_INPUT_ERROR;

          case NOT_ENOUGH_OBJECTS:
          tee_error(od, run_type, overall_line_num,
            E_TOO_FEW_STRUCTURES_FOR_CV, FFDSEL_FAILED);
          return PARSE_INPUT_ERROR;

          default:
          tee_printf(od, M_TOOL_SUCCESS, nesting, command, "FFDSEL");
          tee_flush(od);
          result = reload_weights_loadings(od);
          if (result) {
            tee_error(od, run_type, overall_line_num,
              E_RELOAD_WEIGHTS_LOADINGS, FFDSEL_FAILED);
            return PARSE_INPUT_ERROR;
          }
          od->valid |= FFDSEL_BIT;
        }
      }
      else {
        od->valid |= FFDSEL_BIT;
      }
    }
    else if (!strcasecmp(arg->me[0], "uvepls")) {
      gettimeofday(&start, NULL);
      if (!(od->valid & PLS_BIT)) {
        tee_error(od, run_type, overall_line_num,
          E_NO_PLS_MODEL, UVEPLS_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      od->uvepls.cv_type = LEAVE_ONE_OUT;
      od->uvepls.groups = 1;
      od->uvepls.runs = od->active_object_num;
      od->uvepls.use_srd_groups = 0;
      if (od->voronoi_num) {
        od->uvepls.use_srd_groups = 1;
      }
      od->uvepls.save_ram = 0;
      if ((parameter = get_args(od, "save_ram"))) {
        if (!strncasecmp(parameter, "y", 1)) {
          od->uvepls.save_ram = 1;
        }
      }
      if ((parameter = get_args(od, "use_srd_groups"))) {
        od->uvepls.use_srd_groups = 0;
        if (!strncasecmp(parameter, "y", 1)) {
          if (!(od->valid & SEED_BIT)) {
            tee_error(od, run_type, overall_line_num,
              E_NO_VALID_SELECTION, "SRD", UVEPLS_FAILED);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
          od->uvepls.use_srd_groups = 1;
        }
      }
      od->uvepls.dummy_range_coefficient = 1.0;
      if ((parameter = get_args(od, "dummy_range_coefficient"))) {
        sscanf(parameter, "%lf", &(od->uvepls.dummy_range_coefficient));
        if ((od->uvepls.dummy_range_coefficient < 1.0)
          || (od->uvepls.dummy_range_coefficient > 5.0)) {
          tee_error(od, run_type, overall_line_num,
            "The dummy_range_coefficient parameter "
            "should lie in the range 0.0 - 5.0.\n%s",
            UVEPLS_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
      }
      od->uvepls.dummy_value_coefficient = 1.0e-10;
      if ((parameter = get_args(od, "dummy_value_coefficient"))) {
        sscanf(parameter, "%lf", &(od->uvepls.dummy_value_coefficient));
        if ((od->uvepls.dummy_value_coefficient < 1.0e-20)
          || (od->uvepls.dummy_value_coefficient > 1.0)) {
          tee_error(od, run_type, overall_line_num,
            "The dummy_value_coefficient parameter should lie "
            "in the range 1.0e-20 - 1.0.\n%s",
            UVEPLS_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
      }

      od->uvepls.uve_m = 0;
      if ((parameter = get_args(od, "uve_m"))) {
        if (!strncmp(parameter, "y", 1)) {
          od->uvepls.uve_m = 1;
        }
      }

      od->uvepls.uve_alpha = 0.0;
      if ((parameter = get_args(od, "uve_alpha"))) {
        sscanf(parameter, "%lf", &(od->uvepls.uve_alpha));
        if ((od->uvepls.uve_alpha < 0.0)
          || (od->uvepls.uve_alpha > 100.0)) {
          tee_error(od, run_type, overall_line_num,
            "The quantile for UVE-ALPHA should lie "
            "in the range 0.0 - 100.0.\n%s",
            UVEPLS_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
      }

      od->uvepls.ive = 0;
      if ((parameter = get_args(od, "ive"))) {
        if (!strncasecmp(parameter, "y", 1)) {
          od->uvepls.ive = 1;
        }
      }

      od->uvepls.ive_percent_limit = 100;
      if ((parameter = get_args(od, "ive_percent_limit"))) {
        sscanf(parameter, "%d", &(od->uvepls.ive_percent_limit));
        if ((od->uvepls.ive_percent_limit < 0)
          || (od->uvepls.ive_percent_limit > 100)) {
          tee_error(od, run_type, overall_line_num,
            "The maximum percentage of IVE-eliminated variables "
            "should lie in the range 0 - 100.\n%s",
            UVEPLS_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
      }

      if (!(run_type & DRY_RUN)) {
        if ((parameter = get_args(od, "type"))) {
          if (!strncasecmp(parameter, "lto", 3)) {
            od->uvepls.cv_type = LEAVE_TWO_OUT;
            od->uvepls.groups = 1;
            od->uvepls.runs = 0;
            for (i = 0; i < (od->active_object_num - 1); ++i) {
              for (j = i + 1; j < od->active_object_num; ++j) {
                ++(od->uvepls.runs);
              }
            }
          }
          else if (!strncasecmp(parameter, "lmo", 3)) {
            od->uvepls.cv_type = LEAVE_MANY_OUT;
            result = check_lmo_parameters(od, UVEPLS_FAILED,
              &(od->uvepls.groups), &(od->uvepls.runs),
              run_type, overall_line_num);
            if (result) {
              fail = !(run_type & INTERACTIVE_RUN);
              continue;
            }
          }
          else if (strncasecmp(parameter, "loo", 3)) {
            tee_error(od, run_type, overall_line_num,
              E_ONLY_LOO_LTO_LMO_CV_ALLOWED, UVEPLS_FAILED);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
        }
        pc = od->pc_num;
        if ((parameter = get_args(od, "pc"))) {
          sscanf(parameter, "%d", &pc);
        }
        if (pc < 1) {
          tee_error(od, run_type, overall_line_num,
            E_AT_LEAST_ONE_PC, UVEPLS_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        if (pc > od->pc_num) {
          tee_error(od, run_type, overall_line_num,
            E_TOO_FEW_MANY_FOR_CURRENT_PLS_MODEL,
            "many PCs", UVEPLS_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        od->uvepls.ive_external_sdep = 0;
        if ((parameter = get_args(od, "ive_external_pred"))) {
          if (!strncasecmp(parameter, "y", 1)) {
            if (od->ext_pred_object_num < 4) {
              tee_error(od, run_type, overall_line_num,
                "A minimum of 4 test set compounds are necessary "
                "to use external validation.\n%s",
                UVEPLS_FAILED);
              fail = !(run_type & INTERACTIVE_RUN);
              continue;
            }
            od->uvepls.ive_external_sdep = 1;
          }
        }

        if (od->uvepls.save_ram) {
          if (open_temp_file(od, od->file[TEMP_CV_COEFF], "cv_coeff")) {
            tee_error(od, run_type, overall_line_num,
              E_TEMP_FILE_CANNOT_BE_OPENED_FOR_WRITING,
              od->file[TEMP_CV_COEFF]->name, UVEPLS_FAILED);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
        }
        /*
        create a temporary file where the composition of each
        reduced model is stored
        */
        if (open_temp_file(od, od->file[TEMP_UVEPLS], "uvepls")) {
          tee_error(od, run_type, overall_line_num,
            E_TEMP_FILE_CANNOT_BE_OPENED_FOR_WRITING,
            od->file[TEMP_UVEPLS]->name, UVEPLS_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        ++command;
        tee_printf(od, M_TOOL_INVOKE, nesting, command, "UVEPLS", line_orig);
        tee_flush(od);
        result = uvepls(od, pc);
        if (od->file[TEMP_CV_COEFF]->handle) {
          fclose(od->file[TEMP_CV_COEFF]->handle);
          od->file[TEMP_CV_COEFF]->handle = NULL;
        }
        gettimeofday(&end, NULL);
        elapsed_time(od, &start, &end);
        switch (result) {
          case OUT_OF_MEMORY:
          tee_error(od, run_type, overall_line_num,
            E_OUT_OF_MEMORY, UVEPLS_FAILED);
          return PARSE_INPUT_ERROR;

          case CANNOT_CREATE_THREAD:
          tee_error(od, run_type, overall_line_num,
            E_THREAD_ERROR, "create", od->error_code,
            UVEPLS_FAILED);
          return PARSE_INPUT_ERROR;

          case CANNOT_JOIN_THREAD:
          tee_error(od, run_type, overall_line_num,
            E_THREAD_ERROR, "join", od->error_code,
            UVEPLS_FAILED);
          return PARSE_INPUT_ERROR;

          case NOT_ENOUGH_OBJECTS:
          tee_error(od, run_type, overall_line_num,
            E_TOO_FEW_STRUCTURES_FOR_CV, UVEPLS_FAILED);
          return PARSE_INPUT_ERROR;

          default:
          if (od->file[TEMP_UVEPLS]->handle) {
            fclose(od->file[TEMP_UVEPLS]->handle);
            od->file[TEMP_UVEPLS]->handle = NULL;
          }
          result = reload_weights_loadings(od);
          if (result) {
            tee_error(od, run_type, overall_line_num,
              E_RELOAD_WEIGHTS_LOADINGS, UVEPLS_FAILED);
            return PARSE_INPUT_ERROR;
          }
          tee_printf(od, M_TOOL_SUCCESS, nesting, command, "UVEPLS");
          tee_flush(od);
          od->valid |= UVEPLS_BIT;
        }
      }
      else {
        od->valid |= UVEPLS_BIT;
      }
    }
    else if (!strcasecmp(arg->me[0], "predict")) {
      gettimeofday(&start, NULL);
      if (!(od->valid & PLS_BIT)) {
        tee_error(od, run_type, overall_line_num,
          E_NO_PLS_MODEL, PREDICT_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      if (!(run_type & DRY_RUN)) {
        pc = od->pc_num;
        if ((parameter = get_args(od, "pc"))) {
          sscanf(parameter, "%d", &pc);
        }
        if (pc < 1) {
          tee_error(od, run_type, overall_line_num,
            E_AT_LEAST_ONE_PC, PREDICT_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        if (pc > od->pc_num) {
          tee_error(od, run_type, overall_line_num,
            E_TOO_FEW_MANY_FOR_CURRENT_PLS_MODEL,
            "many PCs", PREDICT_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        if (!(od->ext_pred_object_num)) {
          tee_error(od, run_type, overall_line_num,
            "No objects have been extracted for prediction.\n%s",
            PREDICT_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        result = open_temp_file(od, od->file[TEMP_EXT_PRED], "ext_pred_y");
        if (result) {
          tee_error(od, run_type, overall_line_num,
            E_TEMP_FILE_CANNOT_BE_OPENED_FOR_WRITING,
            od->file[TEMP_EXT_PRED]->name, PREDICT_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        od->file[ASCII_IN]->name[0] = '\0';
        if ((parameter = get_args(od, "file"))) {
          strcpy(od->file[ASCII_IN]->name, parameter);
          if (!(od->file[ASCII_IN]->handle =
            fopen(od->file[ASCII_IN]->name, "wb+"))) {
            tee_error(od, run_type, overall_line_num,
              E_FILE_CANNOT_BE_OPENED_FOR_WRITING,
              od->file[ASCII_IN]->name, PREDICT_FAILED);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
        }
        ++command;
        tee_printf(od, M_TOOL_INVOKE, nesting, command, "PREDICT", line_orig);
        tee_flush(od);
        result = predict(od, pc);
        gettimeofday(&end, NULL);
        elapsed_time(od, &start, &end);
        switch (result) {
          case OUT_OF_MEMORY:
          tee_error(od, run_type, overall_line_num,
            E_OUT_OF_MEMORY, PREDICT_FAILED);
          return PARSE_INPUT_ERROR;

          case CANNOT_READ_TEMP_FILE:
          tee_error(od, run_type, overall_line_num,
            E_ERROR_IN_READING_TEMP_FILE, "TEMP_FIELD", PREDICT_FAILED);
          return PARSE_INPUT_ERROR;

          case CANNOT_WRITE_TEMP_FILE:
          tee_error(od, run_type, overall_line_num,
            E_ERROR_IN_WRITING_TEMP_FILE, "TEMP_FIELD", PREDICT_FAILED);
          return PARSE_INPUT_ERROR;

          default:
          if ((parameter = get_args(od, "scores"))
            && (!strncasecmp(parameter, "X", 1))) {
            print_pls_scores(od, X_SCORES | PREDICT_BIT);
          }
          tee_printf(od, M_TOOL_SUCCESS, nesting, command, "PREDICT");
          tee_flush(od);
          od->valid |= PREDICT_BIT;
        }
        if (od->file[TEMP_EXT_PRED]->handle) {
          fclose(od->file[TEMP_EXT_PRED]->handle);
          od->file[TEMP_EXT_PRED]->handle = NULL;
        }
      }
      else {
        od->valid |= PREDICT_BIT;
      }
      if (od->file[ASCII_IN]->handle) {
        fclose(od->file[ASCII_IN]->handle);
        od->file[ASCII_IN]->handle = NULL;
      }
    }
    else if (!strcasecmp(arg->me[0], "d_optimal")) {
      gettimeofday(&start, NULL);
      if (!(od->valid & PLS_BIT)) {
        tee_error(od, run_type, overall_line_num,
          E_NO_PLS_MODEL, D_OPTIMAL_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      if (!(run_type & DRY_RUN)) {
        candidates = od->mal.x_weights->m;
        design_points = 0;
        synonym = 0;
        pc = od->pc_num;
        if ((parameter = get_args(od, "pc"))) {
          sscanf(parameter, "%d", &pc);
        }
        if (pc < 1) {
          tee_error(od, run_type, overall_line_num,
            E_AT_LEAST_ONE_PC, D_OPTIMAL_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        if (pc > od->pc_num) {
          tee_error(od, run_type, overall_line_num,
            E_TOO_FEW_MANY_FOR_CURRENT_PLS_MODEL,
            "many PCs", D_OPTIMAL_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
      }
      design_points = (int)safe_rint(0.5 * (double)candidates);
      if ((parameter = get_args(od, "percent_remove"))) {
        ++synonym;
        if (!(run_type & DRY_RUN)) {
          sscanf(parameter, "%d", &percent_remove);
          if ((percent_remove < 1) || (percent_remove > 90)) {
            tee_error(od, run_type, overall_line_num,
              "The percent of variables which D_OPTIMAL "
              "will remove should be in the range 1 - 90%%.\n%s",
              D_OPTIMAL_FAILED);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
          design_points = (int)safe_rint((double)(100 - percent_remove)
            / 100.0 * (double)candidates);
        }
      }
      if ((parameter = get_args(od, "design_points"))) {
        ++synonym;
        if (!(run_type & DRY_RUN)) {
          sscanf(parameter, "%d", &design_points);
          if ((design_points < (int)(0.1 * (double)candidates))
            || (design_points > (candidates - 1))) {
            tee_error(od, run_type, overall_line_num,
              "The number of variables which D_OPTIMAL "
              "will keep should be in the range %d - %d.\n%s",
              (int)(0.1 * (double)candidates),
              candidates - 1, D_OPTIMAL_FAILED);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
        }
      }
      if (synonym != 1) {
        tee_error(od, run_type, overall_line_num,
          "Please specify either the percent_remove "
          "or design_points keyword.\n%s",
          D_OPTIMAL_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      type = WEIGHTS;
      if ((parameter = get_args(od, "type"))) {
        if (!strncasecmp(parameter, "loading", 7)) {
          type = LOADINGS;
        }
        else if (strncasecmp(parameter, "weight", 6)) {
          tee_error(od, run_type, overall_line_num,
            "Only \"WEIGHTS\" and \"LOADINGS\" types are allowed.\n%s", 
            D_OPTIMAL_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
      }
      if (!(run_type & DRY_RUN)) {
        ++command;
        tee_printf(od, M_TOOL_INVOKE, nesting, command, "D_OPTIMAL", line_orig);
        tee_flush(od);
        result = d_optimal_select(od, pc,
          design_points, type);
        gettimeofday(&end, NULL);
        elapsed_time(od, &start, &end);
        switch (result) {

          case OUT_OF_MEMORY:
          tee_error(od, run_type, overall_line_num,
            E_OUT_OF_MEMORY, D_OPTIMAL_FAILED);
          return PARSE_INPUT_ERROR;

          case SINGULAR_MATRIX:
          tee_error(od, run_type, overall_line_num,
            "The D-optimal algorithm ran into a singular matrix; "
            "please consider a less aggressive variable reduction.\n%s",
            D_OPTIMAL_FAILED);
          return PARSE_INPUT_ERROR;

          case CANNOT_READ_TEMP_FILE:
          tee_error(od, run_type, overall_line_num,
            E_ERROR_IN_READING_TEMP_FILE, "TEMP_FIELD", D_OPTIMAL_FAILED);
          return PARSE_INPUT_ERROR;

          default:
          tee_printf(od, M_TOOL_SUCCESS, nesting, command, "D_OPTIMAL");
          tee_flush(od);
          od->valid |= D_OPTIMAL_BIT;
        }
      }
      else {
        od->valid |= D_OPTIMAL_BIT;
      }
    }
    else if (!strcasecmp(arg->me[0], "srd")) {
      gettimeofday(&start, NULL);
      if (!(od->valid & PLS_BIT)) {
        tee_error(od, run_type, overall_line_num, E_NO_PLS_MODEL, SRD_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      type = WEIGHTS;
      if ((parameter = get_args(od, "type"))) {
        if (!strncasecmp(parameter, "loading", 7)) {
          type = LOADINGS;
        }
        else if (strncasecmp(parameter, "weight", 6)) {
          tee_error(od, run_type, overall_line_num,
            "Only \"WEIGHTS\" and \"LOADINGS\" "
            "types are allowed.\n%s",
            SRD_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
      }
      if (!(run_type & DRY_RUN)) {
        pc = od->pc_num;
        if ((parameter = get_args(od, "pc"))) {
          sscanf(parameter, "%d", &pc);
        }
        if (pc < 1) {
          tee_error(od, run_type, overall_line_num,
            E_AT_LEAST_ONE_PC, SRD_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        if (pc > od->pc_num) {
          tee_error(od, run_type, overall_line_num,
            E_TOO_FEW_MANY_FOR_CURRENT_PLS_MODEL,
            "many PCs", SRD_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        default_seeds[0] = od->x_vars / 10;
        default_seeds[1] = od->overall_active_x_vars / 2;
        default_seeds[2] = 3000;
        qsort(default_seeds, 3, sizeof(int), compare_integers);
        seeds = default_seeds[0];
        for (i = 0; i < 3; ++i) {
          grid_size[i] = od->grid.end_coord[i]
            - od->grid.start_coord[i];
        }
        grid_diagonal = sqrt(square(grid_size[0])
          + square(grid_size[1]) + square(grid_size[2]));
        if ((parameter = get_args(od, "seeds"))) {
          sscanf(parameter, "%d", &seeds);
          if ((seeds < 1) || (seeds > (od->overall_active_x_vars))) {
            tee_error(od, run_type, overall_line_num,
              "The number of seeds should be in the range 1 - %d.\n%s",
              od->overall_active_x_vars, SRD_FAILED);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
        }
        critical_distance = 1.0;
        if ((parameter = get_args(od, "critical_distance"))) {
          sscanf(parameter, "%lf", &critical_distance);
          if ((critical_distance < 0.0)
            || (critical_distance > grid_diagonal)) {
            tee_error(od, run_type, overall_line_num,
              "The critical distance should be in the range "
              "0.00 - %.2lf angstrom.\n%s",
              grid_diagonal, SRD_FAILED);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
        }
        srd_collapse = 0;
        if ((parameter = get_args(od, "collapse"))) {
          if (!strncasecmp(parameter, "y", 1)) {
            srd_collapse = 1;
            collapse_distance = 2.0;
            if ((parameter = get_args(od, "collapse_distance"))) {
              sscanf(parameter, "%lf", &collapse_distance);
              if ((collapse_distance < 0.0)
                || (collapse_distance > grid_diagonal)) {
                tee_error(od, run_type, overall_line_num,
                  "The collapse distance should be in the range "
                  "0.00 - %.2lf angstrom.\n%s",
                  grid_diagonal, SRD_FAILED);
                fail = !(run_type & INTERACTIVE_RUN);
                continue;
              }
            }
          }
        }
        ++command;
        tee_printf(od, M_TOOL_INVOKE, nesting, command, "SRD", line_orig);
        tee_flush(od);
        result = srd(od, pc, seeds, type, srd_collapse,
          critical_distance, collapse_distance);
        gettimeofday(&end, NULL);
        elapsed_time(od, &start, &end);
        switch (result) {
          case SINGULAR_MATRIX:
          tee_error(od, run_type, overall_line_num,
            "The D-optimal algorithm ran into a singular matrix; "
            "please consider increasing the number of seeds.\n%s",
            SRD_FAILED);
          return PARSE_INPUT_ERROR;

          case OUT_OF_MEMORY:
          tee_error(od, run_type, overall_line_num,
            E_OUT_OF_MEMORY, SRD_FAILED);
          return PARSE_INPUT_ERROR;

          default:
          tee_printf(od, M_TOOL_SUCCESS, nesting, command, "SRD");
          tee_flush(od);
          od->valid |= SEED_BIT;
        }
      }
      else {
        od->valid |= SEED_BIT;
      }
    }
    else if (!strcasecmp(arg->me[0], "plot")) {
      gettimeofday(&start, NULL);
      if (!(parameter = get_args(od, "type"))) {
        tee_error(od, run_type, overall_line_num,
          "Please select which type of plot "
          "you wish to generate.\n%s",
          PLOT_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      type = 0;
      max_pc = od->pc_num;
      if (!strncasecmp(parameter, "recalc_vs_exp", 13)) {
        type = PLS_PLOT | VS_EXP_PLOT | RECALC_VS_EXP;
      }
      else if (!strncasecmp(parameter, "pred_vs_exp", 11)) {
        if (!(od->valid & CV_BIT)) {
          tee_error(od, run_type, overall_line_num,
            E_NO_CV_MODEL, PLOT_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        max_pc = od->cv.pc_num;
        type = PLS_PLOT | VS_EXP_PLOT | PRED_VS_EXP;
      }
      else if (!strncasecmp(parameter, "ext_pred_vs_exp", 15)) {
        if (!(od->valid & PREDICT_BIT)) {
          tee_error(od, run_type, overall_line_num,
            E_NO_EXTERNAL_PREDICTION, PLOT_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        type = PLS_PLOT | VS_EXP_PLOT | EXT_PRED_VS_EXP;
      }
      else if (!strncasecmp(parameter, "scrambled_q2_vs_r2", 18)) {
        if (!(od->valid & SCRAMBLE_BIT)) {
          tee_error(od, run_type, overall_line_num,
            E_NO_SCRAMBLE_ANALYSIS, PLOT_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        type = SCRAMBLED_Q2_VS_R2;
      }
      else if (!strncasecmp(parameter, "scrambled_secv_vs_r2", 20)) {
        if (!(od->valid & SCRAMBLE_BIT)) {
          tee_error(od, run_type, overall_line_num,
            E_NO_SCRAMBLE_ANALYSIS, PLOT_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
        type = SCRAMBLED_SECV_VS_R2;
      }
      else if (!strncasecmp(parameter, "sdec", 4)) {
        type = PLS_PLOT | SDEC_VS_PC_PLOT;
      }
      else if (!strncasecmp(parameter, "r2", 2)) {
        type = PLS_PLOT | R2_VS_PC_PLOT;
      }
      else if (!strncasecmp(parameter, "sdep", 4)) {
        type = PLS_PLOT | SDEP_VS_PC_PLOT;
      }
      else if (!strncasecmp(parameter, "q2", 2)) {
        type = PLS_PLOT | Q2_VS_PC_PLOT;
      }
      else if (!strncasecmp(parameter, "pls_x_vs_y", 10)) {
        type = PLS_PLOT | PLS_X_VS_Y_SCORES_PLOT;
      }
      else if (!strncasecmp(parameter, "pls_loading", 11)) {
        type = PLS_PLOT | LOADINGS_PLOT;
      }
      else if (!strncasecmp(parameter, "pls_weight", 10)) {
        type = PLS_PLOT | WEIGHTS_PLOT;
      }
      else if (!strncasecmp(parameter, "pls_score", 9)) {
        type = PLS_PLOT | SCORES_PLOT;
      }
      else if (!strncasecmp(parameter, "pca_loading", 11)) {
        type = PCA_PLOT | LOADINGS_PLOT;
      }
      else if (!strncasecmp(parameter, "pca_score", 9)) {
        type = PCA_PLOT | SCORES_PLOT;
      }
      if (!type) {
        tee_error(od, run_type, overall_line_num,
          "Only \"RECALC_VS_EXP\", \"PRED_VS_EXP\", \"RECALC_RES\", "
          "\"PRED_RES\", \"SDEC_VS_PC\", \"SDEP_VS_PC\", \"R2_VS_PC\", "
          "\"Q2_VS_PC\", \"PLS_X_VS_Y_SCORES\", "
          "\"PLS_LOADINGS\", \"PLS_WEIGHTS\", \"PLS_SCORES\", "
          "\"PCA_LOADINGS\", \"PCA_SCORES\", "
          "\"SCRAMBLED_Q2_VS_R2\" and \"SCRAMBLED_SECV_VS_R2\" "
          "plots are supported.\n%s",
          PLOT_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      if (type & PLS_PLOT) {
        if (!(od->valid & PLS_BIT)) {
          tee_error(od, run_type, overall_line_num,
            E_NO_PLS_MODEL, PLOT_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
      }
      else if (type & PCA_PLOT) {
        if (!(od->valid & PCA_BIT)) {
          tee_error(od, run_type, overall_line_num,
            E_NO_PCA_ANALYSIS, PLOT_FAILED);
          fail = !(run_type & INTERACTIVE_RUN);
          continue;
        }
      }
      if (!(parameter = get_args(od, "file"))) {
        tee_error(od, run_type, overall_line_num,
          "Please specify a file where plot data should be saved.\n%s",
          PLOT_FAILED);
        fail = !(run_type & INTERACTIVE_RUN);
        continue;
      }
      memset(file_basename, 0, BUF_LEN);
      strncpy(file_basename, parameter, BUF_LEN - 100);
      start_gnuplot = (od->gnuplot.use_gnuplot ? run_type & INTERACTIVE_RUN : 0);
      if (!(run_type & DRY_RUN)) {
        if ((parameter = get_args(od, "residuals"))) {
          if (!strncasecmp(parameter, "y", 1)) {
            type |= RESIDUALS;
          }
        }
        label = 0;
        if ((parameter = get_args(od, "label"))) {
          if (!strncasecmp(parameter, "number", 6)) {
            label = NUMBER_LABEL;
          }
          else if (!strncasecmp(parameter, "id", 2)) {
            label = ID_LABEL;
          }
          else if (!strncasecmp(parameter, "name", 4)) {
            label = NAME_LABEL;
          }
          else if (strncasecmp(parameter, "none", 4)) {
            tee_error(od, run_type, overall_line_num,
              "Only \"NAME\", \"NUMBER\" AND \"NONE\", "
              "keywords are accepted.\n%s", PLOT_FAILED);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
        }
        if ((type & VS_EXP_PLOT) || (type & PLS_X_VS_Y_SCORES_PLOT)) {

          pc = max_pc;
          if ((parameter = get_args(od, "pc"))) {
            sscanf(parameter, "%d", &pc);
          }
          if (pc < 1) {
            tee_error(od, run_type, overall_line_num,
              E_AT_LEAST_ONE_PC, PLOT_FAILED);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
          if (pc > max_pc) {
            tee_error(od, run_type, overall_line_num,
              E_TOO_FEW_MANY_FOR_CURRENT_PLS_MODEL,
              "many PCs", PLOT_FAILED);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
          strcpy(comma_hyphen_list, "all");
          if ((parameter = get_args(od, "y_var_list"))) {
            strcpy(comma_hyphen_list, parameter);
          }
          result = parse_comma_hyphen_list_to_array
            (od, comma_hyphen_list, Y_VAR_LIST);
          switch (result) {
            case OUT_OF_MEMORY:
            tee_error(od, run_type, overall_line_num,
              E_OUT_OF_MEMORY, PLOT_FAILED);
            return PARSE_INPUT_ERROR;

            case INVALID_LIST_RANGE:
            tee_error(od, run_type, overall_line_num,
              E_LIST_PARSING, "Y variables", "PLOT", PLOT_FAILED);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
          result = set(od, Y_VAR_LIST, OPERATE_BIT, 1, SILENT);
          if (result) {
            tee_error(od, run_type, overall_line_num,
              E_Y_VAR_ALLOWED_RANGE, od->y_vars, PLOT_FAILED);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
        }
        else if ((type & PLS_PLOT) || (type & PCA_PLOT)) {
          if (type & PLS_PLOT) {
            pc = od->pc_num;
          }
          else {
            od->file[TEMP_PCA_LOADINGS]->handle =
              fopen(od->file[TEMP_PCA_LOADINGS]->name, "rb");
            if (!(od->file[TEMP_PCA_LOADINGS]->handle)) {
              tee_error(od, run_type, overall_line_num,
                E_ERROR_IN_READING_TEMP_FILE, "TEMP_FIELD", PLOT_FAILED);
              return PARSE_INPUT_ERROR;
            }
            od->file[TEMP_PCA_SCORES]->handle =
              fopen(od->file[TEMP_PCA_SCORES]->name, "rb");
            if (!(od->file[TEMP_PCA_SCORES]->handle)) {
              tee_error(od, run_type, overall_line_num,
                E_ERROR_IN_READING_TEMP_FILE, "TEMP_FIELD", PLOT_FAILED);
              return PARSE_INPUT_ERROR;
            }
            actual_len = fread(&pc, sizeof(int), 1,
              od->file[TEMP_PCA_LOADINGS]->handle);
            if (actual_len != 1) {
              tee_error(od, run_type, overall_line_num,
                E_ERROR_IN_READING_TEMP_FILE, "TEMP_FIELD", PLOT_FAILED);
              return PARSE_INPUT_ERROR;
            }
          }
          if (pc < 2) {
            tee_error(od, run_type, overall_line_num,
              "At least 2 PCs are needed to generate such a plot.\n%s",
              PLOT_FAILED);
            fail = !(run_type & INTERACTIVE_RUN);
            continue;
          }
          pc_axis[0] = 1;
          pc_axis[1] = 2;
          pc_axis[2] = 0;
          if ((parameter = get_args(od, "pc_x"))) {
            sscanf(parameter, "%d", &pc_axis[0]);
            if ((pc_axis[0] < 1) || (pc_axis[0] > pc)) {
              tee_error(od, run_type, overall_line_num,
                E_PC_ALLOWED_RANGE, pc, PLOT_FAILED);
              fail = !(run_type & INTERACTIVE_RUN);
              continue;
            }
          }
          if ((parameter = get_args(od, "pc_y"))) {
            sscanf(parameter, "%d", &pc_axis[1]);
            if ((pc_axis[1] < 1) || (pc_axis[1] > pc)) {
              tee_error(od, run_type, overall_line_num,
                E_PC_ALLOWED_RANGE, pc, PLOT_FAILED);
              fail = !(run_type & INTERACTIVE_RUN);
              continue;
            }
          }
          if ((parameter = get_args(od, "pc_z"))) {
            sscanf(parameter, "%d", &pc_axis[2]);
            if ((pc_axis[2] < 1) || (pc_axis[2] > pc)) {
              tee_error(od, run_type, overall_line_num,
                E_PC_ALLOWED_RANGE, pc, PLOT_FAILED);
              fail = !(run_type & INTERACTIVE_RUN);
              continue;
            }
          }
        }
        ++command;
        tee_printf(od, M_TOOL_INVOKE, nesting, command, "PLOT", line_orig);
        tee_flush(od);
        result = plot(od, file_basename,
          type, label, pc, pc_axis);
        gettimeofday(&end, NULL);
        elapsed_time(od, &start, &end);
        switch (result) {

          case OUT_OF_MEMORY:
          tee_error(od, run_type, overall_line_num,
            E_OUT_OF_MEMORY, SRD_FAILED);
          return PARSE_INPUT_ERROR;

          case CANNOT_READ_TEMP_FILE:
          tee_error(od, run_type, overall_line_num,
            E_ERROR_IN_READING_TEMP_FILE, "TEMP_FIELD", SRD_FAILED);
          return PARSE_INPUT_ERROR;

          case PREMATURE_PLT_DAT_EOF:
          tee_error(od, run_type, overall_line_num,
            E_FILE_CANNOT_BE_OPENED_FOR_WRITING,
            od->file[PLT_DAT_OUT]->name, SRD_FAILED);
          return PARSE_INPUT_ERROR;

          case PREMATURE_PLT_CMD_EOF:
          tee_error(od, run_type, overall_line_num,
            E_FILE_CANNOT_BE_OPENED_FOR_WRITING,
            od->file[PLT_CMD_OUT]->name, SRD_FAILED);
          return PARSE_INPUT_ERROR;

          default:
          tee_printf(od, M_TOOL_SUCCESS, nesting, command, "PLOT");
          tee_flush(od);
          if (start_gnuplot) {
            #ifndef WIN32
            sprintf(buffer,
              "%s %s %s",
              od->gnuplot.gnuplot_exe, GNUPLOT_ARGS,
              od->file[PLT_CMD_OUT]->name);
            pid1 = fork();
            if (!pid1) {
              setpgid(0, 0);
              /*
              we are the child process; close all inherited
              file descriptors
              */
              memset(&rlp, 0, sizeof(struct rlimit));
              if (!getrlimit(RLIMIT_NOFILE, &rlp)) {
                des = (int)(rlp.rlim_cur);
                while (des >= 0) {
                  close(des);
                  --des;
                }
                /*
                open /dev/null for stdin and duplicate
                the file descriptor also for stdout and stderr,
                then start gnuplot and exit
                */
                des = open("/dev/null", O_RDWR);
                if (des != -1) {
                  dup2(des, STDIN_FILENO);
                  dup2(des, STDOUT_FILENO);
                  dup2(des, STDERR_FILENO);
                  system(buffer);
                }
              }
              exit(0);
            }
            #else
            sprintf(buffer, " %s %s",
              od->file[PLT_CMD_OUT]->name,
              GNUPLOT_ARGS);
            ShellExecute(NULL, TEXT("open"), TEXT(od->gnuplot.gnuplot_exe),
              TEXT(buffer), NULL, SW_SHOW);
            #endif
          }
        }
      }
    }
    else {
      tee_error(od, run_type, overall_line_num,
        "%s: unknown command.\n\n", arg->me[0]);
      if (!(run_type & INTERACTIVE_RUN)) {
        tee_printf(od, PACKAGE_NAME" failed.\n");
        fail = 1;
        continue;
      }
    }
  }
  if (fail) {
    return PARSE_INPUT_RECOVERABLE_ERROR;
  }
  if (run_type & DRY_RUN) {
    for (i = 0; i < od->file_num; ++i) {
      if ((i != MAIN_INPUT) && (i != MAIN_OUTPUT)) {
        if (od->file[i]->handle) {
          fclose(od->file[i]->handle);
          od->file[i]->handle = NULL;
        }
      }
    }
    od->valid = 0;
    od->object_num = 0;
    od->field_num = 0;
    rewind(input_stream);
  }
    
  return 0;
}


#define NUM_SYNONYM_LISTS    4

int parse_synonym_lists(O3Data *od, char *tool_name, char *tool_msg,
  int synonym_list, int *list_type, int default_list,
  int run_type, int overall_line_num)
{
  char *list_names[] = {
    "field_list",
    "object_list",
    "id_list",
    "struct_list"
  };
  char *list_elem[] = {
    "fields",
    "objects",
    "object IDs",
    "structure IDs"
  };
  char *parameter = NULL;
  char comma_hyphen_list[BUF_LEN];
  int i;
  int synonym = 0;
  int result = 0;


  memset(comma_hyphen_list, 0, BUF_LEN);
  if (default_list) {
    *list_type = (1 << default_list);
    strcpy(comma_hyphen_list, "all");
  }
  for (i = 0; i < NUM_SYNONYM_LISTS; ++i) {
    if ((synonym_list & (1 << i))
      && (parameter = get_args(od, list_names[i]))) {
      *list_type = (1 << i);
      ++synonym;
      strcpy(comma_hyphen_list, parameter);
    }
  }
  if (((!default_list) && (!synonym)) || (synonym > 1)) {
    tee_error(od, run_type, overall_line_num,
      E_SUPPLY_OBJECT_ID_STRUCT_FIELD_LIST,
      tool_name, tool_msg);
    return PARSE_INPUT_RECOVERABLE_ERROR;
  }
  if ((*list_type & ((1 << OBJECT_LIST) | (1 << ID_LIST)
    | (1 << STRUCT_LIST))) && (!(od->grid.object_num))) {
    tee_error(od, run_type, overall_line_num,
      E_NO_OBJECTS_PRESENT, tool_msg);
    return PARSE_INPUT_RECOVERABLE_ERROR;
  }
  if ((*list_type & (1 << FIELD_LIST)) && (!(od->grid.field_num))) {
    tee_error(od, run_type, overall_line_num,
      E_NO_FIELDS_PRESENT, tool_msg);
    return PARSE_INPUT_RECOVERABLE_ERROR;
  }
  if ((synonym_list & (1 << FROM_FILE)) && (parameter = get_args(od, "file"))) {
    strcpy(od->file[ASCII_IN]->name, parameter);
    if (!(od->file[ASCII_IN]->handle =
      fopen(od->file[ASCII_IN]->name, "rb"))) {
      tee_error(od, run_type, overall_line_num,
        E_FILE_CANNOT_BE_OPENED_FOR_READING,
        od->file[ASCII_IN]->name, tool_msg);
      return PARSE_INPUT_RECOVERABLE_ERROR;
    }
    *list_type |= (1 << FROM_FILE);
  }
  else {
    result = parse_comma_hyphen_list_to_array(od,
      comma_hyphen_list, intlog2(*list_type));
    switch (result) {
      case OUT_OF_MEMORY:
      tee_error(od, run_type, overall_line_num,
        E_OUT_OF_MEMORY, tool_msg);
      return PARSE_INPUT_ERROR;

      case INVALID_LIST_RANGE:
      for (i = 0; i < NUM_SYNONYM_LISTS; ++i) {
        if (*list_type & (1 << i)) {
          tee_error(od, run_type, overall_line_num, E_LIST_PARSING,
            list_elem[i], tool_name, tool_msg);
          break;
        }
      }
      return PARSE_INPUT_RECOVERABLE_ERROR;
    }
    if (!(run_type & DRY_RUN)) {
      result = set(od, intlog2(*list_type), OPERATE_BIT, 1, SILENT);
      if (result) {
        for (i = 0; i < NUM_SYNONYM_LISTS; ++i) {
          if (*list_type & (1 << i)) {
            tee_error(od, run_type, overall_line_num,
              "The specified %s are out of range.\n%s",
              list_elem[i], tool_msg);
            break;
          }
        }
        return PARSE_INPUT_RECOVERABLE_ERROR;
      }
    }
  }
  
  return 0;
}


int check_lmo_parameters(O3Data *od, char *tool_msg, int *groups, int *runs,
  int run_type, int overall_line_num)
{
  char *parameter;
  int active_struct_num;
  
  
  *groups = 5;
  if ((parameter = get_args(od, "groups"))) {
    sscanf(parameter, "%d", groups);
  }
  if (*groups < 2) {
    tee_error(od, run_type, overall_line_num,
      E_TOO_FEW_MANY_FOR_AVAILABLE_DATA, "few groups",
      tool_msg);
    return PARSE_INPUT_RECOVERABLE_ERROR;
  }
  get_attr_struct_ave(od, 0, ACTIVE_BIT, &active_struct_num, NULL);
  if (*groups > active_struct_num) {
    tee_error(od, run_type, overall_line_num,
      E_TOO_FEW_MANY_FOR_AVAILABLE_DATA, "many groups",
      tool_msg);
    return PARSE_INPUT_RECOVERABLE_ERROR;
  }
  *runs = 20;
  if ((parameter = get_args(od, "runs"))) {
    sscanf(parameter, "%d", runs);
  }
  if (*runs < 1) {
    tee_error(od, run_type, overall_line_num,
      "At least one run must be performed "
      "for LMO cross-validation.\n%s",
      CV_FAILED);
    return PARSE_INPUT_RECOVERABLE_ERROR;
  }
  
  return 0;
}
