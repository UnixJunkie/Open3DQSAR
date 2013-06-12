/*

update_jmol.c

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
#include <include/prog_exe_info.h>
#ifndef WIN32
#ifndef INVALID_SOCKET
#define INVALID_SOCKET  -1
#endif
#ifndef SOCKET_ERROR
#define SOCKET_ERROR  -1
#endif
#endif


int send_jmol_command(O3Data *od, char *command)
{
  const char jmol_cmd_start[] =
    "{\"type\":\"command\", \"command\":\"";
  const char jmol_cmd_end[] =
    "\"}" NEWLINE;
  char *buffer;
  int error;
  
  
  if (!(buffer = malloc(strlen(jmol_cmd_start)
    + strlen(command) + strlen(jmol_cmd_end) + 1))) {
    return OUT_OF_MEMORY;
  }
  sprintf(buffer, "%s%s%s", jmol_cmd_start,
    command, jmol_cmd_end);
  error = send(od->jmol.sock, buffer, strlen(buffer), 0);
  free(buffer);
  
  return ((error == -1) ? CANNOT_SEND_JMOL_COMMAND : 0);
}


int update_jmol(O3Data *od)
{
  const char jmol_init[] =
    "background BLACK\n"
    "zap\n"
    "set antialiasDisplay ON\n"
    "set perspectiveDepth OFF\n"
    "set showMultipleBonds OFF\n"
    "set bondmode AND\n"
    "set frank OFF\n"
    "set ambient 50\n"
    "set specpower 40\n"
    "set appendNew true\n";
  const char jmol_after_load[] =
    "select ALL\n"
    "hbonds OFF\n"
    "spin OFF\n"
    "trace OFF\n"
    "ribbons OFF\n"
    "cartoons OFF\n"
    "label OFF\n"
    "monitor OFF\n"
    "wireframe 0.06\n"
    "spacefill OFF\n"
    "color HYDROGEN [204, 255, 255]\n"
    "color CARBON [204, 204, 204]\n"
    "color NITROGEN [0, 0, 255]\n"
    "color OXYGEN [255, 0, 0]\n"
    "color PHOSPHORUS [191, 0, 191]\n"
    "color SULPHUR [230, 198, 64]\n"
    "color FLUORINE [191, 255, 64]\n"
    "color CHLORINE [0, 255, 0]\n"
    "color BROMINE [142, 57, 28]\n"
    "color IODINE [140, 64, 153]\n"
    "select HYDROGEN AND WITHIN(1.15, CARBON)\n"
    "hide SELECTED\n"
    "frame ALL\n"
    "reset\n"
    "set refreshing TRUE\n";
  char buffer[BUF_LEN];
  int first;
  int i;
  int j;
  int n;
  int found;
  int error = 0;
  int attempts;
  FileDescriptor temp_jmol_fd;
  ProgExeInfo prog_exe_info;
  #ifdef WIN32
  WORD wVersionRequested;
  WSADATA wsaData;
  #endif


  memset(buffer, 0, BUF_LEN);
  if (!(od->jmol.use_jmol)) {
    return 0;
  }
  if (!(od->jmol.jmol_pid)) {
    while ((!error) && (od->jmol.port < JMOL_LAST_PORT)) {
      od->jmol.sock = socket(AF_INET, SOCK_STREAM, 0);
      if (od->jmol.sock == INVALID_SOCKET) {
        return 0;
      }
      memset(&(od->jmol.sockaddr), 0, sizeof(od->jmol.sockaddr));
      od->jmol.sockaddr.sin_family = AF_INET;
      od->jmol.sockaddr.sin_addr.s_addr = inet_addr(LOCALHOST_IP);
      od->jmol.sockaddr.sin_port = htons(od->jmol.port);
      error = connect(od->jmol.sock, (struct sockaddr *)&(od->jmol.sockaddr),
        sizeof(od->jmol.sockaddr));
      if (!error) {
        ++(od->jmol.port);
        #ifdef WIN32
        closesocket(od->jmol.sock);
        #else
        close(od->jmol.sock);
        #endif
      }
    }
    if (!error) {
      return 0;
    }
    memset(&prog_exe_info, 0, sizeof(ProgExeInfo));
    strcpy(buffer, od->jmol.jmol_exe);
    get_dirname(buffer);
    prog_exe_info.proc_env = fill_env(od, minimal_env, buffer, 0);
    if (!(prog_exe_info.proc_env)) {
      return OUT_OF_MEMORY;
    }
    prog_exe_info.exedir = buffer;
    sprintf(prog_exe_info.command_line, "%s -J \"sync -%d\"",
      od->jmol.jmol_exe, od->jmol.port);
    ext_program_exe(&prog_exe_info, &error);
    od->jmol.jmol_pid = 1;
    error = -1;
    attempts = 0;
    while (error && (attempts < MAX_ATTEMPTS_JMOL)) {
      error = connect(od->jmol.sock, (struct sockaddr *)&(od->jmol.sockaddr),
        sizeof(od->jmol.sockaddr));
      if (error == -1) {
        #ifndef WIN32
        usleep(1000 * MS_SLEEP_BEFORE_RETRY);
        #else
        Sleep(MS_SLEEP_BEFORE_RETRY);
        #endif
        ++attempts;
      }
    }
  }
  if (error) {
    return 0;
  }
  if (open_temp_file(od, &temp_jmol_fd, "jmol")) {
    return CANNOT_WRITE_TEMP_FILE;
  }
  if (!(od->pel.jmol_old_object_id)) {
    fputs(jmol_init, temp_jmol_fd.handle);
  }
  fputs("select ALL; hide SELECTED; set refreshing FALSE\n", temp_jmol_fd.handle);
  if (od->jmol.grid_box) {
    fputs("delete $grid_box\n"
      "select DUMMY; zap SELECTED\n",
      temp_jmol_fd.handle);
  }
  if (od->grid.nodes[0]) {
    fprintf(temp_jmol_fd.handle,
      "boundbox CORNERS {%.4f, %.4f, %.4f} {%.4f, %.4f, %.4f} OFF; "
      "draw ID grid_box BOUNDBOX MESH NOFILL WIDTH 0.05 COLOR {128, 255, 128}\n"
      "data \"APPEND dummybox\"|2|dummybox|"
      "Du %.4f %.4f %.4f|Du %.4f %.4f %.4f|END \"APPEND dummybox\"\n"
      "select DUMMY; hide SELECTED\n",
      od->grid.start_coord[0], od->grid.start_coord[1], od->grid.start_coord[2],
      od->grid.end_coord[0], od->grid.end_coord[1], od->grid.end_coord[2],
      od->grid.start_coord[0], od->grid.start_coord[1], od->grid.start_coord[2],
      od->grid.end_coord[0], od->grid.end_coord[1], od->grid.end_coord[2]);
    od->jmol.grid_box = 1;
  }
  od->pel.jmol_object_id = int_perm_resize(od->pel.jmol_object_id,
    (od->field_num ? (od->active_object_num + od->ext_pred_object_num)
    : od->grid.object_num));
  if (!(od->pel.jmol_object_id)) {
    if (temp_jmol_fd.handle) {
      fclose(temp_jmol_fd.handle);
      temp_jmol_fd.handle = NULL;
    }
    return OUT_OF_MEMORY;
  }
  for (i = 0, j = 0; i < od->grid.object_num; ++i) {
    if (get_object_attr(od, i, ACTIVE_BIT)
      || get_object_attr(od, i, PREDICT_BIT)) {
      od->pel.jmol_object_id->pe[j] = od->al.mol_info[i]->object_id;
      ++j;
    }
  }
  n = -1;
  if (od->pel.jmol_old_object_id) {
    for (n = 0, i = 0; i < od->pel.jmol_old_object_id->size; ++i) {
      for (j = 0, found = 0; j < od->pel.jmol_object_id->size; ++j) {
        if ((found = (od->pel.jmol_old_object_id->pe[i] ==
          od->pel.jmol_object_id->pe[j]))) {
          break;
        }
      }
      if (!found) {
        ++n;
        fprintf(temp_jmol_fd.handle, "zap %d.1\n",
          od->pel.jmol_old_object_id->pe[i]);
      }
    }
  }
  first = ((n == -1) || (od->pel.jmol_old_object_id
    && (n == od->pel.jmol_old_object_id->size)));
  for (i = 0; i < od->pel.jmol_object_id->size; ++i) {
    found = 0;
    if (od->pel.jmol_old_object_id) {
      for (j = 0; j < od->pel.jmol_old_object_id->size; ++j) {
         if ((found = (od->pel.jmol_object_id->pe[i] ==
          od->pel.jmol_old_object_id->pe[j]))) {
          break;
        }
      }
    }
    if (!found) {
      fprintf(temp_jmol_fd.handle,
        "load %s \"%s%c%04d.mol2\"\n",
        (first ? "" : "APPEND"),
        od->field.mol_dir, SEPARATOR,
        od->al.mol_info[i]->object_id);
      first = 0;
    }
  }
  od->pel.jmol_old_object_id = int_perm_resize
    (od->pel.jmol_old_object_id, od->pel.jmol_object_id->size);
  if (!(od->pel.jmol_old_object_id)) {
    if (temp_jmol_fd.handle) {
      fclose(temp_jmol_fd.handle);
      temp_jmol_fd.handle = NULL;
    }
    return OUT_OF_MEMORY;
  }
  if (od->pel.jmol_object_id->size) {
    memcpy(od->pel.jmol_old_object_id->pe, od->pel.jmol_object_id->pe,
      od->pel.jmol_object_id->size * sizeof(int));
  }
  fputs(jmol_after_load, temp_jmol_fd.handle);
  for (i = 0; i < od->grid.object_num; ++i) {
    if ((od->field_num && (get_object_attr(od, i, ACTIVE_BIT)
      || get_object_attr(od, i, PREDICT_BIT))) || (!(od->field_num))) {
      fprintf(temp_jmol_fd.handle,
        "select CARBON AND %d.1; color %s\n",
        od->al.mol_info[i]->object_id,
        (od->field_num ? (get_object_attr(od, i, ACTIVE_BIT)
        ? JMOL_TRAINING_SET_CARBON : JMOL_TEST_SET_CARBON)
        : JMOL_TRAINING_SET_CARBON));
    }
  }
  if (temp_jmol_fd.handle) {
    fclose(temp_jmol_fd.handle);
    temp_jmol_fd.handle = NULL;
  }
  sprintf(buffer, "script \\\"%s\\\"", temp_jmol_fd.name);
  send_jmol_command(od, buffer);

  return 0;
}
