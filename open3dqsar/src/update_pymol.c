/*

update_pymol.c

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
#include <include/proc_env.h>
#include <include/prog_exe_info.h>


void vertex_xyz(O3Data *od, FILE *handle, int x, int y, int z)
{
  fprintf(handle,
    "VERTEX, %.4f, %.4f, %.4f, ",
    (x ? od->grid.end_coord[0] : od->grid.start_coord[0]),
    (y ? od->grid.end_coord[1] : od->grid.start_coord[1]),
    (z ? od->grid.end_coord[2] : od->grid.start_coord[2]));
}


int update_pymol(O3Data *od)
{
  char pymol_init[] =
    "from pymol.cgo import *\n"
    "from pymol import cmd\n"
    "reinitialize\n"
    "viewport 800, 600\n"
    "set stick_radius = 0.06\n"
    "set_color hydrogen = [0.8, 1.0, 1.0]\n"
    "set_color carbon = [0.8, 0.8, 0.8]\n"
    "set_color nitrogen = [0.0, 0.0, 1.0]\n"
    "set_color oxygen = [1.0, 0.0, 0.0]\n"
    "set_color phosphorus = [0.75, 0.0, 0.75]\n"
    "set_color sulphur = [0.9, 0.775, 0.25]\n"
    "set_color fluorine = [0.75, 1.00, 0.25]\n"
    "set_color chlorine = [0.0, 1.0, 0.0]\n"
    "set_color bromine = [0.555, 0.222, 0.111]\n"
    "set_color iodine = [0.55, 0.25, 0.6]\n";
  char buffer[BUF_LEN];
  int i;
  int j;
  int x;
  int y;
  int z;
  int found;
  int error;
  FileDescriptor temp_pymol_fd;
  ProgExeInfo prog_exe_info;
  #ifdef WIN32
  DWORD n_chr = 0;
  #endif


  memset(buffer, 0, BUF_LEN);
  memset(&temp_pymol_fd, 0, sizeof(FileDescriptor));
  if (!(od->pymol.use_pymol)) {
    return 0;
  }
  if (!(od->pymol.pymol_pid)) {
    memset(&prog_exe_info, 0, sizeof(ProgExeInfo));
    strcpy(buffer, od->pymol.pymol_exe);
    get_dirname(buffer);
    prog_exe_info.proc_env = fill_env(od, minimal_env, buffer, 0);
    if (!(prog_exe_info.proc_env)) {
      return OUT_OF_MEMORY;
    }
    od->pymol.proc_env = prog_exe_info.proc_env;
    prog_exe_info.need_stdin = NEED_STDIN_LEAVE_READ_PIPE_OPEN;
    prog_exe_info.stdout_fd = 0;
    prog_exe_info.stderr_fd = 0;
    prog_exe_info.exedir = buffer;
    prog_exe_info.sep_proc_grp = 1;
    strcpy(prog_exe_info.command_line, od->pymol.pymol_exe);
    ext_program_exe(&prog_exe_info, &error);
    od->pymol.pymol_pid = 1;
    #ifndef WIN32
    if (!error) {
      od->pymol.pymol_pipe = prog_exe_info.pipe_des[1];
      if (!(od->pymol.pymol_handle = fdopen(prog_exe_info.pipe_des[1], "w"))) {
        return FL_CANNOT_CREATE_CHANNELS;
      }
    }
    #else
    if (!error) {
      CloseHandle(prog_exe_info.proc_info.hProcess);
      CloseHandle(prog_exe_info.proc_info.hThread);
      od->pymol.pymol_handle = prog_exe_info.stdin_wr;
    }
    #endif
  }
  if (open_temp_file(od, &temp_pymol_fd, "pymol")) {
    return CANNOT_WRITE_TEMP_FILE;
  }
  if (!(od->pel.pymol_old_object_id)) {
    fputs(pymol_init, temp_pymol_fd.handle);
  }
  if (getcwd(buffer, BUF_LEN)) {
    fprintf(temp_pymol_fd.handle, "cd %s\n", buffer);
  }
  fputs("disable all\n", temp_pymol_fd.handle);
  if (od->pymol.grid_box) {
    fputs("delete grid_box\n", temp_pymol_fd.handle);
  }
  if (od->grid.nodes[0]) {
    fputs("cmd.load_cgo([LINEWIDTH, float(1.0), BEGIN, LINES, "
      "COLOR, float(0.5), float(1.0), float(0.5), ",
      temp_pymol_fd.handle);
    for (x = 0; x <= 1; ++x) {
      for (y = 0; y <= 1; ++y) {
        for (z = 0; z <= 1; ++z) {
          vertex_xyz(od, temp_pymol_fd.handle, x, y, z);
        }
      }
    }
    for (z = 0; z <= 1; ++z) {
      for (y = 0; y <= 1; ++y) {
        for (x = 0; x <= 1; ++x) {
          vertex_xyz(od, temp_pymol_fd.handle, x, y, z);
        }
      }
    }
    for (z = 0; z <= 1; ++z) {
      for (x = 0; x <= 1; ++x) {
        for (y = 0; y <= 1; ++y) {
          vertex_xyz(od, temp_pymol_fd.handle, x, y, z);
        }
      }
    }
    fputs("END], \"grid_box\")\n", temp_pymol_fd.handle);
    od->pymol.grid_box = 1;
  }
  od->pel.pymol_object_id = int_perm_resize(od->pel.pymol_object_id,
    (od->field_num ? (od->active_object_num + od->ext_pred_object_num)
    : od->grid.object_num));
  if (!(od->pel.pymol_object_id)) {
    if (temp_pymol_fd.handle) {
      fclose(temp_pymol_fd.handle);
      temp_pymol_fd.handle = NULL;
    }
    return OUT_OF_MEMORY;
  }
  for (i = 0, j = 0; i < od->grid.object_num; ++i) {
    if (get_object_attr(od, i, ACTIVE_BIT)
      || get_object_attr(od, i, PREDICT_BIT)) {
      od->pel.pymol_object_id->pe[j] = od->al.mol_info[i]->object_id;
      ++j;
    }
  }
  if (od->pel.pymol_old_object_id) {
    for (i = 0; i < od->pel.pymol_old_object_id->size; ++i) {
      for (j = 0, found = 0; j < od->pel.pymol_object_id->size; ++j) {
        if ((found = (od->pel.pymol_old_object_id->pe[i] ==
          od->pel.pymol_object_id->pe[j]))) {
          break;
        }
      }
      if (!found) {
        fprintf(temp_pymol_fd.handle, "delete %04d\n",
          od->pel.pymol_old_object_id->pe[i]);
      }
    }
  }
  for (i = 0; i < od->pel.pymol_object_id->size; ++i) {
    found = 0;
    if (od->pel.pymol_old_object_id) {
      for (j = 0; j < od->pel.pymol_old_object_id->size; ++j) {
         if ((found = (od->pel.pymol_object_id->pe[i] ==
          od->pel.pymol_old_object_id->pe[j]))) {
          break;
        }
      }
    }
    if (!found) {
      fprintf(temp_pymol_fd.handle,
        "load %s%c%04d.mol2, %04d, 0, mol2\n",
         od->field.mol_dir, SEPARATOR,
        od->al.mol_info[i]->object_id,
        od->al.mol_info[i]->object_id);
    }
  }
  od->pel.pymol_old_object_id = int_perm_resize
    (od->pel.pymol_old_object_id, od->pel.pymol_object_id->size);
  if (!(od->pel.pymol_old_object_id)) {
    if (temp_pymol_fd.handle) {
      fclose(temp_pymol_fd.handle);
      temp_pymol_fd.handle = NULL;
    }
    return OUT_OF_MEMORY;
  }
  if (od->pel.pymol_object_id->size) {
    memcpy(od->pel.pymol_old_object_id->pe, od->pel.pymol_object_id->pe,
      od->pel.pymol_object_id->size * sizeof(int));
  }
  for (i = 0; i < od->grid.object_num; ++i) {
    if ((od->field_num && (get_object_attr(od, i, ACTIVE_BIT)
      || get_object_attr(od, i, PREDICT_BIT))) || (!(od->field_num))) {
      fprintf(temp_pymol_fd.handle,
        "color %s, %04d and symbol C\n",
        (od->field_num ? (get_object_attr(od, i, ACTIVE_BIT)
        ? PYMOL_TRAINING_SET_CARBON : PYMOL_TEST_SET_CARBON)
        : PYMOL_TRAINING_SET_CARBON),
        od->al.mol_info[i]->object_id);
    }
  }
  fputs("hide lines\n"
    "show sticks, (not "
    "(hydro and neighbor symbol C))\n"
    "reset\n"
    "enable all\n",
    temp_pymol_fd.handle);
  if (temp_pymol_fd.handle) {
    fclose(temp_pymol_fd.handle);
    temp_pymol_fd.handle = NULL;
  }
  sprintf(buffer, "@%s\n", temp_pymol_fd.name);
  FWRITE_WRAP(od->pymol.pymol_handle, buffer, &n_chr);
  FFLUSH_WRAP(od->pymol.pymol_handle);

  return 0;
}
