/*

grid_box.c

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


int create_box(O3Data *od, GridInfo *temp_grid, double outgap, int from_file)
{
  char buffer[BUF_LEN];
  int i;
  double start_coord[3];
  double end_coord[3];


  memset(buffer, 0, BUF_LEN);
  if (from_file) {
    if (!fgrep(od->file[ASCII_IN]->handle, buffer, "Points")) {
      return PREMATURE_DAT_EOF;
    }
    sscanf(buffer, "%*s %d", &(od->grid.x_vars));
    if (!fgrep(od->file[ASCII_IN]->handle, buffer, "Lower Corner")) {
      return PREMATURE_DAT_EOF;
    }
    sscanf(buffer, "%*s %*s %lf %lf %lf", &start_coord[0],
      &start_coord[1], &start_coord[2]);
    if (!fgrep(od->file[ASCII_IN]->handle, buffer, "Number steps")) {
      return PREMATURE_DAT_EOF;
    }
    sscanf(buffer, "%*s %*s %d %d %d", &(od->grid.nodes[0]),
      &(od->grid.nodes[1]), &(od->grid.nodes[2]));
    if (od->grid.x_vars != od->grid.nodes[0] * od->grid.nodes[1] * od->grid.nodes[2]) {
      return PREMATURE_DAT_EOF;
    }
    if (!fgrep(od->file[ASCII_IN]->handle, buffer, "Step Size")) {
      return PREMATURE_DAT_EOF;
    }
    sscanf(buffer, "%*s %*s %f %f %f", &(od->grid.step[0]),
      &(od->grid.step[1]), &(od->grid.step[2]));
    for (i = 0; i < 3; ++i) {
      end_coord[i] = start_coord[i] + (double)(od->grid.nodes[i] - 1)
        * (double)(od->grid.step[i]);
      od->grid.start_coord[i] = (float)start_coord[i];
      od->grid.end_coord[i] = (float)end_coord[i];
    }
  }
  else {
    od->grid.x_vars = 1;
    for (i = 0; i < 3; ++i) {
      od->grid.nodes[i] = 1;
      od->grid.step[i] = temp_grid->step[i];
      if (outgap < 0.0) {
        start_coord[i] = (double)(temp_grid->start_coord[i]);
        end_coord[i] = (double)(temp_grid->start_coord[i]);
      }
      else {
        temp_grid->start_coord[i] = (float)(od->min_coord[i]);
        temp_grid->end_coord[i] = (float)(od->max_coord[i] + outgap);
        start_coord[i] = (double)((int)safe_rint
          ((double)(temp_grid->start_coord[i]) - outgap));
        end_coord[i] = start_coord[i];
      }
      if (!(temp_grid->nodes[i])) {
        while (((double)(temp_grid->end_coord[i]) - end_coord[i]) > GRID_TOLERANCE) {
          end_coord[i] += (double)(temp_grid->step[i]);
          ++(od->grid.nodes[i]);
        }
      }
      else {
        while (od->grid.nodes[i] < temp_grid->nodes[i]) {
          end_coord[i] += (double)(temp_grid->step[i]);
          ++(od->grid.nodes[i]);
        }
      }
      od->grid.x_vars *= od->grid.nodes[i];
      od->grid.start_coord[i] = (float)start_coord[i];
      od->grid.end_coord[i] = (float)end_coord[i];
    }
  }
  od->x_vars = od->grid.x_vars;
  od->object_num = od->grid.object_num;
  od->active_object_num = od->object_num;
  od->object_pagesize = (od->x_vars * sizeof(float) / od->mmap_pagesize +
    (od->x_vars * sizeof(float) % od->mmap_pagesize ? 1 : 0))
    * od->mmap_pagesize;
  
  return 0;
}


void remove_box(O3Data *od)
{
  int object_num;
  
  
  object_num = od->grid.object_num;
  memset(&(od->grid), 0, sizeof(GridInfo));
  od->grid.object_num = object_num;
  od->grid.struct_num = object_num;
  od->x_vars = 0;
}


void print_grid_coordinates(O3Data *od, GridInfo *grid_info)
{
  int i;
  
  
  for (i = 0; i < 3; ++i) {
    tee_printf(od, "%c start, %c end coordinates:    %.4f, %.4f\n",
      'X' + i, 'X' + i, grid_info->start_coord[i], grid_info->end_coord[i]);
  }
  for (i = 0; i < 3; ++i) {
    tee_printf(od, "%c nodes:                       %d\n",
      'X' + i, grid_info->nodes[i]);
  }
  for (i = 0; i < 3; ++i) {
    tee_printf(od, "%c step:                        %.4f\n",
      'X' + i, grid_info->step[i]);
  }
  tee_printf(od, "Grid points:                   %d\n\n", grid_info->x_vars);
}


void print_grid_comparison(O3Data *od)
{
  tee_printf(od,
    "\n"
    "Existing grid:\n"
    "--------------\n");
  print_grid_coordinates(od, &(od->grid));
  tee_printf(od,
    "Imported grid:\n"
    "--------------\n");
  print_grid_coordinates(od, &(od->newgrid));
}


int match_grids(O3Data *od)
{
  int i;
  
  
  for (i = 0; i < 3; ++i) {
    if (((int)safe_rint((double)(od->newgrid.start_coord[i]) * 1.0e03) != 
      (int)safe_rint((double)(od->grid.start_coord[i]) * 1.0e03))
      || ((int)safe_rint((double)(od->newgrid.end_coord[i]) * 1.0e03) != 
      (int)safe_rint((double)(od->grid.end_coord[i]) * 1.0e03))
      || ((int)safe_rint((double)(od->newgrid.step[i]) * 1.0e03) != 
      (int)safe_rint((double)(od->grid.step[i]) * 1.0e03))
      || (od->newgrid.nodes[i] != od->grid.nodes[i])) {
      return 0;
    }
  }
  
  return 1;
}
