/*

save_dat.c

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


int save_dat(O3Data *od, int file_id)
{
  char buffer[LARGE_BUF_LEN];
  char print_buf[BUF_LEN];
  char o3g_header[TITLE_LEN + 1];
  int binary_int;
  int actual_len;
  int mol_len;
  int to_be_read;
  int not_deleted_field_num;
  int not_deleted_object_num;
  int field_num;
  int object_num;
  int i;
  int j;
  int k;
  uint64_t valid;
  float value;
  FileDescriptor mol_fd;
  fzPtr *dat_out;
  

  memset(&mol_fd, 0, sizeof(FileDescriptor));
  dat_out = (fzPtr *)(od->file[file_id]->handle);
  /*
  as the first 4-byte word of the .dat file
  write (int)1 as an endianness indicator
  */
  binary_int = 1;
  fzwrite(&binary_int, sizeof(int), 1, dat_out);
  /*
  then, write a short header about software version
  */
  memset(buffer, 0, LARGE_BUF_LEN);
  memset(o3g_header, ' ', TITLE_LEN);
  o3g_header[TITLE_LEN] = '\0';
  sprintf(o3g_header, PACKAGE_NAME" v %s", VERSION);
  for (i = 0; i < TITLE_LEN; ++i) {
    if (!o3g_header[i]) {
      o3g_header[i] = ' ';
      break;
    }
  }
  fzwrite(o3g_header, 1, TITLE_LEN, dat_out);
  not_deleted_field_num = 0;
  for (i = 0; i < od->field_num; ++i) {
    if (!get_field_attr(od, i, DELETE_BIT)) {
      ++not_deleted_field_num;
    }
  }
  not_deleted_object_num = 0;
  for (i = 0; i < od->grid.object_num; ++i) {
    if (!get_object_attr(od, i, DELETE_BIT)) {
      ++not_deleted_object_num;
    }
  }
  /*
  then, write 4 int:
  - number of fields;
  - number of objects;
  - number of x_vars.
  - number of y_vars;
  */
  fzwrite(&not_deleted_field_num, sizeof(int), 1, dat_out);
  fzwrite(&not_deleted_object_num, sizeof(int), 1, dat_out);
  fzwrite(&(od->x_vars), sizeof(int), 1, dat_out);
  actual_len = fzwrite(&(od->y_vars), sizeof(int), 1, dat_out);
  if (actual_len != 1) {
    return PREMATURE_DAT_EOF;
  }

  /*
  then, write 6 float:
  - x, y, z start coordinates of 3D-grid;
  - x, y, z end coordinates of 3D-grid;
  */
  for (i = 0; i < 3; ++i) {
    actual_len = fzwrite(&(od->grid.start_coord[i]), sizeof(float), 1, dat_out);
    if (actual_len != 1) {
      return PREMATURE_DAT_EOF;
    }
  }
  for (i = 0; i < 3; ++i) {
    actual_len = fzwrite(&(od->grid.end_coord[i]), sizeof(float), 1, dat_out);
    if (actual_len != 1) {
      return PREMATURE_DAT_EOF;
    }
  }
  
  /*
  then, write 3 int:
  - number of nodes in x, y, z directions;
  */
  for (i = 0; i < 3; ++i) {
    actual_len = fzwrite(&(od->grid.nodes[i]), sizeof(float), 1, dat_out);
    if (actual_len != 1) {
      return PREMATURE_DAT_EOF;
    }
  }

  /*
  then, write 3 float:
  - step size in x, y, z directions;
  */
  for (i = 0; i < 3; ++i) {
    actual_len = fzwrite(&(od->grid.step[i]), sizeof(float), 1, dat_out);
    if (actual_len != 1) {
      return PREMATURE_DAT_EOF;
    }
  }
  
  /*
  then, write an unsigned long:
  - od->valid, after zeroing the PLS, CV, and PREDICT bits
  */
  valid = (uint64_t)(od->valid & (~(PLS_BIT | CV_BIT | PREDICT_BIT)));
  fzwrite(&valid, sizeof(uint64_t), 1, dat_out);
  
  /*
  then, write a number of XData structures equal to active_field_num
  */
  if (od->active_field_num > 0) {
    fzwrite(od->mel.x_data,
      sizeof(XData), od->active_field_num, dat_out);
  }
  field_num = 0;
  for (i = 0; i < od->field_num; ++i) {
    if (get_field_attr(od, i, DELETE_BIT)) {
      continue;
    }
    object_num = 0;
    for (j = 0; j < od->object_num; ++j) {
      if (get_object_attr(od, j, DELETE_BIT)) {
        continue;
      }
      /*
      write out the header to the .dat file
      "HEADER                          "
      "OBJECT_NAME                     "
      */
      memset(print_buf, ' ', 64);
      snprintf(print_buf, 64, "%-32s%-32s", DAT_HEADER,
        od->al.mol_info[j]->object_name);
      fzwrite(print_buf, 1, 64, dat_out);
      /*
      write 3 int:
      - field number
      - object number
      - object ID
      */
      fzwrite(&field_num, sizeof(int), 1, dat_out);
      fzwrite(&object_num, sizeof(int), 1, dat_out);
      actual_len = fzwrite(&(od->al.mol_info[j]->object_id), sizeof(int), 1, dat_out);
      if (actual_len != 1) {
        return PREMATURE_DAT_EOF;
      }
      /*
      write MOL information
      */
      mol_len = 0;
      if (!field_num) {
        sprintf(mol_fd.name, "%s%c%04d.mol", od->field.mol_dir,
          SEPARATOR, od->al.mol_info[j]->object_id);
        mol_fd.handle = fopen(mol_fd.name, "rb");
        if (mol_fd.handle) {
          actual_len = LARGE_BUF_LEN;
          while (actual_len == LARGE_BUF_LEN) {
            actual_len = fread(buffer, 1, LARGE_BUF_LEN, mol_fd.handle);
            mol_len += actual_len;
          }
          rewind(mol_fd.handle);
        }
      }
      actual_len = fzwrite(&mol_len, sizeof(int), 1, dat_out);
      if (actual_len != 1) {
        return PREMATURE_DAT_EOF;
      }
      if (!field_num) {
        if (mol_fd.handle) {
          actual_len = LARGE_BUF_LEN;
          while (actual_len == LARGE_BUF_LEN) {
            mol_len = fread(buffer, 1, LARGE_BUF_LEN, mol_fd.handle);
            if (mol_len > 0) {
              actual_len = fzwrite(buffer, 1, mol_len, dat_out);
              if (actual_len != mol_len) {
                return PREMATURE_DAT_EOF;
              }
            }
            else {
              break;
            }
          }
          fclose(mol_fd.handle);
          mol_fd.handle = NULL;
        }
      }
      /*
      write x_vars
      */
      if (od->save_ram) {
        if (od->mmap_field_num != -1) {
          for (k = 0; k < od->object_num; ++k) {
            #ifndef WIN32
            munmap(od->mel.x_var_array[od->mmap_field_num][k], od->object_pagesize);
            #else
            UnmapViewOfFile(od->mel.x_var_array[od->mmap_field_num][k]);
            CloseHandle(od->al.mol_info[k]->hMapHandle);
            #endif
          }
          od->mmap_field_num = -1;
        }
        if (fseek(od->file[TEMP_FIELD_DATA + i]->handle,
          od->object_pagesize * j, SEEK_SET)) {
          return PREMATURE_DAT_EOF;
        }
        to_be_read = LARGE_BUF_LEN;
        for (k = (sizeof(float) * od->x_vars / LARGE_BUF_LEN); k >= 0; --k) {
          if (!k) {
            to_be_read = sizeof(float) * od->x_vars % LARGE_BUF_LEN;
          }
          actual_len = fread(buffer, 1, to_be_read,
            od->file[TEMP_FIELD_DATA + i]->handle);
          if (actual_len != to_be_read) {
            return PREMATURE_DAT_EOF;
          }
          actual_len = fzwrite(buffer, 1, to_be_read, dat_out);
          if (actual_len != to_be_read) {
            return PREMATURE_DAT_EOF;
          }
        }
      }
      else {
        actual_len = fzwrite(od->mel.x_var_array[i][j],
          sizeof(float), od->x_vars, dat_out);
        if (actual_len != od->x_vars) {
          return PREMATURE_DAT_EOF;
        }
      }
      ++object_num;
    }
    ++field_num;
  }
  if (!(od->field_num)) {
    field_num = -1;
    object_num = 0;
    for (j = 0; j < od->grid.object_num; ++j) {
      if (get_object_attr(od, j, DELETE_BIT)) {
        continue;
      }
      /*
      write out the header to the .dat file
      "HEADER                          "
      "OBJECT_NAME                     "
      */
      memset(print_buf, ' ', 64);
      snprintf(print_buf, 64, "%-32s%-32s", DAT_HEADER,
        od->al.mol_info[j]->object_name);
      fzwrite(print_buf, 1, 64, dat_out);
      /*
      write 3 int:
      - field number
      - object number
      - object ID
      */
      fzwrite(&field_num, sizeof(int), 1, dat_out);
      fzwrite(&object_num, sizeof(int), 1, dat_out);
      actual_len = fzwrite(&(od->al.mol_info[j]->object_id), sizeof(int), 1, dat_out);
      if (actual_len != 1) {
        return PREMATURE_DAT_EOF;
      }
      /*
      write MOL information
      */
      mol_len = 0;
      sprintf(mol_fd.name, "%s%c%04d.mol", od->field.mol_dir,
        SEPARATOR, od->al.mol_info[j]->object_id);
      mol_fd.handle = fopen(mol_fd.name, "rb");
      if (mol_fd.handle) {
        actual_len = LARGE_BUF_LEN;
        while (actual_len == LARGE_BUF_LEN) {
          actual_len = fread(buffer, 1, LARGE_BUF_LEN, mol_fd.handle);
          mol_len += actual_len;
        }
        rewind(mol_fd.handle);
      }
      actual_len = fzwrite(&mol_len, sizeof(int), 1, dat_out);
      if (actual_len != 1) {
        return PREMATURE_DAT_EOF;
      }
      if (mol_fd.handle) {
        actual_len = LARGE_BUF_LEN;
        while (actual_len == LARGE_BUF_LEN) {
          mol_len = fread(buffer, 1, LARGE_BUF_LEN, mol_fd.handle);
          if (mol_len > 0) {
            actual_len = fzwrite(buffer, 1, mol_len, dat_out);
            if (actual_len != mol_len) {
              return PREMATURE_DAT_EOF;
            }
          }
          else {
            break;
          }
        }
        fclose(mol_fd.handle);
        mol_fd.handle = NULL;
      }
      ++object_num;
    }
  }
  /*
  Now save attributes:
  - field_attr
  - object_attr
  - object_weight
  - x_var_attr
  */
  for (i = 0; i < od->field_num; ++i) {
    if (!get_field_attr(od, i, DELETE_BIT)) {
      fzwrite(&(od->mel.field_attr[i]),
        sizeof(uint16_t), 1, dat_out);
    }
  }
  for (i = 0; i < od->grid.object_num; ++i) {
    if (!get_object_attr(od, i, DELETE_BIT)) {
      fzwrite(&(od->mel.object_attr[i]),
        sizeof(uint16_t), 1, dat_out);
    }
  }
  for (i = 0; i < od->grid.object_num; ++i) {
    if (!get_object_attr(od, i, DELETE_BIT)) {
      fzwrite(&(od->mel.object_weight[i]),
        sizeof(double), 1, dat_out);
    }
  }
  for (i = 0; i < od->field_num; ++i) {
    if (!get_field_attr(od, i, DELETE_BIT)) {
      actual_len = fzwrite(od->mel.x_var_attr[i],
      sizeof(uint16_t), od->x_vars, dat_out);
      if (actual_len != od->x_vars) {
        return PREMATURE_DAT_EOF;
      }
    }
  }
  if (IS_O3Q(od)) {
    /*
    if a SRD analysis is present, then save
    relevant data
    */
    if (od->valid & SEED_BIT) {
      fzwrite(&(od->voronoi_num), sizeof(int), 1, dat_out);
      fzwrite(od->mel.voronoi_active, sizeof(int),
        od->voronoi_num, dat_out);
      for (i = 0; i < od->voronoi_num; ++i) {
        fzwrite(&(od->mel.voronoi_fill[i]),
          sizeof(int), 1, dat_out);
        fzwrite(od->al.voronoi_composition[i],
          sizeof(int), od->mel.voronoi_fill[i],
          dat_out);
      }
      for (i = 0; i < od->field_num; ++i) {
        fzwrite(&(od->mel.seed_count[i]),
          sizeof(int), 1, dat_out);
        fzwrite(&(od->mel.seed_count_before_collapse[i]),
          sizeof(int), 1, dat_out);
        fzwrite(&(od->mel.group_zero[i]),
          sizeof(int), 1, dat_out);
        actual_len = fzwrite(od->cimal.voronoi_buf->me[i],
          sizeof(int), od->x_vars, dat_out);
        if (actual_len != od->x_vars) {
          return PREMATURE_DAT_EOF;
        }
      }
    }
  }
  /*
  if dependent variables have been imported,
  then save them too
  */
  if (od->y_vars) {
    /*
    write out the header to the .dat file
    "HEADER                          "
    */
    memset(print_buf, ' ', 32);
    sprintf(print_buf, "%-32s", DAT_HEADER);
    fzwrite(print_buf, 1, 32, dat_out);
    /*
    print out all dependent variable names
    */
    for (i = 0; i < od->y_vars; ++i) {  
      memset(print_buf, ' ', 32);
      snprintf(print_buf, 32, "%-32s", od->cimal.y_var_name->me[i]);
      fzwrite(print_buf, 1, 32, dat_out);
    }
    /*
    then, write a number of doubles equal to y_vars:
    - y_weight_coefficients;
    */
    fzwrite(od->mel.y_data,
      sizeof(YData), od->y_vars, dat_out);
    object_num = 0;
    for (i = 0; i < od->grid.object_num; ++i) {
      if (!get_object_attr(od, i, DELETE_BIT)) {
        /*
        write 1 int:
        - object number
        */
        actual_len = fzwrite(&object_num, sizeof(int), 1, dat_out);
        if (actual_len != 1) {
          return PREMATURE_DAT_EOF;
        }
        /*
        write y_vars
        */
        for (j = 0; j < od->y_vars; ++j) {
          value = (float)get_y_value(od, i, j, 0);
          actual_len = fzwrite(&value, sizeof(float), 1, dat_out);
          if (actual_len != 1) {
            return PREMATURE_DAT_EOF;
          }
        }
        ++object_num;
      }
    }
    /*
    Now save y_var_attr attributes
    */
    actual_len = fzwrite(od->mel.y_var_attr, sizeof(uint16_t),
      od->y_vars, dat_out);
    if (actual_len != od->y_vars) {
      return PREMATURE_DAT_EOF;
    }
  }
  
  if (dat_out) {
    fzclose(dat_out);
  }
  od->file[file_id]->handle = NULL;

  return 0;
}
