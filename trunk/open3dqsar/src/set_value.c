/*

set_value.c

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


#include <include/o3header.h>
#ifdef WIN32
#include <windows.h>
#endif


int check_mmap(O3Data *od, int field_num)
{
  #ifdef WIN32
  char map_name[TITLE_LEN];
  #endif
  int object_num;
  void *mmap_segment;


  if (field_num != od->mmap_field_num) {
    if (od->mmap_field_num != -1) {
      for (object_num = 0; object_num < od->object_num; ++object_num) {
        #ifndef WIN32
        munmap(od->mel.x_var_array[od->mmap_field_num][object_num], od->object_pagesize);
        #else
        UnmapViewOfFile(od->mel.x_var_array[od->mmap_field_num][object_num]);
        CloseHandle(od->al.mol_info[object_num]->hMapHandle);
        #endif
      }
    }
    for (object_num = 0; object_num < od->object_num; ++object_num) {
      #ifndef WIN32
      mmap_segment = mmap(NULL, od->object_pagesize,
        PROT_READ | PROT_WRITE, MAP_SHARED,
        fileno(od->file[TEMP_FIELD_DATA + field_num]->handle),
        od->object_pagesize * object_num);
      if (mmap_segment == MAP_FAILED) {
        return OUT_OF_MEMORY;
      }
      #else
      sprintf(map_name, "Local\\%02d_%04d", field_num + 1, object_num + 1);
      od->al.mol_info[object_num]->hMapHandle = CreateFileMapping
        ((HANDLE)_get_osfhandle(fileno(od->file[TEMP_FIELD_DATA + field_num]->handle)),
        NULL, PAGE_READWRITE, 0, 0, map_name);
      if (od->al.mol_info[object_num]->hMapHandle == NULL) {
        return OUT_OF_MEMORY;
      }
      mmap_segment = MapViewOfFile
        (od->al.mol_info[object_num]->hMapHandle, FILE_MAP_ALL_ACCESS,
        0, (DWORD)(od->object_pagesize * object_num), od->object_pagesize);
      if (mmap_segment == NULL) {
        return OUT_OF_MEMORY;
      }
      #endif
      od->mel.x_var_array[field_num][object_num] = mmap_segment;
      od->mmap_field_num = field_num;
    }
  }
  
  return 0;
}


void sync_field_mmap(O3Data *od)
{
  int field_num;
  
  
  for (field_num = 0; field_num < od->field_num; ++field_num) {
    fflush(od->file[TEMP_FIELD_DATA + field_num]->handle);
  }
}


int set_x_value(O3Data *od, int field_num,
  int object_num, int x_var, double value)
{
  int place;
  int result;
  float float_value;
  
  
  float_value = (float)value;
  place = od->object_pagesize * object_num + x_var * sizeof(float);
  if (od->save_ram) {
    if ((result = check_mmap(od, field_num))) {
      return OUT_OF_MEMORY;
    }
  }
  od->mel.x_var_array
    [field_num][object_num][x_var] = float_value;
  
  return 0;
}


int set_x_value_unbuffered(O3Data *od, int field_num,
  int object_num, int x_var, double value)
{
  int place;
  int actual_len;
  float float_value;
  
  
  float_value = (float)value;
  place = od->object_pagesize * object_num + x_var * sizeof(float);
  if (od->save_ram) {
    if (od->mel.mutex) {
      #ifndef WIN32
      pthread_mutex_lock(od->mel.mutex);
      #else
      WaitForSingleObject(*(od->mel.mutex), INFINITE);
      #endif
    }
    fseek(od->file[TEMP_FIELD_DATA + field_num]->handle, place, SEEK_SET);
    actual_len = fwrite(&float_value, sizeof(float), 1,
      od->file[TEMP_FIELD_DATA + field_num]->handle);
    if (od->mel.mutex) {
      #ifndef WIN32
      pthread_mutex_unlock(od->mel.mutex);
      #else
      ReleaseMutex(*(od->mel.mutex));
      #endif
    }
    if (actual_len != 1) {
      return CANNOT_WRITE_TEMP_FILE;
    }
  }
  else {
    od->mel.x_var_array
      [field_num][object_num][x_var] = float_value;
  }
  
  return 0;
}


void set_y_value(O3Data *od, int object_num, int y_var, double value)
{
  od->mel.y_var_array[od->y_vars
    * object_num + y_var] = (float)value;
}
