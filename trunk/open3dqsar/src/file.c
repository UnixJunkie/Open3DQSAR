/*

file.c

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
#ifdef WIN32
#include <windows.h>
#endif


int dexist(char *dirname)
{
  DIR *handle = NULL;
  
  
  if (dirname && dirname[0]) {
    handle = opendir(dirname);
    if (handle) {
      closedir(handle);
    }
  }
  
  return (handle ? 1 : 0);
}


int fexist(char *filename)
{
  FILE *handle = NULL;
  
  
  if (filename && filename[0]) {
    handle = fopen(filename, "rb");
    if (handle) {
      fclose(handle);
    }
  }
  
  return (handle ? 1 : 0);
}


char *get_basename(char *filename)
{
  int len;
  
  
  if (!filename) {
    return NULL;
  }
  len = strlen(filename);
  while (len && (!ISSEPARATOR(filename[len - 1]))) {
    --len;
  }
  return (&filename[len]);
}


char *get_dirname(char *filename)
{
  int len;
  
  
  if (!filename) {
    return NULL;
  }
  len = strlen(filename);
  while (len && (!ISSEPARATOR(filename[len - 1]))) {
    --len;
  }
  if (!len) {
    filename[0] = '.';
    len = 2;
  }
  filename[len - 1] = '\0';

  return filename;
}


int open_temp_dir(O3Data *od, char *root_dir, char *id_string, char *temp_dir_name)
{
  int result;
  
  
  if (!root_dir) {
    root_dir = od->temp_dir;
  }
  if (!(od->temp_dir)) {
    return CANNOT_WRITE_TEMP_FILE;
  }
  sprintf(temp_dir_name,
    "%s%c%s.%d.%s.XXXXXX", root_dir, SEPARATOR,
    od->package_code, (int)getpid(), id_string);
  errno = 0;
  result = (mkdtemp(temp_dir_name) ? 0 : -1);
  if ((result == -1) && (errno != EEXIST)) {
    return CANNOT_WRITE_TEMP_FILE;
  }

  return 0;
}


int open_perm_dir(O3Data *od, char *root_dir, char *id_string, char *perm_dir_name)
{
  char date_string[MAX_NAME_LEN];
  int result;
  
  
  if (!root_dir) {
    root_dir = od->temp_dir;
  }
  if (!(od->temp_dir)) {
    return CANNOT_WRITE_TEMP_FILE;
  }
  if (!fill_date_string(date_string)) {
    strcpy(date_string, "unknown_time");
  }
  sprintf(perm_dir_name,
    "%s%c%s_%s_%s_XXXXXX", root_dir, SEPARATOR,
    od->package_code, id_string, date_string);
  errno = 0;
  result = (mkdtemp(perm_dir_name) ? 0 : -1);
  if ((result == -1) && (errno != EEXIST)) {
    return CANNOT_WRITE_TEMP_FILE;
  }

  return 0;
}


int open_temp_file(O3Data *od, FileDescriptor *file_descriptor, char *id_string)
{
  int file_des;
  
  
  if (!(od->temp_dir)) {
    return CANNOT_WRITE_TEMP_FILE;
  }
  /*
  if there is one, remove the old temporary file
  connected with this FileDescriptor
  */
  if (file_descriptor->name[0]) {
    remove(file_descriptor->name);
    memset(file_descriptor->name, 0, BUF_LEN);
  }
  sprintf(file_descriptor->name,
    "%s%c%s.%d.%s.XXXXXX", od->temp_dir, SEPARATOR,
    od->package_code, (int)getpid(), id_string);
  errno = 0;
  file_des = mkstemp(file_descriptor->name);
  if ((file_des == -1) && (errno != EEXIST)) {
    return CANNOT_WRITE_TEMP_FILE;
  }
  file_descriptor->handle = fdopen(file_des, "wb+");
  if (!(file_descriptor->handle)) {
    return CANNOT_WRITE_TEMP_FILE;
  }

  return 0;
}


void absolute_path(char *string)
{
  char current_dir[BUF_LEN];
  char copy[BUF_LEN];
  int is_absolute = 0;
  
  
  if (string) {
    #ifdef WIN32
    is_absolute = (isalpha((int)string[0])
      && (string[1] == ':') && strchr("/\\", (int)string[2]));
    #endif
    strcpy(copy, string);
    if (!is_absolute) {
      is_absolute = (strchr("/\\", (int)copy[0]) ? 1 : 0);
    }
    if (!is_absolute) {
      if (getcwd(current_dir, BUF_LEN - 1)) {
        current_dir[BUF_LEN - 1] = '\0';
        sprintf(string, "%s%c%s", current_dir, SEPARATOR,
          ((!strncmp(copy, "./", 2)) || (!strncmp(copy, ".\\", 2))
          ? &copy[2] : copy));
      }
    }
  }
}


int is_in_path(char *program, char *path_to_program)
{
  char *path;
  char *path_copy;
  char *ptr;
  char *context = NULL;
  int found;
  
  
  memset(path_to_program, 0, BUF_LEN);
  if (!(path = getenv(EXE_PATH))) {
    return 0;
  }
  if (!(path_copy = strdup(path))) {
    return 0;
  }
  ptr = strtok_r(path_copy, PATH_SEPARATOR"\n\r\0", &context);
  found = 0;
  while (ptr && (!found)) {
    sprintf(path_to_program, "%s%c%s", ptr, SEPARATOR, program);
    found = fexist(path_to_program);
    ptr = strtok_r(NULL, PATH_SEPARATOR"\n\r\0", &context);
  }
  free(path_copy);
  
  return found;
}


int up_n_levels(char *path, int levels)
{
  int len;
  int result;
  
  
  len = strlen(path);
  while (len && (path[len - 1] == SEPARATOR)) {
    --len;
  }
  while (levels && len) {
    if (path[len - 1] == SEPARATOR) {
      --levels;
      while (len && (path[len - 1] == SEPARATOR)) {
        --len;
      }
    }
    else {
      --len;
    }
  }
  result = ((!levels) && len);
  if (result) {
    path[len] = '\0';
  }
  
  return result;
}


int zipFiletime(char *filename, zip_fileinfo *zfi)
{
  int ret = 0;
  #ifdef WIN32
  FILETIME ft;
  HANDLE handle;
  WIN32_FIND_DATAA fd;

  
  handle = FindFirstFileA(filename, &fd);
  if (handle != INVALID_HANDLE_VALUE) {
    FileTimeToLocalFileTime(&(fd.ftLastWriteTime), &ft);
    FileTimeToDosDateTime(&ft, ((LPWORD)&(zfi->dosDate)) + 1,
      ((LPWORD)&(zfi->dosDate)) + 0);
    FindClose(handle);
    ret = 1;
  }
  #else
  char name[BUF_LEN + 1];
  int len = 0;
  struct stat s;
  struct tm filedate;
  time_t tm_t = 0;


  if (strcmp(filename, "-")) {
    len = strlen(filename);
    if (len > BUF_LEN) {
      len = BUF_LEN;
    }
    strncpy(name, filename, BUF_LEN - 1);
    name[BUF_LEN] = '\0';

    if (name[len - 1] == '/') {
      name[len - 1] = '\0';
    }
    if (stat(name, &s) == 0) {
      tm_t = s.st_mtime;
      ret = 1;
    }
  }
  memset(&filedate, 0, sizeof(struct tm));
  localtime_r(&tm_t, &filedate);

  zfi->tmz_date.tm_sec  = filedate.tm_sec;
  zfi->tmz_date.tm_min  = filedate.tm_min;
  zfi->tmz_date.tm_hour = filedate.tm_hour;
  zfi->tmz_date.tm_mday = filedate.tm_mday;
  zfi->tmz_date.tm_mon  = filedate.tm_mon ;
  zfi->tmz_date.tm_year = filedate.tm_year;
  #endif
  
  return ret;
}


char *get_basename_no_ext(char *filename)
{
  char *basename_no_ext;
  int len;
  int found = 0;
  

  len = strlen(filename);
  while (len && (!ISSEPARATOR(filename[len - 1]))) {
    --len;
  }
  if (!(basename_no_ext = strdup(&filename[len]))) {
    return NULL;
  }
  len = strlen(basename_no_ext);
  while (len && (!(found = (basename_no_ext[len - 1] == '.')))) {
    --len;
  }
  if (found) {
    basename_no_ext[len - 1] = '\0';
  }
  
  return basename_no_ext;
}


zipFile zipOpenWrite(char *filename)
{
  char *basename_no_zip;
  int ret;
  zipFile handle;
  zip_fileinfo zfi;
  
  
  memset(&zfi, 0, sizeof(zip_fileinfo));
  if (!(basename_no_zip = get_basename_no_ext(filename))) {
    return NULL;
  }
  handle = zipOpen(filename, APPEND_STATUS_CREATE);
  if (handle) {
    zipFiletime(filename, &zfi);
    ret = zipOpenNewFileInZip(handle, basename_no_zip, &zfi,
      NULL, 0, NULL, 0, NULL, Z_DEFLATED, Z_DEFAULT_COMPRESSION);
    free(basename_no_zip);
    if (ret) {
      zipClose(handle, NULL);
      handle = NULL;
    }
  }

  return handle;
}


unzFile zipOpenRead(char *filename)
{
  int ret;
  unzFile handle;
  
  
  handle = unzOpen(filename);
  if (handle) {
    ret = unzGoToFirstFile(handle);
    if (!ret) {
      ret = unzOpenCurrentFile(handle);
    }
    if (ret) {
      unzClose(handle);
      handle = NULL;
    }
  }

  return handle;
}


int zipCloseWrite(unzFile handle)
{
  int ret;


  if (!(ret = zipCloseFileInZip(handle))) {
    ret = zipClose(handle, NULL);
  }

  return ret;
}


int zipCloseRead(zipFile handle)
{
  int ret;


  if (!(ret = ((unzCloseCurrentFile(handle) == UNZ_CRCERROR) ? 1 : 0))) {
    ret = unzClose(handle);
  }

  return ret;
}


fzPtr *fzopen(char *filename, char *mode)
{
  char rw_mode[8];
  int len;
  fzPtr *fz_ptr;
  
  
  fz_ptr = (fzPtr *)malloc(sizeof(fzPtr));
  if (!fz_ptr) {
    return NULL;
  }
  memset(fz_ptr, 0, sizeof(fzPtr));
  if (tolower(mode[0]) == 'w') {
    strcpy(rw_mode, "wb");
    fz_ptr->zip_type |= ZIP_MODE_WRITE;
  }
  else {
    strcpy(rw_mode, "rb");
    fz_ptr->zip_type |= ZIP_MODE_READ;
  }
  len = strlen(filename);
  if ((len >= 3) && (!strncasecmp(&filename[len - 3], ".gz", 3))) {
    fz_ptr->zip_type |= GZIP_FILE_HANDLE;
    fz_ptr->gzip_file_handle = gzopen(filename, rw_mode);
  }
  else if ((len >= 4) && (!strncasecmp(&filename[len - 4], ".zip", 4))) {
    fz_ptr->zip_type |= ZIP_FILE_HANDLE;
    if (fz_ptr->zip_type & ZIP_MODE_WRITE) {
      fz_ptr->zip_file_handle = zipOpenWrite(filename);
    }
    else {
      fz_ptr->unz_file_handle = zipOpenRead(filename);
    }
  }
  else {
    fz_ptr->zip_type |= NORMAL_FILE_HANDLE;
    fz_ptr->normal_file_handle = fopen(filename, rw_mode);
  }
  if ((!(fz_ptr->gzip_file_handle))
    && (!(fz_ptr->zip_file_handle))
    && (!(fz_ptr->unz_file_handle))
    && (!(fz_ptr->normal_file_handle))) {
    free(fz_ptr);
    return NULL;
  }
  
  return fz_ptr;
}


int fzclose(fzPtr *fz_ptr)
{
  int ret = 0;
  
  if (fz_ptr) {
    if ((fz_ptr->zip_type & GZIP_FILE_HANDLE)
      && fz_ptr->gzip_file_handle) {
      ret = gzclose(fz_ptr->gzip_file_handle);
    }
    else if ((fz_ptr->zip_type & ZIP_FILE_HANDLE)
      && (fz_ptr->zip_type & ZIP_MODE_WRITE)
      && fz_ptr->zip_file_handle) {
      ret = zipCloseWrite(fz_ptr->zip_file_handle);
    }
    else if ((fz_ptr->zip_type & ZIP_FILE_HANDLE)
      && (fz_ptr->zip_type & ZIP_MODE_READ)
      && fz_ptr->unz_file_handle) {
      ret = zipCloseRead(fz_ptr->unz_file_handle);
      if (fz_ptr->buf) {
        free(fz_ptr->buf);
      }
    }
    else if ((fz_ptr->zip_type & NORMAL_FILE_HANDLE)
      && fz_ptr->normal_file_handle) {
      ret = fclose(fz_ptr->normal_file_handle);
    }
    free(fz_ptr);
  }
  
  return ret;
}


int fzputs(fzPtr *fz_ptr, char *data)
{
  char cr[4];
  int cr_len;
  int ret = 0;
  int real_len;


  if ((!fz_ptr) || (!data)) {
    return 0;
  }
  real_len = strlen(data);
  #ifdef WIN32
  strcpy(cr, "\r\n");
  cr_len = 2;
  #else
  strcpy(cr, "\n");
  cr_len = 1;
  #endif
  if ((fz_ptr->zip_type & GZIP_FILE_HANDLE)
    && fz_ptr->gzip_file_handle) {
    ret = gzwrite(fz_ptr->gzip_file_handle, data, real_len);
    if (ret) {
      ret = gzwrite(fz_ptr->gzip_file_handle, cr, cr_len);
    }
    if (!ret) {
      ret = EOF;
    }
  }
  else if ((fz_ptr->zip_type & ZIP_FILE_HANDLE)
    && (fz_ptr->zip_type & ZIP_MODE_WRITE)
    && fz_ptr->zip_file_handle) {
    ret = (zipWriteInFileInZip
      (fz_ptr->zip_file_handle, data, real_len) ? EOF : 1);
    if (ret) {
      ret = (zipWriteInFileInZip
        (fz_ptr->zip_file_handle, cr, cr_len) ? EOF : 1);
    }
  }
  else if ((fz_ptr->zip_type & NORMAL_FILE_HANDLE)
    && fz_ptr->normal_file_handle) {
    ret = fwrite(data, 1, real_len, fz_ptr->normal_file_handle);
    if (ret == real_len) {
      ret = fwrite(cr, 1, cr_len, fz_ptr->normal_file_handle);
      if (ret != cr_len) {
        ret = EOF;
      }
    }
    else {
      ret = EOF;
    }
  }
    
  return ret;
}


char *fzgets(char *data, int len, fzPtr *fz_ptr)
{
  char *ret = NULL;
  int real_len = 0;
  int i = 0;
  int eol = 0;
  int crfound = 0;


  if ((!fz_ptr) || (!data)) {
    return NULL;
  }
  if (!(fz_ptr->buf)) {
    fz_ptr->buf = malloc(FZ_BUF_LEN);
    if (!(fz_ptr->buf)) {
      return NULL;
    }
  }
  while ((!eol) && (i < len)) {
    if (fz_ptr->pos >= fz_ptr->data_len) {
      if ((fz_ptr->zip_type & GZIP_FILE_HANDLE)
        && fz_ptr->gzip_file_handle) {
        real_len = gzread
          (fz_ptr->gzip_file_handle, fz_ptr->buf, FZ_BUF_LEN);
        if (real_len <= 0) {
          if (real_len == 0) {
            ret = (crfound ? data : NULL);
          }
          return ret;
        }
      }
      else if ((fz_ptr->zip_type & ZIP_FILE_HANDLE)
        && (fz_ptr->zip_type & ZIP_MODE_READ)
        && fz_ptr->unz_file_handle) {
        real_len = unzReadCurrentFile
          (fz_ptr->unz_file_handle, fz_ptr->buf, FZ_BUF_LEN);
        if (real_len <= 0) {
          if (real_len == 0) {
            ret = (crfound ? data : NULL);
          }
          return ret;
        }
      }
      else if ((fz_ptr->zip_type & NORMAL_FILE_HANDLE)
        && fz_ptr->normal_file_handle) {
        real_len = fread(fz_ptr->buf, 1, FZ_BUF_LEN,
          fz_ptr->normal_file_handle);
        if (real_len == 0) {
          if (feof(fz_ptr->normal_file_handle)) {
            ret = (crfound ? data : NULL);
          }
          return ret;
        }
      }
      fz_ptr->data_len = real_len;
      fz_ptr->pos = 0;
    }
    if (fz_ptr->buf[fz_ptr->pos] == '\r') {
      ++crfound;
    }
    eol = ((fz_ptr->buf[fz_ptr->pos] == '\n')
      || (crfound > 1));
    if ((!eol) && (crfound == 1)
      && (fz_ptr->buf[fz_ptr->pos] != '\r')) {
      ++crfound;
    }
    if (eol || (crfound > 1)) {
      if (crfound > 1) {
        --(fz_ptr->pos);
      }
      crfound = 0;
      fz_ptr->buf[fz_ptr->pos] = '\n';
      eol = 1;
    }
    if ((!eol) && (!crfound)) {
      data[i] = fz_ptr->buf[fz_ptr->pos];
      ++i;
    }
    ++(fz_ptr->pos);
  }
    
  return data;
}


int fzwrite(void *data, size_t size, size_t count, fzPtr *fz_ptr)
{
  int ret;
  int real_len;


  ret = count;
  real_len = count * size;
  if ((!fz_ptr) || (!data)) {
    return -1;
  }
  if ((fz_ptr->zip_type & GZIP_FILE_HANDLE)
    && fz_ptr->gzip_file_handle) {
    ret = gzwrite(fz_ptr->gzip_file_handle, data, real_len);
    if (!ret) {
      ret = -1;
    }
  }
  else if ((fz_ptr->zip_type & ZIP_FILE_HANDLE)
    && (fz_ptr->zip_type & ZIP_MODE_WRITE)
    && fz_ptr->zip_file_handle) {
    ret = (zipWriteInFileInZip
      (fz_ptr->zip_file_handle, data, real_len) ? -1 : count);
  }
  else if ((fz_ptr->zip_type & NORMAL_FILE_HANDLE)
    && fz_ptr->normal_file_handle) {
    ret = fwrite(data, size, count, fz_ptr->normal_file_handle);
  }
    
  return ret;
}


int fzread(void *data, size_t size, size_t count, fzPtr *fz_ptr)
{
  int ret = 0;
  int real_len;
  int len;


  if ((!fz_ptr) || (!data)) {
    return 0;
  }
  real_len = count * size;
  if ((fz_ptr->zip_type & GZIP_FILE_HANDLE)
    && fz_ptr->gzip_file_handle) {
    ret = (((len = gzread(fz_ptr->gzip_file_handle, data, real_len)) ==
      real_len) ? count : (len / size));
  }
  else if ((fz_ptr->zip_type & ZIP_FILE_HANDLE)
    && (fz_ptr->zip_type & ZIP_MODE_READ)
    && fz_ptr->unz_file_handle) {
    ret = (((len = unzReadCurrentFile(fz_ptr->unz_file_handle,
      data, real_len)) == real_len) ? count : (len / size));
  }
  else if ((fz_ptr->zip_type & NORMAL_FILE_HANDLE)
    && fz_ptr->normal_file_handle) {
    ret = fread(data, size, count, fz_ptr->normal_file_handle);
  }
    
  return ret;
}


int fzseek(fzPtr *fz_ptr, long int offset, int whence)
{
  char buffer[LARGE_BUF_LEN];
  int i = 0;
  int times;
  int actual_len;
  int ret = 0;
  
  
  if ((!fz_ptr) || (whence != SEEK_CUR) || (offset < 0)) {
    return -1;
  }
  times = offset / LARGE_BUF_LEN + 1;
  while ((!ret) && (i < times)) {
    actual_len = ((i == (times - 1))
      ? (offset % LARGE_BUF_LEN) : LARGE_BUF_LEN);
    ret = ((fzread(buffer, 1, actual_len, fz_ptr) == actual_len) ? 0 : -1);
    ++i;
  }
  
  return ret;
}


void fzrewind(fzPtr *fz_ptr)
{
  if (fz_ptr) {
    if ((fz_ptr->zip_type & GZIP_FILE_HANDLE)
      && fz_ptr->gzip_file_handle) {
      gzrewind(fz_ptr->gzip_file_handle);
    }
    else if ((fz_ptr->zip_type & ZIP_FILE_HANDLE)
      && (fz_ptr->zip_type & ZIP_MODE_READ)
      && fz_ptr->unz_file_handle) {
      if (!unzCloseCurrentFile(fz_ptr->unz_file_handle)) {
        unzOpenCurrentFile(fz_ptr->unz_file_handle);
      }
    }
    else if ((fz_ptr->zip_type & NORMAL_FILE_HANDLE)
      && fz_ptr->normal_file_handle) {
      rewind(fz_ptr->normal_file_handle);
    }
    fz_ptr->data_len = 0;
    fz_ptr->pos = 0;
  }
}
