/*

grid_write.c

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

#include <math.h>
#include <include/o3header.h>
#include <include/interpolation_coefficients.h>
#include <include/interpolated_value.h>


void copy_plane_to_buffer(O3Data *od,
  float *float_xy_mat, float *buf_float_xy_mat)
{
  int x;
  int y;
  
  
  /*
  the original plane is:
  
  #abcdef#
  #ghijkl#
  #mnopqr#
  
  the copy plane is:
  
  aabcdeff
  aabcdeff
  gghijkll
  mmnopqrr
  mmnopqrr
  
  */
  /*
  fill in the inner rectangle, leaving an empty frame
  */
  for (y = 0; y < od->grid.nodes[1]; ++y) {
    for (x = 0; x < od->grid.nodes[0]; ++x) {
      buf_float_xy_mat[(od->grid.nodes[0] + 2)
        * (y + 1) + x + 1] =
        float_xy_mat[(od->grid.nodes[0] + 2)
        * y + x + 1];
    }
  }
  /*
  fill the side columns
  */
  for (y = 1; y < (od->grid.nodes[1] + 1); ++y) {
    buf_float_xy_mat[(od->grid.nodes[0] + 2) * y] =
      buf_float_xy_mat[(od->grid.nodes[0] + 2) * y + 1];
    buf_float_xy_mat[(od->grid.nodes[0] + 2) *
      y + od->grid.nodes[0] + 1] =
      buf_float_xy_mat[(od->grid.nodes[0] + 2)
      * y + od->grid.nodes[0]];
  }
  /*
  fill the first and the last row
  */
  for (x = 0; x < (od->grid.nodes[0] + 2); ++x) {
    buf_float_xy_mat[x] =
      buf_float_xy_mat[od->grid.nodes[0] + 2 + x];
    buf_float_xy_mat[(od->grid.nodes[0] + 2)
      * (od->grid.nodes[1] + 1) + x] =
      buf_float_xy_mat[(od->grid.nodes[0] + 2)
        * od->grid.nodes[1] + x];
  }
}


int write_grid_plane(O3Data *od, FILE *plane_file,
  int z, int interpolate, int swap_endianness)
{
  int i;
  int actual_len;
  int y_node;
  int x = 0;
  int y = 0;
  int xx;
  int yy;
  int zz;
  int n;
  int nx;
  int ny;
  int nz;
  int buf_xy_width;
  int buf_xy_height;
  int out_xy_width;
  int out_xy_height;
  int *int_matrix;
  float *float_xy_mat;
  float *buf_float_xy_mat[4];
  float *out_float_xy_mat = NULL;
  double x_incr;
  double y_incr;
  double z_incr;
  double f_vec[64];
  double ic_vec[64];

  
  float_xy_mat = od->mel.float_xy_mat;
  if (!interpolate) {
    /*
    write meta-data
    */
    int_matrix = (int *)float_xy_mat;
    for (y_node = 0; y_node < od->grid.nodes[1]; ++y_node) {
      int_matrix[(od->grid.nodes[0] + 2) * y_node] =
        (od->grid.nodes[0]) * sizeof(int);
      int_matrix[(od->grid.nodes[0] + 2) 
        * y_node + (od->grid.nodes[0]) + 1] =
        (od->grid.nodes[0]) * sizeof(int);
    }
    /*
    fix endianness and write out to file
    */
    fix_endianness(float_xy_mat, sizeof(float),
      (od->grid.nodes[0] + 2) * od->grid.nodes[1], swap_endianness);
    actual_len = fwrite(float_xy_mat, sizeof(float),
      (od->grid.nodes[0] + 2)
      * od->grid.nodes[1], plane_file);
    if (actual_len != (od->grid.nodes[0] + 2)
      * od->grid.nodes[1]) {
      O3_ERROR_LOCATE(&(od->task));
      return PREMATURE_EOF;
    }
  }
  else {
    buf_xy_width = od->grid.nodes[0] + 2;
    buf_xy_height = od->grid.nodes[1] + 2;
    out_xy_width = od->grid.nodes[0]
      + interpolate * (od->grid.nodes[0] - 1) + 2;
    out_xy_height = od->grid.nodes[1]
      + interpolate * (od->grid.nodes[1] - 1);
    out_float_xy_mat = od->mel.out_float_xy_mat;
    int_matrix = (int *)out_float_xy_mat;
    i = 0;
    while (i < 4) {
      buf_float_xy_mat[i] = od->mel.buf_float_xy_mat[i];
      if (!buf_float_xy_mat[i]) {
        buf_float_xy_mat[i] = (float *)malloc
          (sizeof(float) * buf_xy_width * buf_xy_height);
        if (!buf_float_xy_mat[i]) {
          return OUT_OF_MEMORY;
        }
        od->mel.buf_float_xy_mat[i] = buf_float_xy_mat[i];
        /*
        if it is plane 0, then duplicate it in buf_float_xy_mat[0]
        and buf_float_xy_mat[1], since it lies on the edge.
        */
        if (!i) {
          out_float_xy_mat = (float *)malloc
            (sizeof(float) * out_xy_width * out_xy_height);
          if (!out_float_xy_mat) {
            return OUT_OF_MEMORY;
          }
          od->mel.out_float_xy_mat = out_float_xy_mat;
          memset(out_float_xy_mat, 0,
            sizeof(float) * out_xy_width * out_xy_height);
          copy_plane_to_buffer(od,
            float_xy_mat, buf_float_xy_mat[0]);
        }
        else {
          copy_plane_to_buffer(od,
            float_xy_mat, buf_float_xy_mat[i]);
          break;
        }
      }
      ++i;
    }
    if (i < 4) {
      return 0;
    }
    /*
    if all buffers have already been filled, then
    do the interpolation and write out the interpolated plane
    */
    while (i) {
      z_incr = (double)0;
      for (zz = 0; zz <= interpolate; ++zz, z_incr += safe_rint
        (((double)1 / (double)(interpolate + 1)) * 1.0e05) / 1.0e05) {
        /*
        if this is the last z plane, then do not interpolate
        anymore, or the grid boundaries will be crossed
        */
        if (zz && (i == 1)) {
          break;
        }
        /*
        build the interpolated grid
        */
        for (y = 1; y < od->grid.nodes[1]; ++y) {
          for (x = 1; x < od->grid.nodes[0]; ++x) {
            n = 0;
            for (nz = 0; nz <= 1; ++nz) {
              for (ny = 0; ny <= 1; ++ny) {
                for (nx = 0; nx <= 1; ++nx) {
                  /*
                  calculate function values and its derivatives
                  in the 8 corners of the cube
                  */
                  f_vec[IV_FUNC + n] = (double)buf_float_xy_mat
                    [nz + 1][buf_xy_width * (y + ny) + x + nx];
                  f_vec[IV_DFDX + n] = ((double)buf_float_xy_mat
                    [nz + 1][buf_xy_width * (y + ny) + x + nx + 1]
                    - (double)buf_float_xy_mat
                    [nz + 1][buf_xy_width * (y + ny) + x + nx - 1]) / 2;
                  f_vec[IV_DFDY + n] = ((double)buf_float_xy_mat
                    [nz + 1][buf_xy_width * (y + ny + 1) + x + nx]
                    - (double)buf_float_xy_mat
                    [nz + 1][buf_xy_width * (y + ny - 1) + x + nx]) / 2;
                  f_vec[IV_DFDZ + n] = ((double)buf_float_xy_mat
                    [nz + 2][buf_xy_width * (y + ny) + x + nx]
                    - (double)buf_float_xy_mat
                    [nz][buf_xy_width * (y + ny) + x + nx]) / 2;
                  f_vec[IV_D2FDXDY + n] = (((double)buf_float_xy_mat
                    [nz + 1][buf_xy_width * (y + ny + 1) + x + nx + 1]
                    - (double)buf_float_xy_mat
                    [nz + 1][buf_xy_width * (y + ny - 1) + x + nx + 1])
                    - ((double)buf_float_xy_mat
                    [nz + 1][buf_xy_width * (y + ny + 1) + x + nx - 1]
                    - (double)buf_float_xy_mat
                    [nz + 1][buf_xy_width * (y + ny - 1) + x + nx - 1])) / 4;
                  f_vec[IV_D2FDXDZ + n] = (((double)buf_float_xy_mat
                    [nz + 2][buf_xy_width * (y + ny) + x + nx + 1]
                    - (double)buf_float_xy_mat
                    [nz][buf_xy_width * (y + ny) + x + nx + 1])
                    - ((double)buf_float_xy_mat
                    [nz + 2][buf_xy_width * (y + ny) + x + nx - 1]
                    - (double)buf_float_xy_mat
                    [nz][buf_xy_width * (y + ny) + x + nx - 1])) / 4;
                  f_vec[IV_D2FDYDZ + n] = (((double)buf_float_xy_mat
                    [nz + 2][buf_xy_width * (y + ny + 1) + x + nx]
                    - (double)buf_float_xy_mat
                    [nz][buf_xy_width * (y + ny + 1) + x + nx])
                    - ((double)buf_float_xy_mat
                    [nz + 2][buf_xy_width * (y + ny - 1) + x + nx]
                    - (double)buf_float_xy_mat
                    [nz][buf_xy_width * (y + ny - 1) + x + nx])) / 4;
                  f_vec[IV_D3FDXDYDZ + n] = ((((double)buf_float_xy_mat
                    [nz + 2][buf_xy_width * (y + ny + 1) + x + nx + 1]
                    - (double)buf_float_xy_mat
                    [nz][buf_xy_width * (y + ny + 1) + x + nx + 1])
                    - ((double)buf_float_xy_mat
                    [nz + 2][buf_xy_width * (y + ny - 1) + x + nx + 1]
                    - (double)buf_float_xy_mat
                    [nz][buf_xy_width * (y + ny - 1) + x + nx + 1]))
                    - (((double)buf_float_xy_mat
                    [nz + 2][buf_xy_width * (y + ny + 1) + x + nx - 1]
                    - (double)buf_float_xy_mat
                    [nz][buf_xy_width * (y + ny + 1) + x + nx - 1])
                    - ((double)buf_float_xy_mat
                    [nz + 2][buf_xy_width * (y + ny - 1) + x + nx - 1]
                    - (double)buf_float_xy_mat
                    [nz][buf_xy_width * (y + ny - 1) + x + nx - 1]))) / 8;
                  ++n;
                }
              }
            }
            cblas_dgemv(CblasColMajor, CblasNoTrans,
              64, 64, 1.0, ic_mat, 64, f_vec,
              1, 0.0, ic_vec, 1);
            
            for (yy = 0, y_incr = 0.0; yy <= interpolate; ++yy, y_incr += safe_rint
              (((double)1 / (double)(interpolate + 1)) * 1.0e05) / 1.0e05) {
              for (xx = 0, x_incr = 0.0; xx <= interpolate; ++xx, x_incr += safe_rint
                (((double)1 / (double)(interpolate + 1)) * 1.0e05) / 1.0e05) {
                /*
                if the grid boundary was reached in one direction,
                then do not interpolate in that direction
                or it will be crossed
                */
                if (((x == (od->grid.nodes[0] - 1)) && xx)
                  || ((y == (od->grid.nodes[1] - 1)) && yy)
                  || ((z == (od->grid.nodes[2] - 1)) && zz)) {
                  continue;
                }
                /*
                write the interpolated value in the grid
                */
                out_float_xy_mat[out_xy_width
                  * ((y - 1) * (interpolate + 1) + yy)
                  + (x - 1) * (interpolate + 1) + xx] =
                  (float)interpolated_value
                  (ic_vec, x_incr, y_incr, z_incr);
              }
            }
          }
        }
        /*
        write meta-data
        */
        for (y_node = 0; y_node < (od->grid.nodes[1] + interpolate
          * (od->grid.nodes[1] - 1)); ++y_node) {
          int_matrix[out_xy_width * y_node] =
            (out_xy_width - 2) * sizeof(int);
          int_matrix[out_xy_width * (y_node + 1) - 1] =
            (out_xy_width - 2) * sizeof(int);
        }
        /*
        fix endianness and write out to file
        */
        fix_endianness(out_float_xy_mat, sizeof(float),
          out_xy_width * out_xy_height, swap_endianness);
        actual_len = fwrite(out_float_xy_mat, sizeof(float),
          out_xy_width * out_xy_height, plane_file);
        if (actual_len != (out_xy_width * out_xy_height)) {
          O3_ERROR_LOCATE(&(od->task));
          return PREMATURE_EOF;
        }
      }
      /*
      shift up buffers:
      1->0, 2->1, 3->2
      and fill buffer 3 with the new matrix
      */
      memcpy(buf_float_xy_mat[0], buf_float_xy_mat[1],
        sizeof(float) * buf_xy_width * buf_xy_height);
      memcpy(buf_float_xy_mat[1], buf_float_xy_mat[2],
        sizeof(float) * buf_xy_width * buf_xy_height);
      memcpy(buf_float_xy_mat[2], buf_float_xy_mat[3],
        sizeof(float) * buf_xy_width * buf_xy_height);
      copy_plane_to_buffer(od, float_xy_mat, buf_float_xy_mat[3]);
      /*
      if the last plane has not been reached yet,
      then exit, otherwise, flush buffers
      */
      if (z != (od->grid.nodes[2] - 1)) {
        break;
      }
      --i;
    }
    /*
    if buffers were flushed, memory can be freed
    */
    if (!i) {
      for (i = 0; i < 4; ++i) {
        if (buf_float_xy_mat[i]) {
          free(buf_float_xy_mat[i]);
          od->mel.buf_float_xy_mat[i] = NULL;
        }
      }
      if (out_float_xy_mat) {
        free(out_float_xy_mat);
        od->mel.out_float_xy_mat = NULL;
      }
    }
  }
  
  return 0;  
}


int write_header(O3Data *od, char *header, int format, int interpolate, int swap_endianness)
{
  int i;
  int j;
  int len;
  int actual_len;
  int npts[3];
  int *int_meta_data;
  float *float_meta_data;
  
  
  for (i = 0; i < 3; ++i) {
    npts[i] = od->grid.nodes[i] + interpolate * (od->grid.nodes[i] - 1);
  }
  switch (format) {
    case INSIGHT_FORMAT:
    int_meta_data = (int *)header;
    float_meta_data = (float *)header;
    /*
    write in the header all information present
    in the od structure
    (see INSIGHT .grd format documentation,
    http://www.scripps.edu/rc/softwaredocs/msi/insight2K/insight/Utilities.html)
    */
    /*
    this integer is useful to determine data endianness;
    it describes the length of the title field
    (TITLE_LEN bytes in our case)
    */
    int_meta_data[0] = TITLE_LEN;
    /*
    this integer closes the title field; again
    its value is equal to the length of the title field
    (TITLE_LEN bytes in our case)
    */
    int_meta_data[16] = TITLE_LEN;
    /*
    this indicates the number of meta-data bytes which follow
    before real data come (METADATA_LEN bytes in our case)
    */
    int_meta_data[17] = METADATA_LEN;
    /*
    this means that x is the fastest-varying coordinate
    */
    int_meta_data[18] = 0;
    /*
    this means that data chunks are 4-byte long
    */
    int_meta_data[19] = 4;
    /*
    data type: 0 = floating point, 1 = integer
    */
    int_meta_data[20] = 0;
    /*
    Cell_x_length
    */
    float_meta_data[21] = od->grid.end_coord[0];
    /*
    Cell_y_length
    */
    float_meta_data[22] = od->grid.end_coord[1];
    /*
    Cell_z_length
    */
    float_meta_data[23] = od->grid.end_coord[2];
    /*
    Cell_x_angle
    */
    float_meta_data[24] = 90.0;
    /*
    Cell_y_angle
    */
    float_meta_data[25] = 90.0;
    /*
    Cell_z_angle
    */
    float_meta_data[26] = 90.0;
    /*
    Cell_x_start
    */
    float_meta_data[27] = od->grid.start_coord[0] / od->grid.end_coord[0];
    /*
    Cell_x_end
    */
    float_meta_data[28] = 1.0;
    /*
    Cell_y_start
    */
    float_meta_data[29] = od->grid.start_coord[1] / od->grid.end_coord[1];
    /*
    Cell_y_end
    */
    float_meta_data[30] = 1.0;
    /*
    Cell_z_start
    */
    float_meta_data[31] = od->grid.start_coord[2] / od->grid.end_coord[2];
    /*
    Cell_z_end
    */
    float_meta_data[32] = 1.0;
    /*
    Cell_x_intervals
    */
    int_meta_data[33] = npts[0] - 1;
    /*
    Cell_y_intervals
    */
    int_meta_data[34] = npts[1] - 1;
    /*
    Cell_z_intervals
    */
    int_meta_data[35] = npts[2] - 1;
    /*
    this indicates the number of meta-data bytes which preceded
    before real data come (METADATA_LEN bytes in our case)
    */
    int_meta_data[36] = METADATA_LEN;
    /*
    write the title
    */
    memset(&header[4], ' ', TITLE_LEN);
    memcpy(&header[4], INSIGHT_GRID_TITLE, strlen(INSIGHT_GRID_TITLE));
    /*
    make sure that the title is null-terminated
    */
    header[4 + TITLE_LEN - 1] = '\0';
    /*
    then convert the header numerical data according
    to the requested endianness
    */
    fix_endianness(header, sizeof(int), 1, swap_endianness);
    fix_endianness(header + TITLE_LEN + 4, sizeof(int),
      METADATA_LEN / 4 + 3, swap_endianness);

    actual_len = fwrite(header, 1, 16 + TITLE_LEN + METADATA_LEN,
      od->file[GRD_OUT]->handle);
    if (actual_len != (16 + TITLE_LEN + METADATA_LEN)) {
      O3_ERROR_LOCATE(&(od->task));
      return PREMATURE_EOF;
    }
    break;

    case MOE_FORMAT:
    int_meta_data = (int *)header;
    float_meta_data = (float *)header;
    sprintf(header, "%s%.2f", MOE_ID_TOKEN, MOE_FORMAT_VERSION);
    float_meta_data[3] = MOE_FORMAT_VERSION;
    int_meta_data[4] = 0;
    len = strlen(MOE_GRID_TITLE);
    int_meta_data[5] = len;
    fix_endianness(&header[12], sizeof(int), 3, swap_endianness);
    strcpy(&header[24], MOE_GRID_TITLE);
    len += 24;
    int_meta_data = (int *)(&header[len]);
    int_meta_data[0] = 3;
    len += sizeof(int);
    fix_endianness(int_meta_data, sizeof(int), 1, swap_endianness);
    actual_len = fwrite(header, 1, len, od->file[GRD_OUT]->handle);
    if (actual_len != len) {
      O3_ERROR_LOCATE(&(od->task));
      return PREMATURE_EOF;
    }
    int_meta_data = (int *)header;
    for (i = 0; i < 3; ++i) {
      int_meta_data[0] = npts[i];
      fix_endianness(int_meta_data, sizeof(int), 1, swap_endianness);
      actual_len = fwrite(int_meta_data, sizeof(int), 1, od->file[GRD_OUT]->handle);
      if (actual_len != 1) {
        O3_ERROR_LOCATE(&(od->task));
        return PREMATURE_EOF;
      }
      for (j = 0; j < npts[i]; ++j) {
        float_meta_data[0] = (float)safe_rint((double)j
          * (double)(od->grid.step[i]) / (double)(interpolate + 1)
          * 1.0e04) / 1.0e04 + (double)(od->grid.start_coord[i]);
        fix_endianness(float_meta_data, sizeof(float), 1, swap_endianness);
        actual_len = fwrite(float_meta_data, sizeof(float), 1, od->file[GRD_OUT]->handle);
        if (actual_len != 1) {
          O3_ERROR_LOCATE(&(od->task));
          return PREMATURE_EOF;
        }
      }
    }
    int_meta_data[0] = 1;
    int_meta_data[1] = 0;
    fix_endianness(int_meta_data, sizeof(int), 2, swap_endianness);
    actual_len = fwrite(int_meta_data, sizeof(int), 2, od->file[GRD_OUT]->handle);
    if (actual_len != 2) {
      O3_ERROR_LOCATE(&(od->task));
      return PREMATURE_EOF;
    }
    break;

    case SYBYL_FORMAT:
    fprintf(od->file[GRD_OUT]->handle,
      "%s\n", SYBYL_GRID_TITLE);
    for (i = 0; i < 3; ++i) {
      fprintf(od->file[GRD_OUT]->handle,
        "%12.6f%12.6lf%12d\n", od->grid.start_coord[i],
        (double)((double)(od->grid.step[i]) / (double)(interpolate + 1)), npts[i]);
    }
    fprintf(od->file[GRD_OUT]->handle,
      "%12.1f%12.1f%12.1f\n", 90.0, 90.0, 90.0);
    
    break;

    case MAESTRO_FORMAT:
    fprintf(od->file[GRD_OUT]->handle, "&plot\n");
    fprintf(od->file[GRD_OUT]->handle, "iplot=  %12d\n", 1);
    fprintf(od->file[GRD_OUT]->handle,
      "origin= %12.6f%12.6f%12.6f\n",
      od->grid.start_coord[0] / BOHR_RADIUS,
      od->grid.start_coord[1] / BOHR_RADIUS,
      od->grid.start_coord[2] / BOHR_RADIUS);
    fprintf(od->file[GRD_OUT]->handle,
      "extentx=%12.6f%12.6f%12.6f\n",
      od->grid.step[0] * od->grid.nodes[0] / BOHR_RADIUS, 0.0, 0.0);
    fprintf(od->file[GRD_OUT]->handle,
      "extenty=%12.6f%12.6f%12.6f\n",
      0.0, od->grid.step[1] * od->grid.nodes[1] / BOHR_RADIUS, 0.0);
    fprintf(od->file[GRD_OUT]->handle,
      "extentz=%12.6f%12.6f%12.6f\n",
      0.0, 0.0, od->grid.step[2] * od->grid.nodes[2] / BOHR_RADIUS);
    fprintf(od->file[GRD_OUT]->handle,
      "npts=   %12d%12d%12d\n", npts[0], npts[1], npts[2]);
    fprintf(od->file[GRD_OUT]->handle, "&end\n");
    
    break;

    default:
    for (i = 0; i < 3; ++i) {
      fprintf(od->file[GRD_OUT]->handle,
        "%c start, %c end coordinates:    %.4f, %.4f\n",
        'X' + i, 'X' + i,
        od->grid.start_coord[i],
        od->grid.end_coord[i]);
    }
    fprintf(od->file[GRD_OUT]->handle,
      "X step, Y step, Z step:        %.4lf, %.4lf, %.4lf\n",
      (double)((double)(od->grid.step[0]) / (double)(interpolate + 1)),
      (double)((double)(od->grid.step[1]) / (double)(interpolate + 1)),
      (double)((double)(od->grid.step[2]) / (double)(interpolate + 1)));
    fprintf(od->file[GRD_OUT]->handle,
      "X nodes, Y nodes, Z nodes:     %d, %d, %d\n",
      npts[0], npts[1], npts[2]);
    fprintf(od->file[GRD_OUT]->handle, "\n");
    break;
  }

  return 0;
}


void set_grid_point(O3Data *od, float *float_xy_mat, VarCoord *varcoord, double value)
{
  float_xy_mat[((od->grid.nodes[0]) + 2) * (int)(varcoord->node[1])
    + (int)(varcoord->node[0]) + 1] = (float)value;
}


double interpolated_value(double *ic_vec, double x_incr, double y_incr, double z_incr)
{
  int i;
  int j;
  int k;
  double x_incr_i;
  double y_incr_j;
  double z_incr_k;
  double value;
  
  
  value = 0.0;
  for (i = 0; i < 4; ++i) {
    for (j = 0; j < 4; ++j) {
      for (k = 0; k < 4; ++k) {
        x_incr_i = pow(x_incr, i);
        y_incr_j = pow(y_incr, j);
        z_incr_k = pow(z_incr, k);
        value += ic_vec[i + j * 4 + k * 16]
          * x_incr_i * y_incr_j * z_incr_k;
      }
    }
  }
  
  return value;
}


int grid_write(O3Data *od, char *filename, int pc_num, int type,
  int sign, int format, int label, int interpolate, int requested_endianness)
{
  char header[TITLE_LEN + METADATA_LEN + 16];
  char basename[BUF_LEN];
  char extension[TITLE_LEN];
  char write_mode[TITLE_LEN];
  int i;
  int j;
  int n;
  int x;
  int y;
  int z;
  int npts[3];
  int filename_len;
  int extension_len;
  int outer_loop_limit;
  int field_num;
  int x_var;
  int x_var_count;
  int z_plane;
  int actual_len = 0;
  int temp_len;
  int swap_endianness = 0;
  int result = 0;
  int n_loadings;
  int printed = 0;
  IntPerm *field_list = NULL;
  IntPerm *y_var_list = NULL;
  IntPerm *object_list = NULL;
  unsigned short state = 0;
  unsigned short bit = 0;
  double value;
  float on_value;
  float float_value;
  float *float_xy_mat = NULL;
  double double_value;
  FILE *grd_out = NULL;
  FILE *plane_file = NULL;
  FILE *x_matrix_handle = NULL;
  VarCoord varcoord;
  int group;
  int group_num = 0;
  IntPerm *group_list = NULL;


  for (i = 0; i < 3; ++i) {
    npts[i] = od->grid.nodes[i] + interpolate * (od->grid.nodes[i] - 1);
  }
  memset(basename, 0, BUF_LEN);
  memset(extension, 0, TITLE_LEN);
  switch (format) {
    case INSIGHT_FORMAT:
    strcpy(extension, INSIGHT_GRD_EXTENSION);
    strcpy(write_mode, "wb+");
    break;
    
    case MAESTRO_FORMAT:
    if (open_temp_file(od, od->file[TEMP_GRD], "maestro_plt")) {
      O3_ERROR_STRING(&(od->task), od->file[TEMP_GRD]->name);
      O3_ERROR_LOCATE(&(od->task));
      return CANNOT_WRITE_TEMP_FILE;
    }
    strcpy(extension, MAESTRO_GRD_EXTENSION);
    strcpy(write_mode, "wb+");
    break;
    
    case MOE_FORMAT:
    if (open_temp_file(od, od->file[TEMP_GRD], "moe_grd")) {
      O3_ERROR_STRING(&(od->task), od->file[TEMP_GRD]->name);
      O3_ERROR_LOCATE(&(od->task));
      return CANNOT_WRITE_TEMP_FILE;
    }
    strcpy(extension, MOE_GRD_EXTENSION);
    strcpy(write_mode, "wb+");
    break;
    
    case SYBYL_FORMAT:
    if (open_temp_file(od, od->file[TEMP_GRD], "sybyl_acnt")) {
      O3_ERROR_STRING(&(od->task), od->file[TEMP_GRD]->name);
      O3_ERROR_LOCATE(&(od->task));
      return CANNOT_WRITE_TEMP_FILE;
    }
    strcpy(extension, SYBYL_GRD_EXTENSION);
    strcpy(write_mode, "wb+");
    break;
    
    default:
    if (open_temp_file(od, od->file[TEMP_GRD], "ascii_grd")) {
      O3_ERROR_STRING(&(od->task), od->file[TEMP_GRD]->name);
      O3_ERROR_LOCATE(&(od->task));
      return CANNOT_WRITE_TEMP_FILE;
    }
    strcpy(extension, ASCII_GRD_EXTENSION);
    strcpy(write_mode, "wb+");
    break;
  }
  extension_len = strlen(extension);
  filename_len = strlen(filename);
  /*
  if the filename already includes the extension,
  then remove it
  */
  if (!strncmp(&filename[filename_len - extension_len],
    extension, extension_len)) {
    strncpy(basename, filename, filename_len - extension_len);
    basename[filename_len - extension_len] = '\0';
  }
  else {
    strncpy(basename, filename, BUF_LEN - 32);
    basename[BUF_LEN - 32] = '\0';
  }
  
  field_list = od->pel.numberlist[FIELD_LIST];
  y_var_list = od->pel.numberlist[Y_VAR_LIST];
  object_list = od->pel.numberlist[OBJECT_LIST];
  if ((type == WEIGHTS) || (type == LOADINGS)
    || (type == PCA_LOADINGS) || (type == COEFFICIENTS)) {
    /*
    Read from the TEMP_X_MATRIX file field_num
    and x_var numbers corresponding to active x_vars
    */
    x_matrix_handle = fopen(od->file[TEMP_X_MATRIX]->name, "rb");
    if (!x_matrix_handle) {
      O3_ERROR_STRING(&(od->task), od->file[TEMP_X_MATRIX]->name);
      O3_ERROR_LOCATE(&(od->task));
      return CANNOT_READ_TEMP_FILE;
    }
    od->file[TEMP_X_MATRIX]->handle = x_matrix_handle;
  }
  
  /*
  desume endianness
  of the machine we are running on
  */
  swap_endianness = (int)fabs(machine_type() - requested_endianness);
  outer_loop_limit = 1;
  switch (type) {
    case OBJECT_FIELD:
    outer_loop_limit = object_list->size;
    break;
    
    case COEFFICIENTS:
    outer_loop_limit = y_var_list->size;
    break;
    
    case GROUPS:
    field_list = od->pel.numberlist[FIELD_LIST];
    field_num = field_list->size;
    group_list = od->pel.numberlist[GROUP_LIST];
    group_num = group_list->size;
    if (!group_num) {
      group_num = od->mel.seed_count[field_list->pe[0] - 1];
      group_list = int_perm_resize(group_list, group_num);
      if (!group_list) {
        od->pel.numberlist[GROUP_LIST] = NULL;
        return OUT_OF_MEMORY;
      }
      od->pel.numberlist[GROUP_LIST] = group_list;
      for (i = 0; i < group_num; ++i) {
        group_list->pe[i] = i + 1;
      }
    }
    if (field_num > 1) {
      if ((group_num > 1)  ||
        ((group_num == 1) && (group_list->pe[0] != 0))) {
        return ONE_FIELD_AT_A_TIME;
      }
    }
    else {
      if ((group_list->pe[0] < 0) || (group_list->pe[group_num - 1]
        > od->mel.seed_count[field_list->pe[0] - 1])) {
        return INVALID_LIST_RANGE;
      }
    }
    break;
  }  
  /*
  allocate a temporary x*y float array
  */
  float_xy_mat = (float *)malloc(sizeof(float)
    * (od->grid.nodes[0] + 2) * od->grid.nodes[1]);
  if (!float_xy_mat) {
    return OUT_OF_MEMORY;
  }
  od->mel.float_xy_mat = float_xy_mat;
  memset(float_xy_mat, 0, sizeof(float)
    * (od->grid.nodes[0] + 2) * od->grid.nodes[1]);
  for (j = 0; j < outer_loop_limit; ++j) {
    for (i = 0; i < field_list->size; ++i) {
      if (!get_field_attr(od, field_list->pe[i] - 1, ACTIVE_BIT)) {
        continue;
      }
      if (od->file[TEMP_GRD]->handle) {
        rewind(od->file[TEMP_GRD]->handle);
      }
      if ((type == WEIGHTS) || (type == LOADINGS)
        || (type == PCA_LOADINGS) || (type == COEFFICIENTS)) {
        od->vel.b = double_vec_resize(od->vel.b,
          od->mal.x_weights->m);
        if (!(od->vel.b)) {
          return OUT_OF_MEMORY;
        }
        switch (type) {
          case WEIGHTS:
          if (reload_weights_loadings(od)) {
            O3_ERROR_STRING(&(od->task), od->file[TEMP_WLS]->name);
            O3_ERROR_LOCATE(&(od->task));
            return CANNOT_READ_TEMP_FILE;
          }
          memcpy(od->vel.b->ve,
            &M_PEEK(od->mal.x_weights, 0, pc_num - 1),
            od->mal.x_weights->m * sizeof(double));
          sprintf(od->file[GRD_OUT]->name, "%s%s%02d%s",
            basename, FIELD_EXTENSION, field_list->pe[i], extension);
          break;

          case LOADINGS:
          if (reload_weights_loadings(od)) {
            O3_ERROR_STRING(&(od->task), od->file[TEMP_WLS]->name);
            O3_ERROR_LOCATE(&(od->task));
            return CANNOT_READ_TEMP_FILE;
          }
          memcpy(od->vel.b->ve,
            &M_PEEK(od->mal.x_loadings, 0, pc_num - 1),
            od->mal.x_loadings->m * sizeof(double));
          sprintf(od->file[GRD_OUT]->name, "%s%s%02d%s",
            basename, FIELD_EXTENSION, field_list->pe[i], extension);
          break;

          case COEFFICIENTS:
          if (reload_coefficients(od, pc_num)) {
            O3_ERROR_STRING(&(od->task), od->file[TEMP_PLS_COEFF]->name);
            O3_ERROR_LOCATE(&(od->task));
            return CANNOT_READ_TEMP_FILE;
          }
          memcpy(od->vel.b->ve,
            &M_PEEK(od->mal.b_coefficients, 0,
            y_var_list->pe[j] - 1), od->mal.
            b_coefficients->m * sizeof(double));
          sprintf(od->file[GRD_OUT]->name, "%s%s%02d%s%02d%s",
            basename, FIELD_EXTENSION, field_list->pe[i],
            Y_VAR_EXTENSION, y_var_list->pe[j], extension);
          break;

          case PCA_LOADINGS:
          /*
          skip the pc number we are not interested in
          */
          if (fseek(od->file[TEMP_PCA_LOADINGS]->handle,
            sizeof(int), SEEK_SET)) {
            O3_ERROR_LOCATE(&(od->task));
            O3_ERROR_STRING(&(od->task), od->file[TEMP_PCA_LOADINGS]->name);
            return CANNOT_READ_TEMP_FILE;
          }
          /*
          read the number of loadings present in the file
          */
          actual_len = fread(&n_loadings, sizeof(int), 1,
            od->file[TEMP_PCA_LOADINGS]->handle);
          if (actual_len != 1) {
            O3_ERROR_LOCATE(&(od->task));
            O3_ERROR_STRING(&(od->task), od->file[TEMP_PCA_LOADINGS]->name);
            return CANNOT_READ_TEMP_FILE;
          }
          if (fseek(od->file[TEMP_PCA_LOADINGS]->handle,
            (pc_num - 1) * n_loadings * sizeof(double), SEEK_CUR)) {
            O3_ERROR_STRING(&(od->task), od->file[TEMP_PCA_LOADINGS]->name);
            O3_ERROR_LOCATE(&(od->task));
            return CANNOT_READ_TEMP_FILE;
          }
          actual_len = fread(od->vel.b->ve,
            sizeof(double), n_loadings,
            od->file[TEMP_PCA_LOADINGS]->handle);
          if (actual_len != n_loadings) {
            O3_ERROR_STRING(&(od->task), od->file[TEMP_PCA_LOADINGS]->name);
            O3_ERROR_LOCATE(&(od->task));
            return CANNOT_READ_TEMP_FILE;
          }
          sprintf(od->file[GRD_OUT]->name, "%s%s%02d%s",
            basename, FIELD_EXTENSION, field_list->pe[i], extension);
          break;
        }
        if (!(grd_out = fopen(od->file[GRD_OUT]->name, write_mode))) {
          O3_ERROR_STRING(&(od->task), od->file[GRD_OUT]->name);
          O3_ERROR_LOCATE(&(od->task));
          return PREMATURE_EOF;
        }
        od->file[GRD_OUT]->handle = grd_out;
        if (format == INSIGHT_FORMAT) {
          plane_file = grd_out;
          O3_ERROR_STRING(&(od->task), od->file[GRD_OUT]->name);
        }
        else {
          plane_file = od->file[TEMP_GRD]->handle;
          O3_ERROR_STRING(&(od->task), od->file[TEMP_GRD]->name);
        }
        /*
        write out the header
        */
        result = write_header(od, header, format, interpolate, swap_endianness);
        if (result) {
          O3_ERROR_LOCATE(&(od->task));
          return PREMATURE_EOF;
        }
        n = 0;
        x_var_count = 0;
        z_plane = 0;
        rewind(x_matrix_handle);
        actual_len = 1;
        while (actual_len) {
          fread(&field_num, sizeof(int), 1, x_matrix_handle);
          actual_len = fread(&x_var, sizeof(int), 1, x_matrix_handle);
          if (!actual_len) {
            field_num = field_list->pe[i] - 1;
            x_var = od->x_vars;
          }
          if (field_num != (field_list->pe[i] - 1)) {
            continue;
          }
          /*
          set grid points to zero until an active variable is found
          */
          while (x_var_count != x_var) {
            var_to_xyz(od, x_var_count, &varcoord);
            if ((int)(varcoord.node[2]) > z_plane) {
              result = write_grid_plane(od, plane_file,
                z_plane, interpolate, swap_endianness);
              if (result) {
                O3_ERROR_LOCATE(&(od->task));
                return PREMATURE_EOF;
              }
              ++z_plane;
            }
            set_grid_point(od, float_xy_mat, &varcoord, (double)0);
            ++x_var_count;
          }
          var_to_xyz(od, x_var_count, &varcoord);
          if ((int)(varcoord.node[2]) > z_plane) {
            result = write_grid_plane(od, plane_file,
              z_plane, interpolate, swap_endianness);
            if (result) {
              O3_ERROR_LOCATE(&(od->task));
              return PREMATURE_EOF;
            }
            ++z_plane;
          }
          if (x_var < od->x_vars) {
            value = od->vel.b->ve[n];
            if (((value > 0.0) && (sign == -1))
              || ((value < 0.0) && (sign == 1))) {
              value = 0.0;
            }
            set_grid_point(od, float_xy_mat, &varcoord, value);
          }
          if (!actual_len) {
            break;
          }
          ++n;
          ++x_var_count;
        }
        if ((format == SYBYL_FORMAT) || (format == ASCII_FORMAT) || (format == XYZ_FORMAT)) {
          /*
          x is the fastest varying coordinate,
          followed by y; z is the slowest
          */
          for (z = 0; z < npts[2]; ++z) {
            for (y = 0; y < npts[1]; ++y) {
              for (x = 0; x < npts[0]; ++x) {
                if (fseek(od->file[TEMP_GRD]->handle,
                  ((npts[0] + 2) * (npts[1] * z + y) + x + 1)
                  * sizeof(float), SEEK_SET)) {
                  O3_ERROR_STRING(&(od->task), od->file[TEMP_GRD]->name);
                  O3_ERROR_LOCATE(&(od->task));
                  return CANNOT_READ_TEMP_FILE;
                }
                temp_len = fread(&float_value, sizeof(float), 1,
                   od->file[TEMP_GRD]->handle);
                if (temp_len != 1) {
                  O3_ERROR_STRING(&(od->task), od->file[TEMP_GRD]->name);
                  O3_ERROR_LOCATE(&(od->task));
                  return CANNOT_READ_TEMP_FILE;
                }
                if (format == XYZ_FORMAT) {
                  fprintf(grd_out, "%16.4lf%16.4lf%16.4lf",
                    safe_rint((double)x
                    * (double)(od->grid.step[0]) / (double)(interpolate + 1)
                    * 1.0e04) / 1.0e04 + (double)(od->grid.start_coord[0]),
                    safe_rint((double)y
                    * (double)(od->grid.step[1]) / (double)(interpolate + 1)
                    * 1.0e04) / 1.0e04 + (double)(od->grid.start_coord[1]),
                    safe_rint((double)z
                    * (double)(od->grid.step[2]) / (double)(interpolate + 1)
                    * 1.0e04) / 1.0e04 + (double)(od->grid.start_coord[2]));
                }
                fprintf(grd_out, "%16.6e\n", float_value);
              }
            }
          }
        }
        else if ((format == MAESTRO_FORMAT) || (format == MOE_FORMAT)) {
          /*
          z is the fastest varying coordinate,
          followed by y; x is the slowest
          */
          for (x = 0; x < npts[0]; ++x) {
            for (y = 0; y < npts[1]; ++y) {
              for (z = 0; z < npts[2]; ++z) {
                if (fseek(od->file[TEMP_GRD]->handle,
                  ((npts[0] + 2) * (npts[1] * z + y) + x + 1)
                  * sizeof(float), SEEK_SET)) {
                  O3_ERROR_STRING(&(od->task), od->file[TEMP_GRD]->name);
                  O3_ERROR_LOCATE(&(od->task));
                  return CANNOT_READ_TEMP_FILE;
                }
                temp_len = fread(&float_value, sizeof(float), 1,
                   od->file[TEMP_GRD]->handle);
                if (temp_len != 1) {
                  O3_ERROR_STRING(&(od->task), od->file[TEMP_GRD]->name);
                  O3_ERROR_LOCATE(&(od->task));
                  return CANNOT_READ_TEMP_FILE;
                }
                if (format == MAESTRO_FORMAT) {
                  fprintf(grd_out, "%16.6e\n", float_value);
                }
                else {
                  double_value = (double)float_value;
                  temp_len = fwrite(&double_value, sizeof(double), 1,  grd_out);
                  if (temp_len != 1) {
                    O3_ERROR_STRING(&(od->task), od->file[TEMP_GRD]->name);
                    O3_ERROR_LOCATE(&(od->task));
                    return CANNOT_WRITE_TEMP_FILE;
                  }
                }
              }
            }
          }
        }
        if (format == MOE_FORMAT) {
          n = 0;
          temp_len = fwrite(&n, sizeof(int), 1,  grd_out);
          if (temp_len != 1) {
            O3_ERROR_STRING(&(od->task), od->file[GRD_OUT]->name);
            O3_ERROR_LOCATE(&(od->task));
            return CANNOT_WRITE_TEMP_FILE;
          }
        }
        if (grd_out) {
          fclose(grd_out);
          od->file[GRD_OUT]->handle = NULL;
        }
      }
      else {
        sprintf(od->file[GRD_OUT]->name, "%s%s%02d%s",
          basename, FIELD_EXTENSION, field_list->pe[i], extension);
        on_value = (float)1;
        switch (type) {
          case TWO_LEVEL:
          bit = TWO_LEVEL_BIT;
          state = 1;
          on_value = (float)2;
          break;
          
          case THREE_LEVEL:
          bit = THREE_LEVEL_BIT;
          state = 1;
          on_value = (float)3;
          break;
          
          case FOUR_LEVEL:
          bit = FOUR_LEVEL_BIT;
          state = 1;
          on_value = (float)4;
          break;
          
          case OBJECT_FIELD:
          sprintf(od->file[GRD_OUT]->name, "%s%s%02d%s%02d%s",
            basename, FIELD_EXTENSION, field_list->pe[i],
            OBJECT_EXTENSION, object_list->pe[j], extension);
          break;
          
          case D_OPTIMAL:
          bit = D_OPTIMAL_BIT;
          state = 0;
          break;
          
          case FFDSEL:
          bit = FFDSEL_BIT;
          state = 0;
          break;
          
          case UVEPLS:
          bit = UVEPLS_BIT;
          state = 0;
          break;
          
          case GROUPS:
          bit = GROUP_BIT | SEED_BIT;
          state = 1;
          break;
        }
        if (!(grd_out = fopen(od->file[GRD_OUT]->name, write_mode))) {
          O3_ERROR_STRING(&(od->task), od->file[GRD_OUT]->name);
          O3_ERROR_LOCATE(&(od->task));
          return PREMATURE_EOF;
        }
        od->file[GRD_OUT]->handle = grd_out;
        if (format == INSIGHT_FORMAT) {
          plane_file = grd_out;
        }
        else {
          plane_file = od->file[TEMP_GRD]->handle;
        }
        /*
        write out the header
        */
        result = write_header(od, header, format, interpolate, swap_endianness);
        if (result) {
          O3_ERROR_LOCATE(&(od->task));
          return PREMATURE_EOF;
        }
        z_plane = 0;
        for (x_var_count = 0; x_var_count < od->x_vars; ++x_var_count) {
          var_to_xyz(od, x_var_count, &varcoord);
          if (type == SD) {
            value = get_x_var_buf(od, field_list->pe[i] - 1,
              x_var_count, AVE_BUF);
          }
          else if (type == OBJECT_FIELD) {
            result = get_x_value(od, field_list->pe[i] - 1,
              object_list->pe[j] - 1, x_var_count, &value, 0);
            if (result) {
              return result;
            }
            /*
            if the ACTIVE_BIT is set, then
            values corresponding to inactive variables
            (i.e., below the stddev cutoff) are written
            as zeros, otherwise the original values are
            used
            */
            if (MISSING(value)) {
              if (!(label & LABEL_MISSING_BIT)) {
                value = 0.0;
              }
            }
            else {
              result = get_x_value(od, field_list->pe[i] - 1,
                object_list->pe[j] - 1, x_var_count, &value,
                CUTOFF_BIT);
              if (result) {
                return result;
              }
              if ((!get_x_var_attr(od, field_list->pe[i] - 1, x_var_count, ACTIVE_BIT))
                && (label & (LABEL_INACTIVE_BIT | ZERO_INACTIVE_BIT))) {
                value = ((label & LABEL_INACTIVE_BIT)
                  ? INACTIVE_VALUE : 0.0);
              }
            }
          }
          else {
            if (get_x_var_attr(od, field_list->pe[i] - 1, x_var_count, ACTIVE_BIT)
              && ((get_x_var_attr(od, field_list->pe[i] - 1,
              x_var_count, bit) ? 1 : 0) == state)) {
              if (type == GROUPS) {
                group = get_voronoi_buf(od, i, x_var_count);
                on_value = (group ? (float)group : (float)-1);
              }
              else {
                value = (double)on_value;
              }
            }
            else {
              value = (double)0;
            }
          }
          if (((int)(varcoord.node[2]) > z_plane)
            || (x_var_count == (od->x_vars - 1))) {
            if (x_var_count == (od->x_vars - 1)) {
              set_grid_point(od, float_xy_mat, &varcoord, value);
            }
            result = write_grid_plane(od, plane_file,
              z_plane, interpolate, swap_endianness);
            if (result) {
              O3_ERROR_LOCATE(&(od->task));
              return PREMATURE_EOF;
            }
            ++z_plane;
          }
          set_grid_point(od, float_xy_mat, &varcoord, value);
        }
        if ((format == SYBYL_FORMAT) || (format == ASCII_FORMAT) || (format == XYZ_FORMAT)) {
          for (z = 0; z < npts[2]; ++z) {
            for (y = 0; y < npts[1]; ++y) {
              for (x = 0; x < npts[0]; ++x) {
                if (fseek(od->file[TEMP_GRD]->handle,
                  ((npts[0] + 2) * (npts[1] * z + y) + x + 1)
                  * sizeof(float), SEEK_SET)) {
                  O3_ERROR_LOCATE(&(od->task));
                  return CANNOT_READ_TEMP_FILE;
                }
                temp_len = fread(&float_value, sizeof(float), 1,
                   od->file[TEMP_GRD]->handle);
                if (temp_len != 1) {
                  O3_ERROR_LOCATE(&(od->task));
                  return CANNOT_READ_TEMP_FILE;
                }
                if (format == XYZ_FORMAT) {
                  fprintf(grd_out, "%16.4lf%16.4lf%16.4lf",
                    safe_rint((double)x
                    * (double)(od->grid.step[0]) / (double)(interpolate + 1)
                    * 1.0e04) / 1.0e04 + (double)(od->grid.start_coord[0]),
                    safe_rint((double)y
                    * (double)(od->grid.step[1]) / (double)(interpolate + 1)
                    * 1.0e04) / 1.0e04 + (double)(od->grid.start_coord[1]),
                    safe_rint((double)z
                    * (double)(od->grid.step[2]) / (double)(interpolate + 1)
                    * 1.0e04) / 1.0e04 + (double)(od->grid.start_coord[2]));
                }
                printed = 0;
                if (label & (LABEL_MISSING_BIT | LABEL_INACTIVE_BIT)) {
                  if (MISSING(float_value)) {
                    fprintf(grd_out, "%16s\n", COMFA_MISSING_VALUE);
                    printed = 1;
                  }
                  else if (INACTIVE(float_value)) {
                    fprintf(grd_out, "%16s\n", "INACTIVE");
                    printed = 1;
                  }
                }
                if (!printed) {
                  fprintf(grd_out, "%16.6e\n", float_value);
                }
              }
            }
          }
        }
        else if ((format == MAESTRO_FORMAT) || (format == MOE_FORMAT)) {
          for (x = 0; x < npts[0]; ++x) {
            for (y = 0; y < npts[1]; ++y) {
              for (z = 0; z < npts[2]; ++z) {
                if (fseek(od->file[TEMP_GRD]->handle,
                  ((npts[0] + 2) * (npts[1] * z + y) + x + 1)
                  * sizeof(float), SEEK_SET)) {
                  O3_ERROR_LOCATE(&(od->task));
                  return CANNOT_READ_TEMP_FILE;
                }
                temp_len = fread(&float_value, sizeof(float), 1,
                   od->file[TEMP_GRD]->handle);
                if (temp_len != 1) {
                  O3_ERROR_LOCATE(&(od->task));
                  return CANNOT_READ_TEMP_FILE;
                }
                if (format == MAESTRO_FORMAT) {
                  fprintf(grd_out, "%16.6e\n", float_value);
                }
                else {
                  double_value = (double)float_value;
                  temp_len = fwrite(&double_value, sizeof(double), 1, grd_out);
                  if (temp_len != 1) {
                    O3_ERROR_LOCATE(&(od->task));
                    return CANNOT_WRITE_TEMP_FILE;
                  }
                }
              }
            }
          }
        }
        if (format == MOE_FORMAT) {
          n = 0;
          temp_len = fwrite(&n, sizeof(int), 1,  grd_out);
          if (temp_len != 1) {
            O3_ERROR_LOCATE(&(od->task));
            return CANNOT_WRITE_TEMP_FILE;
          }
        }
        if (grd_out) {
          fclose(grd_out);
          od->file[GRD_OUT]->handle = NULL;
        }
      }
    }
  }
  if (od->file[TEMP_GRD]->handle) {
    fclose(od->file[TEMP_GRD]->handle);
    od->file[TEMP_GRD]->handle = NULL;
  }
  if (od->file[TEMP_PCA_LOADINGS]->handle) {
    fclose(od->file[TEMP_PCA_LOADINGS]->handle);
    od->file[TEMP_PCA_LOADINGS]->handle = NULL;
  }
  if (x_matrix_handle) {
    fclose(x_matrix_handle);
    od->file[TEMP_X_MATRIX]->handle = NULL;
  }
  if (float_xy_mat) {
    free(float_xy_mat);
    od->mel.float_xy_mat = NULL;
  }
  if (od->vel.b) {
    double_vec_free(od->vel.b);
    od->vel.b = NULL;
  }

  return 0;
}
