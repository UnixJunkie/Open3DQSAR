/*

o3header.h

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


#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <stdio.h>
#include <ctype.h>
#include <errno.h>
#include <fcntl.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <unistd.h>
#include <dirent.h>
#include <time.h>
#include <signal.h>
#include <stddef.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/stat.h>
#ifndef WIN32
#include <sys/wait.h>
#include <sys/resource.h>
#include <sys/utsname.h>
#include <sys/mman.h>
#include <pthread.h>
#include <termios.h>
#include <fnmatch.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#else
#include <windef.h>
#include <winsock2.h>
#ifdef _TIMESPEC_DEFINED
#define HAVE_STRUCT_TIMESPEC 1
#endif
#endif
#include <math.h>
#ifdef O3Q
#include "cdflib.h"
#endif
#ifdef HAVE_CLAPACK_H
#ifndef HAVE_MKL_LAPACK_H
#include <clapack.h>
#endif
#endif
#ifdef HAVE_LAPACKE_H
#include <lapacke.h>
#endif
#ifdef HAVE_CBLAS_H
#if ((!defined(HAVE_MKL_CBLAS_H)) && (!defined(HAVE_SUNPERF_H)))
#include <cblas.h>
#endif
#endif
#ifdef HAVE_MKL_H
#include <mkl.h>
#endif
#ifdef HAVE_MKL_CBLAS_H
#include <mkl_cblas.h>
#endif
#ifdef HAVE_MKL_LAPACK_H
#include <mkl_lapack.h>
#endif
#ifdef HAVE_SUNPERF_H
#include <sunperf.h>
#endif
#include <zlib.h>
#ifdef HAVE_MINIZIP_ZIP_H
#include <minizip/zip.h>
#endif
#ifdef HAVE_MINIZIP_UNZIP_H
#include <minizip/unzip.h>
#endif
#include <include/safe_rint.h>

#ifndef PRIO_MIN
#define PRIO_MIN      -20
#endif
#ifndef PRIO_MAX
#define PRIO_MAX      20
#endif
#define O3Q_LITTLE_ENDIAN 0
#define O3Q_BIG_ENDIAN 1
#ifndef INFINITY
#define INFINITY      HUGE_VAL
#endif
#define O3_ERROR_LOCATE(task)    (task)->line = __LINE__; \
          strcpy((task)->file, __FILE__); \
          strcpy((task)->func, __PRETTY_FUNCTION__)
#define O3_ERROR_STRING(task, x)  strcpy((task)->string, x)
#define O3_ERROR_PRINT(task)    if (od->debug) print_debug_info(od, task)
#define IS_O3A(x)      (x->package_code[2] == 'A')
#define IS_O3G(x)      (x->package_code[2] == 'G')
#define IS_O3Q(x)      (x->package_code[2] == 'Q')
#define square(x)      ((x) * (x))
#define absval(x)      ((x) < 0 ? (-(x)) : (x))
#define angle2rad(x)      (M_PI / 180.0 * (double)(x))
#define rad2angle(x)      ((double)(x) / M_PI * 180.0)
#define double_mat_alloc(m, n)    double_mat_resize(NULL, m, n)
#define double_vec_alloc(size)    double_vec_resize(NULL, size)
#define int_perm_alloc(size)    int_perm_resize(NULL, size)
#define set_x_value_xyz(od, field_num, object_num, varcoord, value) \
  set_x_value(od, field_num, object_num, xyz_to_var(od, varcoord), value)
#define set_x_value_xyz_unbuffered(od, field_num, object_num, varcoord, value) \
  set_x_value_unbuffered(od, field_num, object_num, xyz_to_var(od, varcoord), value)
#define M_POKE(double_mat, m, n, value)  (double_mat)->base[(m) + (n) * ((double_mat)->max_m)] = value
#define M_PEEK(double_mat, m, n)  ((double_mat)->base[(m) + (n) * ((double_mat)->max_m)])
#define MISSING(x)      (x > 1.0e36)
#define INACTIVE(x)      (x < -1.0e36)
#ifdef HAVE_LIBMKL
#define LWORK_BLOCK_SIZE    64
#endif
#ifdef HAVE_LIBACCELERATE
#define LWORK_BLOCK_SIZE    64
#endif
#define MAX_STACK      1024

#define BOHR_RADIUS      0.52917715
#define EV_TO_KCAL      23.06035
#define HARTREE_TO_KCAL      627.5095
#define R_KCAL_K_MOL      1.9858775e-03
#define MERSENNE_N      624
#define MERSENNE_M      397
#define DUMMY_COST      100000
#define DUMMY_Y_VAR_VALUE      1.0
#define MAX_H_BINS      20
#define SYMMETRY_H_BINS      5
#define MATRIX_A      0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK      0x80000000UL /* most significant w-r bits */
#define LOWER_MASK      0x7fffffffUL /* least significant r bits */
#define BUF_LEN        1024
#define LARGE_BUF_LEN      8192
#define FZ_BUF_LEN      65536
#define TITLE_LEN      60
#define METADATA_LEN      72
#define SHORT_PROMPT      "> "
#define LOCALHOST_IP      "127.0.0.1"
#define JMOL_FIRST_PORT      49152
#define JMOL_LAST_PORT      65535
#define MOLDEN_INP_EXT      ".mdninp"
#define GAMESS_PUNCH_EXT    ".dat"
#define GAUSSIAN_CUBE_EXT    ".cube"
#define TURBOMOLE_COSMO_EXT    ".cosmo"
#define GRUB_FILENAME      "grub.dat"
#define BABEL_PATH_ENV      "O3_BABEL_PATH"
#define BABEL_DATADIR_ENV    "BABEL_DATADIR"
#define BABEL_LIBDIR_ENV    "BABEL_LIBDIR"
#define TEMP_DIR_ENV      "O3_TEMP_DIR"
#define DAT_HEADER      "HEADER"
#define SDF_DELIMITER      "$$$$"
#define MOL_DELIMITER      "M  END"
#define MD_GRID_PDB_FILE_HEADER    "# PDB FILE IN MD GRID FORMAT GENERATED BY "PACKAGE_NAME_UPPERCASE
#define INSIGHT_GRID_TITLE    "INSIGHT GRID GENERATED BY "PACKAGE_NAME_UPPERCASE
#define INSIGHT_GRD_EXTENSION    ".grd"
#define MOE_GRID_TITLE      "MOE GRID GENERATED BY "PACKAGE_NAME_UPPERCASE
#define MOE_ID_TOKEN      "#grid"
#define MOE_FORMAT_VERSION    2008.03
#define MOE_GRD_EXTENSION    ".mgd"
#define COMFA_MISSING_VALUE    "MISSING"
#define SYBYL_GRID_TITLE    "SYBYL ASCII CONTOUR GENERATED BY "PACKAGE_NAME_UPPERCASE
#define SYBYL_GRD_EXTENSION    ".acnt"
#define CUBE_GRID_TITLE    "CUBE GRID GENERATED BY "PACKAGE_NAME_UPPERCASE
#define CUBE_GRD_EXTENSION    GAUSSIAN_CUBE_EXT
#define MAESTRO_GRD_EXTENSION    ".plt"
#define OPENDX_GRD_EXTENSION    ".dx"
#define ASCII_GRD_EXTENSION    ".agrd"
#define FIELD_EXTENSION      "_fld-"
#define OBJECT_EXTENSION    "_obj-"
#define Y_VAR_EXTENSION      "_y-"
#define PLT_DAT_EXTENSION    ".txt"
#define PLT_CMD_EXTENSION    ".gnuplot"
#define PYMOL_TRAINING_SET_CARBON    "gray80"
#define PYMOL_TEST_SET_CARBON      "brightorange"
#define JMOL_TRAINING_SET_CARBON    "Jmol"
#define JMOL_TEST_SET_CARBON      "orange"
#define FIREFLY_OPTIONS      "-stdext -r -f -p"
#define FIREFLY_OPTIONS      "-stdext -r -f -p"
#define GAMESS_NORMAL_TERMINATION  "GAMESS TERMINATED NORMALLY"
#define FIREFLY_NORMAL_TERMINATION  "FIREFLY TERMINATED NORMALLY"
#define GAUSSIAN_NORMAL_TERMINATION  "Normal termination"
#define TURBOMOLE_NORMAL_TERMINATION  "all done"
#define TINKER_MINIMIZE_TERMINATION  "Normal Termination"
#define TINKER_DYNAMIC_TERMINATION  "Instantaneous Values for Frame saved"
#define TINKER_MMFF94_PRM_FILE    "mmff.prm"
#define TINKER_MMFF94S_PRM_FILE    "mmffs.prm"
#define PHARAO_NO_HYBRID    "--noHybrid"
#define PHARAO_MERGE      "-m"
#define HEADER_FOUND      1
#define EOF_FOUND      2
#define PREMATURE_EOF      10
#define PREMATURE_DAT_EOF    14
#define PREMATURE_DEP_IN_EOF    15
#define PREMATURE_PLT_DAT_EOF    16
#define PREMATURE_PLT_CMD_EOF    17
#define OUT_OF_MEMORY      50
#define PARSE_INPUT_ERROR    60
#define PARSE_INPUT_RECOVERABLE_ERROR  61
#define NOT_ENOUGH_OBJECTS    110
#define TOO_MANY_OBJECTS    111
#define NOT_ENOUGH_Y_VARS    120
#define CANNOT_FIND_Y_VAR_NAME    121
#define WRONG_NUMBER_OF_Y_VARS    130
#define WRONG_NUMBER_OF_X_VARS    140
#define WRONG_NUMBER_OF_FIELDS    141
#define INVALID_LIST_RANGE    150
#define PARSE_NLEVEL_TYPE_ERROR    160
#define CANNOT_WRITE_TEMP_FILE    170
#define CANNOT_READ_TEMP_FILE    180
#define CANNOT_READ_OUT_FILE    181
#define CANNOT_READ_LOG_FILE    182
#define Y_VAR_LOW_SD      190
#define STATUS_INCONSISTENCY    230
#define ONE_FIELD_AT_A_TIME    270
#define STUDENT_T_OUT_OF_BOUNDS    280
#define SINGULAR_MATRIX      290
#define CANNOT_CREATE_THREAD    300
#define CANNOT_JOIN_THREAD    301
#define GRID_NOT_MATCHING    320
#define GRID_NOT_MATCHING_OFF_CENTER  321
#define GRID_NOT_MATCHING_OUT_OF_BOUNDS  322
#define CANNOT_READ_GRID_DATA    323
#define CANNOT_FIND_MO      324
#define OBJECTS_NOT_MATCHING    330
#define TIME_NOT_AVAILABLE    350
#define CANNOT_CHANGE_DIR    370
#define QM_ABNORMAL_TERMINATION    400
#define CANNOT_WRITE_QM_INP_FILE  401
#define CANNOT_READ_QM_FCHK_FILE  404
#define CANNOT_READ_QM_GCUBE_FILE  405
#define OPENBABEL_ERROR      410
#define CANNOT_COPY_GRUB_DAT    420
#define GRID_ERROR      430
#define POSSIBLE_GRID_BUG    431
#define ERROR_EXTRACTING_PHAR    440
#define ERROR_IN_ALIGNMENT    441
#define ERROR_MERGING_FILES    442
#define CANNOT_READ_ORIGINAL_SDF  450
#define CANNOT_WRITE_ALIGNED_SDF  460
#define CANNOT_WRITE_ROTOTRANSED_SDF  461
#define N_ATOM_BOND_MISMATCH    470
#define BABEL_PLUGINS_NOT_FOUND    480
#define BABEL_NOT_WORKING    481
#define BABEL_TOO_OLD      482
#define ERROR_IN_FILTER_EXTRACT_SPLIT  490
#define ERROR_IN_FILTER_INTER    491
#define NOTHING_TO_DO_FILTER    492
#define ERROR_IN_FILTER_INTRA    493
#define NOT_ENOUGH_PARAMETERS    500
#define NOT_ENOUGH_VALUES    510
#define INCORRECT_NUMBER_OF_VALUES  511
#define WRONG_DATA_FORMAT    520
#define WRONG_PARAMETER_NAME    530
#define DUPLICATE_PARAMETER_NAME  531
#define CANNOT_OPEN_DIRECTORY    550
#define CANNOT_CREATE_DIRECTORY    560
#define CS3D_ERROR      570
#define CANNOT_SEND_JMOL_COMMAND      580
#define OUT_FILE_NOT_EMPTY      590
#define LOG_FILE_NOT_EMPTY      591
#define BAD_DX_HEADER      600
#define FL_CANNOT_CREATE_CHANNELS  (1<<0)
#define FL_CANNOT_CREATE_PROCESS  (1<<1)
#define FL_CANNOT_READ_OUT_FILE    (1<<2)
#define FL_CANNOT_READ_SDF_FILE    (1<<2)
#define FL_CANNOT_READ_CONF_FILE  (1<<3)
#define FL_CANNOT_READ_FCHK_FILE  (1<<4)
#define FL_CANNOT_FIND_CONF    (1<<4)
#define FL_CANNOT_READ_INP_FILE    (1<<5)
#define FL_CANNOT_WRITE_INP_FILE  (1<<6)
#define FL_CANNOT_WRITE_SDF_FILE  (1<<6)
#define FL_CANNOT_READ_GCUBE_FILE  (1<<7)
#define FL_ABNORMAL_TERMINATION    (1<<8)
#define FL_OUT_OF_MEMORY    (1<<9)
#define FL_CANNOT_CREATE_SCRDIR    (1<<10)
#define FL_CANNOT_CREATE_MDNDIR    (1<<10)
#define FL_CANNOT_CHDIR      (1<<11)
#define FL_CANNOT_READ_OB_OUTPUT  (1<<12)
#define FL_CANNOT_READ_PHARAO_OUTPUT  (1<<12)
#define FL_PHARAO_ERROR      (1<<13)
#define FL_GRIN_ERROR      (1<<14)
#define FL_ATOM_TYPE_MISMATCH    (1<<14)
#define FL_CANNOT_READ_MOL_FILE    (1<<15)
#define FL_CANNOT_READ_TEMP_FILE  (1<<16)
#define FL_CANNOT_WRITE_TEMP_FILE  (1<<17)
#define FL_UNKNOWN_ATOM_TYPE    (1<<18)
#define MAX_ARG        32
#define MAX_FREE_FORMAT_PARAMETERS  5
#define MAX_DATA_FIELDS      4
#define DATA_OBJECT_NUM      0
#define DATA_BEST_OBJECT_NUM    1
#define DATA_N_CONF      0
#define DATA_N_CONF_OVERALL    1
#define DELETED_CONF      1
#define DELETE_CANDIDATE_CONF    2
#define TEMPLATE_DB      0
#define CANDIDATE_DB      1
#define ANY_DB        2
#define MAX_DB        3
#define BOUND        0
#define WATER        1
#define MAX_STATES      2
#define TEMPLATE_OBJECT_NUM    0
#define TEMPLATE_CONF_NUM    1
#define MOVED_OBJECT_NUM    2
#define MOVED_CONF_NUM      3
#define MAX_ATTEMPTS_FILE    10
#define MAX_ATTEMPTS_JMOL    100
#define O3_MAX_SLOT      10
#define MS_SLEEP_BEFORE_RETRY    100
#define CLEAN_UP_SUFFIX      ".CLEAN_UP_DONT_DELETE_ME"
#define AROMATIC      4
#define AROMATIC_TEMP      8
#define RING_BIT      2
#define SP2_NON_AROMATIC    5
#define SP2_ANY        6
#define INTERACTIVE_RUN      1
#define DRY_RUN        2
#define FULL_MODEL      1
#define CV_MODEL      2
#define LEAVE_ONE_OUT      2
#define LEAVE_TWO_OUT      4
#define LEAVE_MANY_OUT      8
#define PARALLEL_CV      4
#define FFDSEL_FULL_MODEL    8
#define FFDSEL_CV_MODEL      16
#define UVEPLS_FULL_MODEL    32
#define UVEPLS_CV_MODEL      64
#define SCRAMBLE_CV_MODEL    128
#define EXTERNAL_PREDICTION    256
#define SILENT_PLS      512
#define MAX_FILES      36
#define BINARY_IN      0
#define ASCII_IN      1
#define MOLFILE_IN      2
#define DAT_IN        3
#define DAT_OUT        4
#define GRD_OUT        5
#define DEP_IN        6
#define PLT_DAT_OUT      7
#define PLT_CMD_OUT      8
#define PREPINP_OUT      9
#define MAIN_INPUT      10
#define MAIN_OUTPUT      11
#define TEMP_START      12
#define TEMP_OUT      TEMP_START
#define TEMP_LOG      TEMP_START + 1
#define TEMP_WLS      TEMP_START + 2
#define TEMP_PLS_COEFF      TEMP_START + 3
#define TEMP_OBJECT_MATCH    TEMP_START + 4
#define TEMP_PCA_SCORES      TEMP_START + 4
#define TEMP_SORTED_MATCH    TEMP_START + 5
#define TEMP_PCA_LOADINGS    TEMP_START + 5
#define TEMP_CALC      TEMP_START + 6
#define TEMP_PRED      TEMP_START + 7
#define TEMP_EXT_PRED      TEMP_START + 8
#define TEMP_X_MATRIX      TEMP_START + 9
#define TEMP_VORONOI      TEMP_START + 10
#define TEMP_DAT      TEMP_START + 11
#define TEMP_SCRAMBLE      TEMP_START + 12
#define TEMP_CV_COEFF      TEMP_START + 13
#define TEMP_UVEPLS      TEMP_START + 14
#define TEMP_GRD      TEMP_START + 15
#define TEMP_MOLFILE      TEMP_START + 16
#define TEMP_PYMOL      TEMP_START + 17
#define TEMP_BABEL      TEMP_START + 18
#define TEMP_FIELD_DATA      MAX_FILES
#ifdef WIN32
#define NORMAL_INK      BACKGROUND_BLUE | BACKGROUND_GREEN | BACKGROUND_RED | BACKGROUND_INTENSITY | FOREGROUND_BLUE
#ifdef O3Q
#define HIGHLIGHT_INK      BACKGROUND_BLUE | BACKGROUND_GREEN | BACKGROUND_RED | BACKGROUND_INTENSITY | FOREGROUND_RED
#elif O3G
#define HIGHLIGHT_INK      BACKGROUND_BLUE | BACKGROUND_GREEN | BACKGROUND_RED | BACKGROUND_INTENSITY | FOREGROUND_RED | FOREGROUND_GREEN
#else
#define HIGHLIGHT_INK      BACKGROUND_BLUE | BACKGROUND_GREEN | BACKGROUND_RED | BACKGROUND_INTENSITY | FOREGROUND_GREEN
#endif
#define SET_INK(x, y)      if ((x)->hOutput && (!((x)->terminal))) { \
            SetConsoleTextAttribute((x)->hOutput, (DWORD)y); \
            tee_printf(x, "\n"); \
          }
#define GET_TEMPDIR      getenv("TEMP")
#define NEWLINE      "\r\n"
#define SEPARATOR      '\\'
#define PATH_SEPARATOR      ";"
#define APPDATA        "APPDATA"
#define HOMEDIR        PACKAGE_NAME
#define HISTORY_FILE      ".history"
#define DEFAULT_SHELL_CMD      "cmd.exe /c"
#define EXE_PATH      "PATH"
#define BABEL_EXE      "babel.exe"
#define TINKER_ANALYZE_EXE    "analyze.exe"
#define TINKER_MINIMIZE_EXE    "minimize.exe"
#define TINKER_OPTIMIZE_EXE    "optimize.exe"
#define TINKER_DYNAMIC_EXE    "dynamic.exe"
#define GRID_EXE      "grid.exe"
#define GRIN_EXE      "grin.exe"
#define OBENERGY_EXE      "obenergy.exe"
#define FORMCHK_EXE      "formchk.exe"
#define FORMCHK_EXT      ".fck"
#define CUBEGEN_EXE      "cubegen.exe"
#define FIREFLY_EXE      "firefly.exe"
#define G03_EXE        "g03.exe"
#define G09_EXE        "g09.exe"
#define GAMESS_DDIKICK_EXE    "ddikick.exe"
#define GAMESS_MPIEXEC_EXE    "mpiexec.exe"
#define GAMESS_SMPD_EXE      "smpd.exe"
#define TURBOMOLE_DEFINE_EXE    "define.exe"
#define TURBOMOLE_DSCF_EXE    "dscf.exe"
#define TURBOMOLE_RIDFT_EXE    "ridft.exe"
#define CS3D_EXE      "cs3d.exe"
#define PHARAO_EXE      "align-it.exe"
#define GNUPLOT_EXE      "wgnuplot.exe"
#define GNUPLOT_ARGS      "-"
#define PYMOL_EXE      "pymol.bat"
#define PYMOL_ARGS      "-p -q"
#define JMOL_EXE      "jmol.bat"
#define JMOL_ARGS      "-Lo -g800x600"
#define FWRITE_WRAP(h, b, n)    WriteFile(h, b, (DWORD)strlen(b), n, NULL)
#define FFLUSH_WRAP(h)      /**/
#define FCLOSE_WRAP(h)      CloseHandle(h)
#define ISSEPARATOR(c)      ((c == '/') || (c == '\\'))
#else
#ifdef __APPLE__
#define DYN_LIBRARY_PATH    "DYLD_LIBRARY_PATH"
#else
#define DYN_LIBRARY_PATH    "LD_LIBRARY_PATH"
#endif
#ifdef __FreeBSD__
#define DYN_32_LIBRARY_PATH    "LD_32_LIBRARY_PATH"
#endif
#define NORMAL_INK      "\x1b[0;34m"
#ifdef O3Q
#define HIGHLIGHT_INK      "\x1b[0;31m"
#elif O3G
#define HIGHLIGHT_INK      "\x1b[0;33m"
#else
#define HIGHLIGHT_INK      "\x1b[0;32m"
#endif
#define DEFAULT_INK      "\x1b[m"
#define SET_INK(x, y)      if ((!((x)->terminal)) && ((x)->prompt)) { \
            printf("%s\n", y); \
          }
#define GET_TEMPDIR      P_tmpdir
#define NEWLINE      "\n"
#define SEPARATOR      '/'
#define PATH_SEPARATOR      ":"
#define APPDATA        "HOME"
#define HOMEDIR        "."PACKAGE_NAME_LOWERCASE
#define HISTORY_FILE      ".history"
#define DEFAULT_SHELL_CMD      "/usr/bin/env bash -c"
#define EXE_PATH      "PATH"
#define BABEL_EXE      "babel"
#define TINKER_ANALYZE_EXE    "analyze"
#define TINKER_MINIMIZE_EXE    "minimize"
#define TINKER_OPTIMIZE_EXE    "optimize"
#define TINKER_DYNAMIC_EXE    "dynamic"
#define GRID_EXE      "grid"
#define GRIN_EXE      "grin"
#define OBENERGY_EXE      "obenergy"
#define FORMCHK_EXE      "formchk"
#define FORMCHK_EXT      ".fchk"
#define CUBEGEN_EXE      "cubegen"
#define FIREFLY_EXE      "firefly"
#define G03_EXE        "g03"
#define G09_EXE        "g09"
#define GAMESS_DDIKICK_EXE    "ddikick.x"
#define TURBOMOLE_DEFINE_EXE    "define"
#define TURBOMOLE_DSCF_EXE    "dscf"
#define TURBOMOLE_RIDFT_EXE    "ridft"
#define CS3D_EXE      "cs3d"
#define PHARAO_EXE      "align-it"
#define GNUPLOT_EXE      "gnuplot"
#define GNUPLOT_ARGS      "-persist"
#define PYMOL_EXE      "pymol"
#define PYMOL_ARGS      "-p -q"
#define JMOL_EXE      "jmol.sh"
#define JMOL_ARGS      "-Lo -g800x600"
#define FWRITE_WRAP(h, b, n)    fwrite(b, 1, strlen(b), h)
#define FFLUSH_WRAP(h)      fflush(h)
#define FCLOSE_WRAP(h)      fclose(h)
#define ISSEPARATOR(c)      (c == '/')
#endif
#define V2000        1
#define V3000        2
#define MAX_TURBOMOLE_RI_MEMORY    200
#define MAX_TURBOMOLE_SCFITER    80
#define MAX_NAME_LEN      32
#define MAX_FUNC_LEN      64
#define MAX_VAR_BUF      2
#define MAX_THREADS      32
#define MAX_BONDS      10
#define MAX_FF_N      2
#define MAX_FF_PARM      4
#define MAX_FF_TYPE_LEN      16
#define Q2_FIT        0
#define SECV_FIT      1
#define AVE_BUF        0
#define STDDEV_BUF      1
#define FIELD_LIST      0
#define OBJECT_LIST      1
#define ID_LIST        2
#define STRUCT_LIST      3
#define Y_VAR_LIST      4
#define NLEVEL_LIST      5
#define GROUP_LIST      6
#define MAX_LIST      7
#define FROM_FILE      MAX_LIST
#define MOST_ACTIVE_PERCENT    100
#define RT_MAT_SIZE      16
#define RT_VEC_SIZE      4
#define CS3D_DEFAULT_IDELSIG    6
#define CS3D_SKIP_HEADER    22
#define DEFAULT_RANDOM_SEED    0x1234ABCDUL
#define COSMO_DEFAULT_RSOLV    1.30
#define ALMOST_ZERO      1.0e-12
#define SMALL_ENERGY_VALUE    1.0e-06
#define PCA_CONV_THRESHOLD    1.0e-12
#define PLS_CONV_THRESHOLD    1.0e-04
#define MSD_THRESHOLD      1.0e-07
#define ENERGY_THRESHOLD    1.0e-12
#define DEFAULT_MAX_ITER_ALIGN    200
#define MIN_IMPROVEMENT_ITER_ALIGN  0.001
#define DEFAULT_RMSD_ITER_ALIGN    0.001
#define DEFAULT_BEST_PERCENT_ITER_ALIGN  0.95
#define DEFAULT_MIN_PERCENT_ITER_ALIGN  0.60
#define DEFAULT_MAX_FAIL_ALIGN    10
#define MAX_CUTOFF      1.0e35
#define MISSING_VALUE      1.0e37
#define INACTIVE_VALUE      -1.0e37
#define MAX_K_EXCHANGE_ITERATIONS  1000000
#define MAX_DELTA_THRESHOLD    1.0e-06
#define MAX_SDM_ITERATIONS    100
#define MAX_CONF_PER_PHARAO_RUN    1000
#define SDM_THRESHOLD_START    0.7
#define SDM_THRESHOLD_STEP    0.3
#define ALIGN_GOLD_COEFFICIENT    1.2
#define O3_SCORING_FUNCTION_ALPHA  5.0
#define O3_SCORING_FUNCTION_BETA  0.5
#define O3_CHARGE_COEFF      5.0
#define THRESHOLD_DIFF_DISTANCE    0.1
#define THRESHOLD_SIMMETRY_COST    100
#define GRID_TOLERANCE      0.0001
#define CHARGE_WEIGHT      10.0
#define TEMPLATE_ON_ROWS    0
#define TEMPLATE_ON_COLS    1
#define SKIP_MOL_NAME      -1
#define RANDOM_WEIGHTS      1
#define EVEN_WEIGHTS      2
#define ANGLE_STEP      60
#define RANDOM_TRANS_COEFF    5.0
#define MIN_Y_VAR_SD      0.1
#define BLOCK_COMPARE      1
#define O3_COMPRESS_GZIP    1
#define O3_COMPRESS_ZIP      2
#define NEED_STDIN_NORMAL    (1<<0)
#define NEED_STDIN_LEAVE_READ_PIPE_OPEN    (1<<1)
#define NORMAL_FILE_HANDLE    (1<<0)
#define ZIP_FILE_HANDLE      (1<<1)
#define GZIP_FILE_HANDLE    (1<<2)
#define ZIP_MODE_READ      (1<<3)
#define ZIP_MODE_WRITE      (1<<4)
#define VERBOSE_BIT      (1<<0)
#define APPEND_BIT      (1<<1)
#define MATCH_ATOM_TYPES_BIT    (1<<0)
#define MATCH_CONFORMERS_BIT    (1<<1)
#define CENTER_TO_ORIGIN_BIT    (1<<1)
#define ALIGN_RANDOM_BIT    (1<<0)
#define ALIGN_PHARAO_BIT    (1<<1)
#define ALIGN_ATOMBASED_BIT    (1<<2)
#define ALIGN_MIXED_BIT      (1<<3)
#define ALIGN_MULTICONF_TEMPLATE_BIT  (1<<4)
#define ALIGN_KEEP_BEST_TEMPLATE_BIT  (1<<6)
#define ALIGN_MULTICONF_CANDIDATE_BIT  (1<<7)
#define ALIGN_TOGGLE_HYBRID_BIT    (1<<8)
#define ALIGN_TOGGLE_MERGE_BIT    (1<<9)
#define ALIGN_PRINT_RMSD_BIT    (1<<10)
#define ALIGN_TOGGLE_LOOP_BIT    (1<<11)
#define ALIGN_ITERATIVE_TEMPLATE_BIT  (1<<12)
#define FILTER_INTRA_CONF_DB_BIT  (1<<0)
#define FILTER_INTER_CONF_DB_BIT  (1<<1)
#define IMPORT_Y_VARS_BIT    (1<<0)
#define ALLOC_MOL_INFO_BIT    (1<<2)
#define LABEL_INACTIVE_BIT    (1<<0)
#define ZERO_INACTIVE_BIT    (1<<1)
#define LABEL_MISSING_BIT    (1<<2)
#define WEIGHT_BIT      (1<<0)
#define CUTOFF_BIT      (1<<1)
#define CHECK_IF_ACTIVE_BIT    (1<<2)
#define ACTIVE_BIT      (1<<0)
#define PLS_BIT        (1<<0)
#define OPERATE_BIT      (1<<1)
#define CV_BIT        (1<<1)
#define PREDICT_BIT      (1<<2)
#define DELETE_BIT      (1<<3)
#define PCA_BIT        (1<<3)
#define SEED_BIT      (1<<4)
#define GROUP_BIT      (1<<5)
#define D_OPTIMAL_BIT      (1<<6)
#define SEL_INCLUDED_BIT    (1<<7)
#define FFDSEL_BIT      (1<<8)
#define UVEPLS_BIT      (1<<9)
#define TWO_LEVEL_BIT      (1<<10)
#define THREE_LEVEL_BIT      (1<<11)
#define FOUR_LEVEL_BIT      (1<<12)
#define SCRAMBLE_BIT      (1<<13)
#define SDF_BIT        (1<<14)
#define CALC_LEVERAGE_BIT      (1<<0)
#define CALC_FIELD_CONTRIB_BIT      (1<<1)
#define QMD_KEEP_INITIAL    (1<<0)
#define QMD_GBSA      (1<<1)
#define QMD_REMOVE_FOLDER    (1<<2)
#define QMD_ALIGN      (1<<3)
#define QMD_REMOVE_DUPLICATES    (1<<4)
#define QMD_DONT_SUPERPOSE    (1<<5)
#define OBJECT_ASSIGNED      (1<<0)
#define OBJECT_FINISHED      (1<<1)
#define OBJECT_COPIED      (1<<2)
#define OBJECT_ALREADY_DONE    (1<<3)
#define OBJECT_TO_BE_DONE    (1<<4)
#define OBJECT_CHANGED      (1<<5)
#define FFDSEL_DUMMY      (1<<0)
#define FFDSEL_FIXED      (1<<1)
#define FFDSEL_UNCERTAIN    (1<<2)
#define FFDSEL_EXCLUDED      (1<<3)
#define X_SCORES      (1<<0)
#define Y_SCORES      (1<<1)
#define O3_PARAM_NUMERIC    (1<<0)
#define O3_PARAM_STRING      (1<<1)
#define O3_PARAM_FILE      (1<<2)
#define O3_PARAM_DIRECTORY    (1<<3)
#define O3_PARAM_PC      (1<<4)
#define O3_PARAM_DESIGN_POINTS    (1<<5)
#define O3_PARAM_N_CPUS      (1<<6)
#define O3_PARAM_PROBE      (1<<7)
#define O3_PARAM_SCRAMBLE_MAX_BINS  (1<<8)
#define O3_PARAM_SRD_SEEDS    (1<<9)
#define O3_PARAM_REF_Y_VAR    (1<<10)
#define O3_COMPLETION_KEYWORD    (1<<0)
#define O3_COMPLETION_PARAMETER    (1<<1)
#define O3_COMPLETION_VALUE    (1<<2)
#define O3_COMPLETION_QUOTE    (1<<3)
#define NLEVEL_INIT      5
#define NAME_LABEL      1
#define NUMBER_LABEL      2
#define ID_LABEL      4
#define CUSTOM_SCALE      1
#define AUTO_SCALE      2
#define BUW_SCALE      4
#define VDW_FIELD      (1<<0)
#define MM_ELE_FIELD      (1<<1)
#define QM_ELE_FIELD      (1<<2)
#define QM_DEN_FIELD      (1<<3)
#define MD_GRID_FIELD      (1<<4)
#define CS3D_FIELD      (1<<5)
#define PREP_GAUSSIAN_INPUT    (1<<6)
#define PREP_GAMESS_INPUT    (1<<7)
#define PREP_FIREFLY_INPUT    (1<<8)
#define PREP_TURBOMOLE_INPUT    (1<<9)
#define PREP_MOLDEN_INPUT    (1<<10)
#define PREP_SYBYL_INPUT    (1<<11)
#define PREP_MOE_GRID_INPUT    (1<<12)
#define O3_LI_ROWSOL      0
#define O3_LI_COLSOL      1
#define O3_LI_FREE      2
#define O3_LI_COLLIST      3
#define O3_LI_MATCHES      4
#define O3_LI_D        5
#define O3_LI_V        6
#define O3_LI_PRED      7
#define O3_MMFF94      0
#define O3_MD_GRID      1
#define O3_RESTRICTED      'R'
#define O3_UNRESTRICTED      'U'
#define O3_HF        "HF"
#define O3_DFT_B3LYP      "B3LYP"
#define O3_STO_3G      (1<<0)
#define O3_3_21G      (1<<1)
#define O3_SV        (1<<1)
#define O3_6_31G      (1<<2)
#define O3_SVP        (1<<2)
#define O3_6_311G      (1<<3)
#define O3_TZVP        (1<<3)
#define O3_6_311GXX      (1<<4)
#define EMSL_BASIS_SET      (1<<5)
#define MULTIPLY      1
#define DIVIDE        2
#define SUM        3
#define SUBTRACT      4
#define OPPOSITE      5
#define ABSOLUTE_VALUE      6
#define POWER        7
#define LOGARITHM_BASE_10    8
#define LOGARITHM_BASE_E    9
#define LOGARITHM_BASE_N    10
#define HIGHEST_VALUE      1
#define LOWEST_VALUE      2
#define ANY_OBJECT      0
#define MATCHING_OBJECTS    1
#define CONST_DIELECTRIC    0
#define DIST_DEP_DIELECTRIC    1
#define PLS_PLOT      (1<<0)
#define PCA_PLOT      (1<<1)
#define VS_EXP_PLOT      (1<<2)
#define RECALC_VS_EXP      (1<<3)
#define EXT_PRED_VS_EXP      (1<<4)
#define PRED_VS_EXP      (1<<5)
#define SDEC_VS_PC_PLOT      (1<<6)
#define R2_VS_PC_PLOT      (1<<7)
#define SDEP_VS_PC_PLOT      (1<<8)
#define Q2_VS_PC_PLOT      (1<<9)
#define PLS_X_VS_Y_SCORES_PLOT    (1<<10)
#define LOADINGS_PLOT      (1<<11)
#define WEIGHTS_PLOT      (1<<12)
#define SCORES_PLOT      (1<<13)
#define RESIDUALS      (1<<14)
#define SCRAMBLED_Q2_VS_R2    (1<<15)
#define SCRAMBLED_SECV_VS_R2    (1<<16)
#define WEIGHTS        'w'
#define LOADINGS      'l'
#define PCA_LOADINGS      'p'
#define COEFFICIENTS      'c'
#define MEAN_X_COEFFICIENTS      'm'
#define SD_X_COEFFICIENTS      's'
#define TWO_LEVEL      '2'
#define THREE_LEVEL      '3'
#define FOUR_LEVEL      '4'
#define D_OPTIMAL      'd'
#define FFDSEL        'f'
#define FIELD_SD        'F'
#define UVEPLS        'u'
#define SEEDS        's'
#define GROUPS        'g'
#define OBJECT_FIELD      'o'
#define SILENT        0
#define FIELD_LIST      0
#define CUTOFF_MIN      1
#define CUTOFF_MAX      2
#define CUTOFF_QUADRATIC    4
#define ZERO_POS      0
#define ZERO_NEG      1
#define ZERO_ALL      2
#define MAX_SPACES      20
#define INSIGHT_FORMAT      0
#define SYBYL_FORMAT      (1<<0)
#define MAESTRO_FORMAT      (1<<1)
#define MOE_FORMAT      (1<<2)
#define ASCII_FORMAT      (1<<3)
#define XYZ_FORMAT      (1<<4)
#define FORMATTED_CUBE_FORMAT      (1<<5)
#define UNFORMATTED_CUBE_FORMAT      (1<<6)
#define OPENDX_FORMAT      (1<<7)
#define FORMATTED_CUBE_INPUT_FILE  1
#define UNFORMATTED_CUBE_INPUT_FILE  2
#define MOLDEN_INPUT_FILE    3
#define MOE_GRID_INPUT_FILE    4
#define GRID_ASCII_INPUT_FILE    5
#define OPENDX_INPUT_FILE    6
#define COSMO_INPUT_FILE    7
#define GAUSSIAN_UNIT_NUMBER    30
#define DONT_USE_WEIGHTS    0
#define USE_MMFF_WEIGHTS    1
#define USE_CHARGE_WEIGHTS    2


typedef struct CLIArgs CLIArgs;
typedef struct O3Data O3Data;
typedef struct ThreadInfo ThreadInfo;
typedef struct GridInfo GridInfo;
typedef struct CharMat CharMat;
typedef struct IntMat IntMat;
typedef struct MemList MemList;
typedef struct ArrayList ArrayList;
typedef struct CIMatList CIMatList;
typedef struct MatList MatList;
typedef struct PermList PermList;
typedef struct VecList VecList;
typedef struct VarCoord VarCoord;
typedef struct FileDescriptor FileDescriptor;
typedef struct XData XData;
typedef struct YData YData;
typedef struct RegexData RegexData;
typedef struct SeedDistMat SeedDistMat;
typedef struct TaskInfo TaskInfo;
typedef struct FieldInfo FieldInfo;
typedef struct MolInfo MolInfo;
typedef struct BondInfo BondInfo;
typedef struct BondList BondList;
typedef struct AtomInfo AtomInfo;
typedef struct AtomPair AtomPair;
typedef struct ProgExeInfo ProgExeInfo;
typedef struct PyMOLInfo PyMOLInfo;
typedef struct JmolInfo JmolInfo;
typedef struct DoubleMat DoubleMat;
typedef struct DoubleVec DoubleVec;
typedef struct IntPerm IntPerm;
typedef struct NodeInfo NodeInfo;
typedef struct RingInfo RingInfo;
typedef struct RotoTransList RotoTransList;
typedef struct AlignInfo AlignInfo;
typedef struct PharConfInfo PharConfInfo;
typedef struct TemplateInfo TemplateInfo;
typedef struct LAPInfo LAPInfo;
typedef struct QMDInfo QMDInfo;
typedef struct ConfInfo ConfInfo;
typedef struct EnvList EnvList;
typedef struct CationList CationList;
typedef struct FFDSELInfo FFDSELInfo;
typedef struct UVEPLSInfo UVEPLSInfo;
typedef struct ScrambleInfo ScrambleInfo;
typedef struct CVInfo CVInfo;
typedef struct GnuplotInfo GnuplotInfo;
typedef struct fzPtr fzPtr;
#ifdef WIN32
typedef unsigned __int64 uint64_t;
#else
typedef struct EditLineData EditLineData;
struct EditLineData {
  char *prompt;
  char *line;
  int status;
  pthread_t thread_id;
};
#endif

struct fzPtr {
  char *buf;
  int zip_type;
  int pos;
  int data_len;
  FILE *normal_file_handle;
  gzFile gzip_file_handle;
  #if (defined HAVE_LIBMINIZIP) && (defined HAVE_MINIZIP_ZIP_H) && (defined HAVE_MINIZIP_UNZIP_H)
  zipFile zip_file_handle;
  unzFile unz_file_handle;
  #endif
};

struct CationList {
  char **cations;
};

struct EnvList {
  char *name;
  char *ext;
  int flag;
};

struct CLIArgs {
  char prompt;
  int input;
  int output;
  char input_file[BUF_LEN];
  char output_file[BUF_LEN];
};

struct VarCoord {
  int node[3];
  double cart[3];
};

struct FileDescriptor {
  char name[BUF_LEN];
  FILE *handle;
  struct FileDescriptor *next;
};

struct XData {
  int active_x_vars;
  int zero_x_values;
  double min_x_value;
  double max_x_value;
  double min_cutoff;
  double max_cutoff;
  double sdcut_x_var;
  double x_weight_coefficient;
  double temp_x_weight_coefficient;
};

struct YData {
  int zero_y_values;
  double min_y_value;
  double max_y_value;
  double y_weight_coefficient;
  double temp_y_weight_coefficient;
};

struct RegexData {
  int struct_num;
  int conf_num;
};

struct SeedDistMat {
  int seed[2];
  double dist;
};

struct TaskInfo {
  char string[BUF_LEN];
  char file[MAX_FUNC_LEN];
  char func[MAX_FUNC_LEN];
  int code;
  int line;
  int data[MAX_DATA_FIELDS];
};

struct BondInfo {
  int num;
  int order;
};

struct BondList {
  int a[2];
  int order;
};

struct AtomPair {
  int a[2];
  int score;
  int cost;
  double dist;
  double weight;
};

struct RotoTransList {
  int pairs;
  AtomPair *sdm;
};

struct AtomInfo {
  char atom_name[MAX_FF_TYPE_LEN];
  char element[MAX_FF_TYPE_LEN];
  char ring;
  char used;
  int n_bonded;
  int n_arom_bonded;
  int sdf_charge;
  int atom_num;
  int atom_type;
  int tinker_type;
  int match;
  BondInfo bonded[MAX_BONDS];
  int arom_bonded[MAX_BONDS];
  double coord[3];
  double parm[MAX_FF_PARM];
  double charge;
  double formal_charge;
};

struct FieldInfo {
  char theory[MAX_NAME_LEN];
  char mol_dir[BUF_LEN];
  char qm_dir[BUF_LEN];
  char md_grid_dir[BUF_LEN];
  char qm_scratch[BUF_LEN];
  char qm_software[MAX_NAME_LEN];
  char qm_exe[BUF_LEN];
  char qm_exe_path[BUF_LEN];
  #ifdef WIN32
  char mpiexec_exe[BUF_LEN];
  #endif
  char cs3d_exe[BUF_LEN];
  char babel_datadir[BUF_LEN];
  char babel_libdir[BUF_LEN];
  char babel_exe_path[BUF_LEN];
  char md_grid_exe_path[BUF_LEN];
  char force_field;
  char diel_dep;
  char smooth_probe_flag;
  char spin;
  char basis_set;
  char match;
  char compress;
  int type;
  int d_func;
  int p_func;
  int f_func;
  int diff_sp;
  int idelsig;
  int max_n_atoms;
  int max_n_heavy_atoms;
  int max_n_bonds;
  double diel_const;
  double md_grid_cutoff;
  AtomInfo probe;
};

struct AlignInfo {
  char pharao_exe[BUF_LEN];
  char pharao_exe_path[BUF_LEN];
  char template_file[BUF_LEN];
  char candidate_file[BUF_LEN];
  char template_dir[BUF_LEN];
  char candidate_dir[BUF_LEN];
  char template_conf_dir[BUF_LEN];
  char candidate_conf_dir[BUF_LEN];
  char align_dir[BUF_LEN];
  char filter_conf_dir[BUF_LEN];
  char align_scratch[BUF_LEN];
  int type;
  int filter_type;
  int n_tasks;
  int max_iter;
  int max_fail;
  double level;
  double gold;
};

struct PharConfInfo {
  char delete;
  int n_phar_points;
  int object_num;
  int conf_num;
};

struct TemplateInfo {
  int num;
  double score;
};

struct LAPInfo {
  int *array[O3_MAX_SLOT];
  int **cost;
  double **diff;
};

struct NodeInfo {
  NodeInfo *next;
  int atom_id;
  int source;
};
  
struct RingInfo {
  int *atom_id;
  int size;
  int arom;
  int ele;
};
  
struct MolInfo {
  char object_name[MAX_NAME_LEN];
  int sdf_version;
  int object_id;
  int name_count;
  int n_atoms;
  int n_heavy_atoms;
  int n_bonds;
  int done;
  int struct_num;
  int conf_num;
  double score;
  double ln_k;
  double exp_g_minus_ln_k;
  AtomInfo **atom;
  #ifdef WIN32
  HANDLE hMapHandle;
  #endif
};
  
struct QMDInfo {
  char tinker_exe_path[BUF_LEN];
  char tinker_prm_path[BUF_LEN];
  char src[BUF_LEN];
  char dest[BUF_LEN];
  char qmd_dir[BUF_LEN];
  char minimizer[MAX_NAME_LEN];
  int options;
  int runs;
  int min_maxiter;
  double diel_const;
  double min_grad;
  double rmsd;
  double range;
  double temperature;
  double window;
  double time_step;
};
  
struct ConfInfo {
  AtomInfo **atom;
  int n_conf;
  int n_atoms;
  int n_heavy_atoms;
  int **h;
  double *coord;
  double energy;
};

struct PyMOLInfo {
  char pymol_exe[BUF_LEN];
  char use_pymol;
  int grid_box;
  #ifndef WIN32
  char **proc_env;
  pid_t pymol_pid;
  int pymol_pipe;
  FILE *pymol_handle;
  #else
  char *proc_env;
  DWORD pymol_pid;
  HANDLE pymol_handle;
  #endif
};

struct JmolInfo {
  char jmol_exe[BUF_LEN];
  char use_jmol;
  int grid_box;
  int port;
  struct sockaddr_in sockaddr;
  #ifndef WIN32
  char **proc_env;
  pid_t jmol_pid;
  int sock;
  #else
  char *proc_env;
  DWORD jmol_pid;
  SOCKET sock;
  #endif
};

struct FFDSELInfo {
  int ffdsel_included_vars;
  int cv_type;
  int groups;
  int runs;
  int percent_dummies;
  int use_srd_groups;
  int retain_uncertain;
  int fold_over;
  int print_sdep;
  int print_effect;
  int confidence_level;
  int n_dummies;
  int design_x;
  int design_y;
  int power_of_two;
  double combination_variable_ratio;
};

struct CVInfo {
  int pc_num;
  int overall_cv_runs;
  int num_predictions;
  int n_threads;
  double tss;
  void *cv_thread;
};

struct UVEPLSInfo {
  int save_ram;
  int uvepls_included_vars;
  int cv_type;
  int groups;
  int runs;
  int use_srd_groups;
  int uve_m;
  int ive;
  int ive_percent_limit;
  int ive_external_sdep;
  double uve_alpha;
  double dummy_range_coefficient;
  double dummy_value_coefficient;
};

struct ScrambleInfo {
  int cv_type;
  int groups;
  int runs;
  int scramblings;
  int overall_scramblings;
  int max_bins;
  int min_bins;
  int fit_order;
  int print_runs;
  double critical_point;
};

struct GnuplotInfo {
  char gnuplot_exe[BUF_LEN];
  char use_gnuplot;
};

struct DoubleMat {
  int m;
  int n;
  int max_m;
  int max_n;
  double *base;
};

struct DoubleVec {
  int size;
  int max_size;
  double *ve;
};

struct IntPerm {
  int size;
  int max_size;
  int *pe;
};

struct CharMat {
  int m;
  int n;
  char **me;
};

struct IntMat {
  int m;
  int n;
  int **me;
};

struct MemList {
  FileDescriptor **file_descriptor;
  FileDescriptor *source;
  XData *x_data;
  YData *y_data;
  uint16_t *object_attr;
  uint16_t *field_attr;
  uint16_t **x_var_attr;
  uint16_t *y_var_attr;
  char *line;
  char *list_orig;
  char *list_copy;
  int *struct_list;
  unsigned long *random_seed_array;
  double *sum;
  double *object_weight;
  double *weighted_value;
  double *score_temp;
  double *field_contrib;
  float ***x_var_array;
  float *float_xy_mat;
  float *buf_float_xy_mat[4];
  float *out_float_xy_mat;
  float *y_var_array;
  double **x_var_buf[MAX_VAR_BUF];
  double *y_var_buf[MAX_VAR_BUF];
  int *ipiv;
  #ifndef HAVE_LIBSUNPERF
  double *work;
  #endif
  double *s;
  #ifndef WIN32
  pthread_mutex_t *mutex;
  #else
  HANDLE *mutex;
  #endif
  ThreadInfo *thread_info[MAX_THREADS];
  unsigned char *ffdsel_status;
  char *ffdsel_included;
  char *uvepls_included;
  int *struct_per_group;
  int *predicted_object_list;
  int *seed_count;
  int *seed_count_before_collapse;
  int *voronoi_fill;
  int *voronoi_active;
  int *group_zero;
  int *exponents;
  int *binary;
  int *dummy_random_list;
  int *bin_populations;
  int *candidate_pos;
  int *per_object_template;
  int *per_object_template_temp;
};

struct ArrayList {
  char **done_objects;
  int **voronoi_composition;
  double **score_matrix;
  MolInfo **mol_info;
  TemplateInfo **candidate_template_object_list;
  VarCoord **seed_coord;
  SeedDistMat **nearest_mat;
  PharConfInfo **phar_conf_list;
  TaskInfo **task_list;
  RotoTransList **rt_list;
  RegexData **regex_list[MAX_STATES];
};

struct CIMatList {
  CharMat *arg;
  CharMat *y_var_name;
  CharMat *ffd_design_mat;
  IntMat *seed_list;
  IntMat *voronoi_buf;
  IntMat *collapse_checked;
  IntMat **group_composition_list;
};

struct MatList {
  DoubleMat *e_mat;
  DoubleMat *large_e_mat;
  DoubleMat *hat_temp1_mat;
  DoubleMat *hat_temp2_mat;
  DoubleMat *hat_mat;
  DoubleMat *x_weights;
  DoubleMat *large_e_mat_ave;
  DoubleMat *b_coefficients;
  DoubleMat *b_coefficients_ave;
  DoubleMat *b_coefficients_sd;
  DoubleMat *b_coefficients_store;
  DoubleMat *f_mat;
  DoubleMat *pred_f_mat;
  DoubleMat *large_f_mat;
  DoubleMat *ordered_f_mat;
  DoubleMat *scrambled_f_mat;
  DoubleMat *large_f_mat_ave;
  DoubleMat *x_scores;
  DoubleMat *x_loadings;
  DoubleMat *x_loadings_tr;
  DoubleMat *design_mat;
  DoubleMat *temp;
  DoubleMat *y_loadings;
  DoubleMat *y_scores;
  DoubleMat *sorted_candidates_mat;
  DoubleMat *support_mat;
  DoubleMat *x_weights_star;
  DoubleMat *sdep_mat;
  DoubleMat *press;
  DoubleMat *ave_press;
  DoubleMat *cum_ave;
  DoubleMat *pos_ave;
  DoubleMat *neg_ave;
  DoubleMat *y_values;
  DoubleMat *r2_fit_mat;
  DoubleMat *fit_temp;
};

struct VecList {
  DoubleVec *b;
  DoubleVec *u;
  DoubleVec *explained_s2_x;
  DoubleVec *e_mat_ave;
  DoubleVec *e_mat_full_ave;
  DoubleVec *sum1;
  DoubleVec *sumweight;
  DoubleVec *ave;
  DoubleVec *stddev;
  DoubleVec *ss;
  DoubleVec *activity_list;
  DoubleVec *heavy_msd_list;
  DoubleVec *v;
  DoubleVec *c;
  DoubleVec *v_new;
  DoubleVec *ro;
  DoubleVec *y_values_ave;
  DoubleVec *explained_s2_y;
  DoubleVec *dummy_e_mat_full_ave;
  DoubleVec *f_mat_full_ave;
  DoubleVec *f_mat_ave;
  DoubleVec *active_value_ave;
  DoubleVec *ave_sdep;
  DoubleVec *best_sdep;
  DoubleVec *ave_sdec;
  DoubleVec *effect;
  DoubleVec *r2;
  DoubleVec *r2_pred;
  DoubleVec *q2;
  DoubleVec *secv;
  DoubleVec *predictivity_list;
  DoubleVec *xn_vec;
  DoubleVec *temp_xn_vec;
  DoubleVec *dxn_vec;
  DoubleVec *dxi_vec;
  DoubleVec *delta_vec;
  DoubleVec *eff_c;
  DoubleVec *eff_dummy_c;
  DoubleVec *median_vec;
  DoubleVec *fit_vec[2];
  DoubleVec *fit_prod[2];
  DoubleVec *fit_coeff[2];
};

struct PermList {
  IntPerm *out_structs;
  IntPerm *numberlist[MAX_LIST];
  IntPerm *pymol_object_id;
  IntPerm *pymol_old_object_id;
  IntPerm *jmol_object_id;
  IntPerm *jmol_old_object_id;
  IntPerm *activity_rank;
  IntPerm *conf_population[MAX_DB];
  IntPerm *dxn_perm;
  IntPerm *dxi_perm;
  IntPerm *delta_perm;
  IntPerm *ordered_var_list;
  IntPerm *sdep_rank;
  IntPerm *eff_c_rank;
  IntPerm *eff_dummy_c_rank;
  IntPerm *predictivity_list_rank;
  IntPerm *median_perm;
  IntPerm *scrambling_order;
  IntPerm *scrambling_temp;
  IntPerm *best_template_object_list;
  IntPerm *failed_template_object_list;
};

struct GridInfo {
  int x_vars;
  int field_num;
  int object_num;
  int struct_num;
  int nodes[3];
  float start_coord[3];
  float end_coord[3];
  float step[3];
};

struct O3Data {
  char *title;
  char package_code[MAX_NAME_LEN];
  char temp_dir[BUF_LEN];
  char home_dir[BUF_LEN];
  char default_folder[BUF_LEN];
  char prompt;
  char terminal;
  int debug;
  int n_proc;
  int error_code;
  int save_ram;
  int file_num;
  int mti; /* mti==MERSENNE_N+1 means mt[MERSENNE_N] is not initialized */
  int object_num;
  int active_object_num;
  int ext_pred_object_num;
  int field_num;
  int active_field_num;
  int voronoi_num;
  int largest_atom_num;
  int argn;
  int overall_active_x_vars;
  int overall_zero_x_values;
  int overall_zero_y_values;
  int x_vars;
  int y_vars;
  int pc_num;
  int mmap_field_num;
  int mmap_pagesize;
  int object_pagesize;
  uint64_t valid;
  unsigned long random_seed;
  unsigned long mt[MERSENNE_N]; /* the array for the state vector  */
  double max_coord[3];
  double min_coord[3];
  FILE *in;
  FILE *out;
  TaskInfo task;
  FieldInfo field;
  GridInfo grid;
  GridInfo newgrid;
  MemList mel;  
  ArrayList al;  
  MatList mal;
  CIMatList cimal;
  VecList vel;
  PermList pel;
  FileDescriptor **file;
  QMDInfo qmd;
  AlignInfo align;
  PyMOLInfo pymol;
  JmolInfo jmol;
  CVInfo cv;
  FFDSELInfo ffdsel;
  UVEPLSInfo uvepls;
  GnuplotInfo gnuplot;
  ScrambleInfo scramble;
  #ifndef WIN32
  struct termios *user_termios;
  pthread_t thread_id[MAX_THREADS];
  void *thread_result[MAX_THREADS];
  #else
  DWORD dwThreadIdArray[MAX_THREADS];
  HANDLE hThreadArray[MAX_THREADS];
  HANDLE hInput;
  HANDLE hOutput;
  #endif
};

struct ThreadInfo {
  int thread_num;
  int n_calc;
  int start;
  int end;
  int model_type;
  int pc_num;
  int groups;
  int data[MAX_DATA_FIELDS];
  int cannot_write_temp_file;
  FileDescriptor temp_pred;
  FileDescriptor temp_cv_coeff;
  O3Data od;
  O3Data od_comp;
};


void absolute_path(char *string);
int add_to_list(IntPerm **list, int elem);
int align_iterative(O3Data *od);
int align_random(O3Data *od);
int align(O3Data *od);
#ifndef WIN32
void *align_atombased_thread(void *pointer);
void *align_single_pharao_thread(void *pointer);
void *align_multi_pharao_thread(void *pointer);
#else
DWORD align_atombased_thread(void *pointer);
DWORD align_single_pharao_thread(void *pointer);
DWORD align_multi_pharao_thread(void *pointer);
#endif
int alignment_exists(O3Data *od, FileDescriptor *sdf_fd);
char **alloc_array(int n, int size);
CharMat *alloc_char_matrix(CharMat *old_char_mat, int m, int n);
ConfInfo *alloc_conf(int n_atoms);
int alloc_average_mat(O3Data *od, int model_type, int cv_type, int groups, int runs);
int alloc_cv_sdep(O3Data *od, int pc_num, int runs);
int alloc_file_descriptor(O3Data *od, int file_num);
int *alloc_int_array(int *old_ptr, int places);
IntMat *alloc_int_matrix(IntMat *old_int_mat, int m, int n);
int alloc_lap_info(LAPInfo *li, int max_n_atoms);
int alloc_object_attr(O3Data *od, int start);
int alloc_pls(O3Data *od, int x_vars, int pc_num, int model_type);
int prepare_scrambling(O3Data *od);
int alloc_threads(O3Data *od);
int alloc_voronoi(O3Data *od, int places);
int alloc_x_var_array(O3Data *od, int num_fields);
int alloc_y_var_array(O3Data *od);
int autoscale_field(O3Data *od);
int autoscale_y_var(O3Data *od);
int average_x_var(O3Data *od, int field_num);
void average_y_var(O3Data *od);
int bond_in_aromatic_ring(RingInfo **ring, int *a);
int break_sdf_to_mol(O3Data *od, TaskInfo *task, FileDescriptor *from_fd, char *to_dir);
int break_sdf_to_sdf(O3Data *od, TaskInfo *task, FileDescriptor *from_fd, char *to_dir);
int calc_active_vars(O3Data *od, int model_type);
void calc_conf_centroid(ConfInfo *conf, double *centroid);
double calc_delta_ij(O3Data *od, DoubleMat *dispersion_mat, int i, int j);
double calc_dxi(O3Data *od, DoubleMat *dispersion_mat, DoubleMat *candidates_mat, int i);
double calc_dxixj(O3Data *od, DoubleMat *dispersion_mat, int i, int j);
double calc_dxixj(O3Data *od, DoubleMat *dispersion_mat, int i, int j);
int calc_field(O3Data *od, void *thread_func, int prep_or_calc);
void calc_large_mat_ave(O3Data *od, DoubleMat *large_mat, DoubleMat *large_mat_ave, int run);
#ifndef WIN32
void *calc_cosmo_thread(void *pointer);
void *calc_md_grid_thread(void *pointer);
void *calc_mm_thread(void *pointer);
void *calc_qm_thread(void *pointer);
#else
DWORD calc_cosmo_thread(void *pointer);
DWORD calc_md_grid_thread(void *pointer);
DWORD calc_mm_thread(void *pointer);
DWORD calc_qm_thread(void *pointer);
#endif
int call_cs3d_program(O3Data *od);
int call_md_grid_program(O3Data *od);
int call_obenergy(O3Data *od, int force_field);
int calc_p_vectors(O3Data *od, int field_num, int seed_num);
int calc_y_values(O3Data *od, int options);
int check_babel(O3Data *od, char *bin);
int check_bond_type(AtomInfo **atom, int *tinker_types, int *a, int i, BondInfo *bond_info, int value);
int check_conf_db(O3Data *od, char *conf_dir, int type, int *wrong_object_id, int *wrong_conf_num);
int check_define(O3Data *od, char *bin);
int check_duplicate_parameter(int **max_vary, int *parameter);
int check_file_pattern(FILE *handle, char *file_pattern, int object_num);
int check_lmo_parameters(O3Data *od, char *tool_msg, int *groups, int *runs, int run_type, int overall_line_num);
int check_mmap(O3Data *od, int field_num);
int check_pharao(O3Data *od, char *bin);
void *check_readline();
int check_regex_name(char *regex_name, int n_regex);
void close_files(O3Data *od, int from);
int compare(O3Data *od, O3Data *od_comp, int type, int verbose);
#ifndef WIN32
void *compare_thread(void *pointer);
#else
DWORD compare_thread(void *pointer);
#endif
int compare_bond_info(const void *a, const void *b);
int compare_bond_list(const void *a, const void *b);
int compare_conf_energy(const void *a, const void *b);
int compare_corr(const void *a, const void *b);
int compare_dist(const void *a, const void *b);
int compare_integers(const void *a, const void *b);
int compare_n_phar_points(const void *a, const void *b);
int compare_regex_data(const void *a, const void *b);
int compare_score(const void *a, const void *b);
int compare_template_score(const void *a, const void *b);
int compare_seed_dist(const void *a, const void *b);
void compute_conf_h(ConfInfo *conf);
int compute_cost_matrix(LAPInfo *li, ConfInfo *moved_conf, ConfInfo *template_conf, int n_bins, int coeff, int options);
int convert_mol(O3Data *od, char *from_filename, char *to_filename, char *from_ext, char *to_ext, char *flags);
void copy_plane_to_buffer(O3Data *od, float *float_xy_mat, float *buf_float_xy_mat);
int create_box(O3Data *od, GridInfo *temp_grid, double outgap, int from_file);
int create_design_support_matrices(O3Data *od, DoubleMat *candidates_mat, int design_points);
int cutoff(O3Data *od, int type, double cutoff);
int cv(O3Data *od, int suggested_pc_num, int model_type, int cv_type, int groups, int runs);
void double_mat_sort_clean_exit(DoubleVec **vec_in, DoubleVec **vec_short, DoubleVec **vec_temp, IntPerm *perm, int columns);
int d_optimal(O3Data *od, int factors, int design_points, int type);
int d_optimal_select(O3Data *od, int pc_num, int percent_remove, int type);
int d_optimal_select(O3Data *od, int pc_num, int percent_remove, int type);
DoubleMat *double_mat_sort(DoubleMat *mat_in);
void double_mat_free(DoubleMat *double_mat);
void double_vec_free(DoubleVec *double_vec);
DoubleVec *double_vec_resize(DoubleVec *double_vec, int size);
DoubleVec *double_vec_sort(DoubleVec *x, IntPerm *order);
void determine_best_cpu_number(O3Data *od, char *parameter);
int dexist(char *dirname);
void double_mat_free(DoubleMat *double_mat);
DoubleMat *double_mat_resize(DoubleMat *double_mat, int m, int n);
#ifdef HAVE_LIBATLAS
void dsyev_(char *jobz, char *uplo, int *n, double *a, int *lda, double *w, double *work, int *lwork, int *info);
#endif
void elapsed_time(O3Data *od, struct timeval *start, struct timeval *end);
int energy(O3Data *od);
#ifndef WIN32
void *energy_thread(void *pointer);
#else
DWORD energy_thread(void *pointer);
#endif
int exclude(O3Data *od, int type, int ref_field);
void ext_program_wait(ProgExeInfo *prog_exe_info, int pid);
int ext_program_exe(ProgExeInfo *prog_exe_info, int *error);
int fcopy(char *from_filename, char *to_filename, char *mode);
int fexist(char *filename);
int ffdsel(O3Data *od, int pc);
#ifndef WIN32
void *ffdsel_thread(void *pointer);
#else
DWORD ffdsel_thread(void *pointer);
#endif
int fgrep(FILE *handle, char *buffer, char *grep_key);
int fill_atom_info(O3Data *od, TaskInfo *task, AtomInfo **atom, BondList **bond_list, int object_num, char force_field);
int fill_date_string(char *date_string);
int fill_md_grid_types(AtomInfo **atom);
#ifndef WIN32
char **fill_env(O3Data *od, EnvList personalized_env[], char *bin, int object_num);
#else
char *fill_env(O3Data *od, EnvList personalized_env[], char *bin, int object_num);
#endif
int fill_numberlist(O3Data *od, int len, int type);
int fill_tinker_bond_info(O3Data *od, FileDescriptor *inp_fd, AtomInfo **atom, BondList **bond_list, int object_num);
int fill_tinker_types(AtomInfo **atom);
int fill_thread_info(O3Data *od, int n_tasks);
int fill_x_matrix_pca(O3Data *od);
int fill_x_matrix(O3Data *od, int model_type, int use_srd_groups);
int fill_x_matrix_scrambled(O3Data *od);
void fill_x_vector(O3Data *od, int object_num, int row, int model_type, int cv_run);
int fill_y_matrix(O3Data *od);
int fill_y_matrix_scrambled(O3Data *od);
void fill_y_vector(O3Data *od, int object_num, int row, int model_type, int cv_run);
int filter(O3Data *od);
int filter_extract_split_phar_thread(void *pointer);
int filter_inter_thread(void *pointer);
int filter_intra_thread(void *pointer);
int filter_sol_vector(LAPInfo *li, ConfInfo *moved_conf, ConfInfo *template_conf, AtomPair *temp_sdm, AtomPair *sdm);
int find_atom_type(O3Data *od, int nb_pos, AtomInfo *atom);
int find_conformation_in_sdf(FILE *handle_in, FILE *handle_out, int conf_num);
int find_vary_speed(O3Data *od, char *name_list, int **max_vary, int **vary, int *field_num, int *object_num, VarCoord *varcoord);
void fix_endianness(void *chunk, int chunk_len, int word_size, int swap_endianness);
int fmove(char *filename1, char *filename2);
void free_cv_groups(O3Data *od, int runs);
void free_cv_sdep(O3Data *od);
void free_parallel_cv(O3Data *od, ThreadInfo **thread_info, int model_type, int cv_type, int runs);
void free_pls(O3Data *od);
void free_array(void *array);
void free_atom_array(O3Data *od);
void free_char_matrix(CharMat *char_mat);
void free_conf(ConfInfo *conf);
void free_lap_info(LAPInfo *li);
void free_mem(O3Data *od);
void free_node(NodeInfo *fnode, int **path, RingInfo **ring, int n_atoms);
void free_threads(O3Data *od);
void free_x_var_array(O3Data *od);
void free_y_var_array(O3Data *od);
char *get_basename_no_ext(char *filename);
#if (defined HAVE_LIBMINIZIP) && (defined HAVE_MINIZIP_ZIP_H) && (defined HAVE_MINIZIP_UNZIP_H)
int zipFiletime(char *filename, zip_fileinfo *zfi);
zipFile zipOpenWrite(char *filename);
unzFile zipOpenRead(char *filename);
int zipCloseWrite(unzFile handle);
int zipCloseRead(zipFile handle);
#endif
fzPtr *fzopen(char *filename, char *mode);
int fzclose(fzPtr *fz_ptr);
int fzputs(fzPtr *fz_ptr, char *data);
char *fzgets(char *data, int len, fzPtr *fz_ptr);
int fzread(void *data, size_t size, size_t count, fzPtr *fz_ptr);
int fzwrite(void *data, size_t size, size_t count, fzPtr *fz_ptr);
void fzrewind(fzPtr *fz_ptr);
int fzseek(fzPtr *fz_ptr, long int offset, int whence);
unsigned long genrand_int32(O3Data *od);
double genrand_real(O3Data *od);
int get_alignment_score(O3Data *od, FileDescriptor *fd, int object_num, double *score, int *best_template_object_num);
char *get_args(O3Data *od, char *parameter_name);
void get_attr_struct_ave(O3Data *od, int y_var, uint16_t attr, int *attr_struct_num, double *attr_value_ave);
char *get_basename(char *filename);
int get_current_time(char *time_string);
int get_cv_coeff(O3Data *od, int cv_run, int y, int x, double *cv_coeff, int save_ram);
int get_datafile_coord(O3Data *od, FileDescriptor *data_fd, int n_atom, int n_total_atoms, int *cube_word_size, double *data_coord, int datafile_type);
char *get_dirname(char *filename);
int get_double_from_ascii_file(FILE *handle, char *read_buffer, int read_buffer_size, char **context, char **ptr, int *back, double *value);
uint16_t get_field_attr(O3Data *od, int field_num, uint16_t attr);
int get_gridkont_data_points(O3Data *od, int new_model, int replace_object_name, int endianness_switch, int dry_run);
int get_number_of_procs();
int get_n_atoms_bonds(MolInfo *mol_info, FILE *handle, char *buffer);
uint16_t get_object_attr(O3Data *od, int object_num, uint16_t attr);
#ifdef WIN32
BOOL GetOSDisplayString(LPTSTR pszOS, int *page_size);
#endif
void get_system_information(O3Data *od);
int get_voronoi_buf(O3Data *od, int field_num, int x_var);
int get_x_value(O3Data *od, int field_num, int object_num, int x_var, double *value, int flag);
uint16_t get_x_var_attr(O3Data *od, int field_num, int x_var, uint16_t attr);
double get_x_var_buf(O3Data *od, int field_num, int x_var, int buf_num);
double get_y_value(O3Data *od, int object_num, int y_var, int flag);
uint16_t get_y_var_attr(O3Data *od, int y_var, uint16_t attr);
double get_y_var_buf(O3Data *od, int y_var, int buf_num);
char *get_y_var_name(char *buffer, char *y_name);
int grid_write(O3Data *od, char *filename, int pc_num, int type, int sign, int format, int label, int interpolate, int requested_endianness);
int import_dependent(O3Data *od, char *name_list);
int import_free_format(O3Data *od, char *name_list, int skip_header, int *n_values);
int import_grid_ascii(O3Data *od, char *regex_name);
int import_opendx(O3Data *od, char *regex_name);
int import_grid_formatted_cube(O3Data *od, int mo);
int import_grid_unformatted_cube(O3Data *od, int mo);
int import_gridkont(O3Data *od, int replace_object_name);
int import_grid_moe(O3Data *od, char *regex_name);
int import_grid_molden(O3Data *od);
void init_cv_sdep(O3Data *od);
void init_genrand(O3Data *od, unsigned long s);
void init_pls(O3Data *od);
void int_perm_free(IntPerm *int_perm);
IntPerm *int_perm_resize(IntPerm *int_perm, int size);
DoubleMat *int_perm_rows(IntPerm *perm, DoubleMat *double_mat1, DoubleMat *double_mat2);
DoubleVec *int_perm_vec(IntPerm *perm, DoubleVec *double_vec1, DoubleVec *double_vec2);
int intlog2(int n);
int is_aromatic_bond(BondList **bond_list, int a1, int a2);
int is_in_list(IntPerm *list, int elem);
int is_in_path(char *program, char *path_to_program);
int join_aligned_files(O3Data *od, int done_array_pos, char *error_filename);
int join_mol_to_sdf(O3Data *od, TaskInfo *task, FileDescriptor *to_fd, char *from_dir);
int join_thread_files(O3Data *od, ThreadInfo **thread_info);
int k_exchange(O3Data *od, DoubleMat *dispersion_mat);
void lap(LAPInfo *li, int dim);
#ifndef WIN32
void *lmo_cv_thread(void *pointer);
void *loo_cv_thread(void *pointer);
void *lto_cv_thread(void *pointer);
#else
DWORD lmo_cv_thread(void *pointer);
DWORD loo_cv_thread(void *pointer);
DWORD lto_cv_thread(void *pointer);
#endif
int load_dat(O3Data *od, int file_id, int options);
int machine_type();
int match_grids(O3Data *od);
int match_objects_with_datafile(O3Data *od, char *file_pattern, int datafile_type);
#ifndef HAVE_MKDTEMP
char *mkdtemp(char *tmpl);
#endif
#ifndef HAVE_MKSTEMP
int mkstemp(char *tmpl);
#endif
int mol_to_sdf(O3Data *od, int object_num, double actual_value);
int nlevel(O3Data *od);
char *o3_completion_generator(const char *text, int state);
char **o3_completion_matches(const char *text, int start, int end);
#ifdef WIN32
void o3_compentry_free(void *mem);
#endif
char *o3_get_keyword(int *keyword_len);
int open_perm_dir(O3Data *od, char *root_dir, char *id_string, char *perm_dir_name);
int open_temp_dir(O3Data *od, char *root_dir, char *id_string, char *temp_dir_name);
int open_temp_file(O3Data *od, FileDescriptor *file_descriptor, char *id_string);
void overall_msd(AtomPair *sdm, int pairs, ConfInfo *moved_conf, ConfInfo *template_conf, double *heavy_msd);
int parallel_cv(O3Data *od, int x_vars, int suggested_pc_num, int model_type, int cv_type, int groups, int runs);
int parse_comma_hyphen_list_to_array(O3Data *od, char *list, int list_number);
double parse_grid_ascii_line(char *line, char *parsed_line, VarCoord *varcoord);
int parse_input(O3Data *od, FILE *input_stream, int run_type);
int parse_o3_line(char *buffer);
int parse_sdf(O3Data *od, int options, char *name_list);
void parse_sdf_coord_line(int sdf_version, char *buffer, char *element, double *coord, int *charge);
int parse_synonym_lists(O3Data *od, char *tool_name, char *tool_msg, int synonym_list, int *list_type, int default_list, int run_type, int overall_line_num);
int pca(O3Data *od, int pc_num);
double pearson_r(DoubleMat *mat, int *x, int check_missing, double missing);
double perform_operation(int type, double value, double factor);
#ifndef WIN32
void *phar_extract_thread(void *pointer);
#else
DWORD phar_extract_thread(void *pointer);
#endif
int plot(O3Data *od, char *filename, int type, int label, int requested_pc_num, int *pc_axis);
void pls(O3Data *od, int suggested_pc_num, int model_type);
int pred_ext_y_values(O3Data *od, int pc_num, int model_type);
int predict(O3Data *od, int pc_num);
int pred_y_values(O3Data *od, ThreadInfo *ti, int pc_num, int model_type, int cv_run);
int preload_best_templates(O3Data *od, FileDescriptor *fd, int *skip, double *score);
int prepare_cv(O3Data *od, int pc_num, int cv_type, int groups, int runs);
void prepare_design_model(O3Data *od, int design_row);
void prepare_rototrans_matrix(double *rt_mat, double *t_mat1, double *t_mat2, double *rad);
int prep_cosmo_input(O3Data *od, TaskInfo *task, AtomInfo **atom, int object_num);
int prep_cs3d_input(O3Data *od);
void prep_moe_grid_input(O3Data *od);
int prep_molden_input(O3Data *od, int object_num);
int prep_qm_input(O3Data *od, TaskInfo *task, AtomInfo **atom, int object_num);
void prep_sybyl_input(O3Data *od);
int print_calc_values(O3Data *od, int options);
void print_debug_info(O3Data *od, TaskInfo *task);
int print_ext_pred_values(O3Data *od);
void print_grid_comparison(O3Data *od);
void print_grid_coordinates(O3Data *od, GridInfo *grid_info);
int print_pred_values(O3Data *od);
void print_pls_scores(O3Data *od, int options);
int print_variables(O3Data *od, int type);
#ifndef WIN32
void program_signal_handler(int signum);
#else
BOOL program_signal_handler(DWORD fdwCtrlType);
#endif
void pseudo_seed_coord(O3Data *od, int field_num, int *seed);
int qmd(O3Data *od);
#ifndef WIN32
void *qmd_thread(void *pointer);
#else
DWORD qmd_thread(void *pointer);
#endif
int read_dx_header(O3Data *od, FileDescriptor *inp_fd, int object_num);
void read_tinker_xyz_n_atoms_energy(char *line, int *n_atoms, double *energy);
int realloc_x_var_array(O3Data *od, int old_object_num);
int realloc_y_var_array(O3Data *od, int old_object_num);
int reload_coefficients(O3Data *od, int pc_num);
int reload_weights_loadings(O3Data *od);
void remove_box(O3Data *od);
void remove_exe(char *string);
void remove_extension(char *string);
int remove_field(O3Data *od);
int remove_from_list(IntPerm **list, int elem);
void remove_newline(char *string);
int remove_object(O3Data *od);
void remove_recursive(char *filename);
void remove_temp_files(char *basename);
int remove_with_prefix(char *temp_dir_string, char *prefix);
int remove_x_vars(O3Data *od, uint16_t attr);
int remove_y_vars(O3Data *od);
int replace_coord(int sdf_version, char *buffer, double *coord);
void replace_orig_y(O3Data *od);
void reset_user_terminal(O3Data *od);
void restore_orig_y(O3Data *od);
int rms_algorithm(int options, AtomPair *sdm, int pairs, ConfInfo *moved_conf, ConfInfo *template_conf, ConfInfo *fitted_conf, double *rt_mat, double *heavy_msd, double *original_heavy_msd);
int rms_algorithm_multi(O3Data *od, O3Data *od_comp, double *rt_mat, double *heavy_msd);
int rototrans(O3Data *od, char *out_sdf_name, double *trans, double *rot);
int save_dat(O3Data *od, int file_id);
double score_alignment(O3Data *od, ConfInfo *template_conf, ConfInfo *fitted_conf, AtomPair *sdm, int pairs);
int scramble(O3Data *od, int pc_num);
int sdcut(O3Data *od, double threshold);
int sdm_algorithm(AtomPair *sdm, ConfInfo *moved_conf, ConfInfo *template_conf, char **used, int options, double threshold);
int send_jmol_command(O3Data *od, char *command);
int set_sel_included_bit(O3Data *od, int use_srd_groups);
void set_voronoi_buf(O3Data *od, int field_num, int x_var, int voronoi_num);
int srd(O3Data *od, int pc_num, int seed_num, int type, int collapse, double critical_distance, double collapse_distance);
int store_weights_loadings(O3Data *od);
int set(O3Data *od, int type, uint16_t attr, int state, int verbose);
void set_field_attr(O3Data *od, int field_num, uint16_t attr, int onoff);
void set_field_weight(O3Data *od, double weight);
void set_grid_point(O3Data *od, float *float_xy_mat, VarCoord *varcoord, double value);
void set_nice_value(O3Data *od, int nice_value);
void set_object_attr(O3Data *od, int object_num, uint16_t attr, int onoff);
int set_object_weight(O3Data *od, double weight, int list_type, int options);
void set_random_seed(O3Data *od, unsigned long seed);
int set_x_value(O3Data *od, int field_num, int object_num, int x_var, double value);
int set_x_value_unbuffered(O3Data *od, int field_num, int object_num, int x_var, double value);
void set_x_var_attr(O3Data *od, int field_num, int x_var, uint16_t attr, int onoff);
void set_x_var_buf(O3Data *od, int field_num, int x_var, int buf_num, double value);
void set_y_value(O3Data *od, int object_num, int y_var, double value);
void set_y_var_attr(O3Data *od, int y_var, uint16_t attr, int onoff);
void set_y_var_buf(O3Data *od, int y_var, int buf_num, double value);
void set_y_var_weight(O3Data *od, double weight);
void slash_to_backslash(char *string);
double squared_euclidean_distance(double *coord1, double *coord2);
void string_to_lowercase(char *string);
int stddev_x_var(O3Data *od, int field_num);
void stddev_y_var(O3Data *od);
#ifndef HAVE_STRTOK_R
char *strtok_r(char *s1, const char *s2, char **lasts);
#endif
int superpose_conf_lap(LAPInfo *li, ConfInfo *moved_conf, ConfInfo *template_conf, ConfInfo *fitted_conf, ConfInfo *progress_conf, AtomPair *temp_sdm, AtomPair *fitted_sdm, char **used, double *rt_mat, double *heavy_msd, double *original_heavy_msd, int *pairs);
int superpose_conf_syst(ConfInfo *moved_conf, ConfInfo *template_conf, ConfInfo *fitted_conf, ConfInfo *progress_conf, ConfInfo *cand_conf, AtomPair *sdm, AtomPair *local_best_sdm, AtomPair *fitted_sdm, char **used, double *rt_mat, int angle_step, double *heavy_msd, double *original_heavy_msd, int *pairs);
void sync_field_mmap(O3Data *od);
int exe_shell_cmd(O3Data *od, char *command, char *exedir, char *shell);
int tanimoto(O3Data *od, int ref_struct);
void tee_error(O3Data *od, int run_type, int overall_line_num, char *fmt, ...);
void tee_flush(O3Data *od);
void tee_printf(O3Data *od, char *fmt, ...);
int timeval_subtract(struct timeval *result, struct timeval *x, struct timeval *y);
int tinker_analyze(O3Data *od, char *work_dir, char *xyz, int object_num, int conf_num);
int tinker_minimize(O3Data *od, char *work_dir, char *xyz, char *xyz_min, int object_num, int conf_num);
int tinker_dynamic(O3Data *od, char *work_dir, char *xyz, int object_num, int conf_num, unsigned long seed);
int transform(O3Data *od, int type, int operation, double value);
void trim_mean_center_x_matrix_pca(O3Data *od);
void trim_mean_center_matrix(O3Data *od, DoubleMat *large_mat, DoubleMat **mat,
  DoubleVec **mat_ave, int model_type, int active_object_num);
void trim_mean_center_x_matrix_hp(O3Data *od, int model_type, int active_object_num, int run);
void trim_mean_center_y_matrix_hp(O3Data *od, int active_object_num, int run);
int up_n_levels(char *path, int levels);
int update_conf_ln_k(O3Data *od, int model_type, int pc_num, double *ln_k_rmsd, int conv_method);
void update_field_object_attr(O3Data *od, int verbose);
int update_mol(O3Data *od);
int update_pymol(O3Data *od);
int update_jmol(O3Data *od);
int uvepls(O3Data *od, int pc);
int v_intersection(int *v1, int *v2);
int v_union(int *v_union, int *v1, int *v2);
void var_to_xyz(O3Data *od, int x_var, VarCoord *varcoord);
void vertex_xyz(O3Data *od, FILE *handle, int x, int y, int z);
int write_aligned_mol(O3Data *od, O3Data *od_comp, TaskInfo *task, ConfInfo *fitted_conf, int object_num);
void write_ffd_design_matrix_col(O3Data *od, int first_element, int col, int decimal);
int write_grid_plane(O3Data *od, FILE *plane_file, int z_plane, int interpolate, int swap_endianness, float *minVal, float *maxVal);
int write_header(O3Data *od, int object_num, char *header, int format, int interpolate, int swap_endianness);
int write_tinker_energy(FileDescriptor *fd, double energy);
int write_tinker_xyz_bnd(O3Data *od, AtomInfo **atom, BondList **d_list, int n_atoms, int object_num, char *xyz_name, char *bnd_name);
int x_var_buw(O3Data *od);
int xyz_to_var(O3Data *od, VarCoord *varcoord);
int y_var_buw(O3Data *od);
int zero(O3Data *od, int type, double level);
