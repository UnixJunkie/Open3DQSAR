/*

keywords.h

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


typedef struct O3ParameterData O3ParameterData;
typedef struct O3KeywordData O3KeywordData;

struct O3ParameterData {
  uint16_t type;
  char *parameter;
  char *choice[MAX_ARG];
};

struct O3KeywordData {
  char *keyword;
  O3ParameterData parameter_data[MAX_ARG];
};

O3KeywordData keyword_data[] =
{
  {
    "box",
    {
      {
        O3_PARAM_STRING, "mode", {
          "SET",
          "GET",
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "outgap", {
          "5.0",
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "step", {
          "1.0",
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "x_start", {
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "x_end", {
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "x_nodes", {
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "y_start", {
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "y_end", {
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "y_nodes", {
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "z_start", {
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "z_end", {
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "z_nodes", {
          NULL
        }
      }, {
        O3_PARAM_FILE, "file", {
          NULL
        }
      }, {  // this is the terminator
        0, NULL, {
          NULL
        }
      }
    }
  }, {
    "calc_field",
    {
      {
        O3_PARAM_STRING, "type", {
          "MM_ELE",
          "VDW",
          "CS3D",
          "MD_GRID",
          "QM_ELE",
          "QM_DEN",
          NULL
        }
      }, {
        O3_PARAM_STRING, "diel_dep", {
          "CONST",
          "DIST",
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "diel_const", {
          "1.0",
          NULL
        }
      }, {
        O3_PARAM_PROBE, "probe_type", {
          NULL
        }
      }, {
        O3_PARAM_STRING, "smooth_probe", {
          "NO",
          "YES",
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "cutoff", {
          "5.0",
          NULL
        }
      }, {
        O3_PARAM_DIRECTORY, "md_grid_dir", {
          NULL
        }
      }, {
        O3_PARAM_STRING, "theory", {
          "HF",
          "DFT",
          "MP2",
          NULL
        }
      }, {
        O3_PARAM_STRING, "basis_set", {
          "6-31G",
          "3-21G",
          "6-311G",
          "STO-3G",
          "SV",
          "SV(P)",
          "TZVP",
          "EMSL_3-21G",
          "EMSL_6-311G",
          "EMSL_6-311Gxx",
          NULL
        }
      }, {
        O3_PARAM_STRING, "spin", {
          "R",
          "U",
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "d_func", {
          "1",
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "p_func", {
          "0",
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "f_func", {
          "0",
          NULL
        }
      }, {
        O3_PARAM_STRING, "diff_sp", {
          "NO",
          "YES",
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "delsig", {
          "6",
          NULL
        }
      }, {
        O3_PARAM_DIRECTORY, "qm_dir", {
          NULL
        }
      }, {
        O3_PARAM_DIRECTORY, "qm_scratch", {
          NULL
        }
      }, {
        O3_PARAM_FILE, "file", {
          NULL
        }
      }, {  // this is the terminator
        0, NULL, {
          NULL
        }
      }
    }
  }, {
    "cd",
    {
      {
        O3_PARAM_DIRECTORY, "dir", {
          NULL
        }
      }, {  // this is the terminator
        0, NULL, {
          NULL
        }
      }
    }
  }, {
    "chdir",
    {
      {
        O3_PARAM_DIRECTORY, "dir", {
          NULL
        }
      }, {  // this is the terminator
        0, NULL, {
          NULL
        }
      }
    }
  }, {
    "cutoff",
    {
      {
        O3_PARAM_STRING, "type", {
          "MIN",
          "MAX",
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "level", {
          NULL
        }
      }, {
        O3_PARAM_STRING, "smooth", {
          "NONE",
          "QUADRATIC",
          NULL
        }
      }, {
        O3_PARAM_STRING, "field_list", {
          "ALL",
          NULL
        }
      }, {  // this is the terminator
        0, NULL, {
          NULL
        }
      }
    }
  }, {
    "cv",
    {
      {
        O3_PARAM_STRING, "type", {
          "LOO",
          "LTO",
          "LMO",
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "runs", {
          "20",
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "groups", {
          "5",
          NULL
        }
      }, {
        O3_PARAM_PC, "pc", {
          NULL
        }
      }, {
        O3_PARAM_FILE, "file", {
          NULL
        }
      }, {  // this is the terminator
        0, NULL, {
          NULL
        }
      }
    }
  }, {
    "dataset",
    {
      {  // this is the terminator
        0, NULL, {
          NULL
        }
      }
    }
  }, {
    "d_optimal",
    {
      {
        O3_PARAM_PC, "pc", {
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "percent_remove", {
          "50",
          NULL
        }
      }, {
        O3_PARAM_DESIGN_POINTS, "design_points", {
          NULL
        }
      }, {
        O3_PARAM_STRING, "type", {
          "WEIGHTS",
          "LOADINGS",
          NULL
        }
      }, {  // this is the terminator
        0, NULL, {
          NULL
        }
      }
    }
  }, {
    "env",
    {
      {
        O3_PARAM_NUMERIC, "random_seed", {
          NULL
        }
      }, {
        O3_PARAM_DIRECTORY, "temp_dir", {
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "nice", {
          #ifndef WIN32
          "20",
          #else
          "BELOW_NORMAL",
          "ABOVE_NORMAL",
          "HIGH",
          "IDLE",
          "NORMAL",
          "REALTIME",
          #endif
          NULL
        }
      }, {
        O3_PARAM_N_CPUS, "n_cpus", {
          NULL
        }
      }, {
        O3_PARAM_DIRECTORY, "babel_path", {
          NULL
        }
      }, {
        O3_PARAM_DIRECTORY, "md_grid_path", {
          NULL
        }
      }, {
        O3_PARAM_FILE, "qm_engine", {
          NULL
        }
      }, {
        O3_PARAM_FILE, "cs3d", {
          NULL
        }
      }, {
        O3_PARAM_FILE, "pymol", {
          NULL
        }
      }, {
        O3_PARAM_FILE, "gnuplot", {
          NULL
        }
      }, {  // this is the terminator
        0, NULL, {
          NULL
        }
      }
    }
  }, {
    "tanimoto",
    {
      {
        O3_PARAM_NUMERIC, "ref_struct", {
          NULL
        }
      }, {
        O3_PARAM_STRING, "field_list", {
          "ALL",
          NULL
        }
      }, {
        O3_PARAM_STRING, "struct_list", {
          "ALL",
          NULL
        }
      }, {  // this is the terminator
        0, NULL, {
          NULL
        }
      }
    }
  }, {
    "exclude",
    {
      {
        O3_PARAM_STRING, "type", {
          "MATCH",
          "ANY",
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "ref_field", {
          NULL
        }
      }, {
        O3_PARAM_STRING, "field_list", {
          "ALL",
          NULL
        }
      }, {  // this is the terminator
        0, NULL, {
          NULL
        }
      }
    }
  }, {
    "exit",
    {
      {  // this is the terminator
        0, NULL, {
          NULL
        }
      }
    }
  }, {
    "stop",
    {
      {  // this is the terminator
        0, NULL, {
          NULL
        }
      }
    }
  }, {
    "quit",
    {
      {  // this is the terminator
        0, NULL, {
          NULL
        }
      }
    }
  }, {
    "export",
    {
      {
        O3_PARAM_STRING, "type", {
          "WEIGHTS",
          "LOADINGS",
          "PCA_LOADINGS",
          "COEFFICIENTS",
          "MEAN_X_COEFFICIENTS",
          "SD_X_COEFFICIENTS",
          "2_LEVEL",
          "3_LEVEL",
          "4_LEVEL",
          "D_OPTIMAL",
          "FFDSEL",
          "UVEPLS",
          "SRD",
          "FIELD_SD",
          "OBJECT_FIELD",
          NULL
        }
      }, {
        O3_PARAM_PC, "pc", {
          NULL
        }
      }, {
        O3_PARAM_STRING, "y_var_list", {
          "ALL",
          NULL
        }
      }, {
        O3_PARAM_STRING, "object_list", {
          "ALL",
          NULL
        }
      }, {
        O3_PARAM_STRING, "id_list", {
          "ALL",
          NULL
        }
      }, {
        O3_PARAM_STRING, "field_list", {
          "ALL",
          NULL
        }
      }, {
        O3_PARAM_STRING, "inactive", {
          "ORIGINAL",
          "ZERO",
          "LABEL",
          NULL
        }
      }, {
        O3_PARAM_STRING, "missing", {
          "ZERO",
          "LABEL",
          NULL
        }
      }, {
        O3_PARAM_STRING, "endianness", {
          "NATIVE",
          "LITTLE_ENDIAN",
          "BIG_ENDIAN",
          NULL
        }
      }, {
        O3_PARAM_FILE, "file", {
          NULL
        }
      }, {
        O3_PARAM_STRING, "format", {
          "INSIGHT",
          "MAESTRO",
          "MOE",
          "SYBYL",
          "FORMATTED_CUBE",
          "UNFORMATTED_CUBE",
          "OPENDX",
          "ASCII",
          "XYZ",
          NULL
        }
      }, {
        O3_PARAM_STRING, "sign", {
          "BOTH",
          "POS",
          "NEG",
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "interpolate", {
          "0",
          NULL
        }
      }, {  // this is the terminator
        0, NULL, {
          NULL
        }
      }
    }
  }, {
    "ffdsel",
    {
      {
        O3_PARAM_STRING, "type", {
          "LOO",
          "LTO",
          "LMO",
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "runs", {
          "20",
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "groups", {
          "5",
          NULL
        }
      }, {
        O3_PARAM_PC, "pc", {
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "percent_dummies", {
          "20",
          NULL
        }
      }, {
        O3_PARAM_STRING, "use_srd_groups", {
          "NO",
          "YES",
          NULL
        }
      }, {
        O3_PARAM_STRING, "retain_uncertain", {
          "YES",
          "NO",
          NULL
        }
      }, {
        O3_PARAM_STRING, "fold_over", {
          "NO",
          "YES",
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "combination_variable_ratio", {
          "2.0",
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "confidence_level", {
          "99.0",
          NULL
        }
      }, {
        O3_PARAM_STRING, "print_sdep", {
          "NO",
          "YES",
          NULL
        }
      }, {
        O3_PARAM_STRING, "print_effect", {
          "NO",
          "YES",
          NULL
        }
      }, {  // this is the terminator
        0, NULL, {
          NULL
        }
      }
    }
  }, {
    "import",
    {
      {
        O3_PARAM_STRING, "type", {
          "SDF",
          "MOL2",
          "DEPENDENT",
          "GRIDKONT",
          "FORMATTED_CUBE",
          "UNFORMATTED_CUBE",
          "MOLDEN",
          "MOE_GRID",
          "GRID_ASCII",
          "OPENDX",
          "FREE_FORMAT",
          NULL
        }
      }, {
        O3_PARAM_FILE, "file", {
          NULL
        }
      }, {
        O3_PARAM_STRING, "y_var_name", {
          "ALL",
          NULL
        }
      }, {
        O3_PARAM_STRING, "replace_object_name", {
          "YES",
          "NO",
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "mo", {
          "ALL",
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "skip_header", {
          NULL
        }
      }, {
        O3_PARAM_STRING, "data_order", {
          "Z_COORD,Y_COORD,X_COORD,N_OBJECT,N_FIELD",
          NULL
        }
      }, {  // this is the terminator
        0, NULL, {
          NULL
        }
      }
    }
  }, {
    "load",
    {
      {
        O3_PARAM_FILE, "file", {
          NULL
        }
      }, {
        O3_PARAM_STRING, "mode", {
          "APPEND",
          "NORMAL",
          NULL
        }
      }, {  // this is the terminator
        0, NULL, {
          NULL
        }
      }
    }
  }, {
    "nlevel",
    {
      {
        O3_PARAM_NUMERIC, "level", {
          "ALL",
          NULL
        }
      }, {
        O3_PARAM_STRING, "field_list", {
          "ALL",
          NULL
        }
      }, {  // this is the terminator
        0, NULL, {
          NULL
        }
      }
    }
  }, {
    "pca",
    {
      {
        O3_PARAM_NUMERIC, "pc", {
          "5",
          NULL
        }
      }, {
        O3_PARAM_FILE, "file", {
          NULL
        }
      }, {  // this is the terminator
        0, NULL, {
          NULL
        }
      }
    }
  }, {
    "plot",
    {
      {
        O3_PARAM_STRING, "type", {
          "RECALC_VS_EXP",
          "PRED_VS_EXP",
          "EXT_PRED_VS_EXP",
          "PLS_X_VS_Y",
          "SDEC",
          "R2",
          "SDEP",
          "Q2",
          "SCRAMBLED_Q2_VS_R2",
          "SCRAMBLED_SECV_VS_R2",
          "PLS_LOADINGS",
          "PLS_WEIGHTS",
          "PLS_SCORES",
          "PCA_LOADINGS",
          "PCA_SCORES",
          NULL
        }
      }, {
        O3_PARAM_PC, "pc", {
          NULL
        }
      }, {
        O3_PARAM_STRING, "y_var_list", {
          "ALL",
          NULL
        }
      }, {
        O3_PARAM_STRING, "residuals", {
          "NO",
          "YES",
          NULL
        }
      }, {
        O3_PARAM_STRING, "label", {
          "NONE",
          "NUMBER",
          "ID",
          "NAME",
          NULL
        }
      }, {
        O3_PARAM_FILE, "file", {
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "pc_x", {
          "1",
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "pc_y", {
          "2",
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "pc_z", {
          NULL
        }
      }, {  // this is the terminator
        0, NULL, {
          NULL
        }
      }
    }
  }, {
    "pls",
    {
      {
        O3_PARAM_PC, "pc", {
          NULL
        }
      }, {
        O3_PARAM_STRING, "scores", {
          "NONE",
          "X",
          "Y",
          "BOTH",
          NULL
        }
      }, {
        O3_PARAM_STRING, "calc_field_contrib", {
          "NO",
          "YES",
          NULL
        }
      }, {
        O3_PARAM_STRING, "calc_leverage", {
          "NO",
          "YES",
          NULL
        }
      }, {
        O3_PARAM_FILE, "file", {
          NULL
        }
      }, {  // this is the terminator
        0, NULL, {
          NULL
        }
      }
    }
  }, {
    "predict",
    {
      {
        O3_PARAM_PC, "pc", {
          NULL
        }
      }, {
        O3_PARAM_STRING, "scores", {
          "NONE",
          "X",
          NULL
        }
      }, {
        O3_PARAM_FILE, "file", {
          NULL
        }
      }, {  // this is the terminator
        0, NULL, {
          NULL
        }
      }
    }
  }, {
    "prepare",
    {
      {
        O3_PARAM_STRING, "type", {
          "CS3D",
          "QM_ELE",
          "QM_DEN",
          "MOE_GRID",
          "SYBYL",
          NULL
        }
      }, {
        O3_PARAM_STRING, "theory", {
          "HF",
          "DFT",
          "MP2",
          NULL
        }
      }, {
        O3_PARAM_STRING, "basis_set", {
          "6-31G",
          "3-21G",
          "6-311G",
          "STO-3G",
          "SV",
          "SVP",
          "TZVP",
          "EMSL_3-21G",
          "EMSL_6-311G",
          "EMSL_6-311Gxx",
          NULL
        }
      }, {
        O3_PARAM_STRING, "spin", {
          "R",
          "U",
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "d_func", {
          "1",
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "p_func", {
          "0",
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "f_func", {
          "0",
          NULL
        }
      }, {
        O3_PARAM_STRING, "diff_sp", {
          "NO",
          "YES",
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "delsig", {
          "6",
          NULL
        }
      }, {
        O3_PARAM_DIRECTORY, "qm_dir", {
          NULL
        }
      }, {
        O3_PARAM_FILE, "file", {
          NULL
        }
      }, {  // this is the terminator
        0, NULL, {
          NULL
        }
      }
    }
  }, {
    "print",
    {
      {
        O3_PARAM_STRING, "type", {
          "2_LEVEL",
          "3_LEVEL",
          "4_LEVEL",
          "D_OPTIMAL",
          "FFDSEL",
          "UVEPLS",
          "SEEDS",
          "GROUPS",
          NULL
        }
      }, {
        O3_PARAM_STRING, "group_list", {
          "ALL",
          NULL
        }
      }, {
        O3_PARAM_STRING, "field_list", {
          "ALL",
          NULL
        }
      }, {  // this is the terminator
        0, NULL, {
          NULL
        }
      }
    }
  }, {
    "remove_box",
    {
      {  // this is the terminator
        0, NULL, {
          NULL
        }
      }
    }
  }, {
    "remove_field",
    {
      {
        O3_PARAM_STRING, "field_list", {
          NULL
        }
      }, {  // this is the terminator
        0, NULL, {
          NULL
        }
      }
    }
  }, {
    "remove_object",
    {
      {
        O3_PARAM_STRING, "object_list", {
          NULL
        }
      }, {
        O3_PARAM_STRING, "id_list", {
          NULL
        }
      }, {  // this is the terminator
        0, NULL, {
          NULL
        }
      }
    }
  }, {
    "remove_x_vars",
    {
      {
        O3_PARAM_STRING, "type", {
          "NLEVEL",
          "D_OPTIMAL",
          "FFDSEL",
          "UVEPLS",
          "GROUPS",
          NULL
        }
      }, {
        O3_PARAM_STRING, "level", {
          "ALL",
          NULL
        }
      }, {
        O3_PARAM_STRING, "group_list", {
          "0",
          NULL
        }
      }, {
        O3_PARAM_STRING, "field_list", {
          "ALL",
          NULL
        }
      }, {  // this is the terminator
        0, NULL, {
          NULL
        }
      }
    }
  }, {
    "remove_y_vars",
    {
      {
        O3_PARAM_STRING, "y_var_list", {
          "ALL",
          NULL
        }
      }, {  // this is the terminator
        0, NULL, {
          NULL
        }
      }
    }
  }, {
    "rototrans",
    {
      {
        O3_PARAM_NUMERIC, "x_trans", {
          "0.0",
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "y_trans", {
          "0.0",
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "z_trans", {
          "0.0",
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "x_rot", {
          "0.0",
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "y_rot", {
          "0.0",
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "z_rot", {
          "0.0",
          NULL
        }
      }, {
        O3_PARAM_FILE, "file", {
          NULL
        }
      }, {  // this is the terminator
        0, NULL, {
          NULL
        }
      }
    }
  }, {
    "save",
    {
      {
        O3_PARAM_FILE, "file", {
          NULL
        }
      }, {  // this is the terminator
        0, NULL, {
          NULL
        }
      }
    }
  }, {
    "scale_object",
    {
      {
        O3_PARAM_STRING, "object_list", {
          "ALL",
          "FILE",
          NULL
        }
      }, {
        O3_PARAM_STRING, "id_list", {
          "ALL",
          "FILE",
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "weight", {
          "1.0",
          "RANDOM",
          NULL
        }
      }, {
        O3_PARAM_FILE, "file", {
          NULL
        }
      }, {  // this is the terminator
        0, NULL, {
          NULL
        }
      }
    }
  }, {
    "scale_x_vars",
    {
      {
        O3_PARAM_STRING, "type", {
          "CUSTOM",
          "AUTO",
          "BUW",
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "weight", {
          "1.0",
          NULL
        }
      }, {
        O3_PARAM_STRING, "field_list", {
          "ALL",
          NULL
        }
      }, {  // this is the terminator
        0, NULL, {
          NULL
        }
      }
    }
  }, {
    "scale_y_vars",
    {
      {
        O3_PARAM_STRING, "type", {
          "CUSTOM",
          "AUTO",
          "BUW",
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "weight", {
          "1.0",
          NULL
        }
      }, {
        O3_PARAM_STRING, "y_var_list", {
          "ALL",
          NULL
        }
      }, {  // this is the terminator
        0, NULL, {
          NULL
        }
      }
    }
  }, {
    "scramble",
    {
      {
        O3_PARAM_STRING, "type", {
          "LOO",
          "LTO",
          "LMO",
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "runs", {
          "20",
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "groups", {
          "5",
          NULL
        }
      }, {
        O3_PARAM_PC, "pc", {
          NULL
        }
      }, {
        O3_PARAM_SCRAMBLE_MAX_BINS, "max_bins", {
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "min_bins", {
          "2",
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "scramblings", {
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "fit_order", {
          "2",
          "3",
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "critical_point", {
          "0.85",
          NULL
        }
      }, {
        O3_PARAM_STRING, "print_runs", {
          "NO",
          "YES",
          NULL
        }
      }, {  // this is the terminator
        0, NULL, {
          NULL
        }
      }
    }
  }, {
    "sdcut",
    {
      {
        O3_PARAM_NUMERIC, "level", {
          NULL
        }
      }, {
        O3_PARAM_STRING, "field_list", {
          "ALL",
          NULL
        }
      }, {  // this is the terminator
        0, NULL, {
          NULL
        }
      }
    }
  }, {
    "set",
    {
      {
        O3_PARAM_STRING, "field_list", {
          NULL
        }
      }, {
        O3_PARAM_STRING, "object_list", {
          NULL
        }
      }, {
        O3_PARAM_STRING, "id_list", {
          NULL
        }
      }, {
        O3_PARAM_STRING, "attribute", {
          "INCLUDED",
          "EXCLUDED",
          "TRAININGSET",
          "TESTSET",
          NULL
        }
      }, {  // this is the terminator
        0, NULL, {
          NULL
        }
      }
    }
  }, {
    "source",
    {
      {
        O3_PARAM_FILE, "file", {
          NULL
        }
      }, {  // this is the terminator
        0, NULL, {
          NULL
        }
      }
    }
  }, {
    "srd",
    {
      {
        O3_PARAM_PC, "pc", {
          NULL
        }
      }, {
        O3_PARAM_STRING, "type", {
          "WEIGHTS",
          "LOADINGS",
          NULL
        }
      }, {
        O3_PARAM_SRD_SEEDS, "seeds", {
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "critical_distance", {
          "1.0",
          NULL
        }
      }, {
        O3_PARAM_STRING, "collapse", {
          "YES",
          "NO",
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "collapse_distance", {
          "2.0",
          NULL
        }
      }, {  // this is the terminator
        0, NULL, {
          NULL
        }
      }
    }
  }, {
    "system",
    {
      {
        O3_PARAM_STRING, "cmd", {
          NULL
        }
      }, {
        O3_PARAM_DIRECTORY, "exedir", {
          NULL
        }
      }, {
        O3_PARAM_STRING, "shell", {
          NULL
        }
      }, {  // this is the terminator
        0, NULL, {
          NULL
        }
      }
    }
  }, {
    "transform",
    {
      {
        O3_PARAM_STRING, "type", {
          "X",
          "Y",
          NULL
        }
      }, {
        O3_PARAM_STRING, "field_list", {
          "ALL",
          NULL
        }
      }, {
        O3_PARAM_STRING, "y_var_list", {
          "ALL",
          NULL
        }
      }, {
        O3_PARAM_STRING, "operation", {
          "MUL",
          "DIV",
          "SUM",
          "SUB",
          "OPP",
          "ABS",
          "POW",
          "LOG10",
          "LN",
          "LOG_BASE_N",
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "value", {
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "exponent", {
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "base", {
          NULL
        }
      }, {  // this is the terminator
        0, NULL, {
          NULL
        }
      }
    }
  }, {
    "uvepls",
    {
      {
        O3_PARAM_STRING, "type", {
          "LOO",
          "LTO",
          "LMO",
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "runs", {
          "20",
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "groups", {
          "5",
          NULL
        }
      }, {
        O3_PARAM_PC, "pc", {
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "dummy_range_coefficient", {
          "1.0",
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "dummy_value_coefficient", {
          "1.0e-10",
          NULL
        }
      }, {
        O3_PARAM_STRING, "use_srd_groups", {
          "NO",
          "YES",
          NULL
        }
      }, {
        O3_PARAM_STRING, "save_ram", {
          "NO",
          "YES",
          NULL
        }
      }, {
        O3_PARAM_STRING, "uve_m", {
          "NO",
          "YES",
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "uve_alpha", {
          "0.0",
          NULL
        }
      }, {
        O3_PARAM_STRING, "ive", {
          "NO",
          "YES",
          NULL
        }
      }, {
        O3_PARAM_NUMERIC, "ive_percent_limit", {
          "100",
          NULL
        }
      }, {
        O3_PARAM_STRING, "ive_external_pred", {
          "NO",
          "YES",
          NULL
        }
      }, {  // this is the terminator
        0, NULL, {
          NULL
        }
      }
    }
  }, {
    "zero",
    {
      {
        O3_PARAM_NUMERIC, "level", {
          NULL
        }
      }, {
        O3_PARAM_STRING, "sign", {
          "BOTH",
          "POS",
          "NEG",
          NULL
        }
      }, {
        O3_PARAM_STRING, "field_list", {
          "ALL",
          NULL
        }
      }, {  // this is the terminator
        0, NULL, {
          NULL
        }
      }
    }
  }, {  // this is the terminator
    NULL,
    {
      {
        O3_PARAM_STRING, NULL, {
          NULL
        }
      }
    }
  }
};
