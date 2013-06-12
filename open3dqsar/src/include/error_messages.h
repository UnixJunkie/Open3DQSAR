/*

error_messages.h

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


char M_NUMBER_OF_CPUS[] =
  "The number of CPUs used by %s has been set to %d.\n\n";
char M_TOOL_INVOKE[] =
  "BGN COMMAND #%02d.%04d - %s tool was invoked as follows:\n"
  "> %s\n\n";
char M_TOOL_SUCCESS[] =
  "END COMMAND #%02d.%04d - %s tool succeeded.\n";
char M_EXECUTABLE[] =
  "The %s executable has been set to \"%s\".\n\n";
char M_EXECUTABLE_PATH[] =
  "The path to %s tools has been set to \"%s\".\n\n";
char M_INPUT_OUTPUT_LOG_DIR[] =
  "The %s directory where input, output and log "
  "files will be put is:\n\"%s\"\n\n";
char E_CANNOT_CHANGE_DIR[] =
  "Cannot change directory to \"%s\".\n%s";
char E_DIR_NOT_EXISTING[] =
  "The directory \"%s\" does not exist or "
  "cannot be accessed.\n%s";
char E_IMPORT_MOLFILE_FIRST[] =
  "Please import some molecule structures first.\n%s";
char E_IMPORT_COSMOTHERM_FIRST[] =
  "Please import COSMOtherm data first.\n%s";
char E_DEFINE_SPACED_PATH[] =
  "The \"define\" program in your TURBOMOLE installation "
  "seems to be outdated; please contact COSMOlogic to get an "
  "updated version.\n%s";
char E_SPECIFY_FILE_TYPE[] =
  "Please specify a type for "
  "the file you wish to import.\n%s";
char E_SPECIFY_STRUCTURE_FILE[] =
  "Please specify a SDF file containing the "
  "structure(s) you wish to %s.\n%s";
char E_GRID_BOX_FIRST[] =
  "Please set a grid box first.\n%s";
char E_CANNOT_ALTER_GRID_BOX[] =
  "The grid box cannot be changed/removed "
  "once fields are present.\n%s\n\n";
char E_CANNOT_ALTER_OBJECT_COORDINATES[] =
  "The object coordinates cannot be altered "
  "once fields are present.\n%s\n\n";
char E_NO_GRID_BOX_PRESENT[] =
  "No grid box is currently present.\n%s";
char E_AT_LEAST_ONE_PC[] =
  "At least 1 PC must be chosen.\n%s";
char E_POSITIVE_NUMBER[] =
  "The %s must "
  "be a positive number.\n%s";
char E_NO_FIELDS_PRESENT[] =
  "No fields are currently present.\n%s";
char E_NO_Y_VARS_PRESENT[] =
  "No Y variables are currently present.\n%s";
char E_NOT_ENOUGH_OBJECTS[] =
  "Expecting exactly %d files with name \"%s\".\n%s";
char E_OPENBABEL_NOT_WORKING[] =
  "A simple OpenBabel test failed with the following error:\n";
char E_OPENBABEL_DATA_PLUGINS[] =
  "OpenBabel data/plugins could not be found.\n"
  "This issue can be fixed setting the "BABEL_DATADIR_ENV" and "
  BABEL_LIBDIR_ENV" environment "
  "variables to suitable values.\n%s\n";
char E_OPENBABEL_MISSING_OR_TOO_OLD[] =
  "OpenBabel binaries could not be found.\n"
  "Please consider downloading and installing the "
  "openbabel_for_open3dtools package "
  "for your system from the "PACKAGE_NAME" website.\n%s\n";
char E_OPENBABEL_PATH[] =
  "Please set the O3_BABEL_PATH environment variable "
  "or use the \"env babel_path\" keyword to indicate the path "
  "to OpenBabel executables.\n%s";
char E_EXECUTABLE_PATH[] =
  "Cannot find the %s binary. Please make sure "
  "it is in the executable path or copy "
  "it in the %s folder.\n%s";
char E_PROGRAM_ERROR[] =
  "%s reported the error message which follows:\n";
char E_CANNOT_READ_PROGRAM_LOG[] =
  "Cannot read %s log file.\n%s";
char E_CS3D_EXE[] =
  "Please set the O3_CS3D environment variable "
  "or use the \"env cs3d\" keyword to the absolute "
  "path to the CS3D executable.\n%s";
char E_MD_GRID_PATH[] =
  "Please set the O3_MD_GRID_PATH environment variable "
  "or use the \"env md_grid_path\" keyword to indicate the path "
  "to Molecular Discovery GRID executables.\n%s";
char E_MD_GRID_ERROR[] =
  "%s reported the error message which follows.\n%s";
char E_ERROR_IN_READING_OB_OUTPUT[] =
  "Cannot read OpenBabel output file \"%s\".\n%s";
char E_ERROR_IN_READING_MOL_FILE[] =
  "Cannot read MOL file \"%s\".\n%s";
char E_ERROR_IN_READING_SDF_FILE[] =
  "Cannot read SDF file \"%s\".\n%s";
char E_ERROR_IN_WRITING_SDF_FILE[] =
  "Cannot write SDF file \"%s\".\n%s";
char E_TEMP_FILE_CANNOT_BE_OPENED_FOR_WRITING[] =
  "The temporary file \"%s\" cannot be opened for writing.\n%s";
char E_TEMP_DIR_CANNOT_BE_CREATED[] =
  "The folder \"%s\" cannot be opened for writing.\n%s";
char E_FILE_CANNOT_BE_OPENED_FOR_WRITING[] =
  "The file \"%s\" cannot be opened for writing.\n%s";
char E_FILE_CANNOT_BE_OPENED_FOR_READING[] =
  "The file \"%s\" cannot be opened for reading.\n%s";
char E_FILE_CORRUPTED_OR_IN_WRONG_FORMAT[] =
  "The %s file \"%s\" appears to be corrupted "
  "or in the wrong format.\n%s";
char E_OBJECTS_NOT_MATCHING[] =
  "None of objects indicated "
  "by the filename pattern \"%s\" "
  "matches object with ID %d (%s) "
  "already loaded in "PACKAGE_NAME".\n%s";
char E_CANNOT_CREATE_PIPE[] =
  "Error in creating a pipe.\n%s";
char E_CANNOT_CREATE_PROCESS[] =
  "Error in spawning a %s process.\n%s";
char E_ERROR_IN_WRITING_TEMP_FILE[] =
  "Error in writing temporary file \"%s\".\n%s";
char E_ERROR_IN_READING_TEMP_FILE[] =
  "Error in reading temporary file \"%s\".\n%s";
char E_CANNOT_READ_WRITE_QM_FILE[] =
  "Cannot %s QM %s file \"%s\"; please check "
  "permissions on your QM folder.\n%s";
char E_CALCULATION_ERROR[] =
  "An error occurred during one or more %s; "
  "please check your input/output/log files.\n%s";
char E_QM_DIR_CANNOT_BE_CREATED[] =
  "A %s folder for QM calculations cannot be created.\n%s";
char E_UNKNOWN_ATOM_TYPE[] =
  "Unknown %s type.\n%s";
char E_Y_VAR_LOW_SD[] =
  "The SD associated with the y variable(s) is too low.\n%s";
char E_OUT_OF_MEMORY[] =
  "Out of memory.\n%s";
char E_NO_PCA_ANALYSIS[] =
  "No PCA analysis is currently present.\n%s";
char E_TOO_FEW_MANY_FOR_CURRENT_PCA_ANALYSIS[] =
  "Too %s were requested with respect "
  "to the current PCA analysis.\n%s";
char E_TOO_FEW_MANY_FOR_AVAILABLE_DATA[] =
  "Too %s were requested with respect "
  "to available data.\n%s";
char E_TOO_FEW_STRUCTURES_FOR_CV[] =
  "There are not enough structures "
  "to carry out the requested CV.\n%s";
char E_ONLY_CUSTOM_AUTO_BUW_ALLOWED[] =
  "Only \"CUSTOM\", \"AUTO\", \"BUW\""
  "scaling types are allowed.\n%s";
char E_THREAD_ERROR[] =
  "Thread error; the code returned "
  "from pthread_%s() is %d.\n%s";
char E_LIST_PARSING[] =
  "Error while parsing the list of %s "
  "on which %s should operate.\n%s";
char E_ALLOWED_OBJECT_RANGE[] =
  "The allowed object range is 1 - %d.\n%s";
char E_CHECK_OBJECT_RANGE[] =
  "Please check your object ID range.\n%s";
char E_SUPPLY_OBJECT_ID_STRUCT_FIELD_LIST[] =
  "Please supply a single keyword among "
  "\"object_list\", \"id_list\", \"struct_list\" "
  "and eventually \"field_list\" "
  "on which %s should operate.\n%s";
char E_NO_VALID_SELECTION[] =
  "No valid %s selection is present.\n%s";
char E_GROUPS_OTHER_THAN_ZERO[] =
  "If groups other than zero should be %s, "
  "only one field may be specified.\n%s";
char E_NOT_ENOUGH_Y_VARS[] =
  "In the file \"%s\" there must be at "
  "least %d y variable%s.\n%s";
char E_WRONG_NUMBER_OF_Y_VARS[] =
  "In the file \"%s\" for all objects there "
  "must be the same number of y variables.\n%s";
char E_CANNOT_FIND_Y_VAR_NAME[] =
  "In the file \"%s\" one of the requested y variable "
  "names does not exist.\n%s";
char E_Y_VAR_ALLOWED_RANGE[] =
  "The allowed y variable range is 1 - %d.\n%s";
char E_STRUCT_ATTRIBUTE_ONLY[] =
  "The \"%s\" attribute can "
  "only be assigned to a structure_list, "
  "not to an object_list, when conformers "
  "are present in the dataset.\n%s";
char E_ONLY_ONE_WILDCARD[] =
  "The basename should contain only one regex "
  "in the format %%(0n)d (n = 1-9).\n%s";
char E_COSMOTHERM_REGEX[] =
  "Please specify a regex name identifying the files "
  "from which you wish to import COSMOtherm data "
  "for ligands in %s.\n%s";
char E_TWO_WILDCARDS[] =
  "The basename should contain two regex "
  "in the format %%(0n)d (n = 1-9).\n%s";
char E_ERROR_IN_GRID_DATA[] =
  "Error on line %d in reading the grid "
  "defined in the file \"%s\".\n%s";
char E_GRIDKONT_NOT_MATCHING[] =
  "The number of objects indicated "
  "in the GRIDKONT file \"%s\" (%d) "
  "does not match the number of objects "
  "already loaded in "PACKAGE_NAME" (%d).\n%s";
char E_GRID_NOT_MATCHING[] =
  "The grid defined in the %s file \"%s\" is not matching "
  "with the one already defined in "PACKAGE_NAME".\n%s";
char E_GRID_NOT_MATCHING_DATA_POINT[] =
  "Data point number %d\n"
  "(%.4f,%.4f,%.4f)\n"
  "in the file \"%s\"\n"
  "%s with respect to the grid "
  "already defined in "PACKAGE_NAME".\n%s";
char E_CANNOT_FIND_MO[] =
  "Cannot find the required MO "
  "in the file \"%s\".\n%s";
char E_NO_OBJECTS_PRESENT[] =
  "No objects are currently present.\n%s";
char E_POLARIZATION_FUNCTIONS[] =
  "The number of polarization %s functions should lie "
  "in the %s range.\n%s";
char E_TRANSFORM_MISSING_OPERATION[] =
  "Please specify the operation which should be carried out "
  "on %c variables.\n%s";
char E_TRANSFORM_MISSING_VALUE[] =
  "Please indicate the value %c variables should be %s.\n%s";
char E_TRANSFORM_MISSING_EXPONENT[] =
  "Please indicate the exponent %c variables "
  "should be elevated to.\n%s";
char E_TRANSFORM_MISSING_BASE[] =
  "Please indicate the base of the logarithm operation "
  "which should be performed on %c variables.\n%s";
char E_ALLOWED_FIELD_RANGE[] =
  "The allowed field range is 1 - %d.\n%s";
char E_WHILE_SOURCING[] =
  "Error while sourcing \"%s\".\n%s";
char E_NO_PLS_MODEL[] =
  "No PLS model is currently present.\n%s";
char E_NO_CV_MODEL[] =
  "No CV model is currently present.\n%s";
char E_NO_EXTERNAL_PREDICTION[] =
  "No external prediction is currently present.\n%s";
char E_NO_SCRAMBLE_ANALYSIS[] =
  "No SCRAMBLE analysis is currently present.\n%s";
char E_TOO_FEW_MANY_FOR_CURRENT_PLS_MODEL[] =
  "Too %s were requested with respect "
  "to the current PLS model.\n%s";
char E_ONLY_LOO_LTO_LMO_CV_ALLOWED[] =
  "Only \"LOO\", \"LTO\", \"LMO\" CV types "
  "are allowed.\n%s";
char E_ONLY_FULL_LOO_LTO_LMO_ALLOWED[] =
  "Only \"FULL\", \"LOO\", \"LTO\", \"LMO\" "
  "COSMOPLS types are allowed.\n%s";
char E_RELOAD_WEIGHTS_LOADINGS[] =
  "An error occurred trying to restore "
  "PLS weights/loadings.\n%s";
char E_PC_ALLOWED_RANGE[] =
  "The allowed PC range is 1 - %d.\n%s";
char E_MISSING_CV_TYPE[] =
  "Please specify the type of CV which "
  "should be performed.\n%s";
char E_PROGRAM_EXIT[] =
  "\nPress ENTER to leave "PACKAGE_NAME".\n";
char O3_FAILED[] =
  PACKAGE_NAME" failed.\n";
char PARSE_INPUT_FAILED[] =
  "Input file parsing failed.\n";
char IMPORT_FAILED[] =
  "IMPORT failed.\n";
char ROTOTRANS_FAILED[] =
  "ROTOTRANS failed.\n";
char BOX_FAILED[] =
  "BOX failed.\n";
char REMOVE_BOX_FAILED[] =
  "REMOVE_BOX failed.\n";
char CALC_FIELD_FAILED[] =
  "CALC_FIELD failed.\n";
char LOAD_FAILED[] =
  "LOAD failed.\n";
char SAVE_FAILED[] =
  "SAVE failed.\n";
char CUTOFF_FAILED[] =
  "CUTOFF failed.\n";
char TRANSFORM_FAILED[] =
  "TRANSFORM failed.\n";
char ZERO_FAILED[] =
  "ZERO failed.\n";
char EXCLUDE_FAILED[] =
  "EXCLUDE failed.\n";
char SDCUT_FAILED[] =
  "SDCUT failed.\n";
char NLEVEL_FAILED[] =
  "NLEVEL failed.\n";
char PCA_FAILED[] =
  "PCA failed.\n";
char SCALE_OBJECT_FAILED[] =
  "SCALE_OBJECT failed.\n";
char SCALE_X_VARS_FAILED[] =
  "SCALE_X_VARS failed.\n";
char SCALE_Y_VARS_FAILED[] =
  "SCALE_Y_VARS failed.\n";
char SET_FAILED[] =
  "SET failed.\n";
char ENV_FAILED[] =
  "ENV failed.\n";
char REMOVE_X_VARS_FAILED[] =
  "REMOVE_X_VARS failed.\n";
char REMOVE_Y_VARS_FAILED[] =
  "REMOVE_Y_VARS failed.\n";
char REMOVE_FIELD_FAILED[] =
  "REMOVE_FIELD failed.\n";
char REMOVE_OBJECT_FAILED[] =
  "REMOVE_OBJECT failed.\n";
char PRINT_FAILED[] =
  "PRINT failed.\n";
char EXPORT_FAILED[] =
  "EXPORT failed.\n";
char PREPARE_FAILED[] =
  "PREPARE failed.\n";
char CHANGE_DIR_FAILED[] =
  "CHANGE_DIR failed.\n";
char SOURCE_FAILED[] =
  "SOURCE failed.\n";
char COSMOPLS_FAILED[] =
  "COSMOPLS failed.\n";
char PLS_FAILED[] =
  "PLS failed.\n";
char CV_FAILED[] =
  "CV failed.\n";
char SCRAMBLE_FAILED[] =
  "SCRAMBLE failed.\n";
char FFDSEL_FAILED[] =
  "FFDSEL failed.\n";
char UVEPLS_FAILED[] =
  "UVEPLS failed.\n";
char PREDICT_FAILED[] =
  "PREDICT failed.\n";
char D_OPTIMAL_FAILED[] =
  "D_OPTIMAL failed.\n";
char SRD_FAILED[] =
  "SRD failed.\n";
char PLOT_FAILED[] =
  "PLOT failed.\n";
