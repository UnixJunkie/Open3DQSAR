# Makefile.am
# 
# is part of
# 
# Open3DQSAR
# ----------
# 
# An open-source software aimed at high-throughput
# chemometric analysis of molecular interaction fields
# 
# Copyright (C) 2009-2018 Paolo Tosco, Thomas Balle
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
# 
# For further information, please contact:
# 
# Paolo Tosco, PhD
# Dipartimento di Scienza e Tecnologia del Farmaco
# Universita' degli Studi di Torino
# Via Pietro Giuria, 9
# 10125 Torino (Italy)
# Phone:  +39 011 670 7680
# Mobile: +39 348 553 7206
# Fax:    +39 011 670 7687
# E-mail: paolo.tosco@unito.it


lib_LTLIBRARIES = libo3q.la
libo3q_la_SOURCES = \
alloc.c \
autoscale.c \
average.c \
buw.c \
calc_active_vars.c \
calc_field.c \
calc_large_mat_ave.c \
calc_y_values.c \
check_deps.c \
check_regex_name.c \
close_files.c \
cutoff.c \
cv.c \
cv_thread.c \
dcdflib.c \
determine_best_cpu_number.c \
d_optimal.c \
double_mat.c \
double_vec.c \
exclude.c \
ffdsel.c \
file.c \
fill_env.c \
fill_matrix.c \
fill_numberlist.c \
fill_vector.c \
free.c \
get_args.c \
get_attr.c \
get_buf.c \
get_number_of_procs.c \
get_system_information.c \
get_value.c \
grid_box.c \
grid_write.c \
import_dependent.c \
import_free_format.c \
import_grid_ascii.c \
import_opendx.c \
import_grid_formatted_cube.c \
import_grid_moe.c \
import_grid_unformatted_cube.c \
import_gridkont.c \
import_grid_molden.c \
init.c \
int_perm.c \
int_perm_op.c \
load_dat.c \
mersenne_twister.c \
nlevel.c \
parallel_cv.c \
parse_comma_hyphen_list_to_array.c \
parse_sdf.c \
pca.c \
plot.c \
pls.c \
predict.c \
pred_y_values.c \
prepare_cv.c \
print_values.c \
print_variables.c \
qsort_functions.c \
reload_coeff_weights_loadings.c \
remove.c \
safe_rint.c \
save_dat.c \
scramble.c \
sdcut.c \
set_attr.c \
set_buf.c \
set.c \
set_nice_value.c \
set_random_seed.c \
set_sel_included_bit.c \
set_value.c \
set_weight.c \
srd.c \
stddev.c \
store_weights_loadings.c \
tanimoto.c \
tee.c \
time.c \
transform.c \
trim_mean_center_matrix.c \
update_field_object_attr.c \
update_jmol.c \
update_pymol.c \
utils.c \
uvepls.c \
var_to_xyz.c \
zero.c \
include/cdflib.h \
include/basis_set.h \
include/error_messages.h \
include/extern.h \
include/ff_parm.h \
include/interpolated_value.h \
include/interpolation_coefficients.h \
include/ipmpar.h \
include/nice_windows.h \
include/o3header.h \
include/proc_env.h \
include/prog_exe_info.h \
include/rl_runtime.h \
include/safe_rint.h

bin_PROGRAMS = \
open3dqsar
open3dqsar_SOURCES = \
main.c \
include/cdflib.h \
include/basis_set.h \
include/error_messages.h \
include/extern.h \
include/ff_parm.h \
include/interpolated_value.h \
include/interpolation_coefficients.h \
include/ipmpar.h \
include/keywords.h \
include/nice_windows.h \
include/o3header.h \
include/proc_env.h \
include/prog_exe_info.h \
include/rl_runtime.h \
include/safe_rint.h \
parse_o3q_input.c

open3dqsar_LDADD = \
-lo3q \
$(WINDOWS_RES) \
$(LIBLAPACK_ATLAS) \
$(LIBLAPACKE) \
$(LIBLAPACK_F77) \
$(LIBATLAS) \
$(LIBBLAS) \
$(LIBMKL) \
$(LIBSUNPERF) \
$(LIBACCELERATE) \
$(LIBEDIT)
