/*

remove.c

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


int remove_field(O3Data *od)
{
  char buffer[BUF_LEN];
  int i;
  int k;
  int deleted_fields;
  int result;
  

  deleted_fields = 0;
  for (i = 0; i < od->field_num; ++i) {
    if (get_field_attr(od, i, OPERATE_BIT)) {
      set_field_attr(od, i, DELETE_BIT, 1);
      set_field_attr(od, i, ACTIVE_BIT, 0);
      ++deleted_fields;
    }
    else {
      for (k = 0; k < od->x_vars; ++k) {
        set_x_var_attr(od, i, k,
          GROUP_BIT | SEED_BIT | SEL_INCLUDED_BIT
          | FFDSEL_BIT | UVEPLS_BIT | D_OPTIMAL_BIT, 0);
      }
    }
  }
  od->valid &= (SDF_BIT | COSMOTHERM_BIT);
  if (deleted_fields < od->field_num) {
    result = open_temp_dir(od, NULL, "temp_dat", buffer);
    if (result) {
      return CANNOT_WRITE_TEMP_FILE;
    }
    sprintf(od->file[TEMP_DAT]->name, "%s%ctemp.dat", buffer, SEPARATOR);
    od->file[TEMP_DAT]->handle = (FILE *)
      fzopen(od->file[TEMP_DAT]->name, "wb");
    if (!(od->file[TEMP_DAT]->handle)) {
      return CANNOT_WRITE_TEMP_FILE;
    }
    result = save_dat(od, TEMP_DAT);
    if (result) {
      return PREMATURE_DAT_EOF;
    }
    od->file[TEMP_DAT]->handle = (FILE *)
      fzopen(od->file[TEMP_DAT]->name, "rb");
    if (!(od->file[TEMP_DAT]->handle)) {
      return CANNOT_READ_TEMP_FILE;
    }
    if (od->field_num) {
      close_files(od, MAX_FILES);
      free_x_var_array(od);
    }
    result = load_dat(od, TEMP_DAT, SILENT);
    if (result) {
      return result;
    }
    remove_recursive(buffer);
    od->file[TEMP_DAT]->name[0] = '\0';
  }
  else {
    od->valid = SDF_BIT;
    od->field_num = 0;
    od->active_field_num = 0;
    od->active_object_num = od->object_num;
    od->ext_pred_object_num = 0;
    if (IS_O3Q(od)) {
      od->ffdsel.use_srd_groups = 0;
      od->uvepls.use_srd_groups = 0;
      od->voronoi_num = 0;
      free_array(od->al.voronoi_composition);
      od->al.voronoi_composition = NULL;
    }
  }

  return 0;
}


int remove_object(O3Data *od)
{
  char buffer[BUF_LEN];
  int i;
  int j;
  int k;
  int deleted_objects;
  int result;
  

  deleted_objects = 0;
  for (i = 0; i < od->grid.object_num; ++i) {
    if (get_object_attr(od, i, OPERATE_BIT)) {
      set_object_attr(od, i, DELETE_BIT, 1);
      set_object_attr(od, i, ACTIVE_BIT, 0);
      ++deleted_objects;
    }
    else {
      for (j = 0; j < od->field_num; ++j) {
        for (k = 0; k < od->x_vars; ++k) {
          set_x_var_attr(od, j, k, GROUP_BIT | SEED_BIT
            | SEL_INCLUDED_BIT | TWO_LEVEL_BIT
            | THREE_LEVEL_BIT | FOUR_LEVEL_BIT
            | FFDSEL_BIT | D_OPTIMAL_BIT | UVEPLS_BIT, 0);
        }
      }
    }
  }
  od->valid &= (SDF_BIT | COSMOTHERM_BIT);
  if (deleted_objects < od->grid.object_num) {
    result = open_temp_dir(od, NULL, "temp_dat", buffer);
    if (result) {
      return CANNOT_WRITE_TEMP_FILE;
    }
    sprintf(od->file[TEMP_DAT]->name, "%s%ctemp.dat", buffer, SEPARATOR);
    od->file[TEMP_DAT]->handle = (FILE *)
      fzopen(od->file[TEMP_DAT]->name, "wb");
    if (!(od->file[TEMP_DAT]->handle)) {
      return CANNOT_WRITE_TEMP_FILE;
    }
    result = save_dat(od, TEMP_DAT);
    if (result) {
      return PREMATURE_DAT_EOF;
    }
    od->file[TEMP_DAT]->handle = (FILE *)
      fzopen(od->file[TEMP_DAT]->name, "rb");
    if (!(od->file[TEMP_DAT]->handle)) {
      return CANNOT_READ_TEMP_FILE;
    }
    if (od->field_num) {
      close_files(od, MAX_FILES);
      free_x_var_array(od);
    }
    sprintf(od->file[TEMP_MOLFILE]->name, "%s%ctemp_loaded_molfile.sdf",
      od->field.mol_dir, SEPARATOR);
    result = load_dat(od, TEMP_DAT, SILENT);
    memset(od->file[TEMP_MOLFILE]->name, 0, BUF_LEN);
    remove(od->file[TEMP_MOLFILE]->name);
    remove_recursive(buffer);
    if (result) {
      return result;
    }
  }
  else {
    od->valid = 0;
    od->field_num = 0;
    od->active_field_num = 0;
    od->object_num = 0;
    od->grid.object_num = 0;
    od->grid.struct_num = 0;
    od->active_object_num = 0;
    od->ext_pred_object_num = 0;
    if (IS_O3Q(od)) {
      od->ffdsel.use_srd_groups = 0;
      od->uvepls.use_srd_groups = 0;
      od->voronoi_num = 0;
      free_array(od->al.voronoi_composition);
      od->al.voronoi_composition = NULL;
    }
    memset(&(od->grid), 0, sizeof(GridInfo));
    free_y_var_array(od);
  }

  return 0;
}


int remove_x_vars(O3Data *od, unsigned short attr)
{
  unsigned short level_bit = 0;
  int i;
  int j;
  int k;
  int n = 0;
  int len;
  IntPerm *numberlist;
  int group_num;
  int field_num;
  int old_active_x_vars;
  int voronoi_fill;
  int overall_seed_count;
  int *voronoi_composition;
  IntPerm *field_list;
  IntPerm *group_list;


  switch (attr) {
  
    case TWO_LEVEL_BIT:
    numberlist = od->pel.numberlist[NLEVEL_LIST];
    len = numberlist->size;
    if ((numberlist->pe[0] < 2) || (numberlist->pe[len - 1] > 4)) {
      return INVALID_LIST_RANGE;
    }
    for (j = 0; j < len; ++j) {
      if (numberlist->pe[j] == 2) {
        level_bit = TWO_LEVEL_BIT;
        n = 2;
      }
      if (numberlist->pe[j] == 3) {
        level_bit = THREE_LEVEL_BIT;
        n = 3;
      }
      if (numberlist->pe[j] == 4) {
        level_bit = FOUR_LEVEL_BIT;
        n = 4;
      }
      tee_printf(od, "Field          Number of %d-level variables removed\n", n);
      tee_printf(od, "--------------------------------------------------\n");
      for (i = 0; i < od->field_num; ++i) {
        if (get_field_attr(od, i, ACTIVE_BIT)) {
          if (get_field_attr(od, i, OPERATE_BIT)) {
            n = 0;
            for (k = 0; k < od->x_vars; ++k) {
              if (get_x_var_attr(od, i, k, level_bit)) {
                ++n;
                set_x_var_attr(od, i, k, DELETE_BIT, 1);
              }
            }
            tee_printf(od, "%5d%45d\n", i + 1, n);
          }
        }
      }
      tee_printf(od, "--------------------------------------------------\n\n");
    }
    break;
    
    case D_OPTIMAL_BIT:
    tee_printf(od, "\n%5s%36s\n", "Field", "Number of x variables removed");
    tee_printf(od, "-----------------------------------------\n");
    for (i = 0; i < od->field_num; ++i) {
      if (get_field_attr(od, i, ACTIVE_BIT)) {
        if (get_field_attr(od, i, OPERATE_BIT)) {
          n = 0;
          for (k = 0; k < od->x_vars; ++k) {
            if (get_x_var_attr(od, i, k, D_OPTIMAL_BIT)) {
              ++n;
              set_x_var_attr(od, i, k, DELETE_BIT, 1);
            }
          }
          tee_printf(od, "%5d%36d\n", i + 1, n);
        }
      }
    }
    tee_printf(od, "-----------------------------------------\n\n");
    break;
    
    case FFDSEL_BIT:
    tee_printf(od, "\n%5s%36s\n", "Field", "Number of x variables removed");
    tee_printf(od, "-----------------------------------------\n");
    for (i = 0; i < od->field_num; ++i) {
      if (get_field_attr(od, i, ACTIVE_BIT)) {
        if (get_field_attr(od, i, OPERATE_BIT)) {
          old_active_x_vars = od->mel.x_data[i].active_x_vars;
          od->mel.x_data[i].active_x_vars = 0;
          for (k = 0; k < od->x_vars; ++k) {
            set_x_var_attr(od, i, k, DELETE_BIT, 1);
            if (get_x_var_attr(od, i, k, SEL_INCLUDED_BIT)
              && (!get_x_var_attr(od, i, k, FFDSEL_BIT))) {
              ++(od->mel.x_data[i].active_x_vars);
              set_x_var_attr(od, i, k, DELETE_BIT, 0);
            }
          }
          tee_printf(od, "%5d%36d\n", i + 1,
            old_active_x_vars
            - od->mel.x_data[i].active_x_vars);
          od->overall_active_x_vars +=
            (od->mel.x_data[i].active_x_vars - old_active_x_vars);
        }
      }
    }
    tee_printf(od, "-----------------------------------------\n\n");
    break;

    case UVEPLS_BIT:
    tee_printf(od, "\n%5s%36s\n", "Field", "Number of x variables removed");
    tee_printf(od, "-----------------------------------------\n");
    for (i = 0; i < od->field_num; ++i) {
      if (get_field_attr(od, i, ACTIVE_BIT)) {
        if (get_field_attr(od, i, OPERATE_BIT)) {
          old_active_x_vars = od->mel.x_data[i].active_x_vars;
          od->mel.x_data[i].active_x_vars = 0;
          for (k = 0; k < od->x_vars; ++k) {
            set_x_var_attr(od, i, k, DELETE_BIT, 1);
            if (get_x_var_attr(od, i, k, SEL_INCLUDED_BIT)
              && (!get_x_var_attr(od, i, k, UVEPLS_BIT))) {
              ++(od->mel.x_data[i].active_x_vars);
              set_x_var_attr(od, i, k, DELETE_BIT, 0);
            }
          }
          tee_printf(od, "%5d%36d\n", i + 1,
            old_active_x_vars - od->mel.x_data[i].active_x_vars);
          od->overall_active_x_vars +=
            (od->mel.x_data[i].active_x_vars - old_active_x_vars);
        }
      }
    }
    tee_printf(od, "-----------------------------------------\n\n");
    break;

    case SEED_BIT:
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
    overall_seed_count = 0;
    for (i = 0; i < (field_list->pe[0] - 1); ++i) {
      overall_seed_count += od->mel.seed_count[i];
    }
    tee_printf(od, "\n%5s%12s%36s\n", "Field", "Group", "Number of x variables removed");
    tee_printf(od, "-----------------------------------------------------\n");
    for (i = 0; i < od->field_num; ++i) {
      if (get_field_attr(od, i, ACTIVE_BIT)) {
        if (get_field_attr(od, i, OPERATE_BIT)) {
          for (j = 0; j < group_num; ++j) {
            if (group_list->pe[j] == 0) {
              n = 0;
              for (k = 0; k < od->x_vars; ++k) {
                if (!get_x_var_attr(od, i, k, GROUP_BIT | SEED_BIT)) {
                  ++n;
                  set_x_var_attr(od, i, k, DELETE_BIT, 1);
                }
              }
            }
            else {
              voronoi_fill = od->mel.voronoi_fill
                [overall_seed_count + group_list->pe[j] - 1];
              voronoi_composition = od->al.voronoi_composition
                [overall_seed_count + group_list->pe[j] - 1];
              n = 0;
              for (k = 0; k < voronoi_fill; ++k) {
                if (get_x_var_attr(od, i,
                  voronoi_composition[k], SEED_BIT | GROUP_BIT)) {
                  set_x_var_attr(od, i, k, DELETE_BIT, 1);
                  ++n;
                }
              }
            }
            tee_printf(od, "%5d%12d%36d\n", i + 1, group_list->pe[j], n);
          }
        }
      }
    }
    tee_printf(od, "-----------------------------------------------------\n\n");
    break;
  }
  for (i = 0; i < od->field_num; ++i) {
    for (k = 0; k < od->x_vars; ++k) {
      set_x_var_attr(od, i, k, GROUP_BIT | SEED_BIT
      | SEL_INCLUDED_BIT | TWO_LEVEL_BIT | THREE_LEVEL_BIT
      | FOUR_LEVEL_BIT | FFDSEL_BIT | D_OPTIMAL_BIT | UVEPLS_BIT, 0);
    }
  }
  od->valid &= (SDF_BIT | COSMOTHERM_BIT);
  od->ffdsel.use_srd_groups = 0;
  od->uvepls.use_srd_groups = 0;
  od->voronoi_num = 0;
  free_array(od->al.voronoi_composition);
  od->al.voronoi_composition = NULL;

  return 0;
}


int remove_y_vars(O3Data *od)
{
  int i;
  int j;
  int all;
  int result;
  

  result = open_temp_file(od, od->file[DEP_IN], "temp_y_var");
  if (result) {
    return CANNOT_WRITE_TEMP_FILE;
  }
  for (j = 0, all = 1; all && (j < od->y_vars); ++j) {
    all = get_y_var_attr(od, j, OPERATE_BIT);
  }
  if (all) {
    free_y_var_array(od);
    update_field_object_attr(od, VERBOSE_BIT);
    result = calc_active_vars(od, FULL_MODEL);
  }
  else {
    for (i = -1; i < od->grid.object_num; ++i) {
      for (j = 0; j < od->y_vars; ++j) {
        if (!get_y_var_attr(od, j, OPERATE_BIT)) {
          if (i == -1) {
            fprintf(od->file[DEP_IN]->handle, "\t%s",
              od->cimal.y_var_name->me[j]);
          }
          else {
            fprintf(od->file[DEP_IN]->handle, "\t%lf",
              get_y_value(od, i, j, 0));
          }
        }
      }
      fprintf(od->file[DEP_IN]->handle, "\n");
    }
    rewind(od->file[DEP_IN]->handle);
    result = import_dependent(od, "all");
    if (od->file[DEP_IN]->handle) {
      fclose(od->file[DEP_IN]->handle);
      od->file[DEP_IN]->handle = NULL;
    }
  }
  
  return result;          
}
