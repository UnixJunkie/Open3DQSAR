/*

qsort_functions.c

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


int compare_integers(const void *a, const void *b)
{
  const int *da = (const int *)a;
  const int *db = (const int *)b;


  return (*da > *db) - (*da < *db);
}


int compare_bond_info(const void *a, const void *b)
{
  const BondInfo *da = (const BondInfo *)a;
  const BondInfo *db = (const BondInfo *)b;


  return (da->num > db->num) - (da->num < db->num);
}


int compare_seed_dist(const void *a, const void *b)
{
  int result;
  int result2;
  const SeedDistMat **da = (const SeedDistMat **)a;
  const SeedDistMat **db = (const SeedDistMat **)b;


  result = (((*da)->dist - (*db)->dist) > MSD_THRESHOLD)
    - (((*da)->dist - (*db)->dist) < -MSD_THRESHOLD);
  return (result ? result
    : ((result2 = (((*da)->seed[0] > (*db)->seed[0])
    - ((*da)->seed[0] < (*db)->seed[0]))) ? result2
    : ((*da)->seed[1] > (*db)->seed[1])
    - ((*da)->seed[1] < (*db)->seed[1])));
}


int compare_regex_data(const void *a, const void *b)
{
  int result;
  const RegexData **da = (const RegexData **)a;
  const RegexData **db = (const RegexData **)b;


  result = ((*da)->struct_num > (*db)->struct_num)
    - ((*da)->struct_num < (*db)->struct_num);
  return (result ? result
    : (((*da)->conf_num > (*db)->conf_num)
    - ((*da)->conf_num < (*db)->conf_num)));
}


int compare_dist(const void *a, const void *b)
{
  int result;
  const AtomPair *da = (const AtomPair *)a;
  const AtomPair *db = (const AtomPair *)b;


  result = ((da->dist - db->dist) > MSD_THRESHOLD)
    - ((da->dist - db->dist) < -MSD_THRESHOLD);
  return (result ? result : ((da->a[0] == db->a[0])
    ? (da->a[1] > db->a[1]) - (da->a[1] < db->a[1])
    : (da->a[0] > db->a[0]) - (da->a[0] < db->a[0])));
}


int compare_score(const void *a, const void *b)
{
  int result;
  const AtomPair *da = (const AtomPair *)a;
  const AtomPair *db = (const AtomPair *)b;


  result = (da->score > db->score) - (da->score < db->score);
  return (result ? result : (da->cost > db->cost) - (da->cost < db->cost));
}


int compare_conf_energy(const void *a, const void *b)
{
  int result;
  const ConfInfo **da = (const ConfInfo **)a;
  const ConfInfo **db = (const ConfInfo **)b;


  result = (((*da)->energy - (*db)->energy) > ENERGY_THRESHOLD)
    - (((*da)->energy - (*db)->energy) < -ENERGY_THRESHOLD);
  return (result ? result : ((*da)->n_conf > (*db)->n_conf)
    - ((*da)->n_conf < (*db)->n_conf));
}


int compare_n_phar_points(const void *a, const void *b)
{
  int result;
  int result2;
  const PharConfInfo **da = (const PharConfInfo **)a;
  const PharConfInfo **db = (const PharConfInfo **)b;


  result = ((*da)->n_phar_points < (*db)->n_phar_points)
    - ((*da)->n_phar_points > (*db)->n_phar_points);
  return (result ? result
    : ((result2 = (((*da)->object_num > (*db)->object_num)
    - ((*da)->object_num < (*db)->object_num))) ? result2
    : ((*da)->conf_num > (*db)->conf_num)
    - ((*da)->conf_num < (*db)->conf_num)));
}


int compare_template_score(const void *a, const void *b)
{
  int result;
  const TemplateInfo **da = (const TemplateInfo **)a;
  const TemplateInfo **db = (const TemplateInfo **)b;


  result = ((*da)->score < (*db)->score)
    - ((*da)->score > (*db)->score);
  return (result ? result : ((*da)->num > (*db)->num)
    - ((*da)->num < (*db)->num));
}
