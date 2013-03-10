/****************************************************************
 *
 * globals_table.cpp:
 *
 ****************************************************************
 *
 * Copyright 2002-2013
 *	Institute for Theoretical and Applied Physics
 *	University of Stuttgart, D-70550 Stuttgart, Germany
 *	http://potfit.itap.physik.uni-stuttgart.de/
 *
 ****************************************************************
 *
 *   This file is part of potfit.
 *
 *   potfit is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   potfit is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with potfit; if not, see <http://www.gnu.org/licenses/>.
 *
 ****************************************************************/

#include <cstdlib>
#include <cstring>

#include "globals_table.h"

#include "../io.h"
#include "../memory.h"
#include "../settings.h"

using namespace POTFIT_NS;

GlobalsTable::GlobalsTable(POTFIT *ptf, int num) : Pointers(ptf) {
  num_globals = num;
  num_free_params = 0;

  param_name = (char **)malloc(num_globals * sizeof(char *));
  if (NULL == param_name)
    io->error("Could not allocate memory for potential name\n");
  for (int i=0; i<num_globals; i++) {
    param_name[i] = (char *)malloc(20 * sizeof(char));
    if (NULL == param_name[i])
      io->error("Could not allocate memory for parameter names\n");
    strcpy(param_name[i],"\0");
  }

  memory->create(values, num_globals, "global parameter values");
  memory->create(val_min, num_globals, "global parameter minimum values");
  memory->create(val_max, num_globals, "global parameter maximum values");
  memory->create(invar_par, num_globals, "global parameter invariant setting");
  memory->create(idx, num_globals, "global parameter indirect index");
  memory->create(usage, num_globals, "global parameter indirect index");

  global_idx = (int ***)malloc(num_globals * sizeof(int **));

  for (int i=0; i<num_globals; i++) {
    values[i] = 0.0;
    val_min[i] = 0.0;
    val_max[i] = 0.0;
    invar_par[i] = 0;
    idx[i] = 0;
    usage[i] = 0;
  }
}

GlobalsTable::~GlobalsTable() {
  memory->destroy(values);
  memory->destroy(val_min);
  memory->destroy(val_max);
  memory->destroy(invar_par);
  memory->destroy(usage);
}

void GlobalsTable::add_param(int index, const char *name, double val, double min, double max) {
  char msg[255];

  // check for duplicate names
  for (int k = 0; k < index; k++) {
    if (strcmp(name, param_name[k]) == 0) {
      sprintf(msg, "Found duplicate global parameter name!\n");
      io->error("%sParameter #%d (%s) is the same as #%d (%s)\n", msg, index + 1, name, k + 1, param_name[k]);
    }
  }

  // check for invariance and proper value (respect boundaries)
  // parameter will not be optimized if min==max */
  if (min == max) {
    invar_par[index] = 1;
  } else if (min > max) {
    double temp = min;
    min = max;
    max = temp;
  } else if ((val < min) || (val > max)) {
    /* Only print warning if we are optimizing */
    if (settings->opt) {
      if (val < min)
        val = min;
      if (val > max)
        val = max;
      sprintf(msg, "Starting value for global parameter %s (#%d) is\n", name, index + 1);
      sprintf(msg, "%soutside of specified adjustment range.\n",msg);
      io->warning("%sResetting it to %f.\n", msg, index + 1, val);
      if (val == 0)
        io->warning("New value is 0 ! Please be careful about this.\n");
    }
  }

  // add parameter if nothing went wrong
  strcpy(param_name[index],name);
  values[index] = val;
  val_min[index] = min;
  val_max[index] = max;

  if (!invar_par[index])
    idx[num_free_params++] = index;

  return;
}

void GlobalsTable::set_value(int index, double val) {
  values[index] = val;
}

int GlobalsTable::get_index(const char *name) {
  for (int i=0;i<num_globals;i++)
    if (strcmp(name, param_name[i])==0)
      return i;

  return -1;
}

void GlobalsTable::add_usage(int global_index, int pot_index, int param_index) {
  if (++usage[global_index] > 1) {
    global_idx[global_index] = (int **)realloc(global_idx[global_index], usage[global_index] * sizeof(int *));
  } else {
    global_idx[global_index] = (int **)malloc(1 * sizeof(int *));
  }
  global_idx[global_index][usage[global_index] - 1] = (int *)malloc(2 * sizeof(int));
  global_idx[global_index][usage[global_index] - 1][0] = pot_index;
  global_idx[global_index][usage[global_index] - 1][1] = param_index;
}

void GlobalsTable::get_value(int index, double *val) {
  val[0] = values[index];
  val[1] = val_min[index];
  val[2] = val_max[index];
}

void GlobalsTable::get_values(int *number, double *values) {
  *number = num_globals;
}