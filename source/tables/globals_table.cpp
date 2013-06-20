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
#include "../potential.h"
#include "../settings.h"

using namespace POTFIT_NS;

GlobalsTable::GlobalsTable(POTFIT *ptf, int num) :
  Pointers(ptf),
  num_globals(num),
  num_free_globals(num),
  opt_pot_start(0)
{
  param_name.resize(num);
  values.resize(num,0.0);
  val_min.resize(num,0.0);
  val_max.resize(num,0.0);
  invar_par.resize(num,0);
  idx.resize(num,0);
  usage.resize(num,0);

  global_idx = new int**[num_globals];

  return;
}

GlobalsTable::~GlobalsTable() {
  param_name.clear();
  values.clear();
  val_min.clear();
  val_max.clear();
  invar_par.clear();
  idx.clear();
  usage.clear();

  return;
}

void GlobalsTable::add_param(const int &index, const char *name, const double &val, const double &min, const double &max) {
  // check for duplicate names
  for (int k = 0; k < index; k++) {
    if (param_name[k].compare(name) == 0) {
      io->error << "Found duplicate global parameter name!" << std::endl;
      io->error << "Parameter #" << index + 1 << " (" << name << ") is the same as #" << k + 1;
      io->error << " (" << param_name[k] << ")." << std::endl;
      io->pexit(EXIT_FAILURE);
    }
  }

  // add parameter if nothing went wrong
  param_name[index] = name;
  values[index] = val;
  val_min[index] = min;
  val_max[index] = max;

  // check for invariance and proper value (respect boundaries)
  // parameter will not be optimized if min==max
  if (val == min == max) {
    invar_par[index] = 1;
  } else if (min > max) {
    val_min[index] = max;
    val_max[index] = min;
  } else if ((val < min) || (val > max)) {
    // Only print warning if we are optimizing
    if (settings->get_opt()) {
      if (val < min) {
        values[index] = min;
      } else if (val > max) {
        values[index] = max;
      }
      io->warning << "Starting value for global parameter " << name << " (#" << index + 1 << ") is" << std::endl;
      io->warning << "outside of specified adjustment range." << std::endl;
      io->warning << "Resetting it to " << val << "." << std::endl;
      if (val == 0)
        io->warning << "New value is 0 ! Please be careful about this." << std::endl;
    }
  }

  return;
}

void GlobalsTable::check_usage(void) {
  for (int i=0; i<num_globals; i++)
    if (usage[i] == 0)
      num_free_globals--;

  return;
}

int GlobalsTable::get_number_params(void) {
  return num_globals;
}

int GlobalsTable::get_number_free_params(void) {
  return num_free_globals;
}


void GlobalsTable::set_value(const int &index, const double &val) {
  values[index] = val;

  return;
}

void GlobalsTable::set_opt_pot_start(const int &val) {
  opt_pot_start = val;

  return;
}

int GlobalsTable::get_index(const char *name) {
  for (int i=0;i<num_globals;i++)
    if (param_name[i].compare(name)==0)
      return i;

  return -1;
}

void GlobalsTable::add_usage(const int &global_index, const int &pot_index, const int &param_index) {
  if (++usage[global_index] > 1) {
    global_idx[global_index] = (int **)realloc(global_idx[global_index], usage[global_index] * sizeof(int *));
  } else {
    global_idx[global_index] = (int **)malloc(1 * sizeof(int *));
  }
  global_idx[global_index][usage[global_index] - 1] = (int *)malloc(2 * sizeof(int));
  global_idx[global_index][usage[global_index] - 1][0] = pot_index;
  global_idx[global_index][usage[global_index] - 1][1] = param_index;

  return;
}

void GlobalsTable::get_value(const int &index, double *val) {
  val[0] = values[index];
  val[1] = val_min[index];
  val[2] = val_max[index];

  return;
}

void GlobalsTable::update_potentials(void) {
  for (int i=0;i<num_free_globals;i++) {
    values[i] = potential->opt->val_p[opt_pot_start + i];
    for (int j=0;j<usage[i];j++) {
      potential->pots[global_idx[i][j][0]]->set_param(global_idx[i][j][1],values[i]);
    }
  }

  return;
}

void GlobalsTable::write_potential(std::ofstream &outfile) {
  outfile << std::endl << "globals" << "\t" << num_globals << std::endl;
  for (int i=0;i<num_globals;i++) {
    outfile << param_name[i] << "\t" << values[i] << "\t" << val_min[i] << "\t" << val_max[i] << std::endl;
  }

  return;
}

