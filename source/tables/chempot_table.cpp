/****************************************************************
 *
 * chempot_table.cpp:
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

#include "chempot_table.h"

#include "../io.h"
#include "../potential.h"
#include "../settings.h"

using namespace POTFIT_NS;

ChempotTable::ChempotTable(POTFIT *ptf, int ntypes) : Pointers(ptf) {
  num_params = ntypes;
  num_invar_params = 0;

  values = new double[num_params];
  val_min = new double[num_params];
  val_max = new double[num_params];
  invar_par = new int[num_params];
  param_name = new char*[num_params];

  for (int i=0;i<num_params;i++) {
    values[i] = 0.0;
    val_min[i] = 0.0;
    val_max[i] = 0.0;
    invar_par[i] = 0;
    param_name[i] = new char[30];
    strcpy(param_name[i],"\0");
  }

  return;
}

ChempotTable::~ChempotTable() {
  delete [] values;
  delete [] val_min;
  delete [] val_max;
  delete [] invar_par;
  for (int i=0;i<num_params;i++)
    delete [] param_name[i];
  delete [] param_name;
}

void ChempotTable::add_value(int i, const char *name, double val, double min, double max) {
  char msg[255];

  if (min == max) {
    invar_par[i] = 1;
  } else if (min > max) {
    double temp = min;
    min = max;
    max = temp;
  } else if ((val < min) || (val > max)) {
    if (settings->get_opt()) {
      if (val < min)
        val = min;
      if (val > max)
        val = max;
      io->warning << "Starting value for " << name << " is outside of specified adjustment range." << std::endl;
      io->warning << "Resetting it to " << val << "." << std::endl;
      if (val == 0)
        io->warning << "New value is 0 ! Please be careful about this." << std::endl;
    }
  }
  strcpy(param_name[i],name);
  values[i] = val;
  val_min[i] = min;
  val_max[i] = max;

  potential->num_params++;
  if (invar_par[i] == 0)
    potential->num_free_params++;

  return;
}

void ChempotTable::set_value(int i, double val) {
  values[i] = val;
}
