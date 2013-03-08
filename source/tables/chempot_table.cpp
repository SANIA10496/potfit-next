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
#include "../memory.h"
#include "../settings.h"

using namespace POTFIT_NS;

ChempotTable::ChempotTable(POTFIT *ptf, int ntypes) : Pointers(ptf) {
  number = ntypes;
  n_invar = 0;

  memory->create(values,number,"chemical potentials");
  memory->create(val_min,number,"chemical potentials");
  memory->create(val_max,number,"chemical potentials");
  memory->create(invar_par,number,"chemical potentials");

  param_name = (char **)malloc(number * sizeof(char *));
  for (int i=0; i<number; i++) {
    param_name[i] = (char *)malloc(30 * sizeof(char));
    strcpy(param_name[i],"\0");
  }
}

ChempotTable::~ChempotTable() {
}

void ChempotTable::add_value(int i, const char *name, double val, double min, double max) {
  char msg[255];

  if (min == max) {
    invar_par[number] = 1;
// what is this good for? global counter for invar_pars per potential^^
//        apt->invar_par[i][apt->globals]++;
  } else if (min > max) {
    double temp = min;
    min = max;
    max = temp;
  } else if ((val < min) || (val > max)) {
    if (settings->opt) {
      if (val < min)
        val = min;
      if (val > max)
        val = max;
      sprintf(msg, "Starting value for %s #%d is ", name, number + 1);
      sprintf(msg, "%soutside of specified adjustment range.\n",msg);
      io->warning("%sResetting it to %f.\n", msg, number + 1, val);
      if (val == 0)
        io->warning("New value is 0 ! Please be careful about this.\n");
    }
  }
  strcpy(param_name[number],name);
  values[number] = val;
  val_min[number] = min;
  val_max[number] = max;

  return;
}

void ChempotTable::set_value(int i, double val) {
  values[i] = val;
}
