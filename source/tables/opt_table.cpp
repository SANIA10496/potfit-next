/****************************************************************
 *
 * opt_table.cpp:
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

#include "opt_table.h"

using namespace POTFIT_NS;

OptTable::OptTable(POTFIT *ptf, int n) :
  Pointers(ptf),
  val_p(NULL) {
  values = new double[n];
  val_min = new double[n];
  val_max = new double[n];
  for (int i=0;i<n;i++) {
    values[i] = 0.0;
    val_min[i] = 0.0;
    val_max[i] = 0.0;
  }
  val_p = values;
}

OptTable::~OptTable() {
  if (NULL != values) {
    delete [] values;
    delete [] val_min;
    delete [] val_max;
  }
}
