/****************************************************************
 *
 * potential.cpp:
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

#include <mpi.h>

#include "force.h"
#include "interaction.h"
#include "memory.h"
#include "potential.h"

using namespace POTFIT_NS;

Potential::Potential(POTFIT *ptf) : Pointers(ptf) {
  enable_cp = 0;
  format = 0;
  have_grad = 0;
  n_invar_pots = 0;

  gradient = NULL;
  invar_pot = NULL;
}

Potential::~Potential() {
  delete [] gradient;
  delete [] invar_pot;
}

void Potential::init() {
  memory->create(gradient,interaction->force->cols,"gradient");
  memory->create(invar_pot,interaction->force->cols,"invar_pot");
}
