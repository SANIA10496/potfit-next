/****************************************************************
 *
 * force.cpp:
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

#include "interaction.h"
#include "force.h"
#include "forces/list_forces.h"

using namespace POTFIT_NS;

Interaction::Interaction(POTFIT *ptf) :
  Pointers(ptf),
  force(NULL)
{}

Interaction::~Interaction() {
  if (NULL != force)
    delete force;

  return;
}

void Interaction::init(void) {
  force = init_force(type);

  return;
}

double Interaction::calc_forces(void) {
  return force->calc_forces();
}

Force *Interaction::init_force(const std::string &force_type)
{
  if (force_type.compare("none") == 0)
    return NULL;
#define FORCE_TYPE
#define ForceType(key,Class) \
  else if (force_type.compare(#key) == 0) \
    return new Class(ptf);
#include "forces/list_forces.h"
#undef FORCE_TYPE
  return NULL;
}

void Interaction::set_type(const std::string &str) {
  type = str;

  return;
}

std::string Interaction::get_type(void) {
  return type;
}
