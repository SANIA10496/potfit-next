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

Interaction::Interaction(POTFIT *ptf) : Pointers(ptf)
{
  strcpy(type,"\0");
  force = NULL;
}

Interaction::~Interaction()
{
  if (force) delete force;
}

void Interaction::init()
{
  force = init_force(type);
}

Force *Interaction::init_force(const char *force_type)
{
  if (strcmp(force_type,"none") == 0)
    return NULL;
#define FORCE_TYPE
#define ForceType(key,Class) \
  else if (strcmp(force_type,#key) == 0) \
    return new Class(ptf);
#include "forces/list_forces.h"
#undef FORCE_TYPE
}
