/****************************************************************
 *
 * settings.h:
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

#ifndef PTF_SETTINGS_H
#define PTF_SETTINGS_H

#include <cstdio>

#include "pointers.h"

namespace POTFIT_NS {

  class Settings : protected Pointers {
  public:
    Settings(class POTFIT *);
    ~Settings();

    char interaction[255];

    int myid;
    int num_cpus;
    int opt;

    double eweight;
    double sweight;
    double extend;
    double evo_threshold;
    double anneal_temp;
    double apot_punish_value;
    double d_eps;
    double dp_cut;
    double dp_tol;
    double dp_mix;

  private:
  };

}

#endif /* PTF_SETTINGS_H */
