/****************************************************************
 *
 * io.h: misc io routines
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

#ifndef PTF_INPUT_H
#define PTF_INPUT_H

#include <iostream>

#include "pointers.h"
#include "types.h"

namespace POTFIT_NS {

  class Input : protected Pointers {
  friend class IO;
  public:
    Input(class POTFIT *, int, char **);
    ~Input();

    void read_parameter_file();
    void read_potential_file();
    void read_config_file();

  private:
    int get_param(const char *, void *, param_t, int, int);
    void check_params();

    char param_file[255]; 	// parameter file
    char config_file[255];
    char maxchfile[255];
    char startpot[255];

    int curline;
    int enable_maxch_file;
  };

}

#endif /* PTF_INPUT_H */
