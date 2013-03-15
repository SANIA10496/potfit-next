/****************************************************************
 *
 * pointers.h: prototype for all classes with pointers
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

#ifndef PTF_POINTERS_H
#define PTF_POINTERS_H

#include "potfit.h"

namespace POTFIT_NS {

  class Pointers {
  public:
    Pointers(POTFIT *ptr) :
      ptf(ptr),
      io(ptr->io),
      input(ptr->input),
      output(ptr->output),
      interaction(ptr->interaction),
      optimization(ptr->optimization),
      random(ptr->random),
      structures(ptr->structures),
      settings(ptr->settings),
      potential(ptr->potential),
      memory(ptr->memory),
      communication(ptr->communication),
      utils(ptr->utils) {}
    virtual ~Pointers() {}

  protected:
    POTFIT *ptf;
    IO *&io;
    Input *&input;
    Output *&output;
    Interaction *&interaction;
    Optimization *&optimization;
    Random *&random;
    Structures *&structures;
    Settings *&settings;
    Potential *&potential;
    Memory *&memory;
    Communication *&communication;
    Utils *&utils;
  };

}

#endif /* PTF_POINTERS_H */
