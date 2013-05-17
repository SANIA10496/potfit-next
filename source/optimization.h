/****************************************************************
 *
 * optimization.h:
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

#ifndef PTF_OPTIMIZATION_H
#define PTF_OPTIMIZATION_H

//#include <string>
#include <vector>

#include "pointers.h"
#include "opt/base_opt.h"

namespace POTFIT_NS {

  class Optimization : protected Pointers {
  public:
    Optimization(class POTFIT *);
    ~Optimization();

    void run(void);

    void add_algorithm(void);
    BaseOpt *init_algorithm(const char *);

    std::vector<std::string> algorithms;
    std::vector<double *> params;
    int num_algs;

  private:
    BaseOpt *opt;
  };
}

#endif // PTF_OPTIMIZATION_H
