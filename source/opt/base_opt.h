/****************************************************************
 *
 * base_opt.h:
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

#ifndef PTF_BASE_OPT_H
#define PTF_BASE_OPT_H

#include "../pointers.h"

namespace POTFIT_NS {

  class BaseOpt : protected Pointers {
  public:
    BaseOpt(class POTFIT *);
    ~BaseOpt();

    virtual void init(double *) = 0;
    virtual void run(void) = 0;

    double *params;
  };
}

#endif // PTF_BASE_OPT_H
