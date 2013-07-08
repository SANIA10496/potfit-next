/****************************************************************
 *
 * func_sheng_phi2.h:
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

#ifdef FUNCTION_TYPE

FunctionType(sheng_phi2,FuncSHENG_PHI2)

#else

#ifndef PTF_FUNC_SHENG_PHI2_H
#define PTF_FUNC_SHENG_PHI2_H

#include "function.h"

namespace POTFIT_NS {

  class FuncSHENG_PHI2 : public Function {
  public:
    FuncSHENG_PHI2();
    ~FuncSHENG_PHI2();

    int num_params(void);
    void calc(const double &, const std::vector<double> &, double *);
  };
}

#endif // PTF_FUNC_SHENG_PHI2_H
#endif // FUNCTION_TYPE
