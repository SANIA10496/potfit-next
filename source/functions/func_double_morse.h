/****************************************************************
 *
 * func_double_morse.h:
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

FunctionType(double_morse,FuncDOUBLE_MORSE)

#else

#ifndef PTF_FUNC_DOUBLE_MORSE_H
#define PTF_FUNC_DOUBLE_MORSE_H

#include "function.h"

namespace POTFIT_NS {

  class FuncDOUBLE_MORSE : public Function {
  public:
    FuncDOUBLE_MORSE();
    ~FuncDOUBLE_MORSE();

    int num_params(void);
    void calc(const double &, double *, double *);
  };
}

#endif // PTF_FUNC_DOUBLE_MORSE_H

#endif // FUNCTION_TYPE
