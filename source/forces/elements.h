/****************************************************************
 *
 * elements.h:
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

#ifndef PTF_ELEMENTS_H
#define PTF_ELEMENTS_H

#include <string>
#include <vector>

namespace POTFIT_NS {

  class Elements {
  public:
    Elements();
    ~Elements();

    double mass_from_number(const int &);
    double mass_from_name(const std::string &);
    int number_from_name(const std::string &);
  private:
    std::vector<std::string> name;
    std::vector<std::string> short_name;
    std::vector<double> mass;
  };
}

#endif /* PTF_ELEMENTS_H */
