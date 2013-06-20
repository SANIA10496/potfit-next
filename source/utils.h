/****************************************************************
 *
 * utils.h:
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

#ifndef PTF_UTILS_H
#define PTF_UTILS_H

#include <ctime>
#include <string>

#include "pointers.h"
#include "types.h"

namespace POTFIT_NS {

  class Utils : protected Pointers {
  public:
    Utils(class POTFIT *);
    ~Utils();

    char *tolowercase(char *);
    double vect_dist(const vector &, const vector &);

    void start_timer(void);
    void end_timer(void);
    const int timediff(void);

    void set_flagfile(const std::string &);
    const int check_for_flagfile(void);

  private:
    time_t t_begin;
    time_t t_end;

    std::string flagfile;
  };

  void power_1(double &, const double &, const double &);
  void power_m(const int &, double *, const double *, const double *);
}

#endif /* PTF_UTILS_H */
