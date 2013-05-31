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

    void set_myid(int);
    int get_myid(void);

    void set_num_cpus(int);
    int get_num_cpus(void);

    void set_extend(double);
    double get_extend(void);

    void set_opt(int);
    int get_opt(void);

    void set_apot_punish(double);
    double get_apot_punish(void);

    void set_eweight(double);
    double get_eweight(void);

    void set_sweight(double);
    double get_sweight(void);

    void set_deps(double);
    double get_deps(void);

    void set_dp_cut(double);
    double get_dp_cut(void);

    void set_dp_tol(double);
    double get_dp_tol(void);

    void set_dp_mix(double);
    double get_dp_mix(void);

  private:
    int myid;
    int num_cpus;
    int opt;

    double eweight;
    double sweight;
    double extend;
    double evo_threshold;
    double apot_punish_value;
    double d_eps;
    double dp_cut;
    double dp_tol;
    double dp_mix;
  };
}

#endif /* PTF_SETTINGS_H */
