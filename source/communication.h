/****************************************************************
 *
 * communication.h:
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

#ifndef PTF_COMMUNICATION_H
#define PTF_COMMUNICATION_H

#include <mpi.h>

#include "pointers.h"

namespace POTFIT_NS {

  class Communication : protected Pointers {
  public:
    Communication(class POTFIT *);
    ~Communication();

    void init(void);
    void broadcast_params(void);
    double gather_forces(const double &);
    void set_config_per_cpu(void);

    MPI::Datatype MPI_VECTOR;
    MPI::Datatype MPI_STENS;

    int exit_flag;
  private:
    int opt_pot_len;
  };

}

#endif /* PTF_COMMUNICATION_H */
