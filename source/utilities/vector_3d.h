/****************************************************************
 *
 * vector_3d.h
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

#ifndef PTF_VECTOR_3D_H
#define PTF_VECTOR_3D_H

namespace POTFIT_NS {

  template <typename T>
  class Vector3D {
  public:
    Vector3D() :
      x(0), y(0), z(0) {};
    Vector3D(T &a, T &b, T &c) :
      x(a), y(b), z(c) {};

  private:
    T x, y, z;
  };
}

#endif // PTF_VECTOR_3D_H
