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

#include <ostream>

namespace POTFIT_NS {

  template <typename T>
  class Vector3D {
  public:
    Vector3D() :
      _x(0), _y(0), _z(0) {};
    Vector3D(T& a, T& b, T& c) :
      _x(a), _y(b), _z(c) {};
    Vector3D(T a, T b, T c) :
      _x(a), _y(b), _z(c) {};
    Vector3D(Vector3D<T>& v) :
      _x(0), _y(0), _z(0) {
        _x = v.x();
        _y = v.y();
        _z = v.z();
      };

    void assign(T& a, T& b, T& c) {
      _x = a; _y = b; _z = c;
      return;
    }

    const T& x(void) {
      return _x;
    }

    const T& y(void) {
      return _y;
    }

    const T& z(void) {
      return _z;
    };

    std::ostream& OutputToStream(std::ostream& lhs) const {
      lhs << "(" << _x << "," << _y << "," << _z << ")";
      return lhs;
    };

  private:
    T _x, _y, _z;
  };

  template <typename T>
  const Vector3D<T>& operator+(const Vector3D<T>& a, const Vector3D<T>& b) {
    Vector3D<T> c;
  };

  template <typename T>
  const Vector3D<T>& operator+(const T& a, const Vector3D<T>& b) {};

  template <typename T>
  const Vector3D<T>& operator+(const Vector3D<T>& a, const T& b) {};

  template <typename T>
  std::ostream& operator<<(std::ostream& lhs, const Vector3D<T>& rhs) {
    return rhs.OutputToStream(lhs);
  };
}

#endif // PTF_VECTOR_3D_H
