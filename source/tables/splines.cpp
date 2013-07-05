/****************************************************************
 *
 * spline.cpp:
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

#include <cstdlib>
#include <cstdio>

#include "splines.h"

using namespace POTFIT_NS;

Splines::Splines(POTFIT *ptf) :
  Pointers(ptf),
  a(0.0),
  p1(0.0),
  p2(0.0),
  d21(0.0),
  d22(0.0)
{}

Splines::~Splines(void) {
  return;
}

/****************************************************************
 *
 * spline_ed: initializes second derivatives used for spline interpolation
 *            (equidistant x[i])
 *
 ****************************************************************/

void Splines::spline_ed(const double &xstep, double *y, const int &n, const double &yp1, const double &ypn, double *y2)
{
  int   i, k;
  double p, qn, un;

  if (n > u.size())
    u.resize(n, 0.0);

  if (yp1 > 0.99e30)
    y2[0] = u[0] = 0.0;
  else {
    y2[0] = -0.5;
    u[0] = (3.0 / (xstep)) * ((y[1] - y[0]) / (xstep) - yp1);
  }
  for (i = 1; i < n - 1; i++) {
    p = 0.5 * y2[i - 1] + 2.0;
    y2[i] = (-0.5) / p;
    u[i] = (y[i + 1] - y[i]) / xstep - (y[i] - y[i - 1]) / (xstep);
    u[i] = (6.0 * u[i] / (2 * xstep) - 0.5 * u[i - 1]) / p;
  }
  if (ypn > 0.99e30)
    qn = un = 0.0;
  else {
    qn = 0.5;
    un = (3.0 / (xstep)) * (ypn - (y[n - 1] - y[n - 2]) / (xstep));
  }
  y2[n - 1] = (un - qn * u[n - 2]) / (qn * y2[n - 2] + 1.0);
  for (k = n - 2; k >= 0; k--)
    y2[k] = y2[k] * y2[k + 1] + u[k];
}

/****************************************************************
 *
 * spline_ne  : initializes second derivatives used for spline interpolation
 *            (nonequidistant x[i])
 *
 ****************************************************************/

void Splines::spline_ne(double *x, double *y, const int &n, const double &yp1, const double &ypn, double *y2)
{
  int   i, k;
  double p, qn, sig, un;

  if (n > u.size())
    u.resize(n, 0.0);

  if (yp1 > 0.99e30)
    y2[0] = u[0] = 0.0;
  else {
    y2[0] = -0.5;
    u[0] = (3.0 / (x[1] - x[0])) * ((y[1] - y[0]) / (x[1] - x[0]) - yp1);
  }
  for (i = 1; i < n - 1; i++) {
    sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]);
    p = sig * y2[i - 1] + 2.0;
    y2[i] = (sig - 1.0) / p;
    u[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]) - (y[i] - y[i - 1]) / (x[i] - x[i - 1]);
    u[i] = (6.0 * u[i] / (x[i + 1] - x[i - 1]) - sig * u[i - 1]) / p;
  }
  if (ypn > 0.99e30)
    qn = un = 0.0;
  else {
    qn = 0.5;
    un = (3.0 / (x[n - 1] - x[n - 2])) * (ypn - (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]));
  }
  y2[n - 1] = (un - qn * u[n - 2]) / (qn * y2[n - 2] + 1.0);
  for (k = n - 2; k >= 0; k--)
    y2[k] = y2[k] * y2[k + 1] + u[k];
}

//double Splines::splint_dir(double *xi, double *d2tab, const int &k, const double &b, const double &step)
//{
//  a = 1.0 - b;
//  p1 = *(xi + k);
//  d21 = *(d2tab + k);
//  p2 = *(xi + k + 1);
//  d22 = *(d2tab + k + 1);

//  return a * p1 + b * p2 + ((a * a * a - a) * d21 + (b * b * b - b) * d22) * (step * step) / 6.0;
//}

//double Splines::splint_comb_dir(double *xi, double *d2tab, const int &k, const double &b, const double &step, double *grad)
//{
//  a = 1.0 - b;
//  p1 = xi[k];
//  d21 = d2tab[k];
//  p2 = xi[k+1];
//  d22 = d2tab[k+1];

//  *grad = (p2 - p1) / step + ((3 * (b * b) - 1) * d22 - (3 * (a * a) - 1) * d21) * step / 6.0;

//  return a * p1 + b * p2 + ((a * a * a - a) * d21 + (b * b * b - b) * d22) * (step * step) / 6.0;
//}
