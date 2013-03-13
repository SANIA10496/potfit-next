/****************************************************************
 *
 * calcpot_table.h:
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

#ifndef PTF_CALCPOT_TABLE_H
#define PTF_CALCPOT_TABLE_H

#include "table.h"

namespace POTFIT_NS {

  class CalcpotTable : protected Pointers {
  public:
    CalcpotTable(class POTFIT *);
    ~CalcpotTable();

  int   len;			/* total length of the table */
  int   idxlen;			/* number of changeable potential values */
  int   ncols;			/* number of columns */
  double *begin;		/* first value in the table */
  double *end;			/* last value in the table */
  double *step;			/* table increment */
  double *invstep;		/* inverse of increment */
  int  *first;			/* index of first entry */
  int  *last;			/* index of last entry */
  double *xcoord;		/* the x-coordinates of sampling points */
  double *table;		/* the actual data */
  double *d2tab;		/* second derivatives of table data for spline int */
  int  *idx;			/* indirect indexing */
  };
}

#endif /* PTF_CALCPOT_TABLE_H */
