/****************************************************************
 *
 * opt_simann.cpp:
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

#include "opt_simann.h"

#include "../interaction.h"
#include "../io.h"
#include "../output.h"
#include "../potential.h"
#include "../random.h"
#include "../settings.h"
#include "../utils.h"

#define EPS 0.1
#define NEPS 4
#define NSTEP 20
#define NTEMP (3*ndim)
#define STEPVAR 2.0
#define TEMPVAR 0.85
#define KMAX 1000

using namespace POTFIT_NS;

OptSimann::OptSimann(POTFIT *ptf) :
  BaseOpt(ptf),
  auto_T(0),
  T(-1.0)
{}

OptSimann::~OptSimann() {}

void OptSimann::init(std::vector<std::string> &p) {
  /* check for automatic temperature */
  if (atof(p[0].c_str()) < 0) {
    auto_T = 1;
  } else if (atof(p[0].c_str()) > 0) {
    T = atof(p[0].c_str());
  } else {
    io->error << "The starting temperature for Simulated Annealing cannot be 0!" << std::endl;
    io->pexit(EXIT_FAILURE);
  }

  return;
}

void OptSimann::run(void) {
  int ndim = potential->get_num_free_params();
  double *Fvar;			/* backlog of Fn vals */
  double *v;			/* step vector */
  double *xi = potential->opt->values;
  double *xopt, *xi2;		/* optimal value */
  int  *naccept;		/* number of accepted changes in dir */
  double F, Fopt, F2;		/* Fn value */
  int   h = 0, j = 0, k = 0, n = 0, m = 0;	/* counters */
  int   loopagain;		/* loop flag */
//#ifndef APOT
//  double width, height;		/* gaussian bump size */
//#endif /* APOT */
//#ifdef DIPOLE
//  FILE *outfile;
//  char *filename = "Dipole.convergency";
//#endif /* DIPOLE */

  if (T == 0. && auto_T != 1)
    return;			/* don't anneal if starttemp equal zero */

  Fvar = new double[KMAX + 5 + NEPS];	/* Backlog of old F values */
  v = new double[ndim];
  xopt = new double[ndim];
  xi2 = new double[ndim];
  naccept = new int[ndim];

  /* init step vector and optimum vector */
  for (n = 0; n < ndim; n++) {
    v[n] = .1;
    naccept[n] = 0;
    xi2[n] = xi[n];
    xopt[n] = xi[n];
  }
  F = interaction->calc_forces();
  Fopt = F;

  /* determine optimum temperature for annealing */
  if (auto_T) {
    int   u = 10 * ndim;
    int   m1 = 0;
    double dF = 0.;
    double chi = .8;

    io->write << "Determining optimal starting temperature T ..." << std::endl;
    for (int e = 0; e < u; e++) {
      for (n = 0; n < ndim; n++)
        xi2[n] = xi[n];
      h = (int)(random->eqdist() * ndim);
      randomize_parameter(h, xi2, v);
      potential->opt->val_p = xi2;
      F2 = interaction->calc_forces();
      if (F2 <= F) {
        m1++;
      } else {
        dF += (F2 - F);
      }
    }
    io->write << "Did " << u << " steps, " << m1 << " were accepted" << std::endl;
    u -= m1;
    dF /= u;

    T = dF / log(u / (u * chi + (1 - chi) * m1));
    if (isnan(T) || isinf(T)) {
      io->error << "Simann failed because T was " << T << ", please set it manually." << std::endl;
      io->pexit(EXIT_FAILURE);
    }
    if (T < 0)
      T = -T;
    io->write << "Setting T=" << T << std::endl << std::endl;
  }

  io->write << "  k\t\tT\tm\t\tF\t\tFopt" << std::endl;
  io->write << std::setw(3) << 0 << "\t";
  io->write << std::setw(8) << T << "\t";
  io->write << std::setw(3) << 0 << "\t";
  io->write << std::setw(8) << F << "\t";
  io->write << std::setw(8) << Fopt << std::endl;

  for (n = 0; n <= NEPS; n++)
    Fvar[n] = F;

  /* annealing loop */
  do {
    for (m = 0; m < NTEMP; m++) {
      for (j = 0; j < NSTEP; j++) {
        for (h = 0; h < ndim; h++) {
          /* Step #1 */
          for (n = 0; n < ndim; n++) {
            xi2[n] = xi[n];
          }
          randomize_parameter(h, xi2, v);
          potential->opt->val_p = xi2;
          F2 = interaction->calc_forces();
          if (F2 <= F) {	/* accept new point */
            for (n = 0; n < ndim; n++)
              xi[n] = xi2[n];
            F = F2;
            naccept[h]++;
            if (F2 < Fopt) {
              for (n = 0; n < ndim; n++)
                xopt[n] = xi2[n];
              Fopt = F2;
              if (!output->get_tempfile().empty()) {
                potential->update_potentials(0);
                output->write_tempfile();
              }
            }
          } else if (random->eqdist() < (exp((F - F2) / T))) {
            for (n = 0; n < ndim; n++)
              xi[n] = xi2[n];
            F = F2;
            naccept[h]++;
          }
        }
      }

      /* Step adjustment */
      for (n = 0; n < ndim; n++) {
        if (naccept[n] > (0.6 * NSTEP))
          v[n] *= (1 + STEPVAR * ((double)naccept[n] / NSTEP - 0.6) / 0.4);
        else if (naccept[n] < (0.4 * NSTEP))
          v[n] /= (1 + STEPVAR * (0.4 - (double)naccept[n] / NSTEP) / 0.4);
        naccept[n] = 0;
      }

      io->write << std::setw(3) << k << "\t";
      io->write << std::setw(8) << T << "\t";
      io->write << std::setw(3) << m + 1 << "\t";
      io->write << std::setw(8) << F << "\t";
      io->write << std::setw(8) << Fopt << std::endl;

      /* End annealing if break flagfile exists */
      if (utils->check_for_flagfile() == 1) {
        io->write << "Annealing terminated in presence of break flagfile !" << std::endl;
        io->write << "Temperature was " << T << ", returning optimum configuration" << std::endl;
        for (n = 0; n < ndim; n++)
          xi[n] = xopt[n];
        F = Fopt;
        k = KMAX + 1;
        break;
      }

// TODO: rescaling
#if !defined APOT && ( defined EAM || defined ADP ) && !defined NORESCALE
//      /* Check for rescaling... every tenth step */
//      if ((m % 10) == 0) {
//        /* Was rescaling necessary ? */
//        printf("Force before rescaling %f\n", F);
//        if (rescale(&opt_pot, 1., 0) != 0.) {
//          /* wake other threads and sync potentials */
//          F = (*calc_forces) (xi, fxi1, 2);
//        }
//        printf("Force after rescaling %f\n", F);
//      }
#endif /* !APOT && ( EAM || ADP ) && !NORESCALE */
    }

    /*Temp adjustment */
    T *= TEMPVAR;
    k++;
    Fvar[k + NEPS] = F;
    loopagain = 0;
    for (n = 1; n <= NEPS; n++) {
      if (fabs(F - Fvar[k - n + NEPS]) > (EPS * F * 0.01))
        loopagain = 1;
    }
    if (!loopagain && ((F - Fopt) > (EPS * F * 0.01))) {
      for (n = 0; n < ndim; n++)
        xi[n] = xopt[n];
      F = Fopt;
      loopagain = 1;
    }
  } while (k < KMAX && loopagain);
  for (n = 0; n < ndim; n++) {
    xi[n] = xopt[n];
  }

  F = Fopt;
  if (!output->get_tempfile().empty()) {
    potential->update_potentials(0);
    output->write_tempfile();
  }

  delete [] Fvar;
  delete [] v;
  delete [] xopt;
  delete [] naccept;
  delete [] xi2;

  return;
}

void OptSimann::randomize_parameter(int n, double *xi, double *v)
{
  const int idxpot = potential->idxpot[n];
  const int idxpar = potential->idxparam[n];

  // analytic potential - only change single parameter
  if (potential->get_format(idxpot) == 0) {
    double temp, rand;
    int   done = 0, count = 0;
    double min, max;

    min = potential->opt->val_min[n];
    max = potential->opt->val_max[n];

    if (v[n] > max - min)
      v[n] = max - min;

    do {
      temp = xi[n];
      rand = 2.0 * random->eqdist() - 1.;
      temp += (rand * v[n]);
      if (temp >= min && temp <= max)
        done = 1;
      count++;
    } while (!done);
    xi[n] = temp;

    // tabulated potentials - change neighbors as well
  } else if (potential->pots[idxpot]->format == 3 || potential->pots[idxpot]->format == 4) {
    double width = fabs(random->normdist());
    double height = random->normdist() * v[n];

    for (int i = 0; i <= 4. * width; i++) {
      /* using idx avoids moving fixed points */
      if ((idxpar + i <= potential->pots[idxpot]->num_params) && ((n + i) < potential->get_num_free_params())) {
        xi[n + i] += gauss((double)i / width) * height;
      }
    }
    for (int i = 1; i <= 4. * width; i++) {
      if ((idxpar - i) >= 0) {
        xi[n - i] += gauss((double)i / width) * height;
      }
    }
  }

  return;
}
