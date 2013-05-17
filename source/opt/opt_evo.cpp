/****************************************************************
 *
 * opt_evo.cpp:
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

#include <cmath>
#include <cstdlib>
#include <limits>

#include "opt_evo.h"

#include "../interaction.h"
#include "../io.h"
#include "../output.h"
#include "../potential.h"
#include "../random.h"
#include "../settings.h"
#include "../utils.h"

#define D (ndim+2) 		// all free parameters plus F and CR
#define NP 15*D			/* number of total population */

#define JR 0.6			/* jumping rate for opposite algorithm */

/* boundary values for self-adapting parameters */
#define F_LOWER 0.1		/* lower value for F */
#define F_UPPER 0.9		/* upper value for F */
#define TAU_1 0.1		/* probability for changing F */
#define TAU_2 0.1		/* probability for changing CR */

using namespace POTFIT_NS;

OptEvo::OptEvo(POTFIT *ptf) : BaseOpt(ptf) {
  tot_cost = NULL;
  tot_P = NULL;

  return;
}

OptEvo::~OptEvo() {
  int ndim = potential->num_free_params;

  for (int i=0; i<2*NP; i++)
    delete [] tot_P[i];
  delete [] tot_P;
  delete [] tot_cost;

  return;
}

void OptEvo::init(double *p) {
  params = p;

  return;
}

void OptEvo::run(void) {
  int   a, b, c;		/* store randomly picked numbers */
//  int   d, e; 			/* enable this line for more vectors */
  int   i, j, k;		/* counters */
  int   count = 0;		/* counter for loops */
  int maxsteps = 1e8;
  int ndim = potential->num_free_params;
  int   jsteps = 0;
  double avg = 0.;		/* average sum of squares for all configurations */
  double crit = 0.;		/* treshold for stopping criterion */
  double evo_threshold = 1e-6;
  double force = 0.;		/* holds the current sum of squares */
  double jumprate = JR;
  double min = 1e10;		/* current minimum for all configurations */
  double max = 0.;		/* current maximum for all configurations */
  double temp = 0.;		/* temp storage */
  double pmin = 0.;		/* lower bound for parameter */
  double pmax = 0.;		/* upper bound for parameter */
  double *best;			/* best configuration */
  double *cost;			/* cost values for all configurations */
  double *trial;		/* current trial configuration */
  double *xi = potential->opt->values;
  double **x1;			/* current population */
  double **x2;			/* next generation */
  FILE *ff;			/* exit flagfile */

  // parameters for evo are:
  // 0: epsilon for error margin
  // 1: max. number of steps
  // 2: not used
  if (params[0] <= 0.) {
    io->write("Skipping evo algorithm with error margin of %f\n",params[0]);
    return;
  }
  evo_threshold = params[0];
  if (params[1] > 0) {
    maxsteps = params[1];
  }

  /* vector with new configuration */
  trial = new double[D];

  /* allocate memory for all configurations */
  x1 = new double*[NP];
  x2 = new double*[NP];
  best = new double[NP];
  cost = new double[NP];
  if (x1 == NULL || x2 == NULL || trial == NULL || cost == NULL || best == NULL)
    io->error("Could not allocate memory for population vector!\n");
  for (i = 0; i < NP; i++) {
    x1[i] = new double[D];
    x2[i] = new double[D];
    if (x1[i] == NULL || x2[i] == NULL)
      io->error("Could not allocate memory for population vector!\n");
    for (j = 0; j < D; j++) {
      x1[i][j] = 0;
      x2[i][j] = 0;
    }
  }

  io->write("Initializing population ... ");
  fflush(stdout);

  init_population(x1, xi, cost);
  for (i = 0; i < NP; i++) {
    if (cost[i] < min) {
      min = cost[i];
      for (j = 0; j < D; j++)
        best[j] = x1[i][j];
    }
    if (cost[i] > max)
      max = cost[i];
  }
  for (i = 0; i < NP; i++)
    avg += cost[i];
  io->write("done\n");

  crit = max - min;

  io->write("Loops\t\tOptimum\t\tAverage error sum\t\tMax-Min\n");
  io->write("%5d\t\t%15f\t%20f\t\t%.2e\n", count, min, avg / (NP), crit);
  fflush(stdout);

  /* main differential evolution loop */
  while (crit >= evo_threshold && min >= evo_threshold && count < maxsteps) {
    max = 0.;
    /* randomly create new populations */
    for (i = 0; i < NP; i++) {
      /* generate random numbers */
      do
        a = (int)floor(random->eqdist() * NP);
      while (a == i);
      do
        b = (int)floor(random->eqdist() * NP);
      while (b == i || b == a);
      /*      do*/
      /*        c = (int)floor(random->eqdist() * NP);*/
      /*      while (c == i || c == a || c == b);*/
      /*      do*/
      /*        d = (int)floor(random->eqdist() * NP);*/
      /*      while (d == i || d == a || d == b || d == c);*/
      /*      do*/
      /*        e = (int)floor(random->eqdist() * NP);*/
      /*      while (e == i || e == a || e == b || e == c || e == d);*/

      j = (int)floor(random->eqdist() * ndim);

      /* self-adaptive parameters */
      if (random->eqdist() < TAU_1)
        trial[D - 2] = F_LOWER + random->eqdist() * F_UPPER;
      else
        trial[D - 2] = x1[i][D - 2];
      if (random->eqdist() < TAU_2)
        trial[D - 1] = random->eqdist();
      else
        trial[D - 1] = x1[i][D - 1];

      /* create trail vectors with different methods */
      for (k = 1; k <= ndim; k++) {
        if (random->eqdist() < trial[D - 1] || k == j) {
          /* DE/rand/1/exp */
          /*          temp = x1[c][idx[j]] + trial[D - 2] * (x1[a][idx[j]] - x1[b][idx[j]]);*/
          /* DE/best/1/exp */
          temp = best[j] + trial[D - 2] * (x1[a][j] - x1[b][j]);
          /* DE/rand/2/exp */
          /*          temp = x1[e][j] + trial[D-2] * (x1[a][j] + x1[b][j] - x1[c][j] - x1[d][j]);*/
          /* DE/best/2/exp */
          /*          temp = best[j] + trial[D-2] * (x1[a][j] + x1[b][j] - x1[c][j] - x1[d][j]);*/
          /* DE/rand-to-best/1/exp */
          /*          temp = x1[c][j] + (1 - trial[D-2]) * (best[j] - x1[c][j]) +*/
          /*            trial[D-2] * (x1[a][j] - x1[b][j]);*/
          /* DE/rand-to-best/2/exp */
          /*          temp = x1[e][j] + (1 - trial[D-2]) * (best[j] - x1[e][j]) +*/
          /*            trial[D-2] * (x1[a][j] + x1[b][j] - x1[c][j] - x1[d][j]);*/
          pmin = potential->pots[potential->idxpot[j]]->get_val_min(potential->idxparam[j]);
          pmax = potential->pots[potential->idxpot[j]]->get_val_max(potential->idxparam[j]);
          if (temp > pmax) {
            trial[j] = pmax;
          } else if (temp < pmin) {
            trial[j] = pmin;
          } else
            trial[j] = temp;
        } else {
          trial[j] = x1[i][j];
        }
        j = (j + 1) % ndim;
      }

      potential->opt->val_p = trial;
      force = interaction->calc_forces();

      if (force < min) {
        for (j = 0; j < D; j++)
          best[j] = trial[j];
        if (output->tempfile != '\0') {
          for (j = 0; j < ndim; j++)
            potential->opt->values[j] = trial[j];
          output->write_tempfile();
        }
        min = force;
      }
      if (force <= cost[i]) {
        for (j = 0; j < D; j++)
          x2[i][j] = trial[j];
        cost[i] = force;
        if (force > max)
          max = force;
      } else {
        for (j = 0; j < D; j++)
          x2[i][j] = x1[i][j];
        if (cost[i] > max)
          max = cost[i];
      }
    }

    if (random->eqdist() < jumprate) {
      opposite_check(x2, cost, 0);
      jsteps++;
      if (jsteps > 10) {
        jumprate *= 0.9;
        jsteps = 0;
      }
    }

    avg = 0.;
    for (i = 0; i < NP; i++)
      avg += cost[i];
    printf("%5d\t\t%15f\t%20f\t\t%.2e\n", count + 1, min, avg / (NP), max - min);
    fflush(stdout);
    for (i = 0; i < NP; i++)
      for (j = 0; j < D; j++)
        x1[i][j] = x2[i][j];
    count++;

    /* End optimization if break flagfile exists */
    if (utils->flagfile != '\0') {
      ff = fopen(utils->flagfile, "r");
      if (NULL != ff) {
        io->write("\nEvolutionary algorithm terminated ");
        io->write("in presence of break flagfile \"%s\"!\n\n", utils->flagfile);
        fclose(ff);
        remove(utils->flagfile);
        break;
      }
    }

    crit = max - min;
  }

  for (j = 0; j < ndim; j++)
    xi[j] = best[j];

  /* clean up */
  for (i = 0; i < NP; i++) {
    delete [] x1[i];
    delete [] x2[i];
  }
  delete [] x1;
  delete [] x2;
  delete [] trial;
  delete [] cost;
  delete [] best;

  return;
}

void OptEvo::init_population(double **pop, double *xi, double *cost) {
  int   i, j;
  int ndim = potential->num_free_params;
  double temp, max, min, val;

  for (i = 0; i < NP; i++) {
    for (j = 0; j < (D - 2); j++)
      pop[i][j] = xi[j];
    pop[i][D - 2] = F_LOWER + random->eqdist() * F_UPPER;
    pop[i][D - 1] = random->eqdist();
  }
  for (i = 1; i < NP; i++) {
    for (j = 0; j < ndim; j++) {
      val = xi[j];
      if (potential->pots[potential->idxpot[j]]->format == 0) {
        min = potential->pots[potential->idxpot[j]]->get_val_min(potential->idxparam[j]);
        max = potential->pots[potential->idxpot[j]]->get_val_max(potential->idxparam[j]);
        /* initialize with normal distribution */
        temp = random->normdist() / 3.;
        if (fabs(temp) > 1)
          temp /= fabs(temp);
        if (temp > 0)
          pop[i][j] = val + temp * (max - val);
        else
          pop[i][j] = val + temp * (val - min);
      } else {
        min = -10. * val;
        max = 10. * val;
        /* initialize with uniform distribution in [-1:1] */
        temp = random->eqdist();
        pop[i][j] = val + temp * (max - min);
      }
    }
  }
  for (i = 0; i < NP; i++) {
    potential->opt->val_p = pop[i];
    cost[i] = interaction->calc_forces();
  }
  opposite_check(pop, cost, 1);

  return;
}

// TODO: check for apot
void OptEvo::opposite_check(double **P, double *costP, int init) {
  int   i, j;
  int ndim = potential->num_free_params;
  double max, min;
  double minp[ndim], maxp[ndim];

  /* allocate memory if not done yet */
  if (tot_P == NULL) {
    tot_P = new double*[2*NP];
    if (tot_P == NULL)
      io->error("Could not allocate memory for opposition vector!\n");
    for (i = 0; i < 2 * NP; i++) {
      tot_P[i] = new double[D];
      for (j = 0; j < D; j++)
        tot_P[i][j] = 0.;
    }
  }
  if (tot_cost == NULL)
    tot_cost = new double[2*NP];
  for (i = 0; i < 2 * NP; i++)
    tot_cost[i] = 0.;

  if (!init) {
    for (i = 0; i < ndim; i++) {
      minp[i] = std::numeric_limits<double>::infinity();
      maxp[i] = -std::numeric_limits<double>::infinity();
    }
    for (i = 0; i < NP; i++) {
      for (j = 0; j < ndim; j++) {
        if (P[i][j] < minp[j])
          minp[j] = P[i][j];
        if (P[i][j] > maxp[j])
          maxp[j] = P[i][j];
      }
    }
  }

  /* generate opposite population */
  for (i = 0; i < NP; i++)
    for (j = 0; j < D; j++)
      tot_P[i][j] = P[i][j];
  for (i = 0; i < NP; i++) {
    for (j = 0; j < ndim; j++) {
      if (init) {
        min = potential->pots[potential->idxpot[j]]->get_val_min(potential->idxparam[j]);
        max = potential->pots[potential->idxpot[j]]->get_val_max(potential->idxparam[j]);
      } else {
        min = minp[j];
        max = maxp[j];
      }
      if (potential->pots[potential->idxpot[j]]->format == 0)
        tot_P[i + NP][j] = min + max - tot_P[i][j];
      else
        tot_P[i + NP][j] = tot_P[i][j];
    }
    tot_P[i + NP][D - 2] = tot_P[i][D - 2];
    tot_P[i + NP][D - 1] = tot_P[i][D - 1];
  }

  /* calculate cost of opposite population */
  for (i = 0; i < NP; i++)
    tot_cost[i] = costP[i];
  for (i = NP; i < 2 * NP; i++) {
    potential->opt->val_p = tot_P[i];
    tot_cost[i] = interaction->calc_forces();
  }

  /* evaluate the NP best individuals from both populations */
  /* sort with quicksort and return NP best indivuals */
  utils->quicksort(tot_cost, 0, 2 * NP - 1, tot_P);
  for (i = 0; i < NP; i++) {
    for (j = 0; j < D; j++)
      P[i][j] = tot_P[i][j];
    costP[i] = tot_cost[i];
  }

  return;
}
