/****************************************************************
 *
 * config.cpp:
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
#include <cstring>
#include <cstdio>
#include <iostream>

#include "config.h"

#include "../io.h"
#include "../potential.h"
#include "../structures.h"
#include "../templates.h"

using namespace POTFIT_NS;

Config::Config(POTFIT *ptf) : Pointers(ptf) {
  use_forces = 0;
  use_stresses = 0;
  num_atoms = 0;

  coh_energy = 0.0;
  conf_weight = 0.0;
  volume = 0.0;

  stress.xx = 0.0;
  stress.yy = 0.0;
  stress.zz = 0.0;
  stress.xy = 0.0;
  stress.zx = 0.0;
  stress.yz = 0.0;

  box_x.x = box_x.y = box_x.z = 0.0;
  box_y.x = box_y.y = box_y.z = 0.0;
  box_z.x = box_z.y = box_z.z = 0.0;

  return;
}

Config::~Config() {
  for (unsigned i = 0; i < atoms.size(); ++i)
    delete atoms[i];
  atoms.clear();

  num_per_type.clear();
  return;
}

void Config::read(FILE *infile, int *line) {
  char *res, *ptr;
  char buffer[1024];
  int h_eng = 0, h_stress = 0, h_boxx = 0, h_boxy = 0, h_boxz = 0;
  int have_contrib = 0;
  int have_contrib_box = 0;
  int cell_scale[3];
  vector iheight;
  vector tbox_x, tbox_y, tbox_z;
  vector cbox_o;				// origin of box of contrib. atoms
  vector cbox_a, cbox_b, cbox_c;		// box vectors for box of contrib. atoms
  vector temp_vect; 				//
  double temp; 					//
  std::vector<vector> sphere_centers;		// centers of the spheres of contrib. atoms
  std::vector<double> sphere_radii;
  int num_spheres = 0;
  fpos_t filepos;

  res = fgets(buffer, 1024, infile);
  (*line)++;
  if (NULL == res)
    io->error("Unexpected end of config file on line %d.", *line);
  h_eng = h_stress = h_boxx = h_boxy = h_boxz = 0;
  if (res[1] == 'N') {	// Atom number line
    if (sscanf(res + 3, "%d %d", &num_atoms, &use_forces) < 2)
      io->error("Error in atom number specification on line %d\n", *line);
  } else {
    io->error("Error - number of atoms missing on line %d\n", *line);
  }

  for (int i=0; i<structures->ntypes; i++)
    num_per_type.push_back(0);

  do {
    res = fgets(buffer, 1024, infile);
    if ((ptr = strchr(res, '\n')) != NULL)
      *ptr = '\0';
    line++;
    // read the box vectors
    if (res[1] == 'X') {
      if (sscanf(res + 3, "%lf %lf %lf\n", &box_x.x, &box_x.y, &box_x.z) == 3)
        h_boxx++;
      else
        io->error("Error in box vector x, line %d\n", *line);
    } else if (res[1] == 'Y') {
      if (sscanf(res + 3, "%lf %lf %lf\n", &box_y.x, &box_y.y, &box_y.z) == 3)
        h_boxy++;
      else
        io->error("Error in box vector y, line %d\n", line);
    } else if (res[1] == 'Z') {
      if (sscanf(res + 3, "%lf %lf %lf\n", &box_z.x, &box_z.y, &box_z.z) == 3)
        h_boxz++;
      else
        io->error("Error in box vector z, line %d\n", line);
      // read the box of contributing particles
    } else if (strncmp(res + 1, "B_O", 3) == 0) {
      if (1 == have_contrib_box) {
        io->error("There can only be one box of contributing atoms!\n");
      }
      if (sscanf(res + 5, "%lf %lf %lf\n", &cbox_o.x, &cbox_o.y, &cbox_o.z) == 3) {
        have_contrib_box = 1;
        have_contrib++;
      } else
        io->error("Error in box of contributing atoms, line %d\n", line);
    } else if (strncmp(res + 1, "B_A", 3) == 0) {
      if (sscanf(res + 5, "%lf %lf %lf\n", &cbox_a.x, &cbox_a.y, &cbox_a.z) == 3) {
        have_contrib++;
      } else
        io->error("Error in box of contributing atoms, line %d\n", line);
    } else if (strncmp(res + 1, "B_B", 3) == 0) {
      if (sscanf(res + 5, "%lf %lf %lf\n", &cbox_b.x, &cbox_b.y, &cbox_b.z) == 3) {
        have_contrib++;
      } else
        io->error("Error in box of contributing atoms, line %d\n", line);
    } else if (strncmp(res + 1, "B_C", 3) == 0) {
      if (sscanf(res + 5, "%lf %lf %lf\n", &cbox_c.x, &cbox_c.y, &cbox_c.z) == 3) {
        have_contrib++;
      } else
        io->error("Error in box of contributing atoms, line %d\n", line);
    } else if (strncmp(res + 1, "B_S", 3) == 0) {
      if (sscanf(res + 5, "%lf %lf %lf %lf\n", &temp_vect.x, &temp_vect.y, &temp_vect.z, &temp) == 4) {
        num_spheres++;
        sphere_centers.push_back(temp_vect);
        sphere_radii.push_back(temp);
      } else {
        io->error("Error in sphere of contributing atoms, line %d\n", line);
      }
    } else if (res[1] == 'E') {
      if (sscanf(res + 3, "%lf\n", &coh_energy) == 1)
        h_eng++;
      else
        io->error("Error in energy on line %d\n", line);
    } else if (res[1] == 'W') {
      if (sscanf(res + 3, "%lf\n", &conf_weight) != 1)
        io->error("Error in configuration weight on line %d\n", line);
    } else if (res[1] == 'C') {
      // TODO
      // read C line and compare to values from potential file
    }

    // read stress
    else if (res[1] == 'S') {
      if (sscanf(res + 3, "%lf %lf %lf %lf %lf %lf\n", &stress.xx, &stress.yy, &stress.zz, &stress.xy,
                 &stress.yz, &stress.zx) == 6)
        use_stresses = 1;
      else
        io->error("Error in stress tensor on line %d\n", *line);
    }
  } while (res[1] != 'F');

  if (!(h_eng && h_boxx && h_boxy && h_boxz))
    io->error("Incomplete box vectors on line %d!", *line);
  if (have_contrib_box && have_contrib != 4)
    io->error("Incomplete box of contributing atoms on line %d!", *line);

  if (use_stresses)
    structures->using_stresses++;
  if (use_forces)
    structures->using_forces++;

  // compute tbox_j such that SPROD(box_i,tbox_j) == delta_ij
  // first unnormalized
  tbox_x = vec_prod(box_y, box_z);
  tbox_y = vec_prod(box_z, box_x);
  tbox_z = vec_prod(box_x, box_y);

  // volume
  volume = SPROD(box_x, tbox_x);
  if (0 == volume)
    io->error("Box edges are parallel\n");

  // normalization
  tbox_x.x /= volume;
  tbox_x.y /= volume;
  tbox_x.z /= volume;
  tbox_y.x /= volume;
  tbox_y.y /= volume;
  tbox_y.z /= volume;
  tbox_z.x /= volume;
  tbox_z.y /= volume;
  tbox_z.z /= volume;

  // read the atoms
  Atom *atom;
  for (int i = 0; i < num_atoms; i++) {
    atoms.push_back(new Atom(ptf));
    atom = atoms[i];
    if (7 > fscanf(infile, "%d %lf %lf %lf %lf %lf %lf\n", &(atom->type),
                   &(atom->pos.x), &(atom->pos.y), &(atom->pos.z), &(atom->force.x), &(atom->force.y),
                   &(atom->force.z)))
      io->error("Corrupt configuration file on line %d\n", *line + 1);
    line++;
    if (atom->type >= structures->ntypes || atom->type < 0)
      io->error("Corrupt configuration file on line %d: Incorrect atom type (%d)\n", *line, atom->type);
    atom->absforce = sqrt(square(atom->force.x) + square(atom->force.y) + square(atom->force.z));
    if (have_contrib_box || num_spheres != 0) {
      double n_a = 0.0, n_b = 0, n_c = 0.0, r = 0.0;

      atom->contrib = 0;
      if (have_contrib_box) {
        temp_vect.x = atom->pos.x - cbox_o.x;
        temp_vect.y = atom->pos.y - cbox_o.y;
        temp_vect.z = atom->pos.z - cbox_o.z;
        n_a = SPROD(temp_vect, cbox_a) / SPROD(cbox_a,cbox_a);
        n_b = SPROD(temp_vect, cbox_b) / SPROD(cbox_b,cbox_b);
        n_c = SPROD(temp_vect, cbox_c) / SPROD(cbox_c,cbox_c);
        if (n_a >= 0 && n_a <= 1)
          if (n_b >= 0 && n_b <= 1)
            if (n_c >= 0 && n_c <= 1)
              atom->contrib = 1;
      }

      for (int i = 0; i < num_spheres; i++) {
        temp_vect.x = (atom->pos.x - sphere_centers[i].x);
        temp_vect.y = (atom->pos.y - sphere_centers[i].y);
        temp_vect.z = (atom->pos.z - sphere_centers[i].z);
        r = SPROD(temp_vect, temp_vect);
        r = sqrt(r);
        if (r < sphere_radii[i])
          atom->contrib = 1;
      }
    }
    else
      atom->contrib = 1;
    num_per_type[atom->type] += 1;
  }
  // check cell size
  // inverse height in direction
  iheight.x = sqrt(SPROD(tbox_x, tbox_x));
  iheight.y = sqrt(SPROD(tbox_y, tbox_y));
  iheight.z = sqrt(SPROD(tbox_z, tbox_z));

  if ((ceil(potential->rcut_max * iheight.x) > 30000)
      || (ceil(potential->rcut_max * iheight.y) > 30000)
      || (ceil(potential->rcut_max * iheight.z) > 30000))
    io->error("Very bizarre small cell size - aborting");

  cell_scale[0] = (int)ceil(potential->rcut_max * iheight.x);
  cell_scale[1] = (int)ceil(potential->rcut_max * iheight.y);
  cell_scale[2] = (int)ceil(potential->rcut_max * iheight.z);

#ifdef DEBUG
  io->write_debug("Checking cell size for configuration %d:\n", nconf);
  io->write_debug("Box dimensions:\n");
  io->write_debug("     %10.6f %10.6f %10.6f\n", box_x.x, box_x.y, box_x.z);
  io->write_debug("     %10.6f %10.6f %10.6f\n", box_y.x, box_y.y, box_y.z);
  io->write_debug("     %10.6f %10.6f %10.6f\n", box_z.x, box_z.y, box_z.z);
  io->write_debug("Box normals:\n");
  io->write_debug("     %10.6f %10.6f %10.6f\n", tbox_x.x, tbox_x.y, tbox_x.z);
  io->write_debug("     %10.6f %10.6f %10.6f\n", tbox_y.x, tbox_y.y, tbox_y.z);
  io->write_debug("     %10.6f %10.6f %10.6f\n", tbox_z.x, tbox_z.y, tbox_z.z);
  io->write_debug("Box heights:\n");
  io->write_debug("     %10.6f %10.6f %10.6f\n", 1. / iheight.x, 1. / iheight.y, 1. / iheight.z);
  io->write_debug("Potential range:  %f\n", rcutmax);
  io->write_debug("Periodic images needed: %d %d %d\n\n",
          2 * cell_scale[0] + 1, 2 * cell_scale[1] + 1, 2 * cell_scale[2] + 1);
#endif /* DEBUG */

  calc_neighbors();

  return;
}

void Config::calc_neighbors(void) {

  return;
}

vector Config::vec_prod(vector a, vector b) {
  vector ret;
  ret.x = a.y * b.z - a.z * b.y;
  ret.y = a.z * b.x - a.x * b.z;
  ret.z = a.x * b.y - a.y * b.x;

  return ret;
}

