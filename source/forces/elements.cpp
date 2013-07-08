/****************************************************************
 *
 * elements.cpp:
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

#include "elements.h"

using namespace POTFIT_NS;

Elements::Elements() {
// 1 - Hydrogen
  name.push_back("Hydrogen");
  short_name.push_back("H");
  mass.push_back(1.007900);
// 2 - Helium
  name.push_back("Helium");
  short_name.push_back("He");
  mass.push_back(4.002600);
// 3 - Lithium
  name.push_back("Lithium");
  short_name.push_back("Li");
  mass.push_back(6.941000);
// 4 - Beryllium
  name.push_back("Beryllium");
  short_name.push_back("Be");
  mass.push_back(9.012180);
// 5 - Boron
  name.push_back("Boron");
  short_name.push_back("B");
  mass.push_back(10.810000);
// 6 - Carbon
  name.push_back("Carbon");
  short_name.push_back("C");
  mass.push_back(12.011000);
// 7 - Nitrogen
  name.push_back("Nitrogen");
  short_name.push_back("N");
  mass.push_back(14.006740);
// 8 - Oxygen
  name.push_back("Oxygen");
  short_name.push_back("O");
  mass.push_back(15.999400);
// 9 - Fluorine
  name.push_back("Fluorine");
  short_name.push_back("F");
  mass.push_back(18.998403);
// 10 - Neon
  name.push_back("Neon");
  short_name.push_back("Ne");
  mass.push_back(20.179000);
// 11 - Sodium
  name.push_back("Sodium");
  short_name.push_back("Na");
  mass.push_back(22.989770);
// 12 - Magnesium
  name.push_back("Magnesium");
  short_name.push_back("Mg");
  mass.push_back(24.305000);
// 13 - Aluminum
  name.push_back("Aluminum");
  short_name.push_back("Al");
  mass.push_back(26.981540);
// 14 - Silicon
  name.push_back("Silicon");
  short_name.push_back("Si");
  mass.push_back(28.086000);
// 15 - Phosphorus
  name.push_back("Phosphorus");
  short_name.push_back("P");
  mass.push_back(30.973760);
// 16 - Sulfur
  name.push_back("Sulfur");
  short_name.push_back("S");
  mass.push_back(32.060000);
// 17 - Chlorine
  name.push_back("Chlorine");
  short_name.push_back("Cl");
  mass.push_back(35.453000);
// 18 - Argon
  name.push_back("Argon");
  short_name.push_back("Ar");
  mass.push_back(39.948000);
// 19 - Potassium
  name.push_back("Potassium");
  short_name.push_back("K");
  mass.push_back(39.098000);
// 20 - Calcium
  name.push_back("Calcium");
  short_name.push_back("Ca");
  mass.push_back(40.080000);
// 21 - Scandium
  name.push_back("Scandium");
  short_name.push_back("Sc");
  mass.push_back(44.955900);
// 22 - Titanium
  name.push_back("Titanium");
  short_name.push_back("Ti");
  mass.push_back(47.900000);
// 23 - Vanadium
  name.push_back("Vanadium");
  short_name.push_back("V");
  mass.push_back(50.941400);
// 24 - Cromium
  name.push_back("Cromium");
  short_name.push_back("Cr");
  mass.push_back(51.996000);
// 25 - Manganese
  name.push_back("Manganese");
  short_name.push_back("Mn");
  mass.push_back(54.938000);
// 26 - Iron
  name.push_back("Iron");
  short_name.push_back("Fe");
  mass.push_back(55.847000);
// 27 - Cobalt
  name.push_back("Cobalt");
  short_name.push_back("Co");
  mass.push_back(58.933200);
// 28 - Nickel
  name.push_back("Nickel");
  short_name.push_back("Ni");
  mass.push_back(58.700000);
// 29 - Copper
  name.push_back("Copper");
  short_name.push_back("Cu");
  mass.push_back(63.546000);
// 30 - Zinc
  name.push_back("Zinc");
  short_name.push_back("Zn");
  mass.push_back(65.380000);
// 31 - Gallium
  name.push_back("Gallium");
  short_name.push_back("Ga");
  mass.push_back(69.720000);
// 32 - Germanium
  name.push_back("Germanium");
  short_name.push_back("Ge");
  mass.push_back(72.590000);
// 33 - Arsenic
  name.push_back("Arsenic");
  short_name.push_back("As");
  mass.push_back(74.921600);
// 34 - Selenium
  name.push_back("Selenium");
  short_name.push_back("Se");
  mass.push_back(78.960000);
// 35 - Bromine
  name.push_back("Bromine");
  short_name.push_back("Br");
  mass.push_back(79.904000);
// 36 - Krypton
  name.push_back("Krypton");
  short_name.push_back("Kr");
  mass.push_back(83.800000);
// 37 - Rubidium
  name.push_back("Rubidium");
  short_name.push_back("Rb");
  mass.push_back(85.467800);
// 38 - Strontium
  name.push_back("Strontium");
  short_name.push_back("Sr");
  mass.push_back(87.620000);
// 39 - Yttrium
  name.push_back("Yttrium");
  short_name.push_back("Y");
  mass.push_back(88.905900);
// 40 - Zirconium
  name.push_back("Zirconium");
  short_name.push_back("Zr");
  mass.push_back(91.220000);
// 41 - Niobium
  name.push_back("Niobium");
  short_name.push_back("Nb");
  mass.push_back(92.906400);
// 42 - Molybdenum
  name.push_back("Molybdenum");
  short_name.push_back("Mo");
  mass.push_back(95.940000);
// 43 - Technetium
  name.push_back("Technetium");
  short_name.push_back("Tc");
  mass.push_back(97.000000);
// 44 - Ruthenium
  name.push_back("Ruthenium");
  short_name.push_back("Ru");
  mass.push_back(101.070000);
// 45 - Rhodium
  name.push_back("Rhodium");
  short_name.push_back("Rh");
  mass.push_back(102.905500);
// 46 - Palladium
  name.push_back("Palladium");
  short_name.push_back("Pd");
  mass.push_back(106.400000);
// 47 - Silver
  name.push_back("Silver");
  short_name.push_back("Ag");
  mass.push_back(107.868000);
// 48 - Cadmium
  name.push_back("Cadmium");
  short_name.push_back("Cd");
  mass.push_back(112.400000);
// 49 - Indium
  name.push_back("Indium");
  short_name.push_back("In");
  mass.push_back(114.820000);
// 50 - Tin
  name.push_back("Tin");
  short_name.push_back("Sn");
  mass.push_back(118.690000);
// 51 - Antimony
  name.push_back("Antimony");
  short_name.push_back("Sb");
  mass.push_back(121.750000);
// 52 - Tellurium
  name.push_back("Tellurium");
  short_name.push_back("Te");
  mass.push_back(127.600000);
// 53 - Iodine
  name.push_back("Iodine");
  short_name.push_back("I");
  mass.push_back(126.904500);
// 54 - Xenon
  name.push_back("Xenon");
  short_name.push_back("Xe");
  mass.push_back(131.300000);
// 55 - Cesium
  name.push_back("Cesium");
  short_name.push_back("Cs");
  mass.push_back(132.905400);
// 56 - Barium
  name.push_back("Barium");
  short_name.push_back("Ba");
  mass.push_back(137.340000);
// 57 - Lanthanum
  name.push_back("Lanthanum");
  short_name.push_back("La");
  mass.push_back(138.905500);
// 58 - Cerium
  name.push_back("Cerium");
  short_name.push_back("Ce");
  mass.push_back(140.120000);
// 59 - Praseodynium
  name.push_back("Praseodynium");
  short_name.push_back("Pr");
  mass.push_back(140.907700);
// 60 - Neodymium
  name.push_back("Neodymium");
  short_name.push_back("Nd");
  mass.push_back(144.240000);
// 61 - Promethium
  name.push_back("Promethium");
  short_name.push_back("Pm");
  mass.push_back(145.000000);
// 62 - Samarium
  name.push_back("Samarium");
  short_name.push_back("Sm");
  mass.push_back(150.400000);
// 63 - Europium
  name.push_back("Europium");
  short_name.push_back("Eu");
  mass.push_back(151.960000);
// 64 - Gadolinium
  name.push_back("Gadolinium");
  short_name.push_back("Gd");
  mass.push_back(157.250000);
// 65 - Terbium
  name.push_back("Terbium");
  short_name.push_back("Tb");
  mass.push_back(158.925400);
// 66 - Dysprosium
  name.push_back("Dysprosium");
  short_name.push_back("Dy");
  mass.push_back(162.500000);
// 67 - Holmium
  name.push_back("Holmium");
  short_name.push_back("Ho");
  mass.push_back(164.930400);
// 68 - Erbium
  name.push_back("Erbium");
  short_name.push_back("Er");
  mass.push_back(167.260000);
// 69 - Thulium
  name.push_back("Thulium");
  short_name.push_back("Tm");
  mass.push_back(168.934200);
// 70 - Ytterbium
  name.push_back("Ytterbium");
  short_name.push_back("Yb");
  mass.push_back(173.040000);
// 71 - Lutetium
  name.push_back("Lutetium");
  short_name.push_back("Lu");
  mass.push_back(174.970000);
// 72 - Hafnium
  name.push_back("Hafnium");
  short_name.push_back("Hf");
  mass.push_back(178.490000);
// 73 - Tantalum
  name.push_back("Tantalum");
  short_name.push_back("Ta");
  mass.push_back(180.947900);
// 74 - Tungsten
  name.push_back("Tungsten");
  short_name.push_back("W");
  mass.push_back(183.500000);
// 75 - Rhenium
  name.push_back("Rhenium");
  short_name.push_back("Re");
  mass.push_back(186.207000);
// 76 - Osmium
  name.push_back("Osmium");
  short_name.push_back("Os");
  mass.push_back(190.200000);
// 77 - Iridium
  name.push_back("Iridium");
  short_name.push_back("Ir");
  mass.push_back(192.220000);
// 78 - Platinum
  name.push_back("Platinum");
  short_name.push_back("Pt");
  mass.push_back(195.090000);
// 79 - Gold
  name.push_back("Gold");
  short_name.push_back("Au");
  mass.push_back(196.966500);
// 80 - Mercury
  name.push_back("Mercury");
  short_name.push_back("Hg");
  mass.push_back(200.590000);
// 81 - Thallium
  name.push_back("Thallium");
  short_name.push_back("Tl");
  mass.push_back(204.370000);
// 82 - Lead
  name.push_back("Lead");
  short_name.push_back("Pb");
  mass.push_back(207.200000);
// 83 - Bismuth
  name.push_back("Bismuth");
  short_name.push_back("Bi");
  mass.push_back(208.980400);
// 84 - Polonium
  name.push_back("Polonium");
  short_name.push_back("Po");
  mass.push_back(209.000000);
// 85 - Astatine
  name.push_back("Astatine");
  short_name.push_back("At");
  mass.push_back(210.000000);
// 86 - Radon
  name.push_back("Radon");
  short_name.push_back("Rn");
  mass.push_back(222.000000);
// 87 - Francium
  name.push_back("Francium");
  short_name.push_back("Fr");
  mass.push_back(223.000000);
// 88 - Radium
  name.push_back("Radium");
  short_name.push_back("Ra");
  mass.push_back(226.025400);
// 89 - Actinium
  name.push_back("Actinium");
  short_name.push_back("Ac");
  mass.push_back(227.000000);
// 90 - Thorium
  name.push_back("Thorium");
  short_name.push_back("Th");
  mass.push_back(232.038100);
// 91 - Protactinium
  name.push_back("Protactinium");
  short_name.push_back("Pa");
  mass.push_back(231.035900);
// 92 - Uranium
  name.push_back("Uranium");
  short_name.push_back("U");
  mass.push_back(238.029000);
// 93 - Neptunium
  name.push_back("Neptunium");
  short_name.push_back("Np");
  mass.push_back(237.048200);
// 94 - Plutonium
  name.push_back("Plutonium");
  short_name.push_back("Pu");
  mass.push_back(244.000000);
// 95 - Americium
  name.push_back("Americium");
  short_name.push_back("Am");
  mass.push_back(243.000000);
// 96 - Curium
  name.push_back("Curium");
  short_name.push_back("Cm");
  mass.push_back(247.000000);
// 97 - Berkelium
  name.push_back("Berkelium");
  short_name.push_back("Bk");
  mass.push_back(247.000000);
// 98 - Californium
  name.push_back("Californium");
  short_name.push_back("Cf");
  mass.push_back(251.000000);
// 99 - Einsteinium
  name.push_back("Einsteinium");
  short_name.push_back("Es");
  mass.push_back(252.000000);
// 100 - Fermium
  name.push_back("Fermium");
  short_name.push_back("Fm");
  mass.push_back(257.000000);
// 101 - Mendelevium
  name.push_back("Mendelevium");
  short_name.push_back("Md");
  mass.push_back(258.000000);
// 102 - Nobelium
  name.push_back("Nobelium");
  short_name.push_back("No");
  mass.push_back(259.000000);
// 103 - Lawrencium
  name.push_back("Lawrencium");
  short_name.push_back("Lr");
  mass.push_back(262.000000);
// 104 - Rutherfordium
  name.push_back("Rutherfordium");
  short_name.push_back("Rf");
  mass.push_back(261.000000);
// 105 - Dubnium
  name.push_back("Dubnium");
  short_name.push_back("Db");
  mass.push_back(262.000000);
// 106 - Seaborgium
  name.push_back("Seaborgium");
  short_name.push_back("XX");
  mass.push_back(263.000000);
// 107 - Bohrium
  name.push_back("Bohrium");
  short_name.push_back("Bf");
  mass.push_back(262.000000);
// 108 - Hassium
  name.push_back("Hassium");
  short_name.push_back("Hs");
  mass.push_back(265.000000);
// 109 - Meitnerium
  name.push_back("Meitnerium");
  short_name.push_back("Mt");
  mass.push_back(265.000000);

  return;
}

Elements::~Elements() {
  name.clear();
  short_name.clear();
  mass.clear();

  return;
}

double Elements::mass_from_number(const int &num) {
  if (num > name.size() || num < 1) {
    return 0.0;
  } else {
    return mass[num];
  }
}

double Elements::mass_from_name(const std::string &n) {
  int   i = 0;

  if (n.length() < 3) {
    for (i = 1; i < name.size(); i++)
      if (n.compare(short_name[i]) == 0)
	return mass[i];
  } else {
    for (i = 1; i < name.size(); i++)
      if (n.compare(name[i]) == 0)
	return mass[i];
  }

  return 0.0;
}

int Elements::number_from_name(const std::string &n) {
  int   i = 0;

  if (n.length() < 3) {
    for (i = 1; i < name.size(); i++)
      if (n.compare(short_name[i]) == 0)
	return i;
  } else {
    for (i = 1; i < name.size(); i++)
      if (n.compare(name[i]) == 0)
	return i;
  }

  return 0.0;
}
