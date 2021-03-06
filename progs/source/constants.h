/*
 * Copyright 2015 Jan von Cosel
 *
 * This file is part of utility-scripts.
 *
 * utility-scripts is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * utility-scripts is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have recieved a copy of the GNU General Public License
 * along with utility-scripts. If not, see <http://www.gnu.org/licenses/>.
 *
 */


#ifndef CONSTANTS_H_
#define CONSTANTS_H_

const double c0 = 299792458.0;
const double Eh = 4.3597438e-18;
const double a0 = 5.291772083e-11;
const double me = 9.10938188e-31;
const double amu = 1.66053873e-27;
const double planck = 6.62606876e-34;
const double hbar = 1.054571596e-34;
const double amu2au = 1822.88848325;
const double Eh2eV = 27.21138344;

const double ang2a0 = 1.0e10 * a0;
const double fs2au = Eh / hbar * 1.0e-15;

// spectrum conversion factor from
// "Computational Strategies for Spectroscopy",
// chapter 2, eq. (2.57)
const double facabs = 703.301;

#endif /* CONSTANTS_H_ */
