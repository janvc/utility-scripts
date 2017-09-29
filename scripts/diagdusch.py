#!/usr/bin/env python3

# Copyright 2017 Jan von Cosel
#
# This file is part of utility-scripts.
#
# utility-scripts is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# utility-scripts is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have recieved a copy of the GNU General Public License
# along with utility-scripts. If not, see <http://www.gnu.org/licenses/>.
#
#

# this script takes an MCTDH operator file used for a spectrum calculation
# and diagonalizes the harmonic potential energy surface to get the
# normal modes' force constants

import sys
import numpy as np

opFileName = sys.argv[1]

def f2float(numString):
    return float(numString.replace('d', 'e'))

def main():
    opFile = open(opFileName, "r")
    opData = opFile.readlines()
    opFile.close()

    # determine start and end of the parameter-section
    for i in range(len(opData)):
        if opData[i].startswith("parameter-section"):
            paraStart = i+1
        if opData[i].startswith("end-parameter-section"):
            paraEnd = i
    Npara = paraEnd - paraStart

    # determine the number of modes
    Nmodes = 0
    for i in range(Npara):
        if opData[paraStart + i].startswith("    mass_q_"):
            Nmodes += 1

    Fmat = np.zeros((Nmodes, Nmodes))

    # read the projected force constants
    for i in range(Nmodes):
        Fmat[i, i] = f2float(opData[paraStart + 3 * Nmodes + i].split()[2])

    # read the couplings
    refIndex = 4 * Nmodes
    for i in range(Nmodes):
        for j in range(i+1, Nmodes):
            vecIndex = int(i * (Nmodes - ((i + 1) / 2)) + j - i - 1)
            Fmat[i, j] = f2float(opData[paraStart + 4 * Nmodes + vecIndex].split()[2])
            Fmat[j, i] = Fmat[i, j]


    for i in range(Nmodes):
        for j in range(Nmodes):
            sys.stdout.write(" {0:14.7e}".format(Fmat[i, j]))
        sys.stdout.write("\n")

    diagFC, trafo = np.linalg.eig(Fmat)
    idx = diagFC.argsort()[::1]
    diagFC = diagFC[idx]
    print(diagFC)

    # compute the ground state zero-point energy
    zpe0 = 0.0
    for i in range(Nmodes):
        zpe0 += 0.5 * np.sqrt(f2float(opData[paraStart + Nmodes + i].split()[2]))

    # compute the excited state zero-point energy
    zpe1 = 0.0
    for i in range(Nmodes):
        zpe1 += 0.5 * np.sqrt(f2float(opData[paraStart + 2 * Nmodes + i].split()[2]))

    # compute the effective excited state zero-point energy
    zpeE = 0.0
    for i in range(Nmodes):
        zpeE += 0.5 * np.sqrt(diagFC[i])

    print(zpe0, zpe1, zpeE)



if __name__ == "__main__":
    main()

