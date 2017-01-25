#!/usr/bin/env python

# Copyright 2016 Jan von Cosel
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

# create MCTDH inputs for excited vibrational states up to a specific
# energy threshold

import sys
import functools
import operator
import os
import math

def main():
    initInputName = sys.argv[1]
    stateNum = int(sys.argv[2])

    initBaseName = os.path.splitext(initInputName)[0]

    initInputFile = open(initInputName, "r")
    initInputData = initInputFile.readlines()
    initInputFile.close()

    # determine the name of the corresponding operator file:
    for i in range(len(initInputData)):
        if ("opname = " in initInputData[i]):
            operFileName = initInputData[i].split()[2] + ".op"

    operFile = open(operFileName, "r")
    operData = operFile.readlines()
    operFile.close()

    # determine the number of DOFs (and the starting
    # line of the force constant block):
    Nmodes = 0
    for i in range(len(operData)):
        if ("mass_q_" in operData[i]):
            Nmodes += 1
            FCline = i

    print("The system has %d modes\n" % Nmodes)

    # read the mode energies:
    Energies = []
    for i in range(Nmodes):
        fcString = str(operData[FCline + 1 + i].split()[2])
        Energies.append(math.sqrt(float(fcString.replace("d", "e"))))

    # estimate the number of modes that need to be excited.
    # if it turns out that this number doesn't suffice, just
    # use the requested mode index as maximum
    Nexc = Nmodes / 10

    while True:
        Emax = Energies[Nexc - 1]
        maxQuanta = [None] * Nexc
        for i in range(Nexc):
            maxQuanta[i] = int(Emax / Energies[i]) + 1


        maxStates = functools.reduce(operator.mul, maxQuanta, 1)

        if (maxStates < stateNum):
            Nexc += 1
        else:
            break

    print("Considering excitations of the first %d modes\n" % Nexc)
    print("up to an energy of %f Eh\n" % Emax)
    print("Highest excited state: ", maxQuanta)
    print("total number of states: %d\n" % maxStates)


    stateVectors = [[0 for x in range(Nexc)] for y in range(maxStates)]

    # create the excited vibrational states:
    for kp in range(Nexc):
        nl = 1
        nr = 1
        for i in range(kp):
            nl *= maxQuanta[i]
        for i in range(Nexc - 1, kp, -1):
            nr *= maxQuanta[i]
        for i in range(nr):
            for j in range(maxQuanta[kp]):
                for k in range(nl):
                    w = i * nl * maxQuanta[kp] + j * nl + k
                    stateVectors[w][kp] = j

    # determine their energies:
    stateEnergies = []
    for i in range(len(stateVectors)):
        stateEnergies.append(0.0)
        for j in range(Nexc):
            stateEnergies[i] += stateVectors[i][j] * Energies[j]

    combinedVector = zip(stateEnergies, stateVectors)
    combinedVector.sort()

    for i in range(len(combinedVector)):
        print(combinedVector[i][0])

    # create the new input file:
    currentBaseName = initBaseName + "_" + str(stateNum).zfill(3)

    currentInputData = []
    currentInputData.append("#\n# excited state "
                            + str(combinedVector[stateNum - 1][1:])
                            + "\n# energy "
                            + str(combinedVector[stateNum - 1][0]) + " Eh\n#\n")
    for j in range(len(initInputData)):
        currentInputData.append(initInputData[j])

    for j in range(len(currentInputData)):
        if (" name = " in currentInputData[j]):
            currentInputData[j] = "    name = " + currentBaseName + "\n"
        if (" eigenf " in currentInputData[j]):
            mode = int(str(currentInputData[j][10:13])) - 1
            if (mode < Nexc):
                currentInputData[j] = "        q_" + str(mode + 1).zfill(3) + "  eigenf  Eq_" + str(mode + 1).zfill(3) + "  pop = " + str(int(combinedVector[stateNum - 1][1][mode]) + 1) + "\n"
    currentInputName = currentBaseName + ".inp"
    currentInputFile = open(currentInputName, "w")
    for item in currentInputData:
        currentInputFile.write(item)
    currentInputFile.close()


if __name__ == '__main__':
    main()


