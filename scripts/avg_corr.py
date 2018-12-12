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


import sys


nRun = int(sys.argv[1])
nameDir = "ml_prop_001"
nPoints = 1001
nModes = 9
nStates = 4


for i in range(nModes):
    for j in range(nStates):
        avgData = [[5.0 * l, 0.0] for l in range(nPoints)]

        for k in range(nRun):
            currFile = open("run_" + str(k+1) + "/thermcorr_00" + str(i+1) + "_" + str(j+2) + ".dat", "r")
            currData = currFile.readlines()
            currFile.close()

            for l in range(nPoints):
                avgData[l][1] += float(currData[l+8].split()[3])

        for l in range(nPoints):
            avgData[l][1] = avgData[l][1] / nRun

        with open("avgcorr_" + str(i+1) + "_" + str(j+2) + ".dat", "w") as avgFile:
            for l in range(nPoints):
                avgFile.write("{0:15.5e} {1:15.5e}\n".format(avgData[l][0], avgData[l][1]))

