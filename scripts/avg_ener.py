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

avgData = [[0.0 for i in range(nModes)] for j in range(nPoints)]

for i in range(nRun):
    currFile = open("run_" + str(i+1) + "/" + nameDir + "/ener.log", "r")
    currData = currFile.readlines()
    currFile.close()

    for j in range(nPoints):
        for k in range(nModes):
            avgData[j][k] += float(currData[j+1].split()[k+1])

for i in range(nPoints):
    for j in range(nModes):
        avgData[i][j] = avgData[i][j] / nRun

with open("avg_ener.dat", "w") as avgFile:
    for i in range(nPoints):
        for j in range(nModes):
            avgFile.write("{0:15.5e}".format(avgData[i][j]))
        avgFile.write("\n")

