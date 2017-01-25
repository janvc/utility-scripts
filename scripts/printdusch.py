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

# read an FCClasses standard output file and print the duschinsky matrix
# as a matrix, so that it can be plotted with gnuplot as an image

import sys
import re
import numpy as np

fileoutputName = sys.argv[1]
Nmodes = int(sys.argv[2])

if Nmodes % 6 == 0:
    Nblocks = Nmodes / 6
else:
    Nblocks = (Nmodes / 6) + 1

fileoutput = open(fileoutputName, "r")
outData = fileoutput.readlines()
fileoutput.close()

# find out where the LAST Duschinsky matrix begins (2nd one, if molecules are being rotated)
duschStartLine = 0
for i in range(len(outData)):
    if "DUSCHINSKY MATRIX" in outData[i]:
        duschStartLine = i

duschEndLine = duschStartLine + Nblocks * (Nmodes + (Nmodes / 10) + (Nmodes / 40) + 7) + 2

duschData = []
for i in range(duschStartLine, duschEndLine):
    tmpList = re.findall(r"[-+]?\d*\.\d+", outData[i])

    if len(tmpList) > 0:
        duschData.append(tmpList)

for i in range(Nblocks):
    del duschData[i * Nmodes]

duschMat = np.zeros((Nmodes, Nmodes))
for i in range(Nblocks):
    for j in range(len(duschData[i * Nmodes])):
        for k in range(Nmodes):
            duschMat[k, i * 6 + j] = float(duschData[i * Nmodes + k][j])

for i in range(Nmodes):
    for j in range(Nmodes):
        sys.stdout.write(" {0:8.5f}".format(duschMat[i,j]))
    sys.stdout.write("\n")

logDusch = np.zeros((Nmodes, Nmodes))
for i in range(Nmodes):
    for j in range(Nmodes):
        if duschMat[i, j] == 0.0:
            logDusch[i, j] = -12
        else:
            logDusch[i, j] = np.log(np.absolute(duschMat[i, j]))

print("\n")

for i in range(Nmodes):
    for j in range(Nmodes):
        sys.stdout.write(" {0:8.2f}".format(logDusch[i,j]))
    sys.stdout.write("\n")

