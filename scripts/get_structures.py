#!/usr/bin/env python3

# Copyright 2018 Jan von Cosel
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

atomList = ["H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar"]

if len(sys.argv) == 1 or "-h" in sys.argv:
    print("get_structures: extract optimized geometries from a Gaussian relaxed surface scan")
    print("and form new inputs, e.g. for single points")
    print("")
    print("usage: get_structures <scan log file> <head of inputs> <foot of inputs> <basename of inputs>")
    print("")
    sys.exit()

scanFileName = sys.argv[1]
headerFile = sys.argv[2]
footerFile = sys.argv[3]
basename = sys.argv[4]

scanFile = open(scanFileName, "r")
scanData = scanFile.readlines()
scanFile.close()

for line in scanData:
    if "NAtoms=" in line:
        natoms = int(line.split()[1])
        break

# determine the starting line of each optimized input orientation
startLines = []
i = 0
while True:
    line = scanData[i]
    i += 1

    if "Optimized Parameters" in line:
        j = i
        while True:
            line2 = scanData[j]
            j -= 1
            if "Input orientation" in line2:
                startLines.append(j+6)
                break

    if i > len(scanData) - 10:
        break

for i in range(len(startLines)):
    fileName = basename + "_" + str(i+1).zfill(4) + ".com"
    with open(fileName, "w") as f1:
        with open(headerFile, "r") as f2:
            for line in f2.readlines():
                f1.write(line)

        sys.stdout.write(str(natoms) + "\n\n")
        for j in range(natoms):
            struLine = scanData[startLines[i]+j].split()
            writeString = str(struLine[1]) + " " + str(struLine[3]) + " " + str(struLine[4]) + " " + str(struLine[5]) + "\n"
            f1.write(writeString)
            atNum = int(struLine[1])
            if atNum <= 18:
                writeString = atomList[atNum - 1]
            else:
                writeString = "X"
            writeString += " " + str(struLine[3]) + " " + str(struLine[4]) + " " + str(struLine[5]) + "\n"
            sys.stdout.write(writeString)

        with open(footerFile, "r") as f2:
            for line in f2.readlines():
                f1.write(line)

        f1.write("\n\n")



