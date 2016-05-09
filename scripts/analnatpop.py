#!/usr/bin/env python

# Copyright 2015 Jan von Cosel
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

popFileName = "lownatpop"
doSort = ""
if len(sys.argv) == 2:
    doSort = sys.argv[1]

popList = []

with open(popFileName) as popFile:
    for line in popFile:
        if "#" in line:     # the layer/node/mode information
            # rule out the first line:
            if not "lownatpop" in line:
                if "layer" in line:
                    layerList = line.split()
                    popList = [0]*(len(layerList)-2)
                elif "node" in line:
                    nodeList = line.split()
                elif "mode" in line:
                    modeList = line.split()
        else:
            for i in range(len(popList)):
                if float(line.split()[i+1]) > float(popList[i]):
                    popList[i] = line.split()[i+1]

outputList = []
for i in range(len(popList)):
    tmpList = []
    tmpList.append(layerList[i+2])
    tmpList.append(nodeList[i+1])
    tmpList.append(modeList[i+1])
    tmpList.append(float(popList[i]))
    outputList.append(tmpList)

# do not sort if '-s' is given:
if "-s" in doSort:
    sortedList = outputList
else:
    sortedList = sorted(outputList, key=lambda l:l[3], reverse=True)

for i in range(len(sortedList)):
    print("{0:10s} {1:10s} {2:10s} {3:5.5f}".format(sortedList[i][0], sortedList[i][1], sortedList[i][2], sortedList[i][3]))


