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

import os
import sys
import subprocess

mctdhExe = "/home/jvcosel/mctdh85.4/bin/binary/x86_64/mctdh85"

initInputName = sys.argv[1]
natPopThres = float(sys.argv[2])

# create the first real input file:
initBaseName = os.path.splitext(initInputName)[0]

initInputFile = open(initInputName, "r")
initInputData = initInputFile.readlines()
initInputFile.close()

currentNameDir = initBaseName + "_01"
currentInputName = currentNameDir + ".inp"

currentInputData = initInputData
for i in range(len(initInputData)):
    if ("name = " in initInputData[i]):
        currentInputData[i] = "    name = " + currentNameDir + "\n"
        break

currentInputFile = open(currentInputName, "w")
for item in currentInputData:
    currentInputFile.write(item)
currentInputFile.close()

# do the first MCTDH calculation:
mctdhCommand = [mctdhExe, "-mnd", currentInputName]
subprocess.check_call(mctdhCommand)

# analyze the natural populations:
popFileName = currentNameDir + "/lownatpop"
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
sortedList = sorted(outputList, key=lambda l:l[3], reverse=True)
maxPop = float(sortedList[0][3])

# determine the starting point of the ML basis section:
for i in range(len(currentInputData)):
    if "mlbasis-section" in currentInputData[i]:
        mlStart = i
        break

# now comes the big calculation loop:
while maxPop > natPopThres:
    
    # create the new input file and improve the SPFs:
    newInputData = currentInputData
    for i in range(len(sortedList)):
        if sortedList[i][3] > natPopThres:
            currLNr = mlStart + int(sortedList[i][1][5:])
            currLn = currentInputData[currLNr]
            currMd = int(sortedList[i][2][5:])
            currSPFId = currLn.index('>') + (2 * currMd)
            currSPFs = int(currLn[currSPFId])
            newInputData[currLNr] = currentInputData[currLNr][:currSPFId - 1] + " " + str(currSPFs + 1) + currentInputData[currLNr][currSPFId + 1:]
    
    # create the new MCTDH input file:
    currentNumber = int(currentNameDir[-2:])
    newNumber = currentNumber + 1
    newNameDir = initBaseName + "_" + str(newNumber).zfill(2)
    newInputName = newNameDir + ".inp"
    newInputFile = open(newInputName, "w")
    for item in newInputData:
        if "name = " in item and not "opname" in item:
            item = "    name = " + newNameDir + "\n"
        newInputFile.write(item)
    newInputFile.close()
    currentInputName = newInputName
    currentNameDir = newNameDir
    
    # run the new calculation:
    mctdhCommand = [mctdhExe, "-mnd", currentInputName]
    subprocess.check_call(mctdhCommand)
    
    # analyze the natural populations:
    popFileName = currentNameDir + "/lownatpop"
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
    sortedList = sorted(outputList, key=lambda l:l[3], reverse=True)
    maxPop = float(sortedList[0][3])


