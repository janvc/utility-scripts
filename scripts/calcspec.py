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
import time


#
# some global variables
#
mctdhExe = "/home/jvcosel/mctdh85.4/bin/binary/x86_64/mctdh85"
initInputName = sys.argv[1]
natPopThres = float(sys.argv[2])


#
# determine the maximum natural population in an ML-MCTDH-calculation
#
def analNatPop(nameDir):
    popFileName = nameDir + "/lownatpop"
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

    return sortedList

def getCurrentTime(nameDir):
    speedFileName = nameDir + "/speed"

    currentTime = 0.0
    with open(speedFileName) as speedFile:
        speedList = speedFile.readlines()
        if not "#" in speedList[-1]:
            currentTime = float(speedList[-1].split()[0])

    return currentTime


#
# here begins the main function
#
if __name__ == '__main__':

    # create the first real input file:
    initBaseName = os.path.splitext(initInputName)[0]

    initInputFile = open(initInputName, "r")
    initInputData = initInputFile.readlines()
    initInputFile.close()

    currentNameDir = initBaseName + "_001"
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

    # determine the start of the ML basis section:
    for i in range(len(currentInputData)):
        if ("mlbasis-section" in currentInputData[i]):
            mlStart = i
            break

    # determine the length of the propagation:
    tfinal = 0.0
    for i in range(len(initInputData)):
        if ("tfinal" in initInputData[i]):
            tfinal = float(initInputData[i].split()[2])



    #
    # this is the big calculation loop
    #
    converged = False
    while True:

        # run the calculation:
        print "Starting calculation " + str(currentNameDir)
        mctdhCommand = [mctdhExe, "-mnd", currentInputName]

        MCTDHproc = subprocess.Popen(mctdhCommand)

        # wait for five minutes to avoid timing issues:
        time.sleep(300)

        # while the calculation is running,
        # check periodically, if we have violated the convergence criterion:
        while True:

            # analyze the natural populations:
            currentPopList = analNatPop(currentNameDir)
            currentTime = getCurrentTime(currentNameDir)

            print "At " + str(currentTime) + " fs, max population is " + str(currentPopList[0][3])

            violated = False
            if (currentPopList[0][3] > natPopThres):
                violated = True

            # check if the convergence criterion has been violated and stop the calculation:
            if (MCTDHproc.poll() is not None):
                break
            if (violated):
                try:
                    MCTDHproc.terminate()
                except OSError:
                    pass
                break

            # determine the waiting time: twait = 10' + (110'/tfinal) * t
            waitTime = int(600.0 + (6600.0 / tfinal) * currentTime)
            print "Waiting for " + str(waitTime) + " seconds"
            time.sleep(waitTime)



        # check, if the MCTDH calculation finished sucessfully and is converged, or not
        if (not violated):
            converged = True
            break
        else:
            # create a new and improved input file:
            newInputData = currentInputData
            for i in range(len(currentPopList)):
                if (currentPopList[i][3] > natPopThres):
                    currLNr = mlStart + int(currentPopList[i][1][5:])
                    currLn = currentInputData[currLNr]
                    currMd = int(currentPopList[i][2][5:])
                    currSPFId = currLn.index('>') + (2 * currMd)
                    currSPFs = int(currLn[currSPFId])
                    if (currentPopList[i][3] > 2.0 * natPopThres):
                        newInputData[currLNr] = currLn[:currSPFId - 1] + " " + str(currSPFs + 2) + currLn[currSPFId + 1:]
                    else:
                        newInputData[currLNr] = currLn[:currSPFId - 1] + " " + str(currSPFs + 1) + currLn[currSPFId + 1:]

            currentNumber = int(currentNameDir[-3:])
            newNumber = currentNumber + 1
            newNameDir = initBaseName + "_" + str(newNumber).zfill(3)
            newInputName = newNameDir + ".inp"
            newInputFile = open(newInputName, "w")
            for item in newInputData:
                if ("name = " in item and not "opname" in item):
                    item = "    name = " + newNameDir + "\n"
                newInputFile.write(item)
            newInputFile.close()

            currentInputName = newInputName
            currentNameDir = newNameDir
            currentInputData = newInputData


