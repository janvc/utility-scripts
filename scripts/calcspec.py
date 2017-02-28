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
mctdhExe = "/home/jvcosel/mctdh85.5/bin/binary/x86_64/mctdh85"
rdgpopExe = "/home/jvcosel/mctdh85.5/bin/binary/x86_64/rdgpop85"
initInputName = sys.argv[1]
natPopThres = float(sys.argv[2])
optimizeGrid = False
gPopThres = 0.0
if len(sys.argv) == 4:
    optimizeGrid = True
    gPopThres = float(sys.argv[3])


#
# determine the maximum natural population in an ML-MCTDH-calculation
#
def analNatPop(nameDir):
    popFileName = nameDir + "/lownatpop"
    popList = []

    if not os.path.isfile(popFileName):
        emptyList = []
        tmpList = []
        for i in range(4):
            tmpList.append(0)
        emptyList.append(tmpList)
        emptyList.append(tmpList)
        return emptyList

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


def analGpop(nameDir):

    # read the grid population:
    gridPopFile = nameDir + "/gridpop"
    rdgpopCommand = [rdgpopExe, "-f", gridPopFile, "0"]
    stdOut = subprocess.check_output(rdgpopCommand)
    popRawData = stdOut.split('\n')

    # find the line where the actual data starts:
    startLine = 0
    while not "Maximal values" in popRawData[startLine]:
        startLine += 1
    startLine += 2

    # read the grid populations:
    popData = []
    for i in range(startLine, len(popRawData), 1):
        temp = []
        try:
            temp.append(float(popRawData[i].split()[2]))
            temp.append(float(popRawData[i].split()[3]))
            temp.append(float(popRawData[i].split()[5]))
            popData.append(temp)
        except (ValueError, IndexError):
            break

    return popData



def getCurrentTime(nameDir):
    speedFileName = nameDir + "/speed"

    if not os.path.isfile(speedFileName):
        return 0.0

    currentTime = 0.0
    with open(speedFileName) as speedFile:
        speedList = speedFile.readlines()
        if not "#" in speedList[-1]:
            currentTime = float(speedList[-1].split()[0])

    return currentTime

def findMCTDHparameter(dataSet, string):
    for i in range(len(dataSet)):
        if string in dataSet[i]:
            return dataSet[i].split("=")[1]

def findMCTDHline(dataSet, string):
    for i in range(len(dataSet)):
        if string in dataSet[i]:
            return i

def checkNormalTerm(nameDir):
    logFileName = nameDir + "/log"

    with open(logFileName) as logFile:
        for line in logFile:
            if "Program terminated normally." in line:
                return True

    return False

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

    # determine the start of the ML basis section and the pbasis section:
    mlStart = int(findMCTDHline(currentInputData, "mlbasis-section"))
    primStart = int(findMCTDHline(currentInputData, "pbasis-section")) + 1

    # determine the length of the propagation:
    tfinal = float(findMCTDHparameter(initInputData, "tfinal"))


    #
    # this is the big calculation loop
    #
    while True:

        # run the calculation:
        print "Starting calculation " + str(currentNameDir)
        mctdhCommand = [mctdhExe, "-mnd", currentInputName]
        print " Time [fs]   max. natural population    max. grid population    waiting time [s]"

        MCTDHproc = subprocess.Popen(mctdhCommand)

        # wait for five minutes to avoid timing issues:
        time.sleep(300)

        # while the calculation is running,
        # check periodically, if we have violated the convergence criterion:
        while True:

            # analyze the populations:
            currentPopList = analNatPop(currentNameDir)
            gPopList = analGpop(currentNameDir)
            currentTime = getCurrentTime(currentNameDir)

            # determine the waiting time: twait = 10' + (110'/tfinal) * t
            waitTime = int(600.0 + (6600.0 / tfinal) * currentTime)
            print '   {0:.3f}           {1:.5f}                {2:.5f}               {3:5d}'.format(currentTime, currentPopList[0][3], max(max(gPopList)), waitTime)
            sys.stdout.flush()

            violated = False
            if currentPopList[0][3] > natPopThres:
                violated = True
            if optimizeGrid and max(max(gPopList)) > gPopThres:
                violated = True

            # check, if the MCTDH process is still running or has terminated abnormally.
            # If so, issue the 'violated' flag
            if (MCTDHproc.poll() is not None):
                if checkNormalTerm(currentNameDir):
                    if violated:
                        finished = False
                    else:
                        finished = True
                else:
                    finished = False
                break
            else:
                if violated:
                    try:
                        MCTDHproc.terminate()
                    except OSError:
                        pass
                    finished = False
                    break

            time.sleep(waitTime)



        # check, if the MCTDH calculation finished sucessfully and is converged, or not
        if (finished):
            break
        else:
            # create a new and improved input file:
            newInputData = currentInputData

            # first, take care of the SPF basis:
            for i in range(len(currentPopList)):
                if (currentPopList[i][3] > natPopThres):
                    currLNr = mlStart + int(currentPopList[i][1][5:])
                    currLn = currentInputData[currLNr]
                    currMd = int(currentPopList[i][2][5:])

                    currSPFlist = currLn[currLn.index('>') + 1:].split()

                    newInputData[currLNr] = currLn[:currLn.index('>') + 1]
                    for j in range(len(currSPFlist)):
                        if (j == currMd - 1):
                            if (currentPopList[i][3] > 2.0 * natPopThres):
                                newInputData[currLNr] += " " + str(int(currSPFlist[j]) + 2)
                            else:
                                newInputData[currLNr] += " " + str(int(currSPFlist[j]) + 1)
                        else:
                            newInputData[currLNr] += " " + str(int(currSPFlist[j]))
                    newInputData[currLNr] += "\n"

            # now do the primitive basis:
            if optimizeGrid:
                for i in range(len(gPopList)):
                    currLNr = primStart + i
                    currLn = currentInputData[currLNr]
                    currLnData = currLn.split()
                    gridStep = (float(currLnData[5]) - float(currLnData[4])) / float(currLnData[2])
                    newPBline = "    q_" + str(i + 1).zfill(3) + "  ho   "

                    if (gPopList[i][2] > gPopThres):    # basis end
                        newPBline += str(int(currLnData[2]) + 2)
                    elif (gPopList[i][2] > 10.0 * gPopThres):
                        newPBline += str(int(currLnData[2]) + 4)
                    else:
                        newPBline += str(int(currLnData[2]))

                    newPBline += "  xi-xf    "

                    if (gPopList[i][0] > gPopThres):    # grid begin
                        newPBline += "%.1f"%(float(currLnData[4]) - gridStep)
                    elif (gPopList[i][0] > 10.0 * gPopThres):
                        newPBline += "%.1f"%(float(currLnData[4]) - gridStep * 2.0)
                    else:
                        newPBline += "%.1f"%(float(currLnData[4]))

                    newPBline += "   "

                    if (gPopList[i][1] > gPopThres):    # grid end
                        newPBline += "%.1f"%(float(currLnData[5]) + gridStep)
                    elif (gPopList[i][1] > 10.0 * gPopThres):
                        newPBline += "%.1f"%(float(currLnData[5]) + gridStep * 2.0)
                    else:
                        newPBline += "%.1f"%(float(currLnData[5]))

                    newPBline += "\n"
                    newInputData[currLNr] = newPBline

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


