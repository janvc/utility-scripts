#!/usr/bin/env python3

# Copyright 2016-2018 Jan von Cosel
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
import shutil
import copy
import getopt


#
# Define some global variables.
#
mctdhExe = "/home/jvcosel/mctdh85.7/bin/binary/x86_64/mctdh85"
rdgpopExe = "/home/jvcosel/mctdh85.7/bin/binary/x86_64/rdgpop85"


#
# Find the value of an MCTDH-type "key = value"-pair.
#
def findMCTDHparameter(dataSet, string):
    for i in range(len(dataSet)):
        if string in dataSet[i].lower():
            return dataSet[i].split("=")[1]


#
# Determine the line number of the first occurrence of the
# string "string" in the "dataSet".
#
def findMCTDHline(dataSet, string):
    for i in range(len(dataSet)):
        if string in dataSet[i].lower():
            return i


#
# Determine the maximum natural population in an ML-MCTDH-calculation
# in the name directory "nameDir".
#
def analNatPop(nameDir):
    popFileName = nameDir + "/lownatpop"
    popList = []

    # if the lownatpop file does not exist, return an empty list
    if not os.path.isfile(popFileName):
        emptyList = []
        tmpList = []
        for i in range(4):
            tmpList.append(0)
        emptyList.append(tmpList)
        emptyList.append(tmpList)
        return emptyList

    fileHasLines = False
    with open(popFileName) as popFile:
        for line in popFile:
            if "#" in line:     # the layer/node/mode information
                fileHasLines = True
                # Rule out the first line:
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

    # if the file exists but is empty, return an empty list as well
    if not fileHasLines:
        emptyList = []
        tmpList = []
        for i in range(4):
            tmpList.append(0)
        emptyList.append(tmpList)
        emptyList.append(tmpList)
        return emptyList

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


#
# Determine the population of the primitive basis of an (ML-)MCTDH calculation
# in the name directory "nameDir" at the ends of the grid as well as the basis.
#
def analGpop(nameDir):
    gridPopFile = nameDir + "/gridpop"

    if not os.path.isfile(gridPopFile):
        emptyList = []
        tmpList = []
        for i in range(4):
            tmpList.append(0)
        emptyList.append(tmpList)
        emptyList.append(tmpList)
        return emptyList

    # Read the grid population:
    rdgpopCommand = [rdgpopExe, "-f", gridPopFile, "0"]
    stdOut = subprocess.check_output(rdgpopCommand)
    stdOut = stdOut.decode("utf-8")
    popRawData = stdOut.split('\n')

    # Find the line where the actual data starts:
    startLine = 0
    while not "Maximal values" in popRawData[startLine]:
        startLine += 1
    startLine += 2

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


#
# Determine the current propagation time of a running
# MCTDH calculation in the name directory "nameDir".
#
def getCurrentTime(nameDir):
    speedFileName = nameDir + "/speed"

    if not os.path.isfile(speedFileName):
        return 0.0

    currentTime = 0.0
    with open(speedFileName) as speedFile:
        speedList = speedFile.readlines()
        try:
            if not "#" in speedList[-1]:
                currentTime = float(speedList[-1].split()[0])
        # if the speed file exists, but is (still) empty, we would get an index error
        except IndexError:
            return 0.0

    return currentTime


#
# Check, whether the MCTDH calculation in the name directory
# "nameDir" has terminated normally.
#
def checkNormalTerm(nameDir):
    logFileName = nameDir + "/log"

    with open(logFileName) as logFile:
        for line in logFile:
            if "Program terminated normally." in line:
                return True

    return False


#
# Contain the important informaation about the MCTDH calculation.
#
class MCTDHcalculation:
    def __init__(self):
        self.initInputName = sys.argv[1]
        self.baseName = os.path.splitext(self.initInputName)[0]
        self.natPopThres = float(sys.argv[2])
        self.optGrid = False
        self.grdPopThres = 0.0
        self.checkInterval = 60.0
        if len(sys.argv) >= 4:
            if float(sys.argv[3]) > 0.0:
                self.optGrid = True
                self.grdPopThres = float(sys.argv[3])
        if len(sys.argv) == 5:
            self.checkInterval = float(sys.argv[4])

        f = open(self.initInputName, "r")
        self.initInputData = f.readlines()
        f.close()

        print("Processing MCTDH input file ", self.initInputName)
        sys.stdout.write("Threshold for ML basis population: {0:3.1e}\n".format(self.natPopThres))
        sys.stdout.write("Threshold for grid population:     {0:3.1e}\n".format(self.grdPopThres))
        sys.stdout.write("Checking the populations every {0:5.1f} seconds\n".format(self.checkInterval))

        self.TFinLNr = int(findMCTDHline(self.initInputData, "tfinal"))
        self.tfinal = float(findMCTDHparameter(self.initInputData, "tfinal"))
        self.mlbStart = int(findMCTDHline(self.initInputData, "mlbasis-section")) + 1
        self.mlbEnd = int(findMCTDHline(self.initInputData, "end-mlbasis-section")) - 1
        self.grdStart = int(findMCTDHline(self.initInputData, "pbasis-section")) + 1
        self.grdEnd = int(findMCTDHline(self.initInputData, "end-pbasis-section")) - 1
        self.iwfStart = int(findMCTDHline(self.initInputData, "init_wf-section")) + 1
        self.iwfEnd = int(findMCTDHline(self.initInputData, "end-init_wf-section")) - 1

        self.finished = False
        self.crashed = False

def main():

    print("This is calcspec!")

    calc = MCTDHcalculation()
    calcIndex = 1

    tmpNameDir = calc.baseName + "_" + str(calcIndex).zfill(3)
    tmpInputName = tmpNameDir + ".inp"
    startInputData = calc.initInputData[:]
    for i in range(len(startInputData)):
        if ("name" in startInputData[i].lower() and "=" in startInputData[i]):
            s = startInputData[i].split("=")[0] + " = " + tmpNameDir + "\n"
            startInputData[i] = s
            break
    currentInputData = startInputData[:]

    currentInputFile = open(tmpInputName, "w")
    for l in currentInputData:
        currentInputFile.write(l)
    currentInputFile.close()

    while True:
        sys.stdout.write("Starting calculation no {0:4d}\n".format(calcIndex))
        sys.stdout.flush()

        mctdhCommand = [mctdhExe, "-mnd", tmpInputName]
        mctdhProc = subprocess.Popen(mctdhCommand)

        while True:
            time.sleep(calc.checkInterval)
            # Get info about the state of the calculation.
            natPopList = analNatPop(tmpNameDir)
            grdPopList = analGpop(tmpNameDir)
            currTime = getCurrentTime(tmpNameDir)
            sys.stdout.write("Time = {0:7.2f} fs, max. nat. pop. = {1:7.5f}, max. grid pop. = {2:8.6f}\n".format(currTime, natPopList[0][3], max(max(grdPopList))))

            # See if we exceeded the convergence thresholds.
            natViol = False
            grdViol = False
            if natPopList[0][3] > calc.natPopThres:
                natViol = True
            if max(max(grdPopList)) > calc.grdPopThres and calc.optGrid:
                grdViol = True

            # Terminate the calculation if needed.
            if mctdhProc.poll() is not None:    # If true, the process is no longer running.
                if checkNormalTerm(tmpNameDir): # The calculation has finished
                    if natViol or grdViol:      # but we are still not converged yet.
                        calc.finished = False
                    else:
                        calc.finished = True
                else:
                    # This branch is only reached if MCTDH is no longer running
                    # but has not terminated normally, i.e. it crashed...
                    calc.finished = False
                    sys.exit("MCTDH has crashed")
                break
            else:
                if natViol or grdViol:
                    calc.finished = False
                    if os.path.isfile(tmpNameDir + "/stop"):    # Remove the stop file and wait
                        os.remove(tmpNameDir + "/stop")         # for MCTDH to terminate itself.
                        print("removing stop file")
                        while True:
                            # This is a potential infinite loop...
                            if mctdhProc.poll() is not None:
                                break
                            else:
                                time.sleep(10)
                    else:
                        # This is also a weird place: why would MCTDH still be
                        # running without there being a stop file?
                        print("terminating MCTDH")
                        try:
                            mctdhProc.terminate()
                        except OSError:
                            pass
                    break

        if calc.finished:
            break

        # We have just finished the calculation. Since we have removed the stop file, there is one more
        # write step to consider, where the thresholds could have been exceeded. Therefore, run
        # analNatPop and analGrdPop again.
        # If the natural populations have exceeded the threshold, improve the ML basis:
        natPopList = analNatPop(tmpNameDir)
        grdPopList = analGpop(tmpNameDir)
        newInputData = currentInputData[:]
        if natViol:
            print("improving the ML basis.")

            # Increase the number of SPFs for the modes where the convergence criterion
            # was violated. If the natural population is greater than twice the
            # threshold (i.e. the basis is REALLY bad) then add two SPFs.
            for i in range(len(natPopList)):
                if natPopList[i][3] > calc.natPopThres:
                    currLNr = calc.mlbStart + int(natPopList[i][1][5:]) - 1
                    currLn = currentInputData[currLNr]
                    currMode = int(natPopList[i][2][5:])
                    currSPFlist = currLn[currLn.index('>') + 1:].split()

                    newInputData[currLNr] = currLn[:currLn.index('>') + 1]
                    for j in range(len(currSPFlist)):
                        if j == currMode - 1:
                            if natPopList[i][3] > 2.0 * calc.natPopThres:
                                newInputData[currLNr] += " " + str(int(currSPFlist[j]) + 2)
                                sys.stdout.write("line {0:4d}, mode {1:4d}: {2:4d} > {3:4d}\n".format(int(natPopList[i][1][5:]) - 1, currMode, int(currSPFlist[j]), int(currSPFlist[j]) + 2))
                            else:
                                newInputData[currLNr] += " " + str(int(currSPFlist[j]) + 1)
                                sys.stdout.write("line {0:4d}, mode {1:4d}: {2:4d} > {3:4d}\n".format(int(natPopList[i][1][5:]) - 1, currMode, int(currSPFlist[j]), int(currSPFlist[j]) + 1))
                        else:
                            newInputData[currLNr] += " " + str(int(currSPFlist[j]))
                    newInputData[currLNr] += "\n"
                sys.stdout.flush()

        # If the primitive basis needs to be expanded, we need to start the calculation
        # from the beginning.
        if grdViol:

            print("improving the primitive basis.")
            # Improve the primitive basis section.
            for i in range(len(grdPopList)):
                currLNr = calc.grdStart + i
                currLn = currentInputData[currLNr]
                currLnData = currLn.split()
                grdStep = (float(currLnData[5]) - float(currLnData[4])) / float(currLnData[2])

                newPBline = "    " + currLnData[0] + "  ho   "

                if grdPopList[i][2] > calc.grdPopThres: # basis end too high
                    newPBline += str(int(currLnData[2]) + 2)
                    sys.stdout.write("mode {0:4d} add 2 basis functions\n".format(i+1))
                elif grdPopList[i][2] > 10.0 * calc.grdPopThres:
                    newPBline += str(int(currLnData[2]) + 4)
                    sys.stdout.write("mode {0:4d} add 4 basis functions\n".format(i+1))
                else:
                    newPBline += str(int(currLnData[2]))

                newPBline += "  xi-xf    "

                if grdPopList[i][0] > calc.grdPopThres:
                    newPBline += "%.1f"%(float(currLnData[4]) - grdStep)
                    sys.stdout.write("mode {0:4d} move grid begin from {1:6.2f} to {2:6.2f}\n".format(i+1, float(currLnData[4]), float(currLnData[4]) - grdStep))
                elif grdPopList[i][0] > 10.0 * calc.grdPopThres:
                    newPBline += "%.1f"%(float(currLnData[4]) - grdStep * 2.0)
                    sys.stdout.write("mode {0:4d} move grid begin from {1:6.2f} to {2:6.2f}\n".format(i+1, float(currLnData[4]), float(currLnData[4]) - grdStep * 2.0))
                else:
                    newPBline += "%.1f"%(float(currLnData[4]))

                newPBline += "    "

                if grdPopList[i][1] > calc.grdPopThres:
                    newPBline += "%.1f"%(float(currLnData[5]) + grdStep)
                    sys.stdout.write("mode {0:4d} move grid end from {1:6.2f} to {2:6.2f}\n".format(i+1, float(currLnData[5]), float(currLnData[5]) + grdStep))
                elif grdPopList[i][1] > 10.0 * calc.grdPopThres:
                    newPBline += "%.1f"%(float(currLnData[5]) + grdStep * 2.0)
                    sys.stdout.write("mode {0:4d} move grid end from {1:6.2f} to {2:6.2f}\n".format(i+1, float(currLnData[5]), float(currLnData[5]) + grdStep * 2.0))
                else:
                    newPBline += "%.1f"%(float(currLnData[5]))

                newPBline += "\n"
                newInputData[currLNr] = newPBline
                sys.stdout.flush()

            # Take the initial input and substitute the newly generated ml-
            # and primitive bases into it.
            currentInputData = startInputData[:]
            currentInputData[calc.mlbStart:calc.mlbEnd] = newInputData[calc.mlbStart:calc.mlbEnd]
            currentInputData[calc.grdStart:calc.grdEnd] = newInputData[calc.grdStart:calc.grdEnd]

        # create the new input file for the next calculation.
        calcIndex += 1
        newNameDir = calc.baseName + "_" + str(calcIndex).zfill(3)
        newInputName = newNameDir + ".inp"
        newInputFile = open(newInputName, "w")
        for item in newInputData:
            if ("name" in item.lower() and "=" in item and not "opname" in item.lower()):
                s = item.split("=")[0] + " = " + newNameDir + "\n"
                item = s
            newInputFile.write(item)
        newInputFile.close()

        tmpInputName = newInputName
        tmpNameDir = newNameDir
        currentInputData = newInputData[:]

    print("calculation completed.")
    sys.stdout.flush()

    print("Normal termination!")

if __name__ == "__main__":
    main()

