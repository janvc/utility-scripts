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

# create wavefunctions for specific vibrational states and use them to
# calculate cross-correlation functions for IVR analysis
#   ----- RPTWF version -----

import sys
import re
import multiprocessing as mp
import subprocess
import os
import shutil


mctdhExe = "/home/jvcosel/mctdh85.5/bin/binary/x86_64/mctdh85"
corrExe = "/home/jvcosel/mctdh85.5/bin/binary/x86_64/crosscorr85"

def main():
    nameDir = sys.argv[1]
    refSLInput = sys.argv[2]
    refMLInput = sys.argv[3]
    prtThres = float(sys.argv[4])

    refSLFile = open(refSLInput, "r")
    refMLFile = open(refMLInput, "r")
    refSLData = refSLFile.readlines()
    refMLData = refMLFile.readlines()
    refSLFile.close()
    refMLFile.close()

    # get the dimensionality of the system from the propagation:
    with open(nameDir + "/op.log") as opFile:
        for line in opFile:
            if "ndof   =" in line:
                nDof = int(line.split()[2])

    if len(sys.argv) == 6 and sys.argv[5] == "-c":
        create_corrs(nDof, refSLData, refMLData, nameDir, refMLInput)

    # create a gnuplot file to show the correlation functions:
    gnuplotName = "thermcorrelations.plt"
    gnuplotFile = open(gnuplotName, "w")
    gnuplotFile.write("plot 'thermcorr_gs.dat' u 1:4 w l")

    for i in range(nDof):
        for j in range(2,6):
            currCorrName = "thermcorr_" + str(i+1).zfill(3) + "_" + str(j) + ".dat"
            maxCorrVal = 0.0
            with open(currCorrName) as corrFile:
                for line in corrFile:
                    if not "#" in line:
                        if float(line.split()[3]) > maxCorrVal:
                            maxCorrVal = float(line.split()[3])

            if maxCorrVal > prtThres:
                writeString = ", '" + currCorrName + "' u 1:4 w l"
                gnuplotFile.write(writeString)

    gnuplotFile.write("\n")
    gnuplotFile.close()

def create_corrs(nModes, refSLData, refMLData, nameDir, refMLInput):
    pool = mp.Pool(processes=4)
    for i in range(nModes):
        for j in range(2,6):
            pool.apply_async(func=run_calc, args=(i+1, j, refSLData, refMLData, nameDir))

    pool.close()
    pool.join()

    # do the calculation with the global ground state:
    refMLDir = os.path.splitext(refMLInput)[0]
    MLgencommand = [mctdhExe, "-mnd", refMLInput]
    subprocess.check_call(MLgencommand)
    shutil.copy2(refMLDir + "/restart", "calc_gs.rst")
    corrcommand = [corrExe, "-f", nameDir + "/psi", "-o", "thermcorr_gs.dat", "-r", "calc_gs.rst"]
    shutil.rmtree(refMLDir)
    subprocess.check_call(corrcommand)
    os.remove("calc_gs.rst")


def run_calc(mode, state, refSLData, refMLData, psiDir):
    newSLData = refSLData[:]
    newMLData = refMLData[:]
    # get the name-directory for the reference calculations:
    for i in range(len(refSLData)):
        if ("name" in refSLData[i] and "=" in refSLData[i] and not "opname" in refSLData[i]):
            dirLine = i
        if "file" in refSLData[i] and "=" in refSLData[i]:
            excLine = i

    baseName = "thermcalc_" + str(mode).zfill(3) + "_" + str(state)
    corrName = "thermcorr_" + str(mode).zfill(3) + "_" + str(state) + ".dat"
    SLinputFileName = baseName + "_sl.inp"
    MLinputFileName = baseName + "_ml.inp"
    SLinputWF = baseName + "_sl.rst"
    MLinputWF = baseName + "_ml.rst"

    newSLData[dirLine] = "    name = " + baseName + "\n"
    excString = "    operate = excite_" + str(mode).zfill(3) + "\n"
    for i in range(state-1):
        newSLData.insert(excLine + 1,excString)

    SLinputFile = open(SLinputFileName, "w")
    for item in newSLData:
        SLinputFile.write(item)
    SLinputFile.close()

    os.mkdir(baseName)
    SLgencommand = [mctdhExe, "-w", SLinputFileName]
    subprocess.check_call(SLgencommand)
    shutil.copy2(baseName + "/restart", SLinputWF)

    for i in range(len(refMLData)):
        if "file" in refMLData[i] and "=" in refMLData[i]:
            rstLine = i
            break
    newMLData[dirLine] = "    name = " + baseName + "\n"
    newMLData[rstLine] = "    file = " + SLinputWF + "\n"

    MLinputFile = open(MLinputFileName, "w")
    for item in newMLData:
        MLinputFile.write(item)
    MLinputFile.close()

    MLgencommand = [mctdhExe, "-w", MLinputFileName]
    subprocess.check_call(MLgencommand)
    shutil.copy2(baseName + "/restart", MLinputWF)
    shutil.rmtree(baseName)

    corrcommand = [corrExe, "-f", psiDir + "/psi", "-o", corrName, "-r", MLinputWF]
    subprocess.check_call(corrcommand)

    os.remove(SLinputWF)
    os.remove(MLinputWF)
    os.remove(SLinputFileName)
    os.remove(MLinputFileName)

if __name__ == "__main__":
    main()

