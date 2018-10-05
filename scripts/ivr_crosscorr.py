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

import sys
import re
import multiprocessing as mp
import subprocess
import os
import shutil


mctdhExe = "/home/jvcosel/mctdh85.7/bin/binary/x86_64/mctdh85"
corrExe = "/home/jvcosel/mctdh85.7/bin/binary/x86_64/crosscorr85"

def main():
    nameDir = sys.argv[1]
    refInput = sys.argv[2]
    prtThres = float(sys.argv[3])

    refFile = open(refInput, "r")
    refData = refFile.readlines()
    refFile.close()

    # get the dimensionality of the system from the propagation:
    with open(nameDir + "/op.log") as opFile:
        for line in opFile:
            if "ndof   =" in line:
                nDof = int(line.split()[2])

    if len(sys.argv) == 5 and sys.argv[4] == "-c":
        create_corrs(nDof, refData, nameDir, refInput)

    # create a gnuplot file to show the correlation functions:
    gnuplotName = "correlations.plt"
    gnuplotFile = open(gnuplotName, "w")
    gnuplotFile.write("plot 'corr_gs.dat' u 1:4 w l")

    for i in range(nDof):
        for j in range(2,6):
            currCorrName = "corr_" + str(i+1).zfill(3) + "_" + str(j) + ".dat"
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

def create_corrs(nModes, refData, nameDir, refInput):
    pool = mp.Pool(processes=4)
    for i in range(nModes):
        for j in range(2,6):
            pool.apply_async(func=run_calc, args=(i+1, j, refData, nameDir))

    pool.close()
    pool.join()

    # do the calculation with the global ground state:
    refDir = os.path.splitext(refInput)[0]
    gencommand = [mctdhExe, "-mnd", refInput]
    corrcommand = [corrExe, "-f", nameDir + "/psi", "-o", "corr_gs.dat", "-r", "calc_gs.rst"]
    subprocess.check_call(gencommand)
    shutil.copy2(refDir + "/restart", "calc_gs.rst")
    shutil.rmtree(refDir)
    subprocess.check_call(corrcommand)
    os.remove("calc_gs.rst")


def run_calc(mode, state, refData, psiDir):
    newData = refData
    # get the name-directory for the reference calculations:
    for i in range(len(refData)):
        if ("name" in refData[i] and "=" in refData[i] and not "opname" in refData[i]):
            dirLine = i
        if "Eq_" + str(mode).zfill(3) in refData[i]:
            modeLine = i

    baseName = "calc_" + str(mode).zfill(3) + "_" + str(state)
    corrName = "corr_" + str(mode).zfill(3) + "_" + str(state) + ".dat"
    inputFileName = baseName + ".inp"
    inputWF = baseName + ".rst"
    gencommand = [mctdhExe, "-w", inputFileName]
    corrcommand = [corrExe, "-f", psiDir + "/psi", "-o", corrName, "-r", inputWF]

    newData[dirLine] = "    name = " + baseName + "\n"
    newData[modeLine] = "        q_" + str(mode).zfill(3) + "  eigenf  Eq_" + str(mode).zfill(3) + "  pop = " + str(state) + "\n"


    inputFile = open(inputFileName, "w")
    for item in newData:
        inputFile.write(item)
    inputFile.close()

    os.mkdir(baseName)
    subprocess.check_call(gencommand)

    shutil.copy2(baseName + "/restart", inputWF)
    shutil.rmtree(baseName)

    subprocess.check_call(corrcommand)

    os.remove(inputWF)
    os.remove(inputFileName)




if __name__ == "__main__":
    main()

