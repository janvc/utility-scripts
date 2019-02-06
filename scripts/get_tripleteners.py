#!/usr/bin/env python3

import sys
import numpy as np

if len(sys.argv) == 1 or "-h" in sys.argv:
    print("get_eners: extract state energies separated by singlet and triplet from TDDFT single points create with get_structures")
    print("")
    print("usage: get_eners <calc basename> <coord1 start> <coord1 npoints> <coord1 stepsize> <coord2 start> <coord2 npoints> <coord2 stepsize> <nstates>")
    print("")
    print("where 'coord1' is the slow coordinate and 'coord2' is the fast coordinate.")
    sys.exit()

eh2ev = 27.21138344

td_baseName = sys.argv[1]
start1 = float(sys.argv[2])
points1 = int(sys.argv[3])
inc1 = float(sys.argv[4])
start2 = float(sys.argv[5])
points2 = int(sys.argv[6])
inc2 = float(sys.argv[7])
nstates = int(sys.argv[8])

nPoints = points1 * points2

stateData = np.zeros((points1, points2, nstates + 1))
multData = np.zeros((points1, points2, nstates + 1))
minEner = 0.0

#
# extract the data from the log files
#
for i in range(points1):        # 1 is the slow coordinate
    for j in range(points2):    # 2 is the fast coordinate
        index = i * points2 + j

        TDfileName = td_baseName + "_" + str(index + 1).zfill(4) + ".log"
        TDfile = open(TDfileName, "r")
        TDdata = TDfile.readlines()
        TDfile.close()

        for line in TDdata:
            if "SCF Done" in line:
                stateData[i, j, 0] = float(line.split()[4])
                multData[i, j, 0] = 0.0
                if stateData[i, j, 0] < minEner:
                    minEner = stateData[i, j, 0]

            if "Excited State" in line:
                stateNum = int(line.split()[2][:-1])
                stateEner = float(line.split()[4])
                stateMult = float(line.split()[9][7:])
                stateData[i, j, stateNum] = stateEner
                multData[i, j, stateNum] = stateMult

#
# shift all ground state energies by the minimum energy
# and convert them to eV
#
for i in range(points1):
    for j in range(points2):
        stateData[i, j, 0] = (stateData[i, j, 0] - minEner) * eh2ev

#
# shift the excitation energies by the ground state energy
#
for i in range(points1):
    for j in range(points2):
        for k in range(1, nstates + 1):
            stateData[i, j, k] = stateData[i, j, k] + stateData[i, j, 0]

# apparently, gaussian does 2D relaxed scans by scanning the fast coordinate,
# incrementing the slow coordinate and then scanning the fast coordinate in the
# reverse direction, meaning that every other scan along the fast coordinate
# needs to be inverted here.
for i in range(points1):
    for j in range(points2):
        # determine if we are on an increasing or
        # decreasing branch of the fast coordinate
        if i % 2 == 0:  # even -> increasing
            jind = j
        else:   # odd -> decreasing
            jind = points2 - j - 1

        x = start1 + (i * inc1)
        y = start2 + (j * inc2)

        sys.stdout.write("{0:8.3f} {1:8.3f}".format(x, y))

        # print only the singlet states:
        for k in range(nstates + 1):
            if abs(multData[i, jind, k]) < 0.1:
                sys.stdout.write("{0:8.4f}".format(stateData[i, jind, k]))

        sys.stdout.write("\n")

sys.stdout.write("\n\n")

# now print the triplet states:
for i in range(points1):
    for j in range(points2):
        # determine if we are on an increasing or
        # decreasing branch of the fast coordinate
        if i % 2 == 0:  # even -> increasing
            jind = j
        else:   # odd -> decreasing
            jind = points2 - j - 1

        x = start1 + (i * inc1)
        y = start2 + (j * inc2)

        sys.stdout.write("{0:8.3f} {1:8.3f}".format(x, y))

        # print only the singlet states:
        for k in range(nstates + 1):
            if abs(abs(multData[i, jind, k]) - 2.0) < 0.1:
                sys.stdout.write("{0:8.4f}".format(stateData[i, jind, k]))

        sys.stdout.write("\n")

