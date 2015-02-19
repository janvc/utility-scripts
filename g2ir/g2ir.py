#!/usr/bin/env python

# G2IR -- This program will create an IR stick spectrum from a Gaussian09 output file
# with a normal mode calculation
#     Copyright 2015 by Jan von Cosel

import sys  # operating system stuff
import re   # regular expressions


FreqFileName = sys.argv[1]
specOutFileName = sys.argv[2]

# find the last line where the normal modes start
# (if the calculation used 'hpmodes', the modes will be listed twice and we will
# only use the second list - the not-hpmodes one)
StartLine = 0
with open(FreqFileName) as FreqFile:
    for num, line in enumerate(FreqFile, 1):
        if ("Harmonic frequencies" in line):
            StartLine = num

if (StartLine == 0):
    sys.exit("no frequencies found.")

# read the part of the file containing the normal modes
nmData = []
with open(FreqFileName) as FreqFile:
    for i in xrange(StartLine):
        FreqFile.next()
    for line in FreqFile:
        nmData.append(line.strip())

# read the data from the raw array into the frequency and intensity array
# !!!!! WARNING !!!!!  not sure what happens if the number of normal modes
#                      is not divisable by 3 (incomplete lines)
frequencyData = []
intensityData = []
for i in range(len(nmData)):
    if ("Frequencies" in nmData[i]):
        if ("IR Inten" in nmData[i+3]):
            freqLine = nmData[i]
            intLine = nmData[i+3]
            frequencyData.append(float(freqLine.split()[2]))
            frequencyData.append(float(freqLine.split()[3]))
            frequencyData.append(float(freqLine.split()[4]))
            intensityData.append(float(intLine.split()[3]))
            intensityData.append(float(intLine.split()[4]))
            intensityData.append(float(intLine.split()[5]))

specOutFile = open(specOutFileName, "w")
specOutFile.write("# wavenumber [cm-1]   Intensity\n#\n")
for i in range(len(frequencyData)):
    specOutFile.write("{0:14.4f} {1:14.4f}\n".format(frequencyData[i], intensityData[i]))

specOutFile.close()

