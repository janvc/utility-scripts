#!/usr/bin/env python3.4

# FCCAnalyse -- An FCClasses analysis script
# this script will extract the strongest stick excitations and the corresponding assignments in
# terms of vibrational modes from an FCClasses output.
#    Copyright 2015 by Jan von Cosel


import sys
import os.path


# define the relevant filenames
convSpecFileName = "fort.18"
assignmentFileName = "fort.21"
specOutFileName = "spectra.dat"
assNormOutFileName = "assignment-normal.txt"
assExcOutFileName = "assignment-excited.txt"

# find out if the program was called with the argument "-f". If so, read in the (existing)
# assignment file(s) and print out only the N strongest excitations
#if (len(sys.argv) != 1 and len(sys.argv) != 4):
#    sys.exit("wrong number of arguments")

if (len(sys.argv) == 4):
    if (sys.argv[1] == "-f" and sys.argv[2].isdigit() and (sys.argv[3] == "n" or sys.argv[3] == "e")):
        nPrint = int(sys.argv[2])
        # first deal with the normal assignments
        if (sys.argv[3] == "n"):
            if (os.path.isfile(assNormOutFileName)):
                readNormAssFile = open(assNormOutFileName, "r")
                readNormAssData = readNormAssFile.readlines()
                readNormAssFile.close()

                for i in range(2, 2 + nPrint):
                    sys.stdout.write(readNormAssData[i])

        # now the excited assignments
        else:
            if (os.path.isfile(assExcOutFileName)):
                readExcAssFile = open(assExcOutFileName, "r")
                readExcAssData = readExcAssFile.readlines()
                readExcAssFile.close()

                for i in range(2, 2 + nPrint):
                    sys.stdout.write(readExcAssData[i])
    else:
        sys.exit("wrong arguments for extracting")

elif (len(sys.argv) == 1):
    # now do the actual analysis of the FCClasses results

    # read the convoluted spectra
    convSpecFile = open(convSpecFileName, "r")
    convSpecData = convSpecFile.readlines()
    convSpecFile.close()

    # first, find out if there even is an excited spectrum
    excitedExists = False
    for i in range(len(convSpecData)):
        if ("total spectrum for mother state no.           2" in convSpecData[i]):
            excitedExists = True

    # find the line where the first total spectrum begins
    i = 0
    while (not "total spectrum for mother state no.           1" in convSpecData[i]):
        i += 1
    i += 2      # adding 2 will get us to the first line that contains numbers
    normSpecStart = i

    # find the line where the second total spectrum begins
    if (excitedExists):
        i = 0
        while (not "total spectrum for mother state no.           2" in convSpecData[i]):
            i += 1
        i += 2      # adding 2 will get us to the first line that contains numbers
        excSpecStart = i

    # get the number of data points in the spectra
    i = 0
    while (convSpecData[normSpecStart + i].rstrip() != ""):
        i += 1
    nPoints = i

    # get the total spectra and write them to the summary file
    specOutFile = open(specOutFileName, "w")
    specOutFile.write("#   Energy [eV]          Int (0K)            Int (hot)\n#\n")
    i = 0
    for i in range(0, nPoints):
        Data = []
        for t in convSpecData[normSpecStart + i].split():
            Data.append(float(t))
        if (excitedExists):
            for t in convSpecData[excSpecStart + i].split():
                Data.append(float(t))
            specOutFile.write("{0:15.5e} {1:20.9e} {2:20.9e}\n".format(Data[0], Data[1], Data[3]))
        else:
            specOutFile.write("{0:15.5e} {1:20.9e}\n".format(Data[0], Data[1]))
    specOutFile.close()


    # read the assignments
    assignmentFile = open(assignmentFileName, "r")
    assignmentData = assignmentFile.readlines()
    assignmentFile.close()

    # get the properties of the stick excitations for the 0K spectrum
    normStickData = []
    i = 0
    for i in range(len(assignmentData)):
        if (assignmentData[i].rstrip() != ""):
            intermediateData = []
            try:
                intermediateData.append(float(assignmentData[i].split()[0]))
                intermediateData.append(float(assignmentData[i].split()[3]))
                intermediateData.append(float(assignmentData[i].split()[6]))
                normStickData.append(intermediateData)
            except ValueError:
                pass
            if ("STATE1 MOTHER STATE N.             2" in assignmentData[i]):
                excAssStart = i
                break

    # sort the list to get the strongest excitations at the top
    normStickDataSorted = sorted(normStickData, key=lambda l:l[2], reverse=True)
    nPeaks = int(input("Enter the number of excitations to be assigned:\n"))

    # get the assignments for the strongest nPeaks excitations
    assNormOutFile = open(assNormOutFileName, "w")
    assNormOutFile.write("#  Index   Energy [eV]       Intensity         {mode, quanta}\n#\n")

    for i in range(0, nPeaks):
        peakIndex = normStickDataSorted[i][0]
        outputString = "{0:7d}".format(int(peakIndex))

        # we have to treat the 0-0 excitation specially, because the format is inconsistent
        if (peakIndex == 1.0):
            outputString += "{0:15.5e} {1:18.9e}".format(normStickDataSorted[i][1], normStickDataSorted[i][2])
            outputString += "     0-0"
        else:
            # find the line of the excitation
            for j in range(len(assignmentData)):
                if (assignmentData[j].rstrip() != ""):
                    if (assignmentData[j].split()[0] == str(peakIndex)):
                        lineNumber = j
                        break

            # append the oscillators
            oscillators = []
            quanta = []
            oscillatorList = assignmentData[lineNumber + 3].split()
            quantumList = assignmentData[lineNumber + 4].split()

            # find out if we have an Oscillator with n > 100. If we do, the list of oscillators
            # will have fewer elements due to missing whitespace
            nHosc = 0       # number of oscillators > 100
            for j in range(len(oscillatorList)):
                if (oscillatorList[j] == "Osc2="):
                    continue
                else:
                    if ("Osc2=" in oscillatorList[j]):
                        nHosc += 1
                        currentOscillator = int(oscillatorList[j].lstrip("Osc2="))
                        currentQuantum = int(quantumList[j + nHosc])
                    else:
                        currentOscillator = int(oscillatorList[j])
                        currentQuantum = int(quantumList[j + nHosc])
                if (currentOscillator == 0):
                    break
                oscillators.append(currentOscillator)
                quanta.append(currentQuantum)
            outputString += "{0:15.5e} {1:18.9e}     ".format(normStickDataSorted[i][1], normStickDataSorted[i][2])
            for j in range(len(oscillators)):
                if (oscillators[j] < 10):
                    outputString += "{0:1d}^{{{1:1d}}}".format(oscillators[j], quanta[j])
                elif (oscillators[j] < 100):
                    outputString += "{0:2d}^{{{1:1d}}}".format(oscillators[j], quanta[j])
                else:
                    outputString += "{0:3d}^{{{1:1d}}}".format(oscillators[j], quanta[j])

                if (j < len(oscillators) - 1):
                    outputString += ","
        outputString += "\n"
        assNormOutFile.write(outputString)

    assNormOutFile.close()

    # do the same thing again for the stick excitations of the excited spectrum
    if (excitedExists):
        excStickData = []
        for i in range(excAssStart, len(assignmentData)):
            if (assignmentData[i].rstrip() != ""):
                intermediateData = []
                try:
                    intermediateData.append(float(assignmentData[i].split()[0]))
                    intermediateData.append(float(assignmentData[i].split()[3]))
                    intermediateData.append(float(assignmentData[i].split()[6]))
                    excStickData.append(intermediateData)
                except ValueError:
                    pass

        # sort the list to get the strongest excitations at the top
        excStickDataSorted = sorted(excStickData, key=lambda l:l[2], reverse=True)

        # get the assignments for the strongest nPeaks excitations
        assExcOutFile = open(assExcOutFileName, "w")
        assExcOutFile.write("#  Index   Energy [eV]       Intensity         {mode, quanta}\n#\n")

        for i in range(nPeaks):
            peakIndex = excStickDataSorted[i][0]
            outputString = "{0:7d}".format(int(peakIndex))

            # find the line of the excitation
            for j in range(len(assignmentData)):
                if (assignmentData[j].rstrip() != ""):
                    if (assignmentData[j].split()[0] == str(peakIndex)):
                        lineNumber = j
                        break
            # take care of the M-0 transition
            if ("state 2 = GROUND" in assignmentData[lineNumber + 5]):
                outputString += "{0:15.5e} {1:18.9e}".format(excStickDataSorted[i][1], excStickDataSorted[i][2])
                outputString += "     M-0"
            else:
                # append the oscillators
                oscillators = []
                quanta = []
                oscillatorList = assignmentData[lineNumber + 6].split()
                quantumList = assignmentData[lineNumber + 7].split()

                nHosc = 0       # number of oscillators > 100
                for j in range(len(oscillatorList)):
                    if (oscillatorList[j] == "Osc2="):
                        continue
                    else:
                        if ("Osc2=" in oscillatorList[j]):
                            nHosc += 1
                            currentOscillator = int(oscillatorList[j].lstrip("Osc2="))
                            currentQuantum = int(quantumList[j + nHosc])
                        else:
                            currentOscillator = int(oscillatorList[j])
                            currentQuantum = int(quantumList[j + nHosc])
                    if (currentOscillator == 0):
                        break
                    oscillators.append(currentOscillator)
                    quanta.append(currentQuantum)
                outputString += "{0:15.5e} {1:18.9e}     ".format(excStickDataSorted[i][1], excStickDataSorted[i][2])
                for j in range(len(oscillators)):
                    if (oscillators[j] < 10):
                        outputString += "{0:1d}^{{{1:1d}}}".format(oscillators[j], quanta[j])
                    elif (oscillators[j] < 100):
                        outputString += "{0:2d}^{{{1:1d}}}".format(oscillators[j], quanta[j])
                    else:
                        outputString += "{0:3d}^{{{1:1d}}}".format(oscillators[j], quanta[j])

                    if (j < len(oscillators) - 1):
                        outputString += ","
            outputString += "\n"
            assExcOutFile.write(outputString)
        assExcOutFile.close()
else:
    sys.exit("wrong number of arguments")

