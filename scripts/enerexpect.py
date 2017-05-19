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


import subprocess
import matplotlib.pyplot as plt
import numpy as np
import math
import sys


oplogFileName = "op.log"
expectProg = "/home/jvcosel/mctdh85.5/bin/binary/x86_64/expect85"

class ModePlotter:
    def __init__(self):
        # determine number of modes:
        self.Nmodes = 0
        with open(oplogFileName) as oplogFile:
            for line in oplogFile:
                if "ndof   =" in line:
                    self.Nmodes = int(line.split()[2])

        # calculate the expectation values:
        self.values = []
        self.deltas = []
        for i in range(self.Nmodes):
            print('processing mode %d of %d'%(i + 1, self.Nmodes))
            operName = "Eq_" + str(i + 1).zfill(3)
            expectCommand = [expectProg, "-o", "no", operName]
            stdOut = subprocess.check_output(expectCommand)
            tempList = stdOut.split()
            self.Nsteps = len(stdOut.split('\n')) - 16
            self.Xaxis = []
            Ytmp = []
            for j in range(self.Nsteps):
                self.Xaxis.append(float(tempList[31 + 4 * j]))
                Ytmp.append(float(tempList[32 + 4 * j]))
            self.values.append(Ytmp)

            # shift all curves to begin at zero energy
            dYtmp = []
            for j in range(self.Nsteps):
                dYtmp.append(Ytmp[j] - Ytmp[0])
            self.deltas.append(dYtmp)

    def onclick(self, event, modePlot):
        time = event.xdata
        ener = event.ydata

        print('x = %f   y = %f'%(event.xdata, event.ydata))

        # determine the pseudo-eucledian distance of the clicked point to all data points
        # and select the point with the lowest distance:
        dMin = math.sqrt((self.Xaxis[0] - time) * (self.Xaxis[0] - time) + (self.deltas[0][0] - ener) * (self.deltas[0][0] - ener))
        mode = 0
        for i in range(self.Nmodes):
            for j in range(self.Nsteps):
                dist = math.sqrt((self.Xaxis[j] - time) * (self.Xaxis[j] - time) + (self.deltas[i][j] - ener) * (self.deltas[i][j] - ener))
                if dist < dMin:
                    dMin = dist
                    mode = i

        print('mode %d'%(mode + 1))

    def trace(self, modeIndex):
        values = []
        values.append(self.Xaxis)
        values.append(self.deltas[modeIndex])
        return values

    def modes(self):
        return self.Nmodes

    def data(self):
        return self.deltas

    def times(self):
        return self.Xaxis


if __name__ == "__main__":
    mp = ModePlotter()

    timeAxis = mp.times()
    modeData = mp.data()

    # write the energies to the output file:
    dataFile = open("ener.log", "w")
    dataFile.write("# Time [fs]")
    for i in range(len(modeData)):
        dataFile.write("   mode {0:3d}  ".format(i + 1))
    dataFile.write("\n")
    for i in range(len(modeData[0])):
        dataFile.write("{0:12.5e}".format(timeAxis[i]))
        for j in range(len(modeData)):
            dataFile.write("{0:13.5e}".format(modeData[j][i]))
        dataFile.write("\n")
    dataFile.close()


    if len(sys.argv) == 1:
        fig = plt.figure()
        ax = fig.add_subplot(111)

        for i in range(mp.modes()):
            data = mp.trace(i)
            ax.plot(data[0], data[1],)

        fig.canvas.mpl_connect('button_press_event', lambda event: mp.onclick(event, ax))

        plt.show()

