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


# Generate an ASCII file containing the VIPER pulse sequence
# suitable for use in an MCTDH calculation.
# The IR excitation pulse is a single-sided exponential with a gaussian
# of FWHM of about 150 fs on the rising edge. The other two pulses are
# gaussian shaped

import sys
import numpy as np

fs2au = 41.34137333656
Eh2wn = 219474.6312

if len(sys.argv) == 2 and sys.argv[1] == "-h":
    print("Usage: genviper sequence tau1 tau2 tau3 pw1 pw2 pw3 freq1 freq2 freq3 a1 a2 a3 tfinal step outputfile")
    print("with times tauX, pwX, tfinal and step in fs")
    print("and the frequencies freqX in cm-1")
    sys.exit()

a_seq = int(sys.argv[1])
a_tau1 = float(sys.argv[2])
a_tau2 = float(sys.argv[3])
a_tau3 = float(sys.argv[4])
a_pw1 = float(sys.argv[5])
a_pw2 = float(sys.argv[6])
a_pw3 = float(sys.argv[7])
a_freq1 = float(sys.argv[8])
a_freq2 = float(sys.argv[9])
a_freq3 = float(sys.argv[10])
a_a1 = float(sys.argv[11])
a_a2 = float(sys.argv[12])
a_a3 = float(sys.argv[13])
a_tfinal = float(sys.argv[14])
a_step = float(sys.argv[15])
a_filename = sys.argv[16]


nsteps = int(a_tfinal / a_step) + 1
sig_1 = a_pw1 * fs2au / (2.0 * np.sqrt(2.0 * np.log(2.0)))
sig_2 = a_pw2 * fs2au / (2.0 * np.sqrt(2.0 * np.log(2.0)))
sig_3 = a_pw3 * fs2au / (2.0 * np.sqrt(2.0 * np.log(2.0)))
tau_1 = a_tau1 * fs2au
tau_2 = a_tau2 * fs2au
tau_3 = a_tau3 * fs2au
frq_1 = a_freq1 / Eh2wn
frq_2 = a_freq2 / Eh2wn
frq_3 = a_freq3 / Eh2wn
dt = a_step * fs2au
tfinal = a_tfinal * fs2au

# sequence encoding
#
#  no.     pulse 1    pulse 2    pulse 3
#
#   1        x
#   2                   x
#   3        x          x
#   4                              x
#   5        x                     x
#   6                   x          x
#   7        x          x          x
#

outfile = open(a_filename, "w")
for i in range(nsteps):
    time = i * dt
    value = 0.0
    if a_seq == 1 or a_seq == 3 or a_seq == 5 or a_seq == 7:
        if time < tau_1:
            value += np.exp(-np.power(time - tau_1, 2) / (13869679.2)) * np.sin(frq_1 * time)
        else:
            value += np.exp(-(time - tau_1) / sig_1) * np.sin(frq_1 * time)

    if a_seq == 2 or a_seq == 3 or a_seq == 6 or a_seq == 7:
        value += np.exp(-np.power(time - tau_2, 2) / (2.0 * sig_2 * sig_2)) * np.sin(frq_2 * time)

    if a_seq == 4 or a_seq == 5 or a_seq == 6 or a_seq == 7:
        value += np.exp(-np.power(time - tau_3, 2) / (2.0 * sig_3 * sig_3)) * np.sin(frq_3 * time)

    if abs(value) < 1.0e-8:
        value = 0.0

    outfile.write("{0:12.8f} {1:12.8f} {2:12.8f} {3:12.8f}\n".format(time, value, 0.0, value))

outfile.close()







