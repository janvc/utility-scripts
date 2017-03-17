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

import random
import sys
import math

psiFileName = sys.argv[1]
spfLen = int(sys.argv[2])

psiFile = open(psiFileName, "r")
psiData = psiFile.readlines()
psiFile.close()

for i in range(spfLen):
    print(psiData[i], end="")

for i in range(spfLen, len(psiData)):
    rp = float(psiData[i].split(",")[0][2:])
    ip = float(psiData[i].split(",")[1][:-2])
    
    phase = random.random() * 2.0 * math.pi

    rpn = rp * math.cos(phase) - ip * math.sin(phase)
    ipn = rp * math.sin(phase) + ip * math.cos(phase)

    print("(", rpn, ",", ipn, ")")


