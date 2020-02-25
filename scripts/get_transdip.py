#!/usr/bin/env python3

import argparse

helpString = ""
helpString += "Computes the derivatives of the transition dipole moment of an "
helpString += "electronic one-photon transition with respect to the atomic "
helpString += "coordinates based on a numerical excited-state frequency "
helpString += "calculation"

parser = argparse.ArgumentParser(description=helpString)
parser.add_argument("-f", "--filename", nargs=1)
parser.add_argument("-s", "--state", type=int, nargs=1)

args = parser.parse_args()

print(args.filename)
print(args.state)

