# utility-scripts
a bunch of scripts that i have come to write during my work in the theoretical chemistry business

In the FCClasses output, the mode analysis block begins with

    "ASSIGNMENTS AND ANALYSIS OF X WITH RESPECT TO Y"

use fccanalyse to extract the strongest vibronic transitions from an already existing
output file:

    fccanalyse -f <number of transitions> <n/e (normal/excited)>

to plot the stick spectrum with gnuplot, including the oscillator excitation labels, use

    plot "assignment-normal.txt" u 2:3 w imp, "" u 2:3:4 with labels rotate left offset 0,char 1


To use specdiff:
Call the program with the data file and the weighting width as the arguments:

    specdiff <file name> <width in eV>

The file must contain three columns representing the energy axis, the ground state spectrum
and the excited state spectrum, respectively

