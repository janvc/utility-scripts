# utility-scripts
a bunch of scripts that i have come to write during my work in the theoretical chemistry business

use fccanalyse to extract the strongest vibronic transitions from an already existing
output file:

    fccanalyse -f <number of transitions> <n/e (normal/excited)

to plot the stick spectrum with gnuplot, including the oscillator excitation labels, use

    plot "assignment-normal.txt" u 2:3 w imp, "" u 2:3:4 with labels rotate left offset 0,char 1

