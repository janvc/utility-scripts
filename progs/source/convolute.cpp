/*
 * Copyright 2015 Jan von Cosel
 *
 * This file is part of utility-scripts.
 *
 * utility-scripts is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * utility-scripts is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have recieved a copy of the GNU General Public License
 * along with utility-scripts. If not, see <http://www.gnu.org/licenses/>.
 *
 */


#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <boost/program_options.hpp>

/*
 * call this program like:
 *
 * ./convolute <type argument> <input file> <FWHM in eV> <# of points>
 *
 * type argument:     '-f' for fcclasses 'fort.19'-type file
 *                    '-g' for gaussian log file with vibrational modes
 * input file:        name of the input file
 * FWHM:              width of the peaks (eV for fcclasses, cm-1 for gaussian)
 * # of points:       number of points in the output data
 */
int main(int argc, char *argv[])
{
	/*
	 * Define run parameters and process command options
	 */
	std::string intFileName;
	enum type
	{
		gaussian,
		fcclasses
	};
	type specType;
	char typeChar;
	double gamma;
	long nPoints;

	namespace po = boost::program_options;
	po::options_description desc("Allowed options");
	desc.add_options()
		("help,h", "produce this help message")
		("type,t", po::value<char>(&typeChar)->required(), "format of the input file: g - Gaussian, f - FCClasses")
		("input,f", po::value<std::string>(&intFileName)->required(), "the input file")
		("width,w", po::value<double>(&gamma)->required(), "FWHM of the spectrum (eV for FCClasses, cm-1 for Gaussian")
		("npoints,n", po::value<long>(&nPoints)->required(), "Number of points in the spectrum")
	;
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	if (vm.count("help"))
	{
		std::cout << desc << std::endl;
		return 1;
	}
	po::notify(vm);
	if (typeChar == 'f')
		specType = fcclasses;
	else if (typeChar == 'g')
		specType = gaussian;
	else
	{
		std::cerr << "Invalid value for spectrum type.\n\n";
		std::cout << desc << std::endl;
		return 1;
	}

	/*
	 * Open the input file
	 */
	std::ifstream intFile(intFileName, std::ifstream::in);
	if (!intFile.good())
	{
		std::cerr << "ERROR: File " << intFileName << " could not be opened.\n";
		return 2;
	}

	/*
	 * read the peaks from the input file:
	 */
	std::vector<double> peakPositions;
	std::vector<double> peakIntensities;
	double tmpPos, tmpInt;

	if (specType == fcclasses)	// fcclasses type output files
		while (intFile >> tmpPos >> tmpInt)
			if (tmpInt != 0)
			{
				double nextPos, nextInt;
				intFile >> nextPos >> nextInt;

				double actualPosition = (tmpPos + nextPos) / 2.0;

				peakPositions.push_back(actualPosition);
				peakIntensities.push_back(tmpInt);
			}
	if (specType == gaussian)	// gaussian type output files
	{
		int Natoms;
		int Nmodes;
		while(!intFile.eof())
		{
			std::string inputLine;
			std::getline(intFile, inputLine);
			// find out the number of atoms
			if (inputLine.find("NAtoms=") != std::string::npos)
			{
				std::istringstream iss(inputLine.substr(9));
				iss >> Natoms;
				Nmodes = 3 * Natoms - 6;
			}
			if (inputLine.find("Harmonic frequencies") != std::string::npos) // we have found the frequencies
			{
				int ModesRead = 0;
				while (ModesRead < Nmodes)
				{
					std::getline(intFile, inputLine);
					if (inputLine.find("Frequencies") != std::string::npos)
					{
						std::istringstream iss(inputLine);
						double freq;
						std::string dummy1, dummy2;
						iss >> dummy1 >> dummy2;
						while (iss >> freq)
							peakPositions.push_back(freq);
					}
					if (inputLine.find("IR Inten") != std::string::npos)
					{
						std::istringstream iss(inputLine);
						double inten;
						std::string dummy1, dummy2, dummy3;
						iss >> dummy1 >> dummy2 >> dummy3;
						while (iss >> inten)
						{
							peakIntensities.push_back(inten);
							ModesRead++;
						}
					}
				}
				break;
			}
		}
	}


	/*
	 * Print the stick spectrum to stdout
	 */
	std::cout << "# stick spectrum:\n";
	if (specType == fcclasses)
		std::cout << "# energy [eV]  ";
	else
		std::cout << "# energy [cm-1]";
	std::cout << "  intensity [arb. units]\n";
	for (int i = 0; i < static_cast<int>(peakPositions.size()); i++)
		std::cout << std::setprecision(8) << std::setw(18) << std::scientific
                  << peakPositions.at(i) << "    " << peakIntensities.at(i) << std::endl;

	/*
	 * Prepare the frequency axis.
	 * The min and max values are determined as follows:
	 * Take the distance between min and max peak positions and
	 * add half of it to the min and max values respectively:
	 *
	 * df = Pmax - Pmin
	 *
	 * fmin = Pmin - df/2
	 * fmax = Pmax + df/2
	 */
	std::vector<double> freqAxis(nPoints);
	double width = peakPositions.back() - peakPositions.front();
	double fmin = peakPositions.front() - (width / 2.0);
	double fmax = peakPositions.back() + (width / 2.0);
	double delta_f = (fmax - fmin) / (nPoints - 1);
	for (int i = 0; i < static_cast<int>(freqAxis.size()); i++)
		freqAxis.at(i) = fmin + (i * delta_f);

	/*
	 * Create the convoluted spectrum:
	 *
	 * The Lorenzian function is defined as:
	 *
	 *             1                  gamma / 2
	 * L(w) = A * ----  *  --------------------------------
	 *             pi        (w - w0)**2  +  (gamma / 2)**2
	 *
	 * where     A     is the intensity of the peak
	 *           gamma is the peak width (FWHM)
	 *           w0    is the center frequency of the peak
	 */
	std::vector<double> intAxis(nPoints, 0.0);
	for (int i = 0; i < static_cast<int>(peakPositions.size()); i++)
		for (int j = 0; j < nPoints; j++)
		{
			double value = peakIntensities.at(i) * ((0.5 * gamma)
					     / (M_PI * ((freqAxis.at(j) - peakPositions.at(i)) * (freqAxis.at(j) - peakPositions.at(i))
					    		 + ((0.5 * gamma) * (0.5 * gamma)))));
			intAxis.at(j) += value;
		}

	/*
	 * Normalize the spectrum:
	 * First, we calculate the overall integral:
	 */
	double integral = 0.0;
	for (int i = 0; i < nPoints; i++)
		integral += delta_f * intAxis.at(i);
	for (int i = 0; i < nPoints; i++)
		intAxis.at(i) /= integral;

	/*
	 * Print the convoluted spectrum to stdout:
	 */
	std::cout << "\n\n\n# convoluted spectrum:\n";
	for (int i = 0; i < nPoints; i++)
		std::cout << std::setprecision(8) << std::setw(18) << std::scientific
				  << freqAxis.at(i) << "    " << intAxis.at(i) << std::endl;

	return 0;
}
