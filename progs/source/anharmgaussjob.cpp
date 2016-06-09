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


#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <eigen3/Eigen/Eigenvalues>
#include <boost/program_options.hpp>
#include "GaussFchk.h"
#include "constants.h"
#include "utilities.h"
#include "vibrationalanalysis.h"

const double shiftFac = 0.03;

void gram_schmidt(Eigen::MatrixXd &d);

/*
 * Create a Gaussian09 input file for a frequency job with a structure
 * displaced by a specific normal mode
 */
int main(int argc, char *argv[])
{
	/*
	 * Process the command options
	 */
	std::string baseFileName;
	int mode;
	int posInt;
	namespace po = boost::program_options;
	po::options_description desc("Allowed options");
	desc.add_options()
		("help,h", "produce this help message")
		("file,f", po::value<std::string>(&baseFileName)->required(), "the Gaussian formatted checkpoint file to be processed")
		("mode,m", po::value<int>(&mode)->required(), "Index of normal mode")
		("positive,p", po::value<int>(&posInt)->required(), "Are we positive")
		;
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	if (vm.count("help"))
	{
		std::cout << desc << std::endl;
		return 1;
	}
	po::notify(vm);
	bool positive = posInt == 1 ? true : false;

	std::ifstream startFile(baseFileName + ".fchk", std::ifstream::in);
	GaussFchk startFchk(startFile);

	VibrationalAnalysis vibAn(startFchk);

	/*
	 * Create the displaced structure:
	 */
	const double signFac = positive ? 1.0 : -1.0;
	Eigen::VectorXd Xshift = vibAn.X() + vibAn.Lcart().col(mode) * signFac * shiftFac / ang2a0;

	/*
	 * Write the Gaussian input file with the new structure:
	 */
	std::ofstream GaussComFile(baseFileName + "_" + std::to_string(mode) + (positive ? "p" : "n") + ".com");
	GaussComFile << "%chk=checkpointfile" << std::endl;
	GaussComFile << "%mem=memory" << std::endl;
	GaussComFile << "#P method basisset" << std::endl;
	GaussComFile << "# int=ultrafine" << std::endl;
	GaussComFile << "# formcheck" << std::endl;
	GaussComFile << "# freq=(noraman,hpmodes,savenm)" << std::endl;
	GaussComFile << "# nosymm" << std::endl;
	GaussComFile << std::endl;
	GaussComFile << "displacement along mode " << mode << " with factor " << signFac << std::endl;
	GaussComFile << std::endl;
	GaussComFile << "0 1" << std::endl;
	for (int i = 0; i < vibAn.Natoms(); i++)
		GaussComFile << std::setw(3) << int(vibAn.atomicNumbers()(i))
		             << std::fixed << std::setprecision(9)
	                 << std::setw(14) << double(Xshift(3 * i + 0)) * ang2a0
			    	 << std::setw(14) << double(Xshift(3 * i + 1)) * ang2a0
				     << std::setw(14) << double(Xshift(3 * i + 2)) * ang2a0 << std::endl;
	GaussComFile << std::endl << std::endl;
}
