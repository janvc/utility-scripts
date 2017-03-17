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
#include <boost/program_options.hpp>
#include "GaussFchk.h"
#include "constants.h"
#include "vibrationalanalysis.h"


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
    int mode1;
    int mode2;
    double shiftFac1;
    double shiftFac2;
	namespace po = boost::program_options;
	po::options_description desc("Allowed options");
	desc.add_options()
		("help,h", "produce this help message")
		("file,f", po::value<std::string>(&baseFileName)->required(), "the Gaussian formatted checkpoint file to be processed")
        ("mode,m", po::value<int>(&mode1)->required(), "Index of 1st normal mode")
        ("node,n", po::value<int>(&mode2)->required(), "Index of 2nd normal mode")
        ("scale,s", po::value<double>(&shiftFac1)->required(), "Scale factor for the 1st displacement")
        ("tcale,t", po::value<double>(&shiftFac2)->required(), "Scale factor for the 2nd displacement")
		;
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	if (vm.count("help"))
	{
		std::cout << desc << std::endl;
		return 1;
	}
	po::notify(vm);

	std::ifstream startFile(baseFileName, std::ifstream::in);
	GaussFchk startFchk(startFile);

	VibrationalAnalysis vibAn(startFchk);

	/*
	 * Create the displaced structure:
	 */
    Eigen::VectorXd modeVec = Eigen::VectorXd::Zero(vibAn.Nmodes());
    modeVec(mode1 - 1) = shiftFac1;
    modeVec(mode2 - 1) = shiftFac2;
    Eigen::VectorXd Xshift = vibAn.X() + vibAn.Lcart() * modeVec / ang2a0;

	/*
	 * Write the Gaussian input file with the new structure:
	 */
	std::cout << "%chk=" << std::endl;
	std::cout << "%mem=" << std::endl;
	std::cout << "#P method basisset" << std::endl;
	std::cout << std::endl;
    std::cout << "displacement along mode " << mode1 << " and " << mode2 << " by " << shiftFac1 << " and " << shiftFac2 << std::endl;
	std::cout << std::endl;
	std::cout << "0 1" << std::endl;
	for (int i = 0; i < vibAn.Natoms(); i++)
		std::cout << std::setw(3) << int(vibAn.atomicNumbers()(i))
		          << std::fixed << std::setprecision(9)
	              << std::setw(14) << double(Xshift(3 * i + 0)) * ang2a0
				  << std::setw(14) << double(Xshift(3 * i + 1)) * ang2a0
				  << std::setw(14) << double(Xshift(3 * i + 2)) * ang2a0 << std::endl;
	std::cout << std::endl << std::endl;
}
