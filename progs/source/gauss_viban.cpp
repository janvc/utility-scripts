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
#include <iomanip>
#include <cmath>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Eigenvalues>
#include <boost/program_options.hpp>
#include "GaussFchk.h"
#include "utilities.h"
#include "constants.h"
#include "vibrationalanalysis.h"


int main(int argc, char *argv[])
{
	/*
	 * Some output-related settings
	 */
	int digits = 5;				// number of digits to write out
	bool clean = true;			// do we want to insert zeros?
	double threshold = 1.0e-15;	// threshold for inserting zeros

	/*
	 * Process the command options
	 */
	std::string filename;
	namespace po = boost::program_options;
	po::options_description desc("Allowed options");
	desc.add_options()
		("help,h", "produce this help message")
		("file,f", po::value<std::string>(&filename)->required(), "the Gaussian formatted checkpoint file to be processed")
	;
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	if (vm.count("help"))
	{
		std::cout << desc << std::endl;
		return 1;
	}
	po::notify(vm);


	std::ifstream fchkFile(filename, std::ifstream::in);
	VibrationalAnalysis vibAn(fchkFile);

	/*
	 * write the results to stdout
	 */
	std::cout << "Performing vibrational analysis of file '" << filename << "'" << std::endl << std::endl;
	std::cout << "Number of atoms:        " << std::setw(6) << vibAn.Natoms() << std::endl;
	std::cout << "Number of normal modes: " << std::setw(6) << vibAn.Nmodes() << std::endl << std::endl;

	std::cout << "atomic coordinates in bohr:\n";
	WriteVector(vibAn.X(), digits, clean, threshold);

	std::cout << std::endl;

	vibAn.prtMinGeo();

	std::cout << "\natomic masses in u:\n";
	WriteVector(vibAn.masses(), digits, clean, threshold);

	std::cout << "\ntotal mass in au: " << vibAn.totalMass() << std::endl << std::endl;

	std::cout << "non-mass-weighted cartesian hessian\n";
	WriteSymmGaussMatrix(vibAn.Fcart());

	std::cout << "\nmass-weighted cartesian hessian\n";
	WriteSymmGaussMatrix(vibAn.Fmwc());

	std::cout << "\neigenvalues of the mass-weighted cartesian hessian (au and cm-1):\n";
	for (int i = 0; i < vibAn.Ncoords(); i++)
		std::cout << std::scientific << std::setprecision(4)
				  << std::setw(5) << i + 1
				  << std::setw(14) << double(vibAn.mwcFC()(i))
				  << "  " << std::fixed << std::setw(10) << double(vibAn.mwcFrq()(i)) << std::endl;

	std::cout << "\ncenter of mass in bohr:\n";
	WriteVector(vibAn.com(), digits, clean, threshold);

	std::cout << "\ninertia tensor:\n";
	WriteSymmGaussMatrix(vibAn.inert());

	std::cout << "\nmoments of inertia:\n";
	WriteVector(vibAn.moments(), digits, clean, threshold);

	std::cout << "\nprincipal axes:\n";
	WriteGaussMatrix(vibAn.prinAxes());

	std::cout << "\nmetric of the principal axes:\n";
	WriteSymmGaussMatrix(vibAn.prinAxes().transpose() * vibAn.prinAxes());

	std::cout << "\nThe D matrix:\n";
	WriteGaussMatrix(vibAn.D());

	std::cout << "\nmetric of D matrix:\n";
	WriteSymmGaussMatrix(vibAn.D().transpose() * vibAn.D());

	std::cout << "\nThe Nvib x Nvib submatrix of the internal Hessian, according to eq. (6):\n";
	WriteSymmGaussMatrix(vibAn.Fint());

	std::cout << "\neigenvalues of the internal hessian (au and cm-1):\n";
	for (int i = 0; i < vibAn.Nmodes(); i++)
		std::cout << std::scientific << std::setprecision(4)
				  << std::setw(5) << i + 1
				  << std::setw(14) << double(vibAn.intFC()(i))
				  << "  " << std::fixed << std::setw(10) << double(vibAn.intFrq()(i)) << std::endl;

	std::cout << "\neigenvectors of the internal Hessian (internal displacements):\n";
	WriteGaussMatrix(vibAn.Lint());

	std::cout << "\nmetric of the eigenvectors:\n";
	WriteSymmGaussMatrix(vibAn.Lint().transpose() * vibAn.Lint());

	std::cout << "\nmass-weighted cartesian displacements:\n";
	WriteGaussMatrix(vibAn.Lmwc());

	std::cout << "\nmetric of mass-weighted cartesian displacements:\n";
	WriteSymmGaussMatrix(vibAn.Lmwc().transpose() * vibAn.Lmwc());

	std::cout << "\ncartesian displacements:\n";
	WriteGaussMatrix(vibAn.Lcart());

	std::cout << std::endl;

	vibAn.prtMinModes();

	std::cout << "\nmetric of the cartesian displacements:\n";
	WriteSymmGaussMatrix(vibAn.Lcart().transpose() * vibAn.Lcart());

	std::cout << "\nreduced masses of the normal modes:\n";
	WriteVector(vibAn.mu(), digits, clean, threshold);

	return 0;
}
