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

void gram_schmidt(Eigen::MatrixXd &d);

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
	std::cout << "Number of atoms: " << vibAn.Natoms() << std::endl;

	std::cout << "atomic coordinates in bohr:\n";
	WriteVector(vibAn.X(), digits, clean, threshold);

	vibAn.prtMinGeo();

	std::cout << "atomic masses in u:\n";
	WriteVector(vibAn.masses(), digits, clean, threshold);

	std::cout << "total mass in au: " << vibAn.totalMass() << std::endl;

	std::cout << "non-mass-weighted cartesian hessian\n";
	WriteMatrix(vibAn.Fcart(), digits, clean, threshold);

	std::cout << "mass-weighted cartesian hessian\n";
	WriteMatrix(vibAn.Fmwc(), digits, clean, threshold);

	std::cout << "eigenvalues of the mass-weighted cartesian hessian (au and cm-1):\n";
	for (int i = 0; i < vibAn.Ncoords(); i++)
		std::cout << std::scientific << std::setprecision(4) << std::setw(12) << double(vibAn.mwcFC()(i))
				  << "  " << std::fixed << std::setw(9) << double(vibAn.mwcFrq()(i)) << std::endl;

	std::cout << "center of mass in bohr:\n";
	WriteVector(vibAn.com(), digits, clean, threshold);

	std::cout << "inertia tensor:\n";
	WriteMatrix(vibAn.inert(), digits, clean, threshold);

	std::cout << "moments of inertia:\n";
	WriteVector(vibAn.moments(), digits, clean, threshold);

	std::cout << "principal axes:\n";
	WriteMatrix(vibAn.prinAxes(), digits, clean, threshold);

	std::cout << "metric of the principal axes:\n";
	WriteMatrix(vibAn.prinAxes().transpose() * vibAn.prinAxes(), digits, clean, threshold);

	std::cout << "The D matrix:\n";
	WriteMatrix(vibAn.D(), digits, clean, threshold);

	std::cout << "metric of D matrix:\n";
	WriteMatrix(vibAn.D().transpose() * vibAn.D(), digits, clean, threshold);

	std::cout << "The Nvib x Nvib submatrix of the internal Hessian, according to eq. (6):\n";
	WriteMatrix(vibAn.Fint(), digits, clean, threshold);

	std::cout << "eigenvalues of the internal hessian (au and cm-1):\n";
	for (int i = 0; i < vibAn.Nmodes(); i++)
		std::cout << std::scientific << std::setprecision(4) << std::setw(12) << double(vibAn.intFC()(i))
				  << "  " << std::fixed << std::setw(9) << double(vibAn.intFrq()(i)) << std::endl;

	std::cout << "eigenvectors of the internal Hessian (internal displacements):\n";
	WriteMatrix(vibAn.Lint(), digits, clean, threshold);

	std::cout << "metric of the eigenvectors:\n";
	WriteMatrix(vibAn.Lint().transpose() * vibAn.Lint(), digits, clean, threshold);

	std::cout << "mass-weighted cartesian displacements:\n";
	WriteMatrix(vibAn.Lmwc(), digits, clean, threshold);

	std::cout << "metric of mass-weighted cartesian displacements:\n";
	WriteMatrix(vibAn.Lmwc().transpose() * vibAn.Lmwc(), digits, clean, threshold);

	std::cout << "cartesian displacements:\n";
	WriteMatrix(vibAn.Lcart(), digits, clean, threshold);

	vibAn.prtMinModes();

	std::cout << "metric of cartesian displacements:\n";
	WriteMatrix(vibAn.Lcart().transpose() * vibAn.Lcart(), digits, clean, threshold);

	std::cout << "reduced masses of the normal modes:\n";
	WriteVector(vibAn.mu(), digits, clean, threshold);

	Eigen::VectorXd xshift = vibAn.X() + vibAn.Lcart().col(56) * 0.03 / ang2a0;

	std::cout << "\nThe shifted structure:\n";
	std::cout << "------------------------------------------\n";
	for (int i = 0; i < vibAn.Natoms(); i++)
		std::cout << std::fixed << std::setprecision(6)
	<< std::setw(6) << i + 1
	<< std::setw(12) << double(xshift(3*i+0)) * ang2a0
	<< std::setw(12) << double(xshift(3*i+1)) * ang2a0
	<< std::setw(12) << double(xshift(3*i+2)) * ang2a0 << std::endl;
	std::cout << "------------------------------------------\n";

	std::cout << "\nDifference between the original and the shifted structure:\n";
	std::cout << "------------------------------------------\n";
	for (int i = 0; i < vibAn.Natoms(); i++)
		std::cout << std::fixed << std::setprecision(6)
	<< std::setw(6) << i + 1
	<< std::setw(12) << (double(xshift(3*i+0)) - double(vibAn.X()(3*i+0))) * ang2a0
	<< std::setw(12) << (double(xshift(3*i+1)) - double(vibAn.X()(3*i+1))) * ang2a0
	<< std::setw(12) << (double(xshift(3*i+2)) - double(vibAn.X()(3*i+2))) * ang2a0 << std::endl;
	std::cout << "------------------------------------------\n";

	return 0;
}
