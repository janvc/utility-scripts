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
#include <cmath>
#include <boost/program_options.hpp>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Eigenvalues>
#include "GaussFchk.h"
#include "constants.h"
#include "utilities.h"


int main(int argc, char *argv[])
{
	int digits = 5;
	bool clean = true;
	double threshold = 1.0e-15;

	/*
	 * Process command line options:
	 */
	std::string FchkFileName;
	namespace po = boost::program_options;
	po::options_description desc("Allowed options");
	desc.add_options()
					("help,h", "produce this help message")
					("file,f", po::value<std::string>(&FchkFileName)->required(), "the FChk file")
					;
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	if (vm.count("help"))
	{
		std::cout << desc << std::endl;
		return 1;
	}
	po::notify(vm);

	std::ifstream FchkFile(FchkFileName, std::ifstream::in);
	GaussFchk Fchk(FchkFile);

	int Natoms = Fchk.ReadInteger("Number of atoms");
	int Ncoords = 3 * Natoms;
	int Nmodes = Ncoords - 6;
	int Ndata = Ncoords * (Ncoords + 1) / 2;

	Eigen::VectorXd Masses = Fchk.ReadVector("Real atomic weights");
	Eigen::MatrixXd origHess = Fchk.ReadSymmetricMatrix("Cartesian Force Constants");
	Eigen::VectorXd derivData = Fchk.ReadVector("Cartesian 3rd/4th derivatives");

	Eigen::MatrixXd MassInvMat = Eigen::MatrixXd::Zero(Ncoords,Ncoords);
	for (int i = 0; i < Natoms; i++)
	{
		MassInvMat(3*i+0, 3*i+0) = Masses(i);
		MassInvMat(3*i+1, 3*i+1) = Masses(i);
		MassInvMat(3*i+2, 3*i+2) = Masses(i);
	}

	std::vector<Eigen::MatrixXd> Hessians;

	for (int i = 0; i < Nmodes; i++)
	{
		Eigen::MatrixXd tmpHess(Ncoords,Ncoords);
		Eigen::VectorXd tmpVec = derivData.segment(i * Ndata, Ndata);

		for (int j = 0; j < Ncoords; j++)
			for (int k = 0; k <= j; k++)
			{
				tmpHess(j,k) = tmpVec((j*(j+1)/2)+k);
				tmpHess(k,j) = tmpHess(j,k);
			}

		Hessians.push_back(tmpHess);
	}

	Eigen::MatrixXd firstHess = Hessians.at(0);
	Eigen::MatrixXd mwoHess = MassInvMat.transpose() * origHess * MassInvMat;
	Eigen::MatrixXd mwHess = MassInvMat.transpose() * firstHess * MassInvMat;
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> origHessDiag(mwoHess);
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> firstHessDiag(firstHess);


	/*
	 * write the results to stdout
	 */
	std::cout << "Number of atoms: " << Natoms << std::endl;
	std::cout << "Number of numbers in the large vector: " << derivData.size() << std::endl;


	std::cout << "atomic masses in u:\n";
    Utils::WriteVector(Masses, digits, clean, threshold);

	std::cout << "first 100 elements of the big data vector:\n";
    Utils::WriteVector(derivData.segment(0, 100), digits, clean, threshold);
	std::cout << "original cartesian Hessian at the eq. geo.:\n";
    Utils::WriteMatrix(origHess, digits, clean, threshold);
	std::cout << "first cartesian hessian\n";
    Utils::WriteMatrix(firstHess, digits, clean, threshold);

	std::cout << "mass-weighted first cartesian hessian\n";
    Utils::WriteMatrix(mwHess, digits, clean, threshold);

	std::cout << "eigenvalues of the mass-weighted original cartesian hessian (au and cm-1):\n";
	for (int i = 0; i < Ncoords; i++)
	{
		double value = origHessDiag.eigenvalues()(i);
		double wavenumber;
		if (value > 0)
			wavenumber = sqrt(value * Eh / (a0 * a0 * amu)) / (2.0 * M_PI * c0 * 100.0);
		else
			wavenumber = -sqrt(-value * Eh / (a0 * a0 * amu)) / (2.0 * M_PI * c0 * 100.0);
		std::cout << std::scientific << std::setprecision(4) << std::setw(12) << value
				  << "  " << std::fixed << std::setw(9) << wavenumber << std::endl;
	}

	std::cout << "eigenvalues of the mass-weighted first cartesian hessian (au and cm-1):\n";
	for (int i = 0; i < Ncoords; i++)
	{
		double value = firstHessDiag.eigenvalues()(i);
		double wavenumber;
		if (value > 0)
			wavenumber = sqrt(value * Eh / (a0 * a0 * amu)) / (2.0 * M_PI * c0 * 100.0);
		else
			wavenumber = -sqrt(-value * Eh / (a0 * a0 * amu)) / (2.0 * M_PI * c0 * 100.0);
		std::cout << std::scientific << std::setprecision(4) << std::setw(12) << value
				  << "  " << std::fixed << std::setw(9) << wavenumber << std::endl;
	}
}
