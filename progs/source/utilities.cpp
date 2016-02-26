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


#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Eigenvalues>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <boost/filesystem.hpp>
#include "utilities.h"
#include "constants.h"


void WriteMatrix(const Eigen::MatrixXd &mat, const int dig, const bool clean, const double thres)
{
	int width = dig + 8;
	// write header:
	std::cout << "  ";
	for (int i = 0; i < mat.cols(); i++)
		std::cout << std::setw(width) << i+1;
	std::cout << std::endl << "       ";
	for (int i = 0; i < mat.cols(); i++)
		std::cout << std::string(width,'-');
	std::cout << std::endl;

	for (int i = 0; i < mat.rows(); i++)
	{
		// write row number
		std::cout << std::setw(5) << i+1 << " |";
		for (int j = 0; j < mat.cols(); j++)
		{
			if (clean)
			{
				if (std::abs(double(mat(i,j))) > thres)
					std::cout << std::scientific << std::setprecision(dig) << std::setw(width) << mat(i,j);
				else
					std::cout << std::setw(width) << "0";
			}
			else
				std::cout << std::scientific << std::setprecision(dig) << std::setw(width) << mat(i,j);
		}
		std::cout << std::endl;

	}
	std::cout << "       ";
	for (int i = 0; i < mat.cols(); i++)
		std::cout << std::string(width,'-');
	std::cout << std::endl;
}


void WriteMatrixToFile(std::ofstream &stream, const Eigen::MatrixXd &mat, const int dig, const bool clean, const double thres)
{
	int width = dig + 8;
	// write header:
	stream << "  ";
	for (int i = 0; i < mat.cols(); i++)
		stream << std::setw(width) << i+1;
	stream << std::endl << "      /";
	for (int i = 0; i < mat.cols(); i++)
		stream << std::string(width,'-');
	stream << std::endl;

	for (int i = 0; i < mat.rows(); i++)
	{
		// write row number
		stream << std::setw(5) << i+1 << " |";
		for (int j = 0; j < mat.cols(); j++)
		{
			if (clean)
			{
				if (std::abs(double(mat(i,j))) > thres)
					stream << std::scientific << std::setprecision(dig) << std::setw(width) << mat(i,j);
				else
					stream << std::setw(width) << "0";
			}
			else
				stream << std::scientific << std::setprecision(dig) << std::setw(width) << mat(i,j);
		}
		stream << std::endl;
	}
	stream << "       ";
	for (int i = 0; i < mat.cols(); i++)
		stream << std::string(width,'-');
	stream << std::endl;
}

void WriteVector(const Eigen::VectorXd &vec, const int dig, const bool clean, const double thres)
{
	int width = dig + 8;

	for (int i = 0; i < vec.size(); i++)
	{
		std::cout << std::setw(5) << i+1 << " |";
		if (clean)
		{
			if (std::abs(double(vec(i))) > thres)
				std::cout << std::scientific << std::setprecision(dig) << std::setw(width) << double(vec(i)) << std::endl;
			else
				std::cout << std::setw(width) << "0" << std::endl;
		}
		else
			std::cout << std::scientific << std::setprecision(dig) << std::setw(width) << double(vec(i)) << std::endl;
	}
}


void WriteVectorToFile(std::ofstream &stream, const Eigen::VectorXd &vec, const int dig, const bool clean, const double thres)
{
	int width = dig + 8;

	for (int i = 0; i < vec.size(); i++)
	{
		stream << std::setw(5) << i+1 << " |";
		if (clean)
		{
			if (std::abs(double(vec(i))) > thres)
				stream << std::scientific << std::setprecision(dig) << std::setw(width) << double(vec(i)) << std::endl;
			else
				stream << std::setw(width) << "0" << std::endl;
		}
		else
			stream << std::scientific << std::setprecision(dig) << std::setw(width) << double(vec(i)) << std::endl;
	}
}


void WriteFortranNumber(std::ofstream &stream, const double number)
{
	std::ostringstream strstr;
	strstr << std::scientific << std::setprecision(8) << std::setw(15) << std::setfill(' ') << number;
	std::string str = strstr.str();
	std::replace(str.begin(), str.end(), 'e', 'd');
	stream << str;
}

Eigen::Vector3d calc_com(Eigen::VectorXd x, Eigen::VectorXd m)
{
	Eigen::Vector3d com = Eigen::Vector3d::Zero();

	for (int i = 0; i < m.size(); i++)
	{
		com(0) += m(i) * x(3 * i + 0);
		com(1) += m(i) * x(3 * i + 1);
		com(2) += m(i) * x(3 * i + 2);
	}
	com /= m.sum();

	return com;
}

Eigen::Matrix3d calc_inert(Eigen::VectorXd x, Eigen::VectorXd m)
{
	Eigen::Matrix3d inert = Eigen::Matrix3d::Zero();

	Eigen::Vector3d com = calc_com(x, m);

	for (int i = 0; i < m.size(); i++)
	{
		inert(0,0) += m(i) * ((x(3*i+1)-com(1))*(x(3*i+1)-com(1)) + (x(3*i+2)-com(2))*(x(3*i+2)-com(2)));
		inert(1,1) += m(i) * ((x(3*i+0)-com(0))*(x(3*i+0)-com(0)) + (x(3*i+2)-com(2))*(x(3*i+2)-com(2)));
		inert(2,2) += m(i) * ((x(3*i+0)-com(0))*(x(3*i+0)-com(0)) + (x(3*i+1)-com(1))*(x(3*i+1)-com(1)));
		inert(0,1) -= m(i) *  (x(3*i+0)-com(0)) * (x(3*i+1)-com(1));
		inert(0,2) -= m(i) *  (x(3*i+0)-com(0)) * (x(3*i+2)-com(2));
		inert(1,2) -= m(i) *  (x(3*i+1)-com(1)) * (x(3*i+2)-com(2));
	}

	inert(1,0) = inert(0,1);
	inert(2,0) = inert(0,2);
	inert(2,1) = inert(1,2);

	return inert;
}

Eigen::MatrixXd vibrationalAnalysis(const Eigen::VectorXd &positions, const Eigen::VectorXd &masses, const Eigen::MatrixXd &hessian)
{
	if (positions.size() != hessian.cols() || positions.size() != 3 * masses.size())
		return Eigen::MatrixXd();

	int Natoms = masses.size();
	int Ncoords = 3 * Natoms;
	int Nmodes = Ncoords - 6;

	//double TotMass = masses.sum();
	Eigen::VectorXd MassVector(Ncoords);
	for (int i = 0; i <  Natoms; i++)
	{
		MassVector(3 * i + 0) = masses(i);
		MassVector(3 * i + 1) = masses(i);
		MassVector(3 * i + 2) = masses(i);
	}
	Eigen::MatrixXd MassInvMat = Eigen::MatrixXd::Zero(Ncoords, Ncoords);
	for (int i = 0; i < Ncoords; i++)
		MassInvMat(i,i) = 1.0 / sqrt(double(MassVector(i)));

	Eigen::MatrixXd f_mwc = MassInvMat.transpose() * hessian * MassInvMat;

	Eigen::Vector3d com = calc_com(positions, masses);
	Eigen::Matrix3d inert = calc_inert(positions, masses);
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> InertDiag(inert);
	Eigen::Matrix3d rotmat = InertDiag.eigenvectors().transpose();

	Eigen::MatrixXd D(Eigen::MatrixXd::Random(Ncoords, Ncoords));
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < Natoms; j++)
		{
			D(3*j+0, i) = 0;
			D(3*j+1, i) = 0;
			D(3*j+2, i) = 0;
			D(3*j+i, i) = sqrt(double(masses(j)));
		}
	for (int i = 0; i < Natoms; i++)
	{
		Eigen::Vector3d x_int;
		x_int(0) = positions(3*i+0) - com(0);
		x_int(1) = positions(3*i+1) - com(1);
		x_int(2) = positions(3*i+2) - com(2);
		x_int = rotmat * x_int;

		D(3*i+0, 3) = (x_int(1) * rotmat(2,0) - x_int(2) * rotmat(1,0)) * sqrt(double(masses(i)));
		D(3*i+1, 3) = (x_int(1) * rotmat(2,1) - x_int(2) * rotmat(1,1)) * sqrt(double(masses(i)));
		D(3*i+2, 3) = (x_int(1) * rotmat(2,2) - x_int(2) * rotmat(1,2)) * sqrt(double(masses(i)));
		D(3*i+0, 4) = (x_int(2) * rotmat(0,0) - x_int(0) * rotmat(2,0)) * sqrt(double(masses(i)));
		D(3*i+1, 4) = (x_int(2) * rotmat(0,1) - x_int(0) * rotmat(2,1)) * sqrt(double(masses(i)));
		D(3*i+2, 4) = (x_int(2) * rotmat(0,2) - x_int(0) * rotmat(2,2)) * sqrt(double(masses(i)));
		D(3*i+0, 5) = (x_int(0) * rotmat(1,0) - x_int(1) * rotmat(0,0)) * sqrt(double(masses(i)));
		D(3*i+1, 5) = (x_int(0) * rotmat(1,1) - x_int(1) * rotmat(0,1)) * sqrt(double(masses(i)));
		D(3*i+2, 5) = (x_int(0) * rotmat(1,2) - x_int(1) * rotmat(0,2)) * sqrt(double(masses(i)));
	}
	for (int i = 0; i < Ncoords; i++)
		D.col(i) /= Eigen::VectorXd(D.col(i)).norm();
	Eigen::MatrixXd newMat = Eigen::MatrixXd::Zero(Ncoords, Ncoords);
	newMat.block(0, 0, Ncoords, 6) = D.block(0, 0, Ncoords, 6);
	for (int i = 6; i < Ncoords; i++)
	{
		Eigen::VectorXd tempVec = D.col(i);
		for (int j = 0; j < i; j++)
		{
			double factor = tempVec.dot(newMat.col(j));
			tempVec -= factor * newMat.col(j);
		}
		tempVec.normalize();
		newMat.col(i) = tempVec;
	}
	D = newMat;

	Eigen::MatrixXd f_int = D.transpose().block(6, 0, Nmodes, Ncoords) * f_mwc * D.block(0, 6, Ncoords, Nmodes);
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> IntHessDiag(f_int);
	Eigen::VectorXd forceConstants = IntHessDiag.eigenvalues();
	Eigen::MatrixXd normalModes = D.block(0, 6, Ncoords, Nmodes) * IntHessDiag.eigenvectors();

	Eigen::MatrixXd result(Ncoords + 1, Nmodes);
	for (int i = 0; i < Nmodes; i++)
	{
		result(0, i) = forceConstants(i);
		result.block(1, i, Ncoords, 1) = normalModes.col(i);
	}

	// return a matrix containing the force constants in atomic units in the first row
	// and the mass-weighted cartesian normal modes in the remaining rows.
	return result;
}

void createMCTDHfiles(const Eigen::MatrixXd &J, const Eigen::VectorXd &K, const Eigen::VectorXd &f1, const Eigen::VectorXd &f2, const double dE, std::ofstream &logFile)
{
	int Nmodes = f2.size();
	/*
	 * ############################################################################################
	 *                                  calculate the new PES parameters
	 * ############################################################################################
	 *
	 *
	 * Calculate the new force constants.
	 */
	Eigen::VectorXd fp(Nmodes);
	for (int m = 0; m < Nmodes; m++)
	{
		fp(m) = 0.0;
		for (int n = 0; n < Nmodes; n++)
			fp(m) += f2(n) * J(n,m) * J(n,m);
	}


	/*
	 * Calculate the first-order coefficients.
	 */
	Eigen::VectorXd kappa(Nmodes);
	for (int m = 0; m < Nmodes; m++)
	{
		kappa(m) = 0.0;
		for (int n = 0; n < Nmodes; n++)
			kappa(m) += f2(n) * K(n) * J(n,m);
	}


	/*
	 * Calculate the couplings.
	 */
	Eigen::MatrixXd phi(Eigen::MatrixXd::Zero(Nmodes,Nmodes));
	Eigen::MatrixXd phiFull(Eigen::MatrixXd::Zero(Nmodes,Nmodes));
	for (int m = 0; m < Nmodes; m++)
	{
		for (int o = m + 1; o < Nmodes; o++)
		{
			phi(m,o) = 0.0;
			phiFull(m,o) = 0.0;
			for (int n = 0; n < Nmodes; n++)
				phi(m,o) += f2(n) * J(n,m) * J(n,o);
			phi(o,m) = phi(m,o);
			phiFull(m,o) = phi(m,o);
			phiFull(o,m) = phi(m,o);
		}
		phiFull(m,m) = fp(m);
	}


	/*
	 * Calculate the energy shifts.
	 */
	Eigen::VectorXd d(Nmodes);
	for (int i = 0; i < Nmodes; i++)
		d(i) = (0.5 * f2(i) * K(i) * K(i));

	/*
	 * write the coupling matrix to the log file:
	 */
	logFile << "The coupling matrix:\n";
	WriteMatrixToFile(logFile, phi, 3, true);
	logFile << "Coupling matrix phi with the force constants fp on the diagonal:\n";
	WriteMatrixToFile(logFile, phiFull, 3, true);


	/*
	 * ############################################################################################
	 *                                  write the MCTDH files
	 * ############################################################################################
	 *
	 *
	 * Now we can finally write the MCTDH input and operator files :)
	 * First, inquire the desired base-name for the files.
	 */
	std::cout << "Enter the base-name for the MCTDH files to be generated.\nThe files <name>.inp and <name>.op will then be written.\n>>> ";
	std::string basename;
	std::cin >> basename;

	std::string inputFileName = basename + ".inp";
	std::string operatorFileName = basename + ".op";

	if (boost::filesystem::exists(inputFileName) || boost::filesystem::exists(operatorFileName))
	{
		std::cout << "One of the MCTDH files already exists. Should they be overwritten? (Y/N)\n>>> ";
		char answer;
		std::cin >> answer;
		if (answer == 'N' || answer == 'n')
			return;
	}
	std::ofstream inputFile(inputFileName);
	std::ofstream operatorFile(operatorFileName);

	inputFile.precision(1);
	operatorFile.precision(8);

	/*
	 * The run-section
	 */
	inputFile << "run-section\n";
	inputFile << "    name =\n";
	inputFile << "    propagation\n";
	inputFile << "    tfinal =\n";
	inputFile << "    tout =\n";
	inputFile << "    tpsi =\n";
	inputFile << "    psi gridpop auto steps graphviz\n";
	inputFile << "end-run-section\n\n";

	/*
	 * The operator-section
	 */
	inputFile << "operator-section\n";
	inputFile << "    opname = " << basename << std::endl;
	inputFile << "end-operator-section\n\n";

	/*
	 * The mlbasis-section
	 */
	// rearrange the modes in order of decreasing coupling
	Eigen::MatrixXd phi_sort = phi.cwiseAbs();	// use the absolute value of the coupling
	for (int i = 0; i < Nmodes; i++)
		for (int j = i + 1; j < Nmodes; j++)
		{
			phi_sort(i,j) /= double(f2(j) / f2(i));	// divide by the frequency ratio
			phi_sort(j,i) = phi_sort(i,j);
		}
	std::vector<int> sortedModes;
	while (phi_sort.norm() > 0.0)
	{
		Eigen::MatrixXd::Index maxRow, maxCol;
		phi_sort.maxCoeff(&maxRow, &maxCol);
		phi_sort(maxRow, maxCol) = 0.0;

		if (std::find(sortedModes.begin(), sortedModes.end(), maxRow) == sortedModes.end())
			sortedModes.push_back(maxRow);

		if (std::find(sortedModes.begin(), sortedModes.end(), maxCol) == sortedModes.end())
					sortedModes.push_back(maxCol);
	}
	// determine the required number of layers
	int layers = 1;
	while (pow(2.0, layers) < Nmodes)
		layers++;
	layers--;
	// determine the number of nodes in each layer
	std::vector<int> nodesPerLayer(layers);
	nodesPerLayer.at(layers - 1) = Nmodes / 2;
	for (int i = layers - 1; i > 0; i--)
	{
		nodesPerLayer.at(i - 1) = nodesPerLayer.at(i) / 2;
	}
	inputFile << "mlbasis-section\n";
	for (int i = 0; i < Nmodes - 1; i += 2)
	{
		if (sortedModes.size() - i == 3)
			inputFile << "    [q_" << std::setfill('0') << std::setw(3) << sortedModes.at(i+0) + 1
				          << " q_" << std::setfill('0') << std::setw(3) << sortedModes.at(i+1) + 1
						  << " q_" << std::setfill('0') << std::setw(3) << sortedModes.at(i+2) + 1 << "]\n";
		else
			inputFile << "    [q_" << std::setfill('0') << std::setw(3) << sortedModes.at(i+0) + 1
						  << " q_" << std::setfill('0') << std::setw(3) << sortedModes.at(i+1) + 1 << "]\n";
	}
	inputFile << "end-mlbasis-section\n\n";

	/*
	 * The pbasis-section
	 */
	inputFile << "pbasis-section\n";
	for (int i = 0; i < Nmodes; i++)
		inputFile << "    q_" << std::setfill('0') << std::setw(3) << i + 1
				  << "  ho  " << std::setw(3) << std::setfill(' ') << lrint(-0.6 * log(double(fp(i)))) + 5 << "  xi-xf  "
				  //
				  // the basis boundarie are -kappa / fp +- 6.5 / fp**1/4
				  //
				  << std::fixed << std::setfill(' ') << std::setw(8)
				  << -double(kappa(i) / fp(i)) - (std::abs(double(kappa(i) / fp(i))) + 6.5 / (sqrt(2.0) * pow(double(fp(i)), 0.25)))
				  << std::fixed << std::setfill(' ') << std::setw(8)
				  << -double(kappa(i) / fp(i)) + (std::abs(double(kappa(i) / fp(i))) + 6.5 / (sqrt(2.0) * pow(double(fp(i)), 0.25))) << std::endl;
	inputFile << "end-pbasis-section\n\n";

	/*
	 * The integrator section
	 */
	inputFile << "integrator-section\n";
	inputFile << "    vmf\n";
	inputFile << "    abm = 6, 1.0d-7, 0.01d0\n";
	inputFile << "end-integrator-section\n\n";

	/*
	 * The init wf section
	 */
	inputFile << "init_wf-section\n";
	inputFile << "    build\n";
	for (int i = 0; i < Nmodes; i++)
		inputFile << "        q_" << std::setfill('0') << std::setw(3) << i + 1
				  << "  eigenf"
				  << "  Eq_" << std::setfill('0') << std::setw(3) << i + 1
				  << "  pop = 1\n";
	inputFile << "    end-build\n";
	inputFile << "end-init_wf-section\n\n";
	inputFile << "end-input\n\n";


	/*
	 * Now the operator file
	 *
	 * First the op-define section
	 */
	operatorFile << "op_define-section\n";
	operatorFile << "    title\n";
	operatorFile << "        " << basename << std::endl;
	operatorFile << "    end-title\n";
	operatorFile << "end-op_define-section\n\n";

	/*
	 * The parameter section
	 */
	operatorFile << "parameter-section\n";
	// the masses
	for (int i = 0; i < Nmodes; i++)
		operatorFile << "    mass_q_" << std::setfill('0') << std::setw(3) << i + 1
					 << "  =  1.0\n";
	// the ground state force constants
	for (int i = 0; i < Nmodes; i++)
	{
		operatorFile << "    f1_" << std::setfill('0') << std::setw(3) << i + 1 << "      = ";
		WriteFortranNumber(operatorFile, double(f1(i)));
		operatorFile << std::endl;
	}
	// the excited state force constants
	for (int i = 0; i < Nmodes; i++)
	{
		operatorFile << "    f2_" << std::setfill('0') << std::setw(3) << i + 1 << "      = ";
		WriteFortranNumber(operatorFile, double(f2(i)));
		operatorFile << std::endl;
	}
	// the new effective excited state force constants
	for (int i = 0; i < Nmodes; i++)
	{
		operatorFile << "    fp_" << std::setfill('0') << std::setw(3) << i + 1 << "      = ";
		WriteFortranNumber(operatorFile, double(fp(i)));
		operatorFile << std::endl;
	}
	// the couplings
	for (int i = 0; i < Nmodes; i++)
		for (int j = i + 1; j < Nmodes; j++)
		{
			operatorFile << "    phi_" << std::setfill('0') << std::setw(3) << i + 1
						 << "_" << std::setfill('0') << std::setw(3) << j + 1 << " = ";
			WriteFortranNumber(operatorFile, double(phi(i,j)));
			operatorFile << std::endl;
		}
	// the first-order coefficients (shifts)
	for (int i = 0; i < Nmodes; i++)
	{
		operatorFile << "    kappa_" << std::setfill('0') << std::setw(3) << i + 1 << "   = ";
		WriteFortranNumber(operatorFile, double(kappa(i)));
		operatorFile << std::endl;
	}
	// the energy offsets
	for (int i = 0; i < Nmodes; i++)
	{
		operatorFile << "    d_" << std::setfill('0') << std::setw(3) << i + 1 << "       = ";
		WriteFortranNumber(operatorFile, double(d(i)));
		operatorFile << std::endl;
	}
	// the electronic offset minus the ground state ZPE
	double zpe1 = 0.0;
	for (int i = 0; i < Nmodes; i++)
		zpe1 += 0.5 * sqrt(double(f1(i)));
	operatorFile << "    dE          = ";
	WriteFortranNumber(operatorFile, dE - zpe1);
	operatorFile << "\nend-parameter-section\n\n";

	/*
	 * The hamiltonian section
	 */
	operatorFile << "hamiltonian-section";
	for (int i = 0; i < Nmodes; i++)
	{
		if (i % 8 == 0)
			operatorFile << std::endl << "modes";
		operatorFile << " | q_" << std::setfill('0') << std::setw(3) << i + 1;
	}
	operatorFile << std::endl;
	for (int i = 0; i < Nmodes; i++)
		operatorFile << "1.0         |" << i + 1 << " KE\n";
	for (int i = 0; i < Nmodes; i++)
		operatorFile << "0.5*fp_" << std::setfill('0') << std::setw(3) << i + 1
					 << "  |" << i + 1 << " q^2\n";
	for (int i = 0; i < Nmodes; i++)
		for (int j = i + 1; j < Nmodes; j++)
			operatorFile << "phi_" << std::setfill('0') << std::setw(3) << i + 1
						 << "_" << std::setfill('0') << std::setw(3) << j + 1
						 << " |" << i + 1 << " q"
						 << " |" << j + 1 << " q\n";
	for (int i = 0; i < Nmodes; i++)
		operatorFile << "kappa_" << std::setfill('0') << std::setw(3) << i + 1
					 << "   |" << i + 1 << " q\n";
	for (int i = 0; i < Nmodes; i++)
			operatorFile << "d_" << std::setfill('0') << std::setw(3) << i + 1
						 << "       |" << i + 1 << " 1\n";
	operatorFile << "dE          |1 1\n";
	operatorFile << "end-hamiltonian-section\n\n";

	/*
	 * One-dimensional Hamiltonians for the ground state normal modes
	 */
	for (int i = 0; i < Nmodes; i++)
	{
		operatorFile << "hamiltonian-section_Eq_" << std::setfill('0') << std::setw(3) << i + 1 << std::endl;
		operatorFile << "usediag\n";
		operatorFile << "modes      | q_" << std::setfill('0') << std::setw(3) << i + 1 << std::endl;
		operatorFile << "1.0        |1 KE\n";
		operatorFile << "0.5*f1_" << std::setfill('0') << std::setw(3) << i + 1 << " |1 q^2\n";
		operatorFile << "end-hamiltonian-section\n\n";
	}
	operatorFile << "end-operator\n";


	/*
	 * Diagonalize the coupling matrix to get the force constants back
	 */
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> phiSolv(phiFull);
	logFile << "Eigenvalues of the full force constant / coupling matrix:\n";
	WriteVectorToFile(logFile, phiSolv.eigenvalues());

	/*
	 * Solve the linear system of the full coupling matrix and the kappa vector
	 * to get the coordinates of the minimum
	 */
	Eigen::ColPivHouseholderQR<Eigen::MatrixXd> phiLin(phiFull);
	Eigen::VectorXd minima = phiLin.solve(-kappa);
	logFile << "minimum coordinates\n";
	WriteVectorToFile(logFile, minima);


	/*
	 * calculate the potential energy at the minimum
	 */
	double Emin = 0.0;

	// first, quadratic term:
	for (int i = 0; i < Nmodes; i++)
		Emin += 0.5 * fp(i) * minima(i) * minima(i);

	// second, coupling term:
	for (int i = 0; i < Nmodes; i++)
		for (int j = i + 1; j < Nmodes; j++)
			Emin += phi(i,j) * minima(i) * minima(j);

	// third, displacement term:
	for (int i = 0; i < Nmodes; i++)
		Emin += kappa(i) * minima(i);

	// fourth, constant term:
	for (int i = 0; i < Nmodes; i++)
		Emin += d(i);

	logFile << "Energy at minimum: " << Emin << std::endl;


	/*
	 * Calculate the 1st moment of the spectrum analytically
	 */
	double moment = dE;
	for (int i = 0; i < Nmodes; i++)
		moment += 0.25 * (fp(i) - f1(i)) / sqrt(double(f1(i)));

	std::cout << "1st moment of the spectrum: " <<std::setprecision(8) << moment << std::endl;
}


