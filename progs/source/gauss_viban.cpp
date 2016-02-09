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

	/*
	 * Read the basic information about the system from the formatted checkpoint file:
	 *
	 * "x"          is the vector of cartesian molecular coordinates in bohr (length 3N)
	 * "Masses"     contains the atomic masses in amu (length N)
	 * "f_cart"     the non-mass-weighted cartesian hessian (size 3N*3N) **** eq. (1) ****
	 * "GaussModes" the matrix of vibrational diplacements as calculated by Gaussian (size 3N*3N-6)
	 */
	std::ifstream fchkFile(filename, std::ifstream::in);
	GaussFchk fchk(fchkFile);
	int Natoms = fchk.ReadInteger("Number of atoms");
	int Ncoords = 3 * Natoms;
	int Nmodes = Ncoords - 6;
	Eigen::VectorXd x = fchk.ReadVector("Current cartesian coordinates");
	Eigen::VectorXd Masses = fchk.ReadVector("Real atomic weights");
	Eigen::MatrixXd f_cart = fchk.ReadSymmetricMatrix("Cartesian Force Constants");
	Eigen::MatrixXd GaussModes = fchk.ReadMatrix("Vib-Modes", Nmodes, Ncoords);

	/*
	 * Calculate some mass-related quantities:
	 *
	 * "TotMass"    The total mass of the molecule in amu
	 * "MassVector" Vector of masses (in amu) associated with every coordinate: (m1, m1, m1, m2, m2, m2, ... mN, mN, mN) (length 3N)
	 * "MassInvMat" Diagonal matrix of inverse square roots of masses: diag(1/sqrt(mi)) (size 3N*3N)
	 */
	double TotMass = Masses.sum();
	Eigen::VectorXd MassVector(Ncoords);
	for (int i = 0; i <  Natoms; i++)
	{
		MassVector(3 * i + 0) = Masses(i);
		MassVector(3 * i + 1) = Masses(i);
		MassVector(3 * i + 2) = Masses(i);
	}
	Eigen::MatrixXd MassInvMat = Eigen::MatrixXd::Zero(Ncoords, Ncoords);
	for (int i = 0; i < Ncoords; i++)
		MassInvMat(i,i) = 1.0 / sqrt(double(MassVector(i)));

	/*
	 * Calculate mass-weighted coordinates:
	 *
	 * qi = xi * sqrt(mi)
	 */
	Eigen::VectorXd q(Ncoords);
	for (int i = 0; i < Ncoords; i++)
		q(i) = x(i) * sqrt(double(MassVector(i)));

	/*
	 * Calculate and diagonalize the mass-weighted hessian:
	 *
	 * f_mwc(i,j) = f_cart(i,j) / sqrt(mi*mj) = M^T * f_cart * M  **** eq. (2) ****
	 * where M is the matrix of inverse masses
	 *
	 * The eigenvalues are the frequencies without
	 * the translation and rotation projected out
	 */
	Eigen::MatrixXd f_mwc = MassInvMat.transpose() * f_cart * MassInvMat;
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> HessianDiag(f_mwc);

	/*
	 * Calculate the center of mass: **** eq. (3) ****
	 */
	Eigen::Vector3d com = calc_com(x, Masses);

	/*
	 * Calculate and diagonalize the inertia tensor: **** eq. (4) ****
	 */
	Eigen::Matrix3d inert = calc_inert(x, Masses);
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> InertDiag(inert);

	/*
	 * Create the D matrix **** sec. 2.3 ****
	 * First, fill the matrix with random numbers,
	 * then put the six vectors in (translation, rotation)
	 */
	Eigen::MatrixXd D(Eigen::MatrixXd::Random(Ncoords, Ncoords));

	for (int i = 0; i < 3; i++)
		for (int j = 0; j < Natoms; j++)
		{
			D(3*j+0, i) = 0;
			D(3*j+1, i) = 0;
			D(3*j+2, i) = 0;
			D(3*j+i, i) = sqrt(double(Masses(j)));
		}

	for (int i = 0; i < Natoms; i++)
	{
		Eigen::Vector3d x_int;
		x_int(0) = x(3*i+0) - com(0);
		x_int(1) = x(3*i+1) - com(1);
		x_int(2) = x(3*i+2) - com(2);
		Eigen::Matrix3d rotmat = InertDiag.eigenvectors().transpose();
		x_int = rotmat * x_int;

		D(3*i+0, 3) = (x_int(1) * rotmat(2,0) - x_int(2) * rotmat(1,0)) * sqrt(double(Masses(i)));
		D(3*i+1, 3) = (x_int(1) * rotmat(2,1) - x_int(2) * rotmat(1,1)) * sqrt(double(Masses(i)));
		D(3*i+2, 3) = (x_int(1) * rotmat(2,2) - x_int(2) * rotmat(1,2)) * sqrt(double(Masses(i)));
		D(3*i+0, 4) = (x_int(2) * rotmat(0,0) - x_int(0) * rotmat(2,0)) * sqrt(double(Masses(i)));
		D(3*i+1, 4) = (x_int(2) * rotmat(0,1) - x_int(0) * rotmat(2,1)) * sqrt(double(Masses(i)));
		D(3*i+2, 4) = (x_int(2) * rotmat(0,2) - x_int(0) * rotmat(2,2)) * sqrt(double(Masses(i)));
		D(3*i+0, 5) = (x_int(0) * rotmat(1,0) - x_int(1) * rotmat(0,0)) * sqrt(double(Masses(i)));
		D(3*i+1, 5) = (x_int(0) * rotmat(1,1) - x_int(1) * rotmat(0,1)) * sqrt(double(Masses(i)));
		D(3*i+2, 5) = (x_int(0) * rotmat(1,2) - x_int(1) * rotmat(0,2)) * sqrt(double(Masses(i)));
	}

	// normalize the D matrix:
	for (int i = 0; i < Ncoords; i++)
		D.col(i) /= Eigen::VectorXd(D.col(i)).norm();

	// Gram-Schmidt-orthogonalize the D matrix:
	gram_schmidt(D);

	/*
	 * Transform the mass-weighted Hessian to internal coordinates using the D matrix  **** eq. (6) ****
	 * Create the 3N-6 x 3N-6 Block matrix from the lower right corner
	 */
	Eigen::MatrixXd f_int = Eigen::MatrixXd(D.transpose() * f_mwc * D).block(6, 6, Nmodes, Nmodes);

	/*
	 * Diagonalize the internal Hessian **** eq. (7) ****
	 */
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> IntHessDiag(f_int);
	Eigen::MatrixXd L = IntHessDiag.eigenvectors();

	/*
	 * Calculate the mass-weighted cartesian displacements  **** eq. (9) ****
	 * and the non-mass-weighted cartesian displacements    **** eq. (10) ****
	 */
	Eigen::MatrixXd l_mwc = D.block(0, 6, Ncoords, Nmodes) * L;
	Eigen::MatrixXd l_cart = MassInvMat * l_mwc;
	Eigen::VectorXd normalizationFactors(Nmodes);
	Eigen::VectorXd reducedMasses(Nmodes);
	for (int i = 0; i < Nmodes; i++)
	{
		double norm = Eigen::VectorXd(l_cart.col(i)).norm();
		normalizationFactors(i) = 1.0 / norm;
		reducedMasses(i) = normalizationFactors(i) * normalizationFactors(i);
		l_cart.col(i) *= normalizationFactors(i);
	}

	/*
	 * write the results to stdout
	 */
	std::cout << "Number of atoms: " << Natoms << std::endl;

	std::cout << "atomic coordinates in bohr:\n";
	WriteVector(x, digits, clean, threshold);

	std::cout << "atomic masses in u:\n";
	WriteVector(Masses, digits, clean, threshold);

	std::cout << "total mass in au: " << TotMass << std::endl;

	std::cout << "non-mass-weighted cartesian hessian\n";
	WriteMatrix(f_cart, digits, clean, threshold);

	std::cout << "mass-weighted cartesian hessian\n";
	WriteMatrix(f_mwc, digits, clean, threshold);

	std::cout << "eigenvalues of the mass-weighted cartesian hessian (au and cm-1):\n";
	for (int i = 0; i < Ncoords; i++)
	{
		double value = HessianDiag.eigenvalues()(i);
		double wavenumber;
		if (value > 0)
			wavenumber = sqrt(value * Eh / (a0 * a0 * amu)) / (2.0 * M_PI * c0 * 100.0);
		else
			wavenumber = -sqrt(-value * Eh / (a0 * a0 * amu)) / (2.0 * M_PI * c0 * 100.0);
		std::cout << std::scientific << std::setprecision(4) << std::setw(12) << value
				  << "  " << std::fixed << std::setw(9) << wavenumber << std::endl;
	}

	std::cout << "center of mass in bohr:\n";
	WriteVector(com, digits, clean, threshold);

	std::cout << "inertia tensor:\n";
	WriteMatrix(inert, digits, clean, threshold);

	std::cout << "moments of inertia:\n";
	WriteVector(InertDiag.eigenvalues(), digits, clean, threshold);

	std::cout << "principal axes:\n";
	WriteMatrix(InertDiag.eigenvectors(), digits, clean, threshold);

	std::cout << "metric of the principal axes:\n";
	WriteMatrix(InertDiag.eigenvectors().transpose() * InertDiag.eigenvectors(), digits, clean, threshold);

	std::cout << "The D matrix:\n";
	WriteMatrix(D, digits, clean, threshold);

	std::cout << "metric of D matrix:\n";
	WriteMatrix(D.transpose() * D, digits, clean, threshold);

	std::cout << "The Nvib x Nvib submatrix of the internal Hessian, according to eq. (6):\n";
	WriteMatrix(f_int, digits, clean, threshold);

	std::cout << "eigenvalues of the internal hessian (au and cm-1):\n";
	for (int i = 0; i < Nmodes; i++)
	{
		double value = IntHessDiag.eigenvalues()(i);
		double wavenumber;
		if (value > 0)
			wavenumber = sqrt(value * Eh / (a0 * a0 * amu)) / (2.0 * M_PI * c0 * 100.0);
		else
			wavenumber = -sqrt(-value * Eh / (a0 * a0 * amu)) / (2.0 * M_PI * c0 * 100.0);
		std::cout << std::scientific << std::setprecision(4) << std::setw(12) << value
				  << "  " << std::fixed << std::setw(9) << wavenumber << std::endl;
	}

	std::cout << "eigenvectors of the internal Hessian (internal displacements):\n";
	WriteMatrix(L, digits, clean, threshold);

	std::cout << "metric of the eigenvectors:\n";
	WriteMatrix(L.transpose() * L, digits, clean, threshold);

	std::cout << "mass-weighted cartesian displacements:\n";
	WriteMatrix(l_mwc, digits, clean, threshold);

	std::cout << "metric of mass-weighted cartesian displacements:\n";
	WriteMatrix(l_mwc.transpose() * l_mwc, digits, clean, threshold);

	std::cout << "cartesian displacements:\n";
	WriteMatrix(l_cart, digits, clean, threshold);

	std::cout << "metric of cartesian displacements:\n";
	WriteMatrix(l_cart.transpose() * l_cart, digits, clean, threshold);

	std::cout << "normalization factors for the cartesian displacements:\n";
	WriteVector(normalizationFactors, digits, clean, threshold);

	std::cout << "reduced masses of the normal modes:\n";
	WriteVector(reducedMasses, digits, clean, threshold);

	return 0;
}

/*
 * Orthonormalize the D matrix
 */
void gram_schmidt(Eigen::MatrixXd &d)
{
	Eigen::MatrixXd newMat = Eigen::MatrixXd::Zero(d.rows(), d.cols());
	newMat.block(0, 0, d.cols(), 6) = d.block(0, 0, d.cols(), 6);

	for (int i = 6; i < d.cols(); i++)
	{
		Eigen::VectorXd temp_vec = d.col(i);
		for (int j = 0; j < i; j++)
		{
			double factor = temp_vec.dot(newMat.col(j));
			temp_vec -= factor * newMat.col(j);
		}
		temp_vec.normalize();
		newMat.col(i) = temp_vec;
	}

	d = newMat;
}
