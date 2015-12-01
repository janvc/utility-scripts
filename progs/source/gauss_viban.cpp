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
#include "GaussFchk.h"
#include "utilities.h"
#include "constants.h"

Eigen::Vector3d calc_com(Eigen::VectorXd x, Eigen::VectorXd m);
Eigen::Matrix3d calc_inert(Eigen::VectorXd x, Eigen::VectorXd m);
void gram_schmidt(Eigen::MatrixXd &d);

int main()
{
	int digits = 4;
	bool clean = true;

	std::ifstream fchkFile("h2o_freq.fchk", std::ifstream::in);
	GaussFchk fchk(fchkFile);

	// read the data from the formatted checkpoint file:
	int Natoms = fchk.ReadInteger("Number of atoms");
	int Ncoords = 3 * Natoms;
	int Nmodes = Ncoords - 6;

	/*
	 * Read the basic information about the system from the formatted checkpoint file:
	 *
	 * "x"          is the vector of cartesian molecular coordinates in bohr (length 3N)
	 * "Masses"     contains the atomic masses in amu (length N)
	 * "f_cart"     the non-mass-weighted cartesian hessian (size 3N*3N) **** eq. (1) ****
	 * "GaussModes" the matrix of vibrational diplacements as calculated by Gaussian (size 3N*3N-6)
	 */
	Eigen::VectorXd x = fchk.ReadVector("Current cartesian coordinates");
	Eigen::VectorXd Masses = fchk.ReadVector("Real atomic weights");
	Eigen::MatrixXd f_cart = fchk.ReadSymmetricMatrix("Cartesian Force Constants");
	Eigen::MatrixXd GaussModes = fchk.ReadMatrix("Vib-Modes", Nmodes, Ncoords);

	// convert masses to atomic units:
	Masses = Masses * amu2au;
	double TotMass = Masses.sum();

	/*
	 * Calculate some mass-related quantities:
	 *
	 * "MassVector" Vector of masses (in au) associated with every coordinate: (m1, m1, m1, m2, m2, m2, ... mN, mN, mN) (length 3N)
	 * "MassInvMat" Diagonal matrix of inverse square roots of masses: diag(1/sqrt(mi)) (size 3N*3N)
	 */
	Eigen::VectorXd MassVector(Ncoords);
	for (int i = 0; i <  Natoms; i++)
	{
		MassVector(3 * i + 0) = Masses(i);
		MassVector(3 * i + 1) = Masses(i);
		MassVector(3 * i + 2) = Masses(i);
	}
	Eigen::MatrixXd MassInvMat = Eigen::MatrixXd::Zero(Ncoords, Ncoords);
	for (int i = 0; i < Ncoords; i++)
		MassInvMat(i,i) = 1.0 / sqrt(MassVector(i));

	/*
	 * Calculate mass-weighted coordinates:
	 *
	 * qi = xi * sqrt(mi)
	 */
	Eigen::VectorXd q(Ncoords);
	for (int i = 0; i < Ncoords; i++)
		q(i) = x(i) * sqrt(MassVector(i));

	/*
	 * WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE
	 * WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE
	 */
	std::cout << "Number of atoms: " << Natoms << std::endl;
	std::cout << "atomic coordinates in bohr:\n";
	WriteVector(x, digits, clean);
	std::cout << "atomic masses in au:\n";
	WriteVector(Masses, digits, clean);
	std::cout << "total mass in au: " << TotMass << std::endl;
	std::cout << "mass vector:\n";
	WriteVector(MassVector, digits, clean);
	std::cout << "mass-weighted coordinates:\n";
	WriteVector(q, digits, clean);
	std::cout << "inverse mass matrix:\n";
	WriteMatrix(MassInvMat, digits, clean);
	std::cout << "non-mass-weighted cartesian hessian\n";
	WriteMatrix(f_cart, digits, clean);
	std::cout << "Gaussian Modes\n";
	WriteMatrix(GaussModes, digits, clean);
	/*
	 * WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE
	 * WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE
	 */

	/*
	 * Calculate and diagonalize
	 * the mass-weighted hessian:
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
	 * WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE
	 * WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE
	 */
	std::cout << "mass-weighted cartesian hessian:\n";
	WriteMatrix(f_mwc, digits, clean);
	std::cout << "eigenvalues of the mass-weighted cartesian hessian (au and cm-1):\n";
	for (int i = 0; i < Ncoords; i++)
	{
		double value = HessianDiag.eigenvalues()(i);
		double wavenumber = std::sqrt(std::abs(value) * Eh / (a0 * a0 * me)) / (2.0 * M_PI * c0 * 100.0);
		std::cout << std::scientific << std::setprecision(4) << std::setw(12) << value
				  << "  " << std::fixed << std::setw(9) << wavenumber << std::endl;
	}
	/*
	 * WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE
	 * WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE
	 */

	/*
	 * Calculate the center of mass: **** eq. (3) ****
	 */
	Eigen::Vector3d com = calc_com(x, Masses);

	/*
	 * WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE
	 * WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE
	 */
	std::cout << "center of mass in bohr:\n";
	WriteVector(com, digits, clean);
	/*
	 * WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE
	 * WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE
	 */

	// Translate the molecule to the center of mass:
	for (int i = 0; i < Natoms; i++)
	{
		x(3 * i + 0) -= com(0);
		x(3 * i + 1) -= com(1);
		x(3 * i + 2) -= com(2);
	}

	/*
	 * Calculate and diagonalize the inertia tensor: **** eq. (4) ****
	 */
	Eigen::Matrix3d inert = calc_inert(x, Masses);
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> InertDiag(inert);

	/*
	 * WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE
	 * WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE
	 */
	std::cout << "inertia tensor:\n";
	WriteMatrix(inert, digits, clean);
	std::cout << "moments of inertia:\n";
	WriteVector(InertDiag.eigenvalues(), digits, clean);
	std::cout << "principal axes:\n";
	WriteMatrix(InertDiag.eigenvectors(), digits, clean);
	std::cout << "metric of the principal axes:\n";
	WriteMatrix(InertDiag.eigenvectors().transpose() * InertDiag.eigenvectors(), digits, clean);
	/*
	 * WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE
	 * WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE
	 */

	// rotate the molecule to the inertia frame:
	Eigen::VectorXd y(Ncoords);
	Eigen::Matrix3d rotmat;
	rotmat = InertDiag.eigenvectors();
	for (int i = 0; i < Natoms; i++)
	{
		y(3 * i + 0) = rotmat.transpose()(0,0) * x(3 * i + 0)
				     + rotmat.transpose()(0,1) * x(3 * i + 1)
					 + rotmat.transpose()(0,2) * x(3 * i + 2);
		y(3 * i + 1) = rotmat.transpose()(1,0) * x(3 * i + 0)
					 + rotmat.transpose()(1,1) * x(3 * i + 1)
					 + rotmat.transpose()(1,2) * x(3 * i + 2);
		y(3 * i + 2) = rotmat.transpose()(2,0) * x(3 * i + 0)
				     + rotmat.transpose()(2,1) * x(3 * i + 1)
					 + rotmat.transpose()(2,2) * x(3 * i + 2);
	}
	x = y;
	/*
	 * WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE
	 * WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE
	 */
	std::cout << "atomic coordinates after transforming to the inertial frame:\n";
	WriteVector(x, digits, clean);
	/*
	 * WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE
	 * WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE
	 */

	// calculate mass-weighted coordinates again:
	for (int i = 0; i < Ncoords; i++)
		q(i) = x(i) * sqrt(MassVector(i));

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
			D(3*j+i, i) = sqrt(Masses(j));
		}

	for (int i = 0; i < Natoms; i++)
	{
		D(3*i+0, 3) =  0;
		D(3*i+1, 3) = -sqrt(Masses(i)) * (x(3*i+0) * rotmat(2,0)
										+ x(3*i+1) * rotmat(2,1)
										+ x(3*i+2) * rotmat(2,2));
		D(3*i+2, 3) =  sqrt(Masses(i)) * (x(3*i+0) * rotmat(1,0)
										+ x(3*i+1) * rotmat(1,1)
										+ x(3*i+2) * rotmat(1,2));
		D(3*i+0, 4) =  sqrt(Masses(i)) * (x(3*i+0) * rotmat(2,0)
										+ x(3*i+1) * rotmat(2,1)
										+ x(3*i+2) * rotmat(2,2));
		D(3*i+1, 4) =  0;
		D(3*i+2, 4) = -sqrt(Masses(i)) * (x(3*i+0) * rotmat(0,0)
										+ x(3*i+1) * rotmat(0,1)
										+ x(3*i+2) * rotmat(0,2));
		D(3*i+0, 5) = -sqrt(Masses(i)) * (x(3*i+0) * rotmat(1,0)
										+ x(3*i+1) * rotmat(1,1)
										+ x(3*i+2) * rotmat(1,2));
		D(3*i+1, 5) =  sqrt(Masses(i)) * (x(3*i+0) * rotmat(0,0)
										+ x(3*i+1) * rotmat(0,1)
										+ x(3*i+2) * rotmat(0,2));
		D(3*i+2, 5) =  0;
	}

	// normalize the D matrix:
	for (int i = 0; i < Ncoords; i++)
		D.col(i) /= Eigen::VectorXd(D.col(i)).norm();

	// Gram-Schmidt-orthogonalize the D matrix:
	gram_schmidt(D);

	/*
	 * WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE
	 * WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE
	 */
	std::cout << "The D matrix:\n";
	WriteMatrix(D, digits, clean);
	std::cout << "metric of D matrix:\n";
	WriteMatrix(D.transpose() * D, digits, clean);
	/*
	 * WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE
	 * WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE
	 */

	/*
	 * Transform the mass-weighted Hessian to internal coordinates using the D matrix  **** eq. (6) ****
	 */
	Eigen::MatrixXd f_int = D.transpose() * f_mwc * D;

	/*
	 * Diagonalize the internal Hessian **** eq. (7) ****
	 */
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> IntHessDiag(f_int);
	Eigen::MatrixXd L = IntHessDiag.eigenvectors();

	/*
	 * WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE
	 * WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE
	 */
	std::cout << "The Hessian transformed to internal coordinates:\n";
	WriteMatrix(f_int, digits, clean);
	std::cout << "eigenvalues of the internal hessian (au and cm-1):\n";
	for (int i = 0; i < Ncoords; i++)
	{
		double value = IntHessDiag.eigenvalues()(i);
		double wavenumber = std::sqrt(std::abs(value) * Eh / (a0 * a0 * me)) / (2.0 * M_PI * c0 * 100.0);
		std::cout << std::scientific << std::setprecision(4) << std::setw(12) << value
				  << "  " << std::fixed << std::setw(9) << wavenumber << std::endl;
	}
	std::cout << "Eigenvectors of the internal Hessian:\n";
	WriteMatrix(L, digits, clean);
	std::cout << "metric of the eigenvectors:\n";
	WriteMatrix(L.transpose() * L, digits, clean);
	/*
	 * WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE
	 * WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE
	 */

	/*
	 * Calculate the mass-weighted cartesian displacements  **** eq. (9) ****
	 * and the non-mass-weighted cartesian displacements  **** eq. (10) ****
	 */
	Eigen::MatrixXd l_mwc = D * L;
	Eigen::MatrixXd l_cart = MassInvMat * l_mwc;

	/*
	 * WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE
	 * WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE
	 */
	std::cout << "mass-weighted displacements:\n";
	WriteMatrix(l_mwc, digits, clean);
	std::cout << "metric of mass-weighted displacements:\n";
	WriteMatrix(l_mwc.transpose() * l_mwc, digits, clean);
	std::cout << "cartesian displacements:\n";
	WriteMatrix(l_cart, digits, clean);
	std::cout << "metric of cartesian displacements:\n";
	WriteMatrix(l_cart.transpose() * l_cart, digits, clean);
	/*
	 * WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE
	 * WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE WRITE
	 */

	return 0;
}

/*
 * Calculates the center of mass of the molecule
 */
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

/*
 * Calculates the moment of inertia of the molecule
 */
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





