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
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/SVD>
#include <eigen3/Eigen/LU>
#include <eigen3/Eigen/Eigenvalues>
#include <boost/program_options.hpp>
#include "GaussFchk.h"
#include "utilities.h"
#include "constants.h"


int main(int argc, char *argv[])
{
	const int digits = 5;
	const bool clean = true;
	const double threshold = 1.0e-15;


	/*
	 * Process command line options:
	 */
	std::string GS_filename;
	std::string ES_filename;
	namespace po = boost::program_options;
	po::options_description desc("Allowed options");
	desc.add_options()
		("help,h", "produce this help message")
		("gsfile,g", po::value<std::string>(&GS_filename)->required(), "the ground state FChk file")
		("esfile,e", po::value<std::string>(&ES_filename)->required(), "the excited state FChk file")
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
	 * Read information from the FChk files:
	 */
	std::ifstream GS_fchkFile(GS_filename, std::ifstream::in);
	std::ifstream ES_fchkFile(ES_filename, std::ifstream::in);
	GaussFchk GS_fchk(GS_fchkFile);
	GaussFchk ES_fchk(ES_fchkFile);

	int Natoms = GS_fchk.ReadInteger("Number of atoms");
	int Ncoords = 3 * Natoms;
	int Nmodes = Ncoords - 6;

	double Eg = GS_fchk.ReadReal("Total Energy");
	double Ee = ES_fchk.ReadReal("Total Energy");
	double deltaE = Ee - Eg;

	Eigen::VectorXd masses = GS_fchk.ReadVector("Real atomic weights") * amu2au;
	Eigen::VectorXd Xg = GS_fchk.ReadVector("Current cartesian coordinates");
	Eigen::VectorXd Xe = ES_fchk.ReadVector("Current cartesian coordinates");
	Eigen::MatrixXd Hg = GS_fchk.ReadSymmetricMatrix("Cartesian Force Constants");
	Eigen::MatrixXd He = ES_fchk.ReadSymmetricMatrix("Cartesian Force Constants");


	/*
	 * Create mass-weighted coordinates:
	 */
	Eigen::VectorXd Qg(Ncoords);
	Eigen::VectorXd Qe(Ncoords);
	for (int i = 0; i < Natoms; i++)
	{
		Qg(3 * i + 0) = Xg(3 * i + 0) * sqrt(double(masses(i)));
		Qe(3 * i + 0) = Xe(3 * i + 0) * sqrt(double(masses(i)));
		Qg(3 * i + 1) = Xg(3 * i + 1) * sqrt(double(masses(i)));
		Qe(3 * i + 1) = Xe(3 * i + 1) * sqrt(double(masses(i)));
		Qg(3 * i + 2) = Xg(3 * i + 2) * sqrt(double(masses(i)));
		Qe(3 * i + 2) = Xe(3 * i + 2) * sqrt(double(masses(i)));
	}


	/*
	 * Begin writing the log file:
	 */
	std::ofstream logFile("log");

	logFile << "Number of atoms:      " << Natoms << std::endl;
	logFile << " Ground state total energy"
			<< std::setw(16) << std::setprecision(12) << Eg << " Eh\n";
	logFile << "Excited state total energy"
			<< std::setw(16) << std::setprecision(12) << Ee << " Eh\n";
	logFile << "adiabatic excitation energy:   "
			<< std::setprecision(9) << deltaE << " Eh\n";
	logFile << "                             = " << deltaE * Eh2eV << "  eV\n\n";

	logFile << "Ground and excited state cartesian coordinates in [bohr]:\n";
	for (int i = 0; i < Ncoords; i++)
		logFile << std::setw(5) << i + 1
				<< std::setw(17) << std::fixed << double(Xg(i))
				<< std::setw(17) << double(Xe(i)) << std::endl;

	logFile << "Ground and excited state mass-weighted coordinates in [bohr * sqrt(me)]:\n";
	for (int i = 0; i < Ncoords; i++)
		logFile << std::setw(5) << i + 1
				<< std::setw(17) << std::fixed << double(Qg(i))
				<< std::setw(17) << double(Qe(i)) << std::endl;
	logFile << std::endl;

	logFile << "Ground state cartesian Hessian in [Eh / bohr**2]:\n";
	WriteMatrixToFile(logFile, Hg, digits, clean, threshold);
	logFile << "Excited state cartesian Hessian in [Eh / bohr**2]:\n";
	WriteMatrixToFile(logFile, He, digits, clean, threshold);


	/*
	 * Perform vibrational analysis:
	 */
	Eigen::MatrixXd vibAn_g = vibrationalAnalysis(Xg, masses, Hg);
	Eigen::MatrixXd vibAn_e = vibrationalAnalysis(Xe, masses, He);
	Eigen::VectorXd FCg = vibAn_g.row(0);
	Eigen::VectorXd FCe = vibAn_e.row(0);
	Eigen::MatrixXd NMg = vibAn_g.block(1, 0, Ncoords, Nmodes);
	Eigen::MatrixXd NMe = vibAn_e.block(1, 0, Ncoords, Nmodes);

	Eigen::VectorXd freq_g(Nmodes);
	Eigen::VectorXd freq_e(Nmodes);
	for (int i = 0; i < Nmodes; i++)
	{
		double value_g = FCg(i);
		double value_e = FCe(i);
		double wavenumber_g, wavenumber_e;
		if (value_g > 0)
			wavenumber_g = sqrt(value_g * Eh / (a0 * a0 * me)) / (2.0 * M_PI * c0 * 100.0);
		else
			wavenumber_g = -sqrt(-value_g * Eh / (a0 * a0 * me)) / (2.0 * M_PI * c0 * 100.0);
		if (value_e > 0)
			wavenumber_e = sqrt(value_e * Eh / (a0 * a0 * me)) / (2.0 * M_PI * c0 * 100.0);
		else
			wavenumber_e = -sqrt(-value_e * Eh / (a0 * a0 * me)) / (2.0 * M_PI * c0 * 100.0);
		freq_g(i) = wavenumber_g;
		freq_e(i) = wavenumber_e;
	}

	logFile << std::endl;
	logFile << "Ground state mass-weighted force constants in [Eh / a0**2 * me]" << std::endl
			<< "and frequencies in [cm**-1]:\n";
	for (int i = 0; i < Nmodes; i++)
		logFile << std::setw(5) << i + 1
				<< std::setw(20) << std::scientific << std::setprecision(9) << double(FCg(i))
				<< std::setw(13) << std::fixed << std::setprecision(4) << double(freq_g(i)) << std::endl;
	logFile << "Excited state mass-weighted force constants in [Eh / a0**2 * me]" << std::endl
			<< "and frequencies in [cm**-1]:\n";
	for (int i = 0; i < Nmodes; i++)
		logFile << std::setw(5) << i + 1
				<< std::setw(20) << std::scientific << std::setprecision(9) << double(FCe(i))
				<< std::setw(13) << std::fixed << std::setprecision(4) << double(freq_e(i)) << std::endl;

	logFile << std::endl;
	logFile << "Ground state mass-weighted normal modes in [a0 * sqrt(me)]:\n";
	WriteMatrixToFile(logFile, NMg, digits, clean, threshold);
	logFile << "Excited state mass-weighted normal modes in [a0 * sqrt(me)]:\n";
	WriteMatrixToFile(logFile, NMe, digits, clean, threshold);
	logFile << "Metric of the mass-weighted ground state normal modes:\n";
	WriteMatrixToFile(logFile, NMg.transpose() * NMg, digits, clean, threshold);
	logFile << "Metric of the mass-weighted excited state normal modes:\n";
	WriteMatrixToFile(logFile, NMe.transpose() * NMe, digits, clean, threshold);
	logFile << std::endl;


	/*
	 * Perform SVD on the normal modes:
	 */
	Eigen::JacobiSVD<Eigen::MatrixXd> SVDg(NMg, Eigen::ComputeThinU | Eigen::ComputeThinV);
	Eigen::JacobiSVD<Eigen::MatrixXd> SVDe(NMe, Eigen::ComputeThinU | Eigen::ComputeThinV);
	Eigen::MatrixXd NMOg = SVDg.matrixU() * SVDg.matrixV().transpose();
	Eigen::MatrixXd NMOe = SVDe.matrixU() * SVDe.matrixV().transpose();

	logFile << "Ground state mass-weighted normal modes after Loewdin orthogonalization in [a0 * sqrt(me)]:\n";
	WriteMatrixToFile(logFile, NMOg, digits, clean, threshold);
	logFile << "Excited state mass-weighted normal modes after Loewdin orthogonalization in [a0 * sqrt(me)]:\n";
	WriteMatrixToFile(logFile, NMOe, digits, clean, threshold);
	logFile << "Metric of the orthogonalized ground state mass-weighted normal modes:\n";
	WriteMatrixToFile(logFile, NMOg.transpose() * NMOg, digits, clean, threshold);
	logFile << "Metric of the orthogonalized excited state mass-weighted normal modes:\n";
	WriteMatrixToFile(logFile, NMOe.transpose() * NMOe, digits, clean, threshold);
	logFile << "Difference between the original and orthogonalized ground state normal modes:\n";
	WriteMatrixToFile(logFile, NMOg - NMg, digits, clean, threshold);
	logFile << "Difference between the original and orthogonalized excited state normal modes:\n";
	WriteMatrixToFile(logFile, NMOe - NMe, digits, clean, threshold);
	logFile << std::endl;


	/*
	 * calculate the rotation matrix to minimize the RMSD between
	 * the ground- and excited-state structures:
	 */

	// shift the molecules to their centers of mass:
	Eigen::Vector3d COMg(calc_com(Xg, masses));
	Eigen::Vector3d COMe(calc_com(Xe, masses));

	logFile << "Center of mass of the ground and excited state:\n";
	for (int i = 0; i < 3; i++)
		logFile << std::setw(16) << std::scientific << std::setprecision(5) << double(COMg(i))
				<< std::setw(16) << std::scientific << std::setprecision(5) << double(COMe(i)) << std::endl;

	Eigen::VectorXd XSg(Ncoords);
	Eigen::VectorXd XSe(Ncoords);
	for (int i = 0; i < Natoms; i++)
		for (int j = 0; j < 3; j++)
		{
			XSg(3*i+j) = Xg(3*i+j) - COMg(j);
			XSe(3*i+j) = Xe(3*i+j) - COMe(j);
		}
	Eigen::VectorXd QSg(Ncoords);
	Eigen::VectorXd QSe(Ncoords);
	for (int i = 0; i < Natoms; i++)
	{
		QSg(3 * i + 0) = XSg(3 * i + 0) * sqrt(double(masses(i)));
		QSe(3 * i + 0) = XSe(3 * i + 0) * sqrt(double(masses(i)));
		QSg(3 * i + 1) = XSg(3 * i + 1) * sqrt(double(masses(i)));
		QSe(3 * i + 1) = XSe(3 * i + 1) * sqrt(double(masses(i)));
		QSg(3 * i + 2) = XSg(3 * i + 2) * sqrt(double(masses(i)));
		QSe(3 * i + 2) = XSe(3 * i + 2) * sqrt(double(masses(i)));
	}

	COMg = calc_com(XSg, masses);
	COMe = calc_com(XSe, masses);

	logFile << "Center of mass of the ground and excited state after shifting:\n";
	for (int i = 0; i < 3; i++)
		logFile << std::setw(16) << std::scientific << std::setprecision(5) << double(COMg(i))
				<< std::setw(16) << std::scientific << std::setprecision(5) << double(COMe(i)) << std::endl;

	// calculate the original RMSD:
	double RMSD_before = 0;
	for (int i = 0; i < Ncoords; i++)
		RMSD_before += std::abs(double(XSg(i) - XSe(i))) * std::abs(double(XSg(i) - XSe(i)));
	RMSD_before /= Natoms;

	logFile << "original RMSD between ground and excited state: " << RMSD_before << std::endl;

	// calculate the cross-correlation matrix between the two structures:
	Eigen::Matrix3d corr(Eigen::Matrix3d::Zero());
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			for(int k = 0; k < Natoms; k++)
				corr(i,j) += XSg(3*k+i) * XSe(3*k+j);

	logFile << "the cross-correlation matrix between the ground and excited state:\n";
	WriteMatrixToFile(logFile, corr, digits, clean, threshold);
	logFile << "Determinant of the correlation matrix: " << corr.determinant() << std::endl;


	// calculate the quaternion matrix:
	Eigen::Matrix4d F(Eigen::Matrix4d::Zero());
	F(0,0) =  corr(0,0) + corr(1,1) + corr(2,2);
	F(1,1) =  corr(0,0) - corr(1,1) - corr(2,2);
	F(2,2) = -corr(0,0) + corr(1,1) - corr(2,2);
	F(3,3) = -corr(0,0) - corr(1,1) + corr(2,2);
	F(0,1) =  corr(1,2) - corr(2,1);
	F(0,2) =  corr(2,0) - corr(0,2);
	F(0,3) =  corr(0,1) - corr(1,0);
	F(1,2) =  corr(0,1) + corr(1,0);
	F(1,3) =  corr(0,2) + corr(2,0);
	F(2,3) =  corr(1,2) + corr(2,1);
	F(1,0) = F(0,1);
	F(2,0) = F(0,2);
	F(3,0) = F(0,3);
	F(2,1) = F(1,2);
	F(3,1) = F(1,3);
	F(3,2) = F(2,3);

	// diagonalize the quaternion matrix and select the best rotation quaternion:
	// (the one corresponding to the largest absolute eigenvalue)
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix4d> Feig(F);

	Eigen::Vector4d lQuart = abs(double(Feig.eigenvalues()(0))) > abs(double(Feig.eigenvalues()(3)))
							? Feig.eigenvectors().block(0, 0, 4, 1) : Feig.eigenvectors().block(0, 3, 4, 1);

	Eigen::Matrix3d rotmat(Eigen::Matrix3d::Zero());

	rotmat(0,0) = lQuart(0)*lQuart(0) + lQuart(1)*lQuart(1) - lQuart(2)*lQuart(2) - lQuart(3)*lQuart(3);
	rotmat(0,1) = 2.0 * (lQuart(1)*lQuart(2) - lQuart(0)*lQuart(3));
	rotmat(0,2) = 2.0 * (lQuart(1)*lQuart(3) + lQuart(0)*lQuart(2));
	rotmat(1,0) = 2.0 * (lQuart(1)*lQuart(2) + lQuart(0)*lQuart(3));
	rotmat(1,1) = lQuart(0)*lQuart(0) - lQuart(1)*lQuart(1) + lQuart(2)*lQuart(2) - lQuart(3)*lQuart(3);
	rotmat(1,2) = 2.0 * (lQuart(2)*lQuart(3) - lQuart(0)*lQuart(1));
	rotmat(2,0) = 2.0 * (lQuart(1)*lQuart(3) - lQuart(0)*lQuart(2));
	rotmat(2,1) = 2.0 * (lQuart(2)*lQuart(3) + lQuart(0)*lQuart(1));
	rotmat(2,2) = lQuart(0)*lQuart(0) - lQuart(1)*lQuart(1) - lQuart(2)*lQuart(2) + lQuart(3)*lQuart(3);

	Eigen::MatrixXd BigRotMat(Eigen::MatrixXd::Zero(Ncoords, Ncoords));
		for (int i = 0; i < Natoms; i++)
			BigRotMat.block(3*i, 3*i, 3, 3) = rotmat;

	logFile << "Quaternion matrix resulting from the correlation matrix:\n";
	WriteMatrixToFile(logFile, F, digits, clean, threshold);
	logFile << "Eigenvalues of the Quaternion matrix:\n";
	WriteVectorToFile(logFile, Feig.eigenvalues(), digits, clean, threshold);
	logFile << "Eigenvectors of the Quaternion matrix:\n";
	WriteMatrixToFile(logFile, Feig.eigenvectors(), digits, clean, threshold);
	logFile << "The best-rotation Eigen-Quaternion:\n";
	WriteVectorToFile(logFile, lQuart, digits, clean, threshold);
	logFile << "The optimal rotation matrix:\n";
	WriteMatrixToFile(logFile, rotmat, digits, clean, threshold);
	logFile << "The big rotation matrix:\n";
	WriteMatrixToFile(logFile, BigRotMat, digits, clean, threshold);

	// rotate the structures
	Eigen::VectorXd XRg(BigRotMat * XSg);
	Eigen::VectorXd QRg(BigRotMat * QSg);

	double RMSD_after = 0;
	for (int i = 0; i < Ncoords; i++)
		RMSD_after += std::abs(double(XRg(i) - Xe(i))) * std::abs(double(XRg(i) - Xe(i)));
	RMSD_after /= Natoms;

	logFile << "new RMSD after Quaternion rotation: " << RMSD_after << std::endl;

	Eigen::MatrixXd J(Eigen::MatrixXd::Zero(Nmodes, Nmodes));
	Eigen::VectorXd K(Eigen::VectorXd::Zero(Nmodes));

	// check if the rotation actually improved the alignment
	// and calculate the Duschinsky matrix and the displacement vector:
	if (RMSD_after < RMSD_before)
	{
		logFile << "The rotation has reduced the RMSD, so we will use the rotated structures\n";

		Eigen::MatrixXd NMRg(BigRotMat * NMOg);

		logFile << "the rotated ground state mass-weighted normal modes:\n";
		WriteMatrixToFile(logFile, NMRg, digits, clean, threshold);
		logFile << "and their metric:\n";
		WriteMatrixToFile(logFile, NMRg.transpose() * NMRg, digits, clean, threshold);

		J = NMOe.transpose() * NMRg;
		K = NMOe.transpose() * (QRg - QSe);
	}
	else
	{
		logFile << "The rotation did not reduce the RMSD, so we will use the original structures\n";

		J = NMOe.transpose() * NMOg;
		K = NMOe.transpose() * (QSg - QSe);
	}

	logFile << "Difference in equilibrium coordinates:\n";
	WriteVectorToFile(logFile, QSg - QSe, digits, clean, threshold);

	logFile << "Duschinsky Matrix J:\n";
	WriteMatrixToFile(logFile, J, digits, clean, threshold);
	logFile << "Displacement Vector K:\n";
	WriteVectorToFile(logFile, K, digits, clean, threshold);

	logFile << "Metric of the Duschinsky matrix:\n";
	WriteMatrixToFile(logFile, J.transpose() * J, digits, clean, threshold);
	logFile << "Determinant of the Duschinsky matrix: " << J.determinant() << std::endl;


	/*
	 * Write MCTDH files:
	 */
	createMCTDHfiles(J, K, FCg, FCe, deltaE, logFile);

	return 0;
}
