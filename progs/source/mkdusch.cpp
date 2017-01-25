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
#include <boost/filesystem.hpp>
#include "GaussFchk.h"
#include "utilities.h"
#include "constants.h"
#include "vibrationalanalysis.h"


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
	std::string basename;
	double paraThres = 0.0;
	bool noJ = false;
	bool noK = false;
	namespace po = boost::program_options;
	po::options_description desc("Allowed options");
	desc.add_options()
		("help,h", "produce this help message")
		("gsfile,g", po::value<std::string>(&GS_filename)->required(), "the ground state FChk file")
		("esfile,e", po::value<std::string>(&ES_filename)->required(), "the excited state FChk file")
		("mctdh,m", po::value<std::string>(&basename)->required(), "base name of the MCTDH files")
		("threshold,t",po::value<double>(&paraThres), "threshold for including couplings/shifts")
		("noj", "ignore Duschinsky rotation (assume unit matrix)")
		("nok", "ignore displacement vector (assume zero vector)")
	;
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	if (vm.count("help"))
	{
		std::cout << desc << std::endl;
		return 1;
	}
	if (vm.count("noj"))
		noJ = true;
	if (vm.count("nok"))
		noK = true;
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

	VibrationalAnalysis GState(GS_fchk);
	VibrationalAnalysis EState(ES_fchk);


	/*
	 * Begin writing the log file:
	 */
	std::ofstream logFile("log");

	logFile << " Ground state fchk file:    " << GS_filename << std::endl;
	logFile << "Excited state fchk file:    " << ES_filename << std::endl;
	logFile << "MCTDH basename:             " << basename << std::endl;
	logFile << "Mode threshold:             " << paraThres << std::endl << std::endl;

	logFile << "Number of atoms:              " << std::setw(4) << Natoms << std::endl;
	logFile << " Ground state total energy  "
			<< std::setw(16) << std::setprecision(12) << Eg << " Eh\n";
	logFile << "Excited state total energy  "
			<< std::setw(16) << std::setprecision(12) << Ee << " Eh\n";
	logFile << "adiabatic excitation energy:     "
			<< std::setprecision(9) << deltaE << " Eh\n";
	logFile << "                               = " << deltaE * Eh2eV << "  eV\n\n";

	logFile << "Ground and excited state cartesian coordinates in [bohr]:\n";
	for (int i = 0; i < Ncoords; i++)
		logFile << std::setw(5) << i + 1
				<< std::setw(17) << std::fixed << double(GState.X()(i))
				<< std::setw(17) << double(EState.X()(i)) << std::endl;

	logFile << "Ground and excited state mass-weighted coordinates in [bohr * sqrt(me)]:\n";
	for (int i = 0; i < Ncoords; i++)
		logFile << std::setw(5) << i + 1
				<< std::setw(17) << std::fixed << double(GState.Q()(i))
				<< std::setw(17) << double(EState.Q()(i)) << std::endl;
	logFile << std::endl;

	logFile << "Ground state cartesian Hessian in [Eh / bohr**2]:\n";
    Utils::WriteSymmGaussMatrixToFile(logFile, GState.Fcart());
	logFile << "\nExcited state cartesian Hessian in [Eh / bohr**2]:\n";
    Utils::WriteSymmGaussMatrixToFile(logFile, EState.Fcart());

	double zpeGtot = 0.0;
	double zpeEtot = 0.0;
	for (int i = 0; i < Nmodes; i++)
	{
		zpeGtot += 0.5 * sqrt(double(GState.intFC()(i)));
		zpeEtot += 0.5 * sqrt(double(EState.intFC()(i)));
	}

	logFile << std::endl;
	logFile << "Ground state mass-weighted force constants in [Eh / a0**2 * me]" << std::endl
			<< "and frequencies in [cm**-1]:\n";
	for (int i = 0; i < Nmodes; i++)
		logFile << std::setw(5) << i + 1
				<< std::setw(20) << std::scientific << std::setprecision(9) << double(GState.intFC()(i))
				<< std::setw(13) << std::fixed << std::setprecision(4) << double(GState.intFrq()(i)) << std::endl;

	logFile << "\nGround state zero-point energy:"
			<< std::setw(20) << std::scientific << std::setprecision(9) << zpeGtot << " Eh, "
			<< std::setw(10) << std::fixed << std::setprecision(4) << zpeGtot * Eh2eV << " eV"<< std::endl;

	logFile << "\nExcited state mass-weighted force constants in [Eh / a0**2 * me]" << std::endl
			<< "and frequencies in [cm**-1]:\n";
	for (int i = 0; i < Nmodes; i++)
		logFile << std::setw(5) << i + 1
				<< std::setw(20) << std::scientific << std::setprecision(9) << double(EState.intFC()(i))
				<< std::setw(13) << std::fixed << std::setprecision(4) << double(EState.intFrq()(i)) << std::endl;

	logFile << "\nExcited state zero-point energy:"
				<< std::setw(20) << std::scientific << std::setprecision(9) << zpeEtot << " Eh, "
				<< std::setw(10) << std::fixed << std::setprecision(4) << zpeEtot * Eh2eV << " eV"<< std::endl;

	logFile << std::endl;
	logFile << "Ground state mass-weighted normal modes in [a0 * sqrt(me)]:\n";
    Utils::WriteGaussMatrixToFile(logFile, GState.Lmwc());
	logFile << "\nExcited state mass-weighted normal modes in [a0 * sqrt(me)]:\n";
    Utils::WriteGaussMatrixToFile(logFile, EState.Lmwc());
	logFile << "\nMetric of the mass-weighted ground state normal modes:\n";
    Utils::WriteSymmGaussMatrixToFile(logFile, GState.Lmwc().transpose() * GState.Lmwc());
	logFile << "\nMetric of the mass-weighted excited state normal modes:\n";
    Utils::WriteSymmGaussMatrixToFile(logFile, EState.Lmwc().transpose() * EState.Lmwc());
	logFile << std::endl;


	/*
	 * Perform SVD on the normal modes:
	 */
	Eigen::MatrixXd NMg = GState.Lmwc();
	Eigen::MatrixXd NMe = EState.Lmwc();
	Eigen::JacobiSVD<Eigen::MatrixXd> SVDg(NMg, Eigen::ComputeThinU | Eigen::ComputeThinV);
	Eigen::JacobiSVD<Eigen::MatrixXd> SVDe(NMe, Eigen::ComputeThinU | Eigen::ComputeThinV);
	Eigen::MatrixXd NMOg = SVDg.matrixU() * SVDg.matrixV().transpose();
	Eigen::MatrixXd NMOe = SVDe.matrixU() * SVDe.matrixV().transpose();

	logFile << "Ground state mass-weighted normal modes after Loewdin orthogonalization in [a0 * sqrt(me)]:\n";
    Utils::WriteGaussMatrixToFile(logFile, NMOg);
	logFile << "\nExcited state mass-weighted normal modes after Loewdin orthogonalization in [a0 * sqrt(me)]:\n";
    Utils::WriteGaussMatrixToFile(logFile, NMOe);
	logFile << "\nMetric of the orthogonalized ground state mass-weighted normal modes:\n";
    Utils::WriteSymmGaussMatrixToFile(logFile, NMOg.transpose() * NMOg);
	logFile << "\nMetric of the orthogonalized excited state mass-weighted normal modes:\n";
    Utils::WriteSymmGaussMatrixToFile(logFile, NMOe.transpose() * NMOe);
	logFile << "\nDifference between the original and orthogonalized ground state normal modes:\n";
    Utils::WriteGaussMatrixToFile(logFile, NMOg - NMg);
	logFile << "\nDifference between the original and orthogonalized excited state normal modes:\n";
    Utils::WriteGaussMatrixToFile(logFile, NMOe - NMe);
	logFile << std::endl;


	/*
	 * calculate the rotation matrix to minimize the RMSD between
	 * the ground- and excited-state structures:
	 */

	// shift the molecules to their centers of mass:
	Eigen::Vector3d COMg = GState.com();
	Eigen::Vector3d COMe = EState.com();

	logFile << "Center of mass of the ground and excited state:\n";
	logFile << "     ground state    excited state\n";
	for (int i = 0; i < 3; i++)
		logFile << std::setw(16) << std::scientific << std::setprecision(5) << double(COMg(i))
				<< std::setw(16) << std::scientific << std::setprecision(5) << double(COMe(i)) << std::endl;

	Eigen::VectorXd XSg(Ncoords);
	Eigen::VectorXd XSe(Ncoords);
	for (int i = 0; i < Natoms; i++)
		for (int j = 0; j < 3; j++)
		{
			XSg(3*i+j) = GState.X()(3*i+j) - COMg(j);
			XSe(3*i+j) = EState.X()(3*i+j) - COMe(j);
		}
	Eigen::VectorXd QSg(Ncoords);
	Eigen::VectorXd QSe(Ncoords);
	Eigen::VectorXd masses = GState.masses();
	for (int i = 0; i < Natoms; i++)
	{
		QSg(3 * i + 0) = XSg(3 * i + 0) * sqrt(double(masses(i)));
		QSe(3 * i + 0) = XSe(3 * i + 0) * sqrt(double(masses(i)));
		QSg(3 * i + 1) = XSg(3 * i + 1) * sqrt(double(masses(i)));
		QSe(3 * i + 1) = XSe(3 * i + 1) * sqrt(double(masses(i)));
		QSg(3 * i + 2) = XSg(3 * i + 2) * sqrt(double(masses(i)));
		QSe(3 * i + 2) = XSe(3 * i + 2) * sqrt(double(masses(i)));
	}

    COMg = Utils::calc_com(XSg, masses);
    COMe = Utils::calc_com(XSe, masses);

	logFile << "\nCenter of mass of the ground and excited state after shifting:\n";
	logFile << "     ground state    excited state\n";
	for (int i = 0; i < 3; i++)
		logFile << std::setw(16) << std::scientific << std::setprecision(5) << double(COMg(i))
				<< std::setw(16) << std::scientific << std::setprecision(5) << double(COMe(i)) << std::endl;

	// calculate the original RMSD:
	double RMSD_before = 0;
	for (int i = 0; i < Ncoords; i++)
		RMSD_before += std::abs(double(XSg(i) - XSe(i))) * std::abs(double(XSg(i) - XSe(i)));
	RMSD_before /= Natoms;

	logFile << "\noriginal RMSD between ground and excited state:    " << RMSD_before << std::endl << std::endl;

	// calculate the cross-correlation matrix between the two structures:
	Eigen::Matrix3d corr(Eigen::Matrix3d::Zero());
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			for(int k = 0; k < Natoms; k++)
				corr(i,j) += XSg(3*k+i) * XSe(3*k+j);

	logFile << "the cross-correlation matrix between the ground and excited state:\n";
    Utils::WriteGaussMatrixToFile(logFile, corr);
	logFile << "\nDeterminant of the correlation matrix:    " << corr.determinant() << std::endl << std::endl;


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
    Utils::WriteSymmGaussMatrixToFile(logFile, F);
	logFile << "\nEigenvalues of the Quaternion matrix:\n";
    Utils::WriteVectorToFile(logFile, Feig.eigenvalues(), digits, clean, threshold);
	logFile << "\nEigenvectors of the Quaternion matrix:\n";
    Utils::WriteGaussMatrixToFile(logFile, Feig.eigenvectors());
	logFile << "\nThe best-rotation Eigen-Quaternion:\n";
    Utils::WriteVectorToFile(logFile, lQuart, digits, clean, threshold);
	logFile << "\nThe optimal rotation matrix:\n";
    Utils::WriteGaussMatrixToFile(logFile, rotmat);
	logFile << "The big rotation matrix:\n";
    Utils::WriteGaussMatrixToFile(logFile, BigRotMat);

	// rotate the structures
	Eigen::VectorXd XRg(BigRotMat * XSg);
	Eigen::VectorXd QRg(BigRotMat * QSg);

	double RMSD_after = 0;
	for (int i = 0; i < Ncoords; i++)
		RMSD_after += std::abs(double(XRg(i) - XSe(i))) * std::abs(double(XRg(i) - XSe(i)));
	RMSD_after /= Natoms;

	logFile << "\nnew RMSD after Quaternion rotation:    " << RMSD_after << std::endl;

	Eigen::MatrixXd J(Eigen::MatrixXd::Zero(Nmodes, Nmodes));
	Eigen::VectorXd K(Eigen::VectorXd::Zero(Nmodes));

	// check if the rotation actually improved the alignment
	// and calculate the Duschinsky matrix and the displacement vector:
	if (RMSD_after < RMSD_before)
	{
		logFile << "The rotation has reduced the RMSD, so we will use the rotated structures\n";

		Eigen::MatrixXd NMRg(BigRotMat * NMOg);

		logFile << "\nthe rotated ground state mass-weighted normal modes:\n";
        Utils::WriteGaussMatrixToFile(logFile, NMRg);
		logFile << "and their metric:\n";
        Utils::WriteSymmGaussMatrixToFile(logFile, NMRg.transpose() * NMRg);

		J = NMOe.transpose() * NMRg;
		K = NMOe.transpose() * (QRg - QSe);
	}
	else
	{
		logFile << "The rotation did not reduce the RMSD, so we will use the original structures\n";

		J = NMOe.transpose() * NMOg;
		K = NMOe.transpose() * (QSg - QSe);
	}

	logFile << "\nDifference in equilibrium coordinates:\n";
    Utils::WriteVectorToFile(logFile, QSg - QSe, digits, clean, threshold);

	logFile << "\nDuschinsky Matrix J:\n";
    Utils::WriteGaussMatrixToFile(logFile, J);
	logFile << "\nDisplacement Vector K:\n";
    Utils::WriteVectorToFile(logFile, K, digits, clean, threshold);

	/*
	 * Sort normal modes by strongest displacement:
	 */
	std::vector<int> KsortNM(Nmodes);
	Eigen::VectorXd Ktmp = K;
	for (int i = 0; i < Nmodes; i++)
	{
		int Imax = 0;
		double Kmax = 0.0;
		for (int j = 0; j < Nmodes; j++)
		{
			if (std::abs(double(Ktmp(j))) > Kmax)
			{
				Imax = j;
				Kmax = std::abs(double(Ktmp(j)));
			}
		}
		Ktmp(Imax) = 0.0;
		KsortNM.at(i) = Imax;
	}
	logFile << "\nNormal modes sorted by displacement:\n";
	for (int i = 0; i < Nmodes; i++)
		logFile << std::setw(5) << KsortNM.at(i) + 1 << std::setw(20) << double(K(KsortNM.at(i))) << std::endl;

	/*
	 * sort normal modes by strongest rotation:
	 */
	std::vector<int> JsortNM(Nmodes);
	std::vector<double> offDiagRatios(Nmodes);
	for (int i = 0; i < Nmodes; i++)
	{
		double diagElement = double(J(i,i)) * double(J(i,i));
		double offDiags = 0.0;
		for (int j = 0; j < Nmodes; j++)
			if (j != i)
				offDiags += double(J(i,j)) * double(J(i,j));
		offDiagRatios.at(i) = offDiags / diagElement;
	}
	std::vector<double> oDRtmp = offDiagRatios;
	for (int i = 0; i < Nmodes; i++)
	{
		int Imax = 0;
		double Rmax = 0.0;
		for (int j = 0; j < Nmodes; j++)
		{
			if (oDRtmp.at(j) > Rmax)
			{
				Imax = j;
				Rmax = oDRtmp.at(j);
			}
		}
		oDRtmp.at(Imax) = 0.0;
		JsortNM.at(i) = Imax;
	}
	logFile << "\nNormal modes sorted by rotation:\n";
		for (int i = 0; i < Nmodes; i++)
			logFile << std::setw(5) << JsortNM.at(i) + 1 << std::setw(20) << offDiagRatios.at(JsortNM.at(i)) << std::endl;


	logFile << "\nMetric of the Duschinsky matrix:\n";
    Utils::WriteSymmGaussMatrixToFile(logFile, J.transpose() * J);
	logFile << "\nDeterminant of the Duschinsky matrix:     " << J.determinant() << std::endl;

	if (noJ)
	{
		logFile << "\nWe are ignoring the Duschinsky rotation in this calculation, setting J = 1\n";
		J = Eigen::MatrixXd::Identity(Nmodes, Nmodes);
	}
	if (noK)
	{
		logFile << "\nWe are ignoring the displacement in this calculation, setting k = 0\n";
		K = Eigen::VectorXd::Zero(Nmodes);
	}

	/*
	 * Write MCTDH files:
	 */
	Eigen::VectorXd f1 = GState.intFC();
	Eigen::VectorXd f2 = EState.intFC();

	// calculate the new force constants:
	Eigen::VectorXd fp = Eigen::VectorXd::Zero(Nmodes);
	for (int m = 0; m < Nmodes; m++)
		for (int n = 0; n < Nmodes; n++)
			fp(m) += f2(n) * J(n,m) * J(n,m);

	// calculate the first-order coefficients:
	Eigen::VectorXd kappa = Eigen::VectorXd::Zero(Nmodes);
	for (int m = 0; m < Nmodes; m++)
		for (int n = 0; n < Nmodes; n++)
			kappa(m) += f2(n) * K(n) * J(n,m);

	// calculate the couplings:
	Eigen::MatrixXd phi(Eigen::MatrixXd::Zero(Nmodes,Nmodes));
	Eigen::MatrixXd phiFull(Eigen::MatrixXd::Zero(Nmodes,Nmodes));
	for (int m = 0; m < Nmodes; m++)
	{
		for (int o = m + 1; o < Nmodes; o++)
		{
			for (int n = 0; n < Nmodes; n++)
				phi(m,o) += f2(n) * J(n,m) * J(n,o);
			phi(o,m) = phi(m,o);
			phiFull(m,o) = phi(m,o);
			phiFull(o,m) = phi(m,o);
		}
		phiFull(m,m) = fp(m);
	}

	// calculate the energy shifts:
	Eigen::VectorXd d(Nmodes);
	for (int i = 0; i < Nmodes; i++)
		d(i) = 0.5 * f2(i) * K(i) * K(i);

	logFile << "\nExcited state effective force constant matrix:\n";
    Utils::WriteSymmGaussMatrixToFile(logFile, phiFull);

	logFile << "\nShift coefficients:\n";
    Utils::WriteVectorToFile(logFile, kappa, digits, clean, threshold);

	logFile << "\nEnergy offsets:\n";
    Utils::WriteVectorToFile(logFile, d, digits, clean, threshold);

	/*
	 * convert everything to frequency-weighted units:
	 */
	Eigen::VectorXd fp_fw(Nmodes);
	for (int i = 0; i < Nmodes; i++)
		fp_fw(i) = double(fp(i)) / sqrt(double(f1(i)));

	Eigen::MatrixXd phiFull_fw(Nmodes, Nmodes);
	for (int i = 0; i < Nmodes; i++)
		for (int j = 0; j < Nmodes; j++)
			phiFull_fw(i,j) = double(phiFull(i,j)) / (pow(double(f1(i)), 0.25) * pow(double(f1(j)), 0.25));

	Eigen::VectorXd kappa_fw(Nmodes);
	for (int i = 0; i < Nmodes; i++)
		kappa_fw(i) = double(kappa(i)) / pow(double(f1(i)), 0.25);

	logFile << "\nFrequency-weighted effective force constants:\n";
    Utils::WriteVectorToFile(logFile, fp_fw, digits, clean, threshold);

	logFile << "\nFrequency-weighted force constant matrix:\n";
    Utils::WriteSymmGaussMatrixToFile(logFile, phiFull_fw);

	logFile << "\nFrequency-weighted shift coefficients:\n";
    Utils::WriteVectorToFile(logFile, kappa_fw, digits, clean, threshold);

	/*
	 * Sort normal modes by strongest frequency-weighted displacement:
	 */
	std::vector<int> kappa_fw_sortNM(Nmodes);
	Eigen::VectorXd kappatmp = kappa_fw;
	for (int i = 0; i < Nmodes; i++)
	{
		int Imax = 0;
		double Kmax = 0.0;
		for (int j = 0; j < Nmodes; j++)
		{
			if (std::abs(double(kappatmp(j))) > Kmax)
			{
				Imax = j;
				Kmax = std::abs(double(kappatmp(j)));
			}
		}
		kappatmp(Imax) = 0.0;
		kappa_fw_sortNM.at(i) = Imax;
	}
	logFile << "\nNormal modes sorted by frequency-weighted displacement:\n";
	for (int i = 0; i < Nmodes; i++)
		logFile << std::setw(5) << kappa_fw_sortNM.at(i) + 1 << std::setw(20) << double(kappa_fw(kappa_fw_sortNM.at(i))) << std::endl;

	/*
	 * Sort normal modes by strongest frequency-weighted coupling:
	 */
	std::vector<int> phi_fw_sortNM(Nmodes);
	std::vector<double> max_phi_fw(Nmodes);
	for (int i = 0; i < Nmodes; i++)
	{
		double maxPhi = 0.0;
		for (int j = 0; j < Nmodes; j++)
			if (i != j && std::abs(double(phiFull_fw(i,j))) > std::abs(maxPhi))
				maxPhi = double(phiFull_fw(i,j));
		max_phi_fw.at(i) = maxPhi;
	}
	std::vector<double> phiSortTmp = max_phi_fw;
	for (int i = 0; i < Nmodes; i++)
	{
		int Imax = 0;
		double Kmax = 0.0;
		for (int j = 0; j < Nmodes; j++)
		{
			if (std::abs(phiSortTmp.at(j)) > Kmax)
			{
				Imax = j;
				Kmax = std::abs(phiSortTmp.at(j));
			}
		}
		phiSortTmp.at(Imax) = 0.0;
		phi_fw_sortNM.at(i) = Imax;
	}
	logFile << "\nNormal modes sorted by frequency-weighted coupling:\n";
	for (int i = 0; i < Nmodes; i++)
		logFile << std::setw(5) << phi_fw_sortNM.at(i) + 1 << std::setw(20) << max_phi_fw.at(phi_fw_sortNM.at(i)) << std::endl;




	/*
	 * Determine, based on the values of phi and kappa, which modes
	 * should be included in the spectrum calculation:
	 */
	std::vector<bool> isPresent(Nmodes);
	std::vector<int> presentModes;
	int Npresent = 0;
	for (int i = 0; i < Nmodes; i++)
	{
		// find the maximum phi:
		double maxPhi = 0.0;
		for (int j = 0; j < Nmodes; j++)
			if (i != j && double(phiFull(i,j)) > maxPhi)
				maxPhi = std::abs(double(phiFull(i,j)));

		if (maxPhi > paraThres || std::abs(double(kappa(i))) > paraThres)
		{
			isPresent.at(i) = true;
			Npresent++;
			presentModes.push_back(i);
		}
		else
			isPresent.at(i) = false;
	}

	logFile << "\nModes included into the calculation:\n";
	for (int i = 0; i < Nmodes; i++)
	{
		logFile << "  mode" << std::setw(4) << i + 1;
		if (isPresent.at(i))
			logFile << "   YES\n";
		else
			logFile << "   NO\n";
	}
	logFile << std::setw(4) << Npresent << " out of" << std::setw(4) << Nmodes << " included.\n";


	/*
	 * Use the modes actually present in the calculation to calculate the energy at the minimum
	 * and energy offsets:
	 */
	// zero-point energies:
	double zpeGinc = 0.0;
	double zpeEinc = 0.0;
	for (int i = 0; i < Npresent; i++)
	{
		zpeGinc += 0.5 * sqrt(double(f1(presentModes.at(i))));
		zpeEinc += 0.5 * sqrt(double(f2(presentModes.at(i))));
	}
	// minimum coordinates:
	Eigen::ColPivHouseholderQR<Eigen::MatrixXd> phiLin(phiFull);
	Eigen::VectorXd minima = phiLin.solve(-kappa);
	// energy at the minimum
	double Emin = 0.0;
	for (int i = 0; i < Npresent; i++)
	{
		for (int j = 0; j < Npresent; j++)
			Emin += 0.5 * double(phiFull(presentModes.at(i),presentModes.at(j)))
						* double(minima(presentModes.at(i)))
						* double(minima(presentModes.at(j)));							// quadratic term
		Emin += double(kappa(presentModes.at(i))) * double(minima(presentModes.at(i)));	// linear term
		Emin += double(d(presentModes.at(i)));											// offset term
	}


	logFile << "\nVertical energies:\n";
	logFile << "                                              Eh               eV\n";
	logFile << "adiabatic excitation energy:           "
			<< std::setw(18) << std::scientific << std::setprecision(9) << deltaE
			<< std::setw(10) << std::fixed << std::setprecision(4) << deltaE * Eh2eV << std::endl;
	logFile << "ground state zero-point energy:        "
				<< std::setw(18) << std::scientific << std::setprecision(9) << zpeGinc
				<< std::setw(10) << std::fixed << std::setprecision(4) << zpeGinc * Eh2eV << std::endl;
	logFile << "excited state zero-point energy:       "
					<< std::setw(18) << std::scientific << std::setprecision(9) << zpeEinc
					<< std::setw(10) << std::fixed << std::setprecision(4) << zpeEinc * Eh2eV << std::endl;
	logFile << "difference between zero-point energies:"
					<< std::setw(18) << std::scientific << std::setprecision(9) << zpeEinc - zpeGinc
					<< std::setw(10) << std::fixed << std::setprecision(4) << (zpeEinc - zpeGinc) * Eh2eV << std::endl;
	logFile << "energy at the minimum (should be zero):"
						<< std::setw(18) << std::scientific << std::setprecision(9) << Emin
						<< std::setw(10) << std::fixed << std::setprecision(4) << Emin * Eh2eV << std::endl;
	logFile << "shifted minimum (E_ad - zpeG):         "
							<< std::setw(18) << std::scientific << std::setprecision(9) << deltaE - zpeGinc + Emin
							<< std::setw(10) << std::fixed << std::setprecision(4) << (deltaE - zpeGinc + Emin) * Eh2eV << std::endl;


	/*
	 * Start writing the MCTDH input files:
	 */
	std::string inputFileName = basename + ".inp";
	std::string operatorFileName = basename + ".op";

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
	inputFile << "mlbasis-section\n";
	for (int i = 0; i < Npresent - 1; i+= 2)
	{
		if (Npresent - i == 3)
			inputFile << "    [q_" << std::setfill('0') << std::setw(3) << presentModes.at(i) + 1
					  << " q_" << std::setfill('0') << std::setw(3) << presentModes.at(i+1) + 1
					  << " q_" << std::setfill('0') << std::setw(3) << presentModes.at(i+2) + 1 << "]\n";
		else
			inputFile << "    [q_" << std::setfill('0') << std::setw(3) << presentModes.at(i) + 1
					  << " q_" << std::setfill('0') << std::setw(3) << presentModes.at(i+1) + 1 << "]\n";
	}
	inputFile << "end-mlbasis-section\n\n";

	/*
	 * The pbasis-section
	 */
	inputFile << "pbasis-section\n";
	for (int i = 0; i < Nmodes; i++)
		if (isPresent.at(i))
		{
			// determine the number of SPFs:
			int numSPFs = lrint(-0.7 * log(double(fp(i)))) + 6;

			// determine the lower and upper bounds of the grid:
			// the basis boundarie are (originally) -kappa / fp +- 6.5 / fp**1/4
			double lowBound = -double(kappa(i) / fp(i)) - (std::abs(double(kappa(i) / fp(i))) + 7.5 / (sqrt(2.0) * pow(double(fp(i)), 0.25)));
			double uppBound = -double(kappa(i) / fp(i)) + (std::abs(double(kappa(i) / fp(i))) + 7.5 / (sqrt(2.0) * pow(double(fp(i)), 0.25)));

			inputFile << "    q_" << std::setfill('0') << std::setw(3) << i + 1
					<< "  ho  " << std::setw(3) << std::setfill(' ') << numSPFs << "  xi-xf  "
					<< std::fixed << std::setfill(' ') << std::setw(8)
			<< lowBound
			<< std::fixed << std::setfill(' ') << std::setw(8)
			<< uppBound << std::endl;
		}
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
		if (isPresent.at(i))
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
		if (isPresent.at(i))
			operatorFile << "    mass_q_" << std::setfill('0') << std::setw(3) << i + 1
						 << "  =  1.0\n";
	// the ground state force constants
	for (int i = 0; i < Nmodes; i++)
		if (isPresent.at(i))
		{
			operatorFile << "    f1_" << std::setfill('0') << std::setw(3) << i + 1 << "      = ";
            Utils::WriteFortranNumber(operatorFile, double(f1(i)));
			operatorFile << std::endl;
		}
	// the new effective excited state force constants
	for (int i = 0; i < Nmodes; i++)
		if (isPresent.at(i))
		{
			operatorFile << "    fp_" << std::setfill('0') << std::setw(3) << i + 1 << "      = ";
            Utils::WriteFortranNumber(operatorFile, double(fp(i)));
			operatorFile << std::endl;
		}
	// the couplings
	for (int i = 0; i < Nmodes; i++)
		for (int j = i + 1; j < Nmodes; j++)
			if (isPresent.at(i) && isPresent.at(j))
			{
				operatorFile << "    phi_" << std::setfill('0') << std::setw(3) << i + 1
						<< "_" << std::setfill('0') << std::setw(3) << j + 1 << " = ";
                Utils::WriteFortranNumber(operatorFile, double(phi(i,j)));
				operatorFile << std::endl;
			}
	// the first-order coefficients (shifts)
	for (int i = 0; i < Nmodes; i++)
		if (isPresent.at(i))
		{
			operatorFile << "    kappa_" << std::setfill('0') << std::setw(3) << i + 1 << "   = ";
            Utils::WriteFortranNumber(operatorFile, double(kappa(i)));
			operatorFile << std::endl;
		}
	// the energy offsets
	for (int i = 0; i < Nmodes; i++)
		if (isPresent.at(i))
		{
			operatorFile << "    d_" << std::setfill('0') << std::setw(3) << i + 1 << "       = ";
            Utils::WriteFortranNumber(operatorFile, double(d(i)));
			operatorFile << std::endl;
		}
	// the electronic offset minus the ground state ZPE
	double zpe1 = 0.0;
	for (int i = 0; i < Nmodes; i++)
		if (isPresent.at(i))
			zpe1 += 0.5 * sqrt(double(f1(i)));
	operatorFile << "    dE          = ";
    Utils::WriteFortranNumber(operatorFile, deltaE - zpeGinc + Emin);
	operatorFile << "\nend-parameter-section\n\n";

	/*
	 * The hamiltonian section
	 */
	operatorFile << "hamiltonian-section";
	for (int i = 0; i < Npresent; i++)
	{
		if (i % 8 == 0)
			operatorFile << std::endl << "modes";
		operatorFile << " | q_" << std::setfill('0') << std::setw(3) << presentModes.at(i) + 1;
	}
	operatorFile << std::endl;
	for (int i = 0; i < Npresent; i++)
		operatorFile << "1.0         |" << i + 1 << " KE\n";
	for (int i = 0; i < Npresent; i++)
			operatorFile << "0.5*fp_" << std::setfill('0') << std::setw(3) << presentModes.at(i) + 1
					 << "  |" << i + 1 << " q^2\n";
	for (int i = 0; i < Npresent; i++)
		for (int j = i + 1; j < Npresent; j++)
				operatorFile << "phi_" << std::setfill('0') << std::setw(3) << presentModes.at(i) + 1
							 << "_" << std::setfill('0') << std::setw(3) << presentModes.at(j) + 1
							 << " |" << i + 1 << " q"
							 << " |" << j + 1 << " q\n";
	for (int i = 0; i < Npresent; i++)
			operatorFile << "kappa_" << std::setfill('0') << std::setw(3) << presentModes.at(i) + 1
						 << "   |" << i + 1 << " q\n";
	for (int i = 0; i < Npresent; i++)
			operatorFile << "d_" << std::setfill('0') << std::setw(3) << presentModes.at(i) + 1
						 << "       |" << i + 1 << " 1\n";
	operatorFile << "dE          |1 1\n";
	operatorFile << "end-hamiltonian-section\n\n";

	/*
	 * One-dimensional Hamiltonians for the ground state normal modes
	 */
	for (int i = 0; i < Nmodes; i++)
		if (isPresent.at(i))
		{
			operatorFile << "hamiltonian-section_Eq_" << std::setfill('0') << std::setw(3) << i + 1 << std::endl;
			operatorFile << "usediag\n";
			operatorFile << "modes      | q_" << std::setfill('0') << std::setw(3) << i + 1 << std::endl;
			operatorFile << "1.0        |1 KE\n";
			operatorFile << "0.5*f1_" << std::setfill('0') << std::setw(3) << i + 1 << " |1 q^2\n";
			operatorFile << "end-hamiltonian-section\n\n";
		}
	operatorFile << "end-operator\n";

	return 0;
}
