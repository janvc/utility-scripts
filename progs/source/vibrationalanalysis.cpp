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


#include <iomanip>
#include <eigen3/Eigen/Eigenvalues>
#include <armadillo>
#include "vibrationalanalysis.h"
#include "constants.h"
#include "GaussFchk.h"
#include "utilities.h"


VibrationalAnalysis::VibrationalAnalysis(GaussFchk &initFchk)
	: m_fchk(&initFchk)
{
	constructor();
}

VibrationalAnalysis::VibrationalAnalysis(std::ifstream &initStream)
	: m_fchk(new GaussFchk(initStream))
{
	constructor();
}

/**
 * Perform vibrational analysis according to the scheme outlined in
 * the article "Vibrational Analysis in Gaussian" by J. Ochterski
 * The calculation steps are referenced by the corresponding equation numbers.
 */
void VibrationalAnalysis::constructor()
{
	/*
	 * Read the data from the checkpoint file:
	 */
	m_Natoms = m_fchk->ReadInteger("Number of atoms");
	m_Ncoords = 3 * m_Natoms;
	m_Nmodes = m_Ncoords - 6;
	X_min = m_fchk->ReadVector("Current cartesian coordinates");
	atNums = m_fchk->ReadIVector("Atomic numbers");
	m_masses = m_fchk->ReadVector("Real atomic weights") * amu2au;
	Fcart_min = m_fchk->ReadSymmetricMatrix("Cartesian Force Constants");
	GaussModes = m_fchk->ReadMatrix("Vib-Modes", m_Nmodes, m_Ncoords);

	MassVec = Eigen::VectorXd::Zero(m_Ncoords);
	for (int i = 0; i < m_Natoms; i++)
	{
		MassVec(3 * i + 0) = m_masses(i);
		MassVec(3 * i + 1) = m_masses(i);
		MassVec(3 * i + 2) = m_masses(i);
	}

	MassInvMat = Eigen::MatrixXd::Zero(m_Ncoords, m_Ncoords);
	for (int i = 0; i < m_Ncoords; i++)
		MassInvMat(i,i) = 1.0 / sqrt(double(MassVec(i)));

	Q_min = MassInvMat.inverse() * X_min;

	/*
	 * Calculate the mass-weighted Hessian      **** eq. (2) ****
	 */
	Fmwc_min = MassInvMat * Fcart_min * MassInvMat;


	/*
	 * Calculate additional properties:
	 */
	calcCom();
	calcInert();


	/*
	 * Calculate the D-matrix       **** eq. (5) ****
	 */
	Dmat = Eigen::MatrixXd::Random(m_Ncoords, m_Ncoords);

	for (int i = 0; i < 3; i++)
		for (int j = 0; j < m_Natoms; j++)
		{
			Dmat(3 * j + 0, i) = 0.0;
			Dmat(3 * j + 1, i) = 0.0;
			Dmat(3 * j + 2, i) = 0.0;
			Dmat(3 * j + i, i) = sqrt(double(m_masses(j)));
		}

	for (int i = 0; i < m_Natoms; i++)
	{
		Dmat(3 * i + 0, 3) = (X_sr(3 * i + 1) * m_prinAxes.transpose()(2,0) - X_sr(3 * i + 2) * m_prinAxes.transpose()(1,0)) * sqrt(double(m_masses(i)));
		Dmat(3 * i + 1, 3) = (X_sr(3 * i + 1) * m_prinAxes.transpose()(2,1) - X_sr(3 * i + 2) * m_prinAxes.transpose()(1,1)) * sqrt(double(m_masses(i)));
		Dmat(3 * i + 2, 3) = (X_sr(3 * i + 1) * m_prinAxes.transpose()(2,2) - X_sr(3 * i + 2) * m_prinAxes.transpose()(1,2)) * sqrt(double(m_masses(i)));
		Dmat(3 * i + 0, 4) = (X_sr(3 * i + 2) * m_prinAxes.transpose()(0,0) - X_sr(3 * i + 0) * m_prinAxes.transpose()(2,0)) * sqrt(double(m_masses(i)));
		Dmat(3 * i + 1, 4) = (X_sr(3 * i + 2) * m_prinAxes.transpose()(0,1) - X_sr(3 * i + 0) * m_prinAxes.transpose()(2,1)) * sqrt(double(m_masses(i)));
		Dmat(3 * i + 2, 4) = (X_sr(3 * i + 2) * m_prinAxes.transpose()(0,2) - X_sr(3 * i + 0) * m_prinAxes.transpose()(2,2)) * sqrt(double(m_masses(i)));
		Dmat(3 * i + 0, 5) = (X_sr(3 * i + 0) * m_prinAxes.transpose()(1,0) - X_sr(3 * i + 1) * m_prinAxes.transpose()(0,0)) * sqrt(double(m_masses(i)));
		Dmat(3 * i + 1, 5) = (X_sr(3 * i + 0) * m_prinAxes.transpose()(1,1) - X_sr(3 * i + 1) * m_prinAxes.transpose()(0,1)) * sqrt(double(m_masses(i)));
		Dmat(3 * i + 2, 5) = (X_sr(3 * i + 0) * m_prinAxes.transpose()(1,2) - X_sr(3 * i + 1) * m_prinAxes.transpose()(0,2)) * sqrt(double(m_masses(i)));
	}

	for (int i = 0; i < m_Ncoords; i++)
		Dmat.col(i) /= Eigen::VectorXd(Dmat.col(i)).norm();

	// Gram-Schmidt-orthogonalization:
	Eigen::MatrixXd tmpMat = Eigen::MatrixXd::Zero(m_Ncoords, m_Ncoords);
	tmpMat.block(0, 0, m_Ncoords, 6) = Dmat.block(0, 0, m_Ncoords, 6);
	for (int i = 6; i < m_Ncoords; i++)
	{
		Eigen::VectorXd tmpVec = Dmat.col(i);
		for (int j = 0; j < i; j++)
		{
			double factor = tmpVec.dot(tmpMat.col(j));
			tmpVec -= factor * tmpMat.col(j);
		}
		tmpVec.normalize();
		tmpMat.col(i) = tmpVec;
	}
	Dmat = tmpMat;


	/*
	 * Calculate the Hessian in terms of 'internal' coordinates    **** eq. (6) ****
	 * It has dimension 3N-6 x 3N-6
	 */
	Fint_min = Eigen::MatrixXd(Dmat.transpose() * Fmwc_min * Dmat).block(6, 6, m_Nmodes, m_Nmodes);


	/*
	 * Diagonalize the Hessians (full mass-weighted and internal) to obtain the vibrational
	 * frequencies and normal modes
	 */
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigSolv;

	eigSolv.compute(Fmwc_min);
	MwcFrcCon = eigSolv.eigenvalues();

	eigSolv.compute(Fint_min);
	IntFrcCon = eigSolv.eigenvalues();
	Lint_min = eigSolv.eigenvectors();

	// calculate mass-weighted and cartesian displacements:
	Lmwc_min = Dmat.block(0, 6, m_Ncoords, m_Nmodes) * Lint_min;
	Lcrt_min = MassInvMat * Lmwc_min;

	// normalize cartesian displacements (and get their reduced masses along the way):
	RedMasses = Eigen::VectorXd::Zero(m_Nmodes);
	for (int i = 0; i < m_Nmodes; i++)
	{
		double normFac = 1.0 / Eigen::VectorXd(Lcrt_min.col(i)).norm();
		RedMasses(i) = normFac * normFac;
		Lcrt_min.col(i) *= normFac;
	}


	/*
	 * Calculate the frequencies in [cm-1] from the force constants:
	 */
	MwcFreqs = Eigen::VectorXd::Zero(m_Ncoords);
	IntFreqs = Eigen::VectorXd::Zero(m_Nmodes);
	AUFreqs = Eigen::VectorXd::Zero(m_Nmodes);
	for (int i = 0; i < m_Ncoords; i++)
	{
		if (MwcFrcCon(i) >= 0)
			MwcFreqs(i) =  sqrt( double(MwcFrcCon(i)) * Eh / (a0 * a0 * me)) / (2.0 * M_PI * c0 * 100.0);
		else
			MwcFreqs(i) = -sqrt(-double(MwcFrcCon(i)) * Eh / (a0 * a0 * me)) / (2.0 * M_PI * c0 * 100.0);
	}
	for (int i = 0; i < m_Nmodes; i++)
	{
		if (IntFrcCon(i) >= 0)
		{
			AUFreqs(i) = sqrt(double(IntFrcCon(i)));
			IntFreqs(i) =  sqrt( double(IntFrcCon(i)) * Eh / (a0 * a0 * me)) / (2.0 * M_PI * c0 * 100.0);
		}
		else
		{
			AUFreqs(i) = sqrt(-double(IntFrcCon(i)));
			IntFreqs(i) = -sqrt(-double(IntFrcCon(i)) * Eh / (a0 * a0 * me)) / (2.0 * M_PI * c0 * 100.0);
		}
	}
}


VibrationalAnalysis::~VibrationalAnalysis()
{
}

void VibrationalAnalysis::createThirdDerivs(std::string &baseName)
{
	/*
	 * read in the displaced Hessians:
	 */
	for (int i = 0; i < m_Nmodes; i++)
	{
		// start with the positive displacement:
		std::string dirName_p = baseName + "_" + std::to_string(i) + "p";
		std::ifstream dispFchkfile_p(dirName_p + "/" + dirName_p + ".fchk", std::ifstream::in);
		GaussFchk dispFchk_p(dispFchkfile_p);
		Eigen::MatrixXd tmpHess = dispFchk_p.ReadSymmetricMatrix("Cartesian Force Constants");
		dispFchkfile_p.close();

		Fcart_Dp.push_back(tmpHess);
		tmpHess = MassInvMat * tmpHess * MassInvMat;
		Fmwc_Dp.push_back(tmpHess);
		tmpHess = Lmwc_min.transpose() * tmpHess * Lmwc_min;
		Fdiag_Dp.push_back(tmpHess);


		// now do the negative displacement:
		std::string dirName_n = baseName + "_" + std::to_string(i) + "n";
		std::ifstream dispFchkfile_n(dirName_n + "/" + dirName_n + ".fchk", std::ifstream::in);
		GaussFchk dispFchk_n(dispFchkfile_n);
		tmpHess = dispFchk_n.ReadSymmetricMatrix("Cartesian Force Constants");
		dispFchkfile_n.close();

		Fcart_Dn.push_back(tmpHess);
		tmpHess = MassInvMat * tmpHess * MassInvMat;
		Fmwc_Dn.push_back(tmpHess);
		tmpHess = Lmwc_min.transpose() * tmpHess * Lmwc_min;
		Fdiag_Dn.push_back(tmpHess);
	}

	arma::Cube<double> phiTmp(m_Nmodes, m_Nmodes, m_Nmodes);

	for (int i = 0; i < m_Nmodes; i++)
		for (int j = 0; j < m_Nmodes; j++)
			for (int k = 0; k < m_Nmodes; k++)
				phiTmp(i,j,k) = (Fdiag_Dp.at(i)(j,k) - Fdiag_Dn.at(i)(j,k)) / (2.0 * shiftFac);

	rawDerivs = phiTmp;

	for (int i = 0; i < m_Nmodes; i++)
		for (int j = 0; j < m_Nmodes; j++)
			for (int k = 0; k < m_Nmodes; k++)
				phiTmp(i,j,k) = (rawDerivs(i,j,k) + rawDerivs(j,k,i) + rawDerivs(k,i,j)) / 3.0;

	avgDerivs = phiTmp;
}

void VibrationalAnalysis::createAnharmMCTDHoper(const std::string &baseName, const double thres)
{
	std::ofstream operFile(baseName + ".op");
	operFile.precision(8);

	/*
	 * The op-define section
	 */
	operFile << "op_define-section\n";
	operFile << "    title\n";
	operFile << "        " << baseName << std::endl;
	operFile << "    end-title\n";
	operFile << "end-op_define-section\n\n";

	/*
	 * The parameter section
	 */
	operFile << "parameter-section\n";
	// write the masses:
	for (int i = 0; i < m_Nmodes; i++)
		operFile << "    mass_q_" << std::setfill('0') << std::setw(3) << i + 1 << "      =  1.0\n";
	// write the diagonal force constants:
	for (int i = 0; i < m_Nmodes; i++)
	{
		operFile << "    f_" << std::setfill('0') << std::setw(3) << i + 1 << "           = ";
		WriteFortranNumber(operFile, double(IntFrcCon(i)));
		operFile << std::endl;
	}
	// write the third-order couplings:
	for (int i = 0; i < m_Nmodes; i++)
		for (int j = i + 1; j < m_Nmodes; j++)
			for (int k = j + 1; k < m_Nmodes; k++)
				if (std::abs(avgDerivs(i,j,k)) >= thres)
				{
					operFile << "    phi_" << std::setfill('0') << std::setw(3) << i + 1
									<< "_" << std::setfill('0') << std::setw(3) << j + 1
									<< "_" << std::setfill('0') << std::setw(3) << k + 1
									<< " = ";
					WriteFortranNumber(operFile, avgDerivs(i,j,k));
					operFile << std::endl;
				}

	operFile << "end-parameter-section\n\n";

	/*
	 * The hamiltonian section
	 */
	operFile << "hamiltonian-section";
	for (int i = 0; i < m_Nmodes; i++)
	{
		if (i % 8 == 0)
			operFile << std::endl << "modes";
		operFile << " | q_" << std::setfill('0') << std::setw(3) << i + 1;
	}
	operFile << std::endl;

	for (int i = 0; i < m_Nmodes; i++)
		operFile << "1.0             |" << i + 1 << " KE\n";

	for (int i = 0; i < m_Nmodes; i++)
		operFile << "0.5*f_" << std::setfill('0') << std::setw(3) << i + 1 << "       |" << i + 1 << " q^2\n";

	for (int i = 0; i < m_Nmodes; i++)
		for (int j = i + 1; j < m_Nmodes; j++)
			for (int k = j + 1; k < m_Nmodes; k++)
				if (std::abs(avgDerivs(i,j,k)) >= thres)
					operFile << "phi_" << std::setfill('0') << std::setw(3) << i + 1
								<< "_" << std::setfill('0') << std::setw(3) << j + 1
								<< "_" << std::setfill('0') << std::setw(3) << k + 1
								<< " |" << i + 1 << " q"
								<< " |" << j + 1 << " q"
								<< " |" << k + 1 << " q\n";

	operFile << "end-hamiltonian-section\n\n";

	/*
	 * Additional Hamiltonians for the initial condition:
	 */
	for (int i = 0; i < m_Nmodes; i++)
	{
		operFile << "hamiltonian-section_Eq_" << std::setfill('0') << std::setw(3) << i + 1 << std::endl;
		operFile << "usediag\n";
		operFile << "modes | q_" << std::setfill('0') << std::setw(3) << i + 1 << std::endl;
		operFile << "1.0   |1 KE\n";
		operFile << "0.5*f_" << std::setfill('0') << std::setw(3) << i + 1 << " |1 q^2\n";
		operFile << "end-hamiltonian-section\n\n";
	}

	operFile << "end-operator\n";
}

void VibrationalAnalysis::calcKappa(GaussFchk &esFchk)
{
	Eigen::VectorXd esGrad_x = esFchk.ReadVector("Cartesian Gradient");

	Eigen::VectorXd esGrad_q = MassInvMat.inverse() * esGrad_x;

	Eigen::VectorXd esGrad_nm = Lmwc_min.transpose() * esGrad_q;

	WriteVector(esGrad_nm, 6, false);
}

void VibrationalAnalysis::prtMinGeo()
{
	std::cout << "                 Molecular geometry at the minimum:\n";
	std::cout << " ---------------------------------------------------------------------\n";
	std::cout << " Center     Atomic      Atomic             Coordinates (Angstroms)\n";
	std::cout << " Number     Number       Type             X           Y           Z\n";
	std::cout << " ---------------------------------------------------------------------\n";

	for (int i = 0; i < m_Natoms; i++)
		std::cout << std::fixed << std::setprecision(6)
				  << std::setw(7) << i + 1
				  << std::setw(11) << int(atNums(i))
				  << std::setw(12) << 0 << "    "
				  << std::setw(12) << double(X_min(3 * i + 0)) * ang2a0
				  << std::setw(12) << double(X_min(3 * i + 1)) * ang2a0
				  << std::setw(12) << double(X_min(3 * i + 2)) * ang2a0 << std::endl;
	std::cout << " ---------------------------------------------------------------------\n";
}

void VibrationalAnalysis::prtMinModes()
{
	std::cout << "         Cartesian normal modes in Gaussian HPModes style:\n";
	std::cout << " ---------------------------------------------------------------------\n";
	std::cout << " Low frequencies ---" << std::fixed << std::setprecision(4)
			  << std::setw(10) << double(MwcFreqs(0))
			  << std::setw(10) << double(MwcFreqs(1))
			  << std::setw(10) << double(MwcFreqs(2))
			  << std::setw(10) << double(MwcFreqs(3))
			  << std::setw(10) << double(MwcFreqs(4))
			  << std::setw(10) << double(MwcFreqs(5))
			  << std::endl << " Low frequencies ---"
			  << std::setw(10) << double(MwcFreqs(6))
			  << std::setw(10) << double(MwcFreqs(7))
			  << std::setw(10) << double(MwcFreqs(8))
			  << std::endl;

	int remainder = m_Nmodes % 5;
	for (int i = 0; i < m_Nmodes - remainder; i += 5)
	{
		std::cout << "                  "
				  << std::setw(10) << i + 1
				  << std::setw(10) << i + 2
				  << std::setw(10) << i + 3
				  << std::setw(10) << i + 4
				  << std::setw(10) << i + 5 << std::endl;
		std::cout << "       Frequencies --- " << std::fixed << std::setprecision(4)
				  << std::setw(10) << double(IntFreqs(i + 0))
				  << std::setw(10) << double(IntFreqs(i + 1))
				  << std::setw(10) << double(IntFreqs(i + 2))
				  << std::setw(10) << double(IntFreqs(i + 3))
				  << std::setw(10) << double(IntFreqs(i + 4)) << std::endl;
		std::cout << "    Reduced masses --- " << std::fixed << std::setprecision(4)
				  << std::setw(10) << double(RedMasses(i + 0))
				  << std::setw(10) << double(RedMasses(i + 1))
				  << std::setw(10) << double(RedMasses(i + 2))
				  << std::setw(10) << double(RedMasses(i + 3))
				  << std::setw(10) << double(RedMasses(i + 4)) << std::endl;
		std::cout << " Coord Atom Element:\n";

		for (int j = 0; j < m_Ncoords; j++)
			std::cout << std::fixed << std::setprecision(5)
			          << std::setw(4) << (j % 3) + 1
					  << std::setw(6) << (j / 3) + 1
					  << std::setw(6) << int(atNums(j / 3)) << "       "
					  << std::setw(10) << double(Lcrt_min(j, i + 0))
					  << std::setw(10) << double(Lcrt_min(j, i + 1))
					  << std::setw(10) << double(Lcrt_min(j, i + 2))
					  << std::setw(10) << double(Lcrt_min(j, i + 3))
					  << std::setw(10) << double(Lcrt_min(j, i + 4)) << std::endl;
	}

	std::cout << "                  ";
	for (int i = m_Nmodes - remainder; i < m_Nmodes; i++)
		std::cout << std::setw(10) << i + 1;
	std::cout << std::endl << "       Frequencies --- " << std::fixed << std::setprecision(4);
	for (int i = m_Nmodes - remainder; i < m_Nmodes; i++)
		std::cout << std::setw(10) << double(IntFreqs(i));
	std::cout << std::endl << "    Reduced masses --- " << std::fixed << std::setprecision(4);
	for (int i = m_Nmodes - remainder; i < m_Nmodes; i++)
		std::cout << std::setw(10) << double(RedMasses(i));
	std::cout << std::endl << " Coord Atom Element:\n";

	for (int j = 0; j < m_Ncoords; j++)
	{
		std::cout << std::fixed << std::setprecision(5)
				  << std::setw(4) << (j % 3) + 1
				  << std::setw(6) << (j / 3) + 1
				  << std::setw(6) << int(atNums(j / 3)) << "       ";

		for (int i = m_Nmodes - remainder; i < m_Nmodes; i++)
			std::cout << std::setw(10) << double(Lcrt_min(j, i));
		std::cout << std::endl;
	}
	std::cout << " ---------------------------------------------------------------------\n";
}

void VibrationalAnalysis::prtDiagHess()
{
	Eigen::MatrixXd diagHess = Lmwc_min.transpose() * Fmwc_min * Lmwc_min;

	WriteMatrix(diagHess, 12, false);
}

void VibrationalAnalysis::prtFreqs()
{
	std::cout << "Force constants and vibrational frequencies of the normal modes:\n";
	for (int i = 0; i < m_Nmodes; i++)
	{
		std::cout << std::setw(4) << i + 1
				  << std::scientific << std::setprecision(7)
				  << std::setw(15) << double(IntFrcCon(i))
				  << std::setw(15) << double(AUFreqs(i))
				  << std::fixed << std::setprecision(4)
				  << std::setw(12) << double(IntFreqs(i)) << std::endl;
	}
}

void VibrationalAnalysis::prtCubics()
{
	std::cout << "Displaced Hessians, in a different style:\n";

	int remainder = m_Nmodes % 5;
	for (int i = 0; i < m_Nmodes - remainder; i += 5)
	{
		std::cout << "  ";
		for (int j = 0; j < 5; j++)
			std::cout << std::setw(15) << i + j + 1;
		std::cout << std::endl;

		int counter = 0;
		for (int j = 0; j < m_Ncoords; j++)
			for (int k = j; k < m_Ncoords; k++)
			{
				counter++;
				std::cout << std::setw(5) << counter;
				for (int l = 0; l < 5; l++)
				{
					Eigen::MatrixXd tmpHess_p = Fcart_Dp.at(i + l);
					Eigen::MatrixXd tmpHess_n = Fcart_Dn.at(i + l);

					double value = (double(tmpHess_p(j, k)) - double(tmpHess_n(j, k))) / (2.0 * shiftFac);
					std::cout << std::scientific << std::setprecision(7) << std::setw(15) << value;
				}
				std::cout << std::endl;
			}
	}
	for (int i = m_Nmodes - remainder; i < m_Nmodes; i++)
		std::cout << std::setw(15) << i + 1;
	std::cout << std::endl;
	int counter = 0;
	for (int j = 0; j < m_Ncoords; j++)
		for (int k = j; k < m_Ncoords; k++)
		{
			counter++;
			std::cout << std::setw(5) << counter;
			for (int i = m_Nmodes - remainder; i < m_Nmodes; i++)
			{
				Eigen::MatrixXd tmpHess_p = Fcart_Dp.at(i);
				Eigen::MatrixXd tmpHess_n = Fcart_Dn.at(i);

				double value = (double(tmpHess_p(j, k)) - double(tmpHess_n(j, k))) / (2.0 * shiftFac);
				std::cout << std::scientific << std::setprecision(7) << std::setw(15) << value;
			}
			std::cout << std::endl;
		}


	std::cout << "Cubic force constants:\n";

	for (int i = 0; i < m_Nmodes; i++)
		for (int j = i; j < m_Nmodes; j++)
			for (int k = j; k < m_Nmodes; k++)
				std::cout << std::setw(4) << i + 1
						  << std::setw(4) << j + 1
						  << std::setw(4) << k + 1
						  << std::scientific << std::setprecision(7)
						  << std::setw(15) << double(rawDerivs(i,j,k))
						  << std::setw(15) << double(rawDerivs(j,k,i))
						  << std::setw(15) << double(rawDerivs(k,i,j))
						  << std::setw(20) << double(avgDerivs(i,j,k))
						  << std::setw(20) << double(avgDerivs(i,j,k)) / sqrt(double(AUFreqs(i)) * double(AUFreqs(j)) * double(AUFreqs(k)))
						  << std::endl;
}

int VibrationalAnalysis::Natoms() const
{
	return m_Natoms;
}

int VibrationalAnalysis::Ncoords() const
{
	return m_Ncoords;
}

int VibrationalAnalysis::Nmodes() const
{
	return m_Nmodes;
}

double VibrationalAnalysis::totalMass() const
{
	return double(m_masses.sum());
}

Eigen::VectorXd VibrationalAnalysis::X() const
{
	return X_min;
}

Eigen::VectorXd VibrationalAnalysis::Q() const
{
	return Q_min;
}

Eigen::MatrixXd VibrationalAnalysis::Fcart() const
{
	return Fcart_min;
}

Eigen::MatrixXd VibrationalAnalysis::Fmwc() const
{
	return Fmwc_min;
}

Eigen::MatrixXd VibrationalAnalysis::Fint() const
{
	return Fint_min;
}

Eigen::MatrixXd VibrationalAnalysis::Lint() const
{
	return Lint_min;
}

Eigen::MatrixXd VibrationalAnalysis::Lmwc() const
{
	return Lmwc_min;
}

Eigen::MatrixXd VibrationalAnalysis::Lcart() const
{
	return Lcrt_min;
}

Eigen::MatrixXd VibrationalAnalysis::D() const
{
	return Dmat;
}

Eigen::VectorXd VibrationalAnalysis::mwcFC() const
{
	return MwcFrcCon;
}

Eigen::VectorXd VibrationalAnalysis::mwcFrq() const
{
	return MwcFreqs;
}

Eigen::VectorXd VibrationalAnalysis::intFC() const
{
	return IntFrcCon;
}

Eigen::VectorXd VibrationalAnalysis::intFrq() const
{
	return IntFreqs;
}

Eigen::VectorXd VibrationalAnalysis::mu() const
{
	return RedMasses;
}

Eigen::Vector3d VibrationalAnalysis::com() const
{
	return m_com;
}

Eigen::Matrix3d VibrationalAnalysis::inert() const
{
	return m_inert;
}

Eigen::Vector3d VibrationalAnalysis::moments() const
{
	return m_moments;
}

Eigen::Matrix3d VibrationalAnalysis::prinAxes() const
{
	return m_prinAxes;
}

Eigen::VectorXd VibrationalAnalysis::masses() const
{
	return m_masses;
}

Eigen::VectorXi VibrationalAnalysis::atomicNumbers() const
{
	return atNums;
}

void VibrationalAnalysis::calcCom()
{
	m_com = Eigen::Vector3d::Zero();

	for (int i = 0; i < m_Natoms; i++)
	{
		m_com(0) += m_masses(i) * X_min(3 * i + 0);
		m_com(1) += m_masses(i) * X_min(3 * i + 1);
		m_com(2) += m_masses(i) * X_min(3 * i + 2);
	}

	m_com /= m_masses.sum();

	X_s = Eigen::VectorXd::Zero(m_Ncoords);

	for (int i = 0; i < m_Natoms; i++)
	{
		X_s(3 * i + 0) = X_min(3 * i + 0) - m_com(0);
		X_s(3 * i + 1) = X_min(3 * i + 1) - m_com(1);
		X_s(3 * i + 2) = X_min(3 * i + 2) - m_com(2);
	}
	Q_s = MassInvMat.inverse() * X_s;
}

void VibrationalAnalysis::calcInert()
{
	m_inert = Eigen::Matrix3d::Zero();

	for (int i = 0; i < m_Natoms; i++)
	{
		m_inert(0,0) += m_masses(i) * (X_s(3 * i + 1) * X_s(3 * i + 1) + X_s(3 * i + 2) * X_s(3 * i + 2));
		m_inert(1,1) += m_masses(i) * (X_s(3 * i + 0) * X_s(3 * i + 0) + X_s(3 * i + 2) * X_s(3 * i + 2));
		m_inert(2,2) += m_masses(i) * (X_s(3 * i + 0) * X_s(3 * i + 0) + X_s(3 * i + 1) * X_s(3 * i + 1));
		m_inert(0,1) -= m_masses(i) *  X_s(3 * i + 0) * X_s(3 * i + 1);
		m_inert(0,2) -= m_masses(i) *  X_s(3 * i + 0) * X_s(3 * i + 2);
		m_inert(1,2) -= m_masses(i) *  X_s(3 * i + 1) * X_s(3 * i + 2);
	}

	m_inert(1,0) = m_inert(0,1);
	m_inert(2,0) = m_inert(0,2);
	m_inert(2,1) = m_inert(1,2);

	Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> InertDiag(m_inert);
	m_moments = InertDiag.eigenvalues();
	m_prinAxes = InertDiag.eigenvectors();

	BigAxes = Eigen::MatrixXd::Zero(m_Ncoords, m_Ncoords);
	for (int i = 0; i < m_Natoms; i++)
		BigAxes.block(3 * i, 3 * i, 3, 3) = m_prinAxes;

	X_sr = BigAxes.transpose() * X_s;
	Q_sr = MassInvMat.inverse() * X_sr;
}

