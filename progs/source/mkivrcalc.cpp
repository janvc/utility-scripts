/*
 * Copyright 2016 Jan von Cosel
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
#include <fstream>
#include <iomanip>
#include <boost/program_options.hpp>
#include <eigen3/Eigen/SVD>
#include <eigen3/Eigen/Eigenvalues>
#include "vibrationalanalysis.h"
#include "GaussFchk.h"
#include "utilities.h"

int main(int argc, char *argv[])
{
	std::string gsFchkName, esFchkName, AnharmLogName, MCTDHbaseName;
	double deriv_thres;
	bool doExcitedState = false;
    bool freqWeight = false;
	namespace po = boost::program_options;
	po::options_description desc("Allowed options");
	desc.add_options()
		("help,h", "produce this help message")
		("gsfile,g", po::value<std::string>(&gsFchkName)->required(), "the ground state Fchk file")
		("esfile,e", po::value<std::string>(&esFchkName), "the excited state Fchk file")
		("ahlog,l", po::value<std::string>(&AnharmLogName)->required(), "the anharmonic log file")
		("mctdh,m", po::value<std::string>(&MCTDHbaseName)->required(), "the MCTDH base filename")
		("thres,t", po::value<double>(&deriv_thres)->required(), "threshold for the anharmonic couplings")
        ("fw,w", "use frequency weighted coordinates")
	;
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	if (vm.count("help"))
	{
		std::cout << desc << std::endl;
		return 1;
	}
	if (vm.count("esfile"))
		doExcitedState = true;
    if (vm.count("fw"))
        freqWeight = true;

	po::notify(vm);


	std::cout << "Ground state Fchk file:       " << gsFchkName << std::endl;
	std::cout << "Excited state Fchk file:      " << esFchkName << std::endl;
	std::cout << "Anharmonic log file:          " << AnharmLogName << std::endl;
	std::cout << "Basename for MCTDH files:     " << MCTDHbaseName << std::endl;
	std::cout << "Threshold for 3rd/4th derivs: " << deriv_thres << std::endl << std::endl;

    if (freqWeight)
    {
        std::cout << "---------------------------------------------------------------\n";
        std::cout << "           WE ARE USING FREQUENCY-WEIGHTED UNITS!\n";
        std::cout << "---------------------------------------------------------------\n";
    }

	std::ifstream gsFchkFile(gsFchkName, std::ifstream::in);
	GaussFchk gsFchk(gsFchkFile);


	VibrationalAnalysis GState(gsFchk);

	GState.readAnharm(AnharmLogName);


	int Natoms = GState.Natoms();
	int Ncoords = GState.Ncoords();
	int Nmodes = GState.Nmodes();

    Eigen::VectorXd f1(Nmodes);
    arma::Cube<double> phi(Nmodes, Nmodes, Nmodes);
    Eigen::VectorXd diagf(Nmodes);

    if (freqWeight)
        for (int i = 0; i < Nmodes; i++)
        {
            f1(i) = std::sqrt(double(GState.intFC()(i)));
            diagf(i) = double(GState.diag4thDerivs()(i)) / double(GState.intFC()(i));

            for (int j = i; j < Nmodes; j++)
                for (int k = j; k < Nmodes; k++)
                    phi(i,j,k) = GState.thirdDerivs()(i,j,k)
                               / pow(double(GState.intFC()(i))
                                   * double(GState.intFC()(j))
                                   * double(GState.intFC()(k)), 0.25);
        }
    else
    {
        f1 = GState.intFC();
        phi = GState.thirdDerivs();
        diagf = GState.diag4thDerivs();
    }


	std::cout << "--------------------------------------------\n";
	std::cout << "     INFORMATION ABOUT THE GROUND STATE\n";
	std::cout << "--------------------------------------------\n\n";

	std::cout << "Number of atoms:          " << Natoms << std::endl;
	std::cout << "Number of normal modes:   " << Nmodes << std::endl;


	std::cout << "\nmass-weighted force constants in [Eh / a0**2 * me]," << std::endl
			  << "frequencies in [cm**-1] and anharmonicities:\n";
	std::cout << " mode     force constant     frequency     3rd derivative      4th derivative\n";
	for (int i = 0; i < Nmodes; i++)
		std::cout << std::setw(5) << i + 1
				  << std::setw(20) << std::scientific << std::setprecision(9) << double(GState.intFC()(i))
				  << std::setw(13) << std::fixed << std::setprecision(4) << double(GState.intFrq()(i))
				  << std::setw(20) << std::scientific << std::setprecision(9) << phi(i,i,i)
				  << std::setw(20) << std::scientific << std::setprecision(9) << double(diagf(i)) << std::endl;


	/*
	 * Determine the strongest 3rd order couplings for each normal mode:
	 */
	std::vector<double> maxCouplings(Nmodes);
	std::vector<double> avgCouplings(Nmodes);
	std::vector<int> maxI(Nmodes);
	std::vector<int> maxJ(Nmodes);
	std::vector<int> maxK(Nmodes);

	for (int i = 0; i < Nmodes; i++)
	{
		double coupling = 0.0;
		avgCouplings.at(i) = 0.0;
		int imax = 0, jmax = 0, kmax = 0;
		for (int j = 0; j < Nmodes; j++)
			for (int k = 0; k < Nmodes; k++)
				if (!(i == j && i == k))
				{
					avgCouplings.at(i) += std::abs(phi(i,j,k)) / double(Nmodes);
					if (std::abs(phi(i,j,k)) > std::abs(coupling))
					{
						coupling = phi(i,j,k);
						imax = i;
						jmax = j;
						kmax = k;
					}
				}
		maxCouplings.at(i) = coupling;
		maxI.at(i) = imax;
		maxJ.at(i) = jmax;
		maxK.at(i) = kmax;
	}

	std::cout << "\nStrongest 3rd order coupling coefficients for each normal mode:\n";
	std::cout << "    i    j    k       coupling\n";
	for (int i = 0; i < Nmodes; i++)
		std::cout << std::setw(5) << maxI.at(i) + 1
				  << std::setw(5) << maxJ.at(i) + 1
				  << std::setw(5) << maxK.at(i) + 1
				  << std::setw(20) << std::scientific << std::setprecision(9) << maxCouplings.at(i) << std::endl;

	std::cout << "\nAverage coupling of each mode:\n";
	std::cout << " mode  avg. coupling\n";
	for (int i = 0; i < Nmodes; i++)
		std::cout << std::setw(5) << i + 1
				  << std::setw(20) << std::scientific << std::setprecision(9) << avgCouplings.at(i) << std::endl;


	std::ofstream IVRinput(MCTDHbaseName + "_ivr.inp");
	std::ofstream IVRoper(MCTDHbaseName + "_ivr.op");

	IVRinput.precision(1);
	IVRoper.precision(8);

	/*
	 * Create the primitive basis for a pure IVR calculation.
	 * This will be overwritten if we also have an excited state present.
	 */
	std::vector<int> nBasis(Nmodes);
	std::vector<double> lowBound(Nmodes);
	std::vector<double> uppBound(Nmodes);

	for (int i = 0; i < Nmodes; i++)
	{
//        nBasis[i] = lrint(-0.7* (freqWeight ? 2.0 : 1.0) * log(double(f1(i)))) + 6;
        nBasis[i] = lrint(-1.4 * log(double(f1(i)))) + 6;
//        lowBound[i] = -7.5 / (sqrt(2.0) * (freqWeight ? std::sqrt(double(f1(i))) : pow(double(f1(i)), 0.25)));
//        uppBound[i] =  7.5 / (sqrt(2.0) * (freqWeight ? std::sqrt(double(f1(i))) : pow(double(f1(i)), 0.25)));
//        if (freqWeight)
//        {
//            lowBound[i] = -7.5 / (sqrt(2.0) * std::sqrt(double(f1(i))));
//            uppBound[i] =  7.5 / (sqrt(2.0) * std::sqrt(double(f1(i))));
//        }
//        else
//        {
            lowBound[i] = -7.5 / (sqrt(2.0) * pow(double(f1(i)), 0.25));
            uppBound[i] =  7.5 / (sqrt(2.0) * pow(double(f1(i)), 0.25));
//        }
	}


	if (doExcitedState)
	{
		std::cout << "Excited state Fchk file:  " << esFchkName << std::endl;
		std::ifstream esFchkFile(esFchkName, std::ifstream::in);
		GaussFchk esFchk(esFchkFile);
		VibrationalAnalysis EState(esFchk);
		double deltaE = esFchk.ReadReal("Total Energy") - gsFchk.ReadReal("Total Energy");

        Eigen::VectorXd f2(Nmodes);
        if (freqWeight)
            for (int i = 0; i < Nmodes; i++)
                f2(i) = double(EState.intFC()(i)) / std::sqrt(double(f1(i)));
        else
            f2 = EState.intFC();

		std::ofstream SpecInput(MCTDHbaseName + "_spec.inp");
		std::ofstream SpecOper(MCTDHbaseName + "_spec.op");
		SpecInput.precision(1);
		SpecOper.precision(8);


		/*
		 * Perform SVD on the normal modes:
		 */
		Eigen::JacobiSVD<Eigen::MatrixXd> SVDg(GState.Lmwc(), Eigen::ComputeThinU | Eigen::ComputeThinV);
		Eigen::JacobiSVD<Eigen::MatrixXd> SVDe(EState.Lmwc(), Eigen::ComputeThinU | Eigen::ComputeThinV);
		Eigen::MatrixXd NMOg = SVDg.matrixU() * SVDg.matrixV().transpose();
		Eigen::MatrixXd NMOe = SVDe.matrixU() * SVDe.matrixV().transpose();


		/*
		 * calculate the rotation matrix to minimize the RMSD between
		 * the ground- and excited-state structures:
		 */
		Eigen::Vector3d COMg = GState.com();
		Eigen::Vector3d COMe = EState.com();

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

		double RMSD_before = 0;
		for (int i = 0; i < Ncoords; i++)
			RMSD_before += std::abs(double(XSg(i) - XSe(i))) * std::abs(double(XSg(i) - XSe(i)));
		RMSD_before /= Natoms;

		Eigen::Matrix3d corr(Eigen::Matrix3d::Zero());
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				for(int k = 0; k < Natoms; k++)
					corr(i,j) += XSg(3*k+i) * XSe(3*k+j);

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

		Eigen::VectorXd XRg(BigRotMat * XSg);
		Eigen::VectorXd QRg(BigRotMat * QSg);

		double RMSD_after = 0;
		for (int i = 0; i < Ncoords; i++)
			RMSD_after += std::abs(double(XRg(i) - XSe(i))) * std::abs(double(XRg(i) - XSe(i)));
		RMSD_after /= Natoms;

		Eigen::MatrixXd J(Eigen::MatrixXd::Zero(Nmodes, Nmodes));
		Eigen::VectorXd K(Eigen::VectorXd::Zero(Nmodes));

		if (RMSD_after < RMSD_before)
		{
			Eigen::MatrixXd NMRg(BigRotMat * NMOg);

			J = NMOe.transpose() * NMRg;
			K = NMOe.transpose() * (QRg - QSe);
		}
		else
		{
			J = NMOe.transpose() * NMOg;
			K = NMOe.transpose() * (QSg - QSe);
		}


		/*
		 * Calculate the required data for the spectrum calculation hamiltonian
		 */
		Eigen::VectorXd fp = Eigen::VectorXd::Zero(Nmodes);
		for (int m = 0; m < Nmodes; m++)
			for (int n = 0; n < Nmodes; n++)
				fp(m) += f2(n) * J(n,m) * J(n,m);

		Eigen::VectorXd kappa = Eigen::VectorXd::Zero(Nmodes);
		for (int m = 0; m < Nmodes; m++)
			for (int n = 0; n < Nmodes; n++)
				kappa(m) += f2(n) * K(n) * J(n,m);

		Eigen::MatrixXd duschRot = Eigen::MatrixXd::Zero(Nmodes, Nmodes);
		for (int m = 0; m < Nmodes; m++)
			for (int o = m + 1; o < Nmodes; o++)
			{
				for (int n = 0; n < Nmodes; n++)
					duschRot(m,o) += f2(n) * J(n,m) * J(n,o);
				duschRot(o,m) = duschRot(m,o);
			}

		Eigen::VectorXd d(Nmodes);
		for (int i = 0; i < Nmodes; i++)
			d(i) = 0.5 * f2(i) * K(i) * K(i);


		/*
		 * Determine the required primitive basis for the two calculations, since they
		 * need to be equal
		 */
		for (int i = 0; i < Nmodes; i++)
		{
			nBasis[i] = lrint(-0.7 * log(double(fp(i)))) + 6;
			lowBound[i] = -double(kappa(i) / fp(i)) - (std::abs(double(kappa(i) / fp(i))) + 7.5 / (sqrt(2.0) * pow(double(fp(i)), 0.25)));
			uppBound[i] = -double(kappa(i) / fp(i)) + (std::abs(double(kappa(i) / fp(i))) + 7.5 / (sqrt(2.0) * pow(double(fp(i)), 0.25)));
		}


		/*
		 * Create the MCTDH operator file for the second calculation:
		 * the spectrum calculation
		 */
		SpecOper << "op_define-section\n";
		SpecOper << "    title\n";
		SpecOper << "        " << MCTDHbaseName << std::endl;
		SpecOper << "    end-title\n";
		SpecOper << "end-op_define-section\n\n";

		SpecOper << "parameter-section\n";
		for (int i = 0; i < Nmodes; i++)
			SpecOper << "    mass_q_" << std::setfill('0') << std::setw(3) << i + 1
					 << "  =  1.0\n";
		for (int i = 0; i < Nmodes; i++)
		{
			SpecOper << "    f1_" << std::setfill('0') << std::setw(3) << i + 1 << "      = ";
            Utils::WriteFortranNumber(SpecOper, double(f1(i)));
			SpecOper << std::endl;
		}
		for (int i = 0; i < Nmodes; i++)
		{
			SpecOper << "    f2_" << std::setfill('0') << std::setw(3) << i + 1 << "      = ";
            Utils::WriteFortranNumber(SpecOper, double(f2(i)));
			SpecOper << std::endl;
		}
		for (int i = 0; i < Nmodes; i++)
		{
			SpecOper << "    fp_" << std::setfill('0') << std::setw(3) << i + 1 << "      = ";
            Utils::WriteFortranNumber(SpecOper, double(fp(i)));
			SpecOper << std::endl;
		}
		for (int i = 0; i < Nmodes; i++)
			for (int j = i + 1; j < Nmodes; j++)
			{
				SpecOper << "    phi_" << std::setfill('0') << std::setw(3) << i + 1
							 << "_" << std::setfill('0') << std::setw(3) << j + 1 << " = ";
                Utils::WriteFortranNumber(SpecOper, double(duschRot(i,j)));
				SpecOper << std::endl;
			}
		// the first-order coefficients (shifts)
		for (int i = 0; i < Nmodes; i++)
		{
			SpecOper << "    kappa_" << std::setfill('0') << std::setw(3) << i + 1 << "   = ";
            Utils::WriteFortranNumber(SpecOper, double(kappa(i)));
			SpecOper << std::endl;
		}
		// the energy offsets
		for (int i = 0; i < Nmodes; i++)
		{
			SpecOper << "    d_" << std::setfill('0') << std::setw(3) << i + 1 << "       = ";
            Utils::WriteFortranNumber(SpecOper, double(d(i)));
			SpecOper << std::endl;
		}
		double zpe1 = 0.0;
		for (int i = 0; i < Nmodes; i++)
			zpe1 += 0.5 * sqrt(double(f1(i)));
		SpecOper << "    dE          = ";
        Utils::WriteFortranNumber(SpecOper, deltaE - zpe1);
		SpecOper << "\nend-parameter-section\n\n";

		SpecOper << "hamiltonian-section";
		for (int i = 0; i < Nmodes; i++)
		{
			if (i % 8 == 0)
				SpecOper << std::endl << "modes";
			SpecOper << " | q_" << std::setfill('0') << std::setw(3) << i + 1;
		}
		SpecOper << std::endl;
		for (int i = 0; i < Nmodes; i++)
			SpecOper << "1.0         |" << i + 1 << " KE\n";
		for (int i = 0; i < Nmodes; i++)
			SpecOper << "0.5*fp_" << std::setfill('0') << std::setw(3) << i + 1
					 << "  |" << i + 1 << " q^2\n";
		for (int i = 0; i < Nmodes; i++)
			for (int j = i + 1; j < Nmodes; j++)
				SpecOper << "phi_" << std::setfill('0') << std::setw(3) << i + 1
						 << "_" << std::setfill('0') << std::setw(3) << j + 1
						 << " |" << i + 1 << " q"
						 << " |" << j + 1 << " q\n";
		for (int i = 0; i < Nmodes; i++)
			SpecOper << "kappa_" << std::setfill('0') << std::setw(3) << i + 1
					 << "   |" << i + 1 << " q\n";
		for (int i = 0; i < Nmodes; i++)
				SpecOper << "d_" << std::setfill('0') << std::setw(3) << i + 1
						 << "       |" << i + 1 << " 1\n";
		SpecOper << "dE          |1 1\n";
		SpecOper << "end-hamiltonian-section\n\n";

		SpecOper << "end-operator\n\n";


		/*
		 * Create the MCTDH input file for the spectrum calculation
		 */
		SpecInput << "run-section\n";
		SpecInput << "    name =\n";
		SpecInput << "    propagation\n";
		SpecInput << "    tfinal =\n";
		SpecInput << "    tout =\n";
		SpecInput << "    tpsi =\n";
		SpecInput << "    psi gridpop auto steps graphviz\n";
		SpecInput << "end-run-section\n\n";

		SpecInput << "operator-section\n";
		SpecInput << "    opname = " << MCTDHbaseName << "_spec" << std::endl;
		SpecInput << "end-operator-section\n\n";

		SpecInput << "mlbasis-section\n";
		for (int i = 0; i < Nmodes - 1; i += 2)
		{
			if (Nmodes - i == 3)
				SpecInput << "    [q_" << std::setfill('0') << std::setw(3) << i + 1
							  << " q_" << std::setfill('0') << std::setw(3) << i+1 + 1
							  << " q_" << std::setfill('0') << std::setw(3) << i+2 + 1 << "]\n";
			else
				SpecInput << "    [q_" << std::setfill('0') << std::setw(3) << i+0 + 1
							  << " q_" << std::setfill('0') << std::setw(3) << i+1 + 1 << "]\n";
		}
		SpecInput << "end-mlbasis-section\n\n";

		SpecInput << "pbasis-section\n";
		for (int i = 0; i < Nmodes; i++)
			SpecInput << "    q_" << std::setfill('0') << std::setw(3) << i + 1
					  << "  ho  " << std::setw(3) << std::setfill(' ') << nBasis[i] << "  xi-xf  "
					  << std::fixed << std::setfill(' ') << std::setw(8)
					  << lowBound[i]
					  << std::fixed << std::setfill(' ') << std::setw(8)
					  << uppBound[i] << std::endl;
		SpecInput << "end-pbasis-section\n\n";

		SpecInput << "integrator-section\n";
		SpecInput << "    vmf\n";
		SpecInput << "    abm = 6, 1.0d-7, 0.01d0\n";
		SpecInput << "end-integrator-section\n\n";

		SpecInput << "init_wf-section\n";
		SpecInput << "    build\n";
		for (int i = 0; i < Nmodes; i++)
			SpecInput << "        q_" << std::setfill('0') << std::setw(3) << i + 1
					 << "  eigenf"
					 << "  Eq_" << std::setfill('0') << std::setw(3) << i + 1
					 << "  pop = 1\n";
		SpecInput << "    end-build\n";
		SpecInput << "end-init_wf-section\n\n";
		SpecInput << "end-input\n\n";
	}

	/*
	 * Create the MCTDH operator file for the first calculation:
	 * The anharmonic IVR propagation
	 */
	IVRoper << "op_define-section\n";
	IVRoper << "    title\n";
    IVRoper << "        " << MCTDHbaseName << (freqWeight ? " frequency weighted" : "") << std::endl;
	IVRoper << "    end-title\n";
	IVRoper << "end-op_define-section\n\n";

	IVRoper << "parameter-section\n";
	for (int i = 0; i < Nmodes; i++)
		IVRoper << "    mass_q_" << std::setfill('0') << std::setw(3) << i + 1
				<< "      =  1.0\n";
	for (int i = 0; i < Nmodes; i++)
	{
		IVRoper << "    f_" << std::setfill('0') << std::setw(3) << i + 1 << "           = ";
        Utils::WriteFortranNumber(IVRoper, double(f1(i)));
		IVRoper << std::endl;
	}
	IVRoper << "    # cubic couplings:\n";
	IVRoper << "    # the numbers here are equal to 1/6 * phi_ijk\n";
	for (int i = 0; i < Nmodes; i++)
		for (int j = i; j < Nmodes; j++)
			for (int k = j; k < Nmodes; k++)
				if (std::abs(phi(i,j,k)) >= deriv_thres)
				{
					IVRoper << "    phi_" << std::setfill('0') << std::setw(3) << i + 1
								   << "_" << std::setfill('0') << std::setw(3) << j + 1
								   << "_" << std::setfill('0') << std::setw(3) << k + 1
								   << " = ";
                    Utils::WriteFortranNumber(IVRoper, phi(i,j,k) / 6.0);
					IVRoper << std::endl;
				}
	IVRoper << "    # diagonal quartic force constants:\n";
	IVRoper << "    # the numbers here are equal to 1/24 * fq_iiii\n";
	for (int i = 0; i < Nmodes; i++)
		if (std::abs(double(diagf(i))) >= deriv_thres)
		{
			IVRoper << "    fdia_" << std::setfill('0') << std::setw(3) << i + 1 << "        = ";
            Utils::WriteFortranNumber(IVRoper, double(diagf(i)) / 24.0);
			IVRoper << std::endl;
		}
	IVRoper << "end-parameter-section\n\n";

	IVRoper << "hamiltonian-section\n";
	for (int i = 0; i < Nmodes; i++)
	{
		if (i % 8 == 0)
			IVRoper << std::endl << "modes";
		IVRoper << " | q_" << std::setfill('0') << std::setw(3) << i + 1;
	}
	IVRoper << std::endl;

	for (int i = 0; i < Nmodes; i++)
        if (freqWeight)
            IVRoper << "f_" << std::setfill('0') << std::setw(3) << i + 1 << "       |" << i + 1 << " KE\n";
        else
            IVRoper << "1.0       |" << i + 1 << " KE\n";

	for (int i = 0; i < Nmodes; i++)
		IVRoper << "0.5*f_" << std::setfill('0') << std::setw(3) << i + 1 << "       |" << i + 1 << " q^2\n";

	for (int i = 0; i < Nmodes; i++)
		for (int j = i; j < Nmodes; j++)
			for (int k = j; k < Nmodes; k++)
				if (std::abs(phi(i,j,k)) >= deriv_thres)
				{
					if (i == j && i == k)	// diagonal elements
						IVRoper << "phi_" << std::setfill('0') << std::setw(3) << i + 1
								   << "_" << std::setfill('0') << std::setw(3) << i + 1
								   << "_" << std::setfill('0') << std::setw(3) << i + 1
								   << " |" << i + 1 << " q^3\n";
					else if (i == j)
						IVRoper << "phi_" << std::setfill('0') << std::setw(3) << i + 1
								   << "_" << std::setfill('0') << std::setw(3) << i + 1
								   << "_" << std::setfill('0') << std::setw(3) << k + 1
								   << " |" << i + 1 << " q^2"
								   << " |" << k + 1 << " q\n";
					else if (j == k)
						IVRoper << "phi_" << std::setfill('0') << std::setw(3) << i + 1
								   << "_" << std::setfill('0') << std::setw(3) << j + 1
								   << "_" << std::setfill('0') << std::setw(3) << j + 1
								   << " |" << i + 1 << " q"
								   << "   |" << j + 1 << " q^2\n";
					else
						IVRoper << "phi_" << std::setfill('0') << std::setw(3) << i + 1
								   << "_" << std::setfill('0') << std::setw(3) << j + 1
								   << "_" << std::setfill('0') << std::setw(3) << k + 1
								   << " |" << i + 1 << " q"
								   << "   |" << j + 1 << " q"
								   << "   |" << k + 1 << " q\n";
				}

	for (int i = 0; i < Nmodes; i++)
		if (std::abs(double(diagf(i))) >= deriv_thres)
			IVRoper << "fdia_" << std::setfill('0') << std::setw(3) << i + 1
					<< "        |" << i + 1 << " q^4\n";

	IVRoper << "end-hamiltonian-section\n\n";

	for (int i = 0; i < Nmodes; i++)
	{
		IVRoper << "hamiltonian-section_Eq_" << std::setfill('0') << std::setw(3) << i + 1 << std::endl;
		IVRoper << "usediag\n";
		IVRoper << "modes | q_" << std::setfill('0') << std::setw(3) << i + 1 << std::endl;
        if (freqWeight)
            IVRoper << "f_"  << std::setfill('0') << std::setw(3) << i + 1 << "       " << "|1 KE\n";
        else
            IVRoper << "1.0             " << "|1 KE\n";
		IVRoper << "0.5*f_" << std::setfill('0') << std::setw(3) << i + 1 << "       |1 q^2\n";
		if (std::abs(phi(i,i,i)) >= deriv_thres)
			IVRoper << "phi_" << std::setfill('0') << std::setw(3) << i + 1 << "_" << std::setfill('0') << std::setw(3) << i + 1 << "_" << std::setfill('0') << std::setw(3) << i + 1 << " |1 q^3\n";
		if (std::abs(double(diagf(i))) >= deriv_thres)
			IVRoper << "fdia_" << std::setfill('0') << std::setw(3) << i + 1 << "        |1 q^4\n";
		IVRoper << "end-hamiltonian-section\n\n";
	}

	IVRoper << "end-operator\n\n";


	/*
	 * Create the MCTDH input file for the IVR propagation
	 */
	IVRinput << "run-section\n";
	IVRinput << "    name =\n";
	IVRinput << "    propagation\n";
	IVRinput << "    tfinal =\n";
	IVRinput << "    tout =\n";
	IVRinput << "    tpsi =\n";
	IVRinput << "    psi gridpop auto steps graphviz\n";
	IVRinput << "end-run-section\n\n";

	IVRinput << "operator-section\n";
	IVRinput << "    opname = " << MCTDHbaseName << "_ivr" << std::endl;
	IVRinput << "end-operator-section\n\n";

	IVRinput << "mlbasis-section\n";
	for (int i = 0; i < Nmodes - 1; i += 2)
	{
		if (Nmodes - i == 3)
			IVRinput << "    [q_" << std::setfill('0') << std::setw(3) << i + 1
						 << " q_" << std::setfill('0') << std::setw(3) << i+1 + 1
						 << " q_" << std::setfill('0') << std::setw(3) << i+2 + 1 << "]\n";
		else
			IVRinput << "    [q_" << std::setfill('0') << std::setw(3) << i+0 + 1
						 << " q_" << std::setfill('0') << std::setw(3) << i+1 + 1 << "]\n";
	}
	IVRinput << "end-mlbasis-section\n\n";

	IVRinput << "pbasis-section\n";
	for (int i = 0; i < Nmodes; i++)
		IVRinput << "    q_" << std::setfill('0') << std::setw(3) << i + 1
				 << "  ho  " << std::setw(3) << std::setfill(' ') << nBasis[i] << "  xi-xf  "
				 << std::fixed << std::setfill(' ') << std::setw(8)
				 << lowBound[i]
				 << std::fixed << std::setfill(' ') << std::setw(8)
				 << uppBound[i] << std::endl;
	IVRinput << "end-pbasis-section\n\n";

	IVRinput << "integrator-section\n";
	IVRinput << "    vmf\n";
	IVRinput << "    abm = 6, 1.0d-7, 0.01d0\n";
	IVRinput << "end-integrator-section\n\n";

	IVRinput << "init_wf-section\n";
	IVRinput << "    build\n";
	for (int i = 0; i < Nmodes; i++)
		IVRinput << "        q_" << std::setfill('0') << std::setw(3) << i + 1
				 << "  eigenf"
				 << "  Eq_" << std::setfill('0') << std::setw(3) << i + 1
				 << "  pop = 1\n";
	IVRinput << "    end-build\n";
	IVRinput << "end-init_wf-section\n\n";
	IVRinput << "end-input\n\n";
}
