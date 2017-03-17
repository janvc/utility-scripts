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


#include <iostream>
#include <iomanip>
#include <fstream>
#include <boost/program_options.hpp>
#include <eigen3/Eigen/Eigenvalues>
#include "vibrationalanalysis.h"
#include "GaussFchk.h"
#include "utilities.h"
#include "constants.h"
#include "mkduschclass.h"

int main(int argc, char *argv[])
{
	std::string filename_g;
	std::string filename_e;
//    std::string basename;
	namespace po = boost::program_options;
	po::options_description desc("Allowed options");
	desc.add_options()
		("help,h", "produce this help message")
        ("gfile,g", po::value<std::string>(&filename_g)->required(), "the GS Fchk file")
        ("efile,e", po::value<std::string>(&filename_e)->required(), "the ES Fchk file")
//        ("mctdh-basename,m", po::value<std::string>(&basename)->required(), "the ES Fchk file")
//        ("freqweight,f", "user frequency-weighted coordinates")
	;
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	if (vm.count("help"))
	{
		std::cout << desc << std::endl;
		return 1;
	}
	po::notify(vm);

    std::ifstream gsStream(filename_g);
    std::ifstream esStream(filename_e);
    GaussFchk gsfile(gsStream);
    GaussFchk esfile(esStream);

    MkDuschClass mkdc(gsfile, esfile);
    mkdc.calcMCTDHdata();

    Eigen::MatrixXd Fex = mkdc.F2().asDiagonal();
    Eigen::MatrixXd duschMat = mkdc.Jfull();


    Eigen::VectorXd kVec = mkdc.Kfull();

    double zpeG = 0.0;
    for (int i = 0; i < mkdc.Nmodes(); i++)
        zpeG += std::sqrt(double(mkdc.F1()(i))) / 2.0;
    double zpeE = 0.0;
    for (int i = 0; i < mkdc.Nmodes(); i++)
        zpeE += std::sqrt(double(mkdc.F2()(i))) / 2.0;

    std::cout << "ground state zpe: " << zpeG * 219474.6312 << std::endl;
    std::cout << "excited state zpe: " << zpeE * 219474.6312 << std::endl;

/*
 *  First moment!
 */
    double M1 = mkdc.Ead();

    Eigen::MatrixXd Ftilde = duschMat * Fex * duschMat.transpose();
    double Eshift = 0.5 * kVec.transpose() * Ftilde * kVec;

    double Efreq = 0.0;
    for (int i = 0; i < mkdc.Nmodes(); i++)
        Efreq += 0.25 * (Ftilde(i,i) - mkdc.F1()(i)) / std::sqrt(double(mkdc.F1()(i)));

    double Evert = mkdc.Ead() + Eshift;
    M1 = Evert + Efreq;

    std::cout << "adiabatic energy: " << mkdc.Ead() * Eh2eV << std::endl;
    std::cout << "shift:            " << Eshift * Eh2eV << std::endl;
    std::cout << "vertical energy:  " << Evert * Eh2eV << std::endl;
    std::cout << "frequency shift:  " << Efreq * Eh2eV << std::endl;
    std::cout << "total 1st moment: " << M1 * Eh2eV << std::endl;

/*
 *  Second moment!!
 */
    double M2 = Evert * Evert;

    double term1 = 0.0;
    for (int i = 0; i < mkdc.Nmodes(); i++)
        term1 += (Ftilde(i,i) - mkdc.F1()(i)) / std::sqrt(double(mkdc.F1()(i)));
    term1 *= Evert / 2.0;

    double term2 = 0.0;
    for (int i = 0; i < mkdc.Nmodes(); i++)
        for (int j = i + 1; j < mkdc.Nmodes(); j++)
            term2 += (Ftilde(i,i) - mkdc.F1()(i)) * (Ftilde(j,j) - mkdc.F1()(j)) / std::sqrt(double(mkdc.F1()(i)) * double(mkdc.F1()(j)));
    term2 /= 8.0;

    Eigen::VectorXd grad = kVec.transpose() * Ftilde;
    double term3 = 0.0;
    for (int i = 0; i < mkdc.Nmodes(); i++)
        term3 += grad(i) * grad(i) / std::sqrt(double(mkdc.F1()(i)));
    term3 /= 2.0;

    double term4 = 0.0;
    for (int i = 0; i < mkdc.Nmodes(); i++)
        term4 += (Ftilde(i,i) - mkdc.F1()(i)) * (Ftilde(i,i) - mkdc.F1()(i)) / mkdc.F1()(i);
    term4 *= 3.0 / 8.0;

    double term5 = 0.0;
    for (int i = 0; i < mkdc.Nmodes(); i++)
        for (int j = i + 1; j < mkdc.Nmodes(); j++)
            term5 += Ftilde(i,j) * Ftilde(i,j) / std::sqrt(double(mkdc.F1()(i)) * double(mkdc.F1()(j)));
    term5 /= 2.0;

    std::cout << "term 0: " << M2 * Eh2eV * Eh2eV << std::endl;
    std::cout << "term 1: " << term1 * Eh2eV * Eh2eV << std::endl;
    std::cout << "term 2: " << term2 * Eh2eV * Eh2eV << std::endl;
    std::cout << "term 3: " << term3 * Eh2eV * Eh2eV << std::endl;
    std::cout << "term 4: " << term4 * Eh2eV * Eh2eV << std::endl;
    std::cout << "term 5: " << term5 * Eh2eV * Eh2eV << std::endl;
    std::cout << "total 2nd moment: " << (M2 + term1 + term2 + term3 + term4 + term5) * Eh2eV * Eh2eV << std::endl;


    /*
     * B vector
     */
    Eigen::MatrixXd Om = Eigen::MatrixXd::Zero(mkdc.Nmodes(),mkdc.Nmodes());
    Eigen::MatrixXd Omp = Eigen::MatrixXd::Zero(mkdc.Nmodes(),mkdc.Nmodes());
    Eigen::MatrixXd Ompsq = Eigen::MatrixXd::Zero(mkdc.Nmodes(),mkdc.Nmodes());
    for (int i = 0; i < mkdc.Nmodes(); i++)
    {
        Om(i,i) = std::sqrt(double(mkdc.F2()(i)));
        Omp(i,i) = std::sqrt(double(mkdc.F1()(i)));
        Ompsq(i,i) = std::sqrt(double(Omp(i,i)));
    }

    Eigen::MatrixXd Xmat = duschMat.transpose() * Omp * duschMat + Om;
    Eigen::VectorXd vecDelta = kVec.transpose() * Ompsq;

    Eigen::MatrixXd totalMat = Eigen::MatrixXd::Identity(mkdc.Nmodes(),mkdc.Nmodes()) - (Ompsq * duschMat * Xmat.inverse() * duschMat.transpose() * Ompsq);

    Eigen::VectorXd bVec = 2.0 * vecDelta.transpose() * totalMat;

    Utils::WriteVector(vecDelta);
    Utils::WriteMatrix(totalMat);
    Utils::WriteVector(bVec);


	return 0;
}
