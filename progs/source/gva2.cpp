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
	namespace po = boost::program_options;
	po::options_description desc("Allowed options");
	desc.add_options()
		("help,h", "produce this help message")
        ("gfile,g", po::value<std::string>(&filename_g)->required(), "the GS Fchk file")
        ("efile,e", po::value<std::string>(&filename_e)->required(), "the ES Fchk file")
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

    double zpeG = 0.0;
    for (int i = 0; i < mkdc.Nmodes(); i++)
        zpeG += std::sqrt(double(mkdc.F1()(i))) / 2.0;
    double zpeE = 0.0;
    for (int i = 0; i < mkdc.Nmodes(); i++)
        zpeE += std::sqrt(double(mkdc.F2()(i))) / 2.0;

    std::cout << "ground state zpe: " << zpeG * 219474.6312 << std::endl;
    std::cout << "excited state zpe: " << zpeE * 219474.6312 << std::endl;

/*
 *  Calculate the spectral moments and the B-vector:
 */
    // get the required data:
    int Nmodes = mkdc.Nmodes();
    double Ead = mkdc.Ead();
    Eigen::MatrixXd Fex = mkdc.F2().asDiagonal();
    Eigen::MatrixXd duschMat = mkdc.Jfull();
    Eigen::VectorXd kVec = mkdc.Kfull();
    Eigen::MatrixXd Ftilde = duschMat * Fex * duschMat.transpose();
    Eigen::VectorXd gOm(Nmodes);
    Eigen::VectorXd grad = -kVec.transpose() * Ftilde;
    Eigen::MatrixXd Om = Eigen::MatrixXd::Zero(Nmodes,Nmodes);
    Eigen::MatrixXd Omp = Eigen::MatrixXd::Zero(Nmodes,Nmodes);
    Eigen::MatrixXd Ompsq = Eigen::MatrixXd::Zero(Nmodes,Nmodes);
    for (int i = 0; i < mkdc.Nmodes(); i++)
    {
        Om(i,i) = std::sqrt(double(mkdc.F2()(i)));
        Omp(i,i) = std::sqrt(double(mkdc.F1()(i)));
        Ompsq(i,i) = std::sqrt(double(Omp(i,i)));
        gOm(i) = Omp(i,i);
    }

    // calculate the first moment at 0K:
    double Eshift = 0.5 * kVec.transpose() * Ftilde * kVec;
    double Evert = Ead + Eshift;
    double Efreq = 0.0;
    for (int i = 0; i < mkdc.Nmodes(); i++)
        Efreq += 0.25 * (Ftilde(i,i) - (gOm(i) * gOm(i))) / double(gOm(i));
    double M1 = Evert + Efreq;

    std::cout << "adiabatic energy: " << Ead * Eh2eV << std::endl;
    std::cout << "shift:            " << Eshift * Eh2eV << std::endl;
    std::cout << "vertical energy:  " << Evert * Eh2eV << std::endl;
    std::cout << "frequency shift:  " << Efreq * Eh2eV << std::endl;
    std::cout << "total 1st moment: " << M1 * Eh2eV << std::endl;

    // First moment with pre-excitation:
    std::vector<double> M1pre;
    for (int i = 0; i < Nmodes; i++)
    {
        M1pre.push_back(M1 + 0.5 * ((Ftilde(i,i) - (gOm(i) * gOm(i))) / gOm(i)));
    }

    // calculate the second moment at 0K:
    double term1 = 0.0;
    for (int i = 0; i < Nmodes; i++)
        term1 += (Ftilde(i,i) - (gOm(i) * gOm(i))) / gOm(i);
    term1 *= Evert / 2.0;

    double term2 = 0.0;
    for (int i = 0; i < Nmodes; i++)
        for (int j = i + 1; j < Nmodes; j++)
            term2 += (Ftilde(i,i) - (gOm(i) * gOm(i))) * (Ftilde(j,j) - (gOm(j) * gOm(j))) / (gOm(i) * gOm(j));
    term2 /= 8.0;

    double term3 = 0.0;
    for (int i = 0; i < Nmodes; i++)
        term3 += grad(i) * grad(i) / gOm(i);
    term3 /= 2.0;

    double term4 = 0.0;
    for (int i = 0; i < Nmodes; i++)
        term4 += (Ftilde(i,i) - (gOm(i) * gOm(i))) * (Ftilde(i,i) - (gOm(i) * gOm(i))) / (gOm(i) * gOm(i));
    term4 *= 3.0 / 8.0;

    double term5 = 0.0;
    for (int i = 0; i < Nmodes; i++)
        for (int j = i + 1; j < Nmodes; j++)
            term5 += Ftilde(i,j) * Ftilde(i,j) / (gOm(i) * gOm(j));
    term5 /= 2.0;

    double M2 = Evert * Evert + term1 + term2 + term3 + term4 + term5;

    std::cout << "term 0: " << Evert * Evert * Eh2eV * Eh2eV << std::endl;
    std::cout << "term 1: " << term1 * Eh2eV * Eh2eV << std::endl;
    std::cout << "term 2: " << term2 * Eh2eV * Eh2eV << std::endl;
    std::cout << "term 3: " << term3 * Eh2eV * Eh2eV << std::endl;
    std::cout << "term 4: " << term4 * Eh2eV * Eh2eV << std::endl;
    std::cout << "term 5: " << term5 * Eh2eV * Eh2eV << std::endl;
    std::cout << "total 2nd moment: " << M2 * Eh2eV * Eh2eV << std::endl;

    // calculate the spectrum width:
    double width = std::sqrt(M2 - (M1 * M1));
    std::cout << "spectrum width: " << width * Eh2eV << std::endl;

    // Second moment with pre-excitation
    std::vector<double> M2pre;
    std::vector<double> T1pre;
    std::vector<double> T2pre;
    std::vector<double> T3pre;
    std::vector<double> T4pre;
    std::vector<double> T5pre;
    std::vector<double> exWidth;

    for (int k = 0; k < Nmodes; k++)
    {
        // first term:
        double tmp1 = Evert * (Ftilde(k,k) - (gOm(k) * gOm(k))) / gOm(k);
        T1pre.push_back(tmp1);

        // second term:
        double tmp2 = 0.0;
        for (int i = 0; i < Nmodes; i++)
            if (i != k)
                tmp2 += (Ftilde(i,i) - (gOm(i) * gOm(i))) / gOm(i);
        tmp2 *= 0.5 * (Ftilde(k,k) - (gOm(k) * gOm(k))) / gOm(k);
        T2pre.push_back(tmp2);

        // third term:
        double tmp3 = grad(k) * grad(k) / gOm(k);
        T3pre.push_back(tmp3);

        // fourth term:
        double tmp4 = 15.0 * (Ftilde(k,k) - (gOm(k) * gOm(k))) * (Ftilde(k,k) - (gOm(k) * gOm(k))) / (8.0 * gOm(k) * gOm(k));
        T4pre.push_back(tmp4);

        // fifth term:
        double tmp5 = 0.0;
        for (int i = 0; i < Nmodes; i++)
            if (i != k)
                tmp5 += 2.0 * Ftilde(i,k) * Ftilde(i,k) / (gOm(i) * gOm(k));
        T5pre.push_back(tmp5);

        M2pre.push_back(M2 + tmp1 + tmp2 + tmp3 + tmp4 + tmp5);

        exWidth.push_back(std::sqrt(M2pre.back() - (M1pre[k] * M1pre[k])));
    }

    std::cout << "First and second moments after pre-excitation:\n";
    std::cout << "  mode          M1         delta M1       M2 term1       M2 term2       M2 term3       M2 term4       M2 term5"
              << "             M2          width    delta width\n";
    for (int i = 0; i < Nmodes; i++)
        std::cout << std::setw(5) << i + 1
                  << std::setw(15) << M1pre[i] * Eh2eV
                  << std::setw(15) << (M1pre[i] - M1) * Eh2eV
                  << std::setw(15) << T1pre[i] * Eh2eV * Eh2eV
                  << std::setw(15) << T2pre[i] * Eh2eV * Eh2eV
                  << std::setw(15) << T3pre[i] * Eh2eV * Eh2eV
                  << std::setw(15) << T4pre[i] * Eh2eV * Eh2eV
                  << std::setw(15) << T5pre[i] * Eh2eV * Eh2eV
                  << std::setw(15) << M2pre[i] * Eh2eV * Eh2eV
                  << std::setw(15) << exWidth[i] * Eh2eV
                  << std::setw(15) << (exWidth[i] - width) * Eh2eV
                  << std::endl;

    // for testing:
//    duschMat = Eigen::MatrixXd::Identity(Nmodes, Nmodes);
//    Omp = Om;
//    for (int i = 0; i < Nmodes; i++)
//        Ompsq(i,i) = std::sqrt(double(Om(i,i)));
    /*
     * B vector
     */
    Eigen::MatrixXd Xmat = duschMat.transpose() * Omp * duschMat + Om;
    Eigen::VectorXd vecDelta = kVec.transpose() * Ompsq;

    Eigen::MatrixXd totalMat = Ompsq * duschMat * Xmat.inverse() * duschMat.transpose() * Ompsq;

    Eigen::VectorXd bVec = 2.0 * vecDelta.transpose() * (Eigen::MatrixXd::Identity(Nmodes, Nmodes) - totalMat);

    std::cout << "Dimensionless displacements along GS normal modes:\n";
    Utils::WriteVector(vecDelta);

    std::cout << "X matrix:\n";
    Utils::WriteMatrix(Xmat);

    std::cout << "total matrix for B vector:\n";
    Utils::WriteMatrix(totalMat);

    std::cout << "norm of total matrix: " << totalMat.norm() << std::endl;

    std::cout << "B vector:\n";
    Utils::WriteVector(bVec);

    /*
     * sort the B vector to determine the best mode(s) for VIPER:
     */
    std::vector<std::pair<double,int> > pairVec;
    for (int i = 0; i < Nmodes; i++)
        pairVec.push_back(std::make_pair(std::abs(double(bVec(i))), i));

    std::sort(pairVec.begin(), pairVec.end());

    std::cout << "B vector sorted by increasing magnitude\nwith corresponding normal mode\n";
    for (int i = 0; i < Nmodes; i++)
        std::cout << pairVec[i].first << "   " << pairVec[i].second + 1 << std::endl;


	return 0;
}
