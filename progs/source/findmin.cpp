/*
 * Copyright 2018 Jan von Cosel
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
#include <random>
#include <boost/program_options.hpp>
#include "GaussFchk.h"
#include "vibrationalanalysis.h"

int main(int argc, char *argv[])
{
    std::string harmName;
    std::string ahName;
    double deriv_thres;
    int Npoints;
    namespace po = boost::program_options;
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help,h", "produce this help message")
        ("fchk,f", po::value<std::string>(&harmName)->required(), "the ground state harmonic Fchk file")
        ("ahlog,l", po::value<std::string>(&ahName)->required(), "the anharmonic log file")
        ("thres,t", po::value<double>(&deriv_thres)->required(), "threshold for the anharmonic couplings")
        ("npoints,n", po::value<int>(&Npoints)->required(), "number of points to compute")
    ;
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    if (vm.count("help"))
    {
        std::cout << desc << std::endl;
        return 1;
    }
    po::notify(vm);

    std::ifstream harmFchkFile(harmName, std::ifstream::in);
    GaussFchk harmFchk(harmFchkFile);
    VibrationalAnalysis GState(harmFchk);
    GState.readAnharm(ahName);

    int Nmodes = GState.Nmodes();

    Eigen::VectorXd f1(Nmodes);
    arma::Cube<double> phi(Nmodes, Nmodes, Nmodes);
    Eigen::VectorXd diagf(Nmodes);

    f1 = GState.intFC();
    phi = GState.thirdDerivs();
    diagf = GState.diag4thDerivs();

    std::vector<double> Bound(Nmodes);

    for (int i = 0; i < Nmodes; i++)
    {
        Bound[i] =  10.0 / (sqrt(2.0) * pow(double(f1(i)), 0.25));
    }

    double minEner = 0.0;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);

    for (int i = 0; i < Npoints; i++)
    {
//        std::cout << "point " << i + 1 << "\n";

        std::vector<double> positions(Nmodes);

        for (int j = 0; j < Nmodes; j++)
            positions[j] = Bound[j] * (2.0 * dis(gen) - 1.0);

//        for (int j = 0; j < Nmodes; j++)
//            std::cout << j + 1 << "   " << positions[j] << "\n";

        double ener = 0.0;
        for (int j = 0; j < Nmodes; j++)
            ener += 0.5 * double(f1(j)) * positions[j] * positions[j];

        for (int j = 0; j < Nmodes; j++)
            for (int k = 0; k < Nmodes; k++)
                for (int l = 0; l < Nmodes; l++)
                    ener += phi(j, k, l) * positions[j] * positions[k] * positions[l] / 6.0;

        for (int j = 0; j < Nmodes; j++)
            ener += double(diagf(j)) * positions[j] * positions[j] * positions[j] * positions[j] / 24.0;

        if (ener < minEner)
        {
            std::cout << "found new minimum: E = " << ener << " au\n";
            for (int j = 0; j < Nmodes; j++)
                std::cout << "mode " << std::setw(5) << j + 1 << std::setw(10) << positions[j] << "\n";

            minEner = ener;
        }
    }
}
