/*
 * Copyright 2017 Jan von Cosel
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
#include <boost/program_options.hpp>
#include "GaussFchk.h"
#include "constants.h"
#include "vibrationalanalysis.h"

/*
 * This program reads structures and normal modes from two Gaussian frequency calculations
 * to create input files for FCClasses. It does so first by reading the normal mode displacements
 * from the formatted checkpoint files and second by generating them from the Hessian
 * through its own normal mode analysis.
 */
int main(int argc, char *argv[])
{
    std::string gsName;
    std::string esName;
    namespace po = boost::program_options;
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help,h", "produce this help message")
        ("gfile,g", po::value<std::string>(&gsName)->required(), "the GS Fchk file")
        ("efile,e", po::value<std::string>(&esName)->required(), "the ES Fchk file")
    ;
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    if (vm.count("help"))
    {
        std::cout << desc << std::endl;
        return 1;
    }
    po::notify(vm);


    std::ifstream gsStream(gsName);
    std::ifstream esStream(esName);
    GaussFchk gsFile(gsStream);
    GaussFchk esFile(esStream);


    // obtain the data from the checkpoint files:
    int Natoms = gsFile.ReadInteger("Number of atoms");
    int Ncoord = 3 * Natoms;
    int Nmodes = Ncoord - 6;
    Eigen::VectorXd geo0 = gsFile.ReadVector("Current cartesian coordinates");
    Eigen::VectorXd geo1 = esFile.ReadVector("Current cartesian coordinates");
    Eigen::MatrixXd vib0 = gsFile.ReadMatrix("Vib-Modes", Nmodes, Ncoord);
    Eigen::MatrixXd vib1 = esFile.ReadMatrix("Vib-Modes", Nmodes, Ncoord);
    Eigen::VectorXd frq0 = gsFile.ReadVector("Vib-E2");
    Eigen::VectorXd frq1 = esFile.ReadVector("Vib-E2");


    std::cout << "===== the ground state file using Gaussian modes starts now =====\n";
    for (int i = 0; i < Ncoord; i++)
        std::cout << std::scientific << std::setw(20) << std::setprecision(10) << double(geo0(i)) * ang2a0 << std::endl;

    for (int i = 0; i < Ncoord; i++)
        for (int j = 0; j < Nmodes; j++)
            std::cout << std::setw(20) << double(vib0(i,j)) << std::endl;

    for (int i = 0; i < Nmodes; i++)
        std::cout << std::setw(20) << double(frq0(i)) << std::endl;

    std::cout << "===== the excited state file using Gaussian modes starts now =====\n";
    for (int i = 0; i < Ncoord; i++)
        std::cout << std::setw(20) << double(geo1(i)) * ang2a0 << std::endl;

    for (int i = 0; i < Ncoord; i++)
        for (int j = 0; j < Nmodes; j++)
            std::cout << std::setw(20) << double(vib1(i,j)) << std::endl;

    for (int i = 0; i < Nmodes; i++)
        std::cout << std::setw(20) << double(frq1(i)) << std::endl;


    // now obtain the same data again, but using our own vibrational analysis:
    VibrationalAnalysis vibAn0(gsFile);
    VibrationalAnalysis vibAn1(esFile);

    std::cout << "===== the ground state file using self-computed modes starts now =====\n";
    for (int i = 0; i < Ncoord; i++)
        std::cout << std::setw(20) << double(vibAn0.X()(i)) * ang2a0 << std::endl;

    for (int i = 0; i < Ncoord; i++)
        for (int j = 0; j < Nmodes; j++)
            std::cout << std::setw(20) << double(vibAn0.Lcart()(i,j)) << std::endl;

    for (int i = 0; i < Nmodes; i++)
        std::cout << std::setw(20) << double(vibAn0.intFrq()(i)) << std::endl;

    std::cout << "===== the excited state file using self-computed modes starts now =====\n";
    for (int i = 0; i < Ncoord; i++)
        std::cout << std::setw(20) << double(vibAn1.X()(i)) * ang2a0 << std::endl;

    for (int i = 0; i < Ncoord; i++)
        for (int j = 0; j < Nmodes; j++)
            std::cout << std::setw(20) << double(vibAn1.Lcart()(i,j)) << std::endl;

    for (int i = 0; i < Nmodes; i++)
        std::cout << std::setw(20) << double(vibAn1.intFrq()(i)) << std::endl;
}
