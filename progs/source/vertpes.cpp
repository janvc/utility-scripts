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
#include <fstream>
#include <iomanip>
#include <boost/program_options.hpp>
#include <eigen3/Eigen/Dense>
#include "GaussFchk.h"
#include "vibrationalanalysis.h"
#include "utilities.h"
#include "constants.h"

int main(int argc, char *argv[])
{
    std::string filename_g;
    std::string filename_e;
    std::string basename;
    namespace po = boost::program_options;
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help,h", "produce this help message")
        ("gfile,g", po::value<std::string>(&filename_g)->required(), "the GS Fchk file")
        ("efile,e", po::value<std::string>(&filename_e)->required(), "the ES Fchk file")
        ("mfile,m", po::value<std::string>(&basename)->required(), "the MCTDH file base name")
        ("zpe,z", "include the excited state ZPE in the potential")
        ("noinvert,n", "do not invert the sign of negative force constants")
    ;
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    if (vm.count("help"))
    {
        std::cout << desc << std::endl;
        return 1;
    }
    bool eZPE = false;
    if (vm.count("zpe"))
    {
        eZPE = true;
    }
    bool invFC = true;
    if (vm.count("noinvert"))
    {
        invFC = false;
    }
    po::notify(vm);

    std::ifstream gsStream(filename_g);
    std::ifstream esStream(filename_e);
    GaussFchk gsfile(gsStream);
    GaussFchk esfile(esStream);

    VibrationalAnalysis gState(gsfile);

    /*
     * read the data from the FCHK files
     */
    int Nmodes = gState.Nmodes();
    double GSener = gsfile.ReadReal("Total Energy");
    double ESener = esfile.ReadReal("Total Energy");
    Eigen::VectorXd FCgrad = esfile.ReadVector("Cartesian Gradient");
    Eigen::MatrixXd FChess = esfile.ReadSymmetricMatrix("Cartesian Force Constants");
    double deVert = ESener - GSener;
    Eigen::MatrixXd MassMat = gState.MassVector().asDiagonal();
    for (int i = 0; i < gState.Ncoords(); i++)
        MassMat(i, i) = std::sqrt(MassMat(i, i));

    /*
     * compute the normal mode gradient and hessian
     * (cf. PCCP 14, 13549 (2012); eqns 14, 15)
     */
    Eigen::VectorXd NMgrad = gState.Lmwc().transpose() * MassMat.inverse() * FCgrad;
    Eigen::MatrixXd NMhess = gState.Lmwc().transpose() * MassMat.inverse() * FChess * MassMat.inverse() * gState.Lmwc();

    std::cout << "Ground state frequencies [cm-1]:\n";
    for (int i = 0; i < Nmodes; i++)
        std::cout << std::setw(5) << i+1 << std::fixed << std::setprecision(4) << std::setw(15)
                  << gState.intFrq()(i) << " cm-1\n";

    std::cout << "Excited state effective frequencies [cm-1]:\n";
    for (int i = 0; i < Nmodes; i++)
    {
        if (double(NMhess(i,i) < 0.0))
            std::cout << std::setw(5) << i+1 << std::fixed << std::setprecision(4) << std::setw(15)
                      << -std::sqrt(-double(NMhess(i,i)) * Eh / (a0 * a0 * me)) / (2.0 * M_PI * c0 * 100.0)
                      << " cm-1 (imag)\n";
        else
            std::cout << std::setw(5) << i+1 << std::fixed << std::setprecision(4) << std::setw(15)
                      << std::sqrt(double(NMhess(i,i)) * Eh / (a0 * a0 * me)) / (2.0 * M_PI * c0 * 100.0)
                      << " cm-1\n";

    }

    /*
     * compute the minimum coordinates and the minimum energy
     */
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> linsys(NMhess);
    Eigen::VectorXd Q0 = linsys.solve(NMgrad);
    double Emin = deVert;
    Emin += NMgrad.transpose() * Q0;
    Emin += 0.5 * Q0.transpose() * NMhess * Q0;
    std::cout << "Emin = " << Emin * Eh2eV << " eV\n";

    std::cout << "Minimum coordinates:\n";
    Utils::WriteVector(Q0);

    /*
     * write the MCTDH files
     */
    std::ofstream inputFile(basename + ".inp");
    std::ofstream operatorFile(basename + ".op");
    inputFile.precision(1);
    operatorFile.precision(8);

    inputFile << "run-section\n";
    inputFile << "    name =\n";
    inputFile << "    propagation\n";
    inputFile << "    tfinal =\n";
    inputFile << "    tout =\n";
    inputFile << "    tpsi =\n";
    inputFile << "    psi gridpop auto steps graphviz\n";
    inputFile << "end-run-section\n\n";

    inputFile << "operator-section\n";
    inputFile << "    opname = " << basename << std::endl;
    inputFile << "end-operator-section\n\n";

    inputFile << "mlbasis-section\n";
    for (int i = 0; i < Nmodes; i++)
    {
        inputFile << "    [q_" << std::setfill('0') << std::setw(3) << i + 1 << "]\n";
    }
    inputFile << "end-mlbasis-section\n\n";

    inputFile << "pbasis-section\n";
    for (int i = 0; i < Nmodes; i++)
    {
        double fc = std::abs(double(NMhess(i,i)));

        // determine the number of basis functions:
        int numBas = lrint(-0.7 * log(fc)) + 6;

        // determine the lower and upper bounds of the grid:

//        double lowBound = -double(NMgrad(i)) / fc - (std::abs(double(NMgrad(i) / fc)) + 7.5 / (sqrt(2.0) * pow(fc, 0.25)));
//        double uppBound = -double(NMgrad(i)) / fc + (std::abs(double(NMgrad(i) / fc)) + 7.5 / (sqrt(2.0) * pow(fc, 0.25)));

        double lowBound = double(Q0(i)) - double(Q0(i)) - 5.0 / std::pow(double(gState.intFC()(i)), 0.25);
        double uppBound = double(Q0(i)) + double(Q0(i)) + 5.0 / std::pow(double(gState.intFC()(i)), 0.25);

        inputFile << "    q_" << std::setfill('0') << std::setw(3) << i + 1
                  << "  ho  " << std::setw(3) << std::setfill(' ') << numBas << "  xi-xf  "
                  << std::fixed << std::setfill(' ') << std::setw(8)
                  << lowBound
                  << std::fixed << std::setfill(' ') << std::setw(8)
                  << uppBound << std::endl;
    }
    inputFile << "end-pbasis-section\n\n";

    inputFile << "integrator-section\n";
    inputFile << "    vmf\n";
    inputFile << "    abm = 6, 1.0d-7, 0.01d0\n";
    inputFile << "end-integrator-section\n\n";

    inputFile << "init_wf-section\n";
    inputFile << "    build\n";
    for (int i = 0; i < Nmodes; i++)
    {
        inputFile << "        q_" << std::setfill('0') << std::setw(3) << i + 1
                  << "  eigenf"
                  << "  Eq_" << std::setfill('0') << std::setw(3) << i + 1
                  << "  pop = 1\n";
    }
    inputFile << "    end-build\n";
    inputFile << "end-init_wf-section\n\n";
    inputFile << "end-input\n\n";


    operatorFile << "op_define-section\n";
    operatorFile << "    title\n";
    operatorFile << "        " << basename << std::endl;
    operatorFile << "    end-title\n";
    operatorFile << "end-op_define-section\n\n";

    operatorFile << "parameter-section\n";
    // the masses
    for (int i = 0; i < Nmodes; i++)
    {
        operatorFile << "    mass_q_" << std::setfill('0') << std::setw(3) << i + 1
                     << "  =  1.0\n";
    }
    // the ground state force constants
    for (int i = 0; i < Nmodes; i++)
    {
        operatorFile << "    f1_" << std::setfill('0') << std::setw(3) << i + 1 << "      = ";
        Utils::WriteFortranNumber(operatorFile, double(gState.intFC()(i)));
        operatorFile << std::endl;
    }
    // the new effective excited state force constants
    for (int i = 0; i < Nmodes; i++)
    {
        operatorFile << "    fp_" << std::setfill('0') << std::setw(3) << i + 1 << "      = ";
        if (double(NMhess(i,i)) < 0.0)
        {
            if (invFC)
            {
                Utils::WriteFortranNumber(operatorFile, -double(NMhess(i,i)));
                std::cout << "Warning: I had to invert the sign of fp_" << std::setfill('0') << std::setw(3) << i + 1 << " to make it positive.\n";
            }
            else
            {
                Utils::WriteFortranNumber(operatorFile, double(NMhess(i,i)));
                std::cout << "Warning: fp_" << std::setfill('0') << std::setw(3) << i + 1 << " is negative.\n";
            }
        }
        else
        {
            Utils::WriteFortranNumber(operatorFile, double(NMhess(i,i)));
        }
        operatorFile << std::endl;
    }
    // the couplings
    for (int i = 0; i < Nmodes; i++)
    {
        for (int j = i + 1; j < Nmodes; j++)
        {
            operatorFile << "    phi_" << std::setfill('0') << std::setw(3) << i + 1
                         << "_" << std::setfill('0') << std::setw(3) << j + 1 << " = ";
            Utils::WriteFortranNumber(operatorFile, double(NMhess(i,j)));
            operatorFile << std::endl;
        }
    }
    // the first-order coefficients (shifts)
    for (int i = 0; i < Nmodes; i++)
    {
        operatorFile << "    kappa_" << std::setfill('0') << std::setw(3) << i + 1 << "   = ";
        Utils::WriteFortranNumber(operatorFile, double(NMgrad(i)));
        operatorFile << std::endl;
    }
    // the electronic offset minus the ground state ZPE
    double zpe1 = 0.0;
    if (eZPE)
    {
        for (int i = 0; i < Nmodes; i++)
            zpe1 += 0.5 * sqrt(double(gState.intFC()(i)));
    }
    operatorFile << "    dE          = ";
    Utils::WriteFortranNumber(operatorFile, deVert - zpe1);
    operatorFile << "\nend-parameter-section\n\n";

    operatorFile << "hamiltonian-section";
    for (int i = 0; i < Nmodes; i++)
    {
        if (i % 8 == 0)
            operatorFile << std::endl << "modes";
        operatorFile << " | q_" << std::setfill('0') << std::setw(3) << i + 1;
    }
    operatorFile << std::endl;
    for (int i = 0; i < Nmodes; i++)
        operatorFile << "1.0         |" << i + 1 << " KE\n";
    for (int i = 0; i < Nmodes; i++)
        operatorFile << "0.5*fp_" << std::setfill('0') << std::setw(3) << i + 1
                     << "  |" << i + 1 << " q^2\n";
    for (int i = 0; i < Nmodes; i++)
        for (int j = i + 1; j < Nmodes; j++)
            operatorFile << "phi_" << std::setfill('0') << std::setw(3) << i + 1
                         << "_" << std::setfill('0') << std::setw(3) << j + 1
                         << " |" << i + 1 << " q"
                         << " |" << j + 1 << " q\n";
    for (int i = 0; i < Nmodes; i++)
        operatorFile << "kappa_" << std::setfill('0') << std::setw(3) << i + 1
                     << "   |" << i + 1 << " q\n";
    operatorFile << "dE          |1 1\n";
    operatorFile << "end-hamiltonian-section\n\n";

    /*
     * One-dimensional Hamiltonians for the ground state normal modes
     */
    for (int i = 0; i < Nmodes; i++)
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


