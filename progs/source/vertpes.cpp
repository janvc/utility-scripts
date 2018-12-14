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
    int Nprint;
    int Npoints;
    namespace po = boost::program_options;
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help,h", "produce this help message")
        ("gfile,g", po::value<std::string>(&filename_g)->required(), "the GS Fchk file")
        ("efile,e", po::value<std::string>(&filename_e)->required(), "the ES Fchk file")
        ("mfile,m", po::value<std::string>(&basename)->required(), "the MCTDH file base name")
        ("restrict,r", po::value<int>(&Nprint), "number of modes to print to MCTDH files")
        ("points,p", po::value<int>(&Npoints), "number of points to try for minimum search")
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

    std::cout << " Ground state fchk file: " << filename_g << "\n";
    std::cout << "Excited state fchk file: " << filename_e << "\n";
    std::cout << "MCTDH basename:          " << basename << "\n";
    std::cout << "Max. number of modes: ";
    if (vm.count("restrict"))
        std::cout << Nprint;
    else
        std::cout << "no restriction";
    std::cout << "\n\n";

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

    if (! vm.count("restrict"))
    {
        Nprint = Nmodes;
    }

    /*
     * compute the normal mode gradient and hessian
     * (cf. PCCP 14, 13549 (2012); eqns 14, 15)
     */
    Eigen::VectorXd NMgrad = gState.Lmwc().transpose() * MassMat.inverse() * FCgrad;
    Eigen::MatrixXd NMhess = gState.Lmwc().transpose() * MassMat.inverse() * FChess * MassMat.inverse() * gState.Lmwc();

    std::cout << "Ground state mass-weighted force constants [Eh / a0**2 * me]\nand frequencies [cm-1]:\n";
    for (int i = 0; i < Nmodes; i++)
        std::cout << std::setw(5) << i+1
                  << std::scientific << std::setprecision(9) << std::setw(20) << double(gState.intFC()(i))
                  << std::fixed << std::setprecision(4) << std::setw(15) << gState.intFrq()(i) << " cm-1\n";

    double zpeG = 0.0;
    for (int i = 0; i < Nmodes; i++)
        zpeG += 0.5 * std::sqrt(double(gState.intFC()(i)));

    std::cout << "Ground state zero-point-energy:"
              << std::setw(20) << std::scientific << std::setprecision(9) << zpeG << " Eh, "
              << std::setw(10) << std::fixed << std::setprecision(4) << zpeG * Eh2eV << " eV\n\n";

    std::cout << "Vertical excitation energy:    "
              << std::setw(20) << std::scientific << std::setprecision(9) << deVert << " Eh, "
              << std::setw(10) << std::fixed << std::setprecision(4) << deVert * Eh2eV << " eV\n\n";

    std::cout << "Excited state normal mode gradients,  effective force constants and frequencies [cm-1]:\n";
    for (int i = 0; i < Nmodes; i++)
    {
        std::cout << std::setw(5) << i+1
                  << std::scientific << std::setprecision(9) << std::setw(20) << double(NMgrad(i))
                  << std::scientific << std::setprecision(9) << std::setw(20) << double(NMhess(i,i));

        if (double(NMhess(i,i) < 0.0))
            std::cout << std::fixed << std::setprecision(4) << std::setw(15)
                      << -std::sqrt(-double(NMhess(i,i)) * Eh / (a0 * a0 * me)) / (2.0 * M_PI * c0 * 100.0)
                      << " cm-1 (imag)\n";
        else
            std::cout << std::fixed << std::setprecision(4) << std::setw(15)
                      << std::sqrt(double(NMhess(i,i)) * Eh / (a0 * a0 * me)) / (2.0 * M_PI * c0 * 100.0)
                      << " cm-1\n";
    }
    std::cout << "\n";

    /*
     * compute the minimum coordinates and the minimum energy
     */
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> linsys(NMhess);
    Eigen::VectorXd Q0 = linsys.solve(-NMgrad);
    double Emin = deVert;
    Emin += NMgrad.transpose() * Q0;
    Emin += 0.5 * Q0.transpose() * NMhess * Q0;

    std::cout << "Minimum coordinates:\n";
    Utils::WriteVector(Q0);
    std::cout << "Excited state minimum energy:"
              << std::setw(20) << std::scientific << std::setprecision(9) << Emin << " Eh, "
              << std::setw(10) << std::fixed << std::setprecision(4) << Emin * Eh2eV << " eV\n\n";

    /*
     * compute the eigenvalues of the effective force constant matrix
     */
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> estateSolv(NMhess);
    Eigen::VectorXd esFC = estateSolv.eigenvalues();

    std::cout << "Excited state diagonalized force constants and frequencies [cm-1]:\n";
    for (int i = 0; i < Nmodes; i++)
    {
        std::cout << std::setw(5) << i+1
                  << std::scientific << std::setprecision(9) << std::setw(20) << double(esFC(i));

        if (double(esFC(i) < 0.0))
            std::cout << std::fixed << std::setprecision(4) << std::setw(15)
                      << -std::sqrt(-double(esFC(i)) * Eh / (a0 * a0 * me)) / (2.0 * M_PI * c0 * 100.0)
                      << " cm-1 (imag)\n";
        else
            std::cout << std::fixed << std::setprecision(4) << std::setw(15)
                      << std::sqrt(double(esFC(i)) * Eh / (a0 * a0 * me)) / (2.0 * M_PI * c0 * 100.0)
                      << " cm-1\n";
    }

    double zpeE = 0.0;
    for (int i = 0; i < Nmodes; i++)
    {
        if (double(esFC(i) < 0.0))
            zpeE += 0.5 * std::sqrt(-double(esFC(i)));
        else
            zpeE += 0.5 * std::sqrt(double(esFC(i)));
    }

    std::cout << "Excited state zero-point-energy (with negative force constants inverted):"
              << std::setw(20) << std::scientific << std::setprecision(9) << zpeE << " Eh, "
              << std::setw(10) << std::fixed << std::setprecision(4) << zpeE * Eh2eV << " eV\n\n";

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
    for (int i = 0; i < Nprint; i++)
    {
        inputFile << "    [q_" << std::setfill('0') << std::setw(3) << i + 1 << "]\n";
    }
    inputFile << "end-mlbasis-section\n\n";

    inputFile << "pbasis-section\n";
    std::vector<double> lowBounds(Nprint);
    std::vector<double> uppBounds(Nprint);
    for (int i = 0; i < Nprint; i++)
    {
        double fc = std::abs(double(NMhess(i,i)));

        // determine the number of basis functions:
        int numBas = 2 * (lrint(-0.7 * log(fc)) + 6);

        // determine the lower and upper bounds of the grid:

//        double lowBound = -double(NMgrad(i)) / fc - (std::abs(double(NMgrad(i) / fc)) + 7.5 / (sqrt(2.0) * pow(fc, 0.25)));
//        double uppBound = -double(NMgrad(i)) / fc + (std::abs(double(NMgrad(i) / fc)) + 7.5 / (sqrt(2.0) * pow(fc, 0.25)));

        lowBounds[i] = double(Q0(i)) - std::abs(double(Q0(i))) - 7.5 / std::pow(double(gState.intFC()(i)), 0.25);
        uppBounds[i] = double(Q0(i)) + std::abs(double(Q0(i))) + 7.5 / std::pow(double(gState.intFC()(i)), 0.25);

        inputFile << "    q_" << std::setfill('0') << std::setw(3) << i + 1
                  << "  ho  " << std::setw(3) << std::setfill(' ') << numBas << "  xi-xf  "
                  << std::fixed << std::setfill(' ') << std::setw(8)
                  << lowBounds[i]
                  << std::fixed << std::setfill(' ') << std::setw(8)
                  << uppBounds[i] << std::endl;
    }
    inputFile << "end-pbasis-section\n\n";

    inputFile << "integrator-section\n";
    inputFile << "    vmf\n";
    inputFile << "    abm = 6, 1.0d-7, 0.01d0\n";
    inputFile << "end-integrator-section\n\n";

    inputFile << "init_wf-section\n";
    inputFile << "    build\n";
    for (int i = 0; i < Nprint; i++)
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
    for (int i = 0; i < Nprint; i++)
    {
        operatorFile << "    mass_q_" << std::setfill('0') << std::setw(3) << i + 1
                     << "  =  1.0\n";
    }
    // the ground state force constants
    for (int i = 0; i < Nprint; i++)
    {
        operatorFile << "    f1_" << std::setfill('0') << std::setw(3) << i + 1 << "      = ";
        Utils::WriteFortranNumber(operatorFile, double(gState.intFC()(i)));
        operatorFile << std::endl;
    }
    // the new effective excited state force constants
    for (int i = 0; i < Nprint; i++)
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
    for (int i = 0; i < Nprint; i++)
    {
        for (int j = i + 1; j < Nprint; j++)
        {
            operatorFile << "    phi_" << std::setfill('0') << std::setw(3) << i + 1
                         << "_" << std::setfill('0') << std::setw(3) << j + 1 << " = ";
            Utils::WriteFortranNumber(operatorFile, double(NMhess(i,j)));
            operatorFile << std::endl;
        }
    }
    // the first-order coefficients (shifts)
    for (int i = 0; i < Nprint; i++)
    {
        operatorFile << "    kappa_" << std::setfill('0') << std::setw(3) << i + 1 << "   = ";
        Utils::WriteFortranNumber(operatorFile, double(NMgrad(i)));
        operatorFile << std::endl;
    }
    // the electronic offset minus the ground state ZPE
    operatorFile << "    dE          = ";
    if (eZPE)
    {
        Utils::WriteFortranNumber(operatorFile, deVert - zpeG);
    }
    else
        Utils::WriteFortranNumber(operatorFile, deVert);
    operatorFile << "\nend-parameter-section\n\n";

    operatorFile << "hamiltonian-section";
    for (int i = 0; i < Nprint; i++)
    {
        if (i % 8 == 0)
            operatorFile << std::endl << "modes";
        operatorFile << " | q_" << std::setfill('0') << std::setw(3) << i + 1;
    }
    operatorFile << std::endl;
    for (int i = 0; i < Nprint; i++)
        operatorFile << "1.0         |" << i + 1 << " KE\n";
    for (int i = 0; i < Nprint; i++)
        operatorFile << "0.5*fp_" << std::setfill('0') << std::setw(3) << i + 1
                     << "  |" << i + 1 << " q^2\n";
    for (int i = 0; i < Nprint; i++)
        for (int j = i + 1; j < Nprint; j++)
            operatorFile << "phi_" << std::setfill('0') << std::setw(3) << i + 1
                         << "_" << std::setfill('0') << std::setw(3) << j + 1
                         << " |" << i + 1 << " q"
                         << " |" << j + 1 << " q\n";
    for (int i = 0; i < Nprint; i++)
        operatorFile << "kappa_" << std::setfill('0') << std::setw(3) << i + 1
                     << "   |" << i + 1 << " q\n";
    operatorFile << "dE          |1 1\n";
    operatorFile << "end-hamiltonian-section\n\n";

    /*
     * One-dimensional Hamiltonians for the ground state normal modes
     */
    for (int i = 0; i < Nprint; i++)
    {
        operatorFile << "hamiltonian-section_Eq_" << std::setfill('0') << std::setw(3) << i + 1 << std::endl;
        operatorFile << "usediag\n";
        operatorFile << "modes      | q_" << std::setfill('0') << std::setw(3) << i + 1 << std::endl;
        operatorFile << "1.0        |1 KE\n";
        operatorFile << "0.5*f1_" << std::setfill('0') << std::setw(3) << i + 1 << " |1 q^2\n";
        operatorFile << "end-hamiltonian-section\n\n";
    }
    operatorFile << "end-operator\n";


    /*
     * look for minima in the PES
     */
    if (vm.count("points"))
    {
        double minEner = Emin;
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0.0, 1.0);

        for (int i = 0; i < Npoints; i++)
        {
            std::vector<double> positions(Nprint);

            for (int j = 0; j < Nprint; j++)
                positions[j] = lowBounds[j] + dis(gen) * (uppBounds[j] - lowBounds[j]);

            double ener = 0.0;
            for (int j = 0; j < Nprint; j++)
            {
                if (double(NMhess(j,j)) < 0.0)
                    ener -= 0.5 * double(NMhess(j,j)) * positions[j] * positions[j];
                else
                    ener += 0.5 * double(NMhess(j,j)) * positions[j] * positions[j];
            }

            for (int j = 0; j < Nprint; j++)
            {
                for (int k = 0; k < Nprint; k++)
                {
                    ener += double(NMhess(j,k)) * positions[j] * positions[k];
                }
            }

            for (int j = 0; j < Nprint; j++)
                ener += double(NMgrad(j)) * positions[j];

            ener += Emin;

//            std::cout << "point " << i << " ener " << ener << "\n";
//            for (int j = 0; j < Nprint; j++)
//                std::cout << "mode " << j + 1 << " " << positions[j] << "\n";

            if (ener < minEner)
            {
                std::cout << "found new minimum: E = " << ener * Eh2eV << " eV\n";
                for (int j = 0; j < Nprint; j++)
                    std::cout << "mode " << std::setw(5) << j + 1 << std::setw(10) << positions[j] << "\n";

                minEner = ener;
            }
        }
    }

    return 0;
}


