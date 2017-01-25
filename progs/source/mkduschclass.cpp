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


#include <iomanip>
#include <eigen3/Eigen/SVD>
#include <eigen3/Eigen/Eigenvalues>
#include "mkduschclass.h"
#include "utilities.h"

MkDuschClass::MkDuschClass(GaussFchk &gState, GaussFchk &eState)
    : m_GroundState(gState)
    , m_ExcitedState(eState)
{
    // read basic stuff from the Fchk files:
    m_Natoms = gState.ReadInteger("Number of atoms");
    m_Ncoords = 3 * m_Natoms;
    m_Nmodes = m_Ncoords - 6;

    m_gsEner = gState.ReadReal("Total Energy");
    m_esEner = eState.ReadReal("Total Energy");
    m_dE = m_esEner - m_gsEner;

    // calculate the total zero-point energies
    // (including all modes)
    m_zpeGtot = 0.0;
    m_zpeEtot = 0.0;
    for (int i = 0; i < m_Nmodes; i++)
    {
        m_zpeGtot += 0.5 * sqrt(double(m_GroundState.intFC()(i)));
        m_zpeEtot += 0.5 * sqrt(double(m_ExcitedState.intFC()(i)));
    }

    // perform SVD on the normal modes:
    Eigen::JacobiSVD<Eigen::MatrixXd> SVDg(m_GroundState.Lmwc(), Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::JacobiSVD<Eigen::MatrixXd> SVDe(m_ExcitedState.Lmwc(), Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::MatrixXd NMOg = SVDg.matrixU() * SVDg.matrixV().transpose();
    Eigen::MatrixXd NMOe = SVDe.matrixU() * SVDe.matrixV().transpose();

    // minimize the RMSD between the structures
    m_XSg = Eigen::VectorXd(m_Ncoords);
    m_XSe = Eigen::VectorXd(m_Ncoords);
    for (int i = 0; i < m_Natoms; i++)
        for (int j = 0; j < 3; j++)
        {
            m_XSg(3*i+j) = m_GroundState.X()(3*i+j) - m_GroundState.com()(j);
            m_XSe(3*i+j) = m_ExcitedState.X()(3*i+j) - m_ExcitedState.com()(j);
        }
    m_QSg = Eigen::VectorXd(m_Ncoords);
    m_QSe = Eigen::VectorXd(m_Ncoords);
    for (int i = 0; i < m_Natoms; i++)
    {
        m_QSg(3*i+0) = m_XSg(3*i+0) * sqrt(double(m_GroundState.masses()(i)));
        m_QSe(3*i+0) = m_XSe(3*i+0) * sqrt(double(m_GroundState.masses()(i)));
        m_QSg(3*i+1) = m_XSg(3*i+1) * sqrt(double(m_GroundState.masses()(i)));
        m_QSe(3*i+1) = m_XSe(3*i+1) * sqrt(double(m_GroundState.masses()(i)));
        m_QSg(3*i+2) = m_XSg(3*i+2) * sqrt(double(m_GroundState.masses()(i)));
        m_QSe(3*i+2) = m_XSe(3*i+2) * sqrt(double(m_GroundState.masses()(i)));
    }

    m_COMg = Utils::calc_com(m_XSg, m_GroundState.masses());
    m_COMe = Utils::calc_com(m_XSe, m_GroundState.masses());

    m_RMSDb = 0.0;
    for (int i = 0; i < m_Ncoords; i++)
        m_RMSDb += std::abs(double(m_XSg(i) - m_XSe(i))) * std::abs(double(m_XSg(i) - m_XSe(i)));
    m_RMSDb /= m_Natoms;

    m_CorrMat = Eigen::Matrix3d::Zero();
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            for (int k = 0; k < m_Natoms; k++)
                m_CorrMat(i,j) += m_XSg(3*k+i) * m_XSe(3*k+j);

    m_F(0,0) =  m_CorrMat(0,0) + m_CorrMat(1,1) + m_CorrMat(2,2);
    m_F(1,1) =  m_CorrMat(0,0) - m_CorrMat(1,1) - m_CorrMat(2,2);
    m_F(2,2) = -m_CorrMat(0,0) + m_CorrMat(1,1) - m_CorrMat(2,2);
    m_F(3,3) = -m_CorrMat(0,0) - m_CorrMat(1,1) + m_CorrMat(2,2);
    m_F(0,1) =  m_CorrMat(1,2) - m_CorrMat(2,1);
    m_F(0,2) =  m_CorrMat(2,0) - m_CorrMat(0,2);
    m_F(0,3) =  m_CorrMat(0,1) - m_CorrMat(1,0);
    m_F(1,2) =  m_CorrMat(0,1) + m_CorrMat(1,0);
    m_F(1,3) =  m_CorrMat(0,2) + m_CorrMat(2,0);
    m_F(2,3) =  m_CorrMat(1,2) + m_CorrMat(2,1);
    m_F(1,0) =  m_F(0,1);
    m_F(2,0) =  m_F(0,2);
    m_F(3,0) =  m_F(0,3);
    m_F(2,1) =  m_F(1,2);
    m_F(3,1) =  m_F(1,3);
    m_F(3,2) =  m_F(2,3);

    // diagonalize the quaternion matrix and select the best rotation quaternion:
    // (the one corresponding to the largest absolute eigenvalue)
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix4d> Feig(m_F);

    Eigen::Vector4d lQuart = abs(double(Feig.eigenvalues()(0))) > abs(double(Feig.eigenvalues()(3)))
                            ? Feig.eigenvectors().block(0, 0, 4, 1) : Feig.eigenvectors().block(0, 3, 4, 1);

    m_RotMat(0,0) = lQuart(0) * lQuart(0) + lQuart(1) * lQuart(1) + lQuart(2) * lQuart(2) + lQuart(3) * lQuart(3);
    m_RotMat(1,1) = lQuart(0) * lQuart(0) - lQuart(1) * lQuart(1) + lQuart(2) * lQuart(2) - lQuart(3) * lQuart(3);
    m_RotMat(2,2) = lQuart(0) * lQuart(0) - lQuart(1) * lQuart(1) - lQuart(2) * lQuart(2) + lQuart(3) * lQuart(3);
    m_RotMat(0,1) = 2.0 * (lQuart(1) * lQuart(2) - lQuart(0) * lQuart(3));
    m_RotMat(1,0) = 2.0 * (lQuart(1) * lQuart(2) + lQuart(0) * lQuart(3));
    m_RotMat(0,2) = 2.0 * (lQuart(1) * lQuart(3) + lQuart(0) * lQuart(2));
    m_RotMat(2,0) = 2.0 * (lQuart(1) * lQuart(3) - lQuart(0) * lQuart(2));
    m_RotMat(1,2) = 2.0 * (lQuart(2) * lQuart(3) - lQuart(0) * lQuart(1));
    m_RotMat(2,1) = 2.0 * (lQuart(2) * lQuart(3) + lQuart(0) * lQuart(1));

    Eigen::MatrixXd BigRotMat(Eigen::MatrixXd::Zero(m_Ncoords, m_Ncoords));
    for (int i = 0; i < m_Natoms; i++)
        BigRotMat.block(3*i, 3*i, 3, 3) = m_RotMat;

    m_XRg = BigRotMat * m_XSg;
    m_QRg = BigRotMat * m_QSg;

    m_RMSDa = 0.0;
    for (int i = 0; i < m_Ncoords; i++)
        m_RMSDa += std::abs(double(m_XRg(i) - m_XSe(i))) * std::abs(double(m_XRg(i) - m_XSe(i)));
    m_RMSDa /= m_Natoms;

    if (m_RMSDa < m_RMSDb)
    {
        m_Jfull = NMOe.transpose() * BigRotMat * NMOg;
        m_Kfull = NMOe.transpose() * (m_QRg - m_QSe);
    }
    else
    {
        m_Jfull = NMOe.transpose() * NMOg;
        m_Kfull = NMOe.transpose() * (m_QSg - m_QSe);
    }
}

void MkDuschClass::calcMCTDHdata()
{
    m_f1 = m_GroundState.intFC();
    m_f2 = m_ExcitedState.intFC();

    m_fp = Eigen::VectorXd::Zero(m_Nmodes);
    for (int m = 0; m < m_Nmodes; m++)
        for (int n = 0; n < m_Nmodes; n++)
            m_fp(m) += m_f2(n) * m_Jfull(n,m) * m_Jfull(n,m);

    m_kappa = Eigen::VectorXd::Zero(m_Nmodes);
    for (int m = 0; m < m_Nmodes; m++)
        for (int n = 0; n < m_Nmodes; n++)
            m_kappa(m) += m_f2(n) * m_Kfull(n) * m_Jfull(n,m);

    m_phi = Eigen::MatrixXd::Zero(m_Nmodes, m_Nmodes);
    for (int m = 0; m < m_Nmodes; m++)
    {
        for (int o = 0; o < m_Nmodes; o++)
        {
            for (int n = 0; n < m_Nmodes; n++)
                m_phi(m,o) += m_f2(n) * m_Jfull(n,m) * m_Jfull(n,o);
            m_phi(o,m) = m_phi(m,o);
        }
        m_phi(m,m) = m_fp(m);
    }

    m_d = Eigen::VectorXd::Zero(m_Nmodes);
    for (int i = 0; i < m_Nmodes; i++)
        m_d(i) = 0.5 * m_f2(i) * m_Kfull(i) * m_Kfull(i);

    m_fp_fw = Eigen::VectorXd::Zero(m_Nmodes);
    for (int i = 0; i < m_Nmodes; i++)
        m_fp_fw(i) = double(m_fp(i)) / sqrt(double(m_f1(i)));

    m_phi_fw = Eigen::MatrixXd::Zero(m_Nmodes, m_Nmodes);
    for (int i = 0; i < m_Nmodes; i++)
        for (int j = 0; j < m_Nmodes; j++)
            m_phi_fw(i,j) = double(m_phi(i,j)) / (pow(double(m_f1(i)), 0.25) * pow(double(m_f1(j)), 0.25));

    m_kappa_fw = Eigen::VectorXd::Zero(m_Nmodes);
    for (int i = 0; i < m_Nmodes; i++)
        m_kappa_fw(i) = double(m_kappa(i)) / pow(double(m_f1(i)), 0.25);
}

void MkDuschClass::selectModes()
{
    /*
     *  sort the kappa values to get the most strongly shifted modes:
     */
    std::vector<double> kappas(m_kappa_fw.data(), m_kappa_fw.data() + m_kappa_fw.size());
    std::vector<int> KappaInd(kappas.size());
    std::iota(KappaInd.begin(), KappaInd.end(), 0);
    // use a lambda expression to sort:
    std::sort(KappaInd.begin(), KappaInd.end(), [&kappas](int i1, int i2){return std::abs(kappas[i1]) > std::abs(kappas[i2]);});

    Utils::WriteVector(m_kappa_fw);

    for (int i = 0; i < KappaInd.size(); i++)
        std::cout << KappaInd.at(i) + 1 << std::endl;

    /*
     * find the strongest couplings by taking the largest value of every column
     * of the coupling matrix and sorting them using the same scheme as above:
     */
    Eigen::MatrixXd phiAbs = m_phi_fw.cwiseAbs();
    // set the diagonal elements (frequencies) to zero:
    for (int i = 0; i < m_Nmodes; i++)
        phiAbs(i,i) = 0.0;
    std::vector<double> couplings;
    for (int i = 0; i < m_Nmodes; i++)
    {
        Eigen::VectorXd tmpVec = phiAbs.col(i);
        couplings.push_back(tmpVec.maxCoeff());
    }

    for (int i = 0; i < couplings.size(); i++)
        std::cout << i + 1 << "  " << couplings.at(i) << std::endl;

    std::vector<int> PhiInd(couplings.size());
    std::iota(PhiInd.begin(), PhiInd.end(), 0);
    std::sort(PhiInd.begin(), PhiInd.end(), [&couplings](int i1, int i2){return std::abs(couplings[i1]) > std::abs(couplings[i2]);});

    for (int i = 0; i < PhiInd.size(); i++)
        std::cout << PhiInd.at(i) + 1 << std::endl;

    m_Nactive = 10;
    std::vector<int> activeModes;
    std::vector<int> kappaCandidates(KappaInd.begin(), KappaInd.begin() + m_Nactive);
    std::vector<int> phiCandidates(PhiInd.begin(), PhiInd.begin() + m_Nactive);

    for (int i = 0; i < m_Nactive; i++)
    {
        std::cout << kappaCandidates.at(i) + 1 << "   " << kappas.at(kappaCandidates.at(i)) << "   "
                  << phiCandidates.at(i) + 1 << "   " << couplings.at(phiCandidates.at(i)) << std::endl;

        if (std::find(activeModes.begin(), activeModes.end(), kappaCandidates.at(i)) == activeModes.end())
            activeModes.push_back(kappaCandidates.at(i));

        if (std::find(activeModes.begin(), activeModes.end(), phiCandidates.at(i)) == activeModes.end())
            activeModes.push_back(phiCandidates.at(i));
    }

    for (int i = 0; i < m_Nactive; i++)
        std::cout << activeModes.at(i) + 1 << std::endl;
}

void MkDuschClass::createMCTDHfiles(const std::string &basename, const bool freqWeight, const bool noJ, const bool noK) const
{
    std::string inputFileName = basename + ".inp";
    std::string operatorFileName = basename + ".op";

    std::ofstream inputFile(inputFileName);
    std::ofstream operatorFile(operatorFileName);

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
    for (int i = 0; i < m_Nmodes - 1; i+= 2)
    {
        if (m_Nmodes - i == 3)
            inputFile << "    [q_" << std::setfill('0') << std::setw(3) << i + 1
                      << " q_" << std::setfill('0') << std::setw(3) << i + 2
                      << " q_" << std::setfill('0') << std::setw(3) << i + 3 << "]\n";
        else
            inputFile << "    [q_" << std::setfill('0') << std::setw(3) << i + 1
                      << " q_" << std::setfill('0') << std::setw(3) << i + 2 << "]\n";
    }
    inputFile << "end-mlbasis-section\n\n";

    inputFile << "pbasis-section\n";
    for (int i = 0; i < m_Nmodes; i++)
    {
        // determine the number of basis functions:
        int numFuncs;
        if (freqWeight)
            numFuncs = lrint(-0.7 * log(double(m_fp_fw(i)))) + 8;
        else
            numFuncs = lrint(-0.7 * log(double(m_fp(i)))) + 6;

        // determine the lower and upper bounds of the grid:
        // the basis boundarie are (originally) -kappa / fp +- 6.5 / fp**1/4
        double lowBound;
        double uppBound;
        if (freqWeight)
        {
            lowBound = -double(m_kappa_fw(i) / m_fp_fw(i))
//                     - (std::abs(double(m_kappa_fw(i) / m_fp_fw(i))) + 1.4 / (1.0 * pow(double(m_fp_fw(i)), 0.16)));
                     - (std::abs(double(m_kappa_fw(i) / m_fp_fw(i))) + 4.0);
            uppBound = -double(m_kappa_fw(i) / m_fp_fw(i))
//                     + (std::abs(double(m_kappa_fw(i) / m_fp_fw(i))) + 1.4 / (1.0 * pow(double(m_fp_fw(i)), 0.16)));
                     + (std::abs(double(m_kappa_fw(i) / m_fp_fw(i))) + 4.0);
        }
        else
        {
            lowBound = -double(m_kappa(i) / m_fp(i))
                     - (std::abs(double(m_kappa(i) / m_fp(i))) + 7.5 / (sqrt(2.0) * pow(double(m_fp(i)), 0.25)));
            uppBound = -double(m_kappa(i) / m_fp(i))
                     + (std::abs(double(m_kappa(i) / m_fp(i))) + 7.5 / (sqrt(2.0) * pow(double(m_fp(i)), 0.25)));
        }

        inputFile << "    q_" << std::setfill('0') << std::setw(3) << i + 1
                  << "  ho  " << std::setw(3) << std::setfill(' ') << numFuncs << "  xi-xf  "
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
    for (int i = 0; i < m_Nmodes; i++)
        inputFile << "        q_" << std::setfill('0') << std::setw(3) << i + 1
                  << "  eigenf"
                  << "  Eq_" << std::setfill('0') << std::setw(3) << i + 1
                  << "  pop = 1\n";
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
    for (int i = 0; i < m_Nmodes; i++)
        operatorFile << "    mass_q_" << std::setfill('0') << std::setw(3) << i + 1
                     << "  =  1.0\n";
    // the ground state force constants
    for (int i = 0; i < m_Nmodes; i++)
        if (freqWeight)
        {
            operatorFile << "    w1_" << std::setfill('0') << std::setw(3) << i + 1 << "      = ";
            Utils::WriteFortranNumber(operatorFile, std::sqrt(double(m_f1(i))));
            operatorFile << std::endl;
        }
        else
        {
            operatorFile << "    f1_" << std::setfill('0') << std::setw(3) << i + 1 << "      = ";
            Utils::WriteFortranNumber(operatorFile, double(m_f1(i)));
            operatorFile << std::endl;
        }
    // the new effective excited state force constants
    for (int i = 0; i < m_Nmodes; i++)
        if (freqWeight)
        {
            operatorFile << "    fp_" << std::setfill('0') << std::setw(3) << i + 1 << "      = ";
            Utils::WriteFortranNumber(operatorFile, double(m_fp_fw(i)));
            operatorFile << std::endl;
        }
        else
        {
            operatorFile << "    fp_" << std::setfill('0') << std::setw(3) << i + 1 << "      = ";
            Utils::WriteFortranNumber(operatorFile, double(m_fp(i)));
            operatorFile << std::endl;
        }
    // the couplings
    for (int i = 0; i < m_Nmodes; i++)
        for (int j = i + 1; j < m_Nmodes; j++)
            if (freqWeight)
            {
                operatorFile << "    phi_" << std::setfill('0') << std::setw(3) << i + 1
                        << "_" << std::setfill('0') << std::setw(3) << j + 1 << " = ";
                Utils::WriteFortranNumber(operatorFile, double(m_phi_fw(i,j)));
                operatorFile << std::endl;
            }
            else
            {
                operatorFile << "    phi_" << std::setfill('0') << std::setw(3) << i + 1
                             << "_" << std::setfill('0') << std::setw(3) << j + 1 << " = ";
                Utils::WriteFortranNumber(operatorFile, double(m_phi(i,j)));
                operatorFile << std::endl;
            }
    // the first-order coefficients (shifts)
    for (int i = 0; i < m_Nmodes; i++)
        if (freqWeight)
        {
            operatorFile << "    kappa_" << std::setfill('0') << std::setw(3) << i + 1 << "   = ";
            Utils::WriteFortranNumber(operatorFile, double(m_kappa_fw(i)));
            operatorFile << std::endl;
        }
        else
        {
            operatorFile << "    kappa_" << std::setfill('0') << std::setw(3) << i + 1 << "   = ";
            Utils::WriteFortranNumber(operatorFile, double(m_kappa(i)));
            operatorFile << std::endl;
        }
    // the energy offsets
    for (int i = 0; i < m_Nmodes; i++)
    {
        operatorFile << "    d_" << std::setfill('0') << std::setw(3) << i + 1 << "       = ";
        Utils::WriteFortranNumber(operatorFile, double(m_d(i)));
        operatorFile << std::endl;
    }
    // the electronic offset minus the ground state ZPE
//    double zpe1 = 0.0;
//    for (int i = 0; i < m_Nmodes; i++)
//        if (isPresent.at(i))
//            zpe1 += 0.5 * sqrt(double(f1(i)));
    operatorFile << "    dE          = ";
    Utils::WriteFortranNumber(operatorFile, m_dE - m_zpeGtot);
    operatorFile << "\nend-parameter-section\n\n";


    operatorFile << "hamiltonian-section";
    for (int i = 0; i < m_Nmodes; i++)
    {
        if (i % 8 == 0)
            operatorFile << std::endl << "modes";
        operatorFile << " | q_" << std::setfill('0') << std::setw(3) << i + 1;
    }
    operatorFile << std::endl;
    for (int i = 0; i < m_Nmodes; i++)
        if (freqWeight)
            operatorFile << "w1_" << std::setfill('0') << std::setw(3) << i + 1 << "      |" << i + 1 << " KE\n";
        else
            operatorFile << "1.0         |" << i + 1 << " KE\n";
    for (int i = 0; i < m_Nmodes; i++)
            operatorFile << "0.5*fp_" << std::setfill('0') << std::setw(3) << i + 1
                     << "  |" << i + 1 << " q^2\n";
    for (int i = 0; i < m_Nmodes; i++)
        for (int j = i + 1; j < m_Nmodes; j++)
                operatorFile << "phi_" << std::setfill('0') << std::setw(3) << i + 1
                             << "_" << std::setfill('0') << std::setw(3) << j + 1
                             << " |" << i + 1 << " q"
                             << " |" << j + 1 << " q\n";
    for (int i = 0; i < m_Nmodes; i++)
            operatorFile << "kappa_" << std::setfill('0') << std::setw(3) << i + 1
                         << "   |" << i + 1 << " q\n";
    for (int i = 0; i < m_Nmodes; i++)
            operatorFile << "d_" << std::setfill('0') << std::setw(3) << i + 1
                         << "       |" << i + 1 << " 1\n";
    operatorFile << "dE          |1 1\n";
    operatorFile << "end-hamiltonian-section\n\n";

    for (int i = 0; i < m_Nmodes; i++)
    {
        operatorFile << "hamiltonian-section_Eq_" << std::setfill('0') << std::setw(3) << i + 1 << std::endl;
        operatorFile << "usediag\n";
        operatorFile << "modes       | q_" << std::setfill('0') << std::setw(3) << i + 1 << std::endl;
        if (freqWeight)
        {
            operatorFile << "w1_" << std::setfill('0') << std::setw(3) << i + 1 << "      |1 KE\n";
            operatorFile << "0.5*w1_" << std::setfill('0') << std::setw(3) << i + 1 << "  |1 q^2\n";
        }
        else
        {
            operatorFile << "1.0        |1 KE\n";
            operatorFile << "0.5*f1_" << std::setfill('0') << std::setw(3) << i + 1 << "  |1 q^2\n";
        }
        operatorFile << "end-hamiltonian-section\n\n";
    }
    operatorFile << "end-operator\n";
}

