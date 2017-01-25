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


#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <cmath>
#include <boost/filesystem.hpp>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Eigenvalues>
#include "constants.h"
#include "utilities.h"


/*
 * Call this program like:
 *
 * ./mctdh_specinput <dE> <threshold>
 *
 * dE:        the adiabatic electronic excitation energy in eV
 * threshold: the cutoff threshold for including terms in
 *            the Hamiltonian (optional, default zero)
 *
 * The following files must be present in the working directory:
 *
 * 'gs_freqs'                 the ground and excited state normal mode frequencies in cm-1
 * 'es_freqs'
 * 'Displacement_Vector.dat'  the vector K and the matrix J as output from FCClasses
 * 'Duschinsky_Matrix.dat'
 */
int main(int argc, char *argv[])
{
	/*
	 * Test for the command arguments:
	 */
	if (argc < 2)
	{
		std::cerr << "ERROR: wrong number of arguments\n";
		return 2;
	}
	double dE = atof(argv[1]);

	/*
	 * Test for existence of the required files.
	 */
	std::ifstream gsfc("gs_freqs", std::ifstream::in);
	if (!gsfc.good())
	{
		std::cerr << "ERROR: File 'gs_freqs' could not be opened." << std::endl;
		return 1;
	}
	std::ifstream esfc("es_freqs", std::ifstream::in);
	if (!esfc.good())
	{
		std::cerr << "ERROR: File 'es_freqs' could not be opened." << std::endl;
		return 1;
	}
	std::ifstream shift("Displacement_Vector.dat", std::ifstream::in);
	if (!shift.good())
	{
		std::cerr << "ERROR: File 'Displacement_Vector.dat' could not be opened." << std::endl;
		return 1;
	}
	std::ifstream dusch("Duschinsky_Matrix.dat", std::ifstream::in);
	if (!dusch.good())
	{
		std::cerr << "ERROR: File 'Duschinsky_Matrix.dat' could not be opened." << std::endl;
		return 1;
	}


	/*
	 * Get the number of lines (aka normal modes) in the ground state mode file.
	 * Check, if the excited state mode file contains the same number of modes.
	 */
	int Nmodes = std::count(std::istreambuf_iterator<char>(gsfc), std::istreambuf_iterator<char>(), '\n');
	if (std::count(std::istreambuf_iterator<char>(esfc), std::istreambuf_iterator<char>(), '\n') != Nmodes)
	{
		std::cerr << "ERROR: The files 'gs_freqs' and 'es_freqs' do not contain the same number of lines." << std::endl;
		return 1;
	}
	gsfc.seekg(0);	// jump back to start of file
	esfc.seekg(0);


	/*
	 * Read the ground and excited state frequencies
	 * as well as the Displacement Vector and the Duschinsky Matrix.
	 */
	Eigen::VectorXd v1(Nmodes);
	Eigen::VectorXd v2(Nmodes);
	Eigen::VectorXd Korig(Nmodes), K(Nmodes);
	Eigen::MatrixXd Jorig(Nmodes, Nmodes), J(Nmodes, Nmodes);
	for (int i = 0; i < 2; i++)
		shift.ignore(std::numeric_limits<std::streamsize>::max(), '\n');		// skip the first two lines in the K file
	for (int i = 0; i < 5; i++)
			dusch.ignore(std::numeric_limits<std::streamsize>::max(), '\n');	// and the first five lines in the J file
	for (int i = 0; i < Nmodes; i++)
	{
		int dummyIndex, dummyIndex2;
		gsfc >> v1(i);
		esfc >> v2(i);
		shift >> dummyIndex >> Korig(i);
		for (int j = 0; j < Nmodes; j++)
			dusch >> dummyIndex >> dummyIndex2 >> Jorig(j,i);
	}


	/*
	 * Calculate our J/K quantities from the FCClasses ones.
	 */
	J = Jorig.transpose();
	K = -Jorig.transpose() * Korig;


	/*
	 * Calculate the force constants from the frequencies.
	 */
	Eigen::VectorXd f1(Nmodes);
	Eigen::VectorXd f2(Nmodes);
	/*
	 * conversion from wavenumbers in cm-1 to mass-weighted force constants
	 * in atomic units
	 *
	 * fac = 4 pi**2 * c**2 * (100 cm/m)**2 * me**2 * a0**4 / hbar**2
	 */
	double factor = 40000.0 * M_PI * M_PI * c0 * c0 * me * me * a0 * a0 * a0 * a0 / (hbar * hbar);
	for (int i = 0; i < Nmodes; i++)
	{
		f1(i) = v1(i) * v1(i) * factor;
		f2(i) = v2(i) * v2(i) * factor;
	}


	/*
	 * Calculate the zero-point energies of the two states and print them out.
	 */
	double zpe1 = 0.0;
	double zpe2 = 0.0;
	for (int i = 0; i < Nmodes; i++)
	{
		zpe1 += 0.5 * sqrt(f1(i));
		zpe2 += 0.5 * sqrt(f2(i));
	}


	/*
	 * Open the log file and write some statistics:
	 */
	std::ofstream logFile("log");

	logFile << "Electronic transition energy: " << dE << std::endl;
	logFile << "Number of normal modes: " << Nmodes << std::endl;
	logFile << "Ground state force constants:\n";
    Utils::WriteVectorToFile(logFile, f1);
	logFile << "Excited state force constants:\n";
    Utils::WriteVectorToFile(logFile, f2);
	logFile << "Ground state zero-point Energy: " << zpe1 << "Eh\n";
	logFile << "Excited state zero-point Energy: " << zpe2 << "Eh\n";
	logFile << "Original Displacement Vector from FCClasses:\n";
    Utils::WriteVectorToFile(logFile, Korig);
	logFile << "Original Duschinsky Matrix from FCClasses:\n";
    Utils::WriteMatrixToFile(logFile, Jorig);
	logFile << "Displacement Vector:\n";
    Utils::WriteVectorToFile(logFile, K);
	logFile << "Duschinsky Matrix:\n";
    Utils::WriteMatrixToFile(logFile, J, 3, true);
	logFile << "Metric of the Duschinsky Matrix:\n";
    Utils::WriteMatrixToFile(logFile, J.transpose() * J, 3, true);
	logFile << "Frobenius norm of the Duschinsky matrix: " << J.norm() << std::endl;
	logFile << "Determinant of the Duschinsky matrix: " << J.determinant() << std::endl;

	/*
	 * ############################################################################################
	 *                                  calculate the new PES parameters
	 * ############################################################################################
	 *
	 *
	 * Calculate the new force constants.
	 */
	Eigen::VectorXd fp(Nmodes);
	for (int m = 0; m < Nmodes; m++)
	{
		fp(m) = 0.0;
		for (int n = 0; n < Nmodes; n++)
			fp(m) += f2(n) * J(n,m) * J(n,m);
	}


	/*
	 * Calculate the first-order coefficients.
	 */
	Eigen::VectorXd kappa(Nmodes);
	for (int m = 0; m < Nmodes; m++)
	{
		kappa(m) = 0.0;
		for (int n = 0; n < Nmodes; n++)
			kappa(m) += f2(n) * K(n) * J(n,m);
	}


	/*
	 * Calculate the couplings.
	 */
	Eigen::MatrixXd phi(Eigen::MatrixXd::Zero(Nmodes,Nmodes));
	Eigen::MatrixXd phiFull(Eigen::MatrixXd::Zero(Nmodes,Nmodes));
	for (int m = 0; m < Nmodes; m++)
	{
		for (int o = m + 1; o < Nmodes; o++)
		{
			phi(m,o) = 0.0;
			phiFull(m,o) = 0.0;
			for (int n = 0; n < Nmodes; n++)
				phi(m,o) += f2(n) * J(n,m) * J(n,o);
			phiFull(m,o) = phi(m,o);
			phiFull(o,m) = phi(m,o);
		}
		phiFull(m,m) = fp(m);
	}


	/*
	 * Calculate the energy shifts.
	 */
	Eigen::VectorXd d(Nmodes);
	for (int i = 0; i < Nmodes; i++)
		d(i) = (0.5 * f2(i) * K(i) * K(i));

	/*
	 * write the coupling matrix to the log file:
	 */
	logFile << "Coupling matrix phi with the force constants fp on the diagonal:\n";
    Utils::WriteMatrixToFile(logFile, phiFull, 3, true);


	/*
	 * ############################################################################################
	 *                                  write the MCTDH files
	 * ############################################################################################
	 *
	 *
	 * Now we can finally write the MCTDH input and operator files :)
	 * First, inquire the desired base-name for the files.
	 */
	std::cout << "Enter the base-name for the MCTDH files to be generated.\nThe files <name>.inp and <name>.op will then be written.\n>>> ";
	std::string basename;
	std::cin >> basename;

	std::string inputFileName = basename + ".inp";
	std::string operatorFileName = basename + ".op";

	if (boost::filesystem::exists(inputFileName) || boost::filesystem::exists(operatorFileName))
	{
		std::cout << "One of the MCTDH files already exists. Should they be overwritten? (Y/N)\n>>> ";
		char answer;
		std::cin >> answer;
		if (answer == 'N' || answer == 'n')
			return 0;
	}
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
	// rearrange the modes in order of decreasing coupling
	Eigen::MatrixXd phi_sort = phi.cwiseAbs();
	std::vector<int> sortedModes;
	while (phi_sort.norm() > 0.0)
	{
		Eigen::MatrixXd::Index maxRow, maxCol;
		phi_sort.maxCoeff(&maxRow, &maxCol);
		phi_sort(maxRow, maxCol) = 0.0;

		if (std::find(sortedModes.begin(), sortedModes.end(), maxRow) == sortedModes.end())
			sortedModes.push_back(maxRow);

		if (std::find(sortedModes.begin(), sortedModes.end(), maxCol) == sortedModes.end())
					sortedModes.push_back(maxCol);
	}
	// determine the required number of layers
	int layers = 1;
	while (pow(2.0, layers) < Nmodes)
		layers++;
	layers--;
	// determine the number of nodes in each layer
	std::vector<int> nodesPerLayer(layers);
	nodesPerLayer.at(layers - 1) = Nmodes / 2;
	for (int i = layers - 1; i > 0; i--)
	{
		nodesPerLayer.at(i - 1) = nodesPerLayer.at(i) / 2;
	}
	inputFile << "mlbasis-section\n";
	for (int i = 0; i < Nmodes - 1; i += 2)
	{
		if (sortedModes.size() - i == 3)
			inputFile << "    [q_" << std::setfill('0') << std::setw(3) << sortedModes.at(i) + 1
				      << " q_" << std::setfill('0') << std::setw(3) << sortedModes.at(i+1) + 1
			    	  << " q_" << std::setfill('0') << std::setw(3) << sortedModes.at(i+2) + 1 << "]\n";
		else
			inputFile << "    [q_" << std::setfill('0') << std::setw(3) << sortedModes.at(i) + 1
			          << " q_" << std::setfill('0') << std::setw(3) << sortedModes.at(i+1) + 1 << "]\n";
	}
	inputFile << "end-mlbasis-section\n\n";

	/*
	 * The pbasis-section
	 */
	inputFile << "pbasis-section\n";
	for (int i = 0; i < Nmodes; i++)
		inputFile << "    q_" << std::setfill('0') << std::setw(3) << i + 1
				  << "  ho  " << std::setw(3) << std::setfill(' ') << lrint(-0.8 * log(fp(i))) << "  xi-xf  "
				  //
				  // the basis boundarie are -kappa / fp +- 4 / fp**1/4
				  //
				  << std::fixed << std::setfill(' ') << std::setw(8) << -(kappa(i) / fp(i)) - (2.1 / pow(fp(i), 0.305))
				  << std::fixed << std::setfill(' ') << std::setw(8) << -(kappa(i) / fp(i)) + (2.1 / pow(fp(i), 0.305)) << std::endl;
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
		operatorFile << "    mass_q_" << std::setfill('0') << std::setw(3) << i + 1
					 << "  =  1.0\n";
	// the ground state force constants
	for (int i = 0; i < Nmodes; i++)
	{
		operatorFile << "    f1_" << std::setfill('0') << std::setw(3) << i + 1 << "      = ";
        Utils::WriteFortranNumber(operatorFile, f1(i));
		operatorFile << std::endl;
	}
	// the excited state force constants
	for (int i = 0; i < Nmodes; i++)
	{
		operatorFile << "    f2_" << std::setfill('0') << std::setw(3) << i + 1 << "      = ";
        Utils::WriteFortranNumber(operatorFile, f2(i));
		operatorFile << std::endl;
	}
	// the new effective excited state force constants
	for (int i = 0; i < Nmodes; i++)
	{
		operatorFile << "    fp_" << std::setfill('0') << std::setw(3) << i + 1 << "      = ";
        Utils::WriteFortranNumber(operatorFile, fp(i));
		operatorFile << std::endl;
	}
	// the couplings
	for (int i = 0; i < Nmodes; i++)
		for (int j = i + 1; j < Nmodes; j++)
		{
			operatorFile << "    phi_" << std::setfill('0') << std::setw(3) << i + 1
						 << "_" << std::setfill('0') << std::setw(3) << j + 1 << " = ";
            Utils::WriteFortranNumber(operatorFile, phi(i,j));
			operatorFile << std::endl;
		}
	// the first-order coefficients (shifts)
	for (int i = 0; i < Nmodes; i++)
	{
		operatorFile << "    kappa_" << std::setfill('0') << std::setw(3) << i + 1 << "   = ";
        Utils::WriteFortranNumber(operatorFile, kappa(i));
		operatorFile << std::endl;
	}
	// the energy offsets
	for (int i = 0; i < Nmodes; i++)
	{
		operatorFile << "    d_" << std::setfill('0') << std::setw(3) << i + 1 << "       = ";
        Utils::WriteFortranNumber(operatorFile, d(i));
		operatorFile << std::endl;
	}
	// the electronic offset minus the ground state ZPE
	operatorFile << "    dE          = ";
    Utils::WriteFortranNumber(operatorFile, dE / Eh2eV - zpe1);
	operatorFile << "\nend-parameter-section\n\n";

	/*
	 * The hamiltonian section
	 */
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
	for (int i = 0; i < Nmodes; i++)
			operatorFile << "d_" << std::setfill('0') << std::setw(3) << i + 1
						 << "       |" << i + 1 << " 1\n";
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


	/*
	 * Diagonalize the coupling matrix to get the force constants back
	 */
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> phiSolv(phiFull);
	logFile << "Eigenvalues of the full force constant / coupling matrix:\n";
    Utils::WriteVectorToFile(logFile, phiSolv.eigenvalues());

	/*
	 * Solve the linear system of the full coupling matrix and the kappa vector
	 * to get the coordinates of the minimum
	 */
	Eigen::ColPivHouseholderQR<Eigen::MatrixXd> phiLin(phiFull);
	Eigen::VectorXd minima = phiLin.solve(-kappa);
	logFile << "minimum coordinates\n";
    Utils::WriteVectorToFile(logFile, minima);


	/*
	 * calculate the potential energy at the minimum
	 */
	double Emin = 0.0;

	// first, quadratic term:
	for (int i = 0; i < Nmodes; i++)
		Emin += 0.5 * fp(i) * minima(i) * minima(i);

	// second, coupling term:
	for (int i = 0; i < Nmodes; i++)
		for (int j = i + 1; j < Nmodes; j++)
			Emin += phi(i,j) * minima(i) * minima(j);

	// third, displacement term:
	for (int i = 0; i < Nmodes; i++)
		Emin += kappa(i) * minima(i);

	// fourth, constant term:
	for (int i = 0; i < Nmodes; i++)
		Emin += d(i);

	logFile << "Energy at minimum: " << Emin << std::endl;


	/*
	 * Calculate the 1st moment of the spectrum analytically
	 */
	double moment = dE;
	for (int i = 0; i < Nmodes; i++)
		moment += 0.25 * (fp(i) - f1(i)) / sqrt(f1(i));

	std::cout << "1st moment of the spectrum: " <<std::setprecision(8) << moment << std::endl;

	return 0;
}
