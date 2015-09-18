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
#include <fstream>
#include <algorithm>
#include <iterator>
#include <cmath>
//#include <boost/system.hpp>
#include <boost/filesystem.hpp>
#include <eigen3/Eigen/Core>
#include "constants.h"


//void writeVector(const Eigen::VectorXd &vector, const size_t dig = 8, const bool clean = false);
void writeVector(const Eigen::VectorXd &vector);
void writeMatrix(const Eigen::MatrixXd &matrix);
void printLayer(const int layer);


int main()
{
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
	gsfc.seekg(0);
	esfc.seekg(0);


	/*
	 * Read the ground and excited state frequencies
	 * as well as the Displacement Vector and the Duschinsky Matrix.
	 */
	Eigen::VectorXd v1(Nmodes);
	Eigen::VectorXd v2(Nmodes);
	Eigen::VectorXd K(Nmodes);
	Eigen::MatrixXd J(Nmodes, Nmodes);
	for (int i = 0; i < 2; i++)
		shift.ignore(std::numeric_limits<std::streamsize>::max(), '\n');		// skip the first two lines in the K file
	for (int i = 0; i < 5; i++)
			dusch.ignore(std::numeric_limits<std::streamsize>::max(), '\n');	// and the first five lines in the J file
	for (int i = 0; i < Nmodes; i++)
	{
		int dummyIndex, dummyIndex2;
		gsfc >> v1(i);
		esfc >> v2(i);
		shift >> dummyIndex >> K(i);
		for (int j = 0; j < Nmodes; j++)
			dusch >> dummyIndex >> dummyIndex2 >> J(j,i);
	}


	/*
	 * Calculate the force constants from the frequencies.
	 */
	Eigen::VectorXd f1(Nmodes);
	Eigen::VectorXd f2(Nmodes);
	for (int i = 0; i < Nmodes; i++)
	{
		double factor = 40000.0 * M_PI * M_PI * c0 * c0 * a0 * a0 * me / Eh;
		f1(i) = v1(i) * v1(i) * factor;
		f2(i) = v2(i) * v2(i) * factor;
	}


	/*
	 * write the values to stdout.
	 */
/*	std::cout << "v1" << std::endl;
	writeVector(v1);
	//for (size_t i = 0; i < v1.size(); i++)
		//std::cout << v1(i) << std::endl;
	std::cout << "v2" << std::endl;
	writeVector(v2);
	std::cout << "f1" << std::endl;
	writeVector(f1);
	std::cout << "f2" << std::endl;
	writeVector(f2);
	std::cout << "K" << std::endl;
	writeVector(K);
	std::cout << "J" << std::endl;
	writeMatrix(J); */


	/*
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
	Eigen::MatrixXd phi(Nmodes,Nmodes);
	for (int m = 0; m < Nmodes; m++)
	{
		for (int o = m + 1; o < Nmodes; o++)
		{
			phi(m,o) = 0.0;
			for (int n = 0; n < Nmodes; n++)
				phi(m,o) += f2(n) * J(n,m) * J(n,o);
		}
	}


	/*
	 * write the new results.
	 */
/*	std::cout << "fp" << std::endl;
	writeVector(fp);
	std::cout << "kappa" << std::endl;
	writeVector(kappa);
	std::cout << "phi" << std::endl;
	writeMatrix(phi); */


	/*
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
	inputFile << "end-run-section\n";
	inputFile << std::endl;


	/*
	 * The operator-section
	 */
	inputFile << "operator-section\n";
	inputFile << "    opname = " << basename << std::endl;
	inputFile << "end-operator-section\n";
	inputFile << std::endl;


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

/*	for (auto ele : sortedModes)
		inputFile << ele << "\n";

	inputFile << "number of layers:" << layers << "\n";
	inputFile << "nodes per layer:\n";
	for (auto ele : nodesPerLayer)
		inputFile << ele << "\n";

	inputFile << "    0>";
	for (int i = 0; i < nodesPerLayer.at(1); i++)
		inputFile << " 2";
	inputFile << std::endl;

	for (int i = 0; i < nodesPerLayer.at(1); i++) */


	inputFile << "end-mlbasis-section\n";

	return 0;
}


/*
 * Write a vector to stdout with 'dig' decimal places (default 8) and setting entries
 * with a magnitude less than a certain threshold to zero if 'clean' is true
 */
//void writeVector(const Eigen::VectorXd &vector, const size_t dig = 8, const bool clean = false)
void writeVector(const Eigen::VectorXd &vector)
{
	//double threshold = vector.norm() / 1.0e-10;

	for (int i = 0; i < vector.size(); i++)
	{
		std::cout << vector(i);
		std::cout << std::endl;
	}
}


void writeMatrix(const Eigen::MatrixXd &matrix)
{
	for (int i = 0; i < matrix.rows(); i++)
	{
		for (int j = 0; j < matrix.cols(); j++)
		{
			std::cout << " " << matrix(i,j);
		}
		std::cout << std::endl;
	}
}

void printLayer(const int layer, const int layers, const std::vector<int> sortedModes, const int lastMode)
{
	if (layer == layers)
	{

	}
	else
		printLayer(layer + 1, layers, sortedModes, lastMode);
}

