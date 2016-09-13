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
#include <fstream>
#include <iomanip>
#include <complex>
#include <boost/program_options.hpp>
#include "constants.h"

int main(int argc, char *argv[])
{
	std::string fileName;
	double tau;
	namespace po = boost::program_options;
	po::options_description desc("Allowed options");
	desc.add_options()
					("help,h", "produce this help message")
					("file,f", po::value<std::string>(&fileName)->required(), "the autocorrelation file")
					("tau,t", po::value<double>(&tau)->required(), "the damping time")
					("phase,p", "set the initial phase to zero")
					;
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	if (vm.count("help"))
	{
		std::cout << desc << std::endl;
		return 1;
	}
	po::notify(vm);

	// convert the damping time to atomic units:
	tau *= fs2au;

	std::ifstream inputFile(fileName, std::ifstream::in);
	std::vector<std::complex<double> > autoInp;
	std::vector<double> time;

	inputFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

	double dummy1, dummy2, dummy3, dummy4;
	while (inputFile >> dummy1 >> dummy2 >> dummy3 >> dummy4)
	{
		time.push_back(dummy1);
		std::complex<double> dummy_comp;
		dummy_comp.real(dummy2);
		dummy_comp.imag(dummy3);
		autoInp.push_back(dummy_comp);
	}
	inputFile.close();


	// calculate the required times and frequencies:
	int N = int(time.size());
	double dt = (time[1] - time[0]) * fs2au;
	double tmax = dt * double(N - 1);

	std::vector<double> doubleTime;
	for (int i = 0; i < 2 * N - 1; i++)
		doubleTime.push_back(-tmax + double(i) * dt);


	// set the initial phase of the autocorrelation function to zero and re-normalize it:
	if (vm.count("phase"))
	{
		double phase = -std::atan2(autoInp[0].imag(), autoInp[0].real());
		double norm = 1.0 / std::abs(autoInp[0]);
		std::cout << "applying a phase of " << phase << " to the autocorrelation function\n";
		for (int i = 0; i < N; i++)
			autoInp[i] *= std::exp(std::complex<double>(0.0, 1.0) * phase) * norm;
	}


	// construct the mirrored autocorrelation function with cosine weighting and exponential damping:
	std::vector<std::complex<double> > autoDouble(2 * N - 1);
	for (int i = 0; i < N; i++)
	{
		autoDouble.at(i) = std::conj(autoInp.at(N - i - 1));
		autoDouble.at(2 * N - 2 - i) = autoInp.at(N - i - 1);
	}

	for (int i = 0; i < 2 * N - 1; i++)
		autoDouble[i] *= std::cos(M_PI * doubleTime[i] / (2.0 * tmax)) * std::exp(-std::abs(doubleTime[i] / tau));

	// perform discrete fourier transformation:
	// (Numerical Recipes p. 497, eq 12.1.7)
	int shift =-(2 * N - 1) / 2;
	std::vector<double> freq(2 * N - 1);
	std::vector<std::complex<double> > spec(2 * N - 1);
	for (int m = shift; m < 2 * N - 1 + shift; m++)
	{
		freq[m - shift] = double(m) / (double(2 * N - 1) * dt);
		spec[m - shift] = std::complex<double>(0.0, 0.0);
		for (int k = 0; k < 2 * N - 1; k++)
			spec[m - shift] += autoDouble[k] * std::exp(2.0 * M_PI * std::complex<double>(0.0, 1.0) * double(k * m) / double(2 * N - 1));
	}

	// transform the frequency axis to angular frequency in eV
	// and scale the spectrum with the square root of the number of points:
	for (int i = 0; i < 2 * N - 1; i++)
	{
		freq[i] *= 2.0 * M_PI * Eh2eV;
		spec[i] /= std::sqrt(double(2 * N - 1));
	}


	// write the results to the output files:
	std::ofstream autoOut("auto_out.dat");
	std::ofstream specOut("spec_out.dat");

	autoOut << "#       Time [au]         Re(Auto)        Im(Auto)        Abs(Auto)\n";
	specOut << "# Frequency [au]      Re(Spec)         Im(Spec)         Abs(Spec)\n";

	for (int i = 0; i < 2 * N - 1; i++)
	{
		autoOut << std::fixed << std::setprecision(10) << std::setw(20) << doubleTime[i]
				<< std::setw(16) << autoDouble[i].real()
				<< std::setw(16) << autoDouble[i].imag()
				<< std::setw(16) << std::abs(autoDouble[i]) << std::endl;
		specOut << std::fixed << std::setprecision(10) << std::setw(15) << freq[i]
				<< std::setw(17) << spec[i].real()
				<< std::setw(17) << spec[i].imag()
				<< std::setw(17) << std::abs(spec[i]) << std::endl;

	}
}
