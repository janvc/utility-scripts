/*
 * Copyright 2016-2017 Jan von Cosel
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
#include <fftw3.h>
#include "constants.h"


std::vector<std::complex<double> > calc_fftw(std::vector<std::complex<double> > data)
{
    int N = data.size();
    std::vector<std::complex<double> > spec(N);

    fftw_complex *fftw_data;
    fftw_complex *fftw_spec;
    fftw_data = reinterpret_cast<fftw_complex*>(data.data());
    fftw_spec = reinterpret_cast<fftw_complex*>(spec.data());
    fftw_plan plan = fftw_plan_dft_1d(N, fftw_data, fftw_spec, FFTW_BACKWARD, FFTW_ESTIMATE);

    fftw_execute(plan);
    fftw_destroy_plan(plan);

    for (int i = 0; i < N; i++)
        spec[i] /= double(N);

    return spec;
}

int main(int argc, char *argv[])
{
    std::string fileName;
    double dipMom;
    double tau;
    double vibFreq;
    namespace po = boost::program_options;
    po::options_description desc("Allowed options");
    desc.add_options()
                    ("help,h", "produce this help message")
                    ("file,f", po::value<std::string>(&fileName)->required(), "the autocorrelation file")
                    ("tau,t", po::value<double>(&tau)->required(), "the damping time")
                    ("dipole,d", po::value<double>(&dipMom)->required(), "square of the transition dipole moment")
                    ("vibfreq,v", po::value<double>(&vibFreq), "frequency in cm-1 of the excited normal mode")
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
    std::vector<double> inpTime;

    inputFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    // read the input correlation function, scale it with the
    // transition dipole moment and convert the time to au.
    double dummy1, dummy2, dummy3, dummy4;
    while (inputFile >> dummy1 >> dummy2 >> dummy3 >> dummy4)
    {
        inpTime.push_back(dummy1 * fs2au);
        std::complex<double> dummy_comp;
        dummy_comp.real(dummy2 * dipMom);
        dummy_comp.imag(dummy3 * dipMom);
        autoInp.push_back(dummy_comp);
    }
    inputFile.close();

    std::cout << "Autocorrelation file: " << fileName << std::endl;
    std::cout << "damping time: " << tau / fs2au << " fs\n";
    std::cout << "transition dipole moment: " << dipMom << " au\n";
    if (vm.count("vibfreq"))
        std::cout << "offset vibrational frequency: " << vibFreq << " cm-1\n";


    // calculate the required lengths, times and frequencies:
    int Ninp = int(inpTime.size());
    int Nnz = 2 * Ninp - 1;
    int N = 1;
    while (N < Nnz)
        N *= 2;
    double dt = (inpTime[1] - inpTime[0]);
    double dw = 2.0 * M_PI / (N * dt);
    std::vector<double> time(N);
    std::vector<double> freq(N);
    for (int i = 0; i < N; i++)
    {
        time[i] = dt * (i - (N / 2) + 1);
        freq[i] = (i + 1) * Eh2eV * dw;
    }

    std::cout << "input number of points: " << Ninp << std::endl;
    std::cout << "total number of points: " << N << std::endl;
    std::cout << "dt: " << dt << " au, " << dt / fs2au << " fs\n";
    std::cout << "dw: " << dw << " au, " << dw * Eh2eV << " eV\n";
    std::cout << "factor: " << (N * N * dt * dt) / (2.0 * M_PI) << std::endl;

    // compute the integral (i.e. area) of the absolute squared MCTDH correlation function:
    double mctdhInt = 0.0;
    for (int i = 0; i < Ninp - 1; i++)
    {
        double leftVal = std::abs(autoInp[i]) * std::abs(autoInp[i]);
        double rightVal = std::abs(autoInp[i+1]) * std::abs(autoInp[i+1]);

        //mctdhInt += 0.5 * dt * (leftVal + rightVal);
        mctdhInt += 0.5 * (leftVal + rightVal);
    }
    std::cout << "integral of the MCTDH correlation function: " << mctdhInt << std::endl;

    // if requested, set the initial phase of the autocorrelation
    // function to zero and re-normalize it:
    if (vm.count("phase"))
    {
        double phase = -std::atan2(autoInp[0].imag(), autoInp[0].real());
        double norm = 1.0 / std::abs(autoInp[0]);
        std::cout << "applying a phase of " << phase << " to the autocorrelation function\n";
        for (int i = 0; i < N; i++)
            autoInp[i] *= std::exp(std::complex<double>(0.0, 1.0) * phase) * norm * dipMom;
    }

    // if the autocorrelation function results from a pre-excited calculation, apply
    // the corresponding phase to account for the additional energy.
    if (vm.count("vibfreq"))
    {
        for (int i = 0; i < Ninp; i++)
            autoInp[i] *= std::exp(std::complex<double>(0.0, 1.0) * vibFreq * inpTime[i] / 219474.6312);
    }

    // apply the damping to the initial correlation function:
    for (int i = 0; i < Ninp; i++)
        autoInp[i] *= std::exp(-inpTime[i] / tau);


    // mirror the correlation function to make use of the hermiticity.
    // because the point at time 0 is mapped onto itself, the new correlation
    // function is 2N-1 long instead of 2N.
    std::vector<std::complex<double> > autoNz(Nnz);
    for (int i = 0; i < Ninp; i++)
    {
        autoNz[i] = std::conj(autoInp[Ninp - i - 1]);
        autoNz[2 * Ninp - 2 - i] = autoInp[Ninp - i - 1];
    }

    // create the full autocorrelation function by padding the mirrored
    // one with zeros on both sides to obtain a power of two.
    std::vector<std::complex<double> > autoFull(N);
    int Npad = (N - Nnz) / 2;
    for (int i = 0; i < Nnz; i++)
        autoFull[Npad + i] = autoNz[i];


    // compute the integral (i.e. area) of the absolute squared correlation function:
    double corrInt = 0.0;
    for (int i = 0; i < N - 1; i++)
    {
        double leftVal = std::abs(autoFull[i]) * std::abs(autoFull[i]);
        double rightVal = std::abs(autoFull[i+1]) * std::abs(autoFull[i+1]);

        //corrInt += 0.5 * dt * (leftVal + rightVal);
        corrInt += 0.5 * (leftVal + rightVal);
    }
    std::cout << "integral of the correlation function: " << corrInt << std::endl;

    // use FFTW to compute the fourier transform.
    std::vector<std::complex<double> > spec = calc_fftw(autoFull);

    // compute the integral (i.e. area) of the absolute squared initial spectrum:
    double specInt = 0.0;
    for (int i = 0; i < N - 1; i++)
    {
        double leftVal = std::abs(spec[i]) * std::abs(spec[i]);
        double rightVal = std::abs(spec[i+1]) * std::abs(spec[i+1]);

        //specInt += 0.5 * dw * (leftVal + rightVal);
        specInt += 0.5 * (leftVal + rightVal);
    }
    std::cout << "integral of the initial spectrum: " << specInt << std::endl;


    // apply the frequency-dependent prefactor and the phase that cancels the imaginary part
    // of the spectrum. We have an additional factor of 100 here to account for the conversion
    // from meters to centimeters in the transition dipole moment.
    std::vector<std::complex<double> > final_spec(N);
    for (int i = 0; i < N; i++)
        final_spec[i] = spec[i] * freq[i] * dt * std::exp(std::complex<double>(0.0, 1.0) * dw * time[0] * double(i)) * facabs * 100.0;

    // compute the integral (i.e. area) of the absolute squared final spectrum:
    double finalInt = 0.0;
    for (int i = 0; i < N - 1; i++)
    {
        double leftVal = std::abs(final_spec[i]) * std::abs(final_spec[i]);
        double rightVal = std::abs(final_spec[i+1]) * std::abs(final_spec[i+1]);

        finalInt += 0.5 * dw * (leftVal + rightVal);
    }
    std::cout << "integral of the final spectrum: " << finalInt << std::endl;


    // write the results to the output files:
    std::ofstream autoOut("auto_out.dat");
    std::ofstream specOut("spec_out.dat");
    std::ofstream fspecOut("final_spec_out.dat");

    autoOut << "#       Time [fs]         Re(Auto)        Im(Auto)        Abs(Auto)\n";
    specOut << "# Frequency [eV]      Re(Spec)        Im(Spec)        Abs(Spec)\n";
    fspecOut << "# Frequency [eV]         Re(Spec)           Im(Spec)           Abs(Spec)\n";

    for (int i = 0; i < N; i++)
    {
        autoOut << std::fixed << std::setprecision(10) << std::setw(20) << time[i] / fs2au
                << std::setw(16) << autoFull[i].real()
                << std::setw(16) << autoFull[i].imag()
                << std::setw(16) << std::abs(autoFull[i]) << std::endl;
        specOut << std::fixed << std::setprecision(10) << std::setw(16) << freq[i]
                << std::setw(16) << spec[i].real()
                << std::setw(16) << spec[i].imag()
                << std::setw(16) << std::abs(spec[i]) << std::endl;
        fspecOut << std::fixed << std::setprecision(10) << std::setw(16) << freq[i]
                 << std::setw(19) << final_spec[i].real()
                 << std::setw(19) << final_spec[i].imag()
                 << std::setw(19) << std::abs(final_spec[i]) << std::endl;

    }


    // Compute the zero'th moment of the spectrum.
    double M0 = 0.0;
    for (int i = 0; i < N / 2; i++)
        M0 += dw * Eh2eV * (final_spec[i].real() + final_spec[i+1].real()) * 0.5;

    std::cout << "zeroth moment: " << M0 << std::endl;

    // Compute the first moment of the spectrum.
    double M1 = 0.0;
    for (int i = 0; i < N / 2; i++)
        M1 += dw * Eh2eV * freq[i] * (final_spec[i].real() + final_spec[i+1].real()) * 0.5;

    std::cout << "first moment: " << M1 / M0 << std::endl;

    // Compute the second moment of the spectrum.
    double M2 = 0.0;
    for (int i = 0; i < N / 2; i++)
        M2 += dw * Eh2eV * freq[i] * freq[i] * (final_spec[i].real() + final_spec[i+1].real()) * 0.5;

    std::cout << "second moment: " << M2 / M0 << std::endl;
    double specWidth = std::sqrt((M2 / M0) - ((M1 / M0) * (M1 / M0)));
    std::cout << "spectrum width: " << specWidth << std::endl;


}
