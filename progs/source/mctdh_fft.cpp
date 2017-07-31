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
                    ("tau,t", po::value<double>(&tau), "the damping time")
                    ("dipole,d", po::value<double>(&dipMom), "square of the transition dipole moment")
                    ("vibfreq,v", po::value<double>(&vibFreq), "frequency in cm-1 of the excited normal mode")
                    ("phase,p", "set the initial phase to zero")
                    ("fcc", "use a correlation function from FCClasses")
                    ;
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    if (vm.count("help"))
    {
        std::cout << desc << std::endl;
        return 1;
    }
    po::notify(vm);


    /*
     * Do general stuff that is common to MCTDH and FCClasses.
     * Everything is being done in atomic units
     */
    std::ifstream inputFile(fileName, std::ifstream::in);
    std::cout << "Autocorrelation file: " << fileName << std::endl;
    std::ofstream specOut("spec_out.dat");
    std::ofstream lsOut("lineshape_out.dat");
    std::ofstream fspecOut("final_spec_out.dat");

    /*
     * Global if construct that separates the MCTDH case from the FCClasses case
     */
    if (! vm.count("fcc"))  // MCTDH
    {
        if (! vm.count("tau"))
        {
            std::cerr << "Error: Using an MCTDH correlation function but no damping time given. Exiting!\n";
            return -1;
        }
        if (! vm.count("dipole"))
        {
            std::cerr << "Error: Using an MCTDH correlation function but no transition dipole moment given. Exiting!\n";
            return -1;
        }

        std::cout << "we are treating an MCTDH correlation function.\n";

        // convert the damping time to atomic units:
        tau *= fs2au;
        std::cout << "damping time: " << tau / fs2au << " fs\n";

        std::vector<std::complex<double> > autoInp;
        std::vector<double> timeInp;

        // read the input correlation function, scale it with the
        // transition dipole moment and convert the time to au.
        double dummy1, dummy2, dummy3, dummy4;
        inputFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        while (inputFile >> dummy1 >> dummy2 >> dummy3 >> dummy4)
        {
            timeInp.push_back(dummy1 * fs2au);
            std::complex<double> dummy_comp;
            dummy_comp.real(dummy2 * dipMom);
            dummy_comp.imag(dummy3 * dipMom);
            autoInp.push_back(dummy_comp);
        }

        std::cout << "transition dipole moment: " << dipMom << " au\n";

        // calculate the required lengths, times and frequencies:
        int Ninp = int(timeInp.size());
        int Nnz = 2 * Ninp - 1;
        int N = 1;
        while (N < Nnz)
            N *= 2;
        double dt = (timeInp[1] - timeInp[0]);
        double dw = 2.0 * M_PI / (N * dt);
        std::vector<double> time(N);
        std::vector<double> freq(N);
        for (int i = 0; i < N; i++)
        {
            time[i] = dt * (i - (N / 2) + 1);
            freq[i] = (i + 1) * dw;
        }

        std::cout << "input number of points: " << Ninp << std::endl;
        std::cout << "total number of points: " << N << std::endl;
        std::cout << "dt: " << dt << " au, " << dt / fs2au << " fs\n";
        std::cout << "dw: " << dw << " au, " << dw * Eh2eV << " eV\n";

        // compute the integral (i.e. area) of the absolute MCTDH correlation function:
        double mctdhInt = 0.0;
        for (int i = 0; i < Ninp - 1; i++)
            mctdhInt += 0.5 * dt * (std::abs(autoInp[i]) + std::abs(autoInp[i+1]));
        std::cout << "integral of the absolute initial MCTDH correlation function: " << mctdhInt << std::endl;

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
            std::cout << "offset vibrational frequency: " << vibFreq << " cm-1\n";
            for (int i = 0; i < Ninp; i++)
                autoInp[i] *= std::exp(std::complex<double>(0.0, 1.0) * vibFreq * timeInp[i] / 219474.6312);
        }

        // apply the damping to the initial correlation function:
        for (int i = 0; i < Ninp; i++)
            autoInp[i] *= std::exp(-timeInp[i] / tau) * std::cos(M_PI * timeInp[i] / (2.0 * timeInp[Ninp-1]));

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

        // compute the integral (i.e. area) of the absolute correlation function:
        double corrInt = 0.0;
        for (int i = 0; i < N - 1; i++)
            corrInt += 0.5 * dt * (std::abs(autoFull[i]) + std::abs(autoFull[i+1]));
        std::cout << "integral of the absolute full correlation function: " << corrInt << std::endl;

        // use FFTW to compute the fourier transform.
        std::vector<std::complex<double> > spec = calc_fftw(autoFull);

        // compute the integral (i.e. area) of the absolute initial spectrum:
        double specInt = 0.0;
        for (int i = 0; i < N - 1; i++)
            specInt += 0.5 * dw * (std::abs(spec[i]) + std::abs(spec[i+1]));
        std::cout << "integral of the absolute initial spectrum: " << specInt << std::endl;

        // multiply the spectrum by the phase factor to get the approximate continuous
        // fourier transform. Also, divide by 2 pi. This will give us the lineshape.
        std::vector<std::complex<double> > lineshape_spec(N);
        for (int i = 0; i < N; i++)
            lineshape_spec[i] = spec[i] * dt * std::exp(std::complex<double>(0.0, 1.0) * dw * time[0] * double(i)) / (2.0 * M_PI);

        // compute the integral (i.e. area) of the real part of the lineshape.
        double lineshapeInt = 0.0;
        for (int i = 0; i < N - 1; i++)
            lineshapeInt += 0.5 * dw * (lineshape_spec[i].real() + lineshape_spec[i+1].real());
        std::cout << "integral of the real part of the lineshape: " << lineshapeInt << std::endl;

        // apply the frequency-dependent prefactor to get the spectrum in absolute units.
        std::vector<std::complex<double> > final_spec(N);
        for (int i = 0; i < N; i++)
            final_spec[i] = lineshape_spec[i] * freq[i] * facabs;

        // compute the integral (i.e. area) of the real part of the final spectrum:
        double finalInt = 0.0;
        for (int i = 0; i < N - 1; i++)
            finalInt += 0.5 * dw * (final_spec[i].real() + final_spec[i+1].real());
        std::cout << "integral of the real part of the final spectrum: " << finalInt << std::endl;

        // write the results to the output files:
        std::ofstream autoOut("auto_out.dat");

        autoOut << "#       Time [fs]         Re(Auto)        Im(Auto)        Abs(Auto)\n";
        specOut << "# Frequency [eV]      Re(Spec)        Im(Spec)        Abs(Spec)\n";
        lsOut << "# Frequency [eV]      Re(Spec)        Im(Spec)        Abs(Spec)\n";
        fspecOut << "# Frequency [eV]         Re(Spec)           Im(Spec)           Abs(Spec)\n";

        for (int i = 0; i < N; i++)
        {
            autoOut << std::fixed << std::setprecision(10) << std::setw(20) << time[i] / fs2au
                    << std::setw(16) << autoFull[i].real()
                    << std::setw(16) << autoFull[i].imag()
                    << std::setw(16) << std::abs(autoFull[i]) << std::endl;
            specOut << std::fixed << std::setprecision(10) << std::setw(16) << freq[i] * Eh2eV
                    << std::setw(16) << spec[i].real()
                    << std::setw(16) << spec[i].imag()
                    << std::setw(16) << std::abs(spec[i]) << std::endl;
            lsOut << std::fixed << std::setprecision(10) << std::setw(16) << freq[i] * Eh2eV
                    << std::setw(16) << lineshape_spec[i].real()
                    << std::setw(16) << lineshape_spec[i].imag()
                    << std::setw(16) << std::abs(lineshape_spec[i]) << std::endl;
            fspecOut << std::fixed << std::setprecision(10) << std::setw(16) << freq[i] * Eh2eV
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
        M1 /= M0;

        std::cout << "first moment: " << M1 << std::endl;

        // Compute the second moment of the spectrum.
        double M2 = 0.0;
        for (int i = 0; i < N / 2; i++)
            M2 += dw * Eh2eV * freq[i] * freq[i] * (final_spec[i].real() + final_spec[i+1].real()) * 0.5;
        M2 /= M0;

        std::cout << "second moment: " << M2 << std::endl;
        double specWidth = std::sqrt(M2 - (M1 * M1));
        std::cout << "spectrum width: " << specWidth << std::endl;
    }
    else    // FCClasses
    {
        std::cout << "we are treating an FCClasses correlation function.\n";

        std::vector<std::complex<double> > autoInp;
        std::vector<double> timeInp;

        double dummy1, dummy2, dummy3;
        while (inputFile >> dummy1 >> dummy2 >> dummy3)
        {
            timeInp.push_back(dummy1 * fs2au);
            std::complex<double> dummy_comp;
            dummy_comp.real(dummy2);
            dummy_comp.imag(dummy3);
            autoInp.push_back(dummy_comp);
        }

        int N = int(timeInp.size());
        double dt = (timeInp[1] - timeInp[0]);
        double dw = 2.0 * M_PI / (N * dt);
        std::vector<double> freq(N);
        for (int i = 0; i < N; i++)
            freq[i] = i * dw;

        std::cout << "total number of points: " << N << std::endl;
        std::cout << "dt: " << dt << " au, " << dt / fs2au << " fs\n";
        std::cout << "dw: " << dw << " au, " << dw * Eh2eV << " eV\n";

        // compute the integral (i.e. area) of the absolute FCClasses correlation function:
        double fccInt = 0.0;
        for (int i = 0; i < N - 1; i++)
            fccInt += 0.5 * dt * (std::abs(autoInp[i]) + std::abs(autoInp[i+1]));
        std::cout << "integral of the absolute initial FCClasses correlation function: " << fccInt << std::endl;

        // use FFTW to compute the fourier transform.
        std::vector<std::complex<double> > spec = calc_fftw(autoInp);

        // compute the integral (i.e. area) of the absolute initial spectrum:
        double specInt = 0.0;
        for (int i = 0; i < N - 1; i++)
            specInt += 0.5 * dw * (std::abs(spec[i]) + std::abs(spec[i+1]));
        std::cout << "integral of the absolute initial spectrum: " << specInt << std::endl;

        // multiply the spectrum by the phase factor to get the approximate continuous
        // fourier transform. Also, divide by 2 pi. This will give us the lineshape.
        std::vector<std::complex<double> > lineshape_spec(N);
        for (int i = 0; i < N; i++)
            lineshape_spec[i] = spec[i] * dt * std::exp(std::complex<double>(0.0, 1.0) * dw * timeInp[0] * double(i)) / (2.0 * M_PI);

        // compute the integral (i.e. area) of the real part of the lineshape.
        double lineshapeInt = 0.0;
        for (int i = 0; i < N - 1; i++)
            lineshapeInt += 0.5 * dw * (lineshape_spec[i].real() + lineshape_spec[i+1].real());
        std::cout << "integral of the real part of the lineshape: " << lineshapeInt << std::endl;

        // apply the frequency-dependent prefactor to get the spectrum in absolute units.
        std::vector<std::complex<double> > final_spec(N);
        for (int i = 0; i < N; i++)
            final_spec[i] = lineshape_spec[i] * freq[i] * facabs;

        // compute the integral (i.e. area) of the real part of the final spectrum:
        double finalInt = 0.0;
        for (int i = 0; i < N - 1; i++)
            finalInt += 0.5 * dw * (final_spec[i].real() + final_spec[i+1].real());
        std::cout << "integral of the real part of the final spectrum: " << finalInt << std::endl;

        specOut << "# Frequency [eV]      Re(Spec)        Im(Spec)        Abs(Spec)\n";
        lsOut << "# Frequency [eV]      Re(Spec)        Im(Spec)        Abs(Spec)\n";
        fspecOut << "# Frequency [eV]         Re(Spec)           Im(Spec)           Abs(Spec)\n";

        for (int i = 0; i < N; i++)
        {
            specOut << std::fixed << std::setprecision(10) << std::setw(16) << freq[i] * Eh2eV
                    << std::setw(16) << spec[i].real()
                    << std::setw(16) << spec[i].imag()
                    << std::setw(16) << std::abs(spec[i]) << std::endl;
            lsOut << std::fixed << std::setprecision(10) << std::setw(16) << freq[i] * Eh2eV
                    << std::setw(16) << lineshape_spec[i].real()
                    << std::setw(16) << lineshape_spec[i].imag()
                    << std::setw(16) << std::abs(lineshape_spec[i]) << std::endl;
            fspecOut << std::fixed << std::setprecision(10) << std::setw(16) << freq[i] * Eh2eV
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
        M1 /= M0;

        std::cout << "first moment: " << M1 << std::endl;

        // Compute the second moment of the spectrum.
        double M2 = 0.0;
        for (int i = 0; i < N / 2; i++)
            M2 += dw * Eh2eV * freq[i] * freq[i] * (final_spec[i].real() + final_spec[i+1].real()) * 0.5;
        M2 /= M0;

        std::cout << "second moment: " << M2 << std::endl;
        double specWidth = std::sqrt(M2 - (M1 * M1));
        std::cout << "spectrum width: " << specWidth << std::endl;
    }
}
