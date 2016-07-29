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
#include <boost/program_options.hpp>
#include "vibrationalanalysis.h"
#include "GaussFchk.h"

int main(int argc, char *argv[])
{
	std::string filename_g;
	std::string filename_e;
	namespace po = boost::program_options;
	po::options_description desc("Allowed options");
	desc.add_options()
		("help,h", "produce this help message")
		("file,f", po::value<std::string>(&filename_g)->required(), "the Fchk file")
	;
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	if (vm.count("help"))
	{
		std::cout << desc << std::endl;
		return 1;
	}
	po::notify(vm);


	std::string basename = filename_g.substr(0, filename_g.find_last_of("."));
	std::string logName = basename + ".log";

	std::ifstream fchkFile_g(filename_g, std::ifstream::in);


	VibrationalAnalysis vibAn(fchkFile_g);

	vibAn.prtMinGeo();
	vibAn.prtFreqs();
	vibAn.prtMinModes();
	vibAn.readAnharm(logName);


	vibAn.createAnharmMCTDHoper("test_mctdh", 2.0e-9);

	return 0;
}
