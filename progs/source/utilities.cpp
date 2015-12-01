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


#include <eigen3/Eigen/Core>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "utilities.h"


void WriteMatrix(const Eigen::MatrixXd &mat, const int dig, const bool clean, const double thres)
{
	int width = dig + 8;
	// write header:
	std::cout << "  ";
	for (int i = 0; i < mat.cols(); i++)
		std::cout << std::setw(width) << i+1;
	std::cout << std::endl;

	for (int i = 0; i < mat.rows(); i++)
	{
		// write row number
		std::cout << std::setw(5) << i+1;
		for (int j = 0; j < mat.cols(); j++)
		{
			if (clean)
			{
				if (std::abs(mat(i,j)) > thres)
					std::cout << std::scientific << std::setprecision(dig) << std::setw(width) << mat(i,j);
				else
					std::cout << std::setw(width) << "0";
			}
			else
				std::cout << std::scientific << std::setprecision(dig) << std::setw(width) << mat(i,j);
		}
		std::cout << std::endl;
	}
}


void WriteVector(const Eigen::VectorXd &vec, const int dig, const bool clean, const double thres)
{
	int width = dig + 8;

	for (int i = 0; i < vec.size(); i++)
	{
		std::cout << std::setw(5) << i+1;
		if (clean)
		{
			if (std::abs(vec(i)) > thres)
				std::cout << std::scientific << std::setprecision(dig) << std::setw(width) << vec(i) << std::endl;
			else
				std::cout << std::setw(width) << "0" << std::endl;
		}
		else
			std::cout << std::scientific << std::setprecision(dig) << std::setw(width) << vec(i) << std::endl;
	}
}


