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
#include <eigen3/Eigen/Eigenvalues>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include "utilities.h"
#include "constants.h"


void WriteMatrix(const Eigen::MatrixXd &mat, const int dig, const bool clean, const double thres)
{
	int width = dig + 8;
	// write header:
	std::cout << "  ";
	for (int i = 0; i < mat.cols(); i++)
		std::cout << std::setw(width) << i+1;
	std::cout << std::endl << "       ";
	for (int i = 0; i < mat.cols(); i++)
		std::cout << std::string(width,'-');
	std::cout << std::endl;

	for (int i = 0; i < mat.rows(); i++)
	{
		// write row number
		std::cout << std::setw(5) << i+1 << " |";
		for (int j = 0; j < mat.cols(); j++)
		{
			if (clean)
			{
				if (std::abs(double(mat(i,j))) > thres)
					std::cout << std::scientific << std::setprecision(dig) << std::setw(width) << mat(i,j);
				else
					std::cout << std::setw(width) << "0";
			}
			else
				std::cout << std::scientific << std::setprecision(dig) << std::setw(width) << mat(i,j);
		}
		std::cout << std::endl;

	}
	std::cout << "       ";
	for (int i = 0; i < mat.cols(); i++)
		std::cout << std::string(width,'-');
	std::cout << std::endl;
}


void WriteMatrixToFile(std::ofstream &stream, const Eigen::MatrixXd &mat, const int dig, const bool clean, const double thres)
{
	int width = dig + 8;
	// write header:
	stream << "  ";
	for (int i = 0; i < mat.cols(); i++)
		stream << std::setw(width) << i+1;
	stream << std::endl << "      /";
	for (int i = 0; i < mat.cols(); i++)
		stream << std::string(width,'-');
	stream << std::endl;

	for (int i = 0; i < mat.rows(); i++)
	{
		// write row number
		stream << std::setw(5) << i+1 << " |";
		for (int j = 0; j < mat.cols(); j++)
		{
			if (clean)
			{
				if (std::abs(double(mat(i,j))) > thres)
					stream << std::scientific << std::setprecision(dig) << std::setw(width) << mat(i,j);
				else
					stream << std::setw(width) << "0";
			}
			else
				stream << std::scientific << std::setprecision(dig) << std::setw(width) << mat(i,j);
		}
		stream << std::endl;
	}
	stream << "       ";
	for (int i = 0; i < mat.cols(); i++)
		stream << std::string(width,'-');
	stream << std::endl;
}

void WriteVector(const Eigen::VectorXd &vec, const int dig, const bool clean, const double thres)
{
	int width = dig + 8;

	for (int i = 0; i < vec.size(); i++)
	{
		std::cout << std::setw(5) << i+1 << " |";
		if (clean)
		{
			if (std::abs(double(vec(i))) > thres)
				std::cout << std::scientific << std::setprecision(dig) << std::setw(width) << double(vec(i)) << std::endl;
			else
				std::cout << std::setw(width) << "0" << std::endl;
		}
		else
			std::cout << std::scientific << std::setprecision(dig) << std::setw(width) << double(vec(i)) << std::endl;
	}
}


void WriteVectorToFile(std::ofstream &stream, const Eigen::VectorXd &vec, const int dig, const bool clean, const double thres)
{
	int width = dig + 8;

	for (int i = 0; i < vec.size(); i++)
	{
		stream << std::setw(5) << i+1 << " |";
		if (clean)
		{
			if (std::abs(double(vec(i))) > thres)
				stream << std::scientific << std::setprecision(dig) << std::setw(width) << double(vec(i)) << std::endl;
			else
				stream << std::setw(width) << "0" << std::endl;
		}
		else
			stream << std::scientific << std::setprecision(dig) << std::setw(width) << double(vec(i)) << std::endl;
	}
}


void WriteFortranNumber(std::ofstream &stream, const double number)
{
	std::ostringstream strstr;
	strstr << std::scientific << std::setprecision(8) << std::setw(15) << std::setfill(' ') << number;
	std::string str = strstr.str();
	std::replace(str.begin(), str.end(), 'e', 'd');
	stream << str;
}

Eigen::Vector3d calc_com(Eigen::VectorXd x, Eigen::VectorXd m)
{
	Eigen::Vector3d com = Eigen::Vector3d::Zero();

	for (int i = 0; i < m.size(); i++)
	{
		com(0) += m(i) * x(3 * i + 0);
		com(1) += m(i) * x(3 * i + 1);
		com(2) += m(i) * x(3 * i + 2);
	}
	com /= m.sum();

	return com;
}

Eigen::Matrix3d calc_inert(Eigen::VectorXd x, Eigen::VectorXd m)
{
	Eigen::Matrix3d inert = Eigen::Matrix3d::Zero();

	Eigen::Vector3d com = calc_com(x, m);

	for (int i = 0; i < m.size(); i++)
	{
		inert(0,0) += m(i) * ((x(3*i+1)-com(1))*(x(3*i+1)-com(1)) + (x(3*i+2)-com(2))*(x(3*i+2)-com(2)));
		inert(1,1) += m(i) * ((x(3*i+0)-com(0))*(x(3*i+0)-com(0)) + (x(3*i+2)-com(2))*(x(3*i+2)-com(2)));
		inert(2,2) += m(i) * ((x(3*i+0)-com(0))*(x(3*i+0)-com(0)) + (x(3*i+1)-com(1))*(x(3*i+1)-com(1)));
		inert(0,1) -= m(i) *  (x(3*i+0)-com(0)) * (x(3*i+1)-com(1));
		inert(0,2) -= m(i) *  (x(3*i+0)-com(0)) * (x(3*i+2)-com(2));
		inert(1,2) -= m(i) *  (x(3*i+1)-com(1)) * (x(3*i+2)-com(2));
	}

	inert(1,0) = inert(0,1);
	inert(2,0) = inert(0,2);
	inert(2,1) = inert(1,2);

	return inert;
}
