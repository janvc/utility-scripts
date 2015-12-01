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


#include <sstream>
#include "GaussFchk.h"

GaussFchk::GaussFchk(std::ifstream &file)
	: m_FchkFile(&file)
{
}

GaussFchk::~GaussFchk()
{
}

int GaussFchk::ReadInteger(const std::string &searchString) const
{
	m_FchkFile->clear();
	m_FchkFile->seekg(0);

	std::string currentLine;
	int number;
	while (std::getline(*m_FchkFile, currentLine))
		if (currentLine.find(searchString) != std::string::npos && currentLine.at(43) == 'I' && currentLine.at(47) != 'N')
		{
			std::istringstream iss(currentLine.substr(50));
			iss >> number;

			return number;
		}

	return 0;
}

double GaussFchk::ReadReal(const std::string &searchString) const
{
	m_FchkFile->clear();
	m_FchkFile->seekg(0);

	std::string currentLine;
	double number;
	while (std::getline(*m_FchkFile, currentLine))
		if (currentLine.find(searchString) != std::string::npos && currentLine.at(43) == 'R' && currentLine.at(47) != 'N')
		{
			std::istringstream iss(currentLine.substr(50));
			iss >> number;

			return number;
		}

	return 0;
}

Eigen::VectorXd GaussFchk::ReadVector(const std::string &searchString) const
{
	m_FchkFile->clear();
	m_FchkFile->seekg(0);

	std::string currentLine;
	int vecSize;
	while (std::getline(*m_FchkFile, currentLine))
	{
		if (currentLine.find(searchString) != std::string::npos)
		{
			std::istringstream iss(currentLine.substr(50));
			iss >> vecSize;
			int remainder = vecSize % 5;

			Eigen::VectorXd Vec(vecSize);
			for (int i = 0; i < vecSize - remainder; i += 5)
				*m_FchkFile >> Vec(i) >> Vec(i+1) >> Vec(i+2) >> Vec(i+3) >> Vec(i+4);
			for (int i = vecSize - remainder; i < vecSize; i++)
				*m_FchkFile >> Vec(i);

			return Vec;
		}
	}
	return Eigen::VectorXd();
}

Eigen::MatrixXd GaussFchk::ReadMatrix(const std::string &searchString, const size_t columns, const size_t rows) const
{
	m_FchkFile->clear();
	m_FchkFile->seekg(0);

	std::string currentLine;
	int matSize = columns * rows;
	while (std::getline(*m_FchkFile, currentLine))
	{
		if (currentLine.find(searchString) != std::string::npos)
		{
			std::istringstream iss(currentLine.substr(50));
			iss >> matSize;
			if (matSize == columns * rows)
			{
				int remainder = matSize % 5;

				Eigen::VectorXd tmpVec(matSize);
				for (int i = 0; i < matSize - remainder; i += 5)
					*m_FchkFile >> tmpVec(i) >> tmpVec(i+1) >> tmpVec(i+2) >> tmpVec(i+3) >> tmpVec(i+4);
				for (int i = matSize - remainder; i < matSize; i++)
					*m_FchkFile >> tmpVec(i);

				Eigen::MatrixXd Mat(rows, columns);
				for (int i = 0; i < columns; i++)
					for (int j = 0; j < rows; j++)
						Mat(j,i) = tmpVec(i * rows + j);

				return Mat;
			}
		}
	}
	return Eigen::MatrixXd();
}

Eigen::MatrixXd GaussFchk::ReadSymmetricMatrix(const std::string &searchString) const
{
	m_FchkFile->clear();
	m_FchkFile->seekg(0);

	std::string currentLine;
	int nNumbers;
	while (std::getline(*m_FchkFile, currentLine))
	{
		if (currentLine.find(searchString) != std::string::npos)
		{
			std::istringstream iss(currentLine.substr(50));
			iss >> nNumbers;
			int N = -0.5 + sqrt(0.25 + 2 * nNumbers);	// number of rows/columns

			if (N * (N + 1) / 2 == nNumbers)
			{
				int remainder = nNumbers % 5;

				Eigen::VectorXd tmpVec(nNumbers);
				for (int i = 0; i < nNumbers - remainder; i += 5)
					*m_FchkFile >> tmpVec(i) >> tmpVec(i+1) >> tmpVec(i+2) >> tmpVec(i+3) >> tmpVec(i+4);
				for (int i = nNumbers - remainder; i < nNumbers; i++)
					*m_FchkFile >> tmpVec(i);

				Eigen::MatrixXd Mat(N, N);
				for (int i = 0; i < N; i++)
					for (int j = 0; j <= i; j++)
					{
						Mat(i,j) = tmpVec((i*(i+1)/2)+j);
						Mat(j,i) = Mat(i,j);
					}

				return Mat;
			}
		}
	}
	return Eigen::MatrixXd();
}
