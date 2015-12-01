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
#include <string>
#include <fstream>
#include <eigen3/Eigen/Core>

#ifndef SOURCE_GAUSSFCHK_H_
#define SOURCE_GAUSSFCHK_H_

class GaussFchk
{
public:
	GaussFchk(std::ifstream &file);
	virtual ~GaussFchk();
	int ReadInteger(const std::string &searchString) const;
	double ReadReal(const std::string &searchString) const;
	Eigen::VectorXd ReadVector(const std::string &searchString) const;
	Eigen::MatrixXd ReadMatrix(const std::string &searchString, const size_t columns, const size_t rows) const;
	Eigen::MatrixXd ReadSymmetricMatrix(const std::string &searchString) const;

private:
	std::ifstream *m_FchkFile;
};

#endif /* SOURCE_GAUSSFCHK_H_ */
