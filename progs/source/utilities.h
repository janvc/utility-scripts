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


#ifndef SOURCE_UTILITIES_H_
#define SOURCE_UTILITIES_H_

void WriteMatrix(const Eigen::MatrixXd &mat, const int dig = 5, const bool clean = false, const double thres = 1.0e-10);
void WriteMatrixToFile(std::ofstream &stream, const Eigen::MatrixXd &mat, const int dig = 5, const bool clean = false, const double thres = 1.0e-10);

void WriteVector(const Eigen::VectorXd &vec, const int dig = 5, const bool clean = false, const double thres = 1.0e-10);
void WriteVectorToFile(std::ofstream &stream, const Eigen::VectorXd &vec, const int dig = 5, const bool clean = false, const double thres = 1.0e-10);

void WriteFortranNumber(std::ofstream &stream, const double number);

Eigen::Vector3d calc_com(Eigen::VectorXd x, Eigen::VectorXd m);
Eigen::Matrix3d calc_inert(Eigen::VectorXd x, Eigen::VectorXd m);

Eigen::MatrixXd vibrationalAnalysis(const Eigen::VectorXd &positions, const Eigen::VectorXd &masses, const Eigen::MatrixXd &hessian);

void createMCTDHfiles(const Eigen::MatrixXd &J, const Eigen::VectorXd &K, const Eigen::VectorXd &f1, const Eigen::VectorXd &f2, const double dE, std::ofstream &logFile);

#endif /* SOURCE_UTILITIES_H_ */
