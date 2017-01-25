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


#include <fstream>
#include <vector>
#include <eigen3/Eigen/Core>
#include <armadillo>

#ifndef SOURCE_VIBRATIONALANALYSIS_H_
#define SOURCE_VIBRATIONALANALYSIS_H_

class GaussFchk;

class VibrationalAnalysis
{
public:
	VibrationalAnalysis(GaussFchk &initFchk);
	VibrationalAnalysis(std::ifstream &initStream);
	virtual ~VibrationalAnalysis();

	void createThirdDerivs(std::string &baseName);
	void readAnharm(const std::string &GaussLogName);
	void createAnharmMCTDHoper(const std::string &baseName, const double thres);
	void calcKappa(GaussFchk &esFchk);

	void prtMinGeo();
	void prtMinModes();
	void prtDiagHess();
	void prtFreqs();
	void prtCubics();

	int Natoms() const;
	int Ncoords() const;
	int Nmodes() const;

	double totalMass() const;

	Eigen::VectorXd X() const;
    Eigen::VectorXd Xs() const;
    Eigen::VectorXd Xsr() const;
	Eigen::VectorXd Q() const;
    Eigen::VectorXd Qs() const;
    Eigen::VectorXd Qsr() const;
	Eigen::VectorXd masses() const;

	Eigen::VectorXi atomicNumbers() const;

	Eigen::MatrixXd Fcart() const;
	Eigen::MatrixXd Fmwc() const;
	Eigen::MatrixXd Fint() const;

	Eigen::MatrixXd D() const;

	Eigen::MatrixXd Lint() const;
	Eigen::MatrixXd Lmwc() const;
	Eigen::MatrixXd Lcart() const;

	Eigen::VectorXd mwcFC() const;
	Eigen::VectorXd mwcFrq() const;

	Eigen::VectorXd intFC() const;
	Eigen::VectorXd intFrq() const;
	Eigen::VectorXd auFrq() const;

	Eigen::VectorXd mu() const;

	Eigen::Vector3d com() const;
	Eigen::Matrix3d inert() const;
	Eigen::Vector3d moments() const;
	Eigen::Matrix3d prinAxes() const;

	arma::Cube<double> thirdDerivs() const;
	Eigen::VectorXd diag4thDerivs() const;


private:
	void constructor();			// dummy function to setup the class
	void calcCom();				// calculate the center of mass
	void calcInert();			// calculate and diagonalize the inertia tensor

	GaussFchk *m_fchk;				// the formatted checkpoint file
	int m_Natoms;					// number of atoms
	int m_Ncoords;					// number of cartesian coordinates (3N)
	int m_Nmodes;					// number of normal modes (3N-6)
	const double shiftFac = 0.03;	// factor used for the displacement in cubic force consts

	Eigen::Vector3d m_com;		// center of mass
	Eigen::Vector3d m_moments;	// principal moments (eigenvalues of the inertia tensor)

	Eigen::Matrix3d m_inert;	// inertia tensor
	Eigen::Matrix3d m_prinAxes;	// principal axes (eigenvectors of the inertia tensor)

	Eigen::VectorXd m_masses;	// atomic masses
	Eigen::VectorXi atNums;		// atomic number
	Eigen::VectorXd MassVec;	// atomic mass vector (length 3N)
	Eigen::VectorXd X_min;		// molecular geometry at the minimum (cartesian)
	Eigen::VectorXd Q_min;		// molecular geometry at the minimum (mass-weighted cart.)
	Eigen::VectorXd X_s;		// molecular geometry relative to the center of mass
	Eigen::VectorXd Q_s;		// mass-weighted version of the above
	Eigen::VectorXd X_sr;		// molecular geometry in the inertia frame
	Eigen::VectorXd Q_sr;		// mass-weighted
	Eigen::VectorXd MwcFrcCon;	// force constants (full Hessian)
	Eigen::VectorXd MwcFreqs;	// vibrational frequencies [cm-1] (full Hessian)
	Eigen::VectorXd IntFrcCon;	// force constants (internal Hessian)
	Eigen::VectorXd IntFreqs;	// vibrational frequencies [cm-1] (internal Hessian)
	Eigen::VectorXd AUFreqs;	// frequencies of the normal modes in atomic units
	Eigen::VectorXd RedMasses;	// 'reduced masses' of the normal modes

	Eigen::VectorXd fourthDerivs;	// diagonal fourth derivatives

	Eigen::MatrixXd MassInvMat;	// diagonal matrix of inverse square roots of atomic masses
	Eigen::MatrixXd Fcart_min;	// Hessian matrix at the minimum (cartesian)
	Eigen::MatrixXd Fmwc_min;	// Hessian matrix at the minimum (mass-weighted cart.)
	Eigen::MatrixXd Fint_min;	// Hessian matrix at the minimum (internal coordinates)
	Eigen::MatrixXd GaussModes;	// Normal modes calculated by gaussian
	Eigen::MatrixXd BigAxes;	// 3N x 3N rotation matrix based on the principal axes
	Eigen::MatrixXd Dmat;		// the D-matrix to transform
	Eigen::MatrixXd Lint_min;	// Normal modes in 'internal coordinates' at the minimum
	Eigen::MatrixXd Lmwc_min;	// mass-weighted cartesian normal modes at the minimum
	Eigen::MatrixXd Lcrt_min;	// cartesian normal modes at the minimum

	std::vector<Eigen::MatrixXd> Fcart_Dp;		// cartesian Hessians at positively displaced geometries
	std::vector<Eigen::MatrixXd> Fcart_Dn;		// cartesian Hessians at negatively displaced geometries
	std::vector<Eigen::MatrixXd> Fmwc_Dp;		// mass-weighted Hessians at positively displaced geometries
	std::vector<Eigen::MatrixXd> Fmwc_Dn;		// mass-weighted Hessians at negatively displaced geometries
	std::vector<Eigen::MatrixXd> Fdiag_Dp;		// mass-weighted 'diagonalized' Hessians at positively displaced geometries
	std::vector<Eigen::MatrixXd> Fdiag_Dn;		// mass-weighted 'diagonalized' Hessians at negatively displaced geometries

	arma::Cube<double> rawDerivs;	// tensor containing the raw third derivatives (individual terms)
	arma::Cube<double> avgDerivs;	// tensor containing the averaged third derivatives
};

#endif /* SOURCE_VIBRATIONALANALYSIS_H_ */
