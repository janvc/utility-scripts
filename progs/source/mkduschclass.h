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


#include "GaussFchk.h"
#include "vibrationalanalysis.h"

#ifndef MKDUSCHCLASS_H
#define MKDUSCHCLASS_H

class MkDuschClass
{
public:
    MkDuschClass(GaussFchk &gState, GaussFchk &eState);

    void createMCTDHfiles(const std::string &basename, const bool freqWeight = false, const bool noJ = false, const bool noK = false) const;
    void writeLog() const;

    void calcMCTDHdata();
    void selectModes();

private:
    int m_Natoms;
    int m_Ncoords;
    int m_Nmodes;
    int m_Nactive;
    double m_gsEner;
    double m_esEner;
    double m_dE;
    double m_zpeGtot;
    double m_zpeEtot;
    double m_RMSDb;
    double m_RMSDa;
    Eigen::Vector3d m_COMg;     // center of mass after shifting
    Eigen::Vector3d m_COMe;
    Eigen::VectorXd m_XSg;      // shifted coordinates
    Eigen::VectorXd m_XSe;
    Eigen::VectorXd m_QSg;
    Eigen::VectorXd m_QSe;
    Eigen::VectorXd m_XRg;      // rotated coordinates
    Eigen::VectorXd m_QRg;
    Eigen::VectorXd m_Kfull;    // full Displacement vector
    Eigen::Matrix3d m_CorrMat;
    Eigen::Matrix3d m_RotMat;   // rotation matrix to minimize RMSD
    Eigen::Matrix4d m_F;        // quaternion matrix
    Eigen::MatrixXd m_Jfull;    // full Duschinsky matrix
    VibrationalAnalysis m_GroundState;
    VibrationalAnalysis m_ExcitedState;

    Eigen::VectorXd m_f1;       //
    Eigen::VectorXd m_f2;       //
    Eigen::VectorXd m_fp;       // parameters for the potential energy surface
    Eigen::MatrixXd m_phi;      //
    Eigen::VectorXd m_kappa;    //
    Eigen::VectorXd m_d;        //

    Eigen::VectorXd m_fp_fw;    //
    Eigen::MatrixXd m_phi_fw;   // additional frequency-weighted parameters
    Eigen::VectorXd m_kappa_fw; //
};

#endif // MKDUSCHCLASS_H
