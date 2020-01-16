/***************************************************
 *
 * KMC MODEL FOR OPTOELECTRIC PROCESSES
 *
 * Author: Ruben Gerritsen
 *
 * Created on 14-01-2020
 *
 * PBC stores and computes all things relating to 
 * the periodic boundary conditions (PBC).
 *
 **************************************************/

#pragma once
#include <Eigen/Dense>
#include "Site.h"

class PBC {
public:

    PBC(double xDim, double yDim, double zDim) { boxDimension[0] = xDim; boxDimension[1] = yDim; boxDimension[2] = zDim; }

    /* Computes the 3vector dr pointing from v to w corrected for periodic boundary conditions. */
    Eigen::Vector3d dr_PBC_corrected(const Eigen::Vector3d& v, const Eigen::Vector3d& w) const {
        return w.array() - v.array() - floor((w.array() - v.array()) / boxDimension + 0.5) * boxDimension;
    }

    /* Puts a 3vector v back in the simulation box.*/
    Eigen::Vector3d updatePostionPBC(const Eigen::Vector3d& v) const {
        return v.array() - floor(v.array() / boxDimension) * boxDimension;
    }
private:
    Eigen::Array3d boxDimension;
};