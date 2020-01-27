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
    Eigen::Vector3d dr_PBC_corrected(const Eigen::Vector3d& v, const Eigen::Vector3d& w) const {\
        static Eigen::Vector3d res;
        res[0] = w[0] - v[0] - std::floor((w[0] - v[0]) / boxDimension[0] + 0.5) * boxDimension[0];
        res[1] = w[1] - v[1] - std::floor((w[1] - v[1]) / boxDimension[1] + 0.5) * boxDimension[1];
        res[2] = w[2] - v[2] - std::floor((w[2] - v[2]) / boxDimension[2] + 0.5) * boxDimension[2];
        return res;
    }

    /* Puts a 3vector v back in the simulation box.*/
    Eigen::Vector3d updatePostionPBC(const Eigen::Vector3d& v) const {
        return v.array() - floor(v.array() / boxDimension.array()) * boxDimension.array();
    }
private:
    Eigen::Vector3d boxDimension;
};