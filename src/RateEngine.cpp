/***************************************************
 *
 * KMC MODEL FOR OPTOELECTRIC PROCESSES
 *
 * Author: Ruben Gerritsen
 *
 * Created on 14-01-2020
 *
 * RateEngine is the class containing all parameters
 * and functions to compute the rates in the KMC.
 *
 **************************************************/

#include <Eigen/Dense>
#include "RateEngine.h"
#include "PBC.h"
#include "PType.h"


double RateEngine::millerAbrahams(const Site& siteOne, const Site& siteTwo, const PType type) const {
    return 1.0;
}
/*
double RateEngine::millerAbrahams(const Site& siteOne, const Site& siteTwo, const PType type) const {
    Eigen::Vector3d dr = pbc.dr_PBC_corrected(siteTwo.getCoordinates(), siteOne.getCoordinates());
    double dist = dr.norm();
    double deltaE = siteTwo.getEnergy(type) - siteOne.getEnergy(type) + E_Field * charge[type] * dr[0];

    if ( deltaE <= 0) {
        return(v0[type] * std::exp(-2 * alpha[type] * dist));
    }
    else {
        return(v0[type] * std::exp(-2 * alpha[type] * dist - deltaE / kBT));
    }
} */

double RateEngine::millerAbrahamsGEN(const Site& siteOne, const Site& siteTwo, const PType type) const {
    Eigen::Vector3d dr = pbc.dr_PBC_corrected(siteTwo.getCoordinates(), siteOne.getCoordinates());
    double dist = dr.norm();
    double deltaE = siteTwo.getEnergy(type) - siteOne.getEnergy(type) + E_Field * charge[type] * dr[0] - deltaE_binding;

    if (deltaE <= 0) {
        return(v0[type] * std::exp(-2 * alpha[type] * dist));
    }
    else {
        return(v0[type] * std::exp(-2 * alpha[type] * dist - deltaE / kBT));
    }
}

double RateEngine::millerAbrahamsDIS(const Site& siteOne, const Site& siteTwo, const PType type) const {
    Eigen::Vector3d dr = pbc.dr_PBC_corrected(siteTwo.getCoordinates(), siteOne.getCoordinates());
    double dist = dr.norm();
    double deltaE = siteTwo.getEnergy(type) - siteOne.getEnergy(type) + E_Field * charge[type] * dr[0] + deltaE_binding;

    if (deltaE <= 0) {
        return(v0[type] * std::exp(-2 * alpha[type] * dist));
    }
    else {
        return(v0[type] * std::exp(-2 * alpha[type] * dist - deltaE / kBT));
    }
}

double RateEngine::millerAbrahamsCT_DIS(const Site& siteOne, const Site& siteTwo, const PType type) const {
    Eigen::Vector3d dr = pbc.dr_PBC_corrected(siteTwo.getCoordinates(), siteOne.getCoordinates());
    double dist = dr.norm();
    double deltaE = siteTwo.getEnergy(type) - siteOne.getEnergy(type) + E_Field * charge[type] * dr[0] + deltaE_CTbinding;

    if (deltaE <= 0) {
        return(v0[type] * std::exp(-2 * alpha[type] * dist));
    }
    else {
        return(v0[type] * std::exp(-2 * alpha[type] * dist - deltaE / kBT));
    }
}

double RateEngine::millerAbrahamsCT(const Site& siteOne, const Site& siteTwo, const PType type) const {
    Eigen::Vector3d dr = pbc.dr_PBC_corrected(siteTwo.getCoordinates(), siteOne.getCoordinates());
    double dist = dr.norm();
    double deltaE = siteTwo.getEnergy(type) - siteOne.getEnergy(type) + E_Field * charge[type] * dr[0] + deltaE_SingtoCT;

    if (deltaE <= 0) {
        return(v0[type] * std::exp(-2 * alpha[type] * dist));
    }
    else {
        return(v0[type] * std::exp(-2 * alpha[type] * dist - deltaE / kBT));
    }
}

double RateEngine::forster(const Site& siteOne, const Site& siteTwo) const {
    Eigen::Vector3d dr = pbc.dr_PBC_corrected(siteTwo.getCoordinates(), siteOne.getCoordinates());
    double dist = dr.norm();
    double deltaE = siteTwo.getEnergy(PType::sing) - siteOne.getEnergy(PType::sing) + E_Field * charge[PType::sing] * dr[0];

    if (deltaE <= 0) {
        return( v0[PType::sing] * std::pow(R_forster / dist, 6));
    }
    else {
        return( v0[PType::sing] * std::pow( R_forster/ dist, 6) * std::exp( - deltaE / kBT) );
    }
}

double RateEngine::decay(const PType type) const {
    switch (type) {
    case PType::sing:
        return lifeTime_singlet;
    case PType::trip:
        return lifeTime_triplet;
    }
}
