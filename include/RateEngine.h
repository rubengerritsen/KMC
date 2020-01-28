/***************************************************
 *
 * KMC MODEL FOR OPTOELECTRIC PROCESSES
 *
 * Author: Ruben Gerritsen
 *
 * Created on 14-01-2020
 *
 * RateEngine is the class containing all parameters
 * and functions to compute the rates used by the KMC.
 *
 **************************************************/

#pragma once
#include <array>
#include "Site.h"
#include "Particle.h"
#include "PBC.h"
#include "EnumNames.h"

class RateEngine {
public:
    RateEngine(std::array<double, 4> v0, std::array<double, 4> alpha, std::array<double, 4> charge, double E_Field, double kBT, PBC& pbc) : 
        v0(v0), alpha(alpha), charge(charge), E_Field(E_Field), kBT(kBT), pbc(pbc) {};
    double millerAbrahams(const Site& siteOne, const Site& siteTwo, const PType type) const;
    double millerAbrahamsCT(const Site& siteOne, const Site& siteTwo, const PType type) const;
    double millerAbrahamsGEN(const Site& siteOne, const Site& siteTwo, const PType type) const;
    double millerAbrahamsDIS(const Site& siteOne, const Site& siteTwo, const PType type) const;
    double millerAbrahamsCT_DIS(const Site& siteOne, const Site& siteTwo, const PType type) const;
    double forster(const Site& siteOne, const Site& siteTwo) const;
    double decay(const PType type) const;


private:
    std::array<double, 4> v0;
    std::array<double, 4> alpha;
    std::array<double, 4> charge;
    double lifeTime_singlet = 1e-4;
    double lifeTime_triplet = 1e-4;
    double deltaE_SingtoCT = 0.48;
    double deltaE_binding = 1.0;
    double deltaE_CTbinding = 0.48;
    double E_Field;
    double kBT;
    double R_forster = 0.3;
    PBC pbc;


};