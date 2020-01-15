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

class RateEngine {
public:
    double millerAbrahams(const Site& siteOne, const Site& siteTwo, const PType type) const;
    double dexter(const Site& siteOne, const Site& siteTwo) const;


private:
    std::array<double, 4> v0;
    std::array<double, 4> alpha;
    std::array<double, 4> charge;
    double E_Field;
    double kBT;



};