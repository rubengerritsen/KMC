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

#include "RateEngine.h"

double RateEngine::millerAbrahams(const Site& siteOne, const Site& siteTwo, const PType type) const {
    double deltaE = siteTwo.getEnergy(type) - siteOne.getEnergy(type) + E_Field * charge[type] * siteTwo.xDistFrom(siteOne);
    if ( deltaE <= 0) {
        return(v0[type] * exp(-2 * alpha[type] * 3));
    }
    else {
        return(v0[type] * exp(-2 * alpha[type] * 3 - deltaE / kBT));
    }
}
