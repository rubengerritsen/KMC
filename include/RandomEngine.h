/***************************************************
 *
 * KMC MODEL FOR OPTOELECTRIC PROCESSES
 *
 * Author: Ruben Gerritsen
 *
 * Created on 14-01-2020
 *
 * This class contains all random stuff. It contains
 * the distributions for the DOSes of the different 
 * particles, the distrubtions to generate a time step
 * etc.
 *
 **************************************************/
#pragma once
#include <random> 
#include <array>
#include "PType.h"

class RandomEngine {
public:
    RandomEngine(int seed) { rng = std::mt19937_64(seed); }
    void initializeParameters(std::array<double, 4> mu, std::array<double, 4> sigma);
    void setNrOfSites(int nr) { siteDist = std::uniform_int_distribution<int>(0,nr-1); }
    double getDOSEnergy(PType type) { return dos[type](rng); }
    double getUniform01() { return uniform01(rng); }
    int getRandomSite() { return siteDist(rng); }
    double getInterArrivalTime(double rate) { return -(1.0 / rate) * log(uniform01(rng)); }

private:
    std::mt19937_64 rng;
    std::array<std::normal_distribution<double>, 4> dos{std::normal_distribution<double>(0.0,1)};
    std::uniform_real_distribution<double> uniform01 { 0.0, 1.0 };
    std::uniform_int_distribution<int> siteDist{ 0, 10 };
};