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
#include "EnumNames.h"

class RandomEngine {
public:
    RandomEngine(int seed) { rng = std::mt19937_64(seed); }
    void setNrOfSites(int nr) { siteDist = std::uniform_int_distribution<int>(0,nr-1); }
    double getUniform01() { return uniform01(rng); }
    int getRandomSite() { return siteDist(rng); }
    double getInterArrivalTime(double rate) { return -(1.0 / rate) * log(uniform01(rng)); }

private:
    std::mt19937_64 rng;
    std::array<std::normal_distribution<double>, 4> dos;
    std::uniform_real_distribution<double> uniform01 { 0.0, 1.0 };
    std::uniform_int_distribution<int> siteDist{ 0, 10 };
};