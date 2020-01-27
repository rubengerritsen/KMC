/***************************************************
 *
 * KMC MODEL FOR OPTOELECTRIC PROCESSES
 *
 * Author: Ruben Gerritsen
 *
 * Created on 14-01-2020
 *
 * This class contains all functions necessary to
 * perform a single monte-carlo simulation.
 **************************************************/
#pragma once
#include <string>
#include <numeric>
#include <iostream>
#include "RateEngine.h"
#include "PBC.h"
#include "RandomEngine.h"
#include "NextEventList.h"
#include "PType.h"
#include "OutputManager.h"

class KmcRun {
public:
    KmcRun(RateEngine rate_engine, PBC pbc, RandomEngine random_engine, int nrOfSteps, std::array<int,4> qt, std::string siteFile, double sR_CutOff, double lR_CutOff) :
        rate_engine(rate_engine), pbc(pbc), random_engine(random_engine), nrOfSteps(nrOfSteps), nrOfParticlesPerType(qt) ,siteFile(siteFile), sR_cutOff(sR_CutOff), lR_cutOff(lR_CutOff) {
            int totalNrOfParticles;
            totalNrOfParticles = std::accumulate(nrOfParticlesPerType.begin(), nrOfParticlesPerType.end(), totalNrOfParticles);
            next_event_list.initializeListSize(totalNrOfParticles * 100); // create space for at least a 100 events per particles
            std::cout << "Initial number of particles in the simulation: " << totalNrOfParticles << "\n";
        }
    void runSimulation();


private:
    RateEngine rate_engine;
    PBC pbc;
    RandomEngine random_engine;
    NextEventList next_event_list {} ;

    /* Storage for the graph and particles */
    std::vector<Site> siteList;
    std::vector<Particle> particleList;

    int nrOfSteps;
    std::string siteFile;
    std::string outputFile = "./output/occ.txt";

    /* Some additional model parameters */
    double lR_cutOff;
    double sR_cutOff;
    double totalTime = 0.0;
    std::array<int,4> nrOfParticlesPerType;

    /* Helper functions */
    void initializeSites();
    void initializeNeighbours();
    void initializeParticles();
    void computeNextEventRates();
    void executeNextEvent();
};