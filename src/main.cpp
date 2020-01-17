/***************************************************
 * 
 * KMC MODEL FOR OPTOELECTRIC PROCESSES
 *
 * Author: Ruben Gerritsen
 *
 * Created on 14-01-2020
 *
 **************************************************/
#include <iostream>
#include <fstream>
#include <string>
#include <array>
#include "Particle.h"
#include "Site.h"
#include "PBC.h"
#include "RateEngine.h"
#include "PType.h"
#include "RandomEngine.h"
#include "SingleRun.h"


void setupAndExecuteSimulation() {
    int SEED{0},nrOfSteps{ 0 };
    double Xmax{0}, Ymax{ 0 }, Zmax{ 0 };
    double sR_CutOff{ 0 }, lR_CutOff{ 0 };
    std::array<int, 4> qt;
    std::array<double, 4> v0;
    std::array<double, 4> alpha;
    std::array<double, 4> charge;
    std::array<double, 4> DOS_mu;
    std::array<double, 4> DOS_sigma;
    double kBT{ 0 };
    double E_Field{ 0 };


    /* Reading all model parameters */
    std::string paramFile = "./input/modelParameters.txt";
    std::ifstream myfile("./input/sites.txt");
    std::string junk;
    if (myfile.is_open()) {
        myfile >> junk >> SEED;
        myfile >> junk >> nrOfSteps;
        myfile >> junk >> sR_CutOff;
        myfile >> junk >> lR_CutOff;
        myfile >> junk >> Xmax;
        myfile >> junk >> Ymax;
        myfile >> junk >> Zmax;
        for (unsigned int i = 0; i < 4; ++i) {
            myfile >> junk >> qt[i];
            myfile >> junk >> DOS_mu[i];
            myfile >> junk >> DOS_sigma[i];
            myfile >> junk >> charge[i];
            myfile >> junk >> v0[i];
            myfile >> junk >> alpha[i];
        }
        myfile >> junk >> kBT;
        myfile >> junk >> E_Field;
    }
    else {
        std::cout << "Unable to open file: " << paramFile << std::endl;
        std::cout << "Terminating execution." << std::endl;
        exit(EXIT_FAILURE);
    }

    /* Setting up helper objects for the simulation */
    PBC pbc(Xmax, Ymax, Zmax);
    RateEngine rate_engine(v0, alpha, charge, E_Field, kBT, pbc);
    RandomEngine random_engine(SEED);

    /* Execution of the experiment*/
    SingleRun experiment{rate_engine, pbc, random_engine, nrOfSteps, "./input/sites.txt" };
    experiment.runSimulation();
}

int main() {
	
    setupAndExecuteSimulation();

    return 0;
}


