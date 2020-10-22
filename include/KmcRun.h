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
#include "EnumNames.h"
#include "Neighbourlist.h"
#include "NextEventList.h"
#include "OutputManager.h"
#include "RandomEngine.h"
#include "RateEngine2.h"
#include "SimulationOptions.h"
#include "Topology.h"
#include <iostream>
#include <numeric>
#include <string>

class KmcRun {
public:
  KmcRun(const Topology &topol, const Neighbourlist &nbList,
         const SimulationOptions simOptions)
      : topol(topol), nbList(nbList),
        rate_engine(simOptions.EField_x, simOptions.kBT,
                    topol.getRateOptions()),
        simOptions(simOptions), random_engine(simOptions.SEED),
        out(simOptions) {
    int totalNrOfParticles;
    totalNrOfParticles = simOptions.nrOfElectrons + simOptions.nrOfHoles +
                         simOptions.nrOfSinglets;
    // Create space for at least a 100 next events
    next_event_list.initializeListSize(totalNrOfParticles * 100);
  }
  void runSimulation();

private:
  const Topology &topol;
  const Neighbourlist &nbList;
  SimulationOptions simOptions;
  RateEngine2 rate_engine;
  RandomEngine random_engine;
  NextEventList next_event_list{};
  OutputManager out;

  /* Storage for the sites and particles */
  std::vector<Site> siteList;
  std::vector<Particle> particleList;

  /* time */
  double totalTime = 0;

  /* Helper functions */
  void initializeSites();
  void initializeNeighbours();
  void initializeParticles();
  void computeNextEventList();
  bool executeNextEvent();
};