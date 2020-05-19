#pragma once
#include "Particle.h"
#include "SimulationOptions.h"
#include "NextEventList.h"
#include "Neighbourlist.h"
#include "Topology.h"
#include <boost/format.hpp>
#include <vector>

class OutputManager {
public:
  OutputManager(SimulationOptions simOptions)
      : simID(simOptions.simID), simOptions(simOptions) {
    particlePathFile = simOptions.outputPath +
                       (boost::format("%d_particlePaths.txt") % simID).str();
    stateFile =
        simOptions.outputPath + (boost::format("%d_state.txt") % simID).str();
    numberFile =
        simOptions.outputPath + (boost::format("%d_numbers.txt") % simID).str();
    siteFile =
        simOptions.outputPath + (boost::format("%d_sites.txt") % simID).str();
    rateFile = simOptions.outputPath + "rates.txt";
    nextEventFile = simOptions.outputPath + "nextEvents.txt";
  };

  void registerParticlePositions(const std::vector<Particle> &particleList,
                                 double time);
  void printRatesToFile(const Neighbourlist &nblist);
  void registerState(const std::vector<Particle> &particleList, double time);
  void registerNumbers(const std::vector<Particle> &particleList,
                       const Topology &topol, double time);
  void printNextEventList(NextEventList &nextEventList);

private:
  int simID;
  SimulationOptions simOptions;
  std::string particlePathFile;
  std::string stateFile;
  std::string numberFile;
  std::string siteFile;
  std::string rateFile;
  std::string nextEventFile;
  bool firstTimePath = true;
  bool firstTimeState = true;
  bool firstTimeNumber = true;
};