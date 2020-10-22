#include "Particle.h"
#include "SimulationOptions.h"
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
    tofFile = simOptions.outputPath + "tof.txt";
  };
  void registerParticlePositions(const std::vector<Particle> &particleList,
                                 double time);
  void registerState(const std::vector<Particle> &particleList, double time);
  void registerNumbers(const std::vector<Particle> &particleList,
                       const Topology &topol, double time);

  void registerTOF(const double time);
  void endTOF();

private:
  int simID;
  SimulationOptions simOptions;
  std::string particlePathFile;
  std::string stateFile;
  std::string numberFile;
  std::string siteFile;
  std::string tofFile;
  bool firstTimePath = true;
  bool firstTimeState = true;
  bool firstTimeNumber = true;
};