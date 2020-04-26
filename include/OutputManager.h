#include "Particle.h"
#include <boost/format.hpp>
#include <vector>
#include "Topology.h"

class OutputManager {
public:
  OutputManager(int simID) : simID(simID) {
    particlePathFile = (boost::format("%d_particlePaths.txt") % simID).str();
    stateFile = (boost::format("%d_state.txt") % simID).str();
    numberFile = (boost::format("%d_numbers.txt") % simID).str();
  };
  void registerParticlePositions(const std::vector<Particle> &particleList,
                                 double time);
  void registerState(const std::vector<Particle> &particleList, double time);
  void registerNumbers(const std::vector<Particle> &particleList,
                       const Topology &topol, double time);

private:
  int simID;
  std::string particlePathFile;
  std::string stateFile;
  std::string numberFile;
  bool firstTimePath = true;
  bool firstTimeState = true;
  bool firstTimeNumber = true;
};