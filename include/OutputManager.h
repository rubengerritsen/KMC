#include "Particle.h"
#include <boost/format.hpp>
#include <vector>

class OutputManager {
public:
  OutputManager(int simID) : simID(simID){
      particlePathFile = (boost::format("%d_particlePaths.txt") % simID).str();
      stateFile =  (boost::format("%d_state.txt") % simID).str();
  };
  void registerParticlePositions(const std::vector<Particle> &particleList, double time);
  void registerState(const std::vector<Particle> &particleList, double time);

private:
  int simID;
  std::string particlePathFile;
  std::string stateFile;
  bool firstTimePath = true;
  bool firstTimeState = true;
};