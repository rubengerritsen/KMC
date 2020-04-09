#include "OutputManager.h"
#include <fstream>

void OutputManager::registerParticlePositions(
    const std::vector<Particle> &particleList, double time) {

  std::ofstream outFile;
  if (firstTimePath) {
    outFile.open(particlePathFile);
    firstTimePath = false;
  } else {
    outFile.open(particlePathFile, std::fstream::app);
  }

  if (outFile.is_open()) {
    for (auto &part : particleList) {
      outFile << boost::format("%4d ") % part.getLocation();
    }
    outFile << "\n";
  } else {
    std::cout << "Could not open output file for simID: " << simID << "\n";
  }
}

void OutputManager::registerState(
    const std::vector<Particle> &particleList, double time) {

  std::ofstream outFile;
  if (firstTimeState) {
    outFile.open(stateFile);
    firstTimeState = false;
  } else {
    outFile.open(stateFile, std::fstream::app);
  }

  if (outFile.is_open()) {
      outFile << time << " ";
    for (auto &part : particleList) {
      outFile << boost::format("%d:%4d ") % part.getType() % part.getLocation();
    }
    outFile << "\n";
  } else {
    std::cout << "Could not open output file for simID: " << simID << "\n";
  }
}