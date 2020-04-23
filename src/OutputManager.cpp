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

void OutputManager::registerState(const std::vector<Particle> &particleList,
                                  double time) {

  std::ofstream outFile;
  std::ofstream outFileNR;
  if (firstTimeState) {
    outFile.open(stateFile);
    outFileNR.open(numberFile);
    firstTimeState = false;
  } else {
    outFile.open(stateFile, std::fstream::app);
    outFileNR.open(numberFile, std::fstream::app);
  }

  if (outFile.is_open()) {
    outFile << boost::format("%12.5f") % time << " ";
    for (auto &part : particleList) {
      outFile << boost::format("%d:%d:%-4d ") % part.isAlive() %
                     part.getType() % part.getLocation();
    }
    outFile << "\n";
  } else {
    std::cout << "Could not open output file for simID: " << simID << "\n";
  }

  std::array<int, 4> count{0};

  if (outFileNR.is_open()) {
    outFileNR << boost::format("%12.5f") % time << " ";
    for (auto &part : particleList) {
      if (part.isAlive()) {
        count[part.getType()] += 1;
      }
    }

    for (auto &nr : count) {
      outFileNR << boost::format("%5d ") % nr;
    }
    outFileNR << "\n";
  }
}

void OutputManager::registerNumbers(const std::vector<Particle> &particleList,
                                  double time) {

  std::ofstream outFileNR;
  if (firstTimeNumber) {
    outFileNR.open(numberFile);
    firstTimeNumber = false;
  } else {
    outFileNR.open(numberFile, std::fstream::app);
  }

  std::array<int, 4> count{0};

  if (outFileNR.is_open()) {
    outFileNR << boost::format("%12.5f") % time << " ";
    for (auto &part : particleList) {
      if (part.isAlive()) {
        count[part.getType()] += 1;
      }
    }

    for (auto &nr : count) {
      outFileNR << boost::format("%5d ") % nr;
    }
    outFileNR << "\n";
  }
}