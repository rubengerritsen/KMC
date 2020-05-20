#include "OutputManager.h"
#include "Neighbour.h"
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

void OutputManager::printRatesToFile(const Neighbourlist &nbList) {
  std::ofstream outFile;
  outFile.open(rateFile);
  if (outFile.is_open()) {
    outFile << boost::format(
        "rate_e rate_h rate_s rate_h_s rate_e_s rate_s_ct_e rate_s_ct_h rate_ct_e rate_ct_h dist\n");
    outFile << boost::format("shortRangNeighbours\n");
    for (int i = 0; i < nbList.getNumberOfSites(); i++) {
      for (auto &nb : nbList.getSRNeighbours(i)) {
        outFile << boost::format("%.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e\n") %
                       nb.rate_e % nb.rate_h % nb.rate_s % nb.rate_h_s %
                       nb.rate_e_s % nb.rate_s_ct_e % nb.rate_s_ct_h % nb.rate_ct_e % nb.rate_ct_h % nb.dr.norm();
      }
    }
    std::cout << "printed sr neighbours \n";
    outFile << boost::format("longRangNeighbours\n");
  /*  for (int i = 0; i < nbList.getNumberOfSites(); i++) {
      for (auto &nb : nbList.getLRNeighbours(i)) {
        outFile << boost::format("%.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e\n") %
                       nb.rate_e % nb.rate_h % nb.rate_s % nb.rate_h_s %
                       nb.rate_e_s % nb.rate_s_ct_e % nb.rate_s_ct_h % nb.rate_ct_e % nb.rate_ct_h % nb.dr.norm();
      }
    } */
    std::cout << "printed lr neighbours \n";
  } else {
    std::cout << "Failed to open rateFile for simID: " << simID << "\n";
  };
}

void OutputManager::printNextEventList(NextEventList &nextEventList){
  std::ofstream outFile;
  outFile.open(nextEventFile);
  if (outFile.is_open()) {
    const std::vector<Transition> &types = nextEventList.getEventTypes();
    const std::vector<double> &rates = nextEventList.getRates();

    for (int i = 0 ; i<rates.size(); i++){
      outFile << boost::format("%d %.5e\n") % types[i] % rates[i];
    }
  }else{
    std::cout << "Failed to open nextEvent Output file for simID: " << simID << "\n";
  }
}

void OutputManager::registerState(const std::vector<Particle> &particleList,
                                  double time) {

  std::ofstream outFile;
  if (firstTimeState) {
    outFile.open(stateFile);
    firstTimeState = false;
  } else {
    outFile.open(stateFile, std::fstream::app);
  }

  if (outFile.is_open()) {
    outFile << boost::format("%12.5f") % time << " ";
    for (auto &part : particleList) {
      if (part.isAlive()) {
        outFile << boost::format("%d %d    ") % part.getType() %
                       part.getLocation();
        if (part.getType() == PType::CT) {
          outFile << boost::format("%d %d    ") % part.getType() %
                         part.getLocationCTelec();
        }
      }
    }
    outFile << "\n";
  } else {
    std::cout << "Could not open output file for simID: " << simID << "\n";
  }
}

void OutputManager::registerNumbers(const std::vector<Particle> &particleList,
                                    const Topology &topol, double time) {

  std::ofstream outFileNR;
  if (firstTimeNumber) {
    outFileNR.open(numberFile);
    firstTimeNumber = false;
  } else {
    outFileNR.open(numberFile, std::fstream::app);
  }

  std::array<int, 12> count{0};

  if (outFileNR.is_open()) {
    outFileNR << boost::format("%12.5f") % time << " ";
    for (auto &part : particleList) {
      if (part.isAlive()) {
        count[part.getType()] += 1;
        if (topol.getMolType(part.getLocation()) == MType::ben) {
          count[part.getType() + 4] += 1;
        } else if (topol.getMolType(part.getLocation()) == MType::tcne) {
          count[part.getType() + 8] += 1;
        } else {
          std::cout << "Unknown molecule type in outputmanager!\n";
          exit(EXIT_FAILURE);
        }
      }
    }

    for (auto &nr : count) {
      outFileNR << boost::format("%5d ") % nr;
    }
    outFileNR << "\n";
  }
}