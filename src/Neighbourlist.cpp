#include "Neighbourlist.h"
#include "EnumNames.h"
#include "RateEngine2.h"
#include <fstream>
#include <iostream>

void Neighbourlist::setupNeighbours(std::string filename,
                                    const Topology &topol) {
  std::ifstream file(filename);
  if (file.is_open()) {
    int id1, id2;
    double dx, dy, dz;
    double rate;
    Eigen::Vector3d temp_dr;
    Neighbour temp1;
    while (file >> id1 >> id2 >> rate) {
      temp1.dr = topol.getCoordinates(id2) - topol.getCoordinates(id1);
      temp1.nb = id2;
      temp1.rate_e = rate;
      temp1.rate_h = rate;
      Neighbours[id1].push_back(temp1);
    }
  } else {
    std::cout << "Unable to open short range neighbour file: " << filename
              << std::endl;
    std::cout << "Terminating execution." << std::endl;
    exit(EXIT_FAILURE);
  }
}

void Neighbourlist::printNeighboursOf(int site) {
  std::cout << "Neighbours: ";
  for (auto const &nb : Neighbours[site]) {
    std::cout << nb.nb << " ";
  }
  std::cout << "\nsRRates: ";
  for (auto const &nb : Neighbours[site]) {
    std::cout << nb.rate_h << " ";
  }
  std::cout << std::endl;
}
