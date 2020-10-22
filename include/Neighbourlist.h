#pragma once
#include "Neighbour.h"
#include "Topology.h"
#include <Eigen/Dense>
#include <string>
#include <vector>

class Neighbourlist {
public:
  Neighbourlist(int numberOfSites) : numberOfSites{numberOfSites} {
    Neighbours.resize(numberOfSites);
  }

  void setupNeighbours(std::string filename, const Topology &topol);
  void printNeighboursOf(int site);

  const std::vector<Neighbour> &getNeighbours(int site) const {
    return Neighbours[site];
  }

private:
  int numberOfSites;
  std::vector<std::vector<Neighbour>> Neighbours;
};