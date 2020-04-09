#pragma once
#include "Neighbour.h"
#include "Topology.h"
#include <Eigen/Dense>
#include <string>
#include <vector>

class Neighbourlist {
public:
  Neighbourlist(int numberOfSites) : numberOfSites{numberOfSites} {
    sRNeighbours.resize(numberOfSites);
    lRNeighbours.resize(numberOfSites);
  }

  void setupShortRangeNeighbours(std::string filename, const Topology &topol);
  void setupLongRangeNeighbours(std::string filename, const Topology &topol);
  void printNeighboursOf(int site);

  const std::vector<Neighbour> &getSRNeighbours(int site) const {
    return sRNeighbours[site];
  }
  const std::vector<Neighbour> &getLRNeighbours(int site) const {
    return lRNeighbours[site];
  }

private:
  int numberOfSites;
  std::vector<std::vector<Neighbour>> sRNeighbours;
  std::vector<std::vector<Neighbour>> lRNeighbours;
};