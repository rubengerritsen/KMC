#pragma once
#include <Eigen/Dense>
#include <string>
#include <vector>
#include "Topology.h"

class Neighbourlist {
public:
  Neighbourlist(int numberOfSites) : numberOfSites{numberOfSites} {
    sRNeighbours.resize(numberOfSites);
    lRNeighbours.resize(numberOfSites);
    dr_sR_list.resize(numberOfSites);
    dr_lR_list.resize(numberOfSites);
    rate_e.resize(numberOfSites);
    rate_h.resize(numberOfSites);
    rate_s.resize(numberOfSites);
  }

  void setupShortRangeNeighbours(std::string filename, const Topology &topol);
  void setupLongRangeNeighbours(std::string filename, const Topology &topol);
  void printNeighboursOf(int site);

  const std::vector<int> &getSRNeighbours(int site) const {
    return sRNeighbours[site];
  }
  const std::vector<int> &getLRNeighbours(int site) const {
    return lRNeighbours[site];
  }

private:

  
  int numberOfSites;
  std::vector<std::vector<int>> sRNeighbours;
  std::vector<std::vector<int>> lRNeighbours;
  std::vector<std::vector<Eigen::Vector3d>> dr_sR_list;
  std::vector<std::vector<Eigen::Vector3d>> dr_lR_list;
  std::vector<std::vector<double>> rate_e;
  std::vector<std::vector<double>> rate_h;
  std::vector<std::vector<double>> rate_s;
};