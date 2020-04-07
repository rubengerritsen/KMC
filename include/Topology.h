#pragma once
#include <Eigen/Dense>
#include <array>
#include <string>
#include <vector>
#include "EnumNames.h"

class Topology {
public:
  Topology() {}
  void readTopologyFromFile(std::string filename);
  void readReorganisationEnergies(std::string filename);
  /* prints first few items of topology class */
  void printHead(int nrOfItems = 10);
  void printReorganisationEnergies();
  int getNrOfSites() { return moleculeType.size(); }

private:
  std::vector<Eigen::Vector3d> siteLocations;
  std::vector<std::vector<double>> siteEnergies;
  std::vector<MType> moleculeType;
  std::array<std::array<double,3>,2> lambdaXtoN;
  std::array<std::array<double,3>,2> lambdaNtoX;
};