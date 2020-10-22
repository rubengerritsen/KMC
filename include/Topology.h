#pragma once
#include "EnumNames.h"
#include "RateOptions.h"
#include <Eigen/Dense>
#include <array>
#include <string>
#include <vector>

class Topology {
public:
  Topology() = default;
  void readTopologyFromFile(std::string filename);
  void readReorganisationEnergies(std::string filename);
  /* prints first few items of topology class */
  void printHead(int nrOfItems = 10) const;
  void printReorganisationEnergies();
  int getNrOfSites() const { return siteLocations.size(); }

  double getKBT() const { return kBT; }
  double getEField() const { return EField_x; }
  double getLambda(int id1, int id2, PType type) const {
    return lambdaXtoN[moleculeType[id1]][type] +
           lambdaNtoX[moleculeType[id2]][type];
  }
  double getDeltaEnergy(int id1, int id2, PType type) const {
    return siteEnergies[id2][type] - siteEnergies[id1][type];
  }
  MType getMolType(int id) const { return moleculeType[id]; }

  double getEnergy(int id, PType type) const { return siteEnergies[id][type]; }
  double getSingletBindingEnergy(int id) const {
    return siteEnergies[id][PType::elec] - siteEnergies[id][PType::hole] -
           siteEnergies[id][PType::sing];
  }
  Eigen::Vector3d getCoordinates(int id) const { return siteLocations[id]; }

  RateOptions getRateOptions() const { return rOptions; }
  void setRateOptions(RateOptions options) { rOptions = options; }

private:
  std::vector<Eigen::Vector3d> siteLocations;
  std::vector<std::vector<double>> siteEnergies;
  std::vector<MType> moleculeType;
  std::array<std::array<double, 3>, 2> lambdaXtoN;
  std::array<std::array<double, 3>, 2> lambdaNtoX;
  RateOptions rOptions;
  double EField_x = 0;
  double kBT = 0.1;
};