#pragma once
#include <string>
#include <vector>
#include <Eigen/Dense>

class Neighbourlist {
public:
  Neighbourlist(int numberOfSites): numberOfSites{numberOfSites} {
      sRNeighbours.resize(numberOfSites);
      lRNeighbours.resize(numberOfSites);
      dr_sR_list.resize(numberOfSites);
      dr_lR_list.resize(numberOfSites);
      coupling_e.resize(numberOfSites);
      coupling_h.resize(numberOfSites);
      coupling_s.resize(numberOfSites);
  }

  void readShortRangeNeighboursFromFile(std::string filename);
  void readLongRangeNeighboursFromFile(std::string filename);
  void printNeighboursOf(int site);

  void preComputeRates();

  const std::vector<int>& getSRNeighbours(int site) const {return sRNeighbours[site];}
  const std::vector<int>& getLRNeighbours(int site) const {return lRNeighbours[site];}



private:
    int numberOfSites;
    std::vector<std::vector<int>> sRNeighbours;
    std::vector<std::vector<int>> lRNeighbours;
    std::vector<std::vector<Eigen::Vector3d>>  dr_sR_list;
    std::vector<std::vector<Eigen::Vector3d>>  dr_lR_list;
    std::vector<std::vector<double>> coupling_e;
    std::vector<std::vector<double>> coupling_h;
    std::vector<std::vector<double>> coupling_s;
};