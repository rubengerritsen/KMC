#include "Neighbourlist.h"
#include "EnumNames.h"
#include <boost/math/constants/constants.hpp>
#include <fstream>
#include <iostream>

void Neighbourlist::setupShortRangeNeighbours(std::string filename,
                                              const Topology &topol) {
  std::ifstream file(filename);
  if (file.is_open()) {
    int id1, id2;
    double dx, dy, dz;
    double j2_e, j2_h;
    Eigen::Vector3d temp_dr;
    while (file >> id1 >> id2 >> dx >> dy >> dz >> j2_e >> j2_h) {
      temp_dr << dx, dy, dz;
      dr_sR_list[id1].push_back(temp_dr);
      dr_sR_list[id2].push_back(-temp_dr);
      sRNeighbours[id1].push_back(id2);
      sRNeighbours[id2].push_back(id1);
      coupling_e[id1].push_back(j2_e);
      coupling_e[id2].push_back(j2_e);
      coupling_h[id1].push_back(j2_h);
      coupling_h[id2].push_back(j2_h);
    }
  } else {
    std::cout << "Unable to open short range neighbour file: " << filename
              << std::endl;
    std::cout << "Terminating execution." << std::endl;
    exit(EXIT_FAILURE);
  }
}

void Neighbourlist::setupLongRangeNeighbours(std::string filename,
                                             const Topology &topol) {
  std::ifstream file(filename);
  if (file.is_open()) {
    int id1, id2;
    double dx, dy, dz;
    double j2_s;
    Eigen::Vector3d temp_dr;
    while (file >> id1 >> id2 >> dx >> dy >> dz >> j2_s) {
      temp_dr << dx, dy, dz;
      dr_lR_list[id1].push_back(temp_dr);
      dr_lR_list[id2].push_back(-temp_dr);
      lRNeighbours[id1].push_back(id2);
      lRNeighbours[id2].push_back(id1);
      coupling_s[id1].push_back(j2_s);
      coupling_s[id2].push_back(j2_s);
    }
  } else {
    std::cout << "Unable to open long range neighbour file: " << filename
              << std::endl;
    std::cout << "Terminating execution." << std::endl;
    exit(EXIT_FAILURE);
  }
}

void Neighbourlist::printNeighboursOf(int site) {
  std::cout << "sRNeighbours: ";
  for (auto const &nb : sRNeighbours[site]) {
    std::cout << nb << " ";
  }
  std::cout << "\nlRNeighbours: ";
  for (auto const &nb : lRNeighbours[site]) {
    std::cout << nb << " ";
  }
  std::cout << std::endl;
}

double Neighbourlist::marcusRate(double jeff2, double lambda_ij, double kBT,
                                 double Energy1, double Energy2) {
  // note time unit of hbar: ps => hbar = 6.582119514e-4
  return 2 * boost::math::constants::pi<double>() / (6.582119514e-4) * jeff2 /
         std::sqrt(4 * boost::math::constants::pi<double>() * lambda_ij * kBT) * std::exp(- (Energy2 - Energy1 - lambda_ij) * (Energy2 - Energy1 - lambda_ij)  / (4*lambda_ij *kBT));
}