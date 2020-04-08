#include "Neighbourlist.h"
#include "EnumNames.h"
#include "RateEngine2.h"
#include <fstream>
#include <iostream>

void Neighbourlist::setupShortRangeNeighbours(std::string filename,
                                              const Topology &topol) {
  std::ifstream file(filename);
  if (file.is_open()) {
    RateEngine2 rates(topol.getEField(), topol.getKBT());
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
      rate_e[id1].push_back(
          rates.marcusRate(j2_e, topol.getLambda(id1, id2, PType::elec),
                           topol.getDeltaEnergy(id1, id2, PType::elec),
                           temp_dr[0], PType::elec));
      rate_e[id2].push_back(
          rates.marcusRate(j2_e, topol.getLambda(id2, id1, PType::elec),
                           topol.getDeltaEnergy(id2, id1, PType::elec),
                           -temp_dr[0], PType::elec));
      rate_h[id1].push_back(
          rates.marcusRate(j2_h, topol.getLambda(id1, id2, PType::hole),
                           topol.getDeltaEnergy(id1, id2, PType::hole),
                           temp_dr[0], PType::hole));
      rate_h[id2].push_back(
          rates.marcusRate(j2_h, topol.getLambda(id2, id1, PType::hole),
                           topol.getDeltaEnergy(id2, id1, PType::hole),
                           -temp_dr[0], PType::hole));
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
  RateEngine2 rates(topol.getEField(), topol.getKBT());
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
      rate_s[id1].push_back(
          rates.marcusRate(j2_s, topol.getLambda(id1, id2, PType::sing),
                           topol.getDeltaEnergy(id1, id2, PType::sing),
                           temp_dr[0], PType::elec));
      rate_s[id2].push_back(
          rates.marcusRate(j2_s, topol.getLambda(id2, id1, PType::sing),
                           topol.getDeltaEnergy(id2, id1, PType::sing),
                           -temp_dr[0], PType::elec));
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
  std::cout << "\nsRRates: ";
  for (auto const &nb : rate_e[site]) {
    std::cout << nb << " ";
  }
  std::cout << std::endl;
}
