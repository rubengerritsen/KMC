#include "Neighbourlist.h"
#include "EnumNames.h"
#include "RateEngine2.h"
#include <fstream>
#include <iostream>

void Neighbourlist::setupShortRangeNeighbours(std::string filename,
                                              const Topology &topol) {
  std::ifstream file(filename);
  if (file.is_open()) {
    RateEngine2 rates(topol.getEField(), topol.getKBT(), topol.getRateOptions());
    int id1, id2;
    double dx, dy, dz;
    double j2_e, j2_h;
    Eigen::Vector3d temp_dr;
    Neighbour temp1;
    Neighbour temp2;
    while (file >> id1 >> id2 >> dx >> dy >> dz >> j2_e >> j2_h) {
      temp_dr << dx, dy, dz;
      temp1.dr = temp_dr;
      temp2.dr = -temp_dr;
      temp1.nb = id2;
      temp2.nb = id1;
      temp1.rate_e = rates.marcusRate(
          j2_e, topol.getLambda(id1, id2, PType::elec),
          topol.getDeltaEnergy(id1, id2, PType::elec), temp_dr[0], PType::elec);
      temp2.rate_e =
          rates.marcusRate(j2_e, topol.getLambda(id2, id1, PType::elec),
                           topol.getDeltaEnergy(id2, id1, PType::elec),
                           -temp_dr[0], PType::elec);
      temp1.rate_h = rates.marcusRate(
          j2_h, topol.getLambda(id1, id2, PType::hole),
          topol.getDeltaEnergy(id1, id2, PType::hole), temp_dr[0], PType::hole);
      temp2.rate_h =
          rates.marcusRate(j2_h, topol.getLambda(id2, id1, PType::hole),
                           topol.getDeltaEnergy(id2, id1, PType::hole),
                           -temp_dr[0], PType::hole);
      // SINGLET GENERATION RATES
      temp1.rate_h_s = rates.singletGeneration(temp_dr, topol.getDeltaEnergy(id1,id2, PType::hole), topol.getSingletBindingEnergy(id2), PType::hole) ;
      temp2.rate_h_s = rates.singletGeneration(-temp_dr, topol.getDeltaEnergy(id2,id1, PType::hole), topol.getSingletBindingEnergy(id1), PType::hole) ;
      temp1.rate_e_s = rates.singletGeneration(temp_dr, topol.getDeltaEnergy(id1,id2, PType::elec), topol.getSingletBindingEnergy(id2), PType::elec) ;
      temp2.rate_e_s = rates.singletGeneration(-temp_dr, topol.getDeltaEnergy(id2,id1, PType::elec), topol.getSingletBindingEnergy(id1), PType::elec) ;
      // SINGLET DISSOCIATION RATES S->CT
      temp1.rate_s_ct_e = rates.singletDissociationCT(temp_dr, topol.getDeltaEnergy(id1,id2, PType::elec), topol.getSingletBindingEnergy(id2), PType::elec);
      temp2.rate_s_ct_e = rates.singletDissociationCT(-temp_dr, topol.getDeltaEnergy(id2,id1, PType::elec), topol.getSingletBindingEnergy(id1), PType::elec);
      temp1.rate_s_ct_h = rates.singletDissociationCT(temp_dr, topol.getDeltaEnergy(id1,id2, PType::hole), topol.getSingletBindingEnergy(id2), PType::hole);
      temp2.rate_s_ct_h = rates.singletDissociationCT(-temp_dr, topol.getDeltaEnergy(id2,id1, PType::hole), topol.getSingletBindingEnergy(id1), PType::hole);
      // CT DISSOCIATION 
      temp1.rate_ct_e = rates.ctDissociation(temp_dr, topol.getDeltaEnergy(id1, id2, PType::elec), PType::elec);
      temp2.rate_ct_e = rates.ctDissociation(-temp_dr, topol.getDeltaEnergy(id2, id1, PType::elec), PType::elec);
      temp1.rate_ct_h = rates.ctDissociation(temp_dr, topol.getDeltaEnergy(id1, id2, PType::hole), PType::hole);
      temp2.rate_ct_h = rates.ctDissociation(-temp_dr, topol.getDeltaEnergy(id2, id1, PType::hole), PType::hole);

      sRNeighbours[id1].push_back(temp1);
      sRNeighbours[id2].push_back(temp2);
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
  RateEngine2 rates(topol.getEField(), topol.getKBT(), topol.getRateOptions());
  if (file.is_open()) {
    int id1, id2;
    double dx, dy, dz;
    double j2_s;
    Eigen::Vector3d temp_dr;
    Neighbour temp1;
    Neighbour temp2;
    while (file >> id1 >> id2 >> dx >> dy >> dz >> j2_s) {
      temp_dr << dx, dy, dz;
      temp1.nb = id2;
      temp2.nb = id1;
      temp1.dr = temp_dr;
      temp2.dr = -temp_dr;

      temp1.rate_s = rates.marcusRate(
          j2_s, topol.getLambda(id1, id2, PType::sing),
          topol.getDeltaEnergy(id1, id2, PType::sing), temp_dr[0], PType::sing);
      temp2.rate_s =
          rates.marcusRate(j2_s, topol.getLambda(id2, id1, PType::sing),
                           topol.getDeltaEnergy(id2, id1, PType::sing),
                           -temp_dr[0], PType::sing);

      lRNeighbours[id1].push_back(temp1);
      lRNeighbours[id2].push_back(temp2);
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
    std::cout << nb.nb << " ";
  }
  std::cout << "\nsRRates: ";
  for (auto const &nb : sRNeighbours[site]) {
    std::cout << nb.rate_e << " ";
  }
  std::cout << std::endl;
}
