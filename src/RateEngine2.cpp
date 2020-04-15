#include "RateEngine2.h"



double RateEngine2::marcusRate(double jeff2, double lambda_ij, double deltaE,
                               double dx, PType type) {
  // returns hopping rate from i->j.
  double totalEnergy;
  if (type == PType::hole) {
    totalEnergy = deltaE - lambda_ij + EField_x * charge[type] * dx;
  } else {
    totalEnergy = -deltaE - lambda_ij + EField_x * charge[type] * dx;
  }
  return 2 * constants::pi / constants::hbar * jeff2 /
         std::sqrt(4 * constants::pi * lambda_ij * kBT) *
         std::exp(-totalEnergy * totalEnergy / (4 * lambda_ij * kBT));
}

double RateEngine2::singletDecay(double singletEnergy, MType type) {
  return sqrtDielectric * (4.0 / 3.0) * constants::alpha * std::pow(singletEnergy, 3) /
         (constants::c2 * constants::hbar3) * mu2[type]; 
}