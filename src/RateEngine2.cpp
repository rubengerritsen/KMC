#include "RateEngine2.h"
#include <boost/math/constants/constants.hpp>

double RateEngine2::marcusRate(double jeff2, double lambda_ij, double deltaE,
                               double dx, PType type) {
  // returns hopping rate from i->j.
  // note time unit of hbar: ps => hbar = 6.582119514e-4
  double totalEnergy = deltaE - lambda_ij + EField_x * charge[type] * dx;
  return 2 * boost::math::constants::pi<double>() / (6.582119514e-4) * jeff2 /
         std::sqrt(4 * boost::math::constants::pi<double>() * lambda_ij * kBT) *
         std::exp(-totalEnergy * totalEnergy / (4 * lambda_ij * kBT));
}