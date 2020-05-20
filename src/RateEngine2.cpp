#include "RateEngine2.h"

double RateEngine2::marcusRate(double jeff2, double lambda_ij, double deltaE,
                               double dx, PType type) const {
  // returns hopping rate from i->j.
  double totalEnergy;
  if (type == PType::hole) {
    totalEnergy = -deltaE + lambda_ij + EField_x * charge[type] * dx;
  } else {
    totalEnergy = deltaE + lambda_ij + EField_x * charge[type] * dx;
  }
  return 2 * constants::pi / constants::hbar * jeff2 /
         std::sqrt(4 * constants::pi * lambda_ij * kBT) *
         std::exp(-totalEnergy * totalEnergy / (4 * lambda_ij * kBT));
}

double RateEngine2::singletDecay(double singletEnergy, MType type) const {
  return sqrtDielectric * (4.0 / 3.0) * constants::alpha *
         std::pow(singletEnergy, 3) / (constants::c2 * constants::hbar3) *
         mu2[type];
}

double RateEngine2::millerAbrahams(Eigen::Vector3d dr, double deltaE,
                                   PType type) const {
  double totalEnergy;
  if (type == PType::hole) {
    totalEnergy = -deltaE + EField_x * charge[type] * dr[0];
  } else {
    totalEnergy = deltaE + EField_x * charge[type] * dr[0];
  }
  if (totalEnergy > 0) {
    return attempt[type] *
           std::exp(-2.0 * alpha[type] * dr.norm() - totalEnergy / kBT);
  } else {
    return attempt[type] * std::exp(-2.0 * alpha[type] * dr.norm());
  }
}

double RateEngine2::singletGeneration(Eigen::Vector3d dr, double deltaE,
                                      double singletBinding, PType type) const {
  double totalEnergy;
  if (type == PType::hole) {
    totalEnergy = -deltaE + EField_x * charge[type] * dr[0] - singletBinding;
  } else {
    totalEnergy = deltaE + EField_x * charge[type] * dr[0] - singletBinding;
  }
  if (totalEnergy > 0) {
    return attempt[type] *
           std::exp(-2.0 * alpha[type] * dr.norm() - totalEnergy / kBT);
  } else {
    return attempt[type] * std::exp(-2.0 * alpha[type] * dr.norm());
  }
}

double RateEngine2::singletDissociationCT(Eigen::Vector3d dr, double deltaE,
                                          double singletBinding,
                                          PType type) const {
  double totalEnergy;
  double CTBinding =
      constants::hrt2ev * 1.0 / (dielectricConstant * constants::nm2bohr * dr.norm());
  if (type == PType::hole) {
    totalEnergy =
        -deltaE + EField_x * charge[type] * dr[0] + singletBinding - CTBinding;
  } else {
    totalEnergy =
        deltaE + EField_x * charge[type] * dr[0] + singletBinding - CTBinding;
  }
  if (totalEnergy > 0) {
    return attempt[type] *
           std::exp(-2.0 * alpha[type] * dr.norm() - totalEnergy / kBT);
  } else {
    return attempt[type] * std::exp(-2.0 * alpha[type] * dr.norm());
  }
}

double RateEngine2::ctDissociation(Eigen::Vector3d dr, double deltaE,
                                   PType type) const {
  double totalEnergy;
  double CTBinding =
      constants::hrt2ev * 1.0 / (dielectricConstant  * constants::nm2bohr * dr.norm());
  if (type == PType::hole) {
    totalEnergy =
        -deltaE + EField_x * charge[type] * dr[0] + CTBinding;
  } else {
    totalEnergy =
        deltaE + EField_x * charge[type] * dr[0] + CTBinding;
  }
  if (totalEnergy > 0) {
    return attempt[type] *
           std::exp(-2.0 * alpha[type] * dr.norm() - totalEnergy / kBT);
  } else {
    return attempt[type] * std::exp(-2.0 * alpha[type] * dr.norm());
  }
}
