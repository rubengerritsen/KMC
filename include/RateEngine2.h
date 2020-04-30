#pragma once
#include "EnumNames.h"
#include <array>
#include "constants.h"
#include <Eigen/Dense>
#include "RateOptions.h"

class RateEngine2 {
public:
  RateEngine2(double EField_x, double kBT, RateOptions options) : EField_x(EField_x), kBT(kBT) {
    sqrtDielectric = options.sqrtDielectric;
    mu2 = options.mu2;
    charge = options.charge;
    alpha = options.alpha;
    attempt = options.attempt;
    dielectricConstant = options.dieletricConstant;
  }

  double marcusRate(double jeff2, double lambda_ij, double deltaE, double dx,
                    PType type) const;
  double singletDecay(double singletEnergy, MType type) const;
  double millerAbrahams(Eigen::Vector3d dr, double deltaE, PType type) const;
  double singletGeneration(Eigen::Vector3d dr, double deltaE,  double singletBinding, PType type) const;
  double singletDissociationCT(Eigen::Vector3d dr, double deltaE,  double singletBinding, PType type) const;
  double ctDissociation(Eigen::Vector3d dr, double deltaE, PType type) const;

private:
  double kBT;
  double EField_x;
  double sqrtDielectric; 
  double dielectricConstant;
  std::array<double, 2> mu2;
  std::array<double, 4> charge;
  std::array<double, 4> alpha;
  std::array<double, 4> attempt;
};
