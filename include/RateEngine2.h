#pragma once
#include "EnumNames.h"
#include <array>
#include "constants.h"

class RateEngine2 {
public:
  RateEngine2(double EField_x, double kBT) : EField_x(EField_x), kBT(kBT) {
    mu2[MType::ben] = std::pow(0.235601519, 2); 
    mu2[MType::tcne] = std::pow(0.1987388, 2);// bohr2nm * 3.7556 = 0.1987...
  }

  double marcusRate(double jeff2, double lambda_ij, double deltaE, double dx,
                    PType type);
  double singletDecay(double singletEnergy, MType type);

private:
  double kBT;
  double EField_x;
  double sqrtDielectric = std::sqrt(2.28); // value of benzene solution
  std::array<double, 2> mu2;
  std::array<double, 4> charge{-1.0, 1.0, 0.0, 0.0};
};
