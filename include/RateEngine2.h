#pragma once
#include <array>
#include "EnumNames.h"

class RateEngine2 {
public:
  RateEngine2(double EField_x, double kBT) : EField_x(EField_x), kBT(kBT) {}

  double marcusRate(double jeff2, double lambda_ij, double deltaE, double dx, PType type);

private:
  double kBT;
  double EField_x;
  std::array<double,4> charge{-1.0, 1.0, 0.0, 0.0};
};
