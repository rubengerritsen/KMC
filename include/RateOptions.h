#pragma once
#include <array>
struct RateOptions {
    double sqrtDielectric = 1.51;
    std::array<double, 2> mu2{0.055508,0.039497 };
    std::array<double, 4> charge{-1.0, 1.0, 0.0, 0.0};
    std::array<double, 4> alpha{4.69, 5.86, 1.0, 1.0};
    std::array<double, 4> attempt{2.49e4, 5.63e4, 1.0, 1.0};
    double dieletricConstant = 4.0;
};