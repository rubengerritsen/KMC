#pragma once
#include <Eigen/Dense>
struct Neighbour {
    int nb =0;
    double rate_e =0;
    double rate_h =0;
    double rate_s =0;
    Eigen::Vector3d dr;
};
