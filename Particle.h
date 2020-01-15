#pragma once
#include <Eigen/Dense>
#include <iostream>

enum PType { elec = 0, hole = 1, trip = 2, sing = 3 };

class Particle {
public:
	Particle(int loc, PType type) : location{ loc }, type{ type }, energyLevel{ 0 } {};
	PType getType() const { return type; }
	double distanceTravelled() const { return dr_travelled.norm(); }
	int getLocation() const { return location; }

private:
	int location;
	PType type;
	int energyLevel;
	Eigen::Vector3d dr_travelled;

};
