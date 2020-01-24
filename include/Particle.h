/***************************************************
 *
 * KMC MODEL FOR OPTOELECTRIC PROCESSES
 *
 * Author: Ruben Gerritsen
 *
 * Created on 14-01-2020
 *
 * Class to store all data concerning a particle.
 **************************************************/

#pragma once
#include <Eigen/Dense>
#include <iostream>
#include "Site.h"
#include "PType.h"


class Particle {
public:
	Particle(int loc, PType type) : location{ loc }, type{ type }, energyLevel{ 0 } {};
	PType getType() const { return type; }
	/* Creates an exciton with 3:1 ratio for trip vs sing states */
	PType makeExciton(double random01) { return random01 > 0.75 ? type = PType::sing : type = PType::trip; }
	void makeCTState(int holeLocation, int elecLocation) { type = PType::CT; locationCTelec = elecLocation; location = holeLocation; }
	void makeElectron(int loc) { type = PType::elec; location = loc; }
	void makeHole(int loc) { type = PType::hole; location = loc; }
	double distanceTravelled() const { return dr_travelled.norm(); }
	int getLocation() const { return location; }
	int getLocationCTelec() const { return locationCTelec; }
	void jumpTo(int loc, Eigen::Vector3d dr) { location = loc; dr_travelled += dr; }
	void setLocation(int loc) { location = loc; }
	void setLocationCTelec(int loc) { locationCTelec = loc; }

private:
	int location;
	int locationCTelec;
	PType type;
	int energyLevel;
	Eigen::Vector3d dr_travelled;
};
