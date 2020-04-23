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
#include "EnumNames.h"


class Particle {
public:
	Particle(int loc, PType type) : location{ loc }, type{ type }, energyLevel{ 0 } {};
	
	PType getType() const { return type; }
	/* Creates an exciton with 3:1 ratio for trip vs sing states */
	PType makeExciton(double random01) { return type = PType::sing;}
	void makeCTState(int holeLocation, int elecLocation) { type = PType::CT; locationCTelec = elecLocation; location = holeLocation; }
	void makeElectron(int loc) { type = PType::elec; location = loc; }
	void makeHole(int loc) { type = PType::hole; location = loc; }

	void jumpTo(int loc, Eigen::Vector3d dr) { location = loc; dr_travelled += dr; }
	void setLocation(int loc) { location = loc; }
	void setLocationCTelec(int loc) { locationCTelec = loc; }

	void killParticle(double time) { alive = false; timeOfDeath = time; }
	bool isAlive() const { return alive; }

	double distanceTravelled() const { return dr_travelled.norm(); }
	int getLocation() const { return location; }
	int getLocationCTelec() const { return locationCTelec; }


private:
	int location;
	int locationCTelec = 0;
	bool alive = true;
	double timeOfDeath = 0.0;

	PType type;
	int energyLevel;
	Eigen::Vector3d dr_travelled;
};
