/***************************************************
 *
 * KMC MODEL FOR OPTOELECTRIC PROCESSES
 *
 * Author: Ruben Gerritsen
 *
 * Created on 14-01-2020
 *
 * Class to store all data concerning a site.
 * 
 **************************************************/

#pragma once
#include "Particle.h"
#include <array>
#include <iostream>
#include "EnumNames.h"

class Site {
public:
	Site(){}
	bool isOccupied(PType type) const { return occupied[type]; }
	int isOccupiedBy(PType type) const;

	/* Set a site occupied with a certain particle, if the site was already occupied use changeOccupied() instead */
	void setOccupied(PType type, int partID, double totalTime) { occupied[type] = true; startOccupation[type] = totalTime; occupiedBy[type] = partID; }
	/* Change the current occupation to another type of occupation */
	void changeOccupied(PType oldType, PType newType, int partID, double totalTime) { occupied[oldType] = false; occupied[newType] = true; startOccupation[newType] = totalTime; occupiedBy[newType] = partID; }
	void freeSite(PType type, double totalTime) { if(occupied[type]){ occupied[type] = false;} else{std::cout << "Attempt to free a non occupied site." << std::endl;} totalOccupation[type] += (totalTime - startOccupation[type]); }
	double getOccupation(PType type, double totalTime) { return occupied[type] ? totalOccupation[type] += (totalTime - startOccupation[type]) : totalOccupation[type]; }

private:
	std::array<bool, 5> occupied {false};
	std::array<int, 5> occupiedBy{ 0 };
	std::array<double, 5> startOccupation{ 0 };
	std::array<double, 5> totalOccupation{ 0 };
};
