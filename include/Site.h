/***************************************************
 *
 * KMC MODEL FOR OPTOELECTRIC PROCESSES
 *
 * Author: Ruben Gerritsen
 *
 * Created on 14-01-2020
 *
 * Class to store all data concerning a site.
 **************************************************/

#pragma once
#include "Particle.h"
#include <vector>
#include <array>
#include <Eigen/Dense>
#include <iostream>
#include "PBC.h"
#include "PType.h"

class Site {
public:
	Site(Eigen::Vector3d coord, std::vector<double> energies);
	double getEnergy(PType pType) const { return energies[pType]; }
	Eigen::Vector3d getCoordinates() const { return coord; }

	void addSRNeighbour(int nb) { sRNeighbours.push_back(nb); }
	void addLRNeighbour(int nb) { lRNeighbours.push_back(nb); }
	bool isOccupied(PType type) const { return occupied[type]; }
	int isOccupiedBy(PType type) const;
	void setOccupied(PType type, int partID, double totalTime) { occupied[type] = true; startOccupation[type] = totalTime; occupiedBy[type] = partID; }
	void freeSite(PType type, double totalTime) { if(occupied[type]){ occupied[type] = false;} else{std::cout << "Attempt to free a non occupied site." << std::endl;} totalOccupation[type] += (totalTime - startOccupation[type]); }
	double getOccupation(PType type, double totalTime) { return occupied[type] ? totalOccupation[type] += (totalTime - startOccupation[type]) : totalOccupation[type]; }
	const std::vector<int>& getSRNeighbours() const { return sRNeighbours; }
	const std::vector<int>& getLRNeighbours() const { return lRNeighbours; }

private:
	std::vector<double> energies;
	Eigen::Vector3d coord;
	std::vector<int> sRNeighbours; // sR = short Range
	std::vector<int> lRNeighbours; // lR = long Range (for Forster transport)
	std::array<bool, 5> occupied {false};
	std::array<int, 5> occupiedBy{ 0 };
	static PBC pbc;
	std::array<double, 5> startOccupation{ 0 };
	std::array<double, 5> totalOccupation{ 0 };
};

std::ostream& operator<<(std::ostream& os, const Site& st);



