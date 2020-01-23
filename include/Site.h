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
#include "PBC.h"
#include "PType.h"

class Site {
public:
	Site(Eigen::Vector3d coord, std::vector<double> energies);
	double getEnergy(PType pType) const { return energies[pType]; }
	Eigen::Vector3d getCoordinates() const { return coord; }

	void addSRNeighbour(int nb) {sRNeighbours.push_back(nb); }
	void addLRNeighbour(int nb) {lRNeighbours.push_back(nb); }
	void addSRRate(PType type, double rate) { sRRates[type].push_back(rate); }
	void addLRRate(double rate) { lRRates.push_back(rate); }
	bool isOccupied(PType type) const { return occupied[type]; }
	void setOccupied(PType type, double totalTime) { occupied[type] = true; startOccupation[type] = totalTime; }
	void freeSite(PType type, double totalTime) { occupied[type] = false; totalOccupation[type] += (totalTime - startOccupation[type]); }
	void computeTotals();
	double getTotalOutRate(PType type) { return totalRates[type]; }
	int getNextHop(double uniform, PType type) const;
	double getOccupation(PType type, double totalTime) { return occupied[type] ? totalOccupation[type] += (totalTime - startOccupation[type]) : totalOccupation[type]; }
	const std::vector<int>& getSRNeighbours() const { return sRNeighbours; }
	const std::vector<int>& getLRNeighbours() const { return lRNeighbours; }

private:
	std::vector<double> energies;
	Eigen::Vector3d coord;
	std::vector<int> sRNeighbours; // sR = short Range
	std::vector<int> lRNeighbours; // lR = long Range (for F�rster transport)
	std::array<std::vector<double>,3> sRRates;
	std::vector<double> lRRates;
	std::array<bool, 4> occupied {false};
	std::array<double,4> totalRates { 0 };
	static PBC pbc;
	std::array<double, 4> startOccupation{ 0 };
	std::array<double, 4> totalOccupation{ 0 };
};

std::ostream& operator<<(std::ostream& os, const Site& st);



