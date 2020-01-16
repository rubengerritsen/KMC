#pragma once
#include <vector>
#include <array>
#include <Eigen/Dense>
#include "Particle.h"
#include "PBC.h"


enum PType;

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
	void setOccupied(PType type, double time) { occupied[type] = true; startOccupation[type] = time; }
	void freeSite(PType type, double time) { occupied[type] = false; totalOccupation[type] += (time - startOccupation[type]); }
	void computeTotals();
	double getTotalOutRate(PType type) { return totalRates[type]; }
	int getNextHop(double uniform, PType type);
	double getOccupation(PType type, double time) { return occupied[type] ? totalOccupation[type] += (time - startOccupation[type]) : totalOccupation[type]; }

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



