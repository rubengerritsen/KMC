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
	double xDistFrom(Site otherSite) const { return coord[0] - otherSite.getCoordinates()[0]; }
	const Eigen::Vector3d& getCoordinates() const { return coord; }
	//Eigen::Vector3d dr_from(Site& otherSite) const { return pbc.dr_3vector(otherSite.getCoordinates(), coord); }

	void addSRNeighbour(int nb) {sRNeighbours.push_back(nb); }
	void addLRNeighbour(int nb) {lRNeighbours.push_back(nb); }
	void addSRRate(PType type, double rate) { sRRates[type].push_back(rate); }
	void addLRRate(double rate) { lRRates.push_back(rate); }
	bool isOccupied(PType type) const { return occupied[type]; }
	void setOccupied(PType type) { occupied[type] = true; }
	void freeSite(PType type) { occupied[type] = false; }
	void computeTotals();
	std::vector<int> getSRNeighbours() const { return sRNeighbours; }
	std::vector<double> getSRRates(PType type) const { return sRRates[type]; }
	double getTotalOutRate(PType type) { return totalRates[type]; }

private:
	std::vector<double> energies;
	Eigen::Vector3d coord;
	std::vector<int> sRNeighbours; // sR = short Range
	std::vector<int> lRNeighbours; // lR = long Range (for Förster transport)
	std::array<std::vector<double>,3> sRRates;
	std::vector<double> lRRates;
	std::array<bool, 4> occupied {false};
	std::array<double,4> totalRates { 0 };
	static PBC pbc;
};

std::ostream& operator<<(std::ostream& os, const Site& st);



