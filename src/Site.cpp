#pragma once
#include "Site.h"
#include <vector>
#include <iostream>
#include "Particle.h"

Site::Site(Eigen::Vector3d coord, std::vector<double> energies) : coord(coord), energies(energies)
{
	if (energies.size() != 4) {
	std::cout << "Site initialized with wrong number of energies\n";
	}
}

std::ostream& operator<<(std::ostream& os, const Site& st) {
	os << st.getCoordinates();
	return os;
}

void Site::computeTotals() {
	double totalRate = 0.0;
	for (auto& nbRate : sRRates[PType::elec]) {
		totalRate += nbRate;
	}
	totalRates[PType::elec] = totalRate;
	totalRate = 0.0;
	for (auto& nbRate : sRRates[PType::hole]) {
		totalRate += nbRate;
	}
	totalRates[PType::hole] = totalRate;
	totalRate = 0.0;
	for (auto& nbRate : sRRates[PType::trip]) {
		totalRate += nbRate;
	}
	totalRates[PType::trip] = totalRate;
	totalRate = 0.0;
	for (auto& nbRate : lRRates) {
		totalRate += nbRate;
	}
	totalRates[PType::sing] = totalRate;
}