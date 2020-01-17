/***************************************************
 *
 * KMC MODEL FOR OPTOELECTRIC PROCESSES
 *
 * Author: Ruben Gerritsen
 *
 * Created on 14-01-2020
 *
 **************************************************/

#include "Site.h"
#include <vector>
#include <iostream>
#include "Particle.h"
#include "PType.h"

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

/* Takes in a uniform random number [0,1) and returns the site to jump to. */
int Site::getNextHop(double uniform, PType type) {
	double select = uniform * totalRates[type];
	double rateSum = 0;
	if (type == PType::sing) {
		for (unsigned int i=0; i < lRRates.size(); ++i) {
			rateSum += lRRates[i];
			if (rateSum >= select) {
				return lRNeighbours[i];
			}
		}
	}
	else {
		for (unsigned int i=0; i < sRRates[type].size(); ++i) {
			rateSum += sRRates[type][i];
			if (rateSum >= select) {
				return sRNeighbours[i];
			}
		}
	}
	std::cout << "getNextHop was unsuccessful" << std::endl; 
	exit(EXIT_FAILURE);
}