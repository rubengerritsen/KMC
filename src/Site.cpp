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

int Site::isOccupiedBy(PType type) const { 
	if (occupied[type]){
		return occupiedBy[type];
		}else { 
			std::cout << "This site is no longer occupied!" << std::endl;
			return occupiedBy[type];
		}
}