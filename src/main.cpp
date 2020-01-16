/***************************************************
 * 
 * KMC MODEL FOR OPTOELECTRIC PROCESSES
 *
 * Author: Ruben Gerritsen
 *
 * Created on 14-01-2020
 *
 **************************************************/
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <random>
#include <cmath>
#include "Particle.h"
#include "Site.h"
#include "PBC.h"


std::vector<Site> siteList;
std::vector<Particle> particleList;
std::mt19937_64 rng(12345);

PBC pbc(67.821, 67.821, 67.821);

std::normal_distribution<double> elec_DOS(5.0, 2.0);
std::normal_distribution<double> hole_DOS(5.0, 2.0);
std::normal_distribution<double> sing_DOS(5.0, 2.0);
std::normal_distribution<double> trip_DOS(5.0, 2.0); 
std::uniform_int_distribution<int> siteDist(0, 511);
std::uniform_real_distribution<double> uniform(0.0, 1.0);

void initializeSites() {
	std::string filename{ "sites.txt" };
	std::ifstream myfile(filename);
	Site *tempSite;
	Eigen::Vector3d tempCoord;
	std::vector<double> tempEnergies(4);
	if (myfile.is_open()) {
		double x;
		double y;
		double z;
		while (myfile >> x >> y >> z) {
			tempCoord << x, y, z;
			tempEnergies[int(PType::elec)] = elec_DOS(rng);
			tempEnergies[int(PType::hole)] = hole_DOS(rng);
			tempEnergies[int(PType::sing)] = sing_DOS(rng);
			tempEnergies[int(PType::trip)] = trip_DOS(rng);
			tempSite = new Site(tempCoord, tempEnergies);
			siteList.push_back(*tempSite);
		}
		myfile.close();
	}
	else {
		std::cout << "Unable to open file: " << filename << std::endl;
		std::cout << "Terminating execution." << std::endl;
		exit(EXIT_FAILURE);
	}
}

Eigen::Vector3d dr_PBC_corrected(const Eigen::Vector3d& v, const Eigen::Vector3d& w, const Eigen::Vector3d& boxDimension) {
	return v.array() - w.array() - floor((v.array() - w.array()) / boxDimension.array() + 0.5) * boxDimension.array();
}

void initializeNeighboursAndRates() {
	Eigen::Vector3d dr(0,0,0), simBoxDimensions(67.821, 67.821, 67.821);
	double dist;
	double lR_cutOff = 25.0, sR_cutOff = 15.0;

	for (unsigned int i = 0; i < siteList.size(); ++i){
		for (unsigned int j = i + 1; j < siteList.size(); ++j) {
			dr = dr_PBC_corrected(siteList[i].getCoordinates(), siteList[j].getCoordinates(), simBoxDimensions);
			dist = dr.norm();
			if (dist <= lR_cutOff) {
				siteList[i].addLRNeighbour(j);
				siteList[i].addLRRate(dexterRate(dist, dr[0], siteList[i].getEnergy(PType::elec), siteList[j].getEnergy(PType::elec), 1, 0.15, 0.0));
				siteList[j].addLRNeighbour(i);
				siteList[i].addLRRate(dexterRate(dist, dr[0], siteList[j].getEnergy(PType::elec), siteList[i].getEnergy(PType::elec), 1, 0.15, 0.0));
				if (dist <= sR_cutOff) {
					siteList[i].addSRNeighbour(j);
					siteList[i].addSRRate(PType::elec, millerAbrahamsRate(dist, dr[0], siteList[i].getEnergy(PType::elec), siteList[j].getEnergy(PType::elec), 1, 0.15, 0.0)); 
					siteList[i].addSRRate(PType::hole, millerAbrahamsRate(dist, dr[0], siteList[i].getEnergy(PType::hole), siteList[j].getEnergy(PType::hole), 1, 0.15, 0.0));
					siteList[i].addSRRate(PType::trip, millerAbrahamsRate(dist, dr[0], siteList[i].getEnergy(PType::trip), siteList[j].getEnergy(PType::trip), 1, 0.15, 0.0));
					siteList[j].addSRNeighbour(i);
					siteList[i].addSRRate(PType::elec, millerAbrahamsRate(dist, dr[0], siteList[i].getEnergy(PType::elec), siteList[j].getEnergy(PType::elec), 1, 0.15, 0.0));
					siteList[i].addSRRate(PType::hole, millerAbrahamsRate(dist, dr[0], siteList[i].getEnergy(PType::hole), siteList[j].getEnergy(PType::hole), 1, 0.15, 0.0));
					siteList[i].addSRRate(PType::trip, millerAbrahamsRate(dist, dr[0], siteList[i].getEnergy(PType::trip), siteList[j].getEnergy(PType::trip), 1, 0.15, 0.0));
				}
			}
		}
	}
	for (auto& site : siteList) {
		site.computeTotals();
	}
}


void initializeParticles() {
	Particle* tempParticle;
	int location = 0;
	for (int i = 0; i < 10; ++i) {
		location = siteDist(rng);
		while (siteList[location].isOccupied(PType::elec)) { //Get a unique location
			location = siteDist(rng);
		}
		tempParticle = new Particle(location, PType::elec);
		siteList[location].setOccupied(PType::elec);
	}
}

void findNextEvent() {
	double totalRate = 0;
	for (auto& part : particleList) {
		totalRate += siteList[part.getLocation()].getTotalOutRate(PType::elec);
	}
	double select = totalRate * uniform(rng);
	totalRate = 0;
	for (auto& part : particleList) {
		totalRate += siteList[part.getLocation()].getTotalOutRate(PType::elec);
		if (totalRate >= select) {

		}
	}

}




int main() {
	//initializeSites();
	//initializeNeighboursAndRates();
	//initializeParticles();

	Eigen::Vector3d v(1, 2, 1.3), w(1.3, 2.1, -0.3);
	PBC pbc{2, 2, 2};

	std::cout << pbc.dr_3vector(v, w);




    return 0;
}


