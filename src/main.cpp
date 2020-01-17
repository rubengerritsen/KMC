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
#include <fstream>
#include "Particle.h"
#include "Site.h"
#include "PBC.h"
#include "RateEngine.h"
#include "PType.h"


std::vector<Site> siteList;
std::vector<Particle> particleList;
std::mt19937_64 rng(145);

std::array<double, 4> v0 {1, 1, 1, 1};
std::array<double, 4> alpha { 0.000008, 0.1, 0.1, 0.1 };
std::array<double, 4> charge{ -1, -1, -1, -1 };
double E_Field{ 0.0 };
double kBT{ 0.026 };

PBC pbc(67.821, 67.821, 67.821);
RateEngine rate_engine(v0, alpha, charge, E_Field, kBT, pbc);

std::normal_distribution<double> elec_DOS(0, 0.026);
std::normal_distribution<double> hole_DOS(0, 0.026);
std::normal_distribution<double> sing_DOS(0, 0.026);
std::normal_distribution<double> trip_DOS(0, 0.026);
std::uniform_int_distribution<int> siteDist(0, 511);
std::uniform_real_distribution<double> uniform(0.0, 1.0);

double totalTime = 0;

void initializeSites() {
	std::string filename{ "./input/sites.txt" };
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


void initializeNeighboursAndRates() {
	double lR_cutOff = 25.0, sR_cutOff = 15.0;
	Eigen::Vector3d dr;
	double dist = 0;

	for (unsigned int i = 0; i < siteList.size(); ++i){
		for (unsigned int j = i + 1; j < siteList.size(); ++j) {
			dr = pbc.dr_PBC_corrected(siteList[i].getCoordinates(), siteList[j].getCoordinates());
			dist = dr.norm();
			if (dist <= lR_cutOff) {
				siteList[i].addLRNeighbour(j);
				siteList[i].addLRRate(rate_engine.dexter(siteList[i], siteList[j]));
				siteList[j].addLRNeighbour(i);
				siteList[j].addLRRate(rate_engine.dexter(siteList[j], siteList[i]));
				if (dist <= sR_cutOff) {
					siteList[i].addSRNeighbour(j);
					siteList[i].addSRRate(PType::elec, rate_engine.millerAbrahams(siteList[i], siteList[j], PType::elec)); 
					siteList[i].addSRRate(PType::hole, rate_engine.millerAbrahams(siteList[i], siteList[j], PType::hole));
					siteList[i].addSRRate(PType::trip, rate_engine.millerAbrahams(siteList[i], siteList[j], PType::trip));
					siteList[j].addSRNeighbour(i);
					siteList[j].addSRRate(PType::elec, rate_engine.millerAbrahams(siteList[j], siteList[i], PType::elec));
					siteList[j].addSRRate(PType::hole, rate_engine.millerAbrahams(siteList[j], siteList[i], PType::hole));
					siteList[j].addSRRate(PType::trip, rate_engine.millerAbrahams(siteList[j], siteList[i], PType::trip));
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
	for (int i = 0; i < 1; ++i) {
		location = siteDist(rng);
		while (siteList[location].isOccupied(PType::elec)) { //Get a unique location
			location = siteDist(rng);
		}
		tempParticle = new Particle(location, PType::elec);
		particleList.push_back(*tempParticle);
		siteList[location].setOccupied(PType::elec, 0.0);
	}
}

void findAndExecuteNextEvent() {
	double totalRate = 0;
	int newLocation;
	int oldLocation;
	for (auto& part : particleList) {
		totalRate += siteList[part.getLocation()].getTotalOutRate(PType::elec);
	}

	totalTime += -(1.0 / totalRate) * log(uniform(rng));

	double select = totalRate * uniform(rng);
	totalRate = 0;
	for (auto& part : particleList) {
		totalRate += siteList[part.getLocation()].getTotalOutRate(PType::elec);
		if (totalRate >= select) {
			newLocation = siteList[part.getLocation()].getNextHop(uniform(rng), part.getType());
			oldLocation = part.getLocation();
			part.jumpTo(newLocation, pbc.dr_PBC_corrected(siteList[oldLocation].getCoordinates(), siteList[newLocation].getCoordinates()));
			siteList[oldLocation].freeSite(part.getType(), totalTime);
			siteList[newLocation].setOccupied(part.getType(), totalTime);
			break;
		}
	}
}

void printSiteOccupation() {
	std::ofstream outFile;
	std::string outputFile = "./output/occ.txt";
	outFile.open(outputFile);
	if (outFile.is_open()) {
		for (auto& site : siteList) {
			outFile << site.getEnergy(PType::elec) << " " << site.getOccupation(PType::elec, totalTime) << std::endl;
		}
	}
	else {
		std::cout << "Could not open output file: " << outputFile << std::endl;
	}
	outFile.close();
}




int main() {
	int maxSteps = 100000;
	initializeSites();
	initializeNeighboursAndRates();
	initializeParticles();
	for (int step = 0; step < maxSteps; ++step) {
		findAndExecuteNextEvent();
	}
	printSiteOccupation();

    return 0;
}


