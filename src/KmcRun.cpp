/***************************************************
 *
 * KMC MODEL FOR OPTOELECTRIC PROCESSES
 *
 * Author: Ruben Gerritsen
 *
 * Created on 14-01-2020
 *
 **************************************************/

#include "KmcRun.h"
#include <iostream>
#include <chrono>
#include "OutputManager.h"

void KmcRun::runSimulation() {
	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

	initializeSites();
	initializeNeighboursAndRates();
	initializeParticles();

	std::cout << "Initialization and setup done." << std::endl;

	for (int step = 0; step < nrOfSteps; ++step) {
		findAndExecuteNextEvent();

		/* Give some feedback on the progress */
		if ((step + 1) % (nrOfSteps / 100) == 0) {
			std::cout << "\rProgress: " << 100.0 * (step + 1) / (nrOfSteps) << "%" << std::flush;
		}
	}
	std::cout << std::endl;

	OutputManager out;
	out.printSiteOccupations(siteList, totalTime);

	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
	std::cout << "Total simulation time: " << (std::chrono::duration_cast<std::chrono::seconds>(end - begin).count()) << "s" << std::endl;

}

void KmcRun::initializeSites() {
	std::ifstream myfile(siteFile);
	Site* tempSite;
	Eigen::Vector3d tempCoord;
	std::vector<double> tempEnergies(4);
	if (myfile.is_open()) {
		double x;
		double y;
		double z;
		while (myfile >> x >> y >> z) {
			tempCoord << x, y, z;
			tempEnergies[int(PType::elec)] = random_engine.getDOSEnergy(PType::elec);
			tempEnergies[int(PType::hole)] = random_engine.getDOSEnergy(PType::hole);
			tempEnergies[int(PType::sing)] = random_engine.getDOSEnergy(PType::sing);
			tempEnergies[int(PType::trip)] = random_engine.getDOSEnergy(PType::trip);
			tempSite = new Site(tempCoord, tempEnergies);
			siteList.push_back(*tempSite);
		}
		myfile.close();
		random_engine.setNrOfSites(siteList.size());
	}
	else {
		std::cout << "Unable to open file: " << siteFile << std::endl;
		std::cout << "Terminating execution." << std::endl;
		exit(EXIT_FAILURE);
	}
}

void KmcRun::initializeNeighboursAndRates() {
	Eigen::Vector3d dr;
	double dist = 0;

	for (unsigned int i = 0; i < siteList.size(); ++i) {
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

void KmcRun::initializeParticles() {
	Particle* tempParticle;
	int location = 0;
	for (int i = 0; i < 100; ++i) {
		location = random_engine.getRandomSite();
		while (siteList[location].isOccupied(PType::elec)) { //Get a unique location
			location = random_engine.getRandomSite();
		}
		tempParticle = new Particle(location, PType::elec);
		particleList.push_back(*tempParticle);
		siteList[location].setOccupied(PType::elec, 0.0);
	}
}

void KmcRun::findAndExecuteNextEvent() {
	double totalHopRate = 0;
	double totalDecayRate = 0;
	for (const auto& part : particleList) {
		/* First we consider all hopping rates */
		totalHopRate += siteList[part.getLocation()].getTotalOutRate(PType::elec);

		/* Next the decay rates (possibilities of sing/trip to decay). */
		if (part.getType() == PType::sing || part.getType() == PType::trip) {
			totalDecayRate += rate_engine.decay(part.getType());
		}

		/* Finally all transition rates from exciton to CT free charges etc. */
		if (part.getType() == PType::sing) {

		}
		
	}

	double totalRate = totalHopRate + totalDecayRate;

	totalTime += random_engine.getInterArrivalTime(totalRate);

	double select = totalRate * random_engine.getUniform01();
	totalRate = 0;
	int newLocation;
	int oldLocation;
	for (auto& part : particleList) {
		totalRate += siteList[part.getLocation()].getTotalOutRate(PType::elec);
		if (totalRate >= select) {
			newLocation = siteList[part.getLocation()].getNextHop(random_engine.getUniform01(), part.getType());
			if (siteList[newLocation].isOccupied(part.getType())) {
				break; // no jump occurs since the site is occupied. 
			}
			oldLocation = part.getLocation();
			part.jumpTo(newLocation, pbc.dr_PBC_corrected(siteList[oldLocation].getCoordinates(), siteList[newLocation].getCoordinates()));
			siteList[oldLocation].freeSite(part.getType(), totalTime);
			siteList[newLocation].setOccupied(part.getType(), totalTime);
			break;
		}
	}
}

