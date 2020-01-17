#include "SingleRun.h"

void SingleRun::runSimulation() {
	initializeSites();
	initializeNeighboursAndRates();
	initializeParticles();
	for (unsigned int i = 0; i < nrOfSteps; ++i) {
		findAndExecuteNextEvent();
	}

	printSiteOccupation();
}

void SingleRun::initializeSites() {
	std::ifstream myfile("./input/sites.txt");
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

void SingleRun::initializeNeighboursAndRates() {
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

void SingleRun::initializeParticles() {
	Particle* tempParticle;
	int location = 0;
	for (int i = 0; i < 1; ++i) {
		location = random_engine.getRandomSite();
		while (siteList[location].isOccupied(PType::elec)) { //Get a unique location
			location = random_engine.getRandomSite();
		}
		tempParticle = new Particle(location, PType::elec);
		particleList.push_back(*tempParticle);
		siteList[location].setOccupied(PType::elec, 0.0);
	}
}

void SingleRun::findAndExecuteNextEvent() {
	double totalRate = 0;
	int newLocation;
	int oldLocation;
	for (auto& part : particleList) {
		totalRate += siteList[part.getLocation()].getTotalOutRate(PType::elec);
	}

	totalTime += random_engine.getInterArrivalTime(totalRate);

	double select = totalRate * random_engine.getUniform01();
	totalRate = 0;
	for (auto& part : particleList) {
		totalRate += siteList[part.getLocation()].getTotalOutRate(PType::elec);
		if (totalRate >= select) {
			newLocation = siteList[part.getLocation()].getNextHop(random_engine.getUniform01(), part.getType());
			oldLocation = part.getLocation();
			part.jumpTo(newLocation, pbc.dr_PBC_corrected(siteList[oldLocation].getCoordinates(), siteList[newLocation].getCoordinates()));
			siteList[oldLocation].freeSite(part.getType(), totalTime);
			siteList[newLocation].setOccupied(part.getType(), totalTime);
			break;
		}
	}
}

void SingleRun::printSiteOccupation() {
	std::ofstream outFile;
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