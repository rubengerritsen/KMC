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
	initializeNeighbours();
	initializeParticles();

	std::cout << "Initialization and setup done." << std::endl;

	for (int step = 0; step < nrOfSteps; ++step) {
		computeNextEventRates();
		executeNextEvent();

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

void KmcRun::initializeNeighbours() {
	Eigen::Vector3d dr;
	double dist = 0;
	for (unsigned int i = 0; i < siteList.size(); ++i) {
		for (unsigned int j = i + 1; j < siteList.size(); ++j) {
			dr = pbc.dr_PBC_corrected(siteList[i].getCoordinates(), siteList[j].getCoordinates());
			dist = dr.norm();
			if (dist <= lR_cutOff) {
				siteList[i].addLRNeighbour(j);
				siteList[j].addLRNeighbour(i);
				if (dist <= sR_cutOff) {
					siteList[i].addSRNeighbour(j);
					siteList[j].addSRNeighbour(i);
				}
			}
		}
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

void KmcRun::computeNextEventRates() {

	next_event_list.clearNextEventList();

	for (unsigned int i = 0; i < particleList.size(); ++i) {
		Particle& part = particleList[i];
		switch (part.getType()) {
		case PType::elec:
			for (const auto& nb : siteList[part.getLocation()].getSRNeighbours()) {
				if (siteList[nb].isOccupied(PType::elec) || siteList[nb].isOccupied(PType::CT) || siteList[nb].isOccupied(PType::sing) || siteList[nb].isOccupied(PType::trip)) {
					; // nothing happens
				}
				else if (siteList[nb].isOccupied(PType::hole)) { //exciton generation
					next_event_list.pushNextEvent(rate_engine.millerAbrahamsGEN(siteList[part.getLocation()], siteList[nb], part.getType()), Transition::excitonGeneration, i, nb);
				}
				else { // normal hop
					next_event_list.pushNextEvent(rate_engine.millerAbrahams(siteList[part.getLocation()], siteList[nb], part.getType()), Transition::normalhop, i, nb);
				}
			}
			break;
		case PType::hole:
			for (const auto& nb : siteList[part.getLocation()].getSRNeighbours()) {
				if (siteList[nb].isOccupied(PType::hole) || siteList[nb].isOccupied(PType::CT) || siteList[nb].isOccupied(PType::sing) || siteList[nb].isOccupied(PType::trip)) {
					; // nothing happens
				}
				else if (siteList[nb].isOccupied(PType::elec)) { // exciton generation
					next_event_list.pushNextEvent(rate_engine.millerAbrahamsGEN(siteList[part.getLocation()], siteList[nb], part.getType()), Transition::excitonGeneration, i, nb);
				}
				else { // normal hop
					next_event_list.pushNextEvent(rate_engine.millerAbrahams(siteList[part.getLocation()], siteList[nb], part.getType()), Transition::normalhop, i, nb);
				}
			}
			break;
		case PType::sing:
			// it can hop, ...
			for (const auto& nb : siteList[part.getLocation()].getLRNeighbours()) { //Note: long range neighbourlist here
				if (!siteList[nb].isOccupied(PType::elec) && !siteList[nb].isOccupied(PType::hole) &&
					!siteList[nb].isOccupied(PType::CT) && !siteList[nb].isOccupied(PType::sing) && !siteList[nb].isOccupied(PType::trip)) {

					next_event_list.pushNextEvent(rate_engine.forster(siteList[part.getLocation()], siteList[nb]), Transition::normalhop ,i, nb); // normal "forster" hop
				}
			}
			// ... it can decay ...
			next_event_list.pushNextEvent(rate_engine.decay(part.getType()), Transition::decay, i, i);
			// ... or it will dissociate into a CT state.
			for (const auto& nb : siteList[part.getLocation()].getSRNeighbours()) { //Note: short range neighbourlist here
				if (!siteList[nb].isOccupied(PType::elec) && !siteList[nb].isOccupied(PType::hole) &&
					!siteList[nb].isOccupied(PType::CT) && !siteList[nb].isOccupied(PType::sing) && !siteList[nb].isOccupied(PType::trip)) {

					next_event_list.pushNextEvent(rate_engine.millerAbrahamsDIS(siteList[part.getLocation()], siteList[nb], PType::elec), Transition::dissociateElec, i, nb);
					next_event_list.pushNextEvent(rate_engine.millerAbrahamsDIS(siteList[part.getLocation()], siteList[nb], PType::hole), Transition::dissociateHole, i, nb);
				}
			}
			break;
		case PType::trip:
			// it can hop, ...
			for (const auto& nb : siteList[part.getLocation()].getSRNeighbours()) { //Note: short range neighbourlist here
				if (!siteList[nb].isOccupied(PType::elec) && !siteList[nb].isOccupied(PType::hole) &&
					!siteList[nb].isOccupied(PType::CT) && !siteList[nb].isOccupied(PType::sing) && !siteList[nb].isOccupied(PType::trip)) {

					next_event_list.pushNextEvent(rate_engine.millerAbrahams(siteList[part.getLocation()], siteList[nb]), Transition::normalhop, i, nb); // normal hop
				}
			}
			// ... it can decay ...
			next_event_list.pushNextEvent(rate_engine.decay(part.getType()), Transition::decay, i, i);
			// ... or it will dissociate into a CT state.
			for (const auto& nb : siteList[part.getLocation()].getSRNeighbours()) { //Note: short range neighbourlist here
				if (!siteList[nb].isOccupied(PType::elec) && !siteList[nb].isOccupied(PType::hole) &&
					!siteList[nb].isOccupied(PType::CT) && !siteList[nb].isOccupied(PType::sing) && !siteList[nb].isOccupied(PType::trip)) {

					next_event_list.pushNextEvent(rate_engine.millerAbrahamsDIS(siteList[part.getLocation()], siteList[nb], PType::elec), Transition::dissociateElec, i, nb);
					next_event_list.pushNextEvent(rate_engine.millerAbrahamsDIS(siteList[part.getLocation()], siteList[nb], PType::hole), Transition::dissociateHole, i, nb);
				}
			}
			break;
		case PType::CT:
			// it can recombine into an exciton (either the hole follows the electron or vice versa) or ...
			next_event_list.pushNextEvent(rate_engine.millerAbrahamsGEN(siteList[part.getLocationCTelec()], siteList[part.getLocation()], PType::elec ), Transition::excitonGeneration,i, part.getLocation());
			next_event_list.pushNextEvent(rate_engine.millerAbrahamsGEN(siteList[part.getLocation()], siteList[part.getLocationCTelec()], PType::elec), Transition::excitonGeneration, i, part.getLocationCTelec());
			// ... it can separate into free charges
			for (const auto& nb : siteList[part.getLocation()].getSRNeighbours()) {
				if (!siteList[nb].isOccupied(PType::elec) && !siteList[nb].isOccupied(PType::hole) &&
					!siteList[nb].isOccupied(PType::CT) && !siteList[nb].isOccupied(PType::sing) && !siteList[nb].isOccupied(PType::trip)) {

					next_event_list.pushNextEvent(rate_engine.millerAbrahamsCT_DIS(siteList[part.getLocation()], siteList[nb], PType::hole), Transition::dissociateHole, i, nb);
				}
			}
			for (const auto& nb : siteList[part.getLocationCTelec()].getSRNeighbours()) {
				if (!siteList[nb].isOccupied(PType::elec) && !siteList[nb].isOccupied(PType::hole) &&
					!siteList[nb].isOccupied(PType::CT) && !siteList[nb].isOccupied(PType::sing) && !siteList[nb].isOccupied(PType::trip)) {

					next_event_list.pushNextEvent(rate_engine.millerAbrahamsCT_DIS(siteList[part.getLocation()], siteList[nb], PType::elec), Transition::dissociateElec, i, nb);
				}
			}
			break;
		}
	}
}



/* UNDER CONSTRUCTION */
void KmcRun::executeNextEvent() {
	
	totalTime += random_engine.getInterArrivalTime(next_event_list.getTotalRate());

	std::tuple<Transition, int, int> nextEvent = next_event_list.getNextEvent(random_engine.getUniform01());
	
	int partID = std::get<1>(nextEvent);
	Particle& part = particleList[partID];
	auto it = particleList.begin() + partID; // get an iterator at the position of the particle in particleList; used below (move and pop_back)
	
	int newLocation = std::get<2>(nextEvent);
	int oldLocation = part.getLocation();

	switch (std::get<0>(nextEvent)) {
	case Transition::normalhop:
		part.jumpTo(newLocation, pbc.dr_PBC_corrected(siteList[oldLocation].getCoordinates(), siteList[newLocation].getCoordinates()));
		siteList[oldLocation].freeSite(part.getType(), totalTime);
		siteList[newLocation].setOccupied(part.getType(), totalTime);
		break;
	case Transition::decay:
		siteList[oldLocation].freeSite(part.getType(), totalTime);
		*it = std::move(particleList.back()); // A trick to do an efficient delete in a vector (move and pop_back)
		particleList.pop_back();
		break;
	case Transition::excitonGeneration:
		if (part.getType() == PType::elec) {

		}
		else { // it's a hole

		}
		break;

	}


			
}

