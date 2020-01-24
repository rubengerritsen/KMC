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
#include <tuple>

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
		std::cout << "Number of sites in the simulation: " << siteList.size() << "\n";
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
	for (int partID = 0; partID < 100; ++partID) {
		location = random_engine.getRandomSite();
		while (siteList[location].isOccupied(PType::elec)) { //Get a unique location
			location = random_engine.getRandomSite();
		}
		tempParticle = new Particle(location, PType::elec);
		particleList.push_back(*tempParticle);
		siteList[location].setOccupied(PType::elec, partID, 0.0);
	}
}

void KmcRun::computeNextEventRates() {

	next_event_list.resetNextEventList();

	for (unsigned int i = 0; i < particleList.size(); ++i) {
		Particle& part = particleList[i];
		switch (part.getType()) {
		case PType::elec:
			for (const auto& nb : siteList[part.getLocation()].getSRNeighbours()) {
				if (siteList[nb].isOccupied(PType::elec) || siteList[nb].isOccupied(PType::CT) || siteList[nb].isOccupied(PType::sing) || siteList[nb].isOccupied(PType::trip)) {
					; // nothing happens
				}
				else if (siteList[nb].isOccupied(PType::hole)) { //exciton generation
					next_event_list.pushNextEvent(rate_engine.millerAbrahamsGEN(siteList[part.getLocation()], siteList[nb], part.getType()), Transition::excitonFromElec, i, nb);
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
					next_event_list.pushNextEvent(rate_engine.millerAbrahamsGEN(siteList[part.getLocation()], siteList[nb], part.getType()), Transition::excitonFromHole, i, nb);
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

					next_event_list.pushNextEvent(rate_engine.millerAbrahamsDIS(siteList[part.getLocation()], siteList[nb], PType::elec), Transition::singToCTViaElec, i, nb);
					next_event_list.pushNextEvent(rate_engine.millerAbrahamsDIS(siteList[part.getLocation()], siteList[nb], PType::hole), Transition::singToCTViaHole, i, nb);
				}
			}
			break;
		case PType::trip:
			// it can hop, ...
			for (const auto& nb : siteList[part.getLocation()].getSRNeighbours()) { //Note: short range neighbourlist here
				if (!siteList[nb].isOccupied(PType::elec) && !siteList[nb].isOccupied(PType::hole) &&
					!siteList[nb].isOccupied(PType::CT) && !siteList[nb].isOccupied(PType::sing) && !siteList[nb].isOccupied(PType::trip)) {

					next_event_list.pushNextEvent(rate_engine.millerAbrahams(siteList[part.getLocation()], siteList[nb], part.getType()), Transition::normalhop, i, nb); // normal hop
				}
			}
			// ... it can decay ...
			next_event_list.pushNextEvent(rate_engine.decay(part.getType()), Transition::decay, i, i);
			// ... or it will dissociate into a CT state.
			for (const auto& nb : siteList[part.getLocation()].getSRNeighbours()) { //Note: short range neighbourlist here
				if (!siteList[nb].isOccupied(PType::elec) && !siteList[nb].isOccupied(PType::hole) &&
					!siteList[nb].isOccupied(PType::CT) && !siteList[nb].isOccupied(PType::sing) && !siteList[nb].isOccupied(PType::trip)) {

					next_event_list.pushNextEvent(rate_engine.millerAbrahamsDIS(siteList[part.getLocation()], siteList[nb], PType::elec), Transition::tripToCTViaElec, i, nb);
					next_event_list.pushNextEvent(rate_engine.millerAbrahamsDIS(siteList[part.getLocation()], siteList[nb], PType::hole), Transition::tripToCTViaHole, i, nb);
				}
			}
			break;
		case PType::CT:
			// it can recombine into an exciton (either the hole follows the electron or vice versa) or ...
			next_event_list.pushNextEvent(rate_engine.millerAbrahamsGEN(siteList[part.getLocationCTelec()], siteList[part.getLocation()], PType::elec ), Transition::excitonFromElecCT,i, part.getLocation());
			next_event_list.pushNextEvent(rate_engine.millerAbrahamsGEN(siteList[part.getLocation()], siteList[part.getLocationCTelec()], PType::hole), Transition::excitonFromHoleCT, i, part.getLocationCTelec());
			// ... it can separate into free charges
			for (const auto& nb : siteList[part.getLocation()].getSRNeighbours()) {
				if (!siteList[nb].isOccupied(PType::elec) && !siteList[nb].isOccupied(PType::hole) &&
					!siteList[nb].isOccupied(PType::CT) && !siteList[nb].isOccupied(PType::sing) && !siteList[nb].isOccupied(PType::trip)) {

					next_event_list.pushNextEvent(rate_engine.millerAbrahamsCT_DIS(siteList[part.getLocation()], siteList[nb], PType::hole), Transition::CTdisViaHole, i, nb);
				}
			}
			for (const auto& nb : siteList[part.getLocationCTelec()].getSRNeighbours()) {
				if (!siteList[nb].isOccupied(PType::elec) && !siteList[nb].isOccupied(PType::hole) &&
					!siteList[nb].isOccupied(PType::CT) && !siteList[nb].isOccupied(PType::sing) && !siteList[nb].isOccupied(PType::trip)) {

					next_event_list.pushNextEvent(rate_engine.millerAbrahamsCT_DIS(siteList[part.getLocation()], siteList[nb], PType::elec), Transition::CTdisViaElec, i, nb);
				}
			}
			break;
		}
	}
}

void KmcRun::executeNextEvent() {
	
	totalTime += random_engine.getInterArrivalTime(next_event_list.getTotalRate());

	std::tuple<Transition, int, int> nextEvent = next_event_list.getNextEvent(random_engine.getUniform01());
	
	int partID = std::get<1>(nextEvent);
	Particle& part = particleList[partID];
	auto it = particleList.begin() + partID; // get an iterator at the position of the particle in particleList; used below (move and pop_back)
	
	int newLocation = std::get<2>(nextEvent);
	int oldLocation = part.getLocation();

	PType type;
	Particle* tempParticle;

	switch (std::get<0>(nextEvent)) {
	case Transition::normalhop:
		part.jumpTo(newLocation, pbc.dr_PBC_corrected(siteList[oldLocation].getCoordinates(), siteList[newLocation].getCoordinates()));
		siteList[oldLocation].freeSite(part.getType(), totalTime);
		siteList[newLocation].setOccupied(part.getType(), partID, totalTime);
		break;

	case Transition::decay:
		siteList[oldLocation].freeSite(part.getType(), totalTime);
		*it = std::move(particleList.back()); // A trick to do an efficient delete in a vector (move and pop_back)
		particleList.pop_back();
		break;

	case Transition::excitonFromElec:
		siteList[oldLocation].freeSite(PType::elec, totalTime);
		type = particleList[siteList[newLocation].isOccupiedBy(PType::hole)].makeExciton(random_engine.getUniform01());
		siteList[newLocation].setOccupied(type, partID, totalTime);
		*it = std::move(particleList.back()); // A trick to do an efficient delete in a vector (move and pop_back)
		particleList.pop_back();
		break;

	case Transition::excitonFromElecCT:
		siteList[part.getLocationCTelec()].freeSite(PType::CT, totalTime); // free the site of the electron
		type = part.makeExciton(random_engine.getUniform01()); // create an exciton in the place of the hole
		siteList[part.getLocation()].setOccupied(type, partID, totalTime);
		break;

	case Transition::excitonFromHole:
		siteList[oldLocation].freeSite(PType::hole, totalTime);
		type = particleList[siteList[newLocation].isOccupiedBy(PType::elec)].makeExciton(random_engine.getUniform01());
		siteList[newLocation].setOccupied(type, partID, totalTime);
		*it = std::move(particleList.back()); // A trick to do an efficient delete in a vector (move and pop_back)
		particleList.pop_back();
		break;

	case Transition::excitonFromHoleCT:
		siteList[part.getLocation()].freeSite(PType::CT, totalTime); // free the site of the hole
		type = part.makeExciton(random_engine.getUniform01()); // create an exciton in the place of the elec
		siteList[part.getLocationCTelec()].setOccupied(type, partID, totalTime);
		part.setLocation(part.getLocationCTelec());
		break;

	case Transition::singToCTViaElec:
		siteList[newLocation].setOccupied(PType::CT, partID, totalTime); //note that new here represents the new position of the electron
		siteList[oldLocation].setOccupied(PType::CT, partID, totalTime); //note that old here represents the original position of the sing 
		part.makeCTState(oldLocation, newLocation);
		break;

	case Transition::singToCTViaHole:
		siteList[newLocation].setOccupied(PType::CT, partID, totalTime); //note that new here represents the new position of the hole
		siteList[oldLocation].setOccupied(PType::CT, partID, totalTime); //note that old here represents the original position of the sing 
		part.makeCTState(newLocation, oldLocation);
		break;

	case Transition::tripToCTViaElec:
		siteList[newLocation].setOccupied(PType::CT, partID, totalTime); //note that new here represents the new position of the electron
		siteList[oldLocation].setOccupied(PType::CT, partID, totalTime); //note that old here represents the original position of the sing 
		part.makeCTState(oldLocation, newLocation);
		break;

	case Transition::tripToCTViaHole:
		siteList[newLocation].setOccupied(PType::CT, partID, totalTime); //note that new here represents the new position of the hole
		siteList[oldLocation].setOccupied(PType::CT, partID, totalTime); //note that old here represents the original position of the sing 
		part.makeCTState(newLocation, oldLocation);
		break;

	case Transition::CTdisViaElec:
		/* free old sites */
		siteList[part.getLocationCTelec()].freeSite(PType::CT, totalTime); 
		siteList[part.getLocation()].freeSite(PType::CT, totalTime);

		/* create the elec and hole */
		part.makeElectron(newLocation);
		tempParticle = new Particle(oldLocation, PType::hole);
		particleList.push_back(*tempParticle);

		/* Set sites occupied */
		siteList[newLocation].setOccupied(PType::elec, partID, totalTime); 
		siteList[oldLocation].setOccupied(PType::hole, particleList.size() - 1, totalTime);
		break;

	case Transition::CTdisViaHole:
		/* free old sites */
		siteList[part.getLocationCTelec()].freeSite(PType::CT, totalTime);
		siteList[part.getLocation()].freeSite(PType::CT, totalTime);

		/* create the elec and hole */
		part.makeHole(newLocation);
		tempParticle = new Particle(oldLocation, PType::elec);
		particleList.push_back(*tempParticle);

		/* Set sites occupied */
		siteList[newLocation].setOccupied(PType::hole, partID, totalTime);
		siteList[oldLocation].setOccupied(PType::elec, particleList.size() - 1, totalTime);
		break;
	}		
}

