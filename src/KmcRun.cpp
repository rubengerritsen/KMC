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
#include "Neighbour.h"
#include <boost/format.hpp>
#include <chrono>
#include <fstream>
#include <iostream>
#include <tuple>

void KmcRun::runSimulation() {
  std::chrono::steady_clock::time_point begin =
      std::chrono::steady_clock::now();

  initializeSites();
  initializeParticles();

  std::cout << "Initialization and setup done." << std::endl;

  int step = 0;
  bool lastEvent;
  while (totalTime < simOptions.maxTime && step < simOptions.maxStep) {
    std::cout << "\r" << double(step) / simOptions.maxStep << "%" << std::flush;
    computeNextEventList();
    lastEvent = executeNextEvent();
    if (lastEvent) {
      break;
    }
    step = step + 1;
  }
  std::cout << std::endl;
  out.endTOF();

  // OutputManager out;
  // out.printSiteOccupations(siteList, totalTime);
  // out.printParticleInfo(particleList);

  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  std::cout
      << "Total simulation time: "
      << (std::chrono::duration_cast<std::chrono::seconds>(end - begin).count())
      << "s" << std::endl;
  std::cout
      << "\n****************************************************************\n"
      << std::endl;
}

void KmcRun::initializeSites() {
  siteList.resize(topol.getNrOfSites());
  random_engine.setNrOfSites(topol.getNrOfSites());
}

void KmcRun::initializeParticles() {
  Particle *tempParticle;
  int location = 0;
  int partID = 0;
  /* electrons */
  for (int i = 0; i < simOptions.nrOfElectrons; ++i) {
    location = i;
    tempParticle = new Particle(location, PType::elec);
    particleList.push_back(*tempParticle);
    siteList[location].setOccupied(PType::elec, partID, 0.0);
    partID++;
  }
  /* holes */
  for (int i = 0; i < simOptions.nrOfHoles; ++i) {
    location = i;
    tempParticle = new Particle(location, PType::hole);
    particleList.push_back(*tempParticle);
    siteList[location].setOccupied(PType::hole, partID, 0.0);
    partID++;
  }
}

void KmcRun::computeNextEventList() {
  next_event_list.reset();
  Neighbour empty_nb;

  for (unsigned int i = 0; i < particleList.size(); ++i) {
    Particle &part = particleList[i];
    if (part.isAlive()) {
      switch (part.getType()) {
      case PType::elec:
        for (const auto &nb : nbList.getNeighbours(part.getLocation())) {
          if (!siteList[nb.nb].isOccupied(PType::elec) &&
              !siteList[nb.nb].isOccupied(PType::hole)) { // normal hop
            next_event_list.pushNextEvent(nb.rate_e, Transition::normalhop, i,
                                          nb);
          }
        }
        break;
      case PType::hole:
        for (const auto &nb : nbList.getNeighbours(part.getLocation())) {
          if (!siteList[nb.nb].isOccupied(PType::hole) &&
              !siteList[nb.nb].isOccupied(PType::hole)) { // normal hop
            next_event_list.pushNextEvent(nb.rate_h, Transition::normalhop, i,
                                          nb);
          }
        }
        break;

      default:
        std::cout << "Particle type not known, something is seriously wrong in "
                     "the setup files or code!\nterminating program."
                  << std::endl;
        exit(EXIT_FAILURE);
      }
    }
  }
}

bool KmcRun::executeNextEvent() {
  totalTime +=
      random_engine.getInterArrivalTime(next_event_list.getTotalRate());

  std::tuple<Transition, int, Neighbour> nextEvent =
      next_event_list.getNextEvent(random_engine.getUniform01());

  if (std::get<0>(nextEvent) == Transition::done) {
    bool noNextEvents = true;
    return noNextEvents;
  }

  int partID = std::get<1>(nextEvent);
  Particle &part = particleList[partID];

  Neighbour newLocation = std::get<2>(nextEvent);
  int oldLocation = part.getLocation();

  PType type;
  Particle *tempParticle;

  switch (std::get<0>(nextEvent)) {
  case Transition::normalhop:
    if (newLocation.nb == simOptions.sink) {
      siteList[oldLocation].freeSite(part.getType(), totalTime);
      part.killParticle(totalTime);
      out.registerTOF(totalTime);
    } else {
      part.jumpTo(newLocation.nb, newLocation.dr);
      siteList[oldLocation].freeSite(part.getType(), totalTime);
      siteList[newLocation.nb].setOccupied(part.getType(), partID, totalTime);
    }
    break;
  default:
    std::cout
        << "Impossible event for doing time of flight, terminating program"
        << std::endl;
    exit(EXIT_FAILURE);
  }
  return false;
}
