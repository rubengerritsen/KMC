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

  int progress = 0;

  // out.registerParticlePositions(particleList, totalTime);
  // out.registerState(particleList, totalTime);
  out.registerNumbers(particleList, topol, totalTime);

  while (totalTime < simOptions.maxTime) {
    computeNextEventRates();
    executeNextEvent();
    // out.registerParticlePositions(particleList, totalTime);
    // out.registerState(particleList, totalTime);
    out.registerNumbers(particleList, topol, totalTime);
  }
  std::cout << std::endl;

  // OutputManager out;
  // out.printSiteOccupations(siteList, totalTime);
  // out.printParticleInfo(particleList);

  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  std::cout
      << "Total simulation time: "
      << (std::chrono::duration_cast<std::chrono::seconds>(end - begin).count())
      << "s" << std::endl;
}

void KmcRun::initializeSites() {
  siteList.resize(topol.getNrOfSites());
  random_engine.setNrOfSites(topol.getNrOfSites());
}

void KmcRun::initializeParticles() {
  Particle *tempParticle;
  if (simOptions.allSinglets) {
    for (int i = 0; i < topol.getNrOfSites(); ++i) {
      tempParticle = new Particle(i, PType::sing);
      particleList.push_back(*tempParticle);
      siteList[i].setOccupied(PType::sing, i, 0.0);
    }
  } else {
    int location = 0;
    int partID = 0;
    /* electrons */
    for (int i = 0; i < simOptions.nrOfElectrons; ++i) {
      location = random_engine.getRandomSite();
      while (siteList[location].isOccupied(PType::elec)) { // Get a unique
                                                           // location
        location = random_engine.getRandomSite();
      }
      tempParticle = new Particle(location, PType::elec);
      particleList.push_back(*tempParticle);
      siteList[location].setOccupied(PType::elec, partID, 0.0);
      partID++;
    }
    /* holes */
    for (int i = 0; i < simOptions.nrOfHoles; ++i) {
      location = random_engine.getRandomSite();
      while (
          siteList[location].isOccupied(PType::elec) ||
          siteList[location].isOccupied(PType::hole)) { // Get a unique location
        location = random_engine.getRandomSite();
      }
      tempParticle = new Particle(location, PType::hole);
      particleList.push_back(*tempParticle);
      siteList[location].setOccupied(PType::hole, partID, 0.0);
      partID++;
    }
    /* singlets */
    for (int i = 0; i < simOptions.nrOfSinglets; ++i) {
      location = random_engine.getRandomSite();
      while (
          siteList[location].isOccupied(PType::elec) ||
          siteList[location].isOccupied(PType::hole) ||
          siteList[location].isOccupied(PType::sing)) { // Get a unique location
        location = random_engine.getRandomSite();
      }
      tempParticle = new Particle(location, PType::sing);
      particleList.push_back(*tempParticle);
      siteList[location].setOccupied(PType::sing, partID, 0.0);
      partID++;
    }
  }
}

void KmcRun::computeNextEventRates() {
  next_event_list.reset();
  Neighbour empty_nb;

  for (unsigned int i = 0; i < particleList.size(); ++i) {
    Particle &part = particleList[i];
    if (part.isAlive()) {
      switch (part.getType()) {
      case PType::elec:
        for (const auto &nb : nbList.getSRNeighbours(part.getLocation())) {
          if (siteList[nb.nb].isOccupied(PType::elec) ||
              siteList[nb.nb].isOccupied(PType::CT) ||
              siteList[nb.nb].isOccupied(PType::sing)) {
            ; // nothing happens
          } else if (siteList[nb.nb].isOccupied(
                         PType::hole)) { // exciton generation
            next_event_list.pushNextEvent(nb.rate_e_s,
                                          Transition::excitonFromElec, i, nb);
          } else { // normal hop
            next_event_list.pushNextEvent(nb.rate_e, Transition::normalhop, i,
                                          nb);
          }
        }
        break;
      case PType::hole:
        for (const auto &nb : nbList.getSRNeighbours(part.getLocation())) {
          if (siteList[nb.nb].isOccupied(PType::hole) ||
              siteList[nb.nb].isOccupied(PType::CT) ||
              siteList[nb.nb].isOccupied(PType::sing) ||
              siteList[nb.nb].isOccupied(PType::trip)) {
            ; // nothing happens
          } else if (siteList[nb.nb].isOccupied(
                         PType::elec)) { // exciton generation
            next_event_list.pushNextEvent(nb.rate_h_s,
                                          Transition::excitonFromHole, i, nb);
          } else { // normal hop
            next_event_list.pushNextEvent(nb.rate_h, Transition::normalhop, i,
                                          nb);
          }
        }
        break;
      case PType::sing:
        // it can hop, ...
        for (const auto &nb :
             nbList.getLRNeighbours(part.getLocation())) { // Note: long range
                                                           // neighbourlist here
          if (!siteList[nb.nb].isOccupied(PType::elec) &&
              !siteList[nb.nb].isOccupied(PType::hole) &&
              !siteList[nb.nb].isOccupied(PType::CT) &&
              !siteList[nb.nb].isOccupied(PType::sing) &&
              !siteList[nb.nb].isOccupied(PType::trip)) {

            next_event_list.pushNextEvent(nb.rate_s, Transition::normalhop, i,
                                          nb);
          }
        }
        // ... it can decay ...
        next_event_list.pushNextEvent(
            rate_engine.singletDecay(
                topol.getEnergy(part.getLocation(), PType::sing),
                topol.getMolType(part.getLocation())),
            Transition::decay, i, empty_nb);
        // ... or it will dissociate into a CT state.
        for (const auto &nb : nbList.getSRNeighbours(
                 part.getLocation())) { // Note: short range neighbourlist here
          if (!siteList[nb.nb].isOccupied(PType::elec) &&
              !siteList[nb.nb].isOccupied(PType::hole) &&
              !siteList[nb.nb].isOccupied(PType::CT) &&
              !siteList[nb.nb].isOccupied(PType::sing)) {
                // No CT-states with two same molecule types
                if( topol.getMolType(part.getLocation()) != topol.getMolType(nb.nb)){
              next_event_list.pushNextEvent(nb.rate_s_ct_e,
                                            Transition::singToCTViaElec, i, nb);
              next_event_list.pushNextEvent(nb.rate_s_ct_h,
                                            Transition::singToCTViaHole, i, nb);
                }
          }
        }
        break;
      case PType::CT:
        // It can separate into free charges
        for (const auto &nb : nbList.getSRNeighbours(part.getLocation())) {
          if (!siteList[nb.nb].isOccupied(PType::elec) &&
              !siteList[nb.nb].isOccupied(PType::hole) &&
              !siteList[nb.nb].isOccupied(PType::CT) &&
              !siteList[nb.nb].isOccupied(PType::sing)) {

            next_event_list.pushNextEvent(
                rate_engine.ctDissociation(
                    nb.dr,
                    topol.getDeltaEnergy(part.getLocation(), nb.nb,
                                         PType::hole),
                    PType::hole),
                Transition::CTdisViaHole, i, nb);
          }
        }
        for (const auto &nb :
             nbList.getSRNeighbours(part.getLocationCTelec())) {
          if (!siteList[nb.nb].isOccupied(PType::elec) &&
              !siteList[nb.nb].isOccupied(PType::hole) &&
              !siteList[nb.nb].isOccupied(PType::CT) &&
              !siteList[nb.nb].isOccupied(PType::sing)) {

            next_event_list.pushNextEvent(
                rate_engine.ctDissociation(
                    nb.dr,
                    topol.getDeltaEnergy(part.getLocationCTelec(), nb.nb,
                                         PType::elec),
                    PType::elec),
                Transition::CTdisViaElec, i, nb);
          }
        }
        break;
      }
    }
  }
}

void KmcRun::executeNextEvent() {
  totalTime +=
      random_engine.getInterArrivalTime(next_event_list.getTotalRate());

  std::tuple<Transition, int, Neighbour> nextEvent =
      next_event_list.getNextEvent(random_engine.getUniform01());

  int partID = std::get<1>(nextEvent);
  Particle &part = particleList[partID];

  Neighbour newLocation = std::get<2>(nextEvent);
  int oldLocation = part.getLocation();

  PType type;
  Particle *tempParticle;

  switch (std::get<0>(nextEvent)) {
  case Transition::normalhop:
    part.jumpTo(newLocation.nb, newLocation.dr);
    siteList[oldLocation].freeSite(part.getType(), totalTime);
    siteList[newLocation.nb].setOccupied(part.getType(), partID, totalTime);
    break;

  case Transition::decay:
    siteList[oldLocation].freeSite(part.getType(), totalTime);
    part.killParticle(totalTime);
    break;

  case Transition::excitonFromElec:
    siteList[oldLocation].freeSite(PType::elec, totalTime);
    type = particleList[siteList[newLocation.nb].isOccupiedBy(PType::hole)]
               .makeExciton(random_engine.getUniform01());
    siteList[newLocation.nb].changeOccupied(PType::hole, type, partID,
                                            totalTime);
    part.killParticle(totalTime);
    break;

  case Transition::excitonFromElecCT:
    siteList[part.getLocationCTelec()].freeSite(
        PType::CT, totalTime); // free the site of the electron
    type = part.makeExciton(
        random_engine
            .getUniform01()); // create an exciton in the place of the hole
    siteList[part.getLocation()].changeOccupied(PType::CT, type, partID,
                                                totalTime);
    break;

  case Transition::excitonFromHole:
    siteList[oldLocation].freeSite(PType::hole, totalTime);
    type = particleList[siteList[newLocation.nb].isOccupiedBy(PType::elec)]
               .makeExciton(random_engine.getUniform01());
    siteList[newLocation.nb].changeOccupied(PType::elec, type, partID,
                                            totalTime);
    part.killParticle(totalTime);
    break;

  case Transition::excitonFromHoleCT:
    siteList[part.getLocation()].freeSite(
        PType::CT, totalTime); // free the site of the hole
    type = part.makeExciton(
        random_engine
            .getUniform01()); // create an exciton in the place of the elec
    siteList[part.getLocationCTelec()].changeOccupied(PType::CT, type, partID,
                                                      totalTime);
    part.setLocation(part.getLocationCTelec());
    break;

  case Transition::singToCTViaElec:
    siteList[newLocation.nb].setOccupied(
        PType::CT, partID, totalTime); // note that new here represents the new
                                       // position of the electron
    siteList[oldLocation].changeOccupied(
        PType::sing, PType::CT, partID,
        totalTime); // note that old here represents the original position of
                    // the sing
    part.makeCTState(oldLocation, newLocation.nb);
    break;

  case Transition::singToCTViaHole:
    siteList[newLocation.nb].setOccupied(
        PType::CT, partID, totalTime); // note that new here represents the new
                                       // position of the hole
    siteList[oldLocation].changeOccupied(
        PType::sing, PType::CT, partID,
        totalTime); // note that old here represents the original position of
                    // the sing
    part.makeCTState(newLocation.nb, oldLocation);
    break;

  case Transition::tripToCTViaElec:
    siteList[newLocation.nb].setOccupied(
        PType::CT, partID, totalTime); // note that new here represents the new
                                       // position of the electron
    siteList[oldLocation].changeOccupied(
        PType::trip, PType::CT, partID,
        totalTime); // note that old here represents the original position of
                    // the sing
    part.makeCTState(oldLocation, newLocation.nb);
    break;

  case Transition::tripToCTViaHole:
    siteList[newLocation.nb].setOccupied(
        PType::CT, partID, totalTime); // note that new here represents the new
                                       // position of the hole
    siteList[oldLocation].changeOccupied(
        PType::trip, PType::CT, partID,
        totalTime); // note that old here represents the original position of
                    // the sing
    part.makeCTState(newLocation.nb, oldLocation);
    break;

  case Transition::CTdisViaElec:
    /* free old sites */
    siteList[part.getLocationCTelec()].freeSite(PType::CT, totalTime);
    siteList[part.getLocation()].freeSite(PType::CT, totalTime);

    /* create the elec and hole */
    part.makeElectron(newLocation.nb);
    tempParticle = new Particle(oldLocation, PType::hole);
    particleList.push_back(*tempParticle);

    /* Set sites occupied */
    siteList[newLocation.nb].setOccupied(PType::elec, partID, totalTime);
    siteList[oldLocation].setOccupied(PType::hole, particleList.size() - 1,
                                      totalTime);
    break;

  case Transition::CTdisViaHole:
    int oldElecLocation = part.getLocationCTelec();

    /* free old sites */
    siteList[part.getLocationCTelec()].freeSite(PType::CT, totalTime);
    siteList[part.getLocation()].freeSite(PType::CT, totalTime);

    /* create the elec and hole */
    part.makeHole(newLocation.nb);
    tempParticle = new Particle(oldElecLocation, PType::elec);
    particleList.push_back(*tempParticle);

    /* Set sites occupied */
    siteList[newLocation.nb].setOccupied(PType::hole, partID, totalTime);
    siteList[oldElecLocation].setOccupied(PType::elec, particleList.size() - 1,
                                          totalTime);
    break;
  }
}
