/***************************************************
 *
 * KMC MODEL FOR OPTOELECTRIC PROCESSES
 *
 * Author: Ruben Gerritsen
 *
 * Created on 14-01-2020
 *
 * Class to store all data concerning next events
 **************************************************/

#include "NextEventList.h"
#include <iostream>
#include <numeric>

std::tuple<Transition, int, Neighbour>
NextEventList::getNextEvent(double random01) const {
  double cumSum = 0;
  double select = totalRate * random01;

  if (totalRate <= 0.0 || cPos == 0) {
    std::cout << "\nNo next event, program done.\n";
    // We need to return something but in reality we only care about "done"
    return std::tuple<Transition, int, Neighbour>{Transition::done, 0,
                                                  newLocation[0]};
  }

  for (int i = 0; i < cPos;
       ++i) { // Note cPos here is smaller then rateList.size()
    cumSum += rateList[i];
    if (cumSum >= select) {
      return std::tuple<Transition, int, Neighbour>{eventType[i], partList[i],
                                                    newLocation[i]};
    }
  }
}

void NextEventList::pushNextEvent(double rate, Transition eventtype, int part,
                                  Neighbour loc) {
  rateList[cPos] = rate;
  eventType[cPos] = eventtype;
  partList[cPos] = part;
  newLocation[cPos] = loc;
  totalRate += rate;

  ++cPos;
  if (cPos >= maxSize) {
    resizeVectors();
  }
}

void NextEventList::resizeVectors() {
  maxSize = (int)std::floor(maxSize * 1.1);
  rateList.resize(maxSize);
  partList.resize(maxSize);
  newLocation.resize(maxSize);
  eventType.resize(maxSize);
}