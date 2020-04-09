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

std::tuple<Transition, int, Neighbour> NextEventList::getNextEvent(double random01) const {
	double cumSum = 0;
	double select = totalRate * random01;

	for (int i = 0; i < cPos ; ++i) { // Note cPos here is smaller then rateList.size()
		cumSum += rateList[i];
		if (cumSum >= select) {
            return std::tuple<Transition, int, Neighbour> {eventType[i], partList[i], newLocation[i]};
		}
	}
    std::cout << "Next event could not be found\n";
    exit(EXIT_FAILURE);
}

void NextEventList::pushNextEvent(double rate, Transition eventtype, int part, Neighbour loc) {
    rateList[cPos] = rate; 
    eventType[cPos] =eventtype; 
    partList[cPos] = part; 
    newLocation[cPos] = loc; 
    totalRate += rate;

    ++cPos;
    if (cPos >= maxSize) {
        resizeVectors();
    }
}

void NextEventList::resizeVectors() {
    maxSize = (int) std::floor(maxSize * 1.1);
    rateList.resize(maxSize);
    partList.resize(maxSize);
    newLocation.resize(maxSize);
    eventType.resize(maxSize);
}