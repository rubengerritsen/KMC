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


std::tuple<int, int> NextEventList::getNextEvent(double random01) const {
	double cumSum = 0;
	double select = totalRate * random01;

	for (unsigned int i = 0; i < rateList; ++i) {
		cumSum += rateList[i];
		if (cumSum >= select) {
			return std::tuple<Transition, int, int>{eventType[i], partList[i], newLocation[i]}
		}
	}
}