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

#include <vector>
#include <tuple>
#include "PType.h"

class NextEventList {
public:

    NextEventList(int maxSize) : maxSize(maxSize) {
        rateList.resize(maxSize);
        partList.resize(maxSize);
        newLocation.resize(maxSize);
        eventType.resize(maxSize);
    }

    void pushNextEvent(double rate, Transition eventtype, int part, int loc);
    void resetNextEventList() { cPos = 0; totalRate = 0.0; }
    double getTotalRate() const { return totalRate; }
    int getNrOfEvents() const { return cPos; }
    
    std::tuple<Transition, int, int> getNextEvent(double random01) const;

    int size() { return rateList.size(); }

private:
    int maxSize;
    int cPos = 0;
    std::vector<double> rateList;
    std::vector<int> partList;
    std::vector<int> newLocation;
    std::vector<Transition> eventType;
    double totalRate = 0.0;
    void resizeVectors();
};