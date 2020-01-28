/***************************************************
 *
 * KMC MODEL FOR OPTOELECTRIC PROCESSES
 *
 * Author: Ruben Gerritsen
 *
 * Created on 14-01-2020
 *
 * Class to store all data concerning next events
 * Also computes the next event.
 **************************************************/

#include <vector>
#include <tuple>
#include <cmath>
#include "EnumNames.h"

class NextEventList {
public:
    void initializeListSize(int size){
        maxSize = size;
        rateList.resize(size);
        partList.resize(size);
        newLocation.resize(size);
        eventType.resize(size);
    }

    void pushNextEvent(double rate, Transition eventtype, int part, int loc);
    void resetNextEventList() { cPos = 0; totalRate = 0.0; }
    double getTotalRate() const { return totalRate; }
    int getNrOfEvents() const { return cPos; }
    
    std::tuple<Transition, int, int> getNextEvent(double random01) const;

    int size() { return rateList.size(); }

private:
    int maxSize = 10;
    int cPos = 0;
    std::vector<double> rateList;
    std::vector<int> partList;
    std::vector<int> newLocation;
    std::vector<Transition> eventType;
    double totalRate = 0.0;
    void resizeVectors();
};