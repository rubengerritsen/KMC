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
#include "PType.h"

class NextEventList {
public:
    void pushNextEvent(double rate, Transition eventtype, int part, int loc) { 
        rateList.push_back(rate); eventType.push_back(eventtype); partList.push_back(part); newLocation.push_back(loc); totalRate += rate;
    }
    
    void clearNextEventList() { rateList.clear(); eventType.clear(); partList.clear(); newLocation.clear(); totalRate = 0.0; }

    double getTotalRate() const { return totalRate; }
    
    std::tuple<Transition, int, int> getNextEvent(double random01) const;

private:
    std::vector<double> rateList;
    std::vector<int> partList;
    std::vector<int> newLocation;
    std::vector<Transition> eventType;
    double totalRate = 0.0;
};