#include "Topology.h"
#include "Neighbourlist.h"
#include <iostream>
#include <sstream>
#include "Logger.h"

class TestCase {
    public:
    TestCase(const Topology& topol, const Neighbourlist& nbList, int simID): topol(topol), nbList(nbList), simID(simID) {    }

    void printTopTest() const {
        Logger{} << "printing from: " << simID << "\n";
        return topol.printHead(1);}

    private:
    const Topology& topol;
    const Neighbourlist& nbList;
    int simID;
};