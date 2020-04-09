#pragma once
struct SimulationOptions {
    int simID = 0;
    double maxTime = 1000000.0;
    double printAt = 0.1;
    int nrOfElectrons = 10;
    int nrOfHoles = 10;
    int nrOfSinglets = 10;
    double kBT = 0.026;
    double EField_x = 0.0;
    int SEED = 12345;
};