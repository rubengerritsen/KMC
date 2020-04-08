/***************************************************
 *
 * KMC MODEL FOR OPTOELECTRIC PROCESSES
 *
 * Author: Ruben Gerritsen
 *
 * Created on 14-01-2020
 *
 * Class that handles all output.
 * 
 **************************************************/
#pragma once
#include <ctime>
#include <string>
#include <vector>
#include "Site.h"
#include "Particle.h"
#include <fstream>

class OutputManager {
public:
	/* Outputs a file with the site occupations and energies (ln: energy occ).*/
	void printSiteOccupations(std::vector<Site>& siteList, double totalTime);

	/* Prints the current state of all particles to the console */
	void printParticleInfo(std::vector<Particle>&);


private:
	std::string outputPath = "./output/";

};