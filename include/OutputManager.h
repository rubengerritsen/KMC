#include <ctime>
#include <boost/format.hpp>
#include <string>
#include <vector>
#include "Site.h"
#include "Particle.h"
#include <fstream>

class OutputManager {
public:
	/* Outputs a file with the site occupations and energies (ln: energy occ).*/
	void printSiteOccupations(std::vector<Site>& siteList, double totalTime);
	void printParticleInfo(std::vector<Particle>&);


private:
	std::string outputPath = "C:/Users/s134864/source/repos/KMC/KMC/output/";

};