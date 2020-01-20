#include <ctime>
#include <boost/format.hpp>
#include <string>
#include <vector>
#include "Site.h"
#include <fstream>

class OutputManager {
public:
	/* Outputs a file with the site occupations and energies (ln: energy occ).*/
	void printSiteOccupations(std::vector<Site>& siteList, double totalTime);


private:
	std::string outputPath = "./output/";

};