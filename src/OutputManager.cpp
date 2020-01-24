#include "OutputManager.h"
#include <ctime>
#include "PType.h"


void OutputManager::printSiteOccupations(std::vector<Site>& siteList, double totalTime) {

	struct tm * ltm;
	time_t now = time(0);
	ltm = localtime( &now);
	ltm->tm_mon = ltm->tm_mon + 1;
	std::string filename = outputPath + str( boost::format("siteOcc_%02d%2d%2d.txt") % ltm->tm_mon % ltm->tm_mday % ltm->tm_hour ) ;

	std::ofstream outFile;
	outFile.open(filename);
	if (outFile.is_open()) {
		for (auto& site : siteList) {
			outFile << site.getEnergy(PType::elec) << " " << site.getOccupation(PType::elec, totalTime) / totalTime << " "
				<< site.getEnergy(PType::hole) << " " << site.getOccupation(PType::hole, totalTime) / totalTime << " "
				<< site.getEnergy(PType::trip) << " " << site.getOccupation(PType::trip, totalTime) / totalTime << " "
				<< site.getEnergy(PType::sing) << " " << site.getOccupation(PType::sing, totalTime) / totalTime << "\n"; 
		}
		std::cout << "Site occupations were printed to:\n\t"<<  filename << "\n";
	}
	else {
		std::cout << "Could not open output file: " << filename << "\n";
	}
	outFile.close();
}