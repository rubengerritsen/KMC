#include "OutputManager.h"
#include <ctime>
#include "EnumNames.h"


void OutputManager::printSiteOccupations(std::vector<Site>& siteList, double totalTime) {

	struct tm * ltm;
	time_t now = time(0);
	ltm = localtime( &now);
	ltm->tm_mon = ltm->tm_mon + 1;
	std::string filename = "fiets";

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
		std::cout << "Could not open output file: " << filename << std::endl;
	}
	outFile.close();
}

void OutputManager::printParticleInfo(std::vector<Particle>& particleList){
	std::cout << "Alive particles: " << std::endl;
	std::array<int,5> nrPerType {0};
	for (auto& part : particleList){
		if(part.isAlive()) nrPerType[part.getType()] += 1;
	}
	for (auto& elem : nrPerType){
		std::cout << elem << "  " ;
	}
	std::cout << std::endl;

	std::cout << "Dead particles: " << std::endl;
	nrPerType = {0, 0, 0, 0, 0};
	for (auto& part : particleList){
		if(!part.isAlive()) nrPerType[part.getType()] += 1;
	}
	for (auto& elem : nrPerType){
		std::cout << elem << "  " ;
	}
	std::cout << std::endl;

}