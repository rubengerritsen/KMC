#include "Topology.h"
#include "EnumNames.h"
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <fstream>
#include <iostream>

void Topology::readTopologyFromFile(std::string filename) {
  std::ifstream file(filename);
  if (file.is_open()) {
    double x, y, z;
    double eng_e, eng_h, eng_s;
    std::string molName;
    Eigen::Vector3d tempCoord;
    std::vector<double> tempENG(3);
    while (file >> x >> y >> z) {
      tempCoord << x, y, z;
      siteLocations.push_back(tempCoord);
    }
  } else {
    std::cout << "Unable to open topology file: " << filename << std::endl;
    std::cout << "Terminating execution." << std::endl;
    exit(EXIT_FAILURE);
  }
}

void Topology::readReorganisationEnergies(std::string filename) {
  boost::property_tree::ptree pt, lambdas;
  read_xml(filename, pt);
  lambdas = pt.get_child("lambdas");

  /* Obviously this is no long term solution */
  std::string moleculeName = "benzene";
  MType molType = MType::ben;
  lambdaXtoN[molType][PType::elec] = lambdas.get<double>(moleculeName + ".e2n");
  lambdaXtoN[molType][PType::hole] = lambdas.get<double>(moleculeName + ".h2n");
  lambdaXtoN[molType][PType::sing] = lambdas.get<double>(moleculeName + ".s2n");
  lambdaNtoX[molType][PType::elec] = lambdas.get<double>(moleculeName + ".n2e");
  lambdaNtoX[molType][PType::hole] = lambdas.get<double>(moleculeName + ".n2h");
  lambdaNtoX[molType][PType::sing] = lambdas.get<double>(moleculeName + ".n2s");

  moleculeName = "tcne";
  molType = MType::tcne;
  lambdaXtoN[molType][PType::elec] = lambdas.get<double>(moleculeName + ".e2n");
  lambdaXtoN[molType][PType::hole] = lambdas.get<double>(moleculeName + ".h2n");
  lambdaXtoN[molType][PType::sing] = lambdas.get<double>(moleculeName + ".s2n");
  lambdaNtoX[molType][PType::elec] = lambdas.get<double>(moleculeName + ".n2e");
  lambdaNtoX[molType][PType::hole] = lambdas.get<double>(moleculeName + ".n2h");
  lambdaNtoX[molType][PType::sing] = lambdas.get<double>(moleculeName + ".n2s");
}

void Topology::printHead(int nrOfItems) const {
  for (int i = 0; i < nrOfItems; i++) {
    std::cout << siteLocations[i][0] << " " << siteLocations[i][1] << " "
              << siteLocations[i][2] << " " << siteEnergies[i][0] << " "
              << siteEnergies[i][1] << " " << siteEnergies[i][2] << " "
              << moleculeType[i] << std::endl;
  }
}

void Topology::printReorganisationEnergies() {
  std::cout << "some reorganisation energies lambda xton for benzene:\n";
  for (int i = 0; i < 3; i++) {
    std::cout << lambdaXtoN[MType::ben][i] << " ";
  }
  std::cout << std::endl;
}