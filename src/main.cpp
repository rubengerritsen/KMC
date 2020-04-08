/***************************************************
 *
 * KMC MODEL FOR OPTOELECTRIC PROCESSES
 *
 * Author: Ruben Gerritsen
 *
 * Created on 14-01-2020
 *
 **************************************************/
#include "EnumNames.h"
#include "KmcRun.h"
#include "Neighbourlist.h"
#include "PBC.h"
#include "Particle.h"
#include "RandomEngine.h"
#include "RateEngine.h"
#include "Site.h"
#include "Topology.h"
#include "TestCase.h"
#include <array>
#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <fstream>
#include <iostream>
#include <string>
#include <omp.h>


void runSimulation(const Topology& topol, const Neighbourlist& nbList, int simID){
  TestCase test(topol, nbList, simID);
  test.printTopTest();
}

void setupAndExecuteSimulation(int ac, char *av[]) {

  /*********************************************/
  /* First we parse the command line arguments */
  /*********************************************/
  std::string optionFile;
  int nrOfThreads = 1;
  try {
    boost::program_options::options_description desc("Possible options");
    desc.add_options()("help,h", "produce help message")(
        "optionfile,o", boost::program_options::value<std::string>(),
        "path to option file (required)")("threads,t",
                                          boost::program_options::value<int>(),
                                          "number of threads to use");

    boost::program_options::variables_map vm;
    boost::program_options::store(
        boost::program_options::parse_command_line(ac, av, desc), vm);
    boost::program_options::notify(vm);

    if (vm.count("help") || ac == 1) {
      std::cout << "General usage: " << av[0]
                << " -o <path/to/optionfile.xml>\n";
      std::cout << desc << "\n";
    }

    if (vm.count("optionfile")) {
      optionFile = vm["optionfile"].as<std::string>();
      std::cout << "Using optionfile: " << optionFile << ".\n";
    } else {
      std::cout << "No option file specified.\nTrying default option file: ../input/options.xml\n";
      optionFile = "../input/options.xml";
    }

    if (vm.count("threads")) {
      std::cout << "Using " << vm["threads"].as<int>() << " threads.\n";
      nrOfThreads = vm["threads"].as<int>();
    } else {
      std::cout << "Default number of threads (1) will be used.\n";
    }
  } catch (std::exception &e) {
    std::cerr << "error: " << e.what() << "\n";
    exit(EXIT_FAILURE);
  } catch (...) {
    std::cerr << "Exception of unknown type!\n";
    exit(EXIT_FAILURE);
  }

  /*****************************************************************/
  /* Next we read the simulation options and setup the simulation. */
  /*****************************************************************/

  // Parse option file
  boost::property_tree::ptree pt, options;
  read_xml(optionFile, pt);
  options = pt.get_child("options");

  // Load Topology
  Topology topol(options.get<double>("kBT"), options.get<double>("electricField_x"));
  topol.readTopologyFromFile(options.get<std::string>("pathToSites"));
  std::cout << "Loaded " << topol.getNrOfSites() << " sites from "
            << options.get<std::string>("pathToSites") << std::endl;

  topol.readReorganisationEnergies(options.get<std::string>("pathToLambdas"));

  // Load Neighbourlist and precompute hopping rates
  std::cout << "Loading neighbours from file this may take a while ...\n";
  Neighbourlist nbList(topol.getNrOfSites());
  nbList.setupShortRangeNeighbours(options.get<std::string>("pathToShortNB"),
                                   topol);
  nbList.setupLongRangeNeighbours(options.get<std::string>("pathToLongNB"),
                                  topol);
  std::cout << "Loaded neighbours from: "
            << options.get<std::string>("pathToShortNB") << " and "
            << options.get<std::string>("pathToLongNB") << std::endl;

  // Redesign KMC run to account for precomputed hopping

  omp_set_num_threads(nrOfThreads);

  #pragma omp parallel for
  for(int i =0; i < 10; i++){
    runSimulation(topol, nbList, i);
  }


  // Run different instances of the same simulation

  // Combine and merge results

  PBC pbc(options.get<double>("Xmax"), options.get<double>("Ymax"),
          options.get<double>("Zmax"));

}

int main(int ac, char *av[]) {

  setupAndExecuteSimulation(ac, av);

  return 0;
}
