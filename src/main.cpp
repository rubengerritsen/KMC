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
#include "RandomEngine.h"
#include "RateOptions.h"
#include "SimulationOptions.h"
#include "Topology.h"
#include <array>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <sys/wait.h>
#include <unistd.h>

void runSimulation(const Topology &topol, const Neighbourlist &nbList,
                   SimulationOptions simOptions) {
  KmcRun simulation(topol, nbList, simOptions);
  simulation.runSimulation();
}

void setupAndExecuteSimulation(int ac, char *av[]) {

  /*********************************************/
  /* First we parse the command line arguments */
  /*********************************************/
  std::string optionFile;
  try {
    boost::program_options::options_description desc("Possible options");
    desc.add_options()("help,h", "produce help message")(
        "optionfile,o", boost::program_options::value<std::string>(),
        "path to option file (required)");

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
      std::cout << "No option file specified.\nTrying default option file: "
                   "../input/options.xml\n";
      optionFile = "../input/options.xml";
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
  Topology topol;
  topol.readTopologyFromFile(options.get<std::string>("pathToSites"));
  std::cout << "Loaded " << topol.getNrOfSites() << " sites from "
            << options.get<std::string>("pathToSites") << std::endl;

  // Set simulation options
  SimulationOptions simOptions;
  simOptions.maxTime = options.get<double>("maxTime");
  simOptions.nrOfElectrons = options.get<int>("nrOfElectrons");
  simOptions.nrOfHoles = options.get<int>("nrOfHoles");
  simOptions.SEED = options.get<int>("SEED");
  simOptions.maxStep = options.get<int>("maxStep");
  simOptions.sink = options.get<int>("sink");

  int nrOfProcesses = options.get<int>("nrOfProcesses");
  int nrOfRuns = options.get<int>("nrOfRunsPerProcess");

  // Load Neighbourlist and precompute hopping rates
  std::cout << "Loading neighbours from file this may take a while ...\n";
  Neighbourlist nbList(topol.getNrOfSites());
  nbList.setupNeighbours(options.get<std::string>("pathToNBList"), topol);
  std::cout << "Loaded neighbours from: "
            << options.get<std::string>("pathToNBList") << std::endl;

  // setup output folder
  struct tm *ltm;
  time_t now = time(0);
  ltm = localtime(&now);
  ltm->tm_mon = ltm->tm_mon + 1;
  ltm->tm_year = ltm->tm_year + 1900;
  std::string foldername =
      "./" +
      (boost::format("%04d%02d%02d_%02d%02d%02d_") % ltm->tm_year %
       ltm->tm_mon % ltm->tm_mday % ltm->tm_hour % ltm->tm_min % ltm->tm_sec)
          .str() +
      options.get<std::string>("simName");

  boost::filesystem::path dir(foldername);
  if (boost::filesystem::create_directory(dir)) {
    std::cout << "Successfully created output directory. \n";
  } else {
    std::cout
        << "Was not able to create output directory, terminating program\n";
    exit(EXIT_FAILURE);
  }

  boost::filesystem::copy_file(optionFile, foldername + "/options.xml");

  simOptions.outputPath = foldername + "/";

  /**************************************************/
  /* Now we run the actual simulations              */
  /**************************************************/
  // Some C trickery to run parallel simulations
  pid_t pid, wpid;
  int status = 0;

  for (int prs = 0; prs < nrOfProcesses; prs++) {
    pid = fork();
    if (pid == 0) {
      // child
      for (int i = 0; i < nrOfRuns; i++) {
        simOptions.simID = nrOfRuns * prs + i;
        // distort SEED with simID
        simOptions.SEED = simOptions.SEED + simOptions.simID;
        runSimulation(topol, nbList, simOptions);
      }
      exit(EXIT_SUCCESS);
    } else if (pid < 0) {
      // failure
      std::cout << "Fork failed! \n";
      exit(EXIT_FAILURE);
    }
  }
  // This is only executed in the parent process
  while ((wpid = wait(&status)) > 0) { // waiting for all child processes
    ;
  }

  std::cout << "Successfully completed " << nrOfRuns * nrOfProcesses
            << " runs of the simulation.\n";
}

int main(int ac, char *av[]) {

  setupAndExecuteSimulation(ac, av);

  return 0;
}
