#!/bin/bash

# set options for full singlets
xmlstarlet edit -L -u "options/allSinglets" -v "1" ../input/options.xml
xmlstarlet edit -L -u "options/registerState" -v "0" ../input/options.xml
xmlstarlet edit -L -u "options/maxTime" -v "1500" ../input/options.xml
xmlstarlet edit -L -u "options/nrOfRunsPerProcess" -v "1" ../input/options.xml
xmlstarlet edit -L -u "options/nrOfProcesses" -v "25" ../input/options.xml
xmlstarlet edit -L -u "options/nrOfHoles" -v "0" ../input/options.xml
xmlstarlet edit -L -u "options/nrOfElectrons" -v "0" ../input/options.xml
xmlstarlet edit -L -u "options/nrOfSinglets" -v "640" ../input/options.xml
xmlstarlet edit -L -u "options/rates/dielectricConstant" -v " 8.0" ../input/options.xml
# run
./KMC

# set options for full singlets reduced CT binding
xmlstarlet edit -L -u "options/allSinglets" -v "1" ../input/options.xml
xmlstarlet edit -L -u "options/registerState" -v "0" ../input/options.xml
xmlstarlet edit -L -u "options/maxTime" -v "1500" ../input/options.xml
xmlstarlet edit -L -u "options/nrOfRunsPerProcess" -v "1" ../input/options.xml
xmlstarlet edit -L -u "options/nrOfProcesses" -v "25" ../input/options.xml
xmlstarlet edit -L -u "options/nrOfHoles" -v "0" ../input/options.xml
xmlstarlet edit -L -u "options/nrOfElectrons" -v "0" ../input/options.xml
xmlstarlet edit -L -u "options/nrOfSinglets" -v "640" ../input/options.xml
xmlstarlet edit -L -u "options/rates/dielectricConstant" -v " 32.0" ../input/options.xml
# run
./KMC


# set options for 10% singlets
xmlstarlet edit -L -u "options/allSinglets" -v "0" ../input/options.xml
xmlstarlet edit -L -u "options/registerState" -v "0" ../input/options.xml
xmlstarlet edit -L -u "options/maxTime" -v "1500" ../input/options.xml
xmlstarlet edit -L -u "options/nrOfRunsPerProcess" -v "1" ../input/options.xml
xmlstarlet edit -L -u "options/nrOfProcesses" -v "25" ../input/options.xml
xmlstarlet edit -L -u "options/nrOfHoles" -v "0" ../input/options.xml
xmlstarlet edit -L -u "options/nrOfElectrons" -v "0" ../input/options.xml
xmlstarlet edit -L -u "options/nrOfSinglets" -v "640" ../input/options.xml
xmlstarlet edit -L -u "options/rates/dielectricConstant" -v " 8.0" ../input/options.xml
# run
./KMC

# set options for 10% singlets reduced CT binding
xmlstarlet edit -L -u "options/allSinglets" -v "0" ../input/options.xml
xmlstarlet edit -L -u "options/registerState" -v "0" ../input/options.xml
xmlstarlet edit -L -u "options/maxTime" -v "1500" ../input/options.xml
xmlstarlet edit -L -u "options/nrOfRunsPerProcess" -v "1" ../input/options.xml
xmlstarlet edit -L -u "options/nrOfProcesses" -v "25" ../input/options.xml
xmlstarlet edit -L -u "options/nrOfHoles" -v "0" ../input/options.xml
xmlstarlet edit -L -u "options/nrOfElectrons" -v "0" ../input/options.xml
xmlstarlet edit -L -u "options/nrOfSinglets" -v "640" ../input/options.xml
xmlstarlet edit -L -u "options/rates/dielectricConstant" -v " 32.0" ../input/options.xml
# run
./KMC

# set options for 10% singlets reduced CT binding and register full state
xmlstarlet edit -L -u "options/allSinglets" -v "0" ../input/options.xml
xmlstarlet edit -L -u "options/registerState" -v "1" ../input/options.xml
xmlstarlet edit -L -u "options/maxTime" -v "1500" ../input/options.xml
xmlstarlet edit -L -u "options/nrOfRunsPerProcess" -v "1" ../input/options.xml
xmlstarlet edit -L -u "options/nrOfProcesses" -v "5" ../input/options.xml
xmlstarlet edit -L -u "options/nrOfHoles" -v "0" ../input/options.xml
xmlstarlet edit -L -u "options/nrOfElectrons" -v "0" ../input/options.xml
xmlstarlet edit -L -u "options/nrOfSinglets" -v "640" ../input/options.xml
xmlstarlet edit -L -u "options/rates/dielectricConstant" -v " 32.0" ../input/options.xml
# run
./KMC