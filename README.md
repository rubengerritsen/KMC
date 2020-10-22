# KMC (time of flight)
To compile, install and run this you only need to run one command `./run` in the main directory.
You maybe need to make `run` executable with `chomd +x run`.

The executable will be in the build directory, it can be run with the `./run`
command or you can run it from within the build directory and specify the option
file.

```bash
./KMC -o <options>.xml
```

## Parallelism is *NOT* supported atm

## Dependencies
To make it work, make sure you have the Boost and Eigen libraries installed.
These are the same dependencies as VOTCA has. So if VOTCA works, this will work.

## Input files
You need to provide three input files, examples can be found in the `input` directory.

### 1. A file with the sites
It should contain the site locations, every site is put on a new line i.e.

```
<xcoordSite0> <ycoordSite0> <zcoordSite0>
<xcoordSite1> <ycoordSite1> <zcoordSite1>
<xcoordSite2> <ycoordSite2> <zcoordSite2>
.
.
```
Since the coordinates are not used in the actual simulation you can also use fake coordinates, that is what I did in the example file.

### 2. A Neigbour List
A neigbourlist is a file with the neigbours and a rate, if a particle can jump from site 0 to site 5 with rate 0.2 than we would add a line `0 5 0.2` to the neigbourlist file. We need to add another pair for the other way around, so a hop from 5 to 0 could have a different rate, say 0.3, then we add `5 0 0.3`. 

```
<startLocation> <targetLocation> <rate>
<startLocation> <targetLocation> <rate>
.
.
```

### 3. An option file that looks like this

```xml
<?xml version="1.0"?>
<options>
  <sink help="numbering starts at 0!">6</sink>
  <SEED help="SEED for random number generator">12345</SEED>
  <simName help="prefix to output folder, can be used as identifier">sim</simName>
  <maxTime help="picoseconds">10</maxTime>
  <maxStep>100</maxStep>
  <nrOfRunsPerProcess help="total number of runs = nrOfFunsPerProcess*nrOfProcesses">1</nrOfRunsPerProcess>
  <nrOfProcesses>1</nrOfProcesses>
  <Xmax>10</Xmax>
  <Ymax>10</Ymax>
  <Zmax>10</Zmax>
  <nrOfElectrons>0</nrOfElectrons>
  <nrOfHoles>1</nrOfHoles>
  <pathToSites help="path to file with sites and energies">../input/sites.txt</pathToSites>
  <pathToNBList help="path to file with short range neighbours">../input/nbList.txt</pathToNBList>
</options>
```

## What the simulator does
If you specify 3 holes (or electrons), these holes will be placed on the first 3 sites (i.e
site 0, 1 and 2 from the sites.txt file). Then the normal KMC algorithm is run. As soon as a hole reaches the sink (indicated in the option file), the hole is removed from the simulation and the time is registered and printed to a file, this is done until all holes (or electrons) have reached the sink.

## Results
The results are collected in a folder within the build directory, this folder is labeled with the date, time and simulation ID (`simName` in the option file). The file is called `tof.txt` to get a quick idea of the time of flight you can run the script below to get the average time of flight for the first, second, ... particle.

```python
#!/usr/bin/env python3

import numpy as np

data = np.loadtxt("tof.txt")
mean = np.mean(data, axis=0)
std  = np.std(data, axis=0)   
    
np.set_printoptions(formatter={'float': '{: 7.2f}'.format})
print(mean)
print(std)
```