/***************************************************
 *
 * KMC MODEL FOR OPTOELECTRIC PROCESSES
 *
 * Author: Ruben Gerritsen
 *
 * Created on 14-01-2020
 *
 * This file contains the global enum definitions.
 *
 **************************************************/

#ifndef ENUMNAMES
#define ENUMNAMES

enum PType { elec = 0, hole = 1, sing = 2, CT = 3, trip = 4 };

enum MType { ben = 0, tcne = 1 };

enum Transition {
  normalhop = 0,
  decay,
  excitonFromElec,
  excitonFromHole,
  excitonFromElecCT,
  excitonFromHoleCT,
  singToCTViaElec,
  singToCTViaHole,
  tripToCTViaElec,
  tripToCTViaHole,
  CTdisViaHole,
  CTdisViaElec,
  done
};

#endif