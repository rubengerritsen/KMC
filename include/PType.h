/***************************************************
 *
 * KMC MODEL FOR OPTOELECTRIC PROCESSES
 *
 * Author: Ruben Gerritsen
 *
 * Created on 14-01-2020
 *
 * 
 **************************************************/

#ifndef PARTICLETYPE
#define PARTICLETYPE

enum PType { elec = 0, hole = 1, trip = 2, sing = 3, CT = 4 };

enum Transition { normalhop = 0, decay, excitonFromElec, excitonFromHole, excitonFromElecCT, excitonFromHoleCT, 
                    singToCTViaElec, singToCTViaHole, tripToCTViaElec, tripToCTViaHole, CTdisViaHole, CTdisViaElec};

#endif