/***************************************************
 *
 * KMC MODEL FOR OPTOELECTRIC PROCESSES
 *
 * Author: Ruben Gerritsen
 *
 * Created on 14-01-2020
 *
 **************************************************/

#include "Site.h"
#include <iostream>

int Site::isOccupiedBy(PType type) const { 
	if (occupied[type]){
		return occupiedBy[type];
		}else { 
			std::cout << "This site is no longer occupied!" << std::endl;
			return occupiedBy[type];
		}
}