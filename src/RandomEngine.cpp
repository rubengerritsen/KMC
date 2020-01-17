#include "RandomEngine.h"

void RandomEngine::initializeParameters(std::array<double, 4> mu, std::array<double, 4> sigma) {
	for (unsigned int i = 0; i < mu.size(); ++i) {
		dos[i] = std::normal_distribution<double>{ mu[i], sigma[i] };
	}
}