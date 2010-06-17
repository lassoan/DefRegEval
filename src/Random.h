//A wrapper class for C++'s rand() function. Seeded with the system's current time upon instantiation

#ifndef RANDOM
#define RANDOM

#include <stdlib.h>
#include <time.h>

class Random
{
public:
	Random();

	float randomFloat(float, float);
	double randomDouble(double, double);
	int randomInt(int, int);
	bool flipCoin(double);
};


#endif