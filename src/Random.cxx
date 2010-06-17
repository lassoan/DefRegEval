#include "Random.h"

Random::Random()
{
  srand((unsigned int)time(NULL));
  rand();
}


//these functions get a random number in a given interval (between "begin" and "end" inclusively)
double Random::randomDouble(double begin, double end)
{
	double r = rand() / ((double)RAND_MAX+1);
	return (r * (end - begin)) + begin;
}

float Random::randomFloat(float begin, float end)
{
	float r = ((float)rand() / ((float)(RAND_MAX)+(float)(1)) );
	return (r * (end - begin)) + begin;
}

int Random::randomInt(int begin, int end)
{
	return ((rand() % (end - begin + 1)) + begin);
}

//this function flips a coin that lands heads with probability p.
//true = heads, false = tails
bool Random::flipCoin(double p)
{
	return (randomDouble(0, 1) <= p);
}
