//A wrapper class for C++'s rand() function. Seeded with the system's current time upon instantiation

#ifndef RANDOM
#define RANDOM

#include <stdlib.h>
#include <time.h>

class Random
{
public:
	Random::Random()
  {
    srand((unsigned int)time(NULL));
    rand();
  }

  //these functions get a random number in a given interval (between "begin" and "end" inclusively)

  float randomFloat(float begin, float end)
  {
    float r = ((float)rand() / ((float)(RAND_MAX)+(float)(1)) );
    return (r * (end - begin)) + begin;
  }

  double randomDouble(double begin, double end)
  {
    double r = rand() / ((double)RAND_MAX+1);
    return (r * (end - begin)) + begin;
  }

	int randomInt(int begin, int end)
  {
	  return ((rand() % (end - begin + 1)) + begin);
  }

  //this function flips a coin that lands heads with probability p.
  //true = heads, false = tails
  bool flipCoin(double p)
  {
	  return (randomDouble(0, 1) <= p);
  }

};

#endif