#ifndef CONTOUR_STATISTICS_H
#define CONTOUR_STATISTICS_H

//This class stores information about each slice in the contour such as range in
//voxel space (max and min), gravity center, and surface points. It is a helper
//class, used intermittently throughout the main workflow. It is mostly used for
//generating the targeting error heat map.

#include <vector>
#include <math.h>
#include "Random.h"

class ContourStatistics	
{
public:
	struct MinMaxPair
	{
		int min;
		int max;
	};

	struct Point
	{
		int x;
		int y;
		int z;
	};

private:
	int numSlices;
	int minX;
	int maxX;
	int minY;
	int maxY;
	int minZ;
	int maxZ;
	bool minMaxZXYSet;
	bool minMaxYXZSet;
	bool minMaxXYZSet;
	std::vector<int> minXForSlices;
	std::vector<int> maxXForSlices;
	std::vector<int> minYForSlices;
	std::vector<int> maxYForSlices;
	std::vector<int> gravityXForSlices;
	std::vector<int> gravityYForSlices;
	MinMaxPair** minMaxZXY;
	MinMaxPair** minMaxYXZ;
	MinMaxPair** minMaxXYZ;
	std::vector<std::vector<Point>> surfacePoints;

	void quickSortSurfacePointsByAngle(std::vector<Point>&, int);

	void destroy()
	{
		int xRange = maxX - minX + 1;
		int yRange = maxY - minY + 1;

		if(minMaxZXYSet)
		{
			for(int x=0; x<xRange; x++)
			{
				delete [] minMaxZXY[x];
			}
			delete [] minMaxZXY;
		}

		if(minMaxYXZSet)
		{
			for(int x=0; x<xRange; x++)
			{
				delete [] minMaxYXZ[x];
			}
			delete [] minMaxYXZ;
		}

		if(minMaxXYZSet)
		{
			for(int y=0; y<yRange; y++)
			{
				delete [] minMaxXYZ[y];
			}
			delete [] minMaxXYZ;
		}
	}

	void copy(const ContourStatistics& statistics)
	{
		numSlices = statistics.numSlices;
		minX = statistics.minX;
		maxX = statistics.maxX;
		minY = statistics.minY;
		maxY = statistics.maxY;
		minZ = statistics.minZ;
		maxZ = statistics.maxZ;
		minMaxZXYSet = statistics.minMaxZXYSet;
		minMaxYXZSet = statistics.minMaxYXZSet;
		minMaxXYZSet = statistics.minMaxXYZSet;

		minXForSlices = statistics.minXForSlices;
		maxXForSlices = statistics.maxXForSlices;
		minYForSlices = statistics.minYForSlices;
		maxYForSlices = statistics.maxYForSlices;
		gravityXForSlices = statistics.gravityXForSlices;
		gravityYForSlices = statistics.gravityYForSlices;

		surfacePoints = statistics.surfacePoints;

		int xRange = maxX - minX + 1;
		int yRange = maxY - minY + 1;
		int zRange = maxZ - minZ + 1;

		if(minMaxZXYSet)
		{
			minMaxZXY = new MinMaxPair*[xRange];
			for(int x=0; x<xRange; x++)
			{
				minMaxZXY[x] = new MinMaxPair[yRange];

				for(int y=0; y<yRange; y++)
				{
					minMaxZXY[x][y] = statistics.minMaxZXY[x][y];
				}
			}
		}

		if(minMaxYXZSet)
		{
			minMaxYXZ = new MinMaxPair*[xRange];
			for(int x=0; x<xRange; x++)
			{
				minMaxYXZ[x] = new MinMaxPair[zRange];

				for(int z=0; z<zRange; z++)
				{
					minMaxYXZ[x][z] = statistics.minMaxYXZ[x][z];
				}
			}
		}

		if(minMaxXYZSet)
		{
			minMaxXYZ = new MinMaxPair*[yRange];
			for(int y=0; y<yRange; y++)
			{
				minMaxXYZ[y] = new MinMaxPair[zRange];

				for(int z=0; z<zRange; z++)
				{
					minMaxXYZ[y][z] = statistics.minMaxXYZ[y][z];
				}
			}
		}
	}

public:

	ContourStatistics()
	{
		minMaxZXYSet = false;
		minMaxYXZSet = false;
		minMaxXYZSet = false;
	}
	ContourStatistics(const ContourStatistics& statistics)
	{
		copy(statistics);
	}
	~ContourStatistics()
	{
		destroy();
	}

	//getters
	int getNumSlices();
	int getMinX();
	int getMaxX();
	int getMinY();
	int getMaxY();
	int getMinZ();
	int getMaxZ();
	int getMinX(int, int);
	int getMaxX(int, int);
	int getMinY(int, int);
	int getMaxY(int, int);
	int getMinZ(int, int);
	int getMaxZ(int, int);
	int getMinXForSlice(int);
	int getMaxXForSlice(int);
	int getMinYForSlice(int);
	int getMaxYForSlice(int);
	int getGravityXForSlice(int);
	int getGravityYForSlice(int);
	std::vector<Point> getSurfacePointsForSlice(int);

	//setters
	void setNumSlices(int);
	void setMinX(int);
	void setMaxX(int);
	void setMinY(int);
	void setMaxY(int);
	void setMinZ(int);
	void setMaxZ(int);
	void setMinMaxZXY(MinMaxPair**);
	void setMinMaxYXZ(MinMaxPair**);
	void setMinMaxXYZ(MinMaxPair**);
	void setMinXForSlices(std::vector<int>);
	void setMaxXForSlices(std::vector<int>);
	void setMinYForSlices(std::vector<int>);
	void setMaxYForSlices(std::vector<int>);
	void setSurfacePoints(std::vector<std::vector<Point>>);
	void addMinXForSlices(int);
	void addMaxXForSlices(int);
	void addMinYForSlices(int);
	void addMaxYForSlices(int);

	//other functions
	void determineGravityCenterForSlices();
	void operator=(const ContourStatistics&);
	double angleOfSurfacePoint(Point, int);
};

#endif