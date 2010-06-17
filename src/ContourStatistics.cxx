#include "ContourStatistics.h"

////////////////////////////////////GETTERS////////////////////////////////////////

int ContourStatistics::getNumSlices()
{
	return numSlices;
}

int ContourStatistics::getMinX()
{
	return minX;
}

int ContourStatistics::getMaxX()
{
	return maxX;
}

int ContourStatistics::getMinY()
{
	return minY;
}

int ContourStatistics::getMaxY()
{
	return maxY;
}

int ContourStatistics::getMinZ()
{
	return minZ;
}

int ContourStatistics::getMaxZ()
{
	return maxZ;
}

int ContourStatistics::getMinX(int y, int z)
{
	return minMaxXYZ[y - minY][z - minZ].min;
}

int ContourStatistics::getMaxX(int y, int z)
{
	return minMaxXYZ[y - minY][z - minZ].max;
}
int ContourStatistics::getMinY(int x, int z)
{
	return minMaxYXZ[x - minX][z - minZ].min;
}

int ContourStatistics::getMaxY(int x, int z)
{
	return minMaxYXZ[x - minX][z - minZ].max;
}

int ContourStatistics::getMinZ(int x, int y)
{
	return minMaxZXY[x - minX][y - minY].min;
}

int ContourStatistics::getMaxZ(int x, int y)
{
	return minMaxZXY[x - minX][y - minY].max;
}

int ContourStatistics::getMinXForSlice(int i)
{
	return minXForSlices.at(i);
}

int ContourStatistics::getMaxXForSlice(int i)
{
	return maxXForSlices.at(i);
}

int ContourStatistics::getMinYForSlice(int i)
{
	return minYForSlices.at(i);
}

int ContourStatistics::getMaxYForSlice(int i)
{
	return maxYForSlices.at(i);
}

int ContourStatistics::getGravityXForSlice(int i)
{
	return gravityXForSlices.at(i);
}

int ContourStatistics::getGravityYForSlice(int i)
{
	return gravityYForSlices.at(i);
}

std::vector<ContourStatistics::Point> ContourStatistics::getSurfacePointsForSlice(int i)
{
	return surfacePoints.at(i);
}

////////////////////////////////////SETTERS////////////////////////////////////////

void ContourStatistics::setNumSlices(int s)
{
	numSlices = s;
}

void ContourStatistics::setMinX(int x)
{
	minX = x;
}

void ContourStatistics::setMaxX(int x)
{
	maxX = x;
}

void ContourStatistics::setMinY(int y)
{
	minY = y;
}

void ContourStatistics::setMaxY(int y)
{
	maxY = y;
}

void ContourStatistics::setMinZ(int z)
{
	minZ = z;
}

void ContourStatistics::setMaxZ(int z)
{
	maxZ = z;
}

void ContourStatistics::setMinMaxZXY(MinMaxPair** minMax)
{
	minMaxZXYSet = true;
	int xRange = maxX - minX + 1;
	int yRange = maxY - minY + 1;

	minMaxZXY = new MinMaxPair*[xRange];
	for(int x=0; x<xRange; x++)
	{
		minMaxZXY[x] = new MinMaxPair[yRange];

		for(int y=0; y<yRange; y++)
		{
			minMaxZXY[x][y] = minMax[x][y];
		}
	}
}

void ContourStatistics::setMinMaxYXZ(MinMaxPair** minMax)
{
	minMaxYXZSet = true;
	int xRange = maxX - minX + 1;
	int zRange = maxZ - minZ + 1;

	minMaxYXZ = new MinMaxPair*[xRange];
	for(int x=0; x<xRange; x++)
	{
		minMaxYXZ[x] = new MinMaxPair[zRange];

		for(int z=0; z<zRange; z++)
		{
			minMaxYXZ[x][z] = minMax[x][z];
		}
	}
}

void ContourStatistics::setMinMaxXYZ(MinMaxPair** minMax)
{
	minMaxXYZSet = true;
	int yRange = maxY - minY + 1;
	int zRange = maxZ - minZ + 1;

	minMaxXYZ = new MinMaxPair*[yRange];
	for(int y=0; y<yRange; y++)
	{
		minMaxXYZ[y] = new MinMaxPair[zRange];

		for(int z=0; z<zRange; z++)
		{
			minMaxXYZ[y][z] = minMax[y][z];
		}
	}
}

void ContourStatistics::setMinXForSlices(std::vector<int> mins)
{
	minXForSlices = mins;
}

void ContourStatistics::setMaxXForSlices(std::vector<int> maxes)
{
	maxXForSlices = maxes;
}

void ContourStatistics::setMinYForSlices(std::vector<int> mins)
{
	minYForSlices = mins;
}

void ContourStatistics::setMaxYForSlices(std::vector<int> maxes)
{
	maxYForSlices = maxes;
}

void ContourStatistics::addMinXForSlices(int min)
{
	minXForSlices.push_back(min);
}

void ContourStatistics::addMaxXForSlices(int max)
{
	maxXForSlices.push_back(max);
}

void ContourStatistics::addMinYForSlices(int min)
{
	minYForSlices.push_back(min);
}

void ContourStatistics::addMaxYForSlices(int max)
{
	maxYForSlices.push_back(max);
}

void ContourStatistics::setSurfacePoints(std::vector<std::vector<ContourStatistics::Point>> points)
{
	surfacePoints = points;

	for(int z=0; z<numSlices; z++)
	{
		std::vector<Point> surfacePointsForSlice = surfacePoints.at(z);
		quickSortSurfacePointsByAngle(surfacePointsForSlice, z);
		surfacePoints[z] = surfacePointsForSlice;
	}
}

////////////////////////////////////Other Functions////////////////////////////////////////

void ContourStatistics::determineGravityCenterForSlices()
{
	//approximation of the gravity center
	for(int i=0; i<numSlices; i++)
	{
		gravityXForSlices.push_back(minXForSlices[i] + (maxXForSlices[i] - minXForSlices[i])/2);
		gravityYForSlices.push_back(minYForSlices[i] + (maxYForSlices[i] - minYForSlices[i])/2);
	}
}

void ContourStatistics::operator=(const ContourStatistics& statistics)
{
	destroy();
	copy(statistics);
}

//This function determines the angle a surface point makes with respect to
//the gravity center of the slice.
double ContourStatistics::angleOfSurfacePoint(Point surfacePoint, int slice)
{
	float PI = 3.1415926535f;
	int xDif = surfacePoint.x - gravityXForSlices.at(slice);
	int yDif = surfacePoint.y - gravityYForSlices.at(slice);
	float distance = sqrt((float)(xDif*xDif + yDif*yDif));
	float angle = acos((float)abs(xDif)/(float)distance);

	if((xDif > 0) && (yDif >= 0)) //first quadrant
	{
		angle = angle;
	}
	else if((xDif <= 0) && (yDif > 0)) //second quadrant
	{
		angle = PI - angle;
	}
	else if((xDif < 0) && (yDif <= 0)) //third quadrant
	{
		angle = PI + angle;
	}
	else //fourth quadrant
	{
		angle = 2*PI - angle;
	}

	return angle;
}

//This is a randomized quicksort of the surface points for a slice with respect to the angle
//they make with the gravity center. This is used when using surface targeting. The sorted
//points easily allow the removal of the points within the angle ANGLE_NOT_TARGETED specified
//in vtkKWProstateErrorMapRenderingWidget.h
void ContourStatistics::quickSortSurfacePointsByAngle(std::vector<Point>& points, int slice)
{
	if(points.size() <= 1)
	{
		return;
	}
	else if(points.size() == 2)
	{
		Point p1 = points.at(0);
		Point p2 = points.at(1);

		points.clear();

		if(angleOfSurfacePoint(p1, slice) < angleOfSurfacePoint(p2, slice))
		{
			points.push_back(p1);
			points.push_back(p2);
		}
		else
		{
			points.push_back(p2);
			points.push_back(p1);
		}
	}
	else
	{
		Random r;

		int pivotIndex = r.randomInt(0, points.size() - 1);;
		Point pivot = points.at(pivotIndex);
		double pivotAngle = angleOfSurfacePoint(pivot, slice);

		std::vector<Point> points1;
		std::vector<Point> points2;

		for(unsigned int i=0; i<points.size(); i++)
		{
			if(i != pivotIndex)
			{
				if(angleOfSurfacePoint(points.at(i), slice) < pivotAngle)
				{
					points1.push_back(points.at(i));
				}
				else
				{
					points2.push_back(points.at(i));
				}
			}
		}

		quickSortSurfacePointsByAngle(points1, slice);
		quickSortSurfacePointsByAngle(points2, slice);

		points.clear();

		for(unsigned int i=0; i<points1.size(); i++)
		{
			points.push_back(points1.at(i));
		}

		points.push_back(pivot);

		for(unsigned int i=0; i<points2.size(); i++)
		{
			points.push_back(points2.at(i));
		}
	}
}