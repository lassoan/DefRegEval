#ifndef CUBIC_SPLINE_SOLVER
#define CUBIC_SPLINE_SOLVER

//This class solves for the natural cubic spline interpolation between
//a given set of points. The solution is an implementation of:
//http://en.wikipedia.org/wiki/Spline_interpolation#Interpolation_using_natural_cubic_spline

//The solution involves solving a system of equations for the coefficients called "z".
//The coefficients, z, as well as other coefficients, h, are used in a formula to generate
//the spline curves.

class CubicSplineSolver
{
private:
	float* xVals;
	float* yVals;
	float* h;
	float* z; //'h' and 'z' are the coefficients needed to solve for the individual spline curves
	float** equationMatrix; //matrix for solving the system of equations using row reduction
	int numPoints;

	void setupEquations();
	void solveEquations();
	float t(int);

public:
	CubicSplineSolver(std::vector<float>& x, std::vector<float>& y)
	{
		//first remove any duplicate values in the x-coordinates
		std::vector<int> duplicateValueIndeces;

		for(int i=0; i<(x.size() - 1); i++)
		{
			if(x.at(i) == x.at(i + 1))
			{
				duplicateValueIndeces.push_back(i);
			}
		}

		numPoints = x.size() - duplicateValueIndeces.size();;
		xVals = new float[numPoints];
		yVals = new float[numPoints];

		if(duplicateValueIndeces.size() > 0)
		{
			int removes = 0;
			int i = 0;

			while((removes < duplicateValueIndeces.size()) && (i < x.size()))
			{
				if(i == duplicateValueIndeces.at(removes))
				{
					removes++;
				}
				else
				{
					xVals[i - removes] = x.at(i);
					yVals[i - removes] = y.at(i);
				}
				i++;
			}

			if(i < x.size())
			{
				for(int j=i; j<x.size(); j++)
				{
					xVals[j - removes] = x.at(j);
					yVals[j - removes] = y.at(j);
				}
			}
		}
		else
		{
			for(int i=0; i<numPoints; i++)
			{
				xVals[i] = x.at(i);
				yVals[i] = y.at(i);
			}
		}

		//setup the coefficients
		h = new float[numPoints - 1];
		z = new float[numPoints];

		//h is easy to solve for - it's just the
		//difference in the x-component
		for(int i=0; i<(numPoints-1); i++)
		{
			h[i] = xVals[i + 1] - xVals[i];
		}

		//z must be solved for in a system of equations
		for(int i=0; i<numPoints; i++)
		{
			z[i] = 0;
		}

		if(numPoints >= 3)
		{
			if(numPoints == 3) //z can be solved for explicitly
			{
				z[0] = 0;
				z[2] = 0;

				z[1] = (float)t(1)/(float)(2*(h[0] + h[1]));
			}
			else if(numPoints == 4) //z can be solved for explicitly
			{
				z[0] = 0;
				z[3] = 0;

				z[2] = ((float)h[1]/(float)(-4*(h[0]+h[1])*(h[1]+h[2]) + 2*h[1]))*(t(1) - (float)(2*(h[0]+h[1])*t(2))/(float)h[1]);
				z[1] = (float)(t(2) - 2*(h[1]+h[2])*z[2])/(float)h[1];
			}
			else
			{
				//must solve the system of equations
				setupEquations();
				solveEquations();
			}
		}
	}

	~CubicSplineSolver()
	{
		delete [] h;
		delete [] z;
		delete [] xVals;
		delete [] yVals;

		if(numPoints > 4)
		{
			for(int i=0; i<(numPoints - 2); i++)
			{
				delete [] equationMatrix[i];
			}

			delete [] equationMatrix;
		}	
	}

	int getNumPoints();
	Polynomial getCubicSpline(int);
	float evaluate(float);
};

#endif