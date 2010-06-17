#include "vtkKWProstateErrorMapRenderingWidget.h"

int CubicSplineSolver::getNumPoints()
{
	return numPoints;
}

//"t" is a value used to clean up the code in the "setupEquations" function
float CubicSplineSolver::t(int i)
{
	return 6 * ((float)(yVals[i + 1] - yVals[i])/(float)h[i] - (float)(yVals[i] - yVals[i - 1])/(float)h[i - 1]);
}

void CubicSplineSolver::setupEquations()
{
	//place all values in a matrix to solve the system of equations in:
	//http://en.wikipedia.org/wiki/Spline_interpolation#Interpolation_using_natural_cubic_spline

	int n = numPoints - 1;

	equationMatrix = new float*[n - 1];
	for(int i=0; i<(n - 1); i++)
	{
		equationMatrix[i] = new float[n];
	}

	for(int i=0; i<(n - 1); i++)
	{
		for(int j=0; j<n; j++)
		{
			equationMatrix[i][j] = 0;
		}
	}

	equationMatrix[0][0] = 2*(h[0] + h[1]);
	equationMatrix[0][1] = h[1];
	equationMatrix[0][n - 1] = t(1);

	equationMatrix[n - 2][n - 3] = h[n - 2];
	equationMatrix[n - 2][n - 2] = 2*(h[n - 2] + h[n - 1]);
	equationMatrix[n - 2][n - 1] = t(n - 1);

	for(int i=1; i<=(n-3); i++)
	{
		equationMatrix[i][i - 1] = h[i];
		equationMatrix[i][i] = 2*(h[i] + h[i + 1]);
		equationMatrix[i][i + 1] = h[i + 2];
		equationMatrix[i][n - 1] = t(i + 1);
	}
}

void CubicSplineSolver::solveEquations()
{
	int n = numPoints - 1;

	//row reduce
	for(int i=1; i<=(n - 2); i++)
	{
		float a = equationMatrix[i - 1][i - 1];
		float b = equationMatrix[i][i - 1];

		for(int j=i-1; j<n; j++)
		{
			equationMatrix[i][j] = a*equationMatrix[i][j] - b*equationMatrix[i-1][j];
		}

		for(int j=i-1; j<n; j++)
		{
			equationMatrix[i][j] = 4.5*equationMatrix[i][j];
		}
	}

	z[n - 1] = (float)equationMatrix[n-2][n-1]/(float)equationMatrix[n-2][n-2];

	//back substitution
	for(int i=n-3; i>=0; i--)
	{
		z[i + 1] = (float)(equationMatrix[i][n-1] - equationMatrix[i][i + 1]*z[i + 2])/(float)equationMatrix[i][i];
	}
}

//This function uses the coefficients z and h to generate the natural cubic
//spline polynomial at index i.
Polynomial CubicSplineSolver::getCubicSpline(int i)
{
	Polynomial p(3);

	//The coefficients for the polynomial are a closed-form of the polynomial Si(x) in:
	//http://en.wikipedia.org/wiki/Spline_interpolation#Interpolation_using_natural_cubic_spline

	p.setCoefficient(3, (float)(z[i+1]-z[i])/(float)(6*h[i]));
	p.setCoefficient(2, (float)(3*xVals[i+1]*z[i])/(float)(6*h[i]) - (float)(3*xVals[i]*z[i+1])/(float)(6*h[i]));
	p.setCoefficient(1, (float)(3*z[i+1]*xVals[i]*xVals[i])/(float)(6*h[i]) - (float)(3*z[i]*xVals[i+1]*xVals[i+1])/(float)(6*h[i]) 
						 + (float)yVals[i+1]/(float)h[i] - (float)(h[i]*z[i+1])/(float)6 + (float)(h[i]*z[i])/(float)6 - (float)yVals[i]/(float)h[i]);
	p.setCoefficient(0, (float)(z[i]*power(xVals[i+1],3))/(float)(6*h[i]) - (float)(z[i+1]*power(xVals[i],3))/(float)(6*h[i]) +
						(float)(h[i]*z[i+1]*xVals[i])/(float)6 - (float)(yVals[i+1]*xVals[i])/(float)h[i] + (float)(yVals[i]*xVals[i+1])/(float)h[i] -
						(float)(h[i]*xVals[i+1]*z[i])/(float)6);

	return p;
}

//This function evaluates for x using the corresponding spline curve
//in the natural cubic spline interpolation of the initial points
float CubicSplineSolver::evaluate(float x)
{
	//first, find which spline curve corresponds to this
	//x value
	int splineIndex;

	if(x < xVals[0])
	{
		splineIndex = 0;
	}
	else if(x > xVals[numPoints - 1])
	{
		splineIndex = numPoints - 2;
	}
	else
	{
		splineIndex = 0;

		for(int i=0; i<(numPoints - 1); i++)
		{
			if((x >= xVals[i]) && (x < xVals[i + 1]))
			{
				break;
			}
			else
			{
				splineIndex++;
			}
		}
	}

	//next, evaluate that x value with the given spline curve
	Polynomial splineFunction = getCubicSpline(splineIndex);
	return splineFunction.evaluate(x);
}