#include "vtkKWProstateErrorMapRenderingWidget.h"

////////////////////////////////Getters////////////////////////////////////////

int Polynomial::getDegree()
{
	return degree;
}

float Polynomial::getCoefficient(int c)
{
	return coefficients[c];
}

////////////////////////////////Setters////////////////////////////////////////

void Polynomial::setCoefficient(int c, float val)
{
	coefficients[c] = val;
}

////////////////////////////////Main Functions////////////////////////////////////////

//This function evaluates the polynomial at x
float Polynomial::evaluate(float x)
{
	float result = coefficients[0];

	for(int i=1; i<=degree; i++)
	{
		result += coefficients[i]*power(x, i);
	}

	return result;
}

//This function returns the derivative of the polynomial
Polynomial Polynomial::derivative()
{
	if(degree == 0)
	{
		Polynomial p(0);
		p.setCoefficient(0, 0);
		return p;
	}
	else
	{
		Polynomial p(degree - 1);

		for(int i=(degree - 1); i>=0; i--)
		{
			p.setCoefficient(i, (i+1)*coefficients[i+1]);
		}

		return p;
	}
}

//This function returns the antiderivative of the polynomial
Polynomial Polynomial::antiderivative()
{
	Polynomial p(degree + 1);

	for(int i=(degree + 1); i>=1; i--)
	{
		p.setCoefficient(i, (float)coefficients[i-1]/(float)i);
	}

	p.setCoefficient(0, 0);

	return p;
}

Polynomial& Polynomial::operator=(const Polynomial& p)
{
	degree = p.degree;

	delete [] coefficients;
	coefficients = new float[degree + 1];

	for(int i=0; i<(degree + 1); i++)
	{
		coefficients[i] = p.coefficients[i];
	}

	return *this;
}

std::ostream& operator<<(std::ostream& out, const Polynomial& p)
{
	out<<"y = ";

	int numNotZero = 0;

	for(int i=p.degree; i>=0; i--)
	{
		if(p.coefficients[i] != 0)
		{
			if(numNotZero > 0)
			{
				if(p.coefficients[i] > 0)
				{
					out<<" + ";
				}
				else
				{
					out<<" ";
				}
			}

			numNotZero++;

			if(p.coefficients[i] != 1)
			{
				out<<p.coefficients[i];
			}
		}

		if((p.coefficients[i] != 0) && (i != 0))
		{
			out<<"x";
			if(i != 1)
			{
				out<<"^"<<i;
			}
		}	
	}

	return out;
}