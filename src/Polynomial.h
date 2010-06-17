#ifndef POLYNOMIAL
#define POLYNOMIAL

//This class is an abstraction of a polynomial. This is used in the simulator only as
//the cubic polynomials which are solved for by the CubicSplineSolver

class Polynomial
{
friend std::ostream& operator<<(std::ostream&, const Polynomial&);

private:
	int degree; //degree of the polynomial... f(x) = x^3 - 2x + 1 is degree 3
	float* coefficients; //list of coefficients... coefficients[0] is the constant, coefficients[1] is the coefficient of x, ...

public:
	Polynomial(int d)
	{
		degree = d;
		coefficients = new float[degree + 1];

		for(int i=0; i<(degree + 1); i++)
		{
			coefficients[i] = 0;
		}
	}
	Polynomial(const Polynomial& p)
	{
		degree = p.degree;
		coefficients = new float[degree + 1];

		for(int i=0; i<(degree + 1); i++)
		{
			coefficients[i] = p.coefficients[i];
		}
	}
	~Polynomial()
	{
		delete [] coefficients;
	}

	Polynomial& operator=(const Polynomial&);

	//getters
	int getDegree();
	float getCoefficient(int);

	//setters
	void setCoefficient(int, float);

	//main functions
	float evaluate(float);
	Polynomial derivative();
	Polynomial antiderivative();
};

#endif