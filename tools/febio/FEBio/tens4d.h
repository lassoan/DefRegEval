#pragma once

#include "mat3d.h"

//-----------------------------------------------------------------------------
//! class for 4th order tensors with major and minor symmetries.

// Due to the major symmetry we can store this tensor as a 6x6 matrix.
// The tensor is stored in column major order:
//
//     / 0   1   3   6   10   15  \
//     |     2   4   7   11   16  |
//     |         5   8   12   17  |
// A = |             9   13   18  |
//     |                 14   19  |
//     \                      20  /
//
// Note that due to the minor symmetry we only store the upper matrix of this tensor
//

class tens4ds
{
public:
	enum { NNZ = 21 };

	// default constructor
	tens4ds(){}

	tens4ds(double m[6][6])
	{
		d[ 0] = m[0][0];
		d[ 1] = m[0][1]; d[ 2] = m[1][1];
		d[ 3] = m[0][2]; d[ 4] = m[1][2]; d[ 5] = m[2][2];
		d[ 6] = m[0][3]; d[ 7] = m[1][3]; d[ 8] = m[2][3]; d[ 9] = m[3][3];
		d[10] = m[0][4]; d[11] = m[1][4]; d[12] = m[2][4]; d[13] = m[3][4]; d[14] = m[4][4];
		d[15] = m[0][5]; d[16] = m[1][5]; d[17] = m[2][5]; d[18] = m[3][5]; d[19] = m[4][5]; d[20] = m[5][5];
	}

	double& operator () (int i, int j)
	{
		const int m[6] = {0, 1, 3, 6, 10, 15};
		if (i<=j) return d[m[j]+i]; else return d[m[i]+j];
	}

	double operator () (int i, int j) const
	{
		const int m[6] = {0, 1, 3, 6, 10, 15};
		if (i<=j) return d[m[j]+i]; else return d[m[i]+j];
	}

	// arithmetic operators
	tens4ds operator + (const tens4ds& t) const;
	tens4ds operator - (const tens4ds& t) const;
	tens4ds operator * (double g) const;
	tens4ds operator / (double g) const;

	// arithmetic assignment operators
	tens4ds& operator += (const tens4ds& t);
	tens4ds& operator -= (const tens4ds& t);
	tens4ds& operator *= (double g);
	tens4ds& operator /= (double g);

	// intialize to zero
	void zero();

	// extract 6x6 matrix
	void extract(double d[6][6]);

public:
	double d[NNZ];	// stored in column major order
};

//! Check positive definiteness of a 4th-order symmetric tensor
bool IsPositiveDefinite(const tens4ds& t);

// outer (dyadic) products for symmetric matrices
tens4ds dyad1s(const mat3ds& a);
tens4ds dyad1s(const mat3ds& a, const mat3ds& b);
tens4ds dyad4s(const mat3ds& a);
tens4ds dyad4s(const mat3ds& a, const mat3ds& b);

inline tens4ds operator * (const double g, const tens4ds& a) { return a*g; }

// The following file contains the actual definition of the class functions
#include "tens4d.hpp"

