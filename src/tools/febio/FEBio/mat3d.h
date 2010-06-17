// mat3d.h: interface for the mat3d class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_MAT3D_H__35D54918_61A9_4E11_A4B2_6D3D33AB98FE__INCLUDED_)
#define AFX_MAT3D_H__35D54918_61A9_4E11_A4B2_6D3D33AB98FE__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

// we'll need the vector class
#include "vec3d.h"

//-----------------------------------------------------------------------------
// The following classes are defined in this file
class mat3d;	// general 3D matrix of doubles
class mat3ds;	// symmetric 3D matrix of doubles
class mat3da;	// anti-symmetric 3D matrix of doubles
class mat3dd;	// diagonal matrix of doubles

//-----------------------------------------------------------------------------
//! This class describes a diagonal matrix of doubles in 3D

class mat3dd
{
public:
	// default constructor
	mat3dd(){}

	// constructors
	mat3dd(double a);
	mat3dd(double a0, double a1, double a2);

	// assignment operators
	mat3dd& operator = (const mat3dd& m);
	mat3dd& operator = (double a);

	// access operators
	double operator () (int i, int j) const;
	double& diag(int i);
	const double& diag(int i) const;

	// arithmetic operators
	mat3dd operator + (const mat3dd& m) const;
	mat3dd operator - (const mat3dd& m) const;
	mat3dd operator * (const mat3dd& m) const;
	mat3dd operator * (double a) const;
	mat3dd operator / (double a) const;

	// arithmetic operators for mat3ds
	mat3ds operator + (const mat3ds& m) const;
	mat3ds operator - (const mat3ds& m) const;
	mat3ds operator * (const mat3ds& m) const;

	// arithmetic operators for mat3d
	mat3d operator + (const mat3d& m) const;
	mat3d operator - (const mat3d& m) const;
	mat3d operator * (const mat3d& m) const;

	// arithmetic operators for mat3da const;
	mat3d operator + (const mat3da& m) const;
	mat3d operator - (const mat3da& m) const;
	mat3d operator * (const mat3da& m) const;

	// arithmetic assignment operators
	mat3dd& operator += (const mat3dd& m);
	mat3dd& operator -= (const mat3dd& m);
	mat3dd& operator *= (const mat3dd& m);
	mat3dd& operator *= (double a);
	mat3dd& operator /= (double a);

	// matrix-vector multiplication
	vec3d operator * (const vec3d& r) const;

	// trace
	double tr() const;

	// determinant
	double det() const;

protected:
	double	d[3];	// the diagonal elements

	friend class mat3d;
	friend class mat3ds;
	friend class mat3da;
};

inline mat3dd operator * (double a, const mat3dd& d) { return d*a; }

//-----------------------------------------------------------------------------
//! This class describes a symmetric 3D matrix of doubles

class mat3ds
{
protected:
	// This enumeration can be used to remember the order
	// in which the components are stored.
	enum {
		XX = 0,
		XY = 1,
		YY = 2,
		XZ = 3,
		YZ = 4,
		ZZ = 5 };
public:
	// default constructor
	mat3ds(){}

	// constructors
	mat3ds(double xx, double yy, double zz, double xy, double yz, double xz);
	mat3ds(const mat3dd& d);

	// access operators
	double& operator () (int i, int j);
	const double& operator () (int i, int j) const;

	double& xx() { return m[XX]; }
	double& yy() { return m[YY]; }
	double& zz() { return m[ZZ]; }
	double& xy() { return m[XY]; }
	double& yz() { return m[YZ]; }
	double& xz() { return m[XZ]; }

	const double& xx() const { return m[XX]; }
	const double& yy() const { return m[YY]; }
	const double& zz() const { return m[ZZ]; }
	const double& xy() const { return m[XY]; }
	const double& yz() const { return m[YZ]; }
	const double& xz() const { return m[XZ]; }

	// arithmetic operators for mat3dd objects
	mat3ds operator + (const mat3dd& d) const;
	mat3ds operator - (const mat3dd& d) const;
	mat3ds operator * (const mat3dd& d) const;

	// arithmetic operators
	mat3ds operator + (const mat3ds& t) const;
	mat3ds operator - (const mat3ds& t) const;
	mat3ds operator * (const mat3ds& t) const;
	mat3ds operator * (double g) const;
	mat3ds operator / (double g) const;

	// arithmetic assignment operators
	mat3ds& operator += (const mat3ds& t);
	mat3ds& operator -= (const mat3ds& t);
	mat3ds& operator *= (const mat3ds& t);
	mat3ds& operator *= (double g);
	mat3ds& operator /= (double g);

	// arithmetic assignment operators for mat3dd
	mat3ds& operator += (const mat3dd& d);
	mat3ds& operator -= (const mat3dd& d);

	// matrix-vector multiplication
	vec3d operator * (const vec3d& r) const;

	// trace
	double tr() const;

	// determinant
	double det() const;

	// intialize to zero
	void zero();

	// deviator
	mat3ds dev() const;

	// determine eigen values and vectors
	void eigen(double d[3], vec3d r[3] = 0);

protected:
	double m[6];	// stores data in the order xx, xy, yy, xz, yz, zz

	friend class mat3dd;
	friend class mat3d;
};

inline mat3ds operator * (double a, const mat3ds& m) { return m*a; }

//-----------------------------------------------------------------------------
//! This class describes an anti-symmetric 3D matrix of doubles
// TODO: expand this class
class mat3da
{
public:
	// default constructor
	mat3da(){}

	// constructors
	mat3da(double xy, double yz, double xz);

	// calculates the antisymmetric matrix from a vector
	mat3da(const vec3d& a);

	// access operator
	double operator () (int i, int j);

protected:
	double	d[3];	// stores xy, yz, xz

	friend class mat3dd;
	friend class mat3ds;
	friend class mat3d;
};


//-----------------------------------------------------------------------------
//! This class describes a general 3D matrix of doubles
class mat3d
{
public:
	// default constructor
	mat3d() {}

	// constructors
	mat3d(double a00, double a01, double a02,
		  double a10, double a11, double a12,
		  double a20, double a21, double a22);

	mat3d(double m[3][3]);

	mat3d(const mat3dd& m);
	mat3d(const mat3ds& m);
	mat3d(const mat3da& m);

	// assignment operators
	mat3d& operator = (const mat3dd& m);
	mat3d& operator = (const mat3ds& m);

	// access operators
	double& operator () (int i, int j);
	const double& operator () (int i, int j) const;
	double* operator [] (int i);

	// arithmetic operators
	mat3d operator + (const mat3d& m) const;
	mat3d operator - (const mat3d& m) const;
	mat3d operator * (const mat3d& m) const;
	mat3d operator * (double a) const;
	mat3d operator / (double a) const;

	// arithmetic operators for mat3dd
	mat3d operator + (const mat3dd& m) const;
	mat3d operator - (const mat3dd& m) const;
	mat3d operator * (const mat3dd& m) const;

	// arithmetic operators for mat3ds
	mat3d operator + (const mat3ds& m) const;
	mat3d operator - (const mat3ds& m) const;
	mat3d operator * (const mat3ds& m) const;

	// arithmetic assignment operators
	mat3d& operator += (const mat3d& m);
	mat3d& operator -= (const mat3d& m);
	mat3d& operator *= (const mat3d& m);
	mat3d& operator *= (double a);
	mat3d& operator /= (double a);

	// arithmetic assignment operators for mat3dd
	mat3d& operator += (const mat3dd& m);
	mat3d& operator -= (const mat3dd& m);
	mat3d& operator *= (const mat3dd& m);

	// arithmetic assignment operators for mat3ds
	mat3d& operator += (const mat3ds& m);
	mat3d& operator -= (const mat3ds& m);
	mat3d& operator *= (const mat3ds& m);

	// matrix-vector muliplication
	vec3d operator * (const vec3d& r) const;

	// determinant
	double det() const;

	// trace
	double trace() const;

	// zero the matrix
	void zero();

	// make unit matrix
	void unit();

	// return a column vector from the matrix
	vec3d col(int j) const;

	// return the symmetric matrix 0.5*(A+A^T)
	mat3ds sym() const;

	// return the antisymmetric matrix 0.5*(A-A^T)
	mat3da skew() const;

	// calculates the inverse
	mat3d inverse() const;

	// calculates the transpose
	mat3d transpose() const;

	// calculates the transposed inverse
	mat3d transinv() const;

	// calculate the skew-symmetric matrix from a vector
	void skew(const vec3d& v);

protected:
	double d[3][3];	// matrix data

	friend class mat3dd;
	friend class mat3ds;
	friend class mat3da;
};

typedef double MATRIX3[3][3];

// outer product for vectors
inline mat3d operator & (const vec3d& a, const vec3d& b)
{
	return mat3d(a.x*b.x, a.x*b.y, a.x*b.z,
				 a.y*b.x, a.y*b.y, a.y*b.z,
				 a.z*b.x, a.z*b.y, a.z*b.z);
}

inline mat3ds dyad(const vec3d& a)
{
	return mat3ds(a.x*a.x, a.y*a.y, a.z*a.z, a.x*a.y, a.y*a.z, a.x*a.z);
}

// c_ij = a_i*b_j + a_j*b_i
inline mat3ds dyads(const vec3d& a, const vec3d& b)
{
	return mat3ds(2.0*a.x*b.x, 2.0*a.y*b.y, 2.0*a.z*b.z, a.x*b.y + a.y*b.x, a.y*b.z + a.z*b.y, a.x*b.z + a.y*b.z);
}

inline void matrix3_copy(MATRIX3& a, double b[3][3])
{
	a[0][0] = b[0][0]; a[0][1] = b[0][1]; a[0][2] = b[0][2];
	a[1][0] = b[1][0]; a[1][1] = b[1][1]; a[1][2] = b[1][2];
	a[2][0] = b[2][0]; a[2][1] = b[2][1]; a[2][2] = b[2][2];
}

// The following file contains the actual definition of the class functions
#include "mat3d.hpp"

#endif // !defined(AFX_MAT3D_H__35D54918_61A9_4E11_A4B2_6D3D33AB98FE__INCLUDED_)
