// FECoordSysMap.cpp: implementation of the FECoordSysMap class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FECoordSysMap.h"

//////////////////////////////////////////////////////////////////////
// FELocalMap
//////////////////////////////////////////////////////////////////////

FELocalMap::FELocalMap()
{
	m_n[0] = 0;
	m_n[1] = 1;
	m_n[2] = 3;
}

void FELocalMap::SetLocalNodes(int n1, int n2, int n3)
{
	m_n[0] = n1;
	m_n[1] = n2;
	m_n[2] = n3;
}

mat3d FELocalMap::LocalElementCoord(FEElement& el, int n)
{
	vec3d* r0 = el.r0();

	vec3d a, b, c, d;
	mat3d Q;

	a = r0[m_n[1]] - r0[m_n[0]];
	a.unit();

	if (m_n[2] != m_n[1])
	{
		d = r0[m_n[2]] - r0[m_n[0]];
	}
	else
	{
		d = vec3d(0,1,0);
		if (fabs(d*a) > 0.999) d = vec3d(1,0,0);
	}

	c = a^d;
	b = c^a;

	b.unit();
	c.unit();

	Q[0][0] = a.x; Q[0][1] = b.x; Q[0][2] = c.x;
	Q[1][0] = a.y; Q[1][1] = b.y; Q[1][2] = c.y;
	Q[2][0] = a.z; Q[2][1] = b.z; Q[2][2] = c.z;

	return Q;
}

//////////////////////////////////////////////////////////////////////
// FESphericalMap
//////////////////////////////////////////////////////////////////////

mat3d FESphericalMap::LocalElementCoord(FEElement& el, int n)
{
	vec3d* r0 = el.r0();
	double* H = el.H(n);
	vec3d a;
	for (int i=0; i<el.Nodes(); ++i) a += r0[i]*H[i];
	a -= m_c;
	a.unit();

	vec3d d = r0[1] - r0[0];
	d.unit();
	if (fabs(a*d) > .99) 
	{
		d = r0[2] - r0[1];
		d.unit();
	}

	vec3d c = a^d;
	vec3d b = c^a;

	a.unit();
	b.unit();
	c.unit();

	mat3d Q;
	Q[0][0] = a.x; Q[0][1] = b.x; Q[0][2] = c.x;
	Q[1][0] = a.y; Q[1][1] = b.y; Q[1][2] = c.y;
	Q[2][0] = a.z; Q[2][1] = b.z; Q[2][2] = c.z;

	return Q;
}

//////////////////////////////////////////////////////////////////////
// FEVectorMap
//////////////////////////////////////////////////////////////////////

mat3d FEVectorMap::LocalElementCoord(FEElement& el, int n)
{
	vec3d a = m_a;
	vec3d d = m_d;

	vec3d c = a^d;
	vec3d b = c^a;

	a.unit();
	b.unit();
	c.unit();

	mat3d Q;
	Q[0][0] = a.x; Q[0][1] = b.x; Q[0][2] = c.x;
	Q[1][0] = a.y; Q[1][1] = b.y; Q[1][2] = c.y;
	Q[2][0] = a.z; Q[2][1] = b.z; Q[2][2] = c.z;

	return Q;
}

//////////////////////////////////////////////////////////////////////
// FERandom2DMap
//////////////////////////////////////////////////////////////////////

inline float frand(int ix, int iy)
{
	static const float f = 1.f / 1073741824.0f;
	const int seed = 15477431;

	int n = seed + ix + 57*iy;
	n = (n << 13)^n;

	return (1.f - ((n*(n*n*15731 + 789221) + 1376132589) & 0x7fffffff) * f);
}

mat3d FERandom2DMap::LocalElementCoord(FEElement& el, int n)
{
	double* H = el.H(n);
	vec3d* r0 = el.r0();
	vec3d a;
	for (int i=0; i<el.Nodes(); ++i) a += r0[i]*H[i];

	double w = 2.0*3.1415926*frand((int)(100*a.x), (int)(100*a.y));

	a.x = cos(w);
	a.y = sin(w);
	a.z = 0;

	vec3d d = r0[1] - r0[0];
	d.unit();
	if (fabs(a*d) > .99) 
	{
		d = r0[2] - r0[1];
		d.unit();
	}

	vec3d c = a^d;
	vec3d b = c^a;

	a.unit();
	b.unit();
	c.unit();

	mat3d Q;
	Q[0][0] = a.x; Q[0][1] = b.x; Q[0][2] = c.x;
	Q[1][0] = a.y; Q[1][1] = b.y; Q[1][2] = c.y;
	Q[2][0] = a.z; Q[2][1] = b.z; Q[2][2] = c.z;

	return Q;
}
