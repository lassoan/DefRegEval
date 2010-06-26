// FEPoroElastic.cpp: implementation of the FEPoroElastic class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FEPoroElastic.h"


//-----------------------------------------------------------------------------
//! FEPoroElastic constructor

FEPoroElastic::FEPoroElastic()
{
}

//-----------------------------------------------------------------------------
//! The stress of a poro-elastic material is the sum of the fluid pressure
//! and the elastic stress. Note that this function is declared in the base class
//! so you do not have to reimplement it in a derived class, unless additional
//! pressure terms are required.

mat3ds FEPoroElastic::Stress(FEMaterialPoint& mp)
{
	FEPoroElasticMaterialPoint& pt = *mp.ExtractData<FEPoroElasticMaterialPoint>();

	// calculate solid material stress
	mat3ds s = m_pBase->Stress(mp);

	// add fluid pressure
	s.xx() -= pt.m_p;
	s.yy() -= pt.m_p;
	s.zz() -= pt.m_p;

	return s;
}

//-----------------------------------------------------------------------------
//! The tangent is the sum of the elastic tangent plus the fluid tangent. Note
//! that this function is declared in the base class, so you don't have to 
//! reimplement it unless additional tangent components are required.

tens4ds FEPoroElastic::Tangent(FEMaterialPoint& mp)
{
	FEPoroElasticMaterialPoint& pt = *mp.ExtractData<FEPoroElasticMaterialPoint>();

	// call solid tangent routine
	tens4ds c = m_pBase->Tangent(mp);

	// fluid pressure
	double p = pt.m_p;

	// adjust tangent for pressures
	double D[6][6] = {0};
	c.extract(D);

	D[0][0] -= -p;
	D[1][1] -= -p;
	D[2][2] -= -p;

	D[0][1] -= p; D[1][0] -= p;
	D[1][2] -= p; D[2][1] -= p;
	D[0][2] -= p; D[2][0] -= p;

	D[3][3] -= -p;
	D[4][4] -= -p;
	D[5][5] -= -p;

	return tens4ds(D);
}
