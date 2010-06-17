// FEMooneyRivlin.h: interface for the FEMooneyRivlin class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FEMOONEYRIVLIN_H__D75E80FD_0E25_4A12_9539_044C9DC4CB41__INCLUDED_)
#define AFX_FEMOONEYRIVLIN_H__D75E80FD_0E25_4A12_9539_044C9DC4CB41__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "FEMaterial.h"

//-----------------------------------------------------------------------------
//! Mooney-Rivlin material

class FEMooneyRivlin : public FEIncompressibleMaterial
{
public:
	FEMooneyRivlin() {}

public:
	double	c1;	//!< Mooney-Rivlin coefficient C1
	double	c2;	//!< Mooney-Rivlin coefficient C2

public:
	//! calculate stress at material point
	virtual mat3ds Stress(FEMaterialPoint& pt);

	//! calculate tangent stiffness at material point
	virtual tens4ds Tangent(FEMaterialPoint& pt);

	//! data initialization
	void Init();

	// declare as registered
	DECLARE_REGISTERED(FEMooneyRivlin);

	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};

#endif // !defined(AFX_FEMOONEYRIVLIN_H__D75E80FD_0E25_4A12_9539_044C9DC4CB41__INCLUDED_)
