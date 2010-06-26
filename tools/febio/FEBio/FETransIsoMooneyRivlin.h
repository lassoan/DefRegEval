// FETransIsoMooneyRivlin.h: interface for the FETransIsoMooneyRivlin class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FETRANSISOMOONEYRIVLIN_H__E918D89B_4CCD_44B9_9731_19CEC4EDF406__INCLUDED_)
#define AFX_FETRANSISOMOONEYRIVLIN_H__E918D89B_4CCD_44B9_9731_19CEC4EDF406__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "FEMaterial.h"

//-----------------------------------------------------------------------------
//! Transversely Isotropic Mooney-Rivlin material

//! This material has an isotopric Mooney-Rivlin basis and single preferred
//! fiber direction.

class FETransIsoMooneyRivlin: public FETransverselyIsotropic
{
public:
	FETransIsoMooneyRivlin () {}

public:
	double	c1;	//!< Mooney-Rivlin coefficient C1
	double	c2;	//!< Mooney-Rivlin coefficient C2

public:
	//! calculate stress at material point
	virtual mat3ds Stress(FEMaterialPoint& pt);

	//! calculate tangent stiffness at material point
	virtual tens4ds Tangent(FEMaterialPoint& pt);

	// declare as registered
	DECLARE_REGISTERED(FETransIsoMooneyRivlin);

	// declare parameter list
	DECLARE_PARAMETER_LIST();
};


#endif // !defined(AFX_FETRANSISOMOONEYRIVLIN_H__E918D89B_4CCD_44B9_9731_19CEC4EDF406__INCLUDED_)
