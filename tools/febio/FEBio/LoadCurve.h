// LoadCurve.h: interface for the LoadCurve class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_LOADCURVE_H__E2D23AD8_1961_49E7_BE3A_3FC885CED53D__INCLUDED_)
#define AFX_LOADCURVE_H__E2D23AD8_1961_49E7_BE3A_3FC885CED53D__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "FECore/vector.h"
#include "Archive.h"

//-----------------------------------------------------------------------------
//! This class implements the concept of a loadcurve.

//! A loadcurve is basically a discretized function of time versus load,
//! where load can be interpreted differently in different contexts.
//! The loadcurve stores the (time,load)-pairs for fixed points, which
//! are input from the input file.
//! In between timesteps, the loadcurve class interpolates the (time,load)
//! data pairs according to the interpolation function.

class FELoadCurve
{
public:
	//! Load point structure
	struct LOADPOINT
	{
		double time;
		double value;
	};

public:
	//! Interpolation functions
	enum INTFUNC { STEP=0, LINEAR=1 };

public:
	//! default constructor
	FELoadCurve() : m_fnc(LINEAR) {}

	//! destructor
	virtual ~FELoadCurve() {}

	//! creates the data points of the load curve
	void Create(int n);

	//! adds a point to the loadcurve
	void Add(double time, double value);

	//! Clears the loadcurve data
	void Clear() { m_lp.clear(); }

	//! set the time and data value of point i of the load curve
	void SetPoint(int i, double time, double val);

	//! returns the value of the load curve at time
	double Value(double time);
	double Value() { return m_value; }

	//! evaluates the loadcurve at time
	void Evaluate(double time)
	{
		m_value = Value(time);
	}

	//! Set the type of interpolation
	void SetInterpolation(INTFUNC fnc) { m_fnc = fnc; }

	//! returns point i
	LOADPOINT& LoadPoint(int i) { return m_lp[i]; }

	//! finds closest load point
	int FindPoint(double t);

	//! return nr of points
	int Points() { return m_lp.size(); }

	//! see if there is a point at time t
	bool HasPoint(double t);

	//! Serialize data to archive
	void Serialize(Archive& ar);

protected:
	vector<LOADPOINT>	m_lp;	//!< load time values

	double			m_value;	//!< last calculated value

	INTFUNC				m_fnc;	//!< interpolation function
};

typedef FELoadCurve::LOADPOINT LOADPOINT;

#endif // !defined(AFX_LOADCURVE_H__E2D23AD8_1961_49E7_BE3A_3FC885CED53D__INCLUDED_)
