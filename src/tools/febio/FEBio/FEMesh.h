// FEMesh.h: interface for the FEMesh class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FEMESH_H__81ABA97F_AD5F_4F1D_8EE9_95B67EBA448E__INCLUDED_)
#define AFX_FEMESH_H__81ABA97F_AD5F_4F1D_8EE9_95B67EBA448E__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "FEElement.h"
#include "Archive.h"
#include "FENodeElemList.h"

///////////////////////////////////////////////////////////////////////////////
// STRUCT: FE_BOUNDING_BOX
//  This class stores the coordinates of a bounding box
//

struct FE_BOUNDING_BOX
{
	vec3d	r0, r1; // coordinates of opposite corners

	// center of box
	vec3d center() { return (r0+r1)*0.5; }

	// dimensions of box
	double width () { return (r1.x-r0.x); }
	double height() { return (r1.y-r0.y); }
	double depth () { return (r1.z-r0.z); }

	// max dimension
	double radius() 
	{ 
		double w = width();
		double h = height();
		double d = depth();

		if ((w>=d) && (w>=h)) return w;
		if ((h>=w) && (h>=d)) return h;
		
		return d;
	}

	void operator += (vec3d r)
	{
		if (r.x < r0.x) r0.x = r.x;
		if (r.y < r0.y) r0.y = r.y;
		if (r.z < r0.z) r0.z = r.z;
		if (r.x > r1.x) r1.x = r.x;
		if (r.y > r1.y) r1.y = r.y;
		if (r.z > r1.z) r1.z = r.z;
	}

	void inflate(double dx, double dy, double dz)
	{
		r0.x -= dx; r1.x += dx;
		r0.y -= dy; r1.y += dy;
		r0.z -= dz; r1.z += dz;
	}

	// check whether a point is inside or not
	bool IsInside(vec3d r) 
	{ 
		return ((r.x>=r0.x) && (r.y>=r0.y) && (r.z>=r0.z) && (r.x<=r1.x) && (r.y<=r1.y) && (r.z<=r1.z));
	}
};

//-----------------------------------------------------------------------------
//! This class defines a finite element node

//! It stores nodal positions and nodal equations numbers and more.

class FENode
{
public:
	// geometry data
	vec3d	m_r0;	//!< initial position
	vec3d	m_v0;	//!< initial velocity

	vec3d	m_rt;	//!< current position
	vec3d	m_vt;	//!< nodal velocity
	vec3d	m_at;	//!< nodal acceleration

	vec3d	m_rp;	//!< position of node at previous time step
	vec3d	m_vp;	//!< previous velocity
	vec3d	m_ap;	//!< previous acceleration


	// shell data
	vec3d	m_D0;	//!< initial director
	vec3d	m_Dt;	//!< current director

	// poroelasticity-data
	double	m_pt;	//!< current pressure

	// heat-conduction data
	double	m_T;	//!< temperature

	// rigid body data
	int		m_rid;	//!< rigid body number

	bool	m_bshell;	//!< does this node belong to a non-rigid shell element?

	int		m_ID[MAX_NDOFS];	//!< nodal equation numbers
	int		m_BC[MAX_NDOFS];	//!< boundary condition
};

class FEMesh;

//-----------------------------------------------------------------------------
//! Defines a node set

class FENodeSet
{
public:
	FENodeSet(FEMesh* pm) : m_pmesh(pm), m_nID(-1) { m_szname[0] = 0; }

	void create(int n)
	{
		assert(n);
		m_Node.setsize(n);
	}

	int size() { return m_Node.size(); }

	int& operator [] (int i) { return m_Node[i]; }

	void SetID(int n) { m_nID = n; }
	int GetID() { return m_nID; }

	void SetName(const char* sz) { strcpy(m_szname, sz); }
	const char* GetName() { return m_szname; }

protected:
	FEMesh*	m_pmesh;
	vector<int>	m_Node;		//!< list of nodes

	int	m_nID;
	char	m_szname[256];
};

//-----------------------------------------------------------------------------
//! Defines a finite element mesh

//! All the geometry data is stored in this class. 

class FEMesh  
{
public:
	//! constructor
	FEMesh();

	//! copy constructor
	FEMesh(FEMesh& m);

	//! destructor
	virtual ~FEMesh();

	//! assignment operator
	FEMesh& operator = (FEMesh& m);

	//! allocate storage for mesh data
	void Create(int nodes, int elems, int shells, int ntruss);

	//! return number of nodes
	int Nodes() { return m_Node.size(); }

	//! count the number of solid elements
	int SolidElements() { return m_Elem.size(); }

	//! count the number of shell elements
	int ShellElements() { return m_Shell.size(); }

	//! counte the number of truss elements
	int TrussElements() { return m_Truss.size(); }

	//! return reference to a node
	FENode& Node(int i) { return m_Node[i]; }

	//! return reference to a solid element
	FESolidElement& SolidElement(int i) { return m_Elem[i]; }

	//! return reference to a shell element
	FEShellElement& ShellElement(int i) { return m_Shell[i]; }

	//! return a reference to a  truss element
	FETrussElement& TrussElement(int i) { return m_Truss[i]; }

	//! update bounding box
	void UpdateBox();

	//! retrieve the bounding box
	FE_BOUNDING_BOX& GetBoundingBox() { return m_box; }

	//! remove isolated vertices
	int RemoveIsolatedVertices();

	//! Unpack solid element data
	void UnpackElement(FESolidElement& el, unsigned int nflag = FE_UNPACK_ALL);

	//! Unpack shell element data
	void UnpackElement(FEShellElement& el, unsigned int nflag = FE_UNPACK_ALL);

	//! Unpack surface element data
	void UnpackElement(FESurfaceElement& el, unsigned int nflag = FE_UNPACK_ALL);

	//! Unpack truss element data
	void UnpackElement(FETrussElement& el, unsigned int flag = FE_UNPACK_ALL);

	//! Reset the mesh data
	void Reset();

	//! initialize the mesh data
	bool Init();

	//! Calculates an elements volume
	double ElementVolume(FEElement& el);

	//! assign a material ID to the entire mesh
	void SetMatID(int n);

	//! adds a node set to the mesh
	void AddNodeSet(FENodeSet* pns) { m_NodeSet.add(pns); }

	//! Find a nodeset by ID
	FENodeSet* FindNodeSet(int nid);

	//! Find a nodeset by name
	FENodeSet* FindNodeSet(const char* szname);

	//! serialize data to or from a binary archive
	void Serialize(Archive& ar);

	//! Get the face nodes from a given element
	int GetFace(FEElement& el, int n, int nf[4]);

	//! return the nr of faces an element has
	int Faces(FEElement& el);

	//! Finds an element from a given ID
	FEElement* FindElementFromID(int nid);

	FENodeElemList& NodeElementList()
	{
		if (m_NEL.Size() != m_Node.size()) m_NEL.Create(*this);
		return m_NEL;
	}

protected:
	vector<FENode>			m_Node;		//!< FE nodes array

	vector<FESolidElement>	m_Elem;		//!< FE solid element array
	vector<FEShellElement>	m_Shell;	//!< FE shell element array
	vector<FETrussElement>	m_Truss;	//!< FE truss element array

	ptr_vector<FENodeSet>	m_NodeSet;	//!< node sets

	FE_BOUNDING_BOX		m_box;	//!< bounding box

	FENodeElemList	m_NEL;
};

#endif // !defined(AFX_FEMESH_H__81ABA97F_AD5F_4F1D_8EE9_95B67EBA448E__INCLUDED_)
