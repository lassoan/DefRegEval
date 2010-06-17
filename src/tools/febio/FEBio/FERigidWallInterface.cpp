// FERigidWallInterface.cpp: implementation of the FERigidWallInterface class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FERigidWallInterface.h"
#include "fem.h"
#include "FESolidSolver.h"
#include "FENNQuery.h"
#include "log.h"

///////////////////////////////////////////////////////////////////////////////
// FERigidWallSurface
///////////////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
//! Finds the (master) element that contains the projection of a (slave) node

FEElement* FERigidWallSurface::FindMasterSegment(vec3d& x, vec3d& q, vec2d& r, bool& binit_nq, double tol)
{
	// get the mesh
	FEMesh& mesh = *m_pmesh;

	// see if we need to initialize the NQ structure
	if (binit_nq) m_NQ.Init();
	binit_nq = false;

	// let's find the closest master node
	int mn = m_NQ.Find(x);

	// mn is a local index, so get the global node number too
	int m = node[mn];

	// get the nodal position
	vec3d r0 = mesh.Node(m).m_rt;

	// now that we found the closest master node, lets see if we can find 
	// the best master element
	int N;

	// loop over all master elements that contain the node mn
	int nval = m_NEL.Valence(mn);
	FEElement** pe = m_NEL.ElementList(mn);
	for (int j=0; j<nval; ++j)
	{
		// get the master element
		FESurfaceElement& el = dynamic_cast<FESurfaceElement&> (*pe[j]);
		N = el.Nodes();

		// project the node on the element
		r[0] = 0;
		r[1] = 0;
		q = ProjectToSurface(el, x, r[0], r[1]);
		if (IsInsideElement(el, r[0], r[1], tol)) return pe[j];
	}

	// we did not find a master surface
	return 0;
}

//-----------------------------------------------------------------------------
//! Creates a surface for use with a sliding interface. All surface data
//! structures are allocated.
//! Note that it is assumed that the element array is already created
//! and initialized.

void FERigidWallSurface::Init()
{
	int i, j, n;

	// always intialize base class first!
	FESurface::Init();

	// get the number of nodes
	int nn = Nodes();

	// allocate other surface data
	gap.create(nn);		// gap funtion
	nu.create(nn);		// node normal 
	pme.create(nn);		// penetrated master element
	rs.create(nn);		// natural coords of projected slave node on master element
	rsp.create(nn);
	Lm.create(nn);
	M.create(nn);
	Lt.create(nn);
	off.create(nn);
	eps.create(nn);

	// set initial values
	gap.zero();
	pme.set(0);
	Lm.zero();
	off.zero();
	eps.set(1.0);

	// we calculate the gap offset values
	// This value is used to take the shell thickness into account
	// note that we force rigid shells to have zero thickness
	FEMesh& m = *m_pmesh;
	vector<double> tag(m.Nodes());
	tag.zero();
	for (i=0; i<m.ShellElements(); ++i)
	{
		FEShellElement& el = m.ShellElement(i);
		n = el.Nodes();
		for (j=0; j<n; ++j) tag[el.m_node[j]] = 0.5*el.m_h0[j];
	}
	for (i=0; i<nn; ++i) off[i] = tag[node[i]];
}

//-----------------------------------------------------------------------------
//! 

vec3d FERigidWallSurface::traction(int inode)
{
	vec3d t(0,0,0);
	FEElement* pe = pme[inode];
	if (pe)
	{
		FESurfaceElement& el = dynamic_cast<FESurfaceElement&>(*pe);
		double Tn = Lm[inode];
		double T1 = Lt[inode][0];
		double T2 = Lt[inode][1];
		double r = rs[inode][0];
		double s = rs[inode][1];

		vec3d tn = nu[inode]*Tn, tt;
		vec3d e[2];
		ContraBaseVectors0(el, r, s, e);
		tt = e[0]*T1 + e[1]*T2;
		t = tn + tt;
	}

	return t;
}

//-----------------------------------------------------------------------------

void FERigidWallSurface::UpdateNormals()
{
	int i, j, jp1, jm1;
	int N = Nodes();
	int NE = Elements();
	for (i=0; i<N; ++i) nu[i] = vec3d(0,0,0);
	vec3d y[4], e1, e2;

	for (i=0; i<NE; ++i)
	{
		FESurfaceElement& el = Element(i);
		int ne = el.Nodes();
		for (j=0; j<ne; ++j) y[j] = Node(el.m_lnode[j]).m_rt;

		for (j=0; j<ne; ++j)
		{
			jp1 = (j+1)%ne;
			jm1 = (j+ne-1)%ne;

			e1 = y[jp1] - y[j];
			e2 = y[jm1] - y[j];

			nu[el.m_lnode[j]] -= e1 ^ e2;						
		}
	}

	for (i=0; i<N; ++i) nu[i].unit();
}


///////////////////////////////////////////////////////////////////////////////
// FERigidWallInterface
///////////////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
//! constructor
FERigidWallInterface::FERigidWallInterface(FEM* pfem) : FEContactInterface(pfem), m_ss(&pfem->m_mesh)
{
	static int count = 1;
	m_ntype = FE_CONTACT_RIGIDWALL;

	m_mp = 0;

	m_nplc = -1;
	m_pplc = 0;
	m_nID = count++;
};


//-----------------------------------------------------------------------------
//! Initializes the rigid wall interface data

void FERigidWallInterface::Init()
{
	// create the surface
	m_ss.Init();

	// set the penalty load curve
	if (m_nplc >= 0) m_pplc = m_pfem->GetLoadCurve(m_nplc);

	// initialize rigid surface
	m_mp->Init();

	// project slave surface onto master surface
	ProjectSurface(m_ss);
}

//-----------------------------------------------------------------------------
//!  Projects the slave surface onto the master plane

void FERigidWallInterface::ProjectSurface(FERigidWallSurface& ss)
{
	vec3d r, q;

	// surface normal
	vec3d np;

	// loop over all slave nodes
	for (int i=0; i<m_ss.Nodes(); ++i)
	{
		// get the nodal position
		r = m_ss.Node(i).m_rt;

		// project this node onto the plane
		q = m_mp->Project(r);

		// get the local surface normal
		np = m_mp->Normal(q);

		// the slave normal is set to the master element normal
		m_ss.nu[i] = np;
	
		// calculate initial gap
		m_ss.gap[i] = -(np*(r - q)) + m_ss.off[i];
	}
}

//-----------------------------------------------------------------------------
//!  Updates rigid wall data

void FERigidWallInterface::Update()
{
	// project slave surface onto master surface
	ProjectSurface(m_ss);
}

//-----------------------------------------------------------------------------

void FERigidWallInterface::ContactForces(vector<double>& F)
{
	int j, k, m, n;
	int nseln;

	double *Gr, *Gs;

	// jacobian
	double detJ;

	vec3d dxr, dxs;
	vec3d* rt, *r0;
	double* w;

	// normal force
	double tn;

	// element contact force vector
	// note that this assumes that the element to which the integration node
	// connects is a four noded quadrilateral
	vector<double> fe(3);

	// the lm array for this force vector
	vector<int> lm(3);

	// the en array
	vector<int> en(1);

	vector<int> sLM;

	FEMesh& mesh = m_pfem->m_mesh;

	FESolidSolver* psolver = dynamic_cast<FESolidSolver*>(m_pfem->m_pStep->m_psolver);

	// penalty value
	double pen = Penalty(), eps;
	
	// loop over all slave facets
	int ne = m_ss.Elements();
	for (j=0; j<ne; ++j)
	{
		// get the slave element
		FESurfaceElement& sel = m_ss.Element(j);
		mesh.UnpackElement(sel);

		sLM = sel.LM();

		nseln = sel.Nodes();

		rt = sel.rt();
		r0 = sel.r0();
		w = sel.GaussWeights();

		// loop over slave element nodes (which are the integration points as well)
		for (n=0; n<nseln; ++n)
		{
			Gr = sel.Gr(n);
			Gs = sel.Gs(n);

			m = sel.m_lnode[n];

			// see if this node's constraint is active
			// that is, if it has a master element associated with it
//			if (ss.Lm[m] >= 0)
			{
				// calculate jacobian
				dxr = dxs = vec3d(0,0,0);
				for (k=0; k<nseln; ++k)
				{
					dxr.x += Gr[k]*r0[k].x;
					dxr.y += Gr[k]*r0[k].y;
					dxr.z += Gr[k]*r0[k].z;

					dxs.x += Gs[k]*r0[k].x;
					dxs.y += Gs[k]*r0[k].y;
					dxs.z += Gs[k]*r0[k].z;
				}

				detJ = (dxr ^ dxs).norm();

				// get slave node normal force
				eps = pen*m_ss.eps[m];
				tn = m_ss.Lm[m] + eps*m_ss.gap[m];
				tn = MBRACKET(tn);

				// get the slave node normal
				vec3d& nu = m_ss.nu[m];

				// calculate force vector
				fe[0] = detJ*w[n]*tn*nu.x;
				fe[1] = detJ*w[n]*tn*nu.y;
				fe[2] = detJ*w[n]*tn*nu.z;
	
				// fill the lm array
				lm[0] = sLM[n*3  ];
				lm[1] = sLM[n*3+1];
				lm[2] = sLM[n*3+2];

				// fill the en array
				en[0] = sel.m_node[n];

				// assemble into global force vector
				psolver->AssembleResidual(en, lm, fe, F);
			}
		}
	}

}

//-----------------------------------------------------------------------------
//! This function calculates the stiffness contribution for the rigid wall
//! interface.
//! \todo I think there are a couple of stiffness terms missing in this formulation

void FERigidWallInterface::ContactStiffness()
{
	int j, k, l, n, m;
	int nseln, ndof;

	matrix ke;

	vector<int> lm(3);
	vector<int> en(1);

	double *Gr, *Gs, *w;
	vec3d *rt, *r0;

	double detJ, tn;
	vec3d dxr, dxs;
	double N[3];

	double gap, Lm;

	vector<int> sLM;

	FEMesh& mesh = m_pfem->m_mesh;

	FESolidSolver* psolver = dynamic_cast<FESolidSolver*>(m_pfem->m_pStep->m_psolver);

	// penalty value
	double pen = Penalty(), eps;

	// loop over all slave elements
	int ne = m_ss.Elements();
	for (j=0; j<ne; ++j)
	{
		FESurfaceElement& se = m_ss.Element(j);
		mesh.UnpackElement(se);

		sLM = se.LM();

		nseln = se.Nodes();

		r0 = se.r0();
		rt = se.rt();
		w = se.GaussWeights();

		// loop over all integration points (that is nodes)
		for (n=0; n<nseln; ++n)
		{
			Gr = se.Gr(n);
			Gs = se.Gs(n);

			m = se.m_lnode[n];

			// see if this node's constraint is active
			// that is, if it has a master element associated with it
//			if (ss.Lm[m] >= 0)
			{
				// calculate jacobian
				dxr = dxs = vec3d(0,0,0);
				for (k=0; k<nseln; ++k)
				{
					dxr.x += Gr[k]*r0[k].x;
					dxr.y += Gr[k]*r0[k].y;
					dxr.z += Gr[k]*r0[k].z;

					dxs.x += Gs[k]*r0[k].x;
					dxs.y += Gs[k]*r0[k].y;
					dxs.z += Gs[k]*r0[k].z;
				}

				detJ = (dxr ^ dxs).norm();

				// slave gap
				gap = m_ss.gap[m];

				// lagrange multiplier
				Lm = m_ss.Lm[m];

				// get slave node normal force
				eps = pen*m_ss.eps[m];

				tn = m_ss.Lm[m] + eps*m_ss.gap[m];
				tn = MBRACKET(tn);

				// get the slave node normal
				vec3d& nu = m_ss.nu[m];

				// set up the N vector
				N[0] = nu.x;
				N[1] = nu.y;
				N[2] = nu.z;

				ndof = 3;

				// fill stiffness matrix
				// TODO: I don't think this is correct, since
				// if the rigid wall is not a plance, some terms
				// are probably missing
				ke.Create(ndof, ndof);
				for (k=0; k<ndof; ++k)
					for (l=0; l<ndof; ++l)
						ke[k][l] = w[n]*detJ*eps*HEAVYSIDE(Lm+eps*gap)*N[k]*N[l];

				// create lm array
				lm[0] = sLM[n*3  ];
				lm[1] = sLM[n*3+1];
				lm[2] = sLM[n*3+2];

				// create the en array
				en[0] = se.m_node[n];
						
				// assemble stiffness matrix
				psolver->AssembleStiffness(en, lm, ke);
			}
		}
	}
}

//-----------------------------------------------------------------------------

bool FERigidWallInterface::Augment(int naug)
{
	int i;
	double Lm;
	bool bconv = true;

	// penalty value
	double pen = Penalty(), eps;

	// calculate initial norms
	double normL0 = 0;
	for (i=0; i<m_ss.Nodes(); ++i) normL0 += m_ss.Lm[i]*m_ss.Lm[i];
	normL0 = sqrt(normL0);

	// update Lagrange multipliers and calculate current norms
	double normL1 = 0;
	double normgc = 0;
	int N = 0;
	for (i=0; i<m_ss.Nodes(); ++i)
	{
		// update Lagrange multipliers
		eps = pen*m_ss.eps[i];

		Lm = m_ss.Lm[i] + eps*m_ss.gap[i];
		Lm = MBRACKET(Lm);
		normL1 += Lm*Lm;
		if (m_ss.gap[i] > 0)
		{
			normgc += m_ss.gap[i]*m_ss.gap[i];
			++N;
		}
	}	
	if (N==0) N = 1;

	normL1 = sqrt(normL1);
	normgc = sqrt(normgc / N);

	// get the logfile
	Logfile& log = GetLogfile();

	// check convergence of constraints
	log.printf(" rigid wall interface # %d\n", m_nID);
	log.printf("                        CURRENT        REQUIRED\n");
	double pctn = 0;
	if (fabs(normL1) > 1e-10) pctn = fabs((normL1 - normL0)/normL1);
	log.printf("    normal force : %15le %15le\n", pctn, m_atol);
	log.printf("    gap function : %15le       ***\n", normgc);
		
	if (pctn >= m_atol)
	{
		bconv = false;
		for (i=0; i<m_ss.Nodes(); ++i)
		{
			// update Lagrange multipliers
			eps = pen*m_ss.eps[i];

			Lm = m_ss.Lm[i] + eps*m_ss.gap[i];
			m_ss.Lm[i] = MBRACKET(Lm);
		}	
	}

	return bconv;
}

//-----------------------------------------------------------------------------

void FERigidWallInterface::Serialize(Archive &ar)
{
	int k, n, mat;
	if (ar.IsSaving())
	{
		ar << m_eps;
		ar << m_atol;
		ar << m_nplc;

		FERigidWallSurface& s = m_ss;

		int ne = s.Elements();
		ar << ne;

		for (k=0; k<ne; ++k)
		{
			FESurfaceElement& el = s.Element(k);
			ar << el.Type();
			ar << el.GetMatID() << el.m_nID << el.m_nrigid;
			ar << el.m_node;
			ar << el.m_lnode;
		}

		ar << s.gap;
		ar << s.nu;
		ar << s.rs;
		ar << s.Lm;
		ar << s.off;
		
		// plane data
		if (dynamic_cast<FEPlane*>(m_mp))
		{
			FEPlane* pp = dynamic_cast<FEPlane*>(m_mp);
			ar << FE_RIGID_PLANE;
			ar << pp->m_nplc;
			double* a = pp->GetEquation();
			ar << a[0] << a[1] << a[2] << a[3];
		}
		else if (dynamic_cast<FERigidSphere*>(m_mp))
		{
			FERigidSphere* ps = dynamic_cast<FERigidSphere*>(m_mp);
			ar << FE_RIGID_SPHERE;
			ar << ps->m_rc;
			ar << ps->m_R;
			ar << ps->m_nplc[0] << ps->m_nplc[1] << ps->m_nplc[2];
		}
	}
	else
	{
		ar >> m_eps;
		ar >> m_atol;
		ar >> m_nplc;

		if (m_nplc >= 0) m_pplc = m_pfem->GetLoadCurve(m_nplc);

		FERigidWallSurface& s = m_ss;

		int ne=0;
		ar >> ne;
		s.Create(ne);

		for (k=0; k<ne; ++k)
		{
			FESurfaceElement& el = s.Element(k);

			ar >> n;
			el.SetType(n);

			ar >> mat >> el.m_nID >> el.m_nrigid;
			ar >> el.m_node;
			ar >> el.m_lnode;

			el.SetMatID(mat);
		}

		// initialize surface
		s.Init();

		ar >> s.gap;
		ar >> s.nu;
		ar >> s.rs;
		ar >> s.Lm;
		ar >> s.off;

		// plane data
		int ntype;
		ar >> ntype;
		switch (ntype)
		{
		case FE_RIGID_PLANE:
			{
				SetMasterSurface(new FEPlane(m_pfem));
				FEPlane& pl = dynamic_cast<FEPlane&>(*m_mp);
				ar >> pl.m_nplc;
				if (pl.m_nplc >= 0) pl.m_pplc = m_pfem->GetLoadCurve(pl.m_nplc);
				double* a = pl.GetEquation();
				ar >> a[0] >> a[1] >> a[2] >> a[4];
			}
			break;
		case FE_RIGID_SPHERE:
			{
				SetMasterSurface(new FERigidSphere(m_pfem));
				FERigidSphere& s = dynamic_cast<FERigidSphere&>(*m_mp);
				ar >> s.m_rc;
				ar >> s.m_R;
				ar >> s.m_nplc[0] >> s.m_nplc[1] >> s.m_nplc[2];
			}
			break;
		default:
			assert(false);
		}
	}
}
