#include "stdafx.h"
#include "FEPeriodicBoundary.h"
#include "Archive.h"
#include "fem.h"
#include "FESolidSolver.h"
#include "log.h"

//-----------------------------------------------------------------------------
//! Creates a surface for use with a sliding interface. All surface data
//! structures are allocated.
//! Note that it is assumed that the element array is already created
//! and initialized.

void FEPeriodicSurface::Init()
{
	// always intialize base class first!
	FESurface::Init();

	// get the number of nodes
	int nn = Nodes();

	// allocate other surface data
	m_gap.create(nn);		// gap funtion
	m_pme.create(nn);		// penetrated master element
	m_rs.create(nn);		// natural coords of projected slave node on master element
	m_Lm.create(nn);		// Lagrangian multipliers

	// set initial values
	m_gap.zero();
	m_pme.set(0);
	m_Lm.zero();
}

//-----------------------------------------------------------------------------
//! Calculate the center of mass for this surface
//!
vec3d FEPeriodicSurface::CenterOfMass()
{
	vec3d c(0,0,0);
	int N = Nodes();
	for (int i=0; i<N; ++i) c += Node(i).m_r0;
	if (N != 0) c /= (double) N;
	return c;
}

//-----------------------------------------------------------------------------
// FEPeriodicBoundary
//-----------------------------------------------------------------------------

FEPeriodicBoundary::FEPeriodicBoundary(FEM* pfem) : FEContactInterface(pfem), m_ss(&pfem->m_mesh), m_ms(&pfem->m_mesh)
{
	static int count = 1;
	m_ntype = FE_PERIODIC_BOUNDARY;

	m_stol = 0.01;
	m_atol = 0;
	m_eps = 0;
	m_npass = 1;

	m_nID = count++;
}

//-----------------------------------------------------------------------------
void FEPeriodicBoundary::Init()
{
	// create the surfaces
	m_ss.Init();
	m_ms.Init();

	// project slave surface onto master surface
	ProjectSurface(m_ss, m_ms, false);
	ProjectSurface(m_ms, m_ss, false);
}

//-----------------------------------------------------------------------------
//! project surface

void FEPeriodicBoundary::ProjectSurface(FEPeriodicSurface& ss, FEPeriodicSurface& ms, bool bmove)
{
	int i, nm;
	double rs[2];

	// get the slave's center of mass
	vec3d cs = ss.CenterOfMass();

	// get the master's center of mass
	vec3d cm = ms.CenterOfMass();

	// get the relative distance
	vec3d cr = cm - cs;

	// unit vector in direction of cr
	// this will serve as the projection distance
	vec3d cn(cr); cn.unit();

	// loop over all slave nodes
	for (i=0; i<ss.Nodes(); ++i)
	{
		FENode& node = ss.Node(i);

		// get the nodal position
		vec3d r0 = node.m_r0;

		// find the intersection with the master surface
		ss.m_pme[i] = ms.FindIntersection(r0, cn, rs, m_stol, &nm);
		assert(ss.m_pme[i]);

		ss.m_rs[i][0] = rs[0];
		ss.m_rs[i][1] = rs[1];
	}
}

//-----------------------------------------------------------------------------
void FEPeriodicBoundary::Update()
{
	int i, j, ne;
	FESurfaceElement* pme;

	FEMesh& mesh = *m_ss.GetMesh();

	vec3d us, um;
	vec3d umi[4];

	// update gap functions
	for (int np=0; np<m_npass; ++np)
	{
		FEPeriodicSurface& ss = (np == 0? m_ss : m_ms);
		FEPeriodicSurface& ms = (np == 0? m_ms : m_ss);

		int N = ss.Nodes();

		for (i=0; i<N; ++i)
		{
			// calculate the slave displacement
			FENode& node = ss.Node(i);
			us = node.m_rt - node.m_r0;

			// get the master element
			pme = ss.m_pme[i];

			// calculate the master displacement
			ne = pme->Nodes();
			for (j=0; j<ne; ++j)
			{
				FENode& node = ms.Node(pme->m_lnode[j]);
				umi[j] = node.m_rt - node.m_r0;
			}
			um = pme->eval(umi, ss.m_rs[i][0], ss.m_rs[i][1]);

			// calculate gap function
			ss.m_gap[i] = us - um;
		}
	}
}

//-----------------------------------------------------------------------------
void FEPeriodicBoundary::ShallowCopy(FEContactInterface &ci)
{
	FEPeriodicBoundary& si = dynamic_cast<FEPeriodicBoundary&>(ci);
	m_ss.ShallowCopy(si.m_ss);
	m_ms.ShallowCopy(si.m_ms);
}

//-----------------------------------------------------------------------------
void FEPeriodicBoundary::ContactForces(vector<double> &F)
{
	int j, k, l, m, n;
	int nseln, nmeln;

	double *Gr, *Gs;

	// jacobian
	double detJ;

	vec3d dxr, dxs;
	vec3d* rt, *r0;
	double* w;

	// natural coordinates of slave node in master element
	double r, s;

	// contact force
	vec3d tc;

	// shape function values
	double N[4];

	// element contact force vector
	vector<double> fe;

	// the lm array for this force vector
	vector<int> lm;

	// the en array
	vector<int> en;

	vector<int> sLM;
	vector<int> mLM;

	FEMesh& mesh = m_pfem->m_mesh;

	FESolidSolver* psolver = dynamic_cast<FESolidSolver*>(m_pfem->m_pStep->m_psolver);

	for (int np=0; np<m_npass; ++np)
	{
		FEPeriodicSurface& ss = (np == 0? m_ss : m_ms);
		FEPeriodicSurface& ms = (np == 0? m_ms : m_ss);

		// loop over all slave facets
		int ne = ss.Elements();
		for (j=0; j<ne; ++j)
		{
			// get the slave element
			FESurfaceElement& sel = ss.Element(j);
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

				// get slave node contact force
				tc = ss.m_Lm[m] + ss.m_gap[m]*m_eps;

				// get the master element
				FESurfaceElement& mel = dynamic_cast<FESurfaceElement&> (*ss.m_pme[m]);
				mesh.UnpackElement(mel, FE_UNPACK_LM);

				mLM = mel.LM();

				// we must unpack the slave element again
				mesh.UnpackElement(sel);

				nmeln = mel.Nodes();

				// isoparametric coordinates of the projected slave node
				// onto the master element
				r = ss.m_rs[m][0];
				s = ss.m_rs[m][1];

				// get the master shape function values at this slave node
				if (nmeln == 4)
				{
					// quadrilateral
					N[0] = 0.25*(1-r)*(1-s);
					N[1] = 0.25*(1+r)*(1-s);
					N[2] = 0.25*(1+r)*(1+s);
					N[3] = 0.25*(1-r)*(1+s);
				}
				else if (nmeln == 3)
				{
					// triangle
					N[0] = 1 - r - s;
					N[1] = r;
					N[2] = s;
				}

				// calculate force vector
				fe.create(3*(nmeln+1));
				fe[0] = -detJ*w[n]*tc.x;
				fe[1] = -detJ*w[n]*tc.y;
				fe[2] = -detJ*w[n]*tc.z;
				for (l=0; l<nmeln; ++l)
				{
					fe[3*(l+1)  ] = detJ*w[n]*tc.x*N[l];
					fe[3*(l+1)+1] = detJ*w[n]*tc.y*N[l];
					fe[3*(l+1)+2] = detJ*w[n]*tc.z*N[l];
				}

				// fill the lm array
				lm.create(3*(nmeln+1));
				lm[0] = sLM[n*3  ];
				lm[1] = sLM[n*3+1];
				lm[2] = sLM[n*3+2];

				for (l=0; l<nmeln; ++l)
				{
					lm[3*(l+1)  ] = mLM[l*3  ];
					lm[3*(l+1)+1] = mLM[l*3+1];
					lm[3*(l+1)+2] = mLM[l*3+2];
				}

				// fill the en array
				en.create(nmeln+1);
				en[0] = sel.m_node[n];
				for (l=0; l<nmeln; ++l) en[l+1] = mel.m_node[l];

				// assemble into global force vector
				psolver->AssembleResidual(en, lm, fe, F);
			}
		}
	}
}

//-----------------------------------------------------------------------------
void FEPeriodicBoundary::ContactStiffness()
{
	int j, k, l, n, m;
	int nseln, nmeln, ndof;

	matrix ke;

	vector<int> lm(15);
	vector<int> en(5);

	double *Gr, *Gs, *w;
	vec3d *rt, *r0;

	vec3d rtm[4];

	double detJ, r, s;
	vec3d dxr, dxs;
	double H[4];

	vec3d gap, Lm, tc;

	// curvature tensor K
	double K[2][2] = {0};

//	double scale = -0.0035*m_fem.m_mesh.GetBoundingBox().radius();

	vector<int> sLM;
	vector<int> mLM;

	FEMesh& mesh = m_pfem->m_mesh;

	FESolidSolver* psolver = dynamic_cast<FESolidSolver*>(m_pfem->m_pStep->m_psolver);

	for (int np=0; np<m_npass; ++np)
	{
		FEPeriodicSurface& ss = (np == 0? m_ss : m_ms);
		FEPeriodicSurface& ms = (np == 0? m_ms : m_ss);

		// loop over all slave elements
		int ne = ss.Elements();
		for (j=0; j<ne; ++j)
		{
			FESurfaceElement& se = ss.Element(j);
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

				// get the master element
				FESurfaceElement& me = dynamic_cast<FESurfaceElement&> (*ss.m_pme[m]);
				mesh.UnpackElement(me, FE_UNPACK_LM);

				mLM = me.LM();
				nmeln = me.Nodes();

				// get the master element node positions
				for (k=0; k<nmeln; ++k) rtm[k] = me.rt()[k];

				// we must unpack the slave element again
				mesh.UnpackElement(se);

				// slave node natural coordinates in master element
				r = ss.m_rs[m][0];
				s = ss.m_rs[m][1];

				// get slave node normal force
				tc = ss.m_Lm[m] + ss.m_gap[m]*m_eps; //ss.T[m];

				// get the master shape function values at this slave node
				if (nmeln == 4)
				{
					// quadrilateral
					H[0] = 0.25*(1-r)*(1-s);
					H[1] = 0.25*(1+r)*(1-s);
					H[2] = 0.25*(1+r)*(1+s);
					H[3] = 0.25*(1-r)*(1+s);
				}
				else if (nmeln == 3)
				{
					// triangle
					H[0] = 1 - r - s;
					H[1] = r;
					H[2] = s;
				}

				// number of degrees of freedom
				ndof = 3*(1 + nmeln);

				// fill stiffness matrix
				ke.Create(ndof, ndof); ke.zero();
				ke[0][0] = w[n]*detJ*m_eps;
				ke[1][1] = w[n]*detJ*m_eps;
				ke[2][2] = w[n]*detJ*m_eps;
				for (k=0; k<nmeln; ++k)
				{
					ke[0][3+3*k  ] = -w[n]*detJ*m_eps*H[k];
					ke[1][3+3*k+1] = -w[n]*detJ*m_eps*H[k];
					ke[2][3+3*k+2] = -w[n]*detJ*m_eps*H[k];

					ke[3+3*k  ][0] = -w[n]*detJ*m_eps*H[k];
					ke[3+3*k+1][1] = -w[n]*detJ*m_eps*H[k];
					ke[3+3*k+2][2] = -w[n]*detJ*m_eps*H[k];
				}
				for (k=0; k<nmeln; ++k)
					for (l=0; l<nmeln; ++l)
					{
						ke[3+3*k  ][3+3*l  ] = w[n]*detJ*m_eps*H[k]*H[l];
						ke[3+3*k+1][3+3*l+1] = w[n]*detJ*m_eps*H[k]*H[l];
						ke[3+3*k+2][3+3*l+2] = w[n]*detJ*m_eps*H[k]*H[l];
					}

				// create lm array
				lm[0] = sLM[n*3  ];
				lm[1] = sLM[n*3+1];
				lm[2] = sLM[n*3+2];

				for (k=0; k<nmeln; ++k)
				{
					lm[3*(k+1)  ] = mLM[k*3  ];
					lm[3*(k+1)+1] = mLM[k*3+1];
					lm[3*(k+1)+2] = mLM[k*3+2];
				}

				// create the en array
				en.create(nmeln+1);
				en[0] = se.m_node[n];
				for (k=0; k<nmeln; ++k) en[k+1] = me.m_node[k];

				// assemble stiffness matrix
				psolver->AssembleStiffness(en, lm, ke);
			}
		}
	}
}

//-----------------------------------------------------------------------------
bool FEPeriodicBoundary::Augment(int naug)
{
	int i;
	bool bconv = true;

	double g;
	vec3d lm;

	// calculate initial norms
	double normL0 = 0;
	for (i=0; i<m_ss.Nodes(); ++i)
	{
		lm = m_ss.m_Lm[i];
		normL0 += lm*lm;
	}
	for (i=0; i<m_ms.Nodes(); ++i)
	{
		lm = m_ms.m_Lm[i];
		normL0 += lm*lm;
	}
	normL0 = sqrt(normL0);

	// update Lagrange multipliers and calculate current norms
	double normL1 = 0;
	double normgc = 0;
	int N = 0;
	for (i=0; i<m_ss.Nodes(); ++i)
	{
		lm = m_ss.m_Lm[i] + m_ss.m_gap[i]*m_eps;

		normL1 += lm*lm;
		g = m_ss.m_gap[i].norm();
		normgc += g*g;
		++N;
	}
	for (i=0; i<m_ms.Nodes(); ++i)
	{
		lm = m_ms.m_Lm[i] + m_ms.m_gap[i]*m_eps;

		normL1 += lm*lm;
		g = m_ms.m_gap[i].norm();
		normgc += g*g;
		++N;
	}
	if (N == 0) N=1;

	normL1 = sqrt(normL1);
	normgc = sqrt(normgc / N);

	// get the logfile
	Logfile& log = GetLogfile();

	// check convergence of constraints
	log.printf(" tied interface # %d\n", m_nID);
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
			m_ss.m_Lm[i] = m_ss.m_Lm[i] + m_ss.m_gap[i]*m_eps;
		}
		for (i=0; i<m_ms.Nodes(); ++i)
		{
			// update Lagrange multipliers
			m_ms.m_Lm[i] = m_ms.m_Lm[i] + m_ms.m_gap[i]*m_eps;
		}
	}

	return bconv;
}

//-----------------------------------------------------------------------------
void FEPeriodicBoundary::Serialize(Archive &ar)
{
	int j, k, n, mat;

	if (ar.IsSaving())
	{
		ar << m_eps;
		ar << m_atol;

		for (j=0; j<2; ++j)
		{
			FEPeriodicSurface& s = (j==0? m_ss : m_ms);

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

			ar << s.m_gap;
			ar << s.m_rs;
			ar << s.m_Lm;
		}
	}
	else
	{
		ar >> m_eps;
		ar >> m_atol;

		for (j=0; j<2; ++j)
		{
			FEPeriodicSurface& s = (j==0? m_ss : m_ms);

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

			// read the contact data
			// Note that we do this after Init() since this data gets
			// initialized to zero there
			ar >> s.m_gap;
			ar >> s.m_rs;
			ar >> s.m_Lm;
		}
	}
}
