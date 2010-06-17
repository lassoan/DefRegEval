#include "stdafx.h"
#include "FEFacet2FacetSliding.h"
#include "fem.h"
#include "FESolidSolver.h"
#include "log.h"

//-----------------------------------------------------------------------------
// FEFacetSlidingSurface
//-----------------------------------------------------------------------------

void FEFacetSlidingSurface::Init()
{
	// initialize surface data first
	FESurface::Init();

	// count how many integration points we have
	int nint = 0, i;
	for (i=0; i<Elements(); ++i)
	{
		FESurfaceElement& el = Element(i);
		nint += el.GaussPoints();
	}

	// allocate data structures
	m_gap.create(nint);
	m_nu.create(nint);
	m_rs.create(nint);
	m_Lm.create(nint);
	m_pme.create(nint);
	m_eps.create(nint);

	// set intial values
	m_gap.zero();
	m_nu.zero();
	m_pme.set(0);
	m_Lm.zero();
	m_eps.set(1);
}

//-----------------------------------------------------------------------------
void FEFacetSlidingSurface::ShallowCopy(FEFacetSlidingSurface &s)
{
	m_Lm  = s.m_Lm;
	m_gap = s.m_gap;
	m_pme.zero();
}

//-----------------------------------------------------------------------------
//! Finds the (master) element that contains the projection of a (slave) node

FEElement* FEFacetSlidingSurface::FindMasterSegment(vec3d& x, vec3d& q, vec2d& r, bool& binit_nq, double tol)
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
// FEFacet2FacetSliding
//-----------------------------------------------------------------------------

FEFacet2FacetSliding::FEFacet2FacetSliding(FEM* pfem) : FEContactInterface(pfem), m_ss(&pfem->m_mesh), m_ms(&pfem->m_mesh)
{
	m_ntype = FE_FACET2FACET_SLIDING;
	static int ncount = 1;
	m_nID = ncount++;

	// default parameters
	m_epsn = 1.0;
	m_knmult = 1.0;
	m_stol = 0.01;
	m_npass = 1;
	m_bautopen = false;

	m_atol = 0.01;
	m_gtol = 0;
	m_naugmin = 0;
	m_naugmax = 10;
}

//-----------------------------------------------------------------------------
//! Initialization routine
void FEFacet2FacetSliding::Init()
{
	// initialize surface data
	m_ss.Init();
	m_ms.Init();

	// calculate penalty factors
	if (m_bautopen) CalcAutoPenalty(m_ss);

	// project slave surface onto master surface
	ProjectSurface(m_ss, m_ms);

	if (m_npass == 2) 
	{
		ProjectSurface(m_ms, m_ss);
		if (m_bautopen) CalcAutoPenalty(m_ms);
	}
}

//-----------------------------------------------------------------------------
void FEFacet2FacetSliding::CalcAutoPenalty(FEFacetSlidingSurface& s)
{
	// get the mesh
	FEMesh& m = m_pfem->m_mesh;

	// loop over all surface elements
	int ni = 0;
	for (int i=0; i<s.Elements(); ++i)
	{
		// get the surface element
		FESurfaceElement& el = s.Element(i);

		// find the element this face belongs to
		FEElement* pe = m.FindElementFromID(el.m_nelem);
		assert(pe);

		// get the area of the surface element
		double A = s.FaceArea(el);

		// get the volume of the volume element
		double V = m.ElementVolume(*pe);

		// calculate a modulus
		double E = AutoPenalty(el, s);

		// calculate penalty
		double eps = E*A/V;

		// assign to integation points of surface element
		int nint = el.GaussPoints();
		for (int j=0; j<nint; ++j, ++ni) s.m_eps[ni] = eps;
	}
}

//-----------------------------------------------------------------------------
//! In this function we project the integration points to the master surface,
//! calculate the projection's natural coordinates and normal vector
//
void FEFacet2FacetSliding::ProjectSurface(FEFacetSlidingSurface &ss, FEFacetSlidingSurface &ms)
{
	bool bfirst = true;

	// keep a running counter of which integration 
	// point we are currently investigating
	int ni = 0;

	// get the mesh
	FEMesh* pm = ss.GetMesh();

	// loop over all slave elements
	for (int i=0; i<ss.Elements(); ++i)
	{
		// get the slave element
		FESurfaceElement& se = ss.Element(i);
		pm->UnpackElement(se);
		vec3d* re = se.rt();
		int nn = se.Nodes();

		// loop over all its integration points
		int nint = se.GaussPoints();
		for (int j=0; j<nint; ++j, ++ni)
		{
			// calculate the global coordinates of this integration point
			double* H = se.H(j);

			vec3d x(0,0,0), q;
			for (int k=0; k<nn; ++k) x += re[k]*H[k];

			// find the master segment this element belongs to
			ss.m_rs[ni] = vec2d(0,0);
			FESurfaceElement* pme = dynamic_cast<FESurfaceElement*>(ms.FindMasterSegment(x, q, ss.m_rs[ni], bfirst, m_stol));
			ss.m_pme[ni] = pme;

			if (pme)
			{
				double r = ss.m_rs[ni][0];
				double s = ss.m_rs[ni][1];

				FESurfaceElement& mel = *pme;

				// the slave normal is set to the master element normal
				ss.m_nu[ni] = ms.SurfaceNormal(mel, r, s);

				// calculate gap
				ss.m_gap[ni] = -ss.m_nu[ni]*(x - q);
			}
			else
			{
				ss.m_gap[ni] = 0;
				ss.m_Lm[ni] = 0;
			}
		}
	}
}

//-----------------------------------------------------------------------------
void FEFacet2FacetSliding::Update()
{
	// project slave surface to master surface
	ProjectSurface(m_ss, m_ms);
	if (m_npass == 2) ProjectSurface(m_ms, m_ss);
}

//-----------------------------------------------------------------------------
void FEFacet2FacetSliding::ShallowCopy(FEContactInterface &ci)
{
	FEFacet2FacetSliding& si = dynamic_cast<FEFacet2FacetSliding&>(ci);
	m_ss.ShallowCopy(si.m_ss);
	m_ms.ShallowCopy(si.m_ms);
}

//-----------------------------------------------------------------------------
void FEFacet2FacetSliding::ContactForces(vector<double>& F)
{
	int i, j, k;
	vector<int> sLM, mLM, LM, en;
	vector<double> fe;

	// keep a running counter of integration points
	int ni = 0;

	// get the mesh
	FEMesh* pm = m_ss.GetMesh();

	// get the solver
	FESolidSolver* psolver = dynamic_cast<FESolidSolver*>(m_pfem->m_pStep->m_psolver);

	double detJ[4], w[4], *Hs, Hm[4];

	for (int np=0; np<m_npass; ++np)
	{
		FEFacetSlidingSurface& ss = (np == 0? m_ss : m_ms);
		FEFacetSlidingSurface& ms = (np == 0? m_ms : m_ss);

		// loop over all slave elements
		ni = 0;
		for (i=0; i<ss.Elements(); ++i)
		{
			FESurfaceElement& se = ss.Element(i);
			pm->UnpackElement(se);
			int nseln = se.Nodes();
			int nint = se.GaussPoints();

			// copy the LM vector
			sLM = se.LM();

			// nodal coordinates
			vec3d* r0 = se.r0();

			// we calculate all the metrics we need before we
			// calculate the nodal forces
			for (j=0; j<nint; ++j)
			{
				double* Gr = se.Gr(j);
				double* Gs = se.Gs(j);

				// calculate jacobian
				// note that we are integrating over the reference surface
				vec3d dxr, dxs;
				for (k=0; k<nseln; ++k)
				{
					dxr.x += Gr[k]*r0[k].x;
					dxr.y += Gr[k]*r0[k].y;
					dxr.z += Gr[k]*r0[k].z;

					dxs.x += Gs[k]*r0[k].x;
					dxs.y += Gs[k]*r0[k].y;
					dxs.z += Gs[k]*r0[k].z;
				}

				// jacobians
				detJ[j] = (dxr ^ dxs).norm();

				// integration weights
				w[j] = se.GaussWeights()[j];
			}

			// loop over all integration points
			for (int j=0; j<nint; ++j, ++ni)
			{
				// get the master element
				FESurfaceElement* pme = ss.m_pme[ni];
				if (pme)
				{
					FESurfaceElement& me = *pme;
					pm->UnpackElement(me);

					int nmeln = me.Nodes();
					mLM = me.LM();

					// calculate degrees of freedom
					int ndof = 3*(nseln + nmeln);

					// build the LM vector
					LM.create(ndof);
					for (k=0; k<nseln; ++k)
					{
						LM[3*k  ] = sLM[3*k  ];
						LM[3*k+1] = sLM[3*k+1];
						LM[3*k+2] = sLM[3*k+2];
					}

					for (k=0; k<nmeln; ++k)
					{
						LM[3*(k+nseln)  ] = mLM[3*k  ];
						LM[3*(k+nseln)+1] = mLM[3*k+1];
						LM[3*(k+nseln)+2] = mLM[3*k+2];
					}

					// build the en vector
					en.create(nseln+nmeln);
					for (k=0; k<nseln; ++k) en[k] = se.m_node[k];
					for (k=0; k<nmeln; ++k) en[k+nseln] = me.m_node[k];

					// calculate shape functions
					Hs = se.H(j);

					double r = ss.m_rs[ni][0];
					double s = ss.m_rs[ni][1];
					if (me.Nodes() == 4)
					{
						Hm[0] = 0.25*(1-r)*(1-s);
						Hm[1] = 0.25*(1+r)*(1-s);
						Hm[2] = 0.25*(1+r)*(1+s);
						Hm[3] = 0.25*(1-r)*(1+s);
					}

					// get normal vector
					vec3d nu = ss.m_nu[ni];

					// gap function
					double g = ss.m_gap[ni];
					
					// lagrange multiplier
					double Lm = ss.m_Lm[ni];

					// penalty value
					double eps = m_epsn*ss.m_eps[ni];

					// contact traction
					double tn = Lm + eps*g;
					tn = MBRACKET(tn);

					// calculate the force vector
					fe.create(ndof);

					for (k=0; k<nseln; ++k)
					{
						fe[3*k  ] = Hs[k]*nu.x;
						fe[3*k+1] = Hs[k]*nu.y;
						fe[3*k+2] = Hs[k]*nu.z;
					}

					for (k=0; k<nmeln; ++k)
					{
						fe[3*(k+nseln)  ] = -Hm[k]*nu.x;
						fe[3*(k+nseln)+1] = -Hm[k]*nu.y;
						fe[3*(k+nseln)+2] = -Hm[k]*nu.z;
					}

					for (k=0; k<ndof; ++k) fe[k] *= tn*detJ[j]*w[j];

					// assemble the global residual
					psolver->AssembleResidual(en, LM, fe, F);
				}
			}
		}
	}
}

//-----------------------------------------------------------------------------

void FEFacet2FacetSliding::ContactStiffness()
{
	int i, j, k, l;
	vector<int> sLM, mLM, LM, en;
	double N[24], T1[24], T2[24], N1[24] = {0}, N2[24] = {0}, D1[24], D2[24], Nb1[24], Nb2[24];
	matrix ke;

	// keep a running counter of integration points
	int ni = 0;

	// get the mesh
	FEMesh* pm = m_ss.GetMesh();

	// get the solver
	FESolidSolver* psolver = dynamic_cast<FESolidSolver*>(m_pfem->m_pStep->m_psolver);

	// get the logfile
	Logfile& log = GetLogfile();

	// see how many reformations we've had to do so far
	int nref = psolver->m_nref;

	// set higher order stiffness mutliplier
	double knmult = m_knmult;
	if (m_knmult < 0)
	{
		int ni = int(-m_knmult);
		if (nref >= ni)
		{
			knmult = 1; 
			log.printf("Higher order stiffness terms included.\n");
		}
		else knmult = 0;
	}

	double detJ[4], w[4], *Hs, Hm[4], Hmr[4], Hms[4];

	for (int np=0; np < m_npass; ++np)
	{
		FEFacetSlidingSurface& ss = (np == 0? m_ss : m_ms);
		FEFacetSlidingSurface& ms = (np == 0? m_ms : m_ss);

		// loop over all slave elements
		ni = 0;
		for (i=0; i<ss.Elements(); ++i)
		{
			FESurfaceElement& se = ss.Element(i);
			pm->UnpackElement(se);
			int nseln = se.Nodes();
			int nint = se.GaussPoints();

			// copy the LM vector
			sLM = se.LM();

			// nodal coordinates
			vec3d* r0 = se.r0();

			// we calculate all the metrics we need before we
			// calculate the nodal forces
			for (j=0; j<nint; ++j)
			{
				double* Gr = se.Gr(j);
				double* Gs = se.Gs(j);

				// calculate jacobian
				// note that we are integrating over the reference surface
				vec3d dxr, dxs;
				for (k=0; k<nseln; ++k)
				{
					dxr.x += Gr[k]*r0[k].x;
					dxr.y += Gr[k]*r0[k].y;
					dxr.z += Gr[k]*r0[k].z;

					dxs.x += Gs[k]*r0[k].x;
					dxs.y += Gs[k]*r0[k].y;
					dxs.z += Gs[k]*r0[k].z;
				}

				// jacobians
				detJ[j] = (dxr ^ dxs).norm();

				// integration weights
				w[j] = se.GaussWeights()[j];
			}

			// loop over all integration points
			for (int j=0; j<nint; ++j, ++ni)
			{
				// get the master element
				FESurfaceElement* pme = ss.m_pme[ni];
				if (pme)
				{
					FESurfaceElement& me = *pme;
					pm->UnpackElement(me);

					int nmeln = me.Nodes();
					mLM = me.LM();

					// calculate degrees of freedom
					int ndof = 3*(nseln + nmeln);

					// build the LM vector
					LM.create(ndof);
					for (k=0; k<nseln; ++k)
					{
						LM[3*k  ] = sLM[3*k  ];
						LM[3*k+1] = sLM[3*k+1];
						LM[3*k+2] = sLM[3*k+2];
					}

					for (k=0; k<nmeln; ++k)
					{
						LM[3*(k+nseln)  ] = mLM[3*k  ];
						LM[3*(k+nseln)+1] = mLM[3*k+1];
						LM[3*(k+nseln)+2] = mLM[3*k+2];
					}

					// build the en vector
					en.create(nseln+nmeln);
					for (k=0; k<nseln; ++k) en[k] = se.m_node[k];
					for (k=0; k<nmeln; ++k) en[k+nseln] = me.m_node[k];

					// calculate shape functions
					Hs = se.H(j);
					double r = ss.m_rs[ni][0];
					double s = ss.m_rs[ni][1];
					if (me.Nodes() == 4)
					{
						Hm[0] = 0.25*(1-r)*(1-s);
						Hm[1] = 0.25*(1+r)*(1-s);
						Hm[2] = 0.25*(1+r)*(1+s);
						Hm[3] = 0.25*(1-r)*(1+s);
					}
					else
					{
						Hm[0] = 1 - r - s;
						Hm[1] = r;
						Hm[2] = s;
					}

					// get normal vector
					vec3d nu = ss.m_nu[ni];

					// gap function
					double g = ss.m_gap[ni];
					
					// lagrange multiplier
					double Lm = ss.m_Lm[ni];

					// penalty value
					double eps = m_epsn*ss.m_eps[ni];

					// contact traction
					double tn = Lm + eps*g;
					tn = MBRACKET(tn);

					double dtn = eps*HEAVYSIDE(Lm + eps*g);

					// calculate the N-vector
					for (k=0; k<nseln; ++k)
					{
						N[3*k  ] = Hs[k]*nu.x;
						N[3*k+1] = Hs[k]*nu.y;
						N[3*k+2] = Hs[k]*nu.z;
					}

					for (k=0; k<nmeln; ++k)
					{
						N[3*(k+nseln)  ] = -Hm[k]*nu.x;
						N[3*(k+nseln)+1] = -Hm[k]*nu.y;
						N[3*(k+nseln)+2] = -Hm[k]*nu.z;
					}

					// --- N O R M A L   S T I F F N E S S ---

					// create the stiffness matrix
					ke.Create(ndof, ndof);

					// add the first order term (= D(tn)*dg )
					for (k=0; k<ndof; ++k)
						for (l=0; l<ndof; ++l) ke[k][l] = dtn*N[k]*N[l]*detJ[j]*w[j];

					// add the higher order terms (= tn*D(dg) )
					if (knmult > 0)
					{
						// calculate the master shape fncs derivatives
						if (nmeln == 4)
						{
							Hmr[0] = -0.25*(1-s); Hms[0] = -0.25*(1-r);
							Hmr[1] =  0.25*(1-s); Hms[1] = -0.25*(1+r);
							Hmr[2] =  0.25*(1+s); Hms[2] =  0.25*(1+r);
							Hmr[3] = -0.25*(1+s); Hms[3] =  0.25*(1-r);
						}
						else
						{
							Hmr[0] = -1; Hms[0] = -1;
							Hmr[1] =  1; Hms[1] =  0;
							Hmr[2] =  0; Hms[2] =  1;
						}

						// get the master nodes
						vec3d* rt = me.rt();

						// get the tangent vectors
						vec3d tau1(0,0,0), tau2(0,0,0);
						for (k=0; k<nmeln; ++k)
						{
							tau1.x += Hmr[k]*rt[k].x;
							tau1.y += Hmr[k]*rt[k].y;
							tau1.z += Hmr[k]*rt[k].z;
		
							tau2.x += Hms[k]*rt[k].x;
							tau2.y += Hms[k]*rt[k].y;
							tau2.z += Hms[k]*rt[k].z;
						}

						// set up the Ti vectors
						for (k=0; k<nseln; ++k)
						{
							T1[k*3  ] = Hs[k]*tau1.x; T2[k*3  ] = Hs[k]*tau2.x;
							T1[k*3+1] = Hs[k]*tau1.y; T2[k*3+1] = Hs[k]*tau2.y;
							T1[k*3+2] = Hs[k]*tau1.z; T2[k*3+2] = Hs[k]*tau2.z;
						}

						for (k=0; k<nmeln; ++k) 
						{
							T1[(k+nseln)*3  ] = -Hm[k]*tau1.x;
							T1[(k+nseln)*3+1] = -Hm[k]*tau1.y;
							T1[(k+nseln)*3+2] = -Hm[k]*tau1.z;

							T2[(k+nseln)*3  ] = -Hm[k]*tau2.x;
							T2[(k+nseln)*3+1] = -Hm[k]*tau2.y;
							T2[(k+nseln)*3+2] = -Hm[k]*tau2.z;
						}

						// set up the Ni vectors
						for (k=0; k<nmeln; ++k) 
						{
							N1[(k+nseln)*3  ] = -Hmr[k]*nu.x;
							N1[(k+nseln)*3+1] = -Hmr[k]*nu.y;
							N1[(k+nseln)*3+2] = -Hmr[k]*nu.z;

							N2[(k+nseln)*3  ] = -Hms[k]*nu.x;
							N2[(k+nseln)*3+1] = -Hms[k]*nu.y;
							N2[(k+nseln)*3+2] = -Hms[k]*nu.z;
						}

						// calculate metric tensor
						mat2d M;
						M[0][0] = tau1*tau1; M[0][1] = tau1*tau2; 
						M[1][0] = tau2*tau1; M[1][1] = tau2*tau2; 

						// calculate reciprocal metric tensor
						mat2d Mi = M.inverse();

						// calculate curvature tensor
						double K[2][2] = {0};
						const double Grs[4] = {0.25, -0.25, 0.25, -0.25};
						if (nmeln == 4)
						{
							for (k=0; k<nmeln; ++k)
							{
								K[0][1] += (nu*rt[k])*Grs[k];
								K[1][0] += (nu*rt[k])*Grs[k];
							}
						}

						// setup A matrix A = M + gK
						double A[2][2];
						A[0][0] = M[0][0] + g*K[0][0];
						A[0][1] = M[0][1] + g*K[0][1];
						A[1][0] = M[1][0] + g*K[1][0];
						A[1][1] = M[1][1] + g*K[1][1];

						// calculate determinant of A
						double detA = A[0][0]*A[1][1] - A[0][1]*A[1][0];

						// setup Di vectors
						for (k=0; k<ndof; ++k)
						{
							D1[k] = (1/detA)*(A[1][1]*(T1[k]+g*N1[k]) - A[0][1]*(T2[k] + g*N2[k]));
							D2[k] = (1/detA)*(A[0][0]*(T2[k]+g*N2[k]) - A[0][1]*(T1[k] + g*N1[k]));
						}

						// setup Nbi vectors
						for (k=0; k<ndof; ++k)
						{
							Nb1[k] = N1[k] - K[0][1]*D2[k];
							Nb2[k] = N2[k] - K[0][1]*D1[k];
						}

						// add it to the stiffness
						double sum;
						for (k=0; k<ndof; ++k)
							for (l=0; l<ndof; ++l)
							{
								sum = Mi[0][0]*Nb1[k]*Nb1[l]+Mi[0][1]*(Nb1[k]*Nb2[l]+Nb2[k]*Nb1[l])+Mi[1][1]*Nb2[k]*Nb2[l];
								sum *= g;
								sum -= D1[k]*N1[l]+D2[k]*N2[l]+N1[k]*D1[l]+N2[k]*D2[l];
								sum += K[0][1]*(D1[k]*D2[l]+D2[k]*D1[l]);
								sum *= tn*knmult;

								ke[k][l] += sum*detJ[j]*w[j];
							}
					}

					// assemble the global residual
					psolver->AssembleStiffness(en, LM, ke);
				}
			}
		}
	}
}

//-----------------------------------------------------------------------------
bool FEFacet2FacetSliding::Augment(int naug)
{
	// make sure we need to augment
	if (!m_blaugon) return true;

	int i;
	double Ln;
	bool bconv = true;

	static double normg0 = 0;

	// --- c a l c u l a t e   i n i t i a l   n o r m s ---
	// a. normal component
	int NS = m_ss.m_Lm.size();
	int NM = m_ms.m_Lm.size();
	double normL0 = 0;
	for (i=0; i<NS; ++i) normL0 += m_ss.m_Lm[i]*m_ss.m_Lm[i];
	for (i=0; i<NM; ++i) normL0 += m_ms.m_Lm[i]*m_ms.m_Lm[i];
	normL0 = sqrt(normL0);

	// --- c a l c u l a t e   c u r r e n t   n o r m s ---
	// a. normal component
	double normL1 = 0;	// force norm
	double normg1 = 0;	// gap norm
	int N = 0;
	for (i=0; i<NS; ++i)
	{
		// penalty value
		double eps = m_epsn*m_ss.m_eps[i];

		// update Lagrange multipliers
		Ln = m_ss.m_Lm[i] + eps*m_ss.m_gap[i];
		Ln = MBRACKET(Ln);

		normL1 += Ln*Ln;

		if (m_ss.m_gap[i] > 0)
		{
			normg1 += m_ss.m_gap[i]*m_ss.m_gap[i];
			++N;
		}
	}	

	for (i=0; i<NM; ++i)
	{
		// penalty value
		double eps = m_epsn*m_ms.m_eps[i];

		// update Lagrange multipliers
		Ln = m_ms.m_Lm[i] + eps*m_ms.m_gap[i];
		Ln = MBRACKET(Ln);

		normL1 += Ln*Ln;
		if (m_ms.m_gap[i] > 0)
		{
			normg1 += m_ms.m_gap[i]*m_ms.m_gap[i];
			++N;
		}
	}
	if (N == 0) N=1;

	normL1 = sqrt(normL1);
	normg1 = sqrt(normg1 / N);

	if (naug == 0) normg0 = 0;

	// get the logfile
	Logfile& log = GetLogfile();

	// calculate and print convergence norms
	double lnorm = 0, gnorm = 0;
	if (normL1 != 0) lnorm = fabs(normL1 - normL0)/normL1; else lnorm = fabs(normL1 - normL0);
	if (normg1 != 0) gnorm = fabs(normg1 - normg0)/normg1; else gnorm = fabs(normg1 - normg0);

	log.printf(" sliding interface # %d\n", m_nID);
	log.printf("                        CURRENT        REQUIRED\n");
	log.printf("    normal force : %15le", lnorm);
	if (m_atol > 0) log.printf("%15le\n", m_atol); else log.printf("       ***\n");
	log.printf("    gap function : %15le", gnorm);
	if (m_gtol > 0) log.printf("%15le\n", m_gtol); else log.printf("       ***\n");

	// check convergence
	bconv = true;
	if ((m_atol > 0) && (lnorm > m_atol)) bconv = false;
	if ((m_gtol > 0) && (gnorm > m_gtol)) bconv = false;
	if (m_naugmin > naug) bconv = false;
	if (m_naugmax <= naug) bconv = true;
		
	if (bconv == false)
	{
		// we did not converge so update multipliers
		for (i=0; i<NS; ++i)
		{
			// penalty value
			double eps = m_epsn*m_ss.m_eps[i];

			// update Lagrange multipliers
			Ln = m_ss.m_Lm[i] + eps*m_ss.m_gap[i];
			m_ss.m_Lm[i] = MBRACKET(Ln);
		}	

		for (i=0; i<NM; ++i)
		{
			// penalty value
			double eps = m_epsn*m_ms.m_eps[i];

			// update Lagrange multipliers
			Ln = m_ms.m_Lm[i] + eps*m_ms.m_gap[i];
			m_ms.m_Lm[i] = MBRACKET(Ln);
		}
	}

	// store the last gap norm
	normg0 = normg1;

	return bconv;
}

//-----------------------------------------------------------------------------
void FEFacet2FacetSliding::Serialize(Archive &ar)
{

}
