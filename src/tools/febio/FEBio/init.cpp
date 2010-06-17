#include "stdafx.h"
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <time.h>
#include "fem.h"
#include "FEException.h"
#include "FENodeReorder.h"
#include "FERigid.h"
#include "log.h"
#include "FESolidSolver.h"
#include "LSDYNAPlotFile.h"

// Forward declarations
void Hello(FILE* fp);

bool err(const char* sz, ...)
{
	// get a pointer to the argument list
	va_list	args;

	// print to file
	va_start(args, sz);
	vfprintf(stderr, sz, args);
	va_end(args);

	return false;
}

//-----------------------------------------------------------------------------
//! This function performs one-time-initialization stuff. All the different 
//! modules are initialized here as well. This routine also performs some
//! data checks

bool FEM::Init()
{
	int i;

	// get the logfile
	Logfile& log = GetLogfile();

	// Open the logfile
	if (!log.is_valid()) 
	{
		log.open(m_szlog);

		// if we don't want to output anything we only output to the logfile
		if (m_pStep->GetPrintLevel() == FE_PRINT_NEVER) log.SetMode(Logfile::FILE_ONLY);

		// print welcome message to file
		Hello(log);
	}

	// check step data
	for (int i=0; i<m_Step.size(); ++i)
	{
		FEAnalysis& step = m_Step[i];
		if (step.m_ntime <= 0) return err("Invalid number of time steps for analysis step %d", i+1);
		if (step.m_dt0   <= 0) return err("Invalid time step size for analysis step %d", i+1);
		if (step.m_bautostep)
		{
//			if (m_pStep->m_dtmin <= 0) return err("Invalid minimum time step size");
//			if (m_pStep->m_dtmax <= 0) return err("Invalid maximum time step size");
		}
	}

	// evaluate all loadcurves at the initial time
	for (i=0; i<LoadCurves(); ++i) m_LC[i].Evaluate(0);

	// if the analysis is run in plain-strain mode we fix all the z-dofs of all nodes
	if (m_bplane_strain)
	{
		for (int i=0; i<m_mesh.Nodes(); ++i) m_mesh.Node(i).m_ID[2] = -1;
	}

	// find and remove isolated vertices
	int ni = m_mesh.RemoveIsolatedVertices();
	if (ni != 0) 
	{
		if (ni == 1)
			log.printbox("WARNING", "%d isolated vertex removed.", ni);
		else
			log.printbox("WARNING", "%d isolated vertices removed.", ni);
	}

	// initialize rigid body data
	if (InitRigidBodies() == false) return false;

	// initialize poro-elastic data
	if (InitPoro() == false) return false;

	// initialize heat condition data
	if (InitHeat() == false) return false;

	// initialize random number generator
	srand((unsigned) time(NULL));

	// initialize mesh data
	// note that this must be done AFTER the elements have been assigned material point data !
	// this is because the mesh data is reset
	// TODO: perhaps I should not reset the mesh data during the initialization
	if (m_mesh.Init() == false) return false;

	// intialize local coordinate data
	bool bmerr = false;
	for (i=0; i<m_mesh.SolidElements(); ++i)
	{
		// unpack element data
		FESolidElement& el = m_mesh.SolidElement(i);

		try
		{
			m_mesh.UnpackElement(el);
		}
		catch (NegativeJacobian e)
		{
			log.printbox("F A T A L   E R R O R", "A negative jacobian was detected at\n integration point %d of element %d.\nDid you use the right node numbering?", e.m_iel, e.m_ng);
			return false;
		}

		if (dynamic_cast<FESolidSolver*>(m_pStep->m_psolver))
		{
			// get the elements material
			FEElasticMaterial* pme = GetElasticMaterial(el.GetMatID());

			// set the local element coordinates
			if (pme)
			{
				if (pme->m_pmap)
				{
					for (int n=0; n<el.GaussPoints(); ++n)
					{
						FEElasticMaterialPoint& pt = *el.m_State[n]->ExtractData<FEElasticMaterialPoint>();
						pt.Q = pme->m_pmap->LocalElementCoord(el, n);
					}
				}
				else
				{
					if (GetDebugFlag())
					{
						// If we get here, then the element has a user-defined fiber axis
						// we should check to see if it has indeed been specified.
						// TODO: This assumes that pt.Q will not get intialized to
						//		 a valid value. I should find another way for checking since I
						//		 would like pt.Q always to be initialized to a decent value.
						if (dynamic_cast<FETransverselyIsotropic*>(pme))
						{
							FEElasticMaterialPoint& pt = *el.m_State[0]->ExtractData<FEElasticMaterialPoint>();
							mat3d& m = pt.Q;
							if (fabs(m.det() - 1) > 1e-7)
							{
								// this element did not get specified a user-defined fiber direction
								log.printbox("ERROR", "Solid element %d was not assigned a fiber direction.", i+1);
								bmerr = true;
							}
						}
					}
				}
			}
		}
	}
	if (bmerr) return false;

	// now do the same thing for shells
	bmerr = false;
	for (i=0; i<m_mesh.ShellElements(); ++i)
	{
		// unpack element data
		FEShellElement& el = m_mesh.ShellElement(i);

		// get the elements material
		FEElasticMaterial* pme = GetElasticMaterial(el.GetMatID());

		// set the local element coordinates
		if (pme)
		{
			if (pme->m_pmap)
			{
				for (int n=0; n<el.GaussPoints(); ++n)
				{
					FEElasticMaterialPoint& pt = *el.m_State[n]->ExtractData<FEElasticMaterialPoint>();
					pt.Q = pme->m_pmap->LocalElementCoord(el, n);
				}
			}
			else
			{
				if (GetDebugFlag())
				{
					// If we get here, then the element has a user-defined fiber direction
					// we should check to see if it has indeed been specified.
					// TODO: This assumes that pt.Q will not get intialized to
					//		 a valid value. I should find another way for checking since I
					//		 would like pt.Q always to be initialized to a decent value.
					if (dynamic_cast<FETransverselyIsotropic*>(pme))
					{
						FEElasticMaterialPoint& pt = *el.m_State[0]->ExtractData<FEElasticMaterialPoint>();
						mat3d& m = pt.Q;
						if (fabs(m.det() - 1) > 1e-7)
						{
							// this element did not get specified a user-defined fiber direction
							log.printbox("ERROR", "Shell element %d was not assigned a fiber direction.", i+1);
							bmerr = true;
						}
					}
				}
			}
		}
	}
	if (bmerr) return false;

	// initialize material data
	char szmat[256] = "Invalid value for material parameter \"%s\" of material %d";
	for (i=0; i<Materials(); ++i)
	{
		// get the material
		FEMaterial* pmat = GetMaterial(i);

		// initialize material data
		try
		{
			pmat->Init();
		}
		catch (MaterialError e)
		{
			log.printf("Failed initializing material %d (name=\"%s\"):\n", i+1, pmat->GetName());
			log.printf("ERROR: %s\n\n", e.Error());
			return false;
		}
		catch (...)
		{
			log.printf("A fatal error occured during material intialization\n\n");
			return false;
		}

		// set the activation load curve
		FETransverselyIsotropic* pm = dynamic_cast<FETransverselyIsotropic*> (pmat);
		if (pm)
		{
			if (pm->m_fib.m_lcna >= 0) pm->m_fib.m_plc = GetLoadCurve(pm->m_fib.m_lcna);
		}
	}

	// initialize contact data
	if (InitContact() == false) return false;

	// check discrete elements
	for (i=0; i<m_DE.size(); ++i)
	{
		FE_DISCRETE_ELEMENT& el = m_DE[i];
		if (el.n1 <0 || el.n1 >= m_mesh.Nodes() ||
			el.n2 <0 || el.n2 >= m_mesh.Nodes() )
		{
			log.printf("Invalid node number for discrete element %d\n", i+1);
			return false;
		}

		if (el.E < 0)
		{
			log.printbox("INPUT ERROR", "Invalide value for sping constant of discrete element %d", i+1);
			return false;
		}
	}

	// intitialize time
	m_ftime = 0;

	// open plot database file
	if (m_pStep->m_nplot != FE_PLOT_NEVER)
	{
		if (m_plot == 0) m_plot = new LSDYNAPlotFile;

		if (m_plot->Open(*this, m_szplot) == false)
		{
			log.printf("ERROR : Failed creating PLOT database\n");
			return false;
		}
	}

	// Since it is assumed that for the first timestep
	// there are no loads or initial displacements, the case n=0 is skipped.
	// Therefor we can output those results here.
	// Offcourse we should actually check if this is indeed
	// the case, otherwise we should also solve for t=0
	if (m_pStep->m_nplot != FE_PLOT_NEVER) m_plot->Write(*this);

	// do the callback
	DoCallback();

	// Alright, all initialization is done, so let's get busy !
	return true;
}

//-----------------------------------------------------------------------------
//! This function initializes the linear constraint table (LCT). This table
//! contains for each dof the linear constraint it belongs to. (or -1 if it is
//! not constraint)

bool FEM::InitConstraints()
{
	int nlin = m_LinC.size();
	if (nlin == 0) return true;

	int i;

	// set the equation numbers for the linear constraints
	list<FELinearConstraint>::iterator it = m_LinC.begin();
	for (i=0; i<nlin; ++i, ++it)
	{
		FELinearConstraint& lc = *it;
		lc.master.neq = m_mesh.Node(lc.master.node).m_ID[lc.master.bc];

		// make sure the master did not get assigned an equation
		assert(lc.master.neq == -1);

		// set the slave equation numbers
		list<FELinearConstraint::SlaveDOF>::iterator is = lc.slave.begin();
		int nn = lc.slave.size();
		for (int n=0; n<nn; ++n, ++is)
		{
			FELinearConstraint::SlaveDOF& sn = *is;
			sn.neq = m_mesh.Node(sn.node).m_ID[sn.bc];
		}		
	}

	// create the linear constraint table
	m_LCT.create(m_mesh.Nodes()*MAX_NDOFS);

	// initialize it to -1
	m_LCT.set(-1);

	list<FELinearConstraint>::iterator ic = m_LinC.begin();
	for (i=0; i<nlin; ++i, ++ic)
	{
		FELinearConstraint& lc = *ic;
		int n = lc.master.node;
		int m = lc.master.bc;
		
		m_LCT[n*MAX_NDOFS+m] = i;
	}

	// to simplify accessing the linear constraint data
	// we store all pointers in an array
	// TODO: perhaps I should store the linear constraints that way
	// anyways and get rid of the list
	m_LCA.create(nlin);
	ic = m_LinC.begin();
	for (i=0; i<nlin; ++i, ++ic) m_LCA[i] = &(*ic);

	// let's do the aug lag linear constraints
	if (m_LCSet.size() > 0)
	{
		int N = m_LCSet.size();
		list<FELinearConstraintSet*>::iterator it = m_LCSet.begin();
		for (i=0; i<N; ++i, ++it) (*it)->Init();
	}

	return true;
}

//-----------------------------------------------------------------------------
//!	This function initializes the equation system.
//! It is assumed that all free dofs up until now have been given an ID >= 0
//! and the fixed or rigid dofs an ID < 0.
//! After this operation the nodal ID array will contain the equation
//! number assigned to the corresponding degree of freedom. To distinguish
//! between free or unconstrained dofs and constrained ones the following rules
//! apply to the ID array:
//!
//!           /
//!          |  >=  0 --> dof j of node i is a free dof
//! ID[i][j] <  == -1 --> dof j of node i is a fixed (no equation assigned too)
//!          |  <  -1 --> dof j of node i is constrained and has equation nr = -ID[i][j]-2
//!           \
//!

bool FEM::InitEquations()
{
	int i, j, n;

	// initialize nr of equations
	m_neq=0;

	// see if we need to optimize the bandwidth
	if (m_bwopt)
	{
		// reorder the node numbers
		vector<int> P(m_mesh.Nodes());
		FENodeReorder mod;
		mod.Apply(m_mesh, P);

		// set the equation numbers
		for (i=0; i<m_mesh.Nodes(); ++i)
		{
			FENode& node = m_mesh.Node(P[i]);
			for (j=0; j<MAX_NDOFS; ++j)
				if (node.m_ID[j] >= 0) node.m_ID[j] = m_neq++;
		}
	}
	else
	{
		// give all free dofs an equation number
		for (i=0; i<m_mesh.Nodes(); ++i)
		{
			FENode& node = m_mesh.Node(i);
			for (j=0; j<MAX_NDOFS; ++j)
				if (node.m_ID[j] >= 0) node.m_ID[j] = m_neq++;
		}
	}


	// Next, we assign equation numbers to the rigid body degrees of freedom
	m_nreq = m_neq;
	for (i=0; i<m_nrb; ++i)
	{
		FERigidBody& RB = m_RB[i];
		FERigid* pm = dynamic_cast<FERigid*>(GetMaterial(RB.m_mat));
		assert(pm);
		for (j=0; j<6; ++j)
			if (pm->m_bc[j] >= 0)
			{
				RB.m_LM[j] = m_neq++;
			}
			else 
				RB.m_LM[j] = -1;
	}

	// we assign the rigid body equation number to
	// Also make sure that the nodes are NOT constrained!
	for (i=0; i<m_mesh.Nodes(); ++i)
	{
		FENode& node = m_mesh.Node(i);
		if (node.m_rid >= 0)
		{
			FERigidBody& RB = m_RB[node.m_rid];
			node.m_ID[0] = -RB.m_LM[0]-2;
			node.m_ID[1] = -RB.m_LM[1]-2;
			node.m_ID[2] = -RB.m_LM[2]-2;
			node.m_ID[7] = -RB.m_LM[3]-2;
			node.m_ID[8] = -RB.m_LM[4]-2;
			node.m_ID[9] = -RB.m_LM[5]-2;
		}
	}

	// adjust the rigid dofs that are prescribed
	for (i=0; i<m_nrb; ++i)
	{
		FERigidBody& RB = m_RB[i];
		FERigid* pm = dynamic_cast<FERigid*>(GetMaterial(RB.m_mat));
		for (j=0; j<6; ++j)
		{
			n = RB.m_LM[j];
			if (pm->m_bc[j] > 0) RB.m_LM[j] = -n-2;
		}
	}

	// All initialization is done
	return true;
}

//-----------------------------------------------------------------------------
//! This function resets the FEM data so that a new run can be done.
//! This routine is called from the optimization routine.

bool FEM::Reset()
{
	int i;

	// reset mesh data
	m_mesh.Reset();

	// initialize rigid body data
	for (i=0; i<m_nrb; ++i)
	{
		// zero total displacements
		m_RB[i].m_Ut[0] = m_RB[i].m_Up[0] = 0;
		m_RB[i].m_Ut[1] = m_RB[i].m_Up[1] = 0;
		m_RB[i].m_Ut[2] = m_RB[i].m_Up[2] = 0;
		m_RB[i].m_Ut[3] = m_RB[i].m_Up[3] = 0;
		m_RB[i].m_Ut[4] = m_RB[i].m_Up[4] = 0;
		m_RB[i].m_Ut[5] = m_RB[i].m_Up[5] = 0;

		// initialize orientation
		m_RB[i].m_qt = quatd(0, vec3d(0,0,1));

		// initialize center of mass
		m_RB[i].m_rt = m_RB[i].m_r0;

		// reset reaction force and torque
		m_RB[i].m_Fr = vec3d(0,0,0);
		m_RB[i].m_Mr = vec3d(0,0,0);
	}

	// set up rigid joints
	if (m_nrj > 0)
	{
		FERigid* pm;
		for (i=0; i<m_nrj; ++i)
		{
			FERigidJoint& rj = m_RJ[i];
			rj.m_F = vec3d(0,0,0);

			pm = dynamic_cast<FERigid*> (GetMaterial(rj.m_nRBa));
			rj.m_nRBa = pm->m_nRB;

			pm = dynamic_cast<FERigid*> (GetMaterial(rj.m_nRBb));
			rj.m_nRBb = pm->m_nRB;

			rj.m_qa0 = rj.m_q0 - m_RB[ rj.m_nRBa ].m_r0;
			rj.m_qb0 = rj.m_q0 - m_RB[ rj.m_nRBb ].m_r0;
		}
	}

	// set the start time
	m_ftime = 0;

	// set first time step
	m_pStep->m_dt = m_pStep->m_dt0;

	// get the logfile
	Logfile& log = GetLogfile();

	// open plot database file
	if (m_pStep->m_nplot != FE_PLOT_NEVER)
	{
		if (m_plot == 0) m_plot = new LSDYNAPlotFile;

		if (m_plot->Open(*this, m_szplot) == false)
		{
			log.printf("ERROR : Failed creating PLOT database\n");
			return false;
		}
	}

	// Since it is assumed that for the first timestep
	// there are no loads or initial displacements, the case n=0 is skipped.
	// Therefor we can output those results here.
	// Offcourse we should actually check if this is indeed
	// the case, otherwise we should also solve for t=0
	if (m_pStep->m_nplot != FE_PLOT_NEVER) m_plot->Write(*this);
/*
	// reset the log file
	if (!log.is_valid())
	{
		log.open(m_szlog);

		// if we don't want to output anything we only output to the logfile
		if (m_pStep->GetPrintLevel() == FE_PRINT_NEVER) log.SetMode(Logfile::FILE_ONLY);

		// print welcome message to file
		Hello(log);
	}
*/
	// do the callback
	DoCallback();

	// All data is reset successfully
	return true;
}
