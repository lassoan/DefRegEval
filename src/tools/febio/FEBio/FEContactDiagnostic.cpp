// FEContactDiagnostic.cpp: implementation of the FEContactDiagnostic class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FEContactDiagnostic.h"
#include "FENeoHookean.h"
#include "FESolidSolver.h"
#include "log.h"

void print_matrix(Logfile& log, FullMatrix& m)
{
	int i, j;
	int N = m.Size();
	int M = m.Size();

	log.printf("\n    ");
	for (i=0; i<N; ++i) log.printf("%15d ", i);
	log.printf("\n----");
	for (i=0; i<N; ++i) log.printf("----------------", i);

	for (i=0; i<N; ++i)
	{
		log.printf("\n%2d: ", i);
		for (j=0; j<M; ++j)
		{
			log.printf("%15lg ", m(i,j));
		}
	}
	log.printf("\n");
}


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

FEContactDiagnostic::FEContactDiagnostic(FEM& fem) : FEDiagnostic(fem)
{

}

FEContactDiagnostic::~FEContactDiagnostic()
{

}

bool FEContactDiagnostic::Run()
{
	// get the solver
	FESolidSolver& solver = dynamic_cast<FESolidSolver&>(*m_fem.m_pStep->m_psolver);
	solver.Init();

	// make sure contact data is up to data
	m_fem.UpdateContact();

	// create the stiffness matrix
	solver.CreateStiffness(true);

	// get the stiffness matrix
	FEStiffnessMatrix& K = *solver.m_pK;
	SparseMatrix *pA = (SparseMatrix*)K;
	FullMatrix& K0 = dynamic_cast<FullMatrix&>(*pA);

	// build the stiffness matrix
	K0.zero();
	solver.ContactStiffness();
//	solver.StiffnessMatrix();

	// get the logfile
	Logfile& log = GetLogfile();

	print_matrix(log, K0);

	// calculate the derivative of the residual
	FullMatrix K1;
	deriv_residual(K1);

	print_matrix(log, K1);

	// calculate difference matrix
	const int N = 48;
	FullMatrix Kd; Kd.Create(N);
	double kij, kmax = 0, k0;
	int i0=-1, j0=-1;
	for (int i=0; i<N; ++i)
		for (int j=0; j<N; ++j)
		{
			Kd(i,j) = K1(i,j) - K0(i,j);
			k0 = fabs(K0(i,j));
			if (k0 < 1e-15) k0 = 1;
			kij = 100.0*fabs(Kd(i,j))/k0;
			if (kij > kmax) 
			{
				kmax = kij;
				i0 = i;
				j0 = j;
			}
		}

	print_matrix(log, Kd);

	log.printf("\nMaximum difference: %lg%% (at (%d,%d))\n", kmax, i0, j0);

	return (kmax < 1e-3);
}

bool FEContactDiagnostic::Init()
{
	FEM& fem = m_fem;
	FEMesh& mesh = fem.m_mesh;

	// --- create the geometry ---

	// currently we simply assume a two-element contact problem
	// so we create two elements
	const double eps = 0.5;
	mesh.Create(16, 2, 0, 0);
	mesh.Node( 0).m_r0 = vec3d(0,0,0);
	mesh.Node( 1).m_r0 = vec3d(1,0,0);
	mesh.Node( 2).m_r0 = vec3d(1,1,0);
	mesh.Node( 3).m_r0 = vec3d(0,1,0);
	mesh.Node( 4).m_r0 = vec3d(0,0,1);
	mesh.Node( 5).m_r0 = vec3d(1,0,1);
	mesh.Node( 6).m_r0 = vec3d(1,1,1);
	mesh.Node( 7).m_r0 = vec3d(0,1,1);
	mesh.Node( 8).m_r0 = vec3d(0,0,1-eps);
	mesh.Node( 9).m_r0 = vec3d(1,0,1-eps);
	mesh.Node(10).m_r0 = vec3d(1,1,1-eps);
	mesh.Node(11).m_r0 = vec3d(0,1,1-eps);
	mesh.Node(12).m_r0 = vec3d(0,0,2-eps);
	mesh.Node(13).m_r0 = vec3d(1,0,2-eps);
	mesh.Node(14).m_r0 = vec3d(1,1,2-eps);
	mesh.Node(15).m_r0 = vec3d(0,1,2-eps);

	for (int i=0; i<16; ++i)
	{
		FENode& node = mesh.Node(i);
		node.m_rt = node.m_r0;

		// set rigid body id
		node.m_rid = -1;

		// open displacement dofs
		node.m_ID[0] = 0;
		node.m_ID[1] = 0;
		node.m_ID[2] = 0;

		// open rotational dofs
		node.m_ID[3] = 0;
		node.m_ID[4] = 0;
		node.m_ID[5] = 0;

		// open pressure dof
		node.m_ID[6] = 0;

		// close the rigid rotational dofs
		node.m_ID[7] = -1;
		node.m_ID[8] = -1;
		node.m_ID[9] = -1;

		node.m_ID[10] = 0;
	}

	FESolidElement& el0 = mesh.SolidElement(0);
	FESolidElement& el1 = mesh.SolidElement(1);

	el0.SetType(FE_HEX);
	el0.m_nID = 1;
	el0.SetMatID(-1);
	el0.m_node[0] = 0;
	el0.m_node[1] = 1;
	el0.m_node[2] = 2;
	el0.m_node[3] = 3;
	el0.m_node[4] = 4;
	el0.m_node[5] = 5;
	el0.m_node[6] = 6;
	el0.m_node[7] = 7;


	el1.SetType(FE_HEX);
	el0.m_nID = 2;
	el0.SetMatID(-1);
	el1.m_node[0] = 8;
	el1.m_node[1] = 9;
	el1.m_node[2] = 10;
	el1.m_node[3] = 11;
	el1.m_node[4] = 12;
	el1.m_node[5] = 13;
	el1.m_node[6] = 14;
	el1.m_node[7] = 15;

	// --- create a material ---
	FENeoHookean* pm = new FENeoHookean;
	pm->m_E = 1;
	pm->m_v = 0.45;
	fem.AddMaterial(pm);

	mesh.SetMatID(0);

	// --- create the sliding interface ---
	FESlidingInterface* ps = new FESlidingInterface(&fem);
	ps->m_atol = 0.1;
	ps->m_eps = 1;
	ps->m_npass = 1;
	ps->m_nsegup = 0;
	FESlidingSurface& ms = ps->m_ms;
	ms.Create(1);
	ms.Element(0).SetType(FE_NIQUAD);
	ms.Element(0).m_node[0] = 4;
	ms.Element(0).m_node[1] = 5;
	ms.Element(0).m_node[2] = 6;
	ms.Element(0).m_node[3] = 7;
	FESlidingSurface& ss = ps->m_ss;
	ss.Create(1);
	ss.Element(0).SetType(FE_NIQUAD);
	ss.Element(0).m_node[0] = 11;
	ss.Element(0).m_node[1] = 10;
	ss.Element(0).m_node[2] = 9;
	ss.Element(0).m_node[3] = 8;
	fem.m_CI.add(ps);

	// --- set fem data ---
	fem.m_nsolver = LU_SOLVER;	// make sure we have the LU solver

	return FEDiagnostic::Init();
}

void FEContactDiagnostic::deriv_residual(FullMatrix& K)
{
	FEM& fem = m_fem;

	// get the solver
	FESolidSolver& solver = dynamic_cast<FESolidSolver&>(*fem.m_pStep->m_psolver);

	// get the mesh
	FEMesh& mesh = fem.m_mesh;

	fem.UpdateContact();

	// first calculate the initial residual
	vector<double> R0(48); R0.zero();
	solver.ContactForces(R0);
//	solver.Residual(R0);

	// now calculate the perturbed residuals
	K.Create(48);
	int i, j, nj;
	int N = mesh.Nodes();
	double dx = 1e-8;
	vector<double> R1(48);
	for (j=0; j<3*N; ++j)
	{
		FENode& node = mesh.Node(j/3);
		nj = j%3;

		switch (nj)
		{
		case 0: node.m_rt.x += dx; break;
		case 1: node.m_rt.y += dx; break;
		case 2: node.m_rt.z += dx; break;
		}

		solver.UpdateStresses();
		fem.UpdateContact();

		R1.zero();
		solver.ContactForces(R1);
//		solver.Residual(R1);

		switch (nj)
		{
		case 0: node.m_rt.x -= dx; break;
		case 1: node.m_rt.y -= dx; break;
		case 2: node.m_rt.z -= dx; break;
		}

		solver.UpdateStresses();
		fem.UpdateContact();

		for (i=0; i<3*N; ++i) K(i,j) = (1-(R1[i] - R0[i])/dx)-1;
	}
}
