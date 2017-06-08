/*!
 * \page Vortex_in_cell_petsc Vortex in Cell 3D
 *
 * # Vortex in Cell 3D ring with ring-lets # {#vic_ringlets}
 *
 * In this example we solve the Navier-Stokes in the vortex formulation in 3D
 * for an incompressible fluid
 *
 *
 */

#define SE_CLASS1
//#define PRINT_STACKTRACE

//! \cond [def system] \endcond

#include "Grid/grid_dist_id.hpp"
#include "Vector/vector_dist.hpp"
#include "Matrix/SparseMatrix.hpp"
#include "Vector/Vector.hpp"


#include "FiniteDifference/FDScheme.hpp"
#include "FiniteDifference/util/common.hpp"
#include "FiniteDifference/eq.hpp"
#include "Solvers/petsc_solver.hpp"
#include "Solvers/petsc_solver.hpp"
#include "Solvers/umfpack_solver.hpp"
#include "interpolation/mp4_kernel.hpp"
#include "interpolation/interpolation.hpp"

constexpr int x = 0;
constexpr int y = 1;
constexpr int z = 2;

////////////////// Simulation Parameters //////////////////////

typedef grid_dist_id<3,float,aggregate<float[3]>> grid_type;
typedef vector_dist<3,float,aggregate<float[3],float[3],float[3],float[3],float[3]>> particles_type;

// radius of the torus
float ringr1 = 5.0/4.0;
// radius of the core of the torus
float ringr2 = 0.5;
// Reinold number
float tgtre  = 10000.0;
// Noise factor for the ring vorticity on z
float ringnz = 0.00009;

// Kinematic viscosity
float nu = 0.0001535;
float sc = 2.0;

float dt = 0.012;


constexpr unsigned int vorticity = 0;
constexpr unsigned int velocity = 0;
constexpr unsigned int p_vel = 1;
constexpr int rhs_part = 2;
constexpr unsigned int old_vort = 3;
constexpr unsigned int old_pos = 4;

//! \cond [bond def eq] \endcond

#include "Vector/vector_dist.hpp"
#include "data_type/aggregate.hpp"

void init_ring(grid_type & gr, const Box<3,float> & domain)
{
	constexpr int nk = 32;

	float ak[nk];
	float bk[nk];

	for (size_t i = 0 ; i < nk ; i++)
	{
	     ak[i] = rand()/RAND_MAX;
	     bk[i] = rand()/RAND_MAX;
	}

	float gamma = nu * tgtre;
	float rinv2 = 1.0f/(ringr2*ringr2);
	float max_vorticity = gamma*rinv2/M_PI;

	auto it = gr.getDomainIterator();

	while (it.isNext())
	{
		auto key_d = it.get();
		auto key = it.getGKey(key_d);

        float tx = (key.get(x)-2)*gr.spacing(x) + domain.getLow(x);
        float ty = (key.get(y)-2)*gr.spacing(y) + domain.getLow(y);
        float tz = (key.get(z)-2)*gr.spacing(z) + domain.getLow(z);
        float theta1 = atan2((ty-0.5f),(tz-0.5f));

        float noise = 0.0f;
        for (int kk=1 ; kk < nk; kk++)
        	noise = noise + sin(kk*(theta1+2.0f*M_PI*ak[kk])) + cos(kk*(theta1+2.0f*M_PI*bk[kk]));

        float rad1r  = sqrt((ty-2.5f)*(ty-2.5f) + (tz-2.5f)*(tz-2.5f)) - ringr1*(1.0f + ringnz * noise);
        float rad1t = tx - 2.5f;
        float rad1sq = rad1r*rad1r + rad1t*rad1t;
        float radstr = -exp(-rad1sq*rinv2)*rinv2*gamma/M_PI;
        gr.template get<vorticity>(key_d)[x] = 0.0f;
        gr.template get<vorticity>(key_d)[y] = -radstr * cos(theta1);
        gr.template get<vorticity>(key_d)[z] = radstr * sin(theta1);

        theta1 = atan2((ty-0.5f),(tz-0.5f));
        rad1r  = sqrt((ty-0.5f)*(ty-0.5f) + (tz-0.5f)*(tz-0.5f)) + ringr1*(1.0f + ringnz * noise);
        rad1t = tx - 0.5f;
        rad1sq = rad1r*rad1r + rad1t*rad1t;
        float rad1sqTILDA = rad1sq*rinv2;
        radstr = exp(-rad1sqTILDA)*rinv2*gamma/M_PI;
        gr.template get<vorticity>(key_d)[x] = 0.0f;
        gr.template get<vorticity>(key_d)[y] = gr.template get<vorticity>(key_d)[y] + radstr * cos(theta1);
        gr.template get<vorticity>(key_d)[z] = gr.template get<vorticity>(key_d)[z] - radstr * sin(theta1);

		++it;
	}
}

// Specification of the poisson equation for the helmotz-hodge projection
// 3D (dims = 3). The field is a scalar value (nvar = 1), bournary are periodic
// type of the the space is float. Final grid that will store \phi (grid_dist_id<.....>)
// The other indicate to construct a Parallel Matrix and a parallel vector on PETSC
// NORMAL_GRID indicate that is a standard grid (non-staggered)
struct poisson_nn_helm
{
        static const unsigned int dims = 3;
        static const unsigned int nvar = 1;
        static const bool boundary[];
        typedef float stype;
        typedef grid_dist_id<3,float,aggregate<float>> b_grid;
        typedef SparseMatrix<double,int,PETSC_BASE> SparseMatrix_type;
        typedef Vector<double,PETSC_BASE> Vector_type;
        static const int grid_type = NORMAL_GRID;
};

const bool poisson_nn_helm::boundary[] = {PERIODIC,PERIODIC,PERIODIC};


void helmotz_hodge_projection(grid_type & gr, const Box<3,float> & domain)
{
	///////////////////////////////////////////////////////////////

	// Convenient constants
	constexpr unsigned int phi = 0;

	// Constructing the poisson equation
	typedef Field<phi,poisson_nn_helm> phi_f;

	// Poisson equation
	typedef Lap<phi_f,poisson_nn_helm> poisson;

	poisson ps;

	// ghost get
	gr.template ghost_get<vorticity>();

    Ghost<3,long int> g(1);

	// Here we create a distributed grid to store the result of the helmotz projection
	grid_dist_id<3,float,aggregate<float>> gr_hl(gr.getDecomposition(),gr.getGridInfo().getSize(),g);

	// Calculate the divergence of the vortex field

	auto it = gr.getDomainIterator();	

	while (it.isNext())
	{
		auto key = it.get();

        gr_hl.template get<phi>(key) = 1.0f/gr.spacing(x)/2.0f*(gr.template get<vorticity>(key.move(x,1))[x] - gr.template get<vorticity>(key.move(x,-1))[x] ) +
			         	 	 1.0f/gr.spacing(y)/2.0f*(gr.template get<vorticity>(key.move(y,1))[y] - gr.template get<vorticity>(key.move(y,-1))[y] ) +
							 1.0f/gr.spacing(z)/2.0f*(gr.template get<vorticity>(key.move(z,1))[z] - gr.template get<vorticity>(key.move(z,-1))[z] );

		++it;
	}

	Padding<3>pd({1,1,1},{1,1,1});
	Ghost<3,long int> stencil_max(1);

	// Finite difference scheme
	FDScheme<poisson_nn_helm> fd(pd, stencil_max, domain, gr_hl.getGridInfo(), gr_hl);

	// fd.impose(ic_eq(),0.0, EQ_3, {0,0},{sz[0]-2,sz[1]-2},true);

	fd.template impose_dit<phi>(ps,gr_hl,0,gr_hl.getDomainIterator());
	fd.template impose(phi_f(),0,0,{-1,-1,-1},{-1,(long int)gr_hl.size(1),(long int)gr_hl.size(2)});
	fd.template impose(phi_f(),0,0,{(long int)gr_hl.size(0),-1,-1},{(long int)gr_hl.size(0),(long int)gr_hl.size(1),(long int)gr_hl.size(2)});

	fd.template impose(phi_f(),0,0,{0,-1,-1},{(long int)gr_hl.size(0)-1,-1,(long int)gr_hl.size(2)});
	fd.template impose(phi_f(),0,0,{0,(long int)gr_hl.size(1),-1},{(long int)gr_hl.size(0)-1,(long int)gr_hl.size(1),(long int)gr_hl.size(2)});

	fd.template impose(phi_f(),0,0,{0,0,-1},{(long int)gr_hl.size(0)-1,(long int)gr_hl.size(1)-1,-1});
	fd.template impose(phi_f(),0,0,{0,0,(long int)gr_hl.size(2)},{(long int)gr_hl.size(0)-1,(long int)gr_hl.size(1)-1,(long int)gr_hl.size(2)});

	// Create an PETSC solver
	petsc_solver<double> solver;

	solver.setSolver(KSPCG);
	solver.setAbsTol(0.0001);
	solver.setRestart(250);
	solver.setMaxIter(500);

	auto & A = fd.getA();
	auto & trpl = A.getMatrixTriplets();

	for (size_t i = 0 ; i < trpl.size() ; i++)
	{
		std::cout << "Triplet: " << trpl.get(i).row() << "   " << trpl.get(i).col() << "   " << trpl.get(i).value() << std::endl;
	}

	// Give to the solver A and b, return x, the solution
	auto phi_s = solver.solve(fd.getA(),fd.getB());

	// Create an UMFPACK solver
/*	umfpack_solver<double> solver;

	// Give to the solver A and b, return x, the solution
	auto phi_s = solver.solve(fd.getA(),fd.getB());*/

	long int szs[3] = {(long int)gr.getGridInfo().size(0),
			           (long int)gr.getGridInfo().size(1),
					   (long int)gr.getGridInfo().size(2)};

	// copy the solution to grid
	fd.template copy<phi>(phi_s,{0,0,0},szs,gr_hl);

	gr_hl.template ghost_get<phi>();

	// Correct the vorticity to make it divergence free

	auto it2 = gr.getDomainIterator();

	while (it2.isNext())
	{
		auto key = it2.get();

		gr.template get<vorticity>(key)[x] -= 1.0f/2.0f/gr.spacing(x)*(gr_hl.template get<phi>(key.move(x,1))-gr_hl.template get<phi>(key.move(x,-1)));
		gr.template get<vorticity>(key)[y] -= 1.0f/2.0f/gr.spacing(y)*(gr_hl.template get<phi>(key.move(y,1))-gr_hl.template get<phi>(key.move(y,-1)));
		gr.template get<vorticity>(key)[z] -= 1.0f/2.0f/gr.spacing(z)*(gr_hl.template get<phi>(key.move(z,1))-gr_hl.template get<phi>(key.move(z,-1)));

		++it2;
	}


	// Here we check the maximum divergence of the ring
	// We also copy the solution to the original grid
	double max = 0.0;

	gr.template ghost_get<vorticity>();
	auto it3 = gr.getDomainIterator();

	while (it3.isNext())
	{
		auto key = it3.get();

        gr_hl.template get<phi>(key) = 1.0f/gr.spacing(x)/2.0f*(gr.template get<vorticity>(key.move(x,1))[x] - gr.template get<vorticity>(key.move(x,-1))[x] ) +
			         	 	 1.0f/gr.spacing(y)/2.0f*(gr.template get<vorticity>(key.move(y,1))[y] - gr.template get<vorticity>(key.move(y,-1))[y] ) +
							 1.0f/gr.spacing(z)/2.0f*(gr.template get<vorticity>(key.move(z,1))[z] - gr.template get<vorticity>(key.move(z,-1))[z] );

        if (gr_hl.template get<phi>(key) > max)
        	max = gr_hl.template get<phi>(key);

		++it3;
	}

	std::cout << "Maximum divergence of the ring MAX " << max << std::endl;

}

/*! \brief Here we calculate the the velocity from vortex on grid
 *
 *
 */
void comp_vel(Box<3,float> & domain, grid_type & g_vort,grid_type & g_vel, petsc_solver<double>::return_type (& phi_s)[3])
{
	// Convenient constants
	constexpr unsigned int phi = 0;

	// Constructing the poisson equation
	typedef Field<phi,poisson_nn_helm> phi_f;

	// Poisson equation
	typedef Lap<phi_f,poisson_nn_helm> poisson;

	Padding<3>pd({1,1,1},{1,1,1});
	Ghost<3,long int> stencil_max(1);

	Ghost<3,long int> g(1);

	long int szs[3] = {(long int)g_vort.getGridInfo().size(0),
			           (long int)g_vort.getGridInfo().size(1),
					   (long int)g_vort.getGridInfo().size(2)};

	float inv_dx = 0.5f / g_vort.spacing(0);
	float inv_dy = 0.5f / g_vort.spacing(1);
	float inv_dz = 0.5f / g_vort.spacing(2);

	// Here we create a distributed grid to store the result of the helmotz projection
	grid_dist_id<3,float,aggregate<float>> gr_ps(g_vort.getDecomposition(),g_vort.getGridInfo().getSize(),g);
	grid_dist_id<3,float,aggregate<float[3]>> gr_ps_v(g_vort.getDecomposition(),g_vort.getGridInfo().getSize(),g);

	// For each component solve a poisson
	for (size_t i = 0 ; i < 3 ; i++)
	{
		// Copy the vorticity in gr_ps
		auto it2 = gr_ps.getDomainIterator();

		// calculate the velocity from the curl of phi
		while (it2.isNext())
		{
			auto key = it2.get();

			gr_ps.get<vorticity>(key) = g_vort.template get<vorticity>(key)[i];

			++it2;
		}


		// Finite difference scheme
		FDScheme<poisson_nn_helm> fd(pd, stencil_max, domain, gr_ps.getGridInfo(), gr_ps);

		poisson ps;

		fd.template impose_dit<phi>(ps,gr_ps,0,gr_ps.getDomainIterator());

		fd.template impose(phi_f(),0,0,{-1,-1,-1},{-1,(long int)gr_ps.size(1),(long int)gr_ps.size(2)});
		fd.template impose(phi_f(),0,0,{(long int)gr_ps.size(0),-1,-1},{(long int)gr_ps.size(0),(long int)gr_ps.size(1),(long int)gr_ps.size(2)});

		fd.template impose(phi_f(),0,0,{0,-1,-1},{(long int)gr_ps.size(0)-1,-1,(long int)gr_ps.size(2)});
		fd.template impose(phi_f(),0,0,{0,(long int)gr_ps.size(1),-1},{(long int)gr_ps.size(0)-1,(long int)gr_ps.size(1),(long int)gr_ps.size(2)});

		fd.template impose(phi_f(),0,0,{0,0,-1},{(long int)gr_ps.size(0)-1,(long int)gr_ps.size(1)-1,-1});
		fd.template impose(phi_f(),0,0,{0,0,(long int)gr_ps.size(2)},{(long int)gr_ps.size(0)-1,(long int)gr_ps.size(1)-1,(long int)gr_ps.size(2)});

		// Create an PETSC solver
		petsc_solver<double> solver;

		solver.setSolver(KSPCG);
		solver.setAbsTol(0.0001);
		solver.setRestart(250);
		solver.setMaxIter(500);

		std::cout << "Solving component " << i << std::endl;

		PetscBool flg;
		SparseMatrix<double,int,PETSC_BASE> & A = fd.getA();
//		MatIsSymmetric(A.getMat(),0.000001,&flg);

		Vector<double,PETSC_BASE> & b = fd.getB();

		// Give to the solver A and b, return x, the solution
		solver.solve(A,phi_s[i],b);

		// Calculate the residual

		petsc_solver<double>::return_type r(gr_ps.size(),gr_ps.getLocalDomainSize());
		Vec & pr = r.getVec();

		PETSC_SAFE_CALL(MatResidual(A.getMat(),b.getVec(),phi_s[i].getVec(),pr));

		PetscReal ne;
		PETSC_SAFE_CALL(VecNorm(pr,NORM_INFINITY,&ne));

		// Create an UMFPACK solver
/*		umfpack_solver<double> solver;

		// Give to the solver A and b, return x, the solution
		auto x = solver.solve(fd.getA(),fd.getB());*/

		std::cout << "Solved component " << i << "  Error: " << ne << "   Symmetric: " << flg << std::endl;

//		std::cout << "Solved" << std::endl;

		// copy the solution to grid
		fd.template copy<phi>(phi_s[i],{0,0,0},szs,gr_ps);

		auto it3 = gr_ps.getDomainIterator();

		// calculate the velocity from the curl of phi
		while (it3.isNext())
		{
			auto key = it3.get();

			gr_ps_v.get<velocity>(key)[i] = gr_ps.get<phi>(key);

			++it3;
		}
	}

	gr_ps_v.ghost_get<phi>();

	auto it3 = gr_ps_v.getDomainIterator();

	// calculate the velocity from the curl of phi
	while (it3.isNext())
	{
		auto key = it3.get();

		float phi_xy = (gr_ps_v.get<phi>(key.move(y,1))[x] - gr_ps_v.get<phi>(key.move(y,-1))[x])*inv_dy;
		float phi_xz = (gr_ps_v.get<phi>(key.move(z,1))[x] - gr_ps_v.get<phi>(key.move(z,-1))[x])*inv_dz;
		float phi_yx = (gr_ps_v.get<phi>(key.move(x,1))[y] - gr_ps_v.get<phi>(key.move(x,-1))[y])*inv_dx;
		float phi_yz = (gr_ps_v.get<phi>(key.move(z,1))[y] - gr_ps_v.get<phi>(key.move(z,-1))[y])*inv_dz;
		float phi_zx = (gr_ps_v.get<phi>(key.move(x,1))[z] - gr_ps_v.get<phi>(key.move(x,-1))[z])*inv_dx;
		float phi_zy = (gr_ps_v.get<phi>(key.move(y,1))[z] - gr_ps_v.get<phi>(key.move(y,-1))[z])*inv_dy;

		g_vel.template get<velocity>(key)[x] = phi_zy - phi_yz;
		g_vel.template get<velocity>(key)[y] = phi_xz - phi_zx;
		g_vel.template get<velocity>(key)[z] = phi_yx - phi_xy;

		++it3;
	}

	g_vel.template ghost_get<velocity>();

	// We check that curl u = vorticity

	auto it4 = gr_ps_v.getDomainIterator();

	double norm_max = 0.0;

	// calculate the velocity from the curl of phi
	while (it4.isNext())
	{
		auto key = it4.get();

		float phi_xy = (g_vel.template get<phi>(key.move(y,1))[x] - g_vel.template get<phi>(key.move(y,-1))[x])*inv_dy;
		float phi_xz = (g_vel.template get<phi>(key.move(z,1))[x] - g_vel.template get<phi>(key.move(z,-1))[x])*inv_dz;
		float phi_yx = (g_vel.template get<phi>(key.move(x,1))[y] - g_vel.template get<phi>(key.move(x,-1))[y])*inv_dx;
		float phi_yz = (g_vel.template get<phi>(key.move(z,1))[y] - g_vel.template get<phi>(key.move(z,-1))[y])*inv_dz;
		float phi_zx = (g_vel.template get<phi>(key.move(x,1))[z] - g_vel.template get<phi>(key.move(x,-1))[z])*inv_dx;
		float phi_zy = (g_vel.template get<phi>(key.move(y,1))[z] - g_vel.template get<phi>(key.move(y,-1))[z])*inv_dy;

		gr_ps_v.template get<velocity>(key)[x] = phi_zy - phi_yz  + g_vort.template get<vorticity>(key)[x];
		gr_ps_v.template get<velocity>(key)[y] = phi_xz - phi_zx + g_vort.template get<vorticity>(key)[y];
		gr_ps_v.template get<velocity>(key)[z] = phi_yx - phi_xy + g_vort.template get<vorticity>(key)[z];

		if (gr_ps_v.template get<velocity>(key)[x] > norm_max)
			norm_max = gr_ps_v.template get<velocity>(key)[x];

		if (gr_ps_v.template get<velocity>(key)[y] > norm_max)
			norm_max = gr_ps_v.template get<velocity>(key)[y];

		if (gr_ps_v.template get<velocity>(key)[z] > norm_max)
			norm_max = gr_ps_v.template get<velocity>(key)[z];

		++it4;
	}

	std::cout << "Norm max: " << norm_max << std::endl;
}

void remesh(particles_type & vd, grid_type & gr,Box<3,float> & domain)
{
	// Remove all particles
	vd.clear();

	auto git = vd.getGridIterator(gr.getGridInfo().getSize());

	while (git.isNext())
	{
		auto p = git.get();
		auto key = git.get_dist();

		vd.add();

		vd.getLastPos()[x] = gr.spacing(x)*p.get(x) + domain.getLow(x);
		vd.getLastPos()[y] = gr.spacing(y)*p.get(y) + domain.getLow(y);
		vd.getLastPos()[z] = gr.spacing(z)*p.get(z) + domain.getLow(z);

		vd.template getLastProp<vorticity>()[x] = gr.template get<vorticity>(key)[x];
		vd.template getLastProp<vorticity>()[y] = gr.template get<vorticity>(key)[y];
		vd.template getLastProp<vorticity>()[z] = gr.template get<vorticity>(key)[z];

		++git;
	}
}

void remesh(particles_type & vd, grid_type & g_rhs, grid_type & g_vel,Box<3,float> & domain)
{
	constexpr int rhs_g = 0;

	// Remove all particles
	vd.clear();

	auto git = vd.getGridIterator(g_rhs.getGridInfo().getSize());

	while (git.isNext())
	{
		auto p = git.get();
		auto key = git.get_dist();

		vd.add();

		vd.getLastPos()[x] = g_rhs.spacing(x)*p.get(x) + domain.getLow(x);
		vd.getLastPos()[y] = g_rhs.spacing(y)*p.get(y) + domain.getLow(y);
		vd.getLastPos()[z] = g_rhs.spacing(z)*p.get(z) + domain.getLow(z);

		vd.template getLastProp<rhs_part>()[x] = g_rhs.template get<rhs_g>(key)[x];
		vd.template getLastProp<rhs_part>()[y] = g_rhs.template get<rhs_g>(key)[y];
		vd.template getLastProp<rhs_part>()[z] = g_rhs.template get<rhs_g>(key)[z];

		vd.template getLastProp<p_vel>()[x] = g_vel.template get<velocity>(key)[x];
		vd.template getLastProp<p_vel>()[y] = g_vel.template get<velocity>(key)[y];
		vd.template getLastProp<p_vel>()[z] = g_vel.template get<velocity>(key)[z];

		++git;
	}
}

template<unsigned int prp> void set_zero(grid_type & gr)
{
	auto it = gr.getDomainGhostIterator();

	// calculate the velocity from the curl of phi
	while (it.isNext())
	{
		auto key = it.get();

		gr.template get<prp>(key)[0] = 0.0;
		gr.template get<prp>(key)[1] = 0.0;
		gr.template get<prp>(key)[2] = 0.0;

		++it;
	}
}


template<unsigned int prp> void set_zero(particles_type & vd)
{
	auto it = vd.getDomainIterator();

	// calculate the velocity from the curl of phi
	while (it.isNext())
	{
		auto key = it.get();

		vd.template getProp<prp>(key)[0] = 0.0;
		vd.template getProp<prp>(key)[1] = 0.0;
		vd.template getProp<prp>(key)[2] = 0.0;

		++it;
	}
}

//////////////////////////// LAPLACIAN 3D ///////////////////////////////////////////


static char help[] = "Solves 3D Laplacian using multigrid.\n\n";

#include <petscdm.h>
#include <petscdmda.h>
#include <petscksp.h>

extern PetscErrorCode ComputeMatrix(KSP,Mat,Mat,void*);
extern PetscErrorCode ComputeRHS(KSP,Vec,void*);

int solve3D()
{
   KSP            ksp;
   DM             da;
   PetscReal      norm;

   PetscInt    i,j,k,mx,my,mz,xm,ym,zm,xs,ys,zs;
   PetscScalar Hx,Hy,Hz;
   PetscScalar ***array;
   Vec         x,b,r;
   Mat         J;

//   PetscInitialize(&argc,&argv,(char*)0,help);

   KSPCreate(PETSC_COMM_WORLD,&ksp);
   DMDACreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_STAR,256,128,128,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,1,1,0,0,0,&da);
   DMDASetInterpolationType(da, DMDA_Q0);

   KSPSetDM(ksp,da);

   KSPSetComputeRHS(ksp,ComputeRHS,NULL);
   KSPSetComputeOperators(ksp,ComputeMatrix,NULL);
   KSPSetFromOptions(ksp);
   KSPSolve(ksp,NULL,NULL);
   KSPGetSolution(ksp,&x);
   KSPGetRhs(ksp,&b);
   KSPGetOperators(ksp,NULL,&J);
   VecDuplicate(b,&r);

   MatMult(J,x,r);
   VecAXPY(r,-1.0,b);
   VecNorm(r,NORM_2,&norm);
   PetscPrintf(PETSC_COMM_WORLD,"Residual norm %g\n",(double)norm);

   DMDAGetInfo(da, 0, &mx, &my, &mz, 0,0,0,0,0,0,0,0,0);
   Hx   = 1.0 / (PetscReal)(mx);
   Hy   = 1.0 / (PetscReal)(my);
   Hz   = 1.0 / (PetscReal)(mz);
   DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);
   DMDAVecGetArray(da, x, &array);

   for (k=zs; k<zs+zm; k++) {
     for (j=ys; j<ys+ym; j++) {
       for (i=xs; i<xs+xm; i++) {
         array[k][j][i] -=
           PetscCosScalar(2*PETSC_PI*(((PetscReal)i+0.5)*Hx))*
           PetscCosScalar(2*PETSC_PI*(((PetscReal)j+0.5)*Hy))*
           PetscCosScalar(2*PETSC_PI*(((PetscReal)k+0.5)*Hz));
       }
     }
   }
   DMDAVecRestoreArray(da, x, &array);
   VecAssemblyBegin(x);
   VecAssemblyEnd(x);

   VecNorm(x,NORM_INFINITY,&norm);
   PetscPrintf(PETSC_COMM_WORLD,"Error norm %g\n",(double)norm);
   VecNorm(x,NORM_1,&norm);
   PetscPrintf(PETSC_COMM_WORLD,"Error norm %g\n",(double)(norm/((PetscReal)(mx)*(PetscReal)(my)*(PetscReal)(mz))));
   VecNorm(x,NORM_2,&norm);
   PetscPrintf(PETSC_COMM_WORLD,"Error norm %g\n",(double)(norm/((PetscReal)(mx)*(PetscReal)(my)*(PetscReal)(mz))));

   VecDestroy(&r);
   KSPDestroy(&ksp);
   DMDestroy(&da);
   PetscFinalize();
   return 0;
 }

PetscErrorCode ComputeRHS(KSP ksp,Vec b,void *ctx)
{
   PetscInt       i,j,k,mx,my,mz,xm,ym,zm,xs,ys,zs;
   PetscScalar    Hx,Hy,Hz;
   PetscScalar    ***array;
   DM             da;
   MatNullSpace   nullspace;

   KSPGetDM(ksp,&da);
   DMDAGetInfo(da, 0, &mx, &my, &mz, 0,0,0,0,0,0,0,0,0);
   Hx   = 1.0 / (PetscReal)(mx);
   Hy   = 1.0 / (PetscReal)(my);
   Hz   = 1.0 / (PetscReal)(mz);
   DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);
   DMDAVecGetArray(da, b, &array);
   for (k=zs; k<zs+zm; k++) {
     for (j=ys; j<ys+ym; j++) {
       for (i=xs; i<xs+xm; i++) {
         array[k][j][i] = 12 * PETSC_PI * PETSC_PI
                          * PetscCosScalar(2*PETSC_PI*(((PetscReal)i+0.5)*Hx))
                          * PetscCosScalar(2*PETSC_PI*(((PetscReal)j+0.5)*Hy))
                          * PetscCosScalar(2*PETSC_PI*(((PetscReal)k+0.5)*Hz))
                          * Hx * Hy * Hz;
       }
     }
   }
   DMDAVecRestoreArray(da, b, &array);
   VecAssemblyBegin(b);
   VecAssemblyEnd(b);

   /* force right hand side to be consistent for singular matrix */
   /* note this is really a hack, normally the model would provide you with a consistent right handside */

   MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,0,&nullspace);
   MatNullSpaceRemove(nullspace,b);
   MatNullSpaceDestroy(&nullspace);
   return(0);
}


 PetscErrorCode ComputeMatrix(KSP ksp, Mat J,Mat jac, void *ctx)
 {
   PetscInt       i,j,k,mx,my,mz,xm,ym,zm,xs,ys,zs,num, numi, numj, numk;
   PetscScalar    v[7],Hx,Hy,Hz,HyHzdHx,HxHzdHy,HxHydHz;
   MatStencil     row, col[7];
   DM             da;
   MatNullSpace   nullspace;

   KSPGetDM(ksp,&da);
   DMDAGetInfo(da,0,&mx,&my,&mz,0,0,0,0,0,0,0,0,0);
   Hx      = 1.0 / (PetscReal)(mx);
   Hy      = 1.0 / (PetscReal)(my);
   Hz      = 1.0 / (PetscReal)(mz);
   HyHzdHx = Hy*Hz/Hx;
   HxHzdHy = Hx*Hz/Hy;
   HxHydHz = Hx*Hy/Hz;
   DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);
   for (k=zs; k<zs+zm; k++) {
     for (j=ys; j<ys+ym; j++) {
       for (i=xs; i<xs+xm; i++) {
         row.i = i; row.j = j; row.k = k;
         if (i==0 || j==0 || k==0 || i==mx-1 || j==my-1 || k==mz-1) {
           num = 0; numi=0; numj=0; numk=0;
           if (k!=0) {
             v[num]     = -HxHydHz;
             col[num].i = i;
             col[num].j = j;
             col[num].k = k-1;
             num++; numk++;
           }
           if (j!=0) {
             v[num]     = -HxHzdHy;
             col[num].i = i;
             col[num].j = j-1;
             col[num].k = k;
             num++; numj++;
             }
           if (i!=0) {
             v[num]     = -HyHzdHx;
             col[num].i = i-1;
             col[num].j = j;
             col[num].k = k;
             num++; numi++;
           }
           if (i!=mx-1) {
             v[num]     = -HyHzdHx;
             col[num].i = i+1;
             col[num].j = j;
             col[num].k = k;
             num++; numi++;
           }
           if (j!=my-1) {
             v[num]     = -HxHzdHy;
             col[num].i = i;
             col[num].j = j+1;
             col[num].k = k;
             num++; numj++;
           }
           if (k!=mz-1) {
             v[num]     = -HxHydHz;
             col[num].i = i;
             col[num].j = j;
             col[num].k = k+1;
             num++; numk++;
           }
           v[num]     = (PetscReal)(numk)*HxHydHz + (PetscReal)(numj)*HxHzdHy + (PetscReal)(numi)*HyHzdHx;
           col[num].i = i;   col[num].j = j;   col[num].k = k;
           num++;
           MatSetValuesStencil(jac,1,&row,num,col,v,INSERT_VALUES);
         } else {
           v[0] = -HxHydHz;                          col[0].i = i;   col[0].j = j;   col[0].k = k-1;
           v[1] = -HxHzdHy;                          col[1].i = i;   col[1].j = j-1; col[1].k = k;
           v[2] = -HyHzdHx;                          col[2].i = i-1; col[2].j = j;   col[2].k = k;
           v[3] = 2.0*(HyHzdHx + HxHzdHy + HxHydHz); col[3].i = i;   col[3].j = j;   col[3].k = k;
           v[4] = -HyHzdHx;                          col[4].i = i+1; col[4].j = j;   col[4].k = k;
           v[5] = -HxHzdHy;                          col[5].i = i;   col[5].j = j+1; col[5].k = k;
           v[6] = -HxHydHz;                          col[6].i = i;   col[6].j = j;   col[6].k = k+1;
           MatSetValuesStencil(jac,1,&row,7,col,v,INSERT_VALUES);
         }
       }
     }
   }
   MatAssemblyBegin(jac,MAT_FINAL_ASSEMBLY);
   MatAssemblyEnd(jac,MAT_FINAL_ASSEMBLY);
   MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,0,&nullspace);
   MatSetNullSpace(jac,nullspace);
   MatNullSpaceDestroy(&nullspace);
   return(0);
 }

/////////////////////////////////////////////////////////////////////////////////////


// Calculate the right hand side of the vorticity formulation
template<typename grid> void calc_rhs(grid & g_vort, grid & g_vel, grid & g_dwp)
{
	constexpr int rhs = 0;

	float fac1 = 2.0f*nu/(g_vort.spacing(0));
	float fac2 = 2.0f*nu/(g_vort.spacing(1));
	float fac3 = 2.0f*nu/(g_vort.spacing(2));

	float fac4 = 0.5f/(g_vort.spacing(0));
	float fac5 = 0.5f/(g_vort.spacing(1));
	float fac6 = 0.5f/(g_vort.spacing(2));

	auto it = g_dwp.getDomainIterator();

	g_vort.template ghost_get<vorticity>();

	while (it.isNext())
	{
		auto key = it.get();

		g_dwp.template get<rhs>(key)[x] = /*fac1*(g_vort.template get<vorticity>(key.move(x,1))[x]+g_vort.template get<vorticity>(key.move(x,-1))[x])+
				                          fac2*(g_vort.template get<vorticity>(key.move(y,1))[x]+g_vort.template get<vorticity>(key.move(y,-1))[x])+
										  fac3*(g_vort.template get<vorticity>(key.move(z,1))[x]+g_vort.template get<vorticity>(key.move(z,-1))[x])-
										  2.0f*(fac1+fac2+fac3)*g_vort.template get<vorticity>(key)[x]+*/
										  fac4*g_vort.template get<vorticity>(key)[x]*
										  (g_vel.template get<velocity>(key.move(x,1))[x] - g_vel.template get<velocity>(key.move(x,-1))[x])+
										  fac5*g_vort.template get<vorticity>(key)[y]*
										  (g_vel.template get<velocity>(key.move(y,1))[x] - g_vel.template get<velocity>(key.move(y,-1))[x])+
										  fac6*g_vort.template get<vorticity>(key)[z]*
										  (g_vel.template get<velocity>(key.move(z,1))[x] - g_vel.template get<velocity>(key.move(z,-1))[x]);

		g_dwp.template get<rhs>(key)[y] = /*fac1*(g_vort.template get<vorticity>(key.move(x,1))[y]+g_vort.template get<vorticity>(key.move(x,-1))[y])+
				                          fac2*(g_vort.template get<vorticity>(key.move(y,1))[y]+g_vort.template get<vorticity>(key.move(y,-1))[y])+
										  fac3*(g_vort.template get<vorticity>(key.move(z,1))[y]+g_vort.template get<vorticity>(key.move(z,-1))[y])-
										  2.0f*(fac1+fac2+fac3)*g_vort.template get<vorticity>(key)[y]+*/
										  fac4*g_vort.template get<vorticity>(key)[x]*
										  (g_vel.template get<velocity>(key.move(x,1))[y] - g_vel.template get<velocity>(key.move(x,-1))[y])+
										  fac5*g_vort.template get<vorticity>(key)[y]*
										  (g_vel.template get<velocity>(key.move(y,1))[y] - g_vel.template get<velocity>(key.move(y,-1))[y])+
										  fac6*g_vort.template get<vorticity>(key)[z]*
										  (g_vel.template get<velocity>(key.move(z,1))[y] - g_vel.template get<velocity>(key.move(z,-1))[y]);


		g_dwp.template get<rhs>(key)[z] = /*fac1*(g_vort.template get<vorticity>(key.move(x,1))[z]+g_vort.template get<vorticity>(key.move(x,-1))[z])+
				                          fac2*(g_vort.template get<vorticity>(key.move(y,1))[z]+g_vort.template get<vorticity>(key.move(y,-1))[z])+
										  fac3*(g_vort.template get<vorticity>(key.move(z,1))[z]+g_vort.template get<vorticity>(key.move(z,-1))[z])-
										  2.0f*(fac1+fac2+fac3)*g_vort.template get<vorticity>(key)[z]+*/
										  fac4*g_vort.template get<vorticity>(key)[x]*
										  (g_vel.template get<velocity>(key.move(x,1))[z] - g_vel.template get<velocity>(key.move(x,-1))[z])+
										  fac5*g_vort.template get<vorticity>(key)[y]*
										  (g_vel.template get<velocity>(key.move(y,1))[z] - g_vel.template get<velocity>(key.move(y,-1))[z])+
										  fac6*g_vort.template get<vorticity>(key)[z]*
										  (g_vel.template get<velocity>(key.move(z,1))[z] - g_vel.template get<velocity>(key.move(z,-1))[z]);

		++it;
	}
}


void rk_step1(particles_type & particles)
{
	auto it = particles.getDomainIterator();

	while (it.isNext())
	{
		auto key = it.get();

		// Save the old vorticity
		particles.getProp<old_vort>(key)[x] = particles.getProp<vorticity>(key)[x];
		particles.getProp<old_vort>(key)[y] = particles.getProp<vorticity>(key)[y];
		particles.getProp<old_vort>(key)[z] = particles.getProp<vorticity>(key)[z];

		particles.getProp<vorticity>(key)[x] = particles.getProp<vorticity>(key)[x] + 0.5f * dt * particles.getProp<rhs_part>(key)[x];
		particles.getProp<vorticity>(key)[y] = particles.getProp<vorticity>(key)[y] + 0.5f * dt * particles.getProp<rhs_part>(key)[y];
		particles.getProp<vorticity>(key)[z] = particles.getProp<vorticity>(key)[z] + 0.5f * dt * particles.getProp<rhs_part>(key)[z];

		++it;
	}

	auto it2 = particles.getDomainIterator();

	while (it2.isNext())
	{
		auto key = it2.get();

		// Save the old position
		particles.getProp<old_pos>(key)[x] = particles.getPos(key)[x];
		particles.getProp<old_pos>(key)[y] = particles.getPos(key)[y];
		particles.getProp<old_pos>(key)[z] = particles.getPos(key)[z];

		particles.getPos(key)[x] = particles.getPos(key)[x] + 0.5f * dt * particles.getProp<p_vel>(key)[x];
		particles.getPos(key)[y] = particles.getPos(key)[y] + 0.5f * dt * particles.getProp<p_vel>(key)[y];
		particles.getPos(key)[z] = particles.getPos(key)[z] + 0.5f * dt * particles.getProp<p_vel>(key)[z];

		++it2;
	}

	particles.map();
}

void rk_step2(particles_type & particles)
{
	auto it = particles.getDomainIterator();

	while (it.isNext())
	{
		auto key = it.get();

		particles.getProp<vorticity>(key)[x] = particles.getProp<old_vort>(key)[x] + dt * particles.getProp<rhs_part>(key)[x];
		particles.getProp<vorticity>(key)[y] = particles.getProp<old_vort>(key)[y] + dt * particles.getProp<rhs_part>(key)[y];
		particles.getProp<vorticity>(key)[z] = particles.getProp<old_vort>(key)[z] + dt * particles.getProp<rhs_part>(key)[z];

		++it;
	}

	auto it2 = particles.getDomainIterator();

	while (it2.isNext())
	{
		auto key = it2.get();

		particles.getPos(key)[x] = particles.getProp<old_pos>(key)[x] + dt * particles.getProp<p_vel>(key)[x];
		particles.getPos(key)[y] = particles.getProp<old_pos>(key)[y] + dt * particles.getProp<p_vel>(key)[y];
		particles.getPos(key)[z] = particles.getProp<old_pos>(key)[z] + dt * particles.getProp<p_vel>(key)[z];

		++it2;
	}

	particles.map();
}

template<typename grid, typename vector> void do_step(vector & particles,
		                                              grid & g_vort,
													  grid & g_vel,
													  grid & g_dvort,
													  Box<3,float> & domain,
													  interpolate<particles_type,grid_type,mp4_kernel<float>> & inte,
													  petsc_solver<double>::return_type (& phi_s)[3])
{
	constexpr int rhs = 0;

	set_zero<vorticity>(g_vort);
	inte.template p2m<vorticity,vorticity>(particles,g_vort);
	g_vort.template ghost_put<add_,vorticity>();

	// Calculate velocity from vorticity
	comp_vel(domain,g_vort,g_vel,phi_s);

	calc_rhs(g_vort,g_vel,g_dvort);

	g_dvort.template ghost_get<rhs>();
	g_vel.template ghost_get<velocity>();
	set_zero<rhs_part>(particles);
	set_zero<p_vel>(particles);
	inte.template m2p<rhs,rhs_part>(g_dvort,particles);
	inte.template m2p<velocity,p_vel>(g_vel,particles);
}

int main(int argc, char* argv[])
{
	// Initialize
	openfpm_init(&argc,&argv);

	// It store the solution to compute velocity
	// It is used as initial guess every time we call the solver
	petsc_solver<double>::return_type phi_s[3];

	// velocity in the grid is the property 0, pressure is the property 1
	constexpr int velocity = 0;
	constexpr int pressure = 1;

	// Domain, a rectangle
	Box<3,float> domain({0.0,0.0,0.0},{22.0,5.0,5.0});

	// Ghost (Not important in this case but required)
	Ghost<3,long int> g(2);

	// Grid points on x=128 y=64 z=64
	long int sz[] = {8,4,4};
	size_t szu[] = {(size_t)sz[0],(size_t)sz[1],(size_t)sz[2]};

	periodicity<3> bc = {{PERIODIC,PERIODIC,PERIODIC}};

	grid_type g_vort(szu,domain,g,bc);
	grid_type g_vel(g_vort.getDecomposition(),szu,g);
	grid_type g_dvort(g_vort.getDecomposition(),szu,g);
	particles_type particles(g_vort.getDecomposition(),0);

	std::cout << "SPACING " << g_vort.spacing(0) << " " << g_vort.spacing(1) << " " << g_vort.spacing(2) << std::endl;

	init_ring(g_vort,domain);
	helmotz_hodge_projection(g_vort,domain);

	remesh(particles,g_vort,domain);

	// Get the total number of particles
	Vcluster & v_cl = create_vcluster();

	size_t tot_part = particles.size_local();
	v_cl.sum(tot_part);
	v_cl.execute();

	// Now we set phi_s
	phi_s[0].resize(tot_part,particles.size_local());
	phi_s[1].resize(tot_part,particles.size_local());
	phi_s[2].resize(tot_part,particles.size_local());

	std::cout << "SETTING ZERO" << std::endl;

	phi_s[0].setZero();
	phi_s[1].setZero();
	phi_s[2].setZero();

	std::cout << "END setting zero" << std::endl;

	interpolate<particles_type,grid_type,mp4_kernel<float>> inte(particles,g_vort);

	// Time Integration

	for (size_t i = 0 ; i < 1 ; i++)
	{
		do_step(particles,g_vort,g_vel,g_dvort,domain,inte,phi_s);

		rk_step1(particles);

		do_step(particles,g_vort,g_vel,g_dvort,domain,inte,phi_s);

		rk_step2(particles);

		set_zero<vorticity>(g_vort);
		inte.template p2m<vorticity,vorticity>(particles,g_vort);
		g_vort.template ghost_put<add_,vorticity>();

		remesh(particles,g_vort,domain);

		std::cout << "Step " << i << std::endl;

		if (i % 10 == 0)
		{
			particles.write("part_out_" + std::to_string(i),VTK_WRITER | FORMAT_BINARY);

			g_vel.write("grid_velocity_" + std::to_string(i), VTK_WRITER | FORMAT_BINARY);
			g_vort.write("grid_vorticity_" + std::to_string(i), VTK_WRITER | FORMAT_BINARY);
		}

	}

	openfpm_finalize();
}

