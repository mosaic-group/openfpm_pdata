/*! \page Vortex_in_cell_petsc_opt Vortex in Cell 3D (Optimization)
 *
 * # Vortex in Cell 3D ring with ringlets optimization # {#vic_ringlets_optimization}
 *
 * In this example we solve the Navier-Stokes equation in the vortex formulation in 3D
 * for an incompressible fluid. This example
 * has the following changes compared to \ref Vortex_in_cell_petsc
 *
 * * Constructing grid is expensive in particular with a lot of cores. For this
 *   reason we create the grid in the main function rather  than in the **comp_vel**
 *   function and **helmotz_hodge_projection**
 *
 *
 * * Constructing also FDScheme is expensive so we construct it once in the main. We set
 *   the left hand side to the poisson operator, and inside the functions **comp_vel**
 *   and **helmotz_hodge_projection** just write the right hand side with **impose_dit_b**
 *
 * \snippet Numerics/Vortex_in_cell/main_vic_petsc_opt.cpp construct grids
 *
 * \snippet Numerics/Vortex_in_cell/main_vic_petsc_opt.cpp create b
 *
 * Another optimization that we do is to use a Hilbert space-filling curve as sub-sub-domain
 * distribution strategy
 *
 */

#include "config.h"

#ifdef HAVE_PETSC

//#define SE_CLASS1
//#define PRINT_STACKTRACE

#include "interpolation/interpolation.hpp"
#include "Grid/grid_dist_id.hpp"
#include "Vector/vector_dist.hpp"
#include "Matrix/SparseMatrix.hpp"
#include "Vector/Vector.hpp"
#include "FiniteDifference/FDScheme.hpp"
#include "Solvers/petsc_solver.hpp"
#include "interpolation/mp4_kernel.hpp"
#include "Solvers/petsc_solver_AMG_report.hpp"
#include "Decomposition/Distribution/SpaceDistribution.hpp"

constexpr int x = 0;
constexpr int y = 1;
constexpr int z = 2;
constexpr unsigned int phi = 0;

// The type of the grids
typedef grid_dist_id<3,float,aggregate<float[3]>,CartDecomposition<3,float,HeapMemory,memory_traits_lin,SpaceDistribution<3,float>>> grid_type;

typedef grid_dist_id<3,float,aggregate<unsigned short int>,CartDecomposition<3,float,HeapMemory,memory_traits_lin,SpaceDistribution<3,float>>> grid_type_vis;

// The type of the grids
typedef grid_dist_id<3,float,aggregate<float>,CartDecomposition<3,float,HeapMemory,memory_traits_lin,SpaceDistribution<3,float>>> grid_type_s;

// The type of the particles
typedef vector_dist<3,
		            float,
		            aggregate<float[3],float[3],float[3],float[3],float[3]>,
		            CartDecomposition<3,float,HeapMemory,
		                              memory_traits_lin,SpaceDistribution<3,float>>,
		            HeapMemory,
		            memory_traits_lin> particles_type;

typedef vector_dist<3,
					float,
					aggregate<float>,
					CartDecomposition<3,float,HeapMemory,
					                  memory_traits_lin,SpaceDistribution<3,float>>,
					HeapMemory,
					memory_traits_lin> particles_type_s;

// radius of the torus
float ringr1 = 1.0;
// radius of the core of the torus
float sigma = 1.0/3.523;
// Reynold number (If you want to use 7500.0 you have to use a grid 1600x400x400)
//float tgtre  = 7500.0;
float tgtre = 3000.0;

// Kinematic viscosity
float nu = 1.0/tgtre;

// Time step
// float dt = 0.0025 for Re 7500
float dt = 0.0125;

#ifdef TEST_RUN
const unsigned int nsteps = 10;
#else
const unsigned int nsteps = 10001;
#endif

// All the properties by index
constexpr unsigned int vorticity = 0;
constexpr unsigned int velocity = 0;
constexpr unsigned int p_vel = 1;
constexpr int rhs_part = 2;
constexpr unsigned int old_vort = 3;
constexpr unsigned int old_pos = 4;

template<typename grid> void calc_and_print_max_div_and_int(grid & g_vort)
{
	g_vort.template ghost_get<vorticity>();
	auto it5 = g_vort.getDomainIterator();

	double max_vort = 0.0;

	double int_vort[3] = {0.0,0.0,0.0};

	while (it5.isNext())
	{
		auto key = it5.get();

		double tmp;

        tmp = 1.0f/g_vort.spacing(x)/2.0f*(g_vort.template get<vorticity>(key.move(x,1))[x] - g_vort.template get<vorticity>(key.move(x,-1))[x] ) +
			         	 	 1.0f/g_vort.spacing(y)/2.0f*(g_vort.template get<vorticity>(key.move(y,1))[y] - g_vort.template get<vorticity>(key.move(y,-1))[y] ) +
							 1.0f/g_vort.spacing(z)/2.0f*(g_vort.template get<vorticity>(key.move(z,1))[z] - g_vort.template get<vorticity>(key.move(z,-1))[z] );

        int_vort[x] += g_vort.template get<vorticity>(key)[x];
        int_vort[y] += g_vort.template get<vorticity>(key)[y];
        int_vort[z] += g_vort.template get<vorticity>(key)[z];

        if (tmp > max_vort)
        	max_vort = tmp;

		++it5;
	}

	Vcluster<> & v_cl = create_vcluster();
	v_cl.max(max_vort);
	v_cl.sum(int_vort[0]);
	v_cl.sum(int_vort[1]);
	v_cl.sum(int_vort[2]);
	v_cl.execute();

	if (v_cl.getProcessUnitID() == 0)
	{std::cout << "Max div for vorticity " << max_vort << "   Integral: " << int_vort[0] << "  " << int_vort[1] << "   " << int_vort[2] << std::endl;}
}

void init_ring(grid_type & gr, const Box<3,float> & domain)
{
	// To add some noise to the vortex ring we create two random
	// vector
	constexpr int nk = 32;

	float ak[nk];
	float bk[nk];

	for (size_t i = 0 ; i < nk ; i++)
	{
	     ak[i] = rand()/RAND_MAX;
	     bk[i] = rand()/RAND_MAX;
	}

	// We calculate the circuitation gamma
	float gamma = nu * tgtre;
	float rinv2 = 1.0f/(sigma*sigma);
	float max_vorticity = gamma*rinv2/M_PI;

	// We go through the grid to initialize the vortex
	auto it = gr.getDomainIterator();

	while (it.isNext())
	{
		auto key_d = it.get();
		auto key = it.getGKey(key_d);

        float tx = (key.get(x)-2)*gr.spacing(x) + domain.getLow(x);
        float ty = (key.get(y)-2)*gr.spacing(y) + domain.getLow(y);
        float tz = (key.get(z)-2)*gr.spacing(z) + domain.getLow(z);
        float theta1 = atan2((ty-domain.getHigh(1)/2.0f),(tz-domain.getHigh(2)/2.0f));


        float rad1r  = sqrt((ty-domain.getHigh(1)/2.0f)*(ty-domain.getHigh(1)/2.0f) + (tz-domain.getHigh(2)/2.0f)*(tz-domain.getHigh(2)/2.0f)) - ringr1;
        float rad1t = tx - 1.0f;
        float rad1sq = rad1r*rad1r + rad1t*rad1t;
        float radstr = -exp(-rad1sq*rinv2)*rinv2*gamma/M_PI;
        gr.template get<vorticity>(key_d)[x] = 0.0f;
        gr.template get<vorticity>(key_d)[y] = -radstr * cos(theta1);
        gr.template get<vorticity>(key_d)[z] = radstr * sin(theta1);

        // kill the axis term

        float rad1r_  = sqrt((ty-domain.getHigh(1)/2.0f)*(ty-domain.getHigh(1)/2.0f) + (tz-domain.getHigh(2)/2.0f)*(tz-domain.getHigh(2)/2.0f)) + ringr1;
        float rad1sqTILDA = rad1sq*rinv2;
        radstr = exp(-rad1sq*rinv2)*rinv2*gamma/M_PI;
        gr.template get<vorticity>(key_d)[x] = 0.0f;
        gr.template get<vorticity>(key_d)[y] = -radstr * cos(theta1);
        gr.template get<vorticity>(key_d)[z] = radstr * sin(theta1);

		++it;
	}
}

struct poisson_nn_helm
{
		//! 3D Stystem
        static const unsigned int dims = 3;
        //! We are solving for \psi that is a scalar
        static const unsigned int nvar = 1;
        //! Boundary conditions
        static const bool boundary[];
        //! type of the spatial coordinates
        typedef float stype;
        //! grid that store \psi
        typedef grid_type_s b_grid;
        //! Sparse matrix used to sove the linear system (we use PETSC)
        typedef SparseMatrix<double,int,PETSC_BASE> SparseMatrix_type;
        //! Vector to solve the system (PETSC)
        typedef Vector<double,PETSC_BASE> Vector_type;
        //! It is a normal grid
        static const int grid_type = NORMAL_GRID;
};

//! boundary conditions are PERIODIC
const bool poisson_nn_helm::boundary[] = {PERIODIC,PERIODIC,PERIODIC};


void helmotz_hodge_projection(grid_dist_id<3,float,aggregate<float>,CartDecomposition<3,float,HeapMemory,memory_traits_lin,SpaceDistribution<3,float>>> & psi,
							  FDScheme<poisson_nn_helm> & fd,
							  grid_type & gr,
							  const Box<3,float> & domain,
							  petsc_solver<double> & solver,
							  petsc_solver<double>::return_type & x_ ,
							  bool init)
{
	// ghost get
	gr.template ghost_get<vorticity>();

	// ghost size of the psi function
    Ghost<3,long int> g(2);

	// Calculate the divergence of the vortex field
	auto it = gr.getDomainIterator();	

	while (it.isNext())
	{
		auto key = it.get();

        psi.template get<phi>(key) = 1.0f/gr.spacing(x)/2.0f*(gr.template get<vorticity>(key.move(x,1))[x] - gr.template get<vorticity>(key.move(x,-1))[x] ) +
			         	 	 1.0f/gr.spacing(y)/2.0f*(gr.template get<vorticity>(key.move(y,1))[y] - gr.template get<vorticity>(key.move(y,-1))[y] ) +
							 1.0f/gr.spacing(z)/2.0f*(gr.template get<vorticity>(key.move(z,1))[z] - gr.template get<vorticity>(key.move(z,-1))[z] );

		++it;
	}

	calc_and_print_max_div_and_int(gr);

	//! \cond [create b] \endcond

	fd.new_b();
	fd.template impose_dit_b<0>(psi,psi.getDomainIterator());

	//! \cond [create b] \endcond

	timer tm_solve;
	if (init == true)
	{
		// Set the conjugate-gradient as method to solve the linear solver
		solver.setSolver(KSPBCGS);

		// Set the absolute tolerance to determine that the method is converged
		solver.setAbsTol(0.0001);

		// Set the maximum number of iterations
		solver.setMaxIter(500);

#ifdef USE_MULTIGRID

        solver.setPreconditioner(PCHYPRE_BOOMERAMG);
        solver.setPreconditionerAMG_nl(6);
        solver.setPreconditionerAMG_maxit(1);
        solver.setPreconditionerAMG_relax("SOR/Jacobi");
        solver.setPreconditionerAMG_cycleType("V",0,4);
        solver.setPreconditionerAMG_coarsenNodalType(0);
        solver.setPreconditionerAMG_coarsen("HMIS");


#endif

		tm_solve.start();

		// Give to the solver A and b, return x, the solution
		solver.solve(fd.getA(),x_,fd.getB());

		tm_solve.stop();
	}
	else
	{
		tm_solve.start();
		solver.solve(x_,fd.getB());
		tm_solve.stop();
	}

	// copy the solution x to the grid psi
	fd.template copy<phi>(x_,psi);

	psi.template ghost_get<phi>();

	// Correct the vorticity to make it divergence free

	auto it2 = gr.getDomainIterator();

	while (it2.isNext())
	{
		auto key = it2.get();

		gr.template get<vorticity>(key)[x] -= 1.0f/2.0f/psi.spacing(x)*(psi.template get<phi>(key.move(x,1))-psi.template get<phi>(key.move(x,-1)));
		gr.template get<vorticity>(key)[y] -= 1.0f/2.0f/psi.spacing(y)*(psi.template get<phi>(key.move(y,1))-psi.template get<phi>(key.move(y,-1)));
		gr.template get<vorticity>(key)[z] -= 1.0f/2.0f/psi.spacing(z)*(psi.template get<phi>(key.move(z,1))-psi.template get<phi>(key.move(z,-1)));

		++it2;
	}

	calc_and_print_max_div_and_int(gr);
}

void remesh(particles_type & vd, grid_type & gr,Box<3,float> & domain)
{
	// Remove all particles because we reinitialize in a grid like position
	vd.clear();

	// Get a grid iterator
	auto git = vd.getGridIterator(gr.getGridInfo().getSize());

	// For each grid point
	while (git.isNext())
	{
		// Get the grid point in global coordinates (key). p is in local coordinates
		auto p = git.get();
		auto key = git.get_dist();

		// Add the particle
		vd.add();

		// Assign the position to the particle in a grid like way
		vd.getLastPos()[x] = gr.spacing(x)*p.get(x) + domain.getLow(x);
		vd.getLastPos()[y] = gr.spacing(y)*p.get(y) + domain.getLow(y);
		vd.getLastPos()[z] = gr.spacing(z)*p.get(z) + domain.getLow(z);

		// Initialize the vorticity
		vd.template getLastProp<vorticity>()[x] = gr.template get<vorticity>(key)[x];
		vd.template getLastProp<vorticity>()[y] = gr.template get<vorticity>(key)[y];
		vd.template getLastProp<vorticity>()[z] = gr.template get<vorticity>(key)[z];

		// next grid point
		++git;
	}

	// redistribute the particles
	vd.map();
}


void comp_vel(grid_type_s & gr_ps,
		      grid_type & phi_v,
			  FDScheme<poisson_nn_helm> & fd,
		      Box<3,float> & domain,
			  grid_type & g_vort,
			  grid_type & g_vel,
			  petsc_solver<double>::return_type (& phi_s)[3],
			  petsc_solver<double> & solver)
{
	// We calculate and print the maximum divergence of the vorticity
	// and the integral of it
	calc_and_print_max_div_and_int(g_vort);

	// For each component solve a poisson
	for (size_t i = 0 ; i < 3 ; i++)
	{
		//! \cond [solve_poisson_comp] \endcond

		// Copy the vorticity component i in gr_ps
		auto it2 = gr_ps.getDomainIterator();

		// calculate the velocity from the curl of phi
		while (it2.isNext())
		{
			auto key = it2.get();

			// copy
			gr_ps.get<vorticity>(key) = g_vort.template get<vorticity>(key)[i];

			++it2;
		}

		// impose the poisson equation using gr_ps = vorticity for the right-hand-side (on the full domain)
		fd.new_b();
		fd.template impose_dit_b<phi>(gr_ps,gr_ps.getDomainIterator());

		solver.setAbsTol(0.01);

		// Get the vector that represent the right-hand-side
		Vector<double,PETSC_BASE> & b = fd.getB();

		timer tm_solve;
		tm_solve.start();

		// Give to the solver A and b in this case we are giving
		// an intitial guess phi_s calculated in the previous
		// time step
		solver.solve(phi_s[i],b);

		tm_solve.stop();

		// Calculate the residual

		solError serr;
		serr = solver.get_residual_error(phi_s[i],b);

		Vcluster<> & v_cl = create_vcluster();
		if (v_cl.getProcessUnitID() == 0)
		{std::cout << "Solved component " << i << "  Error: " << serr.err_inf << std::endl;}

		// copy the solution to grid
		fd.template copy<phi>(phi_s[i],gr_ps);

		auto it3 = gr_ps.getDomainIterator();

		// calculate the velocity from the curl of phi
		while (it3.isNext())
		{
			auto key = it3.get();

			phi_v.get<velocity>(key)[i] = gr_ps.get<phi>(key);

			++it3;
		}
	}

	phi_v.ghost_get<phi>();

	float inv_dx = 0.5f / g_vort.spacing(0);
	float inv_dy = 0.5f / g_vort.spacing(1);
	float inv_dz = 0.5f / g_vort.spacing(2);

	auto it3 = phi_v.getDomainIterator();

	// calculate the velocity from the curl of phi
	while (it3.isNext())
	{
		auto key = it3.get();

		float phi_xy = (phi_v.get<phi>(key.move(y,1))[x] - phi_v.get<phi>(key.move(y,-1))[x])*inv_dy;
		float phi_xz = (phi_v.get<phi>(key.move(z,1))[x] - phi_v.get<phi>(key.move(z,-1))[x])*inv_dz;
		float phi_yx = (phi_v.get<phi>(key.move(x,1))[y] - phi_v.get<phi>(key.move(x,-1))[y])*inv_dx;
		float phi_yz = (phi_v.get<phi>(key.move(z,1))[y] - phi_v.get<phi>(key.move(z,-1))[y])*inv_dz;
		float phi_zx = (phi_v.get<phi>(key.move(x,1))[z] - phi_v.get<phi>(key.move(x,-1))[z])*inv_dx;
		float phi_zy = (phi_v.get<phi>(key.move(y,1))[z] - phi_v.get<phi>(key.move(y,-1))[z])*inv_dy;

		g_vel.template get<velocity>(key)[x] = phi_zy - phi_yz;
		g_vel.template get<velocity>(key)[y] = phi_xz - phi_zx;
		g_vel.template get<velocity>(key)[z] = phi_yx - phi_xy;

		++it3;
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


template<typename grid> void calc_rhs(grid & g_vort, grid & g_vel, grid & g_dwp)
{
	// usefull constant
	constexpr int rhs = 0;

	// calculate several pre-factors for the stencil finite
	// difference
	float fac1 = 1.0f*nu/(g_vort.spacing(0)*g_vort.spacing(0));
	float fac2 = 1.0f*nu/(g_vort.spacing(1)*g_vort.spacing(1));
	float fac3 = 1.0f*nu/(g_vort.spacing(2)*g_vort.spacing(2));

	float fac4 = 0.5f/(g_vort.spacing(0));
	float fac5 = 0.5f/(g_vort.spacing(1));
	float fac6 = 0.5f/(g_vort.spacing(2));

	auto it = g_dwp.getDomainIterator();

	g_vort.template ghost_get<vorticity>();

	while (it.isNext())
	{
		auto key = it.get();

		g_dwp.template get<rhs>(key)[x] = fac1*(g_vort.template get<vorticity>(key.move(x,1))[x]+g_vort.template get<vorticity>(key.move(x,-1))[x])+
				                          fac2*(g_vort.template get<vorticity>(key.move(y,1))[x]+g_vort.template get<vorticity>(key.move(y,-1))[x])+
										  fac3*(g_vort.template get<vorticity>(key.move(z,1))[x]+g_vort.template get<vorticity>(key.move(z,-1))[x])-
										  2.0f*(fac1+fac2+fac3)*g_vort.template get<vorticity>(key)[x]+
										  fac4*g_vort.template get<vorticity>(key)[x]*
										  (g_vel.template get<velocity>(key.move(x,1))[x] - g_vel.template get<velocity>(key.move(x,-1))[x])+
										  fac5*g_vort.template get<vorticity>(key)[y]*
										  (g_vel.template get<velocity>(key.move(y,1))[x] - g_vel.template get<velocity>(key.move(y,-1))[x])+
										  fac6*g_vort.template get<vorticity>(key)[z]*
										  (g_vel.template get<velocity>(key.move(z,1))[x] - g_vel.template get<velocity>(key.move(z,-1))[x]);

		g_dwp.template get<rhs>(key)[y] = fac1*(g_vort.template get<vorticity>(key.move(x,1))[y]+g_vort.template get<vorticity>(key.move(x,-1))[y])+
				                          fac2*(g_vort.template get<vorticity>(key.move(y,1))[y]+g_vort.template get<vorticity>(key.move(y,-1))[y])+
										  fac3*(g_vort.template get<vorticity>(key.move(z,1))[y]+g_vort.template get<vorticity>(key.move(z,-1))[y])-
										  2.0f*(fac1+fac2+fac3)*g_vort.template get<vorticity>(key)[y]+
										  fac4*g_vort.template get<vorticity>(key)[x]*
										  (g_vel.template get<velocity>(key.move(x,1))[y] - g_vel.template get<velocity>(key.move(x,-1))[y])+
										  fac5*g_vort.template get<vorticity>(key)[y]*
										  (g_vel.template get<velocity>(key.move(y,1))[y] - g_vel.template get<velocity>(key.move(y,-1))[y])+
										  fac6*g_vort.template get<vorticity>(key)[z]*
										  (g_vel.template get<velocity>(key.move(z,1))[y] - g_vel.template get<velocity>(key.move(z,-1))[y]);


		g_dwp.template get<rhs>(key)[z] = fac1*(g_vort.template get<vorticity>(key.move(x,1))[z]+g_vort.template get<vorticity>(key.move(x,-1))[z])+
				                          fac2*(g_vort.template get<vorticity>(key.move(y,1))[z]+g_vort.template get<vorticity>(key.move(y,-1))[z])+
										  fac3*(g_vort.template get<vorticity>(key.move(z,1))[z]+g_vort.template get<vorticity>(key.move(z,-1))[z])-
										  2.0f*(fac1+fac2+fac3)*g_vort.template get<vorticity>(key)[z]+
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

template<typename grid, typename vector> void do_step(grid_type_s & psi,
													  grid_type & phi_v,
													  FDScheme<poisson_nn_helm> & fd,
													  vector & particles,
		                                              grid & g_vort,
													  grid & g_vel,
													  grid & g_dvort,
													  Box<3,float> & domain,
													  interpolate<particles_type,grid_type,mp4_kernel<float>> & inte,
													  petsc_solver<double>::return_type (& phi_s)[3],
													  petsc_solver<double> & solver)
{
	constexpr int rhs = 0;

	set_zero<vorticity>(g_vort);
	inte.template p2m<vorticity,vorticity>(particles,g_vort);

	g_vort.template ghost_put<add_,vorticity>();

	// Calculate velocity from vorticity
	comp_vel(psi,phi_v,fd,domain,g_vort,g_vel,phi_s,solver);
	calc_rhs(g_vort,g_vel,g_dvort);

	g_dvort.template ghost_get<rhs>();
	g_vel.template ghost_get<velocity>();
	set_zero<rhs_part>(particles);
	set_zero<p_vel>(particles);
	inte.template m2p<rhs,rhs_part>(g_dvort,particles);
	inte.template m2p<velocity,p_vel>(g_vel,particles);

	//! \cond [inte_m2p] \endcond
}

template<typename vector, typename grid> void check_point_and_save(vector & particles,
																   grid & g_vort,
																   grid & g_vel,
																   grid & g_dvort,
																   size_t i)
{
	Vcluster<> & v_cl = create_vcluster();

	if (v_cl.getProcessingUnits() < 24)
	{
		particles.write("part_out_" + std::to_string(i),VTK_WRITER | FORMAT_BINARY);
		g_vort.template ghost_get<vorticity>();
		g_vel.template ghost_get<velocity>();
		g_vel.write_frame("grid_velocity",i, VTK_WRITER | FORMAT_BINARY);
		g_vort.write_frame("grid_vorticity",i, VTK_WRITER | FORMAT_BINARY);
	}

	// In order to reduce the size of the saved data we apply a threshold.
	// We only save particles with vorticity higher than 0.1
	particles_type_s part_save(particles.getDecomposition(),0);

	auto it_s = particles.getDomainIterator();

	while (it_s.isNext())
	{
		auto p = it_s.get();

		float vort_magn = sqrt(particles.template getProp<vorticity>(p)[0] * particles.template getProp<vorticity>(p)[0] +
							   particles.template getProp<vorticity>(p)[1] * particles.template getProp<vorticity>(p)[1] +
							   particles.template getProp<vorticity>(p)[2] * particles.template getProp<vorticity>(p)[2]);

		if (vort_magn < 0.1)
		{
			++it_s;
			continue;
		}

		part_save.add();

		part_save.getLastPos()[0] = particles.getPos(p)[0];
		part_save.getLastPos()[1] = particles.getPos(p)[1];
		part_save.getLastPos()[2] = particles.getPos(p)[2];

		part_save.template getLastProp<0>() = vort_magn;

		++it_s;
	}

	part_save.map();

	particles.save("check_point");
}

int main(int argc, char* argv[])
{
	// Initialize
//    openfpm_init(&argc,&argv);

    openfpm_init(&argc,&argv, init_options::in_situ_visualization);
	{
	// Domain, a rectangle
	// For the grid 1600x400x400 use
	// Box<3,float> domain({0.0,0.0,0.0},{22.0,5.57,5.57});
	Box<3,float> domain({0.0,0.0,0.0},{3.57,3.57,3.57});

	// Ghost (Not important in this case but required)
	Ghost<3,long int> g(2);

	// Grid points on x=128 y=64 z=64
	// if we use Re = 7500
	//  long int sz[] = {1600,400,400};
	long int sz[] = {128,128,128};
	size_t szu[] = {(size_t)sz[0],(size_t)sz[1],(size_t)sz[2]};

	periodicity<3> bc = {{PERIODIC,PERIODIC,PERIODIC}};

	//! \cond [construct grids] \endcond

	grid_type g_vort(szu,domain,g,bc);
	grid_type g_vel(g_vort.getDecomposition(),szu,g);
	grid_type g_dvort(g_vort.getDecomposition(),szu,g);
	grid_type_vis g_vis(g_vort.getDecomposition(),szu,g);
	g_vis.visualize();
	particles_type particles(g_vort.getDecomposition(),0);

	// Construct an FDScheme is heavy so we construct it here
	// We also construct distributed temporal grids here because
	// they are expensive

	// Here we create temporal distributed grid, create a distributed grid is expensive we do it once outside
	// And the vectorial phi_v
	grid_type_s psi(g_vort.getDecomposition(),g_vort.getGridInfo().getSize(),g);
	grid_type phi_v(g_vort.getDecomposition(),g_vort.getGridInfo().getSize(),g);

	// In order to create a matrix that represent the poisson equation we have to indicate
	// we have to indicate the maximum extension of the stencil and we we need an extra layer
	// of points in case we have to impose particular boundary conditions.
	Ghost<3,long int> stencil_max(2);

	// We define a field phi_f
	typedef Field<phi,poisson_nn_helm> phi_f;

	// We assemble the poisson equation doing the
	// poisson of the Laplacian of phi using Laplacian
	// central scheme (where the both derivative use central scheme, not
	// forward composed backward like usually)
	typedef Lap<phi_f,poisson_nn_helm,CENTRAL_SYM> poisson;

	FDScheme<poisson_nn_helm> fd(stencil_max, domain, psi);

	fd.template impose_dit<0>(poisson(),psi,psi.getDomainIterator());

	//! \cond [construct grids] \endcond

	// It store the solution to compute velocity
	// It is used as initial guess every time we call the solver
	Vector<double,PETSC_BASE> phi_s[3];
	Vector<double,PETSC_BASE> x_;

	// Parallel object
	Vcluster<> & v_cl = create_vcluster();

	// print out the spacing of the grid
	if (v_cl.getProcessUnitID() == 0)
	{std::cout << "SPACING " << g_vort.spacing(0) << " " << g_vort.spacing(1) << " " << g_vort.spacing(2) << std::endl;}

	// initialize the ring step 1
	init_ring(g_vort,domain);

	x_.resize(g_vort.size(),g_vort.getLocalDomainSize());
	x_.setZero();

	// Create a PETSC solver to get the solution x
	petsc_solver<double> solver;

	// Do the helmotz projection step 2
	helmotz_hodge_projection(psi,fd,g_vort,domain,solver,x_,true);

	// initialize the particle on a mesh step 3
	remesh(particles,g_vort,domain);

	// Get the total number of particles
	size_t tot_part = particles.size_local();
	v_cl.sum(tot_part);
	v_cl.execute();

	// Now we set the size of phi_s
	phi_s[0].resize(g_vort.size(),g_vort.getLocalDomainSize());
	phi_s[1].resize(g_vort.size(),g_vort.getLocalDomainSize());
	phi_s[2].resize(g_vort.size(),g_vort.getLocalDomainSize());

	// And we initially set it to zero
	phi_s[0].setZero();
	phi_s[1].setZero();
	phi_s[2].setZero();

	// We create the interpolation object outside the loop cycles
	// to avoid to recreate it at every cycle. Internally this object
	// create some data-structure to optimize the conversion particle
	// position to sub-domain. If there are no re-balancing it is safe
	// to reuse-it
	interpolate<particles_type,grid_type,mp4_kernel<float>> inte(particles,g_vort);

	// With more than 24 core we rely on the HDF5 checkpoint restart file
	if (v_cl.getProcessingUnits() < 24)
	{g_vort.write("grid_vorticity_init", VTK_WRITER | FORMAT_BINARY);}


	// Before entering the loop we check if we want to restart from
	// a previous simulation
	size_t i = 0;

	// if we have an argument
	if (argc > 1)
	{
		// take that argument
		std::string restart(argv[1]);

		// convert it into number
		i = std::stoi(restart);

		// load the file
		particles.load("check_point_" + std::to_string(i));

		// Print to inform that we are restarting from a
		// particular position
		if (v_cl.getProcessUnitID() == 0)
		{std::cout << "Restarting from " << i << std::endl;}
	}

	// Time Integration
	for ( ; i < nsteps ; i++)
	{
        if (v_cl.getProcessUnitID() == 0)
	    std::cout<<"In the time loop " <<std::endl;
		// do step 4-5-6-7
		do_step(psi,phi_v,fd,particles,g_vort,g_vel,g_dvort,domain,inte,phi_s,solver);

		// do step 8
		rk_step1(particles);

		// do step 9-10-11-12
		do_step(psi,phi_v,fd,particles,g_vort,g_vel,g_dvort,domain,inte,phi_s,solver);

		// do step 13
		rk_step2(particles);

		// so step 14
		set_zero<vorticity>(g_vort);
		inte.template p2m<vorticity,vorticity>(particles,g_vort);
		g_vort.template ghost_put<add_,vorticity>();

		// helmotz-hodge projection
		helmotz_hodge_projection(psi,fd,g_vort,domain,solver,x_,false);

		remesh(particles,g_vort,domain);

		// print the step number
		if (v_cl.getProcessUnitID() == 0)
		{std::cout << "Step " << i << std::endl;}

		// every 100 steps write the output
		if (i % 100 == 0)		{check_point_and_save(particles,g_vort,g_vel,g_dvort,i);}


        //find max and min velocity
        auto it1 = g_vel.getDomainIterator();
		float maxVel = 0.0f;
		float minVel = 65535.0f;

		while(it1.isNext())
        {
		    auto key = it1.get();

            float curVel = (float) sqrt( g_vel.template get<velocity>(key)[0] * g_vel.template get<velocity>(key)[0] +
                                       g_vel.template get<velocity>(key)[1] * g_vel.template get<velocity>(key)[1] +
                                       g_vel.template get<velocity>(key)[2] * g_vel.template get<velocity>(key)[2] );

            if(curVel > maxVel)
            {
                maxVel = curVel;
            }

            if(curVel < minVel)
            {
                minVel = curVel;
            }
            ++it1;
        }

        if (v_cl.getProcessUnitID() == 0)
		std::cout<<"The maximum velocity is "<<maxVel << " and the minimum is " << minVel <<std::endl;

        // calculate the magnitude of velocity
        auto it2 = g_vel.getDomainIterator();
        while (it2.isNext())
		{
			auto key = it2.get();

			float curVel = (float) sqrt( g_vel.template get<velocity>(key)[0] * g_vel.template get<velocity>(key)[0] +
                                                               g_vel.template get<velocity>(key)[1] * g_vel.template get<velocity>(key)[1] +
                                                               g_vel.template get<velocity>(key)[2] * g_vel.template get<velocity>(key)[2] );

			float scaled = (curVel / (maxVel - minVel)) * 65535;
			// copy
			g_vis.get<0>(key) = (unsigned short)(scaled);

			++it2;
		}

	}
	}

	openfpm_finalize();
}

#else

int main(int argc, char* argv[])
{
        return 0;
}

#endif
