/*! \page Vortex_in_cell_petsc Vortex in Cell 3D
 *
 * # Vortex in Cell 3D ring with ringlets # {#vic_ringlets}
 *
 * In this example we solve the Navier-Stokes equation in the vortex formulation in 3D
 * for an incompressible fluid.
 *
 * ## Numerical method ## {#num_vic_mt}
 *
 * In this code we solve the Navier-stokes equation for incompressible fluid in the
 * vorticity formulation. We first recall the Navier-stokes equation in vorticity formulation
 *
 * \f$ \nabla \times u = -w \f$
 *
 * \f$  \frac{Dw}{dt} = (w \cdot \vec \nabla) u + \nu \nabla^{2} w \f$ (5)
 *
 * Where \f$w\f$ is the vorticity and \f$u\f$ is the velocity of the fluid.
 * With high Reynold number \f$Re = \frac{uL}{\nu}\f$ and the term \f$uL\f$ significantly
 * smaller than the reynold number we have that \f$\nu\f$ is small and the term \f$\nu \nabla^{2} w\f$ is
 * negligible. The algorithm can be expressed with the following pseudo code.
 *
 * \verbatim

	1) Initialize the vortex ring on grid
	2) Do an helmotz hodge projection to make the vorticity divergent free
	3) Initialize particles on the same position as the grid or remesh

	while (t < t_end) do
		4) Interpolate mesh vorticity on particles
		5) calculate velocity u from the vorticity w
		6) interpolate velocity u to particles
		7) calculate the right-hand-side on grid and interpolate on particles
		8) move particles accordingly to the velocity
		9) interpolate vorticity w back to grid
	end while

   \endverbatim
 *
 * This pseudo code show how to solve the equation above using euler integration
 * In case of Runge-kutta of order two the pseudo code change into
 *
 *
 * \verbatim

	1) Initialize the vortex ring on grid
	2) Do an helmotz hodge projection to make the vorticity divergent free
	3) Initialize particles on the same position as the grid or remesh

	while (t < t_end) do
		4) Interpolate mesh vorticity on particles
		5) calculate velocity u from the vorticity w
		6) interpolate velocity u to particles
		7) calculate the right-hand-side on grid and interpolate on particles
		8) move particles accordingly to the velocity and save the old position in x_old
		9) interpolate vorticity w back to grid

		10) Interpolate mesh vorticity on particles
		11) calculate velocity u from the vorticity w
		12) interpolate velocity u to particles
		13) calculate the right-hand-side on grid and interpolate on particles
		14) move particles accordingly to the velocity starting from x_old
	end while

   \endverbatim
 *
 * In the following we explain how each step is implemented in the code
 *
 * ## Inclusion ## {#num_vic_inc}
 *
 * This example code need several components. First because is a particle
 * mesh example we have to activate "grid_dist_id.hpp" and "vector_dist_id.hpp".
 * Because we use a finite-difference scheme and linear-algebra to calculate the
 * velocity out of the vorticity, we have to include "FDScheme.hpp" to produce
 * out of the finite difference scheme a matrix that represent the linear-system
 * to solve. "SparseMatrix.hpp" is the Sparse-Matrix that contain the linear
 * system. "Vector.hpp" is the data-structure that contain the solution of the
 * linear system. "petsc_solver.hpp" is the library to invert the linear system.
 * Because we have to interpolate between particles and grid we the to include
 * "interpolate.hpp" as interpolation kernel we use the mp4, so we include the
 * "mp4_kernel.hpp"
 *
 * For convenience we also define the particles type and the grid type and some
 * convenient constants
 *
 * \snippet Numerics/Vortex_in_cell/main_vic_petsc.cpp inclusion
 *
 *
 *
 */

//#define SE_CLASS1
//#define PRINT_STACKTRACE

//! \cond [inclusion] \endcond

#include "Grid/grid_dist_id.hpp"
#include "Vector/vector_dist.hpp"
#include "Matrix/SparseMatrix.hpp"
#include "Vector/Vector.hpp"
#include "FiniteDifference/FDScheme.hpp"
#include "Solvers/petsc_solver.hpp"
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
// Reynold number
float tgtre  = 10000.0;
// Noise factor for the ring vorticity on z
float ringnz = 0.00009;

// Kinematic viscosity
float nu = 0.0001535;

// Time step
float dt = 0.006;


constexpr unsigned int vorticity = 0;
constexpr unsigned int velocity = 0;
constexpr unsigned int p_vel = 1;
constexpr int rhs_part = 2;
constexpr unsigned int old_vort = 3;
constexpr unsigned int old_pos = 4;

//! \cond [inclusion] \endcond

/*! \page Vortex_in_cell_petsc Vortex in Cell 3D
 *
 * # Step 1: Initialization of the vortex ring # {#vic_ring_init}
 *
 * In this function we initialize the vortex ring. The vortex ring is
 * initialized accordingly to these formula.
 *
 * \f$ w(t = 0) =  \frac{\Gamma}{\pi \sigma^{2}} e^{-(s/ \sigma)^2} \f$
 *
 * \f$ s^2 = (z-z_c)^{2} + ((x-x_c)^2 + (y-y_c)^2 - R^2) \f$
 *
 * \f$ \Gamma = \nu Re \f$
 *
 * With this initialization the vortex ring look like the one in figure
 *
 * \image html int_vortex_arrow_small.jpg "Vortex ring initialization the arrow indicate the direction where the vortex point while the colors indicate the magnitude from blue (low) to red (high)"
 *
 * \snippet Numerics/Vortex_in_cell/main_vic_petsc.cpp init_vort
 *
 *
 *
 */

//! \cond [init_vort] \endcond

/*
 * gr is the grid where we are initializing the vortex ring
 * domain is the simulation domain
 *
 */
void init_ring(grid_type & gr, const Box<3,float> & domain)
{
	// To add some noise to the vortex ring we create two random
	// vector
	constexpr int nk = 32;

	float ak[nk];
	float bk[nk];

	for (size_t i = 0 ; i < nk ; i++)
	{
	     ak[i] = 0.0/*rand()/RAND_MAX*/;
	     bk[i] = 0.0/*rand()/RAND_MAX*/;
	}

	// We calculate the circuitation gamma
	float gamma = nu * tgtre;
	float rinv2 = 1.0f/(ringr2*ringr2);
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
        float theta1 = atan2((ty-2.5f),(tz-2.5f));

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

/*        theta1 = atan2((ty-2.5f),(tz-2.5f));
        rad1r  = sqrt((ty-2.5f)*(ty-2.5f) + (tz-2.5f)*(tz-2.5f)) + ringr1*(1.0f + ringnz * noise);
        rad1t = tx - 2.5f;
        rad1sq = rad1r*rad1r + rad1t*rad1t;
        float rad1sqTILDA = rad1sq*rinv2;
        radstr = exp(-rad1sqTILDA)*rinv2*gamma/M_PI;
        gr.template get<vorticity>(key_d)[x] = 0.0f;
        gr.template get<vorticity>(key_d)[y] = gr.template get<vorticity>(key_d)[y] + radstr * cos(theta1);
        gr.template get<vorticity>(key_d)[z] = gr.template get<vorticity>(key_d)[z] - radstr * sin(theta1);*/

		++it;
	}
}

//! \cond [init_vort] \endcond

//! \cond [poisson_syseq] \endcond

// Specification of the poisson equation for the helmotz-hodge projection
// 3D (dims = 3). The field is a scalar value (nvar = 1), bournary are periodic
// type of the the space is float. Final grid that will store \phi, the result (grid_dist_id<.....>)
// The other indicate which kind of Matrix to use to construct the linear system and
// which kind of vector to construct for the right hand side. Here we use a PETSC Sparse Matrix
// and PETSC vector. NORMAL_GRID indicate that is a standard grid (non-staggered)
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

//! \cond [poisson_syseq] \endcond


/*! \page Vortex_in_cell_petsc Vortex in Cell 3D
 *
 * # Step 2: Helmotz-hodge projection # {#vic_hlm_proj}
 *
 * The Helnotz hodge projection is required in order to make the vorticity divergent
 * free. The Helmotz-holde projection work in this way. A field can be divided into
 * a curl-free part and a divergent-free part.
 *
 * \f$ w = w_{rot} + w_{div} \f$
 *
 * with
 *
 * \f$ \vec \nabla \times w_{rot} = 0 \f$
 *
 * \f$  \nabla \cdot w_{div} = 0 \f$
 *
 * To have a vorticity divergent free we have to get the component (3) \f$w_{div} = w - w_{rot}\f$.
 * In particular it hold
 *
 * \f$ \nabla \cdot w = \nabla \cdot w_{rot} \f$
 *
 * Bacause \f$ \vec \nabla \times w_{rot} = 0 \f$ we can introduce a field \f$ \psi \f$
 * such that
 *
 * (2) \f$ w_{rot} = \vec \nabla \psi \f$
 *
 * Doing the  \f$  \nabla \cdot \vec \nabla \psi \f$ we obtain
 *
 * \f$ \nabla \cdot \vec \nabla \psi = \nabla^{2} \psi = \nabla \cdot w_{rot} = \nabla \cdot w \f$
 *
 * so we lead to this equation
 *
 * (1) \f$ \nabla^{2} \psi = \nabla \cdot w  \f$
 *
 * Solving the equation for \f$ \psi \f$ we can obtain \f$ w_{rot} \f$ doing the gradient of \f$ \psi \f$
 * and finally correct \f$ w \f$ obtaining \f$ w_{div} \f$
 *
 * The **helmotz_hodge_projection** function do this correction to the vorticity
 *
 * In particular it solve the equation (1) it calculate \f$ w_{rot} \f$
 * using (2) and correct the vorticity using using (3)
 *
 *
 * ## Poisson equation ##
 *
 * To solve a poisson equation on a grid using finite-difference, we need to create
 * an object that carry information about the system of equations
 *
 * \snippet Numerics/Vortex_in_cell/main_vic_petsc.cpp poisson_syseq
 *
 * Once created this object we can define the equation we are trying to solve
 * the left-hand-side of the equation \f$ \nabla^{2} \psi \f$
 *
 * \snippet Numerics/Vortex_in_cell/main_vic_petsc.cpp poisson_obj_eq
 *
 * Before to construct the linear system we also calculate the divergence of the
 * vorticity \f$ \nabla \cdot w \f$ that will be the right-hand-side
 *
 * \snippet Numerics/Vortex_in_cell/main_vic_petsc.cpp calc_div_vort
 *
 * For statistic purpose we also calculate the integral of the vorticity and
 * the maximum divergence of the vorticity
 *
 * \snippet Numerics/Vortex_in_cell/main_vic_petsc.cpp sanity_div_integral
 *
 * Finally we can create the object FDScheme using the object **poisson_nn_helm**
 * as template variable. In addition to the constructor we have to specify the maximum extension of the stencil, the domain and the
 * grid that will store the result. At this point we can impose an equation to construct
 * our SparseMatrix. In this example we are imposing the poisson equation with right hand
 * side equal to the divergence of vorticity (note: to avoid to create another field we use
 *  \f$ \psi \f$ to preliminary store the divergence of the vorticity). Imposing the
 *  equations produce an invertible SparseMatrix **A** and a right-hand-side Vector **b**.
 *
 * \snippet Numerics/Vortex_in_cell/main_vic_petsc.cpp create_fdscheme
 *
 * Because we need \f$ x = A^{-1}b \f$. We have to invert and solve a linear system.
 * In this case we use the Conjugate-gradient-Method an iterative solver. Such method
 * is controlled by two parameters. One is the tollerance that determine when the
 * method is converged, the second one is the maximum number of iterations to avoid that
 * the method go into infinite loop. After we set the parameters of the solver we can the
 * the solution **x**. Finaly we copy back the solution **x** into the grid \f$ \psi \f$.
 *
 * \snippet Numerics/Vortex_in_cell/main_vic_petsc.cpp div_free_check
 *
 * ### Note ###
 *
 * Because we are solving the poisson equation in periodic boundary conditions the Matrix has
 * determinat equal to zero. This mean that \f$ \psi \f$ has no unique solution (if it has one).
 * In order to recover one we have to ensure that the integral of the righ hand side or vorticity
 * is zero. (In our case is the case). We have to ensure that across time the integral of the
 * vorticity is conserved. (In our case is the case if we consider the \f$ \nu = 0 \f$ and \f$
 * \nabla \cdot w = 0 \f$ we can rewrite (5) in a conservative way \f$  \frac{Dw}{dt} = div(w \otimes v) \f$ ).
 * Is also good to notice that the solution that you get is the one with \f$ \int v  = 0 \f$
 *
 *
 *
 * \snippet Numerics/Vortex_in_cell/main_vic_petsc.cpp solve_petsc
 *
 * ## Correction ## {#vort_correction}
 *
 * After we got our solution for \f$ \psi \f$ we can calculate the correction of the vorticity
 * doing the gradient of \f$ \psi \f$.
 *
 * \snippet Numerics/Vortex_in_cell/main_vic_petsc.cpp vort_correction
 *
 * We also do a sanity check and we control that the vorticity remain
 * divergent-free. Getting the maximum value of the divergence and printing out
 * its value
 *
 *
 */

/*
 * gr vorticity grid where we apply the correction
 * domain simulation domain
 *
 */
void helmotz_hodge_projection(grid_type & gr, const Box<3,float> & domain)
{
	///////////////////////////////////////////////////////////////

	 //! \cond [poisson_obj_eq] \endcond

	// Convenient constant
	constexpr unsigned int phi = 0;

	// We define a field phi_f
	typedef Field<phi,poisson_nn_helm> phi_f;

	// We assemble the poisson equation doing the
	// poisson of the Laplacian of phi using Laplacian
	// central scheme (where the both derivative use central scheme, not
	// forward composed backward like usually)
	typedef Lap<phi_f,poisson_nn_helm,CENTRAL_SYM> poisson;

	//! \cond [poisson_obj_eq] \endcond

	//! \cond [calc_div_vort] \endcond

	// ghost get
	gr.template ghost_get<vorticity>();

	// ghost size of the psi function
    Ghost<3,long int> g(2);

	// Here we create a distributed grid to store the result of the helmotz projection
	grid_dist_id<3,float,aggregate<float>> psi(gr.getDecomposition(),gr.getGridInfo().getSize(),g);

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

	//! \cond [calc_div_vort] \endcond

	//! \cond [sanity_div_integral] \endcond

	// Get iterator across the domain points
	auto it5 = gr.getDomainIterator();

	// maximum vorticity
	double max_vort = 0.0;

	// Total integral
	double int_vort[3] = {0.0,0.0,0.0};

	// For each grid point
	while (it5.isNext())
	{
		auto key = it5.get();

		double tmp;

		// Calculate the vorticity
        tmp = 1.0f/gr.spacing(x)/2.0f*(gr.template get<vorticity>(key.move(x,1))[x] - gr.template get<vorticity>(key.move(x,-1))[x] ) +
			         	 	 1.0f/gr.spacing(y)/2.0f*(gr.template get<vorticity>(key.move(y,1))[y] - gr.template get<vorticity>(key.move(y,-1))[y] ) +
							 1.0f/gr.spacing(z)/2.0f*(gr.template get<vorticity>(key.move(z,1))[z] - gr.template get<vorticity>(key.move(z,-1))[z] );

        // Calculate the integral of the vorticity on each direction
        int_vort[x] += gr.template get<vorticity>(key)[x];
        int_vort[y] += gr.template get<vorticity>(key)[y];
        int_vort[z] += gr.template get<vorticity>(key)[z];

        // Store the maximum vorticity
        if (tmp > max_vort)
        	max_vort = tmp;

		++it5;
	}

	std::cout << "Max div for vorticity " << max_vort << "   Integral: " << int_vort[0] << "  " << int_vort[1] << "   " << int_vort[2] << std::endl;

	//! \cond [sanity_div_integral] \endcond

	//! \cond [create_fdscheme] \endcond

	// In order to create a matrix that represent the poisson equation we have to indicate
	// we have to indicate the maximum extension of the stencil and we we need an extra layer
	// of points in case we have to impose particular boundary conditions.
	Ghost<3,long int> stencil_max(2);

	// Here we define our Finite difference disctetization scheme object
	FDScheme<poisson_nn_helm> fd(stencil_max, domain, psi);

	// impose the poisson equation, using right hand side psi on the full domain (psi.getDomainIterator)
	// the template paramereter instead define which property of phi carry the righ-hand-side
	// in this case phi has only one property, so the property 0 carry the right hand side
	fd.template impose_dit<0>(poisson(),psi,psi.getDomainIterator());

	//! \cond [create_fdscheme] \endcond

	//! \cond [solve_petsc] \endcond

	// Create a PETSC solver to get the solution x
	petsc_solver<double> solver;

	// Set the conjugate-gradient as method to solve the linear solver
	solver.setSolver(KSPCG);

	// Set the absolute tolerance to determine that the method is converged
	solver.setAbsTol(0.001);

	// Set the maximum number of iterations
	solver.setMaxIter(500);

	// Give to the solver A and b, return x, the solution
	auto x_ = solver.solve(fd.getA(),fd.getB());

	// copy the solution x to the grid psi
	fd.template copy<phi>(x_,psi);

	//! \cond [solve_petsc] \endcond

	//! \cond [vort_correction] \endcond

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

	//! \cond [vort_correction] \endcond

	//! \cond [div_free_check] \endcond

	// Here we check the maximum divergence of the ring
	// We also copy the solution to the original grid
	double max = 0.0;

	gr.template ghost_get<vorticity>();
	auto it3 = gr.getDomainIterator();

	while (it3.isNext())
	{
		auto key = it3.get();

        psi.template get<phi>(key) = 1.0f/gr.spacing(x)/2.0f*(gr.template get<vorticity>(key.move(x,1))[x] - gr.template get<vorticity>(key.move(x,-1))[x] ) +
			         	 	 1.0f/gr.spacing(y)/2.0f*(gr.template get<vorticity>(key.move(y,1))[y] - gr.template get<vorticity>(key.move(y,-1))[y] ) +
							 1.0f/gr.spacing(z)/2.0f*(gr.template get<vorticity>(key.move(z,1))[z] - gr.template get<vorticity>(key.move(z,-1))[z] );



        if (psi.template get<phi>(key) > max)
        	max = psi.template get<phi>(key);

		++it3;
	}

	std::cout << "Maximum divergence of the ring MAX " << max << std::endl;

	//! \cond [div_free_check] \endcond

}

/*! \page Vortex_in_cell_petsc Vortex in Cell 3D
 *
 * # Step 3: Remeshing vorticity # {#vic_remesh_vort}
 *
 * After that we initialized the vorticity on the grid, initialize the particles
 * in a grid like position and we interpolate the vorticity on particles. Because
 * of the particles position being in a grid-like position and the symmetry of the
 * interpolation kernels, the re-mesh step simply reduce to initialize the particle
 * in a grid like position and assign the property vorticity of the particles equal to the
 * grid vorticity.
 *
 * \snippet Numerics/Vortex_in_cell/main_vic_petsc.cpp remesh_part
 *
 */

//! \cond [remesh_part] \endcond

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

//! \cond [remesh_part] \endcond

/*! \page Vortex_in_cell_petsc Vortex in Cell 3D
 *
 * # Step 5: Compute velocity from vorticity # {#vic_vel_from_vort}
 *
 * Computing the velocity from vorticity is done in the following way. Given
 *
 * \f$ \vec \nabla \times u = -w \f$
 *
 * We intrododuce the stream line function defined as
 *
 * \f$ \nabla \times \phi = u \f$ (7)
 *
 * \f$ \nabla \cdot \phi = 0 \f$
 *
 * We obtain
 *
 * \f$ \nabla \times \nabla \times \phi = -w = \vec \nabla (\nabla \cdot  \phi) - \nabla^{2} \phi  \f$
 *
 * Because the divergence of \f$ \phi \f$ is defined to be zero we have
 *
 * \f$ \nabla^{2} \phi = w \f$
 *
 * The velocity can be recovered by the equation (7)
 *
 * Putting into code what explained before, we again generate a poisson
 * object
 *
 * \snippet Numerics/Vortex_in_cell/main_vic_petsc.cpp poisson_obj_eq
 *
 * For statistic reason we do a check of the total integral of the vorticity
 * and maximum divergence of the vorticity
 *
 * \snippet Numerics/Vortex_in_cell/main_vic_petsc.cpp sanity_int_div
 *
 * In order to calculate the velocity out of the vorticity, we solve a poisson
 * equation like we did in helmotz-projection equation, but we do it for each
 * component \f$ i \f$ of the vorticity.
 *
 * \snippet Numerics/Vortex_in_cell/main_vic_petsc.cpp solve_poisson_comp
 *
 * We save the \f$ \phi \f$ component into **phi_s**
 *
 * \snippet Numerics/Vortex_in_cell/main_vic_petsc.cpp solve_poisson_comp
 *
 * Once we filled phi_s we can implement (7) doing the curl of **phi_s**
 * and recovering the velocity v
 *
 * \snippet Numerics/Vortex_in_cell/main_vic_petsc.cpp curl_phi_v
 *
 * \snippet Numerics/Vortex_in_cell/main_vic_petsc.cpp solve_poisson_comp
 *
 * \snippet Numerics/Vortex_in_cell/main_vic_petsc.cpp cross_check_curl_u_vort
 *
 */

//! \cond [comp_vel] \endcond

void comp_vel(Box<3,float> & domain, grid_type & g_vort,grid_type & g_vel, petsc_solver<double>::return_type (& phi_s)[3])
{
	//! \cond [poisson_obj_eq] \endcond

	// Convenient constant
	constexpr unsigned int phi = 0;

	// We define a field phi_f
	typedef Field<phi,poisson_nn_helm> phi_f;

	// We assemble the poisson equation doing the
	// poisson of the Laplacian of phi using Laplacian
	// central scheme (where the both derivative use central scheme, not
	// forward composed backward like usually)
	typedef Lap<phi_f,poisson_nn_helm,CENTRAL_SYM> poisson;

	//! \cond [poisson_obj_eq] \endcond

	Ghost<3,long int> g(2);

	// Here we create a distributed grid to store the result of the helmotz projection
	grid_dist_id<3,float,aggregate<float>> gr_ps(g_vort.getDecomposition(),g_vort.getGridInfo().getSize(),g);
	grid_dist_id<3,float,aggregate<float[3]>> phi_v(g_vort.getDecomposition(),g_vort.getGridInfo().getSize(),g);

	// Check one we control the the divergence of vorticity is zero

	//! \cond [sanity_int_div] \endcond

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

	std::cout << "Max div for vorticity " << max_vort << "   Integral: " << int_vort[0] << "  " << int_vort[1] << "   " << int_vort[2] << std::endl;

	//! \cond [sanity_int_div] \endcond

	// For each component solve a poisson
	for (size_t i = 0 ; i < 3 ; i++)
	{
		//! \cond [solve_poisson_comp] \endcond

		// Copy the vorticity in gr_ps
		auto it2 = gr_ps.getDomainIterator();

		// calculate the velocity from the curl of phi
		while (it2.isNext())
		{
			auto key = it2.get();

			gr_ps.get<vorticity>(key) = g_vort.template get<vorticity>(key)[i];

			++it2;
		}

		Ghost<3,long int> stencil_max(2);

		// Finite difference scheme
		FDScheme<poisson_nn_helm> fd(stencil_max, domain, gr_ps);

		poisson ps;

		fd.template impose_dit<phi>(ps,gr_ps,gr_ps.getDomainIterator());

		// Create an PETSC solver
		petsc_solver<double> solver;

		solver.setSolver(KSPCG);
		solver.setAbsTol(0.001);
		solver.setMaxIter(500);

		std::cout << "Solving component " << i << std::endl;

		PetscBool flg;
		SparseMatrix<double,int,PETSC_BASE> & A = fd.getA();

		Vector<double,PETSC_BASE> & b = fd.getB();

		// Give to the solver A and b, return x, the solution
		solver.solve(A,phi_s[i],b);

		//! \cond [solve_poisson_comp] \endcond

		//! \cond [print_residual_copy] \endcond

		// Calculate the residual

		petsc_solver<double>::return_type r(gr_ps.size(),gr_ps.getLocalDomainSize());
		Vec & pr = r.getVec();

		PETSC_SAFE_CALL(MatResidual(A.getMat(),b.getVec(),phi_s[i].getVec(),pr));

		PetscReal ne;
		PETSC_SAFE_CALL(VecNorm(pr,NORM_INFINITY,&ne));


		std::cout << "Solved component " << i << "  Error: " << ne << "   Symmetric: " << flg << std::endl;

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

		//! \cond [print_residual_copy] \endcond
	}

	//! \cond [curl_phi_v] \endcond

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

	//! \cond [curl_phi_v] \endcond

	//! \cond [cross_check_curl_u_vort] \endcond

	g_vel.template ghost_get<velocity>();

	// We check that curl u = vorticity

	auto it4 = phi_v.getDomainIterator();

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

		phi_v.template get<velocity>(key)[x] = phi_zy - phi_yz  + g_vort.template get<vorticity>(key)[x];
		phi_v.template get<velocity>(key)[y] = phi_xz - phi_zx + g_vort.template get<vorticity>(key)[y];
		phi_v.template get<velocity>(key)[z] = phi_yx - phi_xy + g_vort.template get<vorticity>(key)[z];

		if (phi_v.template get<velocity>(key)[x] > norm_max)
			norm_max = phi_v.template get<velocity>(key)[x];

		if (phi_v.template get<velocity>(key)[y] > norm_max)
			norm_max = phi_v.template get<velocity>(key)[y];

		if (phi_v.template get<velocity>(key)[z] > norm_max)
			norm_max = phi_v.template get<velocity>(key)[z];

		++it4;
	}

	std::cout << "Norm max: " << norm_max << std::endl;

	//! \cond [cross_check_curl_u_vort] \endcond
}

//! \cond [comp_vel] \endcond

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

	vd.map();
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
	long int sz[] = {512,128,128};
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

	g_vort.write("grid_vorticity_init", VTK_WRITER | FORMAT_BINARY);


	// Before entering the loop we check if we have to restart a previous simulation

	size_t i = 0;

	if (argc > 1)
	{
		std::string restart(argv[1]);

		i = std::stoi(restart);

		particles.load("check_point_" + std::to_string(i));

		particles.map();

		std::cout << "Restarting from " << i << std::endl;
	}

	// Time Integration

	for ( ; i < 10001 ; i++)
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

		if (i % 100 == 0)
		{
			particles.write("part_out_" + std::to_string(i),VTK_WRITER | FORMAT_BINARY);

			g_vort.ghost_get<vorticity>();
			g_vel.ghost_get<velocity>();
			g_vel.write("grid_velocity_" + std::to_string(i), VTK_WRITER | FORMAT_BINARY);
			g_vort.write("grid_vorticity_" + std::to_string(i), VTK_WRITER | FORMAT_BINARY);

			// Save also for checkpoint restart

			particles.save("check_point_" + std::to_string(i+1));
		}

	}

	openfpm_finalize();
}

