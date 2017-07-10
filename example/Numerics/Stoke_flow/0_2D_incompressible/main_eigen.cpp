/*! \page Stokes_flow Stokes flow
 *
 * \subpage Stokes_0_2D
 * \subpage Stokes_0_2D_petsc
 * \subpage Stokes_1_3D
 * \subpage Stokes_1_3D_petsc
 *
 */

/*!
 * \page Stokes_0_2D Stokes incompressible 2D eigen
 *
 * # Stokes incompressible 2D # {#num_sk_inc_2D}
 *
 * In this example we try to solve the 2D stokes equation for incompressible
 * fluid with Reynold number = 0
 *
 * Such equation require the inversion of a system.  We will
 * show how to produce such system and solve it using Finite differences with
 * staggered grid. The system of equation to solve is the following
 *
 *
 * \f[ \eta \partial^{2}_{x} v_x + \eta \partial^{2}_{y} v_x - \partial_{x} P  = 0 \f]
 * \f[ \eta \partial^{2}_{x} v_y + \eta \partial^{2}_{y} v_y - \partial_{y} P  = 0 \f]
 * \f[ \partial_{x} v_x + \partial_{y} v_y  = 0 \f]
 *
 * \f$ v_x = 0 \quad v_y = 1 \f$ at x = L
 *
 * \f$ v_x = 0 \quad v_y = 0 \f$ otherwise
 *
 * ## General Model properties ## {#num_sk_inc_2D_gmp}
 *
 * In order to do this we have to define a structure that contain the main information
 * about the system
 *
 * \snippet Numerics/Stoke_flow/0_2D_incompressible/main_eigen.cpp def system
 *
 */

//! \cond [def system] \endcond

#include "Grid/grid_dist_id.hpp"
#include "Matrix/SparseMatrix.hpp"
#include "Vector/Vector.hpp"
#include "FiniteDifference/FDScheme.hpp"
#include "FiniteDifference/util/common.hpp"
#include "FiniteDifference/eq.hpp"
#include "Solvers/umfpack_solver.hpp"
#include "Solvers/petsc_solver.hpp"

struct lid_nn
{
	// dimensionaly of the equation (2D problem 3D problem ...)
	static const unsigned int dims = 2;

	// number of fields in the system v_x, v_y, P so a total of 3
	static const unsigned int nvar = 3;

	// boundary conditions PERIODIC OR NON_PERIODIC
	static const bool boundary[];

	// type of space float, double, ...
	typedef float stype;

	// type of base grid, it is the distributed grid that will store the result
	// Note the first property is a 2D vector (velocity), the second is a scalar (Pressure)
	typedef grid_dist_id<2,float,aggregate<float[2],float>,CartDecomposition<2,float>> b_grid;

	// type of SparseMatrix, for the linear system, this parameter is bounded by the solver
	// that you are using, in case of umfpack using <double,int> it is the only possible choice
	// By default SparseMatrix is EIGEN based
	typedef SparseMatrix<double,int> SparseMatrix_type;

	// type of Vector for the linear system, this parameter is bounded by the solver
	// that you are using, in case of umfpack using <double> it is the only possible choice
	typedef Vector<double> Vector_type;

	// Define that the underline grid where we discretize the system of equation is staggered
	static const int grid_type = STAGGERED_GRID;
};

const bool lid_nn::boundary[] = {NON_PERIODIC,NON_PERIODIC};

//! \cond [def system] \endcond


/*!
 *
 * \page Stokes_0_2D Stokes incompressible 2D eigen
 *
 * ## System of equation modeling ## {#num_sk_inc_2D_sem}
 *
 * We model the equations. This part will change in the near future to use template expression
 * parsing and simplify the process of defining equations.
 *
 * \f$  \eta v_x \nabla v_x = eta\_lap\_vx \quad \nabla = \partial^{2}_{x} + \partial^{2}_{y} \f$ Step1
 *
 * \f$  \partial_{x} P = p\_x \f$ Step 2
 *
 * \f$  -p\_x = m\_ p\_ x \f$ Step 3
 *
 * \f$ eta\_lap\_vx - m\_p\_x \f$ Step4. This is the first equation in the system
 *
 * The second equation definition is similar to the first one
 *
 * \f$ \partial^{forward}_{x} v_x = dx\_vx \f$ Step5
 *
 * \f$ \partial^{forward}_{y} v_y = dy\_vy \f$ Step6
 *
 * \f$ dx\_vx + dy\_vy \f$ Step 7. This is the third equation in the system
 *
 *
 * \snippet Numerics/Stoke_flow/0_2D_incompressible/main_eigen.cpp def equation
 *
 */

//! \cond [def equation] \endcond

// Constant Field
struct eta
{
	typedef void const_field;

	static float val()	{return 1.0;}
};

// Convenient constants
constexpr unsigned int v[] = {0,1};
constexpr unsigned int P = 2;
constexpr unsigned int ic = 2;
constexpr int x = 0;
constexpr int y = 1;

// Create field that we have v_x, v_y, P
typedef Field<v[x],lid_nn> v_x;                 // Definition 1 v_x
typedef Field<v[y],lid_nn> v_y;                 // Definition 2 v_y
typedef Field<P,lid_nn> Prs;					// Definition 3 Prs

// Eq1 V_x

typedef mul<eta,Lap<v_x,lid_nn>,lid_nn> eta_lap_vx;       // Step 1
typedef D<x,Prs,lid_nn> p_x;                              // Step 2
typedef minus<p_x,lid_nn> m_p_x;                          // Step 3
typedef sum<eta_lap_vx,m_p_x,lid_nn> vx_eq;               // Step 4

// Eq2 V_y

typedef mul<eta,Lap<v_y,lid_nn>,lid_nn> eta_lap_vy;
typedef D<y,Prs,lid_nn> p_y;
typedef minus<p_y,lid_nn> m_p_y;
typedef sum<eta_lap_vy,m_p_y,lid_nn> vy_eq;

// Eq3 Incompressibility

typedef D<x,v_x,lid_nn,FORWARD> dx_vx;                   // Step 5
typedef D<y,v_y,lid_nn,FORWARD> dy_vy;                   // Step 6
typedef sum<dx_vx,dy_vy,lid_nn> ic_eq;                   // Step 7

//! \cond [def equation] \endcond

/*!
 * \page Stokes_0_2D Stokes incompressible 2D eigen
 *
 * In case of boundary conditions and staggered grid we need a particular set of equations
 * at boundaries. Explain in detail is out of the scope of this example, but to qualitatively
 * have an idea consider the following staggered cell
 *
 	 	 \verbatim

		+--$--+
		|     |
		#  *  #
		|     |
		0--$--+

	  # = velocity(x)
	  $ = velocity(y)
	  * = pressure

		\endverbatim
 *
 * As we can see several properties has different position in the cell.
 * Consider we want to impose \f$v_y = 0\f$ on \f$x=0\f$. Because there are not
 * points at \f$x = 0\f$ we have to interpolate between $ of this cell
 * and $ of the previous cell on y Averaging using the Avg operator. The same apply for
 * \f$v_x\f$ on  \f$y=0\f$. Similar arguments can be done for other boundaries in
 *  order to finally reach a well defined system
 *
 * \snippet Numerics/Stoke_flow/0_2D_incompressible/main_eigen.cpp bond def eq
 *
 */

//! \cond [bond def eq] \endcond

// Equation for boundary conditions

// Directional Avg
typedef Avg<x,v_y,lid_nn> avg_vy;
typedef Avg<y,v_x,lid_nn> avg_vx;

typedef Avg<x,v_y,lid_nn,FORWARD> avg_vy_f;
typedef Avg<y,v_x,lid_nn,FORWARD> avg_vx_f;

// Usefull constants (as MACRO)
#define EQ_1 0
#define EQ_2 1
#define EQ_3 2

//! \cond [bond def eq] \endcond

#include "Vector/vector_dist.hpp"
#include "data_type/aggregate.hpp"

int main(int argc, char* argv[])
{
	/*!
	 * \page Stokes_0_2D Stokes incompressible 2D eigen
	 *
	 * ## Initialization ## {#num_sk_inc_2D_init}
	 *
	 * After model our equation we:
	 * * Initialize the library
	 * * Define some useful constants
	 * * define Ghost size
	 * * Non-periodic boundary conditions
	 * * Padding domain expansion
	 *
	 * Padding and Ghost differ in the fact the padding extend the domain.
	 * Ghost is an extension for each sub-domain
	 *
	 * \snippet Numerics/Stoke_flow/0_2D_incompressible/main_eigen.cpp init
	 *
	 */

	//! \cond [init] \endcond

	// Initialize
	openfpm_init(&argc,&argv);

	// velocity in the grid is the property 0, pressure is the property 1
	constexpr int velocity = 0;
	constexpr int pressure = 1;

	// Domain, a rectangle
	Box<2,float> domain({0.0,0.0},{3.0,1.0});

	// Ghost (Not important in this case but required)
	Ghost<2,float> g(0.01);

	// Grid points on x=256 and y=64
	long int sz[] = {96,32};
	size_t szu[2];
	szu[0] = (size_t)sz[0];
	szu[1] = (size_t)sz[1];

	// We need one more point on the left and down part of the domain
	// This is given by the boundary conditions that we impose.
	//
	Padding<2> pd({1,1},{0,0});

	//! \cond [init] \endcond

	/*!
	 * \page Stokes_0_2D Stokes incompressible 2D eigen
	 *
	 * Distributed grid that store the solution
	 *
	 * \see \ref e0_s_grid_inst
	 *
	 * \snippet Numerics/Stoke_flow/0_2D_incompressible/main_eigen.cpp grid inst
	 *
	 */

	//! \cond [grid inst] \endcond

	grid_dist_id<2,float,aggregate<float[2],float>> g_dist(szu,domain,g);

	//! \cond [grid inst] \endcond

	/*!
	 * \page Stokes_0_2D Stokes incompressible 2D eigen
	 *
	 * Solving the system above require the solution of a system like that
	 *
	 * \f$ Ax = b \quad x = A^{-1}b\f$
	 *
	 * where A is the system the discretize the left hand side of the equations + boundary conditions
	 * and b discretize the right hand size + boundary conditions
	 *
	 * FDScheme is the object that we use to produce the Matrix A and the vector b.
	 * Such object require the maximum extension of the stencil
	 *
	 * \snippet Numerics/Stoke_flow/0_2D_incompressible/main_eigen.cpp fd scheme
	 *
	 */

	//! \cond [fd scheme] \endcond

	// It is the maximum extension of the stencil (order 2 laplacian stencil has extension 1)
	Ghost<2,long int> stencil_max(1);

	// Finite difference scheme
	FDScheme<lid_nn> fd(pd, stencil_max, domain, g_dist);

	//! \cond [fd scheme] \endcond

	/*!
	 * \page Stokes_0_2D Stokes incompressible 2D eigen
	 *
	 * ## Impose the equation on the domain ## {#num_sk_inc_2D_ied}
	 *
	 * Here we impose the system of equation, we start from the incompressibility Eq imposed in the bulk with the
	 * exception of the first point {0,0} and than we set P = 0 in {0,0}, why we are doing this is again
	 * mathematical to have a well defined system, an intuitive explanation is that P and P + c are both
	 * solution for the incompressibility equation, this produce an ill-posed problem to make it well posed
	 * we set one point in this case {0,0} the pressure to a fixed constant for convenience P = 0
	 *
	 * The best way to understand what we are doing is to draw a smaller example like 8x8.
	 * Considering that we have one additional point on the left for padding we have a grid
	 * 9x9. If on each point we have v_x v_y and P unknown we have
	 * 9x9x3 = 243 unknown. In order to fully determine and unique solution we have to
	 * impose 243 condition. The code under impose (in the case of 9x9) between domain
	 * and bulk 243 conditions.
	 *
	 * \snippet Numerics/Stoke_flow/0_2D_incompressible/main_eigen.cpp impose eq dom
	 *
	 *
	 */

	//! \cond [impose eq dom] \endcond

	fd.impose(ic_eq(),0.0, EQ_3, {0,0},{sz[0]-2,sz[1]-2},true);
	fd.impose(Prs(),  0.0, EQ_3, {0,0},{0,0});

	// Here we impose the Eq1 and Eq2
	fd.impose(vx_eq(),0.0, EQ_1, {1,0},{sz[0]-2,sz[1]-2});
	fd.impose(vy_eq(),0.0, EQ_2, {0,1},{sz[0]-2,sz[1]-2});

	// v_x and v_y
	// Imposing B1
	fd.impose(v_x(),0.0, EQ_1, {0,0},{0,sz[1]-2});
	fd.impose(avg_vy_f(),0.0, EQ_2 , {-1,0},{-1,sz[1]-1});
	// Imposing B2
	fd.impose(v_x(),0.0, EQ_1, {sz[0]-1,0},{sz[0]-1,sz[1]-2});
	fd.impose(avg_vy(),1.0, EQ_2,    {sz[0]-1,0},{sz[0]-1,sz[1]-1});

	// Imposing B3
	fd.impose(avg_vx_f(),0.0, EQ_1, {0,-1},{sz[0]-1,-1});
	fd.impose(v_y(), 0.0, EQ_2, {0,0},{sz[0]-2,0});
	// Imposing B4
	fd.impose(avg_vx(),0.0, EQ_1,   {0,sz[1]-1},{sz[0]-1,sz[1]-1});
	fd.impose(v_y(), 0.0, EQ_2, {0,sz[1]-1},{sz[0]-2,sz[1]-1});

	// When we pad the grid, there are points of the grid that are not
	// touched by the previous condition. Mathematically this lead
	// to have too many variables for the conditions that we are imposing.
	// Here we are imposing variables that we do not touch to zero
	//

	// Padding pressure
	fd.impose(Prs(), 0.0, EQ_3, {-1,-1},{sz[0]-1,-1});
	fd.impose(Prs(), 0.0, EQ_3, {-1,sz[1]-1},{sz[0]-1,sz[1]-1});
	fd.impose(Prs(), 0.0, EQ_3, {-1,0},{-1,sz[1]-2});
	fd.impose(Prs(), 0.0, EQ_3, {sz[0]-1,0},{sz[0]-1,sz[1]-2});

	// Impose v_x Padding Impose v_y padding
	fd.impose(v_x(), 0.0, EQ_1, {-1,-1},{-1,sz[1]-1});
	fd.impose(v_y(), 0.0, EQ_2, {-1,-1},{sz[0]-1,-1});

	//! \cond [impose eq dom] \endcond

	/*!
	 * \page Stokes_0_2D Stokes incompressible 2D eigen
	 *
	 * ## Solve the system of equation ## {#num_sk_inc_2D_sse}
	 *
	 * Once we imposed all the equations we can retrieve the Matrix A and the vector b
	 * and pass these two element to the solver. In this example we are using a serial
	 * direct solver Umfpack.
	 *
	 * \snippet Numerics/Stoke_flow/0_2D_incompressible/main_eigen.cpp solver
	 *
	 */

	//! \cond [solver] \endcond

	// Create an UMFPACK solver
	umfpack_solver<double> solver;

	// Give to the solver A and b, return x, the solution
	auto x = solver.solve(fd.getA(),fd.getB());

	//! \cond [solver] \endcond

	/*!
	 * \page Stokes_0_2D Stokes incompressible 2D eigen
	 *
	 * ## Copy the solution on the grid and write on VTK ## {#num_sk_inc_2D_csg}
	 *
	 * Once we have the solution we copy it on the grid
	 *
	 * \snippet Numerics/Stoke_flow/0_2D_incompressible/main_eigen.cpp copy write
	 *
	 */

	//! \cond [copy write] \endcond

	fd.template copy<velocity,pressure>(x,{0,0},{sz[0]-1,sz[1]-1},g_dist);

	g_dist.write("lid_driven_cavity_p_umfpack");

	//! \cond [copy write] \endcond


	/*!
	 * \page Stokes_0_2D Stokes incompressible 2D eigen
	 *
	 * ## Finalize ## {#num_sk_inc_2D_fin}
	 *
	 *  At the very end of the program we have always to de-initialize the library
	 *
	 * \snippet Numerics/Stoke_flow/0_2D_incompressible/main_eigen.cpp fin lib
	 *
	 */

	//! \cond [fin lib] \endcond

	openfpm_finalize();

	//! \cond [fin lib] \endcond


	/*!
	 * \page Stokes_0_2D Stokes incompressible 2D eigen
	 *
	 * # Full code # {#num_sk_inc_2D_code}
	 *
	 * \include Numerics/Stoke_flow/0_2D_incompressible/main_eigen.cpp
	 *
	 */
}
