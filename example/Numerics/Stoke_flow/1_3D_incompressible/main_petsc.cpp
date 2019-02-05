/*!
 * \page Stokes_1_3D_petsc Stokes incompressible 3D petsc
 *
 * # Stokes incompressible 3D # {#num_sk_inc_3D_petsc}
 *
 * In this example we try to solve the 3D stokes equation for incompressible
 * fluid with Reynold number = 0
 *
 * Such equation require the inversion of a system. We will
 * show how to produce such system and solve it using Finite differences with
 * staggered grid. The system of equation to solve is the following
 *
 *
 * \f[ \eta \partial^{2}_{x} v_x + \eta \partial^{2}_{y} v_x + \eta \partial^{2}_{z} v_x - \partial_{x} P  = 0 \f]
 * \f[ \eta \partial^{2}_{x} v_y + \eta \partial^{2}_{y} v_y + \eta \partial^{2}_{z} v_y - \partial_{y} P  = 0 \f]
 * \f[ \eta \partial^{2}_{x} v_z + \eta \partial^{2}_{y} v_z + \eta \partial^{2}_{z} v_z - \partial_{y} P  = 0 \f]
 * \f[ \partial_{x} v_x + \partial_{y} v_y + \partial_{z} v_z = 0 \f]
 *
 * \f$ v_x = 0 \quad v_y = 1 \quad v_z = 1 \f$ at x = L
 *
 * \f$ v_x = 0 \quad v_y = 0 \quad v_z = 0 \f$ otherwise
 *
 * ## General Model properties ## {#num_sk_inc_3D_gmp}
 *
 * In order to do this we have to define a structure that contain the main information
 * about the system
 *
 * \snippet Numerics/Stoke_flow/1_3D_incompressible/main_petsc.cpp def system
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
	// dimensionaly of the equation ( 3D problem ...)
	static const unsigned int dims = 3;

	// number of fields in the system v_x, v_y, v_z, P so a total of 4
	static const unsigned int nvar = 4;

	// boundary conditions PERIODIC OR NON_PERIODIC
	static const bool boundary[];

	// type of space float, double, ...
	typedef float stype;

	// type of base grid, it is the distributed grid that will store the result
	// Note the first property is a 3D vector (velocity), the second is a scalar (Pressure)
	typedef grid_dist_id<3,float,aggregate<float[3],float>,CartDecomposition<3,float>> b_grid;

	// type of SparseMatrix, for the linear system, this parameter is bounded by the solver
	// that you are using, in case of umfpack using <double,int> it is the only possible choice
	// We select PETSC SparseMatrix implementation
	typedef SparseMatrix<double,int,PETSC_BASE> SparseMatrix_type;

	// type of Vector for the linear system, this parameter is bounded by the solver
	// that you are using, in case of umfpack using <double> it is the only possible choice
	// We select PETSC implementation
	typedef Vector<double,PETSC_BASE> Vector_type;

	// Define that the underline grid where we discretize the system of equation is staggered
	static const int grid_type = STAGGERED_GRID;
};

const bool lid_nn::boundary[] = {NON_PERIODIC,NON_PERIODIC,NON_PERIODIC};

//! \cond [def system] \endcond


/*!
 *
 * \page Stokes_1_3D_petsc Stokes incompressible 3D petsc
 *
 * ## System of equation modeling ## {#num_sk_inc_3D_petsc_sem}
 *
 * We model the equations. This part will change in the near future to use template expression
 * parsing and simplify the process of defining equations.
 *
 * \f$  \eta v_x \nabla v_x = eta\_lap\_vx \quad \nabla = \partial^{2}_{x} + \partial^{2}_{y} + \partial^{2}_{z}\f$ Step1
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
 * \f$ \partial^{forward}_{z} v_z = dz\_vz \f$ Step7
 *
 * \f$ dx\_vx + dy\_vy + dz_vz \f$ Step 8. This is the third equation in the system
 *
 *
 * \snippet Numerics/Stoke_flow/1_3D_incompressible/main_petsc.cpp def equation
 *
 */

//! \cond [def equation] \endcond

// Constant Field
struct eta
{
	//! indicate it is a constant field
	typedef void const_field;

	//! Eta is constant  one
	static float val()	{return 1.0;}
};

// Model the equations

constexpr unsigned int v[] = {0,1,2};
constexpr unsigned int P = 3;
constexpr unsigned int ic = 3;
constexpr int x = 0;
constexpr int y = 1;
constexpr int z = 2;

typedef Field<v[x],lid_nn> v_x;
typedef Field<v[y],lid_nn> v_y;
typedef Field<v[z],lid_nn> v_z;
typedef Field<P,lid_nn> Prs;

// Eq1 V_x

typedef mul<eta,Lap<v_x,lid_nn>,lid_nn> eta_lap_vx; // Step1
typedef D<x,Prs,lid_nn> p_x; // Step 2
typedef minus<p_x,lid_nn> m_p_x; // Step3
typedef sum<eta_lap_vx,m_p_x,lid_nn> vx_eq; // Step4

// Eq2 V_y

typedef mul<eta,Lap<v_y,lid_nn>,lid_nn> eta_lap_vy;
typedef D<y,Prs,lid_nn> p_y;
typedef minus<p_y,lid_nn> m_p_y;
typedef sum<eta_lap_vy,m_p_y,lid_nn> vy_eq;

// Eq3 V_z

typedef mul<eta,Lap<v_z,lid_nn>,lid_nn> eta_lap_vz;
typedef D<z,Prs,lid_nn> p_z;
typedef minus<p_z,lid_nn> m_p_z;
typedef sum<eta_lap_vz,m_p_z,lid_nn> vz_eq;

// Eq4 Incompressibility

typedef D<x,v_x,lid_nn,FORWARD> dx_vx; // Step5
typedef D<y,v_y,lid_nn,FORWARD> dy_vy; // Step6
typedef D<z,v_z,lid_nn,FORWARD> dz_vz; // Step 7
typedef sum<dx_vx,dy_vy,dz_vz,lid_nn> ic_eq; // Step8

//! \cond [def equation] \endcond

/*!
 * \page Stokes_1_3D_petsc Stokes incompressible 3D petsc
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
 *  order to finally reach a well defined system. Such observation can be extended In 3D
 *  leading to the following equations
 *
 * \snippet Numerics/Stoke_flow/1_3D_incompressible/main_petsc.cpp bond def eq
 *
 */

//! \cond [bond def eq] \endcond

// Directional Avg
typedef Avg<x,v_y,lid_nn> avg_x_vy;
typedef Avg<z,v_y,lid_nn> avg_z_vy;

typedef Avg<y,v_x,lid_nn> avg_y_vx;
typedef Avg<z,v_x,lid_nn> avg_z_vx;

typedef Avg<y,v_z,lid_nn> avg_y_vz;
typedef Avg<x,v_z,lid_nn> avg_x_vz;

typedef Avg<x,v_y,lid_nn,FORWARD> avg_x_vy_f;
typedef Avg<z,v_y,lid_nn,FORWARD> avg_z_vy_f;

typedef Avg<y,v_x,lid_nn,FORWARD> avg_y_vx_f;
typedef Avg<z,v_x,lid_nn,FORWARD> avg_z_vx_f;

typedef Avg<y,v_z,lid_nn,FORWARD> avg_y_vz_f;
typedef Avg<x,v_z,lid_nn,FORWARD> avg_x_vz_f;


#define EQ_1 0
#define EQ_2 1
#define EQ_3 2
#define EQ_4 3

//! \cond [bond def eq] \endcond

#include "Vector/vector_dist.hpp"
#include "data_type/aggregate.hpp"

int main(int argc, char* argv[])
{
	/*!
	 * \page Stokes_1_3D_petsc Stokes incompressible 3D petsc
	 *
	 * ## Initialization ## {#num_sk_inc_petsc_3D_init}
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
	 * \snippet Numerics/Stoke_flow/0_2D_incompressible/main_petsc.cpp init
	 *
	 */

	//! \cond [init] \endcond

	// Initialize
	openfpm_init(&argc,&argv);
	{
	// velocity in the grid is the property 0, pressure is the property 1
	constexpr int velocity = 0;
	constexpr int pressure = 1;

	// Domain
	Box<3,float> domain({0.0,0.0,0.0},{3.0,1.0,1.0});

	// Ghost (Not important in this case but required)
	Ghost<3,float> g(0.01);

	// Grid points on x=36 and y=12 z=12
	long int sz[] = {36,12,12};
	size_t szu[3];
	szu[0] = (size_t)sz[0];
	szu[1] = (size_t)sz[1];
	szu[2] = (size_t)sz[2];


	// We need one more point on the left and down part of the domain
	// This is given by the boundary conditions that we impose.
	//
	Padding<3> pd({1,1,1},{0,0,0});

	//! \cond [init] \endcond

	/*!
	 * \page Stokes_1_3D_petsc Stokes incompressible 3D petsc
	 *
	 * Distributed grid that store the solution
	 *
	 * \see \ref e0_s_grid_inst
	 *
	 * \snippet Numerics/Stoke_flow/1_3D_incompressible/main_petsc.cpp grid inst
	 *
	 */

	//! \cond [grid inst] \endcond

	grid_dist_id<3,float,aggregate<float[3],float>> g_dist(szu,domain,g);

	//! \cond [grid inst] \endcond

	/*!
	 * \page Stokes_1_3D_petsc Stokes incompressible 3D petsc
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
	 * \snippet Numerics/Stoke_flow/1_3D_incompressible/main_petsc.cpp fd scheme
	 *
	 */

	//! \cond [fd scheme] \endcond

	// It is the maximum extension of the stencil (order 2 laplacian stencil has extension 1)
	Ghost<3,long int> stencil_max(1);

	// Finite difference scheme
	FDScheme<lid_nn> fd(pd, stencil_max, domain, g_dist);

	//! \cond [fd scheme] \endcond

	/*!
	 * \page Stokes_1_3D_petsc Stokes incompressible 3D petsc
	 *
	 * ## Impose the equation on the domain ## {#num_sk_inc_3D_ied}
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
	 * \snippet Numerics/Stoke_flow/1_3D_incompressible/main_petsc.cpp impose eq dom
	 *
	 *
	 */

	//! \cond [impose eq dom] \endcond

	// start and end of the bulk

	fd.impose(ic_eq(),0.0, EQ_4, {0,0,0},{sz[0]-2,sz[1]-2,sz[2]-2},true);
	fd.impose(Prs(),  0.0, EQ_4, {0,0,0},{0,0,0});
	fd.impose(vx_eq(),0.0, EQ_1, {1,0},{sz[0]-2,sz[1]-2,sz[2]-2});
	fd.impose(vy_eq(),0.0, EQ_2, {0,1},{sz[0]-2,sz[1]-2,sz[2]-2});
	fd.impose(vz_eq(),0.0, EQ_3, {0,0,1},{sz[0]-2,sz[1]-2,sz[2]-2});

	// v_x
	// R L (Right,Left)
	fd.impose(v_x(),0.0, EQ_1, {0,0,0},      {0,sz[1]-2,sz[2]-2});
	fd.impose(v_x(),0.0, EQ_1, {sz[0]-1,0,0},{sz[0]-1,sz[1]-2,sz[2]-2});

	// T B (Top,Bottom)
	fd.impose(avg_y_vx_f(),0.0, EQ_1, {0,-1,0},     {sz[0]-1,-1,sz[2]-2});
	fd.impose(avg_y_vx(),0.0, EQ_1,   {0,sz[1]-1,0},{sz[0]-1,sz[1]-1,sz[2]-2});

	// A F (Forward,Backward)
	fd.impose(avg_z_vx_f(),0.0, EQ_1, {0,-1,-1},     {sz[0]-1,sz[1]-1,-1});
	fd.impose(avg_z_vx(),0.0, EQ_1, {0,-1,sz[2]-1},{sz[0]-1,sz[1]-1,sz[2]-1});

	// v_y
	// R L
	fd.impose(avg_x_vy_f(),0.0, EQ_2,  {-1,0,0},     {-1,sz[1]-1,sz[2]-2});
	fd.impose(avg_x_vy(),1.0, EQ_2,    {sz[0]-1,0,0},{sz[0]-1,sz[1]-1,sz[2]-2});

	// T B
	fd.impose(v_y(), 0.0, EQ_2, {0,0,0},      {sz[0]-2,0,sz[2]-2});
	fd.impose(v_y(), 0.0, EQ_2, {0,sz[1]-1,0},{sz[0]-2,sz[1]-1,sz[2]-2});

	// F A
	fd.impose(avg_z_vy(),0.0, EQ_2,   {-1,0,sz[2]-1}, {sz[0]-1,sz[1]-1,sz[2]-1});
	fd.impose(avg_z_vy_f(),0.0, EQ_2, {-1,0,-1},      {sz[0]-1,sz[1]-1,-1});

	// v_z
	// R L
	fd.impose(avg_x_vz_f(),0.0, EQ_3, {-1,0,0},     {-1,sz[1]-2,sz[2]-1});
	fd.impose(avg_x_vz(),1.0, EQ_3,   {sz[0]-1,0,0},{sz[0]-1,sz[1]-2,sz[2]-1});

	// T B
	fd.impose(avg_y_vz(),0.0, EQ_3, {-1,sz[1]-1,0},{sz[0]-1,sz[1]-1,sz[2]-1});
	fd.impose(avg_y_vz_f(),0.0, EQ_3, {-1,-1,0},   {sz[0]-1,-1,sz[2]-1});

	// F A
	fd.impose(v_z(),0.0, EQ_3, {0,0,0},      {sz[0]-2,sz[1]-2,0});
	fd.impose(v_z(),0.0, EQ_3, {0,0,sz[2]-1},{sz[0]-2,sz[1]-2,sz[2]-1});

	// When we pad the grid, there are points of the grid that are not
	// touched by the previous condition. Mathematically this lead
	// to have too many variables for the conditions that we are imposing.
	// Here we are imposing variables that we do not touch to zero
	//

	// L R
	fd.impose(Prs(), 0.0, EQ_4, {-1,-1,-1},{-1,sz[1]-1,sz[2]-1});
	fd.impose(Prs(), 0.0, EQ_4, {sz[0]-1,-1,-1},{sz[0]-1,sz[1]-1,sz[2]-1});

	// T B
	fd.impose(Prs(), 0.0, EQ_4, {0,sz[1]-1,-1}, {sz[0]-2,sz[1]-1,sz[2]-1});
	fd.impose(Prs(), 0.0, EQ_4, {0,-1     ,-1}, {sz[0]-2,-1,     sz[2]-1});

	// F A
	fd.impose(Prs(), 0.0, EQ_4, {0,0,sz[2]-1}, {sz[0]-2,sz[1]-2,sz[2]-1});
	fd.impose(Prs(), 0.0, EQ_4, {0,0,-1},      {sz[0]-2,sz[1]-2,-1});

	// Impose v_x  v_y v_z padding
	fd.impose(v_x(), 0.0, EQ_1, {-1,-1,-1},{-1,sz[1]-1,sz[2]-1});
	fd.impose(v_y(), 0.0, EQ_2, {-1,-1,-1},{sz[0]-1,-1,sz[2]-1});
	fd.impose(v_z(), 0.0, EQ_3, {-1,-1,-1},{sz[0]-1,sz[1]-1,-1});

	//! \cond [impose eq dom] \endcond

	/*!
	 * \page Stokes_1_3D_petsc Stokes incompressible 3D petsc
	 *
	 * ## Solve the system of equation ## {#num_sk_inc_3D_petsc_sse}
	 *
	 * Once we imposed all the equations we can retrieve the Matrix A and the vector b
	 * and pass these two element to the solver. In this example we are using  PETSC solvers
	 *  direct/Iterative solvers. While Umfpack
	 * has only one solver, PETSC wrap several solvers. The function best_solve set the solver in
	 * the modality to try multiple solvers to solve your system. The subsequent call to solve produce a report
	 * of all the solvers tried comparing them in error-convergence and speed. If you do not use
	 * best_solve try to solve your system with the default solver GMRES (That is the most robust iterative solver
	 *  method)
	 *
	 * \snippet Numerics/Stoke_flow/1_3D_incompressible/main_petsc.cpp solver
	 *
	 */

	//! \cond [solver] \endcond

	// Create an PETSC solver
	petsc_solver<double> solver;

	// Warning try many solver and collect statistics require a lot of time
	// To just solve you can comment this line
//	solver.best_solve();

	// Give to the solver A and b, return x, the solution
	auto x = solver.solve(fd.getA(),fd.getB());

	//! \cond [solver] \endcond

	/*!
	 * \page Stokes_1_3D_petsc Stokes incompressible 3D petsc
	 *
	 * ## Copy the solution on the grid and write on VTK ## {#num_sk_inc_3D_petsc_csg}
	 *
	 * Once we have the solution we copy it on the grid
	 *
	 * \snippet Numerics/Stoke_flow/1_3D_incompressible/main_petsc.cpp copy write
	 *
	 */

	//! \cond [copy write] \endcond

	// copy the solution to grid
	fd.template copy<velocity,pressure>(x,{0,0},{sz[0]-1,sz[1]-1,sz[2]-1},g_dist);

	g_dist.write("lid_driven_cavity_p_petsc");

	//! \cond [copy write] \endcond


	/*!
	 * \page Stokes_1_3D_petsc Stokes incompressible 3D petsc
	 *
	 * ## Finalize ## {#num_sk_inc_3D_petsc_fin}
	 *
	 *  At the very end of the program we have always to de-initialize the library
	 *
	 * \snippet Numerics/Stoke_flow/1_3D_incompressible/main_petsc.cpp fin lib
	 *
	 */

	//! \cond [fin lib] \endcond
	}
	openfpm_finalize();

	//! \cond [fin lib] \endcond


	/*!
	 * \page Stokes_1_3D_petsc Stokes incompressible 3D petsc
	 *
	 * # Full code # {#num_sk_inc_3D_petsc_code}
	 *
	 * \include Numerics/Stoke_flow/1_3D_incompressible/main_petsc.cpp
	 *
	 */
}
