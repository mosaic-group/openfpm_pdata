#include "Grid/grid_dist_id.hpp"
#include "data_type/aggregate.hpp"
#include "timer.hpp"

/*!
 *
 * \page Grid_3_gs_3D Gray Scott in 3D
 *
 * # Solving a gray scott-system in 3D # {#e3_gs_gray_scott}
 *
 * This example is just an extension of the 2D Gray scott example.
 * Here we show how to solve a non-linear reaction diffusion system in 3D
 *
 * In figure is the final solution of the problem
 *
 * \htmlonly
 * <img src="http://ppmcore.mpi-cbg.de/web/images/examples/gray_scott_3d/gs_alpha.png"/>
 * \endhtmlonly
 *
 * \see \ref Grid_2_solve_eq
 *
 * \snippet Grid/3_gray_scott/main.cpp constants
 * 
 */

//! \cond [constants] \endcond

constexpr int U = 0;
constexpr int V = 1;

constexpr int x = 0;
constexpr int y = 1;
constexpr int z = 2;

//! \cond [constants] \endcond

void init(grid_dist_id<3,double,aggregate<double,double> > & Old, grid_dist_id<3,double,aggregate<double,double> > & New, Box<3,double> & domain)
{
	auto it = Old.getDomainIterator();

	while (it.isNext())
	{
		// Get the local grid key
		auto key = it.get();

		// Old values U and V
		Old.template get<U>(key) = 1.0;
		Old.template get<V>(key) = 0.0;

		// Old values U and V
		New.template get<U>(key) = 0.0;
		New.template get<V>(key) = 0.0;

		++it;
	}

	long int x_start = Old.size(0)*1.55f/domain.getHigh(0);
	long int y_start = Old.size(1)*1.55f/domain.getHigh(1);
	long int z_start = Old.size(1)*1.55f/domain.getHigh(2);

	long int x_stop = Old.size(0)*1.85f/domain.getHigh(0);
	long int y_stop = Old.size(1)*1.85f/domain.getHigh(1);
	long int z_stop = Old.size(1)*1.85f/domain.getHigh(2);

	grid_key_dx<3> start({x_start,y_start,z_start});
	grid_key_dx<3> stop ({x_stop,y_stop,z_stop});
	auto it_init = Old.getSubDomainIterator(start,stop);

	while (it_init.isNext())
	{
		auto key = it_init.get();

                Old.template get<U>(key) = 0.5 + (((double)std::rand())/RAND_MAX -0.5)/10.0;
                Old.template get<V>(key) = 0.25 + (((double)std::rand())/RAND_MAX -0.5)/20.0;

		++it_init;
	}
}


int main(int argc, char* argv[])
{
	openfpm_init(&argc,&argv);

	// domain
	Box<3,double> domain({0.0,0.0},{2.5,2.5,2.5});
	
	// grid size
        size_t sz[3] = {128,128,128};

	// Define periodicity of the grid
	periodicity<3> bc = {PERIODIC,PERIODIC,PERIODIC};
	
	// Ghost in grid unit
	Ghost<3,long int> g(1);
	
	// deltaT
	double deltaT = 1;

	// Diffusion constant for specie U
	double du = 2*1e-5;

	// Diffusion constant for specie V
	double dv = 1*1e-5;

	// Number of timesteps
        size_t timeSteps = 5000;

	// K and F (Physical constant in the equation)
        double K = 0.053;
        double F = 0.014;

	//! \cond [init lib] \endcond

	/*!
	 * \page Grid_3_gs_3D Gray Scott in 3D
	 *
	 * Here we create 2 distributed grid in 2D Old and New. In particular because we want that
	 * the second grid is distributed across processors in the same way we pass the decomposition
	 * of the Old grid to the New one in the constructor with **Old.getDecomposition()**. Doing this,
	 * we force the two grid to have the same decomposition.
	 *
	 * \snippet Grid/3_gray_scott/main.cpp init grid
	 *
	 */

	//! \cond [init grid] \endcond

	grid_dist_id<3, double, aggregate<double,double>> Old(sz,domain,g,bc);

	// New grid with the decomposition of the old grid
        grid_dist_id<3, double, aggregate<double,double>> New(Old.getDecomposition(),sz,g);

	
	// spacing of the grid on x and y
	double spacing[3] = {Old.spacing(0),Old.spacing(1),Old.spacing(2)};

	init(Old,New,domain);

	// sync the ghost
	size_t count = 0;
	Old.template ghost_get<U,V>();

	// because we assume that spacing[x] == spacing[y] we use formula 2
	// and we calculate the prefactor of Eq 2
	double uFactor = deltaT * du/(spacing[x]*spacing[x]);
	double vFactor = deltaT * dv/(spacing[x]*spacing[x]);

	timer tot_sim;
	tot_sim.start();

	grid_key_dx<3> star_stencil_3D[7] = {{0,0,0},
	                                         {0,0,-1},
											 {0,0,1},
											 {0,-1,0},
											 {0,1,0},
											 {-1,0,0},
											 {1,0,0}};

	for (size_t i = 0; i < timeSteps; ++i)
	{
		if (i % 300 == 0)
			std::cout << "STEP: " << i << std::endl;

		auto it = Old.getDomainIteratorStencil(star_stencil_3D);;

		while (it.isNext())
		{
			// center point
			auto Cp = it.getStencil<0>();

			// plus,minus X,Y,Z
			auto mx = it.getStencil<1>();
			auto px = it.getStencil<2>();
			auto my = it.getStencil<3>();
			auto py = it.getStencil<4>();
			auto mz = it.getStencil<5>();
			auto pz = it.getStencil<6>();

			std::cout << &New.get<U>(Cp) << std::endl;


			// update based on Eq 2
			New.get<U>(Cp) = Old.get<U>(Cp) + uFactor * (
										Old.get<U>(px) +
										Old.get<U>(mx) +
										Old.get<U>(py) +
										Old.get<U>(my) +
										Old.get<U>(pz) +
										Old.get<U>(mz) -
										6.0*Old.get<U>(Cp)) +
										- deltaT * Old.get<U>(Cp) * Old.get<V>(Cp) * Old.get<V>(Cp) +
										- deltaT * F * (Old.get<U>(Cp) - 1.0);


			// update based on Eq 2
			New.get<V>(Cp) = Old.get<V>(Cp) + vFactor * (
										Old.get<V>(px) +
										Old.get<V>(mx) +
										Old.get<V>(py) +
										Old.get<V>(my) +
										Old.get<V>(pz) +
                                        Old.get<V>(mz) -
										6*Old.get<V>(Cp)) +
										deltaT * Old.get<U>(Cp) * Old.get<V>(Cp) * Old.get<V>(Cp) +
										- deltaT * (F+K) * Old.get<V>(Cp);

			// Next point in the grid
			++it;
		}

		// Here we copy New into the old grid in preparation of the new step
		// It would be better to alternate, but using this we can show the usage
		// of the function copy. To note that copy work only on two grid of the same
		// decomposition. If you want to copy also the decomposition, or force to be
		// exactly the same, use Old = New
		Old.copy(New);

		// After copy we synchronize again the ghost part U and V
		Old.ghost_get<U,V>();

		// Every 30 time step we output the configuration for
		// visualization
		if (i % 60 == 0)
		{
			Old.write_frame("output",count,VTK_WRITER | FORMAT_BINARY);
			count++;
		}
	}
	
	tot_sim.stop();
	std::cout << "Total simulation: " << tot_sim.getwct() << std::endl;

	//! \cond [time stepping] \endcond

	/*!
	 * \page Grid_3_gs_3D Gray Scott in 3D
	 *
	 * ## Finalize ##
	 *
	 * Deinitialize the library
	 *
	 * \snippet Grid/3_gray_scott/main.cpp finalize
	 *
	 */

	//! \cond [finalize] \endcond

	openfpm_finalize();

	//! \cond [finalize] \endcond
}
