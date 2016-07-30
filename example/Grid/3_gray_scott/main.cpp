#include "Grid/grid_dist_id.hpp"
#include "data_type/aggregate.hpp"
#include "timer.hpp"

/*
 * ### WIKI 1 ###
 *
 * ## Simple example of grid usage for solving gray scott-system
 *
 * This example show the usage of periodic grid with ghost part given in grid units to solve
 * the following system of equations
 *
 * \frac{\partial u}{\partial t} = D_u \laplacian u -uv^2 + F(1-u)
 * \frac{\partial v}{\partial t} = D_v \laplacian v -uv^2 - (F + k)v
 * 
 * ### WIKI END ###
 * 
 */

constexpr int U = 0;
constexpr int V = 1;


/*
 * ### WIKI 2 ###
 *
 * Initialize the field U
 *
 */
void init(grid_dist_id<2,double,aggregate<double,double> > & Old, grid_dist_id<2,double,aggregate<double,double> > & New, Box<2,double> & domain)
{
	//
	// ### WIKI ###
	//
	// Here we initialize the full domain + boundary conditions
	//
	auto it = Old.getDomainGhostIterator();

	while (it.isNext())
	{
		//
		// ### WIKI 6 ###
		//
		// Get the local grid key
		//
		auto key = it.get();

		// Old values U and V
		Old.template get<U>(key) = 1.0;
		Old.template get<V>(key) = 0.0;

		// Old values U and V
		New.template get<U>(key) = 0.0;
		New.template get<V>(key) = 0.0;

		++it;
	}

	//
	// ### WIKI 20 ###
	//
	// We initialize a box in the center of the domain with different values
	//
	grid_key_dx<2> start({(long int)std::floor(Old.size(0)*1.55f/domain.getHigh(0)),(long int)std::floor(Old.size(1)*1.55f/domain.getHigh(1))});
	grid_key_dx<2> stop ({(long int)std::ceil (Old.size(0)*1.85f/domain.getHigh(0)),(long int)std::ceil (Old.size(1)*1.85f/domain.getHigh(1))});
	auto it_init = Old.getSubDomainIterator(start,stop);

	while (it_init.isNext())
	{
		auto key = it_init.get();

		Old.template get<U>(key) = 0.5 + (((double)std::rand())/RAND_MAX -0.5)/100.0;
		Old.template get<V>(key) = 0.25 + (((double)std::rand())/RAND_MAX -0.5)/200.0;

		++it_init;
	}

	//
	// ### WIKI END ###

}

//
// ### WIKI 4 ###
//
// Usefull constant
//
constexpr int x = 0;
constexpr int y = 1;

int main(int argc, char* argv[])
{
	//
	// ### WIKI 2 ###
	//
	// Initialize the library
	//
	openfpm_init(&argc,&argv);
	
	//
	// ### WIKI 3 ###
	//
	// Create
	// * A 2D box that define the domain
	// * an array of 2 unsigned integer that will define the size of the grid on each dimension
	// * A Ghost object that will define the extension of the ghost part for each sub-domain in grid point unit

	Box<2,double> domain({0.0,0.0},{2.5,2.5});
	size_t sz[2];
	sz[0] = 128;
	sz[1] = 128;
	
	// Define periodicity of the grid
	periodicity<2> bc = {PERIODIC,PERIODIC};
	
	// Ghost in grid unit
	Ghost<2,long int> g(1);
	
	// deltaT
	double deltaT = 1;

	// Diffusion constant for specie U
	double du = 2*1e-5;

	// Diffusion constant for specie V
	double dv = 1*1e-5;

	// Number of timesteps
	size_t timeSteps = 15000;

	// K and F (Physical constant in the equation)
	double K = 0.055;
	double F = 0.03;


	//
	// ### WIKI 4 ###
	//
	// Create a distributed grid in 2D (1° template parameter) space in with double precision (2° template parameter)
	// each grid point contain a scalar (double),
	//
	// Constructor parameters:
	//
	// * sz: size of the grid on each dimension
	// * domain: where the grid is defined
	// * g: ghost extension
	// * bc: boundary conditions
	//
	grid_dist_id<2, double, aggregate<double,double>> Old(sz,domain,g,bc);
	grid_dist_id<2, double, aggregate<double,double>> New(Old.getDecomposition(),sz,domain,g);

	
	// spacing
	double spacing[2] = {Old.spacing(0),Old.spacing(1)};

	//
	// ### WIKI 6 ###
	//
	// Initialize U and fill the boundary conditions
	//

	init(Old,New,domain);

	timer time;
	time.start();

	// ### WIKI 7 ###
	//
	// Do 10000 iteration of Red-Black Gauss-Siedel
	//

	// sync the ghost
	size_t count;
	Old.template ghost_get<0,1>();

	double uFactor = deltaT * du/(spacing[0]*spacing[0]);
	double vFactor = deltaT * dv/(spacing[0]*spacing[0]);

	for (size_t i = 0; i < timeSteps; ++i)
	{
		auto it = Old.getDomainIterator();

		while (it.isNext())
		{
			auto key = it.get();

			New.get<U>(key) = Old.get<U>(key) + uFactor * (
										Old.get<U>(key.move(x,1)) +
										Old.get<U>(key.move(x,-1)) +
										Old.get<U>(key.move(y,1)) +
										Old.get<U>(key.move(y,-1)) +
										-4.0*Old.get<U>(key)) +
										- deltaT * Old.get<U>(key) * Old.get<V>(key) * Old.get<V>(key) +
										- deltaT * F * (Old.get<U>(key) - 1.0);


			New.get<V>(key) = Old.get<V>(key) + vFactor * (
										Old.get<V>(key.move(x,1)) +
										Old.get<V>(key.move(x,-1)) +
										Old.get<V>(key.move(y,1)) +
										Old.get<V>(key.move(y,-1)) -
										4*Old.get<V>(key)) +
										deltaT * Old.get<U>(key) * Old.get<V>(key) * Old.get<V>(key) +
										- deltaT * (F+K) * Old.get<V>(key);

			++it;
		}

		Old.copy(New);
		Old.ghost_get<0,1>();

		if (i % 100 == 0)
		{
			Old.write("output",count);
			count++;
		}
	}
	
	//
	// ### WIKI 14 ###
	//
	// Deinitialize the library
	//
	openfpm_finalize();
}
