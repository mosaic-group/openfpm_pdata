#include "Grid/grid_dist_id.hpp"
#include "data_type/scalar.hpp"
#include "Decomposition/CartDecomposition.hpp"

/*
 * ### WIKI 1 ###
 *
 * ## Simple example
 * 
 * This example show several basic functionalities of the distributed grid
 * 
 * ### WIKI END ###
 * 
 */


int main(int argc, char* argv[])
{
	//
	// ### WIKI 2 ###
	//
	// Initialize the library and several objects 
	//
	init_global_v_cluster(&argc,&argv);
	
	//
	// ### WIKI 3 ###
	//
	// Create several object needed later, in particular
	// * A 3D box that define the domain
	// * an array of 3 unsigned integer that define the size of the grid on each dimension
	// * A Ghost object that will define the extension of the ghost part for each sub-domain in physical units
	
	Box<3,float> domain({0.0,0.0,0.0},{1.0,1.0,1.0});
	size_t sz[3];
	sz[0] = 100;
	sz[1] = 100;
	sz[2] = 100;
	
	// Ghost
	Ghost<3,float> g(0.01);
	
	//
	// ### WIKI 4 ###
	//
	// Create a distributed grid in 3D (1° template parameter) defined in R^3 with float precision (2° template parameter)
	// using a CartesianDecomposition strategy (3° parameter) (the parameter 1° and 2° inside CartDecomposition must match 1° and 2°
	// of grid_dist_id)
	//
	// Constructor parameters:
	//
	// * sz: size of the grid on each dimension
	// * domain: where the grid is defined
	// * g: ghost extension
	//
	grid_dist_id<3, float, scalar<float[3]>, CartDecomposition<3,float>> g_dist(sz,domain,g);
	
	// ### WIKI 5 ###
	//
	// Get an iterator that go throught the point of the domain (No ghost)
	//
	auto dom = g_dist.getDomainIterator();

	// ### WIKI END ###
	
	size_t count = 0;
	
	// Iterate over all the points
	while (dom.isNext())
	{
		//
		// ### WIKI 6 ###
		//
		// Get the local grid key, the local grid key store internaly the sub-domain id (each sub-domain contain a grid)
		// and the local grid point id identified by 2 integers in 2D 3 integer in 3D and so on. These two dinstinc element are
		// available with key.getSub() and key.getKey()
		//
		auto key = dom.get();
		
		//
		// ### WIKI 7 ###
		//
		// Here we convert the local grid position, into global position, key_g contain 3 integers that identify the position
		// of the grid point in global coordinates
		//
		//
		auto key_g = g_dist.getGKey(key);

		//
		// ### WIKI 8 ###
		//
		// we write on the grid point of position (i,j,k) the value i*i + j*j + k*k on the component [0] of the vector
		g_dist.template get<0>(key)[0] = key_g.get(0)*key_g.get(0) + key_g.get(1)*key_g.get(1) + key_g.get(2)*key_g.get(2);
		
		// ### WIKI END ###
		
		// Count the points
		count++;

		//
		// ### WIKI 9 ###
		//
		// next point
		++dom;
		
		// ### WIKI END ###
	}

	//
	// ### WIKI 10 ###
	//
	// Each sub-domain has an extended part, that is materially contained from another processor that in general is not synchronized
	// ghost_get<0> synchronize the property 0 (the vector) in the ghost part
	//
	//
	g_dist.template ghost_get<0>();
	
	//
	// ### WIKI 11 ###
	//
	// count contain the number of points the local processor contain, if we are interested to count the total number across the processor
	// we can use the function add, to sum across processors. First we have to get an instance of Vcluster, queue an operation of add with
	// the variable count and finaly execute. All the operation are asynchronous, execute work like a barrier and ensure that all the 
	// queued operations are executed
	//
	Vcluster & vcl = g_dist.getVC();

	vcl.sum(count);
	vcl.execute();
	
	// only master output
	if (vcl.getProcessUnitID() == 0)
	  std::cout << "Number of points: " << count << "\n";

	//
	// ### WIKI 12 ###
	//
	// Finally we want a nice output to visualize the information stored by the distributed grid
	//
	g_dist.write("output");
	
	//
	// ### WIKI 13 ###
	//
	// For debugging porpouse and demostration we output the decomposition
	//
	g_dist.getDecomposition().write("out_dec");
	
	//
	// ### WIKI 14 ###
	//
	// Deinitialize the library
	//
	delete_global_v_cluster();
}
