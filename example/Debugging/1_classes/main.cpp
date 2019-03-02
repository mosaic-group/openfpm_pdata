/*
 * ### WIKI 1 ###
 *
 * ## Security enhancement example
 *
 * This example show how to see where an allocation or corruption happen offline and online.
 * Every time an error occur, the library output where the detection happen filename and line,
 *  in order to debug, there is an online option and an offline option
 *
 *  * online: put a breakpoint on the indicated line with your preferred debugger
 *  * offline: set ulimit -c unlimited to activate the core dump file and open the core dump with your debugger
 *
 */

#define SE_CLASS1
#define SE_CLASS2
// SE_CLASS2 is unsupported if not used in combination with SE_CLASS2_TRACK_ONLY
#define SE_CLASS2_ONLY_TRACK
#define SE_CLASS3
#define PRINT_STACKTRACE
#define THROW_ON_ERROR
#include "Memleak_check.hpp"
#include "Grid/grid_dist_id.hpp"
#include "Decomposition/CartDecomposition.hpp"
#include "Point_test.hpp"

/*
 * ### WIKI END ###
 */


int main(int argc, char* argv[])
{
	//
	// ### WIKI 2 ###
	//
	// Here we Initialize the library,
	// * message_on_allocation set a message to print when one allocation is reached, the filename and line number can be used to set a breakpoint and analyze the stacktrace.
	// * throw_on_allocation throw when one allocation is reached, producing the termination of the program and a core dump (if no try catch is set-up)
	//
	openfpm_init(&argc,&argv);
	Vcluster<> & v_cl = create_vcluster();

	throw_on_alloc(10);
	// message_on_alloc(10);

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
	//
	grid_dist_id<3, float, aggregate<float[3]>, CartDecomposition<3,float>> * g_dist = new grid_dist_id<3, float, aggregate<float[3]> >(sz,domain,g);

	//
	// ### WIKI 6 ###
	//
	// print allocated structures
	//

	if (v_cl.getProcessUnitID() == 0)
	{print_alloc();}

	//
	// ### WIKI 5 ###
	//
	// delete g_dist
	//

	delete g_dist;

	//
	// ### WIKI 6 ###
	//
	// On purpose we try to access a deleted object
	//

	g_dist->getGridInfo();

	//
	// ### WIKI 13 ###
	//
	// Deinitialize the library
	//
	openfpm_finalize();
}

