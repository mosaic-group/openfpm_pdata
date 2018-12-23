/*
 * ### WIKI 1 ###
 *
 * ## Security enhancement example
 *
 * This example show several basic functionalities of Security Enhancements
 *
 */

#define SE_CLASS1
#define SE_CLASS2
#define SE_CLASS3
#define THROW_ON_ERROR
#include "Memleak_check.hpp"
#include "Vector/vector_dist.hpp"
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
	// With print_unalloc we can check how much memory has been allocated and which structure
	// has been allocated, initially there are not
	//

	std::cout << "Allocated memory before initializing \n";
    print_alloc();
    std::cout << "\n";
    std::cout << "\n";
    std::cout << "\n";

	//
	// ### WIKI 3 ###
	//
	// Here we Initialize the library, than we create a uniform random generator between 0 and 1 to to generate particles
	// randomly in the domain, we create a Box that define our domain, boundary conditions and ghost
	//
	openfpm_init(&argc,&argv);
	Vcluster<> & v_cl = create_vcluster();
	
	typedef Point<2,float> s;

	Box<2,float> box({0.0,0.0},{1.0,1.0});
        size_t bc[2]={NON_PERIODIC,NON_PERIODIC};
	Ghost<2,float> g(0.01);

	//
	// ### WIKI 4 ###
	//
	// Here we ask again for the used memory, as we can see Vcluster and several other structures encapsulated inside
	// Vcluster register itself
	//
	if (v_cl.getProcessUnitID() == 0)
	{
		std::cout << "Allocated memory after initialization \n";
		print_alloc();
		std::cout << "\n";
		std::cout << "\n";
		std::cout << "\n";
	}

	//
	// ### WIKI 5 ###
	//
	// Here we are creating a distributed vector defined by the following parameters
	// 
	// * Dimensionality of the space where the objects live 2D (1° template parameters)
	// * Type of the space, float (2° template parameters)
	// * Information stored by each object (3* template parameters), in this case a Point_test store 4 scalars
	//   1 vector and an asymmetric tensor of rank 2
	// * Strategy used to decompose the space
	// 
	// The Constructor instead require:
	//
	// * Number of particles 4096 in this case
	// * Domain where is defined this structure
	//
	// The following construct a vector where each processor has 4096 / N_proc (N_proc = number of processor)
	// objects with an undefined position in space. This non-space decomposition is also called data-driven
	// decomposition
	//
	{
		vector_dist<2,float, Point_test<float> > vd(4096,box,bc,g);

		//
		// ### WIKI 6 ###
		//
		// we create a key that for sure overflow the local datastructure, 2048 with this program is started with more than 3
		// processors, try and catch are optionals in case you want to recover from a buffer overflow
		//
		try
        {
			vect_dist_key_dx vt(5048);
			auto it = vd.getPos(vt);
        }
		catch (size_t e)
		{
			std::cerr << "Error notification of overflow \n";
		}
	}
	//
	// ### WIKI 7 ###
	//
	// At this point the vector went out of the scope and if destroyed
	// we create, now two of them using new
	//

	vector_dist<2,float, Point_test<float> > * vd1 = new vector_dist<2,float, Point_test<float>, CartDecomposition<2,float> >(4096,box,bc,g);
	vector_dist<2,float, Point_test<float> > * vd2 = new vector_dist<2,float, Point_test<float>, CartDecomposition<2,float> >(4096,box,bc,g);

	//
	// ### WIKI 8 ###
	//
	// we can check that these two structure produce an explicit allocation checking
	// for registered pointers and structures with print_alloc, in the list we see 2 additional
	// entry for distributed vector in yellow, pdata to work use the data structures that register
	// itself in magenta, the same things happen for the real memory allocation from devices in
	// fully green
	//

	if (v_cl.getProcessUnitID() == 0)
	{
		std::cout << "Allocated memory with 2 vectors \n";
		print_alloc();
		std::cout << "\n";
		std::cout << "\n";
		std::cout << "\n";
	}

	//
	// ### WIKI 9 ###
	//
	// we can also ask to the structure to identify their-self in the list
	//

    std::cout << "Vector id: " << vd1->who() << "\n";
    std::cout << "Vector id: " << vd2->who() << "\n";

	//
	// ### WIKI 10 ###
	//
	// delete vd1 and print allocated memory, one distributed vector disappear
	//

	delete vd1;

	if (v_cl.getProcessUnitID() == 0)
	{
		std::cout << "Allocated memory with 1 vector \n";
		print_alloc();
    	std::cout << "\n";
    	std::cout << "\n";
    	std::cout << "\n";
	}

	//
	// ### WIKI 11 ###
	//
	// delete vd2 and print allocated memory, all distributed vector de-register
	//

	delete vd2;

	if (v_cl.getProcessUnitID() == 0)
	{
		std::cout << "Allocated memory with 1 vector \n";
		print_alloc();
		std::cout << "\n";
		std::cout << "\n";
		std::cout << "\n";
	}

	//
	// ### WIKI 12 ###
	//
	// Try to use a deleted object
	//
	try
    {
		vect_dist_key_dx vt(0);
		auto it = vd1->getPos(vt);
    }
	catch (std::exception & e)
	{
		std::cerr << "Error notification of invalid usage of deleted object \n";
	}

	//
	// ### WIKI 13 ###
	//
	// Deinitialize the library
	//
	openfpm_finalize();
}

