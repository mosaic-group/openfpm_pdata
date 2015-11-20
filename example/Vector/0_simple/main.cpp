#include "Vector/vector_dist.hpp"
#include "Decomposition/CartDecomposition.hpp"

/*
 * ### WIKI 1 ###
 *
 * ## Simple example
 * 
 * This example show several basic functionalities of the distributed vector
 * 
 * ### WIKI END ###
 * 
 */

/*
 * ### WIKI 2 ###
 *
 * We define a particle structure it contain 4 scalars one vector with 3 components
 * and a tensor of rank 2 3x3
 *
 * ### WIKI END ###
 *
 */

template<typename T> class Particle
{
public:

	typedef boost::fusion::vector<T,T,T,T,T[3],T[3][3]> type;

	type data;

	static const unsigned int x = 0;
	static const unsigned int y = 1;
	static const unsigned int z = 2;
	static const unsigned int s = 3;
	static const unsigned int v = 4;
	static const unsigned int t = 5;
	static const unsigned int max_prop = 6;
};

int main(int argc, char* argv[])
{
	//
	// ### WIKI 2 ###
	//
	// Here we Initialize the library, than we create a uniform random generator between 0 and 1 to to generate particles
	// randomly in the domain, we create a Box that define our domain
	//
	init_global_v_cluster(&argc,&argv);
	Vcluster & v_cl = *global_v_cluster;
	
	typedef Point<2,float> s;

	// set the seed
	// create the random generator engine
	std::srand(v_cl.getProcessUnitID());
	std::default_random_engine eg;
	std::uniform_real_distribution<float> ud(0.0f, 1.0f);

	Box<2,float> box({0.0,0.0},{1.0,1.0});
	
	//
	// ### WIKI 3 ###
	//
	// Here we are creating a distributed vector defined by the following parameters
	// 
	// * Dimensionality of the space where the objects live 2D (1° template parameters)
	// * Type of the space, float (2° template parameters)
	// * Information stored by each object (3* template parameters), in this case a Point_test store 4 scalars
	//   1 vector and an asymmetric tensor of rank 2
	// * Strategy used to decompose the space
	// 
	// Constructor instead require:
	//
	// * Number of particles 4096 in this case
	// * Domain where is defined this structure
	//
	// The following construct a vector where each processor has 4096 / N_proc (N_proc = number of processor)
	// objects with an undefined position in space. This non-space decomposition is also called data-driven
	// decomposition
	//
	vector_dist<2,float, Particle<float>, CartDecomposition<2,float> > vd(4096,box);

	//
	// ### WIKI 5 ###
	//
	// Get an iterator that go throught the objects, in an undefined position state and define its position
	//
	auto it = vd.getIterator();

	while (it.isNext())
	{
		auto key = it.get();

		vd.template getPos<s::x>(key)[0] = ud(eg);
		vd.template getPos<s::x>(key)[1] = ud(eg);

		++it;
	}

	//
	// ### WIKI 6 ###
	//
	// Once we define the position, we distribute them according to the default decomposition
	// The default decomposition is created even before assigning the position to the object
	// (This will probably change in future)
	//
	vd.map();

	//
	// ### WIKI 7 ###
	//
	// We get the object that store the decomposition, than we iterate again across all the objects, we count them
	// and we confirm that all the particles are local
	//
	size_t cnt = 0;
	auto & ct = vd.getDecomposition();
	it = vd.getIterator();

	while (it.isNext())
	{
		auto key = it.get();

		
		if (ct.isLocal(vd.template getPos<s::x>(key)) == false)
			std::cerr << "Error particle is not local" << "\n";

		cnt++;

		++it;
	}

	//
	// ### WIKI 8 ###
	//
	// cnt contain the number of object the local processor contain, if we are interested to count the total number across the processor
	// we can use the function add, to sum across processors. First we have to get an instance of Vcluster, queue an operation of add with
	// the variable count and finaly execute. All the operations are asynchronous, execute work like a barrier and ensure that all the 
	// queued operations are executed
	//
	v_cl.sum(cnt);
	v_cl.execute();
	
	//
	// ### WIKI 9 ###
	//
	// Output the particle position for each processor
	//

	vd.write("output");

	//
	// ### WIKI 10 ###
	//
	// Deinitialize the library
	//
	delete_global_v_cluster();
}
