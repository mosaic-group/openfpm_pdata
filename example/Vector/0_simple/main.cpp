#include "Vector/vector_dist.hpp"

/*
 * ### WIKI 1 ###
 *
 * ## Simple example
 * 
 * This example show several basic functionalities of the distributed vector, A distributed vector is nothing else than
 * a set of particles in an N dimensional space
 * 
 * ### WIKI END ###
 * 
 */

int main(int argc, char* argv[])
{
	//
	// ### WIKI 2 ###
	//
	// Here we Initialize the library, than we create a uniform random generator between 0 and 1 to to generate particles
	// randomly in the domain, we create a Box that define our domain, boundary conditions, and ghost
	//
	
        // initialize the library
	openfpm_init(&argc,&argv);

	// Here we define our domain a 2D box with internals from 0 to 1.0 for x and y
	Box<2,float> domain({0.0,0.0},{1.0,1.0});

	// Here we define the boundary conditions of our problem
    size_t bc[2]={PERIODIC,PERIODIC};

	// extended boundary around the domain, and the processor domain
	Ghost<2,float> g(0.01);
	
	//
	// ### WIKI 3 ###
	//
	// Here we are creating a distributed vector defined by the following parameters
	//
	// * 2 is the Dimensionality of the space where the objects live
	// * float is the type used for the spatial coordinate of the particles
	// * float,float[3],float[3][3] is the information stored by each particle a scalar float, a vector float[3] and a tensor of rank 2 float[3][3]
	//   the list of properties must be put into an aggregate data astructure aggregate<prop1,prop2,prop3, ... >
	// 
	// vd is the instantiation of the object
	//
	// The Constructor instead require:
	//
	// * Number of particles 4096 in this case
	// * Domain where is defined this structure
	// * bc boundary conditions
	// * g Ghost
	//
	// The following construct a vector where each processor has 4096 / N_proc (N_proc = number of processor)
	// objects with an undefined position in space. This non-space decomposition is also called data-driven
	// decomposition
	//
	vector_dist<2,float, aggregate<float,float[3],float[3][3]> > vd(4096,domain,bc,g);

	// the scalar is the element at position 0 in the aggregate
	const int scalar = 0;

	// the vector is the element at position 1 in the aggregate
	const int vector = 1;

	// the tensor is the element at position 2 in the aggregate
	const int tensor = 2;

	//
	// ### WIKI 5 ###
	//
	// Get an iterator that go through the 4096 particles, in an undefined position state and define its position
	//
	auto it = vd.getDomainIterator();

	while (it.isNext())
	{
		auto key = it.get();

		// we define x, assign a random position between 0.0 and 1.0
		vd.getPos(key)[0] = rand() / RAND_MAX;

		// we define y, assign a random position between 0.0 and 1.0
		vd.getPos(key)[1] = rand() / RAND_MAX;

		// next particle
		++it;
	}

	//
	// ### WIKI 6 ###
	//
	// Once we define the position, we distribute them according to the default space decomposition
	// The default decomposition is created even before assigning the position to the object. It determine
	// which part of space each processor manage
	//
	vd.map();

	//
	// ### WIKI 7 ###
	//
	// We get the object that store the decomposition, than we iterate again across all the objects, we count them
	// and we confirm that all the particles are local
	//
	//Counter we use it later
	size_t cnt = 0;

	// Get the space decomposition
	auto & ct = vd.getDecomposition();

	// Get a particle iterator
	it = vd.getDomainIterator();

	// For each particle ...
	while (it.isNext())
	{
		// ... p
		auto p = it.get();

		// we set the properties of the particle p
		
                // the scalar property
		vd.template getProp<scalar>(p) = 1.0;

		vd.template getProp<vector>(p)[0] = 1.0;
		vd.template getProp<vector>(p)[1] = 1.0;
		vd.template getProp<vector>(p)[2] = 1.0;

		vd.template getProp<tensor>(p)[0][0] = 1.0;
		vd.template getProp<tensor>(p)[0][1] = 1.0;
		vd.template getProp<tensor>(p)[0][2] = 1.0;
		vd.template getProp<tensor>(p)[1][0] = 1.0;
		vd.template getProp<tensor>(p)[1][1] = 1.0;
		vd.template getProp<tensor>(p)[1][2] = 1.0;
		vd.template getProp<tensor>(p)[2][0] = 1.0;
		vd.template getProp<tensor>(p)[2][1] = 1.0;
		vd.template getProp<tensor>(p)[2][2] = 1.0;

		// increment the counter
		cnt++;

		// next particle
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
	
	auto & v_cl = create_vcluster();
	v_cl.sum(cnt);
	v_cl.execute();
	
	//
	// ### WIKI 9 ###
	//
	// Output the particle position for each processor
	//

	vd.write("output",VTK_WRITER);

	//
	// ### WIKI 10 ###
	//
	// Deinitialize the library
	//
	openfpm_finalize();
}
