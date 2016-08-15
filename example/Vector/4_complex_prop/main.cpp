/*!
 *
 * \page Vector_4_complex_prop Vector 4 complex property
 *
 *
 * [TOC]
 *
 *
 * # Vector 4 complex property # {#vector_example_cp}
 *
 *
 * This example show how we can use complex properties in a vector
 *
 */

#include "Vector/vector_dist.hpp"

int main(int argc, char* argv[])
{
	/*!
	 *
	 * \page Vector_4_complex_prop Vector 4 complex property
	 *
	 *
	 * ## Initialization and vector creation ##
	 *
	 * After we initialize the library we can create a vector with complex properties
	 * with the following line
	 *
	 * \snippet Vector/4_complex_prop/main.cpp vect create
	 *
	 * In This this particular case every particle carry a scalar,
	 * a vector in form of float[3], a Point, a list
	 * in form of vector of float and a list of custom structures
	 *
	 * \snippet Vector/4_complex_prop/main.cpp struct A
	 *
	 *
	 */

    // initialize the library
	openfpm_init(&argc,&argv);

	// Here we define our domain a 2D box with internals from 0 to 1.0 for x and y
	Box<2,float> domain({0.0,0.0},{1.0,1.0});

	// Here we define the boundary conditions of our problem
    size_t bc[2]={PERIODIC,PERIODIC};

	// extended boundary around the domain, and the processor domain
	Ghost<2,float> g(0.01);
	
	//! \cond [struct A] \endcond

	// The a custom structure
	struct A
	{
		float p1;
		int p2;

		A() {};

		A(float p1, int p2)
		:p1(p1),p2(p2)
		{}
	};

	//! \cond [struct A] \endcond

	// the scalar is the element at position 0 in the aggregate
	constexpr int scalar = 0;

	// the vector is the element at position 1 in the aggregate
	constexpr int vector = 1;

	// the tensor is the element at position 2 in the aggregate
	constexpr int point = 2;

	// A list1
	constexpr int list = 3;

	// A listA
	constexpr int listA = 4;

	//! \cond [vect create] \endcond

	vector_dist<2,float, aggregate<float,
	                               float[3],
								   Point<3,double>,
								   openfpm::vector<float>,
	                               openfpm::vector<A>>>
	vd(4096,domain,bc,g);

	//! \cond [vect create] \endcond

	/*!
	 *
	 * \page Vector_4_complex_prop Vector 4 complex property
	 *
	 *
	 * ## Assign values to properties ##
	 *
	 * Assign values to properties does not changes, from the simple case. Consider
	 * now that each particle has a list, so when we can get the property listA for particle p
	 * and resize such list with **vd.getProp<listA>(p).resize(...)**. We can add new elements at the
	 * end with **vd.getProp<listA>(p).add(...)** and get some element with **vd.getProp<listA>(p).get(i)**.
	 * More in general from vd.getProp<listA>(p) we can use the full openfpm::vector interface.
	 *
	 * \snippet Vector/4_complex_prop/main.cpp vect assign
	 *
	 */

	//! \cond [vect assign] \endcond

	auto it = vd.getDomainIterator();

	while (it.isNext())
	{
		auto p = it.get();

		// we define x, assign a random position between 0.0 and 1.0
		vd.getPos(p)[0] = (float)rand() / RAND_MAX;

		// we define y, assign a random position between 0.0 and 1.0
		vd.getPos(p)[1] = (float)rand() / RAND_MAX;


		vd.getProp<scalar>(p) = 1.0;

		vd.getProp<vector>(p)[0] = 1.0;
		vd.getProp<vector>(p)[1] = 1.0;
		vd.getProp<vector>(p)[2] = 1.0;

		vd.getProp<point>(p).get(0) = 1.0;
		vd.getProp<point>(p).get(1) = 1.0;
		vd.getProp<point>(p).get(2) = 1.0;

		size_t n_cp = (float)10 * rand()/RAND_MAX;

		vd.getProp<listA>(p).resize(n_cp);

		for (size_t i = 0 ; i < n_cp ; i++)
		{
			vd.getProp<list>(p).add(i + 10);
			vd.getProp<list>(p).add(i + 20);
			vd.getProp<list>(p).add(i + 30);

			vd.getProp<listA>(p).get(i) = A(i+10.0,i+20.0);
			vd.getProp<listA>(p).get(i) = A(i+30.0,i+40.0);
			vd.getProp<listA>(p).get(i) = A(i+50.0,i+60.0);
		}

		// next particle
		++it;
	}

	//! \cond [vect assign] \endcond

	/*!
	 *
	 * \page Vector_4_complex_prop Vector 4 complex property
	 *
	 *
	 * ## Mapping and ghost_get ##
	 *
	 * Particles are redistributed across processors but only the scalar,vector and the point
	 * are communicated (properties 0,1,2). A lot of time complex properties can be recomputed and
	 * communicate them is not a good idea. The same concept also apply for ghost_get
	 *
	 * \note OpenFPM <= 0.5.0 cannot communicate complex properties like a vector or other structure
	 *                        that are not POD object
	 *
	 * \note OpenFPM > 0.5.0 Does not have such limitation
	 *
	 *
	 * \snippet Vector/4_complex_prop/main.cpp vect map ghost
	 *
	 */

	//! \cond [vect map ghost] \endcond

	// Particles are redistribued across the processors but only the scalar,vector, and point properties
	// are transfert
	vd.map_list<scalar,vector,point>();
	
	// Synchronize the ghost
	vd.ghost_get<scalar,vector,point>();

	//! \cond [vect map ghost] \endcond

	/*!
	 *
	 * \page Vector_4_complex_prop Vector 4 complex property
	 *
	 *
	 * ## Output and VTK visualization ##
	 *
	 * Vector with complex properties can be still be visualized, because unknown properties are
	 * automatically excluded
	 *
	 * \snippet Vector/4_complex_prop/main.cpp vtk
	 *
	 */

	//! \cond [vtk] \endcond

	vd.write("particles");

	//! \cond [vtk] \endcond

	/*!
	 * \page Vector_4_complex_prop Vector 4 complex property
	 *
	 * ## Finalize ## {#finalize}
	 *
	 *  At the very end of the program we have always to de-initialize the library
	 *
	 * \snippet Vector/4_complex_prop/main.cpp finalize
	 *
	 */

	//! \cond [finalize] \endcond

	openfpm_finalize();

	//! \cond [finalize] \endcond

	/*!
	 * \page Vector_4_complex_prop Vector 4 complex property
	 *
	 * # Full code # {#code}
	 *
	 * \include Vector/4_complex_prop/main.cpp
	 *
	 */
}
