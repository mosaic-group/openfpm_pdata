/*! \page Vector_4_cp Vector 4 complex properties and serialization
 *
 * \subpage Vector_4_complex_prop
 * \subpage Vector_4_complex_prop_ser
 *
 */

/*!
 *
 * \page Vector_4_complex_prop Vector 4 complex properties
 *
 *
 * [TOC]
 *
 *
 * # Vector 4 complex properties # {#vector_example_cp}
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
	 * \page Vector_4_complex_prop Vector 4 complex properties
	 *
	 *
	 * ## Initialization and vector creation ##
	 *
	 * We first initialize the library and define useful constants
	 *
	 * \see \ref e0_s_init
	 *
	 * \snippet Vector/4_complex_prop/main.cpp lib init
	 *
	 * We also define a custom structure
	 *
	 * \snippet Vector/4_complex_prop/main.cpp struct A
	 *
	 * After we initialize the library we can create a vector with complex properties
	 * with the following line
	 *
	 * \snippet Vector/4_complex_prop/main.cpp vect create
	 *
	 * In this this particular case every particle carry a scalar,
	 * a vector in form of float[3], a Point, a list
	 * in form of vector of float and a list of custom structures, and a vector of vector.
	 * In general particles can have properties of arbitrary complexity.
	 *
	 * \warning For arbitrary complexity mean that we can use any openfpm data structure with and arbitrary nested complexity.
	 *          For example a openfpm::vector<aggregate<grid_cpu<openfpm::vector<aggregate<double,double[3]>>>,openfpm::vector<float>> is valid
	 * \verbatim

	       particle
	          *
	        vector
	         / \
	        /   \
	     grid    vector<float>
	      /\
	     /  \
	double  double[3]

	 * \endverbatim
	 *
	 * Our custom data-structure A is defined below. Note that this data-structure
	 * does not have pointers
	 *
	 * \snippet Vector/4_complex_prop/main.cpp struct A
	 *
	 *
	 * \warning custom data structure are allowed only if they does not have pointer.
	 *          In case they have pointer we have to define how to serialize our data-structure
	 *
	 * \see \ref vector_example_cp_ser
	 *
	 */

	//! \cond [lib init] \endcond

    // initialize the library
	openfpm_init(&argc,&argv);

	// Here we define our domain a 2D box with internals from 0 to 1.0 for x and y
	Box<2,float> domain({0.0,0.0},{1.0,1.0});

	// Here we define the boundary conditions of our problem
    size_t bc[2]={PERIODIC,PERIODIC};

	// extended boundary around the domain, and the processor domain
	Ghost<2,float> g(0.01);
	
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

	// A list of list
	constexpr int listlist = 5;

	//! \cond [lib init] \endcond

	//! \cond [struct A] \endcond

	// The custom structure
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

	//! \cond [vect create] \endcond

	vector_dist<2,float, aggregate<float,
	                               float[3],
								   Point<3,double>,
								   openfpm::vector<float>,
	                               openfpm::vector<A>,
								   openfpm::vector<openfpm::vector<float>>> >
	vd(4096,domain,bc,g);

	//! \cond [vect create] \endcond

	/*!
	 *
	 * \page Vector_4_complex_prop Vector 4 complex properties
	 *
	 *
	 * ## Assign values to properties ##
	 *
	 * Assign values to properties does not changes, from the simple case. Consider
	 * now that each particle has a list, so when we can get the property listA for particle p
	 * and resize such list with **vd.getProp<listA>(p).resize(...)**. We can add new elements at the
	 * end with **vd.getProp<listA>(p).add(...)** and get some element of this list with **vd.getProp<listA>(p).get(i)**.
	 * More in general vd.getProp<listA>(p) return a reference to the openfpm::vector contained by the particle.
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

		size_t n_cp = (float)10.0 * rand()/RAND_MAX;

		vd.getProp<listA>(p).resize(n_cp);

		for (size_t i = 0 ; i < n_cp ; i++)
		{
			vd.getProp<list>(p).add(i + 10);
			vd.getProp<list>(p).add(i + 20);
			vd.getProp<list>(p).add(i + 30);

			vd.getProp<listA>(p).get(i) = A(i+10.0,i+20.0);
		}

		vd.getProp<listlist>(p).resize(2);
		vd.getProp<listlist>(p).get(0).resize(2);
		vd.getProp<listlist>(p).get(1).resize(2);

		vd.getProp<listlist>(p).get(0).get(0) = 1.0;
		vd.getProp<listlist>(p).get(0).get(1) = 2.0;
		vd.getProp<listlist>(p).get(1).get(0) = 3.0;
		vd.getProp<listlist>(p).get(1).get(1) = 4.0;

		// next particle
		++it;
	}

	//! \cond [vect assign] \endcond

	/*!
	 *
	 * \page Vector_4_complex_prop Vector 4 complex properties
	 *
	 *
	 * ## Mapping and ghost_get ##
	 *
	 * Particles are redistributed across processors all properties are communicated but instead of
	 * using map we use **map_list** that we can use to select properties.
	 * A lot of time complex properties can be recomputed and communicate them is not a good idea.
	 * The same concept also apply for ghost_get. In general we choose which properties to communicate
	 *
	 *
	 * \see \ref e0_s_map
	 *
	 * \see \ref e1_part_ghost
	 *
	 * \snippet Vector/4_complex_prop/main.cpp vect map ghost
	 *
	 */

	//! \cond [vect map ghost] \endcond

	// Particles are redistribued across the processors but only the scalar,vector, and point properties
	// are transfert
	vd.map_list<scalar,vector,point,list,listA,listlist>();
	
	// Synchronize the ghost
	vd.ghost_get<scalar,vector,point,listA,listlist>();

	//! \cond [vect map ghost] \endcond

	/*!
	 *
	 * \page Vector_4_complex_prop Vector 4 complex properties
	 *
	 *
	 * ## Output and VTK visualization ##
	 *
	 * Vector with complex properties can be still be visualized, because unknown properties are
	 * automatically excluded
	 *
	 * \see \ref e0_s_vis_vtk
	 *
	 * \snippet Vector/4_complex_prop/main.cpp vtk
	 *
	 */

	//! \cond [vtk] \endcond

	vd.write("particles");

	//! \cond [vtk] \endcond

	/*!
	 *
	 * \page Vector_4_complex_prop Vector 4 complex properties
	 *
	 * ## Print 4 particles in the ghost area ##
	 *
	 * Here we print that the first 4 particles to show that the list of A and the list of list are filled
	 * and the ghosts contain the correct information
	 *
	 * \snippet Vector/4_complex_prop/main.cpp print ghost info
	 *
	 */

	//! \cond [print ghost info] \endcond

	size_t fg = vd.size_local();

	Vcluster & v_cl = create_vcluster();
	if (v_cl.getProcessUnitID() == 0)
	{
		for ( ; fg < vd.size_local()+4 ; fg++)
		{
			std::cout << "List of A" << std::endl;
			for (size_t i = 0 ; i < vd.getProp<listA>(fg).size() ; i++)
				std::cout << "Element: " << i << "   p1=" << vd.getProp<listA>(fg).get(i).p1 << "   p2=" << vd.getProp<listA>(fg).get(i).p2 << std::endl;

			std::cout << "List of list" << std::endl;
			for (size_t i = 0 ; i < vd.getProp<listlist>(fg).size() ; i++)
			{
				for (size_t j = 0 ; j < vd.getProp<listlist>(fg).get(i).size() ; j++)
					std::cout << "Element: " << i << "  " << j << "   " << vd.getProp<listlist>(fg).get(i).get(j) << std::endl;
			}
		}
	}

	//! \cond [print ghost info] \endcond

	/*!
	 * \page Vector_4_complex_prop Vector 4 complex properties
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
	 * \page Vector_4_complex_prop Vector 4 complex properties
	 *
	 * # Full code # {#code}
	 *
	 * \include Vector/4_complex_prop/main.cpp
	 *
	 */
}
