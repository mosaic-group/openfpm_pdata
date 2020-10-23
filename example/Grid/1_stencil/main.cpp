#include "Grid/grid_dist_id.hpp"
#include "data_type/aggregate.hpp"
#include "Decomposition/CartDecomposition.hpp"

/*!
 * \page Grid_1_stencil Stencil example
 *
 *
 * # Stencil example and ghost # {#e1_st}
 *
 * This example show how to move grid_key in order to create a Laplacian stencil,
 * be careful, the function move is convenient, but not the fastest implementation.
 * We also show how to do ghost communications
 *
 */

/*!
 * \page Grid_1_stencil Stencil example
 *
 * Define some convenient constants and types
 *
 * \snippet Grid/1_stencil/main.cpp useful constant
 *
 */

//! \cond [useful constant] \endcond

constexpr size_t x = 0;
constexpr size_t y = 1;
constexpr size_t z = 2;

constexpr size_t A = 0;
constexpr size_t B = 1;

//! \cond [useful constant] \endcond

int main(int argc, char* argv[])
{
	/*!
	 * \page Grid_1_stencil Stencil example
	 *
	 * ## Initialization ## {#e1_st_init}
	 *
	 * Initialize the library and several objects
	 *
	 * \see \ref e0_s_initialization
	 *
	 * \snippet Grid/1_stencil/main.cpp parameters
	 *
	 *
	 */

	//! \cond [parameters] \endcond

	openfpm_init(&argc,&argv);

	// domain
	Box<3,float> domain({0.0,0.0,0.0},{1.0,1.0,1.0});

	// grid sizes
	size_t sz[3] = {100,100,100};

	// ghost extension
	Ghost<3,float> g(0.03);

	//! \cond [parameters] \endcond

	/*!
	 * \page Grid_1_stencil Stencil example
	 *
	 * ## Grid create ## {#e1_st_inst}
	 *
	 * Create a distributed grid in 3D. With typedef we create an alias name for aggregate<float[3],float[3]>.
	 * In practice the type of grid_point == aggregate<float[3],float[3]>
	 *
	 * \see \ref e0_s_grid_inst
	 *
	 * \snippet Grid/1_stencil/main.cpp grid
	 *
	 */

	//! \cond [grid] \endcond

	// a convenient alias for aggregate<...>
	typedef aggregate<float,float> grid_point;

	grid_dist_id<3, float, grid_point> g_dist(sz,domain,g);

	//! \cond [grid] \endcond

	/*!
	 * \page Grid_1_stencil Stencil example
	 *
	 * ## Loop over grid points ## {#e1_s_loop_gp}
	 *
	 * Get an iterator that go through the point of the domain (No ghost)
	 *
	 * \see \ref e0_s_loop_gp
	 *
	 * \snippet Grid/1_stencil/main.cpp iterator
	 * \snippet Grid/1_stencil/main.cpp iterator2
	 *
	 */

	//! \cond [iterator] \endcond

	auto dom = g_dist.getDomainIterator();

	while (dom.isNext())
	{

		//! \cond [iterator] \endcond

		/*!
		 * \page Grid_1_stencil Stencil example
		 *
		 * Inside the cycle we get the local grid key
		 *
		 * \see \ref e0_s_grid_coord
		 *
		 * \snippet Grid/1_stencil/main.cpp local key
		 *
		 */

		//! \cond [local key] \endcond

		auto key = dom.get();

		//! \cond [local key] \endcond

		/*!
		 * \page Grid_1_stencil Stencil example
		 *
		 * We convert the local grid position, into global position, key_g contain 3 integers that identify the position
		 * of the grid point in global coordinates
		 *
		 * \see \ref e0_s_grid_coord
		 *
		 * \snippet Grid/1_stencil/main.cpp global key
		 *
		 */

		//! \cond [global key] \endcond

		auto key_g = g_dist.getGKey(key);

		//! \cond [global key] \endcond

		/*!
		 * \page Grid_1_stencil Stencil example
		 *
		 * we write on the grid point of position (i,j,k) the value i*i + j*j + k*k on the property A.
		 * Mathematically is equivalent to the function
		 *
		 * \f$ f(x,y,z) = x^2 + y^2 + z^2 \f$
		 *
		 * \snippet Grid/1_stencil/main.cpp function
		 *
		 */

		//! \cond [function] \endcond

		g_dist.template get<A>(key) = key_g.get(0)*key_g.get(0) + key_g.get(1)*key_g.get(1) + key_g.get(2)*key_g.get(2);

		//! \cond [function] \endcond

		//! \cond [iterator2] \endcond

		++dom;
	}

	//! \cond [iterator2] \endcond

	/*!
	 * \page Grid_1_stencil Stencil example
	 *
	 * ## Ghost ## {#e1_s_ghost}
	 *
	 * Each sub-domain has an extended part, that is materially contained into another processor.
	 * In general is not synchronized
	 * ghost_get<A> synchronize the property A in the ghost part
	 *
	 * \snippet Grid/1_stencil/main.cpp ghost
	 *
	 */

	//! \cond [ghost] \endcond

	g_dist.template ghost_get<A>();
	
	//! \cond [ghost] \endcond

	/*!
	 * \page Grid_1_stencil Stencil example
	 *
	 * Get again another iterator, iterate across all the domain points, calculating a Laplace stencil. Write the
	 * result on B
	 *
	 * \snippet Grid/1_stencil/main.cpp laplacian
	 *
	 */

	//! \cond [laplacian] \endcond

	auto dom2 = g_dist.getDomainIterator();
	
	while (dom2.isNext())
	{
		auto key = dom2.get();

		// Laplace stencil
		g_dist.template get<B>(key) = g_dist.template get<A>(key.move(x,1)) + g_dist.template get<A>(key.move(x,-1)) +
		                                 g_dist.template get<A>(key.move(y,1)) + g_dist.template get<A>(key.move(y,-1)) +
										 g_dist.template get<A>(key.move(z,1)) + g_dist.template get<A>(key.move(z,-1)) -
										 6*g_dist.template get<A>(key);

		++dom2;
	}

	//! \cond [laplacian] \endcond

	/*!
	 * \page Grid_1_stencil Stencil example
	 *
	 *
	 * Finally we want a nice output to visualize the information stored by the distributed grid
	 *
	 * \see \ref e0_s_VTK_vis
	 *
	 * \snippet Grid/1_stencil/main.cpp output
	 *
	 */

	//! \cond [output] \endcond

	g_dist.write("output");

	//! \cond [output] \endcond

	/*!
	 * \page Grid_1_stencil Stencil example
	 *
	 * Deinitialize the library
	 *
	 * \snippet Grid/1_stencil/main.cpp finalize
	 *
	 */

	//! \cond [finalize] \endcond

	openfpm_finalize();

	//! \cond [finalize] \endcond

	/*!
	 * \page Grid_1_stencil Stencil example
	 *
	 * # Full code # {#code}
	 *
	 * \include Grid/1_stencil/main.cpp
	 *
	 */
}
