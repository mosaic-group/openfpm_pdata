#include "Grid/grid_dist_id.hpp"
#include "data_type/aggregate.hpp"

/*!
 * \page Grid_2_solve_eq Grid 2 solve eq
 *
 * [TOC]
 *
 * # Simple example of grid usage to solve an equation # {#e2_solve_eq}
 *
 * This example show the usage of grid to solve the following equation
 *
 * \f$\frac{\partial^2 u}{\partial^2 x} + \frac{\partial^2 u}{\partial^2 y} = 1\f$
 * 
 * \f$u(x,y) = 0 \f$
 *
 * at the boundary
 * 
 *
 * ## Field initialization ## {#e2_se_finit}
 *
 */

//! \cond [field init] \endcond

void init(grid_dist_id<2,double,aggregate<double> > & g_dist, const size_t (& sz)[2])
{
	/*!
	 * \page Grid_2_solve_eq Grid 2 solve eq
	 *
	 * In order to initialize the field U, first we get an iterator that cover
	 *  domain + Ghost to iterate all the grid points.
	 * Inside the cycle we initialize domain and border (Ghost part to 0.0)
	 *
	 * \see \ref e0_s_loop_gp
	 *
	 */

	//! \cond [iterator] \endcond

	// Get the iterator
	auto it = g_dist.getDomainGhostIterator();

	// For each point in the grid
	while (it.isNext())
	{
		//! \cond [iterator] \endcond

		/*!
		 * \page Grid_2_solve_eq Grid 2 solve eq
		 *
		 * Get the local grid key
		 *
		 * \see \ref e0_s_grid_coord
		 *
		 * \snippet Grid/2_solve_eq/main.cpp local key
		 *
		 */

		//! \cond [local key] \endcond

		auto key = it.get();

		//! \cond [local key] \endcond

		/*!
		 * \page Grid_2_solve_eq Grid 2 solve eq
		 *
		 *
		 * Here we convert the local grid position, into global position.
		 * key_g contain 3 integers that identify the position of the grid point
		 * in global coordinates
		 *
		 * \see \ref e0_s_grid_coord
		 *
		 * \snippet Grid/2_solve_eq/main.cpp global key
		 *
		 */

		//! \cond [global key] \endcond

		auto key_g = g_dist.getGKey(key);

		//! \cond [global key] \endcond

		/*!
		 * \page Grid_2_solve_eq Grid 2 solve eq
		 *
		 * Initialize to 0, domain + boundary
		 *
		 * \snippet Grid/2_solve_eq/main.cpp init field zero
		 *
		 * The full function look like this
		 *
		 * \snippet Grid/2_solve_eq/main.cpp field init
		 *
		 */

		//! \cond [init field zero] \endcond

		if (key_g.get(0) == 0 || key_g.get(0) == sz[0] ||
			key_g.get(1) == 0 || key_g.get(1) == sz[1])
		{
			// Boundary part
			g_dist.template get<0>(key) = 0.0;
		}
		else
		{
			// Internal part
			g_dist.template get<0>(key) = 0.0;
		}

		//! \cond [init field zero] \endcond

		//! \cond [iterator2] \endcond

		++it;

	}

	//! \cond [iterator2] \endcond
}

//! \cond [field init] \endcond

constexpr int x = 0;
constexpr int y = 1;

int main(int argc, char* argv[])
{
	/*!
	 * \page Grid_2_solve_eq Grid 2 solve eq
	 *
	 * ## Initialization ##
	 *
	 * Initialize the library
	 *
	 * \see \ref e0_s_initialization
	 *
	 * \snippet Grid/2_solve_eq/main.cpp ofp_init
	 *
	 */

	//! \cond [ofp_init] \endcond

	openfpm_init(&argc,&argv);
	
	//! \cond [ofp_init] \endcond

	/*!
	 * \page Grid_2_solve_eq Grid 2 solve eq
	 *
	 * ## Grid instantiation and initialization ##
	 *
	 * Create
	 * * A 2D box that define the domain
	 * * an array of 2 unsigned integer that will define the size of the grid on each dimension
	 * * A Ghost object that will define the extension of the ghost part for each sub-domain in grid point unit
	 *
	 * \snippet Grid/2_solve_eq/main.cpp ofp_par
	 *
	 */

	//! \cond [ofp_par] \endcond

	Box<2,double> domain({-1.0,-1.0},{1.0,1.0});
	size_t sz[2] = {64,64};
	
	periodicity<2> bc = {NON_PERIODIC,NON_PERIODIC};
	
	// Ghost in grid unit
	Ghost<2,long int> g(1);
	
	//! \cond [ofp_par] \endcond

	/*!
	 * \page Grid_2_solve_eq Grid 2 solve eq
	 *
	 * Create a distributed grid in 2D (1° template parameter) space in with double precision (2° template parameter)
	 * each grid point contain a scalar (double),
	 *
	 * Constructor parameters:
	 *
	 * * sz: size of the grid on each dimension
	 * * domain: where the grid is defined
	 * * g: ghost extension
	 * * bc: boundary conditions
	 *
	 * \see \ref e0_s_grid_inst
	 *
	 * \snippet Grid/2_solve_eq/main.cpp grid inst
	 *
	 */

	//! \cond [grid inst] \endcond

	grid_dist_id<2, double, aggregate<double>> g_dist(sz,domain,g,bc);
	
	// spacing between two points
	double spacing[2] = {g_dist.spacing(0),g_dist.spacing(1)};

	//! \cond [grid inst] \endcond

	/*!
	 * \page Grid_2_solve_eq Grid 2 solve eq
	 *
	 * Initialize U and fill the boundary conditions
	 *
	 * \see \ref e2_se_field_init
	 *
	 * \snippet Grid/2_solve_eq/main.cpp grid init
	 *
	 */

	//! \cond [grid init] \endcond

	init(g_dist,sz);

	//! \cond [grid init] \endcond

	/*!
	 * \page Grid_2_solve_eq Grid 2 solve eq
	 *
	 * ## %Ghost synchronization ##
	 *
	 * Before the computation read the ghost point we have to guarantee that they have
	 *  updated values.
	 *
	 *  \snippet Grid/2_solve_eq/main.cpp ghost sync
	 *
	 */

	//! \cond [ghost sync] \endcond

	// sync the ghost property 0
	g_dist.template ghost_get<0>();

	//! \cond [ghost sync] \endcond

	/*!
	 * \page Grid_2_solve_eq Grid 2 solve eq
	 *
	 * ## Red-Black alghorithm ##
	 *
	 * Do 10000 iteration of Red-Black Gauss-Siedel iterations
	 *
	 *
	 */

	//! \cond [gs_alg] \endcond

	// flag that indicate if we are processing red or black
	// we start from red
	bool red_black = true;

	// 10000 iterations
	for (size_t i = 0 ; i < 10000 ; i++)
	{
		/*!
		 * \page Grid_2_solve_eq Grid 2 solve eq
		 *
		 * Get an iterator that go through the points of the grid (No ghost)
		 * To compute one iteration.
		 *
		 * \see \ref e0_s_loop_gp
		 * \see \ref e0_s_grid_coord
		 *
		 * \snippet Grid/2_solve_eq/main.cpp gs_it
		 *
		 */

		//! \cond [gs_it] \endcond

		auto dom = g_dist.getDomainIterator();
		
		// Iterate over all the points
		while (dom.isNext())
		{

			// Get the local grid key
			auto key = dom.get();

			// Here we convert the local grid position, into global position, key_g contain 3 integers that identify the position
			// of the grid point in global coordinates
			auto key_g = g_dist.getGKey(key);


			//
			// If we are processing red and is odd jump to the next point
			// If we are processing black and is even jump to the next point
			//
			if (red_black == false && (key_g.get(0) + key_g.get(1)) % 2 == 0)
			{++dom; continue;}
			else if (red_black == true && (key_g.get(0) + key_g.get(1)) % 2 == 1)
			{++dom; continue;}

			//
			// Update the grid values
			//
			// P.S. The keyword template is removed, it is possible only if we are in a function
			// without template parameters (if you are unsure use the keyword template)
			//
			g_dist.get<0>(key) = (g_dist.get<0>(key.move(x,1)) + g_dist.template get<0>(key.move(x,-1)) +
                    			 g_dist.get<0>(key.move(y,1)) + g_dist.template get<0>(key.move(y,-1)) +
								 - 1.0)/4.0;

			//
			// next point (red/black)
			++dom;
		}

		//! \cond [gs_it] \endcond

		/*!
		 * \page Grid_2_solve_eq Grid 2 solve eq
		 *
		 *
		 * Once an iteration is done we have to synchronize the ghosts
		 * to start a new iteration. Consider that we calculated the red points
		 * in the next iteration the red points in the ghost are used in reading.
		 * This mean that the ghost must have the updated values
		 *
		 * \snippet Grid/2_solve_eq/main.cpp ghost sync2
		 *
		 *
		 */

		//! \cond [ghost sync2] \endcond

		g_dist.template ghost_get<0>();

		// switch from red to black or black to red
		red_black = !red_black;

		//! \cond [ghost sync2] \endcond
	}
	
	/*!
	 * \page Grid_2_solve_eq Grid 2 solve eq
	 *
	 * The full Algorithm look like this
	 *
	 * \snippet Grid/2_solve_eq/main.cpp gs_alg
	 */

	//! \cond [gs_alg] \endcond

	/*!
	 * \page Grid_2_solve_eq Grid 2 solve eq
	 *
	 * ## Solution statistic ##
	 *
	 * Once we got the solution we want to check if it really satisfy the equation,
	 * and calculate the error (norm infinity)
	 *
	 * \snippet Grid/2_solve_eq/main.cpp sol stat
	 *
	 */

	//! \cond [sol stat] \endcond

	// It contain the error
	double error = 0.0;

	// Get the iterator
	auto dom = g_dist.getDomainIterator();

	// Iterate over all the points
	while (dom.isNext())
	{
		// same the the grid point and the global grid point
		auto key = dom.get();

		// Calculate the error on each point
		// The error is how much the solution does not respect the equation
		double error_tmp = abs((g_dist.get<0>(key.move(x,1)) + g_dist.get<0>(key.move(x,-1)) +
                		   g_dist.get<0>(key.move(y,1)) + g_dist.get<0>(key.move(y,-1)) +
						   - 4.0*g_dist.get<0>(key)) - 1.0);

		// In the norm infinity the maximum error across all the point is
		// important
		if (error_tmp > error)
			error = error_tmp;

		// next point
		++dom;

	}

	// Get the maximum across processor to calculate the norm infinity of the error
	// Norm infinity of the error is the maximum error across all the grid points
	Vcluster & v_cl = create_vcluster();
	v_cl.max(error);
	v_cl.execute();

	// Only master print the norm infinity
	if (v_cl.getProcessUnitID() == 0)
			std::cout << "Error norm infinity: " << error << std::endl;

	//! \cond [sol stat] \endcond

	/*!
	 * \page Grid_2_solve_eq Grid 2 solve eq
	 *
	 * ## VTK Write and visualization ##
	 *
	 * Finally we want a nice output to visualize the information stored by the distributed grid.
	 * The function write by default produce VTK files. One for each processor that can be visualized
	 * with the programs like paraview
	 *
	 * \see \ref e0_s_VTK_vis
	 *
	 * \snippet Grid/2_solve_eq/main.cpp write
	 *
	 */

	//! \cond [write] \endcond

	g_dist.write("output");

	//! \cond [write] \endcond
	

	/*!
	 * \page Grid_2_solve_eq Grid 2 solve eq
	 *
	 * ## Finalize ##
	 *
	 *  At the very end of the program we have always to de-initialize the library
	 *
	 * \snippet Grid/2_solve_eq/main.cpp finalize
	 *
	 */

	//! \cond [finalize] \endcond

	openfpm_finalize();

	//! \cond [finalize] \endcond

	/*!
	 * \page Grid_2_solve_eq Grid 2 solve eq
	 *
	 * # Full code # {#code}
	 *
	 * \include Grid/2_solve_eq/main.cpp
	 *
	 */
}
