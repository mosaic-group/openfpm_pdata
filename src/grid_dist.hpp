#ifndef COM_UNIT_HPP
#define COM_UNIT_HPP

#include <vector>
#include "map_grid.hpp"

/*! \brief Distributed grid iterator
 *
 * Iterator across the local element of the distributed grid
 *
 *
 *
 */

template<unsigned int dim>
class grid_dist_iterator
{
	//! List if the grids on which iterate

	std::vector<grid_key_dx<dim>> gList;

	public:

	/*! \brief Constructor of the distributed grid
	 *
	 * Constructor of the distributed grid
	 *
	 * \param vector of the local grids
	 *
	 */
	grid_dist_iterator(std::vector<grid_key_dx<dim>> & gk)
	:gList(gk)
	{

	}

	~grid_dist_iterator()
	{

	}
};

/*! \brief This is a distributed grid
 *
 * Implementation of a distributed grid. A distributed grid is a grid distributed
 * across processors
 *
 * \dim Dimensionality of the grid
 * \T type of grid
 * \Decomposition Class that decompose the grid
 * \Mem Is the allocator
 * \device type of base structure is going to store the data
 *
 */

template<unsigned int dim, typename T, typename Decomposition, typename Mem = typename memory_traits_lin< typename T::type >::type, typename device_grid=device_g<dim,T>::cpu >
class grid_dist
{
	//! Local grids
	device_grid * loc_grid;

	//! Space Decomposition
	Decomposition dec;

	//! constructor
	grid_dist(Decomposition dec)
	{
		// ! Create an hyper-cube approximation
		// ! in order to work on grid_dist the decomposition
		// ! has to be a set of hyper-cube

		dec.HyperCube();

		// Get the number of local grid needed

		size_t n_grid = dec.getNHyperCube();

		// create local grids for each hyper-cube

		loc_grid = new device_grid[n_grid];

		// Allocate the grids

		for (size_t i = 0 ; i < n_grid ; i++)
		{
			loc_grid[i].setDimensions(dec.getGridDims());
		}
	}

	//! Get iterator
	getIterator()
	{

	}

	//! Destructor
	~grid_dist()
	{
		// destroy the memory

		delete [] loc_grid;
	}
};

#endif
