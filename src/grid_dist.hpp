#ifndef COM_UNIT_HPP
#define COM_UNIT_HPP

#include <vector>
#include "map_grid.hpp"

/*! \brief Grid key for a distributed grid
 *
 * Grid key for a distributed grid
 *
 */

template<unsigned int dim>
class grid_dist_key_dx
{
	//! grid list counter

	size_t g_c;

	//! Local grid iterator

	grid_key_dx_iterator<dim> a_it;
};

/*! \brief Distributed grid iterator
 *
 * Iterator across the local element of the distributed grid
 *
 */

template<unsigned int dim>
class grid_dist_iterator
{
	//! grid list counter

	size_t g_c;

	//! List if the grids on which iterate

	std::vector<grid_key_dx<dim>> & gList;

	//! Actual iterator

	grid_key_dx_iterator<dim> a_it;

	public:

	/*! \brief Constructor of the distributed grid
	 *
	 * Constructor of the distributed grid
	 *
	 * \param vector of the local grids
	 *
	 */
	grid_dist_iterator(std::vector<grid_key_dx<dim>> & gk)
	:g_c(0),gList(gk)
	{
		// Initialize with the current iterator
		// with the first grid

		a_it = gList[0].getIterator();
	}

	// Destructor
	~grid_dist_iterator()
	{}

	/*! \brief Get the next element
	 *
	 * Get the next element
	 *
	 * \return the next grid_key
	 *
	 */

	grid_key_dx_iterator<dim> operator++()
	{
		a_it++;

		// check if a_it is at the end

		if (a_it.isEnd() == false)
			return *this;
		else
		{
			// switch to the new grid

			g_c++;

			// get the next grid iterator

			a_it = a_it = gList[g_c].getIterator();

			// increment to a valid point

			a_it++;
		}

		return *this;
	}

	/*! \brief Check if there is the next element
	 *
	 * Check if there is the next element
	 *
	 * \return true if there is the next, false otherwise
	 *
	 */

	bool isEnd()
	{
		// If there are no other grid stop

		if (g_c >= gList.size())
			return true;
	}

	/*! \brief Get the actual key
	 *
	 * Get the actual key
	 *
	 * \return the actual key
	 *
	 */
	grid_dist_key_dx<dim> get()
	{
		return a_it;
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

template<unsigned int dim, typename T, typename Decomposition, typename Mem = typename memory_traits_lin< typename T::type >::type, typename device_grid=device_g<dim,T> >
class grid_dist
{
	//! Local grids
	device_grid * loc_grid;

	//! Space Decomposition
	Decomposition dec;

public:

	//! constructor
	grid_dist()
	:loc_grid(NULL)
	{
	}

	/*! \brief get the object that store the decomposition information
	 *
	 * get the object that store the decomposition information
	 *
	 * \return the decomposition object
	 *
	 */

	Decomposition & getDecomposition()
	{
		return dec;
	}

	/*! \brief Decompose the domain
	 *
	 * Decompose the domain
	 *
	 */

	void Decompose()
	{
		// ! Create an hyper-cube approximation.
		// ! In order to work on grid_dist the decomposition
		// ! has to be a set of hyper-cube

		dec.hyperCube();

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
/*	getIterator()
	{

	}*/

	//! Destructor
	~grid_dist()
	{
		// destroy the memory

		delete [] loc_grid;
	}
};

#endif
