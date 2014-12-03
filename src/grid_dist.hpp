#ifndef COM_UNIT_HPP
#define COM_UNIT_HPP

#include <vector>
#include "map_grid.hpp"
#include "VCluster.hpp"
#include "Space/SpaceBox.hpp"

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

	//! List of the grids we are going to iterate
	std::vector<grid_key_dx<dim>> & gList;

	//! Actual iterator
	grid_key_dx_iterator<dim> a_it;

	public:

	/*! \brief Constructor of the distributed grid
	 *
	 * \param gk std::vector of the local grid
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
 * \param dim Dimensionality of the grid
 * \param T type of grid
 * \param Decomposition Class that decompose the grid for example CartDecomposition
 * \param Mem Is the allocator
 * \param device type of base structure is going to store the data
 *
 */

template<unsigned int dim, typename T, typename Decomposition, typename Mem = typename memory_traits_lin< typename T::type >::type, typename device_grid=grid_cpu<dim,T> >
class grid_dist
{
	//! Local grids
	device_grid * loc_grid;

	//! Space Decomposition
	Decomposition dec;

	//! Size of the grid on each dimension
	std::vector<size_t> g_res;

	//! Communicator class

	Vcluster & v_cl;

	/*! \brief Get the grid size
	 *
	 * Get the grid size, given a domain, the resolution on it and another spaceBox
	 * it give the size on all directions of the local grid
	 *
	 * \param sp SpaceBox enclosing the local grid
	 * \param domain Space box enclosing the physical domain or part of it
	 * \param v_size grid size on this physical domain
	 *
	 * \return An std::vector representing the local grid on each dimension
	 *
	 */
	std::vector<size_t> getGridSize(SpaceBox<dim,typename Decomposition::domain_type> & sp, Box<dim,typename Decomposition::domain_type> & domain, std::vector<size_t> & v_size)
	{
		std::vector<size_t> tmp;
		for (size_t d = 0 ; d < dim ; d++)
		{
			//! Get the grid size compared to the domain space and its resolution
			typename Decomposition::domain_type dim_sz = (sp.getHigh(d) - sp.getLow(d)) / ((domain.getHigh(d) - domain.getLow(d)) / v_size[d]) + 0.5;

			// push the size of the local grid
			tmp.push_back(dim_sz);
		}
		return tmp;
	}

public:

	//! constructor
	grid_dist(Vcluster v_cl)
	:loc_grid(NULL),v_cl(v_cl),dec(*global_v_cluster)
	{
	}

	//! constructor
	grid_dist()
	:loc_grid(NULL),v_cl(*global_v_cluster),dec(*global_v_cluster)
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

		// Each processing unit take a space
		// Set the space associated to this process unit

		dec.setSpace(v_cl.getProcessUnitID());

		// Allocate the grids

		for (size_t i = 0 ; i < n_grid ; i++)
		{
			// Get the local hyper-cube

			SpaceBox<dim,typename Decomposition::domain_type> sp = dec.getLocalHyperCube(i);

			// Calculate the local grid size

			std::vector<size_t> l_res = getGridSize(sp,dec.getDomain(),g_res);

			// Set the dimensions of the local grid

			loc_grid[i].setDimensions(l_res);
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
