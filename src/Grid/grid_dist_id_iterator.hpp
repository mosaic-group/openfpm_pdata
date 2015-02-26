/*
 * grid_dist_id_iterator_sub.hpp
 *
 *  Created on: Feb 4, 2015
 *      Author: Pietro Incardona
 */

#ifndef GRID_DIST_ID_ITERATOR_HPP_
#define GRID_DIST_ID_ITERATOR_HPP_

#include "grid_dist_key.hpp"
#include "VCluster.hpp"

/*! \brief Distributed grid iterator
 *
 * Iterator across the local element of the distributed grid
 *
 */

template<unsigned int dim, typename device_grid>
class grid_dist_iterator
{
	//! grid list counter
	size_t g_c;

	//! List of the grids we are going to iterate
	Vcluster_object_array<device_grid> & gList;

	//! Actual iterator
	grid_key_dx_iterator_sub<dim> a_it;

	//! margin of the grid iterator
	size_t m;

	public:

	/*! \brief Constructor of the distributed grid
	 *
	 * \param gk std::vector of the local grid
	 *
	 */
	grid_dist_iterator(Vcluster_object_array<device_grid> & gk, size_t m)
	:g_c(0),gList(gk),m(m)
	{
		// Initialize the current iterator
		// with the first grid

		a_it.reinitialize(gList[0].getDomainIterator());
	}

	// Destructor
	~grid_dist_iterator()
	{
	}

	/*! \brief Get the next element
	 *
	 * \return the next grid_key
	 *
	 */

	grid_dist_iterator<dim,device_grid> operator++()
	{
		++a_it;

		// check if a_it is at the end

		if (a_it.isNext() == true)
			return *this;
		else
		{
			// switch to the new grid

			g_c++;

			// get the next grid iterator

			if (g_c < gList.size())
				a_it = gList[g_c].getDomainIterator();
		}

		return *this;
	}

	/*! \brief Check if there is the next element
	 *
	 * \return true if there is the next, false otherwise
	 *
	 */

	bool isNext()
	{
		// If there are no other grid stop

		if (g_c >= gList.size())
			return false;

		return true;
	}

	/*! \brief Get the actual key
	 *
	 * \return the actual key
	 *
	 */
	grid_dist_key_dx<dim> get()
	{
		return grid_dist_key_dx<dim>(g_c,a_it.get());
	}
};


#endif /* GRID_DIST_ID_ITERATOR_SUB_HPP_ */
