/*
 * grid_dist_id_iterator.hpp
 *
 *  Created on: Feb 4, 2015
 *      Author: Pietro Incardona
 */

#ifndef GRID_DIST_ID_ITERATOR_HPP_
#define GRID_DIST_ID_ITERATOR_HPP_

#include "grid_dist_key.hpp"
#include "Grid/grid.hpp"

/*! \brief Distributed grid iterator
 *
 * Iterator across the local element of the distributed grid
 *
 */

template<unsigned int dim, typename l_grid>
class grid_dist_iterator
{
	//! grid list counter

	size_t g_c;

	//! List of the grids we are going to iterate
	std::vector<l_grid> & gList;

	//! Actual iterator
	grid_key_dx_iterator<dim> a_it;

	public:

	/*! \brief Constructor of the distributed grid
	 *
	 * \param gk std::vector of the local grid
	 *
	 */
	grid_dist_iterator(std::vector<l_grid> & gk)
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


#endif /* GRID_DIST_ID_ITERATOR_HPP_ */
